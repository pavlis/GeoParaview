#include <math.h>
#include <vector>
#include <fstream>
#include <algorithm>
#include <ctime>
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/math/distributions/normal.hpp>
using boost::math::normal;
#include "filter++.h"
#include "slowness.h"
#include "PPWDeconProcessor.h"
#include "SimpleWavelets.h"
#include "SeismicPlot.h"
#include "TraceEditPlot.h"
#include "TimeWindowPicker.h"
//#include "MatlabPlotter.h"
using namespace std;
using namespace SEISPP;
using namespace boost;
/* Copied form http://en.wikibooks.org/wiki/C++_Programming/Libraries/Boost */
double SampleNormal(double mean, double sigma)
{
    /* Generate a random number from the current time */
    static mt19937 rng(static_cast<unsigned> (std::time(0)));
    /* This defines a normal distribution to generate numbers */
    normal_distribution<double> norm_dist(mean,sigma);
    // bind random number generator to distribution to forma function object
    variate_generator<mt19937&, normal_distribution<double> >  
        normal_sampler(rng, norm_dist);
    return(normal_sampler());
}
vector<SlownessVector> gridtovector(RectangularSlownessGrid& ugrid)
{
	vector<SlownessVector> svvec;
	int i,j;
	svvec.reserve((ugrid.nux)*(ugrid.nuy));
	for(i=0;i<ugrid.nux;++i)
		for(j=0;j<ugrid.nuy;++j)
		{
			svvec.push_back(ugrid.slow(i,j));
		}
	return svvec;
}
/* Thgis procedure scans a time 
   series object for maximum value and resets t0 to that value.  
   It then does something which a library routine should not do and that is 
   to rescale the data by the peak amplitude.  Something the decon algorithm here
   needs, but which is not in line with a library function which would do only one thing.
   It returns a copy of the input with t0 changed to that value.

   This routine could be made a library routine but this behaviour should be changed to
   alter the data in place and not do the amplitude change.
   */
TimeSeries ShiftT02MaxAmp(TimeSeries& d)
{
    vector<double>::iterator smax;
    smax=max_element(d.s.begin(),d.s.end());
    int imax=distance(d.s.begin(),smax);
    double dt0=d.time(imax);
    TimeSeries result(d);
    //result.t0-=dt0;
    double scale=1.0/(*smax);
    dscal(result.ns,scale,&(result.s[0]),1);
    //DEBUG
    /*
    cout << "Wavelet max amp at sample number "<<imax<<" = time shift of "<< dt0
        << " New t0="<<result.t0<<" scale  by factor="<<scale<<endl;
        */
    return(result);
}
void BuildFakeData(ThreeComponentEnsemble& d, RectangularSlownessGrid& g,
        Pf *pf,TimeSeries& wavelet,double sigma)
{
    /* get local values of this to avoid repeated confusing get calls */
    int nsta=d.member.size();
    int nsamp=d.member[0].ns;
    double dt=d.member[0].dt;
    int i,j;
    int lag;
    double ux,uy;
    double amp[3];
    Tbl *t;
    t=pfget_tbl(pf,"PlaneWaves");
    if(t==NULL)
    {
        cerr << "Missing PlaveWaves Tbl in parameter file.  Exiting"<<endl;
        exit(-1);
    }
    for(i=0;i<maxtbl(t);++i)
    {
        double dnorth,deast,lagtime;
        int stalag;
        char *line=(char *)gettbl(t,i);
        sscanf(line,"%lf%lf%lf%lf%lf%d",&ux,&uy,amp,amp+1,amp+2,&lag);
        for(j=0;j<nsta;++j)
        {
            /* assume these are set */
            dnorth=d.member[j].get_double("dnorth");
            deast=d.member[j].get_double("deast");
            lagtime=ux*deast+uy*dnorth;
            stalag=SEISPP::nint(lagtime/dt);
            stalag+=lag;
            if(stalag<0 || stalag>=nsamp)
            {
                cerr << "BuildFakeData:  illegal computed = "<<lag<<endl;
                cerr << "Input parameter error.  Encountered trying slowness="
                    <<ux<<", "<<uy<<endl
                    <<"dnorth="<<dnorth<<" deast="<<deast<<" using base lag="
                    <<lag<<endl
                    <<"Computed lag = "<<stalag
                    <<" which exceed range of 0 to "<<nsamp<<endl;
                exit(-1);
            }
            for(int k=0;k<3;++k)d.member[j].u(k,stalag)+=amp[k];
        }
    }
    /* Add gausian noise */
    for(j=0;j<nsta;++j)
    {
        int ns=d.member[j].ns;
        ns=3*ns;
        double *ptr=d.member[j].u.get_address(0,0);
        for(i=0;i<ns;++i)
        {
            double n=SampleNormal(0.0,sigma);
            *ptr=(*ptr)+n;
            ++ptr;
        }
    }

    /* now convolve with the desired wavelet */
    for(j=0;j<nsta;++j){
        d.member[j]=sparse_convolve(wavelet,d.member[j]);
    }
}
auto_ptr<ThreeComponentEnsemble> BuildNoiseEnsemble(ThreeComponentEnsemble& patternd,TimeSeries& wavelet,double sigma)
{
    auto_ptr<ThreeComponentEnsemble> result(new ThreeComponentEnsemble(patternd));
    vector<ThreeComponentSeismogram>::iterator dptr;
    for(dptr=result->member.begin();dptr!=result->member.end();++dptr)
    {
        int ns=dptr->ns;
        ns=3*ns;
        double  *ptr=dptr->u.get_address(0,0);
        for(int k=0;k<ns;++k,ptr++)
        {
            *ptr=SampleNormal(0.0,sigma);
        }
        /* Set t0 and endtime fields to an arbitrary number */
        dptr->t0=100000.0;
        dptr->put("time",dptr->t0);
        dptr->put("endtime",dptr->endtime());
        *dptr=sparse_convolve(wavelet,*dptr);
    }
    return(result);
}
void usage()
{
    cerr << "test_seisw [-pf pffile]"<<endl;
    exit(-1);
}
bool SEISPP::SEISPP_verbose(true);
int main(int argc, char **argv)
{
    Pf *pf;
    string pfname("test_seisw");
    int i;
    for(i=1;i<argc;++i)
    {
        string sarg(argv[i]);
        if(sarg=="-pf")
        {
            ++argc;
            pfname=string(argv[i]);
        }
        else
        {
            usage();
        }
    }
    if(pfread(const_cast<char *>(pfname.c_str()),&pf))
    {
        cerr << "pfread failed on required parameter file=test_pwdecon.pf"
            <<endl;
        exit(-1);
    }
    try {
        Metadata control(pf);
        string stafname=control.get_string("station_geometry_file");
        ifstream stain;
        stain.open(stafname.c_str(),ios::in);
        if(stain.fail())
        {
            cerr << "Open failed on station_geometry_file="
                    << stafname<<endl;
            usage();
        }
        vector<double> dnorth,deast;
        char *inputline=new char[512];
        while(stain.getline(inputline,512))
        {
            stringstream ss(inputline);
            double x,y;
            ss>>x;   ss>> y;
            deast.push_back(x);
            dnorth.push_back(y);
        }
        delete [] inputline;
        int nsta=deast.size();
        int nsamp=control.get_int("trace_length");
        /* Create an empty ensemble to be filled with simulation data.*/
        ThreeComponentEnsemble d(nsta,nsamp);
        for(i=0;i<nsta;++i)d.member[i].u.zero();
        /* global data parameters */
        double dt=control.get_double("dt");
        double t0=control.get_double("t0");
        /* For this test we always place pseudostation at lat,lon=0,0 to avoid
           hassles of geographic conversions.  As a result conversion to lat
           and lon in radians is trivial and done is this loop. */
        const double Rearth(6378.164);
        for(i=0;i<nsta;++i){
            d.member[i].put("dnorth",dnorth[i]);
            d.member[i].put("deast",deast[i]);
            double sta_lat,sta_lon;
            sta_lat=dnorth[i]/Rearth;
            sta_lon=deast[i]/Rearth;
            /* PPWDeconProcessor assumes these were loaded in degrees*/
            sta_lat=deg(sta_lat);
            sta_lon=deg(sta_lon);
            d.member[i].put("sta_lat",sta_lat);
            d.member[i].put("sta_lon",sta_lon);
            d.member[i].put("sta_elev",0.0);
            /* Because we used the crude constructor we need to set these too*/
            d.member[i].t0=t0;
            d.member[i].dt=dt;
            d.member[i].live=true;
            d.member[i].tref=relative;
            /* These are required by some functions*/
            d.member[i].put("samprate",1.0/dt);
            d.member[i].put("time",t0);
            d.member[i].put("nsamp",nsamp);
        }
        RectangularSlownessGrid grid(pf,string("SlownessGridParameters"));
        vector<SlownessVector> ulist=gridtovector(grid);
        /* Go ahead and build the wavelet vector immediately.  For now 
        always normalize by peak amplitude  */
        string wavelet_type=control.get_string("wavelet_type");
        double width_parameter=control.get_double("width_parameter");
        int wavelet_length=control.get_int("wavelet_length");
        WaveletNormalizationMethod nm=PEAK;
        TimeSeries wavelet;
	if(wavelet_type=="gaussian")
        {
            wavelet=gaussian_wavelet(wavelet_length,dt,width_parameter,nm);
        }
	else {
	    if(wavelet_type!="ricker")
	    {
	        cerr << "Warning:  unknown wavelet type="<<wavelet_type<<endl
	                << "Defaulting to ricker wavelet"<<endl;
	    }
            wavelet=ricker_wavelet(wavelet_length,dt,width_parameter,nm);
        }
        //SeismicPlot plotter(control);
	//SeismicPlot plotter;
        //plotter.plot(wavelet);
        /* Set assumed incident wave slowness.  Note this is independent
        of Tbl list of slowness vectors used in BuildFakeData below so 
        be aware */
        double ux0=control.get_double("ux0");
        double uy0=control.get_double("uy0");
        SlownessVector u0(ux0,uy0);
        /* noise generation parameters */
        double sigma=control.get_double("noise_sigma");
        /* This procedure loads plane waves into the ensemble with 
        a given wavelet type.  It does a kind of evil thing passing
        some key things through the pf, but this is a test program so
        who cares. */
        BuildFakeData(d,grid,pf,wavelet,sigma);
        // Start with a basic plot.
        cout << "Trying basic plot function"<<endl;
        SeismicPlot bp(control);
        bp.plot(d);
        cout << "Trying trace edit plotter"<<endl;
        TraceEditPlot tp(control);
        tp.plot(d);
        set<int> killlist=tp.report_kills();
        set<int>::iterator iptr;
        for(iptr=killlist.begin();iptr!=killlist.end();++iptr)
        {
            cout << "Killed trace number "<<*iptr<<endl;
        }
        cout << "Trying TimeWindowPicker"<<endl;
        TimeWindowPicker twp(control);
        twp.plot(d);
        TimeWindow win=twp.get();
        cout << "Window picked is" << win.start << " to " <<win.end<<endl; 
    } catch (SeisppError& serr)
    {
        serr.log_error();
    }
}
