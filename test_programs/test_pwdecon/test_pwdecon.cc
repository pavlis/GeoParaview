#include <math.h>
#include <vector>
#include <fstream>
#include "slowness.h"
#include "PPWDeconProcessor.h"
#include "SeismicPlot.h"
#include <Xm/MainW.h>
//#include "MatlabPlotter.h"
using namespace std;
using namespace SEISPP;
/* May want to put this on in dmatrix.h */
/*
double Linf(dmatrix& A) {
    double *aptr;
    aptr=A.get_address(0,0);
    double Ainf=*aptr;
    ++aptr;
    int na=A.rows() * A.columns();
    for(int i=1;i<na;++i,++aptr){
        if((*aptr)>Ainf) Ainf=(*aptr);
    }
    return(Ainf);
}
*/
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
/* Use the above to compute wavelet*d.  Return result.  */
ThreeComponentSeismogram convolve(TimeSeries& wavelet, 
	ThreeComponentSeismogram& d)
{
	if( (wavelet.tref==absolute) || (d.tref==absolute) )
		throw SeisppError(string("Error (convolve procedure): ")
			+ "both functions to be convolved must have "
			+ "relative time base");
	ThreeComponentSeismogram out3c(d);
        int nw=wavelet.ns;
        double *wptr;
        wptr=&(wavelet.s[0]);
	/* Add a generous padding for out3c*/
	int nsout=d.ns+2*nw;

	out3c.u=dmatrix(3,nsout);
        out3c.u.zero();
        out3c.t0=d.t0-(out3c.dt)*static_cast<double>(wavelet.ns);
	out3c.ns=nsout;
        /* oi is the position of the moving index position in out3c */
        int oi=out3c.sample_number(d.t0);
        /* si is the index to the point where the wavelet is to be inserted. offset by 0 of wavelet*/
        int si=oi-wavelet.sample_number(0.0);
        if(si<0) throw SeisppError("Error computed out3c index is less than 0 ");
        int i,k;
        for(i=0;i<d.ns;++i,++si){
            double *sptr=out3c.u.get_address(0,si);
            double *dptr=d.u.get_address(0,i);
            for(k=0;k<3;++k){
                if((*dptr)!=0.0){
                    daxpy(nw,(*dptr),wptr,1,sptr,3);
                }
                ++sptr;
                ++dptr;
            }
            // REmove me when this is known to work
            if((si+nw)>out3c.ns) SeisppError("Error in convolution:  exit to prevent out3c overrun");
        }
	return out3c;
}

enum NormalizeMethod {PEAK, AREA, NONE};

TimeSeries gaussian_wavelet(int n, double dt, double sigma,
	NormalizeMethod normalize)
{
	TimeSeries gwavelet(n);
	double tlength;
	tlength=dt*static_cast<double>(n-1);
        gwavelet.ns=n;
	gwavelet.dt=dt;
	gwavelet.t0=-tlength/2.0;
	gwavelet.live=true;
	gwavelet.tref=relative;
	gwavelet.put("samprate",1.0/dt);
	gwavelet.put("time",gwavelet.t0);
	gwavelet.put("nsamp",n);
	int i;
	for(i=0;i<n;++i)
	{
		gwavelet.s[i]=exp(-pow(gwavelet.time(i)/sigma,2.0));
	}
	double normfactor;
	switch (normalize)
	{
	case PEAK:
                int peaksample=gwavelet.sample_number(0.0);
		normfactor=gwavelet.s[peaksample];
		break;
	case AREA:
		for(i=0,normfactor=0.0;i<n;++i)
			normfactor += gwavelet.s[i];
		break;
	case NONE:
		return gwavelet;
	}
	for(i=0;i<n;++i) gwavelet.s[i]/=normfactor;
	return gwavelet;
}
TimeSeries ricker_wavelet(int n, double dt, double nu,
	NormalizeMethod normalize)
{
	TimeSeries rwavelet(n);
	double tlength;
	tlength=dt*static_cast<double>(n-1);
        rwavelet.ns=n;
	rwavelet.dt=dt;
	rwavelet.t0=-tlength/2.0;
	rwavelet.live=true;
	rwavelet.tref=relative;
	rwavelet.put("samprate",1.0/dt);
	rwavelet.put("time",rwavelet.t0);
	rwavelet.put("nsamp",n);
	double pinu2=M_PI*M_PI*nu*nu;
	int i;
	for(i=0;i<n;++i)
	{
		double t2=rwavelet.time(i);
		t2=t2*t2;
		rwavelet.s[i]=(1.0 - 2.0*pinu2*t2)
				*exp(-pinu2*t2);
	}
	double normfactor;
	switch (normalize)
	{
	case PEAK:
                int peaksample=rwavelet.sample_number(0.0);
		normfactor=rwavelet.s[peaksample];
		break;
	case AREA:
		for(i=0,normfactor=0.0;i<n;++i)
			normfactor += rwavelet.s[i];
		break;
	case NONE:
		return rwavelet;
	}
	for(i=0;i<n;++i) rwavelet.s[i]/=normfactor;
	return rwavelet;
}
void BuildFakeData(ThreeComponentEnsemble& d, RectangularSlownessGrid& g,
        Pf *pf,TimeSeries& wavelet)
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
            for(int k=0;k<3;++k)d.member[j].u(k,lag)=amp[k];
        }
    }
    /* now convolve with the desired wavelet */
    for(j=0;j<nsta;++j){
        d.member[j]=convolve(wavelet,d.member[j]);
    }
}
void usage()
{
    cerr << "test_pwdecon [-pf pffile]"<<endl;
    exit(-1);
}
bool SEISPP::SEISPP_verbose(true);
int main(int argc, char **argv)
{
    Pf *pf;
    string pfname("test_pwdecon");
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
        NormalizeMethod nm=PEAK;
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
        SeismicPlot plotter(control);
	//SeismicPlot plotter;
        plotter.plot(wavelet);
        /* Set assumed incident wave slowness.  Note this is independent
        of Tbl list of slowness vectors used in BuildFakeData below so 
        be aware */
        double ux0=control.get_double("ux0");
        double uy0=control.get_double("uy0");
        SlownessVector u0(ux0,uy0);
        /* This procedure loads plane waves into the ensemble with 
        a given wavelet type.  It does a kind of evil thing passing
        some key things through the pf, but this is a test program so
        who cares. */
        BuildFakeData(d,grid,pf,wavelet);
        /* Manual test of Widget setup inline here.  Temporary */
        /*
        XtSetLanguageProc(NULL, NULL, NULL);
        XtToolkitThreadInitialize();
        Widget toplevel,main_w;
#define APP_CLASS "seismicplot"
        XtAppContext AppContext;
        toplevel = XtVaAppInitialize(&AppContext, APP_CLASS,NULL,0,
                                &argc,argv, NULL,NULL);
        main_w = XtVaCreateManagedWidget ((char *)"main window",
                xmMainWindowWidgetClass, toplevel,
                XmNscrollBarDisplayPolicy, XmAS_NEEDED,
                XmNscrollingPolicy,        XmAUTOMATIC,
                XmNwidth,1000,
                XmNheight,800,
                NULL);

        
        Arg args[4];
        int nargs=0;
        XtSetArg(args[nargs],(char *) ExmNzoomFactor,100); nargs++;
        XtSetArg(args[nargs],(char *) ExmNdisplayOnly,1); nargs++;
        Widget seisw[2];
        seisw[0]=ExmCreateSeisw(main_w,(char *)"Seisw",args,nargs);
        auto_ptr<TimeSeriesEnsemble>x1=ExtractComponent(d,0);
        TimeSeriesEnsemble *tse=x1.get();
        XtVaSetValues(seisw[0],ExmNseiswMetadata,(XtPointer)(&control),
                ExmNseiswEnsemble,(XtPointer)(tse),NULL);
        XtRealizeWidget(toplevel);
        XtAppMainLoop(AppContext);
        */

        /* END temporary test of X stuff */
        auto_ptr<TimeSeriesEnsemble>x1=ExtractComponent(d,0);
        //plotter.plot(*x1);
/*  DEBUG
cout << "Dump of fake data"<<endl;
for(int im=0;im<d.member.size();++im){
    cout << "Member = "<<im<<endl;
    cout << dynamic_cast<BasicTimeSeries&>(d.member[im])<<endl
    << dynamic_cast<Metadata&>(d.member[im])<<endl;
    cout << "max value of data="<<Linf(d.member[im].u) <<endl;
}
*/
        /* NOTE:  first version of this tests only fixed wavelet data.
        To test variable wavelet will need to add a second call to 
        BuildFakeData with a different wavelet.  Then sum the two ensembles
        to get the mixed wavelet result. */

        PPWDeconProcessor deconoperator(d,u0,ulist,wavelet,control);
        /* run this for 0 0 pseudostation point as initial test. This will
           expand below */
        auto_ptr<ThreeComponentEnsemble> decondat=deconoperator.compute(0.0,0.0);

    } catch (SeisppError& serr)
    {
        serr.log_error();
    }
}
