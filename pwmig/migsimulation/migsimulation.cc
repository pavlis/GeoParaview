#include <fstream>
#include <vector>
#include <set>
#include "stock.h"
#include "seispp.h"
#include "Metadata.h"
#include "SimplePSPrimarySynthetic.h"
#include "filter++.h"
#include "SimpleWavelets.h"
enum SyntheticType {SIMPLE,POINTSOURCE,EXTERN_LAYERED};
using namespace std;
using namespace SEISPP;
bool SEISPP::SEISPP_verbose(false);
SyntheticSeismogram *CreateSimpleGenerator(Metadata& control, Pf *pf)
{
    try {
        Tbl *t;
        t=pfget_tbl(pf,"interfaces");
        vector<double> d,a;
        int i;
        for(i=0;i<maxtbl(t);++i)
        {
            char *line;
            line=(char *)gettbl(t,i);
            double depth,amp;
            sscanf(line,"%lf%lf",&depth,&amp);
            d.push_back(depth);
            a.push_back(amp);
        }
        SimplePSPrimarySynthetic *result = new SimplePSPrimarySynthetic(control,d,a);
        return result;
    } catch (...) {throw;};
}
void build_db_view(DatascopeHandle& dbh, Metadata& control,Pf *pf)
{
    try{
        /* The following is copied exactly from pwstack.  Not very
           elegant, but a way to guarantee the gathers we bring
           in are compatible with pwstack.  If that code changes this
           will need to change too. */
        string dbviewmode(control.get_string("database_view_mode"));
	if(dbviewmode=="dbprocess")
        	dbh=DatascopeHandle(dbh,pf,string("dbprocess_commands"));
	else if(dbviewmode=="use_wfdisc")
	{
		dbh.lookup("arrival");
		list<string> j1,j2;
		j1.push_back("sta");
		j1.push_back("wfdisc.time::wfdisc.endtime");
		j2.push_back("sta");
		j2.push_back("arrival.time");
		dbh.leftjoin("wfdisc",j1,j2);
		dbh.natural_join("assoc");
		dbh.natural_join("origin");
		dbh.natural_join("event");
		dbh.subset("orid==prefor");
		dbh.natural_join("sitechan");
		dbh.natural_join("site");
		if(SEISPP_verbose) cout << "working view size="
			<<dbh.number_tuples()<<endl;
		list<string> sortkeys;
		sortkeys.push_back("evid");
		sortkeys.push_back("sta");
		sortkeys.push_back("chan");
		dbh.sort(sortkeys);
		list<string> gkey;
		gkey.push_back("evid");
		gkey.push_back("sta");
		dbh.group(gkey);
	}
	else if(dbviewmode=="use_wfprocess")
	{
		dbh.lookup("event");
		dbh.natural_join("origin");
		dbh.subset("orid==prefor");
		dbh.natural_join("assoc");
		dbh.natural_join("arrival");
		DatascopeHandle ljhandle(dbh);
		ljhandle.lookup("wfprocess");
		ljhandle.natural_join("sclink");
		ljhandle.natural_join("evlink");
		list<string> jk;
		jk.push_back("evid");
		jk.push_back("sta");
		dbh.join(ljhandle,jk,jk);
		dbh.natural_join("site");
		list<string> sortkeys;
		sortkeys.push_back("evid");
		sortkeys.push_back("sta");
		dbh.sort(sortkeys);
	}
	else
	{
		cerr << "Illegal option for parameter database_view_mode="<<dbviewmode;
		exit(-1);
	}
        list<string> group_keys;
        group_keys.push_back("evid");        
        dbh.group(group_keys);        
        dbh.rewind();
    }catch(...) {throw;};
}

set<int> load_eventset(string fname)
{
    ifstream fin;
    fin.open(fname.c_str(),ios::in);
    if(fin.fail()) throw SeisppError(string("load_eventset:  ")
            +"Cannot open file="+fname);
    int val;
    set<int> result;
    if(SEISPP_verbose) cout << "List of Events To Replicate"<<endl;
    while(fin.good())
    {
        fin>>val;
        result.insert(val);
        if(SEISPP_verbose)cout << val<<endl;
    }
    fin.close();
    return(result);
}

void usage()
{
    cerr << "migsimulation dbin dbout [-evf eventlistfile -pf pffile -V]"<<endl;
    exit(-1);
}
int main(int argc, char **argv)
{
    //WARNING WARNING:  frozen for now.  Change when a new synthetic generator is added 
    SyntheticType syntype(SIMPLE);
    if(argc<3) usage();
    string dbin_name, dbout_name, pfname("migsimulation");
    dbin_name=string(argv[1]);
    dbout_name=string(argv[2]);
    string evlfile("NONE");
    bool check_evlist(false);
    int i;
    for(i=3;i<argc;++i)
    {
        string strarg(argv[i]);
        if(strarg=="-pf")
        {
            ++i;
            if(i==argc) usage();
            pfname=string(argv[i]);
        }
        else if(strarg=="-evf")
        {
            ++i;
            if(i==argc) usage();
            evlfile=string(argv[i]);
            check_evlist=true;
        }
        else if(strarg=="-V")
            SEISPP_verbose=true;
        else
            usage();
    }
    SyntheticSeismogram *synbase;
    try {
        /* For the present freeze this to be css3.0 schema */
        AttributeMap am("css3.0");
        DatascopeHandle dbhin(dbin_name,true);
        DatascopeHandle dbhout(dbout_name,false);
        DatascopeHandle dbhsc(dbhout);
        dbhsc.lookup("sclink");
        DatascopeHandle dbhev(dbhout);
        dbhev.lookup("evlink");
        dbhout.lookup("wfprocess");
        Pf *pf;
        if(pfread(const_cast<char *>(pfname.c_str()),&pf))
        {
            cerr << "pfread failed for pffile="<<pfname<<endl;
            usage();
        }
        Metadata control(pf);
        double tsfull = control.get_double("data_time_window_start");
        double tefull = control.get_double("data_time_window_end");
        TimeWindow data_window(tsfull,tefull);
        /* These should be identical to those defined for pwstack */
        MetadataList station_mdl=pfget_mdlist(pf,"station_metadata");
        MetadataList ensemble_mdl=pfget_mdlist(pf,"ensemble_metadata");
        MetadataList output_mdl=pfget_mdlist(pf,"wfprocess_metadata");
        MetadataList wfdisc_mdl=pfget_mdlist(pf,"wfdisc_metadata");
        /* for now we freeze these names for output */
        vector<string> chanmap;
        chanmap.push_back("E");
        chanmap.push_back("N");
        chanmap.push_back("Z");
        string output_dir=control.get_string("output_directory");
        string dfile_base=control.get_string("output_data_file_base_name");
        /* the next group sets if synthetics are filtered with an
           idealized wavelet or with a standard time invariant filter */
        string wavelet_type=control.get_string("wavelet_type");
        TimeSeries wavelet;
        double datadt=control.get_double("data_sample_interval");
        // these are not needed for TimeInvariantFilter, but are 
        // needed for gaussian or ricker wavelets 
        int nwsamp = control.get_int("wavelet_length");
        double wavelet_width_parameter
            =control.get_double("wavelet_width_parameter");
        WaveletNormalizationMethod nm=PEAK;
        TimeInvariantFilter *filter;
        if(wavelet_type=="filter")
        {
            /* Use a none definition for filter to do nothing to synthetic
               output */
            string filparam=control.get_string("filter");
            filter=new TimeInvariantFilter(filparam);
        }
        else if(wavelet_type=="gaussian")
        {
            wavelet=gaussian_wavelet(nwsamp,datadt,
                    wavelet_width_parameter,nm);
        }
        else if(wavelet_type=="ricker")
        {
            wavelet=ricker_wavelet(nwsamp,datadt,
                    wavelet_width_parameter,nm);
        }
        else
        {
            cerr << "wavelet_type parameter="<<wavelet_type
                << " is not supported."<<endl
                << "Must be one of:  filter, gaussian, or ricker"<<endl;
            usage();
        }
        /* build event (evid) list if requested */
        set<int> eventset;
        if(check_evlist)
        {
            eventset=load_eventset(evlfile);
        }

        /* Construct the synthetic seismogram generator objects */
        switch(syntype){
            /* For now only one option and it is also default */
        case SIMPLE:
        default:
            // Necessary because I wanted to not use a pf for this object 
            synbase=CreateSimpleGenerator(control,pf);
        };
        build_db_view(dbhin,control,pf);
        int nevents=dbhin.number_tuples();
        if(SEISPP_verbose) cout << "Processing "<<nevents<<" simulated events"
                            <<endl;
        int rec;
        for(rec=0;rec<nevents;++rec,++dbhin)
        {
            int evid=dbhin.get_int("evid");
            if(check_evlist)
            {
                if(eventset.find(evid)==eventset.end()) continue;
            }
            if(SEISPP_verbose) cout << "Working on evid="<<evid<<endl;
            auto_ptr<ThreeComponentEnsemble> ensemble(
                  new ThreeComponentEnsemble(dynamic_cast<DatabaseHandle&>(dbhin),
                        station_mdl, ensemble_mdl,am));
            ensemble=ArrivalTimeReference(*ensemble,"arrival.time",
                data_window);
            double slat,slon,sz,otime;
            slat=ensemble->get_double("origin.lat");
            slon=ensemble->get_double("origin.lon");
            slat=rad(slat);
            slon=rad(slon);
            sz=ensemble->get_double("origin.depth");
            otime=ensemble->get_double("origin.time");
            //int evid=ensemble->get_int("evid");
            //For now freeze this as taup calculator iasp91
            Hypocenter hypo(slat,slon,sz,otime,string("tttaup"),string("iasp91"));
            vector<ThreeComponentSeismogram>::iterator dptr;
            for(dptr=ensemble->member.begin();
                dptr!=ensemble->member.end();++dptr)
            {
                double rlat,rlon,relev;
                rlat=dptr->get_double("site.lat");
                rlon=dptr->get_double("site.lon");
                rlat=rad(rlat);
                rlon=rad(rlon);
                relev=dptr->get_double("site.elev");
                ThreeComponentSeismogram syndata(synbase->Compute3C(*dptr,
                            hypo,rlat,rlon,-relev,string("")));
                if(wavelet_type=="filter")
                {
                    filter->apply(syndata);
                }
                else
                {
                    syndata=sparse_convolve(wavelet,syndata);
                }
                /* Data are written to event files in a user 
                   selected directory.  This puts this in the trace
                   headers/metadata area. */
                syndata.put("dir",output_dir);
                stringstream ss;
                ss<<dfile_base<<"_"<<evid;
                syndata.put("dfile",ss.str());
                string sta=dptr->get_string("sta");
                syndata.put("sta",sta);
                /* We have to put syndata back to an absolute time 
                   frame.  Done through arrival.time header field */
                double atime=syndata.get_double("arrival.time");
                syndata.rtoa(atime);
                string algorithm("migsimulation");
                syndata.put("wfprocess.algorithm",algorithm);
                int wfrec=dbsave(syndata,dbhout.db,string("wfprocess"),
                        output_mdl,am);
                int pwfid;
                dbhout.db.record=wfrec;
                pwfid=dbhout.get_int("pwfid");
                /* use datatype for chan here.  May be trouble, but 
                   a starting point */
                string chan=dbhout.get_string("datatype");
                dbhev.append();
                dbhev.put("pwfid",pwfid);
                dbhev.put("evid",evid);
                dbhsc.append();
                dbhsc.put("sta",sta);
                dbhsc.put("chan",chan);
                dbhsc.put("pwfid",pwfid);
                /* now deal with wfdisc.  This assumes
                 dbsave does a lookup on dbhout for wfdisc.
                 Note wfdisc data will be interleaved with wfprocess
                 data all in one large event file*/
                wfrec=dbsave(syndata,dbhout.db,string("wfdisc"),
                        wfdisc_mdl,am,chanmap,true);
            }
        } 
    }catch(SeisppError& serr)
    {
        serr.log_error();
    }
    catch(dmatrix_error& derr)
    {
        derr.log_error();
    }

}

