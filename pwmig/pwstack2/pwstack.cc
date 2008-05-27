#include <string>
#include <sstream>
#include "stock.h"
#include "elog.h"
#include "pf.h"
#include "glputil.h"
#include "gclgrid.h"
#include "Hypocenter.h"
#include "seispp.h"
#include "PwmigFileHandle.h"
#include "pwstack.h"

using namespace std;
using namespace SEISPP;

bool Verbose;

void usage()
{
    cbanner((char *)"$Revision: 1.10 $ $Date: 2008/05/27 20:42:59 $",
        (char *)"dbin [-v -V -pf pfname]",
        (char *)"Gary Pavlis",
        (char *)"Indiana University",
        (char *)"pavlis@indiana.edu") ;
    exit(-1);
}


void load_file_globals(PwmigFileHandle& fh,int evid, double olat, double olon, 
	double odepth, double otime, string gridname)
{
	fh.filehdr.evid=evid;
	fh.filehdr.slat=olat;
	fh.filehdr.slon=olon;
	fh.filehdr.sdepth=odepth;
	fh.filehdr.stime=otime;
	strncpy(fh.filehdr.gridname,gridname.c_str(),16);
}
#ifdef MATLABDEBUG
MatlabProcessor mp(stdout);
#endif

bool SEISPP::SEISPP_verbose(false);
int main(int argc, char **argv)
{
    char *dbname_in, *dbname_out;
    Dbptr db;
    int i,j;
    char *pfin=NULL;
    int sift=0;
    char subset_string[128];
    Pf *pf;
    Pf *pfnext;
    char *ensemble_tag;
    char *dir;

    // This is the input data ensemble
    ThreeComponentEnsemble *din;
    // Tbl tag used in pf to define depth-dependent apeture.
    // frozen here as this constant but passed as a ariable to the
    // appropriate constructor below
    string aperture_tag="depth_dependent_aperture";
    string mdlist_tag="copy_metadata_list";

    ios::sync_with_stdio();
    /* Initialize the error log and write a version notice */
    //elog_init (argc, argv);

    /* usual cracking of command line */
    if(argc < 2) usage();
    dbname_in = argv[1];

    for(i=2;i<argc;++i)
    {
        if(!strcmp(argv[i],"-V"))
            usage();
        else if(!strcmp(argv[i],"-v"))
        {
            Verbose=true;
            SEISPP::SEISPP_verbose=true;
        }
        else if(!strcmp(argv[i],"-pf"))
        {
            ++i;
            if(i>=argc) usage();
            pfin = argv[i];
        }
        else
            usage();
    }
    /* this sets defaults */
    if(pfin == NULL) pfin = strdup("pwstack");

    i = pfread(pfin,&pf);
    if(i != 0) die(1,(char *)"Pfread error\n");
    try
    {

        // This builds the grid of plane wave components
        RectangularSlownessGrid ugrid(pf,"Slowness_Grid_Definition");
        // control parameters on stack process
        double ts,te;
        ts = pfget_double(pf,(char *)"stack_time_start");
        te = pfget_double(pf,(char *)"stack_time_end");
        // the data are windowed around arrivals to this interval
        // normally should be ts:te with sufficient allowance for moveout
        double tsfull, tefull;                    // normally longer than ts and te to allow
        tsfull = pfget_double(pf,"data_time_window_start");
        tefull = pfget_double(pf,"data_time_window_end");
        TimeWindow data_window(tsfull,tefull);
        dir = pfget_string(pf,"waveform_directory");
        if(dir==NULL)
            dir=strdup("./pwstack");
        if(makedir(dir))
            elog_die(0,"Cannot create directory %s\n",dir);
	char *dfilebase=pfget_string(pf,"output_data_file_base_name");
        //
        // station_mdl defines database attributes copied to
        // metadata space of data object.  ensemble_mdl are
        // globals copied to the metadata area for the entire ensemble
        //
        MetadataList station_mdl=pfget_mdlist(pf,"station_metadata");
        MetadataList ensemble_mdl=pfget_mdlist(pf,"ensemble_metadata");
        AttributeMap InputAM("css3.0");
        DepthDependentAperture aperture(pf,aperture_tag);
	/* These will throw an exception if not defined.  
	We use a parameter to turn these off if desired.
	Not very elegant, but functional. */
        TopMute mute(pf,string("Data_Top_Mute"));
	bool enable_data_mute=pfget_boolean(pf,"enable_data_mute");
	if(enable_data_mute)
		mute.enabled=true;
	else
		mute.enabled=false;
        TopMute stackmute(pf,string("Stack_Top_Mute"));
	bool enable_stack_mute=pfget_boolean(pf,"enable_stack_mute");
	if(enable_stack_mute)
		stackmute.enabled=true;
	else
		stackmute.enabled=false;
        /* new parameters April 2007 for coherence attribute calculator*/
        double dtcoh,cohwinlen;
        dtcoh=pfget_double(pf,"coherence_sample_interval");
        cohwinlen=pfget_double(pf,"coherence_average_length");
        /* This database open must create a doubly grouped view to
        // correctly define a three component ensemble grouping
        // normally this will be done by a hidden dbprocess list
        // stored in the master pf directory
        */
        string dbnmi(dbname_in);                  // this temporary seems necessary for g++
        DatascopeHandle dbh(dbnmi,false);
	string dbviewmode(pfget_string(pf,"database_view_mode"));
	if(dbviewmode=="dbprocess")
        	dbh=DatascopeHandle(dbh,pf,string("dbprocess_commands"));
	else if(dbviewmode=="use_wfdisc")
	{
		dbh.lookup("arrival");
		list<string> j1,j2;
		j1.push_back("wfdisc.time::wfdisc.endtime");
		j2.push_back("arrival.time");
		dbh.leftjoin("wfdisc",j1,j2);
		dbh.natural_join("assoc");
		dbh.natural_join("origin");
		dbh.natural_join("event");
		dbh.subset("orid==prefor");
		dbh.natural_join("sitechan");
		dbh.natural_join("site");
		list<string> sortkeys;
		sortkeys.push_back("evid");
		sortkeys.push_back("sta");
		sortkeys.push_back("chan");
		dbh.sort(sortkeys);
	}
	else if(dbviewmode=="use_wfprocess")
	{
		dbh.lookup("event");
		dbh.natural_join("origin");
		dbh.subset("orid==prefor");
		dbh.natural_join("assoc");
		dbh.natural_join("arrival");
		list<string> j1,j2;
		j1.push_back("arrival.time");
		j2.push_back("wfprocess.time::wfprocess.endtime");
		dbh.join("wfprocess",j1,j2);
		dbh.natural_join("sclink");
		dbh.subset("sclink.sta == arrival.sta");
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
	cout << "Processing begins on database " 
		<<  dbname_in << endl
		<<"Working database has "<<dbh.number_tuples() <<" rows."<<endl;
        list<string> group_keys;
        group_keys.push_back("evid");
        dbh.group(group_keys);
        dbh.rewind();
	cout << "This run will process data from "<<dbh.number_tuples()<<" events."<<endl;

        // We need to load the primary GCLgrid that defines the
        // location of pseudostation points.  It is assumed we
        // will build a site and sitechan table in the output database
        // in that define these virtual stations in a different
        // application.  i.e. dbho will eventually need a site table
        // produced separately to be useful.  This program should not
        // need this information though.
        //
        char *stagridname;
        stagridname=pfget_string(pf,"pseudostation_grid_name");
        string grdnm(stagridname);
        GCLgrid stagrid(dbh.db,grdnm);

	char ssbuf[256];
        int rec;
        //DEBUG
        /*
        MatlabProcessor mp(stdout);
        */
        for(rec=0,dbh.rewind();rec<dbh.number_tuples();++rec,++dbh)
        {
            int iret;
            int evid;
            double olat,olon,odepth,otime;
            // ensemble is read once for entire grid in this
            // version of the program.  This assumes newer
            // methods like that under development by Fan
            // will prove better than pseudostation method
            //

            din = new ThreeComponentEnsemble(
                dynamic_cast<DatabaseHandle&>(dbh),
                station_mdl, ensemble_mdl,InputAM);
            //DEBUG
#ifdef MATLABDEBUG
            string chans[3]={"x1","x2","x3"};
            mp.load(*din,chans);
            mp.process(string("wigb(x1);figure;wigb(x2);figure;wigb(x3);"));
            mp.run_interactive();
#endif
            auto_ptr<ThreeComponentEnsemble>
                ensemble=ArrivalTimeReference(*din,"arrival.time",
                data_window);
            // Release this potentially large memory area
            delete din;
            //DEBUG
#ifdef MATLABDEBUG
            mp.load(*ensemble,chans);
            mp.process(string("figure; wigb(x1)"));
            mp.process(string("figure; wigb(x2)"));
            mp.run_interactive();
#endif
            // this should probably be in a try block, but we need to
            // extract it here or we extract it many times later.
            evid=ensemble->get_int("evid");
            olat=ensemble->get_double("origin.lat");
            olon=ensemble->get_double("origin.lon");
            // Database has these in degrees, but we need them in radians here.
            olat=rad(olat);  olon=rad(olon);
            odepth=ensemble->get_double("origin.depth");
            otime=ensemble->get_double("origin.time");
	    /* A way to asssure this is cleared.  May not be necessary */
	    ssbuf[0]='\0';
	    stringstream ss(ssbuf);
	    ss << dir << "/" << string(dfilebase) << "_" << evid;
	    string dfilebase=ss.str();
	    string dfile=dfilebase+DataFileExtension;
	    string coh3cf=dfilebase+Coh3CExtension;
	    string cohf=dfilebase+CohExtension;
	    PwmigFileHandle dfh(dfile,false,false);
	    load_file_globals(dfh,evid,olat,olon,odepth,otime,stagridname);
	    PwmigFileHandle coh3cfh(coh3cf,false,false);
	    load_file_globals(coh3cfh,evid,olat,olon,odepth,otime,stagridname);
	    PwmigFileHandle cohfh(cohf,false,true);
	    load_file_globals(cohfh,evid,olat,olon,odepth,otime,stagridname);
            double lat0,lon0,elev0;
	    cout << "Beginning process for event id = "<<evid<<endl;
            for(i=0;i<stagrid.n1;++i)
                for(j=0;j<stagrid.n2;++j)
            {
                lat0=stagrid.lat(i,j);
                lon0=stagrid.lon(i,j);
                elev0=stagrid.depth(i,j);
                elev0=-elev0;                     // z is positive down, elev is positive up
                // ux0 and uy0 are the incident wavefield's
                // slowness vector.  Stacks are made relative to this
                // and the data are aligned using P wave arrival times
                try
                {
                    //
                    // this is a backdoor method to pass this information
                    // to the main processing program called just below.
                    // perhaps bad form, but it is logical as these
                    // are metadata related to these data by any measure.
                    ensemble->put("ix1",i);
                    ensemble->put("ix2",j);
                    ensemble->put("lat0",lat0);
                    ensemble->put("lon0",lon0);
                    ensemble->put("elev0",elev0);
                    ensemble->put("gridname",stagridname);
                    Hypocenter hypo(olat,olon,odepth,otime,
                        string("tttaup"),string("iasp91"));
                    SlownessVector slow=hypo.pslow(lat0,lon0,elev0);
                    ensemble->put("ux0",slow.ux);
                    ensemble->put("uy0",slow.uy);

                    if(SEISPP_verbose) cout << "Working on grid i,j = " 
					<< i << "," << j << endl;
                    iret=pwstack_ensemble(*ensemble,
                        ugrid,
                        mute,
                        stackmute,
                        ts,
                        te,
                        aperture,
                        dtcoh,
                        cohwinlen,
                        ensemble_mdl,
                        dfh,
			coh3cfh,
			cohfh);
                    if((iret<0) && (Verbose))
                        cout << "Ensemble number "<< rec
                            << " has no data in pseudoarray aperture"<<endl
                            << "Pseudostation grid point indices (i,j)="
                            << "("<<i<<","<<j<<")"<<endl;
                } catch (MetadataError& mderr1)
                {
                    mderr1.log_error();
                    cerr << "Ensemble " << rec << " data skipped" << endl;
                    cerr << "Pseudostation grid point indices (i,j)="
                        << "("<<i<<","<<j<<")"<<endl;
                }
                catch (SeisppError serr)
                {
                    serr.log_error();
                }
            }
        }
    } catch (SeisppError& err)
    {
        err.log_error();
    }
    catch (MetadataError& mderr)
    {
        mderr.log_error();
    }
    // The GCLgrid library db routines throw simple int exceptions
    // This should only be entered if the GCLgrid constructor fails.
    catch (int gclerr)
    {
        // these all put a message to the antelope elog list
        // so we just call the routine that dumps this error one exit.
        elog_die(1,"GCLgrid constructor failed\n");
    }
}
