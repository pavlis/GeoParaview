/* throw away program for experiment on S to P conversion receiver 
functions.  Code borrows heavily from pwmigcoh and uses some 
frozen constructs as a one shot test.  */

#include <stdio.h>
#include <string>
#include <sstream>
#include "stock.h"
#include "elog.h"
#include "pf.h"
#include "glputil.h"
#include "gclgrid.h"
#include "Hypocenter.h"
#include "seispp.h"
#include "pwstack.h"

using namespace std;
using namespace SEISPP;
#ifdef MATLABDEBUG
MatlabProcessor mp(stdout);
#endif
/* used to set required metadata in trace objects.  Crude bu t
necessary for this special read procedure */
void push_ts_metadata(TimeSeries *d)
{
	d->put("samprate",1.0/d->dt);
	d->put("time",d->t0);
	d->put("nsamp",d->ns);
	d->put("components_are_cardinal",true);
}
/* Specialized procedure to read files produced by Michael Landes
from Seismic Handler.  These are ascii files with a frozen structure
not worth describing as it is too volatile.  Procedure 
reads a list of file names and parses them to return a TimeSeriesEnsemble
object ready to be used in later processing.

This should be viewed as an experimental procedure that should be burned
when we are finished with it.*/

TimeSeriesEnsemble *ReadSHAscii(list<string> fnames)
{
	const int BUFSIZE(1024);
	char line[BUFSIZE];
	/* These are frozen in this file format */
	const double dt(0.2);
	const int ns(600);
	const double t0(-10.0);
	TimeSeriesEnsemble *result=new TimeSeriesEnsemble();
	FILE *fp;  /* plain C file handle for stdio */
	list<string>::iterator fnptr;
	char dummy1[10],dummy2[10];
	double lat,lon;
	double rdummy;
	double val;
	int idummy;
	for(fnptr=fnames.begin();fnptr!=fnames.end();++fnptr)
	{
		fp=fopen(fnptr->c_str(),"r");
		if(fp==NULL)
		{
			cerr << "ReadSHAscii:  cannot open input file="
				<< *fnptr<<endl;
			exit(-1);
		}
		/* drop the first header line.  We don't need it. */
		fgets(line,BUFSIZE,fp);
		int count=0;
		TimeSeries *d;
		d=NULL;
		while(fgets(line,BUFSIZE,fp)!=NULL)
		{
			if(line[0]=='>')
			{
				if(count!=0)
				{
					d->ns=d->s.size();
					push_ts_metadata(d);
					d->live=true;
					result->member.push_back(*d);
				}
				if(d!=NULL) delete d;
				d=new TimeSeries();
				d->dt=dt;
				d->t0=t0;
				/* A crude way to skip the first 58 bytes */
				char *stmp;
				stmp=strdup(line+59);
				sscanf(stmp,"%s%lf%s%lf",dummy1,&lat,dummy2,&lon,0);
				free(stmp);
				d->put("site.lat",lat);
				d->put("site.lon",lon);
				d->put("site.elev",0.0);
			}
			else
			{
				sscanf(line,"%lf%d%lg",&rdummy,&idummy,&val);
				d->s.push_back(val);
			}
			++count;
		}
		/* deal with the last member in this file */
		d->ns=d->s.size();
		push_ts_metadata(d);
		d->live=true;
		result->member.push_back(*d);
	}
	return(result);
}
/* Special kludge for this program.  to avoid changing the pwstack
algorithm need to convert the scalar data to 3c.  Here we will do
this by zeroing x1 and x2 and loading data in x3.  Correct for 
S to P conversion */
auto_ptr<ThreeComponentEnsemble> lto3c(TimeSeriesEnsemble& dz)
{
	auto_ptr<ThreeComponentEnsemble> result(new ThreeComponentEnsemble(dynamic_cast<Metadata&>(dz),
			dz.member.size()));;
	int i,j;
	for(i=0;i<dz.member.size();++i)
	{
		TimeSeries x=dz.member[i];
		ThreeComponentSeismogram d3c(x,false);
		d3c.u.zero();
		for(j=0;j<d3c.ns;++j)
			d3c.u(2,j)=dz.member[i].s[j];
		d3c.live=true;
		result->member.push_back(d3c);
	}
	return(result);
}


void usage()
{
    cbanner((char *)"$Revision: 1.3 $ $Date: 2008/12/29 16:17:40 $",
        (char *)"dbout [-pf pfname] < infilelist",
        (char *)"Gary Pavlis",
        (char *)"Indiana University",
        (char *)"pavlis@indiana.edu") ;
    exit(-1);
}

bool SEISPP::SEISPP_verbose(true);

int main(int argc, char **argv)
{
    char *dbname_out;
    Dbptr db;
    int i,j;
    char *pfin=NULL;
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
    dbname_out = argv[1];

    for(i=2;i<argc;++i)
    {
        if(!strcmp(argv[i],"-pf"))
        {
            ++i;
            if(i>=argc) usage();
            pfin = argv[i];
        }
        else
            usage();
    }
    /* this sets defaults */
    if(pfin == NULL) pfin = strdup("S2Pstacker");

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
/*
        tsfull = pfget_double(pf,"data_time_window_start");
        tefull = pfget_double(pf,"data_time_window_end");
        TimeWindow data_window(tsfull,tefull);
*/
        dir = pfget_string(pf,"waveform_directory");
        if(dir==NULL)
            dir=strdup("./pwstack");
        if(makedir(dir))
            elog_die(0,"Cannot create directory %s\n",dir);
        //
        // station_mdl defines database attributes copied to
        // metadata space of data object.  ensemble_mdl are
        // globals copied to the metadata area for the entire ensemble
        //
        MetadataList station_mdl=pfget_mdlist(pf,"station_metadata");
        MetadataList ensemble_mdl=pfget_mdlist(pf,"ensemble_metadata");
        MetadataList stack_mdl=pfget_mdlist(pf,"stack_metadata");
        AttributeMap OutputAM("css3.0");
        AttributeMap InputAM("css3.0");
        DepthDependentAperture aperture(pf,aperture_tag);
        TopMute mute(pf,string("Data_Top_Mute"));
        TopMute stackmute(pf,string("Stack_Top_Mute"));
        /* new parameters April 2007 for coherence attribute calculator*/
        double dtcoh,cohwinlen;
        dtcoh=pfget_double(pf,"coherence_sample_interval");
        cohwinlen=pfget_double(pf,"coherence_average_length");
	/* New parameter special for this code - fold cutoff */
	int minfold=pfget_int(pf,"stack_fold_cutoff");
        string dbnmo(dbname_out);
        DatascopeHandle dbho(dbnmo,false);

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
        GCLgrid stagrid(dbho.db,grdnm);
	char infile[128];

        //for(rec=0,dbh.rewind();rec<dbh.number_tuples();++rec,++dbh)
	// crude hack to make this work.  Force one pass.
	for(int foo=0;foo<1;++foo)
        {
            int iret;
            int evid;
            double olat,olon,odepth,otime;
	
		/* read list of seismic handler output files */
		list<string> shfnames;
		while(scanf("%s",infile)!=EOF)
		{
			shfnames.push_back(string(infile));
		}

	    TimeSeriesEnsemble *dlcomp=ReadSHAscii(shfnames);
/*
            din = new ThreeComponentEnsemble(
                dynamic_cast<DatabaseHandle&>(dbh),
                station_mdl, ensemble_mdl,InputAM);
*/
/*
            auto_ptr<ThreeComponentEnsemble>
                ensemble=ArrivalTimeReference(*din,"arrival.time",
                data_window);
            // Release this potentially large memory area
            delete din;
*/
            auto_ptr<ThreeComponentEnsemble> ensemble=lto3c(*dlcomp);
	    delete dlcomp;
	    ensemble->put("dir",dir);
            //DEBUG
#ifdef MATLABDEBUG
/*
	    string chans[3]={"x1","x2","x3"};
            mp.load(*ensemble,chans);
            mp.process(string("wigb(x3)"));
            mp.run_interactive();
*/
#endif
            // this should probably be in a try block, but we need to
            // extract it here or we extract it many times later.
/* This was pwstack code for setting origin information 
            evid=ensemble->get_int("evid");
            olat=ensemble->get_double("origin.lat");
            olon=ensemble->get_double("origin.lon");
            odepth=ensemble->get_double("origin.depth");
            otime=ensemble->get_double("origin.time");
*/
/* Frozen values here for Bolivar.  Not general, but let's see if this
proves useful before we make it right. */
		evid=1;
		ensemble->put("evid",evid);
		olat=-8.0;  olon=112.5;  // Almost antipode of gclgrid used here
		odepth=0.0;
		otime=0.0;
		ensemble->put("origin.lat",olat);
		ensemble->put("origin.lon",olon);
		ensemble->put("origin.time",otime);
		ensemble->put("origin.depth",odepth);
		/* Apply top mutes */
		ApplyTopMute(*ensemble,mute);
            // Database has these in degrees, but we need them in radians here.
            olat=rad(olat);  olon=rad(olon);
            double lat0,lon0,elev0;
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

                    cerr << "Working on grid i,j = " << i << "," << j << endl;
                    iret=pwstack_ensemble(*ensemble,
                        ugrid,
                        stackmute,
                        ts,
                        te,
                        aperture,
                        dtcoh,
                        cohwinlen,
			minfold,
                        ensemble_mdl,
                        stack_mdl,
                        OutputAM,
                        dir,
                        dbho);
                } catch (MetadataError& mderr1)
                {
                    mderr1.log_error();
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
