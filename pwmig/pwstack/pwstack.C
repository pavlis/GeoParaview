#include <string>
#include "stock.h"
#include "elog.h"
#include "pf.h"
#include "glputil.h"
#include "gclgrid.h"
#include "pwstack.h"

using namespace std;
using namespace SEISPP;

bool Verbose;
void usage()
{
        cbanner((char *)"$Revision: 1.11 $ $Date: 2005/01/16 14:07:37 $",
		(char *)"dbin dbout [-v -V -pf pfname]",
                (char *)"Gary Pavlis",
                (char *)"Indiana University",
                (char *)"pavlis@indiana.edu") ;
	exit(-1);
}
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
	Three_Component_Ensemble *din;
	// This is the output ensemble built on a grid of plane waves
	Three_Component_Ensemble *dout;
	// Tbl tag used in pf to define depth-dependent apeture.
	// frozen here as this constant but passed as a ariable to the 
	// appropriate constructor below
	string aperture_tag="depth_dependent_aperture";
	string mdlist_tag="copy_metadata_list";


	ios::sync_with_stdio();
        /* Initialize the error log and write a version notice */
        elog_init (argc, argv);

        /* usual cracking of command line */
        if(argc < 2) usage();
        dbname_in = argv[1];
        dbname_out = argv[2];

        for(i=2;i<argc;++i)
        {
		if(!strcmp(argv[i],"-V"))
			usage();
		else if(!strcmp(argv[i],"-v"))
			Verbose=true;
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
	try{

	    // This builds the grid of plane wave components
	    Rectangular_Slowness_Grid ugrid(pf,"Slowness_Grid_Definition");
	    // control parameters on stack process
	    double ts,te;
	    ts = pfget_double(pf,(char *)"stack_time_start");
	    te = pfget_double(pf,(char *)"stack_time_end");
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
    	    Metadata_list station_mdl=pfget_mdlist(pf,"station_metadata");
    	    Metadata_list ensemble_mdl=pfget_mdlist(pf,"ensemble_metadata");
    	    Metadata_list stack_mdl=pfget_mdlist(pf,"ensemble_metadata");
	    Attribute_Map am;
    	    Depth_Dependent_Aperture aperture(pf,aperture_tag);
	    Top_Mute mute(pf,string("Data_Top_Mute"));
	    Top_Mute stackmute(pf,string("Stack_Top_Mute"));
	    /* This database open must create a doubly grouped view to
	    // correctly define a three component ensemble grouping 
	    // normally this will be done by a hidden dbprocess list
	    // stored in the master pf directory
	    */
	    string dbnmi(dbname_in);  // this temporary seems necessary for g++
	    Datascope_Handle dbh(dbnmi,false);
	    dbh=Datascope_Handle(dbh.db,pf,string("dbprocess_commands"));
	    // This is the output database which probably should always
	    // be different than the input database.
	    string dbnmo(dbname_out);
	    Datascope_Handle dbho(dbnmo,false);

	    // We need to load the primary GCLgrid that defines the 
	    // location of pseudostation points.  It is assumed we
	    // will build a site and sitechan table in the output database
	    // in that define these virtual stations in a different 
	    // application.  i.e. dbho will eventually need a site table
	    // produced separately to be useful.  This program should not
	    // need this information though.
	    //
	    char *stagridname;
	    stagridname=pfget_string(pf,"Station_GCLgrid_name");
	    string grdnm(stagridname);
	    GCLgrid stagrid(dbh.db,grdnm);

	    int rec;
	    for(rec=0,dbh.rewind();rec<dbh.number_tuples();++rec,++dbh)
	    {
		int iret;
		// ensemble is read once for entire grid in this 
		// version of the program.  This assumes newer 
		// methods like that under development by Fan 
		// will prove better than pseudostation method
		//
		
		din = new Three_Component_Ensemble(
		    dynamic_cast<Database_Handle&>(dbh),
		    station_mdl, ensemble_mdl,am);
		double lat0,lon0,elev0,ux0,uy0;
		for(i=0;i<stagrid.n1;++i) for(j=0;j<stagrid.n2;++j)
		{
			lat0=stagrid.lat(i,j);
			lon0=stagrid.lon(i,j);
			elev0=stagrid.depth(i,j);
			elev0=-elev0;  // z is positive down, elev is positive up
		// ux0 and uy0 are the incident wavefield's 
		// slowness vector.  Stacks are made relative to this
		// and the data are aligned using P wave arrival times
		try {
		    // This code used to contain these two lines which
		    // were then passed to pwstack_ensemble.  We now 
		    // simply treat them as metadata extracted in the
		    // process function for consistency with lat0 and lon0
		    // set below
		    //ux0=din->get_double("ux0");
		    //uy0=din->get_double("uy0");
		    //
		    // this is a backdoor method to pass this information
		    // to the main processing program called just below.
		    // perhaps bad form, but it is logical as these
		    // are metadata related to these data by any measure.
		    din->put_metadata("ix1",i);
		    din->put_metadata("ix2",j);
		    din->put_metadata("lat0",lat0);
		    din->put_metadata("lon0",lon0);
		    din->put_metadata("elev0",elev0);
		    din->put_metadata("gridname",stagridname);
		    iret=pwstack_ensemble(din,
			ugrid,
			mute,
			stackmute,
			ts,
			te,
			aperture,
			station_mdl,
			stack_mdl,
			am,
			dbho);
		    if((iret<0) && (Verbose))
			cout << "Ensemble number "<< rec << " has no data in pseudoarray aperture"<<endl;
		} catch (Metadata_error mderr1)
		{
			mderr1.log_error();
			cerr << "Ensemble " << i << " data skipped" << endl;
		}
		delete din;
	    }
	    }
	} catch (seispp_error err)
	{
		err.log_error();
	}
	catch (seispp_dberror derr)
	{
		derr.log_error();
	}
	catch (Metadata_error mderr)
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

