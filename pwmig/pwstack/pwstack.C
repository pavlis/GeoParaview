#include <string>
#include "stock.h"
#include "elog.h"
#include "pf.h"
#include "glputil.h"
#include "pwstack.h"

using namespace std;
using namespace SEISPP;

bool Verbose;
void usage()
{
        cbanner((char *)"$Revision: 1.9 $ $Date: 2004/09/23 20:36:00 $",
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
	int i;
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

	if(dbopen(dbname_in,"r+",&db))
		die(1,"Cannot open input database");

	// extract all the stuff we need for the main processing
	// routine, pwstack_process, below
	Rectangular_Slowness_Grid ugrid(pf,"Slowness_Grid_Definition");
	double ts,te;
	ts = pfget_double(pf,(char *)"stack_time_start");
	te = pfget_double(pf,(char *)"stack_time_end");
	dir = pfget_string(pf,"waveform_directory");
	if(dir==NULL)
		dir=strdup("./pwstack");
	if(makedir(dir))
		elog_die(0,"Cannot create directory %s\n",dir);
	try{
	   // station_mdl defines database attributes copied to 
	   // metadata space of data object.  ensemble_mdl are 
	   // globals copied to the metadata area for the entire ensemble
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
	    Datascope_Handle dbh(db,pf,
		string("dbprocess_commands"));
	    Datascope_Handle dbho(string(dbname_out),false);
	    for(i=0,dbh.rewind();i<dbh.number_tuples();++i,++dbh)
	    {
		int iret;
		din = new Three_Component_Ensemble(
		    dynamic_cast<Database_Handle&>(dbh),
		    station_mdl, ensemble_mdl,am);
		double lat0,lon0,ux0,uy0;
		lat0=din->get_double("lat0");
		lon0=din->get_double("lon0");
		// ux0 and uy0 are the incident wavefield's 
		// slowness vector.  Stacks are made relative to this
		// and the data are aligned using P wave arrival times
		ux0=din->get_double("ux0");
		uy0=din->get_double("uy0");
		try {
		    iret=pwstack_ensemble(din,
			ugrid,
			mute,
			stackmute,
			lat0,
			lon0,
			ux0,
			uy0,
			ts,
			te,
			aperture,
			station_mdl,
			stack_mdl,
			am,
			dbho);
		    if((iret<0) && (Verbose))
			cout << "Ensemble number "<< i << " has no data in pseudoarray aperture"<<endl;
		} catch (Metadata_error mderr1)
		{
			mderr1.log_error();
			cerr << "Ensemble " << i << " data skipped" << endl;
		}
		delete din;
	    }
	} catch (seispp_error err)
	{
		err.log_error();
	}
	catch (Metadata_error mderr)
	{
		mderr.log_error();
	}
}	

