#include <string>
#include "stock.h"
#include "elog.h"
#include "pf.h"
#include "glputil.h"
#include "pfstream.h"
#include "pwstack.h"
#ifdef MPI_SET
        #include <mpi.h>
#endif

void usage()
{
        cbanner((char *)"$Revision: 1.5 $ $Date: 2003/07/05 14:37:41 $",
		(char *)"pfstreamin pfstreamout [-V -pf pfname]",
                (char *)"Gary Pavlis",
                (char *)"Indiana University",
                (char *)"pavlis@indiana.edu") ;
	exit(-1);
}
int main(int argc, char **argv)
{
	char *pfstreamin, *pfstreamout;
	int i;
	char *pfin=NULL;
	int sift=0;
        char subset_string[128];
	Pf *pf;
	Pf *pfnext;
	Pfstream_handle *pfshi, *pfsho;
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


	        /* Initialize the error log and write a version notice */
        elog_init (argc, argv);

        /* usual cracking of command line */
        if(argc < 3) usage();
        pfstreamin = argv[1];
        pfstreamout = argv[2];

        for(i=2;i<argc;++i)
        {
		if(!strcmp(argv[i],"-V"))
			usage();
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

        /* This utility causes the program to die if required parameters
        are missing */
        check_required_pf(pf);
	// extract all the stuff we need for the main processing
	// routine, pwstack_process, below
	Rectangular_Slowness_Grid ugrid(pf);
	Top_Mute mute(pf);
	Top_Mute stackmute(pf);
	double ts,te;
	ts = pfget_double(pf,(char *)"stack_time_start");
	te = pfget_double(pf,(char *)"stack_time_end");
	ensemble_tag = pfget_string(pf,"pfensemble_tag");
	if(ensemble_tag==NULL)
		elog_die(0,"Missing required parameter pfensemble_tag\n");
	// This constructor can throw an exception we will intentionally
	// ignore.  An uncaught exception will abort the program anyway
	//
	Depth_Dependent_Aperture aperture(pf,aperture_tag);
	list<Metadata_typedef> mdlist;
	mdlist=pfget_mdlist(pf,mdlist_tag);
	dir = pfget_string(pf,"waveform_directory");
	if(dir==NULL)
		dir=strdup("./pwstack");
	if(makedir(dir))
		elog_die(0,"Cannot create directory %s\n",dir);
	

	// Start up reader and writer
	pfshi = pfstream_start_read_thread(pfstreamin);
	pfsho = pfstream_start_write_thread(pfstreamout);
	// this is the object that holds the primary processing input here
	Three_Component_Ensemble *indata_ptr;
	Three_Component_Ensemble& indata=*indata_ptr; //stardard use of reference

	try{
	    while( (indata_ptr = get_next_3c_ensemble(pfshi,ensemble_tag,mdlist))
		!=NULL)
	    {
		double lat0,lon0,ux0,uy0;
		int iret;
		lat0=indata.md.get_double("lat0");
		lon0=indata.md.get_double("lon0");
		ux0=indata.md.get_double("ux0");
		uy0=indata.md.get_double("uy0");
		iret = pwstack_ensemble(indata,
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
			mdlist,
			dir,
			pfsho);
		if(iret)

		// We can release the input data now
		delete indata_ptr;
	    }
	} catch (seispp_error err)
	{
		err.log_error();
		elog_die(0,"Fatal processing error:  cannot fix exception\n");
	}
}	

