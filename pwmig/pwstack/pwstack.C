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
        cbanner((char *)"$Revision: 1.4 $ $Date: 2003/06/01 16:47:57 $",
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
	ts = pfget_double(pf,"stack_time_start");
	te = pfget_double(pf,"stack_time_end");
	try{ 
		Depth_Dependent_Apeture(pf,aperture_tag);
	} catch (string serr)
	{
		elog_die(0,"pf error:  %s\n",serr.c_str());
	}
	list<Metadata_typedef> *mdlptr;
	mdlptr=pfget_mdlist(pf,mdlist_tag);
	list<Metadata_typedef> mdlist& = *mdlptr;

	// Start up reader and writer
	pfshi = pfstream_start_read_thread(pfstreamin);
	pfsho = pfstream_start_write_thread(pfstreamout);


	while ( (pfnext = pfstream_get_next_ensemble(pfshi))!=NULL)
	{
		Three_Component_Ensemble *indata_ptr;
		indata_ptr = pfstream_load_3censemble(pfnext);
		Three_Component_Ensemble& indata=*indata_ptr;
		// may need to crack the data object here with something like
		//din = build_input_ensemble(pfnext);
		pffree(pfnext);

		//dout = pwstack_ensemble(din, pf);
		// We can release the input data now
		//delete din;
		// save output ensemble here 
		//delete dout;
	}
	// Residual cleanup here
}	

