#include <string>
#include "stock.h"
#include "elog.h"
#include "pf.h"
#include "pfstream.h"
#ifdef MPI_SET
        #include <mpi.h>
#endif

void usage()
{
        cbanner("$Revision: 1.2 $ $Date: 2003/03/07 15:43:57 $",
		"pfstreamin pfstreamout [-V -pf pfname]"
                "Gary Pavlis",
                "Indiana University",
                "pavlis@indiana.edu") ;
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
	Pf_ensemble *pfe;
	Pfstream_handle *pfshi, *pfsho;

	// This is the input data ensemble
	Three_Component_Ensemble *din;
	// This is the output ensemble built on a grid of plane waves
	Three_Component_Ensemble *dout;


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
        if(i != 0) die(1,"Pfread error\n");

        /* This utility causes the program to die if required parameters
        are missing */
        check_required_pf(pf);

	// Start up reader and writer
	pfshi = pfstream_start_read_thread(pfstreamin);
	pfsho = pfstream_start_write_thread(pfstreamout);

	// initialize read and write here.  May change so won't do just
	// now

	while ( (pfe = pfstream_get_next_ensemble(pfshi))!=NULL)
	{
		// may need to crack the data object here with something like
		din = build_ensemble(pfe);
		free_Pf_ensemble(pfe);

		dout = pwstack_ensemble(din, pf);
		// We can release the input data now
		delete din;
		// save output ensemble here 
		delete dout;
	}
	// Residual cleanup here
}	

