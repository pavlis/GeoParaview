#include <stdio.h>
#include <strings.h>
#include "stock.h"
#include "elog.h"
#include "db.h"
#include "pwmig.h"
void usage()
{
	elog_die(0,"Usage:  pwstack db [-sift x -pf pf]\n");
}
main(int argc, char **argv)
{
	char *version="1.0 summer 2000\nAuthor:  Gary Pavlis";
	char *dbname;  /* input database name */
	Dbptr db,dbgrp;
	int i;
	char *pfin=NULL;
	int sift=0;
        char subset_string[128];
        Tbl *grpkeys;
        Tbl *sortkeys;
	Tbl *joinkeys;
	Pf *pf;
	char *gridname;
	int ntraces, ngathers;
	GCL2dgrid *grd;

	        /* Initialize the error log and write a version notice */
        elog_init (argc, argv);
        elog_log (0, "%s version %s\n", argv[0], version) ;

        /* usual cracking of command line */
        if(argc < 2) usage(argv[0]);
        dbname = argv[1];

        for(i=2;i<argc;++i)
        {
              if(!strcmp(argv[i],"-pf"))
                {
                        ++i;
                        if(i>=argc) usage(argv[0]);
                        pfin = argv[i];
                }
                else if(!strcmp(argv[i],"-sift"))
                {
                        ++i;
                        if(i>=argc) usage(argv[0]);
                        sift_exp = argv[i];
                        sift = 1;
                }
                else
                        usage(argv[0]);
        }
        /* this sets defaults */
        if(pfin == NULL) pfin = strdup("pwstack");

        i = pfread(pfin,&pf);
        if(i != 0) die(1,"Pfread error\n");

        /* This utility causes the program to die if required parameters
        are missing */
        check_required_pf(pf);

	if(dbopen(dbname,"r+",&db) == dbINVALID)
                die(1,"Unable to open input database %s\n",dbname);
	db = dblookup(db,0,"deconwfdisc",0,0);
	gridname = pfget_string(pf,"pseudostation_gridname");
	sprintf(subset_string,"gridname =~ /%s/",gridname);
	db = dbsubset(db,subset_string,0);
	if(dbj.record == dbINVALID)
		die(1,"dbsubset of %s with expression %s failed\n",
                                dbname, subset_string);
	joinkeys = newtlb(1);
	pushtbl(joinkeys,"sta");
	db = dbjoin(db,dblookup(db,0,"site",0,0),joinkeys,joinkeys,
		  		0,0,0);
	if(db.table == dbINVALID)
		die(1,"deconwfdisc->site join failed");

	sorkeys = newtbl(5);
	pushtbl(sortkeys,"evid");
	pushtbl(sortkeys,"ix1");
	pushtbl(sortkeys,"ix2");
	pushtbl(sortkeys,"sta");
	pushtbl(sortkeys,"chan");
	db = dbsort(db,sortkeys,0,0);
	if(db.record == dbINVALID)
		die(1,"dbsort of deconwfdisc->site view failed\n");

	grpkeys = newtbl(3);
	pushtbl(grpkeys,"evid");
	pushtbl(grpkeys,"ix1");
	pushtbl(grpkeys,"ix2");
	dbgrp=dbgroup(db,grpkeys,0,0);

	dbquery(db,dbRECORD_COUNT,&ntraces);
	dbquery(dbgrp,dbRECORD_COUNT,&ngathers);
	elog_log(0,"Beginning processing of %d seismograms in %d gathers\n",
		ntraces,ngathers);

	/* This routine loads the grid identified by gridname */
	grd = GCL2dgrid_load_db(db,gridname);

	/* Now we loop through each group and process the 
	individual gathers they define */
	for(dbgrp.record=0;dbgrp.record<ngathers;++dbgrp.record)
	{
		int retcode;

		retcode = pwstack_process(dbgrp,grd,pf);
	}
}	

