#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "elog.h" 
#include "stock.h"
#include "arrays.h"
#include "pf.h"
#include "db.h"
#include "coords.h"
#include "location.h"
#include "dbpmel.h"

/* This small function parses an input string made up of
a comma seperated list of integer ranges returning a Tbl
of ints after expanding the full range defined by the 
string.  An example will make this clearer:
if *gstr is "1,2,4-6" The output Tbl list will be
1,2,4,5,6.  No test is made for overlapping integer
ranges.  

Author: Gary Pavlis
Written:  Sept 2000
*/
Tbl *parse_gridlist_string(char *gstr)
{
	Tbl *t;
	char *sp;
	char *sie;
	int is,ie,i;

	t = newtbl(0);
	sp = strtok(gstr,",");
	do
	{
		if( (sie=strchr(sp,'-')) == NULL)
		{
			is = atoi(sp);
			ie = is;
		}
		else
		{
			*sie = '\0';
			++sie;
			is = atoi(sp);
			ie = atoi(sie);
		}
		for(i=is;i<=ie;++i)
		{
			pushtbl(t,(void *)i);
		}
	} while((sp = strtok(NULL,","))!=NULL);
	return(t);
}
void usage()
{
	elog_die(0,"Usage:  dbpmel db gridlist [-sift expression -pf file]\n\twhere gridlist is a comma seperated list of grid points\n");
}
/* main for dbpmel*/	
void main(int argc, char **argv)
{
	char *dbin;  /* Input db name */
	Tbl *gridlist;
	Dbptr db;  /* input db pointer */
	Dbptr dbv;  /* set to view formed by join */
	char *pfin=NULL;  /* input parameter file */
	char *sift_exp;  /* sift expression for subset */
	int sift = 0;  /* default is no sift.  */
	Tbl *sortkeys;
	Tbl *joinkey1, *joinkey2;
	int nevents;
	/* db row variables */
	int nrows, nrows_raw;

	Pf *pf;
	char *version="1.0";
	int i;

	/* Initialize the error log and write a version notice */
	elog_init (argc, argv) ;
	elog_log (0, "%s version %s\n", argv[0], version) ;

	if(argc < 3) usage();
	dbin = argv[1];
	gridlist = parse_gridlist_string(argv[2]);
	
	for(i=3;i<argc;++i)
	{
		if(!strcmp(argv[i],"-pf"))
		{
			++i;
			if(i>=argc) usage();
			pfin = argv[i];
		}
		else if(!strcmp(argv[i],"-sift"))
		{
			++i;
			if(i>=argc) usage();
			sift_exp = argv[i];
			sift = 1;
		}
		else
			usage();
	}
	/* set default this way*/
	if(pfin == NULL) pfin = strdup("dbpmel");
	i = pfread(pfin,&pf);
	if(i != 0) die(1,"Pfread error\n");
	check_required_pf(pf);

	/* Set up main database view.  This is a derived from code
	in the related genloc program called relocate.
	Always join assoc, arrival, and site.  We join site 
	to make sure station table is properly dynamic to account for
	time changes.  With this setup, the stations can even move
	around and this should still work.*/
	if(dbopen(dbin,"r",&db) == dbINVALID) 
		die(1,"Unable to open input database %s\n",dbin);
	/* We call this routine that uses dbprocess driven by the 
	parameter file definition tagged by the Tbl dbpmel_dbview*/
	dbv = dbform_working_view(db,pf,"dbpmel_dbview_definition");
	dbquery(dbv, dbRECORD_COUNT, &nrows);

	/* We can't use dbprocess for this join as the current 
	form does not allow this use of mixed keys*/
	joinkey1 = newtbl(0);
	joinkey2 = newtbl(0);
	pushtbl(joinkey1,"arrival.sta");
	pushtbl(joinkey1,"arrival.time");
	pushtbl(joinkey2,"sta");
	pushtbl(joinkey2,"ondate::offdate");
	dbv = dbjoin ( dbv, dblookup(db,0,"site",0,0),
			&joinkey1,&joinkey2,0,0,0);
	if(dbv.table == dbINVALID)
		die(1,"cluster->event->origin->assoc->arrival->site join failed\n");

	/* Subset using sift_key if requested */
	if(sift)
	{
		dbv = dbsubset(dbv,sift_exp,0);
		if(dbv.record == dbINVALID)
			die(1,"dbsubset of %s with expression %s failed\n",
				dbin, sift_exp);
	}

	/* First we have to run a unique key sort in the following order
	to remove redundant picks made on multiple channels.  We will
	issue a warning if the record count changes.  This was found
	to be a common problem that had to be repaired automatically.*/
	dbquery(dbv, dbRECORD_COUNT, &nrows_raw);
	sortkeys = newtbl(0);
	pushtbl(sortkeys,"gridid");
	pushtbl(sortkeys,"evid");
	pushtbl(sortkeys,"sta");
	pushtbl(sortkeys,"phase");
	dbv = dbsort(dbv,sortkeys,dbUNIQUE,0);
	dbquery(dbv, dbRECORD_COUNT, &nrows);
	if(nrows != nrows_raw)
		complain(0,"Input database has duplicate picks of one or more phases on multiple channels\n\
Which picks will be used here is unpredictable\n\
%d total picks, %d unique\nContinuing\n", nrows_raw, nrows);

	if(dbpmel_process(dbv,gridlist,pf))
	{
		elog_complain(0,"Errors in dbpmel_process\n");
		exit(-1);
	}
	exit(0);
}
