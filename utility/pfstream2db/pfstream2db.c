/* This program reads a pfstream and saves attributes in a relational database
with a dictionary defined by an control parameter file.  

Usage:  pfstream2db file db [-v -V -pf pffile]

file - input pfstream. This can be a fifo or regular file provided it 
	follows the pfstream format.
-v verbose flag
-V exit with usage line 
-pf use alternative pf control file pffile.  

Primary control is driven by a complex parameter file described in a
man page (or at least will be).
Author:  Gary Pavlis
Written:  Octovber 2002
*/
#include <stdio.h>
#include <string.h>
#include "stock.h"
#include "db.h"
#include "elog.h"
#include "pfstream.h"

void usage()
{
	cbanner("1.0",
		"file db [-V -v -pf pffile -sift expression]",
		"Gary L. Pavlis", "Indiana University",
		"pavlis@indiana.edu");
	exit(-1);
}
/* Cracks a description for a group of attributes headed by "tag" 
and surrounded by &Arr{}.  The tagged group is assumed to be a set of
database table names.  Within the Tbl tagged by that list are attribute
map definitions for that table.  A small example will make the format
clearer:

save_by_row &Arr{
arrival &Tbl{
	arid int arid
	phase string phase
}
assoc &Tbl{
	orid int orid
	delta double delta 
}
}

In this example char *tag="save_by_row"

IMPORTANT:  note because we use the same function used for the
inverse operation in db2pfstream, the entries for each table
are in the order:  dbname type pfname
This is somewhat backward from a copy convention one might
guess (it isn't from to it s to from).  

Author:  Gary Pavlis
Written:  October 2002
*/
Arr *build_attribute_list(char *tag,Pf *pfi)
{
	Tbl *t;
	Arr *a;
	Pf *pf;
	Tbl *table_list;
	int i;

	a = newarr(0);
	/* silently return an empty list if the tag is not defined.
	This allows loops to be driven cleanly by range of the arr*/
	if(pfget(pfi,tag,"(void **)&pf) != PFARR) return(ta);
	table_list = pfkeys(pf);
	for(i=0;i<maxtbl(t);++i)
	{
		char *table;
		table = gettbl(table_list,i);
		/* this libpfstream function pushes Attribute_map objects 
		onto the current list */
		t=pfget_Attribute_map(pf,table);
		/* this places the list of Attribute_map pointers into 
		an associative array keyed by the table name */
		setarr(a,table,t);
	}
	return(a);
}
		

}
int PF2DBSVerbose;
main(int argc, char **argv)
{
	Dbptr db;
	/* These lists contain secondary lists containing pointers
	to Attribute_maps.  They are built from input parameter file 
	descriptions */
	Tbl *save_by_row,*save_by_group,*save_by_ensemble;
	Pf_ensemble *pfe;
	

	PF2DBSVerbose=0;

/*crack the command line*/
	if(argc<3) usage();
	elog_init(argc,argv);
	if(dbopen(argv[2],"r+",&db)==dbINVALID)
	{
		elog_complain(0,"dbopen failed on database %s\n",argv[1]);
		usage();
	}
	fp=fopen(argv[1],"r");
	if(fp==NULL)
	{
		elog_complain(0,"Cannot open %s\n",argv[2]);
		usage();
	}
	for(i=3;i<argc;++i)
	{
		if(!strcmp(argv[i],"-pf"))
		{
			++i;
			if(i>argc) usage();
			pfi=argv[i];
		}
		else if(!strcmp(argv[i],"-V"))
			usage();
		else if(!strcmp(argv[i],"-v"))
			PFS2DBVerbose=1;
		else
			usage();
	}
	
	if(pfi==NULL) pfi=strdup("db2pfstream");
	if(pfread(pfi,&pf))
		elog_die(0,"Error reading parameter file %s.pf\n",pfi);

	save_by_row = build_attribute_list("save_by_row",pf);
	save_by_group = build_attribute_list("save_by_group",pf);
	save_by_ensemble = build_attribute_list("save_by_ensemble",pf);

	/* We have to open the input stream without exclusive access
	priveleges to guarantee a fifo is kept open until the last 
	data block is transmitted.  The read routine in the while
	loop below opens and closes the file for each block.  This
	will cause the last block of data to be flushed and lost when
	the writer we are feeding from terminates.  */
	fd = open(streamin,O_RDONLY);

	while((pfi=pfstream_read(streamin))!=NULL)
	{
		int ii;

		pfe = pfget_Pf_ensemble(pfi,"arrivals");
		for(i=0,ii=0;i<pfe->nmembers;++i)
		{
			
