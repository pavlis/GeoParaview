/* This program builds a pfstream from a specified database view and
related list of attributes.  The pfstream is pushed to an output file
name that would normally be a fifo, but could also be a regular file.

Usage:  db2pfstream db file [-n -v -pf pffile -sift expressin]

where 
db is input database
file is the output to send the stream to
-n dryrun, no output
-v verbose flag
-pf use alternate parameter pffile instead of default 
-sift apply additional subset condition to output view
	(used to avoid cluttering up parameter file script driven runs)

Most of the behaviour of this program is driven by the parameter file,
but I wont' document that here and leave it to the man page.

Author:  GAry Pavlis
Written:  September 2002
*/
#include <stdio.h>
#include <string.h>
#include "stock.h"
#include "db.h"
#include "elog.h"
#include "db2pfstream.h"
#include "glputil.h"

void usage()
{
	cbanner("1.0",
		"db2pfstream db file [-n -V -v -pf pffile -sift expression]",
		"Gary L. Pavlis", "Indiana University",
		"pavlis@indiana.edu");
	exit(-1);
}

/* Output an ensemble to stream fp.  Writes header block using pfhead
then writes npfe pieces with separators defined by a simple numbering
scheme.  Blocks is terminated with __EOF___, which pf(5) says is treated
as and end of file for a stream.  
*/
void pfensemble_to_stream(Pf *pfhead,Pf **pfensemble,int npfe,FILE *fp)
{
	int ipf,i;
	pfout(fp,pfhead);
	fprintf(fp,"ensemble &Arr{\n");
	for(i=0,ipf=1;i<npfe;++i,++ipf)
	{
		fprintf(fp,"%6.6d &Arr{\n",ipf);
		pfout(fp,pfensemble[i]);
		fprintf(fp,"}\n");
	}
	fprintf(fp,"__EOF__\n");
}

/*  This function adds a list of attributes stored in one
row of a datascope database to a pf.  What attributes are extracted
and converted is controlled by the Attribute_map list stored
in t.  

Arguments:
	db - row of database to use to extract info to create or
		update pf.
	t - Tbl containing pointers to Attribute_map objects that
		define which attributes to extract and possible renaming in pf
	pf - destination pf object.

Returns:
	normal return is 0 if there are no problems.  The routine checks
	for null fields requested.  An error count is incremented each time
	this occurs.  Thus a nonzero error count indicates a potential 
	problem with missing data in the output pf.  This condition 
	should be trapped if the attribute requested is essential.  This
	is why the main program below groups attributes into required and
	not required.  When verbose mode is on each null causes a 
	message to be posted on the error log.  

Author:  Gary L. Pavlis
Written:  September 2002
*/
int db2pf(Dbptr db, Tbl *t, Pf *pf)
{
	int i;  
	char *nullvalue;
	int iattrib, inull;
	double dattrib, dnull;
	char sattrib[128];
	Attribute_map *m;
	int error_count=0;
	Dbptr dbf;

	for(i=0;i<maxtbl(t);++i)
	{
		m = (Attribute_map *)gettbl(t,i);
		dbf=dblookup(db,0,0,m->dbname,0);
		dbquery(dbf,dbNULL,&nullvalue);
		switch(m->type)
		{
		case DBREAL:
			if(dbgetv(db,0,m->dbname,&dattrib,0)==dbINVALID)
			{
				register_error(0,"dbgetv error for %s\n",
					m->dbname);
				++error_count;
				break;
			}
			dnull=atof(nullvalue);
			if(dattrib==dnull)
			{
				if(DB2PFS_verbose)
					register_error(0,"Attribute %s is null\n",
						m->dbname);
				++error_count;
				break;
			}
			pfput_double(pf,m->oname,dattrib);
			break;
		case DBINT:
			if(dbgetv(db,0,m->dbname,&iattrib,0)==dbINVALID)
			{
				register_error(0,"dbgetv error for %s\n",
					m->dbname);
				++error_count;
				break;
			}
			inull=atoi(nullvalue);
			if(iattrib==inull)
			{
				if(DB2PFS_verbose)
					register_error(0,"Attribute %s is null\n",
						m->dbname);
				++error_count;
				break;
			}
			pfput_int(pf,m->oname,iattrib);
			break;
		case DBSTR:
		default:
			if(dbgetv(db,0,m->dbname,sattrib,0)==dbINVALID)
			{
				register_error(0,"dbgetv error for %s\n",
					m->dbname);
				++error_count;
				break;
			}
			/*conversion not necessary here since dbNULL
			query returns a string for null value*/
			if(!strcmp(sattrib,nullvalue))
			{
				if(DB2PFS_verbose)
					register_error(0,"Attribute %s is null\n",
						m->dbname);
				++error_count;
				break;
			}
			pfput_string(pf,m->oname,sattrib);
		}
	}
	return(error_count);
}
/* see top of file for description of the main program here */
	

int DB2PFS_verbose;
main(int argc, char **argv)
{
	Dbptr db, dbv, dbge, dbgg;
	/*ensemble_mode is primary grouping, grouping_on set if secondary on*/
	int ensemble_mode,grouping_on;
	Tbl *ensemble_keys;
	Tbl *group_keys;
	int i;
	char *pfi=NULL;
	Pf *pf;
	Pf *pfo_head;
	Tbl *process_list;
	int dryrun=0;
	char *sift_exp;
	int sift;
/*Attributes listed in require will abort program if missing
passthrough parameters will generate an error message but not
cause the program to abort.  Both contain pointers to Attribute_map*/
	Tbl *require,*passthrough;
	Attribute_map *amap;
	FILE *fp;
	int ierr;

	DB2PFS_verbose=0;

/*crack the command line*/
	if(argc<3) usage();
	elog_init(argc,argv);
	if(dbopen(argv[1],"r",&db)==dbINVALID)
	{
		elog_complain(0,"dbopen failed on database %s\n",argv[1]);
		usage();
	}
	fp=fopen(argv[2],"r+");
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
			DB2PFS_verbose=1;
		else if(!strcmp(argv[i],"-n"))
			dryrun=1;
		else if(!strcmp(argv[i],"-sift"))
		{
			++i;
			if(i>argc) usage();
			sift=1;
			sift_exp=argv[i];
		}
		else
			usage();
	}
	
	if(pfi==NULL) pfi=strdup("db2pfstream");
	if(pfread(pfi,&pf))
		elog_die(0,"Error reading parameter file %s.pf\n",pfi);
	ensemble_mode = pfget_boolean(pf,"ensemble_mode");
	if(ensemble_mode)
	{
		ensemble_keys=pfget_tbl(pf,"ensemble_keys");
		if(ensemble_keys==NULL)
			elog_die(0,
			 "ensemble_mode is on, but no grouping keys defined\nCheck parameter file\n");
		group_keys=pfget_tbl(pf,"grouping_keys");
		if(group_keys==NULL)
			grouping_on = 0;
		else
			grouping_on = 1;
	}
	if(DB2PFS_verbose && ensemble_mode)
	{
		char sm[128];
		char *key;
		strcpy(sm,"Defining ensemble with keys: ");
		for(i=0;i<maxtbl(ensemble_keys);++i)
		{
			key=gettbl(ensemble_keys,i);
			strcat(sm,key);
		}
		elog_log(0,"%s\n",sm);
		if(grouping_on)
		{
			strcpy(sm,"Grouping ensemble with keys:  ");
			for(i=0;i<maxtbl(group_keys);++i)
			{
				key = gettbl(group_keys,i);
				strcat(sm,key);
			}
			elog_log(0,"%s\n",sm);
		}
	}

	
	/*This function calls dbprocess using a tbl from pf*/
	dbv=dbform_working_view(db,pf,"dbprocess_list");
	if(dbv.record==dbINVALID)
		elog_die(0,"dbprocess failed:  check database and parameter file parameter dbprocess_list\n");
	/*Now we form the groupings if needed */
	if(ensemble_mode)
	{
		dbge = dbgroup(dbv,ensemble_keys,"ensemble",1);
		if(dbge.record == dbINVALID)
			elog_die(0,"dbgroup with ensemble keys failed\n");
		if(grouping_on)
		{
			dbgg = dbgroup(dbv,group_keys,"groups",2);
			if(dbgg.record == dbINVALID)
				elog_die(0,"dbgroup on secondary grouping keys failed\n");
		}
	}

/*This builds the maps of which attributes will be saved */
	require = pfget_Attribute_map(pf,"require");
	if(require==NULL)
		elog_die(0,"Error in parameter file for parameter require\n");
	passthrough = pfget_Attribute_map(pf,"passthrough");
	if(passthrough==NULL)
		elog_die(0,"Error in parameter file for parameter passthrough\n");
	/* Unfortunately, we have two different loops for ensemble mode
	 and single object mode */
	if(ensemble_mode)
	{
		int number_ensembles;
		Dbptr db_bundle_1;
		Dbptr db_bundle_2;
		int is_ensemble, ie_ensemble;
		int isg, ieg;
		int number_groups;

		pfo_head=pfnew(0);
		pfput_tbl(pfo_head,"ensemble_keys",ensemble_keys);
		if(grouping_on) pfput_tbl(pfo_head,"grouping_keys",group_keys);
	

		dbquery(dbge,dbRECORD_COUNT,&number_ensembles);
		if(DB2PFS_verbose) elog_log(0,"database has %d ensembles\n",
					number_ensembles);
		if(grouping_on) 
		{
			dbgg.record=0;
			dbquery(dbgg,dbRECORD_COUNT,&number_groups);
			if(DB2PFS_verbose)
			  elog_log(0,"database has %d subgroups\n",
				number_groups);
		}
		for(dbge.record=0;dbge.record<number_ensembles;
					++dbge.record)
		{
			Pf **pfensemble;
			int nmembers;
			char *grp_records;
			Tbl *grplist;

			dbgetv(dbge,0,"bundle",&db_bundle_1);
			dbget_range(db_bundle_1,&is_ensemble,&ie_ensemble);
			nmembers = ie_ensemble-is_ensemble+1;
			allot(Pf **,pfensemble,nmembers);
			/*Now loop through and load the array of pf's
			with parameters defined by the attribute maps*/
			for(i=0,dbv.record=is_ensemble;i<nmembers;
				++i,++dbv.record)
			{
				pfensemble[i]=pfnew(PFARR);
				ierr=db2pf(dbv,require,pfensemble[i]);
				if(ierr!=0)
				{
				    if(dryrun)
					elog_notify(0,"%d database errors in ensemble %d need to be fixed\n",
						ierr,dbge.record);
				    else
					elog_die(0,"%d database errors in ensemble %d\nCannot continue without required attributes\n",
						ierr,dbge.record);
				}
				ierr=db2pf(dbv,passthrough,pfensemble[i]);
				if(DB2PFS_verbose)
					elog_log(0,"%d errors in extracting passthrough parameters for ensemble %d\n",
						ierr,dbge.record);
			}
			
			if(grouping_on)
			{
				grplist=newtbl(0);
				do {
				     dbgetv(dbgg,0,"bundle",&db_bundle_2,0);
				     dbget_range(db_bundle_2,&isg,&ieg);
				     grp_records=malloc(30);
				     sprintf(grp_records,"%d %d",
					isg-is_ensemble,ieg-is_ensemble);
				     pushtbl(grplist,grp_records);
				     ++dbgg.record;
				}
				while (ieg<ie_ensemble);
				pfput_tbl(pfo_head,"group_records",grplist);
			}

pfout(stdout,pfo_head);
/*
			pfensemble_to_stream(pfo_head,pfensemble,nmembers,fp);
*/
					
			for(i=0;i<nmembers;++i)
			{
				pffree(pfensemble[i]);
			}
			free(pfensemble);
		}
	}
	else
	{
	/*This is the simpler, single block mode */
		int nrecords;
		Pf *pfo;

		dbquery(dbv,dbRECORD_COUNT,&nrecords);
		pfo=pfnew(PFARR);
		for(dbv.record=0;dbv.record<nrecords;++dbv.record)
		{
			ierr=db2pf(dbv,require,pfo);
			if(ierr!=0)
			{
			    if(dryrun)
				elog_notify(0,"%d database errors for record %d need to be fixed\n",
						ierr,dbv.record);
			    else
				elog_die(0,"%d database errors in record %d\nCannot continue without required attributes\n",
						ierr,dbv.record);
			}
			ierr=db2pf(dbv,passthrough,pfo);
			if(DB2PFS_verbose)
				elog_log(0,"%d errors in extracting passthrough parameters for ensemble %d\n",
						ierr,dbv.record);
			pfout(fp,pfo);
			fprintf(fp,"__EOF__\n");
		}
	}
	exit(0);
}
