#include "db.h"
/* This is a general routine that uses a different concept for using
a database to build a processing flow environment based on shell
scripts.  The procedure opens a special table that is used to define
tables that are to be joined to the core table via a set of keys
defined in the database.   The algorithm is pretty careless about 
the join keys and assumes blindly that the defining tables are 
properly defined.  This should be true because these tables would
be expected to be totally static and changed only when a new
algorithm is added to the "package" that uses this concept.  

Arguments:
	db - database pointer of open database to search
	algorithm - we search for the entries related to this processing
		module and only look at rows that match this field

Return:
	Returns a database pointer to view created by the defining
	join.  Returns dbINVALID in the record field
	of the returned dbptr in the event of an error.
*/
Dbptr dbpwfjoin( Dbptr db, char *algorithm)
{
	Dbptr dbva,dbvb,dbaux;
	char keys[3][10];
	int joinlevel;
	char auxtable[10];
	Tbl *joinkeys,*sortkeys;
	char sset[30];
	int nrows;
	int jkey;

	dbva = dblookup(db,0,"pwfjoin",0,0);
	sprintf(sset,"algorithm =~ /%s/",algorithm);
	dbvb = dbsubset(dbva,sset,0);
	/* These frees are probably not really needed because
	these tables will always be small, but it would cause a
	memory leak without them.*/
	dbfree(dbva);
	
	sortkeys = newtbl(2);
	pushtlb(sortkeys,"algorithm");
	pushtbl(sortkeys,"joinlevel");
	dbva = dbsort(dbvb,sortkeys,0,0);
	dbfree(dbvb);
	freetbl(sortkeys,free);

	dbquery(dbva,dbRECORD_COUNT,&nrec);
	if(nrec <= 0)
	{
		elog_notify(0,"Error in pwfjoin table definitions for algorithm %s\nNo rows in join definitions\n",algorithm);
		dbva.record = dbINVALID;
		return(dbva);
	}
	joinkeys = newtbl(0);
	/* The first table is always the core table */
	dbvb = dblookup(db,0,"wfprocess",0,0);
	for(dbva.record=0;dbva.record<nrec;++dbva.record)
	{
		dbgetv(dbva,0,"auxtable",auxtable,
			"jkey1",keys[0],
			"jkey2",keys[1],
			"jkey3",keys[2],
			0);
		dbaux = dblookup(db,0,auxtable,0,0);
		if(dbaux == dbINVALID)
		{
			elog_notify(0,"Table %s not defined in schema and cannot be joined\n",
				auxtable);
			dbva.record = dbINVALID);
			return(dbva);
		}
		for(jkey=0;jkey<3;++jkey)
		{
			/* This needs to be checked in debugging.  
			I'm guessing a null string is returned 
			when an attibute is not defined*/
			if(keys[jkey] == NULL) break;
			pushtbl(joinkeys,keys[jkey]);
		}
		dbvb = dbjoin(dbvb,dbaux,joinkeys,joinkeys,0,0,0);
		freetbl(joinkeys,0);
	}
	return(dbvb);
}


	

