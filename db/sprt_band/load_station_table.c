#include <math.h>
#include <string.h>
#include "stock.h"
#include "arrays.h"
#include "db.h"
#include "location.h"
/* This function scans a site table and builds an associative
array of Station strctures indexed by the station name.

This is simpler version of a routine in libgenloc.
The MTS prefix was required to avoid collisions with libgenloc names
*/

Arr *MTSload_station_table(Dbptr db)
{
	Arr *a;
	int i;
        Station *s;
	/* Values set by dbgetv from row of view */
	double lat, lon, elev;
	char staname[12];
	int nsite;

	a = newarr(0);
	
	db = dblookup(db,0,"site",0,0);
	dbquery(db,dbRECORD_COUNT,&nsite);
	for(db.record=0;db.record<nsite;++db.record)
	{
		if((dbgetv( db, 0,
			"site.sta",staname,
			"site.lat",&lat,
			"site.lon",&lon,
			"site.elev",&elev,
			0)) == dbINVALID) die(1,"load_station_table:  dbgetv error on station table at row %d\n",
				db.record);
		/* This assumes a station doesn't actually moves.  We use
		this because we can then just pass through the full site 
		table adding a new entry to the arr only when a new row
		is found.  Because site tables have ondata:offdate times
		this is necessary */
		s = (Station *)getarr(a,staname);
		if(s==NULL)
		{
			s = (Station *) malloc(sizeof(Station));
                	if(s == NULL) 
				die(1,"load_station_table:  Cannot malloc Station structure\n");
			strcpy(s->name,staname);
			s->lat = lat;
			s->lon = lon;
			s->elev = elev;
			setarr(a,s->name,s);
		}
	}
	return(a);
}
