#include "stock.h"
#include "coords.h"
#include "arrays.h"
#include "db.h"
#include "elog.h"
#include "pf.h"
#include "gclgrid.h"
#include "glputil.h"
void usage()
{
	banner(Program_Name, "$Revision: 1.5 $ $Date: 2001/08/05 08:34:06 $") ;
	elog_die(0,"usage:  %s db [-pf pfname]\n",Program_Name);
}
/* This is an internal definition of an event location object that holds
all the things we need here.
*/
typedef struct evloc {
	double lat, lon, z;
	int evid;
} EVENTlocation;
/* This function loads up the full event catalog of EVENTlocation structures
for the database view pointed to by db.  Note this should be a database
view formed by event->origin subsetted with orid==prefor. 
It returns a Tbl loaded with these pointers.  
*/
Tbl *load_full_catalog(Dbptr db)
{
	Tbl *t;
	int nrecords;
	int i;
	EVENTlocation *e;

	t = newtbl(0);
	dbquery(db,dbRECORD_COUNT,&nrecords);
	for(db.record=0;db.record<nrecords;++db.record)
	{
		allot(EVENTlocation *,e,1);
		dbgetv(db,0,"origin.lat",&(e->lat),
			"origin.lon",&(e->lon),
			"origin.depth",&(e->z),
			"evid",&(e->evid),0);
		e->lat = rad(e->lat);
		e->lon = rad(e->lon);
		pushtbl(t,e);
	}
	return(t);
}
/* Scans a list of EVENTlocation structures stored in all
and builds a list of pointers to EVENTlocation structures
that have a distance < radius and with zmin<=z<zmax.
Depths are assumed to be in km and lat, lon values are 
all assumed to ahve been converted to radians.  
Function returns a fresh list that should be freed when
finished.  freed as freetbl(t,0), that is -- do not call 
free on the pointers */

Tbl *find_close_events(Tbl *all,double lat0,double lon0,
	double zmin, double zmax, double radius)
{
	EVENTlocation *e;
	int i;
	double delta, azimuth;
	Tbl *t;

	t = newtbl(0);

	for(i=0;i<maxtbl(all);++i)
	{
		e = (EVENTlocation *)gettbl(all,i);
		/* it is wise to bypass the distance test if the
		depth test is not satisfied */
		if( ((e->z)>zmax) || ((e->z)<zmin) )continue;

		dist(e->lat,e->lon,lat0,lon0,&delta,&azimuth);
		delta = deg(delta);
		if(deg2km(delta)<radius)pushtbl(t,e);
	}
	return(t);
}


/* This program is a companion to dbpmel.  It reads a gclgrid file and an input
database and creates an css3.0 extension table called cluster that links
events in the database to each grid point.  The recipe used to do this is 
driven by several input parameters.  The basic approach is that a cylindrical
volume is defined around each grid point (well, slightly distorted cyclinder in
a spherical earth).  Such a cylinder is described by it's radius and depth
range.  The depth range is controlled by the input parameter file in 
a straightforward way.  The distance range is dynamic and controlled by
four parameters:  minimum_radius, maximum_radius, radius_step_size, 
and minumum_event_count.  The basic algorithm is that the radius is 
expanded in radius_step_size increments until it either reaches the
maximum_radius or the count exceeds the minimum_event_count.  No entries
in cluster are produced for events for which the count drops below this
minimum at the maximum_radius.  


usage:  cluster db [-pf pffile]

Author:  Gary Pavlis
Written:  July 2001
Revised:  Revision 1.4 changed the algorithm because of a serious 
performance problem in earlier versions.  Initial runs on database
of the order of 8000 events took processing times of the order
of 5 days (estimated, never let one go to completion).  That version
used endless calls to dbsubset using the distance calculator 
that is part of datascope's expression calculator.  To speed the
code the ENTIRE catalog is now loaded into memory the distance
calculation is all done with native binary quantities.  I 
expect at least an order of magnitude increase in speed.
*/

void main(int argc, char **argv)
{
	char *dbin;  /* Input db name */
	Dbptr db,dbv,dbc,dbh;
	Tbl *proctbl;  /* operations passed to dbprocess */
	int i,j,k,ie;
	char *pfin=NULL;
	Pf *pf;

	char *gridname;
	GCL3Dgrid *grd;
	int gridid;
	/* search control parameters */
	double rmin, rmax, dr,dz;
	int minimum_events;
	double gridlat,gridlon,gridr,gridz;
	double zmin,zmax;
	double search_radius,search_radius_km;
	int nrecs;
	int evid;
	double hypocen_lat,hypocen_lon,hypocen_z;
	
	EVENTlocation *e;
	Tbl *allevents,*keepers;

	elog_init(argc,argv);
	elog_notify (0, "$Revision: 1.5 $ $Date: 2001/08/05 08:34:06 $") ;
	if(argc<2) usage();

	dbin = argv[1];

	for(i=2;i<argc;++i)
	{
                if(!strcmp(argv[i],"-pf"))
                {
                        ++i;
                        if(i>=argc) usage();
                        pfin = argv[i];
                }
		else
			usage();
	}
	if(pfin==NULL) pfin=strdup("cluster");
	if(pfread(pfin,&pf)) elog_die(0,"pfread error for pffile=%s\n",pfin);
	rmin=pfget_double(pf,"minimum_radius");
	rmax=pfget_double(pf,"maximum_radius");
	dr = pfget_double(pf,"radius_step_size");
	minimum_events = pfget_int(pf,"minimum_event_count");
	dz = pfget_double(pf,"depth_range");
	dz /= 2.0;

        if(dbopen(dbin,"r+",&db) == dbINVALID)
                die(1,"Unable to open input database %s\n",dbin);

	/* load the gcl grid from the database by name. */
	gridname = pfget_string(pf,"GCLgrid_name");
	if(gridname == NULL) 
		elog_die(0,"Missing required parameter GCLgrid_name\n");
	grd = GCL3Dgrid_load_db(db,gridname);
	if(!(grd->geographic_defined))
		elog_die(0,"Must have geographic components of GCLgrid object with name %s defined for %s to run\n",gridname,Program_Name);

	if(grd == NULL) elog_die(0,"Problems in GCL3Dgrid_load_db\n");
	proctbl = strtbl("dbopen event","dbjoin origin","dbsubset orid==prefor",0);
	dbv = dbprocess(db,proctbl,0);
	dbquery(dbv,dbRECORD_COUNT,&nrecs);
	if(nrecs<=0)
		elog_die(0,"Working view is empty\nCheck event->origin join\n");
	else
		elog_log(0,"Working view has %d events\n",nrecs);
	dbc = dblookup(db,0,"cluster",0,0);
	dbh = dblookup(db,0,"hypocentroid",0,0);

	/* This loads up the full catalog  of events */
	allevents = load_full_catalog(dbv);

	/*3d looping */
	elog_log(0,"Grid point hit counts (lat, long, count)\n");
	for(i=0,gridid=0;i<(grd->n1);++i)
	    for(j=0;j<(grd->n2);++j)
		for(k=0;k<(grd->n3);++k)
		{
		    ++gridid;
		    search_radius_km = rmin;
		    gridz = r_to_depth(grd->r[i][j][k],grd->lat[i][j][k]);
		    zmin = gridz - dz;
		    zmax = gridz + dz;
		    while(search_radius_km<=rmax)
		    {
			search_radius = km2deg(search_radius_km);
			gridz = r_to_depth(grd->r[i][j][k],
					grd->lat[i][j][k]);
			/*note gclgrid stores lat/lon in radians
			while db routines assume degrees */
			/* somewhat arbitrary ceiling */
			if(zmin<-10.0) zmin=-10.0;
			if(zmax < zmin)elog_die(0,"Grid setup problem:  depth floor computed as %lf km, which is above all earth's surface\n",zmax);

			keepers = find_close_events(allevents,
					grd->lat[i][j][k],grd->lon[i][j][k],
					zmin,zmax,search_radius_km);

			nrecs = maxtbl(keepers);
			if(nrecs>=minimum_events) break;

			freetbl(keepers,0);
			search_radius_km += dr;
		    }
		    if(nrecs>=minimum_events) 
		    {
			hypocen_lat = 0.0;
			hypocen_lon = 0.0;
			hypocen_z = 0.0;
			elog_log(0,"%lf %lf %d\n",deg(grd->lat[i][j][k]),
			deg(grd->lon[i][j][k]),nrecs);
			for(ie=0;ie<nrecs;++ie)
			{
			/* We have to be careful about crossing
			the equator or the prime meridian*/
			int lon_is_positive;
			double lat,lon,z;

			e = (EVENTlocation *)gettbl(keepers,ie);
			lat = deg(e->lat);
			lon = deg(e->lon);
			z = (e->z);
			if(ie==0)
			{
				if(lon>0.0)
					lon_is_positive=1;
				else
					lon_is_positive=0;
			}
			if(lon_is_positive)
			{
				if(lon>0.0)
					hypocen_lon += lon;
				else
					hypocen_lon += (360.0+lon);
			}
			else
			{
				if(lon<0.0)
					hypocen_lon += lon;
				else
					hypocen_lon += (lon-360.0);
			}
			hypocen_lat += lat;
			hypocen_z += z;
			if(dbaddv(dbc,0,"gridname",gridname,
				"gridid",gridid,
				"evid",e->evid,
					0)<0) elog_die(0,"Error appending to cluster table for gridid %d and evid %d\n",
					gridid,evid);
			}
			hypocen_lat /= (double)nrecs;
			hypocen_lon /= (double)nrecs;
			hypocen_z /= (double)nrecs;
		    }
		    else
		    {
			hypocen_lat=deg(grd->lat[i][j][k]);
			hypocen_lon=deg(grd->lon[i][j][k]);
			hypocen_z=gridz;
			nrecs=0;
		    }
		    if(dbaddv(dbh,0,"gridname",gridname,
			"gridid",gridid,
			"dlat",deg(grd->lat[i][j][k]),
			"dlon",deg(grd->lon[i][j][k]),
			"depth",gridz,
			"hclat",hypocen_lat,
			"hclon",hypocen_lon,
			"hcdepth",hypocen_z,
			"nass",nrecs,
			"delta",search_radius,
			"ztop",zmin,
			"zbot",zmax,0) < 0) elog_die(0,
			    "Error updating hypocentroid table for grid point %d,%d,%d\n",
					i,j,k);
		}
	freetbl(allevents,free);
}
