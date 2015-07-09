#include <stdlib.h>
#include <values.h>

#include "stock.h"
#include "coords.h"
#include "arrays.h"
#include "db.h"
#include "elog.h"
#include "pf.h"
#include "location.h"
#include "gclgrid.h"
#include "glputil.h"
void usage()
{
	banner(Program_Name, "$Revision: 1.1 $ $Date: 2001/08/07 21:39:56 $") ;
	elog_die(0,"usage:  %s db [-pf pfname]\n",Program_Name);
}

/* This program is a companion to dbpmel.  It builds a simulated catalog 
constructed from a gclgrid.  I expect to progressively add options to 
this program so don't expect it to be pretty.  

Input control comes through a parameter file and results are all written
to an output database.  


usage:  simcat db [-pf pffile]

Author:  Gary Pavlis
Written:  August 2001
If this program looks a lot like cluster it is because it is.  I basically 
hacked cluster.c apart and built this one from it.
*/

void main(int argc, char **argv)
{
	char *dbin;  /* Input db name */
	Dbptr db,dbc,dbh,dbe,dbo,dba,dbas;
	int i,j,k,ie,ia,ip;
	char *pfin=NULL;
	Pf *pf;

	char *gridname;
	GCL3Dgrid *grd;
	int gridid;
	
	/* holds station table in an associative array */
	Arr *s;
	Tbl *skeys;
	Pf *vpf;
	
	/* model definitions and phase handles to them */
	char *vmodel;
	Arr *arr_phase;
	Tbl *phaselist;  /* phases to actually use */
	Phase_handle *phandle;
	char *phase;

	int nsta,nphases,nused;
	int evid;
	double hypocen_lat,hypocen_lon,hypocen_z;
	int nperg,seed;
	double pdel;
	long int rn;
	double ptest;  /* random number between 0 and 1 */
	double time=0.0,dt;
	Arrival a;
	Hypocenter hypo;
	Travel_Time_Function_Output tto;

	elog_init(argc,argv);
	elog_notify (0, "$Revision: 1.1 $ $Date: 2001/08/07 21:39:56 $") ;
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
	if(pfin==NULL) pfin=strdup("simcat");
	if(pfread(pfin,&pf)) elog_die(0,"pfread error for pffile=%s\n",pfin);
	
	nperg = dbget_int(pf,"number_of_events_per_grid_point");
	pdel = dbget_double(pf,"probability_of_deletion");
	seed = dbget_int(pf,"random_number_seed");
	dt = dbget_double(pf,"event_time_offset");
	
	/*This initializes the random number generator */
	srandom((unsigned int)seed);
	
	/* This is essentially the same function used in pmel to load station array */
	s = load_stations(db,pf);
	skeys=keysarr(s);
	nsta = maxtbl(skeys);
	
	/* this builds genloc phase handles */
	vmodel = pfget_string(pf,"travel_time_model");
        if(pfload("GENLOC_MODELS", "tables/genloc", vmodel, &vpf) != 0)
                elog_die(0,"pfload failed reading working velocity model %s\n",
                        vmodel);
        arr_phase = parse_phase_parameter_file(vpf);
        pffree(vpf);

	phaselist = pfget_tbl(pf,"phases_to_use");
	nphases=maxtbl(phaselist);
	for(i=0;i<nphases;++i)
	{
		phase = (char *)gettbl(phaselist,i);
		if(getarr(arr_phase,phase) == NULL)
			elog_die(0,"Phase %s not defined for model %s\n",
				phase,vmodel);
	}

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

	dbc = dblookup(db,0,"cluster",0,0);
	dbh = dblookup(db,0,"hypocentroid",0,0);
	dbe = dblookup(db,0,"event",0,0);
	dba = dblookup(db,0,"arrival",0,0);
	dbo = dblookup(db,0,"origin",0,0);
	dbas=dblookup(db,0,"assoc",0,0);


	/*3d looping */
	elog_log(0,"Grid point hit counts (lat, long, count)\n");
	for(i=0,gridid=0;i<(grd->n1);++i)
	    for(j=0;j<(grd->n2);++j)
		for(k=0;k<(grd->n3);++k)
		{
		    ++gridid;
		    gridz = r_to_depth(grd->r[i][j][k],grd->lat[i][j][k]);
		    hypocen_lat = deg(grd->lat[i][j][k]);
		    hypocen_lon = deg(grd->lon[i][j][k]);
		    
		    hypo.lat = grd->lat[i][j][k];
		    hypo.lon = grd->lon[i][j][k];
		    hypo.z = gridz;
		    hypo.time = time;

		    if(dbaddv(dbh,0,"gridname",gridname,
			"gridid",gridid,
			"dlat",deg(grd->lat[i][j][k]),
			"dlon",deg(grd->lon[i][j][k]),
			"depth",gridz,
			"hclat",hypocen_lat,
			"hclon",hypocen_lon,
			"hcdepth",hypocen_z,
			"nass",nperg,
			"delta",grd->dx,
			"ztop",gridz-(grd->dr)/2.0,
			"zbot",gridz+(grd->dr)/2.0,0) < 0) elog_die(0,
			    "Error updating hypocentroid table for grid point %d,%d,%d\n",
					i,j,k);

		    for(ie=0;ie<nperg;++ie,time+=dt)
		    {
		        evid = dbnextid(db,"evid");
		        orid = dbnextid(db,"orid");
		        dbaddv(dbe,0,"evid",evid,"prefor",orid,"auth","simcat",0);
		        dbaddv(dbc,0,"gridname",gridname,"gridid",gridid,"evid",evid,0);

		        for(ip=0,nused=0;i<nphases;++ip)
		        {
		            phase = gettbl(phaselist.ip);
		            phandle = (Phase_handle *)getarr(arr_phase,phase);
		            a.phase = phandle;
		    	    for(ia=0;ia<nsta;++ia)
		    	    {
		    		if(pdel>0.0)
		    		{
		    			/* MAXLONG is defined in values.h this gives uniform
		    			distribution between 0 and 1 */
		    			rn=random();
		    			prob = ((double)rn)/((double)MAXLONG);
		    			if(prob<=pdel) continue;
		    		}
		    		s = (Station *)gettbl(skeys,i);
		    		a.sta = s;
		    		
		    		tto = calculate_travel_time(a,hypo,RESIDUALS_ONLY);
		    		/* Since this is a simulator this is a fatal error */
		    		if(tto.time == TIME_INVALID)
		    			elog_die(0,"Travel time calculator failed for station %s at grid point %d\n",
		    				s->name,gridid);
		    		a.time = time + tto.time;
		    		arid = dbnextid(db,"arid");
		    		dbaddv(dba,0,"sta",s->name,
		    			"time",a.time,
		    			"arid",arid,
		    			"chan","BHZ",
		    			"deltim",phandle->deltat0,
		    			"iphase",phase,
		    			"auth","simcat",0);
		    		dbaddv(dbas,0,"arid",arid,
		    			"orid",orid,
		    			"sta",s->name,
		    			"phase",phase,
		    			"timeres",0.0,
		    			"timedef","d",0);
		    				
		    		++nused;
		        }
		        dbaddv(dbo,0,"lat",deg(hypo.lat),
		        	"lon",deg(hypo.lon),
		        	"depth",hypo.z,
		        	"time",hypo.time,
		        	"orid",orid,
		        	"evid",evid,
		        	"nass",nused,
		        	"ndef",nused",
		        	"algorithm","simcat",
		        	"auth","simcat",0);
		   }
}
