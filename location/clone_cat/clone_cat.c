#include <stdlib.h>
#include <values.h>
#include "stock.h"
#include "elog.h"
#include "pf.h"
#include "db.h"
#include "tt.h"
#include "location.h"
#include "glputil.h"
void usage()
{
	elog_die(0,"clone_cat db\n\nNote:  overwrites original db\n");
}
/* This program will clone an original database in place by 
substituting theoretical arrival times for measured arrival 
times.  This effectively makes an idealized data set with the
same arrival constellation as the parent data.  Most of it's 
behaviour is controlled by the input parameter file that 
is frozen to be derived from the name of the program (no options).

Author:  Gary Pavlis
Written:  September 2002
*/

int main(int argc,char **argv)
{
	Dbptr db,dba,dbas;
	Pf *pf;
	int nrec;
	Tbl *kpat,*tpat;
	Hook *mhook=NULL;
	Tbl *matches;
	int seed;
	double default_deltim;
	int use_deltim;
	int add_errors;
	Arr *arr_phase;
	char *vmodel;
	Pf *vpf;
	Phase_handle *p;
	Travel_Time_Function_Output t;
	Ray_Endpoints x;


	elog_init(argc,argv);
	if(argc!=2) usage();

	if(pfread("clone_cat",&pf))elog_die(0,"pfread error\n");
	vmodel=pfget_string(pf,"travel_time_model");
	if(vmodel==NULL)
		elog_die(0,
		"travel_time_model parameter not in\
input pf\n");
        if(pfload("GENLOC_MODELS", "tables/genloc", vmodel, &vpf) !=
0)
                elog_die(0,"pfload failed reading working velocity\
model %s\n",
                        vmodel);
        arr_phase = parse_phase_parameter_file(vpf);
        pffree(vpf);

	add_errors = pfget_boolean(pf,"add_measurement_error");
	if(add_errors)
	{
		use_deltim=pfget_boolean(pf,"use_deltim");
		default_deltim=pfget_double(pf,"default_deltim");
		seed = pfget_int(pf,"random_number_seed");
		srandom((unsigned int)seed);
	}
	if(dbopen(argv[1],"r+",&db)==dbINVALID)
		elog_die(0,"Unable to open input database\n");
	dba = dblookup(db,0,"arrival",0,0);
	dbas = dblookup(db,0,"arrival",0,0);
	dbas.record = dbSCRATCH;
	db = dbform_working_view(db,pf,"workingview");
	dbquery(db,dbRECORD_COUNT,&nrec);
	fprintf(stdout,"Working view has %d records\n",nrec);
	kpat=strtbl("arid",0);
	tpat=strtbl("arid",0);
	for(db.record=0;db.record<nrec;++db.record)
	{
		double rlat,rlon,slat,slon;
		double otime;
		double elev,rz,sz;
		double atime,time_error;
		char phase[10];
		int arid;
		double tt;
		Tbl *tlist=NULL;
		int match_row;
		double deltim;
		char sta[10];
		double *sc;
		long int irand;
		double drand;

		if(dbgetv(db,0,"sta",sta,
			"origin.lat",&slat,
			"origin.lon",&slon,
			"origin.depth",&sz,
			"origin.time",&otime,
			"site.lat",&rlat,
			"site.lon",&rlon,
			"site.elev",&elev,
			"phase",phase,
			"deltim",&deltim,
			"arid",&arid,0)==dbINVALID)
		   elog_die(0,"dbgetv error at row %d\n",nrec);
		x.sta=strdup(sta);
		x.slat = slat;
		x.slon = slon;
		x.sz = sz;
		x.rlat=rlat;
		x.rlon=rlon;
		x.rz=-elev;
		/* for a simulation we want to die if this happens*/
		p = (Phase_handle *)getarr(arr_phase,phase);
		if(p==NULL) 
		{
			elog_complain(0,
"Lookup failed for phase %s on station %s\nSkipped\n",
				phase,sta);
			continue;
		}
		t = p->ttcalc(x,phase,RESIDUALS_ONLY);
		if(t.time==TIME_INVALID)
		{
			elog_complain(0,
			   "travel time calculator failure\n");
			continue;
		}
		atime=otime+t.time;
		/*manually get the station corrections and add them*/
		sc = (double *)getarr(p->time_station_corrections,sta);
		if(sc!=NULL) atime += (*sc);

		if(add_errors)
		{
			if(deltim<=0.0) deltim=default_deltim;
			/* drand becomes uniform (-1,1)*/
			irand = random()-MAXLONG/2;
			drand =2.0*((double)irand)
				/((double)MAXLONG);
			time_error=drand*deltim;
			atime += time_error;
		}
		dbputv(dbas,0,"arid",arid,0);
		dba.record = dbALL;
		dbmatches(dbas,dba,&kpat,&tpat,&mhook,&matches);
		if(maxtbl(matches)!=1) 
			elog_die(0,
			 "dbmatches found %d matches for arid %d\n",
				maxtbl(matches),arid);
		dba.record = (int)gettbl(matches,0);
		dbputv(dba,0,"time",atime,0);
	}
	exit(0);
}
