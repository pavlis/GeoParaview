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

#define KMPERDEG 111.19

Hypocenter db_load_initial(Dbptr dbv,int row)
{
	Hypocenter h;
	dbv.record = row;
	if(dbgetv(dbv, 0, 
		"origin.lat", &(h.lat),
		"origin.lon", &(h.lon),
		"origin.depth",&(h.z),
		"origin.time", &(h.time),
		0) == dbINVALID)
			die(1,"relocate:  dbgetv error fetching previous location data\nFailure at line %d of database view\n",row);
	/* This initializes parts of the hypocenter stucture that define
        this as an initial location. */
        h.dz = 0.0;
        h.dx = 0.0;
        h.dy = 0.0;
        h.dt = 0.0;
        h.rms_raw = -1.0;
        h.rms_weighted = -1.0;
        h.interquartile = -1.0;
        h.number_data = 0;
        h.degrees_of_freedom = 0;
	h.lat0 = h.lat;
	h.lon0=h.lon;
	h.z0 = h.z;
	h.t0 = h.time;
	return(h);
}

/* This is a group of procedures to update a set of css3.0 
database tables that are changed in a relocation.  Some information
is copied, and some is updated.  It is kind of ugly, but a necessary
evil with this particular schema. */

/* save_origin adds a single record for each event to the origin
table and returns the orid of this event returned by dbnextid
for the output table.  

Arguments:

dbi - input database pointer of complex view used internally in
	this program 
is, ie - range of rows of dbi containing data for this event
depth_fixed - set true when fixed depth solution is used 
hypos - hypocenter object for relocated solution
dbo - output db pointer.  

Function returns the orid assigned to this event from dbnextid

Author:  Gary L. Pavlis
Written:  February 1997
*/

int save_origin(Dbptr dbi, int is, int ie, int depth_fixed,
	Hypocenter h, Dbptr dbo)
{
	int orid;

	/* These are parameters copied from input db -- names = css names */
	int evid;
	int jdate;
	int grn;
	int srn;
	char etype[8];
	double mb;
	int mbid;
	double ms;
	int msid;
	double ml;
	int mlid;
	/* These obtained from this solution */

	int nass;
	int ndef;
	char dtype[2];
	char algorithm[16]="dbpmel";
	char *auth;
	

	/* these are intentionally left null: ndp,depdp,commid*/


	/* set but obtained directly from hypo structure 
	lat, lon, depth, time */

	dbo = dblookup(dbo,0,"origin",0,0);
	
	/* Grab what we need from the input db for copying.  Note
	that because the joins used here, each of the input db
	records have the same origin entries.  Thus, we just fetch
	stuff from row is.  */
	dbi.record = is;
	if( dbgetv(dbi,0,
		"origin.evid", &evid,
		"origin.jdate", &jdate,
		"origin.grn", &grn,
		"origin.srn", &srn,
		"origin.etype", etype,
		"origin.mb", &mb,
		"origin.mbid", &mbid,
		"origin.ms", &ms,
		"origin.msid", &msid,
		"origin.ml", &ml,
		"origin.mlid", &mlid,
				0) == dbINVALID)
	{
		die(1,"save_origin: dbgetv error reading origin fields of input view at record %d\n",
				is);
	}
	nass = ie - is + 1;
	/* ndef is potentially wrong by this calculation.  We set ndef 
	from the number of degrees of freedom in the solution + npar
	where npar is 3 when depth is fixed and 4 otherwise.  This does
	not match the real definition of ndef, especially when array 
	data are used because each arrival then then has 3 data points
	used to constrain the solution.  These is either an error in 
	my (glp) understanding of the css3.0 schema or a flaw in its 
	design.*/
	if(depth_fixed)
	{
		ndef = h.degrees_of_freedom + 3;
		strcpy(dtype,"r");
	}
	else
	{
		ndef = h.degrees_of_freedom + 4;
		strcpy(dtype,"f");
	}
	auth = cuserid(NULL);
	orid = dbnextid(dbo,"orid");
	if(orid<0) elog_die(0,"dbnextid failure for orid: check write permission on lastid table\n");
	if(dbaddv(dbo,0,
                "lat",h.lat,
                "lon",h.lon,
                "depth",h.z,
                "time",h.time,
                "orid",orid,
                "evid",evid,
                "jdate",jdate,
                "nass",nass,
                "ndef",ndef,
                "grn",grn,
                "srn",srn,
                "etype",etype,
                "dtype",dtype,
                "mb",mb,
                "mbid",mbid,
                "ms",ms,
                "msid",msid,
                "ml",ml,
                "mlid",mlid,
                "algorithm",algorithm,
                "auth",auth,
			0) == dbINVALID)
	{
		die(1,"save_origin: dbaddv error writing orid %d\n",
				orid);
	}
	return(orid);
}
/*This function adds a single row to the origerr table for this event.

arguments:
	orid - orid assign to this solution
	h - Hypocenter object for this solution
	dbo - output db

This version doesn't set much of this large table.  Eventually, it 
needs to include all the error terms, but that function is not 
written yet. 

Author:  Gary L. Pavlis
Written:  February 1997
*/
void save_origerr(int orid, Hypocenter h, double **C, Dbptr dbo)
{
	double sdobs; 
	/* Intentionally ignored: smajax,sminax,strike,sdepth,stime,
					conf,commid */

	dbo = dblookup(dbo,0,"origerr",0,0);

	/* Bad news here is that we can't save the more useful 
	statistics that this program calculates.  css3.0 only allows
	sdobs = sswr/ndgf */

	sdobs = h.rms_raw;
	if(dbaddv(dbo,0,
		"orid",orid,
                "sxx",C[0][0],
                "syy",C[1][1],
                "szz",C[2][2],
                "stt",C[3][3],
                "sxy",C[0][1],
                "sxz",C[0][2],
                "syz",C[1][2],
                "stx",C[0][3],
                "sty",C[1][3],
                "stz",C[2][3],
		"sdobs",sdobs,
			0) == dbINVALID)
	{
		die(1,"save_origerr: dbaddv error writing origerr record for orid %d\n",
				orid);
	}
}

/*This function forms a new assoc table for this solution using 
input view, dbv, from is to ie (rows) as a pattern.  This is 
more complex than it's siblings above because of the need to deal
with residuals.  It gets especially ugly because of the fact that
I insist that slowness vectors use components rather than the
polar form.  

This algorithm is quite different from an earlier version used
in earlier version of genloc.  The need to support S-P times
led me to completely rewrite this function.  I uses a similar
algorithm to the earlier one, but completely different in
detail.  The key difference is that the earlier version worked
from the Tbl of strings returned by ggnloc while this version
uses Residual structures added to the Arrival and Slowness_vector
structures by Dan Quinlan when he used my codes to write
the initial version of dbgenloc.  

arguments:
	dbi - input db view -- event table is same for each row of
		this view.  looks only at the first
	is, ie - range of rows of dbi containing data for this event
	orid - output db orid for this solution
	arrival - Tbl of Arrival structures (res = Residual field
		is used from this only in the S-P case)
	slow - Tbl of Slowness_vector structures similar to arrival.  
	h - Hypocenter object for this solution
	dbo - output db

Author:  GAry Pavlis
Written:  February 1997
Revised:  Sept. 2000 (near total rewrite)
*/
void save_assoc(Dbptr dbi, int is, int ie, int orid, char *vmodel,
	Tbl *arrivals,Tbl *slow,Hypocenter h, Dbptr dbo)
{
	/* These fields are copied from input assoc table */
	int arid; 
	char sta[8];
	char phase[10];
	double belief;
	/* These fields are set here */
	double delta;
	double seaz;
	double esaz;
	double timeres;
	double azres;
	double slores;
	double wgt;
	char timedef[2],slodef[2], azdef[2];
	
	/* intentionally ignored:  emares, commid */


	/* passed through arg list;  orid*/

	/* We use this to produce a keyed arr list of the residual 
	list passed into here as a Tbl */

	Arr *worktarr, *workuarr; 
	int i;
	char key[40]; 
	Arrival *aptr;
	Slowness_vector *uptr;

	double uxresid, uyresid,uxpred,uypred,slowres;
	double stalat, stalon; 
	double duphi;
	int irec;

	dbo = dblookup(dbo,0,"assoc",0,0);

	/* Build an associative array keyed by arid:phase for
	each entry in the arrival and slowness vector Tbls.  arid is
	set with a simple %d string with leading zeros.  This may
	be less efficient than writing a compare function for ints, but
	this is safer and the array indexed will not be that huge anyway*/
	worktarr = newarr(0);
	workuarr = newarr(0);
	for(i=0;i<maxtbl(arrivals);++i)
	{
		char *phase1,*phase2;
		aptr = (Arrival *)gettbl(arrivals,i);
		phase1 = strdup(aptr->phase->name);
		phase2=strchr(phase1,'-');
		if(phase2)
		{
			/* phase names with a - will land here.  We assume
			these are things like S-P.  We duplicate the
			aptr pointers for the two phases.  These 
			are checked later for consistency with original. */
			*phase2='\0';
			++phase2;
			sprintf(key,"%08.8d:%s",aptr->arid,phase1);
			setarr(worktarr,key,aptr);
			if((aptr->arid2)>=0)
			{
				sprintf(key,"%08.8d:%s",aptr->arid2,phase2);
				setarr(worktarr,key,aptr);
			}
			else
			{
				elog_complain(0,"WARNING (save_assoc): minus phase inconsistency.\nParsed %s to be %s minus %s\nArid of %s is the illegal value %d and will be skipped\n",
				aptr->phase->name,phase1,phase2,phase2,
				aptr->arid2);
			}	
		}
		else
		{
			sprintf(key,"%08.8d",aptr->arid);
			setarr(worktarr,key,aptr);
		}
	}
	/* Similar for slowness, but we don't need to fret about phases
	like S-P in here since the concept makes no sense for a slowness
	vector measurements*/
	for(i=0;i<maxtbl(slow);++i)
	{
		uptr = (Slowness_vector *)gettbl(slow,i);
		sprintf(key,"%08.8d",uptr->arid);
		setarr(workuarr,key,uptr);
	}
	/* We now loop through the input db view usingit as a
	pattern for the output db.*/
	for(dbi.record=is;dbi.record < ie;++dbi.record)
	{
		/*We need to set the null values of these variables at top
		of this loop to be sure the db rows are written correctly.
		We do this with an undocumented trick. Setting record to
		dbNULL and calling dbgetv sets the values to their null
		field values */
		dbo.record = dbNULL;
		dbgetv(dbo,0,"timeres",&timeres,
			"azres",&azres,
			"slores",&slowres,0);
		wgt = 0.0;
		strcpy(timedef,"n");
		strcpy(slodef,"n");
		strcpy(azdef,"n");
		if( dbgetv(dbi,0,
      	    		"assoc.arid",&arid,
          		"assoc.sta",sta,
          		"assoc.phase",phase,
          		"assoc.belief",&belief,
			"site.lot",&stalat,
			"site.lon",&stalon,
				0) == dbINVALID)
		{
		  register_error(0,"save_assoc: dbgetv error reading assoc fields of input view at record %d\nRows of assoc will not be copied\n",
				dbi.record);
		  continue;
		}

		dist(rad(h.lat),rad(h.lon),rad(stalat),rad(stalon),
				&delta,&esaz);
		dist(rad(stalat),rad(stalon),rad(h.lat),rad(h.lon),
				&delta,&seaz);
		delta = deg(delta);
		seaz = deg(seaz);
		esaz = deg(esaz);

		sprintf(key,"%08.8d",arid);
		aptr = (Arrival *)getarr(worktarr,key);
		uptr = (Slowness_vector *)getarr(workuarr,key);
		if((aptr==NULL) && (uptr==NULL))
		{
			register_error(0,"arid mismatch:  cannot find arid %d expected to associate with orid %d in internal data list.\naOutput assoc table may be incomplete\n",
				arid,orid);
			continue;
		}
		if(aptr==NULL)
		{
			timeres = 0.0;
			wgt = 0.0;
			irec = -1;
		}
		else
		{
			timeres = aptr->res.raw_residual;
			wgt = (aptr->res.weighted_residual)/timeres;
			irec = dbaddv(dbo,0,"arid",arid,
                        	"orid",orid,
                        	"sta",sta,
                        	"phase",phase,
                        	"belief",belief,
                        	"delta",delta,
                        	"seaz",seaz,
                        	"esaz",esaz,
                        	"timeres",timeres,
                        	"timedef",timedef,
                                "wgt",wgt,
                        	"vmodel",vmodel,0);
			if(irec==dbINVALID)
			{
				complain(0,"dbaddv error on assoc table\nCannot add record %d\n",
					dbi.record);
				continue;
			}	
		}
		/* Slowness is complicated by the fact that css3.0 stores
		the slowness vector in polar form and in s/deg.*/
		if(uptr!=NULL)
		{
			uxresid = uptr->xres.raw_residual;
			uyresid = uptr->yres.raw_residual;
			/* We compute the theoretical slowness vector
			by using umeasured and residuals then compute
			the azimuth error from the dot product*/
			uxpred = (uptr->ux)-uxresid;
			uypred = (uptr->uy)-uyresid;
			duphi = acos((uxpred*(uptr->ux)+uypred*(uptr->uy))
			    /(hypot(uptr->ux,uptr->uy)*hypot(uxpred,uypred)));
			azres=deg(duphi);
			slores = hypot(uxresid,uyresid);
			slores *= KMPERDEG;
			strcpy(azdef,"d");
			strcpy(slodef,"d");
			if(irec<0)
			{
				wgt = (uptr->xres.residual_weight)
					*(uptr->xres.other_weights);
				dbaddv(dbo,0,
                       			"arid",arid,
                        		"orid",orid,
                        		"sta",sta,
                        		"phase",phase,
                        		"belief",belief,
                        		"delta",delta,
                        		"seaz",seaz,
                        		"esaz",esaz,
                       			"timedef",timedef,
                        		"azres",azres,
                       		 	"azdef",azdef,
                        		"slores",slores,
                        		"slodef",slodef,
                        		"wgt",wgt,
                        		"vmodel",vmodel,0);
			}
			else
			{
				dbo.record = irec;
				dbputv(dbo,0,"azres",azres,
                       		 	"azdef",azdef,
                        		"slores",slores,
                        		"slodef",slodef,0);
			}
		}
	}
	freearr(worktarr,0);
	freearr(workuarr,0);
}
/*
This function saves new origin estimates for a group of nevents
hypocenters in the origin table.  It does this by matching evids,
grabbing all auxiliary fields from the origin table for the prefor,
and writing these back with updated hypocenter information.
The matching is done by a dumb linear search of the input evid
array.  This was intentional as the assumption was that it was
faster to do a linear search on a vector of ints than the more
general approach of using dbmatches against the input database.
Since the evid list will always be far smaller than the full
catalog of events this assumption should always be appropriate.
If this proves to be a bottleneck, the solution should be to 
use a faster indexing method on the evid list using an Arr
or some other binary tree based index.  I judged the complication
of adding such an indexing was not justified as I expect nevents
will normally be pretty small.  Further, this function is only
called once for each group of events located simultaneously.
In the time it took to write this I probably could done it right
instead of making excuses.

Arguments
        dbv - database to read and write to.  HOWEVER base dbv given
		MUST be a pointer to the basic working view that
		drives pmel_process.   
	reclist - Tbl of records of group pointers for each seismic
		event of grouping defined in input view dbv.  This
		must be consistent with the view pointer passed as 
		db or chaos will clearly result.  
        nevents - number of events in this group
		(a consistency check is made against the size of
		the reclist Tbl, if they aren't the same it hints
		at big trouble.  This maybe should be a die, but 
		it make it only a complain for now.)
        h - vector of Hypocenter structures containing new
                hypocenter estimates to be added to the database.
        evid - parallel vector of ints that are the css3.0 evid
                ids with each element of h.  That is, evid[i] is
                the event id for the hypocenter stored in h[i].
	ta - array of Tbl pointers of length nevents that contain
		arrival data for each h and evid.  That is, ta is
		a parallel array to h and evid of points to Tbls.
	fixlist - dbpmel's Arr that contains indexed list of
		fixed coordinate information. 
	vmodel - velocity model named (required for assoc table) 
	o - control structure passed to ggnloc - assumed consistent
		with initial calls to build covariance matrix.

Author:  Gary Pavlis
*/
int dbpmel_save_hypos(Dbptr dbv, 
		Tbl *reclist,
		int nevents, 
		Hypocenter *h,
		int *evid,
		Tbl **ta,
		Arr *fixlist,
		Pf *pf)
{
	int i,j,match,is,ie;
	int nrec;
	int current_evid;
	int orid;
	Dbptr dbbdl;
	char *fix;
	Tbl *tu=newtbl(0);  /* Null table always for pmel of slowness data*/
	double **C;
	float emodel[4];
	char *vmodel;
	Location_options o;

	o = parse_options_pf(pf);
	vmodel = pfget_string(pf,"velocity_model_name");
	C = dmatrix(0,3,0,3);

	nrec = maxtbl(reclist);
	if(nevents != nrec)
		elog_complain(0,"dbpmel_save_hypos (WARNING): mismatch in event count\nDatabase defines %d events in working group while process list is %d\nMore errors will probably follow this one\n",
			nevents,nrec);
	for(i=0;i<nrec;++i)
	{
		int fixdepth;
		dbv.record = (int)gettbl(reclist,i);
		dbgetv(dbv,0,"evid",&current_evid,"bundle",&dbbdl,0);
		/* this is the simple linear search discussed above*/
		for(j=0,match=-1;j<nevents;++j)
		{
			if(evid[j]==current_evid)
			{
				match = j;
				break;
			}
		}
		if(match<0)
		{
			elog_notify(0,"Cannot find evid %d in event process list;  cannot save hypocenter for this event\n",
				current_evid);
		}
		else
		{
			dbget_range(dbbdl,&is,&ie);
			fix = get_fixlist(fixlist,current_evid);
			if(fix==NULL)
				fixdepth=0;
			else
			{
				if(strchr(fix,'z')==NULL)
					fixdepth=0;
				else
					fixdepth=1;
			}
			orid = save_origin(dbv,is,ie,fixdepth,h[match],dbv);
			save_assoc(dbv,is,ie,orid,vmodel,
				ta[match],tu,h[match],dbv);
			predicted_errors(h[match],ta[match],tu,o,C,emodel);
			/* This routine should ultimately be more 
			sophisticated to include error estimates defined
			in pavlis(1991).  This is a placeholder that puts
			some basic stuff in origerr. */
			save_origerr(orid,h[match],C,dbv);
		
		}
	}
	free_matrix((char **)C,0,3,0);
	freetbl(tu,free);
        return(0);
}
void dbpmel_save_sc(int gridid, Dbptr db,SCMatrix *s,Pf *pf)
{
	char *pmelrun,*gridname;
	int *iptr;
	int phacol;
	Tbl *phaskeys,stakeys;
	int i,j,icol;
	char *phase,*sta;
	
	/* we assume these have been verified to be in pf through 
	check_required_pf*/
	pmelrun=pfget_string(pf,"pmel_run_name");
	gridname = pfget_string(pf,"gridname");
	
        phskeys = keysarr(s->phase_index);
        stakeys = keysarr(s->sta_index);
        db = dblookup(db,0,"gridscor",0,0);

        for(j=0;j<maxtbl(phskeys);++j)
        {
        	phase = (char *)gettbl(phskeys,j);
		iptr = (int *)getarr(s->phase_index,phase);
		phacol = *iptr;
		for(i=0;i<maxtbl(stakeys);++i)
		{
			sta = (char *)gettbl(stakeys,i);
			iptr = (int *)getarr(s->sta_index,sta);
			icol = phacol + (*iptr);
			/* This trap is probably superfluous, but better to
			be paranoid */
			if((icol>=0) && (icol<s->ncol))
				dbaddv(db,0,"sta",sta,
					"phase",phase,
					"gridname",gridname,
					"gridid",gridid,
					"pmelrun",pmelrun,
					"tsc",s->sc[icol],
					"tscref",s->scref[icol],
					"tscbias",s->scbias[icol],0);
			}
			else
				elog_die(0,"dbpmel_save_sc:  computed index %d for accessing station correction vector is outside range of 1 to %d\nProbably access violation making continuation ill advised\n",
					icol+1,s->ncol);
		}
	}
	freetbl(phskeys,0);
	freetbl(stakeys,0);
}
			
	
