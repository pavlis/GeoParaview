#include <math.h>
#include "stock.h"
#include "arrays.h"
#include "db.h"
#include "pf.h"
#include "elog.h"
#include "dbpmel.h"

/* Simple little function to read the 3D reference model pf 
like the normal phase description in genloc, but nested inside
of the label  
3Dphase_definitions {  all the phase description stuff }
*/
Arr *parse_3D_phase(Pf *pf)
{
	Pf *pf3d;
	Arr *a;

	if(pfget(pf,"3Dphase_definitions",&pf3d) != PFSTRING)
		die(0,"pf error reading 3Dphase_definition block\n");
	a = parse_phase_parameter_file(pf3d);
	return(a);
}

initialize_hypocenter(Hypocenter *h)
{
	h->dx=0.0;
	h->dy=0.0;
	h->dz=0.0;
	h->dt=0.0;
	h->lat=0.0;
	h->lon=0.0;
	h->z = 0.0;
	h->lat0=0.0;
	h->lon0=0.0;
	h->z0 = 0.0;
	h->rms_raw=0.0;
	h->rms_weighted=0.0;
	h->interquartile = 0.0;
	h->number_data = 0;
	h->degrees_of_freedom = 0;
}
/* This small routine loads the hypocentroid from the view pointer.
It will fail if the hypocentroid extension table has not been joined
to the working view.

Arguments:
	dbv - working view pointer
	rec - record to read hypocentroid coordinates from
	h - Hypocenter structure in which hypocentroid is returned.
		(Only space coordinate components are set)

Returns:
0 - aok
-1 - dbgetv error
*/	
int load_hypocentroid(Dbptr dbv,int rec, Hypocenter *h)
{
	double lat,lon,depth;
	dbv.record = rec;
	if(dbgetv(dbv,0,"hypocentroid.dlat",&lat,
		"hypocentroid.dlon",&lon,
		"hypocentroid.depth",&depth,0) == dbINVALID)
	{
		elog_notify(0,"dbgetv error reading hypocentroid coordinates from row %d of working view\n",
			rec);
		h->lat = 0.0;
		h->lon = 0.0;
		h->z = 0.0;
		return(-1);
	}
	h->lat = lat;
	h->lon = lon;
	h->z = depth;
	return(0);
}
#define EVIDGRP "evidgroup"
/* this is the main processing routine for dbpmel.  It takes an input
list of grid point, which are defined by integers gridid that are
passed through gridid_list, and processes data associated with each
of these grid points in sequence.   It uses datascope db manipulations
to form a view of arrivals associated with the current grid point 
and runs pmel against arrivals for that group of events.  Detailed
behaviour comes through a long list of options defined in the 
parameter file.  

Arguments:
	db - input db view of formed by a complex sequence (see below).
	gridid - input Tbl list of gridids stored as pure
		integers in the Tbl.   i.e. they are extracted
		as gridid = (int)gettbl(...)
	pf - parameter space of options.  It is important that
		this pf object be checked with a routine that 
		verifies if key parameters are set using the function
		check_required_pf.  This is necessary because pf
		is accessed deep in this code to set check some
		of the many complex options.  This stops the program
		from aborting at inconvenient points within a 
		processing loop.

The view defined by db is very complicated.  This comment could go
out of data, but as of the date of initial writing the view was
formed by running dbprocess with this sequence:
        dbopen cluster
        dbjoin event
        dbjoin origin
        dbsubset orid==prefor
        dbjoin assoc
        dbjoin arrival
        dbjoin site
	dbjoin hypocentroid
        dbsort gridid evid

This is followed by a join with site using ondate:offdate keys to 
make the join work properly in varying station geometry.  A subset
by a sift key is also applied, but this should only reduce data not
alter the view.  The view is next sorted by the key
gridid:evid:sta:phase with a dbUNIQUE to prevent redundant picks on
multiple channels. 

Author:  Gary Pavlis
Written:  Fall 2000
*/
int dbpmel_process(Dbptr db, Tbl *gridid_list,Pf *pf)
{
	Location_options o;
	Arr *arr_phase;
	char *vmodel;
	/* Holds indexed list of stations with bad clocks problems 
	that should be used only with S-P type timing */
	Arr *badclocks;
	/* Hold station table.   We make the array table empty always. */
	Arr *stations,*arrays=newarr(0);
	int nbcs;
	int i,j,k;
	Tbl *grptbl;
	Dbptr dbevid_grp;  /* dbgroup pointers */
	/* Used by dbmatches */
	static Hook *hook=NULL;
	Tbl *reclist;
	Dbptr dbgs;
	Dbptr dbbundle;
	Dbptr dbcs;
	int evid;
	Tbl **ta;  /* We use an array of Tbls for arrival times to 
			allow a simple loop through an event group.*/
	Hypocenter *h0;  
	int *evid;  /* parallel array to h0 of event ids for each h0[i]*/
	Hypocenter hypocentroid;
	/* The associative array here is loaded from the pf and contains
	information on which events are to be treated as calibration events.
	The fixlist array of character strings is passed to pmel as a 
	parallel of vectors defining events with fixed coordinates */
	Arr *events_to_fix;
	char **fixlist;
	/* This is the structure that encapsulates all the station
	correction elements.  WARNING: it is built up in pieces below,
	and later modifications most handle this carefully.*/
	SCMatrix *smatrix;
	int worksize;
	Tbl *converge;

	initialize_hypocenter(hypocentroid);

	/* DB is now set up correctly, now we turn to the parameter files */
	o = parse_options_pf (pf);
 	arr_phase = parse_phase_parameter_file(pf);
	/* We want to explicitly clear the station corrections in the 
	phase handle list if any were set by the above function.  If not,
	we have a real mess with corrections in corrections.  */
	clear_station_corrections(arr_phase);

	/* Now load the model used to compute bias correction with 
	a (presumed) 3D model.  Not really required 3D, but we call it that
	anyway.*/
	arr_phase_3D = parse_3D_phase(pf);

	/* load the station table and use it to create the station
	indexes for the smatrix structure and setup the result
	vectors stored there.*/
	stations = dbpmel_load_stations(db,pf);
	smatrix = create_SCMatrix(stations,arr_phase);

	/* These routines set up definitions of stations to use 
	S-P type phases with due to "bad clocks".   We use this
	here because we want to either throw away data from stations
	with bad clocks or use the S-P feature.  This sets up
	that feature, but lets the decision of how to handle
	this problem be set in pf.*/
	badclocks = newarr(0);
	if(db_badclock_definitions(db,pf,badclocks))
		elog_notify(0,"Problems in database definitions of bad clocks time periods\n");
	pfget_badclocks(pf,badclocks);
	/* nbcs is used as a boolean below, but the cntarr sets it 
	to a nonzero value if there are any entries in the badclock 
	array.  This is useful for efficiency to bypass the messy
	S-P code unless it is necessary. */
	nbcs = cntarr(badclocks);

	vmodel = pfget_string(pf,"velocity_model_name");

	/* Note pfget_arr returns null if the name is missing.
	To make the result less prone to failure in this case
	we create an empty array in that case an complain.*/
	events_to_fix = pfget_arr(pf,"pmel_calibration_events");
	if(events_to_fix == NULL)
	{
		elog_complain(0,"pmel_calibration_events parameter (Arr) not in parameter space\nAssuming an empty list\nCorrect before the next run of this program\n");
		events_to_fix=newarr(0);
	}
	/* We next create a grouping of the working view by
	the gridid:evid key.  We use dbmatches below to match
	gridids and the record list from dbmatches is the set
	of group pointers for the collection of events matching
	a given gridid.*/
	grptbl = newtbl(2);
	pushtbl(grptbl,"gridid");
	pushtbl(grptbl,"evid");
	dbevid_grp = dbgroup(db,grptbl,EVIDGRP,1);
	if(dbevid_grp.record == dbINVALID)
		die(0,"dbgroup failed on gridid:evid bundling\n");
	dbgs = dblookup(db,0,EVIDGRP,0,0);
	dbgs.record = dbSCRATCH;

	dbcs = dblookup(db,0,"gridstat",0,0);
	
	for(i=0;i<maxtbl(gridlist);++i)
	{
		int gridid;
		int nevents;
		int is,ie;
		int ndata;
		Dbptr dbhypo; 

		gridid = (int)gettbl(gridlist,i);

		dbputv(dbgs,0,"gridid",gridid,0);
		dbmatches(dbgs,dbevid_grp,&grptbl,&grptbl,&hook,
					&reclist);
		nevents = maxtbl(reclist);
		if(nevents<=0)
		{
			fprintf(stdout,"No data for gridid = %d\n",gridid);
			continue;
		}
		allot(Tbl **,ta,nevents);
		allot(Hypocenter *,h0,nevents);
		allot(double *,evid,nevents);
		allot(char **,fixlist,nevents);
		
		/* reclist now contains a collection of record numbers
		for gridid:evid grouped parts of the working view. */
		for(j=0,ndata=0;j<nevents;++j)
		{
			char testevid[12];
			dbevid_grp.record = (int)gettbl(reclist,i);
			dbgetv(dgevid_grp,0,"evid",evid+j,
				"bundle",&dbbundle,0);
			dbget_range(dbbundle,&is,&ie);

			ta[j] = dbload_arrival_table(dbbundle,
                                is,ie,stations, arr_phase);
			if(nbcs)
			{
				if(minus_phase_arrival_edit(ta[j],
					arr_phase,badclocsk))
				{
					elog_notify(0,"Warning (dbpmel_process):  problems in editing arrival table for minus phases in minus_phases_arrival function\n");
				}
			}
			ndata += maxtbl(ta[j]);
			h0[j] = db_load_initial(dbbundle,is);
			/*This sets fixed coordinates on calibration events.
			It depends on getarr returning a NULL pointer if
			the there is no match. */
			sprintf(testevid,"%d",evid[j]);
			fixlist[i] = (char *)getarr(events_to_fix,testevid);
			if(fixlist[i]!=NULL)
				elog_log(0,"Using evid %d as calibration event fixing coordinates as: %s\n",
					evid,fixlist[i]);
			
		}
		/* We assume the hypocentroid has been joined to this
		view and we can just grab the first row of the view pointer
		to get the hypocentroid location for this group. */
		if(load_hypocentroid(dbbundle,is,&hypocentroid))
		{
			elog_complain(0,"Error loading hypocentroid from working view for gridid=%d;  Skipping to next gridid in processing list\n",
				gridid);
			continue;
		}

		/* This function alters the phase handles by setting the
		station corrections to those computed from the difference
		in travel time between the model defined in the arr_phase_3D
		definition and that in arr_phase.  */
		ierr = initialize_station_corrections(arr_phase,arr_phase_3D,
			stations) )
		if(ierr>0)
		{
			elog_notify(0,"%d problems setting path anomaly corrections for gridid=%d\n",
				ierr,gridid);
		}
		else if(ierr<0)
		{
			elog_complain(0,"Cannot compute any path anomaly corrections for gridid=%d\nSkipping to next grid point\n",
				gridid);
			for(k=0;k<nevents;++k) freetbl(ta[i],free);
			free(ta);
			continue;
		}
		/* We make the S work space of size ndata times smatrix->ncol
		to avoid having to could fixed coordinates.  We could work in
		a space as small as ndata-4*nevents if all events had free 
		coordinates,but this would be too small if any coordinates are
		fixed.  Hence, we may waste some memory, but this should not
		be a large factor.  Note the use of realloc which causes this
		workspace to be dynamic from group to group */
		worksize = sizeof(double)*ndata*(smatrix->ncol);
		if((realloc(smatrix->S,worksize))!=worksize)
			die(0,"realloc failure trying to build S matrix of %d bytes\n",
				workspace);
		s->nrow = ndata;
		for(k=0;k<(ndata*(smatrix->ncol));++k) smatrix->S[k]=0.0;

		/* This routine creates an initializes the station
		correction matrix to 0s. note it assumes the nrow and ncol
		components are already set.*/

		/* This computes the set of reference station corrections
		for this group of events */
		if(compute_scref(scmatrix, &hypocentroid,
				arr_phase,arr_phase_3D))
		{
			error call
		}
		/* This is the main processing routine.  It was 
		intentionally built without any db hooks to make it
		more portable */
		converge=pmel(nevents,ta,h0,fixlist,&hypocentroid,&smatrix,pf);
		fprintf(stdout,"Cluster id=%d pmel convergence reason\n",
			gridid);
		fprintf(stdout,"sswr/dgf = %lf, dgf = %d, rawrms = %lf\n",
			smatrix.sswrodgf,smatrix.ndgf,smatrix.rawrms);
		dbpmel_save_hypos(nevents,h0,evid);
		if(dbaddv(dbcs,0,"gridid",gridid,
				"time",now(),
				"sswrodgf",smatrix.sswrodgf,
				"ndgf",smatrix.ndgf,
				"rawrms",smatrix.rawrms,0) == dbINVALID)
		{
			elog_complain(0,"dbaddv error for gridid %d adding to gridstat table\n",
				gridid);
		}
		for(k=0;k<nevents;++k) freetbl(ta[i],free);
		free(ta);
		free(fixlist);
	}
	freearr(events_to_fix,0);
	destroy_SCMatrix(smatrix);
	return(0);
}
