#include <string.h>
#include "db.h"
#include "gclcoords.h"
#include "pwmig.h"
#include "multiwavelet.h"

#define PM_MWARRAY 1
#define PM_MWSS 2
#define PM_MODEL_BASED 3
/* This function uses the analytic formula from Aki and Richards for
particle motion emergence angle for a P wave with a particular slowness
emergent at a half space.  The formula is not just Snell's law because
it include the effects of the P to S conversion at the free surface.

Arguments:
	vp0, vs0 - P and S velocities of the half space = surface 
	ux,uy - slowness vector components in standard geographical
		coordinates (x=+east, y=+north)

Returns a unit vector expressed as spherical coordinate angles
theta and phi.

Author:  G Pavlis
Written:  August 2000
*/
Spherical_Coordinate pm_halfspace_model(
		double vp0,
		double vs0,
		double ux,
		double uy)
{
	Spherical_Coordinate s;
	double sin_i,vpvs2,sini2;

	s.radius = 1.0;
	s.phi = atan2(uy,ux);
	sin_i = hypot(ux,uy)*vp0;
	sini2 = sin_i*sin_i;
	vpvs2 = vp0/vs0;
	vpvs2 *= vpvs2;
	s.theta = atan2( 2.0*sin_i*sqrt(vpvs2 - sini2),vpvs2 - 2.0*sini2);
	return(s);
}
/* This is a parallel function to the above but using array averaged
particle motions determined by the multiwavelet array processing
program (mwap).  The function accesses an extension table created
by mwap called mwpm.  It will call die if this table is not found.
The algorithm used pulls three parameters from the pf parameter space:
fc, the center frequency of the wavelets used in mwap; bankid,
a integer identifier for a family of multiwavelet functions,
and an "array name" that tags the array average estimate in mwpm.  These
are all keys in the mwpm table.  We use a cautious algorithm
with fc (floating point matches are always dangerous) were we call
dbmatches on evid and bankid.  We find the row with the closest fc to
request when multiple matches are found.  This may be excessively 
paranoid, but it should always work.

Arguments:
	pf - parameter file object handle.  As noted above two parameters
		are extracted from pf -- fc, bankid, and array_name
	db - db handle used to lookup mwpm table
	evid - event id of current event being processed 

Return:  always returns a Spherical_Coordinate structure with 
radius set to 1.0.  Angles are in geographical coordinates
(x is + east).  Nonfatal errors are flagged by setting the 
radius to a negative number.  Missing parameters in pf, as always,
will cause the program to die.  Other problems are viewed nonfatal
and should be trapped although errors will be posted on elog.  
Normal procedure should be to recover by calling the model based
function above in the event of nonfatal errors.  

Author:  Gary L. Pavlis
Written:  August 2000
*/

Spherical_Coordinate pm_mwarray(Pf *pf,Dbptr db,int evid)
{
	double fc,fctest,fcmatch;
	int bankid;
	char *array_name;
	Tbl *matchkeys;
	Spherical_Coordinate s;
	Dbptr dbscr;
	Hook *hook=NULL;
	Tbl *mtbl;
	int nmatches;

	s.radius = 1.0;
	s.theta = 0.0;
	s.phi = 0.0;

	bankid = pfget_int(pf,"bankid");
	fc = pfget_double(pf,"wavelet_fc");
	array_name = pfget_string(pf,"array_name");

	db = dblookup(db,0,"mwpm",0,0);
	if(db.record == dbINVALID)
	{
		elog_notify(0,"mwpm table lookup failed -- cannot get polarization\n");
		s.radius = -1.0;
		return(s);
	}
	dbscr = dblookup(db,0,"mwpm",0,dbSCRATCH);
	dbputv(dbscr,0,"sta",array_name,"bankid",bankid,"evid",evid,0);
	matchkeys = strtbl("sta","bankid","evid",0);
	nmatches = dbmatches(dbscr,db,matchkeys,matchkeys,&hook,&mtbl);

	if(nmatches == dbINVALID)
	{
		s.radius = -1.0;
		elog_notify(0,"dbmatches failed with evid=%d, bankid=%d, and array name=%s\n",
			evid, bankid,array_name);

	}
	else
	{
		db.record = (int)gettbl(mtbl,0);
		if(dbgetv(db,0,"fc",&fcmatch,
			"majoraz",&majoraz,
			"majorema",&majorema,0) == dbINVALID)
		{
			s.radius = -1.0;
			elog_notify(0,"dbgetv error on mwpm table for evid %d\n",
				evid);
			return(s);
		}
		if(nmatches > 1)
		{
			int i;
			for(i=1;i<maxtbl(mtbl);++i)
			{
				db.record = (int)gettbl(mtbl,i);
				/* We bypass testing dbgetv here assuming
				that if the above didn't generate an error
				this is unlikely to */
				dbgetv(db,0,"fc",&fctest,0);
				if(fabs(fctest-fc)<fabs(fcmatch-fc))
				{
					fcmatch = fctest;
					dbgetv(db,0,"majoraz",&majoraz,
						"majorema",&majorema,0);
				}
			}
		}
		/* We have to change angles to go from azimuth to geographic
		coordinates and convert to radians */
		s.phi = M_PI_2 - rad(majoraz);
		s.theta = rad(majorema);
	}
	return(s);
}
/* This function defines a new set of angles for a ray based coordinate
system for the slowness vector ux, uy.  It does this in a somewhat 
odd way because particle motions can be strongly skewed from the 
theoretical value.  Thus we use angles relative to to a reference
ray direction not a pure model based estimate.  Said another way, 
we compute from a model the difference in emergence angle and 
azimuth for input slowness vectors uxref,uyref and ux,uy.  These
angles are added to the angles in rayref.  The angle manipulations
are done be assure the results are in the correct range for spherical
coordinates (theta 0 to pi/2 and -pi<phi<pi).
Author:  Gary Pavlis
Written:  August 2000
*/
Spherical_Coordinate compute_ray_coordinates(
	Spherical_Coordinate rayref,
	double vp0,
	double vs0,
	double uxref,
	double uyref,
	double ux,
	double uy)
{
	Spherical_Coordinate sref,snew;
	double dphi, dtheta;

	sref = pm_halfspace_model(vp0,vs0,uxref,uyref);
	snew = pm_halfspace_model(vp0,vs0,ux,uy);
	dtheta = snew.theta - sref.theta;
	dphi = snew.phi = sref.phi;
	rayref.theta += dtheta;
	rayref.phi += dphi;
	/* When the zonal angles goes negative, we assume it means we jumped
	across the z axis.  We fix this by adding 180 degrees to the azimuth
	and making the zonal angle positive*/
	if(rayref.theta < 0.0)
	{
		rayref.theta = - rayref.theta;
		rayref.phi += M_PI;
	}
	/* Make sure the azimuth is between -Pi and PI*/
	while(rayref.phi > M_PI) rayref.phi -= M_PI;
	while(rayref.phi < -M_PI) rayref.phi += M_PI;
	return(rayref);
}

int pwstack_process(Dbptr dbgrp, GCL2Dgrid *grd,Pf *pf)
{
	PWMIGraw_trace_gather *g,*raygath;
	int i;
	char *pmode;
	int pmtype=PM_MODEL_BASED;
	double vp0,vs0;
	double *weights;
	int evid, ix1, ix2;
	double *dnorth, *deast;
	double ux,uy;
	double dux,duy; 
	double *moveout;
	int save_r, save_t, save_l;
	Spherical_Coordinate ray0,raynow;
	PWslowness_grid ugrid;
	int i,j;

	g = dbread_decon_data(dbgrp,pf);
	if(g == NULL) return(-1);

	/* the above routine fails if the following dbgetv fails so
	we don't bother to test for dbgetv failure again */
	dbgetv(dbgrp,0,"evid",&evid,"ix1",&ix1,"ix2",&ix2,0);

	/* We now decide how we want to deal with particle motions.
	Several options are allowed that need to have a lengthy 
	documentation ultimately */
	pmmode = pfget_string(pf,"ray_coordinate_definition");
	if(!(strcmp(pmmode,"mwarray")) pmtype = PM_MWARRAY;
	if(!(strcmp(pmmode,"mwss")) pmtype = PM_MWSS;
	vp0 = pfget_double(pf,"P_surface_velocity");
	vs0 = pfget_double(pf,"S_surface_velocity");

	switch(pmmode)
	{
		case (PM_MWSS):
			elog_notify(0,"Single station mode not currently supported\nReverting to mwarray mode\n");
		case (PM_MWARRAY):
			/* This function reads an array average estimate
			of particle motion from multiwavelet program
			mwap extension tables.  */
			ray0 = pm_mwarray(pf,dbgrp,evid);
			/* Errors cause this to fall into default 
			case with current logic. */
			if(ray0.radius>0.0) break;
			elog_notify(0,"Errors in pm_mwarray -- reverting to model based particle motion method\n);
		case (PM_MODEL_BASED):
		default:
			/* default is a model based estimate */
			ray0 = pm_halfspace_model(vp0,vs0,g->ux,g->uy);
	}


	/* Apply front end mutes to all traces */
	mute_gather(g->nsta, g->nsamp, 
		g->x, g->y, g->z, g->si, pf);

	/* compute a set of weights for each station.  This assumes
	blindly that the grid object is consistent with the indexing
	used to locate the given pseudostation point */
	weights = compute_pseudostation_weights(g->nsta,
		grd->lat[ix1][ix2],grd->lon[ix1][ix2],
		g->lat, g->lon, pf);
	/* This function builds the structure that defines the set of
	slowness vectors that are to be used for stacks */
	ugrid = setup_slowness_grid(pf);

	/* We need dnorth, deast vectors to compute moveout sensibly
	for this program */
	allot(double *,dnorth,g->nsta);
	allot(double *,deast,g->nsta);
	geographic_to_dne(g->nsta,grd->lat[ix1][ix2],grd->lon[ix1][ix2],
		g->lat,g->lon,dnorth,deast);

	/* these booleans determine which components are saved 
	for later migration */
	save_r = pfget_boolean(pf,"save_R_component");
	save_t = pfget_boolean(pf,"save_T_component");
	save_l = pfget_boolean(pf,"save_L_component");

	/* This loop is a bit odd because we use integers to define the
	range, but simultaneously use the for loops initialization and
	increment features to set the ordered pairs (ux,uy).  This
	is safer than using using floating point arithmetic to define
	the range of the loops.  */
	for(i=0,ux=ugrid.uxlow;i<ugrid.nux;++i,ux+=ugrid.dux)
	for(j=0,uy=ugrid.uylow;j<ugrid.nuy;++j,uy+=ugrid.duy)
	{
		float *stack;
		/* We want to rotate to ray coordinates for the current
		slowness vector.  However, we use angles relative to that
		defined as ray0 to allow measured particle motions that
		can differ significantly from that predicted theoretically.
		This complexity is irrelevant if the model based ray
		coordinate defintion is used */
		raynow =  compute_ray_coordinates(ray0,vp0,vs0,
						g->ux,g->uy,ux,uy);

		/* Rotate PW gather to ray coordinates.  Currently a
		single transformation matrix is used, but the logic
		would not change if we allowed station dependent
		transformations.   */
		raygath = rotate_pwgather(g,raynow);
		/* The input gather is assumed prealigned with the slowness
		vector defined by raygath->ux,uy.  We use relative moveouts
		from this base moveout for the actual stacks */
		dux = ux - raygath->ux;
		duy = uy - raygath->uy;
		moveout = compute_pmmoveout(raygath->nsta,deast,dnorth,dux,duy);
		
		/*compute ray coordinate transformation for this ux, uy
		as a relative angle from incident ray direction.  This
		function needs to compute an angle relative to rays
		so we can use delta angles only so we can utilize
		base particle motions if we want to. */
		to_be_written()

		if(save_t)
		{
			stack = weighted_stack_with_moveout(xp, raygath->nsta,
				raygath->nsamp, raygath->si,weights, moveout);
			call dbsave routine with "T" to define component 
			free(stack);
		}
		
		if(save_r)
		{
			stack = weighted_stack_with_moveout(yp, raygath->nsta,
				raygath->nsamp, raygath->si,weights, moveout);
			call dbsave routine with "T" to define component 
			free(stack);
		}
		if(save_l)
		{
			stack = weighted_stack_with_moveout(zp, raygath->nsta,
				raygath->nsamp, raygath->si,weights, moveout);
			call dbsave routine with "T" to define component 
			free(stack);
		}


	}

	free(dnorth);
	free(deast);
	free(weights);
	destroy_PWMIGraw_trace_gather(g);
	destroy_PWMIGraw_trace_gather(raygath);
}
	
