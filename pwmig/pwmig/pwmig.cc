#include <cstdio>
#include <string>
#include "stock.h"
#include "db.h"
#include "seispp.h"
#include "vmodel.h"

void usage()
{
        cbanner((char *)"$Revision: 1.2 $ $Date: 2003/09/09 11:07:39 $",
                (char *)"db pfi pfo [-V -pf pfname]",
                (char *)"Gary Pavlis",
                (char *)"Indiana University",
                (char *)"pavlis@indiana.edu") ;
        exit(-1);
}

/* An extension table called raygrid is used as a cross reference between the 
generic GCLgrid object definitions and a specific form of GCLgrid used in this plane
wave migration I call a ray grid in this implementation.  A ray grid is the distorted
grid shown in the Poppeliers and Pavlis paper figure with the top surface at the earth
and the depth defined by ray paths.  The raygrid table cross-references with the
gclgdisk table.  This is necessary because the gclgdisk table is indexed by a 
grid name string while a ray grid is parameterized by the slowness vector that
it was built from.  The table is indexed by a generic gridid integer and a name
(area) that defines which of possible multiple grids are to be associated with 
this instance.  This implementation assumes a rectangular slowness grid is used
indexed by two integers ix1 and ix2.  This was generalized as I expect that a
future implementation detail might be to use some other griding scheme for slowness
vectors.  

Arguments:
	db - Antelope database pointer to open database. (does not need to be raygrid table)
	area - area key used to find correct entry (see above)
	ix1, ix2 - integer indexes for desired slowness vector.
	u - slowness vector associated with input data.  This is used for a simple consistency
		check.  In the current implementation the ray grids are created independently
		from the plane wave stacking of the data.  If they do not match within a 
		reasonable tolerance the program will issue an elog message and blunder on.
		This is a necessary test because current implementation uses a pf to describe
		the slowness grid but this is used in two different programs and it would
		be too easy to make them inconsistent.

This function will die on several errors that should only happen if the database is screwed up.
In all cases it is dumb to proceed if this is the case so die is appropriate.

Author:  Gary Pavlis
Written:  May 2003
*/
Hook *matchhook=NULL;  // needs to be global to avoid memory leak or inefficiency of recompute
const double DUTEST=0.001;
string get_gridname(Dbptr db, string area, int iux1, int iux2, Slowness_vector u)
{
	Dbptr dbg,dbgs;
	Tbl *matchkey;
	Tbl *matches;
	char gname[20];
	string sout;
	double ux,uy,dumag;

	dbg=dblookup(db,0,"raygrid",0,0);
	if(dbg.record==dbINVALID)
		die(0,"dblookup of raygrid table failed.  Cannot continue\n");
	dbgs=dblookup(db,0,"raygrid",0,0);
	dbgs.record = dbSCRATCH;
	matchkey=strtbl("area","iux1","iux2",0);
	dbputv(dbgs,0,"area",area,"iux1",iux1,"iux2",iux2,0);
	dbmatches(dbgs,dbg,matchkey,matchkey,&matchhook,&matches);
	nmatches=maxtbl(matches);
	if(nmatches<=0)
		die(0,"No match in raygrid table for index (%d,%d) and area %s\n",
			iux1,iux2,area.c_str());
	if(nmatches>1)
	{
		die(0,"Multiple matches in raygrid table for index (%d,%d) and area %s\n",
			iux1,iux2,area.c_str());
	}
	dbg.record=(int)gettbl(matches,0);
	dbgetv(dbg,0,"gridname",gname,"ux",&ux,"uy",&uy,0);
	dumag=hypot(ux-u.ux,uy-u.uy);
	if(dumag>DUTEST) 
		elog_notify(0,"WARNING (get_gridname):  slowness vector mismatch.\nData have (%lf,%lf) but matching database entry is (%lf,%lf)\n",
			u.ux,u.uy,ux,uy);
	sout=gname;
	return(sout);
}
// fetch event information from global metadata object and return in Hypocenter object
Hypocenter get_event_hypocenter(Three_Component_Ensemble *d)
{
	Hypocenter h;
	try {
		h.lat = d->md.get_double("origin.lat");
		h.lon = d->md.get_double("origin.lon");
		h.z = d->md.get_double("origin.depth");
		h.time = d->md.get_double("origin.time");
	} catch (Metadata_error mde)
	{
		mde.log_error();
		die(0,"Missing required origin information in input data stream\n");
	}
	return(h);
}
// small companion to below blindly computes distance to input index point, (i,j,k), from
// previous point (-1 of each)
double raydist(GCLgrid3d *ray,int i0, int j0, int k0)
{
	double dx1,dx2,dx3;
	dx1 = ray->x1[i0][j0][k0] - ray->x1[i0][j0][k0+1];
	dx2 = ray->x2[i0][j0][k0] - ray->x2[i0][j0][k0+1];
	dx3 = ray->x3[i0][j0][k0] - ray->x3[i0][j0][k0+1];
	return( sqrt( dx1*dx1 + dx2*dx2 + dx3*dx3);
}
/* This is a procedural routine to compute a 3d array of S wave travel times for
ray the ray geometry passed through "raygrid" and using the velocities passed through
the Vs grid.  A 1d referen model, Vs1d, is used as a default.  That is, because all
3d models have a finite extent it is highly likely a point will map outside the domain
of the model.  This procedure automatically defaults to the 1d reference model when 
a point is found outside the grid.  

The routine here is procedural by design to avoid the excessive copying that would
inevitably occur if we started copying the monster GCLgrid objects. That is, the 
OOP version of this would create a new Vs field mapped into the raygrid object coordinate
system and then integrate that to get times.  Here we do the same thing without forming
that intermediate monster but just ask for velocities at each point.  
double ***compute_Swave_gridtime(GCLgrid3d *Vs, Velocity_model *Vs1d, GCLgrid3d *raygrid)
{
	double ***Stime;
	int i,j,k;
	int nfailures=0;
	int n30;
	
	// don't trap memory errors here as the function dies if allocs fail
	Stime = create_3dgrid_contiguous(raygrid->n1, raygrid->n2,raygrid->n3);
	Vs->reset_index();	
	n30=raygrid->n3 - 1;  // top of the grid is the last entry

	for(i=0;i<raygrid->n1;++i)
		for(j=0;j<raygrid->n2;++j)
		{
			double t,dt;
			// First fill the top surface values with zeros as the time reference datum.
			Stime[i][j][n30]=0.0;

			for(k=n30-1,t=0.0;k>=0;--k)
			{
				int lret;
				int index[3];

				lret = Vs->lookup(raygrid->x1[i][j][k],
					raygrid->x2[i][j][k],
					raygrid->x3[i][j][k]);
				if(lret==0)
				{
					Vs->get_index(index);
					vel = Vs->interpolate(raygrid->x1[i][j][k],
                                        	raygrid->x2[i][j][k],
                                        	raygrid->x3[i][j][k]);

				}
				else
				{
				// silently revert to the 1D reference model on any failure.
				// We do count how many times we fail, however, an issue a
				// diagnostic if every call failed.
					++nfailures;
					vel = Vs1d->getv(raygrid->depth(i,j,k));
				}
				dt = raydist(raygrid,i,j,k)/vel;
				Stime[i][j][k]=Stime[i][j][k-1] + dt;
			}
		}
					
}
/* small companion to immediately below.  Creates a new ray path truncated to radius r0.
Throws a simple integer exception if r0 is below outside the range of input ray.  -1 if
above the top, 1 if below.
*/
RayPathSphere *truncate_ray_path(RayPathSphere *ray0, double r0) throw(int)
{
	int nn;  // number of points in truncated ray
	int nr;  // saved to point at truncation point
	int i;
	double frac,dr0,dr;

	// this assumes the input ray is oriented so the points start at the surface
	// and work downward.  This is reversed in the companion program that calls
	// this function
	if(r0>ray0->r[0]) throw -1;
	for(nn=0;nn<ray0->npts;++i)if(ray0->r[nn]<=r0) break;
	if(nn>=ray0->npts) throw 1;
	// need to add one point at the end that is interpolated plus one for 0 start
	nr = nn;  // index points to ray0 position just above r0
	nn += 2;  // length of output vector = +2, 1 for extra point, another because pts versus intervals
	RayPathSphere *rayout = new RayPathSphere(nn);
	for(i=0;i<nn-1;++i)
	{
		rayout->r[i] = ray0->r[i];
		rayout->delta[i] = ray0->delta[i];
		rayout->t[i] = ray0->t[i];
	}
	// handle the last point by just taking the fraction of the total distance
	// spanned by the last point to r0
	rayout->r[nr-1]=r0;
	dr=ray0>r[nr-1]-ray0->r[nr];
	dr0=ray0[nr-1]-r0;
	frac = dr0/dr;
	rayout->delta[nn-1]= (ray0->delta[nn])*frac;
	rayout->t[nn-1]= (ray0->t[nn])*frac;
	return(rayout);
}
	
/* Computes ray path in GCLgrid3d Cartesian reference frame for specified receiver point
and slowness vector.  Ray is traced with a 1d reference model passed through vp.  
Only the path is computed, not the travel times because the path will be integrated 
later using a 3d model.  The path is returned in a 3xn (n is variable) array of 
points that begin at the ix1,ix2,ix3 grid position and terminates at datum.  
That is, the points are an ordered sequence from scatter point toward the surface.  

The algorithm works by using the ray computed in standard coordinates (ray0), truncating
the standard ray to depth of requested point, and then applying coordinate transformations
to put that ray in space in the correct positions.  The later is done in 


Arguments:
	raygrid - pointer to GCLgrid3d object to which this path is to be associated.
	ix1,ix2,ix3 - position in ray grid that output ray is to originate from
		(scatter point)
	u0 - incident ray slowness vector
	ray0 - ray path for incident wavefield in standard coordinates 
		(see object definition)  
Note both raygrid and ray0 are passed as pointers for efficiency


Author:  Gary L. Pavlis
Written:  May 2003
*/
dmatrix compute_P_incident_ray(GCLgrid3d *raygrid,
	int ix1,
		int ix2, 
			int ix3,
				Slowness_vector u0;
					RayPathSphere *ray0) 
{
	//convenient shorthands
	GCLgrid3d& grid=raygrid;
	double r0;  // convenient shorthands
	RayPathSphere *ray;
	double azimuth;
	int i,ii;
	double r_ray;
	
	r0 = grid.depth(ix1,ix2,ix3);
	azimuth=u0.baz();
	
	try {
		ray = truncate_ray_path(ray0,r0);
		dmatrix coords(3,ray->npts);
		// this is a little ugly because we have to reverse the orientation
		// of the ray path.  The output order is opposite of that computed for
		// the RayPathSphere object
		for(i=0,ii=ray->npts-1,r_ray=grid.r[ix1][ix2][ix3];i<ray->npts;++i)
		{
			double lat,lon;
			double deltar;
			Geographic_point this_point;
			latlon(grid.lat[ix1][ix2][ix3],grid.lon[ix1][ix2][ix3],
				ray->delta[ii],azimuth,&lat,&lon);
			this_point = grid.gtoc(lat,lon,r_ray);
			coords(0,i)=this_point.x1;
			coords(1,i)=this_point.x2;
			coords(2,i)=this_point.x3;
			// we need to use deltar values because of ellipticity
			// that is the ray object is in a standard coordinate system for
			// a spherical earth but a GCLgrid uses a standard ellipsoid
			if(ii>1)
			{
				deltar = ray->r[ii-1]-ray->r[ii];
				r_ray -= deltar;
			}
			else
			{
				// this only is executed for the last point.  
				// well, the last two passes but the final one has no effect
				// we avoid a closure error from ellipticity by forcing
				// this point to the surface ellipsoid value
				deltar = r0_ellipse(lat) - ray->r[ii];
//DEBUG
if(ii==0) cerr<<"Closure error="<<deltar<<endl;
			}
		}
	
		delete ray;
	} 
	// either of these errors are  programming error and need to cause the program to die
	catch(int ierr);
	{
		if(ierr<0)
			die(0,"ray trace failed for incident ray:  image point radius %lf is above top surface datum\nProgramming error requires a fix\n",
				r0;
		else if(ierr>0)
			die(0,"ray trace failed for incident ray:  image point radius %lf is below the bottom of the reference ray\nThis is a programming error and caller needs to prevent this\n",
				r0);
	}
}
/* Returns a 3x nx-1  matrix of ray path unit vectors pointing in
direction of tangent to the path defined by the 3 x nx  input matrix 
passed as ray.  The algorithm is a simple forward difference scheme 
with no checking for possible roundoff problems.  Because we use forward
differences the output vector at point i is computed from the segment 
immediately forward (x[i+1]-x[i]).

Note the dmatrix is created with new here and needs to be cleared by caller.

Author:  Gary Pavlis
Written:  August 2003
*/

dmatrix *ray_path_tangent(dmatrix& ray)
{
	dmatrix *gptr;
	dmatrix& gamma=*gptr;
	int *sz,nx;
	double dx[3];
	double nrmdx;
	int i,j;

	sz = ray.size();
	// assume number rows = 3.   could test and throw exception, but this
	// function is not viewed as a library function
	nx=sz[1];
	delete [] sz;
	gptr = new dmatrix(3,nx-1);

	for(i=0;i<nx-1;++i)
	{
		dx[0]=ray(0,i+1)-ray(0,i);
		dx[1]=ray(1,i+1)-ray(1,i);
		dx[2]=ray(2,i+1)-ray(2,i);
		nrmdx=dnrm2(3,dx,1);
		dscal(3,1.0/nrmdx,dx,1);
		for(j=0;j<3;++j) gamma(j,i)=dx[j];
	}
	return(gptr);
}
/* computes depth dependent transformation matrix assuming a scattering
is specular.  That is, transformation matrices will yield a form of
L,R,T 
*/

class Ray_Transformation_Operator
{
public:
	int npoints;
	dmatrix U[];
	Ray_Transformation_Operator(int np)
	{
		U=new dmatrix(3,3)[np];
	}
	// constructor for constant azimuth 
	Ray_Transformation_Operator(GCLgrid&g, dmatrix& path, double azimuth);
	// constructor for depth dependent operator
	Ray_Transformation_Operator(GCLgrid&g, dmatrix& path);
	~Ray_Transformation_Operator(){delete [] U;);
	dmatrix *apply(dmatrix& in);
};
dmatrix *Ray_Transformation_Operator::apply(dmatrix& in)
{
	int i,j;
	double *p;
	int *insize=in.size()
	int nrow=insize[0],ncol=insize[1];
	delete insize;
	// This could be trapped as an exception but since this
	// is expected to only be used in this program this is
	// simpler
	if( (nrow!=3) !! (ncol!=in.npoints) )
	{
		cerr << "Coding error:  Ray_Transformation_Operator has "
			<< in.npoints << "elements but matrix passed is "
			<< nrow << "X" << ncol << endl;
		exit(-1);
	}
	dmatrix *optr=new dmatrix(nrow,ncol);
	dmatrix& out=*optr;

	// my dmatrix class doesn't know about matrix vector multiply
	// so I hand code this using an admittedly opaque approach
	// using pointers.  Algorithm works because dmatrix elements
	// are stored in packed storage ala fortran
	for(i=0;i<ncol;++i)
	{
		double work[3];
		work[0]=in(0,i);
		work[1]=in(1,i);
		work[2]=in(2,i);
		p=U[i].get_address(0,0);
		out(0,i)=p[0]*work[0]+p[1]*work[1]+p[2]*work[2];
		out(1,i)=p[3]*work[0]+p[4]*work[1]+p[5]*work[2];
		out(2,i)=p[6]*work[0]+p[7]*work[1]+p[8]*work[2];
	}
	return(optr);
}
// Constructor for simple cast with all matrices going to surface R,T,L coordinates
Ray_Transformation_Operator::Ray_Transformation_Operator(GCLgrid& g, 
	dmatrix& path,double azimuth)
{
	int *insize=path.size();
	int np=insize[1];
	delete [] insize;
	dmatrix U0(3,3);  // final matrix copied to all elements here
	double x[3];  //work vector to compute unit vectors loaded to U0
	Geographic_point x0_geo;

	//first we need the geographic coordinates of the ray emergence point
	x0_geo = g.ctog(path(0,0),path(1,0),path(2,0));
	
	// this holds transformation matrix from local geo TO the GCLgrid 
	// reference frame
	dmatrix Ugeo=ustrans(g,x0_geo.lat,x0_geo.lon);
	
	U = new dmatrix(3,3)[np];
	// Get the L direction from the first pair of points in path
	x[0]= path(0,0) - path(0,1);
	x[1]= path(1,0) - path(1,1);
	x[2]= path(2,0) - path(2,1);
	
// ran out of time here.  See "rotate" function in rotation.C of libseispp 
// to derive transformation matrix for ray coordinates.  This is currently
// potentially backwards in a couple places.  Should be able to just
// use the matrix computed in the rotate function as this is a conversion
to ray coordinates.  The next function needs the complications I was
working away at here.

		
	
	

/* Computes the domega term for the inverse Radon transform inversion formula.
Computed solid angle subtened the region formed by the difference u-du/2 to u+du/2
in both the x and y directions.  The product is the solid angle subtended.
The geometry is a bit ugly because we have to 

ahhh noooooo 

New algorithm comes to mind.  I'll use the ray trace funtion with u computed as
ux+du and uy+dy.

Going back to main to put that code in.

double compute_solid_angle(Slowness_vector u,double du,double *gradP,double *gradS,
				double svel)

/*  Plane wave migration code main.  This main program does most of the work of the migration code.
Basic algorithm is driven by an  input stream but it also reads a number of fixed entities (notably
the ray grid geometrys and velocity models that are precomputed) from a database.  An antelope pf
is used to control processing parameters for the program.

Required pf quantities:

Strings:  
P_velocity_model3d_name
S_velocity_model3d_name
P_velocity_model_name
ensemble_tag = pfstream required tab surrounding metadata to be parsed for this ensemble
incident_slowness_method = 


Input data stream is assumed organized as input 3C trace ensembles.  An interface routine takes
the input apart and returns a Three_Component_Ensemble object.  
The global metadata object for the ensemble must contain the following:

ux, uy - slowness vector used to stack these data
iux, iuy - index to slowness grid that ux,uy were created from.
origin.lat, origin.lon, origin.depth, origin.time - source location of this event (the ensemble is
	assumed to be a common event gather that is to be migrated)


Each three-component trace must define the following:
ix,iy - index to parent pseusostation grid
elev - elevation that can be used to compute a geometric static

*/

Pfstream_handle *pfshi,*pfsho;
int main(int argc, char **argv)
{
        char *pfin=NULL;
        char *dbname;
        Dbptr db,dbv;
	char *pfstream_in,*pfstream_out;;
	char *Pmodel1d_name;
	char *Pmodel3d_name,*Smodel3d_name;
	char *pvfnm="Pvelocity",*svfnm="Svelocity";
	char *pfsin_tag;  // name attached to input ensemble in pfstream
	Three_Component_Ensemble *pwdata;
	int i,j,k;


        ios::sync_to_stdio();
        elog_init(argc,argv);

        if(argc < 4) usage();
        dbname = argv[1];
	pfstream_in = argv[2];
	pfstream_out = argv[3];

        for(i=2;i<argc;++i)
        {
                if(!strcmp(argv[i],"-V"))
                        usage();
                else if(!strcmp(argv[i],"-pf"))
                {
                        ++i;
                        if(i>=argc) usage();
                        pfin = argv[i];
                }
		else
			usage();
	}
	/* this sets defaults */
        if(pfin == NULL) pfin = strdup("pwmigprecompute");
        if(pfread(pfin,&pf)) die(1,(char *)"pfread error\n");

	dbopen(dbname,'r+',&db);
	if(db.record == dbINVALID) die(1,"Cannot open database %s\n",dbname);
	Pmodel3d_name=pfget_string(pf,"P_velocity_model3d_name");
	Smodel3d_name=pfget_string(pf,"S_velocity_model3d_name");
	Pmodel1d_name=pfget_string(pf,"P_velocity_model1d_name");
	Smodel1d_name=pfget_string(pf,"S_velocity_model1d_name");
	pfsin_tag = pfget_string(pf,"ensemble_tag");
	if(Pfmodel3d_name == NULL || Smodel3d_name==NULL 
		|| Pmodel1d_name==NULL || pfsin_tag==NULL)
		die(0,"Missing required parameters\n");

	try {
		GCLscalarfield3d Vp3d(db,Pmodel3d_name,pvfnm);
		GCLscalarfield3d Vs3d(db,Smodel3d_name,svfnm);
		Velocity_Model Vp1d(db,Pmodel_name,pvfnm);
		Velocity_Model Vs1d(db,Smodel_name,svfnm);

		// launch reader here.  interface under development this uses current version
		pfshi = pfstream_start_read_thread(pfstream_in);
		if(pfshi == NULL) 
			die(1,"Problems in starting read thread from stream %s\n",pfstream_in);
		// need a writer here too to push weighted results to stacker
		pfsho = pfstream_start_write_thread(pfstream_out);
		if(pfsho == NULL) 
			die(1,"Problems in starting write thread to stream %s\n",pfstream_out);

// NOTE  THIS ROUTINE MUST PUSH GLOBAL SOME PARAMETERS TO GLOBAL PF FOR THE ENSEMBLE.  SEE
// EXTRACTIONS BELOW
// ANOTHER NOTE TO MYSELF.  NEED A CONSISTENT WAY TO FLAG GAPS THIS NEEDS TO BE PART OF LIBSEISPP
// THIS LOOP WILL BE PARALLELIZED EVENTUALLY

//YET ANOTHER NOTE.  HAVE TO DECIDE HOW TO HANDLE ZERO DATUM.  DO WE DO A GEOMETRIC STATIC
// AS IN REFLECTION PROCESSING?  PROBABLY SHOULD PUT THIS IN PWSTACK.  THAT IS STACK SHOULD
// SHOULD ALWAYS BE SHIFTED TO DATUM.
// Solution for now.  Modified pwstack to put an elev estimate for the pseudostation
// grid point computed as a weighted sum of elevations of contributing stations.  Algorithm
// below should use a vertical ray approximation to this static to datum or it will get
// unnecessarily complex.  A p=0 approximation for an elevation correction is pretty
// good for teleseismic data and is unlikely to be a big problem.
		while( pwdata = fetch_next_pwensemble(pfshi,pfsin_tag) != NULL)
		{
			Slowness_vector ustack;
			Hypocenter hypo;
			string gridname;
			GCLgrid3d *grdptr;
			GCLgrid3d& raygrid;  // convenient shorthand because we keep changing this
			double ***ts;  // S wave travel time grid
			string area;  // key for database when multiple areas are stored in database
			RayPathSphere *rayupdux,*rayupduy;
			int iux1, iux2;   // index positions for slowness grid
			Three_Component_Seismogram::iterator s3cptr;
			int n30;  // convenient since top surface is at n3-1
			// we get a slowness vector for the stacks but we need to 
			// use the event location and compute slowness vectors for
			// each point for the incident wavefield.  
			ustack = get_stack_slowness(pwdata);
			hypo = get_event_hypocenter(pwdata);
			area = pwdata->md.get_string("area");

			// find the GCLgrid that matches the slowness vector for this ensemble
			string get_gridname(db, area, iux1, iux2, ustack);
			gridname = get_gridname(ustack);
			// note this pointer has raygrid as an alias
			grdptr = new GCLgrid3d(db,gridname.c_str());
			n30 = raygrid.n3 - 1;
			// This function returns grid of S wave travel times for this ray-based GCLgrid
			// because we do not try to interpolate this function laterally we store this
			// only as a 3D array of dimension parallel to the parent GCLgrid.  This
			// is a dangerous programming style but desirable to save space
			// THIS NEEDS A CHANGE.  THIS SHOULD BE MOVVED TO THE MAIN LOOP BELOW
			// AND ONLY RETURN A VECTOR OF TIMES.  THIS WILL ALSO ALLOW THE SAME
			// FUNCTION TO BE USED TO COMPUTE P WAVE TIMES
			ts = compute_Swave_gridtime(Vs3d,Vs1d,grdptr);
			//We need to compute this pair of ray paths to compute solid 
			// angles for the integration variable.  One is used to 
			// compute angle subtended by adding du to the x component
			// and the other is for adding du to y component.  
			// The sign of the component is used to choose if the
			// du constant is to be added or subtracted.  We subtract
			// du when u component is > 0 and add otherwise.  This
			// makes sure we don't compute a perturbation that would
			// turn inside the image volume
			Slowness_vector updux=ustack;
			Slowness_vector upduy=ustack;
			if(updux.ux<=0)
				updux.ux+=du;
			else
				updux.uy-=du;
			if(upduy.uy<=0)
				upudy.uy+=du;
			else
				upudy.uy-=du;
			rayupdux=new RayPathSphere(Vs1d,updux.mag(),
				zmax,tmax,dt,"t");
			rayupduy=new RayPathSphere(Vs1d,upduy.mag(),
				zmax,tmax,dt,"t");

			// loop through the ensemble associating each trace with a 
			// point in the raygrid and interpolating to get time to depth.
			// that is the basic loop but there are lots more details

			for(s3cptr=pwdata->tcse.begin();s3cptr!=pwdata->tcse.end();++s3cptr)
			{
				Three_Component_Seismogram& seis3c=*s3cptr;
				dmatrix work(3,seis3c.nsamp);
				Slowness_vector uincident;
				dmatrix Pinc_ray();
				dmatrix gradTs(3,raygrid.n3),gradTp();
				double zmax,tmax;
				RayPathSphere *ray0;
				i = seis3c.md.get_int("ix1");
				j = seis3c.md.get_int("ix2");
				if(i<0) || i>raygrid.n1 || j<0 || j>raygrid.n2)
					elog_die(0,"grid index of data = (%d,%d) is outside range or ray grid\nCheck migration_precompute parameters and/or input stream metadata\n",
						i,j);
				// assumes zero datum to estimate incident wave slowness vector
				// NOTE THERE SHOULD BE AN OPTION HERE.  COMPUTE FROM A MODEL
				// OR USE INPUT STORED IN METADATA
				uincident = hypo.pslow(raygrid.lat[i][j][n30],
						raygrid.lon[i][j][n30],
						0.0);
				// We need a reference ray for the incident wavefield.
				// This ray is chopped up and transformed to GCLgrid coordinates
				// in the integration loop below.  We compute it here to 
				// avoid constantly recomputing it
				zmax = raygrid.depth(i,j,raygrid.r[0]);
				zmax *= 1.2;  // Small upward adjustment to avoid rubber banding
						// problems in 3d media
				tmax=10000.0;  // arbitrary large number.  depend on zmax 
				ray0=new RayPathSphere(Vp1d,uincident.mag(),
								zmax,tmax,dt,"t");

				//  This incident ray is now reused in the loop below
				//  we work up the ray path (0 is the bottom) computing
				//  times and the weight parameters for the inverse 
				//  scattering formulation.  
				for(k=0;k<raygrid.n3;++k)
				{
					double tlag;
					double nup[3],nus[3];
					double gradP[3],gradS[3],gsum[3];
					double nrmgsum;  // norm of gradP + gradS
					Pinc_ray = compute_P_incident_ray(raygrid,i,j,k,
							uincident,ray0);
					// get unit vector directions for incident P and 
					// scattered S and store in nup and nus respectively
					// The 2 point increment is not a mistake for the 
					// P ray.  This is done to avoid a possible
					// roundoff error in the first point that is randomly
					// truncated in the preceding function
					// Note these point in the gradP and gradS directions
					// This requires opposite signs on the finite differences
					// as our ray geometry is defined internally here
					for(l=0;l<3;++l) nup[l]=Pinc_ray(l,2)-Pinc_ray(l,0);
					if(k==raygrid.n3-1)
					{
						nus[0]=raygrid.x1[i][j][k-1]-raygrid.x1[i][j][k];
						nus[1]=raygrid.x2[i][j][k-1]-raygrid.x2[i][j][k];
						nus[2]=raygrid.x3[i][j][k-1]-raygrid.x3[i][j][k];
					}
					else
					{
						nus[0]=raygrid.x1[i][j][k]-raygrid.x1[i][j][k+1];
						nus[1]=raygrid.x2[i][j][k]-raygrid.x2[i][j][k+1];
						nus[2]=raygrid.x3[i][j][k]-raygrid.x3[i][j][k+1];
					}
					// Need to call a function here to get vp0 an vs0 = velocity
					// at this point in space.  needs to be bombproof and give
					// a 1d reference model if not defined.  
					tp0 = ???
					ts0 = ???
					dr3sxv(vp0,nup,gradP);
					dr3sxv(vs0,nus,gradS);
					dr3add(gradP,gradS,gsum);
					// this routine computes the solid angle that defines
					// the variable of integration in the inverse radon
					// transform formula
					domega = compute_solid_angle(ustack,du,gradS,gradP);
					tlag = compute_lag(raygrid,i,j,k,
							Pinc_ray,uincident,Vp3d,Vs3d);
		


				}
				delete ray0;
				delete rayupdux;
				delete rayupduy;
			}

			delete grdptr;
			free_3dgrid_contiguous(ts,raygrid.n1,raygrid.n2);

		}
	}
	catch (GCLgrid_dberror gcldbe)
	{
		gcldbe.log_error();
		die(1,"GCLgrid library fatal error\n");
	}
	catch (Metadata_error mde)
	{
		mde.log_error();
		die(1,"Required processing metadata missing from input stream\n");
	}
}

into pseudocode mode:
Modified August 20, 2003 to not use precompute function

	set up grid stack volume as a GCLvectorfield3d with 3 vectors using pf or metadata object
	open database
	get 3d P velocity model GCLscalarfield from database
	get 3d S velocity model GCLscalarfield from database
	check consistency (require P and S to match exactly, implement with an == operator in gclgrid library) throw exception and exit if not (DONE:  now in libgclgrid )`
	get 1d P velocity model
	get 1d S velocity model
	get parameters for slowness grid object
	build slowness grid object
	start read thread to bring in data through pfstream

	while( data are available get next ensemble )
	{
	    foreach ensemble member
	    {
		build 3c trace object T
		get ux, uy from T
		get lat, lon of pseudostation position
		get ux0, uy0 (incident wave slowness vector)
		trace ray for ux,uy
		trace base ray for ux0,uy0
		match ux,uy to slowness grid (may use index)
		determine updux = u+dux
		determine upduy = u+duy
		trace rays for updux and upduy
		refernce all rays to S3d and P3D
		compute pathintegral for u using S3d
		for i=0 to end of S raypath
			compute P ray path from point i to surface
			compute pathintegral of truncated ray path using P3d
			compute deltax = xr - xp where xp = emergence point of p ray
			compute lag using dt[i] = TS-TP+(?)u0.deltax
			compute gamma_s and gamma_p (unit vectors) (NEW FUNCTION ABOVE)
			compute transformation matrices for each ray point
			compute tau_P and Tau_S
			compute || tau_P - tau_S ||
			compute solid angle domega[i]
			compute weight[i]
		interpolate data onto dt[i] grid (resample operation)
		for i=0 to end of S raypath
			data[i] *= weight[i]
		add metadata required 
		push to output
	    }
	}
		
// some items about object output form:
- needs to be <vector> type of higher order objectg
- object needs lat,lon, r format not GCLgrid cartesian coordinates
- object should have ncomponents 
			
		
			
		
/// older pseudocode
		crack input to build a 3C ensemble object of input data
		get slowness vector from input metadata
		match to databased
		retrieve and create GCLgrid3d object for this slowness
		map S velocities to this grid
		compute 3d travel times in this grid
		
		get slowness vector for incident wavefield
		foreach 3c ensemble member
			compute ray geometry for inident ray 
			get nearest neighbor slowness vectors
			use nearest neighbors to compute solid angle weight vector along ray
			match data to surface GCLgrid point
			compute travel time delay function 
			map data from time to ray GCLgrid for this point (need to be standard coord here)
			compute plane wave migration weights from grad Ts and grad Tp vectors
			compute transformation matrices for each depth point 
			apply depth variable transformation matrices
			apply weights
			sum to image grid
			accumulate weights for standard grid
		}
	}
	save output grid to database
}

Current gaps:

1. how complex do the depth variable transformation matrices need to be?  Are they even necessary
2. How is this best done in MPI?  The algorithm above only works serial.  To run MPI I think it
would be better to just push the data out to a stacking processor program that is not necessarily
connected.  Alternatively, the results could be written to a database table,  but that seems
awkward.   
			




