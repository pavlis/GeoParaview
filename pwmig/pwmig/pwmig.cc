#include <cstdio>
#include <string>
#include "stock.h"
#include "db.h"
#include "seispp.h"
#include "vmodel.h"

void usage()
{
        cbanner((char *)"$Revision: 1.3 $ $Date: 2003/11/28 17:42:44 $",
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
/* Returns a 3x nx  matrix of ray path unit vectors pointing in
direction of tangent to the path defined by the 3 x nx  input matrix 
passed as ray.  The algorithm is a simple forward difference scheme 
with no checking for possible roundoff problems.  Because we use forward
differences the output vector at point i is computed from the segment 
immediately forward (x[i+1]-x[i]).  Note this means the sign convention
is determined solely from the path variables.  Each vector is normalized
to unity.  The last point two points in the output are equal to keep
sizes consistent and is the only choice with a forward difference scheme.

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
	gptr = new dmatrix(3,nx);

	// debug trap
	if(sz[1]<2 || sz[0]!=3)
	{
		cerr << "ray_path_tangent:  input path size "<<sz[0]<<" by " <<sz[1]<2 <<endl;
		exit(1);
	}

	for(i=0;i<nx-1;++i)
	{
		dx[0]=ray(0,i+1)-ray(0,i);
		dx[1]=ray(1,i+1)-ray(1,i);
		dx[2]=ray(2,i+1)-ray(2,i);
		nrmdx=dnrm2(3,dx,1);
		// debug trap.  Won't happen if used correctly. Could to be removed eventually
		if(nrmdx<=0.0)
		{
			cerr<<"ray_path_tangent: coding error two points on path are equal."<<endl;
			exit(1);
		} 
		dscal(3,1.0/nrmdx,dx,1);
		for(j=0;j<3;++j) gamma(j,i)=dx[j];
	}
	// copy the last point
	for(j=0;j<3;++j) gamma(j,nx-1)=gamma(j,nx-2);
	return(gptr);
}
/* similar function to above in concept, but this one computes the set of local
verticals along the ray path in the GCLgrid cartesian system
*/
dmatrix *compute_local_verticals(dmatrix& ray)
{
	dmatrix *vptr;
	dmatrix& verticals=*vptr;
	int *sz,nx;
	Geographic_point x0_geo;
	Cartesian_point x0_c;
	const double DR=100.0;

	double dx[3];
	double nrmdx;
	int i,j;

	sz = ray.size();
	// assume number rows = 3.   could test and throw exception, but this
	// function is not viewed as a library function
	nx=sz[1];
	delete [] sz;
	vptr = new dmatrix(3,nx);
	for(i=0;i<nx;++i)
	{
		x0_geo = g.ctog(path(0,0),path(1,0),path(2,0));
		x0_geo.r -= DR;
		x0_c = g.gtoc(x0_geo);
		verticals(0,i) = (path(0,i) - x0_c.x1)/DR;
		verticals(1,i) = (path(1,i) - x0_c.x2)/DR;
		verticals(2,i) = (path(2,i) - x0_c.x3)/DR;
	}
	return(vptr);
}
/* computes depth dependent transformation matrix assuming a scattering
is specular.  That is, transformation matrices will yield a form of
L,R,T 
*/

class Ray_Transformation_Operator
{
public:
	int npoints;
	Ray_Transformation_Operator(int np)
	{
		for(int i=0;i<np;++i)
		{
			dmatrix *dmtp = new dmatrix(3,3);
			U.push_back(*dmtp);
		}
	}
	// constructor for constant azimuth 
	Ray_Transformation_Operator(GCLgrid&g, dmatrix& path, double azimuth);
	// constructor for depth dependent operator
	Ray_Transformation_Operator(GCLgrid&g, dmatrix& path);
	~Ray_Transformation_Operator(){delete [] U;);
	dmatrix *apply(dmatrix& in);
private:
	vector<dmatrix> U;
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
/* Constructor for simple case with all matrices going to surface R,T,L coordinates.
The primary contents is a vector of dmatrix object holding (in this case) the
same transformation matrix at every point.  

Arguments:
	g - reference GCLgrid object.  Mainly used here as the coordinate reference
		frame for path
	path - 3xNp matrix of points that defines the ray path
	azimuth - angle (radians) of the ray propagation direction 
		at the surface.  Note this is geographical azimuth NOT the
		phi angle in spherical coordinates.
*/
Ray_Transformation_Operator::Ray_Transformation_Operator(GCLgrid& g, 
	dmatrix& path,double azimuth)
{
	int *insize=path.size();
	int np=insize[1];
	delete [] insize;
	dmatrix U0(3,3);  // final matrix copied to all elements here
	double x[3];  //work vector to compute unit vectors loaded to U0
	Geographic_point x0_geo;
	Cartesian_point x0_c;
	double x_vertical[3];  // point to local vertical direction at x_geo
	const double DR=100.0;
	double xdotxv,theta,phi;
	int i;

	//first we need the geographic coordinates of the ray emergence point
	x0_geo = g.ctog(path(0,0),path(1,0),path(2,0));
	x0_geo.r -= DR;
	x0_c = g.gtoc(x0_geo);
	x_vertical[0] = (path(0,0) - x0_c.x1)/DR;
	x_vertical[1] = (path(1,0) - x0_c.x2)/DR;
	x_vertical[2] = (path(2,0) - x0_c.x3)/DR;
	
	// Get the L direction from the first pair of points in path
	x[0]= path(0,0) - path(0,1);
	x[1]= path(1,0) - path(1,1);
	x[2]= path(2,0) - path(2,1);
	nrmx = dnrm2(3,x,1);
	dscal(3,1.0/nrmx,x,1);  // normalize to unit vector
	// need to check for a vertical vector and revert to identity 
	// to avoid roundoff problems
	xdotxv = ddot(3,x,1,x_vertical,1);
	// azimuth is from North while theta in spherical coordinates is measured from x1
	phi = M_PI_2 - azimuth;
	if(fabs(xdotxv)<=DBL_EPSILON)
	{
		// L vertical means theta=0.0
		a=cos(phi);
		b=sin(phi);
		c=1.0;
		d=0.0;
	}
	else
	{
		// dot product is preserved so we can use this here
		theta = acos(xdotxv);
		/* The following is copied from the rotate function in rotation.C of libseispp.
		I did this to allow parallel coding for this and the more general case below.
		It causes a maintenance problem as if there is an error there it will have
		been propagated here. */
	        a = cos(phi);
	        b = sin(phi);
	        c = cos(theta);
	        d = sin(theta);
	}
	U0(0,0) = a*c;
        U0(1,0) = b*c;
        U0(2,0) = d;
        U0(0,1) = -b;
        U0(1,1) = a;
        U0(2,1) = 0.0;
        U0(0,2) = -a*d;
        U0(1,2) = -b*d;
        U0(2,2) = c;
	for(int i=0;i<np;++i)
	{
		// Confusing C++ memory use here.  HOpe this is right
		dmatrix dmtp(3,3);
		dmatrix *dmtp = new dmatrix(3,3);
		*dtmp = U0;
		U.push_back(*dmtp);  
	}
}
/* Now the much more complex algorithm for the case where we assume the ray path
is a specular reflection.  Then we have a different ray transformation matrix 
at each point.  We compute this using this algorithm:
  The L coordinate turns upward.  At the scatter point we get it from the scattered ray path
  direction and assume that is one thing we can anchor on.

  The scattered S has a polarization in the plane of scattering formed by the incident P and
  scattered S.  The complementary component is then readily found by the cross product of the
  P and S paths.  

  We get the scattered S direction (at the scattering point) from a cross product between 
  local L and the local "SH-like" direction derived in the previous step.

  We propagate the components to the surface with a simple directional change.  L remains
  L  and the S components obey Snell's law rules.  

Arguments:
	g-parent GCLgrid that defines coordinate system (only used for coordinate
		transformations here)
	path - 3xnp matrix of points defining the path.  It is assumed path starts at
		the surface and is oriented downward.
	azimuth - angle (radians) of the ray propagation direction 
		at the surface.  Note this is geographical azimuth NOT the
		phi angle in spherical coordinates.
	gamma_P - 3xnp matrix of tangent (unit) vectors for incident wave ray path
		at each point in path.  It is assumed gamma_P vectors point along
		the direction of propagation of the incident wavefield (nominally upward)

The output array of matrices when applied to data will yield output with x1 = generalized
R, x2= generalized T, and x3=L.  

*/
Ray_Transformation_Operator::Ray_Transformation_Operator(GCLgrid& g, 
	dmatrix& path,double azimuth, dmatrix& gamma_P)
{
	int *insize=path.size();
	int np=insize[1];
	delete [] insize;
	Geographic_point x0_geo;
	Cartesian_point x0_c;
	double x_vertical[3];  // point to local vertical direction at x_geo
	const double DR=100.0;
	const double PARALLEL_TEST=FLT_EPSILON;  // a conservative test for zero length double vector
	int i,j;

	// First we need to compute a set of ray path tangents
	dmatrix *tanptr;
	dmatrix &tangents = tanptr;
	tanptr = ray_path_tangent(path);

	// We next need a set of local vertical vectors. 
	dmatrix *vptr;
	dmatrix &local_verticals = *vptr;
	vptr = compute_local_verticals(path);

	// We call the simpler constructor first above and use it to provide
	// transformation matrices to ray coordinates -- the starting point here
	// For old C programmers like me:  This will be deleted when it goes out
	// of scope when the function exits so we don't need a delete below.
	Ray_Transformation_Operator Raytrans0(g,path,azimuth);

	// work down the ray path building up the transformation matrix at each point
	U=new dmatrix(3,3)[np];	
	for(i=0;i<np;++i)
	{
		double Lscatter[3],Tscatter[3],Rscatter[3];
		double nu0;  // unit vector in direction gamma_P
		// Tp, Rp, and Zp for an orthogonal basis for earth coordinates
		// at the scattering point that are standard 1D propagator coordinates
		// That is Zp is local vertical, Rp is Sv director for S ray path,
		// and Tp is Sh.  
		double Zp[3],Rp[3],Tp[3];  // radial and tangential for S ray path 
		dmatrix work(3,3);
		// copy the tangent vector Lscatter from the tangent vectors 
		// matrix for convenience and clarity.  Assume Lscatter is unit vector
		dcopy(3,tangents.get_address(0,i),1,Lscatter,1);
		// do the same for the incident ray direction but here we want to
		// normalize to a unit vector
		dcopy(3,gamma_P.get_address(0,i),1,nu0,1);
		d3norm(nu0);
		// First derive the S ray path propagation coordinate system
		// first copy Zp
		dcopy(3,local_verticals.get_address(0,i),1,Zp,1);
		// Now get radial from Lscatter = S ray path tangent.  Have to 
		// handle case with Lscatter vertical carefully.  In that case, 
		// derive radial from nu0.
		d3cros(Zp,Lscatter,Tp);
		if(dnrm2(3,Tp,1)<PARALLEL_TEST)
                {
			d3cros(Zp,nu0,Tp);
			// excessively paranoid, but could happen
			if(dnrm2(3,Tp,1)<PARALLEL_TEST)
                        {
                                cerr << "Warning:  cannot handle singular geometry"
					<< endl;
				<< "Antipodal event scattered to vertical incidence"
					<< endl
				<< "Using T=x1 and R=x2 of GCL coordinate system"
					<< endl
				<< "Probable discontinuity at normal incidence"<<endl;
				Tp[0]=1.0;  Tp[1]=0.0;  Tp[2]=0.0;
			}
		}
		d3norm(Tp);
		// Radial is now derived as Lscatter X Tp to make right handed cartesian
		d3cros(Lscatter,Tp,1);
		// Normalization d3norm(Rp) not necessary here because Tp and Lscatter
		// are orthogonal unit vectors
		//
		// Now we need to derive the specular reflection SV vector (Rscatter)
		// and SH vector (Tscatter).  A confusion is this is much like the
		// way we computed Rp and Tp above but using different vectors as
		// building blocks.  First, Tscatter is derived from cross product of
		// nu0 and Lscatter.  Sign is important.  We aim for a coordinate system
		// where T,R,L = x1,x2,x3 form a right handed coordinate system
		// This happens to be a nonstandard convention (textbooks would
		// use L,R,T = x,y,z) but I'm deriving too much code from 
		// multiwavelet code that used this convention.  
		//
		d3cros(Lscatter,nu0,Tscatter);  // Tscatter is x1
		// as above, need to handle case when L and T are parallel
		if(dnrm2(3,Tscatter,1)<PARALLEL_TEST)
		{
			// In this case, Tscatter will be the same as Tp
			dcopy(3,Tp,1,Tscatter,1);
		}
		// Do need to normalize
		d3norm(Tscatter);  
		// In a right handed coordinate system x2 = x3 X x1 = L X T
		d3cros(Lscatter,Tscatter,Rscatter);

		//
		/* Now we derive the transformation matrix.  The matrixs computed
		computed here will transform a vector in ray coordinates to a vector
		with in the specular reflection coordinate system used here.
		Note that for the present case of a 1d medium this is a simple
		rotation around x3, but we derive it from dot products.  
		In the symbols of this algorithm it transforms from the 
		Tp,Rp, Lscatter basis to that of Tscatter, Rscatter, Lscatter.
		This works because we assume a simple ray theory propagator to
		project the wavefield to the depth of the scattering point.
		(Note: my working theory notes derived the inverse of this
		transform = transpose of this one.
		*/

		work(0,0) = ddot(3,Tscatter,1,Tp,1);
		work(0,1) = ddot(3,Tscatter,1,Rp,1);
		work(0,2) = 0.0;
		work(1,0) = ddot(3,Rscatter,1,Tp,1);
		work(1,1) = ddot(3,Rscatter,1,Rp,1);
		work(1,2) = 0.0;
		work(2,0) = 0.0;
		work(2,1) = 0.0;
		work(2,2) = 1.0;

		U[i] = work*Raytrans0.U[i];

	}

	delete tangents;
	delete local_verticals;
}
/* This function computes the domega term for the inverse Radon transform inversion
formula.  It works according to this basic algorithm.
1.  If u0 is the slowness vector for the current ray path, compute four neighboring
rays at 1/2 grid spacings:  u0+dux1/2, u0-dux1/2, u0+duy1/2, and u0-duy1/2 where the
1/2 means the 1/2 a grid spacing.  In slowness space this forms a rectangle surrounding
the u0 point.
2.  Compute a ray for each of these points using constructors for a RayPathSphere 
object and the routine GCLgrid_Ray_Project.  Note these are computed with the equal
depth option, not equal time.
3.  Compute gradS for each point on each ray.  
4.  Interpolate the input gradP vectors to depths (the input is on an irregular
grid).  
5.  Compute the gradS+gradP sums to build an array of unit vectors for each ray.
6.  compute solid angles between + and - pairs of finite difference rays using
a vector cross product and a small angle approximation (if a and b are unit vectors
axb = |a||b|sin theta = sin theta ~ theta for small angles )  Solid angle is product
of the two angles.

Returns a vector of domega terms on the grid defined by the input ray path.

Arguments:

	u0 - input slowness vector (note passed as Slowness_vector object)
	dux, duy - base slowness grid increments.
	vmod - 1d velocity model used for ray tracing 
	zmax, dz - ray trace parameters.  trace to depth zmax with depth increment dz
	g - gclgrid for geometry
	ix1, ix2  - index position of receiver point in gclgrid g.
	path - ray path of scattered S in GCLgrid cartesian coordinate system.
	gradP - parallel matrix to path of P time gradient vectors.

Returns an np by 1 dmatrix object containing the domega term.  I chose to use this
instead of the STL vector to reduce baggage.  All I need is the numbers in the vector and the
convenience of the size elements of the dmatrix np (length of output) being passed iwth
the object.  a dmatrix is fine for that and the contents can be manipulated with the BLAS.

Author:  Gary Pavlis
Written:  NOvember 2003
*/
	

double *compute_domega_for_path(Slowness_vector& u0,double dux, double duy,
	Velocity_Model_1d& vmod, GCLgrid3d& g,int ix1, int ix2,
	dmatrix& path, dmatrix& gradP)
{
	Slowness_vector udx1, udx2, udy1, udy2;
	dmatrix *pathdx1, *pathdx2, *pathdy1, *pathdy2;
	dmatrix *gradSdx1, *gradSdx2, *gradSdy1, *gradSdy2;
	int i;
	int npath; // used to define length of valid output 

	udx1 = u0;
	udx1.ux -= dux/2.0;
	udx2=u0;
	dux2.ux += dux/2.0;
	udy1 = u0;
	udy1.uy -= duy/2.0;
	udy2=u0;
	udy2.uy += duy/2.0;
	//
	// Now call the constructors for the basic ray paths for these four rays.  
	// We set the constructor in depth model and ASSUME zmax is small enough that
	// the ray for each u will not turn within vmod.  
	//
	RayPathSphere raydx1(vmod,udx1,zmax,1.0e99,del,"z");
	//
	// Now call the constructors for the basic ray paths for these four rays.  
	// We set the constructor in depth model and ASSUME zmax is small enough that
	// the ray for each u will not turn within vmod.  
	//
	RayPathSphere raydx1(vmod,udx1.mag(),zmax,1.0e99,dz,"z");
	RayPathSphere raydx2(vmod,udx2.mag(),zmax,1.0e99,dz,"z");
	RayPathSphere raydy1(vmod,udy1.mag(),zmax,1.0e99,dz,"z");
	RayPathSphere raydy1(vmod,udy1.mag(),zmax,1.0e99,dz,"z");
	//
	// project these into the GCLgrid coordinate system
	//
	try {
		pathdx1 = GCLgrid_Ray_Project(g,raydx1&, udx1.azimuth(),ix1,ix2);
		pathdx2 = GCLgrid_Ray_Project(g,raydx2&, udx2.azimuth(),ix1,ix2);
		pathdy1 = GCLgrid_Ray_Project(g,raydy1&, udy1.azimuth(),ix1,ix2);
		pathdy2 = GCLgrid_Ray_Project(g,raydy2&, udy2.azimuth(),ix1,ix2);
	}  catch (GCLgrid_error)
	{
		err.log_error();
		die(0,"Unrecoverable error:  requires a bug fix\n");
	}
	//
	// get gradS values by first computing tangents as unit vector and then
	// scaling them by slowness at that depth. Note 1d slowness is used for
	// this calculation, not a 3d value.  This is an approximation that avoid
	// tremendous headaches that could occur otherwise.  It should not be a
	// significant error.
	//
	gradSdx1=ray_path_tangent(pathdx1);
	gradSdx2=ray_path_tangent(pathdx2);
	gradSdy1=ray_path_tangent(pathdy1);
	gradSdy2=ray_path_tangent(pathdy2);
	//
	// Check all these matrices have the same size and truncate the 
	// output if necessary writing a diagnostic in that situation.
	//
	npath = gradSdx1.nc;
	if(npath>gradSdx2.nc) npath=gradSdx2.nc;
	if(npath>gradSdy1.nc) npath=gradSdy1.nc;
	if(npath>gradSdy2.nc) npath=gradSdy2.nc;
	if(npath!= gradSdx1.nc)
	{
		cerr << "compute_domega_for_path: irregular path length for bundle of 4 rays used to compute domega" << endl

			<< "Slowness grid probably defines turning rays in thie model" << endl;
	}
	for(i=0;i<npath;++i)
	{
		double slow;
		double depth;
		// because we used equal depth grid, we can just choose one of these to get depth
		depth = R_EQUATOR - raydx1.r[i];  
		slow = 1.0/vmod.getv(depth);
		dscal(3,slow,gradSdx1.get_address(0,i),1);
		dscal(3,slow,gradSdx2.get_address(0,i),1);
		dscal(3,slow,gradSdy1.get_address(0,i),1);
		dscal(3,slow,gradSdy2.get_address(0,i),1);
	}

	//
	// Now we have interpolate the gradP vectors onto this uniform depth grid.  
	// to do that we use the GCLgrid's geometry function to build a vector of 
	// depths
	//
	int npath_P;
	npath_P = gradP.nc;
	double *z = new double[npath_P];
	dmatrix gradP_regular(3,npath);
	Geographic_point gclp;
	for(i=0;i<npath_P;++i)
	{
		gclp = g.ctog(gradP(0,i), gradP(1,i),gradP(2,i));
		z[i] = r0_ellipse(gclp.lat) - gclp.r;
	}
	INTERPOLATOR1d::linear_vector_irregular_to_regular(z,gradP,
		0.0,dz,gradP_regular);
	delete [] z;   // finished with this so we better get rid of it
	//
	// Now we can compute sum of gradP and gradS to form unit vectors
	// the nu matrices are created with a 1,1 dimension initially to 
	// define the object, but are assumed rebuilt by the assignment operator
	// below.
	//
	dmatrix nudx1(1,1);
	dmatrix nudx2(1,1);
	dmatrix nudy1(1,1);
	dmatrix nudy2(1,1);
	nudx1 = gradP_regular + gradSdx1;
	nudx2 = gradP_regular + gradSdx2;
	nudy1 = gradP_regular + gradSdy1;
	nudy2 = gradP_regular + gradSdy2;
	// normalize using antelope function d3norm
	for(i=0;i<npath;++i)
	{
		d3norm(nudx1.get_address(0,i));
		d3norm(nudx2.get_address(0,i));
		d3norm(nudy1.get_address(0,i));
		d3norm(nudy2.get_address(0,i));
	}
	//
	// Now get domega using cross products between pairs of unit vectors and 
	// small angle approximation (sin theta approx theta when theta is small)
	//
	double *doptr = new dmatrix(npath,1);
	dmatrix domega& = *doptr;
	for(i=0;i<npath;++i)
	{
		double cross[3];
		double dtheta_x, dtheta_y;
		d3cros(nudx1,nudx2,cross);
		dtheta_x = r3mag(cross);
		d3cros(nudy1,nudy2,cross);
		dtheta_y = r3mag(cross);
		domega(i,0) = dtheta_x*dtheta_y;
	}

	// because these are defined as pointers we have to explicitly delete these
	// objects
	delete pathdx1;  delete pathdx2;  delete pathdy1;  delete pathdy2;
	delete gradSdx1;  delete gradSdx2;  delete gradSdy1;  delete gradSdy2;

	return(doptr);  // pointer to domega defined above 
}
	

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
/* This was a design change.  Originally I planned a precompute function.
Decided the work involved was not that severe and we could just recompute
this grid every time.  Now we have to construct this GCLgrid incrementally
That is, the final result is a GCLvectorfield3d object with the projected
data into that field.  
			string get_gridname(db, area, iux1, iux2, ustack);
			gridname = get_gridname(ustack);
			// note this pointer has raygrid as an alias
			grdptr = new GCLgrid3d(db,gridname.c_str());
			n30 = raygrid.n3 - 1;
end area needing replacement*/
			grdptr =  Build_GCLraygrid(parent,ustack,vs1d,zmax,tmax,dt);
			// here we need to call a constructor that builds a 
			// new 2d grid object with seismograms associatged with points



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
			compute gamma_s and gamma_p (unit vectors) and store in dmatrix
			compute transformation matrices for each ray point
			compute tau_P and Tau_S
			compute || tau_P - tau_S ||
			compute solid angle domega[i]
			compute weight[i]
		interpolate data onto rathpath times (resample operation)
		compute vector of transformation matrix for each point on path
		apply transformation matrices to data
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
			




