#include <cstdio>
#include <string>
#include "stock.h"
#include "db.h"
#include "seispp.h"
#include "vmodel.h"

void usage()
{
        cbanner((char *)"$Revision: 1.5 $ $Date: 2004/10/30 14:09:31 $",
                (char *)"db  [-V -pf pfname]",
                (char *)"Gary Pavlis",
                (char *)"Indiana University",
                (char *)"pavlis@indiana.edu") ;
        exit(-1);
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
double ***compute_Swave_gridtime(GCLscalarfield3d& Vs, Velocity_model& Vs1d, GCLgrid3d& raygrid)
{
	double ***Stime;
	int i,j,k;
	int nfailures=0;
	int n30;
	Geographic_point gp;
	Cartesian_point cx;
	
	// don't trap memory errors here as the function dies if allocs fail
	Stime = create_3dgrid_contiguous(raygrid->n1, raygrid->n2,raygrid->n3);
	Vs.reset_index();	
	n30=raygrid.n3 - 1;  // top of the grid is the last entry

	for(i=0;i<raygrid.n1;++i)
		for(j=0;j<raygrid.n2;++j)
		{
			double t,dt;
			// First fill the top surface values with zeros as the time reference datum.
			Stime[i][j][n30]=0.0;

			for(k=n30-1,t=0.0;k>=0;--k)
			{
				int lret;
				int index[3];
				// we have to convert to geo coordinates and then 
				// into reference from for S velocity model 
				gp = raygrid.geo_coordinates(i,j,k);
				cx = Vs.gtoc(gp.lat,gp.lon,gp.r);

				lret = Vs.lookup(cx.x1,cx.x2,cx.x3);
				if(lret==0)
				{
					Vs.get_index(index);
					vel = Vs.interpolate(cx.x1,cx.x2,cx.x3);
				}
				else
				{
				// silently revert to the 1D reference model on any failure.
				// We do count how many times we fail, however, and issue a
				// diagnostic if every call failed.
					++nfailures;
					vel = Vs1d.getv(raygrid.depth(i,j,k));
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

int main(int argc, char **argv)
{
        Dbptr db,dbv;
	string Pmodel1d_name;
	const string pvfnm("Pvelocity"), svfnm("Svelocity");
	Three_Component_Ensemble *pwdata;
	int i,j,k;
	string pfin("pwmig");


        ios::sync_to_stdio();
        elog_init(argc,argv);

        if(argc < 2) usage();
	string dbname(argv[1]);

        for(i=2;i<argc;++i)
        {
                if(!strcmp(argv[i],"-V"))
                        usage();
                else if(!strcmp(argv[i],"-pf"))
                {
                        ++i;
                        if(i>=argc) usage();
                        pfin = string(argv[i]);
                }
		else
			usage();
	}
        if(pfread(pfin.c_str(),&pf)) die(1,(char *)"pfread error\n");
	try {
		Metadata control(pfin);
		string Pmodel3d_name=control.get_string("P_velocity_model3d_name");
		string Smodel3d_name=control.get_string("S_velocity_model3d_name");
		string Pmodel1d_name=control.get_string("P_velocity_model1d_name");
		string Smodel1d_name=control.get_string("S_velocity_model1d_name");
		string parent_grid_name=control.get_string("Parent_GCLgrid_Name");
		// Create a database handle for reading data objects
		dbopen(dbname,'r+',&db);
		if(db.record == dbINVALID) die(1,"Cannot open database %s\n",dbname);
		
		GCLscalarfield3d Vp3d(db,Pmodel3d_name.c_str(),pvfnm.c_str());
		GCLscalarfield3d Vs3d(db,Smodel3d_name.c_str(),svfnm.c_str());
		double ***Stimes;  // holds parallel array of S travel times for raygrid
		Velocity_Model Vp1d(db,Pmodel_name,pvfnm);
		Velocity_Model Vs1d(db,Smodel_name,svfnm);
		GCLgrid parent(db,parent_grid_name.c_str());


//YET ANOTHER NOTE.  HAVE TO DECIDE HOW TO HANDLE ZERO DATUM.  DO WE DO A GEOMETRIC STATIC
// AS IN REFLECTION PROCESSING?  PROBABLY SHOULD PUT THIS IN PWSTACK.  THAT IS STACK SHOULD
// SHOULD ALWAYS BE SHIFTED TO DATUM.
// Solution for now.  Modified pwstack to put an elev estimate for the pseudostation
// grid point computed as a weighted sum of elevations of contributing stations.  Algorithm
// below should use a vertical ray approximation to this static to datum or it will get
// unnecessarily complex.  A p=0 approximation for an elevation correction is pretty
// good for teleseismic data and is unlikely to be a big problem.
		// Load the default attribute map
		Attribute_Map am();
		// Need this set of attributes extraction lists to 
		// load correct metadata into the ensmble and each member
		Metadata_list mdlin=pfget_mdlist(pf,"Station_mdlist");
		Metadata_list mdens=pfget_mdlist(pf,"Ensemble_mdlist");
		Datascope_Handle dbh(db,pf,"dbprocess_commands");
		dbh.rewind();
		for(int record=0;record<dbh.number_tuples();++dbh,++record)
		{
			// first read in the next ensemble
			pwdata = new Three_Component_Ensemble(dynamic_cast<Database_Handle&>
					(dbh),mdlin,mdens,am);
			Hypocenter hypo;
			hypo.lat=pwdata->get_double("origin.lat");
			hypo.lon=pwdata->get_double("origin.lon");
			hypo.z=pwdata->get_double("origin.depth");
			hypo.time=pwdata->get_double("origin.time");
			// The data are assumed grouped in the ensemble by a common
			// slowness vector for the plane wave stack.  We extract this
			// from the ensemble metadata.  
			Slowness_vector ustack;
			// Kind of evil OOP style placing data members inline line is instead
			// of in a constructor.  Prone to error if class structure changes
			ustack.ux=pwdata->get_double("ux");
			ustack.uy=pwdata->get_double("uy");
			// This constucts the slowness vector for the incident ray
			// Uses a model.  Eventually may want an option to use something
			// like the output of mwap.  Note also this uses one slowness vector
			// for entire grid.  This may be a bad approximation for large 
			// aperture arrays like the USArray. 
			Slowness_vector uincident=hypo.pslow(parent.lat0,parent.lon0,0.0);
			string gridname;
			GCLgrid3d *grdptr;
			GCLgrid3d& raygrid=*grdptr;  // convenient shorthand because we keep changing this
			string area;  // key for database when multiple areas are stored in database
			RayPathSphere *rayupdux,*rayupduy;
			int iux1, iux2;   // index positions for slowness grid

			// This builds a 3d grid using ustack slowness rays with
			// 2d grid parent as the surface grid array

			grdptr =  Build_GCLraygrid(parent,ustack,vs1d,zmax,tmax,dt);
			int n30;  // convenient since top surface is at n3-1
			n30 = raygrid.n3 - 1;
			// Form 3 dimensional vector of S wave travel times from the
			// surface to points in raygrid.  This uses a 3d model
			// where define reverting to a 1d model if a ray passes outside
			// the input model
			Stimes=compute_Swave_gridtime(Vs3d, Vs1d,raygrid);

			// loop through the ensemble associating each trace with a 
			// point in the raygrid and interpolating to get time to depth.
			// that is the basic loop but there are lots more details

			for(s3cptr=pwdata->tcse.begin();s3cptr!=pwdata->tcse.end();++s3cptr)
			{
				Three_Component_Seismogram& seis3c=*s3cptr;
				dmatrix work(3,seis3c.nsamp);
				dmatrix Pinc_ray();
				dmatrix gradTs(3,raygrid.n3),gradTp();
				double zmax,tmax;
				RayPathSphere *ray0;
				i = seis3c.get_int("ix1");
				j = seis3c.get_int("ix2");
				if(i<0) || i>raygrid.n1 || j<0 || j>raygrid.n2)
				{
					cerr << "grid index of data = ("
						<< i << ","<< j 
						<< "is outside range of grid. <<endl
						<< "Check grid definitions in database and pf"<<endl;
					exit(-1);
				}
// REWRITE: OCT 24,2004 stopped here
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
			free_3dgrid_contiguous(Stimes,raygrid.n1,raygrid.n2);

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
	catch (seispp_error seer)
	{
		seer.log_error();
	}
}


