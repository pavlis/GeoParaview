#include <cstdio>
#include <string>
#include <sstream>
#include <float.h>
#include <memory>
#include "stock.h"
#include "coords.h"
#include "pf.h"
#include "db.h"
#include "dmatrix.h"
#include "gclgrid.h"
#include "ray1d.h"
#include "ensemble.h"
#include "seispp.h"
#include "interpolator1d.h"
//DEBUG
//#include "vtk_output.h"
#ifdef MATLABDEBUG
#include "MatlabProcessor.h"
#endif
using namespace std;
using namespace SEISPP;
using namespace INTERPOLATOR1D;
#include "pwmig.h"

//DEBUG  Uncomment if not using matlab for debu
#ifdef MATLABDEBUG
MatlabProcessor mp(stdout);
#endif
void usage()
{
        cbanner((char *)"$Revision: 1.17 $ $Date: 2007/04/11 19:36:28 $",
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
	dx1 = ray->x1[i0][j0][k0] - ray->x1[i0][j0][k0-1];
	dx2 = ray->x2[i0][j0][k0] - ray->x2[i0][j0][k0-1];
	dx3 = ray->x3[i0][j0][k0] - ray->x3[i0][j0][k0-1];
	return( sqrt( dx1*dx1 + dx2*dx2 + dx3*dx3));
}
/* computes travel times along a ray path defined by the input dmatrix of Cartesian
coordinate (path) using a 3D model defined by the GCLsclarfield3d variable U3d.  
The path is assumed to be congruent with the ray geometry GCLgrid (raygrid)
but not necessarily congruent with the velocity model field.  Thus we use the
not equal operator to test for conguence between raygride and U3d.  
It is assumed the U3d holds SLOWNESS values in units of s/km.

Because most 3D Earth models have finite extents it is highly likely a ray can 
traverse outside the reference model.  We allow for this here and use the 1d model
(V1d) as a background reference.  That is, when the ray moves outside the 3D refernce
model the 1d model is used.  Note this ASSUMES that the input path starts inside the
3d model box and passes outside somewhere along path.  If the path starts outside the
V3d region the 1d model will be silently used.  i.e.  the 3d model will not be used
if the path runs from outside into the box.  This is because the algorithm used is
a simple one which simply checks the length of the vector returned by the pathintegral
library function against the size of path.  This should be generalized if this were
ever made into a library function.

Note path is intentionally passed by value not by reference as remap 
path can redefine it.
*/

vector<double> compute_gridttime(GCLscalarfield3d& U3d, 
	int ix1, int ix2,
		VelocityModel_1d& V1d, 
			GCLgrid3d& raygrid,
				dmatrix& path)
{

	dmatrix thispath(path);
	// test to be sure the grids are congruent or paths will not mesh
	if(dynamic_cast<BasicGCLgrid&>(U3d)
			!=dynamic_cast<BasicGCLgrid&>(raygrid))
	{
		thispath=remap_path(raygrid,path,U3d);
	}
	vector<double> times=pathintegral(U3d,thispath);
	// This fills the output vector with 1d values if the path wanders outside
	// the 3d grid bounds
	if(times.size()!= thispath.columns())
	{
		double vel,dt;
		int k,kstart;
		// Watch for special case when stations is outside
		// the 3d model completely
		if(times.size()<=0)
			kstart=1;
		else
			kstart=times.size();
		for(k=kstart;k<thispath.columns();++k)
		{
			vel=V1d.getv(raygrid.depth(ix1,ix2,k));
			dt = raydist(&raygrid,ix1,ix2,k)/vel;
			times.push_back(times[k-1]+dt);
		}
	}
	return(times);
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
	int nx;
	double dx[3];
	double nrmdx;
	int i,j;

	nx=ray.columns();
	dmatrix *gptr;
	gptr = new dmatrix(3,nx);
	dmatrix& gamma=*gptr;

	// debug trap
	if(nx<2 || ray.rows()!=3)
	{
		cerr << "ray_path_tangent:  input path size "
			<< ray.rows() 
			<<" by " 
			<< ray.columns() <<endl;
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

/* Computes a dmatrix of gradient S vectors for a ray GCLgrid defining the x3
grid line from x1,x2 grid position ix1,ix2 (convention in this program).
Velocities used to scale the vectors are derived from the 1D velocity model
passed as Vs1d.  This is an approximation but one that should make vitually no
difference in the results for the simplicity and speed it (likely) buys.
*/

dmatrix compute_gradS(GCLgrid3d& raygrid,int ix1, int ix2, VelocityModel_1d& Vs1d)
{
	int k,kk,i;
	dmatrix gradS(3,raygrid.n3);
	double vs;
	double nrm;
	// We start at the surface and work down 
	// gradS vectors should point upward
	// The raygrid lines are oriented from the bottom up so sign is 
	// as below
	for(k=0,kk=raygrid.n3-1;k<raygrid.n3-1;++k,--kk)
	{
		gradS(0,k)=raygrid.x1[ix1][ix2][kk]
				- raygrid.x1[ix1][ix2][kk-1];
		gradS(1,k)=raygrid.x2[ix1][ix2][kk]
				- raygrid.x2[ix1][ix2][kk-1];
		gradS(2,k)=raygrid.x3[ix1][ix2][kk]
				- raygrid.x3[ix1][ix2][kk-1];
		// use the velocity from the upper point
		vs=Vs1d.getv(raygrid.depth(ix1,ix2,kk));
		// gradient is unit vector multiplied by slowness
		nrm=dnrm2(3,gradS.get_address(0,k),1);
		for(i=0;i<3;++i) gradS(i,k)/=(nrm*vs);
	}
	// Final point is just a copy of the second to last as we run out 
	// of points in forward difference
	for(i=0;i<3;++i)gradS(i,raygrid.n3-1)=gradS(i,raygrid.n3-2);
	return(gradS);
}
vector<double> compute_unit_normal_P(GCLgrid3d& raygrid,double x1, double x2, double x3)
{
	// Ray path tangents do not vary rapidly with position.  We just get the nearest
	// neighbor and compute by finite differences
	vector<double>nu; // unit vector here normally pointing generally up and in some direction
	int index[3];
	int err;
	double dx;
	double nrmx;
	double rtest;
	int k;
	Geographic_point geo;

	err=raygrid.lookup(x1,x2,x3);
	switch(err)
	{
		case 1:
		case -1:
			// If the lookup failed we start here.  Since we are using a grid
			// constructed of 1d rays we can cheat and jump to the center of the 
			// grid and try to use something reasonably approximates the corect
			// ray direction
			raygrid.reset_index();
			raygrid.get_index(index);
			geo=raygrid.ctog(x1,x2,x3);
			rtest=geo.r;
			// hunting upward for the first point on the ray through the 
			// origin that is just above the radius defined in geo
			// we set the index position there for the calculation below
			for(k=0;k<raygrid.n3;++k)
			{
				geo=raygrid.geo_coordinates(index[0],index[1],k);
				if(geo.r>rtest) 
				{
					index[2]=k;
					break;
				}
			}
			if(k>=raygrid.n3) index[2]= raygrid.n3 - 2;  // set to one below surface
		case 0:
		case 2:
			raygrid.get_index(index);
			// Standard fix for forward diff.  Back off one if at the top edge
			if(index[2]>=raygrid.n3-1)--index[2];  
			dx=raygrid.x1[index[0]][index[1]][index[2]+1];
			dx-=raygrid.x1[index[0]][index[1]][index[2]];
			nu.push_back(dx);
			dx=raygrid.x2[index[0]][index[1]][index[2]+1];
			dx-=raygrid.x2[index[0]][index[1]][index[2]];
			nu.push_back(dx);
			dx=raygrid.x3[index[0]][index[1]][index[2]+1];
			dx-=raygrid.x3[index[0]][index[1]][index[2]];
			nu.push_back(dx);
			nrmx=dnrm2(3,&(nu[0]),1);
			dscal(3,1.0/nrmx,&(nu[0]),1);
	}
	return(nu);
}

//@{
// Computes a GCLgrid defined field (raygrid) for an incident P wavefield.  
// The result is travel times from a given hypocenter to points inside a GCLgrid3d
// object.  The times are computed from a 1d refernce model to the base of the model
// and then using approximate ray tracing from the base up.  
// Actual algorithm computes times to the surface and then subtracts path integral times
// to get 3d travel times at points in the grid.  The grid geometry is defined by
// incident wave rays traced with a 1d reference model.  That is, they are like the
// S wavefield raygrid used in the main program of pwmig.
// 
// @param pstagrid parent pseudostation grid used to build this grid.
// @param border_pad is defines additional padding used to extend the pseudostation
//     grid on which the data are defined.  This is necessary (essential really) because
//     in pwmig this grid is used for all plane wave components while the S grids are
//     different for each plane wave component.  This should be made large enough that
//     it can contain all scattering points connected to S rays from the pseudostationh
//     grid.  
// @param UP3d P wave slowness defined with a 3d field object
// @param vp1d P wave 1d reference model.  Used to trace traces
// @param h hypocenter of source that created wavefield being imaged
// @param zmax ray trace parameters passed directly to BuildGCLraygrid
// @param tmax ray trace parameters passed directly to BuildGCLraygrid
// @param dt ray trace parameters passed directly to BuildGCLraygrid
//
//@returns auto_ptr object containing the travel time field
//@throws catches all exceptions and passes them along.  Several functions
//      called in this procedure could throw several possible exceptions I don't
//      know well which is the reason for the catch all syntax.
//
//@}
GCLscalarfield3d *ComputeIncidentWaveRaygrid(GCLgrid& pstagrid,
				int border_pad,
					GCLscalarfield3d& UP3d,
						VelocityModel_1d vp1d,
								Hypocenter h,
									double zmax,
									double tmax, 
									double dt)
{
	int i,j,k;
	// Incident wave slowness vector at grid origin
	// Always use 0 elevation datum for this calculation
	// note this isn't always used.  
	SlownessVector u0=h.pslow(pstagrid.lat0,pstagrid.lon0,0.0);
	// First set up the new expanded grid.  The new grid has border_pad cells
	// added on each side of the parent grid (pstagrid).  
	int n1new, n2new;
	int i0new, j0new;
	n1new = pstagrid.n1 + 2*border_pad;
	n2new = pstagrid.n2 + 2*border_pad;
	i0new = pstagrid.i0 + border_pad;
	j0new = pstagrid.j0 + border_pad;
	try {
		// this assumes the pstagrid is a "regular" gclgrid created with this same constructor.
		// Otherwise there may be small disconnect, although it should not really matter
		// provided dx1_nome and dx2_nom approximate the parent well.
		GCLgrid ng2d(n1new, n2new, string("ng2d"),
				pstagrid.lat0,pstagrid.lon0,pstagrid.r0,
				pstagrid.azimuth_y, pstagrid.dx1_nom, pstagrid.dx2_nom, i0new, j0new);
		// Use ng2d to compute a 1d reference model travel time grid
		// Note we explicitly do no include elevation in this calculation assuming this
		// will be handled with statics at another point
		GCLscalarfield Tp0(ng2d);
		for(i=0;i<Tp0.n1;++i)
		{
			for(j=0;j<Tp0.n2;++j)
			{
				Tp0.val[i][j]=h.ptime(Tp0.lat(i,j),Tp0.lon(i,j),0.0);
			}
		}
		auto_ptr<GCLgrid3d> IncidentRaygrid(Build_GCLraygrid(false, ng2d, u0,h,vp1d,zmax,tmax,dt));
		GCLscalarfield3d *Tp=new GCLscalarfield3d(*IncidentRaygrid);
		int kk, kend=Tp->n3 - 1;
		for(i=0;i<Tp->n1;++i)
		{
			for(j=0;j<Tp->n2;++j)
			{
				auto_ptr<dmatrix> path(extract_gridline(*IncidentRaygrid,
							i,j,kend,3,true));
				// this is the pwmig function using 1d if path is outside 3d model
				vector<double>Ptimes=compute_gridttime(UP3d,i,j,vp1d,*IncidentRaygrid,*path);
				for(k=0,kk=kend;k<=kend;++k,--kk)
					Tp->val[i][j][k] = Tp0.val[i][j] - Ptimes[kk];
			}
		}
		return Tp;
	} catch (...) {throw;}
}

// used in function below to flag an error
// has file scope so pwmig main can do a consistent check
const double TTERROR=-9999999999.0;
//@{
// computes lag time for P to S conversions in 3d.
// @param raygrid S plane wave component ray grid
// @param ix1 coordinate 1 index position of scattering point where the lag is to be computed
// @param ix2 coordinate 2 index position of scattering point where the lag is to be computed
// @param ix3 coordinate 3 index position of scattering point where the lag is to be computed
// @param TP scalarfield holding P wave travel times from source.  Assumed congruent with raygrid
// @param h location of source
//
//@}
double compute_scatter_point_Ptime(GCLgrid3d& raygrid, int ix1, int ix2, int ix3,
		GCLscalarfield3d& TP,Hypocenter& h)
{
	//This algorithm ASSUMES raygrid and TP are congruent. I avoids the overhead
	// of constant testing because it will be called a lot.
	// First get the P times to the pseudostation and and image points by 
	// calling the field lookup and interpolators.  
	double Tpx, Tpr;
	int err;
	double x1,x2,x3;
	// First compute P time to image point
	x1=raygrid.x1[ix1][ix2][ix3];
	x2=raygrid.x2[ix1][ix2][ix3];
	x3=raygrid.x3[ix1][ix2][ix3];
	err=TP.lookup(x1,x2,x3);
	switch (err)
	{
		case 0:
		case 2:
			Tpx=TP.interpolate(x1,x2,x3);
			break;
		default:
		// This condition is difficult to recover from.  The image point
		// depth can be too deep for many travel time calculators.  
		// I choose to not throw an exception in this case as exceptions
		// are a slow mechanism and this could happen often.
		// REturn a large negative number instead to signal this error.
		// The factor of 2 allows a direct test below against TTERROR and
		// is failsafe against small negative lags that could potentially arise.
			TP.reset_index();
			return(2.0*TTERROR);
	}
	return(Tpx);
}
double compute_receiver_Ptime(GCLgrid3d& raygrid, int ix1, int ix2,
		GCLscalarfield3d& TP,Hypocenter& h)
{
	double Tpr;  // result
	int err;
	double x1,x2,x3;
	int nsurface=raygrid.n3 - 1;
	x1=raygrid.x1[ix1][ix2][nsurface];
	x2=raygrid.x2[ix1][ix2][nsurface];
	x3=raygrid.x3[ix1][ix2][nsurface];
	err=TP.lookup(x1,x2,x3);
	switch (err)
	{
		case 0:
		// lookup almost always returns 2 in this case because receiver
		// is on top surface of grid.
		case 2:
			Tpr=TP.interpolate(x1,x2,x3);
			break;
		default:
		// Here we can recover from a lookup failure by reverting to the
		// 1d calculator
			TP.reset_index();
			Tpr=h.ptime(raygrid.lat(ix1,ix2,nsurface),
					raygrid.lon(ix1,ix2,nsurface),0.0);
	}
	return(Tpr);
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

	u0 - input slowness vector (note passed as SlownessVector object)
		u0 is base slowness for S raygrid passed as g (below)
	dux, duy - base slowness grid increments.
	vmod - 1d velocity model used for ray tracing 
	zmax, dz - ray trace parameters.  trace to depth zmax with depth increment dz
	g - gclgrid for S ray geometry
	ix1, ix2  - index position of receiver point in gclgrid g.
	gradP - parallel matrix to path of P time gradient vectors.
	zP - vector of depths parallel to gradP 
	(Note algorithm assumes gradP,zP are ordered from surface down)

Returns an np length STL vector object.

Author:  Gary Pavlis
Written:  NOvember 2003
*/
	

vector<double> compute_domega_for_path(SlownessVector& u0,double dux, double duy,
	VelocityModel_1d& vmod, 
		double zmax,
			double dz,
				GCLgrid3d& g,int ix1, int ix2,
					dmatrix& gradP,
						vector<double>& zP)
{
	SlownessVector udx1, udx2, udy1, udy2;
	dmatrix *pathdx1, *pathdx2, *pathdy1, *pathdy2;
	dmatrix *gradSdx1, *gradSdx2, *gradSdy1, *gradSdy2;
	int i;
	int npath; // used to define length of valid output 

	udx1 = u0;
	udx1.ux -= dux/2.0;
	udx2=u0;
	udx2.ux += dux/2.0;
	udy1 = u0;
	udy1.uy -= duy/2.0;
	udy2=u0;
	udy2.uy += duy/2.0;
	//
	// Now call the constructors for the basic ray paths for these four rays.  
	// We set the constructor in depth model and ASSUME zmax is small enough that
	// the ray for each u will not turn within vmod.  
	//
	RayPathSphere raydx1(vmod,udx1.mag(),zmax,1.0e99,dz,"z");
	RayPathSphere raydx2(vmod,udx2.mag(),zmax,1.0e99,dz,"z");
	RayPathSphere raydy1(vmod,udy1.mag(),zmax,1.0e99,dz,"z");
	RayPathSphere raydy2(vmod,udy2.mag(),zmax,1.0e99,dz,"z");
/*
mp.load(raydx1.r,string("rdx1"));
mp.load(raydx2.r,string("rdx2"));
mp.load(raydy1.r,string("rdy1"));
mp.load(raydy2.r,string("rdy2"));
mp.load(raydx1.delta,string("deltax1"));
mp.load(raydx2.delta,string("deltax2"));
mp.load(raydy1.delta,string("deltay1"));
mp.load(raydy2.delta,string("deltay2"));
mp.load(raydx1.t,string("tdx1"));
mp.load(raydx2.t,string("tdx2"));
mp.load(raydy1.t,string("tdy1"));
mp.load(raydy2.t,string("tdy2"));
*/
	//
	// project these into the GCLgrid coordinate system
	//
	try {
		pathdx1 = GCLgrid_Ray_project_down(g,raydx1, udx1.azimuth(),ix1,ix2,
				g.n3-1);
		pathdx2 = GCLgrid_Ray_project_down(g,raydx2, udx2.azimuth(),ix1,ix2,
				g.n3-1);
		pathdy1 = GCLgrid_Ray_project_down(g,raydy1, udy1.azimuth(),ix1,ix2,
				g.n3-1);
		pathdy2 = GCLgrid_Ray_project_down(g,raydy2, udy2.azimuth(),ix1,ix2,
				g.n3-1);
/*
mp.load(*pathdx1,string("pdx1"));
mp.load(*pathdx2,string("pdx2"));
mp.load(*pathdy1,string("pdy1"));
mp.load(*pathdy2,string("pdy2"));
mp.process(string("plotrays(pdx1,pdx2,pdy1,pdy2)"));
*/
	}  catch (GCLgrid_error err)
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
	gradSdx1=ray_path_tangent(*pathdx1);
	gradSdx2=ray_path_tangent(*pathdx2);
	gradSdy1=ray_path_tangent(*pathdy1);
	gradSdy2=ray_path_tangent(*pathdy2);
	delete pathdx1;  delete pathdx2;  delete pathdy1;  delete pathdy2;
	//
	// Check all these matrices have the same size and truncate the 
	// output if necessary writing a diagnostic in that situation.
	//
	npath = gradSdx1->columns();
	if(npath>gradSdx2->columns()) npath=gradSdx2->columns();
	if(npath>gradSdy1->columns()) npath=gradSdy1->columns();
	if(npath>gradSdy2->columns()) npath=gradSdy2->columns();
	if(npath!= gradSdx1->columns())
	{
		cerr << "compute_domega_for_path: irregular path length for bundle of 4 rays used to compute domega" << endl

			<< "Slowness grid probably defines turning rays in thie model" << endl;
	}
	for(i=0;i<npath;++i)
	{
		double slow;
		double depth;
		depth = raydx1.depth(i);
		slow = 1.0/vmod.getv(depth);
		dscal(3,slow,gradSdx1->get_address(0,i),1);
		dscal(3,slow,gradSdx2->get_address(0,i),1);
		dscal(3,slow,gradSdy1->get_address(0,i),1);
		dscal(3,slow,gradSdy2->get_address(0,i),1);
	}
/*
mp.load(gradP,string("gp"));
mp.load(*gradSdx1,string("gsdx1"));
mp.load(*gradSdx2,string("gsdx2"));
mp.load(*gradSdy1,string("gsdy1"));
mp.load(*gradSdy2,string("gsdy2"));
mp.process(string("figure; unitvectorplot(gsdx1);"));
mp.process(string("figure; unitvectorplot(gsdx2);"));
mp.process(string("figure; unitvectorplot(gsdy1);"));
mp.process(string("figure; unitvectorplot(gsdy2);"));
*/
	//
	// Now we interpolate the gradS vector arrays onto the gradP set
	// of depths. This produces a domega congruent with gradS and gradP
	// in the sense that the domega values will correspond to the same
	// node depths as those computed for gradS and gradP.  These can
	// be copied on return into the outputs, without any more interpolation
	//  albeit in the reverse order.
	//
	int npath_P;
	double z0=raydx1.depth(0);  // assume all four rays have same z0
	npath_P = gradP.columns();
	if(npath_P!=zP.size())
		throw SeisppError(string("compute_domega_for_path:  ")
		  + string("vector grid point depths and gradP matrix do not match"));
	dmatrix temp(3,npath_P); // Holds interpolated gradS values before replacement
	linear_vector_regular_to_irregular(z0,dz,*gradSdx1,&(zP[0]),temp);
	dmatrix nudx1(gradP+temp);  // creates vector gradP + gradS
	delete gradSdx1;  // Dont' need this any longer
	// Same for dx2, dy1, and dy2
	linear_vector_regular_to_irregular(z0,dz,*gradSdx2,&(zP[0]),temp);
	dmatrix nudx2(gradP+temp);  
	delete gradSdx2;  
	linear_vector_regular_to_irregular(z0,dz,*gradSdy1,&(zP[0]),temp);
	dmatrix nudy1(gradP+temp);  
	delete gradSdy1;  
	linear_vector_regular_to_irregular(z0,dz,*gradSdy2,&(zP[0]),temp);
	dmatrix nudy2(gradP+temp);  
	delete gradSdy2; 
	// normalize using antelope function d3norm to make these unit vectors
	for(i=0;i<npath_P;++i)
	{
		dr3norm(nudx1.get_address(0,i));
		dr3norm(nudx2.get_address(0,i));
		dr3norm(nudy1.get_address(0,i));
		dr3norm(nudy2.get_address(0,i));
	}
/*
mp.load(nudx1,string("nudx1"));
mp.load(nudx2,string("nudx2"));
mp.load(nudy1,string("nudy1"));
mp.load(nudy2,string("nudy2"));
mp.process(string("figure; unitvectorplot(nudx1);"));
mp.process(string("figure; unitvectorplot(nudx2);"));
mp.process(string("figure; unitvectorplot(nudy1);"));
mp.process(string("figure; unitvectorplot(nudy2);"));
mp.run_interactive();
*/
	//
	// Now get domega using cross products between pairs of unit vectors and 
	// small angle approximation (sin theta approx theta when theta is small)
	//
	vector<double>domega;
	for(i=0;i<npath_P;++i)
	{
		double cross[3];
		double dtheta_x, dtheta_y,domega_i;
		dr3cros(nudx1.get_address(0,i),nudx2.get_address(0,i),cross);
		dtheta_x = dr3mag(cross);
		dr3cros(nudy1.get_address(0,i),nudy2.get_address(0,i),cross);
		dtheta_y = dr3mag(cross);
		// abs shouldn't be necessary really, but better
		// safe than sorry given 0 cost
		domega_i=fabs(dtheta_x*dtheta_y);
		domega.push_back(domega_i);
	}
	return(domega);  
}
/* Much simpler routine to compute weight of this each member of a ray path
 * in generalized radon transform formula.  gradTp and gradTs are gradients
 * for P and S ray directions a defined in the Poppeliers and Pavlis (2003)
 * and Bostock et al. (2001).  Note this function drops the A term that
 * could be computed from geometric spreading because the assumption is
 * the this term is negligible for plane waves
 */
vector<double> compute_weight_for_path(dmatrix& gradTp,dmatrix& gradTs)
{
	const double eightpisq=78.95683521;
	double sum[3],nrmsum;
	int np=gradTp.columns();
	int i,j;
	vector<double> weight(np);

	for(j=0;j<np;++j)
	{
		for(i=0;i<3;++i)
		{
			sum[i]=gradTp(i,j)+gradTs(i,j);
		}
		nrmsum=dnrm2(3,sum,1);
		weight[j]=nrmsum*nrmsum*nrmsum/eightpisq;
	}
	return(weight);
}
/* Computing a running average of a specific length for a 
series of numbers stored in x.  npts is the number of points
in the smoother.  npts/2 on the left and right are padded as the 
constants for the first and last npts that can be averaged.
i.e. there is not endpoint tapering.  Return vector is x
smoothed as requested. 
*/
vector<double> running_average(vector<double>& x, int ns)
{
	int nx=x.size();
	int i,ii,j;
	double avg;
	vector<double> result(nx);

	if(nx<=ns)
	{
		for(i=0,avg=0.0;i<nx;++i)
			avg+=x[i];
		avg=avg/static_cast<double>(nx);
		for(i=0;i<nx;++i) result[i]=avg;
	}
	else
	{
		int npo2=(ns+1)/2;
		double avglen=static_cast<double>(ns);
		for(i=npo2-1;i<nx-(ns-npo2);++i)
		{
			for(j=i,ii=0,avg=0.0;ii<ns;++ii,++j)
				avg+=x[j-npo2+1];
			avg/=avglen;
			result[i]=avg;
		}
		// pad front and back
		for(i=npo2-2;i>=0;--i) result[i]=result[i+1];
		for(i=nx-(ns-npo2);i<nx;++i) result[i]=result[i-1];
	}
	return(result);
}
// builds dfile names for migratione output
string MakeDfileName(string base, int evid)
{
	string dfile;
	ostringstream sbuf(dfile);
	sbuf<<base<<"_"<<evid;
	return(sbuf.str());
}
// 
// We store velocity models externally in the database, but internally
// slowness is more appropriate.  This converts a velocity field to a 
// slowness field in an efficient way by altering the val array in place.
//
void VelocityFieldToSlowness(GCLscalarfield3d& g)
{
	for(int i=0;i<g.n1;++i)
		for(int j=0;j<g.n2;++j)
			for(int k=0;k<g.n3;++k) g.val[i][j][k]=1.0/g.val[i][j][k];
}
SlownessVector slowness_average(ThreeComponentEnsemble *d)
{
	double ux,uy;
	int n;
	n=0;
	ux=0.0;
	uy=0.0;
	for(int i=0;i<d->member.size();++i)
	{
		if(d->member[i].live)
		{
			ux+=d->member[i].get_double("ux");
			uy+=d->member[i].get_double("uy");
			++n;
		}
	}
	SlownessVector result;
	result.ux=ux/static_cast<double>(n);
	result.uy=uy/static_cast<double>(n);
	return(result);
}
/* Small helper retrieves hypocenter for current event.
This was made a procedure to isolate the oddities required
to manipulate the group pointer.  Done using datascope API
as there is no current way to do this in my generalized C++
database api.   I don't know enough about alternatives to
know how to generalize this.
*/
Hypocenter GetHypoFromDB(DatascopeHandle& dbh)
{
	int ierr;
	int is,ie;  // used for start and end record in a bundle
	Hypocenter hypo;
	Dbptr dbb;  // Bundle result
	const string base_emess("GetHypoFromDB:  ");
	// Get bundle pointer for lower level grouping 
	ierr=dbgetv(dbh.db,0,"bundle",&dbb,0);
	if(ierr==dbINVALID)
		throw SeisppError(base_emess
		  + string("dbgetv failed getting bundle pointer"));
	dbget_range(dbb,&is,&ie);
	// Assume one record is as good as another for next step
	dbb.record=is;
	ierr=dbgetv(dbb,0,"origin.lat",&(hypo.lat),
			"origin.lon",&(hypo.lon),
			"origin.depth",&(hypo.z),
			"origin.time",&(hypo.time),
			0);
	if(ierr==dbINVALID)
		throw SeisppError(base_emess
		  + string("dbgetv failed getting hypo parameters"));
	// db stores these in degrees, but we use radians internally
	hypo.lat=rad(hypo.lat);
	hypo.lon=rad(hypo.lon);
	return(hypo);
}
//DEBUG
// funtion using matlab engine to play a movie of a 3d scalar field
/*
void HorizontalSlicer(MatlabProcessor& p,GCLscalarfield3d& grid)
{
	dmatrix A(grid.n1,grid.n2);
	string name("A");
	string command("imagesc(A);pause(1);");
	int i,j,k;
	for(k=0;k<grid.n3;++k)
	{
		for(i=0;i<grid.n1;++i)
			for(j=0;j<grid.n2;++j)
				A(i,j)=grid.val[i][j][k];
		p.load(A,name);
		p.process(command);
	}
}
void SliceX1(MatlabProcessor& p,GCLscalarfield3d& grid)
{
	dmatrix A(grid.n3,grid.n1);
	string name("A");
        string command("imagesc(A);pause(1);");
        int i,j,k;
	for(j=0;j<grid.n2;++j)
	{
		for(i=0;i<grid.n1;++i)
			for(k=0;k<grid.n3;++k)
				A(k,i)=grid.val[i][j][k];
		p.load(A,name);
		p.process(command);
	}
}
*/


int main(int argc, char **argv)
{
        Dbptr db,dbv;
	string Pmodel1d_name;
	const string pvfnm("P"), svfnm("S");
	ThreeComponentEnsemble *pwdata;
	int is,i,j,k,l;
	int kk;
	string pfin("pwmig");
	dmatrix gradTs;
	Pf *pf;
	const string schema("pwmig1.1");
	int border_pad;

        ios::sync_with_stdio();
        elog_init(argc,argv);

        if(argc < 2) usage();
	string dbname(argv[1]);
	string dfile;  // used repeatedly below for data file names

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
        if(pfread(const_cast<char *>(pfin.c_str()),&pf)) die(1,(char *)"pfread error\n");
	try {
		Metadata control(pf);
		string velocity_grid_name=control.get_string("velocity_model_grid");
		string Pmodel3d_name=control.get_string("P_velocity_model3d_name");
		string Smodel3d_name=control.get_string("S_velocity_model3d_name");
		string Pmodel1d_name=control.get_string("P_velocity_model1d_name");
		string Smodel1d_name=control.get_string("S_velocity_model1d_name");
		border_pad = control.get_int("border_padding");
		string parent_grid_name=control.get_string("Parent_GCLgrid_Name");
		string stack_grid_name=control.get_string("stack_grid_name");
		string fielddir=control.get_string("output_field_directory");
		// Make sure this directory exists
		if(makedir(const_cast<char *>(fielddir.c_str())))
		{
			cerr << "Failure trying to create directory for output fields = "
				<< fielddir << endl;
			exit(-1);
		}
		string dfilebase=control.get_string("output_filename_base");
		int evid;  // output file name is build from dfilebase+"_"+evid
		bool use_depth_variable_transformation
			=control.get_bool("use_depth_variable_transformation");
		double zmax=control.get_double("maximum_depth");
		double tmax=control.get_double("maximum_time_lag");
		double dux=control.get_double("slowness_grid_deltau");
		double dt=control.get_double("data_sample_interval");
		double duy=dux;
		double dz=control.get_double("ray_trace_depth_increment");
		//
		// added Jan 2007 after recognition dweight and domega
		// were nearly invariant for a given slowness and 
		// were prone to creating artifacts at velocity discontinuities
		//
		bool rcomp_wt;
		rcomp_wt=control.get_bool("recompute_weight_functions");
		int nwtsmooth=control.get_int("weighting_function_smoother_length");
		bool smooth_wt;
		if(nwtsmooth<=0)
			smooth_wt=false;
		else
			smooth_wt=true;
		// Create a database handle for reading data objects
		dbopen(const_cast<char *>(dbname.c_str()),"r+",&db);
		if(db.record == dbINVALID)
		{
			cerr << "Cannot open database " << dbname <<endl;
			exit(-1);
		}
		// These constructors load a velocity model into a GCLgrid
		GCLscalarfield3d Up3d(db,velocity_grid_name,Pmodel3d_name);
		GCLscalarfield3d Us3d(db,velocity_grid_name,Smodel3d_name);
		// 
		// For this program we need to convert velocities to slowness.
		// This function called on P and S models does this in a procedural
		// way by altering the field variables in place.
		//
		VelocityFieldToSlowness(Up3d);
		VelocityFieldToSlowness(Us3d);
//DEBUG
//MatlabProcessor mp(string("zemlya"));
//MatlabProcessor mp(stdout);
/*
cout << "Display P velocity model horizontal slices"<<endl;
HorizontalSlicer(mp,Up3d);
cout << "Display S velocity model horizontal slices"<<endl;
sleep(5);
HorizontalSlicer(mp,Us3d);
*/
		VelocityModel_1d Vp1d(db,Pmodel1d_name,pvfnm);
		VelocityModel_1d Vs1d(db,Smodel1d_name,svfnm);
		GCLgrid parent(db,const_cast<char*>(parent_grid_name.c_str()) );

		// This loads the image volume assuming this was precomputed with
		// makegclgrid

		GCLvectorfield3d migrated_image(db,stack_grid_name,"",3);
		GCLscalarfield3d omega(db,stack_grid_name,"");
		GCLscalarfield3d weights(db,stack_grid_name,"");
		// 
		if(dynamic_cast<BasicGCLgrid&>(migrated_image)
			!= dynamic_cast<BasicGCLgrid&>(parent))
		{
			cerr << "pwmig(warning):  image grid and pseudostation grid are not congruent"
				<< endl
				<< "  This can increase run time"
				<< endl;
		}
		// Let's make sure these are initialized
		migrated_image.zero();
		omega.zero();
		weights.zero();
		//
		// silently remap the velocity model grids if necessary
		//
		if(dynamic_cast<BasicGCLgrid&>(Up3d)
			!= dynamic_cast<BasicGCLgrid&>(parent))
				remap_grid(dynamic_cast<GCLgrid3d&>(Up3d),
					dynamic_cast<BasicGCLgrid&>(parent));
		if(dynamic_cast<BasicGCLgrid&>(Us3d)
			!= dynamic_cast<BasicGCLgrid&>(parent))
				remap_grid(dynamic_cast<GCLgrid3d&>(Us3d),
					dynamic_cast<BasicGCLgrid&>(parent));
		

//YET ANOTHER NOTE.  HAVE TO DECIDE HOW TO HANDLE ZERO DATUM.  DO WE DO A GEOMETRIC STATIC
// AS IN REFLECTION PROCESSING?  PROBABLY SHOULD PUT THIS IN PWSTACK.  THAT IS STACK SHOULD
// SHOULD ALWAYS BE SHIFTED TO DATUM.
// Solution for now.  Modified pwstack to put an elev estimate for the pseudostation
// grid point computed as a weighted sum of elevations of contributing stations.  Algorithm
// below should use a vertical ray approximation to this static to datum or it will get
// unnecessarily complex.  A p=0 approximation for an elevation correction is pretty
// good for teleseismic data and is unlikely to be a big problem.
		// Load the default attribute map
		AttributeMap am(schema);
		// Need this set of attributes extraction lists to 
		// load correct metadata into the ensmble and each member
		MetadataList mdlin=pfget_mdlist(pf,"Station_mdlist");
		MetadataList mdens=pfget_mdlist(pf,"Ensemble_mdlist");
		// Outer loop over events.  This algorithm assumes
		// dbh is doubly grouped by evid:gridid and then by evid
		DatascopeHandle dbhens(db,pf,"dbprocess_commands");
		list<string> groupkey;
		groupkey.push_back(string("evid"));
		DatascopeHandle dbh(dbhens);
		dbh.group(groupkey);
		dbh.rewind();
		for(int record=0;record<dbh.number_tuples();++dbh,++record)
		{
			// First get the hypo for the current event
			// Made a procedure to insulate form database api dependence
			evid=dbh.get_int("evid");
			// Get the range for this group and set the db pointer
			// to first record.  This uses Datascope API directory
			// and hence is a portability issue for future use
			DBBundle dbb=dbh.get_range();
			dbhens.db.record=dbb.start_record;
			Hypocenter hypo=GetHypoFromDB(dbhens);
			// We build the P wave travel time field once for each
			// event as it is constant for all.  The border_pad parameter
			// is needed because plane waves come from outside the range
			// of the incident rays.  
			auto_ptr<GCLscalarfield3d> TPptr(ComputeIncidentWaveRaygrid(parent,
					border_pad,Up3d,Vp1d,hypo,zmax,tmax,dt));
			// Now loop over plane wave components
			for(dbhens.db.record=dbb.start_record;
				dbhens.db.record<dbb.end_record;++dbhens.db.record)
			{
			// first read in the next ensemble
			pwdata = new ThreeComponentEnsemble(dynamic_cast<DatabaseHandle&>
					(dbhens),mdlin,mdens,am);
#ifdef MATLABDEBUG
/*
string chans[3]={"x1","x2","x3"};
for(i=0;i<pwdata->member.size();++i)
{
int itest,jtest;
double smax=0;
for(jtest=0;jtest<3;++jtest) for(itest=0;itest<pwdata->member[i].ns;++itest) 
		smax=max(smax,pwdata->member[i].u(jtest,itest));
cerr << "smax="<<smax<<endl;
	mp.load(pwdata->member[i],"s");
	mp.process(string("subplot(3,1,1),plot(s(1,:));")
		+string("subplot(3,1,2),plot(s(2,:));")
		+string("subplot(3,1,3),plot(s(3,:));") );
	mp.process(string("pause(1);"));
}
mp.load(*pwdata,chans);
mp.process(string("figure;wigb(x1);figure;wigb(x2);figure;wigb(x3);"));
mp.run_interactive();
*/
#endif
cout << "DEBUG:  processing ensemble "<<record<<endl;
cout << "DEBUG: ensemble size="<<pwdata->member.size()<<endl;
int gridid=pwdata->get_int("gridid");
cout << "DEBUG:  this data is gridid="<<gridid<<endl;

			// A useful warning
			if(dt!=pwdata->member[0].dt)
			{
				cerr << "pwmig(WARNING):  evid="
					<< evid
					<< "sample rate mismatch. "<<endl
					<< "expected dt="<<dt
					<< "found dt="<<pwdata->member[0].dt
					<<endl
					<< "Data dt overrides dt from pf"<<endl;
				dt=pwdata->member[0].dt;
			}
			Ray_Transformation_Operator *troptr;
			// The data are assumed grouped in the ensemble by a common
			// Values do vary and we compute an average
			SlownessVector ustack=slowness_average(pwdata);
cout << "DEBUG: ux,uy="<<ustack.ux<<","<<ustack.uy<<endl;
/*
if(fabs(ustack.uy)>=0.007)
{
	cerr << "Skipping this component for debug"<<endl;
	continue;
}
*/
			string gridname;
			GCLgrid3d *grdptr;
			//
			// This builds a 3d grid using ustack slowness rays with
			// 2d grid parent as the surface grid array
			// use first member of ensemble to get dt
			grdptr =  Build_GCLraygrid(true,parent,ustack,hypo,
				Vs1d,zmax,tmax,dt);
			GCLgrid3d& raygrid=*grdptr;  // convenient shorthand because we keep changing this
			int n30;  // convenient since top surface is at n3-1
			n30 = raygrid.n3 - 1;
cout << "n3="<<raygrid.n3<<endl;
			// create work spaces for accumation of this component
			// This is a large memory model.  I'll use it until it proves
			// intractable
			// Below assumes the val arrays in these fields
			// have been initialized to zero.
			GCLscalarfield3d dweight(raygrid);
			GCLscalarfield3d domega(raygrid);
			GCLvectorfield3d pwdgrid(raygrid,3);

			// loop through the ensemble associating each trace with a 
			// point in the raygrid and interpolating to get time to depth.

			// 
			// this boolean is needed as part of optional setting 
			// recomputing weighting functions for each ray path.
			// When these weighting functions are computed only
			// once for each plane wave component (rcomp_wf false)
			// this seems the only safe and maintainable way to 
			// make sure these functions are always initialized.
			// They are way down in the code with too many 
			// potential errors that could cause a break condition
			// 
			bool weight_functions_set=false;
			// Warning:  n3 must not change in the loop below
			int n3=raygrid.n3;
			vector<double> domega_ij(n3),dweight_ij(n3);
			for(is=0;is<pwdata->member.size();++is)
			{
				// convenient shorthand variables.  ns is data length
				// while n3 is the ray path onto which data are mapped
				int ns=pwdata->member[is].ns;
				double zmaxray,tmaxray;
				vector<double> Stime(n3), SPtime(n3);
				double t0;

				dt=pwdata->member[is].dt;
				// first decide in which grid cell to place this seismogram 
				i = pwdata->member[is].get_int("ix1");
				j = pwdata->member[is].get_int("ix2");
				t0=pwdata->member[is].t0;
				if( (i<0) || (i>raygrid.n1) || (j<0) || (j>raygrid.n2) )
				{
					cerr << "grid index of data = ("
						<< i << ","<< j 
						<< "is outside range of grid." <<endl
						<< "Check grid definitions in database and pf"<<endl;
					exit(-1);
				}
				// We need a reference ray for the incident wavefield.
				// This ray is chopped up and transformed to GCLgrid coordinates
				// in the integration loop below.  We compute it here to 
				// avoid constantly recomputing it
				zmaxray = raygrid.depth(i,j,0);
				zmaxray *= 1.2;  // Small upward adjustment to avoid rubber banding
						// problems in 3d media
				tmaxray=10000.0;  // arbitrary large number.  depend on zmax 

				// Now compute a series of vectorish quantities that 
				// are parallel with the ray path associated with the 
				// current seismogram (i,j).  First we compute the S wave
				// travel time down this ray path (i.e. ordered from surface down)
				dmatrix *pathptr=extract_gridline(raygrid,i,j,raygrid.n3-1,3,true);
				Stime = compute_gridttime(Us3d,i,j,Vs1d,raygrid,*pathptr);
				// Now we compute the gradient in the S ray travel time
				// for each point on the ray.  Could have been done with 
				// less memory use in the loop below, but the simplification
				// it provides seems useful to me
				gradTs = compute_gradS(raygrid,i,j,Vs1d);
//DEBUG
#ifdef MATLABDEBUG
/*
mp.load(gradTs,string("gradts"));
mp.load(&(Stime[0]),Stime.size(),string("stime"));
cout << "Loaded gradts and stime"<<endl;
cout << "plotting gradts"<<endl;
mp.process(string("plot3c(gradts);"));
mp.process(string("pause(1)"));
mp.load(*pathptr,string("spath"));
cout << "plotting S wave path components for (i,j)="<<i<<","<<j;
mp.process(string("plot3c(spath);"));
*/
#endif
				dmatrix gradTp(3,n3);
				dmatrix nup(3,n3);
				vector<double> zP(n3);

				// Now loop over each point on this ray computing the 
				// p wave path and using it to compute the P time and
				// fill the vector gradP matrix
				// Note the double index, one running forward (k) and
				// a second (kk) running backward.  This complication is
				// caused by raygrid having upward oriented curves while
				// the S ray path (pathptr) is oriented downward.  
				bool tcompute_problem=false;
				// Only need to compute P time to surface once for 
				// each ray so we compute it before the loop below
				double Tpr=compute_receiver_Ptime(raygrid,i,j,*TPptr,hypo);
				for(k=0,kk=raygrid.n3-1;k<raygrid.n3;++k,--kk)
				{
					double tlag,Tpx;
					vector<double>nu;
					double vp;

					nu = compute_unit_normal_P(*TPptr,raygrid.x1[i][j][k],
						raygrid.x2[i][j][k], raygrid.x3[i][j][k]);
					for(l=0;l<3;++l) nup(l,k)=nu[l];
					// zP and gradTp are computed in reverse order
					zP[kk]=raygrid.depth(i,j,k);
					vp=Vp1d.getv(zP[kk]);
					for(l=0;l<3;++l) gradTp(l,kk)=nu[l]/vp;
//DEBUG

					// Compute time lag between incident P and scattered S
					Tpx=compute_scatter_point_Ptime(raygrid,i,j,k,*TPptr,hypo);
					tlag=Tpx+Stime[kk]-Tpr;
					if(tlag>=0)
						SPtime[kk]=tlag;
					else
					{
						if(tlag>-dt)
							SPtime[kk]=0.0;
						else
						{
cout << "tcompute_problem flag raised: tlag="<<tlag<<" at kk="<<kk<<endl;
							tcompute_problem=true;
							break;
						}
					}
				}
//DEBUG
/*
cout << "SP lag times\n";
for(l=0;l<n3;++l) cout << SPtime[l] << endl;
mp.load(&(SPtime[0]),SPtime.size(),string("sptime"));
mp.process(string("plot(sptime);"));
cout << "plotting gradTp"<<endl;
mp.load(gradTp,string("gradTp"));
mp.process(string("plot3c(gradtp);pause(1)"));
*/
				// skip this ray if there were travel time computation problems
				if(tcompute_problem)
				{
					cerr << "Warning:  ensemble "<< record 
						<< "grid position index ("
						<< i << ","
						<< j << ") Ray project dropped" << endl
						<< "Failure in computing lag.  "
						<< "Probably need to increase border_pad"<<endl;
					continue;

				}
				// We now interpolate the data with tlag values to map
				// the data from time to an absolute location in space
				// Note carefully SPtime and work sense is inverted from
				// the raygrid.  i.e. they are oriented for a ray from the surface
				// to the bottom.  Below we have to reverse this
				dmatrix work(3,SPtime.size());
				linear_vector_regular_to_irregular(t0,dt,pwdata->member[is].u,
					&(SPtime[0]),work);
#ifdef MATLABDEBUG
//DEBUG
/*
cout << "Input data number of samples="<<pwdata->member[is].ns<<endl;
cout << "DAta before interpolation" << endl;
for(int ifoo=0;ifoo<pwdata->member[is].ns;++ifoo)
	cout << pwdata->member[is].time(ifoo) << " "
		<< pwdata->member[is].u(0,ifoo) << " "
		<< pwdata->member[is].u(1,ifoo) << " "
		<< pwdata->member[is].u(2,ifoo) << endl;
*/
mp.load(pwdata->member[is],string("de"));
mp.load(&(SPtime[0]),SPtime.size(),string("sptime"));
mp.load(work,string("di"));
cerr << "Loaded sptime and interpolated data (work matrix)"<<endl;
mp.process(string("plot6c(de,sptime,di);pause(1.0);"));
#endif

				// Compute the transformation vector for each 
				// ray point in space.  This comes in two flavors.
				// A fast, inexact version and a exact slow version
				//
//DEBUG NOTE:  check that this is correctly oriented.  i.e. order of the work array
// is inverted compared to raygrid 
				if(use_depth_variable_transformation)
				{
					// the slow, exact method
					troptr=new Ray_Transformation_Operator(
						parent,
						*pathptr,
						ustack.azimuth(),
						nup);
				}
				else
				{
					// the fast approximation
					troptr=new Ray_Transformation_Operator(
						parent,
						*pathptr,ustack.azimuth());
				}
//DEBUG
//mp.load(work,string("work"));
//cerr << "work matrix loaded.  Before ray transformation"<endl;
//mp.process(string("plot3c(work);"));
				Ray_Transformation_Operator& trans_operator=*troptr;
				work=trans_operator.apply(work);
/*
mp.load(work,string("work2"));
cerr << "work matrix loaded as work2.  After ray transformaton"<endl;
//mp.process(string("figure; plot3c(work);"));
//mp.run_interactive();
mp.process(string("plot3c6(work,work2)"));
*/
				// done with these now
				delete pathptr;
				delete troptr;

				// This computes domega for a ray path in constant dz mode
				// We would have to interpolate anyway to mesh with 
				// time data choose the faster, coarser dz method here
				//
				if(rcomp_wt || !weight_functions_set)
				{
				  domega_ij=compute_domega_for_path(ustack,dux,duy,
					Vs1d, zmaxray,dz,
					raygrid, i, j, gradTp,zP);
				  dweight_ij=compute_weight_for_path(gradTp,gradTs);
// DEBUG SECTION
/*
cout << "Plotting dweight and domega for (i,j)="<<i<<","<<j<<endl;
mp.load(&(domega_ij[0]),domega_ij.size(),string("domega"));
mp.load(&(dweight_ij[0]),dweight_ij.size(),string("dweight"));
mp.process(string("plotow;figure;"));
*/
				  if(smooth_wt)
				  {
					domega_ij=running_average(domega_ij,nwtsmooth);
					dweight_ij=running_average(dweight_ij,nwtsmooth);
				  }
/*
mp.load(&(domega_ij[0]),domega_ij.size(),string("domega"));
mp.load(&(dweight_ij[0]),dweight_ij.size(),string("dweight"));
mp.process(string("plotow;title 'after smoother';"));
*/
				  weight_functions_set=true;
				}
					
				//
				// copy transformed data to vector field
				// copy weights and domega at same time
				//

				for(k=0,kk=raygrid.n3-1;k<raygrid.n3;++k,--kk)
				{
					for(l=0;l<3;++l)
					{
						// DEBUG
						//pwdgrid.val[i][j][k][l]=work(l,kk);
						pwdgrid.val[i][j][k][l]=work(l,kk)
							*dweight_ij[kk]*domega_ij[kk];
					}
					dweight.val[i][j][k]=dweight_ij[kk];
					domega.val[i][j][k]=domega_ij[kk];
				}
#ifdef MATLABDEBUG
//DEBUG
/*
dmatrix dproj(3,n3);
for(k=0;k<n3;++k)
	for(l=0;l<3;++l) dproj(l,k)=pwdgrid.val[i][j][k][l];
mp.load(dproj,string("d"));
mp.process(string("plot3c(d);pause(0.1);"));
*/
#endif
			}
			delete pwdata;
/*
cout << "plane wave data component" << endl;
bool debugexit;
debugexit=true;
//cout << pwdgrid;
GCLscalarfield3d *sfptr;
sfptr=extract_component(pwdgrid,0);
//output_gcl3d_to_vtksg(*sfptr,"component0.vtk");
SliceX1(mp,*sfptr);
delete sfptr;
sfptr=extract_component(pwdgrid,1);
//output_gcl3d_to_vtksg(*sfptr,"component1.vtk");
SliceX1(mp,*sfptr);
delete sfptr;
sfptr=extract_component(pwdgrid,2);
//output_gcl3d_to_vtksg(*sfptr,"component2.vtk");
SliceX1(mp,*sfptr);
delete sfptr;
*/
			// last but not least, add this component to the stack
			//
			migrated_image += pwdgrid;
			weights += dweight;
			omega += domega;
//DEBUG  save partial sums
/*
                dfile=MakeDfileName(dfilebase+string("_psum"),gridid+1000);
                migrated_image.dbsave(db,"",fielddir,dfile,dfile);
*/


			delete grdptr;
//if(debugexit) exit(1);

		} // Bottom of plane wave component loop
/*
mp.run_interactive();
GCLscalarfield3d *sfptr;
cout << "Display x1 component with matlab"<<endl;
sfptr=extract_component(migrated_image,0);
SliceX1(mp,*sfptr);
delete sfptr;
cout << "Display x2 component with matlab"<<endl;
sfptr=extract_component(migrated_image,1);
SliceX1(mp,*sfptr);
delete sfptr;
cout << "Display x3 component with matlab"<<endl;
sfptr=extract_component(migrated_image,2);
SliceX1(mp,*sfptr);
delete sfptr;
cout << "Display weights with matlab"<<endl;
SliceX1(mp,weights);
cout << "Display omega volume with matlab"<<endl;
SliceX1(mp,omega);
*/


		//
		// Save results
		// Second arg as empty string causes parent grid 
		// to not be saved and only the migrated data are saved
		// In all cases we force the field name in the output to
		// be the same as the output data file name, dfile
		// dfile is generated as base+tag+evid
		//

		dfile=MakeDfileName(dfilebase+string("_data"),evid);
		migrated_image.dbsave(db,"",fielddir,dfile,dfile);
		dfile=MakeDfileName(dfilebase+string("_wt"),evid);
		weights.dbsave(db,"",fielddir,dfile,dfile);
		dfile=MakeDfileName(dfilebase+string("_omega"),evid);
		omega.dbsave(db,"",fielddir,dfile,dfile);
		//
		// zero field values so when we loop back they can
		// be reused
		//
		migrated_image.zero();
		weights.zero();
		omega.zero();
		} // bottom of event loop
	}
	catch (GCLgrid_error gcle)
	{
		gcle.log_error();
		die(1,"GCLgrid library fatal error\n");
	}
	catch (MetadataError mde)
	{
		mde.log_error();
		die(1,"Required processing metadata missing from input stream\n");
	}
	catch (SeisppError seer)
	{
		seer.log_error();
		die(1,"seispp library function error\n");
	}
	catch (int ierr)
	{
		die(1,"Something threw a simple int exception of %d\n",
			ierr);
	}
	catch (...)
	{
		die(1,"Unhandled exception was thrown\n");
	}
}


