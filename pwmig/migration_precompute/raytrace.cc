#include <iostream>
#include <list>
#include "db.h"
#include "vmodel.h"
#include "ray1d.h"
#include "dmatrix.h"
#include "gclgrid.h"

/*
 * Traces a ray and returns result in RayPathSphere object 
 * for velocity model vmod and ray parameter u.  vmod is
 * assumed to have units of s/km and u is assumed to have
 * units of s/km.  Internally s/radians are used, but 
 * to use s/rad externally would be very confusing.
 * When mode=="z" aims for equal z intervals of size del.
 * When mode=="t" aims for equal time intervals of size del.
 * Ray is traced to depth zmax.  
 *
 * Note that we have to handle the common situation of a turning
 * ray.  That is, if zmax is below the turning point of a ray
 * with the requested ray parameter we handle this specially because
 * if we did not the travel time formulae contain sqrt of negative
 * numbers.  This is handled in a simple way here by simply stopping
 * the ray tracing if the next point down goes singular.  Whether 
 * this is important to the caller depends on the use.  Rather than
 * throw an exception, which isn't necessary really, the user is 
 * simply warned that they should test the actual output's maximum
 * depth compared to zmax when this is a problem.
 *
 * A similar feature is implemented when running in equal time
 * step mode.  This could do nasty things too if the requested
 * zmax was very close to a turning point.  It is possible in
 * this situation that the function could try to compute a vector
 * that is unbounded because the dz step length for each dt 
 * gets tiny.  To avoid this we also include the tmax parameter.
 * The ray is truncated when either z>zmax or t>tmax.  If you 
 * are sure the ray won't bottom for the requested p simply make
 * tmax large.
 */

RayPathSphere *Trace_Ray_Sphere(Velocity_Model& vmod,
		double u, double zmax, double tmax,
		double del, string mode)
{
	list <double> ddlist;
	list <double> dzlist;
	list <double> dtlist;
	list <double>::iterator id, iz, it;
	double t,z,x;  // accumulated time, depth and distance
	double r;  // radius now = R0-z
	const double R0=6378.17;  
	double ddelta,dz,dt;
	double vel;
	int npoints;
	double p; 
	int i;

	p = u*R0; // this ray parameter will have units of s/radian

	t=0.0;
	z=0.0;

	while(z<zmax && t<tmax)
	{
		// see Lay and Wallace section on spherical earth
		// travel time claculations for formulae used here
		double eta;  
		double root_term;  
		vel=vmod.getv(z);
		r=R0-z;
		eta=r/vel; //
		if(p>eta)break;  // exit loop to avoid a negative sqrt
		root_term=sqrt(eta*eta-p*p);
		// mode switches between equal depth steps and equal
		// time steps
		if(mode=="z")
		{
			dz=del;
			dt=eta*eta*dz/(r*root_term);
		}
		else
		{
			dt=del;
			dz=r*r*root_term*dt/(eta*eta);
		}
		// Note this is in radians
		ddelta=p*dz/(r*root_term);
		t+=dt;
		z+=dz;
		ddlist.push_back(ddelta);
		dtlist.push_back(dt);
		dzlist.push_back(dz);
	}
	// all the lists have to be the same size so we randomly
	// choose one to compute the size of the output vector
	// Is +1 because this is segments versus points
	npoints=ddlist.size()+1;
	RayPathSphere *path=new RayPathSphere(npoints);
	// Now we build the output structure by simple accumulation
	// The lists will be delted on exit and provided a useful
	// alternative to realloc.  Because this lists are parallel
	// we loop over one and just passively increment the iterators
	// on the others.  This could have been done with a list of a
	// temporary structure, but that's a judgment call on memory
	// use versus clarity.
	path->r[0]=R0;
	path->t[0]=0.0;
	path->delta[0]=0.0;
	for(i=1,it=dtlist.begin(),id=ddlist.begin(),iz=dzlist.begin();
			it!=dtlist.end();++it,++id,++iz,++i) 
	{
		path->r[i]=path->r[i-1] - (*iz);
		path->t[i]=path->t[i-1] + (*it);
		path->delta[i]=path->delta[i-1] + (*id);
	}
	// Don't need to set path->npts assuming constructor set it
	
	return(path);
}
/*  This function takes an input of a ray path in spherical coordinates
 *  defined by a RayPathSphere object.  This is assumed to be spherical
 *  coordinate definitions of the path in a radius-distance (radians)
 *  reference frame.  This path is translated to the origin of the
 *  (input)GLCgrid object and projected along the positive x1 baseline
 *  of this GLCgrid.  Note this tacitly assumes the GCLgrid passed
 *  is a "standard" grid and not a more general GCLgrid object allowed
 *  by the definition of a GCLgrid.  
 *
 *  The output dmatrix contains cartesian vectors in the GCLgrid 
 *  reference frame for this path.  Note this means x2 coordinates
 *  are always machine zeros because the path is projected parallel
 *  to the x1 baseline.  
 *
 *  Note this function maybe should return a pointer to avoid
 *  copying this potentially fairly large object.
 *
 *  Author:  Gary L. Pavlis
*/
/*
dmatrix GCLgrid_Ray_gtoc(GCLgrid& grid, RayPathSphere& path)
{
	dmatrix pathout(3,path.npts);  // the output object
	double azimuth;
	Cartesian_point this_point;

	// GCLgrid stores azimuth y, we want azimuth x (x1)
	azimuth = grid.azimuth_y + M_PI_2;

	for(int i=0;i<path.npts;++i)
	{
		double lat,lon;
		latlon(grid.lat0,grid.lon0,path.r[i],azimuth,&lat,&lon);
		this_point = grid.gtoc(lat,lon,r[i]);
		pathout(0,i)=this_point.x1;
		pathout(1,i)=this_point.x2;
		pathout(2,i)=this_point.x3;
	}
	return(pathout);
}

*/
