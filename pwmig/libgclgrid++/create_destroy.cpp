#include <string.h>
#include "stock.h"
#include "coords.h"
#include "gclgrid.h"
/* These are create and free routines for 2d and 3d GCLgrid objects in plain C.
They are used below for C++ constructors and destructors.  The MOST 
IMPORTANT thing to note about the constructions used here is that I have
implemented these built on contiguous blocks of memory.  This is NOT necessarily
the way one would always build 2D and 3D arrays in C.  This is probably a bad
idea for future maintainability of this code as someone might be tempted to
write their own destructor or add a constructor that used a different memory model.
Future users beware on this point.  */


/* Many gcl functions work with 3d C grids.  C implements multidimensioned
arrays pointer arrays, which makes this kind of ugly.  Anyway, this
function will create the F90 equivalent of  x[0:n1-1][0:n2-1][0:n3-1]
and return a ***double pointer for the 3d array.  The approach 
used here allocates the full block of memory in one single call to 
guarantee the result is in a contiguous block of memory. 
This is useful for reading and writing this type of entity because
then a single call to fread and fwrite can be used.  
Note, however, that the actual order in memory remains the C convention
not the FOTRAN form.  i.e. the rightmost index defines members contingous
in memory.
Author:  GAry Pavlis
Written:  June 2000
*/

double ***create_3dgrid_contiguous(int n1, int n2, int n3)
{
	double ***ptr3d;
	double **ptr2ptr;
	double *ptr;
	int i,j;

	allot(double *,ptr,n1*n2*n3);
	allot(double ***,ptr3d,n1);
	for(i=0;i<n1;++i)
	{
		allot(double **,ptr2ptr,n2);
		ptr3d[i] = ptr2ptr;
	}
	for(i=0;i<n1;++i)
	{
		for(j=0;j<n2;++j)
		{
			ptr3d[i][j] = ptr + n2*n3*i + n3*j;
		}
	}
	return(ptr3d);
}		
/* Free routine for a 3d array.  pointers arrays require exra attention
to avoid a memory leak.   This function assumes standard C indexing
*/
void free_3dgrid_contiguous(double ***x,int n1, int n2)
{
	int i;
	double *ptr;

	/* this clears the work space*/
	ptr = x[0][0];
	free((void *)ptr);

	/* The pointer arrays are still then and have to be freed 
	seperately */
	for(i=0;i<n1;++i)  free((void *)x[i]);
	free((void *)x);
}
/* parallel routines to the above for 2d */
double **create_2dgrid_contiguous(int n1, int n2)
{
	double **ptr2ptr;
	double *ptr;
	int i;

	allot(double *,ptr,n1*n2);
	allot(double **,ptr2ptr,n1);
	for(i=0;i<n1;++i)
		ptr2ptr[i] = ptr + n2*i;
	return(ptr2ptr);
}
void free_2dgrid_contiguous(double **x,int n1)
{
	double *ptr;

	ptr = x[0];
	free((void *)ptr);
	free((void *)x);
}
//
// C++ constructors 
//
GCLgrid::GCLgrid (int n1size, int n2size)
{
	name[0]='\0';
	lat0=0.0; lon0=0.0; r0=0.0;
	azimuth_y=0.0;  dx1_nom=0.0;  dx2_nom=0.0;
	i0=0;  j0=0;
	n1=n1size;
	n2=n2size;
	x1=create_2dgrid_contiguous(n1size,n2size);
	x2=create_2dgrid_contiguous(n1size,n2size);
	x3=create_2dgrid_contiguous(n1size,n2size);
	lat=create_2dgrid_contiguous(n1size,n2size);
	lon=create_2dgrid_contiguous(n1size,n2size);
	r=create_2dgrid_contiguous(n1size,n2size);
	cartesian_defined=1;
	geographic_defined=1;
}

GCLgrid3d::GCLgrid3d (int n1size, int n2size, int n3size)
{
	name[0]='\0';
	lat0=0.0; lon0=0.0; r0=0.0;
	azimuth_y=0.0;  dx1_nom=0.0;  dx2_nom=0.0; dx3_nom=0.0;
	i0=0;  j0=0;  k0=0;
	n1=n1size;
	n2=n2size;
	n3=n3size;
	x1=create_3dgrid_contiguous(n1size,n2size,n3size);
	x2=create_3dgrid_contiguous(n1size,n2size,n3size);
	x3=create_3dgrid_contiguous(n1size,n2size,n3size);
	lat=create_3dgrid_contiguous(n1size,n2size,n3size);
	lon=create_3dgrid_contiguous(n1size,n2size,n3size);
	r=create_3dgrid_contiguous(n1size,n2size,n3size);
	cartesian_defined=1;
	geographic_defined=1;
}
//  Copy constructor is far more complex and here it is
GCLgrid::GCLgrid(const GCLgrid& g)
{
	int i,j;
	//
	//I think the books officially say copying the scalars defaulted
	//I, however, have found depending on subtle standard features
	//dangerous so  I'll do this explicitly.
	//
	strcpy(name, g.name);
	lat0=g.lat0;
	lon0=g.lon0;
	r0=g.r0;
	azimuth_y=g.azimuth_y;
	dx1_nom=g.dx1_nom;
	dx2_nom=g.dx2_nom;
	n1=g.n1;
	n2=g.n2;
	i0=g.i0;
	j0=g.j0;
	xlow=g.xlow;
	xhigh=g.xhigh;
	ylow=g.ylow;
	yhigh=g.yhigh;
	zlow=g.zlow;
	zhigh=g.zhigh;
	x1=create_2dgrid_contiguous(n1,n2);
	x2=create_2dgrid_contiguous(n1,n2);
	x3=create_2dgrid_contiguous(n1,n2);
	lat=create_2dgrid_contiguous(n1,n2);
	lon=create_2dgrid_contiguous(n1,n2);
	r=create_2dgrid_contiguous(n1,n2);
	//
	//I use separate loops for each array here as this is highly
	//optimized on most compilers 
	//
	for(i=0;i<n1;++i)
		for(j=0;j<n2;++j) x1[i][j]=g.x1[i][j];
	for(i=0;i<n1;++i)
		for(j=0;j<n2;++j) x2[i][j]=g.x2[i][j];
	for(i=0;i<n1;++i)
		for(j=0;j<n2;++j) x3[i][j]=g.x3[i][j];
	for(i=0;i<n1;++i)
		for(j=0;j<n2;++j) lat[i][j]=g.lat[i][j];
	for(i=0;i<n1;++i)
		for(j=0;j<n2;++j) lon[i][j]=g.lon[i][j];
	for(i=0;i<n1;++i)
		for(j=0;j<n2;++j) r[i][j]=g.r[i][j];

}
GCLgrid3d::GCLgrid3d(const GCLgrid3d& g)
{
	int i,j,k;
	strcpy(name,g.name);
	lat0=g.lat0;
	lon0=g.lon0;
	r0=g.r0;
	azimuth_y=g.azimuth_y;
	dx1_nom=g.dx1_nom;
	dx2_nom=g.dx2_nom;
	dx3_nom=g.dx3_nom;
	n1=g.n1;
	n2=g.n2;
	n3=g.n3;
	i0=g.i0;
	j0=g.j0;
	xlow=g.xlow;
	xhigh=g.xhigh;
	ylow=g.ylow;
	yhigh=g.yhigh;
	zlow=g.zlow;
	zhigh=g.zhigh;
	x1=create_3dgrid_contiguous(n1,n2,n3);
	x2=create_3dgrid_contiguous(n1,n2,n3);
	x3=create_3dgrid_contiguous(n1,n2,n3);
	lat=create_3dgrid_contiguous(n1,n2,n3);
	lon=create_3dgrid_contiguous(n1,n2,n3);
	r=create_3dgrid_contiguous(n1,n2,n3);
	for(i=0;i<n1;++i)
		for(j=0;j<n2;++j) 
			for(k=0;k<n3;++k) x1[i][j]=g.x1[i][j];
	for(i=0;i<n1;++i)
		for(j=0;j<n2;++j) 
			for(k=0;k<n3;++k) x2[i][j]=g.x2[i][j];
	for(i=0;i<n1;++i)
		for(j=0;j<n2;++j) 
			for(k=0;k<n3;++k) x3[i][j]=g.x3[i][j];
	for(i=0;i<n1;++i)
		for(j=0;j<n2;++j) 
			for(k=0;k<n3;++k) lat[i][j]=g.lat[i][j];
	for(i=0;i<n1;++i)
		for(j=0;j<n2;++j) 
			for(k=0;k<n3;++k) lon[i][j]=g.lon[i][j];
	for(i=0;i<n1;++i)
		for(j=0;j<n2;++j) 
			for(k=0;k<n3;++k) r[i][j]=g.r[i][j];

}
//
//  assignment operator is almost like the copy constructor but 
//  we protect against self assignment.
//
GCLgrid& GCLgrid::operator=(const GCLgrid& g)
{
	if(this != &g)  // avoid self assignment
	{
		int i,j;
		strcpy(name,g.name);
		lat0=g.lat0;
		lon0=g.lon0;
		r0=g.r0;
		azimuth_y=g.azimuth_y;
		dx1_nom=g.dx1_nom;
		dx2_nom=g.dx2_nom;
		n1=g.n1;
		n2=g.n2;
		i0=g.i0;
		j0=g.j0;
		xlow=g.xlow;
		xhigh=g.xhigh;
		ylow=g.ylow;
		yhigh=g.yhigh;
		zlow=g.zlow;
		zhigh=g.zhigh;
		x1=create_2dgrid_contiguous(n1,n2);
		x2=create_2dgrid_contiguous(n1,n2);
		x3=create_2dgrid_contiguous(n1,n2);
		lat=create_2dgrid_contiguous(n1,n2);
		lon=create_2dgrid_contiguous(n1,n2);
		r=create_2dgrid_contiguous(n1,n2);
		//
		//I use separate loops for each array here as this is highly
		//optimized on most compilers 
		//
		for(i=0;i<n1;++i)
			for(j=0;j<n2;++j) x1[i][j]=g.x1[i][j];
		for(i=0;i<n1;++i)
			for(j=0;j<n2;++j) x2[i][j]=g.x2[i][j];
		for(i=0;i<n1;++i)
			for(j=0;j<n2;++j) x3[i][j]=g.x3[i][j];
		for(i=0;i<n1;++i)
			for(j=0;j<n2;++j) lat[i][j]=g.lat[i][j];
		for(i=0;i<n1;++i)
			for(j=0;j<n2;++j) lon[i][j]=g.lon[i][j];
		for(i=0;i<n1;++i)
			for(j=0;j<n2;++j) r[i][j]=g.r[i][j];
	
	}
	return *this;
}
//
//Same for 3D grid
//
GCLgrid3d& GCLgrid3d::operator=(const GCLgrid3d& g)
{
	if(this != &g)  // avoid self assignment
	{
		int i,j,k;
		strcpy(name,g.name);
		lat0=g.lat0;
		lon0=g.lon0;
		r0=g.r0;
		azimuth_y=g.azimuth_y;
		dx1_nom=g.dx1_nom;
		dx2_nom=g.dx2_nom;
		dx3_nom=g.dx3_nom;
		n1=g.n1;
		n2=g.n2;
		n3=g.n3;
		i0=g.i0;
		j0=g.j0;
		xlow=g.xlow;
		xhigh=g.xhigh;
		ylow=g.ylow;
		yhigh=g.yhigh;
		zlow=g.zlow;
		zhigh=g.zhigh;
		x1=create_3dgrid_contiguous(n1,n2,n3);
		x2=create_3dgrid_contiguous(n1,n2,n3);
		x3=create_3dgrid_contiguous(n1,n2,n3);
		lat=create_3dgrid_contiguous(n1,n2,n3);
		lon=create_3dgrid_contiguous(n1,n2,n3);
		r=create_3dgrid_contiguous(n1,n2,n3);
		for(i=0;i<n1;++i)
			for(j=0;j<n2;++j) 
				for(k=0;k<n3;++k) x1[i][j]=g.x1[i][j];
		for(i=0;i<n1;++i)
			for(j=0;j<n2;++j) 
				for(k=0;k<n3;++k) x2[i][j]=g.x2[i][j];
		for(i=0;i<n1;++i)
			for(j=0;j<n2;++j) 
				for(k=0;k<n3;++k) x3[i][j]=g.x3[i][j];
		for(i=0;i<n1;++i)
			for(j=0;j<n2;++j) 
				for(k=0;k<n3;++k) lat[i][j]=g.lat[i][j];
		for(i=0;i<n1;++i)
			for(j=0;j<n2;++j) 
				for(k=0;k<n3;++k) lon[i][j]=g.lon[i][j];
		for(i=0;i<n1;++i)
			for(j=0;j<n2;++j) 
				for(k=0;k<n3;++k) r[i][j]=g.r[i][j];
	
	}
	return *this;
}

/* This small function builds a baseline vector of latitude and
longitude values along a specified azimuth with a prescribed 
origin.

Arguments:
	lat0, lon0 - origin in radians
	phi - angle to rotate the x axis (+ east) baseline.  That is
		the baseline is established as a great circle passing
		through the origin point with an angle phi measured positive
		northward from east.  (counterclock looking from above)
		Note that away from the origin this angle relative to north
		will change.
	dx - distance spacing between points on the baseline (in radians)
	n - number of points along the baseline.
	i0 - position of the origin point along the baseline vector. 
		NOTE:  we use C convention so 0 is the first point in
		the vector.
	lat, lon - hold n vectors of latitude, longitude pairs of the 
		baseline points.  Note that lat[i0], lon[i0] should 
		equal lat0, lon0 respectively.
Author:  G Pavlis
Written:  Aug 2000
*/

void build_baseline(double lat0, double lon0, double phi, 
			double dx,int n,int i0,
			double *lat, double *lon)
{
	double azimuth,az;
	int i;
	double delta;
	
	/* We want azimuth between 0 and 2*pi */
	azimuth = M_PI_2 - phi;
	if(azimuth < 0.0 ) azimuth += (2.0*M_PI);

	for(i=0,delta=dx*((double)-i0);i<n;++i,delta+=dx)
	{
		if(i<i0)
			latlon(lat0,lon0,fabs(delta),azimuth-M_PI,lat+i,lon+i);
		else if(i==i0)
		{
			lat[i] = lat0;
			lon[i] = lon0;
		}
		else
		{
			latlon(lat0,lon0,delta,azimuth,lat+i,lon+i);
		}
	}
}
//
// This is the C++ constructor that uses a pf and coord routines to compute
// a GCLgrid on the fly using the description coming from pf
//
GCLgrid::GCLgrid(int n1in, int n2in, 
	char *namein, 
	double lat0in, double lon0in, double r0in,
	double azimuth_yin, double dx1_nomin, double dx2_nomin, 
	int i0in, int j0in)

{		

	/* pole to baseline */
	double pole_lat, pole_lon;
	int i,j,k;
	double deltax, deltay, delta;  
	double z0,z;
	double x[3],x0[3];
	double rmatrix[9],rmtrans[9];
	double xwork[3],xwork2[3];
	double rotation_angle,azimuth_i;
	double dx1_rad,dx2_rad;

	strcpy(name,namein);
	lat0=lat0in;
	lon0=lon0in;
	lat0 = rad(lat0);
	lon0 = rad(lon0);
	r0 = r0_ellipse(lat0);
	azimuth_y=azimuth_yin;
	rotation_angle=-azimuth_y;
	dx1_nom=dx1_nomin;
	dx2_nom=dx2_nomin;
	dx1_rad = dx1_nom/r0;
	dx2_rad = dx2_nom/r0;
	n1=n1in;
	n2=n2in;
	i0=i0in;
	j0=j0in;

	/* These are the vector of lat and lons used to define the
	baseline the effectively forms the equator of of the rotated
	coordinates passing through the origin*/
	double *baseline_lat=new double[n1in];
	double *baseline_lon=new double[n1in];

	build_baseline(lat0,lon0,
		rotation_angle,dx1_rad,n1,i0,baseline_lat,baseline_lon);
	/* We need to compute the pole to our base line to use as a 
	target for grid lines that are locally perpendicular along
	the grid lines -- like longitude lines at the equator */
	latlon(lat0,lon0,M_PI_2,-rotation_angle,&pole_lat,&pole_lon);

	/* We now allocate memory for the large coordinate arrays
	themselves and set the elements of the structure */
	lat = create_2dgrid_contiguous(n1,n2);
	lon = create_2dgrid_contiguous(n1,n2);
	r = create_2dgrid_contiguous(n1,n2);
	x1 = create_2dgrid_contiguous(n1,n2);
	x2 = create_2dgrid_contiguous(n1,n2);
	x3 = create_2dgrid_contiguous(n1,n2);
	cartesian_defined = 1;
	geographic_defined = 1;

	/* We now complete the latitude/longitude grid by projecting
	lines from the baseline toward the computed pole.  Grid points 
	will be equally spaced along the gcp toward the pole, but the
	lines they form will converge.  This is complicated greatly by
	the grid origin parameters.  We could make this a function,
	but it is simpler to just write it inline.*/
	for(i=0;i<n1;++i)
	{
		dist(baseline_lat[i],baseline_lon[i],pole_lat, pole_lon,
			&delta,&azimuth_i);
		if(azimuth_i < 0.0) azimuth_i += (2.0*M_PI);
		delta = dx2_rad*(-(double)j0);
		for(j=0;j<n2;++j,delta+=dx2_rad)
		{
			if(j<j0)
				latlon(baseline_lat[i],baseline_lon[i],
					fabs(delta),azimuth_i-M_PI,
					lat[i]+j,lon[i]+j);
			else if(j==j0)
			{
				lat[i][j] = baseline_lat[i];
				lon[i][j] = baseline_lon[i];
			}
			else
				latlon(baseline_lat[i],baseline_lon[i],
					fabs(delta),azimuth_i,
					lat[i]+j,lon[i]+j);
			r[i][j] = r0_ellipse(lat[i][j]);
		}
	}

	/* We now compute the Cartesian coordinates of all these points
	with an origin at the earth's center.  After we do compute all 
	the Cartesian coordinates we then do a translation of the 
	coordinates to the grid origin to make the numbers more 
	manageable*/
	for(i=0;i<n1;++i) 
		for(j=0;j<n2;++j)
		{
			dsphcar(lon[i][j],lat[i][j],x);
			x1[i][j] = x[0]*r[i][j];
			x2[i][j] = x[1]*r[i][j];
			x3[i][j] = x[2]*r[i][j];
		}
	/* This computes the radial direction from the latitude and longitude.
	We translate the origin in this direction below */
	dsphcar(lon0,lat0,x0);
	/* We construct the transformation matrix for rotation of coordinates
	from the origin using unit vectors computed as follows:
	column 3= local vertical = copy of x0
	column 2 = local y = constructed from pole to baseline
	column 1 = y cross z 
	This yields a rotation matrix stored in fortran order in rmatrix 
	Note use of ugly pointer arithmetic done to store the matrix this way */
	for(i=0;i<3;++i)rmatrix[i+6]=x0[i];
	dsphcar(pole_lon,pole_lat,rmatrix+3);
	dr3cros(rmatrix+3,rmatrix+6,rmatrix);
	/* We actually need the transpose of this matrix */
	dr3tran(rmatrix,rmtrans);

	/* We now apply the combined translation and rotation 
	change of coordinates */
	r0=r0_ellipse(lat0);
	for(i=0;i<3;++i) x0[i]*=r0;
	for(i=0;i<n1;++i) 
		for(j=0;j<n2;++j)
		{
			xwork[0] = x1[i][j];
			xwork[1] = x2[i][j];
			xwork[2] = x3[i][j];

			dr3sub(xwork,x0,xwork);
			dr3mxv(rmtrans,xwork,xwork2);

			x1[i][j] = xwork2[0];
			x2[i][j] = xwork2[1];
			x3[i][j] = xwork2[2];
		}

	/* We have to compute the extents parameters as the minimum 
	and maximum in each cartesian direction */
	xlow=x1[0][0];
	xhigh=x1[0][0];
	ylow=x2[0][0];
	yhigh=x2[0][0];
	zlow=x3[0][0];
	zhigh=x3[0][0];
	for(i=0;i<n1;++i)
		for(j=0;j<n2;++j)
		{
			xlow = MIN(x1[i][j],xlow);
			xhigh = MAX(x1[i][j],xhigh);
			ylow = MIN(x2[i][j],ylow);
			yhigh = MAX(x2[i][j],yhigh);
			zlow = MIN(x3[i][j],zlow);
			zhigh = MAX(x3[i][j],zhigh);
		}
	xlow = xlow;
	ylow = ylow;
	zlow = zlow;
	xhigh = xhigh;
	yhigh = yhigh;
	zhigh = zhigh;

	delete baseline_lat;
	delete baseline_lon;
}
//
// This is the C++ constructor that uses a pf and coord routines to compute
// a GCLgrid on the fly using the description coming from pf
//
GCLgrid3d::GCLgrid3d(int n1in, int n2in, int n3in,
	char *namein, 
	double lat0in, double lon0in, double r0in,
	double azimuth_yin, double dx1_nomin, double dx2_nomin,double dx3_nomin, 
	int i0in, int j0in)
{		

	/* pole to baseline */
	double pole_lat, pole_lon;
	int i,j,k;
	double deltax, deltay,delta;  
	double z0,z;
	double x[3],x0[3];
	double rmatrix[9],rmtrans[9];
	double xwork[3],xwork2[3];
	double dx1_rad,dx2_rad;
	double rotation_angle, azimuth_i;

	strcpy(name,namein);
	lat0=lat0in;
	lon0=lon0in;
	lat0 = rad(lat0);
	lon0 = rad(lon0);
	r0 = r0_ellipse(lat0);
	azimuth_y=azimuth_yin;
	rotation_angle=-azimuth_y;
	dx1_nom=dx1_nomin;
	dx2_nom=dx2_nomin;
	dx3_nom=dx3_nomin;
	dx1_rad = dx1_nom/r0;
	dx2_rad = dx2_nom/r0;
	n1=n1in;
	n2=n2in;
	n3=n3in;
	i0=i0in;
	j0=j0in;
//
// We force the origin of the r coordinate to be at the bottom of the grid
//
	k0=0;
	z0 = ((double)(n3-1))*dx3_nom;
	/* These are the vector of lat and lons used to define the
	baseline the effectively forms the equator of of the rotated
	coordinates passing through the origin*/
	double *baseline_lat=new double[n1in];
	double *baseline_lon=new double[n1in];

	build_baseline(lat0,lon0,
		rotation_angle,dx1_rad,n1,i0,baseline_lat,baseline_lon);
	/* We need to compute the pole to our base line to use as a 
	target for grid lines that are locally perpendicular along
	the grid lines -- like longitude lines at the equator */
	latlon(lat0,lon0,M_PI_2,-rotation_angle,&pole_lat,&pole_lon);

	/* We now allocate memory for the large coordinate arrays
	themselves and set the elements of the structure */
	cartesian_defined = 1;
	geographic_defined = 1;
	lat = create_3dgrid_contiguous(n1,n2,n3);
	lon = create_3dgrid_contiguous(n1,n2,n3);
	r = create_3dgrid_contiguous(n1,n2,n3);
	x1 = create_3dgrid_contiguous(n1,n2,n3);
	x2 = create_3dgrid_contiguous(n1,n2,n3);
	x3 = create_3dgrid_contiguous(n1,n2,n3);

	r0 = r0-z0;
	azimuth_y = -rotation_angle;

	/* We now complete the latitude/longitude grid by projecting
	lines from the baseline toward the computed pole.  Grid points 
	will be equally spaced along the gcp toward the pole, but the
	lines they form will converge.  This is complicated greatly by
	the grid origin parameters.  We could make this a function,
	but it is simpler to just write it inline. Because this function
	descended from the a C program call makegclgrid, which built both
	a 2d and 3d grid in one pass, I'm stealing that algorithm.  That is
	the top of the grid (at n3-1) is built first, then the remainder is
	filled in with same lat lons but with varying depth.  So, first this
	loop creates the top sheet of the grid */
	int top_of_grid=n3-1;
	for(i=0;i<n1;++i)
	{
		dist(baseline_lat[i],baseline_lon[i],pole_lat, pole_lon,
			&delta,&azimuth_i);
		if(azimuth_i < 0.0) azimuth_i += (2.0*M_PI);
		delta = dx2_rad*(-(double)j0);
		for(j=0;j<n2;++j,delta+=dx2_rad)
		{
			if(j<j0)
				latlon(baseline_lat[i],baseline_lon[i],
					fabs(delta),azimuth_i-M_PI,
					lat[i][j]+top_of_grid,lon[i][j]+top_of_grid);
			else if(j==j0)
			{
				lat[i][j][top_of_grid] = baseline_lat[i];
				lon[i][j][top_of_grid] = baseline_lon[i];
			}
			else
				latlon(baseline_lat[i],baseline_lon[i],
					fabs(delta),azimuth_i,
					lat[i][j]+top_of_grid,lon[i][j]+top_of_grid);
			r[i][j][top_of_grid] = r0_ellipse(lat[i][j][top_of_grid]);
		}
	}
//
//Now we fill in the rest of the grid keeping the lat lon of the top for all
//levels but varying r -- i.e. these are spherical shells.
//
	for(k=0,z=z0;k<n3-1;++k,z-=dx3_nom)
	{
		for(i=0;i<n1;++i)
			for(j=0;j<n2;++j)
			{
				lat[i][j][k]=lat[i][j][top_of_grid];
				lon[i][j][k]=lon[i][j][top_of_grid];
				r[i][j][k] = r0_ellipse(lat[i][j][top_of_grid]) - z;
			}
	}

	/* We now compute the Cartesian coordinates of all these points
	with an origin at the earth's center.  After we do compute all 
	the Cartesian coordinates we then do a translation of the 
	coordinates to the grid origin to make the numbers more 
	manageable.  First we computes the radial direction from the 
	latitude and longitude.
	We translate the origin in this direction below */
	dsphcar(lon0,lat0,x0);
	/* We construct the transformation matrix for rotation of coordinates
	from the origin using unit vectors computed as follows:
	column 3= local vertical = copy of x0
	column 2 = local y = constructed from pole to baseline
	column 1 = y cross z 
	This yields a rotation matrix stored in fortran order in rmatrix 
	Note use of ugly pointer arithmetic done to store the matrix this way */
	for(i=0;i<3;++i)rmatrix[i+6]=x0[i];
	dsphcar(pole_lon,pole_lat,rmatrix+3);
	dr3cros(rmatrix+3,rmatrix+6,rmatrix);
	/* We actually need the transpose of this matrix */
	dr3tran(rmatrix,rmtrans);

	/* We now apply the combined translation and rotation 
	change of coordinates */
	r0=r0_ellipse(lat0);
	for(i=0;i<3;++i) x0[i]*=r0;

	/* Now apply the same transformations to the 3d grid */

	for(i=0;i<n1;++i) 
	    for(j=0;j<n2;++j)
	    {
		for(k=0;k<n3;++k)
		{
			dsphcar(lon[i][j][k],lat[i][j][k],x);
			dr3sxv(r[i][j][k],x,xwork);

			dr3sub(xwork,x0,xwork2);
			dr3mxv(rmtrans,xwork2,xwork);

			x1[i][j][k] = xwork[0];
			x2[i][j][k] = xwork[1];
			x3[i][j][k] = xwork[2];
		}
	    }
	/* We have to compute the extents parameters as the minimum 
	and maximum in each cartesian direction */

	xlow=x1[0][0][0];
	xhigh=x1[0][0][0];
	ylow=x2[0][0][0];
	yhigh=x2[0][0][0];
	zlow=x3[0][0][0];
	zhigh=x3[0][0][0];

	for(i=1;i<n1;++i)
	    for(j=1;j<n2;++j)
		for(k=1;k<n3;++k)
		{
			xlow = MIN(x1[i][j][k],xlow);
			xhigh = MAX(x1[i][j][k],xhigh);
			ylow = MIN(x2[i][j][k],ylow);
			yhigh = MAX(x2[i][j][k],yhigh);
			zlow = MIN(x3[i][j][k],zlow);
			zhigh = MAX(x3[i][j][k],zhigh);
		}

	delete baseline_lat;
	delete baseline_lon;
}
//
//C++ destructors
//
GCLgrid::~GCLgrid()
{
	if(x1!=NULL) free_2dgrid_contiguous(x1,n1);
	if(x2!=NULL) free_2dgrid_contiguous(x2,n1);
	if(x3!=NULL) free_2dgrid_contiguous(x3,n1);
	if(lat!=NULL) free_2dgrid_contiguous(lat,n1);
	if(lon!=NULL) free_2dgrid_contiguous(lon,n1);
	if(r!=NULL) free_2dgrid_contiguous(r,n1);
}
GCLgrid3d::~GCLgrid3d()
{
	if(x1!=NULL) free_3dgrid_contiguous(x1,n1,n2);
	if(x2!=NULL) free_3dgrid_contiguous(x2,n1,n2);
	if(x3!=NULL) free_3dgrid_contiguous(x3,n1,n2);
	if(lat!=NULL) free_3dgrid_contiguous(lat,n1,n2);
	if(lon!=NULL) free_3dgrid_contiguous(lon,n1,n2);
	if(r!=NULL) free_3dgrid_contiguous(r,n1,n2);
}

