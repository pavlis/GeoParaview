#define RADIUS_EARTH 6378.164
/* For plane wave moveout computations a local cartesian coordinate
system is used wrt to a particular origin.  This approximation is
reasonable until the distance of a station from the origin 
becomes a significant fraction of the epicentral distance.  
In pseudostation stacking we can minimize the impact of this approximation
by continually translating the origin to the pseudostation point and 
computing distances to all stations wrt to that point.  Since large
distance stations receive a low weight, time alignments at large
distances become unimportant.  this function is used to computer 
vector distances in this context.

Arguments:
	nsta - number of stations to process = length of lat and lon vectors
			(see below).
	lat0, lon0 - origin to compute dnorth deast from
	lat, lon - vectors of length nsta of station coordinates to be 
		converted.
	dnorth, deast - vectors of length nsta to contain the results 
		(These are geographic dirctions +north and + east
		respectively)

Author:  G Pavlis
Written:  June 2000
*/
void geographic_to_dne(int nsta,double lat0, double lon0,
	double *lat, double *lon, double *dnorth, double *deast)
{
	int i;
	double azimuth;
	double d; 

	for(i=0;i<nsta;++i)
	{
		dist(rad(lat0),rad(lon0),rad(lat[i]),rad(lon[i]),&d,&azimuth);
		d *= RADIUS_EARTH;
		deast = d*sin(azimuth);
		dnorth = d*cos(azimuth);
	}
}
/* this function allocates and returns a vector of weights using
a 2d Gaussian function centered on lat0 and lon0 to implement
pseudostation stacking as described in Neal and Pavlis (1999) grl.

Arguments:
	nsta - number of stations 
	lat0, lon0 - image point in geographical coordinates
	lat, lon - vectors of length nsta containing geographic coordinates
		of the nsta stations
	pf - parameter file handle. (used to define aperture parameters)

The function returns a vector of nsta doubles that are the weights
to use for pseudostation stacking at lat0, lon0.  

This function is slightly inefficient in it's computing distances exactly
using the dist function, but the cost of this is not high enough to 
warrant the cartesian approximation.  In addition, the use of a Pf object
to pass aperture information is inefficient, but was done because it
will make it easy to make a plug compatible replacement module for
experimentation with alternative weighting functions to a Gauusian.

Author:  G Pavlis
Written:  June 2000
*/
double *compute_pseudostation_weights(int nsta, double lat0, double lon0,
	double *lat, double *lon, Pf *pf)
{
	int i;
	double *w;
	double aperture;
	double cutoff;
	double distance, azimuth, denom;

	aperture = pfget_double(pf,"pseudostation_smoother_width");
	cutoff = pfget_double(pf,"pseudostation_distance_cutoff");
	denom = 2.0*aperture*aperture;

	allot(double *,w,nsta);
	for(i=0;i<nsta;++i)
	{
		dist(lat0,lon0,lat[i],lon[i],&distance,&azimuth
		distance *= RADIUS_EARTH;
		if(distance>cutoff)
			w[i] = 0.0;
		else
			w[i] = exp(distance*distance/denom);
	}
	return(w);
}
double *compute_pwmoveout(int nsta,double *deast, double *dnorth,
		double ux, double uy)
{
	double *moveout;
	int i;

	allot(double *,moveout,nsta);
	for(i=0;i<nsta;++i)
		moveout[i] = ux*deast + uy*dnorth;
	return(moveout);
}
/* This function fills the PWslowness_grid structure defined in
pwmig.h from definitions in the parameter file passed through the
handl pf and returns the filled out structure 

Author:  G Pavlis
Written:  August 2000
*/ 
PWslowness_grid setup_slowness_grid(Pf *pf)
{
	PWslowness_grid g;
	g.nx = pfget_int(pf,"ux_number_gridpoints");
	g.ny = pfget_int(pf,"uy_number_gridpoints");
	g.uxlow = pfget_double(pf,"ux_lower_limit");
	g.uxhigh = pfget_double(pf,"ux_upper_limit");
	g.uylow = pfget_double(pf,"uy_lower_limit");
	g.uxhigh = pfget_double(pf,"uy_upper_limit");
	if( (g.uxhigh<g.uxlow) || (g.uyhigh<g.uylow))
		die(0,"Illegal slowness grid specification ux=(%lf,%lf) and uy=(%lf,%lf)\nCheck parameter file\n",
			g.uxlow,g.uxhigh,g.uylow,g.uyhigh);
	g.dux = (g.uxhigh-g.uxlow)/((double)(g.nx - 1));
	g.duy = (g.uyhigh-g.uylow)/((double)(g.ny - 1));
	return(g);
}
