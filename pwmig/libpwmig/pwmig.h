/* This file contains typdefs for structures used interally in 
programs used for the plane wave migration method introduced by
Poppeliers and Pavlis (2000).  */


/* This is a gather object for time aligned seismic data with a fixed
sample rate.  It uses parallel arrays of simple data types to define
the gather.  It is appropriate for either raw, time aligned data or
deconvolved seismic data.  Note time alignment allows us to drop
the time field required in raw data gathers */
typedef struct tracegather_ {
	int nsta;
	int nsamp;
	char gridname[10];
	int ix1, ix2;  /* pseudostation indices */
	double ux,uy;  /* Base moveout slowness vector in km/s */
	double si;  /* sample interval in seconds */
	char **sta;   /* character string vector of station names */
	double *lat, *lon;  /* parallel vectors of latitude and longitude
		for each sta.  e.g. lat[i], lon[i] are location of sta[i]*/
	float **x, **y, **z;   /* trace data in cartesian coordinates 
				defined wrt to transformation matrix */
	double U[9];  /* Transformation matrix of x,y,z coordinates relative
			to geographical standard coordinates (fortran order
			3x3 matrix)  I when in standard coordinates */
} PWMIGraw_trace_gather;
/* A slowness grid in this the pwmig package is rectangular in 
slowness space with uniformly spaced grid points.  There
are nx grid points between uxlow and uxhigh and ny grid points
between uylow and uyhigh.  */
typedef PWslowness_grid_ {
	int nx,ny;
	double uxlow, uxhigh, uylow, uyhigh;
	double dux, duy;
} PWslowness_grid;
