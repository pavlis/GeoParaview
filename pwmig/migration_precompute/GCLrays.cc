#include "gclgrid.h"
#include "vmodel.h"


/* This function takes a 2D GCLgrid and a given ray parameter
 * and builds a GCLgrid3D object defined by rays derived from a
 * give velocity model and the input ray parameter.  This 
 * builds rays at different back azimuths because for a large
 * grid the local slowness vector for a fixed azimuth will point
 * in different directions.  
 *
 *
 * Ok, this needs to be done differently.  In fact it impacts the
 * plane wave code too.  We are going to need to estimate an epicentral
 * distance give a ray parameter at the origin and then provide
 * an algorithm to compute the ray parameter at any point in a GCLgrid
 * derivable from that p.  That is, what p is correct for source at 
 * x and receiver at each of the GCLgrid points.  Call it a virtual
 * source point.  
 * */

new code goes here.  
