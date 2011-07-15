#ifndef _REGIONALCOORDINATES_H_
#define _REGIONALCOORDINATES_H_

#include "gclgrid.h"
/* \brief Define a properly georefrenced regional coordinate system.

   At local and regional scales in the earth it is normally more convenient
to define a set of local, Cartesian coordinates than constantly deal with 
spherical geometry.  This object provides a general solution to this problem
to a set of local coordinates that can implicitly be switched at any time
to geographic coordinates.  The coordinates and the underlying algorithms
are identical to that used in the GCLgrid library. 

Note the geometry of the coordinate system generated is defined only by
two quantities (4 parameters):  (1) point on the globe that is to be the 
origin, and (2) an orientation angle.  This works because the x1-x2 plane is 
always made to be the local horizontal at the origin.  An angle is used to 
define how the standard x1=+East and x2=+North axes are rotated to define the
local coordinate frame.  Many applications will not want to rotate the 
axes, but this is convenient, or example, in the GCLgrid library to define
a volume not oriented by the cardinal directions. 
*/
class RegionalCoordinates
{
public:
    /* Default constructor.  Place holder you should not use. */
    RegionalCoordinates();
    /* Primary constructor.

       The local coordinate system requires only a definition of the point that
    on the earth to use as the origin of the coordinate system and a 
    direction to define one of the axes. This constructor uses those to build
    the object.

    \param olat is the origin latitude in radians
    \param olon is the origin longitude in radians
    \param oradius is the radius (in km) for the origin.  (Use the r0_ellipse
    procedure in libgclgrid to make this the surface with ellipticity).
    \param azn is the angle to rotate x2 from north to define the coordinate 
      system.  azn=0 means x2 will point north.  azn=M_PI/4 points the x2 
      axis at north 45 degrees east.  Thus, azn uses the geographic azimuth
      convention and is in radians.
    */
    RegionalCoordinates(double olat, double olon, double oradius, double azn);
    /*! Standard copy constructor.*/
    RegionalCoordinates(const RegionalCoordinates& parent);
    /*! Standard assignment operator. */
    RegionalCoordinates& operator=(const RegionalCoordinates& parent);
    /* Return cartesian coordinates of a geographic point.
       \param lat is latitude in radians
       \param lon is longitude in radians
       \param r is radius in km.
       */
    Cartesian_point  cartesian(double lat,double lon, double r);
    /* Return cartesian coordinates of a geographic point.
       Point is specified by Geographic_point object. (radian and km units)*/
    Cartesian_point  cartesian(Geographic_point gp);
    /* Return cartesian coordinates of a geographic point.
       Point is passed as a 3 vector with the following content:
       x[0]=lat, x[1]=lon, x[2]=r (units radians and km). */
    Cartesian_point cartesian(double x[3]);
    Geographic_point geographic(double x1,double x2, double x3);
    Geographic_point geographic(double x[3]);
    Geographic_point geographic(Cartesian_point cp);
    Geographic_point origin();
    double aznorth_angle(){return(azimuth_y);};
private:
    double gtoc_rmatrix[3][3];
    double translation_vector[3];
    double lat0;
    double lon0;
    double r0;
    double azimuth_y;
};
#endif
