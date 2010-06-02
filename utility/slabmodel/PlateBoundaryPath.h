#ifndef _PLATEBOUNDARYPATH_H_
#define _PLATEBOUNDARYPATH_H_
#include <float.h>
#include "gclgrid.h"
#include "GeoPath.h"
/*! Object to project a small circle path defining a plate boundary.

  In plate tectonics transform plate boundaries and motion vectors
are defined by small circles (latitude like lines) relative to a pole
of rotation.  This object simplifies the task of creating these
lines.  Curves are parameterized by time.  Normally that would
be in years. */

class PlateBoundaryPath : public GeoPath
{
public:
    PlateBoundaryPath(double pla, double plo, double ola, double olo,double avel);
    PlateBoundaryPath(const PlateBoundaryPath& parent);
    PlateBoundaryPath& operator=(const PlateBoundaryPath& parent);
    /* Return coordinates as a GCLgrid Geographic_point object at distance 
       s (in km) from origin along the small circle */
    Geographic_point position(double s);
    double latitude(double s);
    double longitude(double s);
    Geographic_point origin();
private:
    double plat,plon,olat,olon;
    double ssdelta;  // small circle distance from pole to origin (radians)
    double azimuth0;  // azimuth at pole to point define as origin for path
    double angular_velocity;  //angular velocity in radians
};
#endif
