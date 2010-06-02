#ifndef _PLGEOPATH_H_
#define _PLGEOPATH_H_
#include "GeoPath.h"
#include "RegionalCoordinates.h"
/*! Geographic 3d path specified as piecewise linear path defined by finite set of control points.
This represents the simplest form of a path defind on the earth.  It is constructed from a finite
set of control points linked by great circle paths.  The distance parameter is defined by 
distance along the path (in km) from the origin.  An approximation is used for depth where it
is presumed the total distance between two points is sqrt(distance(lat1,lon1,lat2,lon2)^2+dz^2)
That is, it is assymptotic to L2 when the two control points are close enough to neglect 
great circle path curvature.  
*/
class PLGeoPath : public GeoPath
{
public:
    /*! Construct from a vector of Geographic_point objects.

      This assumes points are specified as radians.  pts[i0] sets
      the origin of the path distance coordinates. */
     PLGeoPath(vector<Geographic_point>& pts,int i0,
             double az0=0.0);
     PLGeoPath(const PLGeoPath& parent);
     Geographic_point position(double sp);
     Cartesian_point position_xyz(double sp);
     double latitude(double sp);
     double longitude(double sp);
     Geographic_point origin();
     /* These return first and last point in the path */
     Geographic_point begin();
     Geographic_point end();
     /* These return the parameter value matching the endpoints */
     double sbegin();
     double send();
     int number_points(){return(lat.size());};
     PLGeoPath& operator=(const PLGeoPath& parent);
     friend ostream& operator<<(ostream& os,PLGeoPath&);
private:
     /*! Store control points in four parallel vector containers */
     vector<double> lat,lon,r,s;
     /*! store this point to allow it to be anywhere on the path and not have
       to interpolate to figure it's location */
     Geographic_point s0;  
     RegionalCoordinates coordxyz;
};
#endif
