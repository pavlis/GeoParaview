#ifndef _PLGEOPATH_H_
#define _PLGEOPATH_H_
#include "GeoCoordError.h"
#include "GeoPath.h"
#include "RegionalCoordinates.h"
/*! Geographic 3d path specified as piecewise linear path defined by finite set of control points.
This represents the simplest form of a path defind on the earth.  It is constructed from a finite
set of control points linked by great circle paths.  The distance parameter is defined by 
distance along the path (in km) from the origin.  An approximation is used for the depth component. It
is presumed the total distance between two points is sqrt(distance(lat1,lon1,lat2,lon2)^2+dz^2)
That is, it is assymptotic to L2 when the two control points are close enough to neglect 
great circle path curvature.  
*/
class PLGeoPath : public GeoPath
{
public:
    /* This must be public to allow the Cartesian_point methods to make
       any sense. */
     RegionalCoordinates coordxyz;
    /*! Default constructor.  Minimal initialization only. */
    PLGeoPath();
    /*! Construct from a vector of Geographic_point objects.

      This constructor uses and order list of points that define the tie points for
      the curve being defined.  The curve created is linear (more precisely great circle
      arc segments with linear depth interpolation) segments defined between the 
      points given to this constructor. 
      \param pts is the ordered list of points used to define this path.
      \param i0 defines which point in the input vector is to be treated as the s=0
        point for the path.  Note this is C convention so the first point is i0=0.
      \param az0 can be used to change the azimuth of the internal coordinate system.  
        Default is 0.0 and there is no reason one should ever need to change this.
      \exception Throws a GeoCoordError exception for a number of potential error conditions.
        */
     PLGeoPath(vector<Geographic_point>& pts,int i0,
             double az0=0.0);
     /*! Standard copy constructor. */
     PLGeoPath(const PLGeoPath& parent);
     /* \brief Return spatial position of point at a given curve parameter.

        This object can be used to define 3D curves within the earth parameterized by 
        distance, sp (in km), from a specified origin.  This method returns the geographical
        coordinates of the curve for a specified distance sp. Note this method always returns
        an answer. If the distance sp is beyond the range of support of control points the curve
        is extrapolated using a linear project of the two points closest to the appropriate
        endpoint.
        \param sp is the distance in km from the origin whose coordinates are requested.
        \return point requested as a Geographic_point object.
        */
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
     /*! Return the position of node i.

       Sometimes we need the actual node locations of a path.  This
       returns the position of node i.

       \param i is the index position of the requested node.
       \exception SeisppError will be thrown if i is out of range.
       */
     Geographic_point node_position(int i);
     /*! Return the Cartesian position of node i.

       Sometimes we need the actual node locations of a path.  This
       returns the position of node i in the internal Cartesian system.

       \param i is the index position of the requested node.
       \exception SeisppError will be thrown if i is out of range.
       */
     Cartesian_point node_position_xyz(int i);
     /*! Return path parameter for node i.

       Sometimes we need to know the path parameter matching a particular
       node point.   

       \param i is the index position of the requested node.
       \exception SeisppError is throw if i is out of range */
     double node_path_parameter(int i);
     PLGeoPath& operator=(const PLGeoPath& parent);
     friend ostream& operator<<(ostream& os,PLGeoPath&);
private:
     /*! Store control points in four parallel vector containers */
     vector<double> lat,lon,r,s;
     /*! store this point to allow it to be anywhere on the path and not have
       to interpolate to figure it's location */
     Geographic_point s0;  
};
#endif
