#ifndef _GEOPOLYGONREGION_H
#define _GEOPOLYGONREGION_H
#include <string>
#include <vector>
#include "gclgrid.h"
#include "GCLgridError.h"
/* \brief Define a polygon region mainly for masking a map area.

   The main purpose of this object is to define a polygonal region 
   on a map and be able to ask if a point is or is not inside this 
   polygon.   This is a basic tool for masking a region.  The primary
   method is thus a boolean function that answers that question yes or 
   no.  

   This object uses the winding number algorithm.  This is useful
   as it works for polygons that are not convex - a restriction of some
   common alternative algorithms.  

   This object uses a simple mapping with longitude=x and latitude=y 
   to define the geometry.  There is a tacit assumption the points are
   specified in degrees as a there are frozen constants that provide
   a incomplete solution to detect a limitation. User should note this 
   algorith will NOT work for polygons spanning a singular longitude 
   point (e.g. date line where longitude changes from 180 to -180) or
   across the poles.   Date line problems are readily fixed by switching 
   to 0 to 360 specification instead of -180 to 180.  Converse applies to 
   crossing the prime meridian.

   Although the input is required to be in degrees be aware that to be
   consistent with everything else in this library the points are
   stored internally in radians.   
   */
class GeoPolygonRegion
{
    public:
        /* \brief Default constructor.

           Explicit initialization to a null polygon. 
           */
        GeoPolygonRegion();
        /* \brief Construct from an ascii file of lon,lat pairs.
          \param fname is a file name
          \exception GCLgridError is thrown if the file open fails
         */
        GeoPolygonRegion(string fname);
        /* \brief Construct from parallel vectors of lat and lon values.

           \param vlat vector of latitude values in degrees.
           \param vlon parallel vector to vlat of longitude values (degrees)

           \exception throws GCLgridError object in one of three conditions:
            (1) vlat and vlon are not the same length, (2) vlat or vlon
            values are outside normal range, and (3) any pair of point have
            a change in either lat or lon angle larger than 90 degrees. The 
            last two are a crude solution to try to detect polygons 
            spanning the date line or across the poles. 
            */
        GeoPolygonRegion(vector<double> vlan, vector<double> vlon);
        /* \brief Construct from vector of geographic point objects.

           \param pts vector of geographic points.  The important
             difference here is that because Geographic_point objects
             are tacitly assumed to be in radians these are ALWAYS 
             converted from radians to degrees internally.

           \exception throws GCLgridError object if latitude and longitude
             values converted to degrees are illegal.  This is an inadequate
             attempt to automatically detect incorrect specification in 
             degrees (Geographic_point data is assumed radians).
            */
        GeoPolygonRegion(vector<Geographic_point> pts);
        /*! Standard copy constructor */
        GeoPolygonRegion(const GeoPolygonRegion& parent);
        /* \brief Method to find if point is inside polygon.

           This method returns true if the point (lat0,lon0) is inside
           the polygon defined by the object.  

           \param lat0 - latitude (degrees) of point to test
           \param lon0 - longitude (degrees) of point to test
           */
        bool is_inside(double lat0, double lon0);
        /*! Standard assignment operator. */
        GeoPolygonRegion& operator=(const GeoPolygonRegion& parent);
        /*! Standard output stream operator.  List points in degrees.*/
        friend ostream& operator << (ostream& os, GeoPolygonRegion& p);
    private:
        int npoints;  // Number of points in polygon
        /* Points stored interleaved x[0],y[0],x[1],y[1],...,y[npoints-1]*/
        vector<double> lon,lat;
};
#endif
