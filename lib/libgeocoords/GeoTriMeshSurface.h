#ifndef _GEOTRIMESHSURFACE_H_
#define _GEOTRIMESHSURFACE_H_
#include <vector>
#include "cgeom.h"
#include "gclgrid.h"
#include "GeoSurface.h"
#include "GeoPolygonRegion.h"

using namespace std;
/* This enum is needed to define coordinate units. It may belong in a more 
   generic include file */
enum GeographicUnits {RADIANS,DEGREES};

/*! Define a surface in geographic coordinates from random points. 

This builds a surface in geographic coordinates (i.e. latitude,
longitude, and radius values) from a random set of point.  It
is implemented using a Delauney trangularization in the 2d space
defined by the lat and lon coordinates of the points with 
radius values defining a function defined by this group of points.
The radius and depth methods automate interpolation of points
inside the convex hull of the set of points used to construct
the object. 

Users of this code should be sure that all data are consistent in 
how they handle lat and lon coordinates.  That is, the object can be 
constructed with lat and lon in degrees or radians, but chaos will 
obviously result if you mix units.  A more subtle restriction is the 
convention used for wrapping the angles.  Latitude is rarely an issue,
but longitude is commonly specified either 0 to 360 degrees or 
-180 to 180 degrees.  Do NOT mix these conventions.   The 
triangularization will still work fine, but the geometry you get
will by crazy.  Finally, you have define what the coordinate units
are when the object is constructed, BUT constructors have a default
that assumes radians.   Finally, because this is based on spherical
coordinates there is the usual warning that it won't work across the 
pole.  Be careful as well if the points cross the longitude wrap 
point.  (i.e. the date line (-180 to 180 convention) or the prime
meridian (0 to 360 convention)).

This is implemented with Kent Lindquist's Computational Geometry library 
in antelope contrib.  
\author Gary L. Pavlis
*/

class GeoTriMeshSurface : public GeoSurface
{
    public:
        /*! Construct from container of Geographic_point data structures. 
         
         \param pts is a vector of Geographic_point objects used to define
            surface.  
         \param units sets units of lat and lon values (default "radians").  
           If set to ANYTHING else will assume degrees. 
         \exception GeoCoordError is thrown if number of points is less than 3.
           This is done since we use triangularization, which is clearly meaningless
           if we have fewer than 3 points 
         */
        GeoTriMeshSurface(vector<Geographic_point> pts,string units="radians");
        /*! Construct from parallel sets of vectors. 

          This is a wrapper routine to allow an alternative specification in
          terms of lat,lon and depth (note GeogrpahicPoint objects are radii) 
          vectors. 
         \param lat is a vector latitude values.
         \param lon is a vector of longitude values
         \param d is a vector of depth values.
         \param units sets units of lat and lon values (default "radians"). 
           If set to ANYTHING else will assume degrees. 
         \exception GeoCoordError object is thrown if three vectors are not the same 
            length.
         */
        GeoTriMeshSurface(vector<double> lat, vector<double> lon, 
                vector<double>d,string units="radians");
        /*! Standard copy constructor. */
        GeoTriMeshSurface(const GeoTriMeshSurface& parent);
        /*! Destructor. 

          In this implementation the grid triangularization is stored 
          internally.  The destructor is nontrivial as it has to 
          release this work space. */
        ~GeoTriMeshSurface();
        /*! Standard assignment operator. */
        GeoTriMeshSurface& operator=(const GeoTriMeshSurface& parent);
        /*! Add a polygonal boundary region.

          This object always is constrained by the convex hull of the set of points.
          Sometimes, however, it is useful to carve out part of the surface 
          particularly when part of the object is not convex.   This is added
          through the GeoPolygonRegion object that is copied and used internally.
          When this is called the boundary is the interior of both the convex
          hull and the polygon defined here. 

          \param poly defines the polygon region. 
          */
        void AddBoundary(const GeoPolygonRegion& poly);
        /*! Return the radius of the surface at lat,lon.  
         
         This method interpolates the surface and returns a radius value
         at a requested lat,lon.  Assumes units of lat, lon are consistent
         with flag used to construct the object.  Returns a nan that
         if point is outside grid.  This should be tested with is_nan 
         procedure in antelope library.*/
        double radius(double lat, double lon);
        /*! Return the depth of the surface at lat,lon.  
         
         This method interpolates the surface and returns a depth value
         at a requested lat,lon.  Assumes units of lat, lon are consistent
         with flag used to construct the object.  Returns a nan that
         if point is outside grid.  This should be tested with is_nan 
         procedure in antelope library.*/
        double depth(double lat, double lon);
        /*! Test if a point is with a region defined by this surface.

          Unless a surface is defined for the entire earth it has bounds.   A
          triangular mesh is implemented here by Delaunay triangularization so 
          this method tests to see if a point is inside the convex hull for 
          the triangulation of that surface. */
        bool is_defined(double lat, double lon);
    private:
        CGPartition *trigrid;
        GeographicUnits units;
        bool polygon_boundary_defined;
        GeoPolygonRegion boundary;
};
#endif
