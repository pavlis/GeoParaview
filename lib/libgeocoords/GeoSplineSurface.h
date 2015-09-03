#ifndef _GEOSPLINESURFACE_H_
#define _GEOSPLINESURFACE_H_

#include "gclgrid.h"
#include "Metadata.h"
#include "cgeom.h"
#include "GeoSurface.h"
#include "GeoPolygonRegion.h"
using namespace std;
using namespace SEISPP;

/*! \brief Geographic coordinate surface defined by bicubic splines.

  Bicubic splines are often used to define surfaces in 2d.  Bicubic splines
  fit a curve through a random collection of points.  This implementation
  can be thought of as C++ interface to the Generic Mapping Tools program
  called surface.  The user is instructed to read the man page for that 
  program, look at the cookbook examples, and refer to the paper
  by Smith and Wessel (1990, Geophysics, 55, pp. 293-305) in order
  to understand the proper application of the algorithm.  
  This applies the method called splines under tension in which the 
  tension parameters can be used to prevent wild oscillations between
  irregularly sampled data.  

  Be aware this implementation is actually based on Kent Lindquist's cgeom 
  library in Antelope.  It would be quite feasible to change the implementation
  to the GMT code, but I chose to use cgeom because I would have had to shred
  and reassemble the C code in gmt while Lindquist's had already built a clean
  C api.  This object is little more than a front end to his implementation.

  A detail of the implementation the user must recognize is that the way points
  are handled that are defined as outside the region of interest.  There are 
  two different approachs.  First, the constructor builds a regular grid 
  between a specified range of specified latitudes and longitudes (a spherical
  rectangle in standard coordinates).  Any point within that "rectangular" area
  should return a valid value when the depth or radius method are called.  
  Because points outside the region of control are often very unreliable with
  spline interpolation, however, this object provides a fast test in the 
  is_defined method using a totally different algorithm.  That is, the input
  points are passed through a Delaunay triangulation routine.  The is_defined
  method can then be viewed as testing to see if the given point is within
  the convex hull defined by the triangulation.

  Finally, note some mixed units documented in the detailed description of
  parameters in this object.  That is, some lat,lon parameters are assumed in
  degrees while other inputs are assumed to be radians.
  */
class GeoSplineSurface : public GeoSurface
{
public:
    /*! Full parameterized constructor for splined surface.

    This constructor allows passing the complete set of parameters to
    build the object as arguments.  This list is long and relate to the
    complexity of this algorithm.  A key limitation of this method is it
    lays down a regular grid of points on a lon,lat grid in the same
    manner as the GMT surface program.  

    \param pts - is the irregular grid of points used to defined the surface.
      The radius component of each point is converted internally to depth 
      using a standard ellipsoid. The latitude and longitude components
      are assumed to be in radians.
    \param minx is the minimum longitude (degrees) of the working rectanglular
      area.
    \param maxx is the maximum longitude (degrees) of the working rectanglular
          area.
    \param miny is the minimum latitude (degrees) of the working rectanglular
      area.
    \param maxy is the maximum latitude (degrees) of the working rectanglular
          area.
    \param dx rectangular grid spacing in longitude
    \param dy rectangular grid spacing in latitude coordinate.
    \param lower allows a floor on an interpolated depth value.  (default
      is -8 km which means this is turned off since that is the highest
      topography on earth.)
    \param upper is a ceiling on returned interpolated depth values.  
      (default is 6371.0 = earth radius which means this is effectively off)
    \param tension_interior bicubic spline tension parameter for
      interior points.  See GMT man page and Smith and Wessel paper
      for guidelines.  (default 0.25)
    \param tension_boundary bicubic spline tension parameter for
      interior points.  See GMT man page and Smith and Wessel paper
      for guidelines.  (default 1.0 which forces the curve defining
      the convex hull to be linear between control points.)
    \param aspect_ratio is used if the regular grid geometry given is
      very different in x and y.  (default 1.0)
    \param overrelaxation is a convergence parameter.  Default of 1.4 is
      known to nearly always work.
    \param convergence specifies the fraction of the range of the input 
      data that defines convergence (default 0.001)
    \param max_iterations is the maximum number of iterations allowed
      before giving up. (default 250)
    \param use_ch sets if the object should apply the convex hull
      of the collection of points.   When true points outside
      the convex hull will be effectively undefined.  (default is false
      because implementation in antelope has known problem with an 
      infinite loop when any 3 points are colinear in the set. )
    */
    GeoSplineSurface(vector<Geographic_point>& pts,
            double minx,
             double maxx,
              double miny,
               double maxy,
                double dx,
                 double dy,
                  double lower=-8.0,
                   double upper=6371.0,
            double tension_interior=0.25,
             double tension_boundary=1.0,
              double aspect_ratio=1.0,
               double overrelaxation=1.4,
                double convergence=0.001,
                 int max_iterations=250,
                  bool use_ch=false);
    /* \brief Construct with control parameters coming from a Metadata object.

       In my library I use the Metadata object heavily to encapsulate a set of 
       adjustable parameters that apply to a particular program.  This interface
       allows one to put the long string of parameters required by the spline
       in tension algorithm into a generic form that can be read to construct
       a Metadata object (currently this means Antelope parameter files, but
       the concept is more general).  The list of keys required are similar, but
       not identical to the arguments to the fully parameterized constructor 
       for this object.  The list of required parameters in argument list order
       are:  grid_longitude_minimum, grid_longitude_maximum, 
       grid_latitude_minimum,grid_latitude_maximum, lower_bound_depth,upper_bound_depth,
       tension_interior, tension_boundary, aspect_ratio, overrelaxation, convergence,
       max_iterations.  The max_iterations parameter is integer and the other are
       all real numbers. 

    \param pts - is the irregular grid of points used to defined the surface.
      The radius component of each point is converted internally to depth 
      using a standard ellipsoid. The latitude and longitude components
      are assumed to be in radians.
    \param inpar is the Metadata object containing the required parameters as listed above.

    \exception GeoCoordError is returned if there are problems.  The most likely is 
      missing Metadata components.  
       */
    GeoSplineSurface(vector<Geographic_point>& pts, Metadata& inpar);
    /*! Standard Copy Constructor */
    GeoSplineSurface(const GeoSplineSurface& parent);
    /*! Standard destructor.

      Nontrivial destructor has to free grid data.
      */
    ~GeoSplineSurface();
    /*! Add a polygonal boundary region.

      Sometimes we need to build a bounding curve that has an 
      arbitrary shape.   This adds that capability.  Default uses
      the convex hull of input points. 

      \param poly defines the polygon region.
      */
    void AddBoundary(const GeoPolygonRegion& poly);
    /*! Return radius of surface.

      This method returns the spherical coordinate radius (in km) of the surface
      being interpolated.  
      \param lat is the latitude (in radian) of the point of interest.
      \param lon is the longitude (in radians) of the point of interest.

      \returns radius in km at the point of interest. 
      \exception GeoCoordError is thrown if the requested point is not within the 
        area of support.  For this reason caller should use the is_defined method 
        before using this method to avoid the overhead of catching many exceptions.
     */
    double radius(double lat, double lon);
    /*! Return depth of surface.

      This method returns the depth to the surface being interpolated (in km)
      at a point of interest.  Depth is relative to the reference ellipsoid.
      \param lat is the latitude (in radian) of the point of interest.
      \param lon is the longitude (in radians) of the point of interest.

      \returns depth in km at the point of interest. 
      \exception GeoCoordError is thrown if the requested point is not within the 
        area of support.  For this reason caller should use the is_defined method 
        before using this method to avoid the overhead of catching many exceptions.
     */
    double depth(double lat, double lon);
    /*! Used to ask if a point is inside a defined surface.

      The interpolation method used in this object does not work outside the area
      of support of data used to defined the surface.  This method should be called
      before any interpolation attempt to effectively ask if the point is known.  
      In this implementation this is necessary to avoid a try/catch block handler,
      which is a far less elegant and far slower solution.
      \param lat is the latitude (in radian) of the point of interest.
      \param lon is the longitude (in radians) of the point of interest.
      */
    bool is_defined(double lat,double lon);
    /*! Standard assignment operator */
    GeoSplineSurface& operator=(const GeoSplineSurface& parent);
    /*! call to enable convex hall editing if not turned on originally.*/
    void enable_convex_hull();
private:
    bool use_convex_hull;
    CGPartition *trigrid;
    CGPointset *pointset;
    CGGrid *gridptr;
    /* extents of regular grid stored in gridptr 
     stored in radians*/
    double latmin,latmax,lonmin,lonmax;
    /* This private method contains common code for constructors.  
       Done for maintainability although it confuses things.  Note
       the fully parameterized constructor does nothing but call
       this routine. */
    void GSSinit(vector<Geographic_point>& pts,
            double minx, double maxx, double miny, double maxy,
            double dx, double dy, double lower, double upper,
            double tension_interior, double tension_boundary, 
            double aspect_ratio, double overrelexation, 
            double convergence, int max_iterations);
    GeoPolygonRegion boundary;
};
#endif
