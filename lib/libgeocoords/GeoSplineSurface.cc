#include <fstream>
#include "SeisppError.h"
#include "GeoCoordError.h"
#include "GeoSplineSurface.h"
using namespace std;
using namespace SEISPP;

GeoSplineSurface::GeoSplineSurface(vector<Geographic_point>& pts,
  double minx,
   double maxx,
    double miny,
     double maxy,
      double dx,
       double dy,
        double lower,
         double upper,
          double tension_interior,
           double tension_boundary,
            double aspect_ratio,
             double overrelaxation,
              double convergence,
               int max_iterations,
                bool use_cg) : boundary()
{
    this->GSSinit(pts,minx,maxx,miny,maxy,dx,dy,lower,upper,
            tension_interior,tension_boundary,aspect_ratio,
            overrelaxation,convergence,max_iterations);
    if(use_cg) 
    {
        trigrid=delaunay_On4(pointset);
        use_convex_hull=true;
    };
}
void GeoSplineSurface::GSSinit(vector<Geographic_point>& pts,
  double minx,
   double maxx,
    double miny,
     double maxy,
      double dx,
       double dy,
        double lower,
         double upper,
          double tension_interior,
           double tension_boundary,
            double aspect_ratio,
             double overrelaxation,
              double convergence,
               int max_iterations) 
{
    latmin=rad(miny);
    latmax=rad(maxy);
    lonmin=rad(minx);
    lonmax=rad(maxx);
    vector<Geographic_point>::iterator ppts;
    CGPointset *pointset=cgpointset_new();
    for(ppts=pts.begin();ppts!=pts.end();++ppts)
    {
        /* Note we store the depth so we have to convert r.  
           Here I use the gclgrid library r0_radius procedure.
           That currently does not have any flexibility for
           the reference ellipsoid. */
        double dep=r0_ellipse(ppts->lat) - ppts->r;
        CGPoint *cgpt=cgpoint_new(ppts->lon,ppts->lat,dep);
        cgpointset_addpoint(pointset,cgpt);
        cgpoint_free(&cgpt);
    }
    CGGrid *reggrd=cggrid_new(rad(minx),rad(maxx),rad(miny),rad(maxy),
            rad(dx),rad(dy));
    /* This may not be necessary, but worth setting explitly. */
    strcpy(reggrd->units,"RADIANS");
    gridptr=cggrid_splinefit_pointset(reggrd,
            pointset,tension_interior,tension_boundary,aspect_ratio,
             overrelaxation,convergence,max_iterations,lower,upper);
}
GeoSplineSurface::GeoSplineSurface(vector<Geographic_point>& pts,
            Metadata& inpar) : boundary()
{
    try {
        double minx=inpar.get_double("grid_longitude_minimum");
        double maxx=inpar.get_double("grid_longitude_maximum");
        double miny=inpar.get_double("grid_latitude_minimum");
        double maxy=inpar.get_double("grid_latitude_maximum");
        double dx=inpar.get_double("grid_longitude_spacing");
        double dy=inpar.get_double("grid_latitude_spacing");
        double lower=inpar.get_double("lower_bound_depth");
        double upper=inpar.get_double("upper_bound_depth");
        double tension_interior=inpar.get_double("tension_interior");
        double tension_boundary=inpar.get_double("tension_boundary");
        double aspect_ratio=inpar.get_double("aspect_ratio");
        double overrelaxation=inpar.get_double("overrelaxation");
        double convergence=inpar.get_double("convergence");
        int max_iterations=inpar.get_int("max_iterations");
        this->GSSinit(pts,minx,maxx,miny,maxy,dx,dy,lower,upper,
            tension_interior,tension_boundary,aspect_ratio,
            overrelaxation,convergence,max_iterations);
        if(inpar.get_bool("use_convex_hull"))
        {
            trigrid=delaunay_On4(pointset);
            use_convex_hull=true;
        }
    } catch (SeisppError& serr) 
    {
        /* Translate metadata get errors into GeoCoordsError
           for consistency with other components of this library */
        throw GeoCoordError(serr.what());
    }
}
GeoSplineSurface::GeoSplineSurface(const GeoSplineSurface& parent)
{
    throw GeoCoordError(string("GeoSplineSurface copy constructor not implemented yet"));
}
GeoSplineSurface& GeoSplineSurface::operator=(const GeoSplineSurface& parent)
{
    throw GeoCoordError(string("GeoSplineSurface operator = not implemented yet"));
}

GeoSplineSurface::~GeoSplineSurface()
{
    cggrid_free(&gridptr);
    cgpointset_free(&pointset);
    cgpartition_free(&trigrid);
}
string depth_error_message(double lat,double lon,string line2)
{
    char ebuf[256];
    stringstream ss(ebuf);
    ss << "GeoSplineSurface::depth method:  "
            << "latitude="<<deg(lat)
            << " longitude="<<deg(lon)<<endl
            << line2<<endl;
    return(ss.str());
}
/* This is the core method here.  radius is only an interface 
   that uses this in this implementation.  An important detail
   of the implementation is that it returns a result for any
   point in the extents defined by the gggrid_spline_pointset 
   (a uniform grid in lat and lon).  Use the is_defined to 
   select points only in the convex hull.  

   Throw an exception if the point is outside the defined region.*/
double GeoSplineSurface::depth(double lat, double lon)
{
    if(lat<latmin || lat>latmax || lon<lonmin || lon>lonmax)
        throw GeoCoordError(depth_error_message(lat,lon,
            string("Point is outside domain of grid defining this surface.")));
    double result;
    result=cggrid_probe(gridptr,lon,lat);
    if(is_nan(result)) 
        throw GeoCoordError(depth_error_message(lat,lon,
            string("lookup method (cggrid_probe) failed for this point.")));
    return(result);
}
double GeoSplineSurface::radius(double lat, double lon)
{
    double r0=r0_ellipse(lat);
    try {
        return(r0-this->depth(lat,lon));
    }catch(...){throw;};
}
void GeoSplineSurface::AddBoundary(const GeoPolygonRegion& poly)
{
    boundary=poly;
}
/* Note the test here is fundamentally different from that used by 
   the depth method.  The depth method compares lat and lon to a rectangular
   region (min,max).  This uses the convex hull defined by the Delaunay
   triangularization stored in trigrid. */
bool GeoSplineSurface::is_defined(double lat,double lon)
{
    if(use_convex_hull)
    {
        CGPolygon *poly=polygon_containing_xy(trigrid,lon,lat);
        if(poly==NULL)
            return false;
    }
    if(boundary.is_inside(lat,lon))
        return true;
    else
        return false;
}
void GeoSplineSurface::enable_convex_hull()
{
    // Do nothing if trigrid is already created
    if(trigrid!=NULL) return;
    trigrid=delaunay_On4(pointset);
    use_convex_hull=true;
}
