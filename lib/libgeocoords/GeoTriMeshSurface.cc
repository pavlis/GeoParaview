#include "coords.h"
#include "gclgrid.h"
#include "GeoCoordError.h"
#include "GeoTriMeshSurface.h"
using namespace std;

GeoTriMeshSurface::GeoTriMeshSurface(vector<Geographic_point> pts,string coordunits)
{
    if(pts.size()<3) 
        throw GeoCoordError(string("GeoTriMeshSurface:")
              + "  insufficent points.  Need at least 3");
    if(coordunits=="radians")
        units=RADIANS;
    else
        units=DEGREES;
    vector<Geographic_point>::iterator ppts;
    CGPointset *pointset=cgpointset_new();
    for(ppts=pts.begin();ppts!=pts.end();++ppts)
    {
        CGPoint *cgpt=cgpoint_new(ppts->lon,ppts->lat,ppts->r);
        cgpointset_addpoint(pointset,cgpt);
        cgpoint_free(&cgpt);
    }
    trigrid=delaunay_On4(pointset);
    cgpointset_free(&pointset);
}
GeoTriMeshSurface::GeoTriMeshSurface(vector<double>lat, vector<double>lon, 
        vector<double> depth,string coordunits)
{
    int npts=lat.size();
    if( (lon.size()!=npts) || (depth.size()!=npts)) throw GeoCoordError(
                "GeoTriMeshSurface constructor:  parallel lat,lon,depth arrays must be same size");
    if(coordunits=="radians")
        units=RADIANS;
    else
        units=DEGREES;
    int i;
    vector<Geographic_point> pts;
    pts.reserve(npts);
    Geographic_point p;
    for(i=0;i<npts;++i)
    {
        p.lat=lat[i];
        p.lon=lon[i];
        if(units==RADIANS)
            p.r=r0_ellipse(lat[i]);
        else
            p.r=r0_ellipse(rad(lat[i]));
        p.r-=depth[i];
        pts.push_back(p);
    }
    try {
        if(units==RADIANS)
            *this=GeoTriMeshSurface(pts);
        else
            *this=GeoTriMeshSurface(pts,string("degrees"));
    } catch(...){throw;};
}
GeoTriMeshSurface::GeoTriMeshSurface(const GeoTriMeshSurface& parent)
{
    throw GeoCoordError(string("GeoTriMeshSurface::copy constructor is not implemented yet"));
    /* This section is impossible to complete until I get info 
       from Kent */
    /*
    void temp;
    trigrid=cgpartition_new(temp);
    trigrid->polygons
    */
    // This section will not work correctly since trigrid is a pointer.
    /*
    trigrid=parent.trigrid;
    */
    units=parent.units;
}
GeoTriMeshSurface::~GeoTriMeshSurface()
{
    cgpartition_free(&trigrid);
}
GeoTriMeshSurface& GeoTriMeshSurface::operator=(const GeoTriMeshSurface& parent)
{
    throw GeoCoordError(string("GeoTriMeshSurface::opertor= is not implemented yet"));
    /*
    if(this!=&parent)
    {
        trigrid=parent.trigrid;
        units=parent.units;
    }
    return(*this);
    */
}
void GeoTriMeshSurface::AddBoundary(const GeoPolygonRegion& poly)
{
    boundary=poly;
}
double GeoTriMeshSurface::radius(double lat, double lon)
{
    /* We need to fetch the raw pointer from the shared_ptr by this method */
    CGPolygon *poly=polygon_containing_xy(trigrid,lon,lat);
    if(poly==NULL) return(a_nan());
    return(cg_triinterp(poly,lon,lat));
}
double GeoTriMeshSurface::depth(double lat, double lon)
{
    double r0;
    if(units==RADIANS)
        r0=r0_ellipse(lat);
    else
        r0=r0_ellipse(rad(lat));
    double r=this->radius(lat,lon);
    return(r0-r);
}
bool GeoTriMeshSurface::is_defined(double lat, double lon)
{
    CGPolygon *poly=polygon_containing_xy(trigrid,lon,lat);
    if(poly==NULL) 
        return false;
    if(boundary.is_inside(lat,lon))
        return true;
    else
        return false;
}
