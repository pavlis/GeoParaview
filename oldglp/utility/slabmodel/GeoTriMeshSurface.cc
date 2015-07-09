#include <boost/smart_ptr.hpp>
#include "coords.h"
#include "gclgrid.h"
#include "SeisppError.h"
#include "GeoTriMeshSurface.h"
using namespace std;
using namespace SEISPP;

GeoTriMeshSurface::GeoTriMeshSurface(vector<Geographic_point> pts,string coordunits)
{
    if(pts.size()<3) 
        throw SeisppError(string("GeoTriMeshSurface:")
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
    trigrid=boost::shared_ptr<CGPartition>(delaunay_On4(pointset));
    cgpointset_free(&pointset);
}
GeoTriMeshSurface::GeoTriMeshSurface(vector<double>lat, vector<double>lon, 
        vector<double> depth,string coordunits)
{
    int npts=lat.size();
    if( (lon.size()!=npts) || (depth.size()!=npts)) throw SeisppError(
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
    trigrid=parent.trigrid;
    units=parent.units;
}
GeoTriMeshSurface& GeoTriMeshSurface::operator=(const GeoTriMeshSurface& parent)
{
    if(this!=&parent)
    {
        trigrid=parent.trigrid;
        units=parent.units;
    }
    return(*this);
}
double GeoTriMeshSurface::radius(double lat, double lon)
{
    /* We need to fetch the raw pointer from the shared_ptr by this method */
    CGPartition *cgpptr=trigrid.get();
    CGPolygon *poly=polygon_containing_xy(cgpptr,lon,lat);
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
