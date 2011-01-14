#include "PLGeoPath.h"
#include "GeoCoordError.h"
#include "interpolator1d.h"
using namespace INTERPOLATOR1D;
PLGeoPath::PLGeoPath() : coordxyz(0.0,0.0,0.0,0.0)
{
    lat.reserve(1);
    lon.reserve(1);
    r.reserve(1);
    s.reserve(1);
    s0.lat=0.0;
    s0.lon=0.0;
    s0.r=6371.0;
}

PLGeoPath::PLGeoPath(vector<Geographic_point>& pts,int i0,double az0)
    : coordxyz(pts[i0].lat,pts[i0].lon,pts[i0].r,az0)
{
    const string base_error("PLGeoPath constructor:  ");
    int npts=pts.size();
    if(i0<0) throw GeoCoordError(base_error
                +"negative origin position index not allowed");
    if(i0>=npts) throw GeoCoordError(base_error
                +"origin position index exceed input array size");
    s0=pts[i0];
    vector<Geographic_point>::iterator ptsptr;
    for(ptsptr=pts.begin();ptsptr!=pts.end();++ptsptr)
    {
        lat.push_back(ptsptr->lat);
        lon.push_back(ptsptr->lon);
        r.push_back(ptsptr->r);
        /* This is a simple way to initialize this container */
        s.push_back(0.0);
    }
    int i;
    const double EarthRadius(6378.164);
    double ddelta,ddkm,az;
    for(i=i0+1;i<npts;++i)
    {
        dist(lat[i-1],lon[i-1],lat[i],lon[i],&ddelta,&az);
        ddkm=ddelta*EarthRadius;
        s[i]=s[i-1]+hypot(ddkm,r[i]-r[i-1]);
    }
    for(i=i0-1;i>=0;--i)
    {
        dist(lat[i+1],lon[i+1],lat[i],lon[i],&ddelta,&az);
        ddkm=ddelta*EarthRadius;
        s[i]=s[i+1]-hypot(ddkm,r[i]-r[i-1]);
    }
}
PLGeoPath::PLGeoPath(const PLGeoPath& parent) : coordxyz(parent.coordxyz)
{
    lat=parent.lat;
    lon=parent.lon;
    r=parent.r;
    s=parent.s;
    s0=parent.s0;
}
PLGeoPath& PLGeoPath::operator=(const PLGeoPath& parent)
{
    if(this!=&parent)
    {
        lat=parent.lat;
        lon=parent.lon;
        r=parent.r;
        s=parent.s;
        s0=parent.s0;
        coordxyz=parent.coordxyz;
    }
    return(*this);
}
/* This method should do a linear extrapolation for points beyond
   the endpoints */
Geographic_point PLGeoPath::position(double sp)
{
    Geographic_point result;
    int npts=s.size();
    int is;
    is=irregular_lookup(sp,&(s[0]),npts);
    if(is<=0)
    {
        result.lat=linear_scalar(sp,s[0],lat[0],s[1],lat[1]);
        result.lon=linear_scalar(sp,s[0],lon[0],s[1],lon[1]);
        result.r=linear_scalar(sp,s[0],r[0],s[1],r[1]);
    }
    else if(is>=(npts-2))
    {
        result.lat=linear_scalar(sp,s[npts-2],lat[npts-2],
                s[npts-1],lat[npts-1]);
        result.lon=linear_scalar(sp,s[npts-2],lon[npts-2],
                s[npts-1],lon[npts-1]);
        result.r=linear_scalar(sp,s[npts-2],r[npts-2],
                s[npts-1],r[npts-1]);
    }
    else
    {
        result.lat=linear_scalar(sp,s[is],lat[is],s[is+1],lat[is+1]);
        result.lon=linear_scalar(sp,s[is],lon[is],s[is+1],lon[is+1]);
        result.r=linear_scalar(sp,s[is],r[is],s[is+1],r[is+1]);
    }
    return(result);
}
Cartesian_point  PLGeoPath::position_xyz(double sp)
{
    Geographic_point gp=this->position(sp);
    return(this->coordxyz.cartesian(gp));
}
double PLGeoPath::latitude(double sp)
{
    Geographic_point gp=this->position(sp);
    return(gp.lat);
}
double PLGeoPath::longitude(double sp)
{
    Geographic_point gp=this->position(sp);
    return(gp.lon);
}

Geographic_point PLGeoPath::origin()
{
    return(s0);
}
Geographic_point PLGeoPath::begin()
{
    Geographic_point result;
    result.lat=lat[0];
    result.lon=lon[0];
    result.r=r[0];
    return(result);
}
Geographic_point PLGeoPath::end()
{
    int np=lat.size();
    --np;
    Geographic_point result;
    result.lat=lat[np];
    result.lon=lon[np];
    result.r=r[np];
    return(result);
}
double PLGeoPath::sbegin()
{
    return(s[0]);
}
double PLGeoPath::send()
{
    return(s[lat.size()-1]);
}
Geographic_point PLGeoPath::node_position(int i)
{
    if(i<0 || i>=(this->number_points()) )
        throw GeoCoordError(string("PLGeoPath::node_position method:  ")
                + "requested point outside range of node points");
    Geographic_point result;
    result.lat=lat[i];
    result.lon=lon[i];
    result.r=r[i];
    return result;
}

ostream& operator<<(ostream& os,PLGeoPath& path)
{
    int i;
    /* assume here all the vectors are the same length */
    int npts=path.lat.size();
    // silently do nothing if npts is le 1
    if(npts>1)
    {
	    for(i=0;i<npts;++i)
	    {
	        os << deg(path.lon[i])<<" "
	            << deg(path.lat[i])<<" "
	            << r0_ellipse(path.lat[i])-path.r[i]<<" "
	            << path.r[i] << endl;
	    }
    }
}

