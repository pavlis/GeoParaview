#include <math.h>
#include "PlateBoundaryPath.h"
#include "SeisppError.h"
PlateBoundaryPath::PlateBoundaryPath(double pla, double plo,
        double ola, double olo,double avel)
{
    plat=pla;
    plon=plo;
    olat=ola;
    olon=olo;
    dist(plat,plon,olat,olon,&ssdelta,&azimuth0);
    angular_velocity=avel;
}
PlateBoundaryPath::PlateBoundaryPath(const PlateBoundaryPath& parent)
{
    plat=parent.plat;
    plon=parent.plon;
    olat=parent.olat;
    olon=parent.olon;
    ssdelta=parent.ssdelta;
    angular_velocity=parent.angular_velocity;
}
PlateBoundaryPath& PlateBoundaryPath::operator=(const PlateBoundaryPath& parent)
{
    if(this!=&parent)
    {
        plat=parent.plat;
        plon=parent.plon;
        olat=parent.olat;
        olon=parent.olon;
        ssdelta=parent.ssdelta;
        angular_velocity=parent.angular_velocity;
    }
    return(*this);
}

Geographic_point PlateBoundaryPath::position(double s)
{
    Geographic_point result;
    if(fabs(s)<FLT_EPSILON) 
    {
        result.lat=olat;
        result.lon=olon;
        result.r=r0_ellipse(olat);
        return(result);
    }
    else
    {
        double daz=s*angular_velocity;
        if(daz>(2.0*M_PI)) throw SEISPP::SeisppError(
            string("PlateBoundaryPath::position method: ")
            + "time parameter passed yields angle greater than 360 degrees");

        daz += azimuth0;
        latlon(plat,plon,ssdelta,daz,&(result.lat),&(result.lon));
        result.r=r0_ellipse(result.lon);
    }
    return result;
}
double PlateBoundaryPath::latitude(double s)
{
    Geographic_point gp=this->position(s);
    return(gp.lat);
}
double PlateBoundaryPath::longitude(double s)
{
    Geographic_point gp=this->position(s);
    return(gp.lon);
}
Geographic_point PlateBoundaryPath::origin()
{
    Geographic_point result;
    result.lat=olat;
    result.lon=olon;
    result.r=r0_ellipse(olat);
    return(result);
}
