#ifndef _GEOSURFACE_H_
#define _GEOSURFACE_H_

class GeoSurface
{
public:
    virtual double radius(double lat, double lon)=0;
    virtual double depth(double lat, double lon)=0;
    virtual bool is_defined(double lat,double lon)=0;
};
#endif
