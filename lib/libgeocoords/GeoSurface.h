#ifndef _GEOSURFACE_H_
#define _GEOSURFACE_H_
/*! Abstract base class for georeferenced surfaces.

  We often want a georeferenced surface.  This is the 
basic interface to such an object in this library.  
This is a pure abstract class. It is intended for use
only with single valued surfaces.  That means, at a given
lat,lon position there either only one or no point on the surface
at that position.
 
 This only defines the interface.  The units of the inputs are not
 defined in this interface, although all current children assume
 lat and lon values in radians and distances in km*/

class GeoSurface
{
public:
    /*! Return the radius of the surface at point lat,lon.*/`
    virtual double radius(double lat, double lon)=0;
    /*! Return the depth of the surface at point lat,lon.*/`
    virtual double depth(double lat, double lon)=0;
    /*! Test to see if the surface is defined at this point.  Returns 
      true of the surface is defined at the given point. */
    virtual bool is_defined(double lat,double lon)=0;
};
#endif
