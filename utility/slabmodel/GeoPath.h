#ifndef _GEOPATH_H_
#define _GEOPATH_H_

#include "gclgrid.h"
/*! Abstract base class for paths in Earth coordinates.

  One of a class of geometric objects it is sometimes useful to describe in the
Earth is a curved path.  This is a general definition of such a curve defined by
a parameteric model with an origin and a generalized coordinate s that defines 
distance along that path.  This is the generic interface with no assumptions 
about the units of s or how the methods are implemented.  It is used to define
a pure abstract base class. */
class GeoPath
{
public:
    virtual Geographic_point position(double s)=0;
    virtual double latitude(double s)=0;
    virtual double longitude(double s)=0;
    virtual Geographic_point origin()=0;
};
#endif
