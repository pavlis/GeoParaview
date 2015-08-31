#ifndef _MULTIMASK_H_
#define _MULTIMASK_H_
#include "GeoPolygonRegion.h"
/*! \brief Object to define a masked region define by multiple
  overlapping polygons.  Polygons define region inside or outside that is
  to be kept. Whether a point is defined as in or out is determined by the 
  logical and of all polygon inside/outside tests. */
class MultiMask
{
    public:
        /*! Construct from an input stream */
        MultiMask(istream& in);
        MultiMask(const MultiMask& parent);
        MultiMask& operator=(const MultiMask& parent);
        bool is_inside(double lat0, double lon0);
    private:
        /* An array of polygons that define the masking geometry. */
        vector<GeoPolygonRegion> polygon;
        /* This is a parallel vector to polygon and defines polarity
           of tests applied to a point.  polygons with true use the 
           is_inside test while polygons defined false use !is_inside.
           */
        vector<bool> inside;   
};
#endif
