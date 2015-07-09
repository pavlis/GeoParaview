#ifndef _PLATEBOUNDARYPATH_H_
#define _PLATEBOUNDARYPATH_H_
#include <float.h>
#include "gclgrid.h"
#include "GeoPath.h"
/*! Object to project a small circle path defining a plate boundary.

  In plate tectonics transform plate boundaries and motion vectors
are defined by small circles (latitude like lines) relative to a pole
of rotation.  This object simplifies the task of creating these
lines.  Curves are parameterized by time.  Normally that would
be in years. */

class PlateBoundaryPath : public GeoPath
{
public:
    PlateBoundaryPath(double pla, double plo, double ola, double olo,double avel);
    PlateBoundaryPath(const PlateBoundaryPath& parent);
    PlateBoundaryPath& operator=(const PlateBoundaryPath& parent);
    /* Return coordinates as a GCLgrid Geographic_point object at elapsed time s
       (in years) assuming origin is time 0 */
    Geographic_point position(double t);
    double latitude(double t);
    double longitude(double t);
    Geographic_point origin();
    /*! Return distance at time t */
    double distance(double t); 
    /*! Return time at arc distance s in km */
    double time(double s);
private:
    // pole and origin coordinates (radians)
    double plat,plon,olat,olon;
    double ssdelta;  // small circle distance from pole to origin (radians)
    double azimuth0;  // azimuth at pole to point define as origin for path
    double angular_velocity;  //angular velocity in radians
};
/*! Object to define flow line paths with variable history.

  In plate tectonics transform boundaries and motion vectors are defined by
  small circles (latitude like lines) relative to a pole of rotation.  Some
  plates have variable time history that is specified as a set of stage poles.
  This object is a variant of PlateBoundaryPath but for plates with stage
  poles. 
  */
class TimeVariablePlateBoundaryPath : public GeoPath
{
public:
    /*! Full constructor from vectors of stage pole data.

      This sets internal data that define this objects behavior using
      a set of parallel vectors.  Note all angles given as inputs to this
      constructor are assumed to be in radians.

      \param spla is a vector of stage pole latitudes
      \param splo is a vector of stage pole longitudes.  Assumed to be same 
        length as spla
      \param times is a vector of the time interval covered by a stage pole.  Note
        carefully this is a time INTERVAL not a time from present.  The stage poles
        are presumed ordered in increasing time from present.  Thus, for example,
        the current pole is presumed to be spla[0],splo[0] and be appropriate until
        dt before present.  UNITS of times is YEARS.  

      \param ang is a vector of angular displacements by stage pole spla[i],
        splo[i] in the increment of time dt[i].  
      \param olat is the latitude of the origin (first point) on a path to 
        be defined by this constructor.  
      \param olon is the longidue of the origin (first point) on a path to 
        be defined by this constructor.  
        */
    TimeVariablePlateBoundaryPath(vector<double>spla, vector<double> splo,
            vector<double> times, vector<double> ang, double ola, double olo);
    TimeVariablePlateBoundaryPath(const TimeVariablePlateBoundaryPath& parent);
    TimeVariablePlateBoundaryPath& operator=(const TimeVariablePlateBoundaryPath& parent);
    /* Return coordinates as a GCLgrid Geographic_point object at time t 
       (in years) */
    Geographic_point position(double t);
    double latitude(double t);
    double longitude(double t);
    Geographic_point origin();
    /* Returns the arc path distance in km for elapsed time t in years */
    double distance(double t);
    /*! Return time at arc distance s in km */
    double time(double s);
private:
    vector<double> splat, splon,dt;
    /* These become vectors to simplify path calculation */
    vector<double> olat,olon;
    /* This vector contains path arc length (radians of the longitude 
       like angle subtended) to each interval point. 
       Thus ss[0]=0, ss[1] is distance to olat[1],olon[1] from olat[0],olon[0],
       etc accumulating for each successive stage pole. */
    vector<double> ss;
    /* This holds small circle distance (radians) from each stage pole to the 
       first point on the path defined by that pole. Same length as
       stage pole vectors */
    vector<double> ssdelta;
    /* Similar vector of azimuth at pole to point defined as origin */
    vector<double> azimuth0;
    /* Angular velocities in radians/year.  Just angles[i]/dt[i],
      but more convenient than storing angles. */
    vector<double>  omegadot;  //angular velocity in radians
    /* Integrated time to origin for each stage pole segment */
    vector<double>  t0;
};
#endif
