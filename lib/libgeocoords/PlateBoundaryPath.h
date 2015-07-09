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
be in years, although the only parameter it interacts with is the 
angular velocity in the constructor.  That is, units of angular 
velocity are radians/time so as long as the time units are consistent
the distances computed will be accurate.  Applications I have used
always use years. 
 
  Note finally that the time parameterization of this object does not
need to correspond to earth history time.  This is possible because
a single rotation pole produces a path that is uniquely defined by 
the pole and any point on the path.  This disconnect from reality can
be used, for example, to draw a trajectory of a flow line for a plate
motion model through a point other data show is on the plate boundary.*/

class PlateBoundaryPath : public GeoPath
{
public:
    /* \brief Primary constructor.

       This is the only constructor for this object.
       \param pla is the latitude of the pole of rotation (radians)
       \param plo is the longitude of the pole of rotation (radians)
       \param ola is the latitude of the origin of the path (radians) where t=0
       \param olo is the longitude of the origin of the path (radians) where t=0
       \param avel is the angular velocity (radians/year)
       */
    PlateBoundaryPath(double pla, double plo, double ola, double olo,double avel);
    /*! Standard copy constructor. */
    PlateBoundaryPath(const PlateBoundaryPath& parent);
    /*! Standard assignment operator. */
    PlateBoundaryPath& operator=(const PlateBoundaryPath& parent);
    /* \brief return geographic coordinates at a requested time. 

       This is the core method of this object.  It returns the position for
       a particular time along the boundary curve.  The Geographic_point object
       is defined in gclgrid.h.

       \param t is the time for which a position is requested.  It can be positive
         or negative since a single pole of rotation has complete symmetry.
       \exception Throws a GeoCoordError exception if angle computed from avel*t 
         is greater than 360 degrees.
         */
    Geographic_point position(double t);
    /* \brief Return the latitude (in radians) of the path at time t. 

       This method is strictly for convenience. It calls the position method and 
       extract and returns the latitude field. 
       \param t is the time for which a position is requested.  It can be positive
         or negative since a single pole of rotation has complete symmetry.
       \exception Throws a GeoCoordError exception if angle computed from avel*t 
         is greater than 360 degrees.
         */
    double latitude(double t);
    /* \brief Return the longitude (in radians) of the path at time t. 

       This method is strictly for convenience. It calls the position method,
       extract the longitude variable, and returns it.
       \param t is the time for which a position is requested.  It can be positive
         or negative since a single pole of rotation has complete symmetry.
       \exception Throws a GeoCoordError exception if angle computed from avel*t 
         is greater than 360 degrees.
         */
    double longitude(double t);
    /*! Return the origin parameterized by t=0. */
    Geographic_point origin();
    /*! Return distance along path (km) at time t in km*/
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
  small circles (latitude like lines) relative to a pole of rotation.  Most
  plates have variable time history that is specified as a set of stage poles.
  This object is a variant of PlateBoundaryPath but for plates with stage
  poles models for longer term motion. 
  */
class TimeVariablePlateBoundaryPath : public GeoPath
{
public:
    /*! Full constructor from vectors of stage pole data.

      This sets internal data that define this objects behavior using
      a set of parallel vectors.  Note all angles given as inputs to this
      constructor are assumed to be in radians. Vectors are assume parallel.
      By that I mean they are the same length and the ith value of one vector
      is linked to the ith value of the other.

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
      \param olon is the longitude of the origin (first point) on a path to 
        be defined by this constructor.  
        */
    TimeVariablePlateBoundaryPath(vector<double>spla, vector<double> splo,
            vector<double> times, vector<double> ang, double ola, double olo);
    /*! Construct from gmt format stage pole data. 

      This constructor builds this object from an ascii file of stage pole data.  The 
      format is the same as that used in GMT programs rotconverter, backprojector, and
      others.  That is stage poles are listed from oldest to youngest with five, 
      white space seperated fields:  longitude, latitude, tstart, tend, angle.  All
      angle quantities are in degrees and times are in Myr.  Note degrees are always
      converted internally to radians in this object. 

      \param fname is the file name to be read
      \param ola latitude of origin of path to be projected (radians)
      \param olo longitude of origin of path to be projected (radians)

      \exception throws a SeisppError exception if the file cannot be opened or if the 
        times are not in decreasing time order. 
        */
    TimeVariablePlateBoundaryPath(string fname,double ola, double olo);

    /*! Standard copy constructor. */
    TimeVariablePlateBoundaryPath(const TimeVariablePlateBoundaryPath& parent);
    /*! Standard assignment operator. */
    TimeVariablePlateBoundaryPath& operator=(const TimeVariablePlateBoundaryPath& parent);
    /* \brief return geographic coordinates at a requested time. 

       This is the core method of this object.  It returns the position for
       a particular time along the boundary curve.  The Geographic_point object
       is defined in gclgrid.h.

       \param t is the time for which a position is requested.  It can be positive
         or negative since a single pole of rotation has complete symmetry.
       \exception Throws a GeoCoordError exception if time is negative.  
         For stage pole data this makes no sense.
         */
    Geographic_point position(double t);
    /* \brief Return the latitude (in radians) of the path at time t. 

       This method is strictly for convenience. It calls the position method,
       extract the latitude field, and returns the result.  It is thus very
       inefficent to use this method if both latitude and longitude are needed.
       \param t is the time for which a position is requested.  It can be positive
         or negative since a single pole of rotation has complete symmetry.
       \exception Throws a GeoCoordError exception if time is negative.  
         For stage pole data this makes no sense.
         */
    double latitude(double t);
    /* \brief Return the longitude (in radians) of the path at time t. 

       This method is strictly for convenience. It calls the position method,
       extract the longitude field, and returns the result.  It is thus very
       inefficent to use this method if both latitude and longitude are needed.
       \param t is the time for which a position is requested.  It can be positive
         or negative since a single pole of rotation has complete symmetry.
       \exception Throws a GeoCoordError exception if time is negative.  
         For stage pole data this makes no sense.
         */
    double longitude(double t);
    /*! Return the origin of the path. */
    Geographic_point origin();
    /*! Returns the arc path distance in km for elapsed time t in years */
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
