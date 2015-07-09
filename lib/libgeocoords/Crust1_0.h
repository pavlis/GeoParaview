#include <string>
#include "VelocityModel_1d.h"
using namespace std;
class Crust1_0
{
public:
    Crust1_0();
    //Crust1_0(const Crust1_0& parent);
    /*! \brief Return an interpolated 1d model at a specified point.

      This method returns a model component as VelocityModel_1d object
      at a specified point.   The model type is defined by the property
      argument.

      \param lat - Latitude of point to get model (degrees)
      \param lon - Longitude of point to get model (degrees)
      \param property - Must be "Pvelocity" or "Vp", "Svelocity" or "Vs", 
         or "Density or "rho""

      \exception - throws a SeisppError for an invalid property or illegal
         point.
      */
    VelocityModel_1d model(double lat, double lon, string property);
    double CrustalThickness(double lat, double lon);
    double MohoDepth(double lat, double lon);
    double SedimentThickness(double lat, double lon);
    double WaterDepth(double lat, double lon);
    double IceThickness(double lat, double lon);
private:
    /* This set of parameters is presently frozen and defined by the 
       data files for crust1.0.   They are set in the constructor 
       and not read because they are defined only by a readme file.*/
    int nlon,nlat,nlayers; // grid dimensions
    double lat0,lon0;  // location of [0][0] point
    double dlon,dlat;  // grid spacing in degrees
    /* These simplify lookup.  They return the integer for binlinear
       interpolation of lon/lat */
    int lookup_lon(double lon);
    int lookup_lat(double lat);
    double ***vp;
    double ***vs;
    double ***rho;
    double ***depth;
    double ***thickness;
    void load_crust1file(double ***x,string fname);
};
/* Generic algorithm to do binlinear interpolator of a gridded 2d data.  

   (tx,ty) defines the point to be interpolated in normalized coordinates.
That is, interpolator is for unit square so tx and ty must be between 0 and
1.  Since this is a low level function we intentionally do not test for
validity of the points for speed, but the caller should handle this. 

The c values are the values of the function to be interpolated at the 
corners.   c00 is lower left, c10 lower right, c01 upper left, and
c11 upper right (normal x,y coordinate convention). 

This template algorithm was adapted from
http://www.scratchapixel.com/lessons/3d-advanced-lessons/interpolation/bilinear-interpolation/
*/
template <class T> T bilinear(const T& tx, const T& ty,
       const T& c00, const T& c10, const T& c01, const T& c11)
{
    return (T(1) - tx) * (T(1) - ty) * c00 
        + tx * (T(1) - ty) * c10 + (T(1) - tx) * ty * c01 + tx * ty * c11;
}
