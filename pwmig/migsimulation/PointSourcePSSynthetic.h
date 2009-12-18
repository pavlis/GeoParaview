#ifndef _POINTSOURCEPSYNTHETIC_H_
#define _POINTSOURCEPSYNTHETIC_H_
#include "SyntheticSeismogram.h"

namespace SEISPP {
using namespace std;
using namespace SEISPP;
class SphericalRayPathArray
{
    public:
        SphericalRayPathArray(VelocityModel1d& vmod,vector<double>& p,
                double zminin, double zmax, double dz);
        SphericalRayPathArray(VelocityModel1d& vmod,double p0, double dp, int np,
                double zminin, double zmax, double dz);
        SphericalRayPathArray(const SphericalRayPathArray& parent);
        pair<double,double> time_and_amp(double delta, double z);
        SphericalRayPathArray& operator=(const SphericalRayPathArray& parent);
    private:
        int nrays;  // size of rays, cached for efficiency
        vector<RayPathSphere> rays;
        vector<double> slow;  // vector of ray parameters (s/km)
        /* zmin is the minimum depth allowed and a reference depth used for amplitude.
           That is, the amplitude for a ray from depth zmin to zero offset is given a unit
           amplitude factor */
        double zmin;  
        int nz;  // nominal depth steps that can be dependent on
        double dz;  // rays are made in equal depth steps, save here.
        dmatrix amps;  // nz x nrays matrix of amplitudes derived from rays
        double zfloor;  // actual depth floor = dz*(nz-1) cached for efficiency
};

class PointSourcePSSynthetic : public SyntheticSeismogram
{
    public:
        PointSourcePSSynthetic(VelocityModel_1d& vmod,Pf *pf);
        //~PointSourcePSSynthetic();
        TimeSeries ComputeScalar(Hypocenter& hypo,
                            double rlat, double rlon, double relev,string type);
        TimeSeries ComputeScalar(const TimeSeries& parent,
                    Hypocenter& hypo,
                            double rlat, double rlon, double relev,string type);
        ThreeComponentSeismogram Compute3C(Hypocenter& hypo,
                 double rlat, double rlon, double relev,string units);
        ThreeComponentSeismogram Compute3C(const ThreeComponentSeismogram& parent,
                 Hypocenter& hypo,
                 double rlat, double rlon, double relev,string units);
    private:
        SphericalRayPathArray rays;
        /* purists might want this to be a single object, but here
           I use two parallel vectors to avoid this unnecessary
           complexity for a private variable.   Could easily change
           implementation of that proved essential */
        vector<Geographic_point> points;
        vector<double> amp;
        /* If a station is located beyond this cutoff distance (radians)
           will not try to compute the point source response.  This is
           essential because we don't handle turning rays and for
           the for which I wrote this is links with the concept of 
           cutoff aperture in pwstack and (in progress) pwdecon*/
        double distance_cutoff;
};
}  // End SEISPP namespace encapsulation
#endif
