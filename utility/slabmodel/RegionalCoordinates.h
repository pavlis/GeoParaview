#ifndef _REGIONALCOORDINATES_H_
#define _REGIONALCOORDINATES_H_

#include "gclgrid.h"
class RegionalCoordinates
{
public:
    RegionalCoordinates();
    RegionalCoordinates(double olat, double olon, double oradius, double azn);
    RegionalCoordinates(const RegionalCoordinates& parent);
    RegionalCoordinates& operator=(const RegionalCoordinates& parent);
    Cartesian_point  cartesian(double lat,double lon, double r);
    Cartesian_point  cartesian(Geographic_point gp);
    Cartesian_point cartesian(double x[3]);
    Geographic_point geographic(double x1,double x2, double x3);
    Geographic_point geographic(double x[3]);
    Geographic_point geographic(Cartesian_point cp);
    Geographic_point origin();
    double aznorth_angle(){return(azimuth_y);};
private:
    double gtoc_rmatrix[3][3];
    double translation_vector[3];
    double lat0;
    double lon0;
    double r0;
    double azimuth_y;
};
#endif
