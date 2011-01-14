#include <math.h>
#include <fstream>
#include <sstream>>
#include "PlateBoundaryPath.h"
#include "GeoCoordError.h"
using namespace std;
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
        if(daz>(2.0*M_PI)) throw GeoCoordError(
            string("PlateBoundaryPath::position method: ")
            + "time parameter passed yields angle greater than 360 degrees");

        daz += azimuth0;
        latlon(plat,plon,ssdelta,daz,&(result.lat),&(result.lon));
        result.r=r0_ellipse(result.lat);
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
double PlateBoundaryPath::distance(double t)
{
    double phi=angular_velocity*t;
    return(phi*ssdelta);
}
double PlateBoundaryPath::time(double s)
{
    double phi=s/ssdelta;
    return(phi/angular_velocity);
}
/*  Now begin section on TimeVariablePlateBoundary object */
TimeVariablePlateBoundaryPath::TimeVariablePlateBoundaryPath(
        vector<double>spla, vector<double> splo,
            vector<double> times, vector<double> ang, double ola, double olo)
        : splat(spla), splon(splo), dt(times)
{
    string base_error("TimeVariablePlateBoundaryPath constructor:  ");
    if( (spla.size()!=splo.size()) || (splo.size()!=times.size())
            || (times.size() != ang.size()) )throw GeoCoordError(
                base_error+" input vector size mismatch");
    /* Push the origin as 0 point of this interval vector */
    olat.push_back(ola);
    olon.push_back(olo);
    int i,npoints;
    npoints=splat.size();
    for(i=0;i<npoints;++i)
        omegadot.push_back(ang[i]/dt[i]);
    for(i=0;i<npoints;++i)
    {
        double lat,lon,r;
        double rto0;  // distance in radians to origin for this segment
        double az;  // temporary azimuth output of dist function
        dist(splat[i],splon[i],olat[i],olon[i],&rto0,&az);
        ssdelta.push_back(rto0);
        azimuth0.push_back(az);
        // Now compute the lat lon of the next segment
        double nextaz,daz;
        daz=dt[i]*omegadot[i];
        ss.push_back(daz);
        nextaz=az+daz;
        latlon(splat[i],splon[i],rto0,nextaz,&lat,&lon);
        olat.push_back(lat);
        olon.push_back(lon);
    }
    /* Clearer to just compute these now.  Note carefully t0 will 
    have one more point than dt.  This is the classic interval
    versus points problem. */
    t0.push_back(0.0);
    double t;
    for(i=0;i<npoints;++i)
    {
        t+=dt[i];
        t0.push_back(t);
    }
}
TimeVariablePlateBoundaryPath::TimeVariablePlateBoundaryPath(string fname,
        double ola, double olo)
{
    const string base_error("TimeVariablePlateBoundaryPath ascii file constructor:  ");
    ifstream infile;
    infile.open(fname.c_str(),ifstream::in);
    if(infile.fail()) throw GeoCoordError(base_error
            + "open failed for file name=" + fname);
    vector<double> ilat,ilon,its,ite,iang;
    char line[512];
    while(infile.getline(line,512))
    {
        stringstream ss(line);
        double dval;
        ss >> dval;
        ilon.push_back(rad(dval));
        ss >> dval;
        ilat.push_back(rad(dval));
        ss >> dval;
        ite.push_back(dval);
        ss >> dval;
        its.push_back(dval);
        ss >> dval;
        iang.push_back(rad(dval));
    }
    infile.close();
    /* Check for input errors.  For first test iang is intentionally tested 
    to handle one line file without enough data. */
    if(iang.size()==0) throw GeoCoordError(base_error
            + "file "+fname+" seems to be empty or be improperly formatted");
    if(ilon.size() != iang.size() ) throw GeoCoordError(base_error
            + "size mistmatch.  File "+fname
            + " probably improperly formatted.");
    for(int i=1;i<ilon.size();++i)
    {
        if(fabs(its[i-1]-ite[i])>0.01) throw GeoCoordError(base_error
                + "Mismatched start and end time intervals. Require matching time intervals.");
    }
    /* We have to reverse the order of the vectors to use the parallel vector constructor
       and to convert times to positive time intervals.  */
    vector<double>irlon,irlat,irdt,irang;
    for(int i=ilon.size()-1;i>0;--i)
    {
        irlon.push_back(ilon[i]);
        irlat.push_back(ilat[i]);
        irang.push_back(iang[i]);
        irdt.push_back(ite[i]-its[i]);
    }
    *this=TimeVariablePlateBoundaryPath(irlat,irlon,irdt,irang,ola,olo);
}


TimeVariablePlateBoundaryPath::TimeVariablePlateBoundaryPath
        (const TimeVariablePlateBoundaryPath& parent)
{
    splat=parent.splat;
    splon=parent.splon;
    dt=parent.dt;
    olat=parent.olat;
    olon=parent.olon;
    ss=parent.ss;
    ssdelta=parent.ssdelta;
    azimuth0=parent.azimuth0;
    omegadot=parent.omegadot;
    t0=parent.t0;
}
TimeVariablePlateBoundaryPath& TimeVariablePlateBoundaryPath::operator=(
        const TimeVariablePlateBoundaryPath& parent)
{
    if(this!=&parent)
    {
        splat=parent.splat;
        splon=parent.splon;
        dt=parent.dt;
        olat=parent.olat;
        olon=parent.olon;
        ss=parent.ss;
        ssdelta=parent.ssdelta;
        azimuth0=parent.azimuth0;
        omegadot=parent.omegadot;
        t0=parent.t0;
    }
    return(*this);
}
Geographic_point TimeVariablePlateBoundaryPath::position(double t)
{
    const string base_error("TimeVariablePlateBoundaryPath::position:  ");
    Geographic_point result;
    if(t<0) throw GeoCoordError(base_error
                + "illegal negative time parameter. Must be nonnegative");
    /* Linear search for first position */
    int npoints=splat.size();
    int i;
    /* Intentionally start at 1 since we will decrement this at end */
    for(i=1;i<npoints;++i)
        if(t0[i]>t) break;
    --i;  
    double dt=t-t0[i];
    double daz=dt*omegadot[i];
    latlon(splat[i],splon[i],ssdelta[i],azimuth0[i]+daz,
            &(result.lat),&(result.lon));
    result.r=r0_ellipse(result.lat);
    return result;
}
double TimeVariablePlateBoundaryPath::latitude(double s)
{
    Geographic_point gp=this->position(s);
    return(gp.lat);
}
double TimeVariablePlateBoundaryPath::longitude(double s)
{
    Geographic_point gp=this->position(s);
    return(gp.lon);
}
Geographic_point TimeVariablePlateBoundaryPath::origin()
{
    Geographic_point result;
    result.lat=olat[0];
    result.lon=olon[0];
    result.r=r0_ellipse(olat[0]);
    return(result);
}
double TimeVariablePlateBoundaryPath::distance(double t)
{
    const string base_error("TimeVariablePlateBoundaryPath::distance:  ");
    if(t<0) throw GeoCoordError(base_error
                + "illegal negative time parameter. Must be nonnegative");
    /* parallel code to position method */
    int npoints=splat.size();
    int i;
    /* Intentionally start at 1 since we will decrement this at end */
    for(i=1;i<npoints;++i)
        if(t0[i]>t) break;
    --i;  
    double dt=t-t0[i];
    double daz=dt*omegadot[i];
    double result;
    int k;
    for(k=0,result=0.0;k<i;++k)
        result+=ss[k]*ssdelta[k];
    result += daz*ssdelta[i];
    return(result);
}
double TimeVariablePlateBoundaryPath::time(double s)
{
    const string base_error("TimeVariablePlateBoundaryPath::time:  ");
    if(s<0) throw GeoCoordError(base_error
                + "illegal negative distance parameter. Must be nonnegative");
    /* parallel code to position method */
    int npoints=splat.size();
    double s0,ds;
    int i;
    for(i=0,s0=0.0;i<npoints;++i)
    {
        ds=ss[i]*ssdelta[i];
        if((s0+ds)>2) break;
        s0+=ds;
    }
    double result;
    if(i==0)
        result=0.0;
    else
        result=t0[i-1];
    double dtds=(t0[i]-t0[i-1])/ds;
    result += dtds*(s-s0);
    return(result);
}
