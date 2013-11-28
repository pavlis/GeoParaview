#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include "seispp.h"
#include "dbpp.h"
#include "gclgrid.h"
#include "dmatrix.h"
using namespace std;
using namespace SEISPP;
/* This routine extracts a depth window from the grid using interpolation
   at interval dz over a range windowlen*dz.  This is not a very 
   efficient way to do this since this will be called many times
   and most of the calculations done here will be the same for 
   all points when using a constant depth slice.  However, I 
   did this because I plan to recycle this code for dipping 
   interfaces later and then this is the necessary way to do this. */
dmatrix ExtractDataWindow(GCLvectorfield3d& g, double lon, double lat, 
        double depth, int windowlen,double dz)
{
    try {
        dmatrix result(3,windowlen);
        result.zero();
        Cartesian_point cp;
        /* compute depth range */
        double zmin,zmax,zhalf; 
        zhalf=static_cast<double>(windowlen-1)*dz/2.0;
        zmin=depth-zhalf;
        /* Silently shift if window hits ceiling */
        if(zmin<0.0) zmin=0.0;
        zmax=zmin+2.0*zhalf;
        /* Assume g is a regular grid with base at constant radius */
        double zfloor=g.depth(0,0,0);
        if(zmax>zfloor)
        {
            if( (zfloor-(2.0*zhalf))<0.0)
            {
                cerr << "Illegal windowlen for dz of this grid"<<endl
                    <<"Window length "<<2.0*zhalf<<" exceed depth range of grid"
                    <<endl;
                exit(-1);
            }
            else
            {
                zmax=zfloor;
                zmin=zfloor-2.0*zhalf;
            }
        }
        int index[3];
        /* Now fill the result matrix filling the matrix from the bottom
           up.  Be warned if you try to plot this it will be a bit 
           upside down compared to standard convention.  Done for speed
           because of a peculiarity of the lookup method. */
        double z;
        Geographic_point gp;
        gp.lon=lon;
        gp.lat=lat;
        gp.r=r0_ellipse(lat)-zmax;
        for(int i=0;i<windowlen;++i)
        {
            cp=g.gtoc(gp);
            int ierr=g.lookup(cp.x1,cp.x2,cp.x3);
            if(ierr!=0) 
            {
                /* Note this returns the count where lookup failed */
                throw i;
            }
            double *vec;
            vec=g.interpolate(cp.x1,cp.x2,cp.x3);
            /*gridstacker sets component 4 to stack count.  If zero 
              it means there is no data for that cell.  Since this 
              is a common situation it cannot be treated as an 
              exception.  Thus we instead return a 1x1 matrix to 
              signal this condition. */
            if(vec[4]<=0.0) 
            {
                delete [] vec;
                return (dmatrix(1,1));
            }
            for(int k=0;k<3;++k) result(k,i)=vec[k];
            delete [] vec;
            gp.r -= dz;
        }
        return(result);
    } catch(...){throw;};
}
double rms(vector<double> d)
{
    vector<double>::iterator dptr;
    double result;
    for(result=0.0,dptr=d.begin();dptr!=d.end();++dptr)
        result+=(*dptr)*(*dptr);
    result=sqrt(result)/static_cast<double>(d.size());
    return(result);
}
double mad(vector<double> d)
{
    vector<double>::iterator dptr;
    vector<double> amp;
    amp.reserve(d.size());
    for(dptr=d.begin();dptr!=d.end();++dptr)
        amp.push_back(fabs(*dptr));
    double result=median<double>(amp);
    return(result);
}
/* These do computations.  All intentionally do not 
   test comp for validity or include an error handler.
   Done purely for efficiency.
   Change if these are transported outside this code.  
   They are also painfully inefficient and repetitious.
   For efficiency it would be preferable to have 
   one routine to extract a component and another to
   compute 3c amplitudes then return all the stats in
   a single data structure.  For this application this
   was considered unnecessary as performance is not
   an issue for current plans.*/
double compute_rms(dmatrix& d,int comp)
{
    vector<double> dv;
    int n=d.columns();
    for(int i=0;i<n;++i) dv.push_back(d(comp,i));
    return(rms(dv));
}
double compute_mad(dmatrix& d,int comp)
{
    vector<double> dv;
    int n=d.columns();
    for(int i=0;i<n;++i) dv.push_back(d(comp,i));
    return(mad(dv));
}
double compute_max(dmatrix& d,int comp)
{
    vector<double> dv;
    vector<double>::iterator maxvalptr;
    int n=d.columns();
    for(int i=0;i<n;++i) dv.push_back(d(comp,i));
    maxvalptr=max_element(dv.begin(),dv.end());
    return(*maxvalptr);
}
double compute_rms3c(dmatrix& d)
{
    vector<double> dv;
    int n=d.columns();
    for(int i=0;i<n;++i) 
    {
        double amp;
        amp=0.0;
        for(int k=0;k<3;++k) 
            amp+=d(k,i)*d(k,i);
        amp=sqrt(amp);
        dv.push_back(amp);
    }
    return(rms(dv));
}
double compute_mad3c(dmatrix& d)
{
    vector<double> dv;
    int n=d.columns();
    for(int i=0;i<n;++i) 
    {
        double amp;
        amp=0.0;
        for(int k=0;k<3;++k) 
            amp+=d(k,i)*d(k,i);
        amp=sqrt(amp);
        dv.push_back(amp);
    }
    return(mad(dv));
}
double compute_max3c(dmatrix& d)
{
    vector<double> dv;
    vector<double>::iterator maxvalptr;
    int n=d.columns();
    for(int i=0;i<n;++i) 
    {
        double amp;
        amp=0.0;
        for(int k=0;k<3;++k) 
            amp+=d(k,i)*d(k,i);
        amp=sqrt(amp);
        dv.push_back(amp);
    }
    maxvalptr=max_element(dv.begin(),dv.end());
    return(*maxvalptr);
}
void usage()
{
    cerr << "horizonstats db gridname fieldname windowsize(km) < inlist"<<endl
        << "gridname:fieldname expected in db.gclfield table"<<endl
        << "Compute with window length windowlen in depth sampled at grid dz/2"<<endl;;
    exit(-1);
}
bool SEISPP::SEISPP_verbose(true);
int main(int argc, char **argv)
{
    if(argc<5) usage();
    ios::sync_with_stdio();
    string dbname(argv[1]);
    string gridname(argv[2]);
    string fieldname(argv[3]);
    int windowsize=atof(argv[4]);
    if(windowsize<0)
    {
        cerr << "Illegal window size = "<<windowsize<<endl;
        usage();
    }
    try{
        DatascopeHandle dbh(dbname,true);
        /* pwmig output is a 5component vector field.  Read it from
           the database */
        GCLvectorfield3d image(dbh,gridname,fieldname,5);
        const double oversampling(2.0);  // oversampling relative to dz
        double dz=image.dx3_nom/oversampling;
        int windowlen=SEISPP::nint(windowsize/dz)+1;
        double plon, plat, pdepth;
        char line[256];
        while(cin.getline(line,256))
        {
            sscanf(line,"%lf%lf%lf",&plon,&plat,&pdepth);
            plon=rad(plon);
            plat=rad(plat);
            try {
                dmatrix d=ExtractDataWindow(image,plon,plat,pdepth,
                        windowlen,dz);
                /* A window with null data in the window has d set
                   as a 1x1 matrix.  Thus this is a test for null
                   data.  We only write results without valid data*/
                if(d.rows()>1)
                {
                    double rmsT=compute_rms(d,0);
                    double madT=compute_mad(d,0);
                    double maxT=compute_max(d,0);
                    maxT=fabs(maxT);
                    double rmsR=compute_rms(d,1);
                    double madR=compute_mad(d,1);
                    double maxR=compute_max(d,1);
                    maxR=fabs(maxR);
                    double rms3d=compute_rms3c(d);
                    double mad3d=compute_mad3c(d);
                    double max3d=compute_max3c(d);
                    cout << deg(plon)<<" "
                        <<deg(plat)<<" "
                        <<pdepth<<" "
                        <<rmsR<<" "
                        <<madR<<" "
                        <<maxR<<" "
                        <<rmsT<<" "
                        <<madT<<" "
                        <<maxT<<" "
                        <<rms3d<<" "
                        <<mad3d<<" "
                        <<max3d<<endl;
                }
            } catch(int ierr)
            {
                cerr << "GCLgrid3d lookup method failed for point "
                    <<line<<endl
                    <<"Lookup failed to find point "<<ierr
                    <<" in window of length "<<windowlen<<endl
                    <<"Point skipped - trying next input point"<<endl;
            };
        }
    } catch (SeisppError& serr)
    {
        serr.log_error();
    }
    catch (GCLgridError& gclerr)
    {
        cerr << gclerr.what()<<endl;
    }
}
