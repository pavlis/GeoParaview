#include <sstream>
#include "stock.h"
#include "arrays.h"
#include "pf.h"
#include "pwstack.h"
using namespace std;
//
// Pf-based constrructor for DepthDependentAperture object.
//
//  Arguments:
//	pf - input parameter file "object"
//	tag - name identifying the Tbl defining the designed input
//
//  This routine throws an exception with a message string if the required
// information from the pf cannot be retrieved.
// This should probably normally be fatal.  A more complex return
// to the error was deemed unnesseary since this is a specialized
// object for pwstack
//
//  Author:  Gary L. Pavlis
//

DepthDependentAperture::DepthDependentAperture(Pf *pf,string tag)
throw (string)
{
    Tbl *tlist;

    tlist = pfget_tbl(pf,(char *)tag.c_str());
    if(tlist==NULL) throw "Aperture definition missing from parameter file";
    npoints = maxtbl(tlist);
    t=new double[npoints];
    aperture=new double[npoints];
    cutoff=new double[npoints];

    for(int i=0;i<maxtbl(tlist);++i)
    {
        char *tline;
        tline = (char *)gettbl(tlist,i);
        stringstream instr(tline);
        instr >> t[i];
        instr >> aperture[i];
        instr >> cutoff[i];
    }
}


DepthDependentAperture::DepthDependentAperture(const DepthDependentAperture& a)
{
    npoints=a.npoints;
    if(npoints>0)
    {
        if(a.t == NULL)
        {
            t=NULL;
        }
        else
        {
            t=new double[npoints];
            dcopy(npoints,a.t,1,t,1);
        }
        if(a.aperture == NULL)
        {
            aperture=NULL;
        }
        else
        {
            aperture=new double[npoints];
            dcopy(npoints,a.aperture,1,aperture,1);
        }
        if(a.cutoff == NULL)
        {
            cutoff=NULL;
        }
        else
        {
            cutoff=new double[npoints];
            dcopy(npoints,a.cutoff,1,cutoff,1);
        }
    }
    else
    {
        t=NULL;
        aperture=NULL;
        cutoff=NULL;
    }
}


double linear_interpolate_irregular(int npts, double *x, double *y,
double xsearch)
{
    int i1,i2;
    double yout;
    //silently return endpoints if outside the range
    if (xsearch<x[0]) return(y[0]);
    if (xsearch>x[npts-1]) return(y[npts-1]);
    for(i2=1;i2<npts;++i2)
        if(x[i2]>xsearch) break;
    i1 = i2-1;
    yout = y[i1] + (xsearch-x[i1])*(y[i2]-y[i1])/(x[i2]-x[i1]);
    return(yout);
}


double DepthDependentAperture::get_aperture(double t0)
{
    return(linear_interpolate_irregular(npoints,t,aperture,t0));
}


double DepthDependentAperture::get_cutoff(double t0)
{
    return(linear_interpolate_irregular(npoints,t,cutoff,t0));
}
