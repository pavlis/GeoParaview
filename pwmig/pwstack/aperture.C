#include <strstream>
#include "stock.h"
#include "arrays.h"
#include "pf.h"
#include "seispp.h"
#include "pwstack.h"
using namespace std;
//
// Loads data defined by a pf into a Depth_Dependent_Aperture object.
// This could be done as a constructor, but having a C object in 
// the pf as an arg to a constructor seems problematic even though
// I've done it elsewhere
//  Arguments:
//	pf - input parameter file "object"
//	tag - name identifying the Tbl defining the designed input
//      d - returned object (assumed created with default constructor.)
//
//  This routine throws an exception with a message string if the required
// information from the pf cannot be retrieved.  Caller should arrane
// for a catch( string s) 
void pfget_aperture_parameters(Pf *pf,string tag,
		Depth_Dependent_Aperture& d) 
{
	Tbl *t;

	t = pfget_tbl(pf,(char *)tag.c_str());
	if(t==NULL) throw "Aperture definition missing from parameter file";
	// Note this assumes d was created with the default constructor
	// that sets the entries null.  To avoid a memory leak we clear these
	// vectors if the entries aren't NULL
	if(d.t != NULL) delete [] d.t;
	if(d.aperture != NULL) delete [] d.aperture;
	if(d.cutoff != NULL) delete [] d.cutoff;
	d.npoints = maxtbl(t);
	d.t=new double[d.npoints];
	d.aperture=new double[d.npoints];
	d.cutoff=new double[d.npoints];

	for(int i=0;i<maxtbl(t);++i)
	{
		char *tline;
		tline = (char *)gettbl(t,i);
		istrstream instr(tline);
		instr >> d.t[i];
		instr >> d.aperture[i];
		instr >> d.cutoff[i];
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
double Depth_Dependent_Aperture::get_aperture(double t0)
{
	return(linear_interpolate_irregular(npoints,t,aperture,t0));
}
double Depth_Dependent_Aperture::get_cutoff(double t0)
{
	return(linear_interpolate_irregular(npoints,t,cutoff,t0));
}

