#include "stock.h"
#include "pf.h"
#include "seispp.h"
//
// Loads data defined by a pf into a Depth_Dependent_Aperture object.
// This could be done as a constructor, but having a C object in 
// the pf as an arg to a constructor seems problematic even though
// I've done it elsewhere
//
void pfget_aperture_parameters(Pf *pf,string tag,
		Depth_Dependent_Aperture& d) 
				throw(pferror pferr)
{
	Tbl *t;

	t = pfget_tbl(pf,(char *)tag);
	if(t==NULL)
	{
		pferr.name=tag;
		pferr.pftype = PFTBL;
		pferr.message = "Aperture definition missing from parameter file";
		throw(pferr);
	}
	// Note this assumes d was created with the default constructor
	// that sets the entries null.  If the parameterized constuctor
	// were called this will create a memory leak.
	d.npoints = maxtbl(t);
	new double d.t=double(d.npoints);
	new double d.aperture=double(d.npoints);
	new double d.cutoff=double(d.npoints);

	for(int i=0;i<maxtbl(t);++i)
	{
		char *tline;
		tline = (char *)gettbl(t,i);
		istream instr(tline);
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

