#include "seispp.h"
class Rectangular_Slowness_Grid
{
public:
	string name;
	double uxlow, uylow;
	double dux, duy;
	int nux, nuy;
	Rectangular_Slowness_Grid(pf *);
	double ux(int i) {return(uxlow+i*dux);};
	double uy(int i) {return(uylow+i*duy);};
}

class Depth_Dependent_Aperture
{
public:
	int npoints;
	double *t;
	double *aperture;
	double *cutoff;
	Depth_Dependent_Aperture()
		{npoints=0;t=NULL;aperture=NULL;cutoff=NULL};
	Depth_Dependent_Aperture(int npt){npoints=npts;
		new double *t=double[npts];
		new double *aperture=double[npts];
		new double *cutoff=double[npts];
	};
	~Depth_Dependent_Aperture(){
		if(t!=NULL)delete [] t; 
		if(aperture!=NULL)delete [] aperture;
		if(cutoff!=NULL)delete [] cutoff;
	};
	double get_aperture(double tnow);
	double get_cutoff(double tnow);
}
class Three_Component_Stack : public Three_Component_Seismogram
{
public:
	int end_top_mute_zone;
	int last_valid_sample;
	int nstack;
	Three_Component_Stack();
	Three_Component_Stack(int nsamp);  // needs a parameterized constructor
}

