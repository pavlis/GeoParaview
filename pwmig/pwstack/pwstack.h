#include "seispp.h"
using namespace SEISPP;

class Depth_Dependent_Aperture
{
public:
	int npoints;
	double *t;
	double *aperture;
	double *cutoff;
	Depth_Dependent_Aperture()
		{npoints=0;t=NULL;aperture=NULL;cutoff=NULL;};
	Depth_Dependent_Aperture(int npt){
		npoints=npt;
		t=new double[npt];
		aperture=new double[npt];
		cutoff=new double[npt];
	};
	~Depth_Dependent_Aperture(){
		if(t!=NULL)delete [] t; 
		if(aperture!=NULL)delete [] aperture;
		if(cutoff!=NULL)delete [] cutoff;
	};
	// pf constructor.  tag is keyword to hold list of values
	Depth_Dependent_Aperture(Pf *pf,string tag) throw(string);
	double get_aperture(double tnow);
	double get_cutoff(double tnow);
};

void geographic_to_dne(double lat0, double lon0, double lat, double lon, 
		double *dnorth, double *deast);
void geographic_to_dne(double lat0, double lon0,
        double lat, double lon, double *dnorth, double *deast);
int compute_pseudostation_weights(int nsta, double *dnorth, double *deast,
                double aperture, double cutoff, double *w);
void compute_pwmoveout(int nsta,double *deast, double *dnorth,
                double ux, double uy, double *moveout);
int pwstack_ensemble(Three_Component_Ensemble *indata,
	Rectangular_Slowness_Grid& ugrid,
	Top_Mute& mute,
	Top_Mute& stackmute,
	double lat0,
	double lon0,
	double ux0,
	double uy0,
	double tstart,
	double tend,
	Depth_Dependent_Aperture& aperture,
	Metadata_list& mdlcopy,
	Metadata_list& mdlout,
	Attribute_Map& am,
	Database_Handle& dbh);
