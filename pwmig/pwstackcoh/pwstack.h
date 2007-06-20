#include "SeisppError.h"
#include "ensemble.h"
#include "mute.h"
#ifdef MATLABDEBUG
#include "MatlabProcessor.h"
#endif
using namespace std;
using namespace SEISPP;
#ifdef MATLABDEBUG
extern MatlabProcessor mp;
#endif

class DepthDependentAperture
{
public:
	int npoints;
	double *t;
	double *aperture;
	double *cutoff;
	DepthDependentAperture()
		{npoints=0;t=NULL;aperture=NULL;cutoff=NULL;};
	DepthDependentAperture(int npt){
		npoints=npt;
		t=new double[npt];
		aperture=new double[npt];
		cutoff=new double[npt];
	};
	DepthDependentAperture(const DepthDependentAperture& parent);
	~DepthDependentAperture(){
		if(t!=NULL)delete [] t; 
		if(aperture!=NULL)delete [] aperture;
		if(cutoff!=NULL)delete [] cutoff;
	};
	// pf constructor.  tag is keyword to hold list of values
	DepthDependentAperture(Pf *pf,string tag) throw(string);
	double get_aperture(double tnow);
	double get_cutoff(double tnow);
};

/* Special internal object for this program used to return
two different coherence estimate: (1) component by component
estimates and (2) vector sum coherence (scalar). 
*/
class Coharray
{
public:
	double t0;
	double dt;
	double ns;
	double winlen;  // width of averaging window in time units
	dmatrix compcoh;
	vector<double> coh;
	Coharray(double t0i,double dti, int nsi,double win)
	{
		t0=t0i;
		dt=dti;
		ns=nsi;
		winlen=win;
		coh.reserve(ns);
		compcoh=dmatrix(3,ns);
		compcoh.zero();
	};
};
void geographic_to_dne(double lat0, double lon0, double lat, double lon, 
		double *dnorth, double *deast);
void geographic_to_dne(double lat0, double lon0,
        double lat, double lon, double *dnorth, double *deast);
int compute_pseudostation_weights(int nsta, double *dnorth, double *deast,
                double aperture, double cutoff, double *w);
void compute_pwmoveout(int nsta,double *deast, double *dnorth,
                double ux, double uy, double *moveout);
int pwstack_ensemble(ThreeComponentEnsemble& indata,
	RectangularSlownessGrid& ugrid,
	TopMute& mute,
	TopMute& stackmute,
	double tstart,
	double tend,
	DepthDependentAperture& aperture,
	double dtcoh,
	double cohtwin,
	MetadataList& mdlcopy,
	MetadataList& mdlout,
	AttributeMap& am,
	string dir,
	DatabaseHandle& dbh);
string virtual_station_name(int ix1, int ix2);
Coharray compute_stack_coherence(vector<dmatrix> d,
	dmatrix& wgt,ThreeComponentSeismogram& stack,
	double dtcoh,double cohavgwin,TopMute& mute);
void save_coh(DatascopeHandle& dbhwf, DatascopeHandle& dbhcoh,
	Coharray& coh, ThreeComponentSeismogram& stack,
	int pwfid, MetadataList& mdl, AttributeMap& am);
//DEBUG
#ifdef MATLABDEBUG
extern MatlabProcessor mp;
#endif
