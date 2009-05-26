#ifndef _PPWSDECONPROCESSOR_H_
#define _PPWSDECONPROCESSOR_H_
#include <vector>
#include "slowness.h"
#include "TimeSeries.h"
#include "seispp.h"

namespace SEISPP {
using namespace std;
using namespace SEISPP;
/*! \brief Data structure for PPWDeconProcessor.

This is a data structure that is an extension of a ThreeComponentSeismogram. 
The main use is efficiency and to simplify rebuilding a weighted gather for 
a new pseudostation. */
class PWDataMember: public ThreeComponentSeismogram
{
public:
	/*! Weighted used to scale original data to produce this object. */
	double weight;
	/*! Location of receiver (extracted from inherited Metadata ) */
	double lat,lon,elev;
	/*! Incident wavefield slowness vector */
	SlownessVector u0;  
	TimeSeries wavelet;
	/* These hold parallel vector of lags for this this plane wave
	stack for each slowness vector */
	vector<int> lag;
	PWDataMember(ThreeComponentSeismogram& draw,
	    TimeSeries& w,
	        SlownessVector& u0in,
			vector<SlownessVector>& rayp,
				double weight,
					double lat0,
						double lon0);
	/*! Standard copy constructor */
	PWDataMember(const PWDataMember& parent);
	/*! Standard assignment operator */
	PWDataMember& operator=(const PWDataMember& parent);
	/*! \brief Primary method.  
	Subtracts a*wavelet shifted b
	lag0 + lag[iu].  Be warned this method does it's work
	with NO error checking.  It assumes a is length 3 and
	more importantly that lag0+lag[iu] will not exceed
	internal bounds.*/
	int remove(int iu, int lag0,double *a);
};
class PWStackMember : public ThreeComponentSeismogram
{
public:
	SlownessVector slow;
	PWStackMember(ThreeComponentSeismogram& din,SlownessVector& sv);
	PWStackMember(const PWStackMember& parent);
	PWStackMember& operator=(const PWStackMember& parent);

	/*! Return the maximum amplitude of this member.

	This is a core method for the iterative method.  Returns
	the maximum amplitude of this 3d (stack) eismogram.  This ALWAYS 
	initiates a recalculation of the amplitudes from the current
	data.  Calls to this method must ALWAY precede calls to 
	the lag_at_max method or the return of that method will 
	be wrong.  The reason is that this method initates the 
	recalculation of the amplitude while lag_at_max simply 
	returns a value. */
	double maxamp();
	int lag_at_max(){return lagmax;};
private:
	double ampmax;
	int lagmax;
	/* This holds an array of 3-c amplitude values.  Update on 
	each call to maxamp() method. */
	vector<double> amp;
};
class PWIndexPosition
{
public:
	PWIndexPosition();
	PWIndexPosition(int iuin, int lag0in, double ampin, double *vin);
	PWIndexPosition(const PWIndexPosition& parent);
	PWIndexPosition& operator=(const PWIndexPosition& parent);
	int iu;
	int lag0;
	double amp;
	double v[3];
};
class PPWDeconProcessor
{
public:
	/*! Slowness vector of incident wavefield */
	SlownessVector u0;
	PPWDeconProcessor(ThreeComponentEnsemble& threecd,
		SlownessVector& u0in,
			vector<SlownessVector>& ulist, 
				TimeSeries& w,
					Metadata& md);
	/*! Constructor for common region gathers.

	Expect that it will be useful to run this on composites from
	common source area (program telecluster).  In this case we take
	a vector of ensembles, stack them without moveouts, stack the
	actual outputs, and result looks just like a single event gather.
	BIG complication of this is that it makes the actual output dependent
	on each station.  As a result the PWmember function has an actual
	output member.  (see above).  For single events  these are all 
	the same. Need to decide how the input gather should be organized.  
	May just need to require evid and sta to defined in trace metadata*/
	PPWDeconProcessor(ThreeComponentEnsemble& threecd,
		SlownessVector& u0in,
			vector<SlownessVector>& ulist,
				vector<TimeSeries>& w,
					Metadata& md);
	/*! Standard copy constructor. */
	PPWDeconProcessor(const PPWDeconProcessor& parent);
	/*! Standard assignment operator.  */
	PPWDeconProcessor& operator=(const PPWDeconProcessor& parent);
				
	/*! Compute method.

	NOTES FOR INITIAL:  convert to documentation when finished.
	
	This method does all the work. It shuld compute a vector of 
	3c seismograms stored in the ensemble.  Each member is one 
	plane wave component for the pseudostation at plat,plon.  
	Note to provide a clean interface into the PwmigFileHandle 
	each member needs each of the following metadata elements set:
	gridid, ix1, ix2, ux0, uy0, ux, uy, elev.  Further ns, dt, and t0
	will be extracted from the BasicTimeSeries attributes.  The
	use of the auto_ptr is a bit painful, but necessary for efficiency.*/
	auto_ptr<ThreeComponentEnsemble> compute(double pslat, double pslon);
private:
	int maxiteration;  /* Limit on number of iterative cycles. */
	int nd; /* d.size() cached for efficiency */
	int nu;  /* slow.size() cached for efficiency */
	int ns_stack;  /* stack ns */
	double t0_stack;  /* Relative time t0 for plane wave stacks */
	double dt_stack;  /* sample interval of stack data */
	double aperture, cutoff;  /* pseudostation sigma and cutoff */
	/* Cache these for efficiency */
	vector<double> stalat,stalon;
	/* This is an nu length vector parallel to slow containing offset
	value for an internal buffer used to accumulate plane wave 
	components.  */
	vector<int> start_offset;
	/* This is used to determine stack buffer size.  It is set
	by constructor as a conservative estimate to allow any 
	feasible plane wave shifts (scaled by maximum aperture and
	range of slowness vectors.  Note because the stack is 3c the
	total stack size is 3*stacksize, but held in a dmatrix object
	created and then destroyed  when the stack method is called*/
	int stacksize;
	/* This object uses a vector of slowness values to be more general
	than the rectangular slowness grid used in pwmig.  This is
	intentional to make this more general. */
	vector<SlownessVector> slow;
	vector<PWDataMember> d;
	vector<TimeSeries> wavelet;
	ThreeComponentEnsemble parent_data;
	/* Plane wave stack components (length nu) */
	vector<PWStackMember> current_stack;
	/*! workhorse internal method to reform stack in place.  i.e
	not copying is done, just components of current_stack are updated*/
	void stack();
	PWIndexPosition maxamp();
	double pseudostation_weight(double lat, double lon,
		double pslat, double pslon, 
			double aperture, double cutoff);
	bool converged();
};
}  /* End SEISPP Namespace encapsulation */
#endif
