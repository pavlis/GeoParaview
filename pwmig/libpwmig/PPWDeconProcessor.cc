#include "coords.h"
#include "PPWDeconProcessor.h"
using namespace std;
using namespace SEISPP;
PWIndexPosition::PWIndexPosition(int iuin, int lag0in, 
	double ampin, double *vin)
{
	iu=iuin;
	lag0=lag0in;
	amp=ampin;
	for(int k=0;k<3;++k)v[k]=vin[k];
}
PWIndexPosition::PWIndexPosition(const PWIndexPosition& parent)
{
	iu=parent.iu;
	lag0=parent.lag0;
	amp=parent.amp;
	for(int k=0;k<3;++k)v[k]=parent.v[k];
}
PWIndexPosition& PWIndexPosition::operator=(const PWIndexPosition& parent)
{
    if(this != &parent)
    {
	iu=parent.iu;
	lag0=parent.lag0;
	amp=parent.amp;
	for(int k=0;k<3;++k)v[k]=parent.v[k];
    }
    return *this;
}

/*! Primary constructor for this object.

Construction of this beast fills in data attributes through 
a combination of initialization and computation.  Be aware raw
data (draw) is simply copied.  This means it is important that the
caller guarantee sufficient padding is present in draw to allow for
all feasible lags of any plane wave in the deconvolution or downstream
problems are likely. 
*/
PWDataMember::PWDataMember(ThreeComponentSeismogram& draw,
	TimeSeries& winp,
	    SlownessVector& u0in,
	        vector<SlownessVector>& rayp,
		    double wt, double lat0, double lon0)
		: ThreeComponentSeismogram(draw), u0(u0in)
{
	weight=wt;
	try {
		/* Scale the data vector by weight */
		dscal(3*ns,weight,u.get_address(0,0),1);
		/* Note the wavelet is not scaled and assumed to 
		have a properly normalizd amplitude */
		wavelet=winp;
		lat=get_double("sta_lat");
		lon=get_double("sta_lon");
		elev=get_double("sta_elev");
		/* Intentionally inline this.  There is a procedure
		in pwstack that does this, but this calculation is so
		simple it seems preferable here to just inline it */
		const double Rearth(6378.164);
		double dnorth,deast,delta,azimuth;
		dist(lat0,lon0,lat,lon,&delta, &azimuth);
		delta *= Rearth;
		deast=delta*sin(azimuth);
		dnorth=delta*cos(azimuth);
		/* Here we compute lag in samples or each slowness vector */
		int np=rayp.size();
		lag.reserve(np);
		int i;
		for(i=0;i<np;++i)
		{
			double ux=u0.ux+rayp[i].ux;
			double uy=u0.uy+rayp[i].uy;
			double tlag=ux*deast+uy*dnorth;
			lag.push_back(nint(tlag/dt));
		}
	}
	catch (MetadataGetError mderr)
	{
		throw mderr;
	}
}

PWDataMember::PWDataMember(const PWDataMember& parent)
  : ThreeComponentSeismogram(dynamic_cast<const ThreeComponentSeismogram&>(parent)),
    u0(parent.u0), lag(parent.lag),wavelet(parent.wavelet)
{
	weight=parent.weight;
	lat=parent.lat;
	lon=parent.lon;
	elev=parent.elev;
}
/* Uses trick described at:  
http://www.fredosaurus.com/notes-cpp/oop-overloading/overloadassign2.html
also in a more coherent/expanded form in 
http://www.acm.org/crossroads/xrds1-4/ovp.html
This essentially invokes the operator= function of the class 
from which this is derived.   In this case this is a ThreeComponentSeismogram.
*/
PWDataMember& PWDataMember::operator=(const PWDataMember& parent)
{
	if(this != &parent)
	{
		this->ThreeComponentSeismogram::operator=(parent);
		/* Now the actual members of this class */
		weight=parent.weight;
		lat=parent.lat;
		lon=parent.lon;
		elev=parent.elev;
		u0=parent.u0;
		wavelet=parent.wavelet;
		lag=parent.lag;
	}
	return *this;
}

int PWDataMember::remove(int iu, int lag0, double *a)
{
	/* Note this method intentionally does NOT check for bounds
	for efficiency.  User must guarantee this to avoid chaos */
	int rlag0=lag[iu]+lag0;	
	double *baseptr=u.get_address(rlag0,0);
	/* could do this with three calls to saxpy, but will inline this
	calculation instead until it proves to be a performance issue.  
	Assume the compiler can make this simple loop pretty efficient */
	int i;
	double scale[3];
	for(i=0;i<3;++i) scale[i]= (-1.0)*a[i]*weight;
	for(i=0;i<wavelet.ns;++i) *(baseptr+3*i) += scale[0]*wavelet[i];
	++baseptr;
	for(i=0;i<wavelet.ns;++i) *(baseptr+3*i) += scale[1]*wavelet[i];
	++baseptr;
	for(i=0;i<wavelet.ns;++i) *(baseptr+3*i) += scale[2]*wavelet[i];
}
PWStackMember::PWStackMember(ThreeComponentSeismogram& din,
	SlownessVector& sv) : ThreeComponentSeismogram(din)
{
	try{
		slow=sv;
		/*We initialize amp vector to 0s to inialize the container.
		Note ns is inherited from BasicTimeSeries*/
		amp.reserve(ns);
		for(int i=0;i<ns;++i)amp.push_back(0.0);
		ampmax=0.0;
		lagmax=0;
	} catch (...) {throw;};
}
/* Use same trick as PWDataMember to copy ThreeComponentSeismogram
attributes */
PWStackMember& PWStackMember::operator=(const PWStackMember& parent)
{
	if(this != &parent)
	{
		this->ThreeComponentSeismogram::operator=(parent);
		slow=parent.slow;
		/* These are the private attributes */
		ampmax=parent.ampmax;
		lagmax=parent.lagmax;
		amp=parent.amp;
	}
	return *this;
}
double PWStackMember::maxamp()
{
	/* Always recompute the amplitudes before searching for
	the maximum */
	for(int i=0;i<ns;++i)
	{
		double *baseptr=u.get_address(0,0);
		/* Intentionally do not use push_back in 
		setting amp components for efficiency. */
		amp[i]=0.0;
		for(int k=0;k<3;++i,baseptr+=3)
				amp[i]+=(*(baseptr))*(*(baseptr));
		amp[i]=sqrt(amp[i]);
	}
	/* This duplicates the PeakAmplitude methods in seispp but
	differences in the algorithm require duplication of small
	code fragment */
	vector<double>::iterator amptr;
	amptr=max_element(amp.begin(),amp.end());
	/* I don't think this is the right name for this algorithm.
	There is an stl algorithm to return the number of elements in
	a container between two iterators */
	lagmax=distance(amp.begin(),amptr);
	return(*amptr);
}

PPWDeconProcessor::PPWDeconProcessor(ThreeComponentEnsemble& threecd,
    SlownessVector& u0in,
	vector<SlownessVector>& ulist,
	    vector<TimeSeries>& w,
		Metadata& md) 
		: u0(u0in),parent_data(threecd), wavelet(w),slow(ulist)
{
    try {
    	maxiteration=md.get_int("maximum_iterations");
	nu=ulist.size();
	if(nd<=0) throw SeisppError(string("PPWDeconProcessor::constructor:  ")
			+ "ensemble passed has no data.");
	dt_stack=parent_data.member[0].dt;
	aperture=md.get_double("pseudostation_array_aperture");
	cutoff=md.get_double("pseudostation_array_cutoff");
	vector<ThreeComponentSeismogram>::iterator dptr;
	stalat.reserve(nd);
	stalon.reserve(nd);
	for(dptr=parent_data.member.begin();
		dptr!=parent_data.member.end();++dptr)
	{
		double lat,lon;
		lat=dptr->get_double("sta_lat");
		lon=dptr->get_double("sta_lon");
		stalat.push_back(rad(lat));
		stalon.push_back(rad(lon));
	}
	/* These are initialized to useless values.  They are set when
	the compute method is called */
	ns_stack=-1;
	t0_stack=0.0;
	start_offset.clear();
	stacksize=-1;
	current_stack.clear();
	d.clear();
	nd=-1;
    } catch (...) {throw;};
	
}
/* This is a near duplicate of the above.  There is, no doubt, a clever
way to utilize the above constructor as part of this one, but in my 
view it is better to suffer this duplication than it is to utilize 
a devious trick and make the code obscure.  HOWEVER, be sure in any 
revisions that errors found in the above propagate to here.  The only
difference here is we build the wavelet ensemble as copies of a single
master */
PPWDeconProcessor::PPWDeconProcessor(ThreeComponentEnsemble& threecd,
    SlownessVector& u0in,
	vector<SlownessVector>& ulist,
	    TimeSeries& w,
		Metadata& md) 
		: u0(u0in), parent_data(threecd), slow(ulist)
{
    try {
    	maxiteration=md.get_int("maximum_iterations");
	nu=ulist.size();
	if(nd<=0) throw SeisppError(string("PPWDeconProcessor::constructor:  ")
			+ "ensemble passed has no data.");
	dt_stack=parent_data.member[0].dt;
	aperture=md.get_double("pseudostation_array_aperture");
	cutoff=md.get_double("pseudostation_array_cutoff");
	vector<ThreeComponentSeismogram>::iterator dptr;
	int npd=parent_data.member.size();
	stalat.reserve(npd);
	stalon.reserve(npd);
	wavelet.reserve(npd);
	for(dptr=parent_data.member.begin();
		dptr!=parent_data.member.end();++dptr)
	{
		double lat,lon;
		lat=dptr->get_double("sta_lat");
		lon=dptr->get_double("sta_lon");
		stalat.push_back(rad(lat));
		stalon.push_back(rad(lon));
		/* This sets all wavelet members to a single waveform*/
		wavelet.push_back(w);
	}
	/* These are initialized to useless values.  They are set when
	the compute method is called */
	ns_stack=-1;
	t0_stack=0.0;
	start_offset.clear();
	stacksize=-1;
	current_stack.clear();
	d.clear();
	nd=-1;
    } catch (...) {throw;};
	
}
PPWDeconProcessor::PPWDeconProcessor(const PPWDeconProcessor& parent)
{
	u0=parent.u0;
	maxiteration=parent.maxiteration;
	nd=parent.nd;
	nu=parent.nu;
	ns_stack=parent.ns_stack;
	t0_stack=parent.t0_stack;
	dt_stack=parent.dt_stack;
	aperture=parent.aperture;
	cutoff=parent.cutoff;
	stalat=parent.stalat;
	stalon=parent.stalon;
	start_offset=parent.start_offset;
	stacksize=parent.stacksize;
	d=parent.d;
	parent_data=parent.parent_data;
	wavelet=parent.wavelet;
	slow=parent.slow;
	current_stack=parent.current_stack;
}
PPWDeconProcessor& PPWDeconProcessor::operator=(const PPWDeconProcessor& parent)
{
}
double PPWDeconProcessor::pseudostation_weight(double lat,double lon,
	double pslat, double pslon, double aperture, double cutoff)
{
	const double Rearth(6378.164);
	double dnorth,deast,delta,azimuth;
	dist(pslat,pslon,lat,lon, &delta, &azimuth);
	delta *= Rearth;
	if(delta>cutoff) 
		return(-1.0);
	else
		return(exp(-(delta*delta)/(2.0*aperture*aperture)));
}
PWIndexPosition PPWDeconProcessor::maxamp()
{
	vector<PWStackMember>::iterator mptr;
	PWIndexPosition result;
	mptr=current_stack.begin();
	result.amp=mptr->maxamp();
	result.lag0=mptr->lag_at_max();
	result.iu=0;
	++mptr;
	for(int i=0;mptr!=current_stack.end();++mptr,++i)
	{
		double testval=mptr->maxamp();
		if(testval>result.amp)
		{
			result.amp=testval;
			result.lag0=mptr->lag_at_max();
			result.iu=i;
		}
	}
	for(int k=0;k<3;++k) result.v[k]
	 		= current_stack[result.iu].u(k,result.lag0);
}
void PPWDeconProcessor::stack()
{
	vector<PWDataMember>::iterator dptr;
	vector<PWStackMember>::iterator sptr;
	vector<SlownessVector>::iterator uptr;
	vector<int>::iterator tlptr;
	/* raw pointers passed to saxpy obtained from get_address method */
	double *drp,*srp;
	int ntosum,i,j,total_lag;
	for(i=0,uptr=slow.begin(),sptr=current_stack.begin();
		sptr!=current_stack.end();++i,++uptr,++sptr)
	{
		sptr->u.zero();
		drp=dptr->u.get_address(0,0);
		ntosum=3*(dptr->ns);
		
		for(dptr=d.begin(),tlptr=start_offset.begin();
			dptr!=d.end();++dptr,tlptr)
		{
			total_lag=dptr->lag[i]+(*tlptr);
			srp=sptr->u.get_address(0,total_lag);
			for(j=0;j<ntosum;++j,srp++,drp++) *srp+=(*drp);
		}
	}
}
/* Private method to appraise if solution has converged.  Intentionally 
made a method because it is not yet clear how to define convergence.
Stubbed initally as we'll just go on count for initial testing.*/
bool PPWDeconProcessor::converged()
{
}
auto_ptr<ThreeComponentEnsemble> PPWDeconProcessor::compute(double pslat,
					double pslon)
{
	try {
		int i,j;
		/* First we have to form the pseudostation gather */
		d.clear();
		double sumwt;
		int npd=parent_data.member.size();
		for(i=0;i<npd;++i)
		{
			double wt;
			double lat,lon;
			wt=pseudostation_weight(stalat[i],stalon[i],pslat,pslon,
				aperture, cutoff);
			if(wt>0.0) 
			{
			  d.push_back(PWDataMember(parent_data.member[i],
			  	wavelet[i],u0,slow,wt,pslat,pslon));
			  sumwt += wt;
			}
		}
		/* This section does two things.  First, it computes a size
		that allows blind sums without error checking.  Second, it initializes
		the seismic data in each element of the container to all zeros. */
		int maxnt0;
		double maxtend;
		vector<int>::iterator lptr;
		vector<PWDataMember>::iterator dptr;
		for(i=0,dptr=d.begin();dptr!=d.end();++i,++dptr)
		{
		    int nt0tmp=dptr->sample_number(0.0);
		    double tendtmp=dptr->endtime();
		    for(j=0,lptr=dptr->lag.begin();
		    	lptr!=dptr->lag.end();++j,++lptr)
		    {
		    	int nt0test=nt0tmp - (*lptr);
			double tendtest=tendtmp+((*lptr)*dt_stack);
			if( (i==0) && (j==0) )
			{
			    maxnt0=nt0test;
			    maxtend=tendtmp;
			}
			else
			{
			    if(maxnt0<nt0test) 
			    	maxnt0=nt0test;
			    if(maxtend<tendtest)
			        maxtend=tendtest;
			}
		    }
		}
		
		/* We add npad on top and bottom to be sure we don't overflow 
		the bounds with computed stack indices.  Do NOT make this too 
		large as it adds an inefficieny. */
		const int npad(2);  
		t0_stack=(maxnt0+npad)*dt_stack;
		ns_stack=nint((maxtend-t0_stack)/dt_stack);
		ns_stack+=npad;
		current_stack.clear();
		current_stack.reserve(nu);
		ThreeComponentSeismogram stack_pattern(ns_stack);
		stack_pattern.live=true;
		stack_pattern.t0=t0_stack;
		stack_pattern.dt=dt_stack;
		stack_pattern.tref=relative;
		// This is ot necessary because the stack method always 0s before stack
		//stack_pattern.zero();
		vector<SlownessVector>::iterator uptr;
		for(uptr=slow.begin();uptr!=slow.end();++uptr)
		{
			current_stack.push_back(PWStackMember(stack_pattern,
				*uptr));
		}
		/* Compute lags in a parallel vector for efficiency */
		start_offset.clear();
		start_offset.reserve(nd);
		for(dptr=d.begin();dptr!=d.end();++dptr)
		{
			int thislag;
			thislag=current_stack[0].sample_number(dptr->t0);
			start_offset.push_back(thislag);
		}

		/* NOTE DURING DEVELOPMENT:  may want to compute coherence
		estimaes.  How to most efficiently do this is less than clear. 
		Probably need a way to access residuals after the iterative
		method finishes and compare it to original data.  This may 
		require caching either raw data with nonzero weights or initial 		stacks.  */
		
		auto_ptr<ThreeComponentEnsemble> result
			= auto_ptr<ThreeComponentEnsemble> (new 
				ThreeComponentEnsemble(nu,ns_stack));
		/* Prep this gather */
		vector<ThreeComponentSeismogram>::iterator mptr;
		for(mptr=result->member.begin(),i=0;
		    mptr!=result->member.end();++mptr,++i)
		{
			/* These need to be set. Defaults are ok for other
			attributes of this object*/
			mptr->ns=ns_stack;
			mptr->t0=t0_stack;
			mptr->dt=dt_stack;
			mptr->tref=relative;
			mptr->u.zero();
			/* These are needed in output for pwmig.  Posting these
			in this way is needed for interface to save results
			with the PwmigFileHandle object. */
			mptr->put("ux0",u0.ux);
			mptr->put("uy0",u0.uy);
			mptr->put("ux",slow[i].ux);
			mptr->put("uy",slow[i].uy);
		}
		/* Main processing loop */
		int itercount;
		bool ctest;
		do {
			this->stack();
			PWIndexPosition ip(this->maxamp());
			for(int k=0;k<3;++k) 
			    result->member[ip.iu].u(k,ip.lag0)+=ip.v[k];
			for(i=0;i<nd;++i)
			{
				d[i].remove(ip.iu,ip.lag0,ip.v);
			}
			++itercount;
			ctest=this->converged();
		} while((itercount<maxiteration) && !ctest);
	} catch(...) {throw;};
}
