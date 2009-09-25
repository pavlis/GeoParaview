#include "coords.h"
#include "PPWDeconProcessor.h"
using namespace std;
using namespace SEISPP;
/* Added for debugging.  Remove when working */
double Linf(dmatrix& A) {
    double *aptr;
    aptr=A.get_address(0,0);
    double Ainf=*aptr;
    ++aptr;
    int na=A.rows() * A.columns();
    for(int i=1;i<na;++i,++aptr){
            if((*aptr)>Ainf) Ainf=(*aptr);
     }
     return(Ainf);
}
PWIndexPosition::PWIndexPosition()
{
    iu=0;
    lag0=0;
    time0=0.0;
    amp=0.0;
    for(int k=0;k<3;++k) v[k]=0.0;
}
PWIndexPosition::PWIndexPosition(int iuin, int lag0in, 
	double t0in, double ampin, double *vin)
{
	iu=iuin;
	lag0=lag0in;
        time0=t0in;
	amp=ampin;
	for(int k=0;k<3;++k)v[k]=vin[k];
}
PWIndexPosition::PWIndexPosition(const PWIndexPosition& parent)
{
	iu=parent.iu;
	lag0=parent.lag0;
	time0=parent.time0;
	amp=parent.amp;
	for(int k=0;k<3;++k)v[k]=parent.v[k];
}
PWIndexPosition& PWIndexPosition::operator=(const PWIndexPosition& parent)
{
    if(this != &parent)
    {
	iu=parent.iu;
	lag0=parent.lag0;
	time0=parent.time0;
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
		/* Note the wavelet is not scaled and assumed to 
		have a properly normalizd amplitude */
		wavelet=winp;
		lat=get_double("sta_lat");
		lon=get_double("sta_lon");
		elev=get_double("sta_elev");
                /* Assume lat,lon stored in Metadata always in degrees */
                lat=rad(lat);
                lon=rad(lon);
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

int PWDataMember::remove(int iu, double t0, double *a)
{
	/* Note this method intentionally does NOT check for bounds
	for efficiency.  User must guarantee this to avoid chaos */
	int rlag0=this->sample_number(t0)+lag[iu]-wavelet.sample_number(0.0);;	
	double *baseptr=u.get_address(0,rlag0);
	int i;
        double scale;
//DEBUG
dmatrix u0=u;
        for(i=0;i<3;++i)
        {
            scale=(-1.0)*a[i]*weight;
            daxpy(wavelet.ns,scale,&(wavelet.s[0]),1,baseptr,3);
            ++baseptr;
        }
//DEBUG
/*
for(i=0;i<u.columns();++i)
    cout << i <<" "
        << u0(0,i)<<" "<<u(0,i)<< " "
        << u0(1,i)<<" "<<u(1,i) << " "
        << u0(2,i)<<" "<<u(2,i)<<endl;
*/
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
PWStackMember::PWStackMember(const PWStackMember& parent) 
    : ThreeComponentSeismogram(parent)
{
	slow=parent.slow;
	/* These are the private attributes */
	ampmax=parent.ampmax;
	lagmax=parent.lagmax;
	amp=parent.amp;
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
	double *baseptr=u.get_address(0,0);
	/* Always recompute the amplitudes before searching for
	the maximum */
	for(int i=0;i<ns;++i)
	{
		/* Intentionally do not use push_back in 
		setting amp components for efficiency. */
		amp[i]=0.0;
		for(int k=0;k<3;++k,baseptr++)
				amp[i]+=(*(baseptr))*(*(baseptr));
		amp[i]=sqrt(amp[i]);
	}
//DEBUG
/*
if(slow.mag()<0.001) {
    cout << "Amplitudes for iu=220)"<<endl; 
    for(int foo=0;foo<ns;++foo) cout <<foo<<" "<< amp[foo]<<endl;
}
*/
	/* This duplicates the PeakAmplitude methods in seispp but
	differences in the algorithm require duplication of small
	code fragment */
	vector<double>::iterator amptr;
	amptr=max_element(amp.begin(),amp.end());
        /* Use the STL distance algorithm to get position of lagmax */
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
        int ndens=threecd.member.size();
	if(ndens<=0) throw SeisppError(string("PPWDeconProcessor::constructor:  ")
			+ "ensemble passed has no data.");
	dt_stack=parent_data.member[0].dt;
	aperture=md.get_double("pseudostation_array_aperture");
	cutoff=md.get_double("pseudostation_array_cutoff");
	stalat.reserve(ndens);
	stalon.reserve(ndens);
	vector<ThreeComponentSeismogram>::iterator dptr;
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
        int ndens=threecd.member.size();
	if(ndens<=0) throw SeisppError(string("PPWDeconProcessor::constructor:  ")
			+ "ensemble passed has no data.");
	dt_stack=parent_data.member[0].dt;
	aperture=md.get_double("pseudostation_array_aperture");
	cutoff=md.get_double("pseudostation_array_cutoff");
	stalat.reserve(ndens);
	stalon.reserve(ndens);
	vector<ThreeComponentSeismogram>::iterator dptr;
	stalat.reserve(ndens);
	stalon.reserve(ndens);
	wavelet.reserve(ndens);
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
    if(this!=&parent)
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
    return *this;
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
	for(int i=1;mptr!=current_stack.end();++mptr,++i)
	{
		double testval=mptr->maxamp();
//DEBUG
//cout << "In PPWDeconProcessor::maxamp:  maxamp for stack"<<i<<"="<<testval<<endl;
		if(testval>result.amp)
		{
			result.amp=testval;
			result.lag0=mptr->lag_at_max();
                        result.time0=mptr->time_at_max();
			result.iu=i;
		}
	}
	for(int k=0;k<3;++k) result.v[k]
	 		= current_stack[result.iu].u(k,result.lag0);
        return(result);
}
void PPWDeconProcessor::stack()
{
	vector<PWDataMember>::iterator dptr;
	vector<PWStackMember>::iterator sptr;
	vector<SlownessVector>::iterator uptr;
	vector<int>::iterator tlptr;
	double *drp,*srp;
	int ntosum,i,j,total_lag;
	for(i=0,uptr=slow.begin(),sptr=current_stack.begin();
		sptr!=current_stack.end();++i,++uptr,++sptr)
	{
		sptr->u.zero();
		srp=sptr->u.get_address(0,0);
		
		for(dptr=d.begin(),tlptr=start_offset.begin();
			dptr!=d.end();++dptr,++tlptr)
		{
			total_lag=dptr->lag[i]+(*tlptr);
                        srp=sptr->u.get_address(0,total_lag);
			drp=dptr->u.get_address(0,0);
                        ntosum=3*(dptr->ns);
//DEBUG
/*
if((total_lag+dptr->ns)>(sptr->ns) ){

  cerr << "Coding error:  "
    << "ntosum="<<ntosum<<" but stack ns ="<< sptr->ns<<endl;
  exit(-1);
}
*/
			for(j=0;j<ntosum;++j,srp++,drp++) *srp+=(*drp);
		}
//cout << "Stack member number "<<i<<" has peak amplitude="<<Linf(sptr->u)<<endl;
	}
}
/* Private method to appraise if solution has converged.  Intentionally 
made a method because it is not yet clear how to define convergence.
Stubbed initally as we'll just go on count for initial testing.*/
bool PPWDeconProcessor::converged(int current_count)
{
    if(current_count<maxiteration)
        return(false);
    else
        return(true);
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
                nd=d.size();
                /* We normalize each PWDataMember by wt/sumwt so we do not
                   need to normalized the stack repeatedly. */
                for(i=0;i<nd;++i){
                    double scale=(d[i].weight)/sumwt;
                    dscal(3*d[i].ns,scale,d[i].u.get_address(0,0),1);
                    d[i].weight=scale;  // need to store normalized weight too
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
		t0_stack=-(maxnt0+npad)*dt_stack;
//cout << "t0_stack="<<t0_stack<<endl;
		ns_stack=nint((maxtend-t0_stack)/dt_stack);
		ns_stack+=npad;
		current_stack.clear();
		current_stack.reserve(nu);
		ThreeComponentSeismogram stack_pattern(ns_stack);
		stack_pattern.live=true;
		stack_pattern.t0=t0_stack;
		stack_pattern.dt=dt_stack;
		stack_pattern.tref=relative;
		// This is not necessary because the stack method always 0s before stack
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
		int itercount(0);
		do {
			this->stack();
			PWIndexPosition ip(this->maxamp());
                            cout << "DEBUG  itercount="<<itercount<<endl
                                << "lag="<<ip.lag0
                                <<" or time="<<ip.time0
                                <<" iu="<<ip.iu
                                << " which is ux,uy="
                                << slow[ip.iu].ux <<" "<<slow[ip.iu].uy
                                << " vector values="
                                << ip.v[0]<<", "
                                << ip.v[1]<<", "
                                << ip.v[2]<<endl;
			for(int k=0;k<3;++k) 
                        {
			    result->member[ip.iu].u(k,ip.lag0)+=ip.v[k];
                        }
			for(i=0;i<nd;++i)
			{
				d[i].remove(ip.iu,ip.time0,ip.v);
			}
			++itercount;
		} while(!converged(itercount));
	} catch(...) {throw;};
}
