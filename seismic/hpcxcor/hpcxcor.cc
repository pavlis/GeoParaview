#include "seispp.h"
#include "Hypocenter.h"
#include "SignalToNoise.h"
#include "ensemble.h"
#include "XcorProcessingEngine.h"
#ifdef MPI_SET
        #include <mpi.h>
#endif
using namespace std;
using namespace SEISPP;
bool SEISPP::SEISPP_verbose(false);
/* This program is a variant of dbxcor and rtxcor.  It is aimed as a batch
processing program to grind through large volumes of 
data without intervention.  This is in contrast to dbxcor which is aimed
at interactive use and rtxcor aimed for being launched in a real time
environment.  

The algorithm is pretty simple.  We construct a XcorProcessingEngine
and feed it data as a series of gathers.  The feeding depends on the 
engine processing mode.  

I intend this program to be workable in an hpc environment allowing
it to run on a many processors simultaneously.  It does this by 
using a procesing queue to drive the processing.  Output is a design
issue that will evolve.  I'm writing this right now at the design 
phase to collect my ideas.

 List of methods used to select one or more traces as reference event
for MultichannelCorrelator */
enum ReferenceTraceMethod {UseAll, MaxSnrDB, MaxSnrComputed, MaxPeakAmplitude};
/* file scope keyword used here */
string dbsnrkey("arrival.snr");
/* This is a helper used by some of the reference trace methods below.
Returns a integer of first live trace.  -1 means no live data were found.
Note this could easily have been templatized for any ensemble type in SEISPP. */
int find_first_live(TimeSeriesEnsemble& d)
{
	int i,ival;
	int ensemblesize=d.member.size();
	for(i=0,ival=-1;i<ensemblesize;++i)
	{
		if(d.member[i].live)
		{
			ival=i;
			break;
		}
	}
	return(ival);
}
/* These are methods tied to ReferenceTraceMethod.  All return list of int
values that are an index position into ensemble member vector.  
Some may return an empty list if the algorithm fails.  Caller needs to handle
this error condition and deal with it.
*/
/* This method returns a trivial list of all member (i=0, ..., d.member.size-1)
or a subset of size MaxToUse, whichever is smaller.  If not really all, 
the list is approximated by skippling through the list at intervals of mod size/MaxToUse 
*/
list<int> UseAll_ChooseReference(TimeSeriesEnsemble& d,int MaxToUse)
{
	if(MaxToUse<1) throw SeisppError(string("UseAll_ChooseReference:  ")
			+"maximum to use argument must be a positive integer");
	list<int> result;
	int ensemblesize=d.member.size();
	int i;
	if(ensemblesize<MaxToUse)
	{
		for(i=0;i<ensemblesize;++i) result.push_back(i);
	}
	else
	{
		int skipinterval=ensemblesize/MaxToUse;
		if(skipinterval<1)skipinterval=1;
		for(i=0;i<MaxToUse;i+=skipinterval)
			result.push_back(i);
	}
	return(result);
}
/* Chooses largest value of posted snr through trace header (metadata)
keyword key.  */
list<int> MaxSnrDB_ChooseReference(TimeSeriesEnsemble& d, string key)
{
	list<int> result;
	vector<double> snrvalues;
	int ensemblesize=d.member.size();
	snrvalues.reserve(ensemblesize);
	int number_fetched=0;
	int i;
	for(i=0;i<ensemblesize;++i)
	{
		double val;
		try {
			val=d.member[i].get_double(key);
			++number_fetched;  
		} catch (MetadataGetError mdge)	
		{
			if(SEISPP_verbose) 
			{
				cerr << "Error MaxSnrDB_ChooseReference: "<<endl;
				mdge.log_error();
			}
			val=-1.0;
		}
		/* We must copy this here to create chaos in the sort because the
		sort method doesn't know this key */
		d.member[i].put(snr_keyword,val);
		snrvalues.push_back(val);
	}
	if(number_fetched<=0)
	{
		// This maybe should be only used if verbose is set, but this seems wise to 
		// me to always log this.
		int ival=find_first_live(d);
		if(ival>=0)
		{
			cerr << "Warning(MaxSNRDB_ChooseReference:  "
			  << "No members of this ensemble had snr posted"
			  <<endl
			  <<"Returning first live data (member="<<ival<<") as default"<<endl;
			result.push_back(ival);
		}
	}
	else
	{
		/* may be a more elegant way to do this than a linear search, but the
		cost here is small */
		int imax=0;
		double snrmax=snrvalues[0];
		for(i=1;i<ensemblesize;++i)
		{
			if(snrvalues[i]>snrmax)
			{
				snrmax=snrvalues[i];
				imax=i;
			}
		}
		result.push_back(imax);
	}
	return(result);
}
/* Returns trace with computed highs snr based in input signal and noise
time windows */
list<int> MaxSNRComputed_ChooseReference(TimeSeriesEnsemble& d, TimeWindow signal, TimeWindow noise)
{
	list<int> result;
	int imax;
	double snrmax;
	/* This template does 99% of the work.  All we do afterward is scan the metadata
	and find the largest value */
	ensemble_SNR_rms<TimeSeriesEnsemble,TimeSeries>(d,signal,noise,snr_keyword);
	int i;
	snrmax=-1.0;
	for(i=find_first_live(d);i<d.member.size();++i)
	{
		double testval;
		if(d.member[i].live)
		{
		 	testval=d.member[i].get_double(snr_keyword);
			if(testval>snrmax)
			{
				snrmax=testval;
				imax=i;
			}
		}
	}
			
	/* Take same philosophy as above if all failed. */
	if(snrmax<0.0)
	{
		int ival=find_first_live(d);
		if(ival>=0)
		{
			cerr << "Warning(MaxSNRComputed_ChooseReference:  "
			  << "SNR compution failed for all "
			  << d.member.size() << " members of this ensemble"
			  <<endl
			  <<"Returning first live data (member="<<ival<<") as default"<<endl;
			result.push_back(ival);
		}
		// Caller must handle this condition of an empty list and should also
		// issue a more specific warning.  In verbose mode this will produce
		// two similar messages.
		else if (SEISPP_verbose)
			cerr << "Warning(MaxSNRComputed_ChooseReference:  "
				<< " no live data in this ensemble"<<endl;
	}
	else
	{
		result.push_back(imax);
	}
	return(result);
}

/* This method uses a peak largest amplitude method to get the largest amplitude 
signal.  The entire signal is scanned because we make use of the template function
in libseispp.  The flow of the algorithm is virtually identical to the rms code
above and was, in fact, produced by editing it to produce this. */
const string ampltude_keyword("peak_amplitude");
list<int> MaxPeakAmplitude_ChooseReference(TimeSeriesEnsemble& d)
{
	list<int> result;
	int imax;
	double amax;
	int i;
	/* amplitude_static_keyword is one of a fixed set of keywords allowed by the
	less than totally generic sort functionals defined in XcorProcessingEngine.h.
	We use it somewhat incorrectly here as a convenience, so beware.*/
	MeasureEnsemblePeakAmplitudes<TimeSeriesEnsemble,TimeSeries>(d,amplitude_static_keyword);
	amax=-1.0;
	for(i=find_first_live(d);i<d.member.size();++i)
	{
		double testval;
		if(d.member[i].live)
		{
		 	testval=d.member[i].get_double(amplitude_static_keyword);
			if(testval>amax)
			{
				amax=testval;
				imax=i;
			}
		}
	}
	/* Take same philosophy as above if all failed. */
	if(amax<0.0)
	{
		int ival=find_first_live(d);
		if(ival>=0)
		{
			cerr << "Warning(MaxPeakAmplitude_ChooseReference:  "
			  << "Peak amplitude compution failed for all "
			  << d.member.size() << " members of this ensemble"
			  <<endl
			  <<"Returning first live data (member="<<ival<<") as default"<<endl;
			result.push_back(ival);
		}
		// Caller must handle this condition of an empty list and should also
		// issue a more specific warning.  In verbose mode this will produce
		// two similar messages.
		else if (SEISPP_verbose)
			cerr << "Warning(MaxPeakAmplitude_ChooseReference:  "
				<< " no live data in this ensemble"<<endl;
	}
	else
	{
		result.push_back(imax);
	}
	return(result);
}
/* This function can be used by any method that posts it's result as a double 
accessible with keyword key from each member.  It returns a list of int 
index values of the original ensemble containing the largest values of 
entries accessible with key.  It works by a bit of a trick.  It takes
the input ensemble and posts the original index positions to the trace
headers (metadata).  It then sorts the ensemble by key selecting the
nmax largest values.  It is VERY INMPORTANT to realize that because the
input ensemble d is is sorted here it is passed by reference and 
NOT by a TimeSeries& (reference) declaration.  This will mysteriously 
stop working if that one character were added.  I seriously considered
an explicit copy internally, but hope these comments are sufficient to
prevent a future error here.  

Warning 2:  key is very restricted right now to match allowed keywords
for ensemble sorting defined in XcorProcessingEngine.h */

list<int> find_largest(TimeSeriesEnsemble d, string key, int nmax)
{
	int i;
	int ensemblesize=d.member.size();
	for(i=0;i<ensemblesize;++i)
		d.member[i].put("index",i);
	if(key==amplitude_static_keyword)
		sort(d.member.begin(),d.member.end(),
			greater_metadata_double<TimeSeries,AMPLITUDE>());
	else if(key==snr_keyword || key==dbsnrkey)
	/* Warning:  this asumes db method copies it's snr value to the snr_keyword
	attribute.  This may do bad things if that is not the case */
		sort(d.member.begin(),d.member.end(),
			greater_metadata_double<TimeSeries,SNR>());
	/* Assume that dead traces get shoved to the end of the sorted output 
	so finding a dead trace is a signal to trucate the output list.
	We silently return an empty list if it ends that way, although it 
	really shouldn't here. */
	list<int> result;
	int ival;
	for(i=0;i<nmax,i<ensemblesize;++i)
	{
		if(!d.member[i].live) break;
		ival=d.member[i].get_int("index");
		result.push_back(ival);
	}	
	return(result);
}

/* List of options for deciding with ensemble is "best" if run in a
mode with multiple reference traces */
enum BestEnsembleMethod {MaxSumWeights, MaxSurvivors};
vector<TimeSeriesEnsemble>::iterator find_best_ensemble(vector<TimeSeriesEnsemble>& trials,
	BestEnsembleMethod bem, Metadata md)
{
	int i;
	vector<TimeSeriesEnsemble>::iterator result,tse;
	/*An odd inconsistency in C is that declarations aren't allowed in an individual 
	case portion of a switch block so we have to declare these here.  We further
	have to extract these from md.  keywords match those in dbxcor for obvious reasons.*/
	double stack_weight_cutoff=md.get_double("stack_weight_cutoff");
	double time_lag_cutoff=md.get_double("time_lag_cutoff");
	double coherence_cutoff=md.get_double("coherence_cutoff");
	double correlation_peak_cutoff=md.get_double("stack_weight_cutoff");
	double stackwt, lag, coherence, correlation;
	vector<int> survivorcounts;
	int count,countmax,ensemblesize;
	double sumwtmax,sumwt;
	switch (bem)
	{
	case MaxSumWeights:
		result=trials.end();
		sumwtmax=-1.0;
		for(tse=trials.begin();tse!=trials.end();++tse)
		{
			ensemblesize=tse->member.size();
			for(i=0,sumwt=0.0;i<ensemblesize;++i)
			{
				if(tse->member[i].live)
				{
					//Intentionally avoid a try block here for efficiency
					stackwt=tse->member[i].get_double(stack_weight_keyword);
					sumwt += stackwt;
				}
			}
			if(sumwt>sumwtmax) result=tse;
		}
		break;
	case MaxSurvivors:
	default:
		result=trials.end();
		countmax=-1;
		for(tse=trials.begin();tse!=trials.end();++tse)
		{
			ensemblesize=tse->member.size();
			for(i=0,count=0;i<ensemblesize;++i)
			{
				if(!tse->member[i].live) continue;
				correlation=tse->member[i].get_double(peakxcor_keyword);
				if(correlation<correlation_peak_cutoff) continue;
				stackwt=tse->member[i].get_double(stack_weight_keyword);
				if(stackwt<stack_weight_cutoff) continue;
				coherence=tse->member[i].get_double(coherence_keyword);
				if(coherence<coherence_cutoff) continue;
				lag=tse->member[i].get_double(moveout_keyword);
				if(fabs(lag)>time_lag_cutoff) continue;
				++count;
			}
			if(count>countmax) result=tse;
		}
			
	}
	return(result);
}


/* This program has a complicated parameter file with options that
may be mutually exclusive.  This procedure encapsulates all required
tests to avoid cluttering up main. 

returns true if parameters pass all tests, false with a posted message to
stderr if there are any issues. 

md is a Metadata representation of a control parameter file */
	
bool validate_parameters(Metadata md,XcorEngineMode mode,string qfile, string origindb)
{
	bool result(true);
	switch(mode)
	{
	case EventGathers:
		break;
	case GenericGathers:
		break;
	case ContinuousDB:
	default:
		if(origindb.size()<=0 || qfile.size()<=0)
		{
			cerr << "Usage error.  "
				<< "The -q and -e arguments are required "
				<< "when running in ContinuousDB processing mode"<<endl;
			result=false;
		}
			
	}
	/* Some required parameters that will cause a downstream abort if not
	set */
	double dtest;
	try {
		dtest=md.get_double("correlation_peak_cutoff");
	} catch (MetadataGetError mdge)
	{
		cerr << "Parameter file error:  missing required parameter "
			<< "correlation_peak_cutoff"<<endl;
		result=false;
	}
	/* this gets a bit repetitious, but unless this list gets real long this is simpler */
	try {
		dtest=md.get_double("stack_weight_cutoff");
	} catch (MetadataGetError mdge)
	{
		cerr << "Parameter file error:  missing required parameter "
			<< "stack_weight_cutof"<<endl;
		result=false;
	}
	try {
		dtest=md.get_double("time_lag_cutoff");
	} catch (MetadataGetError mdge)
	{
		cerr << "Parameter file error:  missing required parameter "
			<< "time_lag_cutoff"<<endl;
		result=false;
	}
	try {
		dtest=md.get_double("coherence_cutoff");
	} catch (MetadataGetError mdge)
	{
		cerr << "Parameter file error:  missing required parameter "
			<< "coherence_cutoff"<<endl;
		result=false;
	}
	
	return result;
}
/* Return this critical state enum value from parameter space */
XcorEngineMode get_engine_mode(Metadata& md)
{
	XcorEngineMode result(ContinuousDB);
	try {
		string pmodestr=md.get_string("processing_mode");
		if(pmodestr=="EventGathers")
			result=EventGathers;
		else if(pmodestr=="GenericGathers")
			result=GenericGathers;
		else if(pmodestr=="ContinuousDB")
			result=ContinuousDB;
		else
		{
			cerr << "Illegal keyword = "<<pmodestr<<" for processing_mode parameter"
				<<endl
				<< "Must be one of:  EventGathers, GenericGathers, or ContinuousDB"
				<<endl;
			exit(-1);
		}
		return(result);
	} catch (MetadataGetError mderr)
	{
		mderr.log_error();
		exit(-1);
	}
}
/* Similar to above to get reference trace method from Metadata */
ReferenceTraceMethod GetReferenceTraceMethod(Metadata& md)
{
	string stest;
	try {
		stest=md.get_string("ReferenceTraceMethod");
	}
	catch (MetadataGetError mderr)
	{
		mderr.log_error();
		exit(-1);
	}
	ReferenceTraceMethod result;
	if(stest=="UseAll")
		result=UseAll;
	else if(stest=="MaxSnrDB" )
		result=MaxSnrDB;
	else if(stest=="MaxSnrComputed")
		result=MaxSnrComputed;
	else if(stest=="MaxPeakAmplitude")
		result=MaxPeakAmplitude;
	else
	{
		cerr<< "Illegal keyword for parameter ReferenceTraceMethod"<<endl
			<<"Must be one of:  UseAll, MaxSnrDB, MaxSnrComputed, or MaxPeakAmplitude"
			<<endl;
		exit(-1);
	}
	return(result);
}
BestEnsembleMethod GetBestEnsembleMethod(Metadata& md)
{
	string stest;
	try {
		stest=md.get_string("BestEnsembleMethod");
	}
	catch (MetadataGetError mderr)
	{
		mderr.log_error();
		exit(-1);
	}
	BestEnsembleMethod result;
	if(stest=="MaxSumWeights")
		result=MaxSumWeights;
	else if(stest=="MaxSurvivors" )
		result=MaxSurvivors;
	else
	{
		cerr<< "Illegal keyword for parameter BestEnsembleMethod"<<endl
			<<"Must be one of:  MaxSumWeights, or MaxSurvivors"
			<<endl;
		exit(-1);
	}
	return(result);
}

Hypocenter load_hypo(DatascopeHandle& dbh,string method, string model)
{
	double lat,lon,depth,otime;
	try {
		lat=dbh.get_double("lat");
		lon=dbh.get_double("lon");
		lat=rad(lat);
		lon=rad(lon);
		depth=dbh.get_double("depth");
		otime=dbh.get_double("origin.time");
		return(Hypocenter(lat,lon,depth,otime,method,model));
	} catch (SeisppError serr)
	{
		serr.log_error();
		exit(-1);
	}
}
int apply_cutoff(TimeSeriesEnsemble& d, Metadata& md)
{
	double stack_weight_cutoff=md.get_double("stack_weight_cutoff");
	double time_lag_cutoff=md.get_double("time_lag_cutoff");
	double coherence_cutoff=md.get_double("coherence_cutoff");
	double correlation_peak_cutoff=md.get_double("stack_weight_cutoff");
	vector<TimeSeries>::iterator memptr;
	int nbad;

	nbad=0;
	for(memptr=d.member.begin();memptr!=d.member.end();++memptr)
	{
		double correlation, coherence, stackwt, lag;
		if(memptr->live)
		{
			/*Intentionally omit try/catch as here we are
			guaranteed these have been set */
			correlation=memptr->get_double(peakxcor_keyword);
			if(correlation<correlation_peak_cutoff)
			{
				memptr->live=false;
				++nbad;
				continue;
			}
			coherence=memptr->get_double(coherence_keyword);
			if(coherence<coherence_cutoff)
			{
				memptr->live=false;
				++nbad;
				continue;
			}
			stackwt=memptr->get_double(stack_weight_keyword);
			if(stackwt<stack_weight_cutoff)
			{
				memptr->live=false;
				++nbad;
				continue;
			}
			lag=memptr->get_double(moveout_keyword);
			if(fabs(lag)>time_lag_cutoff)
			{
				memptr->live=false;
				++nbad;
			}
		}
	}
}
/* This template and related procedures are used to provide a way to 
resort an ensemble in the original order. */
const string index_keyword("member_number");
template <typename T> struct less_integer_header_attribute
		: public binary_function<T,T,bool>
{
	string key;
	less_integer_header_attribute(string keyin)
	{
		key=keyin;
	};
	bool operator () (T a, T b)
	{
		int ia,ib;
		try {
			ia=a.get_int(key);
		}catch(...){return true;};
		try {
			ib=b.get_int(key);
		}catch(...){return false;};
		if(ia<ib)
			return true;
		else
			return false;
	}
};
			
	
void set_index(TimeSeriesEnsemble *tse)
{
	int i;
	for(i=0;i<tse->member.size();++i)
	{
		tse->member[i].put(index_keyword,i);
	}
}


/* This function kills traces in parent ensemble using pass1 pattern.  
The method used is very very dependent on an assumption that parent is
already in index_keyword order.  Thus, pass1 is sorted by index_keyword
and we assume we can do a kill in index_keyword order.  Said another
way we assume that parent.member[i] is the same data as pass1.member[i]
for all i.  
*/
void kill_pass1_bad(TimeSeriesEnsemble& pass1, TimeSeriesEnsemble& parent)
{
	if(pass1.member.size() != parent.member.size()) throw SeisppError(string("kill_pass1_bad:  ")
				+ "pass1 and parent ensemble sizes do not match.");
	// Sort pass1 by index_keyword using the generic integer ordering function
	sort(pass1.member.begin(),pass1.member.end(),less_integer_header_attribute<TimeSeries>(index_keyword));
	int i;
	for(i=0;i<pass1.member.size();++i)
		if(!pass1.member[i].live) parent.member[i].live=false;
	
}

void usage()
{
	cerr << "hpcxcor db [-o dbout -e origindb -q queuefile -pf pffile]"<<endl;
	exit(-1);
}
int main(int argc, char **argv)
{
	int i;
	if(argc<2) usage();
        elog_init(argc,argv);
#ifdef MPI_SET

/*
        Initialize the MPI parallel environment, MPI_Init
        must be called before any other MPI call.
        rank -- rank of a specific process
        mp -- the total number of processes. This whole
                process group is called MPI_COMM_WORLD
*/
        int rank, np;

        MPI_Init(&argc, &argv);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &np);
#endif

	string pfname("hpcxcor");
	string dbname(argv[1]);
	string dbout(dbname);  // same as input by default
	string queuefile("");
	string origindb("");
	for(i=2;i<argc;++i)
	{
		string argtest(argv[i]);
		if(argtest=="-pf")
		{
			++i;
			if(i>=argc)usage();
			pfname=string(argv[i]);
		}
		else if(argtest=="-q")
		{
			++i;
			if(i>=argc)usage();
			queuefile=string(argv[i]);
		}
		else if(argtest=="-o")
		{
			++i;
			if(i>=argc)usage();
			dbout=string(argv[i]);
		}
		else if(argtest=="-e")
		{
			++i;
			if(i>=argc)usage();
			origindb=string(argv[i]);
		}
		else
			usage();
	}
	Pf *pf;
	if(pfread(const_cast<char *>(pfname.c_str()),&pf))
	{
		cerr << "Error reading pf file = "<<pfname<<endl;
		exit(-1);
	}
	try {
		string refkey("reference_trace_used");
		AttributeMap("css3.0");
		Metadata global_md(pf);
		XcorAnalysisSetting asetting(global_md);
		// Window parameters
		double ts,te;
		ts=global_md.get_double("beam_window_start_time");
		te=global_md.get_double("beam_window_end_time");
		TimeWindow beamtw(ts,te);
		asetting.set_beam_tw(beamtw);
		ts=global_md.get_double("robust_window_start_time");
		te=global_md.get_double("robust_window_end_time");
		TimeWindow robusttw(ts,te);
		asetting.set_robust_tw(robusttw);
		ts=global_md.get_double("SNR_signal_window_start_time");
		te=global_md.get_double("SNR_signal_window_end_time");
		TimeWindow SNR_signal_tw(ts,te);
		ts=global_md.get_double("SNR_noise_window_start_time");
		te=global_md.get_double("SNR_noise_window_end_time");
		TimeWindow SNR_noise_tw(ts,te);
		string method(global_md.get_string("TTmethod"));
		string model(global_md.get_string("TTmodel"));
		int MaxTrials=global_md.get_int("MaximumReferenceStationTrials");
		ReferenceTraceMethod refmeth=GetReferenceTraceMethod(global_md);
		BestEnsembleMethod bem=GetBestEnsembleMethod(global_md);
		/* This is an ugly way to get this key attribute from parameter
		space, but stuck with it.  This should really be available to the 
		user from the user interface to XcorProcessingEngine, but I 
		don't care to burden the interface with this so will duplicate
		parsing.  This is, or course, a serious maintenance issue because
		if the keyword is changed once place it will break the other.  We
		do it here because this will exit the program if there are problems
		cracking this and the valid_parameters procedure has total
		dependency on this key attribute for the engine. */
		XcorEngineMode enginemode=get_engine_mode(global_md);
		/* This is here to do some sanity checks on pf an exit if
		there are inconsisties.   This function should return a false
		if parameters are inconsistent.*/
		if(!validate_parameters(global_md,enginemode,queuefile,origindb))
		{
			cerr << "hcpxcor:  parameter inconsistency.  Cannot continue."
				<< "Fix errors noted above and try again"<<endl;
			exit(-3);
		}
		XcorProcessingEngine xpe(pf,asetting,dbname,dbout,queuefile);
		DatabaseHandle *dbhwf = xpe.get_db(string("waveformdb"));
		Hypocenter h;  /* Only used in ContinousDB mode, but have to declare now */
		/* In continuous mode we need to loop through the event database.  We
		do this by building a handle to origindb and a linked queue */
		DatascopeHandle dbhorigin;
		DatascopeProcessingQueue *dpq;  //Must be pointer as there is no opeator = for this
		if(enginemode==ContinuousDB)
		{
			dbhorigin=DatascopeHandle(origindb,true);
			dbhorigin.lookup("event");
			dbhorigin.natural_join("origin");
			dbhorigin.subset("orid==prefor");
			if(dbhorigin.number_tuples()<=0)
			{
				cerr << "event->origin prefor subset is empty.  Check database "
					<< origindb<<endl;
				exit(-1);
			}
			dbhorigin.rewind();
			dpq=new DatascopeProcessingQueue(dbhorigin,queuefile);
		}

		/* Top of loop over gathers.   Algorithm loops through each gather.
		Eventually plan an mpi split here.  i.e. this loop is a perfect
		data-driven, embarassingly parallel process. */
		int gathernumber(0);
		bool keeplooping(true);
		/* We push sorted copies of each ensemble to this vector of ensembles.
		This could get out of control for large ensembles, but we'll try it 
		until it proves problematic. */
		vector<TimeSeriesEnsemble> trials;
		int ensemblenumber=-1;  // -1 because we have to increment it at the top of the loop
		int load_status;  // Return code from load_data used to test for end of data 
		do {
		    ++ensemblenumber;
#ifdef MPI_SET
		    if(gathernumber%np == rank)
		    {
#endif
			
			/* load data.  Unfortunately, different modes do this differently */
			switch(enginemode)
			{
			case EventGathers:
			case GenericGathers:
				/* this method is odd because it marks the outcome of 
				the last gather processed.  Works correctly for first pass
				doing nothing to the queue. */
				load_status=xpe.load_data(*dbhwf,FINISHED);
				if(load_status!=0) keeplooping=false;
				break;
			case ContinuousDB:
				dpq->set_to_current(dbhorigin);
				h=load_hypo(dbhorigin,method,model);
				xpe.load_data(h);
			};
			if(enginemode==ContinuousDB)
			{
				if(!dpq->has_data()) keeplooping=false;
			}
			if(!keeplooping) break;
			TimeSeriesEnsemble *tse=xpe.get_waveforms_gui();
			int ensemblesize=tse->member.size();
			/* Silently skip any ensemble with only one member */
			/* For initial debug, require at least 3 members */
			if(ensemblesize<4) continue;
			/* We have to deal with the critical issue of setting
			the event to be used as the reference trace for the MultichannelCorrelator.
			We handle that here through an enum with different methods called.
			This will make it easier to expand new options that might surface.
			Each procedure returns a list of ensemble members to be used in 
			a succession of trials.  Note for maintenance:  be careful xpe
			does not unintentionally sort tse since we are given raw pointer
			to the original data that will be sorted during analysis, for example. */
			list<int> reference_traces;
			list<int>::iterator rtptr;
			switch(refmeth)
			{
			case UseAll:
				reference_traces=UseAll_ChooseReference(*tse,MaxTrials);
				break;
			case MaxSnrDB:
				reference_traces=MaxSnrDB_ChooseReference(*tse,dbsnrkey);
				if(reference_traces.size()<1) break;  // handle this error below
				if(MaxTrials>1) reference_traces=find_largest(*tse,dbsnrkey,MaxTrials);
				break;
			case MaxSnrComputed:
				reference_traces=MaxSNRComputed_ChooseReference(*tse,SNR_signal_tw, SNR_noise_tw);
				if(reference_traces.size()<1) break;  // handle this error below
				if(MaxTrials>1) reference_traces=find_largest(*tse,snr_keyword,MaxTrials);
				break;
			case MaxPeakAmplitude:
				reference_traces=MaxPeakAmplitude_ChooseReference(*tse);
				if(reference_traces.size()<1) break;  // handle this error below
				if(MaxTrials>1) reference_traces=find_largest(*tse,amplitude_static_keyword,MaxTrials);
			}
			if(reference_traces.size()<=0) 
			{
				cerr << "WARNING:  Reference trace selection method failed for ensemble number "
					<< ensemblenumber <<endl
					<< "Skipping these data"<<endl;
				if(enginemode==ContinuousDB)
				{
					if(!dpq->has_data()) break;
				}
				else
				{
					continue;
				}
			}
			trials.clear();
			for(rtptr=reference_traces.begin(),i=0;
				rtptr!=reference_traces.end();
						++rtptr,++i)
			{
				asetting.set_ref_trace(*rtptr);
				/* necessary because analysis call nearly always scrambles ensemble.
				This is a bit inefficient, but hard to avoid without creating
				havoc. */
				if(i!=0) xpe.restore_original_ensemble();
				xpe.change_analysis_setting(asetting);
				/* Must refresh this pointer since restore_original_ensemble is
				guaranteed to change it */
				tse=xpe.get_waveforms_gui();
				set_index(tse);
				/*It is convenient to post this to the ensemble to avoiding having
				to recompute it for pass2 */
				tse->put(refkey,*rtptr);
				xpe.analyze();
				trials.push_back(*tse);
			}
			/*Select the version considered best */
			vector<TimeSeriesEnsemble>::iterator bestens;
			/* Passing this the global_md is a backdoor way to get 
			any parameters needed by any of the allowed methods. 
			Not best practice, but preferable to modifying the engine. */
			bestens=find_best_ensemble(trials, bem,
				global_md);
			/* This procedure applies cutoff criteria to best ensemble returning
			a count of how many members were marked bad.  We pass these again through
			the global_md object. */
			int ibad=apply_cutoff(*bestens, global_md);
			if(ibad>0 && (ensemblesize-ibad>2))
			{
				/* When data get deleted we always rerun the analysis.  This 
				is more complicated than it might seem as bestens is a copy 
				of original.  We first refresh the master, then edit it against
				bestens to mark the data that failed in the first pass as bad. */
				xpe.restore_original_ensemble();
				tse=xpe.get_waveforms_gui();
				set_index(tse);
				kill_pass1_bad(*bestens, *tse);
				int iref=tse->get_int(refkey);
				asetting.set_ref_trace(iref);
				xpe.change_analysis_setting(asetting);
				xpe.analyze();
			}
			/* We cannot use the save_results methods in the engine 
			because it would be very evil to delete and modify the
			tse pointer.  Hence we largely duplicate code from XcorProcessingEngine
			in the procedure called here.  see above*/
			//save_results(enginemode,xpe,*bestens);

#ifdef MPI_SET
		    }
#endif
		    ++gathernumber;
		    if(enginemode==ContinuousDB) ++dbhorigin;
		    /* Because of inconsistent input methods, we have this ugly construct to 
		    break this loop */
		    switch(enginemode)
		    {
		    case EventGathers:
		    case GenericGathers:
			if(load_status) keeplooping=false;
			break;
		    case ContinuousDB:
			if(!dpq->has_data()) keeplooping=false;
		    };
			
		}
		while(keeplooping);
	}
	catch(MetadataError mde)
	{
		mde.log_error();
	}
	catch(SeisppError spe)
	{
		spe.log_error();
	}
}
