#include "seispp.h"
#include "XcorProcessingEngine.h"
#ifdef MPI_SET
        #include <mpi.h>
#endif
using namespace std;
using namespace SEISPP;
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
*/

/* List of methods used to select one or more traces as reference event
for MultichannelCorrelator */
enum ReferenceTraceMethod {UseAll, MaxSnrDB, MaxSnrComputed, MaxPeakAmplitude};
/* This is a helper used by some of the reference trace methods below.
Returns a integer of first live trace.  -1 means no live data were found.
Note this could easily have been templatized for any ensemble type in SEISPP. */
int find_first_live(TimeSeriesEnsemble& d)
{
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
		result.reserve(ensemblesize);
		for(i=0;i<ensemble_size;++i) result.push_back(i);
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
	srnvalues.reserve(ensemblesize);
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
			if(SeisppVerbose) 
			{
				cerr << "Error MaxSnrDB_ChooseReference: "<<endl;
				mdge.log_error();
			}
			val=-1.0;
		}
		snrvalues.push_back(val);
	}
	if(number_fetched<=0)
	{
		// This maybe should be only used if verbose is set, but this seems wise to 
		// me to always log this.
		ival=find_first_live(d);
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
list<int> MaxSnrComputed_ChooseReference(TimeSeriesEnsemble& d, TimeWindow signal, TimeWindow noise)
{
	list<int> result;
	int imax;
	double snrmax;
	vector<TimeSeries>::iterator iptr;
	/* This loop never tests live boolean because it depends on undocumented fact that
	SNR_rms tests this and returns a -1 immediately if so marked.*/
	imax=0;  // assume we aren't called if d is empty
	iptr=d.member.begin();
	snrmax=SNR_rms(*iptr,signal,noise);
	++iptr;
	for(i=1;i!=d.member.end();++iptr,++i)
	{
		double testval;
		testval=SNR_rms(*iptr,signal,noise);
		if(testval>snrmax)
		{
			snrmax=testval
			imax=i;
		}
	}
	/* Take same philosophy as above if all failed. */
	if(snrmax<0.0)
double SNR_rms(TimeSeries& d, TimeWindow signal, TimeWindow noise)

}
list<int> MaxPeakAmplitude_ChooseReference(TimeSeriesEnsemble& d, TimeWindow signal, TimeWindow noise)
{
double PeakAmplitude(TimeSeries *p)

}
		
		
}
/* List of options for deciding with ensemble is "best" if run in a
mode with multiple reference traces */
enum BestEnsembleMethod {MaxSumWeights, MaxSurvivors};


/* Finds the reference trace as the index position by a
simple linear search for maximum signal to noise ratio. 
May eventually want a more complex recipe.*/
int choose_reference(DatascopeHandle& dbh,
		TimeSeriesEnsemble *d)
{
	string maxsta,sta;
	double maxsnr,snr;
	int i;
	try {
	dbh.rewind();
	for(i=0;i<dbh.number_tuples();++i,++dbh)
	{
		if(i==0)
		{
			maxsnr=dbh.get_double("snr");
			maxsta=dbh.get_string("sta");
		}
		else
		{
			snr=dbh.get_double("snr");
			sta=dbh.get_string("sta");
			if(snr>maxsnr)
			{
				maxsta=sta;
				maxsnr=snr;
			}
		}
	}
	} catch (SeisppDberror se)
	{
		se.log_error();
		cerr << "hpcxcor defaulting to first trace as reference"<<endl;
		return(0);
	}
	// Now search for this station in the ensemble
	for(i=0;i<d->member.size();++i)
	{
		if(d->member[i].live)
		{
			sta=d->member[i].get_string("sta");
			if(sta==maxsta)  return(i);
		}
	}
	cerr << "hpcxcor data mismatch in search for reference trace"<<endl
		<<"arrival maximum snr station="<<sta<<" has no trace data  "
		<<"Default to 0"<< endl;
	return(0);
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
		if(origindb.size()<=0 || qfile.size()<=0)
			cerr << "Usage error.  "
				<< "The -q and -e arguments are required "
				<< "when running in ContinuousDB processing mode"<<endl;
			result=false;
			
	default:
	};
	/* dummy for now, always true */
	return result;
}
/* Return this critical state enum value from parameter space */
XcorEngineMode get_engine_mode(Metadata& md)
{
	XcorEngineMode result;
	try {
		string pmodestr=md.get_string("procesing_mode");
		if(pmodestr=="EventGathers")
			result=EventGathers;
		else if(pmodestr=="GenericGathers")
			result==GenericGathers;
		else if(pmodstr=="ContinuousDB")
			result=ContinuousDB;
		else
		{
			cerr << "Illegal keyword = "<<pmodestr<<" for procesing_mode parameter"
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

void usage()
{
	cerr << "hpcxcor db [-o dbout -e origindb -q queuefile -pf pffile]"<<endl;
	exit(-1);
}
bool SEISPP::SEISPP_verbose(false);
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
		if(argtest=="-q")
		{
			++i
			if(i>=argc)usage();
			queuefile=string(argv[i]);
		}
		if(argtest=="-o")
		{
			++i
			if(i>=argc)usage();
			dbout=string(argv[i]);
		}
		if(argtest=="-e")
		{
			++i
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
		string method(global_md.get_string("TTmethod"));
		string model(global_md.get_string("TTmodel"));
		int MaxTrials=global_md.get_int("MaximumReferenceStationTrials");
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
		if(validate_parameters(global_md,enginemode,queuefile,origindb))
		{
			cerr << "hcpxcor:  parameter inconsistency.  Cannot continue."
				<< "Fix errors noted above and try again"<<endl;
			exit(-3);
		}
		XcorProcessingEngine xpe(pf,asetting,dbname,dbout,queuefile);
		/* Not totally necessary in this context, but a useful test to 
		require anyway*/
		if(!xpe.validate_setting()) 
		{
			cerr << "XcorProcessingEngine::validate_setting returned false"
				<< endl <<"Cannot continue.  Exiting. "<<endl;
			exit(-2);
		}
		/* The engine contains an image of the waveform db view.  We need it.*/
		DatabaseHandle *dbhwf = xpe.get_db(dbname);
		Hypocenter h;  /* Only used in ContinousDB mode, but have to declare now */
		/* In continuous mode we need to loop through the event database.  We
		do this by building a handle to origindb and a linked queue */
		DatascopeHandle dbhorigin;
		if(enginemode==ContinuousDB)
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
			dpq=DatascopeProcessingQueue(dbhorigin,queuefile);
		}

		/* Top of loop over gathers.   Algorithm loops through each gather.
		Eventually plan an mpi split here.  i.e. this loop is a perfect
		data-driven, embarassingly parallel process. &/
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
				break;
			case ContinuousDB:
				dpq.set_to_current(dbhorigin);
				h=load_hypo(dbhorigin);
				xpe.load_data(h);
			};
			TimeSeriesEnsemble *tse=xpe.get_waveforms_gui();
			/* Silently skip any ensemble with only one member */
			if(tse->member.size()<2) 
			{
				if(enginemode==ContinuousDB)
				{
					if(!dpq.has_data()) break;
				}
				continue;
			}
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
enum ReferenceTraceMethod {UseAll, MaxSnrDB, MaxSnrComputed, MaxPeakAmplitude};
			switch(ReferenceTraceMethod)
			{
			case UseAll:
				reference_traces=UseAll_ChooseReference(*tse,MaxTrials);
			case MaxSnrDB:
			case MaxSnrComputed:
			case MaxPeakAmplitude:
			};
			if(reference_traces.size()<=0) 
			{
				cerr << "WARNING:  Reference trace selection method failed for ensemble number "
					<< ensemblenumber <<endl
					<< "Skipping these data"<<endl;
				if(enginemode==ContinuousDB)
				{
					if(!dpq.has_data()) break;
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
				xpe.analyze();
				/* Must refresh this pointer since restore_original_ensemble is
				guaranteed to change it */
				tse=xpe.get_waveforms_gui();
				trials.push_back(*tse);
			}
			/*Select the version considered best and save */
			vector<TimeSeriesEnsemble>::iterator bestens;
			/* Passing this the global_md is a backdoor way to get 
			any parameters needed by any of the allowed methods. 
			Not best practice, but preferable to modifying the engine. */
			bestens=find_best_ensemble(trials, BestEnsembleMethod,
				global_md);
			/* We cannot use the save_results methods in the engine 
			because it would be very evil to delete and modify the
			tse pointer.  Hence we largely duplicate code from XcorProcessingEngine
			in the procedure called here.  see above*/
			save_results(enginemode,xpe,*bestens);
				

#ifdef MPI_SET
		    }
#endif
		    ++gathernumber;
		    if(enginemode==ContinousDB) ++dbhorigin;
		    /* Because of inconsistent input methods, we have this ugly construct to 
		    break this loop */
		    switch(enginemode):
		    {
		    case EventGathers:
		    case GenericGathers:
			if(load_status) keeplooping=false;
			break;
		    case ContinuousDB:
			if(!dpq.has_data()) keeplooping=false;
		    };
			
		}
		while(keeplooping);
		// The current method has a default filter applied
		// immediately after reading.  Here that is the filter.
		// To extend to multiple filters will require 
		// a different approach
		xpe.load_data(h);
		// Now we need to select a reference trace.  We use
		// snr from the tmp database. xpe wants the index
		// position of this trace that we have to search for.
		TimeSeriesEnsemble *ens=xpe.get_waveforms_gui();
cout << "data loaded.  Number of members in this ensemble="
	<<ens->member.size()<<endl;
		dbh.natural_join(string("assoc"));
		dbh.natural_join(string("arrival"));
		int reftrace=choose_reference(dbh,ens);
		asetting.set_ref_trace(reftrace);
cout << "Using reference trace ="<<reftrace<<endl;
		xpe.change_analysis_setting(asetting);
cout << "Running correlator"<<endl;
		MultichannelCorrelator *mce=xpe.analyze();
		orid=dbnextid(dbh.db,"orid");
cout << "Saving results tagged with orid="<<orid<<endl;
		xpe.save_results(evid,orid,h);
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
