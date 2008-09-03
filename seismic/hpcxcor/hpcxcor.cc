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
		int load_status;  // Return code from load_data used to test for end of data 
		do {
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
			case MaxSnrDB:
			case MaxSnrComputed:
			case MaxPeakAmplitude:
			};
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
