#include "Metadata.h"
#include "dbpp.h"
#include "MatlabProcessor.h"
using namespace std;
using namespace SEISPP;
/*! \brief Generic algorithm to load arrival times from a database.

In passive array processing a very common need is to extract time windows
around marked phase pick or a theoretical arrival time.  For measured
arrival times the css database has an awkward link to waveforms that 
causes problems when dealing with continuous data.  This is one element of
a set of functions aimed at dealing with this problem.  

This particular procedure aims to take an input ensemble of data and 
find matches for arrival times in an external database.  For each member
of the ensemble with a matching arrival it posts the arrival time to the
generalized header (Metadata) of the parent object.  To avoid testing 
for valid arrivals in the ensemble after this procedure is run a index
of data with valid arrivals is returned as an stl list container of ints.

\param dat is the input ensemble passed as a reference to the vector 
	container of data objects.  The type of the objects in the container
	is generic except class T MUST be an object that inherits Metadata.
	At time of this writing T could be TimeSeries, ThreeComponentSeismogram,
	or ComplexTimeSeries. 
\param dbh match handle to a Datascope database containing arrival. 
	The algorithm actually only cares that a find using the Metadata
	area of the vector data components will provide a unique match
	into the table the handle references (normally arrival).  This 
	tacitly assumes the handle has been properly constructed and the
	proper attributes have been loaded with the data to provide a unique
	match into arrival.  With css this REQUIRES that each component of
	dat MUST have had orid posted before calling this function.  There
	is no other unambiguous way to match waveforms in this context to
	a unique arrival. 
\param keyword is the attribute name used to posted the arrival time data.

\return list of members of input vector with valid arrival times posted.
	The contents of the list can be traversed and used with operator[]
	to extract valid data from the ensemble.  
*/
template <class T>
	list<int> LoadArrivalTimes(vector<T>& dat, DatascopeMatchHandle& dbh,
		const string keyword)
{
	vector<T>::iterator d;
	int i;
	list<int> data_with_arrivals;
	const string base_error("Warning (LoadArrivalTimes): ");
	for(dat.start(),i=0;d!=dat.end();++d,++i)
	{
		double atime;  //  arrival time.  
		if(d->live)
		{
		// First see if there is an arrival for this
		// station.  If not, skip it. 
			list<int> records
				=dbh.find(dynamic_cast<Metadata&>(*d));
			// if no arrival silently skip data for this station
			if(records.size()<=0) continue;
			if(records.size()>1)
			{
				string sta=d->get_string("sta");
				cerr << base_error 
					<< "found "
					<< records.size()
					<< " arrivals for station "
					<< sta <<endl
					<< "Using first found "
					<< "in database view"<<endl;
			}
			Dbptr db=dbh.db;
			// tricky usage here.  begin() returns
			// an iterator so the * operator gets the
			// value = record number of the match
			db.record=*(records.begin());
			if(dbgetv(db,0,"arrival.time",&atime,0)
				== dbINVALID) 
			{
				string sta=d->get_string("sta");
				cerr << base_error
					<< "dbgetv failed"
					<< " in attempt to obtain"
					<< " arrival time for station"
					<< sta << endl
					<< "Data from this station"
					<< " will be dropped"<<endl;	
			}
			d->put(keyword,atime);
			data_with_arrivals.push_back(i);
		}
	}
	return(data_with_arrivals);
} 
ThreeComponentEnsemble *BuildRegularGather(ThreeComponentEnsemble& raw,
	DatascopeMatchHandle& dbh, 
	ResamplingDefinitions& rdef,
	double target_dt)
{
	list<int>  records;  // return of find function for match handle
	const samprate_tolerance(0.01);  // fractional sample rate tolerance
	vector<ThreeComponentSeismogram>::iterator d;
	nmembers=raw.member.size();
	TimeSeries *x1,*x2,*x3;
	ThreeComponentEnsemble *result;
	result = new ThreeComponentEnsemble(raw);
	// An inefficiency here, but this allow us to discard dead
	// traces and problem data from ensemble as we assemble the
	// new one.
	result->member.clear();
	// Use an iterator as it simplifies the symbols.  could
	// use operator[], but this is faster and actually cleaner here.
	for(d=raw->member.start();d!=raw->member.end();++d)
	{
		double atime;  //  arrival time.  
		if(d->live)
		{
			// First see if there is an arrival for this
			// station.  If not, skip it. 
			try { 
				list<int> records
					=dbh.find(dynamic_cast<Metadata&>(*d));
				// if no arrival skip data for this station
				if(records.size()<=0) continue;
				if(records.size()>1)
				{
					string sta=d->get_string("sta");
					cerr << "Warning:  found"
						<< records.size()
						<< " arrivals for station "
						<< sta <<endl
						<< "Using first found "
						<< "in database view"<<endl;
				}
				Dbptr db=dbh.db;
				// tricky usage here.  begin() returns
				// an iterator so the * operator gets the
				// value = record number of the match
				db.record=*(records.begin());
				if(dbgetv(db,0,"arrival.time",&atime,0)
					== dbINVALID) 
				{
					string sta=d->get_string("sta");
					cerr << "Warning:  dbgetv failed"
						<< " in attempt to obtain"
						<< " arrival time for station"
						<< sta << endl
						<< "Data from this station"
						<< " will be dropped"<<endl;
				}
				d->put(
			} 
			catch (SeisppError serr) 
			{
				serr.log_error();
				string sta=d->get_string("sta");
				cerr << "Data for for station "
					<< sta << " dropped"<<endl;
			} 
				
			// partial clone used to hold result
			ThreeComponentSeismogram new(*d);  
			x1=ExtractComponent(*d,0);
			x2=ExtractComponent(*d,1);
			x3=ExtractComponent(*d,2);
			if( (abs( (d->dt)-target_dt)/target_dt) > samprate_tolerance)
			{
				x1=ResampleTimeSeries(x1,rdef,target_dt,false);
				x2=ResampleTimeSeries(x2,rdef,target_dt,false);
				x3=ResampleTimeSeries(x3,rdef,target_dt,false);
			}
			

}
void usage()
{
	cerr << "subarraydecon db [-b beamdb -s beamsubsetexpression -pf pfname]" << endl;
	exit(-1);
}
int main(int argc, char **argv)
{
	// This is a C++ thing.  Not required on all platforms, but
	// recommended by all the books.  
	ios::sync_with_stdio();
	// This is a standard C construct to crack the command line
	if(argc<2) usage();
	string pfin(argv[0]);
	string dbin(argv[1]);
	string beamdb=dbin;
	bool subsetbeam=false;
	string beam_subset_expression;
        for(int i=2;i<argc;++i)
        {
                if(!strcmp(argv[i],"-V"))
                        usage();
                else if(!strcmp(argv[i],"-pf"))
                {
                        ++i;
                        if(i>=argc) usage();
                        pfin = string(argv[i]);
                }
		else if(!strcmp(argv[i],"-b"))
		{
                        ++i;
                        if(i>=argc) usage();
			beamdb=string(argv[i]);
		}
		else if(!strcmp(argv[i],"-s"))
		{
                        ++i;
                        if(i>=argc) usage();
			subsetbeam=true;
			beam_subset_expression=string(argv[i]);
		}
                else
                        usage();
        }
	// Standard Antelope C construct to read a parameter file
	// Well almost standard.  const_cast is a C++ idiom required
	// due to a collision of constantness with Antelope libraries.
	Pf *pf;
        if(pfread(const_cast<char *>(pfin.c_str()),&pf)) 
	{
		cerr << "pfread error for pf file="<<pfin<<".pf"<<endl;
		exit(-1);
	}
	Metadata control(pf);
	MetadataList mdbeam=pfget_mdlist(pf,"Beam_mdlist");
	MetadataList mdens=pfget_mdlist(pf,"Ensemble_mdlist");
	MetadataList mdtrace=pfget_mdlist(pf,"Trace_mdlist");
	try {
		// Get set of control parameters
		string netname=control.get_string("netname");
		string phase=control.get_string("phase");
		double ts,te;
		ts=control.get_double("data_window_start");
		te=control.get_double("data_window_end");
		TimeWindow datatwin(ts,te);
		double tpad=control.get_double("data_time_pad");
		StationchannelMap stachanmap(pf);
		string schema=control.get_string("AttributeMap_schema");
		// First get all the database components assembled.
		// There are three here:  input data, output data, and beam data
		DatascopeHandle dbh(dbin,false);
		DatascopeHandle dbhbeam(dbh);
		if(dbbeam!=dbin)
			dbhbeam=DatascopeHandle(dbbeam,true);
		DatascopeHandle dbhbv(dbhbeam,pf,"dbprocess_commands_for_beams");
		DatascopeHandle dbho(dbh);
		dbho.lookup(string("wfprocess"));
		// subset the beam database view if required
		if(subsetbeam)
		{
			dbhbv.subset(beam_subset_expression)
			if(dbhbv.number_tuples()<=0) 
			{
				cerr << "No data in beam db view "
					<< "after subset condition = "
					<< beam_subset_expression <<endl
					<< "Exiting"<<endl;
				exit(-1);
			}
		}
		dbhbv.rewind();  // Not essential, but good practice to initialize
		AttributeMap am(schema);  // needed for ensemble constructor
		// Launch matlab and sets up the communication channels
		MatlabProcessor mp(stdout);
		// 3C ensemble are loaded into Matlab as three matrices with these names
		string chans[3]={"x1","x2","x3"};
		// These are pointers to dynamic objects used in the loop
		// below.  
		ThreeComponentEnsemble *rawdata=NULL, *regular_gather=NULL,
				*decondata=NULL;
		SeismicArray *stations=NULL;
		int record;
		for(record=0;record<dbhbv.number_tuples();++record,++dbh)
		{

			// Read the beam trace and load it as beam
			TimeSeries beam(dynamic_cast<DatabaseHandle>(dbhbv),
				mdbeam,am);
			mp.load(beam,string("beam"));
			lat=beam.get_double(string("lat"));
			lon=beam.get_double(string("lon"));
			depth=beam.get_double(string("depth"));
			otime=beam.get_double(string("origin.time"));
			Hypocenter hypo(lat,lon,depth,otime);
			// On the first record we need to load the station
			// geometry object
			if(record==0)
			{
				stations = new SeismicArray(
					dynamic_cast<DatabaseHandle>(dbh),
					hypo.time,netname);
			}
			else
			{
				TimeWindow test(hypo.time,hypo.time+2000.0);
				if(!stations->GeometryIsValid(twin)
				{
					delete stations;
					stations = new SeismicArray(
					dynamic_cast<DatabaseHandle>(dbh),
					hypo.time,netname);
				}
			}
			// Read the raw data using the time window based constructor
			rawdata=array_get_data(*stations,hypo,
				phase,twin,tpad,dynamic_cast<DatabaseHandle&>(dbh),
				stachanmap,mdens,mdtrace,am);
			regular_gather=BuildRegularGather(rawdata, + other args);
			delete rawdata;
			// regular gather is now assumed to contain data with
			// a common start time.  irregular start time is allowed
			// with matlab interface
			mp.load(*regular_gather,chans);
			mp.process(string("beamdecon"));
			dmatrix x1d=mp.retrieve(string("x1d"));
			dmatrix x2d=mp.retrieve(string("x2d"));
			dmatrix x3d=mp.retrieve(string("x3d"));
			decondata=AssembleDeconGather(regular_gather,x1d,x2d,x3d);
			delete regular_gather;
			SaveResult(dbho,decondata,mdlo);
			delete decondata;
		}
	}
	catch (SeisppError serr)
	{
		serr.log_error();
	}
	catch (MetadataGetError mde)
	{
		mde.log_error();
	}
	catch (...)
	{
		cerr << "Something threw an unexpected exception"<<endl;
	}
}


pseudocode for this algorithm:
1) crack command line as usual.  expected usage prog dbdata [dbbeam -s subset]
2) create control Metadata object from pf
3) create DatascopeHandle of input data to be used for time window reads
3a) build DatascopeHandle for output data
4) create DatascopeHandle for beam data
5) subset beam data view 
6 foreach row in beam db
	7) create beam trace
	8) send beam trace to matlab
	9) get time window 
	10) construct ensemble in arrival time reference frame
	11) send ensemble to matlab
	12) retrieve deconvolved ensemble
	13) save results to db as 3c data stored in wfprocess
14) endloop

