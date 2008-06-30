#include "SeisppKeywords.h"
#include "TimeSeries.h"
#include "ThreeComponentSeismogram.h"
#include "ensemble.h"
#include "Metadata.h"
#include "Hypocenter.h"
#include "dbpp.h"
#include "filter++.h"
#include "resample.h"
#include "seispp.h"
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
	Normally this should reference the standard catalog view.  That
	is:  event->origin->assoc->arrival subsetted to orid==prefor.
	The algorithm actually only cares that a find using the Metadata
	area of the vector data components will provide a unique match
	into the table the handle references (normally arrival).  This 
	tacitly assumes the handle has been properly constructed and the
	proper attributes have been loaded with the data to provide a unique
	match into arrival.  With css this REQUIRES that each component of
	dat MUST have had evid and/or orid posted before calling this function.  There
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
	for(d=dat.begin(),i=0;d!=dat.end();++d,++i)
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
			char csta[10];
			if(dbgetv(db,0,"arrival.time",&atime,"sta",csta,0)
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
	double target_dt,
	TimeWindow processing_window)
{
	const string arrival_keyword("arrival.time");
	const samprate_tolerance(0.01);  // fractional sample rate tolerance
	int nmembers=raw.member.size();
	auto_ptr<TimeSeries> x1,x2,x3;
	ThreeComponentEnsemble *result;
	result = new ThreeComponentEnsemble(raw);
	// An inefficiency here, but this allow us to discard dead
	// traces and problem data from ensemble as we assemble the
	// new one.
	result->member.clear();
	result->member.reserve(raw.member.size());
	// Load arrivals from database.  List returned is index into raw of
	// data with valid arrivals loaded
	list<int> data_with_arrivals=LoadArrivalTimes<ThreeComponentSeismogram>
					(raw.member,dbh,arrival_keyword);
	list<int>::iterator index;
	for(index=data_with_arrivals.begin();index!=data_with_arrivals.end();++index)
	{
		ThreeComponentSeismogram d=raw.member[*index];
		if(d.live)
		{
		try {
cout << d.get_string("sta")<<" has arrival time ="
	<<strtime(d.get_double("arrival.time"))<<endl;
			d.rotate_to_standard();	
			// partial clone used to hold result
			ThreeComponentSeismogram d3c(d);  
			x1=auto_ptr<TimeSeries>(ExtractComponent(d,0));
			x2=auto_ptr<TimeSeries>(ExtractComponent(d,1));
			x3=auto_ptr<TimeSeries>(ExtractComponent(d,2));
			// resample if necessary.  Using auto_ptr to avoid temporary pointer
			// and as good practice to avoid memory leaks
			if( (abs( (d.dt)-target_dt)/target_dt) > samprate_tolerance)
			{
				*x1=ResampleTimeSeries(*x1,rdef,target_dt,false);
				*x2=ResampleTimeSeries(*x2,rdef,target_dt,false);
				*x3=ResampleTimeSeries(*x3,rdef,target_dt,false);
			}
			// This procedure returns an auto_ptr.  An inconsistency in
			// SEISPP due to evolutionary development
			x1=ArrivalTimeReference(*x1,arrival_keyword,processing_window);
			x2=ArrivalTimeReference(*x2,arrival_keyword,processing_window);
			x3=ArrivalTimeReference(*x3,arrival_keyword,processing_window);
			// safer than using attribute ns in x1
			// assumes all three components are equal length,
			// which is pretty much guaranteed here
			int ns=x1->s.size();
			d3c.ns=ns;
			d3c.dt=x1->dt;
			d3c.t0=x1->t0;
			d3c.tref=x1->tref;
			d3c.u=dmatrix(3,ns);
			// Using blas here for speed
			dcopy(ns,&(x1->s[0]),1,d3c.u.get_address(0,0),3);
			dcopy(ns,&(x2->s[0]),1,d3c.u.get_address(1,0),3);
			dcopy(ns,&(x3->s[0]),1,d3c.u.get_address(2,0),3);
			result->member.push_back(d3c);
		} catch (SeisppError serr)
		{
			// Minor maintenance issue here.  Frozen name is
			// assumed to be in Metadata area.  Avoiding second
			// handler for possible MetadataError intentionally
			string sta=d.get_string("sta");
			cerr << "Problem assembling 3C seismogram for station "
				<< sta <<endl;
			serr.log_error();
			cerr << "Data for this station dropped"<<endl;
		}
		}
			
	}
	return(result);
}
void PostEvid(ThreeComponentEnsemble *d,int evid)
{
	vector<ThreeComponentSeismogram>::iterator dptr;
	for(dptr=d->member.begin();dptr!=d->member.end();++dptr)
		dptr->put("evid",evid);
}
/*! \brief Builds a standard catalog view from a CSS3.0 database.

Much passive array processing is built around the css3.0 schema.
The css3.0 schema has a standard view that defines the definitive
catalog for a network.  This is formed by the join of:
	event->origin->assoc->arrival
and is usually (as here) subsetted to only keep rows with
orid equal to the "preferred origin" (prefor).  This procedure
builds a handle to this database view.

\param dbh is a handle to the parent Datascope database 
	from which the view is to be constructed.  It need
	only reference the correct database. 
\return new handle that refers to the "standard" catalog view
	described above.
*/
DatascopeHandle StandardCatalogView(DatascopeHandle& dbh)
{
	DatascopeHandle result(dbh);
	dbh.lookup("event");
	dbh.natural_join("origin");
	string ss_to_prefor("orid==prefor");
	dbh.subset(ss_to_prefor);
	dbh.natural_join("assoc");
	dbh.natural_join("arrival");
	return(dbh);
}
void SaveResult(DatascopeHandle& dbh,
	auto_ptr<ThreeComponentEnsemble> gather,
		AttributeMap& amo,
			MetadataList& mdlo)
{
	const string table("wfprocess");
	vector<ThreeComponentSeismogram>::iterator d;
	for(d=gather->member.begin();d!=gather->member.end();++d)
	{
		try {
			dbsave(*d,dbh.db,table,mdlo,amo);
		} catch (SeisppError serr) {
			string sta=d->get_string("sta");
			cerr << "Error saving station "<<sta<<endl;
			serr.log_error();
		}
	}
}
/*! \brief Multiple stage TimeInvariantFilter operator.

Sometimes it is necessary to apply a series of stages of
filtering to equalize different data sets or just to 
simply the process of defining a multiple filter chain.
This object simplifies that process by allowing the definition
of a chain of filters and a set of supplied methods to apply
these filters to data.
*/
class MultiStageFilter
{
public:
	/*! \brief Default constructors.  

	Loads a null filter definition.  That is the default is
	a do nothing operator. */
	MultiStageFilter();
	/*! \brief Construct the filter definitions from a string.

	This is currently the prime constructor.  A string is parsed
	into tokens that describe a series of filters that will be
	chained together.  The strings that define each individual
	filter type are parsed into blocks determined by the separator
	argument.  The parsed strings are currently used to construct
	a set of TimeInvariantFilter processing objects.
	An example helps explain how this would be used.  If we
	passed "DEMEAN; BW 0.5 5 2 5" and define ";" as the separator
	this would yield a two stage filter:  demean followed by a 
	0.5 to 2 Hz bandpass filter.

	\param filterlist contains the set of filter recipes to use.
	\param separator is the string used as a separator between
		the recipes for each filter description.
	*/
	MultiStageFilter(string filterlist,string separator);
	template <class T> void  apply(T& d);
private:
	list<TimeInvariantFilter> stages;
};
MultiStageFilter::MultiStageFilter()
{
	TimeInvariantFilter f(string("none"));
	stages.push_back(f);
}
MultiStageFilter::MultiStageFilter(string filterlist,string separator)
{
    try {
	const string white(" \t\n");
	int current=0,end_current;
	string stmp;
	string filterparam;
	// Strip any leading white space
	stmp=filterlist;
	if((current=stmp.find_first_not_of(white,0)) != 0)
	{
		stmp.erase(0,current);
		current=0;
	}
	int endstmp=stmp.size();
	do {
		end_current=stmp.find_first_of(separator,current);
		if(end_current<0)
		{
			filterparam.assign(stmp,current,endstmp);
			stages.push_back(TimeInvariantFilter(filterparam));
			break;
		}			
		filterparam.assign(stmp,current,end_current-current);
		stages.push_back(TimeInvariantFilter(filterparam));
		current=stmp.find_first_not_of(separator,end_current);
		current=stmp.find_first_not_of(white,current);
	}
	while(current<endstmp && current>=0);
    } catch (...) {throw;};
}
		
	
	
/* For now this only works on ensembles.  */
template <class T> void MultiStageFilter::apply(T& d)
{
    try {
	list<TimeInvariantFilter>::iterator filt;
	for(filt=stages.begin();filt!=stages.end();++filt)
		FilterEnsemble(d,*filt);
	/* Note this template could be made to work with TimeSeries
	 or ThreeComponentSeismogram components individually if 
	we used a typeid check on T and then used this line
	for that type of beast:  filt->apply(d);
	*/
    } catch (...) {throw;};
}
void ApplyFST(ThreeComponentEnsemble& e,Hypocenter& hypo)
{
	double vp0(6.0),vs0(3.5);
	for(int i=0;i<e.member.size();++i)
	{
		double lat,lon,elev;
		lat=e.member[i].get_double("lat");
		lon=e.member[i].get_double("lon");
		elev=e.member[i].get_double("elev");
		SlownessVector u=hypo.pslow(lat,lon,elev);
		e.member[i].free_surface_transformation(u,vp0,vs0);
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
	MetadataList mdtrace=pfget_mdlist(pf,"station_mdlist");
	MetadataList mdlo=pfget_mdlist(pf,"output_mdlist");
	try {
		// Get set of control parameters
		string netname=control.get_string("netname");
		string phase=control.get_string("phase");
		double ts,te;
		ts=control.get_double("data_window_start");
		te=control.get_double("data_window_end");
		TimeWindow datatwin(ts,te);
		double tpad=control.get_double("data_time_pad");
		ts=control.get_double("processing_window_start");
		te=control.get_double("processing_window_end");
		TimeWindow processing_twin(ts,te);
		StationChannelMap stachanmap(pf);
		string schemain=control.get_string("InputAttributeMap");
		string schemaout=control.get_string("OutputAttributeMap");
		double target_dt=control.get_double("target_sample_interval");
		ResamplingDefinitions rdef(pf);
		// First get all the database components assembled.
		// There are three here:  input data, output data, and beam data
		DatascopeHandle dbh(dbin,false);
		/* Build a match handle for arrivals using the standard catalog
		view of the join of event->origin->assoc->arrival */
		AttributeMap am(schemain);  
		AttributeMap amo(schemaout);  
		DatascopeHandle dbcatalog=StandardCatalogView(dbh);
		list<string> matchkeys;
		matchkeys.push_back(string("sta"));
		matchkeys.push_back(string("evid"));
		DatascopeMatchHandle dbhm(dbcatalog,string(""),matchkeys,am);
		/*Now construct a handle to get beam traces */
		DatascopeHandle dbhbeam(dbh);
		if(beamdb!=dbin)
			dbhbeam=DatascopeHandle(beamdb,true);
		DatascopeHandle dbhbv(dbhbeam.db,pf,
			string("dbprocess_commands_for_beams"));
		DatascopeHandle dbho(dbh);
		dbho.lookup(string("wfprocess"));
		// subset the beam database view if required
		if(subsetbeam)
		{
			dbhbv.subset(beam_subset_expression);
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
		// End section on prep for beam data
		//
		// Now launch matlab and sets up the communication channels
		MatlabProcessor mp(stdout);
		// 3C ensemble are loaded into Matlab as three matrices with these names
		string chans[3]={"x1","x2","x3"};
		// These are pointers to dynamic objects used in the loop
		// below.  
		ThreeComponentEnsemble *rawdata=NULL, *regular_gather=NULL,
				*decondata=NULL;
		SeismicArray *stations=NULL;
		/* Assorted variables needed in the loop below */
		double lat,lon,depth,otime;
		int evid;
		int record;
		string filter_param;
		for(record=0;record<dbhbv.number_tuples();++record,++dbhbv)
		{

			// Read the beam trace and load it as beam
			TimeSeries beam(dynamic_cast<DatabaseHandle&>(dbhbv),
				mdbeam,am);
			mp.load(beam,string("beam"));
// DEBUG
mp.process(string("plot(beam);"));
			lat=beam.get_double(string("origin.lat"));
			lon=beam.get_double(string("origin.lon"));
			depth=beam.get_double(string("origin.depth"));
			otime=beam.get_double(string("origin.time"));
			evid=beam.get_int(string("evid"));
			// origin coordinates are degrees in the db,
			// but must be radians for hypocenter object
			lat=rad(lat);
			lon=rad(lon);
			Hypocenter hypo(lat,lon,depth,otime,
				string("tttaup"),string("iasp91"));
			// For now we'll always filter the data consistent
			// with the beam.
			try {
				filter_param=beam.get_string("xcorbeam.filter");
			} catch (MetadataGetError mdge) {
				cerr << "Warning:  evid="<<evid<<endl
					<<"filter attribute not defined"
					<<" using default of DEMEAN"<<endl;
				filter_param=string("DEMEAN");
			}
			MultiStageFilter filt(filter_param,string(";"));
			// On the first record we need to load the station
			// geometry object
			if(record==0)
			{
				stations = new SeismicArray(
					dynamic_cast<DatabaseHandle&>(dbh),
					hypo.time,netname);
			}
			else
			{
				TimeWindow test(hypo.time,hypo.time+2000.0);
				if(!stations->GeometryIsValid(test))
				{
					delete stations;
					stations = new SeismicArray(
					dynamic_cast<DatabaseHandle&>(dbh),
					hypo.time,netname);
				}
			}
			// Read the raw data using the time window based constructor
			rawdata=array_get_data(*stations,hypo,
				phase,datatwin,tpad,dynamic_cast<DatabaseHandle&>(dbh),
				stachanmap,mdens,mdtrace,am);
			// Filter the raw data.  
			filt.apply(*rawdata);
			PostEvid(rawdata,evid);

			regular_gather=BuildRegularGather(*rawdata, dbhm,rdef,target_dt,
					processing_twin);
			delete rawdata;
			ApplyFST(*regular_gather,hypo);
			// regular gather is now assumed to contain data with
			// a common start time.  irregular start time is allowed
			// with matlab interface
			mp.load(*regular_gather,chans);
//DEBUG
			char savstr[100];
			sprintf(savstr,"save %d_event",evid);
			mp.process(savstr);
// Uncomment when beamdecon script is working
//			mp.process(string("beamdecon"));
/*
			auto_ptr<ThreeComponentEnsemble> newgath=MatlabRetrieve(mp,
					chans,*regular_gather);
*/
			delete regular_gather;
//DEBUG
// uncomment when working
//			SaveResult(dbho,newgath,amo,mdlo);
			delete decondata;
		}
	}
	catch (SeisppError serr)
	{
		serr.log_error();
	}
	catch (MetadataGetError mdge)
	{
		mdge.log_error();
	}
	catch (MetadataError mde)
	{
		mde.log_error();
	}
	catch (...)
	{
		cerr << "Something threw an unexpected exception"<<endl;
	}
}

