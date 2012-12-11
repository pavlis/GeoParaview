/* This is a copy of a generic algorithm taken from the source code
   file for extract_events.  That program is in antelope contrib.   
   A copy is here for use in this incomplete program.   Long term this
   file should be placed in libseispp, but for now it is baggage left here.
   
GLP:  Dec 10, 2012
*/
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
using namespace std;
using namespace SEISPP;
template <class T> list<int> LoadArrivalTimes(std::vector<T>& dat,
                DatascopeMatchHandle& dbh,
		    const string keyword)
{
        std::vector<T>::iterator d;
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
