#include "dbpp.h"
#include "VectorStatistics.h"
using namespace std;
using namespace SEISPP;
void usage()
{
cerr << "xcoravg dbarr db2 [-v -pf pffile]"<<endl
	<<" dbarr - db with arrival and first order dbxcor xsaa table"<<endl
	<<" db2 - db with stacked beams from second-order dbxcor run"<<endl;
exit(-1);
}
bool SEISPP::SEISPP_verbose(false);
int main(int argc, char **argv)
{
	if(argc<3) usage();
	string dbarr(argv[1]);
	string db2(argv[2]);
	string pffile("xcoravg");
	int i;
	for(i=3;i<argc;++i)
	{
		string sarg(argv[i]);
		if(sarg=="-v")
			SEISPP_verbose=true;
                else if(sarg=="-pf")
                {
                        ++i;
                        if(i>=argc)usage();
                        pffile=string(argv[i]);
                }
		else
			usage();
	}
        Pf *pf;
        if(pfread((char *)(pffile.c_str()),&pf))
        {
                cerr << "pfread failed for file = "<<pffile<<endl;
                usage();
        }
	try{
		DatascopeHandle dbh2(db2,false);
		dbh2.lookup("xsaa");
		if(dbh2.number_tuples()<=0)
		{
			cerr << "Fatal error:  database table "
				<< db2 <<".xsaa" << " is empty"<<endl;
			usage();
		}
		if(SEISPP_verbose) 
			cout << "Double beam xsaa table has "<<dbh2.number_tuples()
			<<" rows"<<endl;
		DatascopeHandle dbh(dbarr,true);
		DatascopeHandle dbhxcat(dbh);
		dbhxcat.lookup("xsaa");
		list<string> xcatskeys,xcatgrpkeys;
		xcatskeys.push_back("phase");
		xcatskeys.push_back("sta");
		xcatskeys.push_back("evid");
		xcatgrpkeys=xcatskeys;
		xcatskeys.push_back("gridid");
		dbhxcat.sort(xcatskeys);
		if(SEISPP_verbose) 
			cout << "Working view built from xsaa number of rows="
				<< dbhxcat.number_tuples()<<endl;
		dbhxcat.group(xcatgrpkeys);
		if(SEISPP_verbose) 
			cout << "Number of phase:sta:evid groups="
				<< dbhxcat.number_tuples()<<endl;
		/* Build a match handle against dbh2. */
		AttributeMap am("css3.0");
		DatascopeMatchHandle db2mh(dbh2,string("xsaa"),
			xcatskeys,am);
		/* Now build another match handle on arrival for db1*/
		list<string> matchkeys;
		matchkeys.push_back("sta");
		/* This is a bit problematic.  arrival uses iphase and we have
		put phase as an attribute in xsaa */
		matchkeys.push_back("iphase");
		matchkeys.push_back("evid");
		/* We need to use a join of arrival and assoc because 
		arid is not in the xsaa table.  This would be simpler
		if that were the case, but elected to do this to avoid
		maintaining arid as a unique key. This is too often 
		abused by me and others and is subject to errors. 
		Dark side is it complicates this code. */
		dbh.lookup("event");
		dbh.natural_join("origin");
		dbh.subset(string("orid==prefor"));
		dbh.natural_join("assoc");
		dbh.natural_join("arrival");
		if(dbh.number_tuples()<=0) 
		{
			cerr << "event->origin->assoc->origin view is empty."
				<<"  Verify database = "<< dbarr<<endl;
			usage();
		}
		/* Null string is a signal to this constructor to use
		dbh.db as the table to match against. */
		DatascopeMatchHandle arrivalhandle(dbh,string(""),
			matchkeys,am);
		/* Now loop through the group view computing statistics for
		each arrival.  The median of the group is used to set the
		new arrival time estimate */
		dbhxcat.rewind();
		int i;
		int ierr;
		/* This holds the list of matching record numbers used to find
		the right row in arrival */
		list<int> matchlist;
		vector<double> arrival_times;
		cout << "evid sta phase db_arrival_time median_dt mean_dt interquartile mad range ndgf " <<endl;
		/* Hold keys extracted with dbgetv in loop below for matching to db2 */
		int gridid, evid;
		char sta[20],phase[10];
		Metadata matchxsaa;
		for(i=0;i<dbhxcat.number_tuples();++i,++dbhxcat)
		{
			Dbptr dbarrival;  // holds Datascope Dbptr fetched from view below
			Dbptr dbview;  // used for looping through group view 
			DBBundle grp=dbhxcat.get_range();
			/* correct because end is +1 as per dbget_range */
			int narr=grp.end_record-grp.start_record;
			dbview=grp.parent;
			arrival_times.clear();
			for(dbview.record=grp.start_record;dbview.record<grp.end_record;++dbview.record)
			{
				double atime;
				ierr=dbgetv(dbview,0,"time",&atime,
					"sta",sta,"phase",phase,"evid",&evid,"gridid",&gridid,0);
				if(ierr==dbINVALID)
				{
					cerr << "Warning:  dbgetv error.  Skipping record number "
						<< dbview.record<<endl;
				}
				/* if there is only one member of the group, there cannot
				be a match in db2.  Thus, we just post the time in that
				case an skip the more complicated logic in the else block
				that follows */
				else if( (grp.end_record-grp.start_record)==1)
				{
					arrival_times.push_back(atime);
				}
				else
				{
					matchxsaa.put("sta",sta);
					matchxsaa.put("phase",phase);
					matchxsaa.put("evid",evid);
					matchxsaa.put("gridid",gridid);
					list<int> recs=db2mh.find(matchxsaa,false);
					int nrecs=recs.size();
					if(nrecs>1)
					{
						cerr << "Warning:  multiple matches for db2 with "
							<< "key (sta:phase:evid:gridid)="
							<< sta <<":"
							<< phase <<":"
							<< evid <<":"
							<< gridid <<endl
						     << "Using first record found"<<endl;
					}
					if(nrecs!=0)
					{
						double adt;
						list<int>::iterator irptr=recs.begin();
						dbh2.db.record=*irptr;
						dbgetv(dbh2.db,0,"time",&adt,0);
						atime += adt;
						arrival_times.push_back(atime);
					}
				}
			}
			narr=arrival_times.size();
			if(narr<=0)
			{
				cerr << "Warning:  no arrivals retrieved for row "
					<< i << " of sorted and grouped xsaa table"<<endl
					<< "Likely cause is skipped gather in dbxcor pass 2"
					<<endl;;
			}
			else
			{
				string stastr=dbhxcat.get_string("sta");
				string phasestr=dbhxcat.get_string("phase");
				evid=dbhxcat.get_int("evid");
				Metadata finder;
				finder.put("sta",stastr);
				finder.put("evid",evid);
				finder.put("iphase",phasestr);
				matchlist=arrivalhandle.find(finder,false);
				if(matchlist.size()<=0)
					cerr << "dbmatches failed for evid:sta:phase="
						<<evid<<" : "
						<<stastr<<" : "
						<<phasestr<<endl;
				else
				{
					if(matchlist.size()>1)
						cerr << "Warning:  evid:sta:phase="
                                                <<evid<<" : "
                                                <<sta<<" : "
                                                <<phase<<endl
						<<"dbmatches found "<<matchlist.size()
						<< " matching rows.  Using first row of list."
						<<endl;
					dbview=dbh.db;
					dbview.record=*(matchlist.begin());
					ierr=dbgetv(dbview,0,"arrival",&dbarrival,0);
					if(ierr==dbINVALID)
						cerr << "Warning:  evid:sta:phase="
                                                <<evid<<" : "
                                                <<sta<<" : "
                                                <<phase<<endl
						<< "Error fetching Dbptr for arrival table.  "
						<< "Results for this arrival will not be in output database"
						<<endl;
					else
					{
					    double oldatime;
					    dbgetv(dbarrival,0,"time",&oldatime,0);
					    if(narr>1)
					    {
						VectorStatistics<double> atstat(arrival_times);
						double atime=atstat.median();
						
						cout << evid <<" "
							<< stastr<<" "
							<< phasestr <<" "
							<< strtime(oldatime) << " "
							<< oldatime-atime << " "
							<< oldatime-atstat.mean() <<" "
							<< atstat.interquartile() <<" "
							<< atstat.mad(atime) <<" "
							<< atstat.range() <<" "
							<< narr-1 <<endl;
					    }
					    else
					    {
						cout << evid <<" "
							<< stastr <<" "
							<< phasestr <<" "
							<< strtime(oldatime) << " "
							<< oldatime-arrival_times[0] << " "
							<< oldatime-arrival_times[0] << " "
							<< "0.0 "
							<< "0.0 "
							<< "0.0 "
							<< "0" <<endl;
					    }
					}
				}
			}
		}
	}
	catch(SeisppError serr)
	{
		serr.log_error();
	}
}

