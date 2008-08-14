#include "dbpp.h"
#include "VectorStatistics.h"
#include "MatlabProcessor.h"
using namespace std;
using namespace SEISPP;
void usage()
{
cerr << "xcoravg dbin dbout [-v -pf pffile]"<<endl;
exit(-1);
}
bool SEISPP::SEISPP_verbose(false);
int main(int argc, char **argv)
{
	if(argc<3) usage();
	string dbname(argv[1]);
	string dbout(argv[2]);
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
		DatascopeHandle dbhout(dbout,false);
		dbhout.lookup("arrival");
		DatascopeHandle dbh(dbname,true);
		DatascopeHandle dbhxcat(dbh);
		dbhxcat.lookup("xsaa");
		list<string> xcatskeys,xcatgrpkeys;
		xcatskeys.push_back("phase");
		xcatskeys.push_back("sta");
		xcatskeys.push_back("evid");
		xcatgrpkeys=xcatskeys;
		xcatskeys.push_back("gridid");
		dbhxcat.sort(xcatskeys);
		dbhxcat.natural_join("wfprocess");
		/* Tricky, but this groups by phase:sta */
		dbhxcat.group(xcatgrpkeys);
		/* Build a match handle on arrival */
		AttributeMap am("css3.0");
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
				<<"  Verify database = "<< dbname<<endl;
			usage();
		}
		/* Null string is a signal to this constructor to use
		dbh.db as the table to match against. */
		DatascopeMatchHandle arrivalhandle(dbh,string(""),
			matchkeys,am);
                /* set up for matlab plotting of stacked data for each gather.
                Need to load a metadata list and then start matlab */
                MetadataList mdtrace=pfget_mdlist(pf,"Trace_mdlist");
                MetadataList mdens=pfget_mdlist(pf,"Ensemble_mdlist");
                MatlabProcessor mp(stdout);
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
		cout << "evid sta phase arrival_time atmean interquartile mad range ndgf " <<endl;
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
				ierr=dbgetv(dbview,0,"time",&atime,0);
				if(ierr==dbINVALID)
				{
					cerr << "Warning:  dbgetv error.  Skipping record number "
						<< dbview.record<<endl;
				}
				else
				{
					arrival_times.push_back(atime);
				}
			}
			narr=arrival_times.size();
			if(narr<=0)
			{
				cerr << "Warning:  no arrivals retrieved for row "
					<< i << " of sorted and grouped xsaa table"<<endl;
			}
			else
			{
				TimeSeriesEnsemble tse(dynamic_cast<DatabaseHandle&>(dbhxcat),
					mdens,mdtrace,am);
	cout << "DEBUG:  number of members in this ensemble = "<<tse.member.size()<<endl;
				mp.load(tse,string("beams"));
				mp.process("wigb(beams);");
				string sta=dbhxcat.get_string("sta");
				string phase=dbhxcat.get_string("phase");
				int evid=dbhxcat.get_int("evid");
				Metadata finder;
				finder.put("sta",sta);
				finder.put("evid",evid);
				finder.put("iphase",phase);
				matchlist=arrivalhandle.find(finder,false);
				if(matchlist.size()<=0)
					cerr << "dbmatches failed for evid:sta:phase="
						<<evid<<" : "
						<<sta<<" : "
						<<phase<<endl;
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
					    if(narr>1)
					    {
						VectorStatistics<double> atstat(arrival_times);
						double atime=atstat.median();
						
						cout << evid <<" "
							<< sta <<" "
							<< phase <<" "
							<< strtime(atime) <<" "
							<< strtime(atstat.mean()) <<" "
							<< atstat.interquartile() <<" "
							<< atstat.mad(atime) <<" "
							<< atstat.range() <<" "
							<< narr-1 <<endl;
					    }
					    else
					    {
						cout << evid <<" "
							<< sta <<" "
							<< phase <<" "
							<< strtime(arrival_times[0]) <<" "
							<< strtime(arrival_times[0]) <<" "
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

