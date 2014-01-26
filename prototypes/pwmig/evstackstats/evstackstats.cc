#include <iostream>
#include <sstream>
#include "dbpp.h"
#include "VectorStatistics.h"
using namespace std;
using namespace SEISPP;

void usage()
{
	cout << "evstackstats db [-minfold n]" <<endl
		<< "Default minfold=1"<<endl;
	exit(-1);
}

int main(int argc, char **argv)
{
	if(argc<2) usage();
	string dbname(argv[1]);
	int minfold(1);
	int i;
	for(i=0;i<argc;++i)
	{
		string sarg(argv[i]);
		if(sarg=="-minfold")
		{
			++i;
			if(i>=argc) usage();
			minfold=atoi(argv[i]);
			if(minfold<1) cerr << "Illegal -minfold value="
					<< minfold<<endl<<"Must be > 1"<<endl;
			usage();
		}
		else
			usage();
	}
	try {
		DatascopeHandle dbh(dbname,true);
		dbh.lookup("stackstats");
		ostringstream ss;
		ss << "fold>"<<minfold;
		dbh.subset(ss.str());
		list<string> keys;
		keys.push_back("fold");
		dbh.sort(keys);
		dbh.group(keys);
		int nfold=dbh.number_tuples();
		cout << "semblance( min, s1_4, median, s3_4, max)"
			<< " coherence (min, s1_4, median, s3_4, max)"<<endl;
		vector<double> coherence;
		vector<double> semblance;
		for(i=0,dbh.rewind();i<nfold;++i,++dbh);
		{
			DBBundle dbb=dbh.get_range();
			Dbptr db=dbh.db;
			coherence.clear();
			semblance.clear();
			int fold=dbh.get_int("fold");
			for(db.record=dbb.start_record;
				db.record<=dbb.end_record;++db.record)
			{
				double coh,semb;
				if(dbgetv(db,0,"coherence",&coh,
					"semblance",&semb,0)==dbINVALID)
				{
					cerr << "dbgetv failed on record "
						<<db.record<<" skipped this row"
						<<endl;
					continue;
				}
				coherence.push_back(coh);
				semblance.push_back(semb);
			}
			if(coherence.size()<=1) 
			{
				cerr << "Insufficent data to compute statistics"
					<<" at fold="<<fold<<endl;
			}
			else
			{
				VectorStatistics<double> cohstat(coherence);
				VectorStatistics<double> semstat(semblance);
				cout << fold <<" "
					<<cohstat.lower_bound()<< " "
					<< cohstat.q1_4()<< " "
					<< cohstat.median()<< " "
					<< cohstat.q3_4()<< " "
					<< cohstat.upper_bound()<< " "
					<< semstat.lower_bound()<< " "
					<< semstat.q1_4()<< " "
					<< semstat.median()<< " "
					<< semstat.q3_4()<< " "
					<< semstat.upper_bound()<<endl;
			}
		}
	} catch (SeisppError serr)
	{
		serr.log_error();
	}
	catch (MetadataGetError mderr)
	{
		mderr.log_error();
	}
}



