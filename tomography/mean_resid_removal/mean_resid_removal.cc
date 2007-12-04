#include <iostream>
#include <string>
#include <vector>
#include <list>
#include "seispp.h"
#include "dbpp.h"
using namespace std;
using namespace SEISPP;
int main(int argc, char **argv)
{
	if(argc!=2) 
	{
		cerr << "Usage:"<<endl<<"    mean_residual_removal db" <<endl;
		exit(-1);
	}
	string dbname(argv[1]);
	cout << "Teleseismic event mean residual removal program" <<endl
		<< "WARNING:  this program will permanently alter assoc timeres attributes"
			" in database = "<<dbname<<endl
		<< "A backup copy of your database should be made before you run this program."<<endl
		<<"Are you sure you want to continue (type y to continue, anything will exit):";
	string ques;
	cin >> ques;
	if(ques!="y") exit(-1);
	cout << "Enter phase name in arrival table to process:";
	string phase;
	cin >> phase;
	char ss_phase[128];
	sprintf(ss_phase,"arrival.iphase=~/%s/",phase.c_str());
	try {
		DatascopeHandle dbh(dbname,false);
		dbh.lookup("event");
		list<string> sortkeys;
		sortkeys.push_back("evid");
		dbh.sort(sortkeys);
		dbh.natural_join("origin");
		string ss_to_prefor("orid==prefor");
		dbh.subset(ss_to_prefor);
		dbh.natural_join("assoc");
		dbh.natural_join("arrival");
		dbh.subset(string(ss_phase));
		list<string> group_keys;
		group_keys.push_back("evid");
		dbh.group(group_keys);
		dbh.rewind();
		list<int> assocrecs;
		list<int>::iterator assocptr;
		vector<double> residuals;
		int evid;
		cout << "Found "<<dbh.number_tuples()<<" events to process"
			<<endl <<endl;
		cout << "event_id    narrivals mean_removed"<<endl;
		int i,j;
		for(i=0;i<dbh.number_tuples();++i,++dbh)
		{
			DBBundle bundle=dbh.get_range();
			Dbptr dbp=bundle.parent;
			Dbptr dbassoc;
			evid=dbh.get_int("evid");
			assocrecs.clear();
			residuals.clear();
			double timeres;
			for(dbp.record=bundle.start_record;
				dbp.record<=bundle.end_record;++dbp.record)
			{
				int iret;
				iret=dbgetv(dbp,0,"assoc",&dbassoc,
					"assoc.timeres",&timeres,0);
				if(iret==dbINVALID)
				{
					cerr << "dbgetv error at row "
						<< dbp.record<<endl;
				}
				else
				{
					assocrecs.push_back(dbassoc.record);
					residuals.push_back(timeres);
				}
			}
			double tmean;
			int ntimes=residuals.size();
			for(j=0,tmean=0.0;j<ntimes;++j)
			{
				tmean+=residuals[i];
			}
			tmean /= static_cast<double>(ntimes);
			cout << evid <<"  "
				<<ntimes<<"  "
				<<tmean<<endl;
			/* assume dbassoc table value points at assoc so we can reuse it 
			without calling lookup. */
			for(j=0,assocptr=assocrecs.begin();j<ntimes;++j,++assocptr)
			{
				dbassoc.record=(*assocptr);
				timeres=residuals[j]-tmean;
				dbputv(dbassoc,0,"timeres",timeres,0);
			}
		}
	}
	catch (SeisppError serr)
	{
		serr.log_error();
		exit(-1);
	}
}
				

