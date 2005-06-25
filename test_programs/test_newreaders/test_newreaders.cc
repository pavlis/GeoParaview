#include <stdio.h>
#include <list>
#include "coords.h"
#include "seispp.h"
using namespace std;
using namespace SEISPP;
void usage(char *prog)
{
	cerr << prog << ":  "
		<< "db < list_of_time_intervals_file"<<endl;
	exit(-1);
}
int main(int argc, char **argv)
{
	TimeSeriesEnsemble *tse;
	if(argc!=2) usage(argv[0]);
	string dbname(argv[1]);
	ios::sync_with_stdio();
	char tstart[80],tend[80];
	Pf *pf;  
	char *pfname="test";
	if(pfread(pfname,&pf))
	{
		cerr << "Error reading pf file test.pf"<<endl;
	}

	try {
		DatascopeHandle dbh(dbname,false); // open read only
		TimeWindow twin;
		while(fscanf(stdin,"%s%s",tstart,tend)!=EOF)
		{
			int i;
			twin.start=str2epoch(tstart);
			twin.end=str2epoch(tend);
			cout << "Attempting read from "<<dbname
				<< "for time interval"
				<< tstart << " to " << tend << endl;

			tse=new TimeSeriesEnsemble(dynamic_cast<DatabaseHandle&>(dbh),
				twin);
			for(i=0;i<tse->member.size();++i)
			{
				cout << "sta="<<tse->member[i].get_string("sta")<<endl;
				cout << "chan="<<tse->member[i].get_string("chan")<<endl;
				cout << "time="<<tse->member[i].get_double("time")<<endl;
				cout << endl;
			}
			delete tse;

		}
	} catch (...){
		cerr << "Test failed.  Something threw an exception."<<endl;
	}
}
