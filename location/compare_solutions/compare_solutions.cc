#include <string>
#include <ostream>
#include "db.h"
#include "pf.h"
#include "elog.h"
using namespace std;
// This section tests a concept of a Database Handle for a db
// Primary purpose of this code is just a simple comparison of two
// catalogs for the genloc paper
//
class Database_Handle
{
public:
	Database_Handle(string dbname, string pfname, string tag);
	~Database_Handle();
	double get_double(string);
	int get_int(string);
	string get_string(string);
	// should have an interator function that could also be 
	// implemented as a ++ operator.  For now use this
	void next_row(){++db.record;};
	int size();
	void rewind();
private:
	Dbptr db;
};
Database_Handle::Database_Handle(string dbname,string pfname, string tag)
{
	if(dbopen(const_cast<char*>(dbname.c_str()),"r+",&db))
		throw string("Failure in dbopen");
	// Probably needs to be an auto_ptr so pffree can be called
	// Maybe should use a Metadata object
	Pf *pf;
	if(pfread(const_cast<char*>(pfname.c_str()),&pf))
		throw string("Failure in pfread");
	Tbl *process_list;
	process_list = pfget_tbl(pf,const_cast<char*>(tag.c_str()));
	if(process_list==NULL)
		throw string("Tag name = "+tag+" not in parameter file");
	db = dbprocess(db,process_list,0);
	if(db.record == dbINVALID)
		throw string("dbprocess failed");
	pffree(pf);
}
Database_Handle::~Database_Handle()
{
	dbclose(db);
}
double Database_Handle::get_double(string name)
{
	double val;
	if(dbgetv(db,0,name.c_str(),&val,0))
		throw string("dbgetv error extracting attribute "+name);
	return(val);
}
int Database_Handle::get_int(string name)
{
	int val;
	if(dbgetv(db,0,name.c_str(),&val,0))
		throw string("dbgetv error extracting attribute "+name);
	return(val);
}
string Database_Handle::get_string(string name)
{
	char val[128];
	if(dbgetv(db,0,name.c_str(),&val,0))
		throw string("dbgetv error extracting attribute "+name);
	return(string(val));
}
int Database_Handle::size()
{
	int ierr;
	int nrec;
	ierr = dbquery(db,dbRECORD_COUNT,&nrec);
	if(ierr<0)
		throw string("dbquery of record count failed");
	return(nrec);
}
void Database_Handle::rewind()
{
	db.record=0;
}
	
int main(int argc, char **argv)
{
	int i;

	if(argc<2)
	{
		cerr << "Usage:  compare_solutions db" << endl;
		exit(-1);
	}
	try
	{
		Database_Handle dbh(string(argv[1]),
			string(argv[0]),string("dbprocess_table"));
		double lat1,lat2,lon1,lon2,z1,z2,t1,t2;
		string alg1,alg2;
		int evid1,evid2;
		double sdobs1,sdobs2;
		cerr << "database has " << dbh.size() << "rows";

		dbh.rewind();
		for(i=0;i<dbh.size();++i,dbh.next_row())
		{
			lat1=dbh.get_double("lat");
			lon1=dbh.get_double("lon");
			z1=dbh.get_double("depth");
			t1=dbh.get_double("time");
			alg1=dbh.get_string("algorithm");
			evid1=dbh.get_int("evid");
			sdobs1=dbh.get_double("sdobs");
			++i;
			dbh.next_row();
			lat2=dbh.get_double("lat");
			lon2=dbh.get_double("lon");
			z2=dbh.get_double("depth");
			t2=dbh.get_double("time");
			alg2=dbh.get_string("algorithm");
			evid2=dbh.get_int("evid");
			sdobs2=dbh.get_double("sdobs");
			if(evid1!=evid2)
			{
				cerr << "Data irregularity at row "<<i<<endl
					<< "evid mismatch between "
					<<evid1 << " and " << evid2 <<endl;
				exit(-2);
			}
			double ssq;
			ssq = 111.0*(lat1-lat2)*(lat1-lat2);
			ssq += 90.9*(lon1-lon2)*(lon1-lon2);
			ssq += (z1-z2)*(z1-z2);
			ssq += 4.5*(t1-t2)*(t1-t2);
			cout << lat1-lat2 << " "
				<< lon1-lon2 << " "
				<< z1-z2 << " "
				<< t1-t2 << " "
				<< sqrt(ssq) << " "
				<< sdobs1-sdobs2 << " "
				<< alg1 <<"-"<<alg2
				<< endl;
		}
	}
	catch (string serr)
	{
		clear_register(1);
		cerr << "Abort message: " << serr << endl;
	}
}
