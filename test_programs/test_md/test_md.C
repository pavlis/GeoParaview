#include "stock.h"
#include "pf.h"
#include "db.h"
#include "metadata.h"
#include "seispp.h"
using namespace SEISPP;
// test routine for metadata functions
int main(int argc, char **argv)
{
	string pfname="test_md";;
	string amname="Attribute_Map";
	int ierr;

	Pf *pf;
	
	//if(pfread(const_cast<char *>(pfname.c_str),&pf))
	if(pfread("test_md",&pf))
	{
		cerr << "Pfread error"<<endl;
		exit(-1);
	}
	try {
		cout << "Trying Attribute_map pf constructor"<<endl;
		SEISPP::Attribute_map am(pf,amname);
		cout << "Trying copy constructor and = operator"<<endl;
		SEISPP::Attribute_map am2(am);
		am=am2;
		cout << "Trying pfget_mdlist function with tag mdlist" << endl;
		Metadata_list mdl=pfget_mdlist(pf,"mdlist");
		cout << "Trying to build Metadata objects" << endl;
		cout << "Trying default constructor" << endl;
		Metadata mdplain;
		cout << mdplain;
		cout << "Trying string constructor with inline data" << endl;
		string smdtest("x 4.5\ny 8.0\ni 27\nMYPAR 2.5\ncpar xyz\n");
		Metadata mds(smdtest);
		cout << mds;
		cout <<"Trying db constructor in testdb"<<endl;
		Dbptr db;
		char *dbname="testdb";
		ierr=dbopen(dbname,"r+",&db);
		db=dblookup(db,0,"wfdisc",0,0);
		db.record = 0;
		Metadata *mddb = new Metadata(db,mdl,am2);
		cout << *mddb;
	}
	catch (Metadata_error mess)
	{
		mess.log_error();
		exit(-1);
	}
	
}
