#include "stock.h"
#include "pf.h"
#include "db.h"
#include "metadata.h"
#include "seispp.h"
using namespace SEISPP;
// test routine for metadata functions
bool SEISPP::SEISPP_verbose(true);
int main(int argc, char **argv)
{
	//string pfname("test_md");
	char *pfname="test_md";
	int ierr;

	Pf *pf;
	
	//if(pfread(const_cast<char *>(pfname.c_str),&pf))
	if(pfread(pfname,&pf))
	{
		cerr << "Pfread error"<<endl;
		exit(-1);
	}
	try {
		cout << "Trying nested pf constructor" << endl;
		Metadata mdnest(pf,string("test_nested_tag"));
		cout<<"test_double="<<mdnest.get_double("test_double")<<endl;
		cout<<"test_int="<<mdnest.get_int("test_int")<<endl;
		cout<<"test_string="<<mdnest.get_string("test_string")<<endl;
                cout << "Trying file constructor on test_md.pf"<<endl;
                Metadata mdfc(string("test_md.pf"));
                cout << mdfc;
		cout << "Trying AttributeMap default constructor"<<endl;
		AttributeMap amd;
		cout << "Default map has "<< amd.attributes.size()
				<< "members" << endl;
		cout << "Trying AttributeMap pf constructor"<<endl;
		SEISPP::AttributeMap am("css3.0");
		cout << "Trying copy constructor and = operator"<<endl;
		SEISPP::AttributeMap am2(am);
		am=am2;
		cout << "Trying pfget_mdlist function with tag mdlist" << endl;
		MetadataList mdl=pfget_mdlist(pf,"mdlist");
		cout << "Trying to copy mdlist" << endl;
		MetadataList mdl2;
		mdl2=mdl;
		cout << "Original and copy of mdlist read" << endl;
		MetadataList::iterator mdi,mdi2;
		for(mdi=mdl.begin(),mdi2=mdl2.begin();mdi!=mdl.end();++mdi,++mdi2)
		{
			cout << mdi->tag 
			<< "," 
			<< mdi2->tag 
			<< "," 
			<<mdi->mdt
			<< "," 
			<<mdi2->mdt<<endl;
		}
		cout << "Testing delete hypothesis on list object " << endl;
		MetadataList *mdl3=new MetadataList(mdl);
		delete mdl3;
		for(mdi=mdl.begin();mdi!=mdl.end();++mdi)
		{
			cout << mdi->tag 
			<< "," 
			<<mdi->mdt<<endl;
		}

		cout << "Trying to build Metadata objects" << endl;
		cout << "Trying default constructor" << endl;
		Metadata mdplain;
		cout << mdplain;
		cout << "Trying string constructor with inline data" << endl;
		string smdtest("x 4.5\ny 8.0\ni 27\nMYPAR 2.5\ncpar xyz\n");
		Metadata mds(smdtest,"string");
		cout << mds;
		cout <<"Trying db constructor in testdb"<<endl;
		Dbptr db;
		string dbname("testdb");
		DatascopeHandle dbh(dbname,pfname,string("dbprocess_list"),false);
		dbh.rewind();
		Metadata *mddb = new Metadata(dynamic_cast<DatabaseHandle&>(dbh),
					mdl,am2);
		cout << *mddb;

	}
        catch (SeisppError& sess)
	{
		sess.log_error();
                exit(-1);
	}
	
}
