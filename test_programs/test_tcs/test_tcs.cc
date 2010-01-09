#include "seispp.h"
using namespace std;
using namespace SEISPP;
int main(int argc, char **argv)
{
    string pfname("test_tcs");
    string amname="Attribute_Map";
    int i;
    int ierr;

    Pf *pf;
    ios::sync_with_stdio();

    if(pfread(const_cast<char *>(pfname.c_str()),&pf))
    {
           cerr << "Pfread error"<<endl;
           exit(-1);
    }
   try
   {
	int i,j;
	cout << "Testing simple constructors" << endl;
	ThreeComponentSeismogram s1();
	ThreeComponentSeismogram s2(100);
	// Initialize the contents of s2 and apply a rotation 
	// matrix
	cout << "trying rotation" << endl;
	for(i=0;i<3;++i)
		for(j=0;j<100;++j)
		{
			if(i==2) s2.u(i,j)=1.0;
			else
				s2.u(i,j)=0.0;
		}
	Spherical_Coordinate sc;
	sc.theta=0.0;
	sc.theta=M_PI_4;
	s2.rotate(sc);
	cout << "one sample of (0,0,1) rotated by theta=45 deg:"
		<< s2.u(0,0) << ", "
		<< s2.u(1,0) << ", "
		<< s2.u(2,0) << endl;
	cout << "Restoring to standard" << endl;
	s2.rotate_to_standard();
	cout << "One sample of restored (0,0,1): "
		<< s2.u(0,0) << ", "
		<< s2.u(1,0) << ", "
		<< s2.u(2,0) << endl;
	ThreeComponentSeismogram s3(100);
	for(i=0;i<3;++i)
		for(j=0;j<100;++j)
		{
			if(i==0) s3.u(i,j)=1.0;
			else
				s3.u(i,j)=0.0;
		}
	// set some other vectors used in test below
	s3.u(0,1)=1.0;
	s3.u(1,1)=1.0;
	s3.u(2,1)=1.0;

	s3.u(0,2)=1.0;
	s3.u(1,2)=1.0;
	s3.u(2,2)=0.0;

	s3.u(0,3)=0.0;
	s3.u(1,3)=0.0;
	s3.u(2,3)=1.0;

	cout << "Trying unit vector rotation to (1,1,1) with (1,0,0) data"<<endl;
	double nu[3]={sqrt(3.0)/3.0,sqrt(3.0)/3.0,sqrt(3.0)/3.0};
	s3.rotate(nu);
	cout << "one sample of (1,0,0) rotated to (1,1,1)"
		<< s3.u(0,0) << ", "
		<< s3.u(1,0) << ", "
		<< s3.u(2,0) << endl;
	cout << "one sample of (1,1,1) rotated to (1,1,1)"
		<< s3.u(0,1) << ", "
		<< s3.u(1,1) << ", "
		<< s3.u(2,1) << endl;
	cout << "Restoring to standard" << endl;
	s3.rotate_to_standard();
	cout << "One sample of restored (1,0,0): "
		<< s3.u(0,0) << ", "
		<< s3.u(1,0) << ", "
		<< s3.u(2,0) << endl;
	cout << "One sample of restored (1,1,1): "
		<< s3.u(0,1) << ", "
		<< s3.u(1,1) << ", "
		<< s3.u(2,1) << endl;
	cout << "Trying spherical coordinate transformation to rotate phi=45 degrees"
	<< endl;
	sc.phi=M_PI_4;
	sc.theta=0.0;
	s3.rotate(sc);
	cout << "one sample of rotated (1,0,0)"
		<< s3.u(0,0) << ", "
		<< s3.u(1,0) << ", "
		<< s3.u(2,0) << endl;
	cout << "one sample of rotated (1,1,0)"
		<< s3.u(0,2) << ", "
		<< s3.u(1,2) << ", "
		<< s3.u(2,2) << endl;
	cout << "one sample of rotated (0,0,1)"
		<< s3.u(0,3) << ", "
		<< s3.u(1,3) << ", "
		<< s3.u(2,3) << endl;
	s3.rotate_to_standard();
	cout << "Applying nonorthogonal transformation"<<endl;
	double a[3][3];
	a[0][0]=1.0;  a[0][1]=1.0;  a[0][2]=1.0;
	a[1][0]=-1.0;  a[1][1]=1.0;  a[1][2]=1.0;
	a[2][0]=0.0;  a[2][1]=-1.0;  a[2][2]=0.0;
	s3.apply_transformation_matrix(a);
	cout << "one sample of transformed (1,0,0)"
		<< s3.u(0,0) << ", "
		<< s3.u(1,0) << ", "
		<< s3.u(2,0) << endl;
	cout << "one sample of transformed (1,1,0)"
		<< s3.u(0,2) << ", "
		<< s3.u(1,2) << ", "
		<< s3.u(2,2) << endl;
	cout << "one sample of transformed (0,0,1)"
		<< s3.u(0,3) << ", "
		<< s3.u(1,3) << ", "
		<< s3.u(2,3) << endl;
	cout << "Testing back conversion to standard coordinates";
	s3.rotate_to_standard();
	cout << "One sample of restored (1,0,0): "
		<< s3.u(0,0) << ", "
		<< s3.u(1,0) << ", "
		<< s3.u(2,0) << endl;
	cout << "One sample of restored (1,1,1): "
		<< s3.u(0,1) << ", "
		<< s3.u(1,1) << ", "
		<< s3.u(2,1) << endl;
	cout << "One sample of restored (1,1,0): "
		<< s3.u(0,2) << ", "
		<< s3.u(1,2) << ", "
		<< s3.u(2,2) << endl;
	cout << "One sample of restored (0,0,1): "
		<< s3.u(0,3) << ", "
		<< s3.u(1,3) << ", "
		<< s3.u(2,3) << endl;

	cout << "Testing multiple, accumulated transformations" << endl;
	s3.rotate(nu);
	s3.apply_transformation_matrix(a);
	s3.rotate(sc);
	s3.rotate_to_standard();
	cout << "One sample of restored (1,0,0): "
		<< s3.u(0,0) << ", "
		<< s3.u(1,0) << ", "
		<< s3.u(2,0) << endl;
	cout << "One sample of restored (1,1,1): "
		<< s3.u(0,1) << ", "
		<< s3.u(1,1) << ", "
		<< s3.u(2,1) << endl;
	cout << "One sample of restored (1,1,0): "
		<< s3.u(0,2) << ", "
		<< s3.u(1,2) << ", "
		<< s3.u(2,2) << endl;
	cout << "One sample of restored (0,0,1): "
		<< s3.u(0,3) << ", "
		<< s3.u(1,3) << ", "
		<< s3.u(2,3) << endl;

	cout << "Trying Db constructor" << endl;
	int ierr;
	// dbsort and dbgroup have to have been done correctly for
	// this to have a hope of working.
	Dbptr db;
	Pf *pf;
	string tag("test_tcs_view");
	if(dbopen("t3ctest","r+",&db))
	{
		cerr << "dbopen of testdb failed" << endl;
		exit(-1);
	}
	if(pfread(const_cast<char *>(pfname.c_str()),&pf))
	{
		cerr << "pfread of "<< pfname << " failed" << endl;
		exit(-1);
	}
	Datascope_Handle dbh(db,pf,tag);
	list<string>groupkeys;
	groupkeys.push_back(string("sta"));
	cout << "Table number before group = "<<dbh.db.table<<endl;
	dbh.group(groupkeys);
	cout << "Table number after group = "<<dbh.db.table<<endl;
	Attribute_Map am(pf,amname);
	Metadata_list mdl=pfget_mdlist(pf,"metadata_input");
	Metadata_list mdlout=pfget_mdlist(pf,"metadata_output");
	Metadata_list mdlcp=pfget_mdlist(pf,"metadata_copy");
	dbh.db.record=0;
	ThreeComponentSeismogram s4(dynamic_cast<Database_Handle&>(dbh),mdlcp,am);
	cout << "Metadata for 3c seismogram constructed from database"<<endl;
	cout << dynamic_cast<Metadata&>(s4);
	cout << "Transformation matrix of this seismogram"<<endl;
	cout << s4.tmatrix[0][0]   
		<< ", "	<<  s4.tmatrix[0][1] 
		<< ", "	<<  s4.tmatrix[0][2]<<endl;
	cout << s4.tmatrix[1][0] 
		<< ", "	<<  s4.tmatrix[1][1] 
		<< ", "	<<  s4.tmatrix[1][2]<<endl;
	cout << s4.tmatrix[2][0] 
		<< ", "	<<  s4.tmatrix[2][1] 
		<< ", "	<<  s4.tmatrix[2][2]<<endl;
	// output test
	Datascope_Handle dbho("testout",false);
	cout << "Testing dbsave function"<<endl;
	s4.put_metadata("directory",".");
	s4.put_metadata("file","testdata_out.w");
	dbsave(s4,dbho.db,"wfdisc",mdlout,am);
	cout << "Successful write to testout database"<<endl;
	cout << "Trying to read another 3c object from database"<<endl;
	++dbh;
	s4=ThreeComponentSeismogram(dynamic_cast<Database_Handle&>(dbh),mdlcp,am);
	s4.put_metadata("directory",".");
	cout << "Metadata for next 3c trace read"<<endl;
	cout << dynamic_cast<Metadata&>(s4);
	cout << "Transformation matrix of this seismogram"<<endl;
	cout << s4.tmatrix[0][0]   
		<< ", "	<<  s4.tmatrix[0][1] 
		<< ", "	<<  s4.tmatrix[0][2]<<endl;
	cout << s4.tmatrix[1][0] 
		<< ", "	<<  s4.tmatrix[1][1] 
		<< ", "	<<  s4.tmatrix[1][2]<<endl;
	cout << s4.tmatrix[2][0] 
		<< ", "	<<  s4.tmatrix[2][1] 
		<< ", "	<<  s4.tmatrix[2][2]<<endl;
	cout << "Testing dbsave function with chanmap arguments"<<endl;
	vector<string> chanmap;
	chanmap.reserve(3);
	chanmap.push_back("BHE");
	chanmap.push_back("BHN");
	chanmap.push_back("BHZ");
	bool oswitch;
	oswitch=true;
	dbsave(s4,dbho.db,"wfdisc",mdlout,am,chanmap,oswitch);
	cout << "Testing free_surface_transformation algorithm"<<endl;
	s4.rotate_to_standard();
	Slowness_vector uvec;
	uvec.ux=0.17085;  // cos(-20deg)/5.5
	uvec.uy=-0.062185; // sin(-20deg)/5.5
	//s4.free_surface_transformation(uvec,8.0,3.5); // test exception
	s4.free_surface_transformation(uvec,5.0,3.5);
	vector<string>chanmap2;
	chanmap2.push_back("T");
	chanmap2.push_back("R");
	chanmap2.push_back("L");
	dbsave(s4,dbho.db,"wfdisc",mdlout,am,chanmap2,false);
    }
    catch (seispp_dberror sderr)
    {
	sderr.log_error();
    }
    catch (seispp_error serr)
    {
	serr.log_error();
    }
}
