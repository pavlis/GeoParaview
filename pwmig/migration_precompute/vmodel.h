#include <string>
using namespace std;

class Velocity_Model_1d_error
{
public:
	string message;
	virtual void log_error(){cerr<<"Velocity_Model_1d object error"<<endl;};
};

class Velocity_Model_1d_dberror : public Velocity_Model_1d_error
{
public:
	string name;
	Velocity_Model_1d_dberror(string modname,string mess){
		name=modname;  message=mess;};
	virtual void log_error(){
		cerr<<"Database error accessing velocity model mod1d table "<<name<<endl;
		cerr<<message;
	};
};

class Velocity_Model_1d
{
public:
	int nlayers;
	double *z,*v,*grad;
	Velocity_Model_1d(){z=NULL;v=NULL;grad=NULL;nlayers=0;};
	Velocity_Model_1d(int n){nlayers=n;
		z=new double[nlayers]; 
		v=new double[nlayers];
		grad=new double[nlayers];};
	Velocity_Model_1d(Dbptr db,string name, string property)
		throw(Velocity_Model_1d_dberror);
	~Velocity_Model_1d()
	{if(z!=NULL)delete[]z; if(v!=NULL)delete[]v; if(grad!=NULL)delete[]grad;};
	double getv(double zin);
};

