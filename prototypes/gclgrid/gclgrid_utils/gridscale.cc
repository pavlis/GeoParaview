#include "dbpp.h"
#include "Metadata.h"
#include "seispp.h"
#include "gclgrid.h"
#include "VelocityModel_1d.h"
using namespace std;
using namespace SEISPP;

void usage()
{
	cerr << "gridscale db grida:fielda * c = gridb:fieldb [-vector -outdir outdir]"<<endl
		<< "where grid and field are gclgrid and gclfield name tags respectively."<<endl
		<< " c is a constant by which fielda is multiplied to produce fieldb"<<endl
		<< "Caution:  gridb may be different from grida, but gridb will be a clone of grida."<<endl
		<< "          This is most useful getting around byte swap limitations on fields."<<endl
		<< "Options:"<<endl
		<< " -vector fields are assumed vector fields (default assumes scalar field)"<<endl
		<< " -outdir set output directory for result"<<endl;
	exit(-1);
}
string splitgridname(string s)
{
	string separator(":,");
	int end_gridstring;
	end_gridstring=s.find_first_of(separator,0);
	if(end_gridstring<0) throw SeisppError(
		string("splitgrid function:  cannot parse argument=")
		+ s
		+string("\nArgument must have a : or , separating grid and field name") );
	string result;
	result.assign(s,0,end_gridstring);
	return(result);
}
string splitfieldname(string s)
{
	string separator(":,");
	int end_gridstring;
	end_gridstring=s.find_first_of(separator,0);
	if(end_gridstring<0) throw SeisppError(
		string("splitfield function:  cannot parse argument=")
		+ s
		+string("\nArgument must have a : or , separating grid and field name") );
	string result;
	result.assign(s,end_gridstring+1,s.length());
	return(result);
}	
bool SEISPP::SEISPP_verbose(true);
int main(int argc, char **argv)
{
	int i;
	ios::sync_with_stdio();
	if(argc<7) usage();

	bool vectormode(false);
	string dbname(argv[1]);
	if(string(argv[3])!="*")usage();
	if(string(argv[5])!="=") usage();
	string outdir("edited_fields");
	// Arg list can override pf contents in this logic
	for(i=7;i<argc;++i)
	{
		string argstr=string(argv[i]);
		if(argstr=="-vector")
			vectormode=true;
		else if(argstr=="-outdir")
		{
			++i;
			outdir=string(argv[i]);
		}
		else
			cerr << "Unknown argument = "<<argstr<<endl;
			usage();
	}
	try {
	    DatascopeHandle dbh(dbname,false);
	    dbh.lookup(string("gclgdisk"));
	    DatascopeHandle dbhg(dbh);
	    dbhg.lookup(string("gclfield"));
	    Dbptr db=dbh.db;
	    Dbptr dbgrd=dbhg.db;
	    /* Split up argv words into grid and field names */
	    string grida,fielda,gridb,fieldb;
	    grida=splitgridname(string(argv[2]));
	    gridb=splitgridname(string(argv[6]));
	    if(grida!=gridb)
	    {
		cerr << "grida = "<<grida << " and gridb="<<gridb<<endl
			<< "gridb will be a clone of a but with a different name tag only"
			<< endl;
	    }
	    fielda=splitfieldname(string(argv[2]));
	    fieldb=splitfieldname(string(argv[6]));
	    double c=atof(argv[4]);
	    string dfileout=fieldb;
		if(vectormode)
		{
			GCLvectorfield3d ga(db,grida,fielda);
			ga*=c;
			ga.name=gridb;
			ga.dbsave(dbgrd,gridb,outdir,fieldb,dfileout);
		}
		else
		{
			GCLscalarfield3d ga(db,grida,fielda);
			ga*=c;
			ga.name=gridb;
			ga.dbsave(dbgrd,gridb,outdir,fieldb,dfileout);
		}
	}
	catch (int ierr)
	{
		cerr << "GCLgrid library function threw error code="
			<< ierr << endl;
	}
	catch (...)
	{
		cerr << "Something threw an exception"<<endl;
	}
}
