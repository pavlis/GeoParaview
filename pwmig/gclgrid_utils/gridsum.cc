#include "dbpp.h"
#include "Metadata.h"
#include "seispp.h"
#include "gclgrid.h"
#include "VelocityModel_1d.h"
using namespace std;
using namespace SEISPP;

void usage()
{
	cerr << "gridsum db grida:fielda +[-] gridb:fieldb = gridc:fieldc [-vector -1d vmodel -1dmodtype [PS] -2dmode -outdir dir]"<<endl;
	cerr << " Note:"<<endl
		<< "(a) -1d switch projects 1d model vmodel to 3d model b"
		<<endl
		<< "  The model vmodel must be defined in mod1d table of df"
		<<endl
		<<"(b) + or - and = must be distinct arguments (white space separator)"
		<<endl;
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
	bool twodmode(false);
	bool signswitch(false);
	bool b1d(false);
	string dbname(argv[1]);
	if(string(argv[3])=="-")
		signswitch=true;
	else if(string(argv[3])!="+") 
		usage();
	if(string(argv[5]) != "=")
		usage();
	string vmod1dname;
	string outdir("edited_fields");
	string Vproperty("Pvelocity");
	// Arg list can override pf contents in this logic
	for(i=7;i<argc;++i)
	{
		string argstr=string(argv[i]);
		if(argstr=="-vector")
			vectormode=true;
		else if(argstr=="-2dmode")
			twodmode=true;
		else if(argstr=="-1d")
		{
			++i;
			b1d=true;
			vmod1dname=string(argv[i]);
		}
		else if(argstr=="-outdir")
		{
			++i;
			outdir=string(argv[i]);
		}
		else if(argstr=="-1dmodtype")
		{
			++i;
			argstr=string(argv[i]);
			if(argstr=="P")
				Vproperty=string("Pvelocity");
			else if(argstr=="S")
                                Vproperty=string("Svelocity");
			else
			{
				cerr << "Error parsing velocity model property"
					<< " must be P or S"<<endl;
				usage();
			}

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
	    string grida,fielda,gridb,fieldb,gridc,fieldc;
	    grida=splitgridname(string(argv[2]));
	    gridb=splitgridname(string(argv[4]));
	    gridc=splitgridname(string(argv[6]));
	    fielda=splitfieldname(string(argv[2]));
	    fieldb=splitfieldname(string(argv[4]));
	    fieldc=splitfieldname(string(argv[6]));
	    string dfileout=fielda+"_"+fieldb+"_sum";
cout << "Field a =" << grida << " " <<fielda <<endl;
cout << "Field b =" << gridb << " " <<fieldb <<endl;
cout << "Field c =" << gridc << " " <<fieldc <<endl;
	    if(twodmode)
	    {
		cerr << "-2dmode is not yet implemented.  Exiting."<<endl;
		usage();
	    }
	    else
	    {
		if(vectormode)
		{
			GCLvectorfield3d ga(db,grida,fielda);
			GCLvectorfield3d gb(db,gridb,fieldb);
			if(signswitch) gb*=-1.0;
			ga+=gb;
			ga.name=gridc;
			ga.dbsave(dbgrd,string(""),outdir,fieldc,dfileout);
		}
		else
		{
			GCLscalarfield3d ga(db,grida,fielda);
			GCLscalarfield3d *gb;
			if(b1d)
			{
				VelocityModel_1d vmod(db,vmod1dname,Vproperty);
				gb=new GCLscalarfield3d(ga);
				initialize_1Dscalar(*gb,vmod.v,vmod.z);
			}
			else
			{
				gb=new GCLscalarfield3d(db,gridb,fieldb);
			}
			if(signswitch) (*gb) *= -1.0;
			ga+=(*gb);
			ga.name=gridc;
			ga.dbsave(dbgrd,string(""),outdir,fieldc,dfileout);
		}
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
