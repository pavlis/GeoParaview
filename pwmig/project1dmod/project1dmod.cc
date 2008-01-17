#include <vector>
#include "seispp.h"
#include "dbpp.h"
#include "gclgrid.h"
#include "VelocityModel_1d.h"
using namespace std;
using namespace SEISPP;

void usage()
{
	cerr << "project1dmod db gridname vmodname [-field fieldname -mt modtype]"<<endl
	<< "Default fieldname=vmodname, modtype=P"<<endl;
	exit(-1);
}
int main(int argc, char **argv)
{
	if(argc<4) usage();
	string dbname(argv[1]);
	string gridname(argv[2]);
	string vmodname(argv[3]);
	string modtype("P");
	string fieldname=vmodname;
	int i;
	for(i=4;i<argc;++i)
	{
		string argstr(argv[i]);
		if(argstr=="-field")
		{
			++i;
			fieldname=string(argv[i]);
		}
		else if(argstr=="-mt")
		{
			++i;
			modtype=string(argv[i]);
		}
		else
		{
			cerr << "Illegal argument="<<argstr<<endl;
			usage();
		}
	}
	try {
		DatascopeHandle dbh(dbname,false);
		DatascopeHandle dbhmod(dbh);
		dbhmod.lookup("mod1d");
		DatascopeHandle dbhgrid(dbh);
		dbhgrid.lookup("gclgdisk");
		GCLgrid3d grid(dbhgrid.db,gridname);
		GCLscalarfield3d mod(grid);
		VelocityModel_1d Vmod(dbhmod.db,vmodname,modtype);
		/* the following works only when the grid has constant
		depth slices for index 3 constant */
		vector<double> zi,vi;
		for(i=mod.n3-1;i>=0;--i)
		{
			double z;
			z=mod.depth(0,0,i);
			zi.push_back(z);
			vi.push_back(Vmod.getv(z));
		}
		initialize_1Dscalar(mod,vi,zi);
		/* for now this is a frozen directory name */
		string fielddir("vmodels");
		string gclgdir("");  /* null as save is not needed*/
		/* Use the fieldname as the file name.  Restrictive,
		but this code may never leave the neighborhood. */
		mod.dbsave(dbhgrid.db,gclgdir,fielddir,fieldname,fieldname);
	} catch (SeisppError serr)
	{
		serr.log_error();
	}
	catch (int ierr)
	{
		cerr << "Some GCLgrid routine threw error number="<<ierr<<endl;
	}

	catch (GCLgrid_error gclerr)
	{
		gclerr.log_error();
	}
}
