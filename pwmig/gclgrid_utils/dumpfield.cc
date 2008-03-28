#include "seispp.h"
#include "dbpp.h"
#include "gclgrid.h"
using namespace std;
using namespace SEISPP;
void usage()
{
	cerr << "dumpfield db gridname fieldname [-component n -3d]"<<endl
		<< " outputs lon,lat,z,value to stdout"<<endl
		<< " -component sets vector mode writing component n"<<endl;
	exit(-1);
}
bool SEISPP::SEISPP_verbose(true);
int main(int argc, char **argv)
{
	if(argc < 4) usage();
	string dbname(argv[1]);
	string gridname(argv[2]);
	string fieldname(argv[3]);
	bool threedmode(false);
	bool vectormode(false);
	int component;
	int i,j,k;
	for(i=4;i<argc;++i)
	{
		string argstr(argv[i]);
		if(argstr=="-3d") 
			threedmode=true;
		else if(argstr=="-component")
		{
			vectormode = true;
			++i;
			component = atof(argv[i]);
		}
		else
			usage();
	}
	try {
		DatascopeHandle dbh(dbname,true);
		if(threedmode)
		{
		   if(vectormode)
		   {
			GCLvectorfield3d field(dbh.db,gridname,fieldname);
			for(i=0;i<field.n1;++i)
			  for(j=0;j<field.n2;++j)
			     for(k=0;k<field.n3;++k)
				cout << deg(field.lon(i,j,k)) <<" "
					<<deg(field.lat(i,j,k))<<" "
					<<field.depth(i,j,k)<<" "
					<<field.val[i][j][k][component]<<endl;
		   }
		   else
		   {
			GCLscalarfield3d field(dbh.db,gridname,fieldname);
			for(i=0;i<field.n1;++i)
			  for(j=0;j<field.n2;++j)
			     for(k=0;k<field.n3;++k)
				cout << deg(field.lon(i,j,k)) <<" "
					<<deg(field.lat(i,j,k))<<" "
					<<field.depth(i,j,k)<<" "
					<<field.val[i][j][k]<<endl;
		    }
			
		}
		else
		{
			GCLgrid grid(dbh.db,gridname);
			for(int j=0;j<grid.n2;++j)
			    for(int i=0;i<grid.n1;++i)
				cout << deg(grid.lon(i,j)) <<" "<<deg(grid.lat(i,j))<<endl;
		}
	} catch (int ierr)
	{
		cerr << "GCLgrid constructor threw exception ="<<ierr<<endl;
		usage();
	}
	catch (SeisppError serr)
	{
		serr.log_error();
		usage();
	}
}
