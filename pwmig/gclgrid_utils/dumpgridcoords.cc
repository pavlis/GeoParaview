#include "seispp.h"
#include "dbpp.h"
#include "gclgrid.h"
using namespace std;
using namespace SEISPP;
void usage()
{
	cerr << "dumpgridcoords db gridname [-3c]"<<endl
		<< " outputs lon,lat pairs to stdout"<<endl;
	exit(-1);
}
int main(int argc, char **argv)
{
	if(argc < 3) usage();
	string dbname(argv[1]);
	string gridname(argv[2]);
	bool threecmode(false);
	if(argc==4) 
	{
		string sw3c(argv[3]);
		if(sw3c!="-3c") usage();
		threecmode=true;
	}
	try {
		DatascopeHandle dbh(dbname,true);
		if(threecmode)
		{
			GCLgrid3d grid(dbh.db,gridname);
			int k=grid.n3 - 1;
			for(int j=0;j<grid.n2;++j)
			    for(int i=0;i<grid.n1;++i)
				cout << deg(grid.lon(i,j,k)) <<" "
					<<deg(grid.lat(i,j,k))<<endl;
			
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
