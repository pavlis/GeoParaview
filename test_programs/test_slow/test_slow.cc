#include "coords.h"
#include "stock.h"
#include "pf.h"
#include "seispp.h"
using namespace SEISPP;
int main(int argc, char **argv)
{
	Pf *pf;
	pfread("test_slow",&pf);
	string tag("Slowness_Grid_Example");
	Rectangular_Slowness_Grid g(pf,tag);
	for(int i=0;i<g.nux;++i)
		for(int j=0;j<g.nuy;++j)
		{
			Slowness_vector u=g.slow(i,j);
			cout << g.ux(i) << " "
				<< u.ux << " "
				<< g.uy(j) << " "
				<< u.uy << " "
				<< deg(u.azimuth()) << " "
				<< deg(u.baz()) << endl;
		}
}
