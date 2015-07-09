#include <iostream>
#include "coords.h"
#include "stock.h"
#include "pf.h"
#include "slowness.h"
using namespace SEISPP;
bool SEISPP::SEISPP_verbose(true);
int main(int argc, char **argv)
{
	Pf *pf;
	pfread("test_slow",&pf);
	string tag("Slowness_Grid_Example");
	RectangularSlownessGrid g(pf,tag);
	for(int i=0;i<g.nux;++i)
		for(int j=0;j<g.nuy;++j)
		{
			SlownessVector u=g.slow(i,j);
			cout << g.ux(i) << " "
				<< u.ux << " "
				<< g.uy(j) << " "
				<< u.uy << " "
				<< deg(u.azimuth()) << " "
				<< deg(u.baz()) << endl;
		}
                cout << "Testing arithmetic operators "<<endl;
                SlownessVector a=g.slow(0,0);
                SlownessVector b=g.slow(1,1);
                cout << "a ux,uy="<<a.ux<<" "<<a.uy<<endl;
                cout << "b ux,uy="<<b.ux<<" "<<b.uy<<endl;
                SlownessVector c;
                c=a+b;
                cout << "a+b ux,uy"<<c.ux<<" "<<c.uy<<endl;
                c=a-b;
                cout << "a-b ux,uy"<<c.ux<<" "<<c.uy<<endl;
}
