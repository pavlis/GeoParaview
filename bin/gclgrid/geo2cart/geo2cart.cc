#include <iostream>
#include <vector>
#include <list>
#include "seispp.h"
#include "gclgrid.h"
extern "C" {
    extern void treex3_(double *, int *, double *, int *, double *);
}
using namespace std;
using namespace SEISPP;
/* This is used to convert exported GMT boundary lines into gclgrid coordinate
by Yinzhi Wang @ Indiana University
   */

void usage()
{
	cerr << "geo2cart gridname "
		<< "[-v] < infile"<<endl
	      << "output is stdout"<<endl;

	exit(-1);
}
bool SEISPP::SEISPP_verbose(false);
int main(int argc, char **argv)
{
	int i,j;
	ios::sync_with_stdio();
	if(argc<2) usage();
	string gridname(argv[1]);
	for(i=2;i<argc;++i)
	{
		string argstr=string(argv[i]);
		if(argstr=="-v")
			SEISPP_verbose=true;
		else
		{
			cerr << "Unknown argument = "<<argstr<<endl;
			usage();
		}
	}
	try 
	{
	    GCLgrid g(gridname);
        double lat,lon;
        string temp_line,line;
        while(getline(cin,temp_line))
        {
         	stringstream line(temp_line);
			if(line.str()[0]=='#') continue;
			if(line.str()[0]=='>')
				cout<<line.str()<<endl;
			else
			{
		       	line>>lon;
		       	line>>lat;
		       	Cartesian_point pot=g.gtoc(rad(lat),rad(lon),g.r0);
		       	if(g.lookup(rad(lat),rad(lon))==0)
		       	{
		       		int indx[2];
		       		g.get_index(indx);
		       		i=indx[0];
		       		j=indx[1];
		       		dmatrix J(3,3),Jinv(3,3);
		       		dvector dxraw(3),dxunit(3);
		       		int three(3);
		       		double det;
		       		double dxi[3],dxj[3],dxk[3];
		       		
		       		dxi[0] = (g.x1[i+1][j]) - (g.x1[i][j]);
					dxi[1] = (g.x2[i+1][j]) - (g.x2[i][j]);
					dxi[2] = (g.x3[i+1][j]) - (g.x3[i][j]);
					dxj[0] = (g.x1[i][j+1]) - (g.x1[i][j]);
					dxj[1] = (g.x2[i][j+1]) - (g.x2[i][j]);
					dxj[2] = (g.x3[i][j+1]) - (g.x3[i][j]);
					
					dxraw(0) = pot.x1 - (g.x1[i][j]);
					dxraw(1) = pot.x2 - (g.x2[i][j]);
					dxraw(2) = pot.x3 - (g.x3[i][j]);
					for(int ii=0;ii<3;ii++)
					{
						J(ii,0)=dxi[ii];
						J(ii,1)=dxj[ii];
						J(ii,2)=1;
					}
					treex3_(J.get_address(0,0),&three,Jinv.get_address(0,0),&three,&det);
					dxunit=Jinv*dxraw;
					
					cout<<(indx[0]+dxunit(0))*100<<"	"<<(indx[1]+dxunit(1))*100<<"	"<<0<<endl;
				}
    		}
        }
	}
	catch (SeisppError& serr)
	{
		serr.log_error();
	}
	catch (GCLgridError& serr)
	{
		cerr<<serr.what()<<endl;
	}
}
