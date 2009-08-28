#include <stdio.h>
#include <string>
#include "stock.h"
#include "coords.h"
#include "pf.h"
#include "db.h"
#include "gclgrid.h"
using namespace std;

int main(int argc, char **argv)
{
	ios::sync_with_stdio();
	int nx=281;
	int ny=281;
	const int nt=400;
	int nz=nt;
	const double dt(0.2);
	const double t0bias=10000000.0; // Done to make absolute times something besides 0.0 
	// we make these constant for now
	const string dir(".");
	const string dfile("response.bin");
	const string datatype("t4");
	const double calib(1.0);
	int decfac=1;
	if(argc==2)
	{
		decfac=atoi(argv[1]);
		cout << "decfac set to "<<decfac<<endl;
	}

	cout << "Calling dbopen"<<endl;
	Dbptr db;
	dbopen("simdata","r+",&db);

	// for gclgrid definition
	// note r0 being negative forces grid plane to be on reference 
	// ellipsoid at sea level
	const double lat0(rad(0.0)),lon0(rad(45.0)),r0(-1.0); 
	const double az(0.0);
	const double dx1n(4.0), dx2n(4.0),dx3n(1.0);
	int iorigin=140;
	int jorigin=140;
	string gridname("fdsim");
	string gclgdir("grids");
	string fielddir=gclgdir;
	string fieldname("subduct");
	string dfileout("subduct.dat");
	if(decfac>1) 
	{
		char appendage[20];
		sprintf(appendage,"_%1.1d",decfac);
		string a(appendage);
		gridname=gridname+a;
		fieldname=fieldname+a;
		dfileout=dfileout+a;
	}
	cout << "Creating parent grid"<<endl;
	int nxgrid,nygrid,nzgrid;
	nxgrid=nx/decfac;  nygrid=ny/decfac;  nzgrid=nz; 
	iorigin=iorigin/decfac;  jorigin=jorigin/decfac;
	GCLgrid3d *parent=new GCLgrid3d(nxgrid,nygrid,nzgrid,gridname,
		lat0,lon0,r0,az,dx1n,dx2n,dx3n,iorigin,jorigin);
	cout << "Creating vector field of size" 
		<<nxgrid << "X"
		<<nygrid << "X"
		<<nzgrid << endl;
	GCLvectorfield3d grid(*parent,3);
	delete parent;
	cout << "Past startup.  Trying to open input file "<<endl;
	int ix,iy,iz,component;
	int foff=0;
	const int doff(1600);
	string fname=dir+string("/")+dfile;
	FILE *fp=fopen(fname.c_str(),"r");
	float rawdata[nt];

	cout << "Starting read data loop"<<endl;
	for(component=0;component<3;++component)
	{
		for(int i=0,ix=0;ix<nx;++ix,i<nxgrid)
		{
			for(int j=0,iy=0;iy<ny;++iy,j<nygrid)
			{
				fread(rawdata,sizeof(float),nt,fp);
				if( (ix%decfac == 0) && (iy%decfac==0) )
				{
				    for(iz=0;iz<nt;++iz) 
					grid.val[i][j][iz][component]=rawdata[iz];
				}
				if(iy%decfac == 0) 
				{
					cout << "Finished input col "<<iy<<
					"= output row "<< j<< endl;
					++j;
				}
			}
			if(ix%decfac==0) 
			{
				cout << "Finished input row "<<ix<<
				"= output row "<< i<< endl;
				++i;
			}
		}
		cout << "Finished component "<< component<<endl;
	}
	cout << "Calling dbsave"<<endl;
	grid.dbsave(db,gclgdir,fielddir,fieldname,dfileout);
}
