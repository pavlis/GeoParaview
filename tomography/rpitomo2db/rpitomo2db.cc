#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include "stock.h"
#include "pf.h"
#include "metadata.h"
#include "gclgrid.h"
using namespace SEISPP;
// Ugly FORTRAN programs 
extern "C" {
void constants_();
void lltoutm_(int*, double *, double *,double *, double *,char *,int *,int);
void utmtoll_(int*, double *, double *,char *,double *, double *,int *,int);
}
void usage()
{
	cerr << "rpitomo2db tomofile db [-vpvs x ]\n";
	exit(-1);
}
int main(int argc, char **argv)
{
	int nx,ny,nz;
	FILE *fp;
	string fname;
	int ntotal,ntest;
	double lat,lon;
	double dx,dy,dz;
	double dxutm,dyutm;
	char *tomofile,*dbname;
	double vpvs=0.707;
	bool save_S=false;

	ios::sync_with_stdio();
	elog_init(argc,argv);
	if(argc<3) usage();
	tomofile=argv[1];
	dbname=argv[2];
	for(int iarg=0;iarg<argc;++iarg)
	{
		string sarg(argv[iarg]);
		if(sarg=="-vpvs")
		{
			++iarg;
			vpvs=atof(argv[iarg]);
			save_S=true;
		}
		else
			usage();
	}
	
	
	Pf *pf;
	char *pfname="rpitomo2db";
	if(pfread(pfname,&pf))
		elog_die(1,"pfread error\n");
	Dbptr db;
	if(dbopen(dbname,"r+",&db)==dbINVALID)
	{
		cerr << "Cannot open output database "<<dbname<<endl;
		exit(-1);
	}
	
	try {
		Metadata control(pf);
		double northing,easting;
		double northing0,easting0;
		double lat0,lon0;
		int zonenumber;
		string sin;
		char utmzone[2];		
		nx=control.get_int("nx");
		ny=control.get_int("ny");
		nz=control.get_int("nz");
		dx=control.get_double("dx");
		dy=control.get_double("dy");
		dxutm=dx*1000.0;
		dyutm=dy*1000.0;
		dz=control.get_double("dz");
		zonenumber=control.get_int("utm_zone_number");
		sin=control.get_string("utm_zone_letter_code");
		easting0=control.get_double("easting0");
		northing0=control.get_double("northing0");
		lat0=control.get_double("grid3d_lat0");
		lon0=control.get_double("grid3d_lon0");
		string gridname=control.get_string("output_gridname_base");
		string modelbasename=control.get_string("model_base_name");
		string fielddir=control.get_string("field_directory");
		string gclgdir=control.get_string("grid_directory");
		// convert these to radians
		lat0=rad(lat0);
		lon0=rad(lon0);
		utmzone[0]=sin[0];  utmzone[1]='\0';
		
		
		FILE *fp=fopen(tomofile,"r");
		if(fp==NULL)
		{
			cerr << "Open of file "<<tomofile
				<<"failed.  Exiting\n";
			exit(-1);
		}
		int ntotal=nx*ny*nz;
		float *buffer=new float[ntotal];
		ntest = fread(buffer,sizeof(float),ntotal,fp);
		if(ntest!=ntotal)
		{
			cerr << "Read error:  read "<<ntest
			   <<" numbers but expected" <<ntotal << endl;
			exit(-1);
		}
		// Steve's program seems to use this 
		int ReferenceEllipsoid=24;
		// necessary to initialize utm conversion function
		constants_();			
		
		// Get the latitude and longitude of origin and 
		// echo these to output for certification of how sensible this
		// result is
		double x,y,z;
		x=0.0;  y=0.0;
/*
		double lat0,lon0;
		utmtoll_(&ReferenceEllipsoid,&y,&x,utmzone,&lat0,&lon0,
			&zonenumber,strlen(utmzone));
		// the above fortran function puts out lat0 and lon0 
		// in degrees.  We need radians
		lat0=rad(lat0);  lon0=rad(lon0);
*/
		// We create the GCLgrid using this full constructor
		// but then reset all the points to use Steve's coordinate
		// system
		int nx0,ny0;  
		nx0=nx/2;
		ny0=ny/2;
		GCLgrid3d *g=new GCLgrid3d(nx,ny,nz,gridname,
			lat0,lon0,r0_ellipse(lat0),0.0,dx,dy,dz,nx0,ny0);
		// use this grid to create an empty scalar field that will hold
		// the velocity model we are going to input.
		GCLscalarfield3d vel(*g);
		// don't need g any longer
		delete g;
		// Now we reset the coordinates to the actual values using the same utm conversion
		// program Steve uses internally in his tomography program.
		// A GCLgrid3d scans x and y like Steve does (scans on x lines, moving to 
		// positive y) but the depth direction is reversed.  Compilcates things some, but
		// not horribly.  
		int i,j,k,kk;
		double x_lowerleft = easting0 - ((double)nx0)*dxutm;
		double y_lowerleft = northing0 - ((double)ny0)*dyutm;
		int offset=0;
		
		for(k=0,kk=vel.n3-1;k<vel.n3;++k,--kk)
			  for(j=0;j<vel.n2;++j) 
				for(i=0;i<vel.n1;++i,++offset)

		{
			Geographic_point geo;
			Cartesian_point cart;
			double ztrue, vtrue,vflattened;
			x = x_lowerleft + ((double)i)*dxutm;
			y = y_lowerleft + ((double)j)*dyutm;
			z = dz*(static_cast<double>(kk));
			utmtoll_(&ReferenceEllipsoid,&y,&x,utmzone,&lat,&lon,
				&zonenumber,strlen(utmzone));
			geo.lat=rad(lat);
			geo.lon=rad(lon);
			// Roecker uses utm coordinate in units of m.  GCLgrid uses
			// units of km.
			vflattened=static_cast<double>(buffer[offset]);
			// Steve's program uses the flattening transformation which
			// tweeks both the depth and velocity to approximately map spherical
			// shells to a parallelpiped.  We have to undo that transformation
			// when we build the gclgrid as it uses true geometry.
			ztrue=uflatz(z);
			vel.val[i][j][kk]=uflatvel(vflattened,z);
			// Steve's program does not use ellipticity, but we add it here.
			// could use a spherical case, but it will mesh cleaner with 
			// the gclgrid library if we do it this way.
			geo.r=r0_ellipse(rad(lat))-ztrue;
			cart=vel.gtoc(geo);
			vel.x1[i][j][k]=cart.x1;
			vel.x2[i][j][k]=cart.x2;
			vel.x3[i][j][k]=cart.x3;
		}
		// save field to database
		string modelname;
		modelname=modelbasename+"_P";
		// Use the model name for the output data file name
		dbsave(db,gclgdir,fielddir,modelname,modelname);
		if(save_S)
		{
			modelname=modelbasename+"_S";
			vel *= vpvs;
			dbsave(db,"",fielddir,modelname,modelname);
		}
	}
	catch (Metadata_get_error mderr)
	{
		mderr.log_error();
		exit(-1);
	}
	catch (int dberr)
	{
		elog_die(1,"dbsave threw exception code %d\n",dberr);
	}
	catch (...)
	{
		cerr << "Something threw an exception that wasn't handled." << endl;
	}
}
			
