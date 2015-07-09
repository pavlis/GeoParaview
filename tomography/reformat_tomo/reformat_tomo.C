#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <sunmath.h>
#include "gclgrid.h"
// Ugly FORTRAN programs 
extern "C" {
void constants_();
void lltoutm_(int*, double *, double *,double *, double *,char *,int *,int);
void utmtoll_(int*, double *, double *,char *,double *, double *,int *,int);
}
int main(int argc, char **argv)
{
	int nx,ny,nz;
	FILE *fp;
	string fname;
	int ntotal,ntest;
	double lat,lon;

	ios::sync_with_stdio();
	cout << "Enter nx,ny, and nz: ";
	cin >> nx;
	cin >> ny;
	cin >> nz;
	cout << "Creating grid of "<<nx<<"X"<<ny<<"X"<<nz<<endl;
	int nx0,ny0;
	nx0=nx/2;
	ny0=ny/2;
	cout << "Assuming origin is at grid position (x,y)=("
		<< nx0 <<","<<ny0<<")"<<endl;
	cout << "Enter dx,dy,dz:";
	double dx,dy,dz;
	cin >> dx;  cin >> dy;  cin >> dz;
	double z_to_output;
	cout << "Enter depth of layer to output (Note this is nearest depth; cannot interpolate):";
	cin >> z_to_output;
	int iz=nint(z_to_output/dz);
	cout << "Will output grid at layer index = " << iz+1<<endl;
	cout << "Enter binary tomography output file to convert:";
	cin >> fname;
	fp = fopen(fname.c_str(),"r");
	if(fp==NULL) 
	{
		cerr << "Open of file "<<fname<<"failed.  Exiting\n";
		exit(-1);
	}
	ntotal=nx*ny*nz;
	float *buffer=new float[ntotal];
	ntest = fread(buffer,sizeof(float),ntotal,fp);
	if(ntest!=ntotal)
	{
		cerr << "Read error:  read "<<ntest<<"numbers but expected" << 
			ntotal << endl;
		exit(-2);
	}
	// For Indiana tomography (jeidi wu) origin northing=4427876.9194520
	// and easting = 542679.94500658.  Zone number is 16 and utmzone="T"
	// This is lat=40.0 and lon = -86.5
	double northing,easting;
	double northing0,easting0;
	integer zonenumber;
	string sin;
	char utmzone[2];
	cout << "Enter UTM Zone number for origin: ";
	cin >> zonenumber;
	cout << "Enter UTM zone letter code (depends on latitude):";
	cin >> sin;
	cout << "Enter UTM easting and northing coordinates for grid origin:";
	cin >> easting0;
	cin >> northing0;
	utmzone[0]=sin[0];  utmzone[1]='\0';
	// Steve's program seems to use this 
	int ReferenceEllipsoid=24;
	// necessary to initialize utm conversion function
	constants_();
	//utmtoll_(&ReferenceEllipsoid,&northing0,&easting0,utmzone,&lat,&lon,&zonenumber,strlen(utmzone));
	ofstream onp;
	string fnout;
	cout << "Enter file name for ascii output" << endl;
	cin >> fnout;
	onp.open(fnout.c_str(),ios::out);
	// Have them in memory now.  Just a matter of rearranging them.
	// We'll push data to the vels array instead of writing out
	// buffer directly
	/* To get coordinates right note this about Steve's grid geometry:
		1.  0 point in the array we just read is at the surface (z=0) in the
			southwest corner of the grid
		2.  Fastest scan (i.e. most rapidly varying index) is west to east.
		3.  next is north south
		4.  final index is layers going deeper in the earth as we get deeper
			into the file 
	*/


	int i,j,k,kk,ioff;
	int offset;
	dmatrix layer(nx,ny); // Will output only one layer at a time
	offset =iz*nx*ny;
	for(i=0,ioff=0;i<nx;++i)
	{
		for(j=0;j<ny;++j,ioff++)
		{
			layer(i,j)=(double)buffer[offset+ioff];
		}
	}
	// remove the background to convert to difference 
	// from reference.  Use the mean as the reference.  
	// This is not ideal, but should get us off the dock.
	/* not included for now 
	double mean;
	for(j=0,mean=0.0;j<ny;++j)
		for(i=0;i<nx;++i) mean += layer(i,j);
	mean /= ((double)(nx*ny));
	for(j=0;j<ny;++j)
		for(i=0;i<nx;++i) layer(i,j) -= mean;
*/
	double x_lowerleft, y_lowerleft;
	double x,y,z;
	x_lowerleft = easting0 - ((double)nx0)*dx;
	y_lowerleft = northing0 - ((double)ny0)*dy;
	for(j=ny-1;j>=0;--j)
		for(i=0;i<nx;++i) 
		{
			x = x_lowerleft + ((double)i)*dx;
			y = y_lowerleft + ((double)j)*dy;
			z = ((double)iz)*dz;
			utmtoll_(&ReferenceEllipsoid,&y,&x,utmzone,&lat,&lon,&zonenumber,strlen(utmzone));
			onp << lat <<" "
				<< lon <<" "
				<< z <<" "
				<<  layer(i,j) << endl;
		}

	onp.close();
}
			
