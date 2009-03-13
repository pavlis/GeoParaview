#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include "stock.h"
#include "pf.h"
#include "seispp.h"
#include "Metadata.h"
#include "gclgrid.h"
using namespace SEISPP;
void usage()
{
	cerr << "xyz2gcl db gridname fieldname [-head n -pf pffile] < infile "<<endl;
	exit(-1);
}
bool SEISPP::SEISPP_verbose(true);
double distance(Cartesian_point& cp, GCLscalarfield3d& f, 
	int i, int j, int k)
{
	double ssqd(0.0);
	double dx=cp.x1-f.x1[i][j][k];
	ssqd+=dx*dx;
	dx=cp.x2-f.x2[i][j][k];
	ssqd+=dx*dx;
	dx=cp.x3-f.x3[i][j][k];
	ssqd+=dx*dx;
	return(sqrt(ssqd));
}
int main(int argc, char **argv)
{
	int nx,ny,nz;
	ios::sync_with_stdio();
	elog_init(argc,argv);
	if(argc<4) usage();

	string dbname(argv[1]);
	string gridname(argv[2]);
	string fieldname(argv[3]);
	string pffile("xyz2gcl");
	int nheaderlines(0);
	for(int iarg=4;iarg<argc;++iarg)
	{
		string sarg(argv[iarg]);
		if(sarg=="-pf")
		{
			++iarg;
			pffile=string(argv[iarg]);
		}
		if(sarg=="-head")
		{
			++iarg;
			nheaderlines=atoi(argv[iarg]);
			if(SEISPP_verbose) cout << "Skipping "
				<< nheaderlines <<" heading lines on input"
				<<endl;
		}
		else
			usage();
	}
	
	Pf *pf;
	if(pfread(const_cast<char *>(pffile.c_str()),&pf))
		elog_die(1,"pfread error\n");
	Dbptr db;
	if(dbopen(const_cast<char *>(dbname.c_str()),"r+",&db)==dbINVALID)
	{
		cerr << "Cannot open output database "<<dbname<<endl;
		exit(-1);
	}
	
	try {
		Metadata control(pf);
		int nx,ny,nz;
		double lat0,lon0,z0,r0,azimuth_y;
		nx=control.get_int("nx");
		ny=control.get_int("ny");
		nz=control.get_int("nz");
		lat0=control.get_double("lat0");
		lon0=control.get_double("lon0");
		z0=control.get_double("z0");
		r0=r0_ellipse(rad(lat0))-z0;
		azimuth_y=control.get_double("azimuth_y");
		string fielddir=control.get_string("field_directory");
		string gclgdir=control.get_string("grid_directory");
		// convert these to radians
		lat0=rad(lat0);
		lon0=rad(lon0);
		int ntotal=nx*ny*nz;

		GCLscalarfield3d field(nx,ny,nz);
		/* We need to set these immediately so we georeferencing
		works correctly */
		field.lat0=lat0;
		field.lon0=lon0;
		field.r0=r0;
		field.azimuth_y=rad(azimuth_y);
		field.set_transformation_matrix();
		/* Loop over input data in an assumed order.  Exit
		if count does not match ntotal.  Assume that if order
		is different for some data there will eventually be
		a filter to reorder to this fixed standard */
		char line[256];
		/* We buffer up the whole model.  Not memory 
		efficient, but assume this is not an issue */
		vector<double> latin,lonin,zin,modin;
		latin.reserve(ntotal);
		lonin.reserve(ntotal);
		zin.reserve(ntotal);
		modin.reserve(ntotal);
		int iskip(0);
		while(cin.getline(line,256))
		{
			double valin;
			if(iskip<nheaderlines)
			{
				++iskip;
				continue;
			}
			stringstream ss(line);
			ss >> valin;
			latin.push_back(valin);
			ss >> valin;
			lonin.push_back(valin);
			ss >> valin;
			zin.push_back(valin);
			ss >> valin;
			modin.push_back(valin);
		}
		if(latin.size()!=ntotal)
		{
			cerr << "xyz2gcl:  input data size mismatch"<<endl
			  << "parameter file declares nx,ny,nz="
			  << nx <<", "<<ny<<", "<<nz<<endl
			  <<"Expected to read "<<nx*ny*nz <<" lines"<<endl
			  <<"Actual input count="<<latin.size();
			exit(-1);
		}
		/* convert lat and lons to radians */
		int i,j,k;
		for(i=0;i<ntotal;++i)
		{
			latin[i]=rad(latin[i]);
			lonin[i]=rad(lonin[i]);
		}
		/* Now load coordinates and field data */
		Cartesian_point cp;
		for(int ipt=0;ipt<ntotal;++ipt)
		{
			int stride1,stride2;
			/*
			stride1=ny*nz;
			stride2=nz;
			i=ipt/stride1;
			j=(ipt/stride2)%ny;
			k = nz - (ipt%nz) - 1;
			*/
			stride1=nx*ny;
			stride2=nx;
			k = nz - ipt/stride1 - 1;
			j=(ipt/stride2)%ny;
			i= ipt%nx;
			double r=r0_ellipse(latin[ipt]);
			r-=zin[ipt];
			cp=field.gtoc(latin[ipt],
				lonin[ipt],r);
//DEBUG
/*
cout << "Geographic:  "<< deg(latin[ipt])
  <<", "<<deg(lonin[ipt])
  <<", "<<r<<endl
  <<"Cartesian:  "<<cp.x1
  <<", "<<cp.x2
  <<", "<<cp.x3<<endl;
  */
			field.x1[i][j][k]=cp.x1;
			field.x2[i][j][k]=cp.x2;
			field.x3[i][j][k]=cp.x3;
			field.val[i][j][k]=modin[ipt];
		}
		/* Now we have to set the other attributes in the grid */
		field.compute_extents();
		field.dx1_nom=fabs(field.x1high-field.x1low)
				/ static_cast<double>(nx-1);
		field.dx2_nom=fabs(field.x1high-field.x1low)
				/ static_cast<double>(ny-1);
		field.dx3_nom=fabs(field.x3high-field.x3low)
				/ static_cast<double>(nz-1);
		cp=field.gtoc(lat0,lon0,r0);
		/* These best set to center before trying this method */
		field.i0=nx/2;
		field.j0=ny/2;
		field.k0=nz/2;
		int originindex[3];
		if(field.lookup(cp.x1,cp.x2,cp.x3) )
		{
			/* Lookup can fail with rectilinear grids.  
			When that happens we revert to a brute
			force search method.  A warning is issued
			because such a model means future trouble*/
			cerr << "Warning:  lookup method failed in finding"
				<< " origin point in this grid"<<endl
				<< "Searching to set origin indices"
				<<endl;
			originindex[0]=0;
			originindex[1]=0;
			originindex[2]=0;
			double dmin,dtest;
			dmin=distance(cp,field,0,0,0);
			for(i=0;i<nx;++i)
			  for(j=0;j<ny;++j)
			    for(k=1;k<nz;++k)
			    {
			    	dtest=distance(cp,field,i,j,k);
				if(dtest<dmin)
				{
					dmin=dtest;
					originindex[0]=i;
					originindex[1]=j;
					originindex[2]=k;
				}
			    }
		}
		else
		{
			field.get_index(originindex);
		}
		field.i0=originindex[0];
		field.j0=originindex[1];
		field.k0=originindex[2];
		field.name=gridname;

		/* save and we're done */
		field.dbsave(db,gclgdir,fielddir,fieldname,fieldname);
	} catch (SeisppError serr)
	{
		serr.log_error();
	}
	catch (int ierr)
	{
		cerr << "GCLgrid error:  something threw error number "
			<< ierr<<endl;
	}
}
