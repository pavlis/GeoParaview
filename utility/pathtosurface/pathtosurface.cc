/* This program takes a path set as a collection of lon,lat pairs and 
projects this path in a specified direction and at a specified dip into
the subsurface.  The result is output as a GCLgrid object (2D).
*/
#include <string>
#include <iostream>
#include "stock.h"
#include "pf.h"
#include "coords.h"
#include "Metadata.h"
#include "seispp.h"
#include "gclgrid.h"
#include "dbpp.h"

using namespace std;
using namespace SEISPP;

void usage()
{
	cout << "pathtosurface db name [-pf pfname -V] "<<endl;
	exit(-1);
}
void copy_master(GCLgrid3d& gin, GCLgrid& gout,string outname)
{
	gout.name=outname;
	gout.lat0=gin.lat0;
	gout.lon0=gin.lon0;
	gout.r0=gin.r0;
	gout.azimuth_y=gin.azimuth_y;
	gout.lat0=gin.lat0;
        for(int i=0;i<3;++i)
	{
	         for(int j=0;j<3;++j) 
		 	gout.gtoc_rmatrix[i][j]=gin.gtoc_rmatrix[i][j];
	         gout.translation_vector[i]=gin.translation_vector[i];
	}
}
bool SEISPP::SEISPP_verbose(false);
int main(int argc, char **argv)
{
	if(argc<3) usage();
	int i,j,k;
	string dbname(argv[1]);
	string gridname(argv[2]);
	string pfname("pathtosurface");
	for(i=3;i<argc;++i)
	{
		string sarg(argv[i]);
		if(sarg=="-pf")
		{
			++i;
			if(i>=argc) usage();
			pfname=string(argv[i]);
		}
		else if(sarg=="-V")
			SEISPP_verbose=true;
		else
			usage();
	}
	Pf *pf;
	if(pfread(const_cast<char *>(pfname.c_str()),&pf))
	{
		cerr << "pfread error"<<endl;
		usage();
	}
	try {
		Metadata control(pf);
		double dip=control.get_double("projection_dip_angle");
		double az=control.get_double("projection_strike_azimuth");
		double dz=control.get_double("projection_depth_increment");
		int nz=control.get_int("number_depth_increments");
		string outputgridname=control.get_string("output_grid_name");
		string outdir=control.get_string("output_directory");
		DatascopeHandle dbh(dbname,false);
		GCLgrid3d grid(dbh.db,gridname);
		/* Load up lon,lat points */
		vector<double> pathlat,pathlon;
		if(SEISPP_verbose) cout << "Longitude Latitude"<<endl;
		while(!cin.eof())
		{
			double lat,lon;
			cin >> lon;
			cin >> lat;
			if(SEISPP_verbose)
			  cout << lon << " " <<lat<<endl;
			pathlon.push_back(lon);
			pathlat.push_back(lat);
		}
		int nx=pathlon.size();
		GCLgrid gout(nx,nz);
		/* There should be a standard method to do this. It is 
		currently lacking in the gclgrid library.*/
		copy_master(grid,gout,outputgridname);

		Geographic_point gp,gp0;
		Cartesian_point cp;
		double delta,ddelta,ddz;
		double lattmp,lontmp;
		ddelta=dz*cos(rad(dip));
		ddz=dz*sin(rad(dip));
		for(i=0;i<nx;++i)
		{
			gp0.lat=rad(pathlat[i]);
			gp0.lon=rad(pathlon[i]);
			gp0.r=r0_ellipse(gp0.lat);
			for(j=0;j<nz;++j)
			{
				gp=gp0;
				if(dip==90.0)
				{
					gp.r=gp0.r-dz*static_cast<double>(j);
				}
				else
				{
					gp.r=gp0.r-ddz*static_cast<double>(j);
					delta=ddelta*static_cast<double>(j);
					latlon(gp.lat,gp.lon,delta,az,
						&lattmp,&lontmp);
					gp.lat=lattmp;
					gp.lon=lontmp;
				}
				cp=gout.gtoc(gp);
				gout.x1[i][j]=cp.x1;
				gout.x2[i][j]=cp.x2;
				gout.x3[i][j]=cp.x3;
			}
		}

		gout.compute_extents();
		/* a couple leftover attributes need to be set */
		gout.dx1_nom=fabs(gout.x1high-gout.x1low)
					/static_cast<double>(nx);
		gout.dx2_nom=dz;
		gout.dbsave(dbh.db,outdir);
	} 
	catch (int ierr)
	{
		cerr << "GCLgrid error:  something threw error number "
			<< ierr<<endl;
	}
	catch (SeisppError serr)
	{
		serr.log_error();
	}
}

