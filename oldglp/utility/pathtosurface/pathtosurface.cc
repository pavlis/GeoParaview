/* This program takes a path set as a collection of lon,lat pairs and 
projects this path in a specified direction and at a specified dip into
the subsurface.  The result is output as a GCLgrid object (2D).

Older versions of this code before Sept 2011 had a dumb bug with 
the dip trig function inverted.  (multiply instead of divide).
This was retained in comments to fix any results using this. 
The only one I remember for sure is the early surface I built for
the SE Caribbean fold belt and the north andean block.  

New version here does a smarter thing anyway and requires specifying
the distance increment which is more stable because one cannot accidentally
do a divide by zero mistake.  

Also added option of a negative offset special for STEEP01 surface
projections.
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
	cout << "pathtosurface db outname [-pf pfname -V] < pathfile "<<endl
		<<"with db as input database and outname being the output gridname for this surface"<<endl
		<<"pathfile is list of lon-lat pairs"  <<endl
		<<"WARNING:  outname must be unique"<<endl;
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
	string outputgridname(argv[2]);
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
                az=rad(az);
		//double dz=control.get_double("projection_depth_increment");
                double dz,dr;
		// Note this isn't really a depth for dips other than 90
                dr=control.get_double("projection_depth_increment");
                /*Normally this parameter should be zero. Positive moves
                  start of surface coordinates in azimuth direction.  Use
                  negative offset in opposite direction */
                double origin_offset=control.get_double("origin_offset");
                /* this combination must be trapped to avoid divide by
                   zero errors.  constant is a bit arbitrary */
                if(fabs(dip)>89.0 && (fabs(origin_offset)>0.1))
                {
                    cerr << "Illegal parameter combination"<<endl
                        << "dip near 90 degrees = "<<dip
                        << "not allowed with nonzero value of origin_offset"
                        <<endl;
                    exit(-1);
                }
		int nr=control.get_int("number_depth_increments");
		string gridname=control.get_string("reference_grid_name");
		string outdir=control.get_string("output_directory");
		DatascopeHandle dbh(dbname,false);
		GCLgrid3d grid(dbh,gridname);
		/* Load up lon,lat points */
		vector<double> pathlat,pathlon,pathz;
		if(SEISPP_verbose) cout << "Longitude Latitude Depth"<<endl;
		while(!cin.eof())
		{
			double lat,lon,z;
			cin >> lon;
			cin >> lat;
                        cin >> z;
			if(SEISPP_verbose)
			  cout << lon << " " <<lat<<" "<<z<<endl;
			pathlon.push_back(lon);
			pathlat.push_back(lat);
			pathz.push_back(z);
		}
		int nx=pathlon.size();
		GCLgrid gout(nx,nr);
		/* There should be a standard method to do this. It is 
		currently lacking in the gclgrid library.*/
		copy_master(grid,gout,outputgridname);

		Geographic_point gp,gp0;
                double z0;
		Cartesian_point cp;
		double delta,ddelta,ddz;
		double lattmp,lontmp;
                /* This is value for equator.  */
                const double kmperdeg(110.57);
                /*  This was the old bug 
		ddelta=dz*cos(rad(dip))/kmperdeg;
		ddz=dz*sin(rad(dip));
                */
                ddelta=dr*cos(rad(dip))/kmperdeg;
                ddz=dr*sin(rad(dip));
		for(i=0;i<nx;++i)
		{
			gp0.lat=rad(pathlat[i]);
			gp0.lon=rad(pathlon[i]);
			gp0.r=r0_ellipse(gp0.lat)-pathz[i];
                        z0=pathz[i];
                        /* Safe test assumes units of km */
                        if(fabs(origin_offset)>1.0)
                        {
                            /* This is an approximation but assume 
                               one would not normally project large 
                               enought distances for this to matter */
                            double delta_offset;
                            delta_offset=fabs(origin_offset)*cos(rad(dip));
                            delta_offset/=kmperdeg;
                            delta_offset=rad(delta_offset);
                            double z_offset=origin_offset*sin(rad(dip));
                            z0=pathz[i]-z_offset;
                            if(origin_offset>0.0)
                            {
                                latlon(gp0.lat,gp0.lon,delta_offset,az,
                                        &lattmp,&lontmp);
                                z0=pathz[i]+z_offset;
                            }
                            else
                            {
                                latlon(gp0.lat,gp0.lon,delta_offset,az+M_PI,
                                        &lattmp,&lontmp);
                                z0=pathz[i]-z_offset;
                                if(z0<0.0)
                                {
                                  cerr << "Warning:  path through point "
                                    "number "<<i<<" projects above sea level"
                                    <<endl
                                    <<"Origin depth="<<z0<<endl;
                                }
                            }
                            gp0.lat=lattmp;
                            gp0.lon=lontmp;
                            gp0.r=r0_ellipse(gp0.lat)-z0;

                        }
			for(j=0;j<nr;++j)
			{
				if(fabs(dip)>89.99)
				{
					gp.r=gp0.r-dr*static_cast<double>(j);
                                        gp.lat=gp0.lat;
                                        gp.lon=gp0.lon;
				}
				else
				{
					delta=rad(ddelta*static_cast<double>(j));
					latlon(gp0.lat,gp0.lon,delta,az,
						&lattmp,&lontmp);
					gp.lat=lattmp;
					gp.lon=lontmp;
                                        gp.r=r0_ellipse(gp.lat)
                                            -z0 - ddz*static_cast<double>(j);
				}
                                //DEBUG
                                cout <<i<<" "<<j<<" "
                                    <<deg(gp.lat)<<" "
                                    <<deg(gp.lon)<<" "
                                    <<r0_ellipse(gp.lat)-gp.r<<endl;
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
		gout.save(dbh,outdir);
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

