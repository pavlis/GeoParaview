#include <stdio.h>
#include <string>
#include "stock.h"
#include "coords.h"
#include "pf.h"
#include "db.h"
#include "gclgrid.h"
#include "Hypocenter.h"
#include "seispp.h"
using namespace std;
using namespace SEISPP;
string virtual_station_name(int ix1, int ix2)
{
        char name[8];
        sprintf(name,"%03.3d%03.3d",ix1,ix2);
        return(string(name));
}

bool SEISPP::SEISPP_verbose(true);
int main(int argc, char **argv)
{
	ios::sync_with_stdio();
	const int nx(281);
	const int ny(281);
	const int nt(400);
	const double dt(0.2);
	const double t0bias=10000000.0; // Done to make absolute times something besides 0.0 
	// we make these constant for now
	const string dir(".");
	const string dfile("response.bin");
	const string datatype("t4");
	const double calib(1.0);
	const double atime(9999257.25);

	// for gclgrid definition
	// note r0 being negative forces grid plane to be on reference 
	// ellipsoid at sea level
	const double lat0(rad(0.0)),lon0(rad(45.0)),r0(-1.0); 
	const double az(0.0);
	const double dx1n(4.0), dx2n(4.0);
	const int iorigin(140),jorigin(140);
	string gridname("fdsim");
	GCLgrid grid(nx,ny,gridname,lat0,lon0,r0,az,dx1n,dx2n,iorigin,jorigin);
	// Eventually need a 3d grid here to build model using uflatz and uflatvel
	// to convert to spherical geometry

	int ix,iy,component;
	int i;
	// Needed for time calculations
	double starttime;
	const double dt0(0.0);  // Time from start of data to P arrival
	Hypocenter hypo;
	hypo.lat=-grid.lat0+rad(1.0);
	hypo.lon=-grid.lon0;
	hypo.z=0.0;
	hypo.time=0.0;
	starttime=dt0-hypo.ptime(grid.lat0,grid.lon0,0.0);
	// evil way to do this.  Generally considered bad form in OOP
	hypo.time=t0bias+starttime;
	// Now push increment starttime by t0bias
	starttime+= t0bias;
	double endtime=starttime+400.0*dt;
	string staname;
	string chan;
	int arid;
	Dbptr db;
	Dbptr dbsite;
	Dbptr dbsitechan;
	Dbptr dbarrival;
	Dbptr dbassoc;
	dbopen("simdata","r+",&db);
	db=dblookup(db,0,"wfdisc",0,0);
	dbsite=dblookup(db,0,"site",0,0);
	dbsitechan=dblookup(db,0,"sitechan",0,0);
	dbarrival=dblookup(db,0,"arrival",0,0);
	dbassoc=dblookup(db,0,"assoc",0,0);
	//
	// save the grid
	//
	grid.dbsave(db,string("grids"));
	int foff=0;
	const int doff(1600);
	int chanid=0;

	for(component=1;component<=3;++component)
	{
		switch (component)
		{
		case 1:
			chan="E";
			break;
		case 2:
			chan="N";
			break;
		case 3:
		default:
			chan="Z";
		}
		for(ix=0;ix<nx;++ix)
		{
			for(iy=0;iy<ny;++iy,foff+=doff)
			{
			/* Decimate by 5 and use only range 70 to 210.  
			This reduces data volume to approaching real data.*/
			if( (ix%5 == 0) && (iy%5 == 0) 
				&&(ix>=70) &&(ix<=210) 
				&&(iy>=70) &&(iy<=210))
			{
				++chanid;
				staname=virtual_station_name(ix,iy);
				dbaddv(db,0,
					"sta",staname.c_str(),
					"chan",chan.c_str(),
					"time",starttime,
					"endtime",endtime,
					"chanid",chanid,
					"nsamp",nt,
					"samprate",5.0,
					"calib",1.0,
					"datatype","t4",
					"dir",dir.c_str(),
					"dfile",dfile.c_str(),
					"foff",foff,0);
				dbaddv(dbsite,0,
					"sta",staname.c_str(),
					"ondate",1960001,
					"lat",deg(grid.lat(ix,iy)),
					"lon",deg(grid.lon(ix,iy)),
					"elev",0.0,0);
				switch (component)
				{
				case 1:
				  dbaddv(dbsitechan,0,
					"sta",staname.c_str(),
					"chan",chan.c_str(),
					"chanid",chanid,
					"hang",0.0,
					"vang",90.0,
					"ondate",1970001,
					"edepth",0.0,0);
				   break;
				case 2:
				  dbaddv(dbsitechan,0,
					"sta",staname.c_str(),
					"chan",chan.c_str(),
					"chanid",chanid,
					"hang",90.0,
					"vang",90.0,
					"ondate",1970001,
					"edepth",0.0,0);
				   break;
				case 3:
				  dbaddv(dbsitechan,0,
					"sta",staname.c_str(),
					"chan",chan.c_str(),
					"chanid",chanid,
					"hang",0.0,
					"vang",0.0,
					"ondate",1970001,
					"edepth",0.0,0);
				   dbarrival.record=dbaddv(dbarrival,0,
					"sta",staname.c_str(),
					"chan",chan.c_str(),
					"chanid",chanid,
					"time",atime,
					"iphase","P",
					"deltim",0.1,0);
				   dbgetv(dbarrival,0,"arid",&arid,0);
				   dbaddv(dbassoc,0,
					"arid",arid,
					"orid",1,
					"sta",staname.c_str(),
					"phase","P",0);

				   break;
				}

			}
			}
			cout << "Finished row "<<ix<<endl;
		}
	}
	db=dblookup(db,0,"event",0,0);
	dbaddv(db,0,"evid",1,"prefor",1,0);
	db=dblookup(db,0,"origin",0,0);
	dbaddv(db,0,"orid",1,"evid",1,
		"lat",-deg(hypo.lat)+1.0,
		"lon",-deg(hypo.lon),
		"depth",hypo.z, 
		"time",hypo.time,0);
}
				
				
	
	

