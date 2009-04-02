/* This program builds a path on the surface of the earth defined by a pole
of rotation, point of origin, and distance range from origin point.  The arc
produced as output can be drawn on a map as the location of a transform
fault through the origin parallel to the plate transport direction at 
that position.  */
#include <iostream>
#include <string>
#include "coords.h"
#include "Metadata.h"
#include "seispp.h"
using namespace std;
using namespace SEISPP;
void usage()
{
	cerr << "tfplateboundary [-pf pffile]"<<endl
		<< "Outputs lon,lat pairs to stdout"<<endl;
	exit(-1);
}
bool SEISPP::SEISPP_verbose(true);
int main(int argc, char **argv)
{
	int i;
	string pfname("tfplateboundary");
	if(argc>1)
	{
		for(i=1;i<argc;++i)
		{
			string sarg(argv[i]);
			if(sarg=="-pf")
			{
				++i;
				if(i>=argc) usage();
				pfname=string(argv[i]);
			}
			else
				usage();
		}
	}
	Pf *pf;
	if(pfread(const_cast<char *>(pfname.c_str()),&pf))
	{
		cerr << "pfread failed for pf file="<<pfname<<endl;
		usage();
	}
	try {
	Metadata control(pf);
	double plat,plon,olat,olon;
	plat=control.get_double("pole_latitude");
	plon=control.get_double("pole_longitude");
	olat=control.get_double("origin_latitude");
	olon=control.get_double("origin_longitude");
	/* These define line endpoints as distance in km from origin.
	Convention is that positive is clockwise rotation around pole.*/
	double lowendpoint,highendpoint;
	lowendpoint=control.get_double("start_distance_from_origin");
	highendpoint=control.get_double("end_distance_from_origin");
	int npath=control.get_int("number_of_points");
	/* convert all lat,lon values to radians */
	plat=rad(plat);
	plon=rad(plon);
	olat=rad(olat);
	olon=rad(olon);
	/* We need to compute the gcp distance from pole to origin */
	double delta, azimuth;
	dist(plat,plon,olat,olon,&delta,&azimuth);
	const double REARTH=6378.17;
	double reffective=REARTH*sin(delta);
	double daz,azrange;
	azrange=(highendpoint-lowendpoint)/reffective;
	daz=azrange/static_cast<double>(npath-1);
	/* compute the starting point azimuth from pole */
	double az0=lowendpoint/reffective;
	azimuth+=az0;
	/* Now compute the points and put them to output */
	for(i=0;i<npath;++i,azimuth+=daz)
	{
		double lat,lon;
		latlon(plat,plon,delta,azimuth,&lat,&lon);
		cout << deg(lon) << " "<<deg(lat)<<endl;
	}
	} catch (MetadataGetError mderr)
	{
		mderr.log_error();
		exit(-1);
	}
}
		

