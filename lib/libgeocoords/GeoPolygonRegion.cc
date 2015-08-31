#include <fstream>
#include <sstream>
#include "GeoPolygonRegion.h"
#include <math.h>
#define deg(r)    ((r) * 180.0/M_PI)
using namespace std;

/* These routines were taken from antelope contrib libttregions */

/* Reimplementation of winding-number algorithm
 Godkin and Pulli, BSSA v.74 1845-1948, 1984.

 K. Lindquist
 Geophysical Institute
 University of Alaska
 1999
*/

double
signed_crossing_number( double *segment )
{
	double	x1, y1, x2, y2;
	double	slope;
	double 	xintercept;
	double	direction;

	x1 = segment[0];
	y1 = segment[1];
	x2 = segment[2];
	y2 = segment[3];

	if( y1 * y2 > 0 )
	{
		/* No crossing--both points on same side of y axis */
		return 0;
	}

	/* Horizontal lines: */
	if( y1 == 0 && y2 == 0 )
	{
		if( x1 * x2 > 0 )
		{
			/* No crossing */
			return 0;
		}
		else
		{
			/* Segment crosses or touches origin */
			return 2;
		}
	}

	if( x2 == x1 )
	{
		/* Vertical lines */
		xintercept = x1;
	}
	else 
	{
		slope = ( y2 - y1 ) / (x2 - x1 );
		xintercept = x1 - y1 / slope;
	}

	if( y2 > y1 )
	{
		direction = 1;
	}
	else 
	{
		direction = -1;
	}	

	if( xintercept > 0 )
	{
		return 0;

	}
	else if( xintercept == 0 ) 
	{
		
		return 2;
	}
	else if( y1 == 0 || y2 == 0 )
	{

		return 0.5 * direction;
	} else
	{
		return direction;
	}
}

double *
segment_number( int seg_index, double *shifted_polygon, int ncoords )
{
	int     nsegs;
	double	*segment;

	nsegs = ncoords / 2;

	if( seg_index > nsegs )
	{
		fprintf( stderr, "Bad segment number\n" );
		return NULL;
	}

	segment = (double *) malloc( 4 * sizeof( double ) );

	if( seg_index == nsegs )
	{
		segment[0] = shifted_polygon[ncoords-2];
		segment[1] = shifted_polygon[ncoords-1];
		segment[2] = shifted_polygon[0];
		segment[3] = shifted_polygon[1];
	} else {
		segment[0] = shifted_polygon[(seg_index-1)*2];
		segment[1] = shifted_polygon[(seg_index-1)*2+1];
		segment[2] = shifted_polygon[(seg_index-1)*2+2];
		segment[3] = shifted_polygon[(seg_index-1)*2+3];
	}

	return segment;
}

double
winding_number( double *shifted_polygon, int ncoords )
{
	double	winding_num = 0;
	double	*segment;
	int	nsegs;
	int	seg_index;

	nsegs = ncoords / 2;
	
	for( seg_index = 1; seg_index <= nsegs; seg_index++ ) 
	{
		segment = segment_number( seg_index, shifted_polygon, ncoords );
		winding_num += signed_crossing_number( segment );
		free( segment );
	}

	return winding_num;
}

/* This is a common error test for problem points */
bool points_not_valid(vector<double> vlat, vector<double> vlon)
{
    int npts=vlat.size();  // assume lat and lon same
    /* First look for points outside range of -180 to 360 for longitude */
    int i;
    for(i=0;i<npts;++i)
    {
        if(vlon[i]>360.0) return true;
        if(vlon[i]<-180.0) return true;
        if(vlat[i]>90.0) return true;
        if(vlat[i]<-90.0) return true;
    }
    /* Look for large jumps */
    const double maxchange(90.0);
    for(i=1;i<npts;++i)
    {
        if(fabs(vlat[i]-vlat[i-1])>maxchange) return true;
        if(fabs(vlon[i]-vlon[i-1])>maxchange) return true;
    }
    return false;
}
GeoPolygonRegion::GeoPolygonRegion()
{
    npoints=0;
}
GeoPolygonRegion::GeoPolygonRegion(string fname)
{
    string base_error("GeoPolygonRegion file constructor:  ");
    ifstream in;
    in.open(fname.c_str(),ios::in);
    if(in.fail())
        throw GCLgridError(base_error
                + "Open failed for file="+fname);
    char inputline[256];
    double latin,lonin;
    lat.clear();   lon.clear();
    while(in.getline(inputline,256))
    {
        stringstream ss(inputline);
        if(inputline[0]=='#') continue;
        if(inputline[0]=='>') break;
        ss >> lonin;
        ss >> latin;
        lat.push_back(latin);
        lon.push_back(lonin);
    }
    if(points_not_valid(lat,lon))throw GCLgridError(base_error
            + "illegal lat or lon values in polygon\n"
            +"Does not support longitude singularity or across poles");
    /* Add a new point if the first and last points are not within this
       (frozen) constant distance.  Crude approach should be made more
       generic if this code is ever released */
    const double deltadegmax(0.001);
    double deltadeg,d2deg;
    npoints=lat.size();
    deltadeg=lat[npoints-1]-lat[0];
    d2deg=deltadeg*deltadeg;
    deltadeg=lon[npoints-1]-lon[0];
    d2deg+=deltadeg*deltadeg;
    if(sqrt(d2deg)>deltadegmax)
    {
        lat.push_back(lat[0]);
        lon.push_back(lon[0]);
    }
    /* Last but not least convert all lat and lon values to radians */
    npoints=lat.size();
    for(int i=0;i<npoints;++i)
    {
        lat[i]=rad(lat[i]);
        lon[i]=rad(lon[i]);
    }
}
GeoPolygonRegion::GeoPolygonRegion(vector<double> vlat, vector<double> vlon)
    : lat(vlat),lon(vlon)
{
    const string base_error("GeoPolygonRegion constructor:  ");
    if(vlat.size()!=vlon.size())
    {
        throw GCLgridError(base_error
                + "lat and lon input vectors have different sizes");
    }
    if(points_not_valid(lat,lon))
        throw GCLgridError(base_error + "illegal lat or lon values in polygon\n"
            +   "Does not support longitude singularity or across poles");
}
GeoPolygonRegion::GeoPolygonRegion(vector<Geographic_point> gpts)
{
    int i;
    const string base_error("GeoPolygonRegion vector<Geographic_point> constructor");
    int npts=gpts.size();
    lat.reserve(npts);
    lon.reserve(npts);
    /* Temporaries for error check */
    vector<double> dlat,dlon;
    dlat.reserve(npts);
    dlon.reserve(npts);
    for(i=0;i<npts;++i)
    {
        /* Geographic points are normally in radians so we create
           the dlat and dlon (degrees) version to allow use of points_not_valid 
           test */
        dlat.push_back(deg(gpts[i].lat));
        dlon.push_back(deg(gpts[i].lon));
        lat.push_back(gpts[i].lat);
        lon.push_back(gpts[i].lon);
    }
    if(points_not_valid(dlat,dlon))
        throw GCLgridError(base_error + "illegal lat or lon values in polygon\n"
                "Does not support longitude singularity or across poles");
}


bool GeoPolygonRegion::is_inside(double lat0, double lon0)
{
    vector<double> shifted_polygon;
    /* In Kent's function these are interleaved so we make it that way */
    shifted_polygon.reserve(npoints*2);
    int i;
    double slat,slon;
    for(i=0;i<npoints;++i)
    {
        /* lon is x and y is lat which for winding means this order */
        slon=lon[i]-lon0;
        slat=lat[i]-lat0;
        shifted_polygon.push_back(slon);
        shifted_polygon.push_back(slat);
    }
    int wn = winding_number(&(shifted_polygon[0]),2*npoints);
    if(wn!=0)
        return(true);
    else
        return(false);
}
GeoPolygonRegion::GeoPolygonRegion(const GeoPolygonRegion& p)
{
    npoints=p.npoints;
    lat=p.lat;
    lon=p.lon;
}
GeoPolygonRegion& GeoPolygonRegion::operator=(const GeoPolygonRegion& p)
{
    if(this!=&p)
    {
        npoints=p.npoints;
        lat=p.lat;
        lon=p.lon;
    }
    return(*this);
}

ostream& operator << (ostream& os, GeoPolygonRegion& p)
{
    int i;
    for(i=0;i<p.npoints;++i)
        os << deg(p.lon[i])<<" "<<deg(p.lat[i])<<endl;
    return(os);
}
