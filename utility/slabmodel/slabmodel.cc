#include <iostream>
#include <fstream>
#include <string>
#include "coords.h"
#include "Metadata.h"
#include "seispp.h"
#include "interpolator1d.h"
#include "GeoTriMeshSurface.h"
#include "PLGeoPath.h"
#include "PlateBoundaryPath.h"
/* this is used several times here.  Beware if you ever cut and paste from this file */
const double REARTH(6378.17);
using namespace std;
using namespace SEISPP;
using namespace INTERPOLATOR1D;
const string prog("slabmodel");
typedef list<Geographic_point> Profile;
typedef vector<Profile> ProfileEnsemble ;
vector<Geographic_point> load_geopointdata(string fname)
{
    const string base_error("load_geopointdata:  ");
    FILE *fp=fopen(fname.c_str(),"r");
    if(fp==NULL) throw SeisppError(base_error
            + "Open failed on file="+fname);
    /* Assume input is decimal degrees lon,lon,depth */
    Geographic_point gp;
    vector<Geographic_point> points;
    double dlat,dlon,depth;
    while(fscanf(fp,"%lf%lf%lf",&dlat,&dlon,&depth)!=EOF)
    {
        if(dlat<-90.0 || dlat>90.0) throw SeisppError(base_error
                + "Inconsistent latitude value must be -90 to 90");
        if(dlon<-180.0 || dlon>360.0) throw SeisppError(base_error
                + "Inconsistent latitude value must be -180 to 360");
        gp.lat=rad(dlat);
        gp.lon=rad(dlon);
        gp.r=r0_ellipse(gp.lat)-depth;
        points.push_back(gp);
    }
    /* PLGeoPath allows a more general origin but here we intentionally 
       force it to first point in list.  Allow azimuth to default to 0 */
    return(points);
}
PLGeoPath resample_PLGeoPath(PLGeoPath& raw,double ds)
{
    if(ds<=0.0) 
        throw SeisppError("resample_PLGeoPath was passed a negative resample interval");
    double s0=raw.sbegin();
    double smax=raw.send();
    Geographic_point gp;
    gp=raw.origin();
    vector<Geographic_point> newpts;
    newpts.push_back(gp);
    double s;
    s=s0+ds;
    while(s<smax)
    {
        gp=raw.position(s);
        newpts.push_back(gp);
        s+=ds;
    }
    return(PLGeoPath(newpts,0));
}

void usage()
{
	cerr << prog <<" [-pf pffile -V]"<<endl
		<< "outputs grid coordinates to stdout";
	exit(-1);
}
bool SEISPP::SEISPP_verbose(false);
int main(int argc, char **argv)
{
	int i;
        ios::sync_with_stdio();
	string pfname(prog);
	for(i=1;i<argc;++i)
	{
		string sarg(argv[i]);
		if(sarg=="-pf")
		{
			++i;
			if(i>=argc) usage();
			pfname=string(argv[i]);
		}
                if(sarg=="-V")
                {
                    SEISPP_verbose=true;
                }
		else
			usage();
	}
	Pf *pf;
	if(pfread(const_cast<char *>(pfname.c_str()),&pf))
	{
		cerr << "pfread failed for pf file="<<pfname<<endl;
		usage();
	}
	try {
	    Metadata control(pf);
	    double plat=control.get_double("pole_latitude");
	    double plon=control.get_double("pole_longitude");
            /* This must have units of radians/year*/
            double angularvelocity=control.get_double("slab_angular_velocity");
            /* The next two  must have time in Mya*/
            double timesampleinterval=control.get_double("time_sample_interval");
            /* Time to run slab motion */
            double modeltime=control.get_double("model_elapsed_time");

	    plat=rad(plat);
	    plon=rad(plon);
            string trenchlinefile=control.get_string("trench_line_filename");
            /* This loads the raw digitized trench position data.  
               could be deleted after resampling, but assume this a 
               trivial memory issue */
            vector<Geographic_point> rawtrenchpoints
                                =load_geopointdata(trenchlinefile);
            PLGeoPath *rawtrenchpath=new PLGeoPath(rawtrenchpoints,0);
            double trench_path_sample_interval=control.get_double(
                    "trench_path_sample_interval");
            PLGeoPath trenchposition(resample_PLGeoPath(*rawtrenchpath,
                    trench_path_sample_interval));
            delete rawtrenchpath;
            /* now load and build the surface to be used to build the 
               model.   Note for simplicity assume the trench path data is
               part of the input set of points */
            string slabdata_filename=control.get_string("slabdata_filename");
            vector<Geographic_point> slabdata=load_geopointdata(slabdata_filename);
            GeoTriMeshSurface slabtrisurf(slabdata);
            /* Now create one path for each point on resampled trench path.*/
            int npoints=SEISPP::nint(modeltime/timesampleinterval)+1;
            int npaths=trenchposition.number_points();
            double Reffective;  // distance in km between pole and origin point 
            for(i=0;i<npaths;++i)
            {
                double s=trench_path_sample_interval*static_cast<double>(i);
                Geographic_point gp=trenchposition.position(s);
                PlateBoundaryPath pbpath(plat,plon,gp.lat,gp.lon,
                        angularvelocity);
                double aztmp;
                dist(plat,plon,gp.lat,gp.lon,&Reffective,&aztmp);
                Reffective=deg2km(deg(Reffective));
                const int oversampling(50);
                double sdt=timesampleinterval/static_cast<double>(oversampling);
                int j;
                vector<Geographic_point> oversampledpath;
                for(j=0;j<(npoints*oversampling+1);++j)
                {
                    gp=pbpath.position(static_cast<double>(j)*sdt);
                    gp.r=slabtrisurf.radius(gp.lat,gp.lon);
                    // Allow the first point to be skipped but issue warning
                    if(is_nan(gp.r)) 
                    {
                        if(j==0)
                        {
                            cerr << "Warning:  origin point not inside convex hull. "<<endl
                                << "A skew of about "<<1.0/static_cast<double>(oversampling)
                                << " of the time sampling rate will be present"<<endl;
                            continue;
                        }
                        else
                            break;
                    }
                    oversampledpath.push_back(gp);
                }
                PLGeoPath plop(oversampledpath,0);
                double sampleintervalkm=Reffective*timesampleinterval*angularvelocity;
                PLGeoPath finalpath=resample_PLGeoPath(plop,sampleintervalkm);
                cout << finalpath;
                cout <<">"<<endl;
                oversampledpath.clear();
            }
        }

        catch (SeisppError& serr)
        {
            serr.log_error();
            exit(-2);
        }
        catch (int ierr)
        {
            cerr << "GCLgrid error:  something threw error code="
                <<ierr<<endl;
            exit(-2);
        }
}
