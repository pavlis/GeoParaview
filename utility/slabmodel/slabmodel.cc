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
    char line[128];
    //while(fscanf(fp,"%lf%lf%lf",&dlat,&dlon,&depth)!=EOF)
    while(fgets(line,128,fp)!=NULL)
    {
        if(line[0]=='#') continue;
        sscanf(line,"%lf%lf%lf",&dlon,&dlat,&depth);
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
PLGeoPath timesample_PLGeoPath(PLGeoPath& raw, 
        vector<double>t,vector<double> s,double dt)
{
    if(t.size()!=s.size()) throw SeisppError(string("timesample_PLGeoPath:  ")
            + "time and distance vector sizes do not match");
    double s0=raw.sbegin();
    double smax=raw.send();
    Geographic_point gp;
    vector<Geographic_point> newpts;
    int nt=t.size();
    if(nt<=1) throw SeisppError(string("timesample_PLGeoPath:  ")
            + "empty vectors for time and distance for path");
    double t0,tmax;
    tmax=t[nt-1];
    //Assume first point is time t=0;
    gp=raw.origin();
    newpts.push_back(gp);
    t0=dt;
    while(t0<tmax)
    {
        // A hideously inefficient linear search, but don't expect
        // nt to ever be huge
        int i;
        for(i=1;i<nt;++i)
            if(t[i]>t0) break;
        double dsdt=(s[i]-s[i-1])/(t[i]-t[i-1]);
        double sp=s[i-1]+(t0-t[i-1])*dsdt;
        gp=raw.position(sp);
        newpts.push_back(gp);
        t0+=dt;
    }
    return(PLGeoPath(newpts,0));
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
GeoPath *BuildPBPObject(double olat, double olon, Pf *pf)
{
    Tbl *t;
    string tname("pole_data");
    t=pfget_tbl(pf,const_cast<char *>(tname.c_str()));
    vector<double>spla,splo,dt,ang;
    /* Read data in gmt rotconverter and backtracker format for stage pole
       data.  (This program requires stage pole input).  That format requires
       stage poles in reverse time order.  Hence we have to read the stage pole
       data in and then reverse the order to produce a path oriented from
       time 0 to some time in the past. */
    for(int i=0;i<maxtbl(t);++i)
    {
        char *line;
        line=(char *)gettbl(t,i);
        stringstream ss(line);
        double lat,lon,timestart,timeend,phi;
        ss>>lon;
        ss>>lat;
        ss>>timeend;
        ss>>timestart;
        ss>>phi;
        lat=rad(lat);
        lon=rad(lon);
        phi=rad(phi);
        spla.insert(spla.begin(),lat);
        splo.insert(splo.begin(),lon);
        dt.insert(dt.begin(),(timeend-timestart)*1000000.0);
        /* These are stage pole rotation angles */
        ang.insert(ang.begin(),phi);
    }
    int npoles=spla.size();
    if(npoles==1)
    {
        PlateBoundaryPath *pbpath=new PlateBoundaryPath(spla[0],splo[0],
                olat,olon,ang[0]/dt[0]);
        return pbpath;
    }
    else
    {
        TimeVariablePlateBoundaryPath *pbpath=new TimeVariablePlateBoundaryPath(
                spla,splo,dt,ang,olat,olon);
        return pbpath;
    }
}
/*  Compute and return great circle distance between two points
    defined by two Geographic_point objects. Returned distance 
    is in radians. */
double GeoDistance(Geographic_point gp0, Geographic_point gp1)
{
    double delta,az;
    dist(gp0.lat,gp0.lon,gp1.lat,gp1.lon,&delta,&az);
    return(delta);
}

/* This pair of functions are used to computed 3d distance and 
   a corrected time increment for paths.  There is duplication of
   the distance calculation, which is not very efficient, but 
   acceptable as long as as this isn't used enormous numbers
   of times.

In both procedures gp0 is first point and gp1 is second.  Significant
only in how horizontal distance is computed as r*delta.  */
double distance_increment(Geographic_point gp0, Geographic_point gp1, 
        double dt)
{
    double delta;
    delta=GeoDistance(gp0,gp1);
    delta*=gp0.r;
    double ds=hypot(delta,fabs(gp1.r-gp0.r));
    return(ds);
}
double adjustedtime(Geographic_point gp0, Geographic_point gp1, 
        double dt)
{
    double delta;
    delta=GeoDistance(gp0,gp1);
    delta*=r0_ellipse(gp0.lat);
    double ds=hypot(delta,fabs(gp1.r-gp0.r));
    return(dt*ds/delta);
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
	    //double plat=control.get_double("pole_latitude");
	    //double plon=control.get_double("pole_longitude");
            /* This must have units of radians/year*/
            //double angularvelocity=control.get_double("slab_angular_velocity");
            /* The next two  must have time in Mya*/
            double timesampleinterval=control.get_double("time_sample_interval");
            /* Time to run slab motion */
            double modeltime=control.get_double("model_elapsed_time");

	    //plat=rad(plat);
	    //plon=rad(plon);
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
            if(SEISPP_verbose) cerr << "Output grid npaths="<<npaths<<endl;
            //double Reffective;  // distance in km between pole and origin point 
            GeoPath *pbpath;
            for(i=0;i<npaths;++i)
            {
                double s=trench_path_sample_interval*static_cast<double>(i);
                Geographic_point gp=trenchposition.position(s);
                pbpath=BuildPBPObject(gp.lat,gp.lon, pf);
                /*
                PlateBoundaryPath pbpath(plat,plon,gp.lat,gp.lon,
                        angularvelocity);
                double aztmp;
                dist(plat,plon,gp.lat,gp.lon,&Reffective,&aztmp);
                Reffective=deg2km(deg(Reffective));
                        */
                const int oversampling(50);
                double sdt=timesampleinterval/static_cast<double>(oversampling);
                int j;
                vector<Geographic_point> oversampledpath;
                /* These are needed if we have to extrapolate the last valid point */
                Geographic_point lastgp;
                /* This vector holds path times corrected for radial distance.
                   This is necessary to create a time grid.*/
                vector<double> corrected_time;
                /* Parallel vector of arc distances */
                vector<double> path_s;
                // set a zero for first of each of these 
                path_s.push_back(0.0);
                corrected_time.push_back(0.0);
                double current_time,current_s;
                double dzdx;
                if(SEISPP_verbose) cerr <<"Oversample path length="
                    <<npoints*oversampling<<endl;
                for(j=0,current_time=0.0,current_s=0.0;
                        j<(npoints*oversampling+1);++j)
                {
                    gp=pbpath->position(static_cast<double>(j)*sdt);
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
                        else if(j<2)
                            break;
                        else
                        {
                            // This computes change in depth (sign switch)
                            dzdx=oversampledpath[j-2].r
                                        - oversampledpath[j-1].r;
                            /* compute distance between j-1 and j-2 to compute dip */
                            double ddelta=GeoDistance(oversampledpath[j-2],oversampledpath[j-1]);
                            dzdx/= ddelta;  // note mixed units.  dz in km, ddelta in radians
                            if(SEISPP_verbose) cerr << "Extending path with dip "
                                <<dzdx<<" from point number "<<j<<endl;
                            // now extend the path to npoints
                            int jlast=j-1;
                            int jj;
                            double s0=path_s[jlast];
                            double r0=oversampledpath[jlast].r;
                            for(jj=j;jj<(npoints*oversampling+1);++jj)
                            {
                                double s;
                                s=static_cast<double>(jj)*sdt;
                                gp=pbpath->position(static_cast<double>(jj)*sdt);
                                // recycle ddelta for same context here
                                ddelta=GeoDistance(lastgp,gp);
                                gp.r=lastgp.r-dzdx*ddelta;
                                oversampledpath.push_back(gp);
                                current_time += adjustedtime(lastgp,gp,sdt);
                                current_s += distance_increment(lastgp,gp,sdt);
                                corrected_time.push_back(current_time);
                                path_s.push_back(current_s);
                                lastgp=gp;
                            }
                            break;
                        }
                    }
                    if(j>0)
                    {
                        current_time += adjustedtime(lastgp,gp,sdt);
                        current_s += distance_increment(lastgp,gp,sdt);
                    }
                    corrected_time.push_back(current_time);
                    path_s.push_back(current_s);
                    oversampledpath.push_back(gp);
                    lastgp=gp;
                }
                if(j>1)
                {
                    //DEBUG
                    /*
                    for(int kd=0;kd<oversampledpath.size();++kd)
                        cout << path_s[kd]<<" "
                            << corrected_time[kd]<<" "
                            << deg(oversampledpath[kd].lat)<<" "
                            << deg(oversampledpath[kd].lon)<<" "
                            << oversampledpath[kd].r<<endl;
                            */
                    PLGeoPath plop(oversampledpath,0);
                    PLGeoPath finalpath=timesample_PLGeoPath(plop,
                            corrected_time,path_s,timesampleinterval);
                    cout << finalpath;
                    cout <<">"<<endl;
                }
                else
                {
                    cerr << "Warning:  Path has zero length "
                        <<"for path point number "<<i<<endl;
                }
                oversampledpath.clear();
                corrected_time.clear();
                path_s.clear();
                delete pbpath;
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
