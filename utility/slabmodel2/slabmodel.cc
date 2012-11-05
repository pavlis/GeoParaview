#include <iostream>
#include <fstream>
#include <string>
#include "coords.h"
#include "Metadata.h"
#include "seispp.h"
#include "interpolator1d.h"
#include "GeoPolygonRegion.h"
#include "GeoTriMeshSurface.h"
#include "GeoSplineSurface.h"
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
bool pointsmatch(Geographic_point p1, Geographic_point p2)
{
    /* Roughly float_epsilon.  Very conservative, could cause
    problems if really close points, but for expected use of
    this program that would be ridiculous anyway */
    const double cutoff(0.0000001);
    double test;
    test=fabs((p1.lat-p2.lat)/p1.lat);
    if(test>cutoff) return false;
    test=fabs((p1.lon-p2.lon)/p1.lon);
    if(test>cutoff) return false;
    return true;
}

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
    int count=0;
    Geographic_point lastpoint;
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
        if(count==0)
            points.push_back(gp);
        else
        {
            /* This is necessary because it causes all the interpolators
               to fail entering an infinite loop */
            if(pointsmatch(gp,lastpoint))
            {
                cerr << "Warning:  duplicate point at "
                    << deg(gp.lat) <<" " << deg(gp.lon)
                    << " at line "<<count<<" dropped"<<endl;
            }
            else
            {
                points.push_back(gp);
            }
        }
        lastpoint=gp;
        ++count;
    }
    /* PLGeoPath allows a more general origin but here we intentionally 
       force it to first point in list.  Allow azimuth to default to 0 */
    return(points);
}
PLGeoPath timesample_PLGeoPath(PLGeoPath& raw, 
        vector<double>t,vector<double> s,double dt,double endtime)
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
    /* terminate on either end of the path or the passed endtime 
       argument*/
    double etest=endtime+(0.01*dt);  //allow some slop for this test
    while(t0<tmax && t0<=etest)
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
	cerr << prog <<" [-tjp -pf pffile -v]"<<endl
		<< "outputs grid coordinates to stdout"<<endl
                << "-tjp option outputs triple junction path"
                << " (error if triple junction mode off)"<<endl;;
	exit(-1);
}
bool SEISPP::SEISPP_verbose(false);
int main(int argc, char **argv)
{
	int i;
        ios::sync_with_stdio();
	string pfname(prog);
        bool output_triple_junction_path(false);
	for(i=1;i<argc;++i)
	{
		string sarg(argv[i]);
		if(sarg=="-pf")
		{
			++i;
			if(i>=argc) usage();
			pfname=string(argv[i]);
		}
                else if(sarg=="-v")
                {
                    SEISPP_verbose=true;
                }
                else if(sarg=="-tjp")
                {
                    output_triple_junction_path=true;
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
            /* Number of points to oversample defining paths.
               Most relevant to spline methods.   Makes larger 
               with splnes and coarse sampling.  Make smaller for 
               triangulation and/or short time sampling */
            int oversampling=control.get_int("oversampling_count");
            double trench_path_sample_interval=control.get_double(
                        "trench_path_sample_interval");
            double maxdip,mindip;
            maxdip=control.get_double("maximum_extension_dip");
            mindip=control.get_double("minimum_extension_dip");
            bool use_local_dip=control.get_bool("use_local_dip");
            bool triple_mode=control.get_bool("triple_junction_mode");
            if(!triple_mode && output_triple_junction_path)
            {
                cerr << "Illegal option combination.  -tjp option cannot "
                    << "be used unles triple_junction_mode is true"<<endl;
                usage();
            }
            PLGeoPath zerotimecurve;
            /* When using triple junction mode derive the zero time line from
               plate boundary path.  Otherwise interpolate digitized line data.
               This if block probably should be a procedure from a style standpoint. */
            if(triple_mode)
            {
                string spfile=control.get_string("stage_pole_file");
                double lat0,lon0;
                lat0=control.get_double("zero_line_origin_latitude");
                lon0=control.get_double("zero_line_origin_longitude");
                lat0=rad(lat0);
                lon0=rad(lon0);
                double z0=control.get_double("zero_line_depth");
                TimeVariablePlateBoundaryPath tvpzp(spfile,lat0,lon0);
                double zldt=control.get_double("zero_line_time_sample_interval");
                vector<Geographic_point> zeroraw;
                double tzl(0.0);
                double modtimemya=modeltime/1000000.0;
                while(tzl<modtimemya)
                {
                    Geographic_point zgp;
                    zgp=tvpzp.position(tzl);
                    /* The path is at the reference ellipsoid.  Subtract the refence
                       depth from each point */
                    zgp.r-=z0;
//DEBUG
//cerr << deg(zgp.lat)<<" " << deg(zgp.lon) << " " << zgp.r 
//    <<" "<<tzl <<endl;
                    zeroraw.push_back(zgp);
                    tzl+=zldt;
                }
                zerotimecurve=PLGeoPath(zeroraw,0);
                if(output_triple_junction_path)
                {
                    cout << zerotimecurve;
                    exit(1);
                }
//DEBUG
//cerr << "PLGeoPath version of zero time curve"<<endl;
//cout << zerotimecurve;
            }
            else
            {
                string trenchlinefile=control.get_string("trench_line_filename");
                vector<Geographic_point> rawtrenchpoints
                                    =load_geopointdata(trenchlinefile);
                PLGeoPath *rawtrenchpath=new PLGeoPath(rawtrenchpoints,0);
                zerotimecurve=resample_PLGeoPath(*rawtrenchpath,
                        trench_path_sample_interval);
                delete rawtrenchpath;
            }

            /* now load and build the surface to be used to build the 
               model.   Note for simplicity assume the trench path data is
               part of the input set of points */
            string slabdata_filename=control.get_string("slabdata_filename");
            vector<Geographic_point> slabdata=load_geopointdata(slabdata_filename);
            /* This oddity is needed to deal with polymophic geosurf */
            GeoSplineSurface *gss;
            GeoTriMeshSurface *gts;
            GeoSurface *geosurf;
            bool spline_surface=control.get_bool("use_bicubic_spline");
            if(spline_surface)
            {
                gss=new GeoSplineSurface(slabdata,control);
                geosurf=dynamic_cast<GeoSurface *>(gss);
            }
            else
            {
                gts=new GeoTriMeshSurface(slabdata);
                geosurf=dynamic_cast<GeoSurface *>(gts);
            }
            /* This adds a bounding curve that may not be convex */
            bool use_bounding_polygon=control.get_bool("use_bounding_polygon");
            if(use_bounding_polygon)
            {
                string boundsfile=control.get_string("bounding_polygon_file");
                GeoPolygonRegion bounds(boundsfile);
                geosurf->AddBoundary(bounds);
            }
            /* Now create one path for each point on resampled trench path.*/
            int npoints=SEISPP::nint(modeltime/timesampleinterval)+1;
            int npaths=zerotimecurve.number_points();
            if(SEISPP_verbose) cerr << "Output grid npaths="<<npaths<<endl;
            //double Reffective;  // distance in km between pole and origin point 
            bool extendpaths=control.get_bool("extendpaths");
            GeoPath *pbpath;
            for(i=0;i<npaths;++i)
            {
                double s;
                Geographic_point gp;
                if(triple_mode)
                {
                    gp=zerotimecurve.node_position(i);
                }
                else
                {
                    s=trench_path_sample_interval*static_cast<double>(i);
                    gp=zerotimecurve.position(s);
                }
//DEBUG
//cerr << "Starting position:  "<<deg(gp.lat)<<" "<<deg(gp.lon)<<endl;
                pbpath=BuildPBPObject(gp.lat,gp.lon, pf);
                /*
                PlateBoundaryPath pbpath(plat,plon,gp.lat,gp.lon,
                        angularvelocity);
                double aztmp;
                dist(plat,plon,gp.lat,gp.lon,&Reffective,&aztmp);
                Reffective=deg2km(deg(Reffective));
                        */
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
                double current_time(0.0),current_s(0.0),ds;
                double dzdx;
                if(SEISPP_verbose) cerr <<"Oversample path length="
                    <<npoints*oversampling<<endl;
                /* Two counters used here.  j is the live index for
                 vectors.  jloop is the total loop counter
                 They are different only if origin outside convex
                 hull error (jloop=0 conditional) is entered */
                for(int jloop=0,j=0;jloop<(npoints*oversampling+1);++j,++jloop)
                {
                    gp=pbpath->position(static_cast<double>(jloop)*sdt);
                    if(geosurf->is_defined(gp.lat,gp.lon))
                    {
                        gp.r=geosurf->radius(gp.lat,gp.lon);
                    }
                    else
                    {
                    // Allow the first point to be skipped but issue warning
                        if(jloop==0)
                        {
                            cerr << "Warning:  origin point not inside convex hull. "<<endl
                                << "A skew of about "<<1.0/static_cast<double>(oversampling)
                                << " of the time sampling rate will be present"<<endl;
                            lastgp=gp;
                            --j;
                            continue;
                        }
                        else if(j<2)
                            break;
                        else if(extendpaths)
                        {

                            // This computes change in depth (sign switch)
                            dzdx=oversampledpath[j-2].r
                                        - oversampledpath[j-1].r;
                            /* compute distance between j-1 and j-2 to compute dip */
                            double ddelta=GeoDistance(oversampledpath[j-2],oversampledpath[j-1]);
                            dzdx/= ddelta;  // note mixed units.  dz in km, ddelta in radians
                            double dipdeg;
                            double R0(oversampledpath[j-1].r);
                            dipdeg=dzdx/R0;
                            dipdeg=deg(atan(dipdeg));
                            /* These could be  precomputed for efficiency but better to leave
                               it here for clarity */
                            if(dipdeg>maxdip)
                            {
                                dzdx=tan(rad(maxdip))*R0;
                                dipdeg=maxdip;
                            }
                            else if(dipdeg<mindip)
                            {
                                dzdx=tan(rad(mindip))*R0;
                                dipdeg=mindip;
                            }
                            int jlast=j-1;
                            if(SEISPP_verbose) 
                                cerr << "Extending path "<<i<<" with dip "
                                <<dipdeg<<" from point number "<<j<<endl
                                <<"Position = " 
                                <<deg(oversampledpath[jlast].lat)<<", "
                                <<deg(oversampledpath[jlast].lon)<<endl;
                            // now extend the path to npoints
                            int jj;
                            double s0=path_s[jlast];
                            double r0=oversampledpath[jlast].r;
                            lastgp=oversampledpath[jlast-1];
                            for(jj=j;jj<(npoints*oversampling+1);++jj)
                            {
                                double s;
                                s=static_cast<double>(jj)*sdt;
                                gp=pbpath->position(s);
                                // recycle ddelta for same context here
                                ddelta=GeoDistance(lastgp,gp);
                                /* This option corrects delta for 
                                   shrinking length with depth */
                                if(use_local_dip) ddelta *= lastgp.r/R0;
                                gp.r=lastgp.r-dzdx*ddelta;
                    //DEBUG
                                /*
                    cout << "Extended: "
                        << deg(gp.lon)<<" "
                        << deg(gp.lat)<<" "
                        <<r0_ellipse(gp.lat)-gp.r<<endl;
                        */
                                oversampledpath.push_back(gp);
                                current_time += adjustedtime(lastgp,gp,sdt);
                                ds = distance_increment(lastgp,gp,sdt);
                                current_s += ds;
                                corrected_time.push_back(current_time);
                                path_s.push_back(current_s);
                                lastgp=gp;
                            }
                            break;
                        }
                        else
                            break;
                    }
                    if(j>0)
                    {
                        current_time += adjustedtime(lastgp,gp,sdt);
                        ds = distance_increment(lastgp,gp,sdt);
                        current_s += ds;
                    }
                    corrected_time.push_back(current_time);
                    path_s.push_back(current_s);
                    //DEBUG
                    /*
                    cout << "Interpolated: "
                        << deg(gp.lon)<<" "
                        << deg(gp.lat)<<" "
                        <<r0_ellipse(gp.lat)-gp.r<<endl;
                        */

                    oversampledpath.push_back(gp);
                    lastgp=gp;
                }
                if(oversampledpath.size()>1)
                {
                    //DEBUG
                    /*
                    for(int kd=0;kd<oversampledpath.size();++kd)
                        cout << path_s[kd]<<" "
                            << corrected_time[kd]<<" "
                            << deg(oversampledpath[kd].lat)<<" "
                            << deg(oversampledpath[kd].lon)<<" "
                            << oversampledpath[kd].r
                            <<" "
                            << r0_ellipse(oversampledpath[kd].lat)
                              - oversampledpath[kd].r<<endl;
                              */
                    PLGeoPath plop(oversampledpath,0);
                    PLGeoPath finalpath=timesample_PLGeoPath(plop,
                            corrected_time,path_s,timesampleinterval,
                            modeltime);
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
        catch (GeoCoordError& gerr)
        {
            cerr << gerr.what()<<endl;
            exit(-3);
        }
        catch (int ierr)
        {
            cerr << "GCLgrid error:  something threw error code="
                <<ierr<<endl;
            exit(-2);
        }
}
