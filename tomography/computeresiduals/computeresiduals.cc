#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <list>
#include "seispp.h"
#include "VelocityModel_1d.h"
#include "ray1d.h"
#include "dbpp.h"
using namespace std;
using namespace SEISPP;
/* use this to load data values and as a container to hold 
   computed values */
class ArrivalRecord
{
    public:
        long evid;
        string sta;
        double lat;   //station latitude (radians)
        double lon;  // station longitude (radians)
        // elevation intentionally ignored since assumed corrected
        // here with crustal correction program
        double time;
        double deltim;
        double timeres;  // residual stored in assoc
        // These are not loaded from db but computed later
        // This allows passing around event data as a vector of these objects
        double predtime;  // time computed using Hypocenter calculator (1d absolute t)
        SlownessVector u;  // computed slowness vector similar to predtime
        double t1d;  // ray trace time through model for 1d reference model
        double t3d;  // ray trace dtime for 3d model using ray computed for 1d mod
        double rawresid;  //computed as time-predtime: generally poorest estimae
        double resid3d; // computed residual with 3d model (timeres-t3d)
        // convenient to define this to output one record 
        friend ostream& operator<<(ostream& os, ArrivalRecord& ar);
};
ostream& operator<<(ostream& os,ArrivalRecord& ar)
{
    os << ar.evid <<" "
       << ar.sta <<" "
       << deg(ar.lat) <<" "
       << deg(ar.lon) <<" "
       << strtime(ar.time) <<" "
       << ar.deltim <<" "
       << ar.timeres <<" "
       << ar.predtime <<" "
       << ar.u.ux <<" "
       << ar.u.uy <<" "
       << ar.t1d <<" "
       << ar.t3d <<" "
       << ar.rawresid <<" "
       << ar.resid3d <<endl;
    return os;
}
/* This procedure is similar to libseispp routines called GCLgrid_Ray_project
   but is a bit more general.  That is the older code used coordinates of the
   grid as is generally needed for the plane wave migation code.  This 
   returns a path in the cartesian reference frame of grid with the input path
   projected according to the direction defined by the slowness vector u. 
   The path will originate from the point lat0,lon0, (both radians), and elev
   (in km).
   */
dmatrix RayPathProject(double lat0, double lon0, double elev0,
        GCLgrid3d& grid,RayPathSphere& path,SlownessVector& u)
{
    double r0=r0_ellipse(lat0)+elev0;
    Cartesian_point this_point;
    try {
        dmatrix pathout(3,path.npts);
        double az=u.baz();
        for(int i=0;i<path.npts;++i)
        {
            double lat,lon;
            latlon(lat0,lon0,path.delta[i],az,&lat,&lon);
            double radius=r0_ellipse(lat)-(path.r[0]-path.r[i]);
            this_point = grid.gtoc(lat,lon,radius);
            pathout(0,i)=this_point.x1;
            pathout(1,i)=this_point.x2;
            pathout(2,i)=this_point.x3;
        }
        return pathout;
    }catch(...){throw;};
}
void compute_residuals(ArrivalRecord& a,VelocityModel_1d& vm1d, 
        GCLscalarfield3d& vm3d,Hypocenter& h,string phase)
{
    // First computer the raw residual as iasp91 residual
    a.predtime=h.phasetime(a.lat,a.lon,0.0,phase);
    a.rawresid=a.time-a.predtime-h.time;
    a.u=h.phaseslow(a.lat,a.lon,0.0,phase);
    // assume origin point is max depth
    double zmax=vm3d.depth(0,0,0);
    double tmax(9999.0);  // effectively off since running always in constant z mode
    // This if frozen as 100 m for now.  Not at all general but should give 
    // more than sufficient approximation
    double dz(0.1);
    RayPathSphere ray(vm1d,a.u.mag(),zmax,tmax,dz,string("z"));
    a.t1d=ray.t[ray.t.size()-1];
    dmatrix path=RayPathProject(a.lat,a.lon,0.0,vm3d,ray,a.u);
    vector<double> t3dv;
    try {
        t3dv=pathintegral(vm3d,path);
    }catch(GCLgrid_error& gerr)
    {
        throw gerr;
    }
    int nt3dv=t3dv.size();
    if(nt3dv!=path.columns())
        throw SeisppError(string("compute_residuals procedure:  ")
            + "pathintegral path was truncated.  Cannot computed 3d time");
    a.t3d=t3dv[nt3dv-1];
    a.resid3d=a.timeres-a.t3d;
}
void remove_mean_x3(GCLscalarfield3d& f)
{
	double mean;
	int i,j,k;
	double nval=static_cast<double>((f.n1)*(f.n2));
	cout << "Mean removal values"<<endl
		<< "grid-index     mean"<<endl;
	for(k=0;k<f.n3;++k)
	{
		mean=0.0;
		for(i=0;i<f.n1;++i)
			for(j=0;j<f.n2;++j)
				mean += f.val[i][j][k];
		mean /= nval;
		cout << k << "      "<<mean<<" "<<1.0/mean<<endl;
		for(i=0;i<f.n1;++i)
			for(j=0;j<f.n2;++j)
				f.val[i][j][k] = f.val[i][j][k]-mean;
	}
}
void VelocityFieldToSlowness(GCLscalarfield3d& g)
{
   for(int i=0;i<g.n1;++i)
      for(int j=0;j<g.n2;++j)
       for(int k=0;k<g.n3;++k) g.val[i][j][k]=1.0/g.val[i][j][k];
}
double mean_removal(list<ArrivalRecord>& arlist)
{
    list<ArrivalRecord>::iterator arptr;
    double sumres;
    for(sumres=0.0,arptr=arlist.begin();arptr!=arlist.end();++arptr)
        sumres += arptr->timeres;
    double mean=sumres/static_cast<double>(arlist.size());
    for(arptr=arlist.begin();arptr!=arlist.end();++arptr)
    {
        arptr->timeres -= mean;
        arptr->resid3d -= mean;
        arptr->rawresid -= mean;
    }
    return(mean);
}
bool SEISPP::SEISPP_verbose(false);
void usage()
{
    cerr << "computeresiduals db outfile"<<endl;
    exit(-1);
}
int main(int argc, char **argv)
{
	if(argc!=3) usage();
	string dbname(argv[1]);
        string outfile(argv[2]);
        cout << "Phase residual computation program"<<endl
            << "Writing output table to file="<<outfile<<endl;
        ofstream ostrm;
        ostrm.open(argv[2],ios::out);
	cout << "Enter phase name in arrival table to process:";
	string phase;
	cin >> phase;
	char ss_phase[128];
	sprintf(ss_phase,"arrival.iphase=~/%s/",phase.c_str());
	try {
		DatascopeHandle dbh(dbname,false);
                cout << "Enter velocity model name to load from database:";
                string model1dname;
                cin >> model1dname;
                // This is not at all general
                string property=phase;
                VelocityModel_1d vmod1d(dbh.db,model1dname,property);
                cout << "Enter gridname and model (fieldname) for 3D model: ";
                string gridname,fieldname;
                cin >> gridname;
                cin >> fieldname;
                cout << "Trying to load fieldname="<<fieldname
                    << " associated with grid name="<<gridname<<endl;
                GCLscalarfield3d vmod3d(dbh,gridname,fieldname);
                cout << "Converting velocity to slowness.  Assuming units are km/s"<<endl;
                VelocityFieldToSlowness(vmod3d);
                cout << "Is this absolute velocity or difference?  Type y if model is absolute: ";
                string ques;
                cin >> ques;
                if(ques=="y")
                {
                    cout << "Assuming model is absolute velocity. Removing mean slowness at each depth level."
                        <<endl;
                    remove_mean_x3(vmod3d);
                }
                bool remove_mean_residual;
                cout << "Remove mean residual frome arrivals for each event?  (enter y to do this):";
                cin >> ques;
                if(ques=="y")
                    remove_mean_residual=true;
                else
                    remove_mean_residual=false;
                /*
                cout << "slowness perturbation field data"<<endl
                    << vmod3d <<endl;
                    */
                /* Now build the database view for the arrival data*/
                cout << "Building database view for arrival data"<<endl;
		dbh.lookup("event");
		list<string> sortkeys;
		sortkeys.push_back("evid");
		dbh.sort(sortkeys);
		dbh.natural_join("origin");
		string ss_to_prefor("orid==prefor");
		dbh.subset(ss_to_prefor);
		dbh.natural_join("assoc");
		dbh.natural_join("arrival");
		dbh.subset(string(ss_phase));
                dbh.natural_join("site");
		list<string> group_keys;
		group_keys.push_back("evid");
		dbh.group(group_keys);
		dbh.rewind();
		int evid;
		cout << "Found "<<dbh.number_tuples()<<" events to process"
			<<endl <<endl;
                list<ArrivalRecord> arrecs;
                list<ArrivalRecord>::iterator arptr;
		int i,j;
		for(i=0;i<dbh.number_tuples();++i,++dbh)
		{
			DBBundle bundle=dbh.get_range();
			Dbptr dbp=bundle.parent;
			Dbptr dbassoc;
			evid=dbh.get_int("evid");
                        double hlat,hlon,hz,htime;
                        arrecs.clear();

			for(dbp.record=bundle.start_record,j=0;
				dbp.record<bundle.end_record;++dbp.record,++j)
			{
				int iret;
                                Hypocenter h;
                                if(j==0)
                                {
                                    iret=dbgetv(dbp,0,"origin.lat",&hlat,
                                        "origin.lon",&hlon,
                                        "origin.time",&htime,
                                        "origin.depth",&hz,NULL);
                                    h=Hypocenter(rad(hlat),rad(hlon),hz,htime,
                                            string("tttaup"),string("iasp91"));
                                }
                                ArrivalRecord ar;
                                ar.evid=evid;
                                char sta[10];
                                double atime;
                                iret=dbgetv(dbp,0,"sta",sta,
                                        "time",&(ar.time),
                                        "deltim",&(ar.deltim),
                                        "site.lat",&(ar.lat),
                                        "site.lon",&(ar.lon),
                                        "assoc.timeres",&atime,NULL);
				if(iret==dbINVALID)
				{
					cerr << "dbgetv error at row "
						<< dbp.record<<endl;
				}
				else
				{
                                    ar.timeres=atime;
                                    ar.sta=string(sta);
                                    ar.lat=rad(ar.lat);
                                    ar.lon=rad(ar.lon);
                                    try {
                                    compute_residuals(ar,vmod1d,vmod3d,h,phase);
                                    } catch (SeisppError& serr)
                                    {
                                        cerr << "Problems with sta="<<sta<<" for evid="<<evid<<endl;
                                        serr.log_error();
                                        continue;
                                    }
                                    arrecs.push_back(ar);
				}
			}
                        if(remove_mean_residual) 
                        {
                            double mean = mean_removal(arrecs);
                            cout << "Evid="<<evid<<" removed mean value of="<<mean
                                << " from "<<arrecs.size()<<" arrivals for this event"<<endl;
                        }
                        for(arptr=arrecs.begin();arptr!=arrecs.end();++arptr)
                        {
                            ostrm << *arptr;
                        }
                        
		}
	}
	catch (SeisppError& serr)
	{
		serr.log_error();
		exit(-1);
	}
        catch (int ierr)
        {
            cerr << "Something threw int exception ="<<ierr<<endl
                << "Most likely this came from a GCLgrid library module"<<endl;
            exit(-1);
        }
}
				

