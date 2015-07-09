#include <iostream>
#include <fstream>
#include <string>
#include "coords.h"
#include "Metadata.h"
#include "seispp.h"
#include "interpolator1d.h"
/* this is used several times here.  Beware if you ever cut and paste from this file */
const double REARTH(6378.17);
using namespace std;
using namespace SEISPP;
using namespace INTERPOLATOR1D;
const string prog("slabproject");
typedef list<Geographic_point> Profile;
typedef vector<Profile> ProfileEnsemble ;
/*! Object to project a small circle path defining a plate boundary.

  In plate tectonics transform plate boundaries and motion vectors
are defined by small circles (latitude like lines) relative to a pole
of rotation.  This object simplifies the task of creating these
lines. */

class PlateBoundaryProjector
{
public:
    PlateBoundaryProjector(double pla, double plo, double ola, double olo);
    /* Return coordinates as a GCLgrid Geographic_point object at distance 
       s (in km) from origin along the small circle */
    Geographic_point position(double s);
private:
    double plat,plon,olat,olon;
    double ssdelta;  // small circle distance from pole to origin (radians)
    double azimuth0;  // azimuth at pole to point define as origin for path
    /* Small circle effective earth radius - stored for convenience */
    double Reffective;  
};
PlateBoundaryProjector::PlateBoundaryProjector(double pla, double plo,
        double ola, double olo)
{
    plat=pla;
    plon=plo;
    olat=ola;
    olon=olo;
    dist(plat,plon,olat,olon,&ssdelta,&azimuth0);
    Reffective=REARTH*sin(ssdelta);
}
Geographic_point PlateBoundaryProjector::position(double s)
{
    double daz=s/Reffective;
    daz += azimuth0;
    Geographic_point result;
    latlon(plat,plon,ssdelta,daz,&(result.lat),&(result.lon));
    result.r=r0_ellipse(result.lon);
    return result;
}
/* pair of simple helper functions to convert between distances in radians
   and distances in km */
double delta2km(double drad)
{
    double delkm=drad*REARTH;
    return(delkm);
}
double km2delta(double dkm)
{
    double drad=dkm/REARTH;
    return(drad);
}

/* This is the file read routine for this program.  Inputs args are:
   fname - file name to be read
   plat, plon - plate pole latitude and longitude (radians)
   dsi - distance sample interval (km).  Important sign convention.
        Will throw an exception if sign of dsi and distance values
        in input differ 
   reverse_profiles - if true the profile data will have order reversed
   on return.  This is necessary sometimes to keep coordinate system
   of slab grid right handled.  
*/
ProfileEnsemble LoadProfiles(string fname,double plat, double plon,
        double dsi, bool reverse_profiles)
{
    FILE *fp;
    const string base_error("LoadProfiles procedure:  ");
    if(fname=="stdin")
    	fp=stdin;
    else
    {
        fp=fopen(fname.c_str(),"r");
        if(fp==NULL)
            throw SeisppError(base_error + "Open file on file="+fname);
    }
    ProfileEnsemble result;
    Profile member;
    int npoints;
    double olat,olon;
    double dist,depth,odepth;
    int pfnumber(1);
    while(fscanf(fp,"%d %lf %lf %lf",&npoints,&olat,&olon,&odepth)!=EOF)
    {
	cout << "Profile number "<<pfnumber<<" distance, depth pairs (km)"<<endl;
	cout << "Origin latitude="<<olat<<" Origin longitude="<<olon<<endl;
	++pfnumber;
        olat=rad(olat);  olon=rad(olon);
        PlateBoundaryProjector projector(plat,plon,olat,olon);
        member.clear();
        /* simple way to push the origin as the first point in profile*/
        Geographic_point nextpoint=projector.position(0.0);;
        nextpoint.r -= odepth;
        member.push_back(nextpoint);
	vector<double> dist_irregular, depth_irregular;
	dist_irregular.reserve(npoints);
	depth_irregular.reserve(npoints);
        for(int i=0;i<npoints;++i)
        {
            if(fscanf(fp,"%lf%lf",&dist,&depth)==EOF) throw SeisppError(base_error
                    + "Read error.  EOF encountered prematurely.\n"
                    +string("Probably formatting error in file=")+fname);
//DEBUG
cout << "dist,depth read="<<dist<<" "<<depth<<endl;
            if(dist*dsi<0.0) throw SeisppError(base_error
                    + "Sign inconsistency of profile sample interval and profile distance data\n"
                    + string("Both must be of same sign"));
	    dist_irregular.push_back(dist);
	    depth_irregular.push_back(depth);
	}
	double *depth_regular;
	/* Found needed to oversample in distance/depth to allow large distance output 
	sampling to not have irregular depths not matching target points */
	double distance_multiplier(50.0);
	double ddfine=dsi/50.0;
	int nfine=static_cast<int>(dist_irregular[npoints-1]/ddfine);
	depth_regular=new double[nfine];
	linear_scalar_irregular_to_regular(npoints,&(dist_irregular[0]),&(depth_irregular[0]),
		nfine,0.0,ddfine,depth_regular);
	/* depth_regular now contains depths resampled at fixed distances of ddfine.  Now build
	grid lines at fixed radial distance intervals honoring the original irregular points
	within a tolerance less than ddfine */
        double lastdepth=odepth;
        double lastdist=0.0;
        double currentdist=0.0;
        double currentdepth=odepth;
	double anchordepth,anchordist;  //last point read, not computed
	anchordepth=odepth;  anchordist=0.0;
        for(int i=0;i<npoints;++i)
        {
	    dist=dist_irregular[i];
	    depth=depth_irregular[i];
	    double theta=atan2((depth-anchordepth),(dist-anchordist));
	    anchordepth=depth;
	    anchordist=dist;
	    double ddz=dsi*sin(theta);
	    double ddx=dsi*cos(theta);
            do {
                currentdist+=ddx;
                double ddist=fabs(currentdist-lastdist);
                currentdepth+=ddz;
                nextpoint=projector.position(currentdist);
                nextpoint.r -= currentdepth;
		//cout << currentdist <<"   "<<currentdepth<<endl;
                member.push_back(nextpoint);
            }while(fabs(currentdist)<fabs(dist));
            lastdist=currentdist;
            lastdepth=currentdepth;
        }
        result.push_back(member);
    }
    if(reverse_profiles)
    {
        ProfileEnsemble::reverse_iterator pitr;
        ProfileEnsemble resrev;
        resrev.reserve(result.size());
        for(pitr=result.rbegin();pitr!=result.rend();++pitr)
        {
            resrev.push_back(*pitr);
        }
        return resrev;
    }
    else
        return result;
}
int MinProfileLength(ProfileEnsemble& profs)
{
    int result=0;
    vector<Profile>::iterator pitr;
    pitr=profs.begin();
    result=pitr->size();
    if(profs.size()==1) return(result);
    ++pitr;
    for(;pitr!=profs.end();++pitr)
    {
        int thissize=pitr->size();
        if(thissize<result)result=thissize;
    }
    return(result);
}

/* Returns a 3xsize(p) matrix of Cartesian points in the coordinates
   of g computed from the geographic points stored in p.  
   */
dmatrix ConvertProfile(Profile& p, BasicGCLgrid& g)
{
    Geographic_point gp;
    Cartesian_point cp;
    int n2=p.size();
    dmatrix result(3,n2);
    result.zero();
    Profile::iterator pitr;
    int j;
    for(j=0,pitr=p.begin();pitr!=p.end();++pitr,++j)
    {
        gp=(*pitr);
        cp=g.gtoc(gp);
        result(0,j)=cp.x1;
        result(1,j)=cp.x2;
        result(2,j)=cp.x3;
    }
    return(result);
}
void copy_gridline(GCLgrid& g, dmatrix& p, int i)
{
    if(g.n2>p.columns())
        throw SeisppError(string("copy_gridline:  coding error")
                +"grid size and profile size are inconsistent");
    int j;
    for(j=0;j<g.n2;++j)
    {
        g.x1[i][j]=p(0,j);
        g.x2[i][j]=p(1,j);
        g.x3[i][j]=p(2,j);
    }
}
/* Helper for RemoveStrain.  Returns Geographic_point for a point
   in a path */
Geographic_point ConvertPoint(GCLgrid& g, dmatrix& p, int i)
{
    Cartesian_point cp;
    cp.x1=p(0,i);
    cp.x2=p(1,i);
    cp.x3=p(2,i);
    return(g.ctog(cp));
}
/* overloaded version for doing inverse of above */
double *ConvertPoint(GCLgrid& g, Geographic_point gp)
{
    Cartesian_point cp;
    cp=g.gtoc(gp);
    double *result=new double[3];
    result[0]=cp.x1;
    result[1]=cp.x2;
    result[2]=cp.x3;
    return(result);
}
/* Helper for RemoveStrain. Strips row 3(4) used to hold
computed longitudinal strain */
dmatrix copy_profile(dmatrix& old)
{
	int np=old.columns();
	if(old.rows()!=4) throw SeisppError(
		string("copy_profile:  coding error.  old nrows not 4"));
	dmatrix result(3,np);
	int i,j;
	for(j=0;j<np;++j)
		for(i=0;i<3;++i) result(i,j)=old(i,j);
	return result;
}

/* grid is GCLgrid used for coordiante conversion, 
anchor is reference path and path is the path to be deformed.
maxdip is the maximum allowed dip (in radians).

Return has coordinates in 0,1,2 and longitudinal strain in 3. */
dmatrix RemoveStrain(GCLgrid& grid, dmatrix& anchor, dmatrix& path,
	double maxdip)
{
//DEBUG
//cout << "anchor path"<<endl<<anchor<<endl;
//cout << "path to deform"<<endl<<path<<endl;
   int np=grid.n2;
   if((path.columns()<np) || (anchor.columns()<np)) 
       throw SeisppError(string("RemoveStrain:  ")
           + "size mismatch on pair of paths to process.  Coding error");
   dmatrix result(4,np);
   double tanmaxdip(tan(fabs(maxdip)));
   double v[3];
   dcopy(3,path.get_address(0,0),1,v,1);
   daxpy(3,-1.0,anchor.get_address(0,0),1,v,1);
   /* v now contains vector between first point in each path.
      Preserve that length */
   double r=dnrm2(3,v,1);
cout << "r="<<r<<endl;
   /* These hold geographic point equivalents of anchor and path points
      respectively */
   Geographic_point gpa, gpp;
   /* copy the first point directly with no change */
   int j;
   for(j=0;j<3;++j) result(j,0)=path(j,0);
   /* main loop */
   for(j=1;j<np;++j)
   {
       double delta,azimuth;
       gpa=ConvertPoint(grid,anchor,j);
       gpp=ConvertPoint(grid,path,j);
       dist(gpa.lat,gpa.lon,gpp.lat,gpp.lon,&delta,&azimuth);
       double distkm=delta2km(delta);
       double az=r0_ellipse(gpa.lat)-gpa.r;
       double pz=r0_ellipse(gpp.lat)-gpp.r;
       double dz=fabs(pz-az);
       double newdist,lstrain;
       double dip=atan2(dz,distkm);
cout << "Computed dip="<<deg(dip)<<endl;
       if(fabs(dip)>maxdip)
       {
       	    newdist=dz/tanmaxdip;
	    lstrain=(sqrt(newdist*newdist+dz*dz)-r)/r;
cout << "In high dip block. strain="<<lstrain<<endl;
	    result(3,j)=lstrain;
	}
	else
	{
	    newdist=sqrt(r*r-dz*dz);
	    result(3,j)=0.0;
	}
       double newdelta=km2delta(newdist);
       latlon(gpa.lat,gpa.lon,newdelta,azimuth,&(gpp.lat),&(gpp.lon));
       double *vptr;
       vptr=ConvertPoint(grid,gpp);
       dcopy(3,vptr,1,result.get_address(0,j),1);
       delete [] vptr;
   }
//DEBUG
cout << "deformed path="<<endl<<result<<endl;
   return result;
}

    
/* This procedure actually builds the grid from profile data.
   It works outward from the mastergridline profile preserving line
   lengths between adjacentprofiles. 

   arguments 
   grid - GCLgrid modified in place to build slab geometry
   master - integer defining which profile is master
   profiles - profile data interpolated to uniform grid in 
    downdip paths.
*/
void BuildNoStrainGrid(GCLscalarfield& grid, int master, ProfileEnsemble& profiles,
	double maxdip)
{
  try {
    /* get local copies of these for readability*/
    int n1,n2;
    n1=grid.n1;
    n2=grid.n2;
    /* get reference to grid used below.  Relic of earlier version
    that used only a grid instead o a field holding longitudinal strain */
    /* It is convenient to store cartesian coordinates of a profile in
       these 3xn2 matrices */
    dmatrix masterprofile,lastprofile,rawprofile;
    /* first construct the master grid profile line */
    masterprofile=ConvertProfile(profiles[master],grid);
    copy_gridline(grid,masterprofile,master);
    /* Start distorting profiles above (higher j index) master */
    lastprofile=masterprofile;
    int j;
    for(j=master+1;j<n1;++j)
    {
        rawprofile=ConvertProfile(profiles[j],grid);
        dmatrix newprofile=RemoveStrain(grid,lastprofile,rawprofile,maxdip);
        copy_gridline(grid,newprofile,j);
	/* Need this to strip row 3(4) */
	lastprofile=copy_profile(newprofile);
    }
    /* Now do the same working down from master */
    for(j=master-1;j>=0;--j)
    {
        rawprofile=ConvertProfile(profiles[j],grid);
        dmatrix newprofile=RemoveStrain(grid,lastprofile,rawprofile,maxdip);
        copy_gridline(grid,newprofile,j);
	lastprofile=copy_profile(newprofile);
    }
  } catch (...){throw;};
}
/* This alternate procedure just builds a grid from profile data 
   without any changes.  That is, the profiles are simply interpolated
   to equally spaced points.  This is useful, for example, to build
   a 3d moho depth map 

   arguments 
   grid - GCLgrid modified in place to build desired grid
   profiles - profile data interpolated to uniform grid in 
    downdip paths.
*/
void BuildPlainGrid(GCLgrid& grid, ProfileEnsemble& profiles)
{
  try {
    int n1,n2;
    n1=grid.n1;
    n2=grid.n2;
    /* As above, store profiles in this 3xn2 matrix */
    dmatrix newprofile;
    int j;
    for(j=0;j<n1;++j)
    {
        newprofile=ConvertProfile(profiles[j],grid);
        copy_gridline(grid,newprofile,j);
    }
  } catch (...){throw;};
}

void usage()
{
	cerr << prog <<" db gridname [-pf pffile -V]"<<endl
		<< "outputs gridname to database db in gclgdisk table"<<endl;
	exit(-1);
}
bool SEISPP::SEISPP_verbose(false);
int main(int argc, char **argv)
{
	int i;
        ios::sync_with_stdio();
	string pfname(prog);
        if(argc<3) usage();
        string dbname(argv[1]);
        string gridname(argv[2]);
	for(i=3;i<argc;++i)
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
        Dbptr db;
        if(dbopen(const_cast<char *>(dbname.c_str()),"r+",&db))
        {
            cerr << "dbopen filed on database="<<dbname<<endl;
            usage();
        }
        db=dblookup(db,0,"gclgdisk",0,0);
        if(db.record==dbINVALID)
        {
            cerr << "dblookup failed for gclgdisk table for database="
                <<dbname<<endl;
            usage();
        }
	try {
	    Metadata control(pf);
            string savedir=control.get_string("output_grid_directory");
	    double plat=control.get_double("pole_latitude");
	    double plon=control.get_double("pole_longitude");
	    plat=rad(plat);
	    plon=rad(plon);
            string profile_fname=control.get_string("profile_data_file");
            double profsi=control.get_double("profile_sample_interval");
            bool reverse_profiles=control.get_bool("reverse_profile_data_order");
            bool remove_strain=control.get_bool("remove_transverse_strain");
            ProfileEnsemble RawProfiles=LoadProfiles(profile_fname,
                    plat,plon,profsi,reverse_profiles);
            int mastergridline=control.get_int("master_grid_line");
            if(mastergridline>=RawProfiles.size())
            {
                cerr << "master_grid_line parameter = "<<mastergridline
                    <<" is not consistent data file"<<endl
                    <<"Data file "<<profile_fname<<" has only "
                    <<RawProfiles.size()<<" profiles defined"<<endl;
                exit(-2);
            }
	    /* critical parameter.  When dip exceeds this longitudinal
	    strain will be allowed and computed.  Ignored unless using
	    no strain option */
	    double maxdip=control.get_double("maximum_allowed_dip");
	    maxdip=rad(maxdip);
	    string outfieldname=control.get_string("output_strain_fieldname");
            if(SEISPP_verbose)
            {
                cout << "Parameter file input"<<endl
                    << control<<endl;
                if(remove_strain) cout << "Profile order will be reversed"<<endl;
                if(remove_strain)
                    cout << "Applying no transverse strain algorithm"<<endl;
                else
                    cout << "Will output interpolated input with no distortion"
                        <<endl;
            }
            Profile master=RawProfiles[mastergridline];
            /* We set the origin of the grid created here to the anchor point
               (first point) for the master profile */
            Profile::iterator mgptr=master.begin();
            int n2grid=MinProfileLength(RawProfiles);
            /* We build the GCLgrid coordinates using the origin of 
            the mastergridline (0) as the coordinate origin.  We also
            set the nominal values of dx1 and dx2 as profsi. Not always
            right, but ok for this application.  The azimuth is also set
            to zero for convenience.  Can use the remap option of the
            vtk converter to make this work right in visualizations */
            GCLgrid grid(RawProfiles.size(),n2grid,gridname,
                    mgptr->lat,mgptr->lon,mgptr->r,0.0,profsi,profsi,
                    mastergridline,0);
            if(remove_strain)
	    {
		GCLscalarfield straingrid(grid);
                BuildNoStrainGrid(straingrid,mastergridline,RawProfiles,
			maxdip);
		straingrid.compute_extents();
		/* always put field and grid data in same dir */
		//straingrid.dbsave(db,savedir,savedir,outfieldname,outfieldname);
	    }
            else
	    {
                BuildPlainGrid(grid,RawProfiles);
                grid.compute_extents();
                grid.dbsave(db,savedir);
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
