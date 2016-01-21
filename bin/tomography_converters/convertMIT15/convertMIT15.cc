#include "VelocityModel_1d.h"
#include "gclgrid.h"
using namespace std;
using namespace SEISPP;
/* this is a special program to convert Scott Burdick's downloaded P tomography
   model to a vtk file.  It has lots of frozen constants built in to the 
   program.  This was produced by modification of code originally written to 
   import the 2010 and 2011 models.   The later models distribute a matlab
   code to interpolate the irregular grid actually used in the mit tomography
   code.   This program was tested against output of the matlab m file
   reformatMIT2015.m. 
   */

void usage()
{
    cout << "convert2gcl db [-slow]"<<endl
        << "-slow outputs delta slowness - default is raw percent velocity change"<<endl
        << "  Note slowness version requires file ak135.mod1d"
        <<endl;
    exit(-1);
}
bool SEISPP::SEISPP_verbose(true);
int main(int argc, char **argv)
{
    try {
    int i,j,k,kk;
    if(argc<2) usage();
    string dbname(argv[1]);
    bool convert_to_slowness(false);
    cout << "Will write results to Antelope db = "<< dbname<<endl;
    DatascopeHandle dbh(dbname,false);
    for(i=2;i<argc;++i)
    {
        string sarg(argv[i]);
        if(sarg=="-slow")
            convert_to_slowness=true;
        else
            usage();
    }
    /* First we load the ak135 reference model from mod1d file.
    We do this always - not elegant but makes the code simpler */
    VelocityModel_1d  vmod(string("ak135.mod1d"),string("mod1d"),
            string("P"));
    double lat0(10.0);
    double lon0(-150.0);
    lat0=rad(lat0);
    lon0=rad(lon0);
    double r0;  // Computed from lat0 using r0_ellipse.  Earlier was wrong
    /* The 2015 model reformat runs 10 to 70 and -150 to -50 at intervals of
       0.25 degrees.   4 per degree means lat is 60*4+1 and lon is
       100*4+1 */
    const int n1vp(401),n2vp(241),n3vp(22);

    //const int n1vp(121),n2vp(201),n3vp(21);
    // 1 is along lon lines, 2 is along lat lines and 3 is vertical
    // oriented from bottom up to be close to gclgrid convention
    const double dx1nom(27.5),dx2nom(24.0),dx3nom(50.0);
    const int i0(0),j0(0);
    /* Scott tabulates the model scans with one line per horizontal position
       with intervals of 50 km ranging to 1000 km.  I will assume first scan
       is zero depth so base is actually 950.  May be off by 25 km and require
       a change.   After that format does scans on longitude lines with 
       longitude changing last.  Here n1 will longitude and n2 latitude
       so have to rearrange order in input loop below. sizes came from
       running nl on ascii file.  Note scan order is: lon,lat,z where
       the last (z) is the fastest varying.
    
    I recreated the same order for the matlab code. */

    /* This grid will not be congruent with receiver function image grid.
       gclfield2vtk must have remapping defined for this to be registered right.
        Here I use simple origin at the first grid point in the input scan order. */
    r0=r0_ellipse(lat0); 
    GCLgrid3d *grid=new GCLgrid3d(n1vp,n2vp,n3vp,string("MIT15"),
            lat0,lon0,r0,0.0,
            dx1nom,dx2nom,dx3nom,i0,j0);

    /* This clones the grid but does not set the field variables*/
    GCLscalarfield3d dVP(*grid);
    delete grid;
    double lat,lon,r;
    Geographic_point gp;
    Cartesian_point cp;
    int n3toread=n3vp-1;  // we have to pad an extra point at surface 
    for(i=0;i<n1vp;++i)
    {
        for(j=0;j<n2vp;++j)
        {
            cin >> lon;
            cin >> lat;
            cout << i << " " << j << " " << lon << " " << lat <<endl;
            lat=rad(lat);
            lon=rad(lon);
            double rsurface=r0_ellipse(lat);
            double depth,v0;
            for(k=n3toread-1,kk=0;kk<n3toread;--k,++kk)
            {
                depth=dx3nom*kk;
                r=rsurface-depth;
                cp=dVP.gtoc(lat,lon,r);
                dVP.x1[i][j][k]=cp.x1;
                dVP.x2[i][j][k]=cp.x2;
                dVP.x3[i][j][k]=cp.x3;
                cin >> dVP.val[i][j][k];
                if(convert_to_slowness)
                {
                    v0=vmod.getv(depth);
                    dVP.val[i][j][k] = -dVP.val[i][j][k]/v0/(1+dVP.val[i][j][k]);
                }
            }
            // copy top value to top cell 
            dVP.val[i][j][n3vp-1]=dVP.val[i][j][n3vp-2];
            r=rsurface+dx3nom;
            cp=dVP.gtoc(lat,lon,r);
            dVP.x1[i][j][n3vp-1]=cp.x1;
            dVP.x2[i][j][n3vp-1]=cp.x2;
            dVP.x3[i][j][n3vp-1]=cp.x3;
        }
    }
    string griddir("grids"),modeldir("vmodels");
    if(convert_to_slowness)
    {
        dVP.save(dbh,griddir,modeldir,string("dUP"),string("MIT15dUP"));
    }
    else
    {
        /* save the dvp data */
        dVP.save(dbh,griddir,modeldir,string("dVP3"),string("MIT11dVP"));
        /* Now we convert all values to absolute velocities using the ak135
           velocity model sorted in the vmod object */
        double v1d,depth;
        for(j=0;j<n2vp;++j)
        {
            for(i=0;i<n1vp;++i)
            {
                for(k=0;k<n3vp;++k)
                {
                    depth=dVP.depth(i,j,k);
                    v1d=vmod.getv(depth);
                /* model is percent perturbation of velocity */
                dVP.val[i][j][k]=v1d+(v1d*dVP.val[i][j][k]);
                }
            }
        }
        dVP.save(dbh,string(""),modeldir,string("MIT11VP"),string("MIT11VP"));
    }
    } catch (int ierr)
    {
        cerr << "dbsave threw error code "<< ierr<<endl
            << "Probably no output"<<endl;
    }
    catch (SeisppError& serr)
    {
        serr.log_error();
    }
}
