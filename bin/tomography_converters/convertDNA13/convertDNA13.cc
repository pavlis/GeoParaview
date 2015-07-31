#include <iostream>
#include <fstream>
#include "VelocityModel_1d.h"
#include "gclgrid.h"
using namespace std;
using namespace SEISPP;
/* this is a special program to convert the Berkeley Dynamic North America model
   aka 2010 to a gcl grid that can be converted to a vtk file for inclusion in the TA scene.
   It has lots of frozen constants built in to the 
   program.  Not originally intended for release, but useful enough to be put
   into pwmig release.  
   */
 

void usage()
{
    cout << "convertucb -dVp|-dVSV|-dVSH|-dVSJ [-dU]< DNA13_data_file"<<endl
        <<"   Output is file name determined by switch.  "<<endl
        << "   -dU sets output to delta slowness (requires file ak135.mod1d)"
        <<endl << "    (default is tabulated percent velocity change)"
        <<endl;
    exit(-1);
}
bool SEISPP::SEISPP_verbose(true);
int main(int argc, char **argv)
{
    try {
    int i,j,k,l,kk;
    if(argc<2)usage();
    /* Default to dVp*/
    int dcol_to_use;  //0 is first data field after r,lat,lon
    string outfile;
    bool convert_to_slowness(false);
    bool model_is_S(true);
    /* Fancy, but need to look for dU flag first so we can 
       set names to dU when appropriate */
    for(i=1;i<argc;++i)
    {
        string sarg(argv[i]);
        if(sarg=="-dU") convert_to_slowness=true;
    }
    for(i=1;i<argc;++i)
    {
        string sarg(argv[i]);
        if(sarg=="-dVp") 
        {
            dcol_to_use=0;
            if(convert_to_slowness)
                outfile="DNA13_dUP";
            else
                outfile="DNA13_dVP";
            model_is_S=false;
        }
        else if(sarg=="-dVSV")
        {
            dcol_to_use=1;
            if(convert_to_slowness)
                outfile="DNA13_dUSV";
            else
                outfile="DNA13_dVSV";
        }
        else if(sarg=="-dVSH")
        {
            dcol_to_use=2;
            if(convert_to_slowness)
                outfile="DNA13_dUSH";
            else
                outfile="DNA13_dVSH";
        }
        else if(sarg=="-dVSJ")
        {
            dcol_to_use=3;
            if(convert_to_slowness)
                outfile="DNA13_dUSJ";
            else
                outfile="DNA13_dVSJ";
        }
        else if(sarg=="-dU")
            cout << "Converting data to slowness perturbations from ak135"<<endl;
        else
            usage();
    }
    double lat0(25.0);
    double lon0(-126.0);
    /* gclgrid requires these in radians */
    lat0=rad(lat0);
    lon0=rad(lon0);
    double r0=r0_ellipse(lat0);
    /* These are frozen constants for DNA10 */
    const int n1vp(109),n2vp(55),n3vp(129);
    const double dx1nom(25.0),dx2nom(25.0),dx3nom(25.0);
    const int i0(0),j0(0);
    VelocityModel_1d *V1d;
    if(convert_to_slowness)
    {
        try{
            if(model_is_S)
                V1d=new VelocityModel_1d(string("ak135.mod1d"),string("mod1d"),
                        string("S"));
            else
                V1d=new VelocityModel_1d(string("ak135.mod1d"),string("mod1d"),
                        string("P"));
        }catch(SeisppError& serr)
        {
            serr.log_error();
            cerr << "You must have the file ak135.mod1d in your working directory"<<endl
                << "This file can be produced with the program convert_ak135"
                <<endl;
            exit(-1);
        }
    }

    /* Build a template for the grid that will define this model */
    GCLgrid3d *grid=new GCLgrid3d(n1vp,n2vp,n3vp,string("UCBtomo13"),
            lat0,lon0,r0,0.0,
            dx1nom,dx2nom,dx3nom,i0,j0);

    /* This clones the grid but does not set the field variables*/
    GCLscalarfield3d dVS(*grid);
    delete grid;
    double lat,lon,r;
    double dep,dvsin;
    int nVcol(4);
    double dVinValues[4];
    double v0;   // 1d velocity value - used only for slowness conversion
    Geographic_point gp;
    Cartesian_point cp;
    char linebuffer[256];
    /* Skip the header line */
    cin.getline(linebuffer,256);
    /* This model scans in order fairly conguenty with GCLgrid geometry but
       a different scan order thatn natural.  Order is longitude fastest, 
       latitude second, and depth last.  All are arranged right handed, at least,
       so there are no inverted indices.  Most unusual is depth moves from bottom up,
       which is more rational in my view.*/
            for(i=0;i<n1vp;++i)
        for(j=0;j<n2vp;++j)
    for(k=0;k<n3vp;++k)
    /*  REMOVE ME.  Retained for initial hacking.
    for(k=n3vp-1,kk=0;kk<n3vp;--k,++kk)
      for(i=0;i<n1vp;++i)
        for(j=n2vp-1;j>=0;--j)
        */
        {
            cin >> r;
            cin >> lat;
            cin >> lon;
            for(l=0;l<nVcol;++l) cin >> dVinValues[l];
            dvsin=dVinValues[dcol_to_use];
            //cin >> dummy;  // needed because format has a column of zeros here 
            lat=rad(lat);
            lon=rad(lon);
            double rsurface=r0_ellipse(lat);
            /*This model uses a spherical approximation but everything
              we've done uses reference ellipsoid.  Have to convert or this 
              will look funny.*/
            dep=6371.0-r;
            r=rsurface-dep;
            cp=dVS.gtoc(lat,lon,r);
            dVS.x1[i][j][k]=cp.x1;
            dVS.x2[i][j][k]=cp.x2;
            dVS.x3[i][j][k]=cp.x3;
            if(convert_to_slowness)
            {
                v0=V1d->getv(dep);
                /* conversion is -1.0*dv/v^2.  Original in percent to 
                require a 1/100 scaling */
                dvsin=-0.01*dvsin/v0;
            }

            dVS.val[i][j][k]=dvsin;
        }
        //string griddir("grids"),modeldir("vmodels");
        /* save the dvs data. arg2 is null string to prevent save error for duplicate
        grid */
            /*
        dVS.save(dynamic_cast<DatabaseHandle&>(dbh),
                griddir,modeldir,string("dVP13"),string("dVPDNA13"));
                */
        // Now save to file 
        string modeldir("DNA13_vmodels");
        dVS.save(outfile,modeldir);
        /* Now we convert all values to absolute velocities using the ak135
           velocity model sorted in the vmod object */
//        double v1d,depth;
//        for(j=0;j<n2vp;++j)
//        {
//            for(i=0;i<n1vp;++i)
//            {
//                for(k=0;k<n3vp;++k)
//                {
//                    depth=dVP.depth(i,j,k);
//                    v1d=vmod.getv(depth);
//                /* model is percent perturbation of velocity */
//                dVP.val[i][j][k]=v1d+(0.01*v1d*dVP.val[i][j][k]);
//                }
//            }
//        }
//        dVP.dbsave(db,string(""),modeldir,string("VPsig08"),string("VPsig08"));
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
