#include "VelocityModel_1d.h"
#include "gclgrid.h"
using namespace std;
using namespace SEISPP;
/* this is a program to convert any global tomography
   models to a GCLgrid file.  It assumes a regular spaced grid of lat 
   and lon. The depth interval is not necessarily constant.
   This was produced from modification of convertMIT15.    
   */

void usage()
{
    cout <<"Usage:"<<endl
        <<"convertLatLonGrid outputname model_name P|S N_lat N_lon N_dep "
        <<"[-d -r -velocity -slow -lti lat_interval -lni lon_interval -di dep_interval] "<<endl<<endl
        <<"-d write to a database with name outputname. The default is write to GCLgrid file."<<endl
        <<"-r specify a regional grid. The default is global grid. "<<endl
        <<"-velocity assumes input to be actual velocity and convert it to velocity change in percent relative to AK135. "
            <<"The default does not do such convertion. "<<endl
        <<"-slow outputs delta slowness - default is raw velocity change in percent"<<endl
        <<"-lti -lni -di defines the nominal intervals of the grid in degree and km. "<<endl
        <<endl;
    exit(-1);
}
bool SEISPP::SEISPP_verbose(true);
int main(int argc, char **argv)
{
    try {
    int i,j,k,kk;
    if(argc<7) usage();
    string outputname(argv[1]);
    string model_name(argv[2]);
    //bool convert_to_slowness(false);
    bool write_to_db(false);
    bool is_regional(false);
    bool convert_to_dv(false);
    bool convert_to_slowness(false);
    double lat_interval(24.0*4);
    double lon_interval(27.5*4);
    double dep_interval(50.0);
    for(i=7;i<argc;++i)
    {
        string sarg(argv[i]);
        if(sarg=="-d")
            write_to_db=true;
        else if(sarg=="-r")
            is_regional=true;
        else if(sarg=="-velocity")
            convert_to_dv=true;
        else if(sarg=="-slow")
            convert_to_slowness=true;
        else if(sarg=="-lti")
        {
            ++i;
            if(i>=argc) usage();
            lat_interval*=atof(argv[i]);
        }
        else if(sarg=="-lni")
        {
            ++i;
            if(i>=argc) usage();
            lon_interval*=atof(argv[i]);
        }
        else if(sarg=="-di")
        {
            ++i;
            if(i>=argc) usage();
            dep_interval=atof(argv[i]);
        }
        else
            usage();
    }
    bool isPm;
    string velocity_com;
    if(argv[3]==string("P"))
    {
        isPm=true;
        velocity_com="P";
    }
    else if(argv[3]==string("S"))
    {
        isPm=false;
        velocity_com="S";
    }
    else
    {
        usage();
    }
    /* First we load the ak135 reference model from mod1d file.
    We do this always - not elegant but makes the code simpler */
    /*VelocityModel_1d  vmod(string("ak135.mod1d"),string("mod1d"),
            string("P"));*/
    double lat0(90.0);
    double lon0(0.0);
    lat0=rad(lat0);
    lon0=rad(lon0);
    double r0;  // Computed from lat0 using r0_ellipse.  
    /* The model grid runs -90 to 90 and -179 to 180 at intervals of
       1 degrees. (-180 to 180 is used to fix the ugly seam) means lat is 181 and lon is
       360 */
    int N_lat(atoi(argv[4])); 
    int N_lon(atoi(argv[5]));
    int N_dep(atoi(argv[6]));
    cout << "N_lat = "<< N_lat << ", N_lon = "<< N_lon << ", N_dep = "<< N_dep <<endl;
    VelocityModel_1d *V1d=NULL;
    if(convert_to_slowness || convert_to_dv)
    {
        try{
            if(isPm)
                V1d=new VelocityModel_1d(string("ak135.mod1d"),string("mod1d"),
                        string("P"));
            else
                V1d=new VelocityModel_1d(string("ak135.mod1d"),string("mod1d"),
                        string("S"));
        }catch(SeisppError& serr)
        {
            serr.log_error();
            cerr << "You must have the file ak135.mod1d in your working directory"<<endl
                << "This file can be produced with the program convert_ak135"
                <<endl;
            exit(-1);
        }
    }
    
    //int n1vp(361),n2vp(181),n3vp((depmax-depmin)/depint+1);
    int n1vp(N_lon),n2vp(N_lat),n3vp(N_dep);
    if(!is_regional)
        n1vp++;
    cout << "n2 = "<< n2vp << ", n1 = "<< n1vp << ", n3 = "<< n3vp <<endl;
    // 1 is along lon lines, 2 is along lat lines and 3 is vertical
    // oriented from bottom up to be close to gclgrid convention
    double dx1nom(lon_interval),dx2nom(lat_interval),dx3nom(dep_interval);
    const int i0(0),j0(0);
    /* Here n1 is longitude and n2 latitude
       so have to rearrange order in input loop below. 
       Note scan order is: z,lat,lon where
       the z is the fastest varying, then lon, then lat. */

    /* This grid will not be congruent with receiver function image grid.
       gclfield2vtk must have remapping defined for this to be registered right.
        Here I use simple origin at the first grid point in the input scan order. */
    r0=r0_ellipse(lat0);
    GCLgrid3d *grid=new GCLgrid3d(n1vp,n2vp,n3vp,model_name,
            lat0,lon0,r0,0.0,
            dx1nom,dx2nom,dx3nom,i0,j0);

    /* This clones the grid but does not set the field variables*/
    GCLscalarfield3d dVP(*grid);
    delete grid;
    double lat,lon,r;
    Geographic_point gp;
    Cartesian_point cp;
    int n3toread=n3vp;  // we have to pad an extra point at surface 
    for(j=0;j<n2vp;++j)
    {
        if(is_regional)
            i=0;
        else
            i=1;
        for(;i<n1vp;++i)
        {
            for(k=n3toread-1,kk=0;kk<n3toread;--k,++kk)
            {
                double depth,v0;
                cin >> depth;
                cin >> lat;
                cin >> lon;
                /*if(lon==180)
                {
                    cin >> depth;
                    cin >> depth;
                    cin >> lat;
                    cin >> lon;
                }*/
                //if(kk==0)
                //    cout << i << " " << j << " " << lon << " " << lat <<endl;
                lat=rad(lat);
                lon=rad(lon);
                double rsurface=r0_ellipse(lat);
                //depth=dx3nom*kk;
                //DEBUG
                //cout<<"Depth: "<<depth<<endl;
                r=rsurface-depth;
                cp=dVP.gtoc(lat,lon,r);
                dVP.x1[i][j][k]=cp.x1;
                dVP.x2[i][j][k]=cp.x2;
                dVP.x3[i][j][k]=cp.x3;
                cin >> dVP.val[i][j][k];
                if(convert_to_dv)
                {
                    v0=V1d->getv(depth);
                    if(dVP.val[i][j][k] == 9999)
                        dVP.val[i][j][k] = 0.0;
                    else
                        dVP.val[i][j][k] = (dVP.val[i][j][k]-v0)/v0*100.0;
                }
                if(convert_to_slowness)
                {
                    v0=V1d->getv(depth);
                    dVP.val[i][j][k] = -0.01*dVP.val[i][j][k]/v0/(1+0.01*dVP.val[i][j][k]);
                }
                if(!is_regional)
                {
                    if(i==n1vp-1)
                    {
                        cp=dVP.gtoc(lat,lon-rad(360.0),r);
                        dVP.x1[0][j][k]=cp.x1;
                        dVP.x2[0][j][k]=cp.x2;
                        dVP.x3[0][j][k]=cp.x3;
                        dVP.val[0][j][k]=dVP.val[i][j][k];
                    }
                }
                /*if(convert_to_slowness)
                {
                    v0=vmod.getv(depth);
                    dVP.val[i][j][k] = -dVP.val[i][j][k]/v0/(1+dVP.val[i][j][k]);
                }*/
                // copy top value to top cell
                /*if(k==n3vp-2)
                {
                    dVP.val[i][j][n3vp-1]=dVP.val[i][j][n3vp-2];
                    r=rsurface+dx3nom;
                    cp=dVP.gtoc(lat,lon,r);
                    dVP.x1[i][j][n3vp-1]=cp.x1;
                    dVP.x2[i][j][n3vp-1]=cp.x2;
                    dVP.x3[i][j][n3vp-1]=cp.x3;
                }*/
            }
        }
    }
    string griddir("grids"),modeldir("vmodels");
    /*if(convert_to_slowness)
    {
        dVP.save(dbh,griddir,modeldir,string("dUP"),string("MIT15dUP"));
    }*/
    /*else
    {*/
        /* save the dvp data */
    
    if(write_to_db)
    {
        cout << "Will write results to Antelope db = "<< outputname<<endl;
        DatascopeHandle dbh(outputname,false);
        if(convert_to_slowness)
            dVP.save(dbh,griddir,modeldir,string("dU")+velocity_com,model_name+string("dU")+velocity_com);
        else
            dVP.save(dbh,griddir,modeldir,string("dV")+velocity_com,model_name+string("dV")+velocity_com);
    }
    if(convert_to_slowness)
        dVP.save(outputname+string("_dU")+velocity_com,string("./"));
    else
        dVP.save(outputname+string("_dV")+velocity_com,string("./"));
        /* Now we convert all values to absolute velocities using the ak135
           velocity model sorted in the vmod object */
        //double v1d,depth;
        //for(j=0;j<n2vp;++j)
        //{
        //    for(i=0;i<n1vp;++i)
        //    {
        //        for(k=0;k<n3vp;++k)
        //        {
        //            depth=dVP.depth(i,j,k);
        //            v1d=vmod.getv(depth);
                /* model is percent perturbation of velocity */
        //        dVP.val[i][j][k]=v1d+(v1d*dVP.val[i][j][k]);
        //        }
        //    }
        //}
        //dVP.save(dbh,string(""),modeldir,string("MIT11VP"),string("MIT11VP"));
    //}
    if(V1d)
        delete V1d;
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
