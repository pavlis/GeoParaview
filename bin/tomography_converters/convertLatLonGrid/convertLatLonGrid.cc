#include "VelocityModel_1d.h"
#include "gclgrid.h"
using namespace std;
using namespace SEISPP;
/* this is a program to convert any global tomography
   models to a GCLgrid file.  It assumes a regular spaced grid of lat from -90 to 90 
   and lon from -179 to 180 with 1 degree spacing. The depth interval should be an 
   argument defined. This was produced by modification of convertMIT15.    
   */

void usage()
{
    cout << "convertLatLonGrid db model_name P|S depth_min depth_max depth_interval "<<endl
        <<endl;
    exit(-1);
}
bool SEISPP::SEISPP_verbose(true);
int main(int argc, char **argv)
{
    try {
    int i,j,k,kk;
    if(argc<4) usage();
    string dbname(argv[1]);
    string model_name(argv[2]);
    bool convert_to_slowness(false);
    cout << "Will write results to Antelope db = "<< dbname<<endl;
    DatascopeHandle dbh(dbname,false);
    /*for(i=2;i<argc;++i)
    {
        string sarg(argv[i]);
        if(sarg=="-slow")
            convert_to_slowness=true;
        else
            usage();
    }*/
    bool isPm;
    string velocity_com;
    if(argv[3]==string("P"))
    {
    	isPm=true;
    	velocity_com="P";
    }
    else
    {
       	isPm=false;
    	velocity_com="S";
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
    double depmin(atof(argv[4])); 
    double depmax(atof(argv[5]));
    double depint(atof(argv[6]));
    cout << "depmin = "<< depmin << "depmax = "<< depmax << "depint = "<< depint <<endl;
    
    int n1vp(361),n2vp(181),n3vp((depmax-depmin)/depint+1);

    // 1 is along lon lines, 2 is along lat lines and 3 is vertical
    // oriented from bottom up to be close to gclgrid convention
    double dx1nom(27.5*4),dx2nom(24.0*4),dx3nom(depint);
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
		for(i=1;i<n1vp;++i)
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
		        cout << i << " " << j << " " << lon << " " << lat <<endl;
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
	            if(i==n1vp-1)
	            {
			        cp=dVP.gtoc(lat,-lon,r);
			        dVP.x1[0][j][k]=cp.x1;
			        dVP.x2[0][j][k]=cp.x2;
			        dVP.x3[0][j][k]=cp.x3;
	            	dVP.val[0][j][k]=dVP.val[i][j][k];
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
        dVP.save(dbh,griddir,modeldir,string("dV")+velocity_com,model_name+string("dV")+velocity_com);
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
