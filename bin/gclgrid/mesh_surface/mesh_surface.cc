#include "seispp.h"
#include "gclgrid.h"
#include "GeoTriMeshSurface.h"
#include "GeoSplineSurface.h"
/* Use externally to check for undefined depths.   This is an elevation at about 10 
   earth radii */
const double UndefinedDepth(-99999.9);
void usage()
{
    cerr << "meshsurface ingrid outf [-pf pffile]< pointfile" <<endl;
    cerr << " outf is GCLscalarfield with depth as value"<<endl;
    exit(-1);
}
bool SEISPP::SEISPP_verbose(true);
int main(int argc, char **argv)
{
    if(argc<3) usage();
    string ingridname(argv[1]);
    string outfieldname(argv[2]);
    string pffile("mesh_surface");
    int i,j;
    for(i=3;i<argc;++i)
    {
        string sarg(argv[i]);
        if(sarg=="-pf")
        {
            ++i;
            if(i>argc) usage();
            pffile=string(argv[i]);
        }
    }
    Pf *pf;
    if(pfread(const_cast<char *>(pffile.c_str()),&pf))
    {
        cerr << "pfread file on pffile="<<pffile<<endl;
        usage();
    }
    try {
        cout << "Running meshsurface program."<<endl
            << "Loading GCLgrid "<<ingridname<<" to define surface meshing"
            <<endl;
        GCLgrid mesh(ingridname);
        /* f is our output field.  This clones the mesh geometry but 
           sets all val entries to 0.0*/
        GCLscalarfield f(mesh);
        /* Control points are read from input.*/
        double latin,lonin,zin;
        vector<Geographic_point> gp;
        int ncps(0);
        char linebuf[512];
        while(cin.getline(linebuf,512))
        {
            stringstream ss(linebuf);
            ss >> lonin; ss>>latin;  ss>>zin;
            Geographic_point gpin;
            gpin.lat=rad(latin);
            gpin.lon=rad(lonin);
            gpin.r=r0_ellipse(gpin.lat)-zin;
            //DEBUG
            cout << deg(gpin.lon) <<" "<<deg(gpin.lat)<<" "<< r0_ellipse(gpin.lat) - gpin.r<<endl;
            gp.push_back(gpin);
            ++ncps;
        }
        cout << "Read "<<ncps<<" control points from stdin"<<endl;
        Metadata md(pf);
        GeoSurface *bptr;  // base pointer passed around for polymophic behaviour
        string surface_type=md.get_string("surface_fitting_method");
        if(surface_type=="DelaunayTriangularization")
        {
            GeoTriMeshSurface *gtmsptr=new GeoTriMeshSurface(gp);
            bptr=dynamic_cast<GeoSurface*>(gtmsptr);
        }
        else if(surface_type=="BicubicSpline")
        {
            GeoSplineSurface *gssptr=new GeoSplineSurface(gp,md);
            bptr=dynamic_cast<GeoSurface*>(gssptr);
        }
        else
        {
            cerr << "Unrecognized value for parameter surface_fitting_method = "
                << surface_type<<endl
                << "Must be either DelaunayTriangularization OR BicubicSpline"
                << "Edit parameter file and try again"<<endl;
            exit(-1);
        }
        int nlive=0;
        for(i=0;i<mesh.n1;++i)
            for(j=0;j<mesh.n2;++j)
            {
                cout << "working on grid point i="<<i<<" j="<<j<<endl;
                double mlat(mesh.lat(i,j));
                double mlon(mesh.lon(i,j));
                double mr(mesh.r(i,j));
                if(bptr->is_defined(mlat,mlon))
                {
                    //DEBUG
                    //cout << deg(mlon) << " "<<deg(mlat)<<" "<<bptr->depth(mlat,mlon)<<endl;
                    f.val[i][j]=bptr->depth(mlat,mlon);
                    mr-=f.val[i][j];
                    ++nlive;
                    Cartesian_point cp=mesh.gtoc(mlat,mlon,mr);
                    f.x1[i][j]=cp.x1;
                    f.x2[i][j]=cp.x2;
                    f.x3[i][j]=cp.x3;
                }
                else
                    // Do not distort grid for undefined points - leaves them at original r
                    f.val[i][j]=UndefinedDepth;
            }
        f.save(outfieldname,string("."));
        cout << "Saved result to field file with root name="<<outfieldname<<endl
            << "There are "<<nlive<<" live points of total in grid of "
            <<f.n1*f.n2<<" points"<<endl;
    }catch(std::exception& stex)
    {
        cerr << stex.what();
        exit(-1);
    }
}


