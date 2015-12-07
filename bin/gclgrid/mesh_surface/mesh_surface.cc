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
        Metadata md(pf);
        /* f is our output field.  This clones the mesh geometry but 
           sets all val entries to 0.0*/
        GCLscalarfield f(mesh);
        string surface_type=md.get_string("surface_fitting_method");
        /* Major switch to control what is posted to field.  When this
           boolean is false (default) only 3 columns are read and the
           scalar is set to depth.   When true an alternative value 
           in column 4 becomes the field variable. */
        bool data_has_attribute=md.get_bool("data_has_attribute");
        /* Control points are read from input.*/
        double latin,lonin,zin;
        vector<Geographic_point> gp;
        /* These duplicate gp vector above, but radius is 
           replaced by radius less the attribute.  Depends 
           upon the spline function converting back, which 
           the current version does.  Very big maintenance
           issue. */
        vector<Geographic_point> gp_attr;
        int ncps(0);
        char linebuf[512];
        double val;
        while(cin.getline(linebuf,512))
        {
            stringstream ss(linebuf);
            ss >> lonin; ss>>latin;  ss>>zin;
            Geographic_point gpin;
            gpin.lat=rad(latin);
            gpin.lon=rad(lonin);
            gpin.r=r0_ellipse(gpin.lat)-zin;
            //DEBUG
            //cout << deg(gpin.lon) <<" "<<deg(gpin.lat)<<" "<< r0_ellipse(gpin.lat) - gpin.r<<endl;
            gp.push_back(gpin);
            if(data_has_attribute)
            {
                ss >> val;
                gpin.r=r0_ellipse(gpin.lat)-val;
                gp_attr.push_back(gpin);
            }
            ++ncps;
        }
        cout << "Read "<<ncps<<" control points from stdin"<<endl;
        //base pointer passed around for polymophic behaviour
        GeoSurface *bptr; 
        // Comparable when using an attribute added to surface
        GeoSurface *abptr;
        if(surface_type=="DelaunayTriangularization")
        {
            GeoTriMeshSurface *gtmsptr=new GeoTriMeshSurface(gp);
            bptr=dynamic_cast<GeoSurface*>(gtmsptr);
            if(data_has_attribute)
            {
                GeoTriMeshSurface *gatmpptr
                    =new GeoTriMeshSurface(gp_attr);
                abptr=dynamic_cast<GeoSurface*>(gatmpptr);
            }
        }
        else if(surface_type=="BicubicSpline")
        {
            GeoSplineSurface *gssptr=new GeoSplineSurface(gp,md);
            bptr=dynamic_cast<GeoSurface*>(gssptr);
            if(data_has_attribute)
            {
                GeoSplineSurface *gatmpptr
                    = new GeoSplineSurface(gp_attr,md);
                abptr=dynamic_cast<GeoSurface*>(gatmpptr);
            }
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
                    if(data_has_attribute)
                        f.val[i][j]=abptr->depth(mlat,mlon);
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
