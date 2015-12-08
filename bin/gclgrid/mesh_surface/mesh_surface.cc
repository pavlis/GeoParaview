#include <iostream>
#include <fstream>
#include <sstream>
#include "seispp.h"
#include "gclgrid.h"
#include "GeoTriMeshSurface.h"
#include "GeoSplineSurface.h"
/* Use externally to check for undefined depths.   This is an elevation at about 10 
   earth radii */
const double UndefinedDepth(-99999.9);
void usage()
{
    cerr << "meshsurface ingrid outf [-af fname -pf pffile]< pointfile" <<endl;
    cerr << " outf is GCLscalarfield with depth as value"<<endl
        <<  " use -af to build a parallel file of attributes on surface"
        <<endl;;
    exit(-1);
}
bool SEISPP::SEISPP_verbose(true);
int main(int argc, char **argv)
{
    if(argc<3) usage();
    string ingridname(argv[1]);
    string outfieldname(argv[2]);
    string pffile("mesh_surface");
    bool data_has_attribute(false);
    string attribute_data_file;
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
        else if(sarg=="-af")
        {
            ++i;
            if(i>argc) usage();
            attribute_data_file=string(argv[i]);
            data_has_attribute=true;
        }
        else
            usage();
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
        /* Control points are read from input.*/
        double latin,lonin,zin;
        vector<Geographic_point> gp;
        /* Similar vector of data for attribute data*/
        vector<Geographic_point> gp_attr;
        int ncps(0);
        char linebuf[512];
        double val;
        Geographic_point gpin;
        while(cin.getline(linebuf,512))
        {
            stringstream ss(linebuf);
            ss >> lonin; ss>>latin;  ss>>zin;
            gpin.lat=rad(latin);
            gpin.lon=rad(lonin);
            gpin.r=r0_ellipse(gpin.lat)-zin;
            gp.push_back(gpin);
            //DEBUG
            //cout << lonin<<" "<<latin<<" "<<zin<<endl;
            ++ncps;
        }
        cout << "Read "<<ncps<<" control points from stdin"<<endl;
        if(data_has_attribute)
        {
            ifstream afin(attribute_data_file.c_str(),ios::in);
            if(afin.fail())
            {
                cerr << "Open failed on attribute file defined by -af "
                    << "flag="<<attribute_data_file<<endl
                    << "Fatal error"<<endl;
                exit(-1);
            }
            ncps=0;
            while(afin.getline(linebuf,512))
            {
                stringstream ss(linebuf);
                ss >> lonin; ss>>latin;  ss>>zin;
                gpin.lat=rad(latin);
                gpin.lon=rad(lonin);
                /* Note zin in this case is not that at all but 
                   the attribute value.  Reusing the variable only*/
                gpin.r=r0_ellipse(gpin.lat)-zin;
                gp_attr.push_back(gpin);
                ++ncps;
            }
            cout << "Read "<<ncps<<" attribute points from file "
                <<attribute_data_file<<endl;
        }
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
                double mr=r0_ellipse(mlat);
                if(bptr->is_defined(mlat,mlon))
                {
                    f.val[i][j]=bptr->depth(mlat,mlon);
                    mr-=f.val[i][j];
                    ++nlive;
                    Cartesian_point cp=mesh.gtoc(mlat,mlon,mr);
                    f.x1[i][j]=cp.x1;
                    f.x2[i][j]=cp.x2;
                    f.x3[i][j]=cp.x3;
                    if(data_has_attribute)
                        f.val[i][j]=abptr->depth(mlat,mlon);
                    //DEBUG
                    cout << deg(mlon) << " "<<deg(f.lon(i,j))<<" "
                        <<deg(mlat)<<" "<<deg(f.lat(i,j))
                        <<" "<<f.depth(i,j)<<" "<<f.val[i][j]<<endl;
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
