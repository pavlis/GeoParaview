#include <stdio.h>
#include <fstream>
#include <vector>
#include "RegionalCoordinates.h"
#include "seispp.h"
#include "Metadata.h"
#include "LatLong-UTMconversion.h"
using namespace std;
using namespace SEISPP;
// Warning:  the vector index will be off by one (C instead of FORTRAN order)
typedef vector<Cartesian_point> PointSet;

PointSet ExtractPoints(FILE *fp)
{
    PointSet utmpoints;
    string tag;
    int node,lastnode;
    double x,y,z;
    char line[256];
    int linenumber(1);
    while(fgets(line,256,fp)!=NULL)
    {
        char chtag[20];
        Cartesian_point pt;
        sscanf(line,"%s",chtag);
        tag=string(chtag);
        /*Silently skip any line without VRTX as the first token */
        if(tag=="VRTX")
        {
            sscanf(line,"%s%d%lf%lf%lf",chtag,&node,&x,&y,&z);
            if(linenumber==1)
                lastnode=node-1;
            else if(node!=(lastnode+1))
            {
                cerr << "Data error in input file"<<endl
                    << "Require points (VRTX tag) values start at 1 and "
                    << "increment by 1"<<endl
                    << "Found node number "<<lastnode<<" followed by "<<node
                    <<endl
                    << "You will need to manually repair this file"
                    <<endl;
                exit(-3);
            }
            //DEBUG
            //cout << node<<" "<<x<<" "<<y<<" "<<z<<endl;
            ++linenumber;
            lastnode=node;
            pt.x1=x;
            pt.x2=y;
            pt.x3=z;
            utmpoints.push_back(pt);
        }
    }
    return(utmpoints);
}
typedef vector<vector<int> > Triangles;

Triangles ExtractTriangles(FILE *fp)
{
    Triangles tri;
    string tag;
    int n1,n2,n3;
    rewind(fp);
    char line[256];
    while(fgets(line,256,fp)!=NULL)
    {
        char chtag[20];
        sscanf(line,"%s",chtag);
        tag=string(chtag);
        /*Silently skip any line without TRGL as the first token */
        if(tag=="TRGL")
        {
            sscanf(line,"%s%d%d%d",chtag,&n1,&n2,&n3);
            vector<int> current_tri;
            current_tri.push_back(n1);
            current_tri.push_back(n2);
            current_tri.push_back(n3);
            tri.push_back(current_tri);
        }
    }
    return(tri);
}
void usage()
{
    cerr << "gocad2vtk infile outfile [-pf pffile]"<<endl;
    cerr << "  converts gocad triangle surface data to vtk file"<<endl;
    exit(-1);
}


bool SEISPP::SEISPP_verbose(true);

int main(int argc, char **argv)
{
    ios::sync_with_stdio();
    if(argc<3) usage();
    string infile(argv[1]);
    string outfile(argv[2]);
    string pfname("gocad2vtk");
    int i,j;
    for(i=3;i<argc;++i)
    {
        string sarg(argv[i]);
        if(sarg=="-pf")
        {
            ++i;
            if(i>=argc) usage();
            pfname=string(argv[i]);
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
    FILE *fpin;
    fpin=fopen(infile.c_str(),"r");
    if(fpin==NULL)
    {
        cerr << "Cannot open input gocad file="<<infile<<endl;
        usage();
    }
    rewind(fpin);   // not really necessary
    try{
        ofstream outstrm(outfile.c_str(),ios::out);
        Metadata control(pf);
        /* These are parameters used to define regional coordinate system
           used for paraview - descendent of gclgrid library. 
           The names are consistent with vtk_gcl_converter */
        double lat,lon,radius,az;
        lat=control.get_double("origin_latitude");
        lon=control.get_double("origin_longitude");
        radius=control.get_double("origin_radius");
        az=control.get_double("azimuth_y_axis"); 
        RegionalCoordinates vtkcoords(rad(lat),rad(lon),radius,rad(az));
        /* UTM parameters */
        string UTMzone=control.get_string("UTM_zone");
        /*Freeze as WGS-84 */
        int RefEllipsoid(23);
        /* Echo all nonessential data to output for quality control */
        cout << "gocad2vtk conversion program"<<endl
            << "This program assumes coordinates are UTM"
            << "The following are all header lines and other lines "
            << "that this program will ignore:"<<endl;
        char line[256];
        char tag[20];
        long foff(0);
        while(fgets(line,256,fpin)!=NULL)
        {
            sscanf(line,"%s",tag);
            // Save foff as we need to seek back when the 
            // vertex keyword is found
            long current=ftell(fpin);
            string stag(tag);
            /* Assume VRTX data start the data section */
            if(stag=="VRTX") break;
            cout << line;
            foff=current;
        }
        string ques;
        cout << "///////////End of Header//////////"<<endl
            << "Ok to continue? (y or n): ";
        cin >> ques;
        if(ques!="y") 
        {
            cout << "Exiting"<<endl;
            exit(-2);
        }
        fseek(fpin,foff,SEEK_SET);
        PointSet utmpoints=ExtractPoints(fpin);
        if(utmpoints.size()<3)
        {
            cerr <<"gocad2vtk problem parsing points from file "<<infile<<endl
                << "Number of points found in file ="<<utmpoints.size()<<endl;
            exit(-1);
        }
        Triangles t=ExtractTriangles(fpin);
        if(t.size()<1)
        {
            cerr << "gocad2vtk problem parsing triangle definitions from file "
                <<infile<<endl
                <<"Found no triangles defined with key TRGL"<<endl;
        }
        int npoints,ntriangles;
        npoints=utmpoints.size();
        ntriangles=t.size();
        cout << "Finished parsing files"<<endl
            << "Number of points read ="<<npoints<<endl
            << "Number of triangle definitions read = "<<ntriangles<<endl;
        /* Next step is to write header material to vtk file */
        outstrm << "# vtk DataFile Version 2.0"<<endl
            << "gocad2vtk conversion"<<endl
            << "ASCII"<<endl
            << "DATASET POLYDATA"<<endl;
        /* Write point data.  This requires a double coordinate conversion
           utm2to dd and then dd to vtk coordinates */
        outstrm << "POINTS " << npoints<<" float"<<endl;
        PointSet::iterator psitr;
        for(psitr=utmpoints.begin();psitr!=utmpoints.end();++psitr)
        {
            double lat,lon;
            UTMtoLL(RefEllipsoid,psitr->x2,psitr->x1,UTMzone.c_str(),lat,lon);
            //DEBUG
            /*
            cout << "Conversion (x,y,lon,lat):  "<<psitr->x1 <<" "<<psitr->x2
                        <<" "<<lon<<" " <<lat<<endl;
                        */
            // convert these to radians 
            lat=rad(lat);
            lon=rad(lon);
            double r0=r0_ellipse(lat);
            /* Assume z values are units of m and positive up.  For
               paraview we need to convert this to radius in km */
            double r=r0;
            r+=((psitr->x3)/1000.0);
            Cartesian_point vtkpoint=vtkcoords.cartesian(lat,lon,r);
            outstrm << vtkpoint.x1 <<" "
                << vtkpoint.x2 << " "
                << vtkpoint.x3 << endl;
        }
        // second number is number of items to read. 
        outstrm << "POLYGONS "<<ntriangles<<" "<<ntriangles*4<<endl;
        for(i=0;i<ntriangles;++i)
        {
            outstrm << "3";
            /* GOCAD index is fortran like starting at 1.  vtk format 
               uses C order starting at 0.  Hence, the -1 factor in nodes 
               number output */
            for(j=0;j<3;++j) outstrm << " "<< t[i][j]-1;
            outstrm << endl;
        }
    }catch(std::exception& stdexcept)
    {
        cerr << "Something threw a low level exception"<<endl
            << "Message posted:  "
            << stdexcept.what()<<endl;
        exit(-1);
    }
    catch(SeisppError& serr)
    {
        serr.log_error();
        exit(-1);
    }
}
