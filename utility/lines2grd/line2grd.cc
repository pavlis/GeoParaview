#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include "coords.h"
/* This is a largely throw away program aimed to take the output of slabmodel
   and reformat it in a way that it can be used as input to grdtodb.  This 
   can be used to build a regular grid geometry.  

   The algorithm is pretty simple.  We simply eat up all the data, look for the
   smallest length line and write the result out as effectively a matrix of 
   points 

   See usage line but note the outname builds an outname.xyz file and an
   outname.pf file required by grdtodb.  
   */
using namespace std;
void usage()
{
    cerr << "line2grd outname [-g gridname -d outdir]< infile"<<endl
        << "  Note outname is root name.  Two output files created with "
        << "extension .xyz and .pf"<<endl
        << "default gridname=slabsurface, dir=slabmodelgrids"<<endl;;
    exit(-1);
}
int main(int argc, char **argv)
{
    if(argc<2) usage();
    string basename(argv[1]);
    string gridname("slabsurface");
    string outdir("slabmodelgrids");
    int i;
    for(i=2;i<argc;++i)
    {
        string sarg(argv[i]);
        if(sarg=="-g")
        {
            ++i;
            gridname=string(argv[i]);
        }
        else if(sarg=="-d")
        {
            ++i;
            outdir=string(argv[i]);
        }
        else
            usage();
    }
    string datfile,pffile;
    datfile=basename+".xyz";  pffile=basename+".pf";
    fstream dfout,pfout;
    dfout.open(datfile.c_str(), ios::out);
    pfout.open(pffile.c_str(),ios::out);
    if(dfout.fail())
    {
        cerr << "Cannot open output file="<< datfile<<endl;
        usage();
    }
    if(pfout.fail())
    {
        cerr << "Cannot open output file pf def file="<<pffile<<endl;
        usage();
    }
    vector<double> lat,lon,depth,r;
    vector<int> start,end;
    bool newsegment;
    int segsize,nseg;
    char line[256];
    /* Initialize these variables at top of the main read loop */
    i=0;
    segsize=0;
    nseg=0;
    newsegment=true;
    int n1=0;   // n1 is determined from min length line segment
    int n2=0;  // n2 is count of the number of line segments
    while(!cin.eof())
    {
        double latin,lonin,zin,rin;
        cin.getline(line,256);
        if(line[0]=='#') continue;
        if(cin.eof())
        {
            if(segsize==0) 
                break;
            else
            {
                line[0]='>';
                cerr << "Warning:  appears > was missing from last "
                    << "line segment.  Check file"<<endl;
            }
        }
        else if(line[0]=='>')
        {
cout<<"Segment end at count i="<<i<<" This segment size="<<segsize<<endl;
            if(segsize>0)
            {
                end.push_back(i);
                if(n1==0)
                    n1=segsize;
                else
                    n1=min(n1,segsize);
                segsize=0;
                newsegment=true;
                ++n2;
            }
        }
        else
        {
            if(newsegment) 
            {
                start.push_back(i);
                newsegment=false;
            }
            sscanf(line,"%lf%lf%lf%lf",&lonin,&latin,&zin,&rin);
            lat.push_back(latin);
            lon.push_back(lonin);
            depth.push_back(zin);
            r.push_back(rin);
            ++segsize;
            ++i;
        }
    }
    pfout << "n1 " << n1<<endl;
    pfout << "n2 " << n2<<endl;
    pfout << "iorigin 0 "<<endl
        << "jorigin 0" <<endl
        <<"lat0 " << lat[0]<<endl
        <<"lon0 " << lon[0]<<endl
        <<"r0 " << r[0]<<endl
        << "gridname "<<gridname <<endl
        << "grid_directory "<< outdir <<endl
        << "vertical_scale_factor -1.0"<<endl
        << "vertical_offset 0.0" << endl;
    // get nominal spacing from first two point.
    // not important enough that these be accurate to do more
    double delta,azimuth;
    dist(rad(lat[0]),rad(lon[0]),rad(lat[1]),rad(lon[1]),&delta,&azimuth);
    double dxnom=delta*6371.0;
    pfout << "dx1n "<<dxnom<<endl;
    pfout << "dx2n "<<dxnom<<endl;
    int j,jj;
    for(i=0;i<start.size();++i)
    {
        for(j=start[i],jj=0;jj<n1;++jj,++j)
        {
            dfout <<lon[j]<<" "<<lat[j]<<" "<<depth[j]<<endl;
        }
    }
}
        



