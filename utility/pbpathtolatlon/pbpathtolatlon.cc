#include <stdio.h>
#include "coords.h"
#include "PlateBoundaryPath.h"
#include "SeisppError.h"
using namespace SEISPP;
void usage()
{
    cerr << "pbpathtolatlon polelat polelon < slabprojectfile " <<endl;
    exit(-1);
}
int main(int argc, char **argv)
{
    if(argc!=3) usage();
    try {
    double plat,plon;
    plat=atof(argv[1]);
    plon=atof(argv[2]);
    plat=rad(plat);
    plon=rad(plon);
    int npoints;
    double olat,olon;
    double distin,depth,odepth;
    double angvel=1.0;  // made unity to simify conversion to distance
    while(fscanf(stdin,"%d %lf %lf %lf",&npoints,&olat,&olon,&odepth)!=EOF)
    {
        // Echo the first point immediately before converting it 
        cout << olat <<" "<<olon<<" "<<odepth<<endl;
        olat=rad(olat);  olon=rad(olon);
        double Reffective,aztmp;
        dist(plat,plon,olat,olon,&Reffective,&aztmp);
        Reffective=deg2km(deg(Reffective));
        double dist2s=1.0/(Reffective*angvel);
        PlateBoundaryPath path(plat,plon,olat,olon,angvel);
        Geographic_point nextpoint;
        for(int i=0;i<npoints;++i)
        {
            if(fscanf(stdin,"%lf%lf",&distin,&depth)==EOF) 
            {
                cerr << "Read error encountered prematurely.  Check count"
                    <<endl;
                exit(-1);
            }
            nextpoint=path.position(distin*dist2s);
            cout << deg(nextpoint.lat) <<" " 
                << deg(nextpoint.lon) << " "
                << depth << endl;
	}
    }
    }catch(SeisppError& serr)
    {
        serr.log_error();
    }
}
