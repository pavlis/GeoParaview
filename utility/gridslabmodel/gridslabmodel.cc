#include <iostream>
#include <vector>
#include "stock.h"
#include "coords.h"
#include "dbpp.h"
#include "gclgrid.h"
#include "seispp.h"
#include "SeisppError.h"
using namespace std;
using namespace SEISPP;
void usage()
{
    cerr << "gridslabmodel db patterngrid outputgrid < slabmodeldata "<<endl
        <<"  patterngrid must be a 3d grid object defined in gclgdisk table"<<endl;
    exit(-1);
}
/* This should be in gclgrid library, but was copied from another program */
void copy_master(GCLgrid3d& gin, GCLgrid& gout,string outname)
{
        gout.name=outname;
        gout.lat0=gin.lat0;
        gout.lon0=gin.lon0;
        gout.r0=gin.r0;
        gout.azimuth_y=gin.azimuth_y;
        gout.lat0=gin.lat0;
        for(int i=0;i<3;++i)
        {
                 for(int j=0;j<3;++j) 
                        gout.gtoc_rmatrix[i][j]=gin.gtoc_rmatrix[i][j];
                 gout.translation_vector[i]=gin.translation_vector[i];
        }
}
bool SEISPP::SEISPP_verbose(true);
int main(int argc, char **argv)
{
    if(argc!=4) usage();
    string dbname(argv[1]);
    string patterngrid(argv[2]);
    string outgrid(argv[3]);
    int i,j;
    bool use2d(false);
    int n1,n1this,n2;
    //n2 is count of second dimension
    n1=0;   n2=0;   n1this=0;
    char line[512];
    vector<double> lat, lon, depth;
    while(cin.getline(line,512))
    {
        double lattmp,lontmp,dtmp;
        if(line[0]=='>') 
        {
            if(n2==0)
                n1=n1this;
            else
            {
                if(n1!=n1this)
                {
                    cerr << "Irregular scan size.  First grid line length="<<n1<<endl
                        << "Current grid line ("<<n2<<") length="<<n1this<<endl
                        << "Cannot continue."<<endl;
                    exit(-1);
                }
            }
            ++n2;
            continue;
        }
        else
        {
            stringstream ss(line);
            ss >> lontmp;
            ss >> lattmp;
            ss >> dtmp;
            lat.push_back(rad(lattmp));
            lon.push_back(rad(lontmp));
            depth.push_back(dtmp);
            ++n1this;
        }
    }
    //DEBUG
    cerr << "n1="<<n1<<" n2="<<n2<<endl;
    try {
        DatascopeHandle dbh(dbname,false);
        dbh.lookup("gclgdisk");
        GCLgrid3d pattern(dbh.db,patterngrid);
        GCLgrid gout(n1,n2);
        copy_master(pattern,gout,outgrid);
        Geographic_point gp;
        Cartesian_point cp;
        int igp;
        for(i=0,igp=0;i<n1;++i)
            for(j=0;j<n2;++j)
            {
                gp.lat=lat[igp];
                gp.lon=lat[igp];
                gp.r=r0_ellipse(gp.lat)-depth[igp];
                cp=gout.gtoc(gp);
                gout.x1[i][j]=cp.x1;
                gout.x2[i][j]=cp.x2;
                gout.x3[i][j]=cp.x3;
                ++igp;
            }
        gout.compute_extents();
        gout.dx1_nom=fabs(gout.x1high-gout.x1low)/static_cast<double>(n1);
        gout.dx2_nom=fabs(gout.x2high-gout.x2low)/static_cast<double>(n2);
        // for now this is frozen
        string outdir("slabmodelgrids");
        gout.dbsave(dbh.db,outdir);
    }catch (int ierr)
    {
        cerr << "GCLgrid error:  something threw error number "
            << ierr <<endl;
    }
    catch (SeisppError& serr)
    {
        serr.log_error();
    }
}


