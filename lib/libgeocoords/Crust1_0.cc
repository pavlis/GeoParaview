#include <iostream>
#include <fstream>
#include "stock.h"
#include "Crust1_0.h"
#include "gclgrid.h"
using namespace std;
using namespace SEISPP;
/* this helper is a key component of the main constructor */
void Crust1_0::load_crust1file(double ***y,string fname)
{
    int i,j,k;
    try {
        ifstream in;
        in.open(fname.c_str(),std::ifstream::in);
        /* Note we invert the latitude order to make grid start 
           at the south pole.  The data file start at 89.5 degrees
           north.  This makes the origin at -89.5 degrees.*/
        for(j=nlat-1;j>0;--j)
            for(i=0;i<nlon;++i)
                for(k=0;k<nlayers;++k)
                    in >> y[i][j][k];
        in.close();
    }catch(...){throw;};
}
Crust1_0::Crust1_0()
{
    /* Could not how to initialize these in the .h file without generating a warning
       with gcc. These would be more appropriate as const parameteters */
    nlon=360;
    nlat=180;
    nlayers=9;
    dlon=1.0;
    dlat=1.0;
    /* Note the data in the file actually start at N 89.5.   I reverse the order of storage
       to mesh more easily with the bilinear interpolator and prevent the confusing (to me 
       anyway) upside down y axis. */
    lat0=-89.5;
    lon0=-179.5;
    const string base_error("Crust1_0 constructor:  ");
    char *vpfile=datapath(NULL,"tables/crust1.0","crust1","vp");
    if(vpfile==NULL) 
        throw SeisppError(base_error+"datapath failed search for crust1.vp file");
    char *vsfile=datapath(NULL,"tables/crust1.0","crust1","vs");
    /* This and the next 2 similar tests are very unlikely to every be
       thrown.  If any of these throw they will generate a small memory 
       leak, but the extra coding required to fix that is not justified.*/
    if(vsfile==NULL) 
        throw SeisppError(base_error+"datapath failed search for crust1.vs file");
    char *rhofile=datapath(NULL,"tables/crust1.0","crust1","rho");
    if(rhofile==NULL) 
        throw SeisppError(base_error+"datapath failed search for crust1.rho file");
    char *bndsfile=datapath(NULL,"tables/crust1.0","crust1","bnds");
    if(bndsfile==NULL) 
        throw SeisppError(base_error+"datapath failed search for crust1.bnds file");
    /* crust1.0 has a frozen structure with nothing in the file at all 
       linked to critical parameters.  In particular, the size of the grid
       is frozen and not in the file.   Here we deal with that using
       const int values in the private area of the object.   We use
       these now to create the work areas. */
    try{
        vp=create_3dgrid_contiguous(nlon,nlat,nlayers);
        if(vp==NULL) throw SeisppError(base_error
                + "alloc failed for vp grid array");
        vs=create_3dgrid_contiguous(nlon,nlat,nlayers);
        if(vs==NULL) throw SeisppError(base_error
                + "alloc failed for vs grid array");
        rho=create_3dgrid_contiguous(nlon,nlat,nlayers);
        if(rho==NULL) throw SeisppError(base_error
                + "alloc failed for rho grid array");
        depth=create_3dgrid_contiguous(nlon,nlat,nlayers);
        if(depth==NULL) throw SeisppError(base_error
                + "alloc failed for depth grid array");
        /* This will have an extra component, but easier to leave
           to allow interpolator to work with this array like
           the other arrays */
        thickness=create_3dgrid_contiguous(nlon,nlat,nlayers);
        if(thickness==NULL) throw SeisppError(base_error
                + "alloc failed for thickness grid array");
    /* Load the model files using the private method load_crust1file */
        load_crust1file(vp,vpfile);
        load_crust1file(vs,vsfile);
        load_crust1file(rho,rhofile);
        load_crust1file(depth,bndsfile);
        /* For efficiency we compute thicknesses at each node as
           there are frequent zero thickness layers in this beast.
           Note the difference in count of 1 for intervals 
           compared to points.  Thickness references the layer 
           below each top.   Also note depth is a bit of a 
           misnomer as the boundaries are tabulated with positive
           up (Moho always has a negative depth value).*/
        for(int i=0;i<nlon;++i)
            for(int j=0;j<nlat;++j)
            {
                for(int k=0;k<8;++k)
                    thickness[i][j][k]=depth[i][j][k]
                                        - depth[i][j][k+1];
                /* set to large value to allow a test for
                   zero thickness to not incorrectly skip */
                thickness[i][j][8]=6000.0;
            }

        free(vpfile);
        free(vsfile);
        free(rhofile);
        free(bndsfile);
    }catch(...){throw;};
}
int Crust1_0::lookup_lon(double lon)
{
    /* Force lon to be between -180 and 180.   A bit excessive 
    here deal with absurd requests */
    int i;
    const int maxtries(10);
    i=0;
    double maxlon=lon0+static_cast<double>(nlon-1)*dlon;
    while(maxlon>180.0 && i<maxtries)
    {
        lon -= 360.0;
        ++i;
    }
    i=0;
    while(lon<lon0 && i<maxtries)
    {
        lon += 360.0;
        ++i;
    }
    int result=(int) (lon-lon0)/dlon;
    return result;
}
int Crust1_0::lookup_lat(double lat)
{
    int result=(int) (lat-lat0)/dlat;
    /* nlat-1 used because we require finding point at lower left of 
       grid square to interpolate.*/
    if( (result<0) || (result>(nlat-1)) ) 
        throw SeisppError(string("Crust1_0::lookup_lat: illegal latitude"));
    return result;
}
/* Use bilinear interpolator for model components.  Complication 
   is that Crust1.0 has zero entries in places that we have to 
   work around.   To do that we check for any zeros in the 4 corners
   and revert to nearest neighbor if any are zero.  Will return
   zero if zero node is closest.

    Note dlon and dlat are normalized delta latitude and delta longitude
    That is, range is 0 to 1 */
double *interpolate_model_grid(double ***y,int i, int j, 
        double dlon, double dlat)
{
    /* The nearest neighbor algorithm requires this test to avoid
       disaster */
    if((dlon<0.0) || (dlon>1.01) || (dlat<0.0) || (dlat>1.01))
        throw SeisppError(string("Crust1_0:  interpolate_model_grid")
              + "illegal delta latitude and/or delta longitude values passed");
    double *result=new double[9];
    /* We blindly assume the i and j values are valid */
    int k,ii,jj;
    bool has_a_zero(false);
    for(k=0;k<9;++k)
    {
        for(ii=i;ii<i+2;++ii)
            for(jj=j;jj<j+2;++jj)
                if(fabs(y[ii][jj][k])<0.01)
                {
                    has_a_zero=true;
                    break;
                }

        if(has_a_zero)
        {
            result[k]=y[i+SEISPP::nint(dlon)][j+SEISPP::nint(dlat)][k];
        }
        else
        {
            result[k]=bilinear<double>(dlon,dlat,y[i][j][k],y[i+1][j][k],
                y[i][j+1][k],y[i+1][j+1][k]);
        }
    }
    return result;
}
VelocityModel_1d Crust1_0::model(double lat, double lon, string property)
{
    try {
        int i,j;
        i=lookup_lon(lon);
        j=lookup_lat(lat);
        double deltalat,deltalon;
        deltalat=(lat - (lat0+dlat*static_cast<double>(j)))/dlat;
        deltalon=(lon - (lon0+dlon*static_cast<double>(i)))/dlon;
        double *dz,*y;
        dz=interpolate_model_grid(thickness,i,j,deltalon,deltalat);
        double *layertops;
        layertops=interpolate_model_grid(depth,i,j,deltalon,deltalat);
        if(property=="Pvelocity" || property=="Vp")
            y=interpolate_model_grid(vp,i,j,deltalon,deltalat);
        else if(property=="Svelocity" || property=="Vp")
            y=interpolate_model_grid(vs,i,j,deltalon,deltalat);
        else if(property=="Density" || property=="rho")
            y=interpolate_model_grid(rho,i,j,deltalon,deltalat);
        VelocityModel_1d result(9);
        int k;
        for(k=0;k<9;++k)
        {
            /* Skip any zero thickness layers and those with zero
               propeties.   (tested against a small value).
               This happens frequently in crust1.0.
               Magic numbers assume units of km and seconds.
               */
            if(((fabs(dz[k])>0.001)) && (fabs(y[k]>0.01)))
            {
                result.v.push_back(y[k]);
                /* Reverse sign of layertops because the VelocityModel_1d
                   object uses depth positive while the crust1.0 boundary
                   file tabulates negative depths (elevations) */
                result.z.push_back(-layertops[k]);
                result.grad.push_back(0.0);
            }
        }
        //DEBUG
        /*
        for(k=0;k<result.v.size();++k)
            cout << k <<" "<< result.v[k]<<" "<<result.z[k]<<endl;
            */
        result.nlayers=result.v.size();
        delete [] dz;
        delete [] y;
        delete [] layertops;
        return result;
    }catch(...){throw;};
}
double Crust1_0::CrustalThickness(double lat, double lon)
{
    try {
        int i,j;
        i=lookup_lon(lon);
        j=lookup_lat(lat);
        double deltalat,deltalon;
        deltalat=(lat - (lat0+dlat*static_cast<double>(j)))/dlat;
        deltalon=(lon - (lon0+dlon*static_cast<double>(i)))/dlon;
        double *tops;
        tops=interpolate_model_grid(depth,i,j,deltalon,deltalat);
        double dz=tops[0]-tops[8];
        delete [] tops;
        return(dz);
    }catch(...){throw;};
}
double Crust1_0::MohoDepth(double lat, double lon)
{
    try {
        int i,j;
        i=lookup_lon(lon);
        j=lookup_lat(lat);
        double deltalat,deltalon;
        deltalat=(lat - (lat0+dlat*static_cast<double>(j)))/dlat;
        deltalon=(lon - (lon0+dlon*static_cast<double>(i)))/dlon;
        double *tops;
        tops=interpolate_model_grid(depth,i,j,deltalon,deltalat);
        double z=-tops[8];
        delete [] tops;;
        return(z);
    }catch(...){throw;};
}
double Crust1_0::SedimentThickness(double lat, double lon)
{
    try {
        int i,j;
        i=lookup_lon(lon);
        j=lookup_lat(lat);
        double deltalat,deltalon;
        deltalat=(lat - (lat0+dlat*static_cast<double>(j)))/dlat;
        deltalon=(lon - (lon0+dlon*static_cast<double>(i)))/dlon;
        double *tops;
        tops=interpolate_model_grid(depth,i,j,deltalon,deltalat);
        double dz=tops[2]-tops[5];
        delete [] tops;
        return(dz);
    }catch(...){throw;};
}
double Crust1_0::WaterDepth(double lat, double lon)
{
    try {
        int i,j;
        i=lookup_lon(lon);
        j=lookup_lat(lat);
        double deltalat,deltalon;
        deltalat=(lat - (lat0+dlat*static_cast<double>(j)))/dlat;
        deltalon=(lon - (lon0+dlon*static_cast<double>(i)))/dlon;
        double *tops;
        tops=interpolate_model_grid(depth,i,j,deltalon,deltalat);
        double z=tops[0]-tops[1];
        delete [] tops;
        return(z);
    }catch(...){throw;};
}
double Crust1_0::IceThickness(double lat, double lon)
{
    try {
        int i,j;
        i=lookup_lon(lon);
        j=lookup_lat(lat);
        double deltalat,deltalon;
        deltalat=(lat - (lat0+dlat*static_cast<double>(j)))/dlat;
        deltalon=(lon - (lon0+dlon*static_cast<double>(i)))/dlon;
        double *tops;
        tops=interpolate_model_grid(depth,i,j,deltalon,deltalat);
        /* It appears we have to use the second field or we get bogus
           ice thickness values in oceanic crust. */
        double dz=tops[1] - tops[2];
        delete [] tops;;
        return(dz);
    }catch(...){throw;};
}
