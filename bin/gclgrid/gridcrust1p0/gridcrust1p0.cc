#include <iostream>
#include "stock.h"
#include "dbpp.h"
#include "Metadata.h"
#include "VelocityModel_1d.h"
#include "Crust1_0.h"
GCLgrid *makegrid(Metadata& p)
{
    /* The parameters used here are identical to makegclgrid 
       for 2d grids.  */
    string gridname=p.get_string("gridname");
    double delta_x1,delta_x2;
    delta_x1=p.get_double("delta_x1");
    delta_x2=p.get_double("delta_x2");
    double az=p.get_double("x_axis_azimuth");
    double lat0,lon0;
    lat0=p.get_double("origin_latitude");
    lon0=p.get_double("origin_longitude");
    /* These need to be converted to radians for the constructor*/
    az=rad(az-90.0);  // constructor want az x2
    lat0=rad(lat0);
    lon0=rad(lon0);
    int n1,n2;
    n1=p.get_int("n1");
    n2=p.get_int("n2");
    /* This grid is created like makegclgrid as at the surface.
    Here we drop the origin offset variables ( 0,0 of last two
    arguments) because grid will be remapped after creation. */
    GCLgrid *g=new GCLgrid(n1,n2,gridname,lat0,lon0,r0_ellipse(lat0),
            az,delta_x1,delta_x2,0,0);
    return(g);
}

void usage()
{
    cerr << "gridcrust1p0 db [-g gridname -dbout -pf pffile]"<<endl
        << " Use -g to use gridname (default created from pf input"<<endl
        << " -dbout will save results to db (default is files)"<<endl;
    exit(-1);
}
bool SEISPP::SEISPP_verbose(false);
int main(int argc, char **argv)
{
    if(argc<2) usage();
    string pffile("gridcrust1p0.pf");
    string dbname(argv[1]);
    bool load_grid_from_db(false);
    bool save_to_db(false);
    string gridname;
    int i,j;
    for(i=2;i<argc;++i)
    {
        string sarg(argv[i]);
        if(sarg=="-pf")
        {
            ++i;
            if(i>=argc) usage();
            pffile=string(argv[i]);
        }
        else if(sarg=="-g")
        {
            ++i;
            if(i>=argc) usage();
            gridname=string(argv[i]);
            load_grid_from_db=true;
        }
        else if(sarg=="-dbout")
            save_to_db=true;
        else
            usage();
    }
    Pf *pf;
    if(pfread(const_cast<char *>(pffile.c_str()),&pf))
    {
        cerr << "pfread failed on pffile="<<pffile<<endl;
        exit(-1);
    }
    try{
        Metadata control(pf);
        DatascopeHandle dbh(dbname,false);
        Crust1_0 crust1p0;
        GCLgrid *g;
        if(load_grid_from_db)
            g=new GCLgrid(dbh,gridname);
        else
            g=makegrid(control);
        /* We always define the coordinate system through this mechanism
           and then use the remap_grid procedure to do the transformation */
        double lat0,lon0,r0;
        lat0=control.get_double("remap_origin_latitude");
        lon0=control.get_double("remap_origin_longitude");
        lat0=rad(lat0);
        lon0=rad(lon0);
        r0=control.get_double("remap_origin_radius");
        double az=control.get_double("remap_azimuth_x");
        az=az-90.0;
        az=rad(az);
        remap_grid(g,lat0,lon0,r0,az);
        /* Now we construct a series of GCLgrids from the Crust1.0 data.
           We write a set of fields for limited horizons. Start with
         the Moho.*/
        string base_name=control.get_string("root_output_name");
        const int dbnamelimit(15),max_ext_size(5);;
        int maxnamesize=dbnamelimit-max_ext_size;
        /* list name size if using db to avoid problems */
        if(save_to_db)
        {
            if(base_name.size()>maxnamesize)
            {
                base_name=base_name.substr(0,maxnamesize);
                cerr << "WARNING:  base name for output fields truncated to "<<base_name<<endl
                    << "String found in pf file would have overflowed fieldname field in db"<<endl;
            }
        }
        string outdir=control.get_string("output_directory");
        GCLscalarfield f(*g);
        /* We can clear g now */
        delete g;
        for(j=0;j<f.n2;++j)
            for(i=0;i<f.n1;++i)
            {
                Geographic_point gp=f.geo_coordinates(i,j);
                /* Crust1_0 object uses degrees not radians as used in gclgrid library*/
                double dlat=deg(gp.lat);
                double dlon=deg(gp.lon);
                double r0=r0_ellipse(gp.lat);
                gp.r = r0 - crust1p0.CrustalThickness(dlat,dlon);
                /* change the point's coordinates.  This algorithm effectively projects
                grid to surface and then projects down.   Not important when generated
                internally, but potentially an issue if grid is read from db */
                Cartesian_point cp=f.gtoc(gp);
                f.x1[i][j]=cp.x1;
                f.x2[i][j]=cp.x2;
                f.x3[i][j]=cp.x3;
                VelocityModel_1d vp=crust1p0.model(dlat,dlon,string("Pvelocity"));
                /* The object interface strips zero thickness layers in crust1.0 so
                   we have to use the last layer velocity for moho velocity */
                f.val[i][j]=vp.v[vp.nlayers-1];
            }
        string fieldname;
        fieldname=base_name+"_MVp";
        if(save_to_db)
            f.save(dbh,outdir,outdir,fieldname,fieldname);
        else
            f.save(fieldname,outdir);
        /* For Vs and rho we only need to extract alternative data so this loop 
           is simpler - no change to coordinates. Painfully repetitious code here, but
           turnin this into a procedure would only confuse this*/
        for(j=0;j<f.n2;++j)
            for(i=0;i<f.n1;++i)
            {
                Geographic_point gp=f.geo_coordinates(i,j);
                double dlat=deg(gp.lat);
                double dlon=deg(gp.lon);
                VelocityModel_1d vs=crust1p0.model(dlat,dlon,string("Svelocity"));
                f.val[i][j]=vs.v[vs.nlayers-1];
            }
        fieldname=base_name+"_MVs";
        if(save_to_db)
            f.save(dbh,outdir,outdir,fieldname,fieldname);
        else
            f.save(fieldname,outdir);
        for(j=0;j<f.n2;++j)
            for(i=0;i<f.n1;++i)
            {
                Geographic_point gp=f.geo_coordinates(i,j);
                double dlat=deg(gp.lat);
                double dlon=deg(gp.lon);
                VelocityModel_1d rho=crust1p0.model(dlat,dlon,string("Density"));
                f.val[i][j]=rho.v[rho.nlayers-1];
            }
        fieldname=base_name+"_Mrho";
        if(save_to_db)
            f.save(dbh,outdir,outdir,fieldname,fieldname);
        else
            f.save(fieldname,outdir);
        /* Finally do depth as a field.   Output depth in km */
        for(j=0;j<f.n2;++j)
            for(i=0;i<f.n1;++i)
            {
                Geographic_point gp=f.geo_coordinates(i,j);
                double dlat=deg(gp.lat);
                double dlon=deg(gp.lon);
                f.val[i][j]=crust1p0.CrustalThickness(dlat,dlon);
            }
        fieldname=base_name+"_depthmoho";
        if(save_to_db)
            f.save(dbh,outdir,outdir,fieldname,fieldname);
        else
            f.save(fieldname,outdir);
    }catch(std::exception err)
    {
        cerr << err.what()<<endl;
    }
    catch(SeisppError& serr)
    {
        serr.log_error();
    }
}

