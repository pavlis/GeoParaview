#include <stdlib.h>
#include <float.h>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include "dmatrix.h"
#include "gclgrid.h"
#include "seispp.h"
#include "Metadata.h"
using namespace std;
using namespace SEISPP;
void usage()
{
  cerr << "ModEM2gclfield infile outputfieldname [-pf pffile]"<<endl;
  exit(-1);
}
/* Subroutine like procedure to adjust grid depths to irregular 
 * depths stored in z.   Note we assume z is reversed relative to 
 * GCLgrid3d convention (i.e. z increases with index - grid does the
 * reverse */
void adjust_depths(GCLgrid3d& g, vector<double> z)
{
    // Sanity check
    int nz=z.size();
    if(nz!=g.n3)
    {
        cerr << "adjust_depths:   fatal error.  Size mismatch"<<endl
            << "grid n3="<<g.n3<<" but depth vector length passed ="
            << nz<<endl;
        exit(-1);
    }
    Geographic_point gp;
    Cartesian_point cp;
    int k,kk;
    int k0=g.n3-1;  // We assume this point is at r0 - constructor does this
    int i,j;
    /* I do not comprehend why this is necessary. Inheritance should
     * handle this but the gcc compiler I'm using does not allow me
     * to use the derived class */
    BasicGCLgrid *bgg;
    bgg=dynamic_cast<BasicGCLgrid*>(&g);
    for(i=0;i<g.n1;++i)
        for(j=0;j<g.n2;++j)
        {
            gp=bgg->ctog(g.x1[i][j][k0],g.x2[i][j][k0],g.x3[i][j][k0]);
            //gp=g.ctog(g.x1[i][j][k0],g.x2[i][j][k0],g.x3[i][j][k0]);
            for(kk=0,k=g.n3-1;kk<nz;++kk,--k)
            {
                Geographic_point gpnow(gp);
                gpnow.r -= z[kk];
                cp=g.gtoc(gpnow);
                g.x1[i][j][k]=cp.x1;
                g.x2[i][j][k]=cp.x2;
                g.x3[i][j][k]=cp.x3;
            }
        }
}

/* Return a vector of accumulated distances (converted to km) from
the input line.   reverse boolean inverts the order as seems necessary
for this format */
vector<double> parse_coord_axis_line(char *line, int nx, bool reverse)
{
  vector<double> xv;
  xv.reserve(nx);
  stringstream ss(line);
  double x,delta;
  int i;
  for(i=0,x=0.0;i<nx;++i)
  {
    ss >> delta;
    delta *= 0.001;
    xv.push_back(x);
    x += delta;
  }
  if(reverse)
  {
    vector<double> xvr;
    xvr.reserve(nx);
    /* x is currently set to the accumulated sum of cell sizes but it is
    off b one cell so we have to correct by delta */
    x -= delta;
    for(i=0;i<nx;++i)
    {
      xvr.push_back(x-xv[i]);
    }
    return xvr;
  }
  else
    return xv;
}
vector<double> crop_axis(vector<double>& x, int lskip, int rskip)
{
  vector<double> xc;
  int i;
  int nx=x.size();
  for(i=0;i<nx;++i)
  {
    if((i>=lskip) &&(i<(nx-rskip))) xc.push_back(x[i]);
  }
  return xc;
}
/* read a block of data from din as a matrix.  Crop the matrix to nx1cXnx2c */
dmatrix read_and_trim_rho_values(ifstream& din, int nx1, int n1skip, int nx1c,
    int nx2, int n2skip, int nx2c)
{
    try{
      dmatrix result(nx1c,nx2c);
      double rho;
      int i,j,ii,jj;
      for(i=0,ii=0;i<nx1;++i)
      {
        for(j=0,jj=0;j<nx2;++j)
        {
          din >> rho;
          // This allows us to read and skip the bottom borders areas
          if((ii>=nx1c) ||(jj>=nx2c)) continue;
          result(ii,jj)=rho;
          if(j>=n2skip)++jj;
        }
        if(i>=n1skip)++ii;
      }
      return(result);
    }catch(...){throw;};
}
void test_axis(vector<double> x, double dx)
{
  /* use a FLT_EPSILON test for equality */
  int i;
  for(i=1;i<x.size();++i)
  {
    double dxi=fabs(x[i]-x[i-1]);
    double test=(dxi-dx)/dx;
    if(test>FLT_EPSILON)
    {
      cerr << "Grid axis error - only currently support regular grids with irregular boundaries"<<endl
        << "Expected constant dx="<<dx<<" found computed interval of "<<dxi
        << " at index "<<i<<endl
        << "Fatal error:  check parameters or improve this code"<<endl;
      exit(-1);
    }
  }
}
int find_smallest(vector<double>& x)
{
    int i,i0;
    int nx=x.size();
    double xmin;
    xmin=fabs(x[0]);
    for(i0=0,i=1;i<nx;++i)
    {
      if(fabs(x[i])<xmin)
      {
          xmin=fabs(x[i]);
          i0=i;
      }
    }
    return i0;
}
      
bool SEISPP::SEISPP_verbose(true);
int main(int argc, char **argv)
{
  int i,j,k;
  if(argc<3)usage();
  string infile(argv[1]);
  string outfieldname(argv[2]);
  string pffile("ModEM2gclfield");
  for(i=3;i<argc;++i)
  {
    string sarg(argv[i]);
    if(sarg=="-pf")
    {
      ++i;
      if(i>=argc)usage();
      pffile=string(argv[i]);
    }
    else
    {
      usage();
    }
  }
  Pf *pf;
  if(pfread(const_cast<char*>(pffile.c_str()),&pf))
  {
    cerr << "pfread failed on pffile="<<pffile<<endl;
    usage();
  }
  try{
    Metadata control(pf);
    /* these are parameters related to grid geometry.   We read make sure
    these are defined first before doing any real work. */
    /* These are model grid origin - assume we can remap them to a different
    origin using gclfield2vtk */
    double lat0,lon0;
    lat0=control.get_double("origin_latitude");
    lon0=control.get_double("origin_longitude");
    lat0=rad(lat0);
    lon0=rad(lon0);
    /*GCLgrid constructor wants a readius for origin.  Here we freeze it to
    surface, although I think the library will ignore this anyway */
    double r0;
    r0=r0_ellipse(lat0);
    string gridname=control.get_string("gridname");
    string outdir=control.get_string("output_directory");
    int n1lskip,n1rskip,n2lskip,n2rskip;
    n1lskip=control.get_int("number_x1_border_left");
    n1rskip=control.get_int("number_x1_border_right");
    n2lskip=control.get_int("number_x2_border_left");
    n2rskip=control.get_int("number_x2_border_right");
    double dx1expected,dx2expected; //Expected regular grid section spacing in m
    dx1expected=control.get_double("dx1");
    dx2expected=control.get_double("dx2");
    ifstream din;
    din.open(infile.c_str(),ios::in);
    if(din.is_open())
    {
      cout << "Reading model data from file "<<infile;
    }
    else
    {
      cerr << "Open failed on input file "<<infile<<endl;
      usage();
    }
    /* We need a fairly large input line buffer because this file format can have
    potentially very long lines */
    const int INBUFSIZE(8192);
    char linebuf[INBUFSIZE];
    din.getline(linebuf,INBUFSIZE);
    const string firstline(" # 3D MT model written by ModEM in WS format");
    if(firstline!=linebuf)
    {
      cerr << "File format error"<<endl
        << "First line read from file is"<<endl
        << linebuf<<endl<<"Expected"<<endl<<firstline<<endl;
      usage();
    }
    din.getline(linebuf,INBUFSIZE);
    stringstream ss(linebuf);
    int nx1,nx2,nx3;
    ss >> nx1;  ss >> nx2; ss >> nx3;
    int ntmp;
    string test;
    ss >> ntmp;  ss >> test;
    if(test!="LOGE")
    {
      cerr << "Expected to find LOGE tag as 5th token in line 2 of input file"
      <<endl <<"Exiting - assumes something is wrong. Check this file"<<endl;
      usage();
    }
    din.getline(linebuf,INBUFSIZE);
    vector<double> x1_all=parse_coord_axis_line(linebuf,nx1,false);
    din.getline(linebuf,INBUFSIZE);
    vector<double> x2_all=parse_coord_axis_line(linebuf,nx2,true);
    din.getline(linebuf,INBUFSIZE);
    /* Parse the z section straight up as we don't have a border issue there */
    stringstream ssz;
    ssz.str(linebuf);
    vector<double> x3v;
    x3v.reserve(nx3);
    double z(0.0);
    cout << "Grid depths (km)- this is allowed to be nonuniform:"<<endl;
    for(i=0;i<nx3;++i)
    {
      double dz;
      ssz >> dz;
      cout << z <<endl;
      x3v.push_back(z);
      z += (dz*0.001);
    }
    cout<<endl;
    /* When I confirm the format need a procedure here to define cropping */
    int nx1c,nx2c;  // cropped sizes
    vector<double> x1v,x2v;
    x1v=crop_axis(x1_all,n1lskip,n1rskip);
    x2v=crop_axis(x2_all,n2lskip,n2rskip);
    /* this little procedure checks the axis for irregular sampling and
    will abort if any points differ from the expected*/
    double dx1,dx2;
    test_axis(x1v,dx1expected);
    dx1=dx1expected;
    cout << "x1 axis checks out with output uniform grid size="<<dx1<<endl;
    test_axis(x2v,dx2expected);
    dx2=dx2expected;
    cout << "x2 axis checks out with output uniform grid size="<<dx2<<endl;
    nx1c=x1v.size();
    nx2c=x2v.size();
    vector<dmatrix> rho;
    rho.reserve(nx3);
    for(k=0;k<nx3;++k)
    {
      din.getline(linebuf,INBUFSIZE);
      /* as the name says this trims boundary nodes using left skip values */
      rho.push_back(read_and_trim_rho_values(din,
        nx1,n1lskip,nx1c,nx2,n2lskip,nx2c));
      /* Each layer is seprated by a blank line we have to skip */
      din.getline(linebuf,INBUFSIZE);
    }
    /* Next read the coordinate offsets */
    double x10,x20,x30;
    din >> x10;
    din >> x20;
    din >> x30;
    x10*=0.001;  x20*=0.001; x30*=0.001;   // convert to km
    /* this is azimuth of grid - I think we can support this although
    I'm not sure the angles are consistent.  */
    double azgrid;
    din>>azgrid;
    if(azgrid!=0.0)
    {
      cerr << "Warning:   read an azimuth of "<<azgrid<<" from the model file"<<endl
        << "This may yield unexpected results - untested features"<<endl;
    }
    /* The x0 values are a vector from the origin to the upper right
    corner of the grid - weird as that means we add those numbers to
    coordinates to get distance from an origin */
    for(i=0;i<nx1c;++i) x1v[i]+=x10;
    for(i=0;i<nx2c;++i) x2v[i]+=x20;
    for(i=0;i<nx3;++i) x3v[i]+=x30;
    /* We hunt for the closest grid point to the origin and uses that to
    set the i0,j0 offsets for the GCLgrid below */
    int i0,j0;
    i0=find_smallest(x1v);
    j0=find_smallest(x2v);
    cout << "Using grid origin a index i0="<<i0
      <<" which has offset from file origin of "<<x1v[i0]<<endl
      <<"Index j0="<<j0<<" which has an offset from file origin of "
      <<x2v[j0]<<endl;
    /* At this point rho_raw has an image of the file stored as
    matrices with the borders cropped.   Each element of the vector is a
    constant depth layer.   The final step is to translate the grid order
    into a GCLgrid3D and then build the scalarfield object.   First build
    the framework GCLgrid.   dx3nominal is set to dx1 arbitrarily - this is a
    irregular mesh in depth so wrong anyway.  Note also nx1c and nx2c are 
    reversed because of the weird geomtry of this grid.*/
    GCLgrid3d grid(nx2c,nx1c,nx3,gridname,lat0,lon0,r0,M_PI_2+azgrid,
      dx1,dx2,dx1,i0,j0);
    /* Clone this grid to build teh scalarfield object - a large memory model
    but these models are tiny by modern standards */
    GCLscalarfield3d f(grid);
    /* Now will just stuff the rho values in the field slots.  This is really
    weird indexing because we have to invert the indices i and j and particularly
    odd way and run the k index backward. */
    int ii,jj,kk;
    for(i=0,jj=0;i<nx1c;++i,++jj)
    {
      for(j=0,ii=0;j<nx2c;++j,++ii)
      {
        for(k=0,kk=nx3-1;k<nx3;++k,--kk)
        {
          f.val[ii][jj][kk]=rho[k](i,j);
          /* Data are stored as log e - convert top physical units */
          f.val[ii][jj][kk]=exp(f.val[ii][jj][kk]);
        }
      }
    }
    /* Last thing we need to do is adjust the depths.   The constructor
     * builds a regular depth grid.   This procedure adjusts the grid
     * geometry for irregular depth */
    adjust_depths(f,x3v);
    cout << "Writing data to directory="
       <<outdir<<" with base name="<<outfieldname<<endl;
    f.save(outfieldname,outdir);
  }catch(SeisppError& serr)
  {
    serr.log_error();
    exit(-2);
  }
  catch(std::exception& stexc)
  {
    cerr << stexc.what()<<endl;
    exit(-2);
  }
}
