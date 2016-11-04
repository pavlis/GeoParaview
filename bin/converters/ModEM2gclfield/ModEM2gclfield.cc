#include <stdlib.h>
#include <float.h>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include "seispp.h"
#include "Metadata.h"
#include "gclgrid.h"
using namespace std;
using namespace SEISPP;
void usage()
{
  cerr << "ModEM2gclfield infile outputfieldname [-pf pffile]"<<endl;
  exit(-1);
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
    double dxi=x[i]-x[i-1];
    double test=fabs((dxi-dx)/dx);
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
    const string firstline("# 3D MT model written by ModEM in WS format");
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
    ss.str(linebuf);
    vector<double> x3v;
    x3v.reserve(nx3);
    double z(0.0);
    cout << "Grid depths (km)- this is allowed to be nonuniform:"<<endl;
    for(i=0;i<nx3;++i)
    {
      double dz;
      ss >> dz;
      cout << z <<" ";
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
    cout << "x1 axis checks out with output uniform grid size="<<dx1<<endl;
    dx1=dx1expected;
    test_axis(x2v,dx2expected);
    cout << "x2 axis checks out with output uniform grid size="<<dx2<<endl;
    dx2=dx2expected;
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
    double x1min,x2min;
    for(i=1,x1min=fabs(x1v[0]),i0=0;i<nx1c;++i)
    {
      if(fabs(x1v[i])<x1min)
      {
        x1min=x1v[i];
        i0=i;
      }
    }
    for(i=1,x2min=fabs(x2v[0]),j0=0;i<nx2c;++i)
    {
      if(fabs(x1v[i])<x1min)
      {
        x2min=x2v[i];
        j0=i;
      }
    }
    cout << "Using grid origin a index i0="<<i0
      <<" which has offset from file origin of "<<x1min<<endl
      <<"Index j0="<<j0<<" which has an offset from file origin of "<<x2min<<endl;
    /* At this point rho_raw has an image of the file stored as
    matrices with the borders cropped.   Each element of the vector is a
    constant depth layer.   The final step is to translate the grid order
    into a GCLgrid3D and then build the scalarfield object.   First build
    the framework GCLgrid.   dx3nominal is set to dx1 arbitrarily - this is a
    irregular mesh in depth so wrong anyway.*/
    GCLgrid3d grid(nx1c,nx2c,nx3,gridname,lat0,lon0,r0,M_PI_2+azgrid,
      dx1,dx2,dx1,i0,j0);
    /* Clone this grid to build teh scalarfield object - a large memory model
    but these models are tiny by modern standards */
    GCLscalarfield3d f(grid);
    /* Now will just stuff the rho values in the field slots.  This is really
    weird indexing because we have to invert the indices i and j and particularly
    odd way and run the k index backward. */
    int ii,jj,kk;
    for(i=0,ii=0;i<nx1c;++i,++ii)
    {
      for(j=0,jj=nx2c-1;j<nx2c;++j,--jj)
      {
        for(k=0,kk=nx3-1;k<nx3;++k,--kk)
        {
          f.val[ii][jj][kk]=rho[k](i,j);
          /* Data are stored as log e - convert top physical units */
          f.val[ii][jj][kk]=exp(f.val[ii][jj][kk]);
        }
      }
    }
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
