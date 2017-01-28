#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "Metadata.h"
#include "gclgrid.h"
void usage()
{
  cerr << "ArcGISgrd2gcl infile outfilebase [-zone xx -gridname gn]" <<endl
    << "Converts txt format regular grid from ArcGIS export to 2D GCLfield"<<endl
    << "infile is the file created by ArcGIS"<<endl
    << "outfilebase is the base name for the converted output gclfield data"
    <<endl
    << "Output is a masked field with depth as the scalar variable"<<endl
    << "Use -zone to change the UTM zone code to xx"<<endl
    << "Use -gridname to change base grid name (default PCStructure)"
    <<endl;
  exit(-1);
}
int main(int argc, char **argv)
{
  int i,j;
  if(argc<3)usage();
  string infile(argv[1]);
  string outbase(argv[2]);
  string utmzone("XXX");   // Need to get this from Dante email
  string gridname("PCStructure");
  for(i=3;i<argc;++i)
  {
    string sarg(argv[i]);
    if(sarg=="-zone")
    {
      ++i;
      if(i>=argc)usage();
      utmzone=string(argv[i]);
    }
    else f(sarg=="-zone")
    {
      ++i;
      if(i>=argc)usage();
      gridname=string(argv[i]);
    }
    else
    {
      usage();
    }
  }
  ifstream ifs;
  const int number_header_lines(6);
  /* This was written for a specific file.  Unclear if arcgis adds other
  key,value pairs than the 6 the file I'm converting uses.  If I knew there
  were other posibilities I would make number_header_lines a command line
  argument.  In any case, we blindly read them and convert them to a
  metadata object.*/
  ifs.open(infile.c_str(),ios::in);
  if(ifs.fail())
  {
    cerr << "ArcGISgrd2gcl:  cannot open input file "<<infile<<endl;
    usage();
  }
  string headerlines;
  for(i=0;i<number_header_lines;++i)
  {
    string linebuf;
    getline(ifs,linebuf);
    headerlines=headerlines+"\n"+linebuf;
  }
  Metadata header(headerlines);
  int nx1=header.get_int("ncols");
  int nx2=header.get_int("nrows");
  double xll=header.get_double("xllcorner");
  double yll=header.get_double("yllcorner");
  double cellsize=header.get_double("cellsize");
  double NODATAvalue=header.get_double("NODATA_value");
  /* To avoid a float equal compare we derive a floor from NODATAvalue and
  test instead for any value smaller than the (negative) NODATAvalue.*/
  double NDfloor=NODATAvalue+1.0;
  /* A sanity check*/
  if(NDfloor>0.0)
  {
    cerr << "ArcGISgrd2gcl:  parsed NODATA_value="<<NODATAvalue
       << " is nonsense"<<endl<<"Expected to be a large negative number"<<endl
       << " Assumptions of this program are violated - bug fix required"
       <<endl
       << "Cannot continue"<<endl;
    exit(-1);
  }
  /* Large memory algorithm - eat up the entire file and make sure the
  size matches what is expected from size parameters */
  vector<double> dvector;
  dvector.reserve(nx1*nx2);
  do
  {
    double dval;
    ifs >> dval;
    if(ifs.fail())break;
    dvector.push_back(dval);
  }while(1);
  if(dvector.size() != (nx1*nx2))
  {
    cerr << "Size mismatch in data read from file="<<infile<<endl
      << "Number of values read from file="<<dvector.size()<<endl
      << "Header parameters give ncol="<<nx1<<" and nrows="<<nx2
      << " which means number expected="<<nx1*nx2<<endl;
    exit(-1);
  }
  /* We have to first build a regular GCLgrid before we apply the mask.
  We assume all the dvector values are in units of m.  horizontal coordinates
  are assumed utm in m.  Dante thinks the file starts in the upper left corner
  but I am going guess it is actually the lower left given the parameter
  names */
  double lat0,lon0,r0;   //frozen origin here to surface r0 at lower left corner
  /* 23 in this call freezes WGS-84 - not at all elegant*/
  UTMtoLL(23,yll,xll,utmzone.c_str(),&lat0,&lon0);
  lat0 *= deg2rad;
  lon0 *= deg2rad;
  r0=r0_ellipse(lat0);
  GCLgrid g(nx1,nx2,gridname,lat0,lon0,r0,0.0,cellsize*0.001,cellsize*0.001,0,0);
  GCLscalarfield f(g);
  double easting,northing,elev;
  double lat,lon,r;
  Cartesian_point cp;
  int iv;
  for(i=0,iv=0;i<nx1;++i)
  {
    for(j=0;j<nx2;++j)
    {
      easting=xll+cellsize*((double)i);
      northing=yll+cellsize*((double)j);
      elev=dvector[iv];
      UTMtoLL(23,northing,easting,utmzone.c_str(),&lat,&lon);
      lat*=deg2rad;
      lon*=deg2rad;
      /*Zero masked values */
      if(dvector[iv]<NDfloor)
      {
        elev=0;
        f.val[i][j]=0.0;
      }
      else
      {
        elev=dvector[iv]/1000.0;
        f.val[i][j]=elev;
      }
      r=r0_ellipse(lat)-elev;
      cp=g.gtoc(lat,lon,r);
      f.x1[i][j]=cp.x1;
      f.x2[i][j]=cp.x2;
      f.x3[i][j]=cp.x3;
      ++iv;
    }
  }
  g.compute_extents();
  GCLMask gmask(g);
  /* Now apply the mask */
  for(i=0,iv=0;i<nx1;++i)
  {
    for(j=0;j<nx2;++j)
    {
      if(dvector[iv]<NDfloor)
        gmask.mask_point(i,j);
      ++iv;
    }
  }
  GCLMaskedScalarField msf(f,gmask);
  msf.save(outfilebase);
}
