#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "seispp.h"
#include "Metadata.h"
#include "gclgrid.h"
#include "GCLMasked.h"
#include "LatLong-UTMconversion.h"
void usage()
{
  cerr << "ArcGISgrd2gcl infile outfilebase [-zone xx "<<endl
      <<"-gridname gn -offset xx -scale xx -gmt gmtout -v]" <<endl
    << "Converts txt format regular grid from ArcGIS export to 2D GCLfield"<<endl
    << "infile is the file created by ArcGIS"<<endl
    << "outfilebase is the base name for the converted output gclfield data"
    <<endl
    << "Output is a masked field with depth as the scalar variable"<<endl
    << "Use -zone to change the UTM zone code to xx"<<endl
    << "Use -gridname to change base grid name (default PCStructure)"<<endl
    << "Use -offset subtract a dc value from all nonnull data (default 0)"<<endl
    << "Use -scale to change scale (default is ft to m conversion)"<<endl
    << "Note:  offset is applied before scale (z=(zraw-offet)*scale)"<<endl
    << "Use -gmt to write file gmtout for import into gmt with xyz2grd"<<endl
    << "(Note xyz2grd does not require a fill for NULL values - just skipped in output)"<<endl
    << "Use -v for verbose (default is silent)"<<endl
    <<endl;
  exit(-1);
}
bool SEISPP::SEISPP_verbose(false);
int main(int argc, char **argv)
{
  int i,j;
  if(argc<3)usage();
  string infile(argv[1]);
  string outbase(argv[2]);
  /* This is correct for great unconformity map.  Letter is superflous. */
  string utmzone("14S");   
  string gridname("PCStructure");
  double zoffset(50000.0);
  double rawz2elev(0.0003048);
  bool write_gmt_output(false);
  string gmtoutputfile;
  for(i=3;i<argc;++i)
  {
    string sarg(argv[i]);
    if(sarg=="-zone")
    {
      ++i;
      if(i>=argc)usage();
      utmzone=string(argv[i]);
    }
    else if(sarg=="-gridname")
    {
      ++i;
      if(i>=argc)usage();
      gridname=string(argv[i]);
    }
    else if(sarg=="-offset")
    {
      ++i;
      if(i>=argc)usage();
      zoffset=atof(argv[i]);
    }
    else if(sarg=="-scale")
    {
      ++i;
      if(i>=argc)usage();
      rawz2elev=atof(argv[i]);
    }
    else if(sarg=="-gmt")
    {
      ++i;
      if(i>=argc)usage();
      write_gmt_output=true;
      gmtoutputfile=string(argv[i]);
    }
    else if(sarg=="-v")
        SEISPP_verbose=true;
    else
    {
      usage();
    }
  }
  if(SEISPP_verbose)
  {
      cout << "ArcGISgrd2gcl converting data in file "<<infile<<endl
          << "Writing results to GCLMaskedScalarField file with base name "
          << outbase<<endl
          << "Data offset="<<zoffset<<endl
          << "Scale factor to apply to elevation values="<<rawz2elev<<endl
          << "Assuming data are UTM coordinates with zone="<<utmzone<<endl
          << "Name tag on output field data="<<gridname<<endl;
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
  ofstream gmtofs;
  if(write_gmt_output)
  {
    gmtofs.open(gmtoutputfile.c_str(),ios::out);
    if(gmtofs.fail())
    {
      cerr << "Cannot open gmt export file="<<gmtoutputfile<<endl;
      usage();
    }
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
  int NODATAvalue=header.get_int("NODATA_value");
  if(SEISPP_verbose)
  {
      cout << "Header data read from input file:"<<endl;
      cout << header<<endl;
  }
  /* A sanity check*/
  if(NODATAvalue>0)
  {
    cerr << "ArcGISgrd2gcl:  parsed NODATA_value="<<NODATAvalue
       << " is nonsense"<<endl<<"Expected to be a large negative number"<<endl
       << " Assumptions of this program are violated - bug fix required"
       <<endl
       << "Cannot continue"<<endl;
    exit(-1);
  }
  double NODATAtest;
  NODATAtest=(double)NODATAvalue;
  NODATAtest += 1.0;  //Make test 1 foot above null value to avoid float equal test
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
  and online sources confirm this. */
  double lat0,lon0,r0;   //frozen origin here to surface r0 at lower left corner
  /* 23 in this call freezes WGS-84 - not at all elegant*/
  UTMtoLL(23,yll,xll,utmzone.c_str(),lat0,lon0);
  lat0 *= deg2rad;
  lon0 *= deg2rad;
  r0=r0_ellipse(lat0);
  GCLgrid g(nx1,nx2,gridname,lat0,lon0,r0,0.0,cellsize*0.001,cellsize*0.001,0,0);
  GCLscalarfield f(g);
  double easting,northing,elev;
  double lat,lon,r;
  Cartesian_point cp;
  int iv;
  /* Scan order is TV order - start upper left corner, scan left to right,
   * next row down, repeat till bottom.  jj is used for down counting. */
  int jj;
  int number_not_null(0);
  for(j=0,iv=0,jj=nx2-1;j<nx2;++j,--jj)
  {
    for(i=0;i<nx1;++i)
    {
      easting=xll+cellsize*((double)i);
      northing=yll+cellsize*((double)jj);
      UTMtoLL(23,northing,easting,utmzone.c_str(),lat,lon);
      /*Zero masked values */
      if(dvector[iv]<=NODATAtest)
      {
        elev=0;
        f.val[i][jj]=0.0;
      }
      else
      {
        elev=(dvector[iv] - zoffset)*rawz2elev;
        f.val[i][jj]=elev;
        ++number_not_null;
        if(write_gmt_output)
        {
            gmtofs<<lat<<" "<<lon<<" "<<elev<<endl;
        }
      }
      lat*=deg2rad;
      lon*=deg2rad;
      r=r0_ellipse(lat)-elev;
      cp=g.gtoc(lat,lon,r);
      f.x1[i][jj]=cp.x1;
      f.x2[i][jj]=cp.x2;
      f.x3[i][jj]=cp.x3;
      ++iv;
    }
  }
  f.compute_extents();
  if(write_gmt_output) gmtofs.close();
  //f.save(outbase,string("."));
  /* Now apply the mask - this is wrong and needs to be fixed 
   * once I establish the scan order*/
  GCLMask gmask(g);
  for(j=0,iv=0,jj=nx2-1;j<nx2;++j,--jj)
  {
    for(i=0;i<nx1;++i)
    {
      if(dvector[iv]<=NODATAtest)
      {
        gmask.mask_point(i,jj);
      }
      else
      {
          gmask.enable_point(i,jj);
      }
      ++iv;
    }
  }
  GCLMaskedScalarField msf(f,gmask);
  msf.save(outbase);
  if(SEISPP_verbose)
  {
    cout << "Saved results to "<<outbase<<" with "<<number_not_null
      << " set cells of "<<nx1*nx2<<" total grid points"<<endl;
  }
}
