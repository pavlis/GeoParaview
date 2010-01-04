#include <fstream>
#include "stock.h"
#include "dbpp.h"
#include "Metadata.h"
#include "seispp.h"
#include "gclgrid.h"
using namespace std;
using namespace SEISPP;

void usage()
{
	cerr << "extract_section db gridname fieldname "
		<< "[-pf pffile -o outfile -vectordata]"<<endl
	      << "(Default outfile is stdout)"<<endl;

	exit(-1);
}
double *extract_vertical(GCLscalarfield3d *g,Geographic_point p,double z0, double dz, int nz)
{
	double *result=new double[nz];
	Cartesian_point pc;
	/* p is not a reference so we alter it in the loop below */
	double z;
	double r0=p.r;
	int i;
	for(i=0;i<nz;++i)
	{
		p.r=r0-z0-dz*static_cast<double>(i);
		pc=g->gtoc(p);
		if(g->lookup(pc.x1,pc.x2,pc.x3)==0)
		{
			result[i]=g->interpolate(pc.x1,pc.x2,pc.x3);
		}
		else
		{
			result[i]=0.0;
		}
	}
	return(result);
}
bool SEISPP::SEISPP_verbose(true);
int main(int argc, char **argv)
{
	int i,j;
	ios::sync_with_stdio();
	if(argc<4) usage();
	string dbname(argv[1]);
	string gridname(argv[2]);
	string fieldname(argv[3]);
	string pfname("extract_section");
	bool vectordata(false);
	ofstream outstrm;
	bool out_to_other(false);

	for(i=4;i<argc;++i)
	{
		string argstr=string(argv[i]);
		if(argstr=="-pf")
		{
			++i;
			pfname=string(argv[i]);
		}
		else if(argstr=="-vectordata")
			vectordata=true;
		else if(argstr=="-o")
		{
			++i;
			out_to_other=true;
			string outfile(argv[i]);
			try {
				outstrm.open(outfile.c_str());
			} catch (ios::failure& var)
			{
				cerr << "Open failure on outfile="<<outfile
					<<endl
					<< "System message:  "
					<< var.what() <<endl;
				usage();
			}
		}
		else
		{
			cerr << "Unknown argument = "<<argstr<<endl;
			usage();
		}
	}
	Pf *pf;
	if(pfread(const_cast<char *>(pfname.c_str()),&pf)) 
	{
		cerr << "pfread failed"<<endl;
		usage();
	}
	/* format switch.  For now just ascii formats, but segy planned */
	enum OutputFormat {Dmatrix, GMT};
	char *formatname=pfget_string(pf,"output_format");
	string stmp(formatname);
	OutputFormat odform;
	if(formatname==NULL)
	{
		cerr << "Required parameter output_format is not in "
			<< "parameter file"<<endl;
		usage();
	}
	else if(stmp=="GMT" || stmp=="gmt")
		odform=GMT;
	else
		odform=Dmatrix;
		
	Tbl *pointlist;
	pointlist=pfget_tbl(pf,"section_points");
	if(pointlist==NULL)
	{
		cerr << "pf error on Tbl section_points"
			<<endl
			<<"Check parameter file"<<endl;
		exit(-1);
	}
	list<Geographic_point> section_points;
	for(i=0;i<maxtbl(pointlist);++i)
	{
		char *line;
		double lat,lon;
		line=(char *)gettbl(pointlist,i);
		sscanf(line,"%lf%lf",&lat,&lon);
		Geographic_point p;
		p.lat=rad(lat);
		p.lon=rad(lon);
		p.r=r0_ellipse(p.lat);
		section_points.push_back(p);
	}
	cerr << "Read section path defined by "
		<< section_points.size()
		<< " control points"<<endl;
	try {
		DatascopeHandle dbh(dbname,true);
		dbh.lookup(string("gclgdisk"));
		DatascopeHandle dbhg(dbh);
		dbhg.lookup(string("gclfield"));
		Dbptr db=dbh.db;
		Dbptr dbgrd=dbhg.db;
		
		Metadata control(pf);
		double dx=control.get_double("path_sample_interval");
		double dz=control.get_double("depth_sample_interval");
		double zmax=control.get_double("maximum_depth");
		double zmin=control.get_double("minimum_depth");
		bool output_latlon=control.get_bool("geographic_output");
		int nz=static_cast<int>((zmax-zmin)/dz) + 1;
		if(nz<=0)
		{
			cerr << "Illegal maximum_depth or minimum_depth "
				<< "parameter specified."<<endl
				<< "Check parameter file"<<endl;
			exit(-1);
		}
		/* Now build the full path of control points as equally spaced 
		as possible.  Using the latlon function from antelope. */
		double del,ddel,az,lat1,lon1,lat2,lon2;
		double delkm,r0,dr;
		int n;
		Geographic_point p={0.0,0.0,0.0};
		list<Geographic_point> path;
		list<Geographic_point>::iterator sptr;
		int ipath;
		for(ipath=0,sptr=section_points.begin();
			ipath<(section_points.size() - 1);++ipath)
		{
			if(sptr==section_points.begin())
			{
				lat1=sptr->lat;
				lon1=sptr->lon;
				r0=sptr->r;
			}
			else
			{
				lat1=p.lat;
				lon1=p.lon;
				r0=p.r;
			}
			++sptr;
			lat2=sptr->lat;
			lon2=sptr->lon;
			dr=(sptr->r)-r0;
			dist(lat1,lon1,lat2,lon2,&del,&az);
			delkm=del*r0;
			n=static_cast<int>(delkm/dx) + 1;
			n=SEISPP::nint(delkm/dx) + 1;
			ddel=dx/r0;
			dr=dr*dx/delkm;
			for(i=0;i<n;++i)
			{
				del=ddel*static_cast<double>(i);
				latlon(lat1,lon1,del,az,&lat2,&lon2);
				p.lat=lat2;
				p.lon=lon2;
				p.r=r0+dr*static_cast<double>(i);
				path.push_back(p);
			}
		} 
		
		GCLscalarfield3d *g;
		if(vectordata)
		{
			int component=control.get_int("component_number");
			GCLvectorfield3d *gvec=new GCLvectorfield3d(db,gridname,fieldname);
			g=extract_component(*gvec,component);
			delete gvec;
		}
		else
		{
			g=new GCLscalarfield3d(db,gridname,fieldname);
		}
		/* This will hold the output.  */
		dmatrix section(nz,path.size());
		section.zero();  /* this allows us to ignore edge effects*/
		double *seis;
		for(sptr=path.begin(),j=0;sptr!=path.end();++sptr,++j)
		{
			seis=extract_vertical(g,*sptr,zmin,dz,nz);
			for(i=0;i<nz;++i) 
			{
				section(i,j)=seis[i]; 
			}
			delete [] seis;
		}
		switch (odform)
		{
		case GMT:
			for(j=0,sptr=path.begin();j<section.columns() && sptr!=path.end();++j,++sptr)
			{
				for(i=0;i<nz;++i)
				{
				    if(out_to_other)
				    {
					if(output_latlon)
					    outstrm << deg(sptr->lat)
						<< " "<<deg(sptr->lon)
						<< " "<<zmin+i*dz
						<< " "<<section(i,j)<<endl;
					else
					    outstrm << j*dx
						<< " "<<zmin+i*dz
						<< " "<<section(i,j)<<endl;
				    }
				    else
				    {
					if(output_latlon)
					    cout << deg(sptr->lat)
						<< " "<<deg(sptr->lon)
						<< " "<<zmin+i*dz
						<< " "<<section(i,j)<<endl;
					else
					    cout << j*dx
						<< " "<<zmin+i*dz
						<< " "<<section(i,j)<<endl;
				    }
				}
				if(out_to_other)
					outstrm<<">"<<endl;
				else
					cout << ">"<<endl;
			}
			break;

		case Dmatrix:
		default:
			if(out_to_other)
				outstrm << section;
			else
				cout << section;
		}
		delete g;
	}
	catch (SeisppError& serr)
	{
		serr.log_error();
	}
}
