#include <iostream>
#include <vector>
#include <list>
#include "gclgrid.h"
#include "dbpp.h"
#include "seispp.h"
using namespace std;
using namespace SEISPP;
enum ObjectType{POINTS,LINES,POLYGONS};
void usage()
{
	cerr << "gcl_vtk_converter db gridname [ [-lines|polygons] -2d -noz -clip] < infile "<<endl
		<< "Converted data written to stdout.  Default format points"<<endl;
	exit(-1);
}
bool in_bounding_box(BasicGCLgrid& g,Cartesian_point cp)
{
	if(cp.x1>g.x1high) return(false);
	if(cp.x1<g.x1low) return(false);
	if(cp.x2>g.x2high) return(false);
	if(cp.x2<g.x2low) return(false);
	if(cp.x3>g.x3high) return(false);
	if(cp.x3<g.x3low) return(false);
	return(true);
}

bool SEISPP::SEISPP_verbose(false);
int main(int argc, char **argv)
{
	if(argc<3) usage();
	string dbname(argv[1]);
	string gridname(argv[2]);
	GCLgrid *g2d;
	GCLgrid3d *g3d;
	BasicGCLgrid *grid;
	bool use2dgrid(false);
        ObjectType objtype(POINTS);
	bool noz(false);
	bool clip(false);
	int i;
	for(i=3;i<argc;++i)
	{
		string testarg(argv[i]);
		if(testarg=="-lines")
			objtype=LINES;
                else if(testarg=="-polygons")
                        objtype=POLYGONS;
		else if(testarg=="-2d")
			use2dgrid=true;
		else if(testarg=="-noz")
			noz=true;
		else if(testarg=="-clip")
			clip=true;
		else
			usage();
	}
	try {
		DatascopeHandle dbh(dbname,true);
		dbh.lookup("gclgdisk");
		if(use2dgrid)
		{
			g2d=new GCLgrid(dbh.db,gridname);
			grid=dynamic_cast<BasicGCLgrid *>(g2d);
		}
		else
		{
			g3d=new GCLgrid3d(dbh.db,gridname);
			grid=dynamic_cast<BasicGCLgrid *>(g3d);
		}
		vector<Cartesian_point> cplist;
		Geographic_point gp;
		Cartesian_point cp;
		double lat,lon,depth,radius;
		/* First write common data for either points or lines for file head */
		cout << "# vtk DataFile Version 2.0"<<endl
                  << argv[0] << " conversion using grid="<<grid->name <<endl
                  << "ASCII" <<endl
                  << "DATASET POLYDATA"<<endl;
		typedef list<int> LineLinks;
		LineLinks thissegment;
		list<LineLinks> segments;
		char line[256];
		int i=0;
                switch(objtype)
                {
                case LINES:
                case POLYGONS:
			while(!cin.eof())
			{
				cin.getline(line,256);
				if(cin.eof())
				{
					/*this is an odd logic construct, but
					we make this look like a > line to force
					the last segment to be pushed to the segments
					list.  We make sure the size isn't zero */
					if(thissegment.size()>1)
					{
						line[0]='>';
						line[1]='\0';
					}
				}
				if(line[0]=='#') continue;
				if(line[0]=='>')
				{
					if(thissegment.size()>1)
					{
						segments.push_back(thissegment);
						thissegment.clear();
					}
					else if(thissegment.size()==1)
					{
						/* Handle special case possible when clipping.
						This is if a line wanders in and out of the bounding
						box it is possible to end up with a stray point we
						need to clear like this. */
						thissegment.clear();
					}
				}
				else
				{
					if(noz)
					{
						sscanf(line,"%lf%lf",&lon,&lat);
						depth=0.0;  // inefficient but clearer
					}
					else
					{
						sscanf(line,"%lf%lf%lf",&lon,&lat,&depth);
					}
					gp.lat=rad(lat);
					gp.lon=rad(lon);
					gp.r=r0_ellipse(gp.lat)-depth;
					cp=grid->gtoc(gp);
					if(clip)
					{
						if(in_bounding_box(*grid,cp))
						{
							cplist.push_back(cp);
							thissegment.push_back(i);
							++i;
						}
						else
						{
							/* truncate segment if it
							wanders outside box. 
							drop if only one point */

							if(thissegment.size()>1)
							   segments.push_back(thissegment);
							thissegment.clear();
						}
					}
					else
					{
						cplist.push_back(cp);
						thissegment.push_back(i);
						++i;
					}
				}
			} 
                        break;
                case POINTS:
                default:
		    while(cin.getline(line,256))
		    {
		    	if(noz)
			{
				sscanf(line,"%lf%lf",&lon,&lat);
				depth=0.0;
			}
			else
			{
				sscanf(line,"%lf%lf%lf",&lon,&lat,&depth);
			}
			gp.lat=rad(lat);
			gp.lon=rad(lon);
			gp.r=r0_ellipse(gp.lat)-depth;
			cp=grid->gtoc(gp);
			if(clip)
			{
				if(in_bounding_box(*grid,cp))
					cplist.push_back(cp);
			}
			else
				cplist.push_back(cp);
                    }
		}

		cout << "POINTS " << cplist.size() << " float"<<endl;
		vector<Cartesian_point>::iterator cpptr;
		for(cpptr=cplist.begin();cpptr!=cplist.end();++cpptr)
		{
			cout << cpptr->x1 << " "
				<< cpptr->x2<<" "
				<< cpptr->x3<<" "<<endl;
		}
		int sizepar;
		list<LineLinks>::iterator llptr;
		LineLinks::iterator lineptr;
                switch(objtype)
		{
                case POLYGONS:
                case LINES:
                default:
			/* For this stupid format we have to count the total number
			of items in segments and number of points for the "size" parameter
			of the LINES keyword line. */
			sizepar=0;
			for(llptr=segments.begin();llptr!=segments.end();++llptr)
			{
				for(lineptr=llptr->begin();lineptr!=llptr->end();++lineptr)
					++sizepar;
			}
			sizepar+=segments.size();
                        switch(objtype)
                        {
                            case POLYGONS:
                                cout << "POLYGONS ";
                                break;
                            case LINES:
                            default:
			        cout <<  "LINES "; 
                        }
                        cout << segments.size() << " " <<sizepar<<endl;
			for(llptr=segments.begin();llptr!=segments.end();++llptr)
			{
				cout << llptr->size() << " ";
				for(lineptr=llptr->begin();lineptr!=llptr->end();++lineptr)
				{
					cout << *lineptr<<" ";
				}
				cout << endl;
			}
		}
	} catch (SeisppError serr)
	{
		serr.log_error();
		exit(-1);
	}
	catch (int ierr)
	{
		cerr << "GCLgrid error.  Something threw error code="<<ierr<<endl;
		exit(-1);
	}
}
