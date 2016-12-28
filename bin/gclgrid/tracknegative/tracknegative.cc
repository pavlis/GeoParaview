#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <list>
#include <vector>
#include <algorithm>
#include <set>
#include <map>
#include <tuple>
#include <stdio.h>
#include <float.h>
#include "stock.h"
#include "pf.h"
#include "perf.h"
#include "seispp.h"
#include "Metadata.h"
#include "gclgrid.h"
#include "GCLMasked.h"
#include "vtk_output.h"

using namespace SEISPP;

int vtk_output_GCLgrid(GCLgrid& g, ofstream& out,string name_tag);
int vtk_output_GCLgrid(GCLscalarfield& g, ofstream& out,string name_tag);
int vtk_output_GCLgrid(GCLMaskedScalarField& g, ofstream& out,string name_tag);

//this is used to find the most negative (smallest) points in a vector 
pair<tuple<int,int>,tuple<int,double,int>> find_largest(vector<tuple<int,int,int,double>> points)
{
	int third=0;
	tuple<int,int,int,double> result=points.front();
	for (int i=0;i<points.size();i++)
		if(std::get<3>(points[i]) < std::get<3>(result))
			result=points[i];
	pair<tuple<int,int>,tuple<int,double,int>> newresult(tie(get<0>(result),get<1>(result)),tie(get<2>(result),get<3>(result),third));
	return(newresult);
}
void trackonestep(map<tuple<int,int>,tuple<int,double,int>> &surface, vector<tuple<int,int,int,double>> &vectorsort, int v_limit, int step_limit, int &step)
{
	std::map<tuple<int,int>,tuple<int,double,int>>::iterator itset;
	std::vector<tuple<int,int,int,double>>::iterator itvec;
	int original_size=surface.size();
	//cout<<"starting size is "<<original_size<<endl;
	cout<<"starting vector size is "<<vectorsort.size()<<endl;
	map<tuple<int,int>,tuple<int,double,int>> surfacecopy(surface);
	for (itset=surfacecopy.begin(); itset!=surfacecopy.end(); ++itset)
	{
		if(!get<2>(itset->second))
		{
			for(int x1=std::get<0>(itset->first)-1;x1<=std::get<0>(itset->first)+1;x1++)
				for(int x2=std::get<1>(itset->first)-1;x2<=std::get<1>(itset->first)+1;x2++)
				{
					//cout<<"x1, x2 are "<<x1<<","<<x2<<endl;
					if(x1==std::get<0>(itset->first) && x2==std::get<1>(itset->first))
						continue;
					vector<tuple<int,int,int,double>> points;
					for(int i=0;i<vectorsort.size();i++)
					{
						if(std::get<0>(vectorsort[i])==x1 && std::get<1>(vectorsort[i])==x2 && abs(std::get<2>(vectorsort[i])-std::get<0>(itset->second))<=v_limit)
							points.push_back(vectorsort[i]);
					}
					if(points.size()>0)
					{
						//cout<<"Continue tracking with "<<points.size()<<" points."<<endl;
						surface.insert(find_largest(points));
					}
				}
			get<2>(surface[itset->first])=1;
		}
	}
	cout<<"Surface updated."<<endl;
	//surfacecopy.~map();
	//int con=0;
	/*for (itvec=vectorsort.begin(); itvec!=vectorsort.end(); )
	{
		try{
			//map<tuple<int,int>,tuple<int,double>>::iterator itfind;
			itset=surface.find(tie(get<0>(*itvec),get<1>(*itvec)));
			if(itset != surface.end() && get<0>(itset->second)==get<2>(*itvec))
				vectorsort.erase(itvec);
			else
				++itvec;
		}catch (...)
		{
		    cerr << "Error: cannot remove element!"<<endl;
		}
		//con++;
		//cout<<"Erase No."<<con<<endl;	
	}*/
	cout<<"original set size is "<<original_size<<endl;
	cout<<"new set size is "<<surface.size()<<endl;
	unsigned i = 0;
	while ( i < vectorsort.size() ) 
	{
		itset=surface.find(tie(get<0>(vectorsort[i]),get<1>(vectorsort[i])));
		if(itset != surface.end() && get<0>(itset->second)==get<2>(vectorsort[i]))
		{
			//cout<<"found"<<endl;
			vectorsort.erase( vectorsort.begin() + i );
			//cout<<"done"<<endl;
		} 
		else 
		{
			//if(vectorsort.size()==299545)
			//	cout<<i<<endl;
			++i;
		}
	}
	if(step==step_limit)
	{
		cout<<"Tracking done at step "<<step<<endl;
		return;
	}
	else
		step++;
	if(original_size<surface.size())
	{
		cout<<"new vector size is "<<vectorsort.size()<<endl;
		trackonestep(surface, vectorsort,v_limit,step_limit,step);
	}
}
void usage()
{
	cerr << "tracknegative db|file outfile [-i -g gridname -f fieldname -r "
		<< "-pf pffile] -V" << endl;
	exit(-1);
}
bool SEISPP::SEISPP_verbose(true);
int main(int argc, char **argv)
{
	int i;

        bool dbmode(true);  // default is db input mode
        ios::sync_with_stdio();
	if(argc<3) usage();
	string dbname(argv[1]);
        string infile(argv[1]);  // redundant but easier to do this way
	string outfile(argv[2]);
	string argstr;
	string pffile("tracknegative.pf");
	const string nodef("NOT DEFINED");
	string gridname(nodef);
	string fieldname(nodef);
	bool remap(false);
	for(i=3;i<argc;++i)
	{
		argstr=string(argv[i]);
		if(argstr=="-pf")
		{
			++i;
			if(i>=argc)usage();
			pffile=string(argv[i]);
		}
                else if(argstr=="-i")
                    dbmode=false;
		else if(argstr=="-V")
		{
                    usage();
		}
		else if(argstr=="-g")
		{
			++i;
			if(i>=argc)usage();
			gridname=string(argv[i]);
		}
		else if(argstr=="-f")
		{
			++i;
			if(i>=argc)usage();
			fieldname=string(argv[i]);
		}
		else if(argstr=="-r")
		{
			remap=true;
		}
		else
		{
			usage();
		}
	}
	try {   // used in multiple blocks below 
        	PfStyleMetadata control=pfread(pffile);
                DatascopeHandle dbh;
                if(dbmode)
                    dbh=DatascopeHandle(dbname,true);

		if(gridname==nodef)
			gridname=control.get_string("gridname");
		if(fieldname==nodef)
			fieldname=control.get_string("fieldname");
		string fieldtype=control.get_string("fieldtype");
		cout << "Converting field = "<<fieldname
			<< " associated with gridname="<<gridname<<endl;
		cout << "Assuming grid type is "<<fieldtype<<endl;
                int nv_expected(5);  // defaulted for pwmig
                int v_target(2);
                if(fieldtype=="vector3d")
                {
                    nv_expected=control.get_int("nv_expected");
                    v_target=control.get_int("v_target"); //vector component to be processed
                }
                string scalars_tag=control.get_string("scalars_name_tag");
                double sdepth=control.get_double("starting_depth");
                double edepth=control.get_double("ending_depth");
                int maxn=control.get_int("max_number_of_surface");
		double upper_limit=control.get_double("amplitude_upper_limit");
		double lower_limit=control.get_double("amplitude_lower_limit");
                int minpoints=control.get_int("minimum_number_of_points_in_a_surface");
		int v_limit=control.get_int("vertical_limit_around_points");
		int step_limit=control.get_int("maximum_number_of_steps_for_tracking");
		/* Slightly odd logic here, but this allows remap off
		to be the default.  pf switch is ignored this way if
		the -r flag was used */
		if(!remap) remap=control.get_bool("remap_grid");
		BasicGCLgrid *rgptr;
                rgptr=NULL;
		if(remap)
		{
			cout << "Remapping enabled"<<endl;
                        /* We create a small temporary GCLgrid to allow
                           use of remap_grid procedure.   */
                        double lat0,lon0,r0,azm;
                        lat0=control.get_double("latitude_origin");
                        lon0=control.get_double("longitude_origin");
                        r0=control.get_double("radius_origin");
                        azm=control.get_double("azimuth_y_axis");
                        cout << "Origin lat,lon,r="
                            <<lat0<<", "<<lon0<<", "<<r0<<endl;
                        cout << "Coordinate system y axis rotation="
                            <<azm<<endl;
                        lat0=rad(lat0);
                        lon0=rad(lon0);
                        azm=rad(azm);
                        GCLgrid *gtmp=new GCLgrid(2,2,string("reference"),
                                lat0,lon0,r0,azm,
                                1.0,1.0,0,0);
                        rgptr=dynamic_cast<BasicGCLgrid*>(gtmp);

		}

                GCLscalarfield3d field;
		if(fieldtype=="scalar3d") 
		{
                        if(dbmode)
			    field=GCLscalarfield3d(dbh,gridname,fieldname);
                        else
                            field=GCLscalarfield3d(infile);
			if(remap)
			{
				// Used to make this optional.  force
				//if(field!=(*rgptr))
				remap_grid(dynamic_cast<GCLgrid3d&>(field),
						*rgptr);
			}
		}
		else if(fieldtype=="vector3d")
		{
                        GCLvectorfield3d vfield;
                        if(dbmode)
			    vfield=GCLvectorfield3d(dbh,gridname,
                                    fieldname,nv_expected);
                        else
                            vfield=GCLvectorfield3d(infile);
			if(remap)
			{
				//if(vfield!=(*rgptr))
					remap_grid(dynamic_cast<GCLgrid3d&>(vfield),
						*rgptr);
			}
			GCLscalarfield3d *sfptr;
			sfptr = extract_component(vfield,v_target);
			field=*sfptr;
			delete sfptr;
		}
		else
		{
			cerr << "Unsupported fieldtype = "<<fieldtype<<endl
				<< "Exiting with no output\n" << endl;
		}
		struct sortingclass {
			bool operator() (tuple<int,int,int,double> i,tuple<int,int,int,double> j) { return (std::get<3>(i)<std::get<3>(j));}
		} custom_sort;
		vector<tuple<int,int,int,double>> vectorsort;
		for(int i=0;i<field.n1;i++)
			for(int j=0;j<field.n2;j++)
				for(int k=0;k<field.n3;k++)
				{
					if(field.depth(i,j,k)>sdepth && field.depth(i,j,k)<edepth && field.val[i][j][k]<upper_limit && field.val[i][j][k]>lower_limit)
					{
						tuple<int,int,int,double> gridpoint;
						std::get<0>(gridpoint)=i;
						std::get<1>(gridpoint)=j;
						std::get<2>(gridpoint)=k;
						std::get<3>(gridpoint)=field.val[i][j][k];
						vectorsort.push_back(gridpoint);
					}
				}
		std::sort (vectorsort.begin(), vectorsort.end(), custom_sort);
		
		//vector<GCLMaskedScalarField> horizons;
		for(int ii=0;ii<maxn && vectorsort.size()>0;ii++)
		{
			tuple<int,int,int,double> seed=vectorsort.front();
			map<tuple<int,int>,tuple<int,double,int>> surface;
			int third=0;
			surface.insert(pair<tuple<int,int>,tuple<int,double,int>>(tie(get<0>(seed),get<1>(seed)),tie(get<2>(seed),get<3>(seed),third)));
			int step=0;
			trackonestep(surface,vectorsort,v_limit,step_limit,step);
			cout<<"Done with surface tracking "<<ii<<" with "<<surface.size()<<"points."<<endl;
			if(surface.size()<minpoints)
			{
				ii--;
				continue;
			}
			
			GCLgrid grid2d(field.n1, field.n2, 
				field.name,
				field.lat0, field.lon0, field.r0,
				field.azimuth_y, field.dx1_nom, field.dx2_nom, 
				field.i0, field.j0);
			
			GCLMaskedGrid tempgrid(grid2d);
			GCLMaskedScalarField newhorizon(tempgrid);
			
			for(int i=0;i<newhorizon.n1;i++)
				for(int j=0;j<newhorizon.n2;j++)
					newhorizon.mask_point(i, j);
			
			for(std::map<tuple<int,int>,tuple<int,double,int>>::iterator itset=surface.begin(); itset!=surface.end(); ++itset)
			{
				newhorizon.enable_point(get<0>(itset->first), get<1>(itset->first));
				newhorizon.x1[get<0>(itset->first)][get<1>(itset->first)]=field.x1[get<0>(itset->first)][get<1>(itset->first)][get<0>(itset->second)];
				newhorizon.x2[get<0>(itset->first)][get<1>(itset->first)]=field.x2[get<0>(itset->first)][get<1>(itset->first)][get<0>(itset->second)];
				newhorizon.x3[get<0>(itset->first)][get<1>(itset->first)]=field.x3[get<0>(itset->first)][get<1>(itset->first)][get<0>(itset->second)];
				//newhorizon.val[get<0>(itset->first)][get<1>(itset->first)]=get<1>(itset->second);
				newhorizon.val[get<0>(itset->first)][get<1>(itset->first)]=field.depth(get<0>(itset->first),get<1>(itset->first),get<0>(itset->second));
			}
			
			//horizons.push_back(newhorizon);
			cout<<"Done with horizon "<<ii<<" with "<<surface.size()<<"points."<<endl;

			stringstream ss;
			ss << outfile <<"_"<<ii<<".vtk";
		        ofstream out;
		        out.open(ss.str().c_str(),ios::out);
			int npoly=vtk_output_GCLgrid(newhorizon,out,scalars_tag);
		        out.close();
		}
		
		//for(int ii=0;ii<horizons.size();ii++)
		//{
                //}

	}
	catch (MetadataGetError& mderr)
	{
		mderr.log_error();
		exit(-1);
	}
        catch(GCLgridError& gerr)
        {
            cerr << "GCLgridError thrown"<<endl<<gerr.what();
        }
        catch(std::exception& sexcp)
        {
            cerr << sexcp.what();
        }
        catch(SeisppError& serr)
        {
            serr.log_error();
        }
	catch (...)
	{
	    cerr << "Something threw an unhandled exception"<<endl;
	}
}
