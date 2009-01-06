#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdio.h>
#include "stock.h"
#include "elog.h"
#include "pf.h"
#include "seispp.h"
#include "Metadata.h"
#include "gclgrid.h"
#include "vtk_output.h"

using namespace SEISPP;

int vtk_output_GCLgrid(GCLgrid& g, string ofile);
/* Removes mean for constant x3 slices.  Important for
proper display of tomography models showing absolute velocities */
void remove_mean_x3(GCLscalarfield3d& f)
{
	double mean;
	int i,j,k;
	double nval=static_cast<double>((f.n1)*(f.n2));
	cout << "Mean removal values"<<endl
		<< "grid-index     mean"<<endl;
	for(k=0;k<f.n3;++k)
	{
		mean=0.0;
		for(i=0;i<f.n1;++i)
			for(j=0;j<f.n2;++j)
				mean += f.val[i][j][k];
		mean /= nval;
		cout << k << "      "<<mean<<endl;
		for(i=0;i<f.n1;++i)
			for(j=0;j<f.n2;++j)
				f.val[i][j][k] = f.val[i][j][k]-mean;
	}
}
		

void usage()
{
	cerr << "gclfield2vtk db outfile [-g gridname -f fieldname -r remapgridname -pf pffile] -V" << endl;
	exit(-1);
}
bool SEISPP::SEISPP_verbose(true);
int main(int argc, char **argv)
{
	int i;

        ios::sync_with_stdio();
	elog_init(argc,argv);
	if(argc<3) usage();
	string dbname(argv[1]);
	string outfile(argv[2]);
	string argstr;
	string pffile("gclfield2vtk");
	const string nodef("NOT DEFINED");
	string gridname(nodef);
	string fieldname(nodef);
	string remapgridname(nodef);
	for(i=3;i<argc;++i)
	{
		argstr=string(argv[i]);
		if(argstr=="-pf")
		{
			++i;
			if(i>=argc)usage();
			pffile=string(argv[i]);
		}
		else if(argstr=="-V")
		{
		    cbanner("1.0",
			"gclfield2vtk db outfile [-pf pffile]",
			"Gary L. Pavlis",
			"Indiana University",
			"pavlis@indiana.edu");
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
			++i;
			if(i>=argc)usage();
			remapgridname=string(argv[i]);
		}
		else
		{
			usage();
		}
	}
        Pf *pf;
        if(pfread(const_cast<char *>(pffile.c_str()),&pf))
                elog_die(1,"pfread error\n");

	try {
        	Metadata control(pf);
		Dbptr db;
		if(dbopen(const_cast<char *>(dbname.c_str()),
				"r+",&db)==dbINVALID)
			elog_die(1,"dbopen failure\n");

		if(gridname==nodef)
			gridname=control.get_string("gridname");
		if(fieldname==nodef)
			fieldname=control.get_string("fieldname");
		string fieldtype=control.get_string("fieldtype");
		cout << "Converting field = "<<fieldname
			<< " associated with gridname="<<gridname<<endl;
		cout << "Assuming grid type is "<<fieldtype<<endl;
		bool remap=control.get_bool("remap_grid");
		BasicGCLgrid *rgptr;
		if(remap)
		{
			if(remapgridname==nodef)
				remapgridname=control.get_string("remapgridname");
			cout << "Remapping emabled.  "
				<<"Will use coordinate system of grid "
				<< remapgridname<<endl;
			int griddim=control.get_int("remapgrid_number_dimensions");
			if(griddim==3)
			{
				GCLgrid3d *gtmp=new GCLgrid3d(db,
					remapgridname);
				rgptr=dynamic_cast<BasicGCLgrid*>(gtmp);
			}
			else if(griddim==2)
			{
				GCLgrid *gtmp=new GCLgrid(db,
					remapgridname);
				rgptr=dynamic_cast<BasicGCLgrid*>(gtmp);
			}
			else
			{
				cerr << "Illegal remap grid dimension="
					<<griddim<<endl
					<<"Fix parameter file"<<endl;
			}
		}
			
		bool SaveAsVectorField;
		int VectorComponent;
		bool rmeanx3=control.get_bool("remove_mean_x3_slices");
		if(fieldtype=="scalar3d")
		{
			GCLscalarfield3d field(db,gridname,fieldname);
			if(remap)
			{
				if(field!=(*rgptr))
					remap_grid(dynamic_cast<GCLgrid3d&>(field),
						*rgptr);
			}
			if(rmeanx3) remove_mean_x3(field);
			output_gcl3d_to_vtksg(field,
				const_cast<string>(outfile),false,true);
		}
		else if(fieldtype=="vector3d")
		{
			SaveAsVectorField=control.get_bool("save_as_vector_field");
			if(SaveAsVectorField)
			{
				cerr << "Vector field output not yet supported but planned"
					<< endl;
			}
			if(rmeanx3)
			{
				cerr << "remove_mean_x3_slices set true:  "
					<< "ignored for vector field="
					<< fieldname<<endl;
			}
			GCLvectorfield3d vfield(db,gridname,fieldname);
			GCLscalarfield3d *sfptr;
			if(remap)
			{
				if(vfield!=(*rgptr))
					remap_grid(dynamic_cast<GCLgrid3d&>(vfield),
						*rgptr);
			}
			for(i=0;i<vfield.nv;++i)
			{
				sfptr = extract_component(vfield,i);
				stringstream ss;
				ss << outfile <<"_"<<i<<".vtk";
				output_gcl3d_to_vtksg(*sfptr,
					ss.str(),false,false);
				delete sfptr;
			}
		}
		else if(fieldtype=="grid2d")
		{
			int npoly;
			GCLgrid g(db,gridname);
			if(remap)
			{
				if(g!=(*rgptr))
					remap_grid(g,*rgptr);
			}
			npoly=vtk_output_GCLgrid(g,outfile);
			cout << "Wrote "<<npoly<<" polygons to output file"<<endl;
		}
		else if(fieldtype=="scalar2d")
		{
			int npoly;
			GCLscalarfield field(db,gridname,fieldname);
			if(remap)
			{
				if(field!=(*rgptr))
					remap_grid(dynamic_cast<GCLgrid&>(field),
						*rgptr);
			}
			vtk_output_GCLgrid(dynamic_cast<GCLgrid&>(field),outfile);
			cout << "Wrote "<<npoly<<" polygons to output file"<<endl;
		}	
		else
		{
			cerr << "Unsupported fieldtype = "<<fieldtype<<endl
				<< "Exiting with no output\n" << endl;
		}

	}
	catch (MetadataGetError mderr)
	{
		mderr.log_error();
		exit(-1);
	}
	catch(int ierr)
	{
		elog_die(1,"GCLgrid library error\n");
	}
	catch (...)
	{
		elog_die(1,"Something threw and unhandled excepton\n");
	}
}
