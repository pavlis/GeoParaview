#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include "stock.h"
#include "elog.h"
#include "pf.h"
#include "metadata.h"
#include "gclgrid.h"
#include "vtk_output.h"

using namespace SEISPP;

void usage()
{
	cerr << "gclfield2vtk db outfile [-pf pffile]" << endl;
	exit(-1);
}
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

		string gridname=control.get_string("gridname");
		string fieldname=control.get_string("fieldname");
		string fieldtype=control.get_string("fieldtype");
		bool SaveAsVectorField;
		int VectorComponent;
		if(fieldtype=="scalar3d")
		{
			GCLscalarfield3d field(db,gridname,fieldname);
			output_gcl3d_to_vtksg(field,
				const_cast<string>(outfile));
		}
		else if(fieldtype=="vector3d")
		{
			SaveAsVectorField=control.get_bool("save_as_vector_field");
			if(SaveAsVectorField)
			{
				cerr << "Vector field output not yet supported but planned"
					<< endl;
			}
			GCLvectorfield3d vfield(db,gridname,fieldname);
			GCLscalarfield3d *sfptr;
			sfptr = extract_component(vfield,0);
			string out0=outfile+"_0";
			output_gcl3d_to_vtksg(*sfptr,
				const_cast<string>(out0));
			delete sfptr;
			sfptr = extract_component(vfield,1);
			string out1=outfile+"_1";
			output_gcl3d_to_vtksg(*sfptr,
				const_cast<string>(out1));
			delete sfptr;
			sfptr = extract_component(vfield,2);
			string out2=outfile+"_2";
			output_gcl3d_to_vtksg(*sfptr,
				const_cast<string>(out2));
			delete sfptr;
		}
		else
		{
			cerr << "Unsupported fieldtype = "<<fieldtype<<endl
				<< "Exiting with no output\n" << endl;
		}

	}
	catch (Metadata_get_error mderr)
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
