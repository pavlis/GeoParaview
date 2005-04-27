#include "gclgrid.h"

void usage()
{
	cerr << "dbremapgrid db [-g gridname -f fieldname -p pattern -og outgrid -of outfield -vector -2D]" << endl;
	exit(-1);
}
int main(int argc, char **argv)
{
	int i;
	ios::sync_with_stdio();
	if(argc<2) usage();

	Pf *pf;
	char *pffile="dbremapgrid";
	if(pfread(pffile,&pf)
	{
		cerr << "pfread failed on file dbremapgrid.pf"<<endl;
		exit(-1);
	}
	string gridname,fieldname,pattern;
	string outgrd,outfld;
	bool vectormode;
	bool twodmode;
	// These variables are not alterable from command line arguments
	string griddir,fielddir,newgridname,fdfileout;
	try {
		Metadata control(pf);
		gridname=control.get_string("input_gridname");
		fieldname=control.get_string("input_fieldname");
		pattern=control.get_string("pattern_grid");
		outgrd=control.get_string("output_gridname");
		outfld=control.get_string("output_fieldname");
		vectormode=control.get_bool("vectormode");
		twodmode=control.get_bool("2Dmode");
		griddir=control.get_string("grid_directory");
		fielddir=control.get_string("field_directory");
		newgridname=control.get_string("new_grid_name");
		fdfileout=control.get_string("output_field_dfile");
	} catch (Metadata_error mderr)
	{
		mderr.log_error();
		exit(-1);
	}
	string dbname(argc[1]);
	string argstr;
	// Arg list can override pf contents in this logic
	for(i=2;i<argc;++i)
	{
		argstr=string(argv[i]);
		if(argstr=="-vector")
			vectormode=true;
		else if(argstr=="-2D")
			twodmode=true;
		else if(argstr=="-g")
		{
			++i;
			gridname=string(argv[i]);
		}
		else if(argstr=="-f")
		{
			++i;
			fieldname=string(argv[i]);
		}
		else if(argstr=="-p")
		{
			++i;
			pattern=string(argv[i]);
		}
		else if(argstr=="-og")
		{
			++i;
			outgrd=string(argv[i]);
		}
		else if(argstr=="-of")
		{
			++i;
			outfld=string(argv[i]);
		}
		else
			cerr << "Unknown argument = "<<argstr<<endl;
			usage();
	}
	Dbptr db;
	if(dbopen((const_cast<char *>(dbname.c_str()),"r+",&db)
			== dbINVALID)  
	{
		cerr << "Open failure on database "<<dbname<<endl;
		usage();
	}
	/* The following block could be done more elegantly with a template
	// In all cases the new grid is saved, but only the field variables
	// are not changed.  As a result we only need a new entry in the field
	// table to connect these values to the alternate grid geometry.
	*/
	Dbptr dbgrd;
	// Don't need to trap an error here as if this lookup fails 
	// the contructors below will also fail.
	dbgrd=dblookup(db,0,(char *)"gclfield",0,0);
	try {
	    if(twodmode)
	    {
		GCLgrid pat(db,pattern);
		if(vectormode)
		{
			GCLvectorfield g(db,gridname,fieldname);
			GCLgrid newg=remap_grid(dynamic_cast<GCLgrid&>(g),
					dynamic_cast<BasicGCLgrid&>(pat));
			newg.name=newgridname;
			string sstring;
			char nvs[20];
			sprintf(nvs,"%d",g.nv);
			sstring = "gridname =~ /"+gridname+"/ && fieldname =~ /"
					+ fieldname "/ && nv=="
					+ string(nvs) + "/"; 
			Dbptr dbss=dbsubset(dbgrd,sstring.c_str(),0);
			int nrec;
			dbquery(dbss,dbRECORD_COUNT,&nrec);
                	if(nrec <= 0)
			{
				cerr <<"Failure forming subset view to get input field parameters"<<endl;
				cerr << "Subset string passed: "<<sstring;
				exit(-1);
			}
			// This uses the feature that dbget with 0 stuffs this record
			// into the scratch record
			dbget(dbss,0);
			dbgrd.record=dbSCRATCH;
			dbputv(dbgrd,0,"gridname",newgridname.c_str(),0);
			dbput(dbgrd,0);
			newg.dbsave(db,griddir);

		}
		else
		{
			GCLscalarfield g(db,gridname,fieldname);
			GCLgrid newg=remap_grid(dynamic_cast<GCLgrid&>(g),
					dynamic_cast<BasicGCLgrid&>(pat));
			newg.name=newgridname;
			sstring = "gridname =~ /"+gridname+"/ && fieldname =~ /"
					+ fieldname "/";
			Dbptr dbss=dbsubset(dbgrd,sstring.c_str(),0);
			int nrec;
			dbquery(dbss,dbRECORD_COUNT,&nrec);
                	if(nrec <= 0)
			{
				cerr <<"Failure forming subset view to get input field parameters"<<endl;
				cerr << "Subset string passed: "<<sstring;
				exit(-1);
			}
			// This uses the feature that dbget with 0 stuffs this record
			// into the scratch record
			dbget(dbss,0);
			dbgrd.record=dbSCRATCH;
			dbputv(dbgrd,0,"gridname",newgridname.c_str(),0);
			dbput(dbgrd,0);
			newg.dbsave(db,griddir);
		}
				
	    }
	    else
	    {
		GCLgrid3d pat(db,pattern);
		if(vectormode)
		{
			GCLvectorfield3d g(db,gridname,fieldname);
			GCLgrid newg=remap_grid(dynamic_cast<GCLgrid&>(g),
					dynamic_cast<BasicGCLgrid&>(pat));
			newg.name=newgridname;
			sstring = "gridname =~ /"+gridname+"/ && fieldname =~ /"
					+ fieldname "/ && nv=="
					+ string(nvs) + "/"; 
			Dbptr dbss=dbsubset(dbgrd,sstring.c_str(),0);
			int nrec;
			dbquery(dbss,dbRECORD_COUNT,&nrec);
                	if(nrec <= 0)
			{
				cerr <<"Failure forming subset view to get input field parameters"<<endl;
				cerr << "Subset string passed: "<<sstring;
				exit(-1);
			}
			// This uses the feature that dbget with 0 stuffs this record
			// into the scratch record
			dbget(dbss,0);
			dbgrd.record=dbSCRATCH;
			dbputv(dbgrd,0,"gridname",newgridname.c_str(),0);
			dbput(dbgrd,0);
			newg.dbsave(db,griddir);
		}
		else
		{
			GCLscalarfield3d g(db,gridname,fieldname);
			GCLgrid newg=remap_grid(dynamic_cast<GCLgrid&>(g),
					dynamic_cast<BasicGCLgrid&>(pat));
			newg.name=newgridname;
			sstring = "gridname =~ /"+gridname+"/ && fieldname =~ /"
					+ fieldname "/";
			Dbptr dbss=dbsubset(dbgrd,sstring.c_str(),0);
			int nrec;
			dbquery(dbss,dbRECORD_COUNT,&nrec);
                	if(nrec <= 0)
			{
				cerr <<"Failure forming subset view to get input field parameters"<<endl;
				cerr << "Subset string passed: "<<sstring;
				exit(-1);
			}
			// This uses the feature that dbget with 0 stuffs this record
			// into the scratch record
			dbget(dbss,0);
			dbgrd.record=dbSCRATCH;
			dbputv(dbgrd,0,"gridname",newgridname.c_str(),0);
			dbput(dbgrd,0);
			newg.dbsave(db,griddir);
		}
	    }
	}
	catch (...)
	{
		cerr << "Something threw an exception"<<endl;
	}
}
