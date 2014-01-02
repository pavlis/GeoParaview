#include "seispp.h"
#include "dbpp.h"
#include "Metadata.h"
#include "gclgrid.h"
#include "SeisppError.h"
using namespace std;
using namespace SEISPP;
/* this program is designed to do a remap operation on GCLgrid objects.
   All use an important trick that must be recognize and is a potential maintenance issue.
   That is, db storage of GCLgrid objects splits the grid and field data storage into 
   to files.   A remap operation really only changes the grid and not the numbers in
   the field.  Hence only the grid is altered and a new entry is created in the
   gclfield table.   To make this doubly confusing it actually would NOT have been 
   necessary to actually make a new file for the modified grid because GCLgrid external
   storage uses lat, lon, r coordinates.   However, I chose to not do that as I thought
   that was really pushing my luck.   
   */

void usage()
{
	cerr << "dbremapgrid db [-g gridname -f fieldname -p pattern -og outgrid -vector -2D]" << endl;
	exit(-1);
}
bool SEISPP::SEISPP_verbose(false);
int main(int argc, char **argv)
{
	int i;
	ios::sync_with_stdio();
	if(argc<2) usage();
        char dir[128],dfile[128]; // these are larger than necessary
        long foff;

	Pf *pf;
	string pffile("dbremapgrid");
	if(pfread(const_cast<char *>(pffile.c_str()),&pf))
	{
		cerr << "pfread failed on file dbremapgrid.pf"<<endl;
		exit(-1);
	}
	string gridname,fieldname,pattern;
	bool vectormode;
	bool twodmode;
	bool use2dpat;
	// These variables are not alterable from command line arguments
	string griddir,newgridname;
	try {
		Metadata control(pf);
		gridname=control.get_string("input_gridname");
		fieldname=control.get_string("input_fieldname");
		pattern=control.get_string("pattern_grid");
		use2dpat=control.get_bool("use_2d_grid_as_pattern");
		vectormode=control.get_bool("vectormode");
		twodmode=control.get_bool("2Dmode");
		griddir=control.get_string("grid_directory");
		newgridname=control.get_string("new_grid_name");
	} catch (SeisppError& mderr)
	{
		mderr.log_error();
		exit(-1);
	}
	string dbname(argv[1]);
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
			newgridname=string(argv[i]);
		}
		else
			cerr << "Unknown argument = "<<argstr<<endl;
			usage();
	}
	try {
	    DatascopeHandle dbh(dbname,false);
	    dbh.lookup(string("gclgdisk"));
	    DatascopeHandle dbhg(dbh);
	    dbhg.lookup(string("gclfield"));
	    Dbptr db=dbh.db;
	    Dbptr dbgrd=dbhg.db;
	    string sstring;
	    BasicGCLgrid *bpat;
	    if(use2dpat) 
	    {
		GCLgrid *pat=new GCLgrid(dbh,pattern);
		bpat=dynamic_cast<BasicGCLgrid*>(pat);
	    }
	    else
	    {	
		GCLgrid3d *pat=new GCLgrid3d(dbh,pattern);
		bpat=dynamic_cast<BasicGCLgrid*>(pat);
	    }
	    if(twodmode)
	    {
		if(vectormode)
		{
			GCLvectorfield g(dbh,gridname,fieldname);
			remap_grid(dynamic_cast<GCLgrid&>(g),*bpat);
			g.name=newgridname;
			string sstring;
			char nvs[20];
			sprintf(nvs,"%d",g.nv);
			sstring = string("gridname =~ /")
					+gridname
					+string("/ && fieldname =~ /")
					+ fieldname 
					+ string("/ && nv==")
					+ string(nvs); 
			Dbptr dbss=dbsubset(dbgrd,const_cast<char *>(sstring.c_str()),NULL);
			long nrec;
			dbquery(dbss,dbRECORD_COUNT,&nrec);
                	if(nrec <= 0)
			{
				cerr <<"Failure forming subset view to get input field parameters"<<endl;
				cerr << "Subset string passed: "<<sstring;
				exit(-1);
			}
                        dbss.record=0;
                        dbgetv(dbss,0,"dir",dir,"dfile",dfile,"foff",&foff,NULL);
			dbaddv(dbgrd,0,"gridname",newgridname.c_str(),"fieldname",fieldname.c_str(),
                                "dimensions",2L,
                                "nv",g.nv,
                                "dir",dir,
                                "dfile",dfile,
                                "foff",foff,NULL);
                        dynamic_cast<GCLgrid&>(g).save(dbh,griddir);

		}
		else
		{
			GCLscalarfield g(dbh,gridname,fieldname);
			remap_grid(dynamic_cast<GCLgrid&>(g),*bpat);
			g.name=newgridname;
			sstring = string("gridname =~ /")
					+gridname
					+string("/ && fieldname =~ /")
					+ fieldname 
					+ string("/");
			Dbptr dbss=dbsubset(dbgrd,const_cast<char *>(sstring.c_str()),NULL);
			long nrec;
			dbquery(dbss,dbRECORD_COUNT,&nrec);
                	if(nrec <= 0)
			{
				cerr <<"Failure forming subset view to get input field parameters"<<endl;
				cerr << "Subset string passed: "<<sstring;
				exit(-1);
			}
                        dbss.record=0;
                        dbgetv(dbss,0,"dir",dir,"dfile",dfile,"foff",&foff,NULL);
			dbaddv(dbgrd,0,"gridname",newgridname.c_str(),"fieldname",fieldname.c_str(),
                                "dimensions",2L,
                                "nv",1L,
                                "dir",dir,
                                "dfile",dfile,
                                "foff",foff,NULL);
                        dynamic_cast<GCLgrid&>(g).save(dbh,griddir);
		}
				
	    }
	    else
	    {
		if(vectormode)
		{
			GCLvectorfield3d g(dbh,gridname,fieldname);
			remap_grid(dynamic_cast<GCLgrid3d&>(g),*bpat);
			g.name=newgridname;
			char nvs[20];
			sprintf(nvs,"%d",g.nv);
			sstring = string("gridname =~ /")
				+ gridname
				+string("/ && fieldname =~ /")
				+ fieldname
				+ string("/ && nv==")
				+ string(nvs); 
			Dbptr dbss=dbsubset(dbgrd,const_cast<char *>(sstring.c_str()),NULL);
			long nrec;
			dbquery(dbss,dbRECORD_COUNT,&nrec);
                	if(nrec <= 0)
			{
				cerr <<"Failure forming subset view to get input field parameters"<<endl;
				cerr << "Subset string passed: "<<sstring;
				exit(-1);
			}
                        dbss.record=0;
                        dbgetv(dbss,0,"dir",dir,"dfile",dfile,"foff",&foff,NULL);
			dbaddv(dbgrd,0,"gridname",newgridname.c_str(),"fieldname",fieldname.c_str(),
                                "dimensions",3L,
                                "nv",g.nv,
                                "dir",dir,
                                "dfile",dfile,
                                "foff",foff,NULL);
                        dynamic_cast<GCLgrid3d&>(g).save(dbh,griddir);
		}
		else
		{
			GCLscalarfield3d g(dbh,gridname,fieldname);
			remap_grid(dynamic_cast<GCLgrid3d&>(g),*bpat);
			g.name=newgridname;
			sstring = string("gridname =~ /")
				+ gridname
				+ string("/ && fieldname =~ /")
				+ fieldname 
				+ string("/");
			Dbptr dbss=dbsubset(dbgrd,const_cast<char *>(sstring.c_str()),NULL);
			long nrec;
			dbquery(dbss,dbRECORD_COUNT,&nrec);
                	if(nrec <= 0)
			{
				cerr <<"Failure forming subset view to get input field parameters"<<endl;
				cerr << "Subset string passed: "<<sstring;
				exit(-1);
			}
                        dbss.record=0;
                        dbgetv(dbss,0,"dir",dir,"dfile",dfile,"foff",&foff,NULL);
			dbaddv(dbgrd,0,"gridname",newgridname.c_str(),"fieldname",fieldname.c_str(),
                                "dimensions",3L,
                                "nv",1L,
                                "dir",dir,
                                "dfile",dfile,
                                "foff",foff,NULL);
                        dynamic_cast<GCLgrid3d&>(g).save(dbh,griddir);
		}
	    }
	}
	catch (SeisppError& serr)
	{
		serr.log_error();
	}
        catch(GCLgridError& gerr)
        {
            cerr << gerr.what()<<endl;
        }
        catch (std::exception& serr)
        {
            cerr << serr.what()<<endl;
        }
	catch (...)
	{
		cerr << "Something threw an exception"<<endl;
	}
}
