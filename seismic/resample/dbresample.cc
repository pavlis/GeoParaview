#include<string>
#include "pf.h"
#include "decimate.h"
usage()
{
	cerr<<"dbresample db [-pf pffile -sift expression -notrim]"
	exit(-1);
}

int main(int argc, char **argv)
{
	Dbptr db;  // database pointer for input 
	Pf *pf;
	Dbptr db;
	string pffile("resample");
	bool trim=true;
	string sift;
	int number_db_rows;
	Time_Series tracein();
	Time_Series traceout();

	// Standard argument parsing in C++ form

	elog_init(argc,argv);
	if(argc<2)
		usage();
	string dbname(argv[1]);
	for(int i=0;i<argc;++i)
	{
		string test;
		test=string(argv[i]);
		if(test=="-pf")
		{
			++i;
			pffile=string(argv[i]);
		}
		else if(test=="-notrim")
			trim = false;
		}
		else if(test == "-sift")
		{
			++i;
			sift=string(argv[i]);
		else
		{
			cerr<<"Illegal argument = "<<argv[i]<<endl;
			usage();
		}
		if(i>=argc) usage();
	}	
	if(pfread(pfname.c_str(),&pf))
		die(0,static_cast<char *>("Failure reading parameter file "+pffile));
	
	// This object is created from a parameter file and defines what operations are to 
	// to be performed on what sample rates 
	try{ Resample_Definitions rsampdef(pf);}
	catch (seispp_error spe)
	{
		spe.log_error();
		exit(-1);
	}
	dtout = pfget_double(pf,"output_sample_rate");
	if(dbopen(dbname.c_str(),"r",&db))
		die(0,static_cast<char *>("dbopen failed on database "+dbname));

	db = dbform_working_view(db,pf,"dbprocess_definition");
	
	dbquery(db,dbRECORD_COUNT,&number_rows_db);
	if(number_rows_db<=0)
	{
		cerr << "Working view has no rows.  Check dbprocess_definition parameter" << endl;
		exit(-1);
	}
	for(dbrecord=0;db.record<number_rows_db;++db.record)
	{
		tracein = Time_Series(db,pf);  // Note pf has to contain a list of attributes to be cloned in output
		traceout = Resample_Time_Series(tracein,rsampdef,dtout,trim);
		dbsave(tsout,pf);
	}
}
