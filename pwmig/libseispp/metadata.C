#include <iostream.h>
#include <string>
using namespace std;
#include "stock.h"
#include "pf.h"
#include "metadata.h"
//
// needed until antelope 4.5
//
Pf *pfdup(Pf*);
//
// This file is used to create a shared object library.  This 
// approach may not be universally applicable, but it is known to work
// in solaris.  When a library is first openned the _init() function is
// executed.  Here we utilized this to initialized a global parameter 
// space for the library used as a default template for metadata
//
static Pf *Metadata_defaults_pf;
void _init()
{
	char *mddname=(char *)"metadata_defaults";
	char *name_from_env;
	char *name;
	int ierr;
	// Get name from the environment
	// revert to default pf name if not defined 
	name_from_env = getenv("METADATA_DEFAULTS");
	if(name_from_env==NULL)
		name=mddname;
	else
		name=name_from_env;

	ierr = pfread(name,&Metadata_defaults_pf);
	if(ierr)
	{
		cerr << "Metadata_default initialization failed\npfread failed for" << name <<".pf\n";
		exit(-1);
	}
}
//
// The default constructor uses the defaults as a template which
// allows new parameters to replace defaults in place.
//
Metadata::Metadata()
{
	pf=pfdup(Metadata_defaults_pf);
};
Metadata::~Metadata()
{
	pffree(pf);
};
Metadata::Metadata(const Metadata& md)
{
	char *pfstr;  
	int ierr;

	// This approach uses the default as a template then compiles
	// the contents of md.pf on top of the defaults.  This is the
	// way in the current pf library to do this.  An ugly solution
	// in terms of memory because of the conversion to and from a string
        pf=pfdup(Metadata_defaults_pf);
	pfstr = pf2string(md.pf);
	ierr = pfcompile(pfstr,&pf);
	// Just post this error and continue.  It should never happen, but
	// at least we handle the return from pfcompile correctly.
	if(ierr) cout << "Metadata copy constructor failure in pfcompile\n";
	free(pfstr);
}
Metadata& Metadata::operator=(const Metadata& md)
{
	if(this != &md)
	{
		pf = pfdup(md.pf);
	}
	return *this;
}
//
// These functions get and convert values
//
double Metadata::get_double(string s) throw(int)
{
	void *result;
	double val;
	int pftype_return;
	char *char_val;

	pftype_return = pfget(pf,(char *)s.c_str(),&result);
	if(pftype_return != PFSTRING) throw pftype_return;
	char_val=pfget_string(pf,(char *)s.c_str());
	val = atof(char_val);
	return(val);
}
int Metadata::get_int(string s)throw(int)
{
	void *result;
	int val;
	int pftype_return;
	char *char_val;
	pftype_return = pfget(pf,(char *)s.c_str(),&result);
	if(pftype_return != PFSTRING) throw pftype_return;
	char_val=pfget_string(pf,(char *)s.c_str());
	val = atoi(char_val);
	return(val);
}
string Metadata::get_string(string s)throw(int)
{
	void *result;
	string val;
	int pftype_return;
	char *char_val;
	pftype_return = pfget(pf,(char *)s.c_str(),&result);
	if(pftype_return != PFSTRING) throw pftype_return;
	char_val=pfget_string(pf,(char *)s.c_str());
	val = char_val; //= is overloaded for this case so is a simple assign
	return(val);
}
bool Metadata::get_bool(string s)throw(int)
{
	void *result;
	int val;
	int pftype_return;
	//
	// Boolean will default to false if not found in the pf
	// this is consistent with documentation and is rational behaviour
	//
	val = pfget_boolean(pf,(char *)s.c_str());
	if(val) 
		return(true);
	else 
		return(false);
}
Tbl *Metadata::get_list(string s)throw(int)
{
	void *result;
	int val;
	int pftype_return;
	Tbl *t;
	pftype_return = pfget(pf,(char *)s.c_str(),&result);
	if(pftype_return != PFTBL) throw pftype_return;
	t = pfget_tbl(pf,(char *)s.c_str());
	return(t);
}
Arr *Metadata::get_map(string s)throw(int)
{
	void *result;
	int val;
	int pftype_return;
	Arr *a;
	pftype_return = pfget(pf,(char *)s.c_str(),&result);
	if(pftype_return != PFARR) throw pftype_return;
	a = pfget_arr(pf,(char *)s.c_str());
	return(a);
}

//
// Functions to put things into metadata object
//
void Metadata::put_metadata(string name, double val)
{
	pfput_double(pf,(char *)name.c_str(),val);
}
void Metadata::put_metadata(string name, int val)
{
	pfput_int(pf,(char *)name.c_str(),val);
}
void Metadata::put_metadata(string name, string val)
{
	pfput_string(pf,(char *)name.c_str(),(char *)val.c_str());
}
// for C style strings, we should not depend on the compiler
void Metadata::put_metadata(string name, char *val)
{
	pfput_string(pf,(char *)name.c_str(),val);
}
void Metadata::put_metadata(string name, bool val)
{
	if(val)
		pfput_boolean(pf,(char *)name.c_str(),1);
	else
		pfput_boolean(pf,(char *)name.c_str(),0);
}
void Metadata::load_metadata(string mdin) throw(int)
{
	int ierr;
	// We might think this was needed:  if(pf!=NULL) pffree(pf);
	// it is not because pfcompile checks for *pf==NULL and
	// assumes it is valid and to be updated if nonzero.
	// This is a handy way to deal with defaults
	ierr = pfcompile((char *)mdin.c_str(),&pf);
	if(ierr!=0) throw ierr;
}
//
// Include an explicit error handler with this library
// Handles errors thrown by get routines above
//
void get_metadata_error_handler(string name, int ecode)
{
	string base_message="get_metadata_by_name:  type conflict for parameter name=";
	switch (ecode)
	{
	case PFFILE:
		cerr << base_message << name << "\nParsed as a pffile\n";
		break;
	case PFTBL:
		cerr << base_message << name << "\nParsed as pftbl, not currently supported\n";
		break;
	case PFARR:
		cerr << base_message << name << "\nParsed as pfarr, not currently supported\n";
		break;
	case PFINVALID:
		cerr << base_message << name << "\nParse failed and returned pfinvalid\nNo Metadata loaded\n";
	default:
		cerr << base_message << name << "\nUnrecognized error code\nFatal:  exiting\n";
		exit(-1);
	}
}
void load_metadata_error_handler(int ecode)
{
	elog_notify(0,(char *)"pfcompile failed in load_metadata\nWill continue but metadata may be incomplete\n");
}

