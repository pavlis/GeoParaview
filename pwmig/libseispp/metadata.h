#ifndef _METADATA_H_
#define _METADATA_H_
#include <string>
using namespace std;
//
// This follows a method of inherited objects as a way to build
// exception handlers described in books by Stroustrup on C++
//
// This is the base class that is bare bones
//
class Metadata_error
{
public:
	string message;
	virtual void log_error(){elog_log(0,(char *)"Pf error: %s", (char *)message);}
};
class Metadata_get_error : public Metadata_error
{

	string mdtype;
	string name;
	Metadata_get_error(string mtd, string n, string mess)
		{mdtype=mdt; name=n;  message=mess; };
	virtual void log_error()
	 {elog_log(0,(char *)"Error attempting to extract parameter %s of type %s\n%s\n",
			(char *)name, (char *)mdtype,(char *), (char *)message);}  
};
class Metadata_load_error : public Metadata_error
{
	int error_code;
	Metadata_load_error(int ierr){error_code=ierr;};
	virtual void log_error()
	 {elog_log(0,(char *)"pfcompile failed in load_metadata\nReturned code=%d\n",
			error_code);};
};


class Metadata
{
public:
        Metadata();
	Metadata(Metadata&);
	Metadata(const Metadata&);
	Metadata& operator=(const Metadata& );
        ~Metadata();

        double get_double(string) throw(Metadata_get_error);
        int get_int(string) throw(Metadata_get_error);
        string get_string(string) throw(Metadata_get_error);
        bool get_bool(string) throw(Metadata_get_error);
	Tbl *get_list(string) throw(Metadata_get_error);  // Antelope specific
	Arr *get_map(string) throw(Metadata_get_error);  // Antelope specific
	//These put new metadata in
	//These use the same name but depend on overloading
        void put_metadata(string,double);
        void put_metadata(string,int);
        void put_metadata(string,string); // C++ basic string 
        void put_metadata(string,char *);  // for plain C strings
        void put_metadata(string,bool);
	void put_metadata(string,Arr *);  // antelope map
	void put_metadata(string,Tbl *);  // antelope list
	// This compiles a large string with pfcompile to load a full pf
        void load_metadata(string) throw(Metadata_get_error);
private:
        Pf *pf;  // Antelope's pf handle
};
//
//This object is used for selective copy
//
enum MDtype {REAL INTEGER STRING LIST MAP BOOLEAN}
typedef struct Metadata_typedef {
	string tag;
	MDtype mdt;
} Metadata_typedef;

//
// Helpers
//
void copy_selected_metadata(Metadata& mdin, Metadata& mdout, 
		list<Metadata_typedef>mdlist& mdlist);

#endif
