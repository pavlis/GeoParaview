#ifndef _METADATA_H_
#define _METADATA_H_
class Metadata
{
public:
        Metadata();
	Metadata(Metadata&);
	Metadata(const Metadata&);
	Metadata& operator=(const Metadata& );
        ~Metadata();

        double get_double(string) throw(int);
        int get_int(string) throw(int);
        string get_string(string) throw(int);
        bool get_bool(string) throw(int);
	Tbl *get_list(string) throw(int);  // Antelope specific
	Arr *get_map(string) throw(int);  // Antelope specific
	//These put new metadata in
        void put_metadata(string,double);
        void put_metadata(string,int);
        void put_metadata(string,string);
        void put_metadata(string,char *);
        void put_metadata(string,bool);
	// This compiles a large string with pfcompile to load a full pf
        void load_metadata(string) throw(int);
private:
        Pf *pf;  // Antelope's pf handle
};
#endif
