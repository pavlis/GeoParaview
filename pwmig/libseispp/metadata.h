#ifndef _METADATA_H_
#define _METADATA_H_
class Metadata
{
public:
        Metadata();
	Metadata(Metadata&);
	Metadata operator=(const Metadata& );
        ~Metadata();
        //Each of these functions returns standard type metadata
        //The list could be expanded to support generic list
        //and map containers (Antelope &Tbl and &Arr in a pf file)
        //but this not necessary for initial development by author
        double get_metadata_by_name(string) throw(int);
        int get_metadata_by_name(string) throw(int);
        string get_metadata_by_name(string) throw(int);
        bool get_metadata_by_name(string) throw(int);
	Tbl *get_metadata_by_name(string) throw(int);  // Antelope specific
	//These put new metadata in
        void put_metadata_by_name(string,double);
        void put_metadata_by_name(string,int);
        void put_metadata_by_name(string,string);
        void put_metadata_by_name(string,char *);
        void put_metadata_by_name(string,bool);
	// This compiles a large string with pfcompile to load a full pf
        void load_metadata(string) throw(int);
private:
        Pf *pf;  // Antelope's pf handle
};
#endif
