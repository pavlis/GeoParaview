#include "sh2db.h"
#include <string>
#include <sstream>
#include <map>
class SHHdrAttribute
{
public:
	string key;
	string shtype;
	int qfnum;
	string shname;
	string table;
	string outname;
	MDtype mdt;
	SHHdrAttribute(string line);
};
SHHdrAttribute::SHHdrAttribute(string line)
{
	istringstream ss(line);
	ss >> key;
	ss >> shtype;
	ss >> qfnum;
	ss >> shname;
	ss >> table;
	if(table=="NULL") table="";
	ss >> outname;
	string mdtname;
	ss >> mdtname;
	if(mdtname=="int") 
		mdt=MDint;
	else if(mdtname=="real")
		mdt=MDreal;
	else if(mdtname=="string")
		mdt=MDstring;
	else
		md=MDinvalid;
}
map<string,SHHdrAttribute> load_shhdrmap(Pf *pf)
{
	map<string,SHHdrAttribute>  result;
	Tbl *t=pfget_tbl(pf,"SHHdrAttributes");
	if(t==NULL)
	{
		cerr << "Parameter file is missing required Tbl SHHdrAttributes"
			<<endl;
		exit(-1);
	}
	int i;
	for(i=0;i<maxtbl(t,i);++i)
	{
		char *s;
		s=(char *) gettbl(t,i);
		SHHdrAttribute sattr(s);
		result[sattr.key]=sattr;
	}
	return(result);
}
		
	
int main(int argc, char **argv)
{
	map<string,SHHdrAttribute> shhdrmap;
	shhdrmap=load_shhdrmap(pf);
}
