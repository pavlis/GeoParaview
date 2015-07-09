#include <string>
#include <map>
#include <iostream>
#include <typeinfo>
#include "stock.h"
#include "pf.h"
#include "seispp.h"
#include "AttributeMap.h" // needed for MDtype definition only
using namespace std;
using namespace SEISPP;
class HeaderAttributeProperties
{
public:
	HeaderAttributeProperties();
	HeaderAttributeProperties(string nm, MDtype mdtin, int nb, int o);
	HeaderAttributeProperties(char *line);
	string name;
	MDtype mdt;
	int nbytes;
	size_t offset;
};
HeaderAttributeProperties::HeaderAttributeProperties()
{
	name=string("INVALID");
	mdt=MDinvalid;
	nbytes=0;
	offset=0;
}
HeaderAttributeProperties::HeaderAttributeProperties(string nm, MDtype mdtin, int nb, int o)
{
	name=nm;
	mdt=mdtin;
	nbytes=nb;
	offset=(size_t) o;
}

HeaderAttributeProperties::HeaderAttributeProperties(char *line)
{
	char cnm[64];
	char tpnm[64];
	int nb,off;
	sscanf(line,"%s%s%d%d",cnm,tpnm,&nb,&off);
	name=string(cnm);
	string tname(tpnm);
	if( (tname==string("REAL")) || tname==string("real") )
		mdt=MDreal;
	if( (tname==string("INT")) || tname==string("int") )
		mdt=MDint;
	if( (tname==string("STRING")) || tname==string("string") )
		mdt=MDstring;
	if( (tname==string("BOOLEAN")) || tname==string("boolean") )
		mdt=MDboolean;
	nbytes=nb;
	offset=(size_t) off;
}

	
class HeaderMap
{
public:
	HeaderMap();
	HeaderMap(Pf *pf, string tag);
	template <class T> T get(string name,unsigned char *h);
	string get_string(string name, unsigned char *h);
	template <class T> void put(string name, T value,unsigned char *h);
	void put_string(string name,string value,unsigned char *h);
	void put_string(string name,char *value,unsigned char *h);
	int size(){return(headersize);};
private:
	int headersize;
	map<string,HeaderAttributeProperties> attributes;
};

template <class T> T HeaderMap::get(string name,unsigned char *headerdata)
{
	map<string,HeaderAttributeProperties>::iterator ithm;
	ithm=attributes.find(name);
	if(ithm==attributes.end())
		throw SeisppError(string("HeaderMap::get attribute name=")
			+ name
			+ string("is not defined for this header type"));
	int nb=ithm->second.nbytes;
	if(sizeof(T)!=nb)
		throw SeisppError(string("HeaderMap::get attribute name=")
			+ name 
			+ string("type size mismatch") );
	++nb;  // Need to add one to allow for trailing null
	size_t o=ithm->second.offset;
	T *val;
	val=reinterpret_cast<T*>(headerdata+o);
	return(*val);
}
template <class T> void HeaderMap::put(string name, T value,unsigned char *headerdata)
{
	map<string,HeaderAttributeProperties>::iterator ithm;
	ithm=attributes.find(name);
	if(ithm==attributes.end())
		throw SeisppError(string("HeaderMap::put attribute named=")
			+ name
			+ string("is not defined for this header type"));
	int nb=ithm->second.nbytes;
	if(sizeof(T)!=nb)
		throw SeisppError(string("HeaderMap::put attribute name=")
			+ name 
			+ string("type size mismatch") );
	// This needs to eventually have an error handler for invalid 
	//MDtype tnm(ithm->second.mdt);
	size_t o=ithm->second.offset;
	T *hptr;
	unsigned char *a;
	a=headerdata+o;
	hptr=reinterpret_cast<T*>(a);
	*hptr=value;
}
string HeaderMap::get_string(string name, unsigned char *h)
{
	map<string,HeaderAttributeProperties>::iterator ithm;
	ithm=attributes.find(name);
	if(ithm==attributes.end())
		throw SeisppError(string("HeaderMap::get attribute name=")
			+ name
			+ string("is not defined for this header type"));
	int nb=ithm->second.nbytes;
	nb=nb+1;
	char *s=new char[nb];
	char *hs=reinterpret_cast<char *>(h);
	for(int i=0;i<nb-1;++i)
		s[i]=hs[i];
	h[nb]='\0';
	string result(s);
	delete s;
	return(result);
}
void HeaderMap::put_string(string name,char *s, unsigned char *h)
{
	map<string,HeaderAttributeProperties>::iterator ithm;
	ithm=attributes.find(name);
	if(ithm==attributes.end())
		throw SeisppError(string("HeaderMap::put attribute named=")
			+ name
			+ string("is not defined for this header type"));
	int nb=ithm->second.nbytes;
	int i;
	for(i=0;i<nb;++i) h[i]='\0';
	// reduce the copy length if s is short
	if(strlen(s)<nb) nb=strlen(s);
	char *scpy=reinterpret_cast<char *>(h);
	for(i=0;i<nb;++i) scpy[i]=s[i];
}
void HeaderMap::put_string(string name,string s, unsigned char *h)
{
	this->put_string(name,s.c_str(),h);
}
HeaderMap::HeaderMap()
{
	headersize=0;
}
HeaderMap::HeaderMap(Pf *pfin,string tag)
{
	Pf *pf;
	// tag is a header type name nested with an Arr
	// Use pfget to grab the one we want.
	if(pfget(pfin,const_cast<char *>(tag.c_str()),(void **)&pf) != PFARR)
		throw SeisppError(
			string("Parameter file has no definition for header named ")
			+ tag);
	Tbl *t;
	t=pfget_tbl(pf,"attributes");
	if(t==NULL)
		throw SeisppError(
			string("Parameter file is missing required Tbl parameter=")
			+ tag);

	for(int i=0;i<maxtbl(t);++i)
	{
		char *line;
		line=(char *)gettbl(t,i);
		HeaderAttributeProperties hap(line);
		attributes[hap.name]=hap;
	}
	freetbl(t,0);
	headersize=pfget_int(pf,"headersize");
}
int main(int argc, char **argv)
{
	Pf *pf;
	if(pfread("test",&pf))
	{
		cerr << "Error reading pf file=test"<<endl;
		exit(-1);
	}
	try {
	string tag("testheader");
	int headsize=pfget_int(pf,"headersize");
	unsigned char *hbuf=new unsigned char[headsize];
	HeaderMap h(pf,tag);
	h.put<double>(string("dtest"),4.0,hbuf);
	double x=h.get<double>(string("dtest"),hbuf);
	cout << "dtest="<<x<<endl;
	h.put<int>(string("itest"),16,hbuf);
	h.put<float>(string("rtest"),555.0,hbuf);
	cout << "rtest="<<h.get<float>(string("rtest"),hbuf)<<endl;
	cout << "itest="<<h.get<int>(string("itest"),hbuf)<<endl;
	string s("teststr");
	h.put_string(string("stest"),s,hbuf);
	cout << "stest="<<h.get_string(string("stest"),hbuf)<<endl;
	} catch (SeisppError serr)
	{
		serr.log_error();
	}
}
