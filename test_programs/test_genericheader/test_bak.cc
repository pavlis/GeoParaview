#include <string>
#include <map>
#include <iostream>
#include <typeinfo>
#include "stock.h"
#include "pf.h"
#include "seispp.h"
using namespace std;
using namespace SEISPP;
class HeaderAttributeProperties
{
public:
	HeaderAttributeProperties();
	HeaderAttributeProperties(string nm, string tname, int nb, int o);
	HeaderAttributeProperties(char *line);
	string name;
	string tname;
	int nbytes;
	size_t offset;
};
HeaderAttributeProperties::HeaderAttributeProperties()
{
	name=string("NONE");
	tname=string("string");
	nbytes=0;
	offset=0;
}
HeaderAttributeProperties::HeaderAttributeProperties(string nm, string tn, int nb, int o)
{
	name=nm;
	tname=tn;
	nbytes=nb;
	offset=(size_t) o;
}

HeaderAttributeProperties::HeaderAttributeProperties(char *line)
{
	char cnm[64];
	char tpnm[64];
	int nb,off;
	sscanf(line,"%s%s%n%n",cnm,tpnm,&nb,&off);
	name=string(cnm);
	tname=string(tpnm);
	nbytes=nb;
	offset=(size_t) off;
}

	
class HeaderMap
{
public:
	HeaderMap();
	HeaderMap(Pf *pf, string tag);
	template <class T> T get(string name,unsigned char *h);
	template <class T> void put(string name, T value,unsigned char *h);
	void put(string name,string value,unsigned char *h);
	void put(string name,char *value,unsigned char *h);
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
		throw SeisppError(string("HeaderMap::get attribute named=")
			+ name
			+ string("is not defined for this header type"));
	int nb=ithm->second.nbytes;
	++nb;  // Need to add one to allow for trailing null
	size_t o=ithm->second.offset;
	//if( (typeid(T)==typeid(basic_string)) || (typeid(T)==typeid(string)) )
	if (typeid(T)==typeid(string)) 
	{
		char *buf=new char[nb];
		for(int i=0;i<nb-1;++i) buf[i]=*(reinterpret_cast<char *>(headerdata+o+i));
		buf[nb]='\0';
		T result(buf);
		delete [] buf;
		return(result);
	}
	else
	{
		T *val;
		if(nb!=sizeof(T)|| typeid(T)!=typeid(tname))
			throw SeisppError(string("HeaderMap::get() type mismatch")
				+ string(" for tname =")+typeid(T).name() );
		val=reinterpret_cast<T*>(headerdata+o);
		return(*val);
	}
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
	string tnm(ithm->second.tname);
	if(nb!=sizeof(T)|| typeid(T)!=typeid(tnm.c_str()))
		throw SeisppError(string("HeaderMap::put() type mismatch for attribute=")
				+ name
				+string("\nRequested type=")+typeid(T).name()
				+string(" but header expects ")+tnm 
				+ string(" for attribute=")+name);
	size_t o=ithm->second.offset;
	T *hptr;
	unsigned char *a;
	a=headerdata+o;
	hptr=reinterpret_cast<T*>(a);
	*hptr=value;
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
}
