#include "seispp.h"
#include "filter++.h"
using namespace SEISPP;
namespace SEISPP
{

Time_Invariant_Filter::Time_Invariant_Filter(string fspec)
{
	filter_spec=fspec;
	char ftkey[20];
	sscanf(fspec.c_str(),"%s",ftkey);
	string sftkey(ftkey);
	if(sftkey=="BW")
	{
		sscanf(fspec.c_str(),"%s%lf%d%lf%d",ftkey,
			&f1,&npole1,&f2,&npole2);
		if( npole1==0 && f1<=0.0 && npole2==0 && f2<=0.0)
			throw seispp_error(string("Illegal filter definition = ")
				+ fspec);
		else if(npole1==0 || f1<=0.0)
			type = lowpass;
		else if(npole2==0 || f2<=0.0)
			type = highpass;
		else
			type = bandpass;
	}
	else
	{
		npole1 = 0;  npole2=0;
		f1=0.0;  f2=0.0;
	}
}
Time_Invariant_Filter::Time_Invariant_Filter(double flow, int npl, 
	double fhigh, int nph)
{
	type = bandpass;
	npole1=npl;  npole2=nph;
	f1=flow;  f2=fhigh;
	char s[50];
	sprintf(s,"BW %lf %d %lf %d",f1,npole1,f2,npole2);
	filter_spec=string(s);
}
// copy constructor -- perhaps unnecessary, but books all say never
// depend on the default
Time_Invariant_Filter::Time_Invariant_Filter(const Time_Invariant_Filter& tin)
{
	npole1=tin.npole1;
	npole2=tin.npole2;
	f1=tin.f1;
	f2=tin.f2;
	filter_spec=tin.filter_spec;
}
// same for assignment operator here
Time_Invariant_Filter& Time_Invariant_Filter::operator=(const Time_Invariant_Filter& tin)
{
	if(&tin!=this)
	{
		npole1=tin.npole1;
		npole2=tin.npole2;
		f1=tin.f1;
		f2=tin.f2;
		filter_spec=tin.filter_spec;
	}
	return(*this);
}

	
double Time_Invariant_Filter::fmin()
{
	switch (type)
	{
	case highpass:
	case bandpass:
		return(f1);
	case lowpass:
	default:
		throw seispp_error(
			string("fmin not defined for filter type")
				+ this->filter_spec);
	}
}
double Time_Invariant_Filter::fmax()
{
	switch (type)
	{
	case lowpass:
	case bandpass:
		return(f2);
	case highpass:
	default:
		throw seispp_error(
			string("fmax not defined for filter type")
				+ this->filter_spec);
	}
}
int Time_Invariant_Filter::fmin_poles()
{
	switch (type)
	{
	case highpass:
	case lowpass:
	case bandpass:
		return(npole1);
	default:
		throw seispp_error(
			string("number poles for fmin corner not defined for filter type")
				+ this->filter_spec);
	}
}
int Time_Invariant_Filter::fmax_poles()
{
	switch (type)
	{
	case highpass:
	case lowpass:
	case bandpass:
		return(npole2);
	default:
		throw seispp_error(
			string("number poles for fmax corner not defined for filter type")
				+ this->filter_spec);
	}
}
string Time_Invariant_Filter::type_description()
{
	string retstr;
	switch (type)
	{
	case highpass:
		retstr=string("highpass: ") + filter_spec;
		break;
	case lowpass:
		retstr=string("lowpass: ")+filter_spec;
		break;
	case bandpass:
		retstr=string("bandpass: ")+filter_spec;
		break;
	case WAA:
		retstr=string("Wood-Anderson");
		break;
	case WAV:
		retstr=string("Wood-Anderson Velocity");
		break;
	case WAD:
		retstr=string("Wood-Anderson Displacement");
		break;
	case DIF:
		retstr=string("First Derivative");
		break;
	case DIF2:
		retstr=string("Second Derivative");
		break;
	case INT:
		retstr=string("Integrator");
		break;
	case INT2:
		retstr=string("Double Integrator");
		break;
	case DEMEAN:
		retstr=string("Mean Removal filter");
		break;
	default:
		retstr=string("Filter type is undefined");
	}
	return(retstr);
}
// apply to simple float vector of length ns and sample rate dt
void Time_Invariant_Filter::apply(int ns, float *s,double dt)
{
	if(trfilter_segs(1,&ns,&dt,&s,const_cast<char*>(filter_spec.c_str()))<0)
			throw seispp_error(string("Error in trfilter_segs"));
}
// same as above for array of doubles
void Time_Invariant_Filter::apply(int ns, double *s,double dt)
{
	int i;
	float *d=new float[ns];
	for(i=0;i<ns;++i) d[i]=static_cast<float>(s[i]);
	if(trfilter_segs(1,&ns,&dt,&d,const_cast<char*>(filter_spec.c_str()))<0)
			throw seispp_error(string("Error in trfilter_segs"));
	for(i=0;i<ns;++i) s[i]=static_cast<double>(d[i]);
	delete [] d;
}
void Time_Invariant_Filter::apply(Time_Series& ts)
{
	int i;
	float *d=new float[ts.ns];
	for(i=0;i<ts.ns;++i) d[i]=static_cast<float>(ts.s[i]);
	if(trfilter_segs(1,&(ts.ns),&(ts.dt),&d,const_cast<char*>(filter_spec.c_str()))<0)
			throw seispp_error(string("Error in trfilter_segs"));
	for(i=0;i<ts.ns;++i) ts.s[i]=static_cast<double>(d[i]);
	delete [] d;
}
void Time_Invariant_Filter::apply(Three_Component_Seismogram& ts)
{
	int i,j;
	float *d=new float[ts.ns];
	for(j=0;j<3;++j)
	{
		for(int i=0;i<ts.ns;++i) d[i]=static_cast<float>(ts.u(j,i));
		if(trfilter_segs(1,&(ts.ns),&(ts.dt),&d,const_cast<char*>(filter_spec.c_str()))<0)
			throw seispp_error(string("Error in trfilter_segs"));
		for(i=0;i<ts.ns;++i) ts.u(j,i)=static_cast<double>(d[i]);
	}
	delete [] d;
}
void Time_Invariant_Filter::apply(Dbptr tr)
{
	if(trfilter(tr,const_cast<char*>(filter_spec.c_str()))<0)
		throw seispp_error(string("Error in trfilter"));
}

// helpers for ensembles.  There is probably a way to do this with templates,
// but it isn't that much code
void Filter_Ensemble(Time_Series_Ensemble& ensemble,Time_Invariant_Filter& filter)
{
	try 
	{
	    for(int i=0;i<ensemble.tse.size();++i)
		filter.apply(ensemble.tse[i]);
	} catch (...) {throw;};
}
void Filter_Ensemble(Three_Component_Ensemble& ensemble,Time_Invariant_Filter& filter)
{
	try 
	{
	    for(int i=0;i<ensemble.tcse.size();++i)
		filter.apply(ensemble.tcse[i]);
	} catch (...) {throw;};
}

}  // end namespace encapsulation for SEISPP
