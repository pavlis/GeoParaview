#include "seispp.h"
#include "filter++.h"
using namespace SEISPP;

Time_Invariant_Filter::Time_Invariant_Filter(string fspec)
{
	filter_spec=fspec;
	char ftkey[20];
	sscanf(fspec.c_str(),"%s",ftkey);
	string sftkey;
	if(sftkey=="BW")
	{


}
Time_Invariant_Filter::fmin()
{
	switch (type)
	{
	case highpass:
	case lowpass:
	case bandpass:
		return(f1);
	default:
		throw seispp_error(
			string("fmin function not defined for filter type")
				+ this->type_description);
}
	
