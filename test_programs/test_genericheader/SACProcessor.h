#ifndef _SAC_PROCESSOR_H_
#define _SAC_PROCESSOR_H_

#include <ifstream>
#include "ExternalProcessor.h"
#include "TimeSeries.h"
#include "ThreeComponentSeismogram.h"
#include "ensemble.h"

namespace SEISPP 
{
using namespace std;
using namespace SEISPP;

class SACProcessor : public ExternalProcessor
{
public:
	SACProcessor();
	SACProcessor(string fname);
	~SACProcessor();
	/*! \brief Load a TimeSeries object.
	
	\param d TimeSeries data to be loaded for processing.
	\param name file name to use to store these data on /tmp
	*/  
	void load(TimeSeries& d,string name,map<string,string>formap);
	void load(ThreeComponentSeismogram& d,string name);
	void load(TimeSeriesEnsemble& tse,string name);
	void load(ThreeComponentEnsemble& tce,string name[3]);
	void process(string text);
	void process(char *s);
	void process(ifstream& strm); 
	TimeSeries retrieve(string name,map<string,string> invmap);
private:
	ofstream sac_command_input;
	pid_t sac_pid;
	HeaderMap hdr;
};

} // End SEISPP namespace encapsulation
#endif
