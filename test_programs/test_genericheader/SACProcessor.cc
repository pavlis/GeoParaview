#include <stdio.h>
#include <iostream>
#include <string>
#include <sys/types.h>
#include <unistd.h>
#include "SACprocessor.h"
namespace SEISPP
{

using namespace std;
const string TmpFileBase("/tmp/SEISPP");
SACProcessor::SACProcessor()
{
	ios::sync_with_stdio();
	string base_error("SACProcessor default constructor: ");
	string default_fifo("/tmp/SACProcessor.fifo");
	if(mkfifo(default_fifo.c_str(),S_IRWXU))
		throw SeisppError(base_error
			+ string("cannot open fifo communication channel=")
			+ default_fifo);
	pid_t pID=fork();
	if (pID == 0) // This is for child process
	{
		execl("sac","<",default_fifo.c_str(),0);
		throw SeisppError(base_error
			+ string("execl failed"));
	}
	else if (pID < 0)  
	{
		throw SeisppError( base_error
			+ string("fork failed"));
	}
	else
	{
		sac_pid=pID;
		// parent needs to open the fifo and store it's definition
		// in the private area
		try {
			sac_command_input.open(default_fifo,ios::out);
		} catch (ios::failure& var)
		{
			string mess=var.what();
			throw SeisppError(base_error
				+ string("parent process cannot open fifo communication channel=")
				+ default_fifo
				+ string("\nofstream open method threw this error: ")
				+ mess);
		}
}
SACProcessor::~SACProcessor()
{
	sac_command_input << "exit";
	sac_command_input.close();
}
void SACProcessor::load(TimeSeries& d,string name)
{
    try {
	this->load(&(d.s[0]),d.ns,name);
    } catch (SeisppError serr) throw;
}
} // End namespace SEISPP declaration
