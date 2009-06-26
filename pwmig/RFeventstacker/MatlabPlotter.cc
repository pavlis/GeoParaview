#include <sstream>
#include "MatlabPlotter.h"
MatlabPlotter::MatlabPlotter() : MatlabProcessor()
{
	if(SEISPP_verbose) cout << "MatlabProcessor created successfully"<<endl;
}
void MatlabPlotter::wigbplot(ThreeComponentEnsemble& d,string chans[3],bool plot_image);
{
    try {
	string buf[256];
	ostringstream ss(buf);
	string command;
	this->load(d,chans);
	ss.clear();
	ss<<"figure("<<fn3ce[0]<<");"<<endl;
	command=ss.str();
	this->process(command);
	command=wigb("+chan[0]+string(");");
	this->process(command);
	if(plot_image)
	{
		command="imagesc("+chan[0]+string(");");
		this->process(command);
	}
	command="title '"+chan[0]+"'";
	this->process(command);

	ss.clear();
	ss<<"figure("<<fn3ce[1]<<");"<<endl;
	command=ss.str();
	this->process(command);
	command=wigb("+chan[1]+string(");");
	this->process(command);
	if(plot_image)
	{
		command="imagesc("+chan[1]+string(");");
		this->process(command);
	}
	command="title '"+chan[1]+"'";
	this->process(command);

	ss.clear();
	ss<<"figure("<<fn3ce[2]<<");"<<endl;
	command=ss.str();
	this->process(command);
	command=wigb("+chan[2]+string(");");
	this->process(command);
	if(plot_image)
	{
		command="imagesc("+chan[2]+string(");");
		this->process(command);
	}
	command="title '"+chan[2]+"'";
	this->process(command);
	
    }
    catch(...){throw;};
}
void MatlabPlotter::wigbplot(TimeSeriesEnsemble& d)
{
    try {
	this->load(d,string("d"));
	string buf[256];
	ostringstream ss(buf);
	string command;
	ss.clear();
	ss<<"figure("<<fntse<<");"<<endl;
	command=ss.str();
	this->process(command);
	command="wigb(d);";
	this->process(command);
	if(plot_image)
	{
		command="imagesc(d);";
		this->process(command);
	}
	command="title 'TimeSeriesEnsemble data'";
	this->process(command);
    }catch(...){throw;};
}
void MatlabPlotter::plot(ThreeComponentSeismogram& d)
{
	try {
		this->load(d,string("d"));
		string buf[256];
		ostringstream ss(buf);
		string command;
		ss.clear();
		ss<<"figure("<<fntcs<<");"<<endl;
		command=ss.str();
		this->process(command);
		command=string("subplot(3,1,1), plot(d(:,1));")
			+ string("subplot(3,1,2),plot(d(:,2));")
			+ string("subplot(3,1,3),plot(d(:,3));");
		this->process(command);
		command="title 'ThreeComponentSeismogram data';";
		this->process(command);
	}catch(...){throw;};
}
void MatlabPlotter::plot(TimeSeries& d)
{
	try {
		this->load(d,string("d"));
		string buf[256];
		ostringstream ss(buf);
		string command;
		ss.clear();
		ss<<"figure("<<fntcs<<");"<<endl;
		command=ss.str();
		this->process(command);
		command="plot(d);";
		this->process(command);
		command="title 'TimeSeries data';";
		this->process(command);
	}catch(...){throw;};
}
