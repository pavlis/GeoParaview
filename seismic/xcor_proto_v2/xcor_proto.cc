// This is a crude driver for version 2 of the array cross
// correlation software.  Version 1 used a driver written by
// Peng Wang that is currently broken from some unresolved xlib
// problems.  This version was produced by glp by simplifying
// the plot library and using a simple, crude command driven
// interface.  This version has this usage:
//    dbxcor db [-pf pffile] < hypo.list
//
#include <iostream>
#include "filter++.h"
#include "AnalysisSetting.h"
#include "XcorProcessingEngine.h"

using namespace std;
using namespace SEISPP;
void usage()
{
	cerr << "dbxcor dbin dbout hypofile [-pf pffile] " <<endl;
	exit(-1);
}
void show_command_menu()
{
	cout << "Enter:"<<endl
		<<"  c - change picks" <<endl
		<<"  a - run analysis" <<endl
		<<"  f - run secondary filter" <<endl
		<<"  x - show cross-correlation functions" <<endl
		<<"  b - show computed beam" <<endl
		<<"  p - pick arrival on beam" <<endl
		<<"  s - skip this event" << endl
		<<"  e - edit mode (pick traces to delete)" <<endl
		<<"  d - done with this event" <<endl			
		<<"  q - quit processing" <<endl;			
}
enum CommandChoice {ChangePick,Analyze,Filter,ShowCorrelations,ShowBeam,
	PickBeam, SkipEvent, Edit, Done,Quit,
	Bad};
CommandChoice get_command()
{
	char sc;
	cin >> sc;
	string s(" ");
	s[0]=sc;
	s[1]='\0';
	if(s=="c")
		return ChangePick;
	else if(s=="a")
		return Analyze;
	else if(s=="f")
		return Filter;
	else if(s=="x")
		return ShowCorrelations;
	else if(s=="b")
		return ShowBeam;
	else if(s=="p")
		return PickBeam;
	else if(s=="s")
		return SkipEvent;
	else if(s=="e")
		return Edit;
	else if(s=="d")
		return Done;
	else if(s=="q")
		return Quit;
	else
		return Bad;
}
int main(int argc, char **argv)
{
	int i;
	if(argc<4) usage();
	string waveform_db_name(argv[1]);
	string result_db_name(argv[2]);
	string hypofile(argv[3]);
	string pfname("dbxcor");
	for(i=4;i<argc;++i)
	{
		string argtest(argv[i]);
		if(argtest=="-pf")
		{
			++i;
			pfname=string(argv[i]);
		}
		else
		{
			usage();
		}
	}
	Pf *pf;
	if(pfread(const_cast<char *>(pfname.c_str()),&pf))
	{
		cerr << "Error reading pf file = "<<pfname<<endl;
		exit(-1);
	}
	try {
	Metadata global_md(pf);
	AnalysisSetting asetting(global_md);
	XcorProcessingEngine xpe(pf,asetting,waveform_db_name,result_db_name);
	double lat,lon,depth,otime;
	const string method("tttaup");
	const string model("iasp91");
	MultichannelCorrelator *mcc;
	char yrmonday[15],day[6],hrminsec[20];
	FILE *fp=fopen(hypofile.c_str(),"r");
	// for now we just count this.  Eventually it will be a 
	// required input.
	int orid=0,evid=0;
	while(fscanf(fp,"%lf%lf%lf%s%s%s",&lat,&lon,&depth,
			yrmonday,day,hrminsec)!=EOF)
	{
		string datestring=string(yrmonday)+" "+string(hrminsec);
		otime=str2epoch(const_cast<char*>(datestring.c_str()));
		Hypocenter h(rad(lat),rad(lon),depth,otime,method,model);
		cout << "Loading data for event: "
			<< lat<<","
			<< lon<<","
			<< depth <<","
			<< datestring<<endl;
		xpe.load_data(h);
		cout << "Data loaded"
			<<endl
			<<"Displaying filtered data"<<endl;
		++orid;
		++evid;
		xpe.do_picks();
		CommandChoice choice;
		PointPick phase_pick;
		// This allows skipping pick on beam.  Defaults 0.  
		// Doesn't really do anything in this prototype, but placed
		// here to emphasize the need for the capability in the 
		// final interface.  glp:  1/13/2006
		double tshift=0.0;
		do {
			string filt;
			char line[80];
			show_command_menu();
			choice=get_command();

			switch (choice)
			{
			case ChangePick:
				xpe.do_picks();
				break;
			case Analyze:
				cout << "Starting analysis processing; stand by"<<endl;
				mcc=xpe.analyze();
				break;
			case Filter:
				cin.getline(line,80);
				filt=string(line);
				try {
					xpe.filter_data(TimeInvariantFilter(filt));
				} catch (SeisppError serr)
				{
					cerr << "filter_data failed using "
					  << filt <<" try again"<<endl;
					serr.log_error();
				}
				break;
			case ShowCorrelations:
				xpe.display_correlations();
				break;
			case ShowBeam:
				xpe.display_beam();
				break;
			case PickBeam:
				phase_pick=xpe.pick_beam();
				xpe.shift_arrivals(phase_pick.time);
				break;

			case Done:
				try {
					xpe.save_results(evid,orid);
				} catch (SeisppError serr)
				{
					cerr << "Problems in save_results routine"<<endl;
					serr.log_error();
					cerr << "Try again or exit"<<endl;
				}
				break;
			case Edit:
				xpe.edit_data();
			case SkipEvent:
				choice=Done;
				break;
			case Quit:
				exit(0);
			default:
				cout << "Illegal/unknown command"<<endl;
			}
		}
		while (choice != Done);
	}
	}
	catch (SeisppError serr)
	{
		serr.log_error();
		cerr << "Fatal error:  exiting"<<endl;
	}
}
