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
		<<"  a - run analysis" <<endl
		<<"  b - show computed beam" <<endl
		<<"  c - change cutoff criteria"<<endl
		<<"  C - pick cutoff interactively"<<endl
		<<"  d - done with this event" <<endl			
		<<"  e - mark traces to delete"<<endl
		<<"  F - run custom filter"<<endl
		<<"  f - run SP filter" <<endl
		<<"  P - change picks" <<endl
		<<"  p - pick arrival on beam" <<endl
		<<"  r - restore original ensemble" <<endl
		<<"  s - skip this event" << endl
		<<"  x - show cross-correlation functions" <<endl
		<<"  Q - quit without saving" <<endl
		<<"  q - quit processing saving current event" <<endl;			
}
enum CommandChoice {ChangePick,Analyze,Edit,Filter,ShowCorrelations,ShowBeam,
	PickBeam, SkipEvent, SortCriteria, Cutoff, Done,Quit,
	CustomFilter,QuitNoSave,RestoreOriginal,
	Bad};
CommandChoice get_command()
{
	char sc[128];
	cin.getline(sc,128);
	string s(sc);
	if(s=="P")
		return ChangePick;
	else if(s=="a")
		return Analyze;
	else if(s=="f")
		return Filter;
	else if(s=="e")
		return Edit;
	else if(s=="x")
		return ShowCorrelations;
	else if(s=="b")
		return ShowBeam;
	else if(s=="p")
		return PickBeam;
	else if(s=="s")
		return SkipEvent;
	else if(s=="c")
		return SortCriteria;
	else if(s=="C")
		return Cutoff;
	else if(s=="d")
		return Done;
	else if(s=="q")
		return Quit;
	else if(s=="F")
		return CustomFilter;
	else if(s=="Q")
		return QuitNoSave;
	else if(s=="r")
		return RestoreOriginal;
	else
		return Bad;
}
SortOrder ask_for_sort_order()
{
	cout << "Enter one of:"<<endl
		<<" c - coherence"<<endl
		<<" x - peak cross correlation"<<endl
		<<" w - stack weight"<<endl
		<<"->";
	string ans;
	cin>>ans;
	if(ans=="c")
		return COHERENCE;
	else if(ans=="x") 
		return CORRELATION_PEAK;
	else if(ans=="w")
		return WEIGHT;
	else
	{
		cerr << "Illegal response = "<<ans
			<<" defaulting to coherence"<<endl;
		return COHERENCE;
	}
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
	AnalysisSetting asetting_default(asetting);
	XcorProcessingEngine xpe(pf,asetting,waveform_db_name,result_db_name);
	double lat,lon,depth,otime;
	const string method("tttaup");
	const string model("iasp91");
	MultichannelCorrelator *mcc;
	char yrmonday[15],day[6],hrminsec[20];
	FILE *fp=fopen(hypofile.c_str(),"r");
	// for now we just count this.  Eventually it will be a 
	// required input.
	int orid,evid;
	while(fscanf(fp,"%d%d%lf%lf%lf%s%s%s",&evid,&orid,
			&lat,&lon,&depth,
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
		asetting.set_ref_trace(xpe.pick_one_trace());
		xpe.change_analysis_setting(asetting);
		CommandChoice choice;
		PointPick phase_pick;
		// This allows skipping pick on beam.  Defaults 0.  
		// Doesn't really do anything in this prototype, but placed
		// here to emphasize the need for the capability in the 
		// final interface.  glp:  1/13/2006
		double tshift=0.0;
		bool beam_picked=false;
		do {
			string filt;
			char line[80];
			show_command_menu();
			choice=get_command();

			switch (choice)
			{
			case ChangePick:
				xpe.do_all_picks();
				break;
			case Analyze:
				cout << "Starting analysis processing; stand by"<<endl;
				mcc=xpe.analyze();
				beam_picked=false;
				break;
			case Edit:
				xpe.edit_data();
				break;

			case CustomFilter:
			case Filter:
				if(choice==Filter)
					filt=string("BW 0.5 5 2 5");
				else
				{
					cout << "Enter custom filter string (Antelope convention):";
					char buffer[128];
					cin.getline(buffer,128);
					filt=string(buffer);
				}
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
				beam_picked=true;
				break;
			case SortCriteria:
				asetting.result_sort_order=ask_for_sort_order();
				xpe.change_analysis_setting(asetting);
				break;
				
			case Cutoff:
				xpe.sort_ensemble();
				xpe.pick_cutoff();
				break;
			case RestoreOriginal:
				xpe.restore_original_ensemble();
				break;
			case Done:
			case Quit:
				if(!beam_picked)
				{
					phase_pick=xpe.pick_beam();
					xpe.shift_arrivals(phase_pick.time);
				}

				try {
					xpe.save_results(evid,orid);
				} catch (SeisppError serr)
				{
					cerr << "Problems in save_results routine"<<endl;
					serr.log_error();
					cerr << "Try again or exit"<<endl;
				}
				if(choice==Quit) exit(0);
				break;
			case QuitNoSave:
				exit(0);
			case SkipEvent:
				choice=Done;
				break;
			default:
				cout << "Illegal/unknown command"<<endl;
			}
		}
		while (choice != Done);
		xpe.change_analysis_setting(asetting_default);
	}
	}
	catch (SeisppError serr)
	{
		serr.log_error();
		cerr << "Fatal error:  exiting"<<endl;
	}
}
