#ifndef _SEISPP_H_
#define _SEISPP_H_
#include <std>
#include "metadata.h"


enum Time_Reference_Type {absolute relative};
class Time_Series
{
public:
	Metadata md;
	double dt,t0;
	in ns;
	Time_Reference_Type tref;
	double *s;
	load_metadata(string); //Loads Metadata object by compiling a string
	Time_Series() {s=NULL;dt=0.0; t0=0.0;};
	Time_Series(int nsin) {dt=0.0; t0=0.0; s = new double(nsin);};
	~Time_series() {  if(s!=NULL) delete s; };
};

class Three_Component_Seismogram : public Time_Series
{
public:
	Metadata md;
	//Using a Time_Series object adds overhead and memory use 
	//by duplicating metadata, but adds generality
	Time_series x1, x2, x3;
	Three_Component_Seismogram();
	Three_Component_Seismogram(int nsamp);
	~Three_Component_Seismogram() { delete x1; delete x2; delete x3; };
};
class Time_Series_Ensemble
{
public:  
	Metadata md;
	vector <Time_Series> ts;

	Time_Series_Ensemble();
	Time_Series_Ensemble(int nts) {
	// This function is bad form, but I see no alernative
	// It's purpose is to initialize the namespace for the 
	// global metadata list.  It should not be advertised
	add_to_ensemble_mdlist(char *s){pushtbl(mdlist,s);
private:
	// mdlist contains a list of names to be pulled from input pf 
	// and copied to global Metadata area for the ensemble
	// Note not typing is needed as this is a string copy
	Tbl *mdlist=newtbl(0);
};
