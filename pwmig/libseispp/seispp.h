#ifndef _SEISPP_H_
#define _SEISPP_H_
#include <vector>
#include "metadata.h"


enum Time_Reference_Type {absolute,relative};
class Time_Series
{
public:
	Metadata md;
	double dt,t0;
	int ns;
	Time_Reference_Type tref;
	double *s;
	Time_Series() {s=NULL;dt=0.0; t0=0.0;};
	Time_Series(int nsin) {dt=0.0; t0=0.0; s = new double(nsin);};
	~Time_Series() {  if(s!=NULL) delete s; };
};

class Three_Component_Seismogram : public Time_Series
{
public:
	Metadata md;
	//Using a Time_Series object adds overhead and memory use 
	//by duplicating metadata, but adds generality
	Time_Series x1, x2, x3;
	Three_Component_Seismogram();
	Three_Component_Seismogram(int nsamp);
	~Three_Component_Seismogram();
};
class Time_Series_Ensemble
{
public:  
	Metadata md;
	vector <Time_Series> ts;

	Time_Series_Ensemble();
	Time_Series_Ensemble(int nts);
	Time_Series_Ensemble(Pf_ensemble *pfe);
	// This function is bad form, but I see no alernative
	// It's purpose is to initialize the namespace for the 
	// global metadata list.  It should not be advertised
	void add_to_ensemble_mdlist(char *s){pushtbl(mdlist,s);};
private:
	// mdlist contains a list of names to be pulled from input pf 
	// and copied to global Metadata area for the ensemble
	// Note no typing is needed as this is a string copy
	Tbl *mdlist;
};

class Three_Component_Ensemble
{
public:
	Metadata md;
	vector <Three_Component_Seismogram> tcs;

	Three_Component_Ensemble();
	Three_Component_Ensemble(int nts);
	Three_Component_Ensemble(Pf_ensemble *pfe);
	void add_to_ensemble_mdlist(char *s){pushtbl(mdlist,s);};
private:
	Tbl *mdlist;
};
#endif
