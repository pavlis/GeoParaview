#ifndef _DECIMATE_H_
#define_DECIMATE_H_
#include <map>
#include "pf.h"

class Decimator
{
public:  
	// [start,end] define sample rate interval = exact
	double start;
	double end;
	double decfac;
	double exact;
	Tbl *declist; // Holds actual decimator defintion
	Decimator(){start=0.0;end=0.0;decfac=1.0;exact=0.0;declist=newtbl(0);};
};

class Dec_Cmp
{
public:
	bool operator()(const Decimator r1, const Decimator r2) const
	{return(r1.end<r2.start);};
};

class Decimation_Definitions
{
	map<Decimator,Dec_Cmp> decmap;
	Decimation_Definitions(string pffile);
	Decimation_Definitions(Pf *pf);
	~Decimation_Definitions();
};

// function prototypes that use the objects defined above

Time_Series *Decimate_Time_Series(Time_Series& ts, Decimation_Definitions& dd);


#endif
