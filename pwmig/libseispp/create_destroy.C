#include "seispp.h"
//
// constructors for the Time_Series object are defined inline
// in seispp.h
//
// here is the parameterized constructor for a 3c seismogram
//
Three_Component_Seismogram::Three_Component_Seismogram(int nsin) : :
{
	for(int i=0; i<3;++i) 
		x[i].s = new double(nsin);
	
}
//No destructor is needed for current Three_Component_Seismogram object
// because the Time_Series destructor and default deletions should
// handle it fine.
//
// Ensemble constructors.  Both just create blank trace or 3c trace objects
// and push them into a vector container.
//
Time_Series_Ensemble::Time_Series_Ensemble(int nts, int nsamp)
{
	for(int i=0; i<nts; ++i)
	{
		Time_Series *ts = new Time_Series(nsamp);
		tse.push_back(ts);
	}
}
Three_Component_Ensemble::Three_Component_Ensemble(int nsta, int nsamp)
{
	for(int i=0;i<nsta;++i)
	{
		Three_Component_Seismogram *tcs 
			= new Three_Component_Seismogram(nsamp);
		tcse.push_back(tcs);
	}
}
