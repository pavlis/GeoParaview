#include "seispp.h"
//
// constructors for the Time_Series object are defined inline
// in seispp.h
//
// Default constructor for Three_Component_Seismogram could be 
// done inline in seispp.h, but it is complication enough I put
// it here
//
Three_Component_Seismogram::Three_Component_Seismogram()
{
	dt=0.0;
	t0=0.0;
	ns=0;
	tref=absolute;
	components_are_orthogonal=true;
	components_are_cardinal=true;
	for(int i=0;i<3;++i)
		for(int j=0;j<3;++j)
			if(i==j) 
				tmatrix[i][i]=1.0;
			else
				tmatrix[i][j]=0.0;
}

//
// here is the parameterized constructor for a 3c seismogram
//
Three_Component_Seismogram::Three_Component_Seismogram(int nsin) 
{
	dt=0.0;
	t0=0.0;
	ns=0;
	tref=absolute;
	components_are_orthogonal=true;
	components_are_cardinal=true;
	for(int i=0;i<3;++i)
		for(int j=0;j<3;++j)
			if(i==j) 
				tmatrix[i][i]=1.0;
			else
				tmatrix[i][j]=0.0;
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
