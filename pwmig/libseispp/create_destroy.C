#include "seispp.h"
//
// constructors for the Time_Series object are defined inline
// in seispp.h.  The only exception is this, which is the copy
// constructor
//
Time_Series::Time_Series(const Time_Series& tsi)
{
	int i;
        dt=tsi.dt;
        t0=tsi.t0;
        ns=tsi.ns;
        tref=tsi.tref;
        md= Metadata(tsi.md);
        s=new double[ns];
        for(i=0;i<tsi.ns;++i) s[i]=tsi.s[i];
}
Time_Series::Time_Series(const Time_Series *tsi)
{
	int i;

        t0=tsi->t0;
        ns=tsi->ns;
        dt=tsi->dt;
        tref=tsi->tref;
        md= Metadata(tsi->md);
        s=new double[ns];
        for(i=0;i<tsi->ns;++i) s[i]=tsi->s[i];
}
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
Three_Component_Seismogram::Three_Component_Seismogram
			(const Three_Component_Seismogram& t3c)
{
	int i,j;
	md=Metadata(t3c.md);
        dt=t3c.dt;
        t0=t3c.t0;
        ns=t3c.ns;
        tref=t3c.tref;
        components_are_orthogonal=t3c.components_are_orthogonal;
        components_are_cardinal=t3c.components_are_cardinal;
	for(i=0;i<3;++i)
		for(j=0;j<3;++j) tmatrix[i][j]=t3c.tmatrix[i][j];
	for(i=0;i<3;++i) x[i]=Time_Series((t3c.x)+i);
}
//No destructor is needed for current Three_Component_Seismogram object
// because the Time_Series destructor and default deletions should
// handle it fine.  

//
// Ensemble constructors.  Both just create blank trace or 3c trace objects
// and push them into a vector container.
//
Time_Series_Ensemble::Time_Series_Ensemble(int nensemble, int nsamples)
{
	for(int i=0; i<nensemble; ++i)
	{
		Time_Series *ts = new Time_Series(nsamples);
		tse.push_back(*ts);
	}
}
Three_Component_Ensemble::Three_Component_Ensemble(int nstations, int nsamples)
{
	for(int i=0;i<nstations;++i)
	{
		Three_Component_Seismogram *tcs 
			= new Three_Component_Seismogram(nsamples);
		tcse.push_back(*tcs);
	}
}
