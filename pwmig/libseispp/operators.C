#include "seispp.h"
Time_Series& Time_Series::operator=(const Time_Series& tsi)
{
	if(this!=&tsi)
	{
		dt=tsi.dt;
		t0=tsi.t0;
		ns=tsi.ns;
		tref=tsi.tref;
		md=tsi.md;  //Metadata assignment operator must create destroy
		if(s!=NULL) delete [] s;
		s=new double[ns];
		for(int i=0;i<tsi.ns;++i) s[i]=tsi.s[i];
	}
	return(*this);
}			
