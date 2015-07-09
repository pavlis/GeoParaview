#include "stock.h"
#include "pf.h"
#include "seispp.h"
#include "SeismicPlot.h"
int main(int argc, char **argv)
{
	int i,j;
	try{
	// Build an empty 10 element ensemble
	TimeSeriesEnsemble ensemble(10,500);
	for(i=0;i<10;++i)
	{
		TimeSeries d(500);
		for(j=0;j<500;++j) d.s[j]=0.0;
		//for(j=0;j<500;++j) d.s.push_back(0.0);
		d.s[10+i]=5.0;
		d.s[11+i]=7.0;
		d.s[12+i]=9.0;
		d.s[13+i]=15.0;
		d.s[14+i]=25.0;
		d.s[15+i]=5.0;
		d.s[20+i]=-5.0;
		d.s[21+i]=-7.0;
		d.s[22+i]=-9.0;
		d.s[23+i]=-15.0;
		d.s[24+i]=-25.0;
		d.s[25+i]=-5.0;
		d.dt=0.1;
		d.t0=0.0;
		d.ns=500;
		d.live=true;
		d.tref=relative;
		ensemble.member[i]=d;
	}
	Pf *pf;
	pfread("SeismicPlot",&pf);
	Metadata md(pf);
	SeismicPlot plot(ensemble,md);
	plot.draw();
	} catch (SeisppError err)
	{
		err.log_error();
		exit(-1);
	}
}
	
	
