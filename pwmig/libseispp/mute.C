#include "seispp.h"

void apply_top_mute(Time_Series &ts,Top_Mute& mute)
{
	int i,i2;
	double t;
	if(ts.t0>mute.t1) return;
	for(i=0;i<ts.ns;++i)
	{
		t = ts.t0 +(ts.dt)*((double)i);
		if(t>mute.t0e) break;
		ts.s[i]=0.0;
	}
	if(i>=ns) return;
	for(i2=i;i2<ts.ns;++i2)
	{
		double weight;
		t = ts.t0 +(ts.dt)*((double)i2);
		if(t>mute.t1s) return;
		weight = (t-ts.t0)/(mute.t1s-mute.t0e);
		ts.s[i2]*=weight;
	}
}
void apply_top_mute_ensemble(Time_Series_Ensemble& t, Top_Mute& mute)
{
	vector iterator i;

	for(i=t.begin();i<t.end();++i)
	{
		apply_top_mute(t.tse[i],mute);
	}
}
void apply_top_mute_3c_ensemble(Three_Component_Ensemble &t3c, Top_Mute& mute)
{
	vector iterator i;
	int j;
	for(i=t3c.begin();i<t3c.end();++i)
	{
		for(j=0;j<2;++j) apply_top_mute(t3c.tcse[i][j],mute);
	}
}


