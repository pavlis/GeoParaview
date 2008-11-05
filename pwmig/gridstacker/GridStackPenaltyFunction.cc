#include "perf.h"
#include "Metadata.h"
#include "GridStackPenaltyFunction.h"
GridStackPenaltyFunction::GridStackPenaltyFunction(Metadata& md)
{
   try {
	string pfuncname=md.get_string("PenaltyFunction");
	if(pfuncname=="simple_coherence")
	{
		gpt=COH;
		cohpow=1.0;
	}
	else if(pfuncname=="coherence_power_function")
	{
		gpt=COHPOW;
		cohpow=md.get_double("power_factor");
	}
	else if(pfuncname=="dbxcor_loss_function")
	{
		gpt=DBXCOR;
		cohpow=md.get_double("power_factor");
	}
	else
	{
		cerr << "GridStackPenaltyFunction constructor:  "
			<< "Unknown name for PenaltyFuntion="
			<< pfuncname<<endl;
		exit(-1);
	}
	
    }
    catch (MetadataGetError mderr)
    {
	mderr.log_error();
	exit(-1);
    }
}
double GridStackPenaltyFunction::weight(int nd, double *d, double *d0)
{
	double sumsqr,sumsqd,r;
	int i;
	double coh,wt;

	for(i=0,sumsqr=0.0;i<nd;++i) 
	{
		r=d[i]-d0[i];
		sumsqr += r*r;
	}
	for(i=0,sumsqd=0.0;i<nd;++i) 
	{
		sumsqd += (d[i])*(d[i]);
	}
	if(sumsqr>=sumsqd)
		coh=0.0;
	else
		coh=1.0-sqrt(sumsqr/sumsqd);
	switch(gpt)
	{
	case COH:
		wt=coh;
		break;
	case DBXCOR:
		wt=ddot(nd,d0,1,d,1);
		wt=wt/sqrt(sumsqd*sumsqr);
		wt=pow(wt,cohpow);
		break;
	case COHPOW:
	default:
		wt=pow(coh,cohpow);
		break;
	}
	return(wt);
}
