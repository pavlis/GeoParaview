#include "SeismicPick.h"
SeismicPick::SeismicPick() : BasicPick()
{
	time=0.0; 
	amplitude=0.0;
	trace_number=-1;
	amplitude_scale=1.0;
	time_scale=1.0;
}

SeismicPick::SeismicPick(const SeismicPick& p) : BasicPick(p)
{
	time=p.time;
	amplitude=p.amplitude;
	trace_number=p.trace_number;
	amplitude_scale=p.amplitude_scale;
}
SeismicPick& SeismicPick::operator=(const SeismicPick& p)
{
	if(this!=&p)
	{
		time=p.time;
		amplitude=p.amplitude;
		trace_number=p.trace_number;
		amplitude_scale=p.amplitude_scale;
	}
	return(*this);
}

