#ifndef SEISMICPICK_H
#define SEISMICPICK_H

class BasicPick
{
public:
	int screenX,screenY;
	BasicPick(){screenX=-1; screenY=-1;};
	BasicPick(int x,int y){screenX=x;screenY=y;};
};

class SeismicPick : BasicPick
{
public:
	double time;
	double amplitude;
	int trace_number;
	double amplitude_scale;
	double time_scale;
	SeismicPick();
	SeismicPick(const SeismicPick& p);
	SeismicPick& operator=(const SeismicPick& p);
};

#endif
	
