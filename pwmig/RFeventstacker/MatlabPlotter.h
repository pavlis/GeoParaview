#ifndef _MATLAB_PLOTTER_H_
#define _MATLAB_PLOTTER_H_

#include "MatlabProcessor.h"
class MatlabPlotter : public MatlabProcessor
{
public:
	/*! Primary constructor.

	Construction in this case is mostly calling the default constructor for the MatlabProcessor
	which is the engine for this derived class.  MatlabProcessor has other options, but this
	constructor assumes matlab is local on this host and the only interaction will be through
	this object. */
	MatlabPlotter();
	/*! Make a wiggle trace shaded area plot of all three components of a 3c gather.

	This method plots three components using 3 windows (one for each component) 
	using a wiggle trace variable area plot and an optional image overlay.
	Multiple calls to this procedure will overwrite each window.  
	\param d is data to plot
	\param channame is a list of channel names to be posted on each plot by matlab title command.
		These are the only hint given as to which window is which component
	\param plot_image if true imagesc is used to superimpose a image background on the wtva plots.

	\exception SeisppError may be thrown in one of several instances.  
i	This method assumes matlab can find an m file called "xwigb.m" in it's search path
	for m-files.  If this is not found you will get an error.  Further if component number
	is invalid you will get an exception.
	*/  
	void wigbplot(ThreeComponentEnsemble& d,int component,string channame);
	/*! Plot a scalar ensemble with a wiggle trace variable area plot.

	This method plots scalar, TimeSeriesEnsemble data with a wiggle trace variable area 
	display.   
	void wigbplot(TimeSeriesEnsemble& d);
	void plot(ThreeComponentSeismogram& d);
	void plot(TimeSeries& d);
private:
	const int fn3ce[3](1,2,3);
	const int fntse(4);
	const int fnts(5);
	const int fntcs(6);
};
#endif
