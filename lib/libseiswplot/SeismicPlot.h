#ifndef _SEISMICPLOT_H_
#define _SEISMICPLOT_H_
#include <pthread.h>
#include "ensemble.h"
#include "Seisw.h"
#include "BasicSeisPlot.h"
#include "Metadata.h"
#include <Xm/MainW.h>
#include <Xm/Form.h>
#include <Xm/PushB.h>
namespace SEISPP
{
using namespace std;
using namespace SEISPP;
/*! This defines a standardized interface to a plot for the SEISPP library.

   This class is an attempt to abstract the basic plot functions required
   to plot and work with seismic data.  It inherits an abstract base class
   called BasicSeismicPlot which is a pure virtual base.  It also inherits
   Metadata as a convenient way to change plot parameters.  An anomaly that
   makes this a nonideal interface is that in the current implementation
   this requires the refresh method to force the new parameters to take
   effect.  It would be better if any parameter changes caused an immediate
   redraw, but the current version based on Xt and the Seisw Motif widget
   just can't do that.  

   Recognize that this implementation is designed only to plot data, not interact with
   it.  If you need that functionality you will need to either write your own GUI 
   or extend this one.  */
class SeismicPlot : public BasicSeisPlot, public Metadata
{
    public:
        /*! Construct a plot handle using default parameters. 

          Defaults for this object are stored in an Antelope a parameter file calld
          SeismicPlot.pf stored in $ANTELOPE/data/pf.  This constructor uses pure
          defaults for parameters. */
        SeismicPlot();
        /* Constructor using Metadata object to define plot parameters.

           The plot geometry is defined by a large group of parameters passed to the
           plotting widget through a Metadata object.  
           \param md is the Metadata object that is assumed to contain the (fairly large) 
                number of parameters needed to define the plot geometry. 
            \exception Will throw a SeisppError object if plot parameters are missing.
        */
        SeismicPlot(Metadata& md);
        /*! Plot an ensemble of data.

          This method plots an ensemble of scalar data in a single window frame. 

          \param d is the ensemble to be plotted. 
          */
        void plot(TimeSeriesEnsemble& d);
        /*! Plot a three component ensemble of data.

          Three-component data could be plotted a variety of ways, but here I choose to plot the 
          data in three parallel window.  Each component is displayed in a window with a label for 
          (C order ) component 0, 1, and 2. 

          \param d is data to be plotted.
          */
        void plot(ThreeComponentEnsemble& d);
        /*! Plot a single time series.

          This is essentially a special case of a TimeSeriesEnsemble with one member and it is
          treated that way.  i.e. you get a plot of the data with one trace in the window.
          */
        void plot(TimeSeries& d);
        /*! Plot a single three-component seismogram.

          A three component seismogram is plotted in the window 0 (the one labelled component 0) 
          on a common time base.  The alternative would have been to display each component in a 
          different window, which would not make sense for common use for this type of plot.  

          \param d is the seismogram to be plotted.
          */
        void plot(ThreeComponentSeismogram& d);
        /*! Make new plot parameters active.

          This plot object inherits another a Metadata object.  You can thus use the 
          get and put methods in Metadata to change parameters internally to this object's
          internal data definitions.  The design of the Widget used to implement this 
          class requires an explicit push of the entire Metadata object through an 
          Xwindows call.  Rather than override the get and put methods I judged it 
          preferable to use this model.  That is, if you want to change an plot
          parameters use the Metadata put methods to change what you want.  When you
          have put all the new parameters you want call this (refresh) method to 
          have them new parameters pushed to the Widget.  
          */
        void refresh();
        /*! Used internally to break the event loop.  

          This really shouldn't be int he public interface, but I couldn't figure out
          a way to get rid of it. */
        void ExitDisplay(){
            EventLoopIsActive=false; 
        };
        /*! Used internally by event loop thread to test for termination.

          his really shouldn't be int he public interface, but I couldn't figure out
                    a way to get rid of it. */
        bool WindowIsActive(){
            return EventLoopIsActive;
        };
        /*! X application context.

          This plotting is based on X and I could not figure out a way to not put
          this in the public interface and get threading to work right.  It really
          doesn't belong here and a caller should never ever touch this. */
        XtAppContext AppContext;
    private:
        bool EventLoopIsActive;  // Set true when an event loop thread is running
        /* These are used to be safe.  They may not be necessary, but 
           because normally these are initialized from main this fakery
           may be needed */
        int argc;  
        char *argv[1];
        Widget toplevel,main_w;
        Widget continue_button;
        Widget rshell[3];  // pop up windows for 3c data comp 2 and 3
        Widget seisw[3];   // This is the seismic display widget
        /* Private method that launches event handler thread */
        void launch_Xevent_thread_handler();
};
}
#endif
