#ifndef _SEISMICPLOT_H_
#define _SEISMICPLOT_H_
#include <pthread.h>
#include "ensemble.h"
#include "Seisw.h"
#include "SeismicPick.h"
#include "BasicSeisPlot.h"
#include "Metadata.h"
#include <Xm/MainW.h>
namespace SEISPP
{
using namespace std;
using namespace SEISPP;
/* This defines a standardized interface to a plot for the SEISPP library.

   This class is an attempt to abstract the basic plot functions required
   to plot and work with seismic data.  It inherits an abstract base class
   called BasicSeismicPlot which is a pure virtual base.  It also inherits
   Metadata as a convenient way to change plot parameters.  An anomaly that
   makes this a nonideal interface is that in the current implementation
   this requires the refresh method to force the new parameters to take
   effect.  It would be better if any parameter changes caused an immediate
   redraw, but the current version based on Xt and the Seisw Motif widget
   just can't do that.  */
class SeismicPlot : public BasicSeisPlot, public Metadata
{
    public:
        SeismicPlot();
        SeismicPlot(Metadata& md);
        void plot(TimeSeriesEnsemble& d);
        void plot(ThreeComponentEnsemble& d);
        void plot(TimeSeries& d);
        void plot(ThreeComponentSeismogram& d);
        double PickTime();
        TimeWindow PickWindow();
        void refresh();
    private:
        /* These are used to be safe.  They may not be necessary, but 
           because normally these are initialized from main this fakery
           may be needed */
        int argc;  
        char *argv[1];
        XtAppContext AppContext;
        Widget toplevel,main_w;
        Widget rshell[2];  // pop up windows for 3c data comp 2 and 3
        /* To allow three componnts we have up to three seismic display 
           widgets in three windows.  Only the first is created by
           constructors. Others are created only when a call to plot
           a 3c ensmble is called.  This is controlled by the nwindows
           integer.  Set initially to 1, then incremented to 3 when 
           a 3c ensemble is to be plotted.
         */
        int nwindows;
        Widget seisw[3];   // This is the seismic display widget
        /*! This is the id number returned by pthread_create for X event
          handler thread */
        int event_handler_thread_id;
        /* Private method that launches event handler thread */
        void launch_Xevent_thread_handler();
};
}
extern "C" {
/*! Routine run by event handler thread */
void *Xevent_loop(void *);
}
#endif
