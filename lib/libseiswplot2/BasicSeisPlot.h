#ifndef _BASICSEISPLOT_H
#define _BASICSEISPLOT_H
#include "ensemble.h"
namespace SEISPP
{
using namespace std;
using namespace SEISPP;
/*! Abstract base class for family of seismic plotting objects.

  This provides a pure virtual base for seismic plotting centered
  around the SEISPP library of data objects.  
  
  \author Gary L Pavlis
  */

class BasicSeisPlot
{
    public:
        virtual void plot(TimeSeriesEnsemble& d)=0;
        virtual void plot(ThreeComponentEnsemble& d)=0;
        virtual void plot(TimeSeries& d)=0;
        virtual void plot(ThreeComponentSeismogram& d)=0;
};
}
#endif
