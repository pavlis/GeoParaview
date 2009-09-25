#ifndef _BASICSEISPLOT_H
#define _BASICSEISPLOT_H
#include "ensemble.h"
namespace SEISPP
{
using namespace std;
using namespace SEISPP;

class BasicSeisPlot
{
    public:
        virtual void plot(TimeSeriesEnsemble& d)=0;
        virtual void plot(ThreeComponentEnsemble& d)=0;
        virtual void plot(TimeSeries& d)=0;
        virtual void plot(ThreeComponentSeismogram& d)=0;
        virtual double PickTime()=0;
        virtual TimeWindow PickWindow()=0;
};
}
#endif
