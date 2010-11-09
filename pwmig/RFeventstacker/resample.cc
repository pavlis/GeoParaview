#include <memory>
#include "ThreeComponentSeismogram.h"
#include "ensemble.h"
#include "resample.h"
using namespace std;
using namespace SEISPP;
ThreeComponentSeismogram Resample(ThreeComponentSeismogram& d,
        ResamplingDefinitions& rd, double dtout, bool trim)
{
    TimeSeries *xi;
    try{
        vector<TimeSeries> x;
        int i;
        for(i=0;i<3;++i)
        {
            xi=ExtractComponent(d,i);
            x.push_back(ResampleTimeSeries(*xi,rd,dtout,trim));
            delete xi;
        }
        ThreeComponentSeismogram result(x,0);
        return(result);
    } catch(...){throw;};
}
auto_ptr<ThreeComponentEnsemble> Resample(ThreeComponentEnsemble& d,
        ResamplingDefinitions& rd, double dtout, bool trim)
{
    try{
        /* We call this constructor for efficiency.  It clones the 
           gather metadata and reserves slots for the resampled data.*/
        auto_ptr<ThreeComponentEnsemble> ddec;
        ddec = auto_ptr<ThreeComponentEnsemble>(new ThreeComponentEnsemble(
                        dynamic_cast<Metadata&>(d),d.member.size()));
        vector<ThreeComponentSeismogram>::iterator dptr;
        for(dptr=d.member.begin();dptr!=d.member.end();++dptr)
        {
            /* For efficiency avoid calling resampler if sample intervals match.*/
            if((dptr->dt)==dtout)
                ddec->member.push_back(*dptr);
            else
            {            
                ThreeComponentSeismogram dmd(Resample(*dptr,rd,dtout,trim));
                ddec->member.push_back(dmd);
            }
        }
        return ddec;
    } catch(...){throw;};
}
