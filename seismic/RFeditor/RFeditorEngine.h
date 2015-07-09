#include "Metadata.h"
#include "TimeSeries.h"
#include "ensemble.h"
#include "TraceEditPlot.h"
using namespace std;
using namespace SEISPP;
const string evidkey("eventid");   // not evid to avoid collision
class RFeditorEngine
{
    public:
        RFeditorEngine(Metadata& params);
        /* \brief Loads ensembles of RFs and enable interactive editing.  

           Load a pair of ensembles for radial and transverse for 
           editing.   Edit procedure first makes the radial window 
           active.  Pick bad traces and exit.  Find matching transverse
           data and kills them before plotting the transverse data
           and making that window live for edits.

           \param r ensemble of radial receiver functions
           \param t ensemble of transverse receiver functions
           
           \return Number of traces loaded.
           */
        int edit(TimeSeriesEnsemble& r, TimeSeriesEnsemble& d);
        TimeWindow report_robust_twin(){return(robust_twin);};
        void set_robust_twin(TimeWindow tw){robust_twin=tw;};
    private:
        bool cutoff_editing;
        TraceEditPlot Rwindow;
        TraceEditPlot Twindow;
        TimeWindow robust_twin;
        /* private method kills T if R marked or R if T marked bad */
        void kill_sibling(); 
};
/* This is the sort order function used for stack sorting */
template <class T> struct less_stackwt
                : public binary_function<T,T,bool> 
{
        bool operator()(T x, T y)
        {
                double valx,valy;
                try {
                        valx=x.get_double(SEISPP::stack_weight_keyword);
                }catch(...) {return true;};
                try {
                        valy=y.get_double(SEISPP::stack_weight_keyword);
                }catch(...) {return false;};

                if(valx<valy)
                        return true;
                else
                        return false;
        }
};
