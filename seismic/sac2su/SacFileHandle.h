#include "GenericFileHandle.h"
using namespace std;
using namespace SEISPP;
/* \brief Special GenericFileHandle for SAC files.

   This handle can be used to read a SAC file.  SAC has no
   concept of ensembles so it has only a method to read 
   a single seismogram.   At this time it also intentionally
   has not write method as I refuse to perpetuate such an 
   abomination.  There are also plenty of other tools to 
   write sac files. 
   */

class SacFileHandle : public GenericFileHandle
{
public:
    SacFileHandle(string filename,AttributeCrossReference& namemap,
                list<string> mdlist,list<string> ensk, list<string>ensmdl);
    TimeSeries GetNextSeismogram();
};
