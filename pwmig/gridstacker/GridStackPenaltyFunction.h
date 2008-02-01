#include "Metadata.h"
using namespace std;
using namespace SEISPP;
enum GridPenaltyType {COH, COHPOW, DBXCOR};
/* Interface to a penalty function to compute a weight grid data.
This is a very simple interface locked to this program.  The only real
purpose is to provide a way to easily add alternative penalty functions
without changing the rest of the code. */
class GridStackPenaltyFunction 
{
public:
        /* Constructor from metadata */
	GridStackPenaltyFunction(Metadata& def);
        /* compute a weight for data d using d0 as an estimate of center (stack). 
        This simple method is sufficient here or now since we use contiguous memory
        arrays of equal length stored in the scratch file.  May need a generalzation
        to a pair of field variables, but this will be far more efficient.*/
        double weight(int nd, double *d, double *d0);
private:
	GridPenaltyType gpt;
	double cohpow;
};
