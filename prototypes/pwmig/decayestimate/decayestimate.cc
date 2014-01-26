/* This small program reads an output image from pwmig and puts out 
   average amplitudes as a function of depth.  Output is a pair of
   depth, amplitude pairs.  The "average" is literally a median.
   Actually think I'll output depth, median, q1/4, q3/4 since these
   come for free through use of the VectorStatistics template */
#include <string>
#include <vector>
#include "dbpp.h"
#include "gclgrid.h"
#include "VectorStatistics.h"
void usage()
{
    cerr << "decayestimate db grid field"<<endl;
    exit(-1);
}
bool SEISPP::SEISPP_verbose(true);
int main(int argc, char **argv)
{
    /* This may need to be a parameter but for use with pwmig
       and current scaling this should be workable to test for a zero value*/
    const double zerotest(1.0e-12);
    if(argc!=4) usage();
    string dbname(argv[1]);
    string gridname(argv[2]);
    string fieldname(argv[3]);
    int i,j,k,l;
    try {
        DatascopeHandle dbh(dbname,true);
        GCLvectorfield3d f(dbh,gridname,fieldname);
        vector<double> rawamp;
        /*count backward because grid top is at n3-1*/
        for(k=f.n3-1;k>0;--k)
        {
            double amp;
            rawamp.clear();
            for(i=0;i<f.n1;++i)
            {
                for(j=0;j<f.n2;++j)
                {
                    amp=0.0;
                    for(l=0;l<3;++l) 
                        amp+=(f.val[i][j][k][l])*(f.val[i][j][k][l]);
                    amp=sqrt(amp);
                    if(amp>zerotest)rawamp.push_back(amp);
                }
            }
            /* Don't compute statistics if count is below 4 */
            if(rawamp.size()>4)
            {
                VectorStatistics<double> astats(rawamp);
                /* assumes a standard image grid with top surface
                   at depth of 0 */
                double depth=f.depth(0,0,k);
                cout << depth <<" "
                    << astats.median() << " "
                    << astats.q1_4() << " "
                    << astats.q3_4() << " "
                    << astats.lower_bound() << " "
                    << astats.upper_bound() << " "
                    << rawamp.size()
                    <<endl;;
            }
        }
    }
    catch(SeisppError& serr)
    {
        serr.log_error();
    }
    catch(std::exception& err)
    {
        cerr << err.what();
    }
}





