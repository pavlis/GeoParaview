#include "stock.h"
#include "Metadata.h"
#include "TraceEditPlot.h"
bool SEISPP::SEISPP_verbose(true);
int main(int argc, char **argv)
{
    try {
    Pf *pf;
    int iret=pfread("testplotter",&pf);
    if(iret!=0)
    {
        cerr << "pfread error"<<endl;
        exit(-1);
    }
    Metadata control(pf);
    cout << "Read pf file.  Calling plotter\n"<<endl;
    TraceEditPlot win(control);
    cout << "Creating fake data ensemble"<<endl;
    TimeSeriesEnsemble tse(10,1000);
    for(int i=0;i<10;++i)
    {
        int j;
        tse.member[i].s.clear();
        for(j=0;j<1000;++j) tse.member[i].s.push_back(0.0);
        tse.member[i].ns=1000;
        tse.member[i].dt=0.1;
        for(j=0;j<25;++j)
        {
            tse.member[i].s[100+2*i+j]=(double)(j+1);
        }
        tse.member[i].live=true;
        //cout << "member "<<i<<endl<<tse.member[i]<<endl;
    }
    cout << tse;
    cout << "Calling plot for fake ensemble"<<endl;
    win.plot(tse);
    cout << "Success.... exiting"<<endl;
    }catch(SeisppError& serr)
    {
        serr.log_error();
    }
}
