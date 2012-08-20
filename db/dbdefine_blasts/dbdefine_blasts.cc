/* This program provides a simple way to force marking events in 
   a specified time window as explosions.  The algorithm is simple.
   The hour field of the origin time is compared to the internally 
   defined range and if between the hours defined the origin attribute
   etype is set ex.   Note this is inclusive (or equal comparison).   
   */
#include <iostream>
#include "dbpp.h"
using namespace std;
using namespace SEISPP;
void usage()
{
    cerr << "dbdefine_blasts db [-s starttime]"<<endl
        <<"    -s option ignores events before starttime"<<endl;;
    exit(-1);
}
bool SEISPP::SEISPP_verbose(true);
int main(int argc, char **argv)
{
    /* Events with hours in this range are those marked automatically ex*/
    const double starthr(15), endhr(23);;
    if(argc<2) usage();
    string dbname(argv[1]);
    int i;
    double starttime(0);
    for(i=2;i<argc;++i)
    {
        string sarg(argv[i]);
        if(sarg=="-s")
        {
            ++i;
            if(i>=argc)usage();
            starttime=str2epoch(argv[i]);
            cout << "dbdefine_blasts considering only events after time "
                << strtime(starttime)<<endl;
        }
        else
            usage();
    }
    try {
        DatascopeHandle dbh(dbname,false);
        dbh.lookup("origin");
        long nrows=dbh.number_tuples();
        if(nrows<=0)
        {
            cerr <<"Origin table for database "
                << dbname <<" is empty"<<endl;
            usage();
        }
        dbh.rewind();
        long rec,nex_set;
        for(rec=0,nex_set=0;rec<nrows;++rec,++dbh)
        {
            double otime;
            otime=dbh.get_double("time");
            char *hour;
            hour=epoch2str(otime,"%H");
            int hr=atoi(hour);
            if(otime>starttime && hr>=starthr && hr<=endhr)
            {
                dbh.put("etype","ex");
                ++nex_set;
            }
            free(hour);
        }
        cout << "dbdefine_blasts set etype field of "
            << nex_set<<" origin rows to explosion tag=ex"<<endl;
    }catch(SeisppError &serr)
    {
        serr.log_error();
        exit(-1);
    }
}
