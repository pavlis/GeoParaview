#include <list>
#include <string>
#include "stock.h"
#include "dbpp.h"
using namespace std;
using namespace SEISPP;
void usage()
{
    cerr << "db2artm db [-s station_subset_string" <<endl
        << "  Writes rpi tomo arrival data file to stdout"<<endl
        << "  station_subset_string is passed directly to dbsubset"<<endl;;
    exit(-1);
}
bool SEISPP::SEISPP_verbose(false);
int main(int argc, char **argv)
{
    if(argc<2) usage();
    string dbname(argv[1]);
    string ss_string;
    bool subset_enabled(false);
    int i;
    for(i=2;i<argc;++i)
    {
        string sarg(argv[i]);
        if(sarg=="-s")
        {
            ++i;
            if(i>=argc) usage();
            ss_string=argv[i];
            subset_enabled=true;
        }
        else
            usage();
    }
    try {
        ios::sync_with_stdio();
        cerr << "Starting db2artm program "<<endl;
        DatascopeHandle dbh(dbname,true);
        dbh.lookup(string("event"));
        dbh.natural_join(string("origin"));
        dbh.subset(string("orid==prefor"));
        dbh.natural_join(string("assoc"));
        dbh.natural_join(string("arrival"));
        dbh.natural_join(string("site"));
        if(subset_enabled)
        {
            cerr << "Subsetting "<<dbname<<" using condition="
                <<ss_string<<endl
                << "input database before subset number of rows="
                << dbh.number_tuples()<<endl;
            dbh.subset(ss_string);
            cerr<< "input database after subset number of rows="
                << dbh.number_tuples()<<endl;
        }
        list<string> sortkeys;
        sortkeys.push_back(string("evid"));
        cerr << "input database view number of rows="
                    <<dbh.number_tuples()<<endl;
        dbh.group(sortkeys);
        int nevents=dbh.number_tuples();
        cerr << "Number of event groups for output="<<nevents<<endl;
        cerr << "Starting to write output"<<endl;
        long iev;
        for(dbh.rewind(),iev=0;iev<nevents;++iev,++dbh)
        {
            DBBundle grp=dbh.get_range();
            Dbptr db=grp.parent;
            db.record=grp.start_record;
            // Get the hypocenter from the first record - in all rows of group
            double lat,lon,depth,time;
            dbgetv(db,0,"origin.time",&time,
                    "lat",&lat,
                    "lon",&lon,
                    "depth",&depth,NULL);
            char *timestr;
            string tform("%Y %j %H %M %S.%s");
            timestr=epoch2str(time,const_cast<char *>(tform.c_str()));
            fprintf(stdout,"%s %10.5lf %10.5lf%10.5lf\n",
                    timestr,lat,lon,depth);
            free(timestr);
            string atform("%Y %j %H %M %S.%s");
            long irec;
            for(irec=grp.start_record;irec<grp.end_record;++irec)
            {
                double atime,deltim;
                char iphase[8]; // perhaps should use assoc.phase 
                char sta[10];
                db.record=irec;
                dbgetv(db,0,"sta",sta,
                        "arrival.time",&atime,
                        "deltim",&deltim,
                        "iphase",iphase,NULL);
                timestr=epoch2str(atime,const_cast<char *>(atform.c_str()));
                fprintf(stdout,"%-6s%s %-2s%15.3lf\n",
                        sta,timestr,iphase,deltim);
                free(timestr);
            }
            fprintf(stdout,"\n");  // need blank line sentinel
        }
    } catch (SeisppError& serr)
    {
        serr.log_error();
    }
}



