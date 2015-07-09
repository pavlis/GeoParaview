/* Use once program to fix data problem in STEEP teleseismic processing 
by Lindsay Kenyon created by a bug in dbxcor.  dbxcor was storing 
residual as an absolute time in arrival.  Thus all this program does is
associate each arrival with it's original arid (match on orid field of
arrival) and add the predicted time from iaspei91.  New arrival table is 
written to output. Assoc requires not changes.

Model is to verify this is correct and have the new arrival table replace
the original.  
*/
#include <list>
#include "coords.h"
#include "dbpp.h"
#include "Hypocenter.h"
using namespace std;
using namespace SEISPP;
bool SEISPP::SEISPP_verbose(true);
void usage()
{
    cerr << "fixresiduals dbin dbout"<<endl;
    exit(-1);
}
int main(int argc,char **argv)
{
    if(argc!=3) usage();
    try {
        DatascopeHandle dbin(string(argv[1]),false);
        DatascopeHandle dbout(string(argv[2]),false);
        dbout.lookup("arrival");
        AttributeMap am(string("css3.0"));
        list<string> matchlist;
        matchlist.push_back(string("orid"));
        DatascopeMatchHandle dborigin(dbin,string("origin"),matchlist,am);
        DatascopeHandle dbarrival(dbin);
        dbarrival.lookup("arrival");
        dbarrival.natural_join("assoc");
        cout << "Size of arrival table="<<dbarrival.number_tuples()<<endl;
        list<string> jkey;
        jkey.push_back(string("sta"));
        dbarrival.join(string("site"),jkey,jkey);
        dbarrival.rewind();  //not really necessary
        cout << "Size after join with site using sta key="<<dbarrival.number_tuples()<<endl;
        /* We will use the raw pointer to add.  dbpp method is 
           very awkward in this case */
        Dbptr dboa;
        /* Loop through arrival table.  Use a very inefficient 
           algorithm to create a Hypocenter object row by row, but
           this is a run once program on a small database so not important*/
        Metadata m;
        dborigin.lookup("origin");
        Dbptr dbmo=dborigin.db;
        long i;
        for(i=0;i<dbarrival.number_tuples();++i,++dbarrival)
        {
            double lat,lon,depth,otime;
            long orid=dbarrival.get_long("orid");
            m.put("orid",orid);
            list<long> recs=dborigin.find(m);
            if(recs.size()==0)
                cerr << "No match for orid="<<orid<<endl;
            else if(recs.size()>1)
                cerr << "Warning multiple orid matches for orid="<<orid
                    <<endl
                    <<"This should not happen since orid must be unique in origin"
                    <<endl
                    <<"Exiting.  Fix the code or the db"<<endl;
            else
            {
                long dbgret;  //used for debug  intentionally ignored as throwaway program
                dbmo.record=*(recs.begin());
                dbgret=dbgetv(dbmo,0,"lat",&lat,
                        "lon",&lon,
                        "depth",&depth,
                        "origin.time",&otime,NULL);
                Hypocenter h(rad(lat),rad(lon),depth,otime,
                        string("tttaup"),string("iasp91"));
                lat=dbarrival.get_double("lat");
                lon=dbarrival.get_double("lon");
                double resid=dbarrival.get_double("arrival.time");
                double elev(0.0); // Intentionally ignore elevation
                /* All data being fixed a p phases so use the simple
                   interface to ptime */
                double predtime=h.ptime(rad(lat),rad(lon),elev); 
                double atime=h.time+predtime+resid;
                /* To complete later.  Load all set attributes
                   and copy to output */
                char sta[10],chan[10],iphase[10],auth[10];
                long arid,jdate;
                double snr,deltim;
                dbgret=dbgetv(dbarrival.db,0,"sta",sta,"chan",chan,
                        "iphase",iphase,"auth",auth,
                        "arid",&arid,"jdate",&jdate,
                        "snr",&snr,"deltim",&deltim,NULL);
                Dbptr dbassoc;
                dbgret=dbgetv(dbarrival.db,0,"assoc",&dbassoc,NULL);
                dbgret=dbputv(dbassoc,0,"timeres",resid,NULL);
                dbgret=dbaddv(dbout.db,0,"sta",sta,"chan",chan,
                        "time",atime,
                        "iphase",iphase,"auth",auth,
                        "arid",arid,"jdate",jdate,
                        "snr",snr,"deltim",deltim,NULL);
            }
        }
    }
    catch (SeisppError& serr)
    {
        serr.log_error();
    }
    catch (...)
    {
        cerr << "Something threw and unexpected exception"<<endl;
    }
}


            

