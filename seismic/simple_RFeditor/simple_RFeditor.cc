#include <string>
#include <fstream>
#include "seispp.h"
#include "SeisppKeywords.h"
#include "TimeWindow.h"
#include "dbpp.h"
#include "Metadata.h"
using namespace std;
using namespace SEISPP;
void usage()
{
    cerr << "simple_RFeditor dbin dbout [-cut cutoff] "<<endl
        << "Default cutoff is 6 events (12 RFs) per station"<<endl;
    exit(-1);
}
/* returns a list of keepers.  What it actually looks for is
   overlapping waveforms.   This algorithm assumes sort order
   is sta:time:chan grouped by sta and the DBBundle (grp) 
   defines rows for a single station referenced by 
   Dbptr db.   

   Algorithm also works only for frozen channel codes with EARS RF
   data (i.e. ITR and ITT)

   The output vector has all overlapping waveforms tagged false and 
   all other marked true.  This is appropriate for this application
   because overlapping events are big trouble for RF estimates.

   */
vector<bool> scan_for_dups(Dbptr db,DBBundle& grp)
{
    vector<bool> keepers;
    const string chan0("ITR"); // mod 2 divisible count
    const string chan1("ITT");
    char thischan[10];
    TimeWindow anchor,testwin;
    int i;
    bool this_ok,last_ok;
    for(i=0,db.record=grp.start_record;db.record<grp.end_record;
                        ++db.record,++i)
    {
        if(i<2)
        {
            dbgetv(db,0L,"time",&anchor.start,"endtime",&anchor.end,NULL);
            last_ok=true;
            this_ok=true;
        }
        else
        {
            dbgetv(db,0L,"chan",thischan,NULL);
            /* ITT should land here or it is an error */
            if(i%2)
            {
                if(string(thischan)!=chan1)
                {
                    cerr << "Data irregularity:   Found "
                        << thischan << " Expecting channel code "<<chan1<<endl;
                    cerr << "Marking this record bad.  Attempting to recover"
                        <<endl;
                    /* Incrementing i is a simple way to try to recover */
                    ++i;
                    this_ok=false;
                }
                /* Intentional do not even test time window for chan1 */
                else
                    this_ok=last_ok;

            }
            else
            {
                if(string(thischan)!=chan0)
                {
                    cerr << "Data irregularity:   Found "
                        << thischan << " Expecting channel code "<<chan0<<endl;
                    cerr << "Marking this record bad.  Attempting to recover"
                        <<endl;
                    this_ok=false;
                    /* Incrementing i is a simple way to try to recover */
                    ++i;
                }
                else
                {
                    dbgetv(db,0L,
                            "time",&testwin.start,
                            "endtime",&testwin.end,NULL);
                    if(testwin.start<anchor.end)
                    {
                        /* We need to back up an kill 2 previous 
                           entries.  This is not foolproof and this 
                           could cause irregularities if this follows
                           a channel irregularity */

                        keepers[i-2]=false;
                        keepers[i-1]=false;
                        this_ok=false;
                        last_ok=false;
                    }
                    else
                    {
                        anchor=testwin;
                        this_ok=true;
                        last_ok=true;
                    }
                }
            }
        }
        //DEBUG
        //cout << strtime(anchor.start) <<" is "<<this_ok<<endl;
        keepers.push_back(this_ok);
    }
    return(keepers);
}




bool SEISPP::SEISPP_verbose(false);
int main(int argc, char **argv)
{
    if(argc<3) usage();
    string dbin_name(argv[1]);
    string dbout_name(argv[2]);
    int cutoff(12);
    int i,j;
    for(i=3;i<argc;++i)
    {
        string sarg(argv[i]);
        /* quotes needed around subsets because sometimes
           station names start with numbers that confuse datascope */
        const string quote("\"");
        if(sarg=="-cut")
        {
            ++i;
            if(i>=argc) usage();
            cutoff=atoi(argv[i]);
            cutoff=2*cutoff;   // assume 2 RFs per station
        }
        else
            usage();
    }
    try {
        ofstream logfile;
        logfile.open("simple_RFeditor.log",ios::app);
        if(logfile.fail())
        {
            cerr << "Cannot open log file simple_RFeditor.log in append mode"
                <<endl << "Check permissions and try again"<<endl;
            exit(-1);
        }
        cout << "Starting simple_RFeditor on database "<<dbin_name<<endl;
        logfile << "simple_RFeditor run on " << dbin_name << " at time "
            << strtime(now()) <<endl;
        AttributeMap am("css3.0");
        /* Open the in and out database */
        DatascopeHandle dbin(dbin_name,true);
        dbin.lookup("wfdisc");
        if(dbin.number_tuples()<=0)
        {
            cerr << "No rows in input wfdisc for database "<<dbin_name<<endl;
            exit(-1);
        }
        logfile << "Editing data from database "<<dbin_name<<endl
            <<"Number of rows in full wfdisc originally openned="
            <<dbin.number_tuples()<<endl;
        DatascopeHandle dbout(dbout_name,false);
        /* First check that the output wfdisc is empty */
        dbout.lookup("wfdisc");
        logfile << "Writing results to "<<dbin_name<<".wfdisc"<<endl
            <<"Number of existing rows in the input wfdisc table is "
            <<dbout.number_tuples()<<endl;
        list<string> sortkeys, groupkeys;
        sortkeys.push_back("sta");
        sortkeys.push_back("time");
        sortkeys.push_back("chan");
        groupkeys.push_back("sta");
        dbin.sort(sortkeys);
        dbin.group(groupkeys);
        cout << "Number of ensembles to process (grouped by sta)="
            << dbin.number_tuples()<<endl;
        logfile << "Number of ensembles to process (grouped by sta)="
            << dbin.number_tuples()<<endl;
        long nsta=dbin.number_tuples();
        /* This is used to delete overlapping waveforms.  Because of 
           sort order this will always delete later version of overlaps. */
        vector<bool> keepers;
        for(i=0,dbin.rewind();i<nsta;++dbin,++i)
        {
            DBBundle grp=dbin.get_range();
            char rec[512];
            long reccount=grp.end_record-grp.start_record;
            Dbptr db;
            db=grp.parent;
            char sta[10];
            dbgetv(db,0L,"sta",sta,NULL);
            cout << "Station "<<sta<<":  ";
            logfile << "Station "<<sta<<":  ";
            if(reccount<=cutoff)
            {
                cout << "Dropped.   Record count="<<reccount<<endl;
                logfile << "Dropped.   Record count="<<reccount<<endl;
            }
            else
            {
                cout << "Writing "<<reccount<<" wfdisc rows"<<endl;
                logfile << "Writing "<<reccount<<" wfdisc rows"<<endl;
                db=grp.parent;
                keepers=scan_for_dups(db,grp);
                for(j=0,db.record=grp.start_record;db.record<grp.end_record;
                        ++db.record,++j)
                {
                    /* The list of attributes is not complete here.  
                       This is restricted to wfdisc entries known to 
                       be set by EARS data.   Do not carry forward 
                       blindly ever. */
                    char chan[10],datatype[4];
                    char dir[70],dfile[40];
                    double time,endtime,samprate;
                    long wfid,jdate,nsamp;
                    dbgetv(db,0L,"sta",sta,"chan",chan,"time",&time,
                            "endtime",&endtime,"wfid",&wfid,
                            "samprate",&samprate,"nsamp",&nsamp,
                            "jdate",&jdate,"dir",dir,"dfile",dfile,
                            "datatype",datatype,NULL);
                    if(keepers[j])
                    {
                        dbaddv(dbout.db,0L,"sta",sta,"chan",chan,"time",time,
                            "endtime",endtime,"wfid",wfid,
                            "samprate",samprate,"nsamp",nsamp,
                            "jdate",jdate,"dir",dir,"dfile",dfile,
                            "datatype",datatype,NULL);
                    }
                    else
                        cout << "Skipping overlapping wavefor  for "
                            << sta << "  with waveform segment "
                            << strtime(time)<<"->"<<strtime(endtime)<<endl;
                }
            }
        }
    }catch(exception& err)
    {
        cerr << "Fatal error:  error message follows."<<endl
            << err.what()<<endl;
    }
}
