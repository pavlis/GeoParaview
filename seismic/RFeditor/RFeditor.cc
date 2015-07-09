#include <string>
#include <fstream>
#include "seispp.h"
#include "SeisppKeywords.h"
#include "dbpp.h"
#include "Metadata.h"
#include "RFeditorEngine.h"
using namespace std;
using namespace SEISPP;
void usage()
{
    cerr << "RFeditor dbin dbout [-pf pffile] [-laststa xx || -ss subset_condition]"<<endl;
    exit(-1);
}
/* Simple algorithm sets an int metadata item (evidkey defined in 
   RFeditorEngine) by a simple counts starting at 1.  
   The algorith is very simple.  It assumes the data were already sorted
   by sta:time:chan.  The approach then is to work through the ensemble
   and whenever t0 of successive traces are within one sample dt group 
   them together.   
   */
void set_eventids(TimeSeriesEnsemble& d)
{
    vector<TimeSeries>::iterator dptr;
    int i,evid;
    double t0test(0.0);
    /* Start with evid 0 to match vector index when not sorted.
       confusing otherwise */
    for(evid=0,i=0,dptr=d.member.begin();dptr!=d.member.end();
            ++i,++dptr)
    {
        /*
        cout << "Member number "<<i<< "Metadata contens"<<endl
            << dynamic_cast<Metadata &>(*dptr)<<endl;
            */
        if(i==0)
        {
            dptr->put(evidkey,evid);
            t0test=dptr->t0;
        }
        else
        {
            /* elegantly obscure C syntax - test if t0 of this member
               is with dt 0 test value */
            if(fabs(t0test-(dptr->t0))<dptr->dt)
                dptr->put(evidkey,evid);
            else
            {
                ++evid;
                t0test=dptr->t0;
                dptr->put(evidkey,evid);
            }
        }
        //DEBUG
        cout << "Member number "<<i<<" assigned evid "<<evid<<endl;
    }
}
TimeSeriesEnsemble extract_by_chan(TimeSeriesEnsemble& d,string chankey)
{
    try {
        int nd=d.member.size();
        int ns0=d.member[0].s.size();
        if(nd==0) throw SeisppError(string("extract_by_chan:  ")
                + "input ensemble has no data");
        /* Inefficient to copy, but this way we get metadata copied.  Use
           until proven to be a barrier */
        TimeSeriesEnsemble result(d);
        result.member.clear();
        vector<TimeSeries>::iterator dptr;
        for(dptr=d.member.begin();dptr!=d.member.end();++dptr)
        {
            string chan=dptr->get_string("chan");
            size_t pos=chan.find(chankey);
            if(pos!=chan.npos)
                result.member.push_back(*dptr);
        }
        return(result);
    }catch(...){throw;};
}
int save_to_db(TimeSeriesEnsemble& r, TimeSeriesEnsemble& t,
        MetadataList& mdl, AttributeMap& am, DatascopeHandle& dbh)
{
    try {
        vector<TimeSeries>::iterator dptr;
        int nlive(0);
        const string table("wfdisc");
        for(dptr=r.member.begin();dptr!=r.member.end();++dptr)
        {
            string sta=dptr->get_string("sta");
            if(dptr->live && (sta!="stack"))
            {
                int rec;
                rec=dbsave(*dptr,dbh.db,table,mdl,am);
//DEBUG
cerr << "Sta="<<dptr->get_string("sta")<<endl;
                ++nlive;
            }
        }
//DEBUG
cerr << "Number radial saved="<<nlive<<endl;
        for(dptr=t.member.begin();dptr!=t.member.end();++dptr)
        {
            if(dptr->live)
            {
                int rec;
                rec=dbsave(*dptr,dbh.db,table,mdl,am);
                ++nlive;
            }
        }
        return(nlive);
    }catch(...){throw;};
}
bool SEISPP::SEISPP_verbose(false);
int main(int argc, char **argv)
{
    if(argc<3) usage();
    string dbin_name(argv[1]);
    string dbout_name(argv[2]);
    string pfname("RFeditor");
    bool apply_subset(false);
    string subset_condition;
    int i;
    for(i=3;i<argc;++i)
    {
        string sarg(argv[i]);
        /* quotes needed around subsets because sometimes
           station names start with numbers that confuse datascope */
        const string quote("\"");
        if(sarg=="-pf")
        {
            ++i;
            if(i>=argc) usage();
            pfname=string(argv[i]);
        }
        else if(sarg=="-ss")
        {
            ++i;
            if(i>=argc) usage();
            apply_subset=true;
            subset_condition=string(argv[i]);
        }
        else if(sarg=="-laststa")
        {
            ++i;
            if(i>=argc) usage();
            apply_subset=true;
            subset_condition=string("sta > ")+quote
                +string(argv[i])+quote;
        }
        else
            usage();
    }
    Pf *pf;
    if(pfread(const_cast<char *>(pfname.c_str()),&pf))
    {
        cerr << "pfread failed on pf file named "<<pfname<<endl;
        exit(-1);
    }
    try {
        Metadata control(pf);
        MetadataList mdlens=pfget_mdlist(pf,"ensemble_mdl");
        MetadataList mdl=pfget_mdlist(pf,"trace_mdl");
        MetadataList mdlout=pfget_mdlist(pf,"output_mdl");
        string rchan=control.get_string("radial_channel_key");
        string tchan=control.get_string("transverse_channel_key");
        int minrfcutoff=control.get_int("minimum_number_receiver_functions");
        /* Open a log file in append mode */
        ofstream logfile;
        logfile.open("RFeditor.log",ios::app);
        if(logfile.fail())
        {
            cerr << "Cannot open log file RFeditor.log in append mode"
                <<endl << "Check permissions and try again"<<endl;
            exit(-1);
        }
        cout << "Starting RFeditor on database "<<dbin_name<<endl;
        logfile << "RFeditor run on " << dbin_name << " at time "
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
            <<"Number of existing rows in the wfdisc table is "
            <<dbout.number_tuples()<<endl;
        /* Prep the input wfdisc table*/
        /* WARNING:  not adequate.  This needs to be changed as we 
           need evid.  May be able to fake this by searching for matching
           start times  */
        if(apply_subset)
        {
            cout << "Applying subset condition="<<subset_condition<<endl;
            dbin.subset(subset_condition);
            cout << "Subset view number of rows="<<dbin.number_tuples()<<endl;
        }
        list<string> sortkeys, groupkeys;
        sortkeys.push_back("sta");
        sortkeys.push_back("time");
        sortkeys.push_back("chan");
        groupkeys.push_back("sta");
        dbin.sort(sortkeys);
        dbin.group(groupkeys);
        cout << "Number of ensembles to process (grouped by sta)="
            << dbin.number_tuples()<<endl;
        /* This launches the editing windows */
        RFeditorEngine rfe(control);
        int nsta=dbin.number_tuples();
        TimeSeriesEnsemble radial,transverse;
        for(i=0,dbin.rewind();i<nsta;++dbin)
        {
            string sta;
            cout << "Calling data reader for ensemble number "<<i<<endl;
            TimeSeriesEnsemble dall(dbin,mdl,mdlens,am);
            sta=dall.get_string("sta");
            logfile << "Read "<<dall.member.size()<<" RF traces for station "
                << sta <<endl;
            cout << "Read "<<dall.member.size()<<" RF traces for station "
                << sta <<endl;
            int nstacount=dall.member.size()/2;
            if(nstacount<minrfcutoff)
            {
                cout << "Station "<<sta<<" dropped.   Count below miminum of "
                    << minrfcutoff<<endl;
                logfile << "Station "<<sta
                    <<" dropped.   Count below miminum of "
                    << minrfcutoff<<endl;
                continue;
            }
            /* This routine sets a metadata item eventid based
               on start time only.  Not a bulletproof approach but one 
               that should work for RF data */
            set_eventids(dall);
            /* use start time as 0.  Also set moveout keyword
            to allow stacking in engine */
            vector<TimeSeries>::iterator im;
            for(im=dall.member.begin();im!=dall.member.end();++im)
            {
                im->ator(im->t0);
                im->put(moveout_keyword,0.0);
            }
            /* Could do this with the database, but I chose this 
               algorithm because I think it will be more robust.
               Main reason is I can use only a string fragment for
               a match instead of demanding a full match */

            try {
                radial=extract_by_chan(dall,rchan);
                cout << "Found "<<radial.member.size()<<" radial component RFs"<<endl;
                transverse=extract_by_chan(dall,tchan);
                cout << "Found "<<radial.member.size()<<" transverse component RFs"<<endl;
            }catch(SeisppError& serr)
            {
                cerr << "Problems in extract_by_chan.  Message "
                    <<serr.what()<<endl
                    <<"Skipping ensemble for station = "<<sta
                    <<endl;
                continue;
            }
            /* Note this is a bit fragile.  Requires sta be in mdlens*/
            cout << "Read "<<dall.member.size()
                <<" receiver functions for station"
                <<sta<<endl
                << "Loading data into editor."<<endl
                << "When read edit radial first then transverse"
                <<endl;
            rfe.edit(radial,transverse);
            for(im=radial.member.begin();im!=radial.member.end();++im)
            {
                double t0=im->get_double("time");
                im->rtoa(t0);
                /* Fragile way to handle this, but skip stack traces */
                if((im->get_string("sta"))=="stack") continue;
                if(!(im->live))
                    logfile << "Deleting " <<sta<<":"<< im->get_string("chan")
                        << " for time "<<strtime(t0)<<endl;
            }
            for(im=transverse.member.begin();im!=transverse.member.end();++im)
            {
                double t0=im->get_double("time");
                im->rtoa(t0);
                if(!(im->live))
                    logfile << "Deleting " <<sta<<":"<< im->get_string("chan")
                        << " for time "<<strtime(t0)<<endl;
            }
            int nsaved=save_to_db(radial,transverse,mdlout,am,dbout);
            logfile << "Saved "<<nsaved<<" RFs for station "<<sta<<endl;
            cout << "Saved "<<nsaved<<" RFs for station "<<sta<<endl;
        }
    }catch(SeisppError& serr)
    {
        serr.log_error();
    }
    catch(exception& stdexcept)
    {
        cerr << "Exception thrown:  "<<stdexcept.what()<<endl;
    }
}


