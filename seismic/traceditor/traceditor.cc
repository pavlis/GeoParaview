#include <string>
#include <fstream>
#include <list>
#include "stock.h"
#include "AttributeCrossReference.h"
#include "GenericFileHandle.h"
#include "TraceEditPlot.h"
using namespace std;
using namespace SEISPP;
bool SEISPP::SEISPP_verbose(false);
string pftbl2string(Pf *pf, const char* tag)
{
    Tbl *t;
    t=pfget_tbl(pf,const_cast<char *>(tag));
    if(t==NULL) throw SeisppError(string("pftbl2string:  ")
            + "pf file is missing Tbl with tag = "
            +tag);
    string result;
    for(int i=0;i<maxtbl(t);++i)
    {
        char *line=(char *)gettbl(t,i);
        result += line;
        result += "\n";
    }
    return result;
}
list<string> pftbl2list(Pf *pf, const char* tag)
{
    Tbl *t;
    t=pfget_tbl(pf,const_cast<char *>(tag));
    if(t==NULL) throw SeisppError(string("pftbl2list:  ")
            + "pf file is missing Tbl with tag = "
            +tag);
    list<string> result;
    for(int i=0;i<maxtbl(t);++i)
    {
        char *line=(char *)gettbl(t,i);
        result.push_back(string(line));
    }
    return result;
}
void usage()
{
    cerr << "traceditor infile outfile [-pf pffile -keepkills -v]"<<endl
        <<"  -keepkills options copies killed traces but marks them dead"
        <<endl;
    exit(-1);
}


int main(int argc, char **argv)
{
    /* First step as always is to crack the command line */
    if(argc<3) usage();
    string infile(argv[1]);
    string outfile(argv[2]);
    string pffile("traceditor");
    bool nowritedead(true);
    int i;
    for(i=3;i<argc;++i)
    {
        string sarg(argv[i]);
        if(sarg=="-pf")
        {
            ++i;
            if(i>=argc) usage();
            sarg=string(argv[i]);
            pffile=sarg;
        }
        else if(sarg=="-keepkills")
        {
            nowritedead=false;
        }
        else if(sarg=="-v")
        {
            SEISPP_verbose=true;
        }
        else
            usage();
    }
    /* This code uses an Antelope parameter file to define
       a large set of parameters.  A large fraction, in this case,
       are defaulted particularly for the seismic display widget */
    Pf *pf;
    if(pfread(const_cast<char *>(pffile.c_str()),&pf))
    {
        cerr << "pfread failed for pffile="<<pffile<<endl;
        usage();
    }
    /* Enter C++ area where we need to handle a range of exceptions */
    try{
        /* Extract and build a set of entities required by the
           generic file handles*/
        string xrefdata=pftbl2string(pf,"xrefdata");
        AttributeCrossReference xref(xrefdata);
        list<string> orderkeys,endslist,tmdlist;
        orderkeys=pftbl2list(pf,"orderkeys");
        endslist=pftbl2list(pf,"ensemble_metadata_list");
        tmdlist=pftbl2list(pf,"trace_metadata_list");
        Metadata control(pf);
        string traceformat=control.get_string("traceformat");
        double yscale=control.get_double("yaxis_scaling");
        if(SEISPP_verbose)
        {
            cout << "Constructing file handle for file="<<infile
                << " with format name="<< traceformat<<endl;
        }
        /* This constructs the input file handle */
        GenericFileHandle handle(infile,traceformat,
                xref,orderkeys,endslist,tmdlist,true,string("nsamp"),
                string("dt"),1.0,true);
        /* Force the output file to be the same format as input. 
           Format conversion is a different problem */
        GenericFileHandle outhandle(outfile,traceformat,
                xref,orderkeys,endslist,tmdlist,false,string("nsamp"),
                string("dt"),1.0,true,nowritedead);
        /* Launch plot window */
        TraceEditPlot editor(control);
        /* Use an auto_ptr instead of defining a TimeSeriesEnsemble
           for efficiency.  Define in the interface to GetNextEnsemble */
        auto_ptr<TimeSeriesEnsemble> din;
        int ensemblenumber(1);
        do {
            cout << "Loading ensemble number "<<ensemblenumber<<endl;
            din=handle.GetNextEnsemble();
            cout << "Number of seismograms read = "<<din->member.size()<<endl
                << "Calling plot gui. "<<endl;
            editor.plot(*din);
            set<int> kills=editor.report_kills();
            set<int>::iterator ikptr;
            for(ikptr=kills.begin();ikptr!=kills.end();++ikptr)
            {
                int ikill=*ikptr;
                din->member[ikill].live=false;
            }
            outhandle.put(*din);
            ++ensemblenumber;
        }while(!handle.eof());
    } catch (SeisppError& serr)
    {
        serr.log_error();
        exit(-1);
    }
    catch (...)
    {
        cerr << "Something threw an unexpected exception type"<<endl;
        exit(-1);
    }
}
