#include <string>
#include <fstream>
#include <list>
#include "stock.h"
#include "AttributeCrossReference.h"
#include "GenericFileHandle.h"
#include "MatlabProcessor.h"
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
/* These small procedures encapsulates the instructions sent
   to matlab to do the picking operation. This is useful to 
   make sure all this is one place */
const string ensemblename("ensemble"),backgroundimagename("bi");
const string mfilename("seispicker");
string onefile_matlab_instructions()
{
    string mi;
    mi=string("[xp,yp]=")
        + mfilename
        + "(" + ensemblename +"," + ensemblename + ",xp,yp);";
    return mi;
}
string twofile_matlab_instructions()
{
    string mi;
    mi=string("[xp,yp]=")
        + mfilename
        + "(" + ensemblename +"," + backgroundimagename + ",xp,yp);";
    return mi;
}
void usage()
{
    cerr << "horizonpicker infile outfile [-b bifile -pf pffile -v]"<<endl
        <<endl<<"where bifile will be used as background image"
       << endl << "Note bifile must be in same format as infile"
        <<endl;
    exit(-1);
}


int main(int argc, char **argv)
{
    if(argc<3) usage();
    cout << "horizonpicker program:  uses matlab as plot engine"<<endl;
    string infile(argv[1]);
    string outfile(argv[2]);
    string pffile("horizonpicker");
    bool usingbifile;
    string bifile("NONE");
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
        else if(sarg=="-b")
        {
            ++i;
            if(i>=argc) usage();
            usingbifile=true;
            sarg=string(argv[i]);
            bifile=sarg;
        }
        else if(sarg=="-v")
        {
            SEISPP_verbose=true;
        }
        else
            usage();
    }
    Pf *pf;
    if(pfread(const_cast<char *>(pffile.c_str()),&pf))
    {
        cerr << "pfread failed for pffile="<<pffile<<endl;
        usage();
    }
    try{
        ofstream outstrm(outfile.c_str(),ios::out);
        string xrefdata=pftbl2string(pf,"xrefdata");
        AttributeCrossReference xref(xrefdata);
        list<string> orderkeys,endslist,tmdlist;
        orderkeys=pftbl2list(pf,"orderkeys");
        endslist=pftbl2list(pf,"ensemble_metadata_list");
        tmdlist=pftbl2list(pf,"trace_metadata_list");
        Metadata control(pf);
        string traceformat=control.get_string("traceformat");
        string logfilename=control.get_string("logfile");
        double yscale=control.get_double("yaxis_scaling");
        MetadataList mdout=pfget_mdlist(pf,"output_file_metadata_list");
        FILE *fp=fopen(logfilename.c_str(),"w");
        if(SEISPP_verbose)
        {
            cout << "Constructing file handle for file="<<infile
                << " with format name="<< traceformat<<endl;
        }
        GenericFileHandle handle(infile,traceformat,
                xref,orderkeys,endslist,tmdlist,true,string("nsamp"),
                string("dt"),1.0,true);
        GenericFileHandle *bifhandle;
        if(usingbifile)
        {
            cout << "Using data from file="<<bifile
                <<" as background image for wiggle plot"<<endl;
            bifhandle=new GenericFileHandle(bifile,traceformat,
                xref,orderkeys,endslist,tmdlist,true,string("nsamp"),
                string("dt"),1.0,true);
        }
        if(SEISPP_verbose) cout << "Launching matlab"<<endl;
        MatlabProcessor mp(fp);
        // Some initialization.  Warning:  frozen to variable in m file
        string prep("clim=[-1,1];xp=[0.0,0.0];yp=[0.0,0.0];");
        mp.process(prep);
        auto_ptr<TimeSeriesEnsemble> din,bidata;
        vector<TimeSeries>::iterator mptr;
        int ensemblenumber(1);
        string runpickinstructions;
        if(usingbifile)
            runpickinstructions=twofile_matlab_instructions();
        else
            runpickinstructions=onefile_matlab_instructions();
        string ans;
        do {
            cout << "Loading ensemble number "<<ensemblenumber<<endl;
            din=handle.GetNextEnsemble();
            if(usingbifile) bidata=bifhandle->GetNextEnsemble();
            cout << "Number of members read = "<<din->member.size()<<endl
                << "Sending data to matlab for plotting and picking"<<endl
                << "Select points with mb1.  Hit enter to exit pick loop."<<endl;
            // This is the symbolic name assigned to the data matrix in matlab
            const string matlabname("ensemble"),matlabbiname("bi"); 
            string ques;
            mp.load(*din,matlabname);
            if(usingbifile) mp.load(*bidata,matlabbiname);
            do {
                mp.process(runpickinstructions);
                cout << "Picks ok? Enter n to repeat, d to drop, y to save: ";
                cin >> ques;
            }while(ques=="n" || ques=="no");
            if(ques=="d") continue;
            vector<double> x=mp.retrieve_vector(string("xp"));
            vector<double> y=mp.retrieve_vector(string("yp"));
            cout << "Saving "<<x.size()<< " picks."<<endl;
            vector<double>::iterator xptr,yptr;
            for(xptr=x.begin(),yptr=y.begin();xptr!=x.end()
                    && yptr!=y.end();++xptr,++yptr)
            {
                i=static_cast<int>(*xptr);
                if(i<0) 
                {
                    cerr << "Warning pick "<< *xptr << " yields negative index"
                        << endl << "Set to 0"<<endl;
                    i=0;
                }
                else if(i>=din->member.size())
                {
                    cerr << "Warning pick "<< *xptr << " overflows array size ="
                        << din->member.size() <<" reset to last trace in gather"
                        <<endl;
                    i=din->member.size() - 1;
                }
                MetadataList::iterator mptr;
                for(mptr=mdout.begin();
                        mptr!=mdout.end();++mptr)
                {
                    switch(mptr->mdt)
                    {
                        case MDreal:
                            outstrm << din->member[i].get_double(mptr->tag);
                            break;
                        case MDint:
                            outstrm << din->member[i].get_int(mptr->tag);
                            break;
                        case MDstring:
                            outstrm << din->member[i].get_string(mptr->tag);
                            break;
                        case MDboolean:
                            outstrm << din->member[i].get_bool(mptr->tag);
                            break;
                        default:
                            cerr << "Error:  Metadata with tag = "
                                << mptr->tag << " has unknown type number "
                                << mptr->mdt<<endl;
                    }
                    outstrm << " ";
                }
                outstrm <<" "<< *yptr <<endl;
            }
            ++ensemblenumber;
        }while(!handle.eof());
    } catch (SeisppError& serr)
    {
        serr.log_error();
        exit(-1);
    }
}
