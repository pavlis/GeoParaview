#include <string>
#include <iostream>
#include <list>
#include "stock.h"
#include "SacFileHandle.h"
using namespace std;
using namespace SEISPP;
bool SEISPP::SEISPP_verbose(true);

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
    cerr << "sac2su [-o outfile -pf pffile] < list_of_sac_files"<<endl
        << "     default outfile is sacfiles.su"<<endl
        << "     list_of_sac_files is white space and/or lines of sac file names"
        <<endl;;
    exit(-1);
}
/* Helpers to handle sac defaults. Return true when not default value */
bool not_default_real(double x)
{
    if(fabs(x+12345.0)<0.1)
        return false;
    else
        return true;
}

/* This handles special header/metadata entries that deal with 
   clashes in concept that have to be handled specially.   
   */
void SacToSUMetadataTranslation(TimeSeries& d)
{
    try {
        /* Force live if we get here */
        d.live=true;
        /* Constant for all.  Defines all as seismic data */
        d.put("trid",1);
        /* This defines a 1000 multiplier and sets decimal degrees for
           coordinates */
        d.put("counit",3);
        d.put("scalco",1000.0);
        d.put("scalel",1000.0);
        double rlat,rlon,elev,edep;
        rlat=d.get_double("stla");
        rlon=d.get_double("stlo");
        elev=d.get_double("stel");
        edep=d.get_double("stdp");
        if(not_default_real(rlat))
        {
            d.put("ry",rlat*1000.0);
        }
        else
            d.put("ry",0.0);
        if(not_default_real(rlon))
        {
            d.put("rx",rlon*1000.0);
        }
        else
            d.put("rx",0.0);
        
        if(not_default_real(elev))
        {
            double etot=elev;
            if(not_default_real(edep))
                etot-=edep;
            d.put("relev",etot*1000.0);
        }
        else
            d.put("relev",0.0);
        double slat,slon;
        slat=d.get_double("evla");
        slon=d.get_double("evlo");
        elev=d.get_double("evel");
        edep=d.get_double("evdp");
        if(not_default_real(slat))
        {
            d.put("sy",slat*1000.0);
        }
        else
            d.put("sy",0.0);
        if(not_default_real(slon))
        {
            d.put("sx",slon*1000.0);
        }
        else
            d.put("sx",0.0);
        /* this is a bit weird because we use sdepth of segy
           but this is computed from edep and elev of sac */
        if(not_default_real(edep) && not_default_real(elev))
        {
            double etot=elev;
            if(not_default_real(edep))
                etot-=edep;
            d.put("sdepth",etot*1000.0);
        }
        else
        {
            d.put("sdepth",0.0);
        }
        double dist=d.get_double("dist");
        /* Not really correct, but we'll just copy this as is.
           Wrong because sac normally has this in degrees while
           su assumes m or ft */
        d.put("offset",dist);
        /* a very bad way to handle a type collision here.  
           su wants the sec field to be real but sac has an int
           for this field in it's header.   */
        int sec=d.get_int("nzsec");
        int msec=d.get_int("msec");
        double fsec=static_cast<double>(sec)+static_cast<double>(msec)/1000.0;
        d.put("sec",fsec);
    }catch(...){throw;};
}
int main(int argc, char **argv)
{
    string outfile("sacfiles.su");
    string pffile("sac2su");
    int i;
    for(i=1;i<argc;++i)
    {
        string sarg(argv[i]);
        if(sarg=="-o")
        {
            ++i;
            if(i==argc) usage();
            outfile=string(argv[i]);
        }
        else if(sarg=="-pf")
        {
            ++i;
            if(i==argc) usage();
            pffile=string(argv[i]);
        }
        else
            usage();
    }

    Pf *pf;
    if(pfread(const_cast<char *>(pffile.c_str()),&pf))
    {
        cerr << "pfread failed for sac2su.pf"<<endl;
        exit(-1);
    }
    int trace_number(1);
    try{
        string inxref=pftbl2string(pf,"SAC_metadata_cross_reference");
        AttributeCrossReference sacxref(inxref);
        string outxref=pftbl2string(pf,"SU_metadata_cross_reference");
        AttributeCrossReference suxref(outxref);
        list<string> tmdlist=pftbl2list(pf,"output_metadata_list");
        list<string> sacmdlist=pftbl2list(pf,"SAC_metadata_list");
        /* output is global so we create it now.  Sac input is 
           commonly many many files */
        /* These simply need to be declared.  They are ignored 
           since the SU file is output */
        list<string> orderkeys,endslist;
        // force initialization
        orderkeys.clear();
        endslist.clear();
        GenericFileHandle suout(outfile,string("SEGYfloat"),
                suxref,orderkeys,endslist,tmdlist,false,string("nsamp"),
                string("dt"),1000000.0,true);
        /* Now loop over the arg list skipping -o option */
        do {
            string infile;
            cin >> infile;
            if(cin.eof()) break;
            cout << "Processing sac file "<<infile<<endl;
            /* We enclose this section in a try/catch block to 
               allow us to drop data from files that cannot be
               openned.  This should, perhaps, be optional. */
            try {
                /* orderkeys and endslist are empty.  SAC has no such concept.
                This is baggage because I couldn't figure out how to build the
                constructor without these as baggage.  Should clean this up if
                I ever elect to release this to the world. */
                SacFileHandle sacin(infile,sacxref,sacmdlist,orderkeys,endslist);
                TimeSeries d=sacin.GetNextSeismogram();
                SacToSUMetadataTranslation(d);
                /* These are counters in segy.  We handle these crudly 
                   for now */
                d.put("fldr",1);
                d.put("tracl",trace_number);
                d.put("tracf",trace_number);
                d.put("tracr",trace_number);
                //cout << dynamic_cast<Metadata&>(d);
                suout.put(d);
                ++trace_number;
            } catch(SeisppError& serr)
            {
                serr.log_error();
                cerr << "Data for file="<<infile
                    << " was not copied to output file."<<endl;
            }
        } while(!cin.eof());
    }
    catch (SeisppError& serr)
    {
        serr.log_error();
        exit(1);
    }
}
