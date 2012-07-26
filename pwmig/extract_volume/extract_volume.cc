#include <fstream>
#include "stock.h"
#include "dbpp.h"
#include "TimeSeries.h"
#include "Metadata.h"
#include "seispp.h"
#include "gclgrid.h"
#include "FixedFormatTrace.h"
#include "GenericFileHandle.h"
/* These are functions in pfutils.cc.  Eventually they will probably land
   in Metadata.h */
string pftbl2string(Pf *pf, const char* tag);
list<string> pftbl2list(Pf *pf, const char* tag);
using namespace std;
using namespace SEISPP;
/* This is a less general version of a similar program called extract_section.
   This program blindly dumps a GCLvectorfield writing the output as seismic
   unix data made to look like 3D seismic data.  It does this by blindly 
   extracting data as original grid point data.  Output is in x1, x2 order with
   each seismogram having 0 at the surface (reflection convention).  

   Header data to write is frozen and will evolve - writing this as a design
   document
   */

void usage()
{
	cerr << "extract_volume db gridname fieldname "
		<< "[-o outfile -pf pffile -v]"<<endl
	      << "Default outfile is called extract_volume.su"<<endl
              << "Output is always in old SU format"<<endl;

	exit(-1);
}
/* Extract i,j grid line (k=0 to g.n3) from grid as a 3C seismogram.   
The k line is reversed because of the way pwmig writes the grid and 
we extract the first 3 components.   This routine does not error checking.
The assumption is the function is locked to this special code and is
not generic.
*/
ThreeComponentSeismogram ExtractFromGrid(GCLvectorfield3d& g, int i, int j)
{
    ThreeComponentSeismogram result(g.n3);
    int k,kk,l;
    for(k=0,kk=g.n3-1;k<g.n3;++k,--kk)
        for(l=0;l<3;++l) result.u(l,k)=g.val[i][j][kk][l];
    /* Required metadata */
    result.put("nsamp",g.n3);
    /* make dt units of meters */
    result.put("dt",g.dx3_nom*1000.0);
    result.ns=g.n3;
    /* This pair doesn't include 1000 factor because segy oddity of
       microsection dt */
    result.dt=g.dx3_nom;
    result.live=true;
    result.put("samprate",1.0/g.dx3_nom);
    result.tref=relative;
    result.t0=0.0;
    result.put("time",0.0);
    /* SU specific */
    /* use tracl to define lines - i fast becomes inline */
    result.put("tracl",i);
    result.put("tracr",1);
    result.put("fldr",1);
    /* This makes the data look like cdp stacked data */
    result.put("cdpt",1);
    result.put("cdp",j*g.n1+i);
    /*segy code for live seismic data. SU allows variable codes
    set later */
    result.put("trid",1);
    /* degree multiplier for lat lon */
    result.put("scalco",10000.0);
    double lat,lon;
    lat=g.lat(i,j,g.n3-1);
    lon=g.lon(i,j,g.n3-1);
    lat=deg(lat);   lon=deg(lon);
    result.put("ry",lat);
    result.put("rx",lon);
    result.put("duse",1);
    return(result);
}
bool SEISPP::SEISPP_verbose(false);
int main(int argc, char **argv)
{
	int i,j;
	ios::sync_with_stdio();
	if(argc<4) usage();
	string dbname(argv[1]);
	string gridname(argv[2]);
	string fieldname(argv[3]);
        string pffile("extract_volume");
	ofstream outstrm;
	bool out_to_other(false);
	string outfile("extract_volume.su");;
	for(i=4;i<argc;++i)
	{
		string argstr=string(argv[i]);
		if(argstr=="-o")
		{
		    ++i;
                    if(i>=argc) usage();
		    out_to_other=true;
		    outfile=string(argv[i]);
		}
                else if(argstr=="-pf")
                {
                    ++i;
                    if(i>=argc) usage();
                    pffile=string(argv[i]);
                }
		else if(argstr=="-v")
			SEISPP_verbose=true;
		else
		{
			cerr << "Unknown argument = "<<argstr<<endl;
			usage();
		}
	}
        Pf *pf;
        if(pfread(const_cast<char *>(pffile.c_str()),&pf))
        {
            cerr << "pfread failed for pf file="<<pffile<<endl;
            exit(-1);
        }
	try {
            /*may eventually want this */
            //Metadata control(pf);
            DatascopeHandle dbh(dbname,true);
            dbh.lookup(string("gclgdisk"));
	    DatascopeHandle dbhg(dbh);
	    dbhg.lookup(string("gclfield"));
	    Dbptr db=dbh.db;
	    Dbptr dbgrd=dbhg.db;
	    GCLvectorfield3d g(dbh,gridname,fieldname,5);
            string xrefstr=pftbl2string(pf,"metadata_cross_reference");
            AttributeCrossReference outxref(xrefstr);
            list<string> tmdlist=pftbl2list(pf,"output_metadata_list");
            /* We create empty lists for these for now since we 
               support only SU output */
            list<string> orderkeys,endslist;
            orderkeys.clear();   endslist.clear();
            /* For now support only SU format.  Could be generalized */
            GenericFileHandle suout(outfile,string("SEGYfloat"),
                    outxref,orderkeys,endslist,tmdlist,false,
                    string("nsamp"),string("dt"),1000000.0,false);
            /* Now loop over grid.  */
            int i,j;
            ThreeComponentSeismogram d;
            for(j=0;j<g.n2;++j)
                for(i=0;i<g.n1;++i)
                {
                    d=ExtractFromGrid(g,i,j);
                    /* We have to extract and write each component
                       as scalar becaue segy/su are not explicitly 
                       vector formats.  write will not allow this. */
                    for(int k=0;k<3;++k)
                    {
                        TimeSeries *comp=ExtractComponent(d,k);
                        /* These are special codes for SU */
                        switch (k)
                        {
                            case (0):
                                comp->put("trid",13);
                                break;
                            case (1):
                                comp->put("trid",14);
                                break;
                            case (2):
                            default:
                                comp->put("trid",12);
                        }
                        //DEBUG
                        cout << *(dynamic_cast<Metadata*>(comp))<<endl;
                        suout.put(*comp);
                        delete comp;
                    }
                }
	}

	catch (SeisppError& serr)
	{
		serr.log_error();
	}
}
