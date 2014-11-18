#include <list>
#include "coords.h"
#include "seispp.h"
#include "SphericalCoordinate.h"
#include "dbpp.h"
#include "Hypocenter.h"
#include "ThreeComponentSeismogram.h"
#include "PfStyleMetadata.h"
DatascopeHandle BuildCatalogView(DatascopeHandle& dbh)
{
	DatascopeHandle result(dbh);
	dbh.lookup("event");
        if(SEISPP_verbose)
            cerr << "Number of rows in event table="<<dbh.number_tuples()<<endl;
	list<string> sortkeys;
	sortkeys.push_back("evid");
	dbh.sort(sortkeys);
	dbh.natural_join("origin");
	string ss_to_prefor("orid==prefor");
	dbh.subset(ss_to_prefor);
        if(SEISPP_verbose)
            cerr << "Number of rows after join origin with orid==prefor"
                << dbh.number_tuples()<<endl;
	dbh.natural_join("assoc");
        if(SEISPP_verbose)
            cerr << "Number of rows after joint with assoc = "
                << dbh.number_tuples()<<endl;
	dbh.natural_join("arrival");
        if(SEISPP_verbose)
            cerr << "Final catalog view size (joining arrival)="
                << dbh.number_tuples()<<endl;
	return(dbh);
}
DatascopeHandle BuildWaveformView(DatascopeHandle& dbh)
{
    DatascopeHandle result(dbh);
    dbh.lookup("wfdisc");
    if(SEISPP_verbose)
            cerr << "BuildWaveformView:  "
                << "Number of rows in wfdisc="<<dbh.number_tuples()<<endl;
    dbh.natural_join("site");
    if(SEISPP_verbose)
            cerr << "BuildWaveformView:  "
                << "Number of rows after natural join of site table="
                <<dbh.number_tuples()<<endl;
    dbh.natural_join("sitechan");
    if(SEISPP_verbose)
            cerr << "BuildWaveformView:  "
                <<"Number of rows after joining sitechan="
                <<dbh.number_tuples()<<endl;
    return(dbh);
}
/* This builds working view by a join of the catalog view (dbcv)
   and the waveform view (dbwf) constructed using the functions
   above.  The view is then sorted and grouped for processing to
   build three-component seismograms. */
DatascopeHandle BuildWorkingView(DatascopeHandle& dbcv,
        DatascopeHandle& dbwv)
{
    try {
        list<string> jk1,jk2,sortkeys,groupkeys;
        DatascopeHandle result(dbcv);
        jk1.push_back(string("arrival.time"));
        jk1.push_back(string("sta"));
        jk2.push_back(string("wfdisc.time::wfdisc.endtime"));
        jk2.push_back(string("sta"));
        result.join(dbwv,jk1,jk2);
        if(SEISPP_verbose)
            cerr << "BuildWorkingView:  "
                << "Working view total final number of rows ="
                << result.number_tuples()<<endl;
        /*This is needed for 3c bundles */
        sortkeys.push_back(string("origin.time"));
        sortkeys.push_back(string("sta"));
        sortkeys.push_back(string("chan"));
        result.sort(sortkeys);
        groupkeys.push_back(string("sta"));
        result.group(groupkeys);
        return(result);
    }
    catch(...){throw;};
}
Hypocenter extract_hypocenter(ThreeComponentSeismogram& d)
{
    double slat,slon,depth,otime;
    try {
        slat=d.get_double("origin.lat");
        slon=d.get_double("origin.lon");
        depth=d.get_double("origin.depth");
        otime=d.get_double("origin.time");
        /* Freeze method and model as tttaup and iasp91.  Not elegant
           but minor importance for now */
        return(Hypocenter(rad(slat),rad(slon),depth,otime,
                    string("tttaup"),string("iasp91")));
    }catch(SeisppError& serr)
    {
        throw serr;
    }
}

        
void usage()
{
    cerr << "RFrotate dbin dbout [-lqt|-fst -v -pf pffile]"<<endl;
    exit(-1);
}
bool SEISPP::SEISPP_verbose(false);
int main(int argc, char **argv)
{
    /* parse the argument list*/
    if(argc<3) usage();
    string dbin_name(argv[1]);
    string dbout_name(argv[2]);
    cout << "RFrotate:   Reading from database "
        << dbin_name<<endl
        << "Writing results to database " <<dbout_name<<endl;
    int i;
    string rotation_type("rtz");
    string pfname("RFrotate");
    for(i=3;i<argc;++i)
    {
        string sarg(argv[i]);
        if(sarg=="-lqt")
            rotation_type=sarg;
        else if(sarg=="-fst")
            rotation_type=sarg;
        else if(sarg=="-v")
            SEISPP_verbose=true;
        else if(sarg=="-pf")
        {
            ++i;
            if(i>=argc) usage();
            sarg=string(argv[i]);
            pfname=sarg;
        }
        else
            usage();
    }
    try {
        PfStyleMetadata md=pfread(pfname);
        /* We need this list of attributes to load from working
           view we pull from the parameter file.  */
        MetadataList tracemdl=get_mdlist(md,string("trace_mdl"));
        MetadataList outmdl=get_mdlist(md,string("output_wfdisc_mdl"));
        list<string> output_channels=md.get_tbl("output_channel_codes");
        AttributeMap am("css3.0");
        double vp0,vs0;
        vp0=md.get_double("vp0");
        vs0=md.get_double("vs0");

        /* First prepare all the database handles */
        DatascopeHandle dbin(dbin_name,true);
        DatascopeHandle dbout(dbout_name,false);
        dbin=BuildCatalogView(dbin);
        DatascopeHandle dbwf(dbin);
        dbwf=BuildWaveformView(dbwf);
        dbin=BuildWorkingView(dbin,dbwf);
        dbout.lookup("wfdisc");
        /* We loop over the group pointer to get one 3c seismogram 
           at a time */
        int nrows=dbin.number_tuples();
        for(i=0;i<nrows;++i,++dbin)
        {
          try{
            ThreeComponentSeismogram d(dbin,tracemdl,am);
            /* This may not be essential, but better safe than sorry */
            d.rotate_to_standard();
            Hypocenter h(extract_hypocenter(d));
            double rlat,rlon,relev;
            rlat=d.get_double("site.lat");
            rlon=d.get_double("site.lon");
            relev=d.get_double("site.elev");
            rlat=rad(rlat);  rlon=rad(rlon);
            /* For now assume P wave phase for slowness vector when
               the fst is used */
            SlownessVector u=h.pslow(rlat,rlon,relev);
            double az=u.azimuth();
            SphericalCoordinate sc_for_rotation;
            if(rotation_type=="fst")
            {
                d.free_surface_transformation(u,vp0,vs0);
            }
            else if(rotation_type=="lqt")
            {
                sc_for_rotation=PMHalfspaceModel(vp0,vs0,u.ux,u.uy);
                d.rotate(sc_for_rotation);
            }
            else
            {
                sc_for_rotation.radius=1.0;
                sc_for_rotation.theta=0.0;
                sc_for_rotation.phi=M_PI_2-az;
                d.rotate(sc_for_rotation);
            }
            /* Now we unbundle the data and write it to output db*/
            TimeSeries *component;
            list<string>::iterator iptr;
            for(i=0,iptr=output_channels.begin();i<3;++i,++iptr) 
            {
                component=ExtractComponent(d,i);
                component->put("chan",*iptr);
                dbsave(*component,dbout.db,string("wfdisc"),outmdl,am);
                delete component;
            }
          }catch(SeisppError& serr)
          {
              cerr << "Something threw a SeisppError exception in main processing"
                  <<endl<<"Error message follows:"<<endl;
              serr.log_error();
          }
        }
    }catch(SeisppError& serr)
    {
        cerr << "Fatal error before processing started"<<endl;
        serr.log_error();
    }
    catch(...)
    {
        cerr << "Something threw an unexpected exception before start of"
            <<" processing"<<endl<<"No data processed"<<endl;
    }
}
