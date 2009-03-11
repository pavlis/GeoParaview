#include <exception>
#include "seispp.h"
#include "stack.h"
#include "Hypocenter.h"
#include "ThreeComponentSeismogram.h"
#include "ensemble.h"
#include "Metadata.h"
#include "dbpp.h"
using namespace std;
using namespace SEISPP;
/* This procedure does not actually seem to be necessary.  Think we
can just build a view containing the hypocentroid and we should be fine.  
Will retain this, however, until I'm sure that is so */

/***************************
EventCatalog LoadHypocentroids(DatascopeHandle dbh,
    string gridname,
	string firstotime, double otoffset,ttmethod,ttmodel)
{
    const string base_error("LoadHypocentroid:  ");
    try {
	double ot0=str2epoch(firstotime.c_str*());
	double hclat,hclon,hcdepth,otimei(ot0);
	MetadataList mdl;
	const Metadata_typedef mdaux={string("gridid"),MDint};
	mdl.push_back(mdaux);
	EventCatalog result(mdl);
	dbh.lookup("hypocentroid");
	dbh.subset(string("gridname=~/")+gridname+"/";
	dbh.rewind();
	int nrec=dbh.number_tuples();
	if(nrec<=0) throw SeisppError(base_error
			+ "hypocentroid table has no data for table"
			+gridname);
	int i,gridid,nass;
	Metadata aux;
	for(i=0;i<nrec;++i,++dbh,otime+=otoffset)
	{
		lat=dbh.get_double("hclat");
		lon=dbh.get_double("hclon");
		depth=dbh.get_double("hcdepth");
		gridid=dbh.get_int("gridid");
		nass=dbh.get_int("nass");
		aux.put("gridid",gridid);
		aux.put("nass",nass);
		Hypocenter h(lat,lon,depth,otime,ttmethod,ttmodel);	
		if(!result.add(h,aux))
		{
			throw SeisppError(base_error
				+ "hypocentroid has duplicates.  fix db");
		}
	}
	// loop over db, load hypo data, set otime, then call "add" method
	// of EventCatalog -- one for each hypo
    } catch (...) {throw;}
}
**********************/

Hypocenter MakeHypocentroid(ThreeComponentEnsemble& d,double otime,
	string ttmethod, string ttmodel)
{
    try {
    	double lat,lon,depth;
	lat=d.get_double("hclat");
	lon=d.get_double("hclon");
	depth=d.get_double("hcdepth");
	Hypocenter h(rad(lat),rad(lon),depth,otime,ttmethod,ttmodel);
	return(h);
    } catch (...){throw;};
}

DatascopeHandle BuildWaveformView(DatascopeHandle& dbhin)
{
	try {
		DatascopeHandle dbh(dbhin);
		dbh.lookup("hypocentroid");
		dbh.natural_join("cluster");
		dbh.natural_join("event");
		dbh.natural_join("origin");
		dbh.subset("orid==prefor");
		dbh.natural_join("assoc");
		dbh.natural_join("arrival");
		DatascopeHandle ljhandle(dbh);
		ljhandle.lookup("wfprocess");
		ljhandle.natural_join("sclink");
		ljhandle.natural_join("evlink");
		list<string> jk;
		jk.push_back("evid");
		jk.push_back("sta");
		dbh.join(ljhandle,jk,jk);
		dbh.natural_join("site");
		list<string> sortkeys;
		sortkeys.push_back("gridid");
		sortkeys.push_back("sta");
		sortkeys.push_back("evid");
		dbh.sort(sortkeys);
		list<string> group_keys;
		group_keys.push_back("gridid");
		group_keys.push_back("sta");
		dbh.group(group_keys);
		dbh.rewind();
		return(dbh);
	} catch (...) {throw;};
}
/* exit with a message if a table has any existing content.
Required for this program.  Bad bad form, but moves us forward.
Ultimately this will be replace by output to a file or an HPC db 
that is not relational */
void exit_if_exists(DatascopeHandle& dbh,char *table)
{
	int nrec=dbh.number_tuples();
	if(nrec>0)
	{
	    cerr << "RFeventstacker:  output database table "
		<< table << " is not empty"<<endl
		<< "This program requires an empty output db"<<endl
		<< "Change output name or remove old tables"<<endl;
	    exit(-1);
	}
}
/* These two procedures are special for this program.  They write a minimal
amount of information ot arrival and assoc. */
long int save_arrival(DatascopeHandle& dbh, ThreeComponentSeismogram& d,
		string arrivalchan, string phase)
{
    try{
	double atime=d.get_double("arrival.time");
	string sta=d.get_string("sta");
	long jday=yearday(atime);
	dbh.append();
	dbh.put("sta",sta);
	dbh.put("time",atime);
	long arid;
	arid = dbnextid(dbh.db,"arid");
	dbh.put("arid",(int)arid);
	dbh.put("jdate",(int)jday);
	dbh.put("chan",arrivalchan);
	dbh.put("iphase",phase);
	return(arid);
    }catch (...) {throw;};
}
void save_assoc(DatascopeHandle& dbh,ThreeComponentSeismogram& d, 
		int orid, int arid, string phase, Hypocenter& h)
{
    try{
	string sta=d.get_string("sta");
	double stalat,stalon,staelev;
	stalat=d.get_double("site.lat");
	stalon=d.get_double("site.lon");
	staelev=d.get_double("site.elev");
	double delta=h.distance(stalat,stalon);
	double seaz=h.seaz(stalat,stalon);
	double esaz=h.esaz(stalat,stalon);
	string timedef("d");
	dbh.append();
	dbh.put("arid",arid);
	dbh.put("orid",orid);
	dbh.put("sta",sta);
	dbh.put("phase",phase);
	dbh.put("delta",delta);
	dbh.put("seaz",seaz);
	dbh.put("esaz",esaz);
	dbh.put("timedef",timedef);
	dbh.put("timeres",0.0);
    }catch (...) {throw;};
}
void usage()
{
	cerr << "RFeventstacker dbin dbout [-pf pfname]" << endl
		<< "dbout must be empty"<<endl;;
	exit(-1);
}
bool SEISPP::SEISPP_verbose(true);
int main(int argc, char **argv)
{
	ios::sync_with_stdio();
	if(argc<3) usage();
	string pfin(argv[0]);
	string dbinname(argv[1]);
	string dboutname(argv[2]);
        for(int i=3;i<argc;++i)
        {
                if(!strcmp(argv[i],"-V"))
                        usage();
                else if(!strcmp(argv[i],"-pf"))
                {
                        ++i;
                        if(i>=argc) usage();
                        pfin = string(argv[i]);
                }
                else
                        usage();
        }
	Pf *pf;
        if(pfread(const_cast<char *>(pfin.c_str()),&pf)) 
	{
		cerr << "pfread error for pf file="<<pfin<<".pf"<<endl;
		exit(-1);
	}
	MetadataList mdlin=pfget_mdlist(pf,"Station_mdlist");
	MetadataList mdens=pfget_mdlist(pf,"Ensemble_mdlist");
	MetadataList mdwfprocess=pfget_mdlist(pf,"wfprocess_mdlist");
	MetadataList mdwfdisc=pfget_mdlist(pf,"wfdisc_mdlist");
	try {
		Metadata control(pf);
		/* This perhaps should be a parameter, but frozen for now*/
		const string atkey("arrival.time");
		double ts,te;
		ts=control.get_double("stack_starttime");
		te=control.get_double("stack_endtime");
		TimeWindow twin(ts,te);
		string phase=control.get_string("phase_for_alignment");
		string arrivalchan=control.get_string("arrival_chan");
		/* All data will be written to this directory and to 
		one file set by dfile.  Intentional to improve performance
		in hpc systems*/
		string dir=control.get_string("output_dir");
		string dfile=control.get_string("output_dfile");
		DatascopeHandle dbh(string(dbinname),true);
		DatascopeHandle dbhwf=BuildWaveformView(dbh);
		dbhwf.rewind();
		DatascopeHandle dbho(string(dboutname),false);
		DatascopeHandle dborigin(dbho);
		dborigin.lookup("origin");
		/* This procedure will exit the program if any of these
		tables have existing content */
		exit_if_exists(dborigin,"origin");
		DatascopeHandle dbevent(dbho);
		dbevent.lookup("event");
		exit_if_exists(dbevent,"event");
		DatascopeHandle dbassoc(dbho);
		dbassoc.lookup("assoc");
		exit_if_exists(dbassoc,"assoc");
		DatascopeHandle dbarrival(dbho);
		dbarrival.lookup("arrival");
		exit_if_exists(dbarrival,"arrival");
		DatascopeHandle dbsclink(dbho);
		dbsclink.lookup("sclink");
		exit_if_exists(dbsclink,"sclink");
		DatascopeHandle dbevlink(dbho);
		dbsclink.lookup("evlink");
		exit_if_exists(dbevlink,"evlink");
		DatascopeHandle dbwfprocess(dbho);
		dbwfprocess.lookup("wfprocess");
		exit_if_exists(dbwfprocess,"wfprocess");
		DatascopeHandle dbwfdisc(dbho);
		dbwfdisc.lookup("wfdisc");
		exit_if_exists(dbwfdisc,"wfdisc");

		const string schema("css3.0");
		AttributeMap am(schema);
		/* This builds a fake catalog of hypocentroids used to 
		guarantee proper associations */
		string firstotime,ttmethod,ttmodel;
		firstotime=control.get_string("first_fake_origin_time");
		double otoffset=control.get_double("fake_origin_time_offsets");
		ttmethod=control.get_string("ttmethod");
		ttmodel=control.get_string("ttmodel");
		double lat,lon,depth,otime;
		otime=str2epoch(const_cast<char *>(firstotime.c_str()));
		/*
		EventCatalog hypocentroids;
		hypocentroids=LoadHypocentroids(dbh,gridname,
			firstotime,otoffset,ttmethod,ttmodel);
		*/
		ThreeComponentEnsemble *pwdataraw;
		auto_ptr<TimeSeriesEnsemble> x1,x2,x3;  // hold components
		vector<string> chanmap;
		chanmap.push_back("E"); 
		chanmap.push_back("N"); 
		chanmap.push_back("Z"); 
		vector<TimeSeries> stack3c;
		stack3c.reserve(3);
		int nrec=dbhwf.number_tuples();
		int evid,orid,arid;
		int record;
		int gridid,lastgridid;
		Hypocenter hcen;
		for(record=0,evid=1,orid=1;record<nrec;++record,++dbh)
		{
			pwdataraw=new ThreeComponentEnsemble(dynamic_cast<DatabaseHandle&>(dbhwf),
				mdlin,mdens,am);
			gridid=pwdataraw->get_int("gridid");
			if(record==0) 
			{
				lastgridid=gridid;
				hcen=MakeHypocentroid(*pwdataraw,otime,
					ttmethod,ttmodel);
			}
			if( (gridid!=lastgridid) || (record==(nrec-1)) )
			{
				dbevent.append();
				dbevent.put("evid",evid);
				dbevent.put("prefor",orid);
				dborigin.put("lat",deg(hcen.lat));
				dborigin.put("lon",deg(hcen.lon));
				dborigin.put("depth",hcen.z);
				dborigin.put("time",hcen.time);
				/* Now we load the next hypo data */
				otime+=otoffset;
				hcen=MakeHypocentroid(*pwdataraw,otime,
					ttmethod,ttmodel);
				++evid;
				++orid;
			}
			string sta=pwdataraw->get_string("sta");
			auto_ptr<ThreeComponentEnsemble> 
			  pwdata = ArrivalTimeReference(*pwdataraw,atkey,twin);
			delete pwdataraw;
			x1=ExtractComponent(*pwdata,0);
			x2=ExtractComponent(*pwdata,1);
			x3=ExtractComponent(*pwdata,2);
			Stack s1(*x1,twin);
			Stack s2(*x2,twin);
			Stack s3(*x3,twin);
			stack3c.clear();
			stack3c.push_back(s1.stack);
			stack3c.push_back(s2.stack);
			stack3c.push_back(s3.stack);
			ThreeComponentSeismogram result(stack3c);
			double stalat=result.get_double("site.lat");
			double stalon=result.get_double("site.lon");
			double staelev=result.get_double("site.elev");
			double ttime=hcen.phasetime(stalat,stalon,staelev,
							phase);
			result.put("dir",dir);
			result.put("dfile",dfile);
			/* result is in relative time.  This makes 
			result arrival time reference = predicted time  for
			this phase */
			double atime=hcen.time+ttime;
			result.rtoa(atime);
			result.put("arrival.time",atime);
			/* Save data to both wfdisc and wfprocess */
			int irec;
			irec=dbsave(result,dbwfdisc.db,string("wfdisc"),
				mdwfdisc,am,chanmap,true);
			irec=dbsave(result,dbwfprocess.db,string("wfprocess"),
				mdwfprocess,am);

			dbwfprocess.db.record=irec;
			int pwfid=dbwfprocess.get_int("pwfid");
			dbevlink.put("evid",evid);
			dbevlink.put("pwfid",pwfid);
			sta=result.get_string("sta");
			dbsclink.put("sta",sta);
			dbsclink.put("chan","3C");
			dbsclink.put("pwfid",pwfid);
			/* We also to have to save arrival and assoc rows */
			long int  arid=save_arrival(dbarrival,result,
					arrivalchan,phase);
			save_assoc(dbassoc,result,orid,arid,phase,hcen);

		}
	}
	catch (SeisppError serr)
	{
		serr.log_error();
	}
	catch (MetadataGetError mgerr)
	{
		mgerr.log_error();
	}
	catch (exception& e)
	{
		cerr << "stdlib excecption caught in main.  Error message: "
			<< e.what()<<endl;
	}
	catch (...)
	{
		cerr << "Something threw an unexpected exception"<<endl;
	}
}
