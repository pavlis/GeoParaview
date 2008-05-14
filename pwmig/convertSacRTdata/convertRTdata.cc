/* Crude program to beat SAC radial receiver function into input format
for pwmig.  That means 3c data objects.  See design file in this directory
for more info. */
#include <math.h>
#include <sstream>
#include "coords.h"
#include "stock.h"
#include "seispp.h"
#include "dbpp.h"
#include "ThreeComponentSeismogram.h"
#include "TimeSeries.h"
#include "Hypocenter.h"
using namespace std;
using namespace SEISPP;
void usage()
{
	cerr << "convertRTdata dbin dbout [-pf pffile]"<<endl;
	exit(-1);
}
bool SEISPP::SEISPP_verbose(true);
int main(int argc, char **argv)
{
	if(argc<3) usage();
	string dbin(argv[1]);
	string dbout(argv[2]);
	const string outtable("wfprocess");
	cout << "convertRTdata:  read from db="<<dbin<<" writing results to db="<<dbout<<endl;
	string pffile("convertRTdata");
	int i;
	for(i=3;i<argc;++i)
	{
		string argstr(argv[i]);
		if(argstr=="-pf") 
		{
			++i;
			if(i>=argc) usage();
			pffile=string(argv[i]);
		}
		else
			usage();
	}
	Pf *pf;
	if(pfread(const_cast<char *>(pffile.c_str()),&pf))
	{
		cerr << "pfread error"<<endl;
		usage();
	}
	MetadataList mdl=pfget_mdlist(pf,"trace_mdl");
	MetadataList mdlout=pfget_mdlist(pf,"output_mdl");
	
		
	try {
		AttributeMap am("css3.0");
		DatascopeHandle dbhi(dbin,true);
		/* This set of db manipulations is a bit odd, 
		but seems to work with these data. */
		/*
		dbhi.lookup("wfdisc");
		dbhi.natural_join("arrival");
		string sstr("iphase=~/P/");
		dbhi.subset(sstr);
		dbhi.natural_join("site");
		list<string> jk;
		jk.push_back("arid");
		dbhi.leftjoin("assoc",jk,jk);
		jk.clear();
		jk.push_back("orid");
		dbhi.leftjoin("origin",jk,jk);
		jk.clear();
		jk.push_back("evid");
		dbhi.leftjoin("event",jk,jk);
		sstr="orid==prefor";
		dbhi.subset(sstr);
		*/

		dbhi.lookup("event");
		dbhi.natural_join("origin");
		dbhi.natural_join("assoc");
		dbhi.natural_join("arrival");
		string sstr("iphase=~/P/");
		dbhi.subset(sstr);
		list<string> jk1,jk2;
		jk1.push_back(string("arrival.time"));
		jk1.push_back(string("sta"));
		jk2.push_back(string("wfdisc.time::wfdisc.endtime"));
		jk2.push_back(string("sta"));
		//dbhi.join("wfdisc",jk1,jk2);
		//dbhi.natural_join("wfdisc");
		dbhi.leftjoin("wfdisc",jk2,jk1);
		dbhi.natural_join("site");

		cout << "Working input view has "<<dbhi.number_tuples()
			<<" rows"<<endl;
		DatascopeHandle dbho(dbout,false);
		dbho.lookup("wfprocess");
		/* Need to put stuff in evlink and sclink for pwmig */
		DatascopeHandle dbhosc(dbho);
		dbhosc.lookup("sclink");
		DatascopeHandle dbhoev(dbho);
		dbhoev.lookup("evlink");
		dbhi.rewind();
		int j,k;
		double slat,slon,sdepth,stime;
		string method("tttaup"),model("iasp91");
		double stalat,stalon,staelev;
		double seaz,az;
		double a,b;
		double atime;
		int evid,pwfid;
		string sta;
		/* frozen for now */
		string outdir("wf3c");
		string dfilebase("RFdata");
		for(i=0;i<dbhi.number_tuples();++i,++dbhi)
		{
			TimeSeries d(dynamic_cast<DatabaseHandle&>(dbhi),mdl,am);
			/* We need to compute an equivalent transformation matrix for
			these data from the hypocenter information.  We'll let this exit
			with an exception if origin information is not present. */
			slat=d.get_double("source_lat");
			slon=d.get_double("source_lon");
			evid=d.get_int("evid");
			cout << "Processing event "<<evid
				<< " located at (lat,lon)=("<<slat
				<<", "<<slon<<") "<<endl;

			slat=rad(slat);
			slon=rad(slon);
			sdepth=d.get_double("source_depth");
			stime=d.get_double("source_time");
			Hypocenter h(slat,slon,sdepth,stime,method,model);
			sta=d.get_string("sta");
			stalat=d.get_double("sta_lat");
			stalon=d.get_double("sta_lon");
			staelev=d.get_double("sta_elev");
			stalat=rad(stalat);
			stalon=rad(stalon);
			/* We have to fudge the start time because
			Fenglin's rf code seems to not keep absolute
			time.  Well, it might be in sac2db, but in
			any case the arrival.time of P is the same
			as the origin time.  Thus, to fix this start 
			time we find the time difference between atime
			and start time of the trace, compute P travel time,
			and then correct t0 to make P time match theoretical
			time exactly. */
			atime=d.get_double("atime");
			double ptime=h.ptime(stalat,stalon,staelev);	
			double oldt0=d.t0;
			double pdt=atime-oldt0;
			cout << sta << " original atime - t0 ="<<pdt<<endl;
			cout << "origin time and atime difference="<< stime-atime<<endl;
			d.t0=stime+ptime-pdt;
			cout << "t0 of this seismogram changed by "
				<< d.t0-oldt0 
				<< " seconds"<<endl;
			
			seaz=h.seaz(stalat,stalon);
			/* radial is defined as propagation direction, but seaz is back azimuth */
			az=seaz+M_PI;
			if(az>(2*M_PI)) az-=M_PI;
			cout <<" station="<<sta
				<<" computed radial azimuth="<<deg(az)<<endl;
			/* Now we have to set the transformation matrix from az.  Convention here
			is we put R into x2, so the rotation angle turns out to be -az when you 
			work through the difference between rotation from the x1 axis and an
			azimuth.  Hence, this is the correction tranformation matrix */
			a=cos(-az);
			b=sin(-az);
			/* by loading these in the metadata for d they will be cloned
			and used to construct the transformation matrix for d3c. */
			d.put("U11",a);
			d.put("U21",-b);
			d.put("U31",0.0);
			d.put("U12",b);
			d.put("U22",a);
			d.put("U32",0.0);
			d.put("U13",0.0);
			d.put("U23",0.0);
			d.put("U33",1.0);
			/* this uses the Metadata from d to create a template or 3c data.*/
			ThreeComponentSeismogram d3c(dynamic_cast<Metadata&>(d),false);
			/* Zero components 1 and 3 and put data from d into component 2*/
			for(k=0;k<d.ns;++k)
				for(j=0;j<3;++j) d3c.u(j,k)=0.0;
			for(k=0;k<d.ns;++k) d3c.u(1,k)=d.s[k];
			// always mark data live
			d3c.live=true;
			d3c.rotate_to_standard();
			/* We need to post these things or we're screwed */
			d3c.put("dir",outdir);
			char buf[128];
			stringstream ss(buf);
			ss<<dfilebase<<i;
			d3c.put("dfile",ss.str());
			d3c.put("wfprocess.algorithm","convertRF");
			d3c.put("timetype","a");
			int rec=dbsave(d3c,dbho.db,outtable,mdlout,am);
			if(rec<0)
			{
				cerr << "dbsave failed working on input record="<<i<<endl;
				exit(-1);
			}
			/* We have to set up these tables because we are using wfprocess
			to store 3c objects.*/
			dbho.db.record=rec;
			pwfid=dbho.get_int("pwfid");
			rec=dbhoev.append();
			dbhoev.put("pwfid",pwfid);
			dbhoev.put("evid",evid);
			rec=dbhosc.append();
			dbhosc.put("sta",sta);
			/* Ignore the fact this has endian implications.  At present datatype
			is used to handle this. */
			dbhosc.put("chan","3C");
			dbhosc.put("pwfid",pwfid);

		}
	
	}
	catch (SeisppError serr)
	{
		serr.log_error();
	}
	catch (MetadataGetError mderr)
	{
		mderr.log_error();
	}
	catch (...)
	{
		cerr << "Something threw an unexpected exception"<<endl;
	}
}
