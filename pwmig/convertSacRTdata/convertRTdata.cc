/* Crude program to beat SAC radial receiver function into input format
for pwmig.  That means 3c data objects.  See design file in this directory
for more info. */
#include <math.h>
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
	for(i=0;i<argc;++i)
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
		dbhi.lookup("event");
		dbhi.natural_join("origin");
		dbhi.natural_join("assoc");
		dbhi.natural_join("arrival");
		dbhi.natural_join("wfdisc");
		dbhi.natural_join("site");
		DatascopeHandle dbho(dbout,false);
		dbho.lookup("wfprocess");
		/* Need to put stuff in evlink and sclink for pwmig */
		DatascopeHandle dbhosc(dbho);
		dbhosc.lookup("sclink");
		DatascopeHandle dbhoev(dbho);
		dbhosc.lookup("evlink");
		dbhi.rewind();
		int j,k;
		double slat,slon,sdepth,stime;
		string method("tttaup"),model("iasp91");
		double stalat,stalon;
		double seaz,az;
		double a,b;
		int evid,pwfid;
		string sta;
		for(i=0;i<dbhi.number_tuples();++i,++dbhi)
		{
			TimeSeries d(dynamic_cast<DatabaseHandle&>(dbhi),mdl,am);
			/* this uses the Metadata from d to create a template or 3c data.*/
			ThreeComponentSeismogram d3c(dynamic_cast<Metadata&>(d),false);
			/* Zero components 1 and 3 and put data from d into component 2*/
			for(k=0;k<d.ns;++k)
				for(j=0;j<3;++j) d3c.u(j,k)=0.0;
			for(k=0;k<d.ns;++k) d3c.u(1,k)=d.s[k];
			/* We need to compute an equivalent transformation matrix for
			these data from the hypocenter information.  We'll let this exit
			with an exception if origin information is not present. */
			slat=d.get_double("source_lat");
			slon=d.get_double("source_lon");
			evid=d.get_int("evid");
			cout << "Processing event "<<evid
				<< "located at (lat,lon)=("<<slat<<", "<<slon<<") ";
			slat=rad(slat);
			slon=rad(slon);
			sdepth=d.get_double("source_depth");
			stime=d.get_double("source_time");
			Hypocenter h(slat,slon,sdepth,stime,method,model);
			sta=d.get_string("sta");
			cout <<" station="<<sta;
			stalat=d.get_double("sta_lat");
			stalon=d.get_double("sta_lon");
			stalat=rad(stalat);
			stalon=rad(stalon);
			seaz=h.seaz(stalat,stalon);
			/* radial is defined as propagation direction, but seaz is back azimuth */
			az=seaz+M_PI;
			if(az>(2*M_PI)) az-=M_PI;
			cout <<" Computed radial azimuth"<<deg(az)<<endl;
			/* Now we have to set the transformation matrix from az.  Convention here
			is we put R into x2, so the rotation angle turns out to be -az when you 
			work through the difference between rotation from the x1 axis and an
			azimuth.  Hence, this is the correction tranformation matrix */
			a=cos(-az);
			b=sin(-az);
			d3c.tmatrix[0][0]=a;
			d3c.tmatrix[1][0]=-b;
			d3c.tmatrix[2][0]=0.0;
			d3c.tmatrix[0][1]=b;
			d3c.tmatrix[1][1]=a;
			d3c.tmatrix[2][1]=0.0;
			d3c.tmatrix[0][2]=0.0;
			d3c.tmatrix[1][2]=0.0;
			d3c.tmatrix[2][2]=1.0;
			d3c.rotate_to_standard();
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
