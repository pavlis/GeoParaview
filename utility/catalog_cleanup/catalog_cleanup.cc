#include <stdlib.h>
#include <string.h>
#include "coords.h"
#include "dbpp.h"
#include "Hypocenter.h"
using namespace std;
using namespace SEISPP;
void usage()
{
	cerr << "catalog_clean dbin dbout tolerance\n";
	exit(-1);
}
Hypocenter load_hypo(Dbptr db)
{
	string method("tttaup");
	string model("iasp91");
	double lat,lon,depth,otime;
	dbgetv(db,0,"lat",&lat,
		"lon",&lon,
		"depth",&depth,
		"time",&otime,0);
	return(Hypocenter(rad(lat),rad(lon),depth,otime,method,model));
}

double space_time_difference(Hypocenter& h1, Hypocenter& h2)
{
	const double velocity(6.0);
	const double Re(6371.0);
	double dist=h1.distance(h2.lat,h2.lon);
	dist*=Re;
	double sqrdist=dist*dist;
	double dz,dt;
	dz=fabs(h1.z-h2.z);
	dt=fabs(h1.time-h2.time);
	dt*=velocity;
	sqrdist += dz*dz + dt*dt;
	return(sqrt(sqrdist));
}
	
	
int main(int argc, char **argv)
{
	Hypocenter currenthypo,previoushypo;
	char dbrow[256];
	int orid;
	if(argc!=4)usage();
	string dbin(argv[1]);
	string dbout(argv[2]);
	double tolerance=atof(argv[3]);
	cout << "catalog_clean:  Writing to output db="<<dbout<<endl
		<<"Using input from origin table of db="<<dbin
		<< "Space time tolerance value="<<tolerance<<endl;
	int nrecin=0,nrecout=0;
	DatascopeHandle dbhi(dbin,true);
	DatascopeHandle dbho(dbout,false);
	DatascopeHandle dbhoe(dbho);
	dbho.lookup("origin");
	dbhoe.lookup("event");
	dbhi.lookup("origin");
	list<string> sortkeys;
	sortkeys.push_back("time");
	dbhi.sort(sortkeys);
	dbhi.rewind();
	if(dbhi.number_tuples()<=0)
	{
		cerr << "Origin table in input database has no rows"<<endl;
		usage();
	}
	dbhi.rewind();
	Dbptr parent;
	bool firstpass;
	for(nrecin=0,firstpass=true;nrecin<dbhi.number_tuples();
			++dbhi,++nrecin)
	{
		int irec;
		currenthypo=load_hypo(dbhi.db);
		if(firstpass || space_time_difference(currenthypo,
				previoushypo)>tolerance)
		{
			/* because dbhi.db is a view we have to fetch
			the parent table this way to allow cloning with
			dbget/dbadd */
			dbgetv(dbhi.db,0,"origin",&parent,0);
			dbget(parent,dbrow);
			irec=dbadd(dbho.db,dbrow);
			dbgetv(dbhi.db,0,"orid",&orid,0);
			dbaddv(dbhoe.db,0,"evid",orid,"prefor",orid,0);
			dbho.db.record=irec;
			dbputv(dbho.db,0,"evid",orid,0);
			previoushypo=currenthypo;
			++nrecout;
			if(firstpass) firstpass=false;
		}
	}
	cout << "Input db had "<<nrecin<<" rows"<<endl
		<<"Wrote "<<nrecout<<" origin and event rows to output"<<endl;
}
