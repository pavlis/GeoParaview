#include <map>
#include <sstream>
#include "StationVariableVelocityModel.h"
#include <string.h>
#include "stock.h"
#include "tt.h"
#include "coords.h"
#include "elog.h"
#include "seispp.h"
#include "dbpp.h"
#include "Hypocenter.h"

using namespace std;
using namespace SEISPP;

DatascopeHandle BuildCatalogView(DatascopeHandle& dbh)
{
	DatascopeHandle result(dbh);
	dbh.lookup("event");
	list<string> sortkeys;
	sortkeys.push_back("evid");
	dbh.sort(sortkeys);
	dbh.natural_join("origin");
	string ss_to_prefor("orid==prefor");
	dbh.subset(ss_to_prefor);
	dbh.natural_join("assoc");
	dbh.natural_join("arrival");
	dbh.natural_join("site");
	list<string> group_keys;
	group_keys.push_back("evid");
	dbh.group(group_keys);
	return(dbh);
}
/* Primary seismology routine in this program.  It computes the delay time,
commonly called tau, for model vmodel using input slowness slow.  This is computed by simple 
block integration from depth -elev to datum at interval dz.  Be aware this simple quadrature scheme
means dz should be made relatively small to avoid errors near velocity discontinuities. 
Note the integration always uses a velocity of the average at top and bottum of each
dz interval.  
*/
double ComputeDelayTime(VelocityModel_1d vmodel, double slow, 
		double elev,double dz, double datum)
{
	const string turningerror("ComputeDelayTimes:  ray became evanescent.");
	const string turnerr2(" Ray parameter passed is too large");
	vector <double> dtau;
	double z=-elev;
	double u;
	// if stations is below datum silently return 0
	if(z>datum) return(0.0);
	double vlast=vmodel.getv(z);
	double v,vbar,dt;
	dtau.push_back(0.0);
	z+=dz;
	while(z<datum)
	{
		v=vmodel.getv(z);
		vbar=(v+vlast)/2.0;
		u=1/vbar;
		if(u<slow) 
			throw SeisppError(turningerror+turnerr2);
		dt=dz*sqrt(u*u-slow*slow);
		dtau.push_back(dt);
		vlast=v;
		z+=dz;
	}
	/* handle the last point */
	if(z>=datum)
	{
		double ddz=datum-(z-dz);
		z=datum;
		vbar=(v+vlast)/2.0;
		u=1/vbar;
		if(u<slow) 
			throw SeisppError(turningerror+turnerr2);
		dt=ddz*sqrt(u*u-slow*slow);
		dtau.push_back(dt);
	}
	/*integrate by simple sum*/
	double result;
	int i;
	/* Intentionaly stated at 1 because first point is 
	always set 0.0*/
	for(i=1,result=0.0;i<dtau.size();++i) result+=dtau[i];
	return(result);
}
list<string> pfget_vdef(Pf *pf)
{
	list<string> vlist;
	Tbl *t;
	t=pfget_tbl(pf,"vmodel_map");
	int i;
	char *line;
	for(i=0;i<maxtbl(t);++i)
	{
		line=(char *)gettbl(t,i);
		vlist.push_back(string(line));
	}
	return(vlist);
}
void usage()
{
	cerr << "crustalcorrection db [-f -v -pf pffile]"<<endl;
	exit(-1);
}
bool SEISPP::SEISPP_verbose(false);
int main(int argc, char **argv)
{
	int i;
	if(argc<2) usage();
	string db(argv[1]);
	string pfname("crustalcorrection");
	bool unix_filter_mode(false);
	for(i=2;i<argc;++i)
	{
		string sarg(argv[i]);
		if(sarg=="-pf")
		{
			++i;
			pfname=string(argv[i]);
		}
		else if(sarg=="-f")
			unix_filter_mode=true;
		else if(sarg=="-v")
			SEISPP_verbose=true;
		else 
			usage();
	}
	Pf *pf;
	if(pfread(const_cast<char *>(pfname.c_str()), &pf) )
	{
		cerr << "crustal correction:  pfread error on pf="<<pfname<<endl;
		usage();
	}
	Metadata control(pf);
	/* The rest of this code is in a try block to catch exceptions from SEISPP routines */
	try {
		double dz=control.get_double("dz");
		double datum=control.get_double("datum");
		string reference_model_name=control.get_string("reference_model");
		string modtype=control.get_string("model_type"); 
		if( ! ((modtype=="P") || (modtype=="S") ) )
		{
			cerr << "Illegal model_type parameter = "<<modtype<<endl
				<< "Must be either P or S"<<endl;
			exit(-2);
		}
		string default_vname=control.get_string("default_velocity_model_name");
		/* mixed use of pf and metadata a bit confusing here. Done because
		the Tbl is easier to handle this way.  This provides the list of 
		phases to be processed.  Allows handling multiple phases in one pass. */
		Tbl *t;
		t=pfget_tbl(pf,"phase_list");
		if(t==NULL)
		{
			cerr << "Missing required parameter phase_list"<<endl;
			exit(-2);
		}
		list<string> phaselist;
		map<string,string> phasemodmap;
		for(i=0;i<maxtbl(t);++i)
		{
			char *ph;
			string thisphase,thismodtype;
			ph=(char *)gettbl(t,i);
			stringstream ss(ph);
			ss>>thisphase;
			ss>>thismodtype;
			phasemodmap[thisphase]=thismodtype;
		}
		/* We test this many times so best to avoid repeated
		calls to this method */
		map<string,string>::iterator pmmend=phasemodmap.end();
		/* We need to not apply this program to arrivals less than 
		some reasonable cutoff to avoid sqrt of negative numbers
		from turning rays.  Easier to do this than add such a check
		to the delay time calculation. */
		double deltamin=control.get_double("minimum_delta");
		/* We use a fixed model/method to recompute predicted arrivals
		times.  This guaranteeds the computed residuals are consistent
		with a particular model.  This is better than a pure correction.*/
		string method=control.get_string("TTmethod");
		string model=control.get_string("TTmodel");
		/* Now open the database.  */
		DatascopeHandle dbhin(db,false);
		DatascopeHandle dbhmod(dbhin);
		dbhmod.lookup("mod1d");
		VelocityModel_1d reference_model(dbhmod.db,reference_model_name,
			modtype);
		/* Load a list of velocity models assigned to each station */
		list<string> vdef=pfget_vdef(pf);
		StationVariableVelocityModel cmodels(dbhmod,vdef,default_vname,modtype);
		double tau,tauref;
		string phase;
		if(unix_filter_mode)
		{
			double delta,depth,slowness,elev;
			string sta;
			while(!cin.eof())
			{
				cin >> phase;
				cin >> sta;
				cin >> elev;
				cin >> delta;
				cin >> depth;
				if(delta<deltamin)
				{
					cerr << "Delta="<<delta
						<< " is less than minimum allowed="
						<< deltamin<<endl;
					continue;
				}
				if(phasemodmap.find(phase)==pmmend)
				{
					cerr << "phase "
						<< phase
						<< " is not on allowed phase list."
						<< "Fix parameter file to change allowed phase list"<<endl;
				}
				slowness=phase_slowness(const_cast<char *>(phase.c_str()),
									delta,depth);
				VelocityModel_1d vmod=cmodels.getmod(sta);
				tauref=ComputeDelayTime(reference_model,slowness,0.0,dz,datum);
				tau=ComputeDelayTime(vmod,slowness,elev,dz,datum);
				cout << phase << " "<<sta<<" "<<tau-tauref<<endl;
			}
		}
		else
		{
			/* This builds a standard catalog view sorted and grouped by evid.
			We then work through the group event by event.*/
			DatascopeHandle dbhcat=BuildCatalogView(dbhin);
			dbhcat.rewind();
			for(i=0;i<dbhcat.number_tuples();++i,++dbhcat)
			{
				DBBundle bundle=dbhcat.get_range();
				Dbptr dbp=bundle.parent;
				double lat,lon,depth,atime,time;
				double elev,slowness;
				double delta;
				double slowmag,predtime;
				double resid;
				int evid;
				char assocphase[10],csta[8];
				dbp.record=bundle.start_record;
				dbgetv(dbp,0,"origin.lat",&lat,
					"origin.lon",&lon,
					"origin.depth",&depth,
					"origin.time",&time,
					"evid",&evid,0);
				Hypocenter hypo(rad(lat),rad(lon),depth,time,method,model);
				Dbptr dbassoc; // need this to fetch the assoc row
				for(dbp.record=bundle.start_record;
					dbp.record<bundle.end_record;++dbp.record)
				{
					/* Using the straight datascope interface here instead
					of the seispp version for efficiency and because we are
					using the Dbptr as the loop counter.*/
					int iret;
					iret=dbgetv(dbp,0,"assoc",&dbassoc,
						"site.lat",&lat,
						"site.lon",&lon,
						"site.elev",&elev,
						"assoc.phase",assocphase,
						"assoc.sta",csta,
						"assoc.delta",&delta,
						"arrival.time",&atime,
						0);
					phase=string(assocphase);
					if(iret==dbINVALID)
					{
						cerr << "dbgetv error at row "
							<< dbp.record
							<< " of working db view"
							<< endl;
						elog_die(1,"Error from antelope");
					}
					if(phasemodmap.find(assocphase)==pmmend)
					{
					   if(SEISPP_verbose)
					   {
						cerr << "skipping data for phase="
							<< phase<<endl
							<<"Not in allowed phase list.  Check parameter file."<< endl;
					   }
					}
					else
					{
						if( (delta<deltamin)
							&& SEISPP_verbose)
						{
							cerr << "sta="
								<< csta
								<<" phase="
								<< assocphase
								<<" evid="
								<< evid<<endl
							<< "delta<deltamin residual calculation skipped"
							<<endl;

						}
						else
						{
							SlownessVector slowness=hypo.phaseslow(rad(lat),
								rad(lon), elev, phase);
							slowmag=slowness.mag();
							predtime=hypo.phasetime(rad(lat),
								rad(lon), elev, phase);
							string sta(csta);
							VelocityModel_1d vmod=cmodels.getmod(sta);
							// Note elev to reference is always set to 0 here 
							tauref=ComputeDelayTime(reference_model,slowness.mag(),0.0,dz,datum);
							tau=ComputeDelayTime(vmod,slowness.mag(),elev,dz,datum);
							resid=atime-hypo.time
								-predtime-(tau-tauref);
							string vmodel=model+":"+sta;
							dbputv(dbassoc,0,"timeres",resid,"vmodel",vmodel.c_str(),0);
							if(SEISPP_verbose)
							{
								cout <<  "sta="
                                                                << csta
                                                                <<" phase="
                                                                << assocphase
                                                                <<" evid="
                                                                << evid
								<< " dt="
								<< tau-tauref
								<<endl;
							}
						}
					}
				}
			}
		}
	
	} catch (SeisppError serr)
	{
		serr.log_error();
		exit(-2);
	}
}

