#include <iostream>
#include <fstream>
#include <string>
#include "stock.h"
#include "pf.h"
#include "db.h"
#include "vmodel.h"
#include "ray1d.h"
int main(int argc, char **argv)
{

	double zmax,dz,p;
	Pf *pf;
	char *pfname="time_to_depth";
	Dbptr db;
	char *dbdir,*dbname,*modname;
	RayPathSphere *Ppath, *Spath;
	int i;
	const double R0=6378.17;
	double dtbin;
        double ttest;

	if(pfread(pfname,&pf))
	{
		cerr << "pfread of "<<pfname<<" failed\n";
		exit(-1);
	}
	dbdir=pfget_string(pf,"velocity_model_directory");
	dbname=pfget_string(pf,"velocity_model_database");
	modname=pfget_string(pf,"velocity_model_name");
	if(dbname==NULL || modname==NULL)
	{
		cerr << "parameter file error"<<endl;
		exit(-1);
	}
	string s1,s2,fullpath,mod;
	s1=dbdir;  s2=dbname;
	fullpath=s1+"/"+dbname;
	if(dbopen((char *)fullpath.c_str(),"r",&db))
	{
		cerr<<"dbopen failure for model database "<<fullpath<<endl;
		exit(-1);
	}
	cout << "Reading P and S velocities for model "<<modname<<endl;
	mod=modname;
	try
	{
		Velocity_Model Vp(db,mod,"P");
		Velocity_Model Vs(db,mod,"S");

		cout << "Enter maximum depth in km" << endl;
		cin >> zmax;
		cout << "Enter desired depth step (km):";
		cin >> dz;
		cout << "Enter ray parameter (s/km): ";
		cin >> p;
		Ppath = Trace_Ray_Sphere(Vp,p,zmax,1000.0,dz,"z");
		Spath = Trace_Ray_Sphere(Vs,p,zmax,1000.0,dz,"z");
		cout << "P path"<<endl;
		for(i=0;i<Ppath->npts;++i)
			cout << Ppath->delta[i]<<","<<Ppath->r[i]<<","<<Ppath->t[i]<<endl;
		cout << "S path"<<endl;
		for(i=0;i<Spath->npts;++i)
			cout << Spath->delta[i]<<","<<Spath->r[i]<<","<<Spath->t[i]<<endl;
		double *dtsmp=new double[Spath->npts];
		if(Ppath->npts != Spath->npts)
		{
			cout << "ray mismatch in points: "<<Ppath->npts
				<< "points for P ray and " <<Spath->npts
				<< "for for S ray" << endl;
			exit(-1);
		}
		for(i=0;i<Ppath->npts;++i) dtsmp[i]=Spath->t[i]-Ppath->t[i];
                cout << "depth P-S time curve"<<endl;
		for(i=0;i<Ppath->npts;++i) 
		{
			double delta;
			delta=Ppath->delta[i]-Spath->delta[i];
			delta *= R0;
			dtsmp[i]+=p*delta;
			cout << dz*((double)i)<<","<<dtsmp[i]<<endl;
		}
		cout << "Enter desired time sample interval"<<endl;
		cin >> dtbin;
                i=0;
		ofstream lagsout("lags.out",ios::out);
		ofstream angles("angles.out",ios::out);
                cout << "S-P time lag depth curve" << endl;
                do
		{
			ttest=dtbin*((double)i);
                        for(int j=0;j<Ppath->npts;++j)
			{
				double depth,vel,theta;
				if(dtsmp[j]>=ttest)
				{
					depth=R0-Ppath->r[j];
					vel = Vp.getv(depth);
					theta=180*asin(p*vel)/M_PI;
					cout << ttest <<" "<<dtsmp[j]<<" "
						<< depth << endl;
					angles << theta << endl;
					lagsout << Spath->t[j] << endl;
					break;
				}
			}
                        ++i;
		}while(ttest<(R0-Ppath->r[Ppath->npts]));
	} 
	catch(Velocity_Model_error(vme))
	{
		vme.log_error();
		exit(-1);
	}

}
