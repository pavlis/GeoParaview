#include <iostream>
#include <fstream>
#include <list>
#include <vector>

#include "perf.h"
#include "coords.h"
#include "seispp.h"
#include "dmatrix.h"
Time_Series *Read_SAC_ascii(string);
using namespace std;
// damped least squares solver companion to dgesvd
// similar to one in libgenloc
//
void lminverse_solver(dmatrix &U,
	double *s,
	dmatrix &Vt,
	double *b,
	int m,
	int n,
	double damp,
	double *x)
{
	int j;
	double *work=new double[n];
	double smax;
	int ns;

	ns=min(m,n);

	// Use relative damping so search for max svalue
	for(smax=0.0,j=0;j<ns;++j)
		if(s[j]>smax) smax=s[j];
	damp *=smax;

	for(j=0;j<ns;++j)
	{
		x[j]=ddot(m,U.get_address(0,j),1,b,1);
		x[j] *= s[j]/(s[j]*s[j] + damp*damp);
	}
	dcopy(ns,x,1,work,1);
	// This blindly assumes Vt is nxn exactly with singular
	// vectors in the rows  
	for(j=0;j<ns;++j)
		x[j] = ddot(n,Vt.get_address(j,0),n,work,1);
	delete [] work;
}

int main(int argc, char **argv)
{

  	string listfile;
	
	list<string> zfiles, rfiles;
	list<string>::const_iterator izf,irf;
	vector<Time_Series> z,r;
	vector<Time_Series>::iterator iz,ir;
	int i;
	Time_Series *s;
	Time_Series *sptr;
	double theta;
	int n,m;


	cout << "Enter file containing list of SAC ascii files for Z and R"<< endl;
	cin >> listfile;
	ifstream inlist(listfile.c_str(), ios::in);
	if(!inlist) 
	{
		cerr << "Cannot open file "<<listfile<<endl;
		exit(1);
	}
	// Read the list of files and put them into a list container
	do {
		string *zp=new string;
                string *rp=new string;
		string& zf=*zp, rf=*rp;
		inlist >> zf;
                if(inlist.eof()) break;
		inlist >> rf;
		// Note this reverses the list.  This is useful
		// for synthetics from herrmann's programs as I've
		// run them at this time.

		zfiles.push_front(zf);
		rfiles.push_front(rf);
	} while (! inlist.eof());
	
	for(izf=zfiles.begin(),i=0;izf!=zfiles.end();++izf,++i)
	{
		try{
			s=Read_SAC_ascii(*izf);
		} catch (SAC_data_error err)
		{
			err.log_error();
			cerr<< "Attempting to read " << *izf << endl;
			exit(1);
		}
		sptr = new Time_Series(*s);
                z.push_back(sptr);
		delete s;
	}
	for(irf=rfiles.begin(),i=0;irf!=rfiles.end();++irf,++i)
	{
		try{
			s=Read_SAC_ascii(*irf);
		} catch (SAC_data_error err)
		{
			err.log_error();
			cerr<< "Attempting to read " << *irf << endl;
			exit(1);
		}
                sptr = new Time_Series(*s);
                r.push_back(sptr);
                delete s;
	}		
	if(r.size() != z.size())
	{
		cerr << "Radial and vertical vectors are not equal"<<endl;
		cerr << "radial size="<<r.size();
		cerr << "vertical size="<<z.size();
		exit(1);
	}
	n=r.size();
	m=z[0].ns;
	cout << "Output matrix is of size "<< m<<" by " << n << endl;
	
	cout << "Enter theta angle(degrees) to rotate for ray coordinates" << endl;
	cin >> theta;
	theta = rad(theta);
	// scan for irregularies
	for(i=0;i<n;++i)
	{
		if(z[i].ns != m) 
		{
			cerr << "Wrong length of "<<z[i].ns
					<<"for vertical trace "<<i<<endl;
			exit(1);
		}
		if(r[i].ns != m) 
		{
			cerr << "Wrong trace length of "<<r[i].ns
					<<"for radial trace "<<i<<endl;
			exit(1);
		}
	}
	
	dmatrix A(m,n,0);
	double *ptr;
	double *zwork=new double[m];
	double *rwork=new double[m];
	double *rhs=new double[m];
	double wz,wr,damp;
	
	// rotate and load matrix
	wz = cos(theta);
	wr = sin(theta);
	for(i=0;i<n;++i)
	{
		dcopy(m,z[i].s,1,zwork,1);
		dcopy(m,r[i].s,1,rwork,1);
		dscal(m,wz,zwork,1);
		ptr=A.get_address(i,0);
		dcopy(m,zwork,1,ptr,1);
		daxpy(m,wr,rwork,1,ptr,1);
	}
	// put the radial in the right hand side.  We'll use rwork for that
	dcopy(m,z[0].s,1,zwork,1);
	dcopy(m,r[0].s,1,rwork,1);
	wz=-sin(theta);
	wr=cos(theta);
	dscal(m,wz,zwork,1);
	dcopy(m,zwork,1,rhs,1);
	daxpy(m,wr,rwork,1,rhs,1);
	// a very slow way to solve this for now, but bombproof
	double *svalue=new double[n];
	double *x=new double[n];
	dmatrix U(m,m,0);
	dmatrix Vt(n,n,0);
	int info;
	string test;
	do {
		string test;
		dmatrix& uref=U,vref=Vt; // kind of odd, but I see now
					// way around using these refs
		cout << "enter damping constant" << endl;
		cin >> damp;
		dgesvd('o','a',m,n,A.get_address(0,0),m,svalue,
			U.get_address(0,0),m,Vt.get_address(0,0),n,&info);
		lminverse_solver(uref,svalue,vref,rhs,m,n,damp,x);
		for(i=0;i<n;++i) cout << x[i] << endl;
		cout<<"Try again?  (y or n):";
		cin >> test;
	} while(test=="y");
}
	

