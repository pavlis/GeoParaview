#include <iostream>
#include <fstream>
#include <list>
#include <vector>
#include <sunmath.h>

#include "perf.h"
#include "coords.h"
#include "seispp.h"
#include "dmatrix.h"
Time_Series *Read_SAC_ascii(string);
using namespace std;
// damped least squares solver companion to dgesvd
// similar to one in libgenloc
//
/*  This funciton is currently commented out because the program
was seg faulting and I decided it was easier to just use matlab
for what it did.

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
cout <<"damping constant="<<damp<<endl;
cout <<"largest singular value="<<smax<<endl;
cout <<"smallest singular value="<<s[n-1]<<endl;
	dcopy(ns,x,1,work,1);

	for(j=0;j<ns;++j)
	{
		work[j]=ddot(m,U.get_address(0,j),1,b,1);
		work[j] *= s[j]/(s[j]*s[j] + damp*damp);
	}
	// This blindly assumes Vt is nxn exactly with singular
	// vectors in the rows  
	for(j=0;j<n;++j) x[j]=0.0;
	for(j=0;j<ns;++j)
		daxpy(n,work[j],Vt.get_address(j,0),n,x,1);
	delete [] work;
}
*/

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
	int n,m,mprime;


	cout << "Enter file containing list of SAC ascii files for Z "<< endl;
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
		inlist >> *zp;
                if(inlist.eof()) break;
		zfiles.push_back(*zp);
	} while (! inlist.eof());
	inlist.close();
	cout << "Found " << zfiles.size() << " files" << endl;

	cout << "Enter file containing list of SAC ascii files for R "<< endl;
	cin >> listfile;
	inlist.open(listfile.c_str(), ios::in);
	if(!inlist) 
	{
		cerr << "Cannot open file "<<listfile<<endl;
		exit(1);
	}
	do {
                string *rp=new string;
		inlist >> *rp;
                if(inlist.eof()) break;
		rfiles.push_back(*rp);
	} while (! inlist.eof());
	inlist.close();
	cout << "Found " << rfiles.size() << " files" << endl;
	if(rfiles.size()!=zfiles.size())
	{
		cerr<<"Abort:  z and r lists are not equal size\n";
		exit(-1);
	}
	cout <<"Enter file containing time lags in seconds"<<endl;
	cin>>listfile;
	inlist.open(listfile.c_str(), ios::in);
	if(!inlist) 
	{
		cerr << "Cannot open file "<<listfile<<endl;
		exit(1);
	}
	cout << "Enter file containing P wave polarization angles" << endl;
	cin>>listfile;
	ifstream anglein(listfile.c_str(), ios::in);
	if(!anglein)
	{
		cerr << "Cannot open file "<<listfile<<endl;
                exit(1);
        }

	n=rfiles.size();
	double *lags=new double[n];
	double *theta=new double[n];
	for(i=0;i<n;++i)
	{
		inlist >> lags[i];
		anglein >> theta[i];
		theta[i]=rad(theta[i]);
		if(inlist.eof() || anglein.eof())
		{
			cerr << "EOF encountered at line " << i <<endl;
			cerr << "Expected to find " << n << "lines in file"<<endl;
			exit(-1);
		}
	}
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
cout << "done reading verticals"<<endl;
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
	m=z[0].ns+3*n;
	mprime=z[0].ns;
	cout << "Output matrix is of size "<< m<<" by " << n << endl;
	
	// scan for irregularies
	for(i=0;i<n;++i)
	{
		if(z[i].ns != z[0].ns) 
		{
			cerr << "WARNING:  Wrong length of "<<z[i].ns
					<<"for vertical trace "<<i<<endl;
		}
		if(r[i].ns != r[0].ns) 
		{
			cerr << "WARNING:  Wrong trace length of "<<r[i].ns
					<<"for radial trace "<<i<<endl;
		}
	}
	
	dmatrix A(m,n);
	dmatrix Az(m,n);  // conventional z convolution matrix
	A.zero();
	Az.zero();
	double *ptr;
/* debug */
/*
for(i=0;i<n;++i)
{
ptr=A.get_address(0,i);
dcopy(mprime,z[i].s,1,ptr,1);
}
ofstream mout("amatrix.out");
mout << A;
mout.close();
*/
	
	double *zwork=new double[m];
	double *rwork=new double[m];
	double *rhs=new double[m];
	double wz,wr,damp;
	
	for(i=0;i<n;++i)
	{
		int i0;
		int m_to_copy;
		wz = cos(theta[i]);
		wr = sin(theta[i]);
		dcopy(mprime,z[i].s,1,zwork,1);
		dcopy(mprime,r[i].s,1,rwork,1);
// debug
/*
if(i==0)
for(int j=0;j<m;++j) cout << zwork[j] << "," << rwork[j] << endl;
*/
		dscal(mprime,wz,zwork,1);
		i0=nint((lags[i]-lags[0])/z[i].dt);
cout << "lag(s)="<<lags[i]<<" yields i0=" << i0 << endl;
// temporary to fudge
//i0-=i;
//cout << i0<<endl;
		ptr=A.get_address(i0,i);
		if(i0+mprime<m)
			m_to_copy=mprime;
		else
		{
			m_to_copy = m-i0;
			cerr << "warning:  column "<<i<<" only copied "
				<< m_to_copy << " of " << mprime 
				<<" rows"<<endl;
		}
		dcopy(m_to_copy,zwork,1,ptr,1);
		daxpy(m_to_copy,wr,rwork,1,ptr,1);
		// This loads he surface l component as convolution
		dcopy(mprime,A.get_address(0,0),1,Az.get_address(i,i),1);
	}
	// put the radial in the right hand side.  We'll use rwork for that
	dcopy(mprime,z[0].s,1,zwork,1);
	dcopy(mprime,r[0].s,1,rwork,1);
ofstream mout("amatrix.out");
mout << A;
mout.close();
mout.open("azmatrix.out");
mout << Az;
mout.close();
	wz=-sin(theta[0]);
	wr=cos(theta[0]);
	dscal(mprime,wz,zwork,1);
	dcopy(mprime,zwork,1,rhs,1);
	daxpy(mprime,wr,rwork,1,rhs,1);
	for(i=mprime;i<m;++i) rhs[i]=0.0;  // this isn't initialized without this
	ofstream rhsout("rhs.out",ios::out);
	for(i=0;i<m;++i) rhsout <<rhs[i] << ","<<A(i,0)<<","<<r[0].s[i]<<","<<z[0].s[i]<<endl;
	rhsout.close();
/*

	// a very slow way to solve this for now, but bombproof
	double *svalue=new double[n];
	double *x=new double[n];
	dmatrix U(m,m);  // placeholder because I use overwrite mode for dgdsvd
	dmatrix Vt(n,n);
	dmatrix& uref=U,vref=Vt; // kind of odd, but I see no
				// way around using these refs
	int info;
	char test[10];
	cout << "Enter s for syndecon solution or c for convolution solution"<<endl;
	cin >> test;
	if(!strcmp(test,"c"))
	{
		cout << "Computing conventional convolutions solution" << endl;
		dgesvd('a','a',m,n,Az.get_address(0,0),m,svalue,
                        U.get_address(0,0),m,Vt.get_address(0,0),n,&info);
	}
	else
	{
		cout << "Computing matrix operator solution" << endl;
		dgesvd('a','a',m,n,A.get_address(0,0),m,svalue,
			U.get_address(0,0),m,Vt.get_address(0,0),n,&info);
	}
	if(info!=0) cerr << "dgesvd returned error code " << info<<endl;
	cout << "singular values"<<endl;
	for(i=0;i<n;++i)cout <<svalue[i]<<endl;
	do {
		cout << "enter damping constant" << endl;
		cin >> damp;
		lminverse_solver(uref,svalue,vref,rhs,m,n,damp,x);
		ofstream sout("syndecon.out",ios::out);
		for(i=0;i<n;++i) sout << x[i] << endl;
		sout.close();
		cout<<"Try again?  (y or n):";
		cin >> test;
	} while(strcmp(test,"n"));
*/
	return(0);
}
	

