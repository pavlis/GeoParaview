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
/*  This funciton is currently commented out because the program
was seg faulting and I decided it was easier to just use matlab
for what it did.
*/

void lminverse_solver(dmatrix& U,
	double *s,
	dmatrix& Vt,
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


/* computes the delay time variable, tau, for a 1d velocity model vmod
for a source a depth zsource and a surface receiver.  (depth=0). 
The variable p is the slowness of the ray used to computed tau.

Integration uses a simple sum making this equivalent to reducing the
model to a set of constant velocity layers at the node points. The
variable dz defines the depth spacing of the grid that is effectively
generated to do this integration.  Hence, dz should normally be made
relatively small to avoid roundoff errors.
*/
double compute_tau(Velocity_Model_1d& vmod, 
	double zsource, 
		double p, 
			double dz) 
{
	double z=0.0,ddz;
	double psq,u;
	double tau=0.0;

	if(zsource==0.0) return(0.0);
	psq = p*p;

	do
	{
		u = 1.0/vmod.getv(z);
		if(u<p)throw(-1);
		if(zsource<=(z+dz))
			ddz = zsource-z;
		else
			ddz =dz;
		tau += sqrt(u*u-psq)*ddz;
		z += dz;
	} while(z<zsource);
	return(tau);
}
	
int main(int argc, char **argv)
{

  	string listfile,vmodfile;
	
	list<string> zfiles, rfiles;
	list<string>::const_iterator izf,irf;
	vector<Time_Series> z,r;
	vector<Time_Series>::iterator iz,ir;
	int i;
	Time_Series *s;
	Time_Series *sptr;
	int n,m,mprime;
	double p;  // ray parameter


	cout << "Enter model file used to generate synthetics"<< endl;
	cin >> vmodfile;
	Velocity_Model_1d vmod(vmodfile.c_str(),"rbh","P");
	Velocity_Model_1d vmodS(vmodfile.c_str(),"rbh","S");
	cout << "Enter ray parameter (s/km) used for synthetics"<<endl;
	cin >> p;
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
	ifstream inlistr(listfile.c_str(), ios::in);
	if(!inlistr) 
	{
		cerr << "Cannot open file "<<listfile<<endl;
		exit(1);
	}
	do {
                string *rp=new string;
		inlistr >> *rp;
                if(inlistr.eof()) break;
		rfiles.push_back(*rp);
	} while (! inlistr.eof());
	inlistr.close();
	cout << "Found " << rfiles.size() << " files" << endl;
	if(rfiles.size()!=zfiles.size())
	{
		cerr<<"Abort:  z and r lists are not equal size\n";
		exit(-1);
	}

	n=rfiles.size();

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
	double theta0;
	
	for(i=0;i<n;++i)
	{
		int lag;
		double tlag,tauP,tauS;
		int m_to_copy;
		double z0;
		const double dz=0.1;
		double theta,v;

		
		z0 = z[i].md.get_double("STEL");
		v = vmod.getv(z0);
		theta=asin(p*v);
		if(i==0)theta0=theta;  // assumes i=0 is surface
		

		wz = cos(theta);
		wr = sin(theta);
		dcopy(mprime,z[i].s,1,zwork,1);
		dcopy(mprime,r[i].s,1,rwork,1);
// debug
/*
if(i==0)
for(int j=0;j<m;++j) cout << zwork[j] << "," << rwork[j] << endl;
*/
		tauP = compute_tau(vmod,z0,p,dz);
		tauS = compute_tau(vmodS,z0,p,dz);
		tlag = tauS;
		lag = nint(tlag/z[i].dt);
// temporary to fudge
//i0-=i;
//cout << i0<<endl;
		ptr=A.get_address(i,i);
		if(lag+mprime<m)
			m_to_copy=mprime;
		else
		{
			m_to_copy = m-lag;
			cerr << "warning:  column "<<i<<" only copied "
				<< m_to_copy << " of " << mprime 
				<<" rows"<<endl;
		}
		// this loads and simultaneously rotates synthetics for this column
		dscal(mprime,wz,zwork,1);
		dcopy(m_to_copy,zwork+lag,1,ptr,1);
		daxpy(m_to_copy,wr,rwork,1,ptr+lag,1);
		// This loads the surface l component as convolution
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
	wz=-sin(theta0);
	wr=cos(theta0);
	dscal(mprime,wz,zwork,1);
	dcopy(mprime,zwork,1,rhs,1);
	daxpy(mprime,wr,rwork,1,rhs,1);
	for(i=mprime;i<m;++i) rhs[i]=0.0;  // this isn't initialized without this
	ofstream rhsout("rhs.out",ios::out);
	for(i=0;i<m;++i) rhsout <<rhs[i] << ","<<A(i,0)<<","<<r[0].s[i]<<","<<z[0].s[i]<<endl;
	rhsout.close();

	// a very slow way to solve this for now, but bombproof
	double *svalue=new double[n];
	double *x=new double[n];
	dmatrix *uptr,*vptr;
	uptr = new dmatrix(m,m);
	vptr = new dmatrix(n,n);
	dmatrix& U=*uptr;
	dmatrix& Vt=*vptr;
	
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
		lminverse_solver(U,svalue,Vt,rhs,m,n,damp,x);
		ofstream sout("syndecon.out",ios::out);
		for(i=0;i<n;++i) sout << x[i] << endl;
		sout.close();
		cout<<"Try again?  (y or n):";
		cin >> test;
	} while(strcmp(test,"n"));
	return(0);
}
	

