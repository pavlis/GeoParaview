#include <iostream>
#include <math.h>
#include "perf.h"
#include "dmatrix.h"

using namespace std;
dmatrix::dmatrix(int nr, int nc)
{
  nrr=nr;
  ncc=nc;
  length=nrr*ncc;
  if(length<1)
      {
      length=1;
      nrr=ncc;
      }
  ary=new double[length];
}

dmatrix::dmatrix(dmatrix& other)
  {
  nrr=other.nrr;
  ncc=other.ncc;
  length=other.length;
  ary=new double[length];
  memcpy(ary,other.ary, length*sizeof(double));
  }

dmatrix::~dmatrix()
{
delete ary;
}

double &dmatrix::operator()(int rowindex, int colindex)
{
  int out_of_range=0;
  if (rowindex>=nrr) out_of_range=1;
  if (rowindex<0) out_of_range=1;
  if (colindex>=ncc) out_of_range=1;
  if (colindex<0) out_of_range=1;
  if (out_of_range)
	throw dmatrix_index_error(nrr,ncc,rowindex,colindex);
  return (ary[rowindex+colindex*nrr]);
}
//
// subtle difference here.  This one returns a pointer to the 
// requested element
//
double* dmatrix::get_address(int rowindex, int colindex)
{
  double *ptr;
  int out_of_range=0;
  if (rowindex>=nrr) out_of_range=1;
  if (rowindex<0) out_of_range=1;
  if (colindex>=ncc) out_of_range=1;
  if (colindex<0) out_of_range=1;
  if (out_of_range)
        throw dmatrix_index_error(nrr,ncc,rowindex,colindex);
  ptr = ary + rowindex + nrr*colindex;
  return(ptr);
}

void dmatrix::operator=(dmatrix& other)
{
if(&other==this) return;
ncc=other.ncc;
nrr=other.nrr;
length=other.length;
delete ary;
ary= new double[length];
memcpy(ary,other.ary, length*sizeof(double));
}

void dmatrix::operator+=(dmatrix& other)
 {
int i;
  if ((nrr!=other.nrr)||(length!=other.length))
	throw dmatrix_size_error(nrr, ncc, other.nrr, other.length);
for(i=0;i<length;i++)
  ary[i]+=other.ary[i];
 }

void dmatrix::operator-=(dmatrix& other)
 {
int i;
  if ((nrr!=other.nrr)||(length!=other.length))
	throw dmatrix_size_error(nrr, ncc, other.nrr, other.length);
for(i=0;i<length;i++)
  ary[i]-=other.ary[i];
 }

dmatrix operator+(dmatrix &x1, dmatrix &x2)
  {
int i;
  if ((x1.nrr!=x2.nrr)||(x1.length!=x2.length))
	throw dmatrix_size_error(x1.nrr, x1.ncc, x2.nrr, x2.length);
 dmatrix tempmat(x1.nrr,x1.ncc);
  for(i=0;i<x1.length;i++) tempmat.ary[i]=x1.ary[i]+x2.ary[i];
return tempmat;
}

dmatrix operator-(dmatrix &x1, dmatrix &x2)
  {
int i;
  if ((x1.nrr!=x2.nrr)||(x1.length!=x2.length))
	throw dmatrix_size_error(x1.nrr, x1.ncc, x2.nrr, x2.length);
  dmatrix tempmat(x1.nrr,x1.ncc);
  for(i=0;i<x1.length;i++) tempmat.ary[i]=x1.ary[i]-x2.ary[i];
return tempmat;
}

//dmatrix operator*(dmatrix& x1,dmatrix& b)
dmatrix operator*(dmatrix& x1,dmatrix& b)
	{
	int i,j,k;
	double xval,bval;
	if(x1.ncc!=b.nrr)
		throw dmatrix_size_error(x1.nrr, x1.ncc, b.nrr, b.length);
	dmatrix prod(x1.nrr,b.ncc);
	for(i=0;i<x1.nrr;i++)
	{
		for(j=0;j<b.ncc;j++)
		{
			prod(i,j)=ddot(x1.ncc,x1.get_address(i,0),x1.nrr,
				b.get_address(0,j),1);
		}
	}
	return prod;
	}

dmatrix operator*(double& x, dmatrix &zx)
  {
  int i;
  dmatrix tempmat(zx.nrr,zx.ncc);
  for(i=0;i<zx.length;i++) tempmat.ary[i]=x*zx.ary[i];
  return tempmat;
  }

dmatrix operator/(dmatrix &zx, double& x)
  {
int i;
  dmatrix tempmat(zx.nrr,zx.ncc);
  for(i=0;i<zx.length;i++) tempmat.ary[i]=zx.ary[i]/x;
return tempmat;
  }  


dmatrix inv(dmatrix& x1)
{
 int i,icol,irow,j,k,jq,iq,*indxc,*indxr,*ipiv;
 double big,dum,pivinv,thold;

if(x1.nrr>x1.ncc)       //for pseudoinverses
  {
  // this produces a small memory leak
  dmatrix x1T(1,1);
  x1T=tr(x1);
  return (inv(x1T*x1))*x1T;
  }
if(x1.nrr<x1.ncc)   // underdetermined crude pseudoinverse
  {
  dmatrix x1T(1,1);
  x1T=tr(x1);
  return (x1T*(inv(x1*x1T)));
  }
 dmatrix a(x1.ncc,x1.ncc);
 a=x1;
 a.ncc=a.nrr=x1.ncc;
 indxc=new int[1+a.nrr];
 indxr=new int[1+a.nrr];
 ipiv=new int[1+a.nrr];
 for (j=1;j<=a.nrr;j++) ipiv[j]=0;
 for (i=1;i<=a.nrr;i++)
   {
    big=0.0;
    for (j=1;j<=a.nrr;j++)
      if (ipiv[j] !=1)
	 for (k=1;k<=a.nrr;k++)
	   {
	   if (ipiv[k]==0)
	    {
	    if (fabs(a(j,k))>=big)
	      {big=fabs(a(j,k));
	       irow=j;
	       icol=k;
	       }
	   } else if (ipiv[k]>1)
		throw dmatrix_inv_error(x1.ncc, x1.nrr);
}
 ++(ipiv[icol]);

 if (irow!=icol)
    {
    for (jq=1;jq<=a.nrr;jq++)   //interchange a(irow,jq) and a(icol,jq)
	{
	thold=a(irow,jq);
	a(irow,jq)=a(icol,jq);
	a(icol,jq)=thold;
	}
    }
 indxr[i]=irow;
 indxc[i]=icol;
 if (a(icol,icol)==0.0)
	throw dmatrix_inv_error(x1.ncc, x1.nrr);
 pivinv=1.0/a(icol,icol);
 a(icol,icol)=1.0;
 for (jq=1;jq<=a.nrr;jq++) a(icol,jq) *= pivinv;
 for (iq=1;iq<=a.nrr;iq++)
     if(iq!=icol)
	{
	dum=a(iq,icol);
	a(iq,icol)=0.0;
	for (jq=1;jq<=a.nrr;jq++) a(iq,jq) -= a(icol,jq)*dum;
	}
 }
 for (jq=a.nrr;jq>=1;jq--)
    {
     if (indxr[jq]!=indxc[jq])
	 for(k=1;k<=a.nrr;k++)
	   {             //interchange a(k,indxr[jq]) and a(k,indxc[jq])
	   thold=a(k,indxr[jq]);
	   a(k,indxr[jq])=a(k,indxc[jq]);
	   a(k,indxc[jq])=thold;
	   }
    }

delete ipiv;
delete indxr;
delete indxc;
a.nrr=x1.nrr;
a.ncc=x1.ncc;
return a;
 }

dmatrix tr(dmatrix& x1)
{
int i,j;
dmatrix temp(x1.ncc,x1.nrr);
for(i=0; i<x1.nrr; i++)
   for(j=0; j<x1.ncc;j++)
	temp(j,i)=x1(i,j);
return temp;
}


ostream& operator<<(ostream& os, dmatrix& x1)
  {
  int i,j;
  for(i=0;i<x1.nrr;i++)
  {
  for(j=0;j<x1.ncc;j++) os << x1(i,j) <<" ";
  os<<"\n";
  }
  return os;
  }

istream& operator>>(istream& is, dmatrix& x1)
  {
  int i,j;
  for(i=0;i<x1.nrr;i++)
  {
  for(j=0;j<x1.ncc;j++) is >> x1(i,j);
  }
  return is;
  }

void dmatrix::zero()
{
	for(int j=0;j<length;++j) ary[j]=0.0;
}
