#include <iostream>
#include "dmatrix.h"

using namespace std;
dmatrix::dmatrix(int nr, int nc, int stval)
{
  nrr=nr;
  ncc=nc;
  matoff=stval;
  length=(nrr-matoff+1)*(ncc-matoff+1);
  if(length<1)
      {
      length=1;
      nrr=ncc=matoff;
      }
  ary=new double[length];
}

dmatrix::dmatrix(const dmatrix& other)
  {
  nrr=other.nrr;
  ncc=other.ncc;
  matoff=other.matoff;
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
  if (rowindex>nrr) out_of_range=1;
  if (rowindex<matoff) out_of_range=1;
  if (colindex>ncc) out_of_range=1;
  if (colindex<matoff) out_of_range=1;
  if (out_of_range)
	throw dmatrix_index_error(nrr,ncc,rowindex,colindex);
// old, stored in column order
//  return (ary[colindex-matoff+(ncc-matoff+1)*(rowindex-matoff)]);
  return (ary[rowindex-matoff+(nrr-matoff+1)*(colindex-matoff)]);
}
//
// subtle difference here.  This one returns a pointer to the 
// requested element
//
double* dmatrix::get_address(int rowindex, int colindex)
{
  double *ptr;
  int out_of_range=0;
  if (rowindex>nrr) out_of_range=1;
  if (rowindex<matoff) out_of_range=1;
  if (colindex>ncc) out_of_range=1;
  if (colindex<matoff) out_of_range=1;
  if (out_of_range)
        throw dmatrix_index_error(nrr,ncc,rowindex,colindex);
  ptr = ary + rowindex-matoff+(nrr-matoff+1)*(colindex-matoff);
  return(ptr);
}

void dmatrix::operator=(const dmatrix& other)
{
if(&other==this) return;
ncc=other.ncc;
nrr=other.nrr;
matoff=other.matoff;
length=other.length;
delete ary;
ary= new double[length];
memcpy(ary,other.ary, length*sizeof(double));
}

void dmatrix::operator+=(const dmatrix& other)
 {
int i;
  if ((nrr!=other.nrr)||(length!=other.length))
	throw dmatrix_size_error(nrr, ncc, other.nrr, other.length);
for(i=0;i<length;i++)
  ary[i]+=other.ary[i];
 }

void dmatrix::operator-=(const dmatrix& other)
 {
int i;
  if ((nrr!=other.nrr)||(length!=other.length))
	throw dmatrix_size_error(nrr, ncc, other.nrr, other.length);
for(i=0;i<length;i++)
  ary[i]-=other.ary[i];
 }

dmatrix operator+(const dmatrix &x1, const dmatrix &x2)
  {
int i;
  if ((x1.nrr!=x2.nrr)||(x1.length!=x2.length))
	throw dmatrix_size_error(x1.nrr, x1.ncc, x2.nrr, x2.length);
 dmatrix tempmat(x1.nrr,x1.ncc,x1.matoff);
  for(i=0;i<x1.length;i++) tempmat.ary[i]=x1.ary[i]+x2.ary[i];
return tempmat;
}

dmatrix operator-(const dmatrix &x1, const dmatrix &x2)
  {
int i;
  if ((x1.nrr!=x2.nrr)||(x1.length!=x2.length))
	throw dmatrix_size_error(x1.nrr, x1.ncc, x2.nrr, x2.length);
  dmatrix tempmat(x1.nrr,x1.ncc,x1.matoff);
  for(i=0;i<x1.length;i++) tempmat.ary[i]=x1.ary[i]-x2.ary[i];
return tempmat;
}

dmatrix operator*(const dmatrix& x1,const dmatrix& b)
	{
	int i,j,k;
	double xval,bval;
	if((x1.ncc!=b.nrr)||(x1.matoff!=b.matoff))
		throw dmatrix_size_error(x1.nrr, x1.ncc, b.nrr, b.length);
	dmatrix prod(x1.nrr,b.ncc,x1.matoff);
	for(i=x1.matoff;i<=x1.nrr;i++)for(j=x1.matoff;j<=b.ncc;j++)
		{prod(i,j)=0.0;
		for(k=x1.matoff;k<=x1.ncc;k++)
		   {
		   xval=x1.ary[k-x1.matoff+(x1.ncc-x1.matoff+1)*(i-x1.matoff)];
		   bval=b.ary[j-b.matoff+(b.ncc-b.matoff+1)*(k-b.matoff)];
		   prod(i,j) += xval*bval;
		   }
		}
	return prod;
	}

dmatrix operator*(const double& x, const dmatrix &zx)
  {
int i;
  dmatrix tempmat(zx.nrr,zx.ncc,zx.matoff);
  for(i=0;i<zx.length;i++) tempmat.ary[i]=x*zx.ary[i];
return tempmat;
  }

dmatrix operator/(const dmatrix &zx, const double& x)
  {
int i;
  dmatrix tempmat(zx.nrr,zx.ncc,zx.matoff);
  for(i=0;i<zx.length;i++) tempmat.ary[i]=zx.ary[i]/x;
return tempmat;
  }  


dmatrix inv(const dmatrix& x1)
{
 int i,icol,irow,j,k,jq,iq,*indxc,*indxr,*ipiv;
 double big,dum,pivinv,thold;

if(x1.nrr>x1.ncc)       //for pseudoinverses
  {
  dmatrix x1T(1,1,1);
  x1T=tr(x1);
  return (inv(x1T*x1))*x1T;
  }
if(x1.nrr<x1.ncc)   // underdetermined crude pseudoinverse
  {
  dmatrix x1T(1,1,1);
  x1T=tr(x1);
  return (x1T*(inv(x1*x1T)));
  }
 dmatrix a(x1.ncc,x1.ncc,x1.matoff);
 a=x1;
 a.matoff=1;
 a.ncc=a.nrr=x1.ncc-x1.matoff+1;
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
a.matoff=x1.matoff;  //set matrix offset and size to match originals
a.nrr=x1.nrr;
a.ncc=x1.ncc;
return a;
 }

dmatrix tr(const dmatrix& x1)
{
int i,j;
dmatrix temp(x1.ncc,x1.nrr,x1.matoff);
for(i=x1.matoff; i<=x1.nrr; i++)
   for(j=x1.matoff; j<=x1.ncc;j++)
temp.ary[i-temp.matoff+(temp.ncc-temp.matoff+1)*(j-temp.matoff)]=
   x1.ary[j-x1.matoff+(x1.ncc-x1.matoff+1)*(i-x1.matoff)];
return temp;
}


ostream& operator<<(ostream& os, dmatrix& x1)
  {
  int i,j;
  for(i=x1.matoff;i<=x1.nrr;i++)
  {
  for(j=x1.matoff;j<=x1.ncc;j++) os << x1(i,j) <<" ";
  os<<"\n";
  }
  return os;
  }

istream& operator>>(istream& is, dmatrix& x1)
  {
  int i,j;
  for(i=x1.matoff;i<=x1.nrr;i++)
  {
  for(j=x1.matoff;j<=x1.ncc;j++) is >> x1(i,j);
  }
  return is;
  }

