/*###################################################################
	Parallel poroelastice wave propagation in 3D Version 1.0
	Last modified: April, 2005
  
	Author:	Dong-Hoon Sheen
		Ph.D. student at Seoul National University
	       	Research associate at Indiana University
	
	Reference: 
	Parallel implementation of a velocity-stress
	    staggered-grid finite-difference method 
	    for poroelastic wave propagation
	    Sheen, D.-H., Tuncay, K., Baag, C.-E. and Ortoleva, P. J.
###################################################################*/

#include "basic_mpi.h"

void sherror(char *fmt, ...)
{
    va_list args;
/*
    if( Myrank == 0 )
    {
    */
	if(EOF == fflush(stdout))
	    fprintf(stderr, "\nerror: fflush failed on stdout");
	va_start(args,fmt);
	vfprintf(stderr, fmt, args);
	vprintf(fmt, args);
	va_end(args);
	/*
    }
    */
    ErrorFlag = 1;
}

void shwarn(char *fmt, ...)
{
    va_list args;

    if( Myrank == 0 )
    {
	if(EOF == fflush(stdout))
	    fprintf(stderr, "\nerror: fflush failed on stdout");
	va_start(args,fmt);
	vfprintf(stderr, fmt, args);
	va_end(args);
    }
}

void ErrorCheck(char *errmsg)
{
    if( get_error_mpi(ErrorFlag) )
    {
	if( Myrank == 0 ) sherror(errmsg);
	MPI_Finalize();
	exit(EXIT_FAILURE);
    }
    return;
}

int Is_inside(int x,int z,int left,int top,int right,int bottom)
{
    if( x < left || x > right || z < top  || z > bottom ) return 0;
    return 1;
}

int Is_inside3D(int x,int y,int z,int left,int front,int top,int right,int back,int bottom)
{
    if( x < left || x > right || y < front || y > back || z < top  || z > bottom ) return 0;
    return 1;
}

FILE *echkfopen(char *fname,char *attr,void (*efnc)(char *fmt,...))
{
    FILE *fp;
    fp = fopen(fname,attr);
    if( (fp == NULL) && ((*efnc) != NULL) ) (*efnc)("%s File open failed\n",fname);
    return fp;
}

float *vector(int low,int high)
{
    float *array;

    if( ErrorFlag ) return NULL;
    if( high - low +1 < 1 )
    {
	sherror("Check Memory allocation range\n");
	return NULL;
    }
    array = (float *)malloc((high-low+1)*sizeof(float));
    if( array == NULL )
    {
	sherror("Memory allocation failure in vector\n");
	return NULL;
    }

    return (array-low);
}

void freevec(float **array,int low,int high)
{
    if( high - low +1 < 1 )
    {
	sherror("Check Memory allocation range\n");
	return;
    }
    free(*array+low);
    *array += low;
    *array = NULL;
}

float **matrix(int nrl, int nrh, int ncl, int nch)
{
    int i,nrow=nrh-nrl+1,ncol=nch-ncl+1;
    float **m;

    if( ErrorFlag ) return NULL;
    if( nrow < 1 || ncol < 1 )
    {
	sherror("Check Memory allocation range\n");
	return NULL;
    }

    m=(float **) malloc((size_t) ((nrow+1)*sizeof(float*)));
    if (!m)
    {
	sherror("allocation failure 1 in matrix\n");
	return NULL;
    }
    m += 1;
    m -= nrl;

    m[nrl]=(float *) malloc((size_t)((nrow*ncol+1)*sizeof(float)));
    if (!m[nrl]) 
    {
	sherror("allocation failure 2 in matrix\n");
	return NULL;
    }
    memset(m[nrl],(int)'\0',(size_t)((nrow*ncol+1)*sizeof(float)));

    m[nrl] += 1;
    m[nrl] -= ncl;
        
    for (i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

    return m;
}

void freemat(float ***m, int nrl, int nrh, int ncl, int nch)
{
    free((char *) ((*m)[nrl]+ncl-1));
    free((char *) ((*m)+nrl-1));
    *m = NULL;
}

float ***cubic(int nrl, int nrh, int ncl, int nch, int ndl, int ndh)
{
    long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
    float ***t;

    if( ErrorFlag ) return NULL;
    if( nrow < 1 || ncol < 1 || ndep < 1)
    {
	sherror("Check Memory allocation range\n");
	return NULL;
    }

    t=(float ***) malloc((size_t) ((nrow+1)*sizeof(float**)));
    if (!t)
    { 
	sherror("allocation failure 1 in cubic\n");
	return NULL;
    }
    t += 1;
    t -= nrl;

    t[nrl]=(float **) malloc((size_t)((nrow*ncol+1)*sizeof(float*)));
    if (!t[nrl]) 
    { 
	sherror("allocation failure 2 in cubic\n");
	return NULL;
    }
    t[nrl] += 1;
    t[nrl] -= ncl;
        
    t[nrl][ncl]=(float *) malloc((size_t)((nrow*ncol*ndep+1)*sizeof(float)));
    if (!t[nrl][ncl]) 
    {
	sherror("allocation failure 3 in cubic\n");
	return NULL;
    }
    memset(t[nrl][ncl],(int)'\0',(size_t)((nrow*ncol*ndep+1)*sizeof(float)));
    t[nrl][ncl] += 1;
    t[nrl][ncl] -= ndl;
        
    for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
    for(i=nrl+1;i<=nrh;i++)
    {
	t[i]=t[i-1]+ncol;
	t[i][ncl]=t[i-1][ncl]+ncol*ndep;
	for (j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
    }
        
    return t;
}


void freecub(float ****t, int nrl, int nrh, int ncl, int nch, int ndl, int ndh)
{
    free((char *) ((*t)[nrl][ncl]+ndl-1));
    free((char *) ((*t)[nrl]+ncl-1));
    free((char *) ((*t)+nrl-1));
    *t = NULL;
}

float ****quadri(int l1,int h1,int l2,int h2,int l3,int h3,int l4,int h4)
{
    long i,n1=h1-l1+1,n2=h2-l2+1,n3=h3-l3+1,n4=h4-l4+1;
    float ****qd;

    if( ErrorFlag ) return NULL;
    if( n1 < 1 || n2 < 1 || n3 < 1 || n4 < 1 )
    {
	sherror("Check Memory allocation range in quad\n");
	return NULL;
    }

    qd = (float ****)malloc( (size_t)(n1+1)*sizeof(float***));
    if( !qd )
    {
	sherror("Memory allocation failure in quad 1: %d - %d\n",l1,h1);
	return NULL;
    }
    qd += 1;
    qd -= l1;

    for(i=l1;i<=h1;i++)
    {
	qd[i] = (float ***)cubic(l2,h2,l3,h3,l4,h4);
	if( !qd[i] )
	{
	    sherror("Memory allocation failure in quad 2(%d): (%d - %d),(%d - %d),(%d - %d)\n",i,l2,h2,l3,h3,l4,h4);
	    return NULL;
	}
    }
    return qd;
}

void freequadri(float *****qd,int l1,int h1,int l2,int h2,int l3,int h3,int l4,int h4)
{
    int i;
    for(i=l1;i<=h1;i++)
	freecub( *qd+l1, l2,h2,l3,h3,l4,h4);

    free((char *) ((*qd)+l1-1));
    *qd = NULL;
}
