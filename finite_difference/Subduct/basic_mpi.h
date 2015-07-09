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

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <mpi.h>

#ifndef BASICMPI
#define BASICMPI

#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) < (b) ? (a) : (b))

int Myrank;
int ErrorFlag;
#else
extern int Myrank;
extern int ErrorFlag;
#endif

void sherror(char *fmt, ...);
void shwarn(char *fmt, ...);
void ErrorCheck(char *errmsg);
int Is_inside(int x,int z,int left,int top,int right,int bottom);
int Is_inside3D(int x,int y,int z,int left,int front,int top,int right,int back,int bottom);
FILE *echkfopen(char *fname,char *attr,void (*efnc)(char *fmt,...));
float *vector(int low,int high);
void freevec(float **array,int low,int high);
float **matrix(int lr,int hr,int lc,int hc);
void freemat(float ***mat,int lr,int hr,int lc,int hc);
float ***cubic(int l1,int h1,int l2,int h2,int l3,int h3);
void freecub(float ****cub,int l1,int h1,int l2,int h2,int l3,int h3);
float ****quadri(int l1,int h1,int l2,int h2,int l3,int h3,int l4,int h4);
void freequadri(float *****quadri,int l1,int h1,int l2,int h2,int l3,int h3,int l4,int h4);
