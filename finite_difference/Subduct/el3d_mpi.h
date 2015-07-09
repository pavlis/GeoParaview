/*###################################################################
	Parallel poroelastice wave propagation in 3D Version 1.0
	April, 2005
  
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
#include <math.h>
#include <mpi.h>

#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) < (b) ? (a) : (b))

#include "mdlpara.h"
#include "read_setup.h"
#include "initialize.h"
#include "update.h"

#ifndef PMWAVE
#define PMWAVE

MDLs elmdls;
MDLv elmdlv;

char Gfname[BUFSIZ];
FILE *Gfp;

#endif

#ifndef BASICMPI
#define BASICMPI
int Myrank;
int ErrorFlag;
#else
extern int Myrank;
extern int ErrorFlag;
#endif

void find_opt(int *m1,int *m2,int *m3,int number);
