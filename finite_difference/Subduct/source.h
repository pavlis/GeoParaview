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

#include "mdlpara.h"
#include "basic_mpi.h"

#ifndef SOURCE
#define SOURCE

#define SRC_CMPNT_ALL	0
#define SRC_CMPNT_X	1
#define SRC_CMPNT_Y	2
#define SRC_CMPNT_Z	3
#define SRC_CMPNT_XZ	4
#define SRC_CMPNT_TXZ	5

extern MDLs elmdls;
extern MDLv elmdlv;

void Make_source(void);
void Free_source(void);
void ricker1(float *src);
void gaussian1(float *src);

#endif
