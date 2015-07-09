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

#ifndef MODEL
#define MODEL

extern MDLs elmdls;
extern MDLv elmdlv;

extern char Gfname[BUFSIZ];
extern FILE *Gfp;

void Make_2layer(void);
void Make_3D_Model(void);
void MakePREM(void);

#endif
