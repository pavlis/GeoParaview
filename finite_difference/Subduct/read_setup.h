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

#include "mdlpara.h"
#include "basic_mpi.h"
#include "getpar.h"
#include "source.h"
#include "abc.h"

#ifndef READSETUP
#define READSETUP

extern char Gfname[BUFSIZ];
extern FILE *Gfp;

extern MDLs elmdls;
extern MDLv elmdlv;

int Read_setup(char *fname);
void validate_input(void);
#endif
