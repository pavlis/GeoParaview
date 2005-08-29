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
#include "abc.h"
#include "surface.h"
#include "output.h"
#include "communicate_mpi.h"

#ifndef UPDATE
#define UPDATE

extern MDLs elmdls;
extern MDLv elmdlv;

void Update(void);

#endif
