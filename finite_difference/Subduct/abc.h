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

#ifndef PML_ABC
#define PML_ABC

/* components of pml region */
#define xUxt  pmlp[0]
#define yUxt  pmlp[1]
#define zUxt  pmlp[2]

#define xUyt  pmlp[3]
#define yUyt  pmlp[4]
#define zUyt  pmlp[5]

#define xUzt  pmlp[6]
#define yUzt  pmlp[7]
#define zUzt  pmlp[8]

#define xTxx  pmlp[9]
#define yTxx  pmlp[10]
#define zTxx  pmlp[11]

#define xTyy  pmlp[12]
#define yTyy  pmlp[13]
#define zTyy  pmlp[14]

#define xTzz  pmlp[15]
#define yTzz  pmlp[16]
#define zTzz  pmlp[17]

#define xTxy  pmlp[18]
#define yTxy  pmlp[19]

#define xTxz  pmlp[20]
#define zTxz  pmlp[21]

#define yTyz  pmlp[22]
#define zTyz  pmlp[23]

extern MDLs elmdls;
extern MDLv elmdlv;

void Alloc_PML_layer(void);
void Make_PML_layer(void);

void Apply_PML_stress(void);
void Apply_PML_velocity(void);
void PML_stress(float ****pmlp,int sx,int ex, int sy,int ey, int sz,int ez);
void PML_velocity(float ****pmlp,int sx,int ex, int sy,int ey, int sz,int ez);
#endif
