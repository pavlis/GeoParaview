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
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <mpi.h>

#include "mdlpara.h"
#include "basic_mpi.h"

#ifndef OUTPUT
#define OUTPUT

extern MDLs elmdls;
extern MDLv elmdlv;

extern char Gfname[BUFSIZ];
extern FILE *Gfp;

#endif

void Print_data(int iter);
void Write_mdl_parameter(void);
void write_yz(float ****outdata,int cmp,int gnum);
void write_xy(float ****outdata,int cmp,int gnum);
void write_xz(float ****outdata,int cmp,int gnum);
void out_func(int gnum,char *plane,int iter,void (*write_func)(float ****outdata,int cmp,int gnum));
void Save_seismogram(int iter);
void Save_seismogram_step(int iter);

