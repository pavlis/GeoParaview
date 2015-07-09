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
#include <mpi.h>

#include "basic_mpi.h"
#include "mdlpara.h"
#include "source.h"
#include "surface.h"
#include "model.h"
#include "output.h"
#include "abc.h"
#include "communicate_mpi.h"

#ifndef InitializeModel
#define InitializeModel

extern MDLs elmdls;
extern MDLv elmdlv;

extern int Myrank;
extern int ErrorFlag;

extern char Gfname[BUFSIZ];
extern FILE *Gfp;

void Initialize(int n1,int n2,int n3);
void Print_Configuration(void);
void Alloc_Mem_Func(void);
void Set_receiver(void);
void Set_source(void);
void Init_model_parameter(void);
int  check_pml_region(int sx,int ex,int sy,int ey,int sz,int ez,int *cord);
void decompose_domain(int cord,int ngrid);
void decompose_domain2(int cord,int ngrid,int *sizes);
void set_priority(int nx,int ny,int nz, int *p1,int *p2,int *p3);
void set_rank_connection(int n1,int n2,int n3,int i,int j,int k);
void set_receiver(void);
void set_boundary(int nx,int ny,int nz,int npml);
void set_snapshot_location(void);

#endif

