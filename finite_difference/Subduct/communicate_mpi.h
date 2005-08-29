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
#include <string.h>
#include <mpi.h>

#include "mdlpara.h"
#include "basic_mpi.h"

#ifndef Communicate_MPI
#define Communicate_MPI

extern MDLs elmdls;
extern MDLv elmdlv;

extern char Gfname[BUFSIZ];
extern FILE *Gfp;

int get_error_mpi(int evalue);

void mem2buf_yz(float ****org,float *buf,int ix,int np);
void buf2mem_yz(float ****org,float *buf,int ix,int np);
void mem2buf_xz(float ****org,float *buf,int jy,int np);
void buf2mem_xz(float ****org,float *buf,int jy,int np);
void mem2buf_xy(float ****org,float *buf,int kz,int np);
void buf2mem_xy(float ****org,float *buf,int kz,int np);
void comm_basic(float ****value,float *buf,int np,int direction,int area,
	void (*mem2buf)(float ****org,float *buf,int ix,int np),
	void (*buf2mem)(float ****org,float *buf,int ix,int np));
void communicate(float ****value, float *buf, int np);
void comm_distribute_setup(int mrank);
void comm_write_info(char *tbuf, char *fname, int type, int csize);
void comm_write_seismogram(void);
#endif



