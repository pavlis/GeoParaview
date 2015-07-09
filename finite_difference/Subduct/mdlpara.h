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

/*===================================================================
 Model Parameter, Source, Absorbing Boundary Condition define/declare 
===================================================================*/

#ifndef MDLPARAMETER
#define MDLPARAMETER

/* time derivative of displacement */
#define Uxt  elmdlv.vel[0]  
#define Uyt  elmdlv.vel[1]  
#define Uzt  elmdlv.vel[2]

/* time derivative of stress and pore fluid pressure */
#define Txx  elmdlv.stress[0]
#define Tyy  elmdlv.stress[1]
#define Tzz  elmdlv.stress[2]
#define Txy  elmdlv.stress[3]
#define Txz  elmdlv.stress[4]
#define Tyz  elmdlv.stress[5]

/* media parameter */
#define M_rho     elmdlv.mpara[0]  
#define M_lambda  elmdlv.mpara[1]
#define M_mu      elmdlv.mpara[2]

/* Snapshot component */
#define SN_X 1
#define SN_Y 2
#define SN_Z 4
/* Receiver component */
#define RC_X 1
#define RC_Y 2
#define RC_Z 4
#define RC_D 8
#define RC_CX 16
#define RC_CY 32
#define RC_CZ 64

struct mdl_static_struct{
	int rank[4];
	int nProcs[4];
	int leftProc[3];
	int rightProc[3];
	int start[3];
	int end[3];

	int nx;
	int ny;
	int nz;
	int px;
	int py;
	int pz;
	int nt;
	int npml;
	int model;
	int interimflag;
	int sample;
	int outstep;
	int rtype;

	int res_type;
	int res_file_type;
	int src_type;
	int src_func_type;
       	int src_cmpnt;
	int shstart;
	int shstep;
	int shcmpnt;
	int shnloc;
	int nshot;

	int ishot;
	int nOutPara,OutParaInfo[8];

        int nres_cmpnt;
	int nrcmpnt;
	int rcmpnt;
	int freesurface;

	float rho;
	float lambda;
	float mu;
	float Vp;
	float Vs;

	float dh;
	float dt;
	float fc;
	float tc;
        float am;
        float wm;
        float npow;
        float apow;
};
typedef struct mdl_static_struct MDLs;

struct mdl_variable_struct{
	int *shxsec;
	int *shgnum;
	int *rcmpnt_list;
	int **cord_pml;

	int **sx;
	int **sy; 
	int **sz;
	int **sloc;
       	int isrc;
	int *srcx;
	int *srcy;
       	int *srcz;
       	int *src_loc;

	int **gx;
	int **gy;
	int **gz;
       	int **rloc;
       	int irec;
       	int *recx;
       	int *recy;
       	int *recz;
       	int *rec_loc;

	int *tnrec;
	int *nsrc;
	int *nrec;

	char *mfile;
	char *wd;
	char *oinfo;
	char *res_name;

	float *src;
	float ***srctfunc;
	float ***seismogrm;
	float ****vel;
	float ****stress;
	float ****mpara;
	float ****pml1;
	float ****pml2;
	float ****pml3;
	float ****pml4;
	float ****pml5;
	float ****pml6;
	float **awx,**awy,**awz;

	void (*src_func)(int iteration);
    	void (*save_seismogram)(int iteration); 
	void (*apply_pml_stress)(void);
	void (*apply_pml_velocity)(void);
	void (*apply_surface_stress)(void);
	void (*apply_surface_velocity)(void);
};
typedef struct mdl_variable_struct MDLv;
#endif

