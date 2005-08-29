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

#include "source.h"
void vel_sfc(int iter)
{
    int isrc,sx,sy,sz,istf,inc;
    if( elmdls.src_func_type >= 4 ) inc = 1;
    else                           inc = 0;

    if( elmdls.src_type == 1 )
    {
	for(isrc=0,istf=0;isrc<elmdlv.isrc;isrc++)
	{
	    sx = elmdlv.srcx[isrc];
	    sy = elmdlv.srcy[isrc];
	    sz = elmdlv.srcz[isrc];
	    Uxt[sx][sy][sz]   += elmdlv.srctfunc[0][istf][iter];
	    Uyt[sx][sy][sz]   += elmdlv.srctfunc[1][istf][iter];
	    Uzt[sx][sy][sz]   += elmdlv.srctfunc[2][istf][iter];
	    istf += inc;
	}
    }
    else if( elmdlv.srcz[0] >= elmdls.start[2] && elmdlv.srcz[0] <= elmdls.end[2] )
    {
	sz =elmdlv.srcz[0];
	istf=0;
	for(sx=elmdls.start[0];sx<=elmdls.end[0];sx++)
	{
	    for(sy=elmdls.start[1];sy<=elmdls.end[1];sy++)
	    {
		Uzt[sx][sy][sz]   += elmdlv.srctfunc[2][istf][iter];
	    }
	}
    }
}

void stress_sfc(int iter)
{
    int isrc,sx,sy,sz,istf,inc;
    if( elmdls.src_func_type >= 4 ) inc = 1;
    else                           inc = 0;

    if( elmdls.src_type == 2 )
    {
	for(isrc=0,istf=0;isrc<elmdlv.isrc;isrc++)
	{
	    sx = elmdlv.srcx[isrc];
	    sy = elmdlv.srcy[isrc];
	    sz = elmdlv.srcz[isrc];
	    Txx[sx][sy][sz] += elmdlv.srctfunc[0][istf][iter];
	    Tyy[sx][sy][sz] += elmdlv.srctfunc[1][istf][iter];
	    Tzz[sx][sy][sz] += elmdlv.srctfunc[2][istf][iter];
	    istf += inc;
	}
    }
    else if( elmdlv.srcz[0] >= elmdls.start[2] && elmdlv.srcz[0] <= elmdls.end[2] )
    {
	sz =elmdlv.srcz[0];
	istf=0;
	for(sx=elmdls.start[0];sx<=elmdls.end[0];sx++)
	{
	    for(sy=elmdls.start[1];sy<=elmdls.end[1];sy++)
	    {
		Tzz[sx][sy][sz] += elmdlv.srctfunc[2][istf][iter];
	    }
	}
    }
}

void null_src_func(int iter)
{
	return;
}
void Make_source(void)
{
    int isrc;

    if( elmdlv.isrc == 0 && (elmdls.src_type == 1 || elmdls.src_type == 2) )
    {
	elmdlv.src_func = (*null_src_func);
	elmdlv.src = NULL;
	elmdlv.srctfunc = NULL;
	return;
    }

    switch(elmdls.src_func_type)
    {
	case 1:
		elmdlv.src = (float *)vector(0,elmdls.nt-1);
		elmdlv.srctfunc = (float ***)cubic(0,3, 0, 0, 0, elmdls.nt-1);
		ricker1(elmdlv.src);
		memcpy(elmdlv.srctfunc[0][0],elmdlv.src,elmdls.nt*sizeof(float));
		memcpy(elmdlv.srctfunc[1][0],elmdlv.src,elmdls.nt*sizeof(float));
		memcpy(elmdlv.srctfunc[2][0],elmdlv.src,elmdls.nt*sizeof(float));
		memcpy(elmdlv.srctfunc[3][0],elmdlv.src,elmdls.nt*sizeof(float));
		break;
	case 2:
		elmdlv.src = (float *)vector(0,elmdls.nt-1);
		elmdlv.srctfunc = (float ***)cubic(0,3, 0, 0, 0, elmdls.nt-1);
		gaussian1(elmdlv.src);
		memcpy(elmdlv.srctfunc[0][0],elmdlv.src,elmdls.nt*sizeof(float));
		memcpy(elmdlv.srctfunc[1][0],elmdlv.src,elmdls.nt*sizeof(float));
		memcpy(elmdlv.srctfunc[2][0],elmdlv.src,elmdls.nt*sizeof(float));
		memcpy(elmdlv.srctfunc[3][0],elmdlv.src,elmdls.nt*sizeof(float));
		break;
    }

    switch(elmdls.src_type)
    {
	case 1:
	case 3:
		elmdlv.src_func = (*vel_sfc);
		break;

	case 2: 
	case 4: 
		elmdlv.src_func = (*stress_sfc);
		break;
	default:
		shwarn("Source type is incorrect(%d).\n Using stress source\n",elmdls.src_type);
		elmdlv.src_func = (*stress_sfc);
    }

    switch(elmdls.src_cmpnt)
    {
	case SRC_CMPNT_X: 
		memset(elmdlv.srctfunc[1][0],(int)'\0',sizeof(float)*elmdls.nt);
		memset(elmdlv.srctfunc[2][0],(int)'\0',sizeof(float)*elmdls.nt);
		break;
	case SRC_CMPNT_Y: 
		memset(elmdlv.srctfunc[0][0],(int)'\0',sizeof(float)*elmdls.nt);
		memset(elmdlv.srctfunc[2][0],(int)'\0',sizeof(float)*elmdls.nt);
		break;
	case SRC_CMPNT_Z: 
		memset(elmdlv.srctfunc[0][0],(int)'\0',sizeof(float)*elmdls.nt);
		memset(elmdlv.srctfunc[1][0],(int)'\0',sizeof(float)*elmdls.nt);
		break;
	case SRC_CMPNT_XZ: 
		memset(elmdlv.srctfunc[2][0],(int)'\0',sizeof(float)*elmdls.nt);
		break;
	case SRC_CMPNT_TXZ: 
		memset(elmdlv.srctfunc[0][0],(int)'\0',sizeof(float)*elmdls.nt);
		memset(elmdlv.srctfunc[1][0],(int)'\0',sizeof(float)*elmdls.nt);
		break;
    }
}

void Free_source(void)
{
    if( elmdlv.isrc == 0 ) return;

    switch(elmdls.src_func_type)
    {
	case 1:
	case 2:
	case 3:
		freevec(&elmdlv.src,0,elmdls.nt-1);
		freecub(&elmdlv.srctfunc,0,3, 0,0, 0,elmdls.nt-1);
		break;
	case 4:
	case 90:
		freecub(&elmdlv.srctfunc,0,2,0,elmdlv.isrc-1,0,elmdls.nt-1);
		break;
    }

    elmdlv.srctfunc=NULL;
}

void ricker1(float *src)
{
    int i;
    float t,tmp;
    for(t=0.,i=0;i<=elmdls.nt;i++)
    {
	tmp = elmdls.fc*(t-elmdls.tc)*M_PI;
	tmp = tmp*tmp;
	src[i] = (1.-2.*tmp)*exp(-tmp)*elmdls.dt*1.e14;
	t += elmdls.dt;
    }
}

void gaussian1(float *src)
{
    int i;
    float t,tmp;
    for(t=0.,i=0;i<=elmdls.nt;i++)
    {
	tmp = elmdls.fc*(t-elmdls.tc)*M_PI;
	tmp = tmp*tmp;
	src[i] = (t-elmdls.tc)*exp(-tmp)*elmdls.dt*1.e14;
	t += elmdls.dt;
    }
}
