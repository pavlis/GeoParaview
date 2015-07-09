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

#include "model.h"

void MakePREM(void)
{
    int i,j,k,nx,ny,nz,layer,Nlayer;
    float rho,lambda,mu,vp,vs,dh;
    float R[]     = {2600., 2600., 2900., 2900., 3380.76, 3374.71, 3374.71, 3367.10, 3359.50, 3435.78, 3489.51, 3543.25, 3723.78, 3849.80 };
    float VP[]    = {5800., 5800., 6800., 6800., 8110.61, 8076.88, 8076.88, 8033.70, 7989.70, 8558.96, 8732.09, 8905.22, 9133.97, 9645.88 };
    float VS[]    = {3200., 3200., 3900., 3900., 4490.94, 4469.53, 4469.53, 4443.61, 4418.85, 4643.91, 4706.90, 4769.89, 4932.59, 5224.28 };
    float Depth[] = {6368., 6356., 6356., 6346.6,6346.6,  6291.0,  6291.0,  6221.0,  6151.0,  6151.0,  6061.0,  5971.0,  5971.0,  5871.0  };

    Nlayer = 14;
    nx = elmdls.nx;
    ny = elmdls.ny;
    nz = elmdls.nz;
    dh = elmdls.dh*0.001;

    for(i=1;i<Nlayer;i++) Depth[i] = Depth[0] - Depth[i];
    Depth[0] = 0.;

    for(layer=0;layer<Nlayer;layer++)
    {
	if( dh*(elmdls.start[2]-1) <= Depth[layer] ) break;
    	else layer++;
    }

    for(k=elmdls.start[2];k<=elmdls.end[2]+1;k++)
    {
	if( dh*(k-1) >= Depth[layer] && k != elmdls.start[2]) layer++;
	if( layer == Nlayer-1 ) layer--;

	if( Depth[layer] == Depth[layer-1] )
	{
	    rho =R[layer];
	    vp = VP[layer];
	    vs = VS[layer];
	}
        else if( Depth[Nlayer-1] < dh*(k-1) )
       	{
	    rho =R[Nlayer-1];
	    vp = VP[Nlayer-1];
	    vs = VS[Nlayer-1];
	}
	else
	{
	    rho = (R[layer-1]-R[layer])/(Depth[layer-1]-Depth[layer])*(dh*(k-1)-Depth[layer]) + R[layer];
	    vp = (VP[layer-1]-VP[layer])/(Depth[layer-1]-Depth[layer])*(dh*(k-1)-Depth[layer]) + VP[layer];
	    vs = (VS[layer-1]-VS[layer])/(Depth[layer-1]-Depth[layer])*(dh*(k-1)-Depth[layer]) + VS[layer];
	}
	mu     = rho*vs*vs;
	lambda = rho*vp*vp - 2.*mu;
	for(i=elmdls.start[0];i<=elmdls.end[0]+1;i++)
	    for(j=elmdls.start[1];j<=elmdls.end[1]+1;j++)
	    {
		   elmdlv.mpara[0][i][j][k] = rho;
		   elmdlv.mpara[1][i][j][k] = lambda;
		   elmdlv.mpara[2][i][j][k] = mu;
	    }
    }
}

void MakePREMsubduct(void)
{
    int i,j,k,nx,ny,nz,layer,Nlayer;
    float rho,lambda,mu,vp,vs,dh;
    float R[]     = {2600., 2600., 2900., 2900., 3380.76, 3374.71, 3374.71, 3367.10, 3359.50, 3435.78, 3489.51, 3543.25, 3723.78, 3849.80 };
    float VP[]    = {5800., 5800., 6800., 6800., 8110.61, 8076.88, 8076.88, 8033.70, 7989.70, 8558.96, 8732.09, 8905.22, 9133.97, 9645.88 };
    float VS[]    = {3200., 3200., 3900., 3900., 4490.94, 4469.53, 4469.53, 4443.61, 4418.85, 4643.91, 4706.90, 4769.89, 4932.59, 5224.28 };
    float Depth[] = {6368., 6356., 6356., 6346.6,6346.6,  6291.0,  6291.0,  6221.0,  6151.0,  6151.0,  6061.0,  5971.0,  5971.0,  5871.0  };
    float submu,sublambda,subratio,hlength,vlength,thick,angle;
    int ihlength,ivlength,ibegin,iend;

    Nlayer = 14;
    nx = elmdls.nx;
    ny = elmdls.ny;
    nz = elmdls.nz;
    dh = elmdls.dh*0.001;

    for(i=1;i<Nlayer;i++) Depth[i] = Depth[0] - Depth[i];
    Depth[0] = 0.;

    thick = Depth[3]; /* 21.4 km */
    vlength = 600.;   /* subducting slab goes down to 250 km depth */
    subratio = 1.1;   /* velocity ratio of subducting slab to surrounding area */
    angle = 45.;
    angle = angle/180.*M_PI;
    hlength = thick / sin(angle);
    ihlength = hlength / dh;
    ivlength = vlength / dh;
    if( ivlength > elmdls.nz ) ivlength=elmdls.nz;

    for(layer=0;layer<Nlayer;layer++)
    {
	if( dh*(elmdls.start[2]-1) <= Depth[layer] ) break;
    	else layer++;
    }

    for(k=elmdls.start[2];k<=elmdls.end[2]+1;k++)
    {
	if( dh*(k-1) >= Depth[layer] && k != elmdls.start[2]) layer++;
	if( layer == Nlayer-1 ) layer--;

	if( (Depth[layer] == Depth[layer-1]) || k == 1 )
	{
	    rho =R[layer];
	    vp = VP[layer];
	    vs = VS[layer];
	}
        else if( Depth[Nlayer-1] < dh*(k-1) )
       	{
	    rho =R[Nlayer-1];
	    vp = VP[Nlayer-1];
	    vs = VS[Nlayer-1];
	}
	else
	{
	    rho = (R[layer-1]-R[layer])/(Depth[layer-1]-Depth[layer])*(dh*(k-1)-Depth[layer]) + R[layer];
	    vp = (VP[layer-1]-VP[layer])/(Depth[layer-1]-Depth[layer])*(dh*(k-1)-Depth[layer]) + VP[layer];
	    vs = (VS[layer-1]-VS[layer])/(Depth[layer-1]-Depth[layer])*(dh*(k-1)-Depth[layer]) + VS[layer];
	}
	mu     = rho*vs*vs;
	lambda = rho*vp*vp - 2.*mu;

       if( !(lambda >=  3.e+10 || lambda <= 1.5e+11) ) printf("layer:%4d, k:%4d, r:%g, mu:%g, lambda:%g, vs:%g, vp:%g\n",layer,k,rho,mu,lambda,vs,vp);
	if( dh*(k-1) > thick && dh*(k-1) < thick+vlength )
	{
	    ibegin = elmdls.nx/2 - (dh*(k-1.) - thick) / tan( angle )/ dh;
	    iend = ibegin + ihlength;
	    if( ibegin < 1 ) ibegin = 1;
	    if( iend < 1 ) ibegin = 0;
	    submu = rho*vs*vs*subratio*subratio;
	    sublambda = rho*vp*vp*subratio*subratio - 2.*submu;
	}
	else ibegin = -1;

	for(i=elmdls.start[0];i<=elmdls.end[0]+1;i++)
	    for(j=elmdls.start[1];j<=elmdls.end[1]+1;j++)
	    {
		   elmdlv.mpara[0][i][j][k] = rho;
		   if( ibegin > 0 && j <= elmdls.ny/2 && i >= ibegin && i <= iend )
		   {
			elmdlv.mpara[1][i][j][k] = sublambda;
			elmdlv.mpara[2][i][j][k] = submu;
		   }
		   else
		   {
			elmdlv.mpara[1][i][j][k] = lambda;
			elmdlv.mpara[2][i][j][k] = mu;
		   }
	    }
    }
}

void Make_2layer(void)
{
    int i,j,k,nx,ny,nz;
    int layer;
    float rho,lambda,mu,vp,vs;

    nx = elmdls.nx;
    ny = elmdls.ny;
    nz = elmdls.nz;

    layer = nz/2;

/* gas in sandstone */
    rho  = 1.885e3;
/* 45 Hz */
    vs   = 0.99201e3;
    vp   = 1.4997e3;

    mu     = rho*vs*vs;
    lambda = rho*vp*vp - 2.*mu;

    for(i=elmdls.start[0];i<=elmdls.end[0]+1;i++)
	for(j=elmdls.start[1];j<=elmdls.end[1]+1;j++)
	    for(k=elmdls.start[2];k<=elmdls.end[2]+1;k++)
	    {
		if( k <= layer ) 
		{
		   elmdlv.mpara[0][i][j][k] = rho;
		   elmdlv.mpara[1][i][j][k] = lambda;
		   elmdlv.mpara[2][i][j][k] = mu;
                }
	    }

/* water in sandstone */
    rho     = 2.155e3;
/* 45 Hz */
    vs   = 0.92779e3;
    vp   = 2.2049e3;

    for(i=elmdls.start[0];i<=elmdls.end[0]+1;i++)
	for(j=elmdls.start[1];j<=elmdls.end[1]+1;j++)
	    for(k=elmdls.start[2];k<=elmdls.end[2]+1;k++)
	    {
		if( k > layer )
		{
		   elmdlv.mpara[0][i][j][k] = rho;
		   elmdlv.mpara[1][i][j][k] = lambda;
		   elmdlv.mpara[2][i][j][k] = mu;
		}
	    }
}

void Make_3D_Model(void)
{
    int i,j,k,nx,ny,nz,j1,j2,j3,j4;
    int z1,z2,z3,z4,xp,yp;
    float rho,lambda,mu,vp,vs;
    float y,xcoef,ycoef,coef;

    nx = elmdls.nx;
    ny = elmdls.ny;
    nz = elmdls.nz;

    /* nx = 500 */
    /* ny = 500 */
    /* nz = 250 */
    z1 = 70;
    z2 = 110;
    z3 = 140;
    z4 = 170;

    xp = 150;
    yp = ny/2;

    xcoef = 200.;
    ycoef = 200.;
    coef = 5.;

    xcoef = 1./xcoef;

    /* water in shale */
    rho     = 2.4685e3;
    vs   = 0.82987e3;
    vp   = 2.0521e3;

    mu     = rho*vs*vs;
    lambda = rho*vp*vp - 2.*mu;

    for(i=elmdls.start[0];i<=elmdls.end[0]+1;i++)
	for(j=elmdls.start[1];j<=elmdls.end[1]+1;j++)
	    for(k=elmdls.start[2];k<=elmdls.end[2]+1;k++)
	    {
		M_rho[i][j][k]    = rho;
		M_lambda[i][j][k] = lambda;
		M_mu[i][j][k]     = mu;
	    }

    /* gas in sandstone */
    rho  = 1.885e3;
    vs   = 0.992e3;
    vp   = 1.4997e3;

    mu     = rho*vs*vs;
    lambda = rho*vp*vp - 2.*mu;

    for(k=z1;k<z2;k++)
    {
	for(i=1;i<=nx+1;i++)
	{
	    y = k-z1-xcoef*(i-xp)*(i-xp);
	    if( i > xp ) y = k-z1-xcoef/coef*(i-xp)*(i-xp);
	    if( y < 0. ) continue;
	    j1 = -sqrt( ycoef*y ) + yp;
	    j2 =  sqrt( ycoef*y ) + yp;
	    if( j1 < 0 ) j1 = 0;
	    if( j2 < 0 ) j2 = -1;
	    if( j1 > ny ) j1 = ny;
	    if( j2 > ny ) j2 = ny;
	    for(j=j1;j<j2;j++)
	    {
		if( i >= elmdls.start[0] && i <= elmdls.end[0]+1 &&
		    j >= elmdls.start[1] && j <= elmdls.end[1]+1 &&
		    k >= elmdls.start[2] && k <= elmdls.end[2]+1 )
		{
		    M_rho[i][j][k]    = rho;
		    M_lambda[i][j][k] = lambda;
		    M_mu[i][j][k]     = mu;
		}
	    }
	}
    }

    /* oil in sandstone */
    rho     = 2.119e3;
    vs   = 0.93564e3;
    vp   = 2.0829e3;

    mu     = rho*vs*vs;
    lambda = rho*vp*vp - 2.*mu;

    for(k=z2;k<z3;k++)
    {
	for(i=1;i<=nx;i++)
	{
	    y = k-z1-xcoef*(i-xp)*(i-xp);
	    if( i > xp ) y = k-z1-xcoef/coef*(i-xp)*(i-xp);
	    if( y < 0. ) continue;
	    j1 = -sqrt( ycoef*y ) + yp;
	    j2 =  sqrt( ycoef*y ) + yp;
	    if( j1 < 0 ) j1 = 0;
	    if( j2 < 0 ) j2 = -1;
	    if( j1 > ny ) j1 = ny;
	    if( j2 > ny ) j2 = ny;
	    for(j=j1;j<=j2;j++)
	    {
		if( i >= elmdls.start[0] && i <= elmdls.end[0]+1 &&
		    j >= elmdls.start[1] && j <= elmdls.end[1]+1 &&
		    k >= elmdls.start[2] && k <= elmdls.end[2]+1 )
		{
		    M_rho[i][j][k]    = rho;
		    M_lambda[i][j][k] = lambda;
		    M_mu[i][j][k]     = mu;
		}
	    }
	}
    }

/* water in sandstone */
    rho     = 2.155e3;
    vs   = 0.92779e3;
    vp   = 2.2049e3;

    mu     = rho*vs*vs;
    lambda = rho*vp*vp - 2.*mu;

    for(k=z3;k<=nz;k++)
    {
	for(i=1;i<=nx;i++)
	{
	    y = k-z1-xcoef*(i-xp)*(i-xp);
	    if( i > xp ) y = k-z1-xcoef/coef*(i-xp)*(i-xp);
	    if( y < 0. ) continue;
	    j1 = -sqrt( ycoef*y ) + yp;
	    j4 =  sqrt( ycoef*y ) + yp;
	    if( k > z4 )
	    {
	         y = k-z4-xcoef*(i-xp)*(i-xp);
		 if( i > xp ) y = k-z4-xcoef/coef*(i-xp)*(i-xp);
		 if( y > 0. )
		 {
		    j2 = -sqrt( ycoef*y ) + yp;
		    j3 =  sqrt( ycoef*y ) + yp;
		 }
		 else
		 {
		    j2 = j4;
		    j3 = j4+1;
		 }
	    }
	    else
	    {
		    j2 = j4;
		    j3 = j4+1;
	    }
	    if( j1 < 0 ) j1 = 0;
	    if( j2 < 0 ) j2 = -1;
	    if( j1 > ny ) j1 = ny;
	    if( j2 > ny ) j2 = ny;
	    if( j3 > ny ) j3 = ny;
	    if( j4 > ny ) j4 = ny;
	    for(j=j1;j<=j2;j++)
	    {
		if( i >= elmdls.start[0] && i <= elmdls.end[0]+1 &&
		    j >= elmdls.start[1] && j <= elmdls.end[1]+1 &&
		    k >= elmdls.start[2] && k <= elmdls.end[2]+1 )
		{
		    M_rho[i][j][k]    = rho;
		    M_lambda[i][j][k] = lambda;
		    M_mu[i][j][k]     = mu;
		}
	    }
	    for(j=j3;j<=j4;j++)
	    {
		if( i >= elmdls.start[0] && i <= elmdls.end[0]+1 &&
		    j >= elmdls.start[1] && j <= elmdls.end[1]+1 &&
		    k >= elmdls.start[2] && k <= elmdls.end[2]+1 )
		{
		    M_rho[i][j][k]    = rho;
		    M_lambda[i][j][k] = lambda;
		    M_mu[i][j][k]     = mu;
		}
	    }
	}
    }
}
