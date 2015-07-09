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

#include "abc.h"

void Alloc_PML_layer(void)
{
    int nx,ny,nz,gridz,flag,npml;

    nx   = elmdls.nx;
    ny   = elmdls.ny;
    nz   = elmdls.nz;
    npml = elmdls.npml;

    if( elmdls.freesurface )	gridz = 3;
    else		  	gridz = 1;

   flag = 0;
   /* Area 1           (1,nx; 1,npml; gridz,nz) ================================= */
   if( check_pml_region(1,nx, 1,npml, gridz,nz, elmdlv.cord_pml[0] ) )
   {
	elmdlv.pml1 = (float ****)quadri(0,23,elmdlv.cord_pml[0][0],	elmdlv.cord_pml[0][1], 
					     elmdlv.cord_pml[0][2],	elmdlv.cord_pml[0][3],
				   	     elmdlv.cord_pml[0][4],	elmdlv.cord_pml[0][5]);
	flag = 1;
   }
   /* Area 2           (1,nx; ny-npml+1,ny; gridz,nz) =========================== */
   if( check_pml_region(1,nx, ny-npml+1,ny, gridz,nz, elmdlv.cord_pml[1] ) )
   {
	elmdlv.pml2 = (float ****)quadri(0,23,elmdlv.cord_pml[1][0],	elmdlv.cord_pml[1][1], 
					     elmdlv.cord_pml[1][2],	elmdlv.cord_pml[1][3],
				   	     elmdlv.cord_pml[1][4],	elmdlv.cord_pml[1][5]);
	flag = 1;
   }
   /* Area 3           (1,npml; npml+1,ny-npml; gridz,nz) ======================= */
   if( check_pml_region(1,npml, npml+1,ny-npml, gridz,nz, elmdlv.cord_pml[2] ) )
   {
	elmdlv.pml3 = (float ****)quadri(0,23,elmdlv.cord_pml[2][0],	elmdlv.cord_pml[2][1], 
					     elmdlv.cord_pml[2][2],	elmdlv.cord_pml[2][3],
				   	     elmdlv.cord_pml[2][4],	elmdlv.cord_pml[2][5]);
	flag = 1;
   }
   /* Area 4           (nx-npml+1,nx; npml+1,ny-npml; gridz,nz) ================= */
   if( check_pml_region(nx-npml+1,nx, npml+1,ny-npml, gridz,nz, elmdlv.cord_pml[3] ) )
   {
	elmdlv.pml4 = (float ****)quadri(0,23,elmdlv.cord_pml[3][0],	elmdlv.cord_pml[3][1], 
					     elmdlv.cord_pml[3][2],	elmdlv.cord_pml[3][3],
				   	     elmdlv.cord_pml[3][4],	elmdlv.cord_pml[3][5]);
	flag = 1;
   }
   /* Area 5           (npml+1,nx-npml; npml+1,ny-npml; nz-npml+1,nz) ======= */
   if( check_pml_region(npml+1,nx-npml, npml+1,ny-npml, nz-npml+1,nz, elmdlv.cord_pml[4] ) )
   {
	elmdlv.pml5 = (float ****)quadri(0,23,elmdlv.cord_pml[4][0],	elmdlv.cord_pml[4][1], 
					     elmdlv.cord_pml[4][2],	elmdlv.cord_pml[4][3],
				   	     elmdlv.cord_pml[4][4],	elmdlv.cord_pml[4][5]);
	flag = 1;
   }
   /* No Free Sufrace Boundary ============================================== */
   /* Area 6           (npml+1,nx-npml; npml+1,ny-npml; 1,npml) ============= */
   if( !elmdls.freesurface )
   {
      if( check_pml_region(npml+1,nx-npml, npml+1,ny-npml, 1,npml, 	elmdlv.cord_pml[5] ) )
      {
	elmdlv.pml6 = (float ****)quadri(0,23,elmdlv.cord_pml[5][0],	elmdlv.cord_pml[5][1], 
					     elmdlv.cord_pml[5][2],	elmdlv.cord_pml[5][3],
				   	     elmdlv.cord_pml[5][4],	elmdlv.cord_pml[5][5]);
	flag = 1;
      }
   }

   if( flag )
   {
	elmdlv.apply_pml_stress   = (*Apply_PML_stress);
	elmdlv.apply_pml_velocity = (*Apply_PML_velocity);
   }
   else
   {
	elmdlv.apply_pml_stress   = NULL;
	elmdlv.apply_pml_velocity = NULL;
   }
}

void Make_PML_layer(void)
{
    int i,z;
    float a,w;
    if( elmdls.npml > 0 )
    {
	elmdlv.awx = (float **)matrix(0,1,1,elmdls.nx);
	elmdlv.awy = (float **)matrix(0,1,1,elmdls.ny);
	elmdlv.awz = (float **)matrix(0,1,1,elmdls.nz);

	for(i=elmdls.npml+1;i<=elmdls.nx-elmdls.npml;i++)
	{
	    elmdlv.awx[0][i] = 1.;
	    elmdlv.awx[1][i] = 0.;
	}
	for(i=elmdls.npml+1;i<=elmdls.ny-elmdls.npml;i++)
	{
	    elmdlv.awy[0][i] = 1.;
	    elmdlv.awy[1][i] = 0.;
	}
	if( elmdls.freesurface ) 	z = 1;
	else				z = elmdls.npml+1;
	for(i=z;i<=elmdls.nz-elmdls.npml;i++)
	{
	    elmdlv.awz[0][i] = 1.;
	    elmdlv.awz[1][i] = 0.;
	}

	for(i=1;i<=elmdls.npml;i++)
	{
	    a=1.+elmdls.am*pow((double)(elmdls.npml-i+1)/elmdls.npml,elmdls.npow);
	    w=2.*M_PI*elmdls.wm*elmdls.fc*pow((double)(elmdls.npml-i+1)/elmdls.npml,elmdls.npow+elmdls.apow);

	    elmdlv.awx[0][i] = a;
	    elmdlv.awx[1][i] = w;
	    elmdlv.awy[0][i] = a;
	    elmdlv.awy[1][i] = w;
	    if( !elmdls.freesurface )
	    {
		elmdlv.awz[0][i] = a;
		elmdlv.awz[1][i] = w;
	    }

	    elmdlv.awx[0][elmdls.nx-i+1] = a;
	    elmdlv.awx[1][elmdls.nx-i+1] = w;
	    elmdlv.awy[0][elmdls.ny-i+1] = a;
	    elmdlv.awy[1][elmdls.ny-i+1] = w;
	    elmdlv.awz[0][elmdls.nz-i+1] = a;
	    elmdlv.awz[1][elmdls.nz-i+1] = w;
	}
	Alloc_PML_layer();
    }
}

void Apply_PML_stress(void)
{
   if( elmdlv.cord_pml[0][0] != -1 )
	PML_stress(elmdlv.pml1,
		elmdlv.cord_pml[0][0],elmdlv.cord_pml[0][1],
		elmdlv.cord_pml[0][2],elmdlv.cord_pml[0][3],
		elmdlv.cord_pml[0][4],elmdlv.cord_pml[0][5]);
   if( elmdlv.cord_pml[1][0] != -1 )
	PML_stress(elmdlv.pml2,
		elmdlv.cord_pml[1][0],elmdlv.cord_pml[1][1],
		elmdlv.cord_pml[1][2],elmdlv.cord_pml[1][3],
		elmdlv.cord_pml[1][4],elmdlv.cord_pml[1][5]);
   if( elmdlv.cord_pml[2][0] != -1 )
	PML_stress(elmdlv.pml3,
		elmdlv.cord_pml[2][0],elmdlv.cord_pml[2][1],
		elmdlv.cord_pml[2][2],elmdlv.cord_pml[2][3],
		elmdlv.cord_pml[2][4],elmdlv.cord_pml[2][5]);
   if( elmdlv.cord_pml[3][0] != -1 )
	PML_stress(elmdlv.pml4,
		elmdlv.cord_pml[3][0],elmdlv.cord_pml[3][1],
		elmdlv.cord_pml[3][2],elmdlv.cord_pml[3][3],
		elmdlv.cord_pml[3][4],elmdlv.cord_pml[3][5]);
   if( elmdlv.cord_pml[4][0] != -1 )
	PML_stress(elmdlv.pml5,
		elmdlv.cord_pml[4][0],elmdlv.cord_pml[4][1],
		elmdlv.cord_pml[4][2],elmdlv.cord_pml[4][3],
		elmdlv.cord_pml[4][4],elmdlv.cord_pml[4][5]);

   if( elmdls.freesurface ) return;
   if( elmdlv.cord_pml[5][0] != -1 )
	PML_stress(elmdlv.pml6,
		elmdlv.cord_pml[5][0],elmdlv.cord_pml[5][1],
		elmdlv.cord_pml[5][2],elmdlv.cord_pml[5][3],
		elmdlv.cord_pml[5][4],elmdlv.cord_pml[5][5]);
}

void Apply_PML_velocity(void)
{
   if( elmdlv.cord_pml[0][0] != -1 )
	PML_velocity(elmdlv.pml1,
		elmdlv.cord_pml[0][0],elmdlv.cord_pml[0][1],
		elmdlv.cord_pml[0][2],elmdlv.cord_pml[0][3],
		elmdlv.cord_pml[0][4],elmdlv.cord_pml[0][5]);
   if( elmdlv.cord_pml[1][0] != -1 )
	PML_velocity(elmdlv.pml2,
		elmdlv.cord_pml[1][0],elmdlv.cord_pml[1][1],
		elmdlv.cord_pml[1][2],elmdlv.cord_pml[1][3],
		elmdlv.cord_pml[1][4],elmdlv.cord_pml[1][5]);
   if( elmdlv.cord_pml[2][0] != -1 )
	PML_velocity(elmdlv.pml3,
		elmdlv.cord_pml[2][0],elmdlv.cord_pml[2][1],
		elmdlv.cord_pml[2][2],elmdlv.cord_pml[2][3],
		elmdlv.cord_pml[2][4],elmdlv.cord_pml[2][5]);
   if( elmdlv.cord_pml[3][0] != -1 )
	PML_velocity(elmdlv.pml4,
		elmdlv.cord_pml[3][0],elmdlv.cord_pml[3][1],
		elmdlv.cord_pml[3][2],elmdlv.cord_pml[3][3],
		elmdlv.cord_pml[3][4],elmdlv.cord_pml[3][5]);
   if( elmdlv.cord_pml[4][0] != -1 )
	PML_velocity(elmdlv.pml5,
		elmdlv.cord_pml[4][0],elmdlv.cord_pml[4][1],
		elmdlv.cord_pml[4][2],elmdlv.cord_pml[4][3],
		elmdlv.cord_pml[4][4],elmdlv.cord_pml[4][5]);

   if( elmdls.freesurface ) return;
   if( elmdlv.cord_pml[5][0] != -1 )
	PML_velocity(elmdlv.pml6,
		elmdlv.cord_pml[5][0],elmdlv.cord_pml[5][1],
		elmdlv.cord_pml[5][2],elmdlv.cord_pml[5][3],
		elmdlv.cord_pml[5][4],elmdlv.cord_pml[5][5]);
}

void PML_stress(float ****pmlp,int sx,int ex, int sy,int ey, int sz,int ez)
{
    int i,j,k;
    static float c0 = 9./8., c1=1./24.;
    float dt=elmdls.dt,dh=1./elmdls.dh;
    float ax,ay,az,wx,wy,wz;
    float fx1,fx2,fy1,fy2,fz1,fz2;
    float mu,lambda;
    float Dx,Dy,Dz;

    for(i=sx;i<=ex;i++)
    {
	ax = elmdlv.awx[0][i];
	wx = elmdlv.awx[1][i];
       	fx1 = (2.*ax-wx*dt)/(2.*ax+wx*dt);
       	fx2 = 2.*dh*dt/(2.*ax+wx*dt);
	for(j=sy;j<=ey;j++)
	{
	    ay = elmdlv.awy[0][j];
	    wy = elmdlv.awy[1][j];
       	    fy1 = (2.*ay-wy*dt)/(2.*ay+wy*dt);
       	    fy2 = 2.*dh*dt/(2.*ay+wy*dt);
	    for(k=sz;k<=ez;k++)
	    {
	        az = elmdlv.awz[0][k];
	        wz = elmdlv.awz[1][k];
       	        fz1 = (2.*az-wz*dt)/(2.*az+wz*dt);
       	        fz2 = 2.*dh*dt/(2.*az+wz*dt);

		lambda = M_lambda[i][j][k];
		mu = lambda + 2.*M_mu[i][j][k];

		Dx = c0*(Uxt[i][j][k]-Uxt[i-1][j][k])-c1*(Uxt[i+1][j][k]-Uxt[i-2][j][k]);
		xTxx[i][j][k]=fx1*xTxx[i][j][k]+fx2*mu*Dx;
		xTyy[i][j][k]=fx1*xTyy[i][j][k]+fx2*lambda*Dx;
		xTzz[i][j][k]=fx1*xTzz[i][j][k]+fx2*lambda*Dx;

		Dy = c0*(Uyt[i][j][k]-Uyt[i][j-1][k])-c1*(Uyt[i][j+1][k]-Uyt[i][j-2][k]);
		yTxx[i][j][k]=fy1*yTxx[i][j][k]+fy2*lambda*Dy;
		yTyy[i][j][k]=fy1*yTyy[i][j][k]+fy2*mu*Dy;
		yTzz[i][j][k]=fy1*yTzz[i][j][k]+fy2*lambda*Dy;
				
		Dz = c0*(Uzt[i][j][k]-Uzt[i][j][k-1])-c1*(Uzt[i][j][k+1]-Uzt[i][j][k-2]);
		zTxx[i][j][k]=fz1*zTxx[i][j][k]+fz2*lambda*Dz;
		zTyy[i][j][k]=fz1*zTyy[i][j][k]+fz2*lambda*Dz;
		zTzz[i][j][k]=fz1*zTzz[i][j][k]+fz2*mu*Dz;
				
		mu = 4./(1./M_mu[i][j][k]+1./M_mu[i][j+1][k]+1./M_mu[i+1][j][k]+1./M_mu[i+1][j+1][k]);
		Dx = (c0*(Uyt[i+1][j][k]-Uyt[i][j][k])-c1*(Uyt[i+2][j][k]-Uyt[i-1][j][k]));
		Dy = (c0*(Uxt[i][j+1][k]-Uxt[i][j][k])-c1*(Uxt[i][j+2][k]-Uxt[i][j-1][k]));
		xTxy[i][j][k]=fx1*xTxy[i][j][k]+fx2*mu*Dx;
		yTxy[i][j][k]=fy1*yTxy[i][j][k]+fy2*mu*Dy;

		mu = 4./(1./M_mu[i][j][k]+1./M_mu[i+1][j][k]+1./M_mu[i][j][k+1]+1./M_mu[i+1][j][k+1]);
		Dx = (c0*(Uzt[i+1][j][k]-Uzt[i][j][k])-c1*(Uzt[i+2][j][k]-Uzt[i-1][j][k]));
		Dz = (c0*(Uxt[i][j][k+1]-Uxt[i][j][k])-c1*(Uxt[i][j][k+2]-Uxt[i][j][k-1]));
		xTxz[i][j][k]=fx1*xTxz[i][j][k]+fx2*mu*Dx;
		zTxz[i][j][k]=fz1*zTxz[i][j][k]+fz2*mu*Dz;

		mu = 4./(1./M_mu[i][j][k]+1./M_mu[i][j+1][k]+1./M_mu[i][j][k+1]+1./M_mu[i][j+1][k+1]);
		Dy = (c0*(Uzt[i][j+1][k]-Uzt[i][j][k])-c1*(Uzt[i][j+2][k]-Uzt[i][j-1][k]));
		Dz = (c0*(Uyt[i][j][k+1]-Uyt[i][j][k])-c1*(Uyt[i][j][k+2]-Uyt[i][j][k-1]));
		yTyz[i][j][k]=fy1*yTyz[i][j][k]+fy2*mu*Dy;
		zTyz[i][j][k]=fz1*zTyz[i][j][k]+fz2*mu*Dz;

		Txx[i][j][k] = xTxx[i][j][k]+yTxx[i][j][k]+zTxx[i][j][k];
		Tyy[i][j][k] = xTyy[i][j][k]+yTyy[i][j][k]+zTyy[i][j][k];
		Tzz[i][j][k] = xTzz[i][j][k]+yTzz[i][j][k]+zTzz[i][j][k];
		Txy[i][j][k] = xTxy[i][j][k]+yTxy[i][j][k];
		Txz[i][j][k] = xTxz[i][j][k]              +zTxz[i][j][k];
		Tyz[i][j][k] =               yTyz[i][j][k]+zTyz[i][j][k];
	    }
	}
    }
}


void PML_velocity(float ****pmlp,int sx,int ex, int sy,int ey, int sz,int ez)
{
    int i,j,k;
    static float c0 = 9./8., c1=1./24.;
    float dt=elmdls.dt,dh=1./elmdls.dh;
    float ax,ay,az,wx,wy,wz;
    float fx1,fx2,fy1,fy2,fz1,fz2;
    float buoy_x,buoy_y,buoy_z;
    float Dx,Dy,Dz;

    for(i=sx;i<=ex;i++)
    {
	ax = elmdlv.awx[0][i];
	wx = elmdlv.awx[1][i];
       	fx1 = (2.*ax-wx*dt)/(2.*ax+wx*dt);
       	fx2 = 2.*dh*dt/(2.*ax+wx*dt);
	for(j=sy;j<=ey;j++)
	{
	    ay = elmdlv.awy[0][j];
	    wy = elmdlv.awy[1][j];
       	    fy1 = (2.*ay-wy*dt)/(2.*ay+wy*dt);
       	    fy2 = 2.*dh*dt/(2.*ay+wy*dt);
	    for(k=sz;k<=ez;k++)
	    {
	        az = elmdlv.awz[0][k];
	        wz = elmdlv.awz[1][k];
       	        fz1 = (2.*az-wz*dt)/(2.*az+wz*dt);
       	        fz2 = 2.*dh*dt/(2.*az+wz*dt);

		buoy_x = (1./M_rho[i][j][k] +1./M_rho[i+1][j][k])*0.5;
		buoy_y = (1./M_rho[i][j][k] +1./M_rho[i][j+1][k])*0.5;
		buoy_z = (1./M_rho[i][j][k] +1./M_rho[i][j][k+1])*0.5;

		Dx = (c0*(Txx[i+1][j][k]-Txx[i][j][k])  -c1*(Txx[i+2][j][k]-Txx[i-1][j][k]));
		Dy = (c0*(Txy[i][j][k]  -Txy[i-1][j][k])-c1*(Txy[i+1][j][k]-Txy[i-2][j][k]));
		Dz = (c0*(Txz[i][j][k]  -Txz[i-1][j][k])-c1*(Txz[i+1][j][k]-Txz[i-2][j][k]));

		xUxt[i][j][k]=fx1*xUxt[i][j][k]+fx2*buoy_x*Dx;
		xUyt[i][j][k]=fx1*xUyt[i][j][k]+fx2*buoy_y*Dy;
		xUzt[i][j][k]=fx1*xUzt[i][j][k]+fx2*buoy_z*Dz;

		Dx = (c0*(Txy[i][j][k]  -Txy[i][j-1][k])-c1*(Txy[i][j+1][k]-Txy[i][j-2][k]));
		Dy = (c0*(Tyy[i][j+1][k]-Tyy[i][j][k])  -c1*(Tyy[i][j+2][k]-Tyy[i][j-1][k]));
		Dz = (c0*(Tyz[i][j][k]  -Tyz[i][j-1][k])-c1*(Tyz[i][j+1][k]-Tyz[i][j-2][k]));

		yUxt[i][j][k]=fy1*yUxt[i][j][k]+fy2*buoy_x*Dx;
		yUyt[i][j][k]=fy1*yUyt[i][j][k]+fy2*buoy_y*Dy;
		yUzt[i][j][k]=fy1*yUzt[i][j][k]+fy2*buoy_z*Dz;

		Dx = (c0*(Txz[i][j][k]  -Txz[i][j][k-1])-c1*(Txz[i][j][k+1]-Txz[i][j][k-2]));
		Dy = (c0*(Tyz[i][j][k]  -Tyz[i][j][k-1])-c1*(Tyz[i][j][k+1]-Tyz[i][j][k-2]));
		Dz = (c0*(Tzz[i][j][k+1]-Tzz[i][j][k])  -c1*(Tzz[i][j][k+2]-Tzz[i][j][k-1]));

		zUxt[i][j][k]=fz1*zUxt[i][j][k]+fz2*buoy_x*Dx;
		zUyt[i][j][k]=fz1*zUyt[i][j][k]+fz2*buoy_y*Dy;
		zUzt[i][j][k]=fz1*zUzt[i][j][k]+fz2*buoy_z*Dz;

                Uxt[i][j][k] = xUxt[i][j][k]+yUxt[i][j][k]+zUxt[i][j][k];
		Uyt[i][j][k] = xUyt[i][j][k]+yUyt[i][j][k]+zUyt[i][j][k];
		Uzt[i][j][k] = xUzt[i][j][k]+yUzt[i][j][k]+zUzt[i][j][k];
	    }
	}
    }
}

/* Clatyon & Engquist (1977): apply A1 to velocity */
/*
void Apply_ABC(void)
{
	register int i,j;
	float lambda,mu,buoy;
	static int nx,nz;
	static float dt,dh;

	nx = elmdl.nx;
	nz = elmdl.nz;
	dt = elmdl.dt;
	dh = elmdl.dh;

	for(i=1,j=nz+1;i<=nx;i++)
	{
	  lambda = M_lambda[i][j];
	  mu      = M_mu[i][j];
	  buoy = 1./M_rho[i][j];
	  Uxt[i][j]-=sqrt(mu*buoy)*(Uxt[i][j]-Uxt[i][j-1])*dt/dh;
	  Uzt[i][j]-=sqrt((lambda+2.*mu)*buoy)*(Uzt[i][j]-Uzt[i][j-1])*dt/dh;
	}

	for(i=0,j=1;j<=nz;j++)
	{
	  lambda = M_lambda[1][j];
	  mu      = M_mu[1][j];
	  buoy = 1./M_rho[1][j];
	  Uxt[i][j]+=sqrt((lambda+2.*mu)*buoy)*(Uxt[i+1][j]-Uxt[i][j])*dt/dh;
	  Uzt[i][j]+=sqrt(mu*buoy)*(Uzt[i+1][j]-Uzt[i][j])*dt/dh;
	}

	for(i=nx+1,j=1;j<=nz;j++)
	{
	  lambda = M_lambda[i][j];
	  mu      = M_mu[i][j];
	  buoy = 1./M_rho[i][j];
	  Uxt[i][j]-=sqrt((lambda+2.*mu)*buoy)*(Uxt[i][j]-Uxt[i-1][j])*dt/dh;
	  Uzt[i][j]-=sqrt(mu*buoy)*(Uzt[i][j]-Uzt[i-1][j])*dt/dh;
	}

	i=0;
	lambda = M_lambda[1][j];
	mu      = M_mu[1][j];
	buoy = 1./M_rho[1][j];
	lambda = (sqrt(mu*buoy)+sqrt((lambda+2.*mu)*buoy))*0.5;
	Uxt[i][j]-=lambda*(Uxt[i][j]-Uxt[i][j-1])*dt/dh;
	Uzt[i][j]-=lambda*(Uzt[i][j]-Uzt[i][j-1])*dt/dh;

	i=nx+1;
	lambda = M_lambda[i][j];
	mu      = M_mu[i][j];
	buoy = 1./M_rho[i][j];
	lambda = (sqrt(mu*buoy)+sqrt((lambda+2.*mu)*buoy))*0.5;
	Uxt[i][j]-=lambda*(Uxt[i][j]-Uxt[i][j-1])*dt/dh;
	Uzt[i][j]-=lambda*(Uzt[i][j]-Uzt[i][j-1])*dt/dh;
}
*/

