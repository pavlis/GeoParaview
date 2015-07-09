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

#include "surface.h"

void Apply_Surface_stress(void)
{
   int i,j,k;
   float tlu,tl,Dx,Dy,Dz,Dp;
   float lambda,mu;
   float dt,dh,c0,c1,deno;

   dt = elmdls.dt;
   dh = 1./elmdls.dh;
   c0 = 9./8.;
   c1 = 1./24.;

   for(i=elmdls.start[0];i<=elmdls.end[0];i++)
   {
	for(j=elmdls.start[1];j<=elmdls.end[1];j++)
	{
              k=1;
	      lambda = M_lambda[i][j][k];
	      mu      = M_mu[i][j][k];

	      tlu 	= dt*(lambda+2.*mu);
	      tl	= dt*lambda;
	      deno    = lambda+2.*mu;
	      	      
	      Tzz[i][j][1] = 0.;
	      Dx = (c0*(Uxt[i][j][1]-Uxt[i-1][j][1])-c1*(Uxt[i+1][j][1]-Uxt[i-2][j][1]))*dh;
	      Dy = (c0*(Uyt[i][j][1]-Uyt[i][j-1][1])-c1*(Uyt[i][j+1][1]-Uyt[i][j-2][1]))*dh;
	      Dz = -lambda*(Dx+Dy)/deno;
	      Dp = Dz;
	      Txx[i][j][1] += tlu*Dx + tl*(Dy+Dz);
	      Tyy[i][j][1] += tlu*Dy + tl*(Dx+Dz);
	      	      
	      mu = 4./(1./M_mu[i][j][1]+1./M_mu[i+1][j][1]+1./M_mu[i][j+1][1]+1./M_mu[i+1][j+1][1]);
	      Dx = (c0*(Uyt[i+1][j][1]-Uyt[i][j][1])-c1*(Uyt[i+2][j][1]-Uyt[i-1][j][1]))*dh;
	      Dy = (c0*(Uxt[i][j+1][1]-Uxt[i][j][1])-c1*(Uxt[i][j+2][1]-Uxt[i][j-1][1]))*dh;
	      Txy[i][j][1]+=dt*mu*(Dx+Dy);
	      	      
	      mu = 4./(1./M_mu[i][j][1]+1./M_mu[i][j][2]+1./M_mu[i+1][j][1]+1./M_mu[i+1][j][2]);
	      Dx = (c0*(Uzt[i+1][j][1]-Uzt[i][j][1])-c1*(Uzt[i+2][j][1]-Uzt[i-1][j][1]))*dh;
	      Dz = (-22.*Uxt[i][j][1]+17.*Uxt[i][j][2]+9.*Uxt[i][j][3]-5.*Uxt[i][j][4]+Uxt[i][j][5])/24.*dh;
	      Txz[i][j][1]+=dt*mu*(Dx+Dz);
	      	      
	      mu = 4./(1./M_mu[i][j][1]+1./M_mu[i][j][2]+1./M_mu[i][j+1][1]+1./M_mu[i][j+1][2]);
	      Dy = (c0*(Uzt[i][j+1][1]-Uzt[i][j][1])-c1*(Uzt[i][j+2][1]-Uzt[i][j-1][1]))*dh;
	      Dz = (-22.*Uyt[i][j][1]+17.*Uyt[i][j][2]+9.*Uyt[i][j][3]-5.*Uyt[i][j][4]+Uyt[i][j][5])/24.*dh;
	      Tyz[i][j][1]+=dt*mu*(Dy+Dz);

	      k=2; 
	      lambda = M_lambda[i][j][k];
	      mu      = M_mu[i][j][k];
	      	      
	      tlu 	= dt*(lambda+2.*mu);
	      tl	= dt*lambda;
	      Dz = -Dp/22.+(-577./528.*Uzt[i][j][1]+201./176.*Uzt[i][j][2]-9./176.*Uzt[i][j][3]+Uzt[i][j][4]/528.)*dh;
	      Dx = (c0*(Uxt[i][j][2]-Uxt[i-1][j][2])-c1*(Uxt[i+1][j][2]-Uxt[i-2][j][2]))*dh;
	      Dy = (c0*(Uyt[i][j][2]-Uyt[i][j-1][2])-c1*(Uyt[i][j+1][2]-Uyt[i][j-2][2]))*dh;
	      Txx[i][j][2] += tlu*Dx + tl*(Dz+Dy);
	      Tyy[i][j][2] += tlu*Dy + tl*(Dx+Dz);
	      Tzz[i][j][2] += tlu*Dz + tl*(Dx+Dy);
	      	      
	      mu = 4./(1./M_mu[i][j][k]+1./M_mu[i][j+1][k]+
	      1./M_mu[i+1][j][k]+1./M_mu[i+1][j+1][k]);
	      Dx = (c0*(Uxt[i][j+1][k]-Uxt[i][j][k])-c1*(Uxt[i][j+2][k]-Uxt[i][j-1][k]))*dh;
	      Dy = (c0*(Uyt[i+1][j][k]-Uyt[i][j][k])-c1*(Uyt[i+2][j][k]-Uyt[i-1][j][k]))*dh;
	      Txy[i][j][k]+=dt*mu*(Dx+Dy);

	      mu = 4./(1./M_mu[i][j][k]+1./M_mu[i+1][j][k]+
				1./M_mu[i][j][k+1]+1./M_mu[i+1][j][k+1]);
	      Dx = (c0*(Uxt[i][j][k+1]-Uxt[i][j][k])-c1*(Uxt[i][j][k+2]-Uxt[i][j][k-1]))*dh;
	      Dz = (c0*(Uzt[i+1][j][k]-Uzt[i][j][k])-c1*(Uzt[i+2][j][k]-Uzt[i-1][j][k]))*dh;
	      Txz[i][j][k]+=dt*mu*(Dx+Dz);

	      mu = 4./(1./M_mu[i][j][k]+1./M_mu[i][j+1][k]+
				1./M_mu[i][j][k+1]+1./M_mu[i][j+1][k+1]);
	      Dy = (c0*(Uyt[i][j][k+1]-Uyt[i][j][k])-c1*(Uyt[i][j][k+2]-Uyt[i][j][k-1]))*dh;
	      Dz = (c0*(Uzt[i][j+1][k]-Uzt[i][j][k])-c1*(Uzt[i][j+2][k]-Uzt[i][j-1][k]))*dh;
	      Tyz[i][j][k]+=dt*mu*(Dy+Dz);
	}
   }
}

void Apply_Surface_velocity(void)
{
   int i,j,k;
   float Dx,Dy,Dz;
   float tbx,tby,tbz;
   float dt,dh,c0,c1;

   dt = elmdls.dt;
   dh = 1./elmdls.dh;
   c0 = 9./8.;
   c1 = 1./24.;

   for(i=elmdls.start[0];i<=elmdls.end[0];i++)
   {
	for(j=elmdls.start[1];j<=elmdls.end[1];j++)
	{
                k=1;
                tbx = (1./M_rho[i][j][k]+1./M_rho[i+1][j][k])*0.5*dt;
                tby = (1./M_rho[i][j][k]+1./M_rho[i][j+1][k])*0.5*dt;
                tbz = (1./M_rho[i][j][k]+1./M_rho[i][j][k+1])*0.5*dt;

		Dx = (c0*(Txx[i+1][j][k]-Txx[i][j][k])  -c1*(Txx[i+2][j][k]-Txx[i-1][j][k]))*dh;
		Dy = (c0*(Txy[i][j][k]  -Txy[i][j-1][k])-c1*(Txy[i][j+1][k]-Txy[i][j-2][k]))*dh;
		Dz = (c0*(Txz[i][j][k]  +Txz[i][j][k])  -c1*(Txz[i][j][k+1]+Txz[i][j][k+1]))*dh;
/*
		Dz = (35./8.*Txz[i][j][1]-35./24.*Txz[i][j][2]+21./40.*Txz[i][j][3]-5./56.*Txz[i][j][4])*dh;
*/
                Uxt[i][j][k]+= tbx*(Dx+Dy+Dz);
				
		Dx = (c0*(Txy[i][j][k]  -Txy[i-1][j][k])-c1*(Txy[i+1][j][k]-Txy[i-2][j][k]))*dh;
		Dy = (c0*(Tyy[i][j+1][k]-Tyy[i][j][k])  -c1*(Tyy[i][j+2][k]-Tyy[i][j-1][k]))*dh;
		Dz = (c0*(Tyz[i][j][k]  +Tyz[i][j][k])  -c1*(Tyz[i][j][k+1]+Tyz[i][j][k+1]))*dh;
/*
		Dz = (35./8.*Tyz[i][j][1]-35./24.*Tyz[i][j][2]+21./40.*Tyz[i][j][3]-5./56.*Tyz[i][j][4])*dh;
*/
                Uyt[i][j][k]+= tby*(Dx+Dy+Dz);

		Dx = (c0*(Txz[i][j][k]  -Txz[i-1][j][k])-c1*(Txz[i+1][j][k]-Txz[i-2][j][k]))*dh;
		Dy = (c0*(Tyz[i][j][k]  -Tyz[i][j-1][k])-c1*(Tyz[i][j+1][k]-Tyz[i][j-2][k]))*dh; 
                Dz = (c0*(Tzz[i][j][k+1]-Tzz[i][j][k])  -c1*(Tzz[i][j][k+2]+Tzz[i][j][k+1]))*dh;
/*
		Dz = (-22.*Tzz[i][j][1]+17.*Tzz[i][j][2]+9.*Tzz[i][j][3]-5.*Tzz[i][j][4]+Tzz[i][j][5])/24.*dh;
		Dp = (-22.*Pfp[i][j][1]+17.*Pfp[i][j][2]+9.*Pfp[i][j][3]-5.*Pfp[i][j][4]+Pfp[i][j][5])/24.*dh;
*/
                Uzt[i][j][k]+=tbz*(Dx+Dy+Dz);

		k=2;

                tbx = (1./M_rho[i][j][k]+1./M_rho[i+1][j][k])*0.5*dt;
                tby = (1./M_rho[i][j][k]+1./M_rho[i][j+1][k])*0.5*dt;
                tbz = (1./M_rho[i][j][k]+1./M_rho[i][j][k+1])*0.5*dt;

		Dx = (c0*(Txx[i+1][j][k]-Txx[i][j][k])  -c1*(Txx[i+2][j][k]-Txx[i-1][j][k]))*dh;
		Dy = (c0*(Txy[i][j][k]  -Txy[i][j-1][k])-c1*(Txy[i][j+1][k]-Txy[i][j-2][k]))*dh;
		Dz = (c0*(Txz[i][j][k]  -Txz[i][j][k-1])-c1*(Txz[i][j][k+1]+Txz[i][j][k-1]))*dh;
/*		
		Dz = (-31./24.*Txz[i][j][1]+29./24.*Txz[i][j][2]-3./40.*Txz[i][j][3]+Txz[i][j][4]/168.)*dh; 
*/
                Uxt[i][j][k]+= tbx*(Dx+Dy+Dz);

		Dx = (c0*(Txy[i][j][k]  -Txy[i-1][j][k])-c1*(Txy[i+1][j][k]-Txy[i-2][j][k]))*dh;
		Dy = (c0*(Tyy[i][j+1][k]-Tyy[i][j][k])  -c1*(Tyy[i][j+2][k]-Tyy[i][j-1][k]))*dh;
		Dz = (c0*(Tyz[i][j][k]  -Tyz[i][j][k-1])-c1*(Tyz[i][j][k+1]+Tyz[i][j][k-1]))*dh;
/*
		Dz = (-31./24.*Tyz[i][j][1]+29./24.*Tyz[i][j][2]-3./40.*Tyz[i][j][3]+Tyz[i][j][4]/168.)*dh;
*/
                Uyt[i][j][k]+= tby*(Dx+Dy+Dz);

		Dx = (c0*(Txz[i][j][k]  -Txz[i-1][j][k])-c1*(Txz[i+1][j][k]-Txz[i-2][j][k]))*dh;
		Dy = (c0*(Tyz[i][j][k]  -Tyz[i][j-1][k])-c1*(Tyz[i][j+1][k]-Tyz[i][j-2][k]))*dh;
		Dz = (c0*(Tzz[i][j][k+1]-Tzz[i][j][k])  -c1*(Tzz[i][j][k+2]-Tzz[i][j][k-1]))*dh;
                Uzt[i][j][k]+=tbz*(Dx+Dy+Dz);
	}
   }
}
