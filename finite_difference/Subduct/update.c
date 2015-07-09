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

#include "update.h"

void Update(void)
{
   int iter,i,j,k;
   int sxr,exr,syr,eyr,szr,ezr,nt,npml;
   float tlu,tl,Dx,Dy,Dz;
   float tbx,tby,tbz,lambda,mu;
   float dt,dh,c0,c1;
   float *tbuf;

   sxr = elmdls.start[0];
   exr = elmdls.end[0];
   syr = elmdls.start[1];
   eyr = elmdls.end[1];
   szr = elmdls.start[2];
   ezr = elmdls.end[2];

   nt = elmdls.nt;
   dt = elmdls.dt;
   dh = 1./elmdls.dh;
   npml = elmdls.npml;
   c0 = 9./8./elmdls.dh;
   c1 = 1./elmdls.dh/24.;

   if( elmdlv.cord_pml[0][0] != -1 ) syr = elmdlv.cord_pml[0][3]+1;
   if( elmdlv.cord_pml[1][0] != -1 ) eyr = elmdlv.cord_pml[1][2]-1;
   if( elmdlv.cord_pml[2][0] != -1 ) sxr = elmdlv.cord_pml[2][1]+1;
   if( elmdlv.cord_pml[3][0] != -1 ) exr = elmdlv.cord_pml[3][0]-1;
   if( elmdlv.cord_pml[4][0] != -1 ) ezr = elmdlv.cord_pml[4][4]-1;
   if( elmdls.freesurface )
	szr = MAX(3, szr);
   else
	if( elmdlv.cord_pml[5][0] != -1 ) szr = elmdlv.cord_pml[5][5]+1;

   k = MAX(elmdls.px*elmdls.py,MAX(elmdls.px*elmdls.pz,elmdls.py*elmdls.pz));

   tbuf = (float *)malloc(2*7*k*sizeof(float));

   for(iter=0;iter<nt;iter++)
   {
	for(i=sxr;i<=exr;i++)
	{
	   for(j=syr;j<=eyr;j++)
	   {
	      for(k=szr;k<=ezr;k++)
	      {
		   tbx = (1./M_rho[i][j][k]+1./M_rho[i+1][j][k])*0.5*dt;
		   tby = (1./M_rho[i][j][k]+1./M_rho[i][j+1][k])*0.5*dt;
		   tbz = (1./M_rho[i][j][k]+1./M_rho[i][j][k+1])*0.5*dt;
		   Dx = (c0*(Txx[i+1][j][k]-Txx[i][j][k])  -c1*(Txx[i+2][j][k]-Txx[i-1][j][k]));
		   Dy = (c0*(Txy[i][j][k]  -Txy[i][j-1][k])-c1*(Txy[i][j+1][k]-Txy[i][j-2][k]));
		   Dz = (c0*(Txz[i][j][k]  -Txz[i][j][k-1])-c1*(Txz[i][j][k+1]-Txz[i][j][k-2]));
		   Uxt[i][j][k] += tbx*(Dx+Dy+Dz); 
				
		   Dx = (c0*(Txy[i][j][k]  -Txy[i-1][j][k])-c1*(Txy[i+1][j][k]-Txy[i-2][j][k]));
		   Dy = (c0*(Tyy[i][j+1][k]-Tyy[i][j][k])  -c1*(Tyy[i][j+2][k]-Tyy[i][j-1][k]));
		   Dz = (c0*(Tyz[i][j][k]  -Tyz[i][j][k-1])-c1*(Tyz[i][j][k+1]-Tyz[i][j][k-2]));
		   Uyt[i][j][k] += tby*(Dx+Dy+Dz);

		   Dx = (c0*(Txz[i][j][k]  -Txz[i-1][j][k])-c1*(Txz[i+1][j][k]-Txz[i-2][j][k]));
		   Dy = (c0*(Tyz[i][j][k]  -Tyz[i][j-1][k])-c1*(Tyz[i][j+1][k]-Tyz[i][j-2][k]));
		   Dz = (c0*(Tzz[i][j][k+1]-Tzz[i][j][k])  -c1*(Tzz[i][j][k+2]-Tzz[i][j][k-1]));
		   Uzt[i][j][k] +=tbz*(Dx+Dy+Dz);
	      }
	   }
	}
	if( elmdlv.apply_pml_velocity != NULL ) elmdlv.apply_pml_velocity();  
	if( elmdlv.apply_surface_velocity != NULL ) elmdlv.apply_surface_velocity();  
	if( elmdls.src_type % 2 == 1 && elmdlv.src != NULL ) elmdlv.src_func(iter);

	communicate(elmdlv.vel,tbuf,3);

	if( elmdlv.save_seismogram ) elmdlv.save_seismogram(iter);

	for(i=sxr;i<=exr;i++)
	{
	   for(j=syr;j<=eyr;j++)
	   {
	     for(k=szr;k<=ezr;k++)
	     {
	       	lambda  = M_lambda[i][j][k];
	       	mu      = M_mu[i][j][k];

		tlu 	= dt*(lambda+2.*mu);
		tl	= dt*lambda;

		Dx = (c0*(Uxt[i][j][k]-Uxt[i-1][j][k])-c1*(Uxt[i+1][j][k]-Uxt[i-2][j][k]));
		Dy = (c0*(Uyt[i][j][k]-Uyt[i][j-1][k])-c1*(Uyt[i][j+1][k]-Uyt[i][j-2][k]));
		Dz = (c0*(Uzt[i][j][k]-Uzt[i][j][k-1])-c1*(Uzt[i][j][k+1]-Uzt[i][j][k-2]));

		Txx[i][j][k]+=tlu*Dx+tl*(Dy+Dz);
		Tyy[i][j][k]+=tlu*Dy+tl*(Dx+Dz);
		Tzz[i][j][k]+=tlu*Dz+tl*(Dx+Dy);

	       	mu = 4./(1./M_mu[i][j][k]+1./M_mu[i][j+1][k]+
				1./M_mu[i+1][j][k]+1./M_mu[i+1][j+1][k]);
		Dx = (c0*(Uxt[i][j+1][k]-Uxt[i][j][k])-c1*(Uxt[i][j+2][k]-Uxt[i][j-1][k]));
		Dy = (c0*(Uyt[i+1][j][k]-Uyt[i][j][k])-c1*(Uyt[i+2][j][k]-Uyt[i-1][j][k]));
		Txy[i][j][k]+=dt*mu*(Dx+Dy);

	       	mu = 4./(1./M_mu[i][j][k]+1./M_mu[i+1][j][k]+
				1./M_mu[i][j][k+1]+1./M_mu[i+1][j][k+1]);
		Dx = (c0*(Uxt[i][j][k+1]-Uxt[i][j][k])-c1*(Uxt[i][j][k+2]-Uxt[i][j][k-1]));
		Dz = (c0*(Uzt[i+1][j][k]-Uzt[i][j][k])-c1*(Uzt[i+2][j][k]-Uzt[i-1][j][k]));
		Txz[i][j][k]+=dt*mu*(Dx+Dz);

	       	mu = 4./(1./M_mu[i][j][k]+1./M_mu[i][j+1][k]+
				1./M_mu[i][j][k+1]+1./M_mu[i][j+1][k+1]);
		Dy = (c0*(Uyt[i][j][k+1]-Uyt[i][j][k])-c1*(Uyt[i][j][k+2]-Uyt[i][j][k-1]));
		Dz = (c0*(Uzt[i][j+1][k]-Uzt[i][j][k])-c1*(Uzt[i][j+2][k]-Uzt[i][j-1][k]));
		Tyz[i][j][k]+=dt*mu*(Dy+Dz);
	     }
	   }
       	}
	if( elmdlv.apply_pml_stress != NULL ) elmdlv.apply_pml_stress();  
	if( elmdlv.apply_surface_stress != NULL ) elmdlv.apply_surface_stress();  

	if( elmdls.src_type % 2 == 0  && elmdlv.src != NULL) elmdlv.src_func(iter);
	communicate(elmdlv.stress,tbuf,6);

	if( elmdls.interimflag ) Print_data(iter+1);
   }
   free(tbuf);
}
