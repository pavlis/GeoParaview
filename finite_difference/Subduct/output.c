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

#include "output.h"

void Print_data(int iter)
{
   int i,j;

   if( elmdlv.src != NULL && elmdls.rank[3] == 0 )
/*   if( elmdlv.src != NULL ) */
   {
	sprintf(Gfname,"%s/source.dat",elmdlv.wd);
	if( (Gfp=echkfopen(Gfname,"a",NULL))!= NULL )
	{
	   fprintf(Gfp,"%g %g\n",elmdls.dt*iter,elmdlv.src[iter]);
	   fclose(Gfp);
	}
   }

   if( (iter+elmdls.shstep-elmdls.shstart)%elmdls.shstep==0 && elmdls.shnloc>0  && iter>=elmdls.shstart ) 
   {
	for(i=0;i<elmdls.shnloc;i++)
	{
	   switch( elmdlv.shxsec[i] )
	   {
		case 0:
			out_func(elmdlv.shgnum[i],"yz",iter,(*write_yz));
			break;
		case 1:
			out_func(elmdlv.shgnum[i],"xz",iter,(*write_xz));
			break;
		case 2:
			out_func(elmdlv.shgnum[i],"xy",iter,(*write_xy));
			break;
	   }
	}
   }
}
void write_yz(float ****outdata,int cmp,int gnum)
{
   int j;
   if( (Gfp=echkfopen(Gfname,"wb",shwarn))!= NULL )
   {
	for(j=elmdls.start[1];j<=elmdls.end[1];j++)
	   fwrite(outdata[cmp][gnum][j]+elmdls.start[2],sizeof(float),elmdls.pz,Gfp);
	fclose(Gfp);
   }
}
void write_xy(float ****outdata,int cmp,int gnum)
{
   int i,j;
   if( (Gfp=echkfopen(Gfname,"wb",shwarn))!= NULL )
   {
	 for(i=elmdls.start[0];i<=elmdls.end[0];i++)
	    for(j=elmdls.start[1];j<=elmdls.end[1];j++)
		 fwrite(outdata[cmp][i][j]+gnum,sizeof(float),1,Gfp);
	fclose(Gfp);
   }
}
void write_xz(float ****outdata,int cmp,int gnum)
{
   int i;
   if( (Gfp=echkfopen(Gfname,"wb",shwarn))!= NULL )
   {
	 for(i=elmdls.start[0];i<=elmdls.end[0];i++)
	   fwrite(outdata[cmp][i][gnum]+elmdls.start[2],sizeof(float),elmdls.pz,Gfp);
	fclose(Gfp);
   }
}

void out_func(int gnum,char *plane,int iter,void (*write_func)(float ****outdata,int cmp,int gnum))
{
   int i;
   char *para_name[] = {"rho","lambda","mu"};

   if( iter == -100 )
   {
	for(i=0; i<elmdls.nOutPara; i++)
	{
	    sprintf(Gfname,"%s/%s_%s%d_%d.bin",elmdlv.wd,para_name[elmdls.OutParaInfo[i]],plane,gnum,elmdls.rank[3]);
	    (*write_func)(elmdlv.mpara,elmdls.OutParaInfo[i],gnum);
	}
	return;
   }
   if( elmdls.shcmpnt & SN_X )
   {
	sprintf(Gfname,"%s/vx_%s%d_%d_%d.bin",elmdlv.wd,plane,gnum,elmdls.rank[3],iter);
	(*write_func)(elmdlv.vel,0,gnum);
   }
   if( elmdls.shcmpnt & SN_Y )
   {
	sprintf(Gfname,"%s/vy_%s%d_%d_%d.bin",elmdlv.wd,plane,gnum,elmdls.rank[3],iter);
	(*write_func)(elmdlv.vel,1,gnum);
   }
   if( elmdls.shcmpnt & SN_Z )
   {
	sprintf(Gfname,"%s/vz_%s%d_%d_%d.bin",elmdlv.wd,plane,gnum,elmdls.rank[3],iter);
	(*write_func)(elmdlv.vel,2,gnum);
   }
}

void Write_mdl_parameter(void)
{
    int i;

    for(i=0;i<elmdls.shnloc;i++)
    {
	switch( elmdlv.shxsec[i] )
	{
	    case 0:
		out_func(elmdlv.shgnum[i],"yz",-100,(*write_yz));
		break;
	    case 1:
		out_func(elmdlv.shgnum[i],"xz",-100,(*write_xz));
		break;
	    case 2:
		out_func(elmdlv.shgnum[i],"xy",-100,(*write_xy));
		break;
	}
    }
}

float value_seismogram(int i,int j,int k,int cmpnt)
{
   float c0,c1,value;
   c0 = 9./8./elmdls.dh;
   c1 = 1./elmdls.dh/24.;

   value = 0.;
   switch(cmpnt)
   {
	case RC_X:
		value = Uxt[i][j][k];
		break;
	case RC_Y:
		value = Uyt[i][j][k];
		break;
	case RC_D:
		value  = c0*(Uxt[i][j][k]-Uxt[i-1][j][k])-c1*(Uxt[i+1][j][k]-Uxt[i-2][j][k]);
		value += c0*(Uyt[i][j][k]-Uyt[i][j-1][k])-c1*(Uyt[i][j+1][k]-Uyt[i][j-2][k]);
		value += c0*(Uzt[i][j][k]-Uzt[i][j][k-1])-c1*(Uzt[i][j][k+1]-Uzt[i][j][k-2]);
		break;
	case RC_CX:
		value  = c0*(Uzt[i][j+1][k]-Uzt[i][j][k])-c1*(Uzt[i][j+2][k]-Uzt[i][j-1][k]);
		value -= c0*(Uyt[i][j][k+1]-Uyt[i][j][k])-c1*(Uyt[i][j][k+2]-Uyt[i][j][k-1]);
		break;
	case RC_CY:
		value  = c0*(Uxt[i][j][k+1]-Uxt[i][j][k])-c1*(Uxt[i][j][k+2]-Uxt[i][j][k-1]);
		value -= c0*(Uzt[i+1][j][k]-Uzt[i][j][k])-c1*(Uzt[i+2][j][k]-Uzt[i-1][j][k]);
		break;
	case RC_CZ:
		value  = c0*(Uyt[i+1][j][k]-Uyt[i][j][k])-c1*(Uyt[i+2][j][k]-Uyt[i-1][j][k]);
		value -= c0*(Uxt[i][j+1][k]-Uxt[i][j][k])-c1*(Uxt[i][j+2][k]-Uxt[i][j-1][k]);
		break;
	case RC_Z:
	default:
		value = Uzt[i][j][k];
		break;
   }
   return value;
}

void Save_seismogram(int iter)
{
   int ir,i,j,k,icmpnt;

   if( iter%elmdls.sample != 0 ) return;

   iter /= elmdls.sample;
   for(ir=0;ir<elmdlv.irec;ir++)
   {
        i = elmdlv.recx[ir];
        j = elmdlv.recy[ir];
        k = elmdlv.recz[ir];
	for(icmpnt=0;icmpnt<elmdls.nrcmpnt;icmpnt++)
            elmdlv.seismogrm[icmpnt][ir][iter] = value_seismogram(i,j,k,elmdlv.rcmpnt_list[icmpnt]);
   }
}

void Save_seismogram_step(int iter)
{
   static int istep=0;
   int ir,i,j,k,icmpnt;

   if( iter == 0 )
   {
	sprintf(Gfname,"%s/seismic_response%d.bin",elmdlv.wd,elmdls.rank[3]);
        if( (Gfp = echkfopen(Gfname,"w",shwarn)) != NULL ) fclose(Gfp);
   }

   if( iter%elmdls.sample != 0 )
   {
	if( (iter!=elmdls.nt-1) || istep == 0  ) return;
	
if( elmdls.rank[3] == 0 ) printf("   istep:%5d <- %5d\n",istep,iter);
	sprintf(Gfname,"%s/seismic_response%d.bin",elmdlv.wd,elmdls.rank[3]);
        if( (Gfp = echkfopen(Gfname,"a",shwarn)) != NULL )
        {
	    for(icmpnt=0; icmpnt<elmdls.nrcmpnt; icmpnt++)
		for(i=0;i<elmdlv.nrec[elmdls.ishot];i++)
			    fwrite(elmdlv.seismogrm[icmpnt][i],sizeof(float),istep,Gfp);
	    fclose(Gfp);
	}
	return;
   }

   for(ir=0;ir<elmdlv.irec;ir++)
   {
        i = elmdlv.recx[ir];
        j = elmdlv.recy[ir];
        k = elmdlv.recz[ir];
	for(icmpnt=0;icmpnt<elmdls.nrcmpnt;icmpnt++)
            elmdlv.seismogrm[icmpnt][ir][istep] = value_seismogram(i,j,k,elmdlv.rcmpnt_list[icmpnt]);
   }
   istep++;
   if( (istep==elmdls.outstep) || (iter==elmdls.nt-1) )
   {
	sprintf(Gfname,"%s/seismic_response%d.bin",elmdlv.wd,elmdls.rank[3]);

if( elmdls.rank[3] == 0 ) printf("istep:%5d <- %5d\n",istep,iter);
        if( (Gfp = echkfopen(Gfname,"a",shwarn)) != NULL )
        {
	    for(icmpnt=0; icmpnt<elmdls.nrcmpnt; icmpnt++)
		for(i=0;i<elmdlv.nrec[elmdls.ishot];i++)
			    fwrite(elmdlv.seismogrm[icmpnt][i],sizeof(float),istep,Gfp);
	    fclose(Gfp);
	}
	istep = 0;
   }
}
