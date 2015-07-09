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

#include "initialize.h"

void Initialize(int n1,int n2,int n3)
{
   int i,j,k,l,nx,ny,nz,nrec,npml;
   float rho,lambda,mu;
   int sizes[] = {120,130,130,70};
   float ratio;

   nx = elmdls.nx;
   ny = elmdls.ny;
   nz = elmdls.nz;
   npml = elmdls.npml;

   set_priority(nx,ny,nz,&i,&j,&k);

   set_rank_connection(n1,n2,n3,i,j,k);

   decompose_domain(0,nx);
   decompose_domain(1,ny);
/* For Gary
   decompose_domain2(2,nz,sizes);
*/
   decompose_domain(2,nz);

   elmdls.px = elmdls.end[0] - elmdls.start[0] + 1;
   elmdls.py = elmdls.end[1] - elmdls.start[1] + 1;
   elmdls.pz = elmdls.end[2] - elmdls.start[2] + 1;

   set_boundary(nx,ny,nz,npml);

   set_snapshot_location();

   if( elmdls.interimflag )
   {
	if( elmdls.rank[3]==0 )
	{
	    sprintf(Gfname,"%s/source.dat",elmdlv.wd);
	    if( (Gfp=echkfopen(Gfname,"wt",shwarn))!= NULL ) fclose(Gfp);
	    sprintf(Gfname,"%s/%s",elmdlv.wd,elmdlv.oinfo);
	    if( (Gfp=echkfopen(Gfname,"wt",shwarn))!= NULL ) fclose(Gfp);
	}
   }
}

void Initialize_fp(void)
{
    elmdlv.isrc = elmdlv.nsrc[ elmdls.ishot ];
    elmdlv.srcx = elmdlv.sx[ elmdls.ishot ];
    elmdlv.srcy = elmdlv.sy[ elmdls.ishot ];
    elmdlv.srcz = elmdlv.sz[ elmdls.ishot ];

    elmdlv.irec = elmdlv.nrec[ elmdls.ishot ];
    elmdlv.recx = elmdlv.gx[ elmdls.ishot ];
    elmdlv.recy = elmdlv.gy[ elmdls.ishot ];
    elmdlv.recz = elmdlv.gz[ elmdls.ishot ];
    elmdlv.rec_loc = elmdlv.rloc[ elmdls.ishot ];

    Make_source();

    if( elmdlv.irec > 0 )
    {
	if( elmdls.outstep == 0 )
	{
	    elmdlv.seismogrm = (float ***)cubic(0,elmdls.nrcmpnt-1,0,elmdlv.irec-1,0,elmdls.nt/elmdls.sample-1);
	    elmdlv.save_seismogram = (*Save_seismogram);
	}
	else
	{
	    elmdlv.seismogrm = (float ***)cubic(0,elmdls.nrcmpnt-1,0,elmdlv.irec-1,0,elmdls.outstep-1);
	    elmdlv.save_seismogram = (*Save_seismogram_step);
	}
    }
    else elmdlv.save_seismogram = NULL;
}

int getnumber(int *ilist,int nlist)
{
    int i,j,icheck,num;

    num = 1;

    for(i=1;i<nlist;i++)
    {
	for(j=i-1,icheck=0; j>=0; j--)
	{
	    if( ilist[i] == ilist[j] ) 
	    {
		icheck=1;
		break;
	    }
	}
	
	if( icheck == 0 ) num++;
    }
    return num;
}

void Print_Configuration(void)
{
   int i,j,csize,rcheck;
   int xnum,ynum,znum;
   char *tbuf;
   
   csize = 450*elmdls.nProcs[3];
   tbuf = (char *)calloc(csize,sizeof(char));

   if( elmdls.rank[3] == 0 ) 
   {
	sprintf(tbuf,"%sModel=%3d,%3d,%3d,%3d,%3d,%3d #nx,ny,nz,nt,sampling\n",tbuf, elmdls.nx,elmdls.ny,elmdls.nz,elmdls.nt,elmdls.sample,elmdls.outstep); 
	sprintf(tbuf,"%snProcs=%3d,%3d,%3d,%3d\n",tbuf, elmdls.nProcs[3],
			elmdls.nProcs[0],elmdls.nProcs[1],elmdls.nProcs[2]);
	sprintf(tbuf,"%ssx=",tbuf);
	for(i=0;i<elmdls.nshot;i++) 
	{
	    sprintf(tbuf,"%s%3d",tbuf, elmdlv.sx[i][0]);
	    for(j=1;j<elmdlv.nsrc[i];j++) sprintf(tbuf,"%s,%3d",tbuf, elmdlv.sx[i][j]);
	    if( i+1 < elmdls.nshot ) sprintf(tbuf,"%s,\t",tbuf);
	}
	sprintf(tbuf,"%s\n",tbuf);
	
	sprintf(tbuf,"%ssy=",tbuf);
	for(i=0;i<elmdls.nshot;i++) 
	{
	    sprintf(tbuf,"%s%3d",tbuf, elmdlv.sy[i][0]);
	    for(j=1;j<elmdlv.nsrc[i];j++) sprintf(tbuf,"%s,%3d",tbuf, elmdlv.sy[i][j]);
	    if( i+1 < elmdls.nshot ) sprintf(tbuf,"%s,\t",tbuf);
	}
	sprintf(tbuf,"%s\n",tbuf);
	
	sprintf(tbuf,"%ssz=",tbuf);
	for(i=0;i<elmdls.nshot;i++) 
	{
	    sprintf(tbuf,"%s%3d",tbuf, elmdlv.sz[i][0]);
	    for(j=1;j<elmdlv.nsrc[i];j++) sprintf(tbuf,"%s,%3d",tbuf, elmdlv.sz[i][j]);
	    if( i+1 < elmdls.nshot ) sprintf(tbuf,"%s,\t",tbuf);
	}

	sprintf(tbuf,"%s\nReceiver number=%d",tbuf,elmdlv.tnrec[0]);
	for(i=1;i<elmdls.nshot;i++) sprintf(tbuf,"%s, %d",tbuf,elmdlv.tnrec[i]);
	sprintf(tbuf,"%s\n",tbuf);
	
	if( elmdls.nOutPara > 0 )
       	{
	    sprintf(tbuf,"%sParameterOut=",tbuf);
	    for(i=0;i<elmdls.nOutPara-1;i++)
		sprintf(tbuf,"%s%d,",tbuf,elmdls.OutParaInfo[i]);
	    sprintf(tbuf,"%s%d\n",tbuf,elmdls.OutParaInfo[i]);
	}
	sprintf(tbuf,"%sShotInfo=%5d,%5d\n",tbuf,elmdls.shstart,elmdls.shstep);
	sprintf(tbuf,"%sShotComponent=%3d #%d(X),%d(Y),%d(Z) component\n",tbuf,elmdls.shcmpnt,SN_X,SN_Y,SN_Z);
	sprintf(tbuf,"%sRecvComponent=%3d\n",tbuf,elmdls.nrcmpnt); 
   }

   sprintf(tbuf,"%sProcInfo%-4d=%4d,%4d,%4d,%4d,%4d,%4d,%4d,%4d,%4d\n",tbuf,
	elmdls.rank[3],elmdls.start[0],elmdls.end[0],
	              elmdls.start[1],elmdls.end[1],
		      elmdls.start[2],elmdls.end[2],
	              elmdls.rank[0],elmdls.rank[1],elmdls.rank[2]);
/*
   if( (elmdls.shnloc>0) && (elmdls.rank[3]==0) ) 
   {
	sprintf(tbuf,"ShotLoc=%d_%d",elmdlv.shxsec[0],elmdlv.shgnum[0]);
	for(i=1;i<elmdls.shnloc;i++)
	{
		sprintf(tbuf,"%s,%d_%d",tbuf,elmdlv.shxsec[i],elmdlv.shgnum[i]);
	}
	sprintf(tbuf,"%s\n",tbuf);
   }
*/
   if( elmdls.shnloc > 0 ) 
   {
	sprintf(tbuf,"%sShotLoc%d=%d_%d",tbuf,elmdls.rank[3],elmdlv.shxsec[0],elmdlv.shgnum[0]);
	for(i=1;i<elmdls.shnloc;i++)
	{
		sprintf(tbuf,"%s,%d_%d",tbuf,elmdlv.shxsec[i],elmdlv.shgnum[i]);
	}
	sprintf(tbuf,"%s\n",tbuf);
   }

   if( elmdls.rtype != 0 )
   {
	for(i=0,rcheck=0;i<elmdls.nshot;i++) if( elmdlv.nrec[i] > 0 ) rcheck++;
	if( rcheck > 0 )
	{
	    for(i=0;i<elmdls.nshot;i++) if( elmdlv.nrec[i] > 0 ) break;

	    xnum = getnumber(elmdlv.gx[i],elmdlv.nrec[i]);
	    ynum = getnumber(elmdlv.gy[i],elmdlv.nrec[i]);
	    znum = getnumber(elmdlv.gz[i],elmdlv.nrec[i]);
	    sprintf(tbuf,"%sRecv%d=%d,%d,%d,%d,%d",tbuf,elmdls.rank[3],i+1,elmdlv.nrec[i],xnum,ynum,znum);
	    i++;
	    for(;i<elmdls.nshot;i++)
		if( elmdlv.nrec[i] > 0 )
		{
		    xnum = getnumber(elmdlv.gx[i],elmdlv.nrec[i]);
		    ynum = getnumber(elmdlv.gy[i],elmdlv.nrec[i]);
		    znum = getnumber(elmdlv.gz[i],elmdlv.nrec[i]);

		    sprintf(tbuf,"%s,%d,%d,%d,%d,%d",tbuf,elmdlv.nrec[i],xnum,ynum,znum);
		}
	    sprintf(tbuf,"%s\n",tbuf);
	}
   }


   sprintf(Gfname,"%s/%s",elmdlv.wd,elmdlv.oinfo);
   comm_write_info(tbuf,Gfname,1,csize);

   free(tbuf); 
}

void Alloc_Mem_Func(void)
{
   elmdlv.vel    = (float ****)quadri(0,2,elmdls.start[0]-2,elmdls.end[0]+2,elmdls.start[1]-2,elmdls.end[1]+2,elmdls.start[2]-2,elmdls.end[2]+2);
   elmdlv.stress = (float ****)quadri(0,5,elmdls.start[0]-2,elmdls.end[0]+2,elmdls.start[1]-2,elmdls.end[1]+2,elmdls.start[2]-2,elmdls.end[2]+2);
   elmdlv.mpara  = (float ****)quadri(0,2,elmdls.start[0]  ,elmdls.end[0]+1,elmdls.start[1]  ,elmdls.end[1]+1,elmdls.start[2]  ,elmdls.end[2]+1);
}

void Init_model_parameter(void)
{
   int i,j,k;

   switch( elmdls.model )
   {
     case 2:
	Make_2layer();
	break;
     case 8:
	Make_3D_Model();
	break;
     case 11:
       	MakePREM();
	break;
     case 12:
       	MakePREMsubduct();
	break;
     case 1:
     default:
	if( elmdls.rho == 0. )
	{
	    elmdls.rho= 1.8e3;
	    elmdls.Vs = 0.9e3;
	    elmdls.Vp = 1.5e3;
	}
	elmdls.mu     = elmdls.Vs*elmdls.Vs*elmdls.rho;
	elmdls.lambda = elmdls.Vp*elmdls.Vp*elmdls.rho - 2.*elmdls.mu;
	for(i=elmdls.start[0];i<=elmdls.end[0]+1;i++)
	{
	   for(j=elmdls.start[1];j<=elmdls.end[1]+1;j++)
	   {
	      for(k=elmdls.start[2];k<=elmdls.end[2]+1;k++)
	      {
		   elmdlv.mpara[0][i][j][k] = elmdls.rho;
		   elmdlv.mpara[1][i][j][k] = elmdls.lambda;
		   elmdlv.mpara[2][i][j][k] = elmdls.mu;
	      }
	   }
	}
	break;
   }

   for(i=elmdls.start[0];i<=elmdls.end[0]+1;i++)
   {
	for(j=elmdls.start[1];j<=elmdls.end[1]+1;j++)
	{
	    for(k=elmdls.start[2];k<=elmdls.end[2]+1;k++)
	    {
		    if( M_rho[i][j][k] * M_lambda[i][j][k] * M_mu[i][j][k] == 0. )
			    sherror("Check model parameters\n");
	    }
	}
   }
}

int check_pml_region(int sx,int ex,int sy,int ey,int sz,int ez,int *cord)
{
   int x1,x2,y1,y2,z1,z2;

   cord[0] = -1;

   x1 = MAX(sx,elmdls.start[0]);
   x2 = MIN(ex,elmdls.end[0]);
   if( x1 > x2 ) return 0;

   y1 = MAX(sy,elmdls.start[1]);
   y2 = MIN(ey,elmdls.end[1]);
   if( y1 > y2 ) return 0;

   z1 = MAX(sz,elmdls.start[2]);
   z2 = MIN(ez,elmdls.end[2]);
   if( z1 > z2 ) return 0;

   cord[0] = x1; cord[1] = x2;
   cord[2] = y1; cord[3] = y2;
   cord[4] = z1; cord[5] = z2;

   return 1;
}

/*BBB
void set_priority(int nx,int ny,int nz, int *p1,int *p2)
*/
void set_priority(int nx,int ny,int nz, int *p1,int *p2,int *p3)
{
#define XD 0
#define YD 1
#define ZD 2
   if( nx < ny )
   {
	if( nz <= nx )
	{
	   *p1 = YD;
	   *p2 = XD;
	   *p3 = ZD;
	}
	else
	{
	   if( nz <= ny )
	   {
		*p1 = YD;
		*p2 = ZD;
		*p3 = XD;
	   }
	   else
	   {
		*p1 = ZD;
		*p2 = YD;
		*p3 = XD;
	   }
	}
   }
   else
   {
	if( nz <= ny )
	{
	   *p1 = XD;
	   *p2 = YD;
	   *p3 = ZD;
	}
	else
	{
	   if( nz <= nx )
	   {
		*p1 = XD;
		*p2 = ZD;
		*p3 = YD;
	   }
	   else
	   {
		*p1 = ZD;
		*p2 = XD;
		*p3 = YD;
	   }
	}
   }
#undef XD
#undef YD
#undef ZD
}


/*BBB
void set_rank_connection(int n,int k,int j)
*/
void set_rank_connection(int n1,int n2,int n3,int i,int j,int k)
{
/* divide model by # of Proc */
/*BBB
   elmdls.nProcs[0] = n/3;
   elmdls.nProcs[1] = n/3;
   elmdls.nProcs[2] = n/3;
   if( n%3 == 1 ) elmdls.nProcs[k]+=1;
   else if( n%3 == 2 ) 
   {
	elmdls.nProcs[k] += 1;
	elmdls.nProcs[j] += 1;
   }
   elmdls.nProcs[0] = pow(2.,elmdls.nProcs[0]);
   elmdls.nProcs[1] = pow(2.,elmdls.nProcs[1]);
   elmdls.nProcs[2] = pow(2.,elmdls.nProcs[2]);
*/
   elmdls.nProcs[i] = n3;
   elmdls.nProcs[j] = n2;
   elmdls.nProcs[k] = n1;

   elmdls.rank[0] =   elmdls.rank[3]%elmdls.nProcs[0];
   elmdls.rank[1]  = (elmdls.rank[3]/elmdls.nProcs[0])%elmdls.nProcs[1];
   elmdls.rank[2]  = (elmdls.rank[3]/elmdls.nProcs[0]/elmdls.nProcs[1])%elmdls.nProcs[2];

   elmdls.leftProc[0]  = (elmdls.rank[0] == 0) ?                 MPI_PROC_NULL : elmdls.rank[3]-1;
   elmdls.rightProc[0] = (elmdls.rank[0] == elmdls.nProcs[0]-1) ? MPI_PROC_NULL : elmdls.rank[3]+1;

   elmdls.leftProc[1]  = (elmdls.rank[1] == 0) ?                 MPI_PROC_NULL : elmdls.rank[3]-elmdls.nProcs[0];
   elmdls.rightProc[1] = (elmdls.rank[1] == elmdls.nProcs[1]-1) ? MPI_PROC_NULL : elmdls.rank[3]+elmdls.nProcs[0];

   elmdls.leftProc[2]  = (elmdls.rank[2] == 0) ?                 MPI_PROC_NULL : elmdls.rank[3]-elmdls.nProcs[0]*elmdls.nProcs[1];
   elmdls.rightProc[2] = (elmdls.rank[2] == elmdls.nProcs[2]-1) ? MPI_PROC_NULL : elmdls.rank[3]+elmdls.nProcs[0]*elmdls.nProcs[1];
}

void decompose_domain(int cord,int ngrid)
{
   int igrid;

   igrid = ngrid / elmdls.nProcs[cord];
   if( elmdls.rank[cord] < (ngrid%elmdls.nProcs[cord]) )
   {
	igrid++;
	elmdls.start[cord] = elmdls.rank[cord]*igrid + 1;
   }
   else
	elmdls.start[cord] = elmdls.rank[cord]*igrid + 1 + (ngrid%elmdls.nProcs[cord]);
   elmdls.end[cord] = elmdls.start[cord] + igrid -1;
}

void decompose_domain2(int cord,int ngrid,int *sizes)
{
   int i,igrid;

   for(i=0,igrid=0;i<elmdls.rank[cord];i++) igrid += sizes[i];

   elmdls.start[cord] = igrid + 1;
   elmdls.end[cord] = igrid+sizes[i]; 
}

void Set_source(void)
{
    int ishot,i,j;
    int *pn,*px,*py,*pz;

    pn	= (int *)calloc(elmdls.nshot,sizeof(int));
    elmdlv.sloc	= (int **)calloc(elmdls.nshot,sizeof(int*));

    if( elmdls.src_type == 3 || elmdls.src_type == 4 )
    {
	for(ishot=0;ishot<elmdls.nshot;ishot++)
	{
	    pn[ishot] = 1;
	    px 	= (int *)calloc(1,sizeof(int));
	    py 	= (int *)calloc(1,sizeof(int));
	    pz 	= (int *)calloc(1,sizeof(int));
	    elmdlv.sloc[ishot]	= (int *)calloc(1,sizeof(int));

	    px[0] = elmdlv.sx[ishot][0];
	    py[0] = elmdlv.sy[ishot][0];
	    pz[0] = elmdlv.sz[ishot][0];
	    elmdlv.sloc[ishot][0] = 0;

	    free(elmdlv.sx[ishot]);
	    free(elmdlv.sy[ishot]);
	    free(elmdlv.sz[ishot]);
	    elmdlv.sx[ishot] = px;
	    elmdlv.sy[ishot] = py;
	    elmdlv.sz[ishot] = pz;
	}
       	free(elmdlv.nsrc);
       	elmdlv.nsrc = pn;
	return;
    }

    for(ishot=0;ishot<elmdls.nshot;ishot++)
    {
	for(i=0,j=0;i<elmdlv.nsrc[ishot];i++)
	    if( Is_inside3D(elmdlv.sx[ishot][i],elmdlv.sy[ishot][i],elmdlv.sz[ishot][i],
			elmdls.start[0],elmdls.start[1],elmdls.start[2],
			elmdls.end[0],elmdls.end[1],elmdls.end[2]) )
		j++;

	pn[ishot] = j;
	if( j > 0 )
	{
	    px 	= (int *)calloc(j,sizeof(int));
	    py 	= (int *)calloc(j,sizeof(int));
	    pz 	= (int *)calloc(j,sizeof(int));
	    elmdlv.sloc[ishot]	= (int *)calloc(j,sizeof(int));

	    for(i=0,j=0;i<elmdlv.nsrc[ishot];i++)
	    {
		if( Is_inside3D(elmdlv.sx[ishot][i],elmdlv.sy[ishot][i],elmdlv.sz[ishot][i],
			elmdls.start[0],elmdls.start[1],elmdls.start[2],
			elmdls.end[0],elmdls.end[1],elmdls.end[2]) )
		{
		    px[j] = elmdlv.sx[ishot][i];
		    py[j] = elmdlv.sy[ishot][i];
		    pz[j] = elmdlv.sz[ishot][i];
		    elmdlv.sloc[ishot][j] = i;
		    j++;
		}
	    }
	    free(elmdlv.sx[ishot]);
	    free(elmdlv.sy[ishot]);
	    free(elmdlv.sz[ishot]);
	    elmdlv.sx[ishot] = px;
	    elmdlv.sy[ishot] = py;
	    elmdlv.sz[ishot] = pz;
	}
	else 
	{
	    free(elmdlv.sx[ishot]);
	    free(elmdlv.sy[ishot]);
	    free(elmdlv.sz[ishot]);
	    elmdlv.sx[ishot] = NULL;
	    elmdlv.sy[ishot] = NULL;
	    elmdlv.sz[ishot] = NULL;
	}
    }
    free(elmdlv.nsrc);
    elmdlv.nsrc = pn;
}

void Set_receiver(void)
{
    int ishot,i,j;
    int *pn,*px,*py,*pz;

    pn	= (int *)calloc(elmdls.nshot,sizeof(int));
    elmdlv.rloc	= (int **)calloc(elmdls.nshot,sizeof(int*));

    for(ishot=0;ishot<elmdls.nshot;ishot++)
    {
	for(i=0,j=0;i<elmdlv.nrec[ishot];i++)
	    if( Is_inside3D(elmdlv.gx[ishot][i],elmdlv.gy[ishot][i],elmdlv.gz[ishot][i],
			elmdls.start[0],elmdls.start[1],elmdls.start[2],
			elmdls.end[0],elmdls.end[1],elmdls.end[2]) )
		j++;

	pn[ishot] = j;
	if( j > 0 )
	{
	    px 	= (int *)calloc(j,sizeof(int));
	    py 	= (int *)calloc(j,sizeof(int));
	    pz 	= (int *)calloc(j,sizeof(int));
	    elmdlv.rloc[ishot]	= (int *)calloc(j,sizeof(int));

	    for(i=0,j=0;i<elmdlv.nrec[ishot];i++)
	    {
		if( Is_inside3D(elmdlv.gx[ishot][i],elmdlv.gy[ishot][i],elmdlv.gz[ishot][i],
			elmdls.start[0],elmdls.start[1],elmdls.start[2],
			elmdls.end[0],elmdls.end[1],elmdls.end[2]) )
		{
		    px[j] = elmdlv.gx[ishot][i];
		    py[j] = elmdlv.gy[ishot][i];
		    pz[j] = elmdlv.gz[ishot][i];
		    elmdlv.rloc[ishot][j] = i;
		    j++;
		}
	    }
	    free(elmdlv.gx[ishot]);
	    free(elmdlv.gy[ishot]);
	    free(elmdlv.gz[ishot]);
	    elmdlv.gx[ishot] = px;
	    elmdlv.gy[ishot] = py;
	    elmdlv.gz[ishot] = pz;
	}
	else 
	{
	    free(elmdlv.gx[ishot]);
	    free(elmdlv.gy[ishot]);
	    free(elmdlv.gz[ishot]);
	    elmdlv.gx[ishot] = NULL;
	    elmdlv.gy[ishot] = NULL;
	    elmdlv.gz[ishot] = NULL;
	}
    }
    free(elmdlv.nrec);
    elmdlv.nrec = pn;
}

void set_boundary(int nx,int ny,int nz,int npml)
{
    int i,nsides,l,sz;

    if( elmdls.freesurface ) 	nsides = 5, sz=3;
    else				nsides = 6, sz=1;

    elmdlv.cord_pml = (int **)calloc(nsides,sizeof(int *));
    for(i=0;i<nsides;i++)
	   elmdlv.cord_pml[i] = (int *)calloc(6,sizeof(int));

    if( !elmdls.rank[2] && elmdls.freesurface )
    {
	elmdlv.apply_surface_stress   = (*Apply_Surface_stress);
	elmdlv.apply_surface_velocity = (*Apply_Surface_velocity);
    }
    else
    {
	elmdlv.apply_surface_stress   = NULL;
	elmdlv.apply_surface_velocity = NULL;
    }
}

void set_snapshot_location(void)
{
   int i,j,k;
   int *xsction_cord,*xsction_grid;

   /* Prepare Snapshot location */
   /*
   for(i=0,j=0;i<elmdls.shnloc;i++)
   {
	if( elmdls.start[ elmdlv.shxsec[i] ] <= elmdlv.shgnum[i] && 
	    elmdls.end  [ elmdlv.shxsec[i] ] >= elmdlv.shgnum[i])
		   j++;
   }
   k = j;
   if( k > 0 )
   {
        xsction_cord = (int *)malloc(k*sizeof(int));
        xsction_grid = (int *)malloc(k*sizeof(int));
   	for(i=0,j=0;i<elmdls.shnloc;i++)
	   if( elmdls.start[ elmdlv.shxsec[i] ] <= elmdlv.shgnum[i] && 
	       elmdls.end  [ elmdlv.shxsec[i] ] >= elmdlv.shgnum[i])
	   {
		   xsction_cord[j] = elmdlv.shxsec[i];
		   xsction_grid[j] = elmdlv.shgnum[i];
		   j++;
	   }
	free(elmdlv.shxsec);
	free(elmdlv.shgnum);
	elmdlv.shxsec = xsction_cord;
	elmdlv.shgnum = xsction_grid;
   }
   else if( elmdls.shnloc > 0 )
   {
	free(elmdlv.shxsec);
	free(elmdlv.shgnum);
	elmdlv.shxsec = NULL;
	elmdlv.shgnum = NULL;
   }
   elmdls.shnloc = k;
   */
   	for(i=0,j=0;i<elmdls.shnloc;i++)
	   if( elmdls.start[ elmdlv.shxsec[i] ] > elmdlv.shgnum[i] || 
	       elmdls.end  [ elmdlv.shxsec[i] ] < elmdlv.shgnum[i])
	   {
		   elmdlv.shxsec[i] = 5;
		   elmdlv.shgnum[i] = 0;
	   }
}
