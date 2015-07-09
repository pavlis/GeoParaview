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

#include "communicate_mpi.h"

int get_error_mpi(int evalue)
{
    int sum=0;

    MPI_Allreduce(&evalue, &sum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    return sum;
}

void mem2buf_yz(float ****org,float *buf,int ix,int np)
{
   register int j,k,nz,ip;
   nz = elmdls.pz;

   for(ip=0,k=0;ip<np;ip++)
   {
     for(j=elmdls.start[1];j<=elmdls.end[1];j++)
     {
        memcpy(buf+nz*k,     org[ip][ix][j]+elmdls.start[2],   sizeof(float)*nz);
        memcpy(buf+nz*(k+1), org[ip][ix+1][j]+elmdls.start[2], sizeof(float)*nz);
        k+=2;
     }
   }
}

void buf2mem_yz(float ****org,float *buf,int ix,int np)
{
   register int j,k,nz,ip;
   nz = elmdls.pz;

   for(ip=0,k=0;ip<np;ip++)
   {
     for(j=elmdls.start[1];j<=elmdls.end[1];j++)
     {
        memcpy(org[ip][ix][j]+elmdls.start[2],     buf+nz*k,        sizeof(float)*nz);
        memcpy(org[ip][ix+1][j]+elmdls.start[2],   buf+nz*(k+1),    sizeof(float)*nz);
        k+=2;
     }
   }
}

void mem2buf_xz(float ****org,float *buf,int jy,int np)
{
   register int i,k,nz,ip;
   nz = elmdls.pz;

   for(ip=0,k=0;ip<np;ip++)
   {
     for(i=elmdls.start[0];i<=elmdls.end[0];i++)
     {
        memcpy(buf+nz*k,     org[ip][i][jy]+elmdls.start[2],   sizeof(float)*nz);
        memcpy(buf+nz*(k+1), org[ip][i][jy+1]+elmdls.start[2], sizeof(float)*nz);
        k+=2;
     }
   }
}

void buf2mem_xz(float ****org,float *buf,int jy,int np)
{
   register int i,k,nz,ip;
   nz = elmdls.pz;

   for(ip=0,k=0;ip<np;ip++)
   {
     for(i=elmdls.start[0];i<=elmdls.end[0];i++)
     {
        memcpy(org[ip][i][jy]+elmdls.start[2],     buf+nz*k,        sizeof(float)*nz);
        memcpy(org[ip][i][jy+1]+elmdls.start[2],   buf+nz*(k+1),    sizeof(float)*nz);
        k+=2;
     }
   }
}

void mem2buf_xy(float ****org,float *buf,int kz,int np)
{
   register int i,j,k,ip;

   for(ip=0,k=0;ip<np;ip++)
   {
     for(i=elmdls.start[0];i<=elmdls.end[0];i++)
     {
       for(j=elmdls.start[1];j<=elmdls.end[1];j++)
       {
	  buf[k++] = org[ip][i][j][kz];
	  buf[k++] = org[ip][i][j][kz+1];
       }
     }
   }
}

void buf2mem_xy(float ****org,float *buf,int kz,int np)
{
   register int i,j,k,ip;

   for(ip=0,k=0;ip<np;ip++)
   {
     for(i=elmdls.start[0];i<=elmdls.end[0];i++)
     {
       for(j=elmdls.start[1];j<=elmdls.end[1];j++)
       {
	  org[ip][i][j][kz]   = buf[k++];
	  org[ip][i][j][kz+1] = buf[k++];
       }
     }
   }
}

void comm_basic(float ****value,float *buf,int np,int direction,int area,
	void (*mem2buf)(float ****org,float *buf,int ix,int np),
	void (*buf2mem)(float ****org,float *buf,int ix,int np))
{
   MPI_Status status;
   if( (elmdls.rank[direction]%2) == 1 )
   {
                (*mem2buf)(value,buf,elmdls.start[direction],  np);
                MPI_Send(buf, 2*np*area, MPI_FLOAT, elmdls.leftProc[direction], 0, MPI_COMM_WORLD );
                MPI_Recv(buf, 2*np*area, MPI_FLOAT, elmdls.leftProc[direction], 0, MPI_COMM_WORLD, &status );
                (*buf2mem)(value,buf,elmdls.start[direction]-2,np);

	if( elmdls.rank[direction] != elmdls.nProcs[direction]-1 )
	{
                (*mem2buf)(value,buf,elmdls.end[direction]-1,np);
                MPI_Send(buf, 2*np*area, MPI_FLOAT, elmdls.rightProc[direction], 0, MPI_COMM_WORLD );
                MPI_Recv(buf, 2*np*area, MPI_FLOAT, elmdls.rightProc[direction], 0, MPI_COMM_WORLD, &status );
                (*buf2mem)(value,buf,elmdls.end[direction]+1,np);
	}
   }
   else
   {
	if( elmdls.rank[direction] != elmdls.nProcs[direction]-1 )
	{
                MPI_Recv(buf, 2*np*area, MPI_FLOAT, elmdls.rightProc[direction], 0, MPI_COMM_WORLD, &status );
                (*buf2mem)(value,buf,elmdls.end[direction]+1,np);
                (*mem2buf)(value,buf,elmdls.end[direction]-1,np);
                MPI_Send(buf, 2*np*area, MPI_FLOAT, elmdls.rightProc[direction], 0, MPI_COMM_WORLD );
	}

	if( elmdls.rank[direction] != 0 )
	{
                MPI_Recv(buf, 2*np*area, MPI_FLOAT, elmdls.leftProc[direction], 0, MPI_COMM_WORLD, &status );
                (*buf2mem)(value,buf,elmdls.start[direction]-2,np);
                (*mem2buf)(value,buf,elmdls.start[direction],np);
                MPI_Send(buf, 2*np*area, MPI_FLOAT, elmdls.leftProc[direction], 0, MPI_COMM_WORLD );
	}
   }
}

void communicate(float ****value,float *buf,int np)
{
/* x direction : yz plane */
   comm_basic(value,buf,np,0,elmdls.py*elmdls.pz,(*mem2buf_yz),(*buf2mem_yz));

/* y direction : xz plane */
   comm_basic(value,buf,np,1,elmdls.px*elmdls.pz,(*mem2buf_xz),(*buf2mem_xz));

/* z direction : xy plane */
   comm_basic(value,buf,np,2,elmdls.px*elmdls.py,(*mem2buf_xy),(*buf2mem_xy));
}

void comm_distribute_setup(int mrank)
{
    int			i,j,k;
    char                *ctemp;
    int                 slength,nsrcrec;
    int                 blocklen[2];
    int                 *itemp;
    MPI_Aint            disp[2];
    MPI_Datatype        type[2];
    MPI_Datatype        NewType;

    type[0] = MPI_INT;
    type[1] = MPI_FLOAT;

/*
integer in MDLs
float in MDLs
 */
    blocklen[0] = 57;
    blocklen[1] = 13;

    MPI_Address(&elmdls.rank[0],disp);
    MPI_Address(&elmdls.rho,disp+1);

    disp[1] -= disp[0];
    disp[0] = 0;

    MPI_Type_struct(2, blocklen, disp, type, &NewType);
    MPI_Type_commit(&NewType);

    if( mrank == 0 )
    {
	if( elmdlv.oinfo != NULL ) elmdls.leftProc[0] = strlen(elmdlv.oinfo);
	else elmdls.leftProc[0] = 0;
	if( elmdlv.wd != NULL )    elmdls.leftProc[1] = strlen(elmdlv.wd);
	else elmdls.leftProc[1] = 0;
	if( elmdlv.mfile != NULL ) elmdls.rightProc[0] = strlen(elmdlv.mfile);
	else elmdls.rightProc[0] = 0;
    }

/* Broadcast elmdls */
    MPI_Bcast(&elmdls,1,NewType,0,MPI_COMM_WORLD);

    MPI_Type_free(&NewType);

    elmdls.rank[3] = mrank;
    itemp = (int *)calloc(2*elmdls.nshot,sizeof(int));
    if( mrank == 0 )
    {
	for(i=0;i<elmdls.nshot;i++)
	{
	    itemp[2*i]   = elmdlv.nsrc[i];
	    itemp[2*i+1] = elmdlv.nrec[i];
	}
	MPI_Bcast(itemp,2*elmdls.nshot,MPI_INT,0,MPI_COMM_WORLD);
    }
    else
    {
	elmdlv.nsrc = (int *)calloc(elmdls.nshot, sizeof(int));
	elmdlv.nrec = (int *)calloc(elmdls.nshot, sizeof(int));
	MPI_Bcast(itemp,2*elmdls.nshot,MPI_INT,0,MPI_COMM_WORLD);
	for(i=0;i<elmdls.nshot;i++)
	{
	    elmdlv.nsrc[i]=itemp[2*i];
	    elmdlv.nrec[i]=itemp[2*i+1];
	}
    }
    free(itemp);

    for(i=0,nsrcrec=0;i<elmdls.nshot;i++)
	nsrcrec += elmdlv.nsrc[i] + elmdlv.nrec[i];

    itemp = (int *)calloc(nsrcrec*3+elmdls.nrcmpnt+2*elmdls.shnloc,sizeof(int));

    slength = elmdls.leftProc[0]+elmdls.leftProc[1]+elmdls.rightProc[0]+1;
    ctemp = (char *)calloc(slength,sizeof(char));

    if( mrank == 0 )
    {
        for(i=0,k=0;i<elmdls.nshot;i++)
	{
	    for(j=0;j<elmdlv.nsrc[i];j++)
	    {
		itemp[k++] = elmdlv.sx[i][j];
		itemp[k++] = elmdlv.sy[i][j];
		itemp[k++] = elmdlv.sz[i][j];
	    }
	    for(j=0;j<elmdlv.nrec[i];j++)
	    {
		itemp[k++] = elmdlv.gx[i][j];
		itemp[k++] = elmdlv.gy[i][j];
		itemp[k++] = elmdlv.gz[i][j];
	    }
	}
	for(i=0;i<elmdls.nrcmpnt;i++) itemp[k++] = elmdlv.rcmpnt_list[i];
	for(i=0;i<elmdls.shnloc;i++)
	{
	    itemp[k++] = elmdlv.shxsec[i];
	    itemp[k++] = elmdlv.shgnum[i];
	}
        if( elmdlv.mfile ) sprintf(ctemp,"%s%s%s",elmdlv.oinfo,elmdlv.wd,elmdlv.mfile);
	else		   sprintf(ctemp,"%s%s",elmdlv.oinfo,elmdlv.wd);
    }

/* Broadcast elmdlv */
    MPI_Bcast(itemp,nsrcrec*3+elmdls.nrcmpnt+2*elmdls.shnloc,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(ctemp,slength,MPI_CHAR,0,MPI_COMM_WORLD);

    if( mrank != 0 )
    {
	elmdlv.sx = (int **)calloc(elmdls.nshot,sizeof(int*));
	elmdlv.sy = (int **)calloc(elmdls.nshot,sizeof(int*));
	elmdlv.sz = (int **)calloc(elmdls.nshot,sizeof(int*));
	elmdlv.gx = (int **)calloc(elmdls.nshot,sizeof(int*));
	elmdlv.gy = (int **)calloc(elmdls.nshot,sizeof(int*));
	elmdlv.gz = (int **)calloc(elmdls.nshot,sizeof(int*));

	elmdlv.rcmpnt_list = (int *)calloc(elmdls.nrcmpnt,sizeof(int));
	elmdlv.shxsec      = (int *)calloc(elmdls.shnloc,sizeof(int));
	elmdlv.shgnum      = (int *)calloc(elmdls.shnloc,sizeof(int));
        for(i=0,k=0;i<elmdls.nshot;i++)
	{
	    elmdlv.sx[i] = (int *)calloc(elmdlv.nsrc[i],sizeof(int));
	    elmdlv.sy[i] = (int *)calloc(elmdlv.nsrc[i],sizeof(int));
	    elmdlv.sz[i] = (int *)calloc(elmdlv.nsrc[i],sizeof(int));
	    elmdlv.gx[i] = (int *)calloc(elmdlv.nrec[i],sizeof(int));
	    elmdlv.gy[i] = (int *)calloc(elmdlv.nrec[i],sizeof(int));
	    elmdlv.gz[i] = (int *)calloc(elmdlv.nrec[i],sizeof(int));
	    for(j=0;j<elmdlv.nsrc[i];j++)
	    {
		elmdlv.sx[i][j]=itemp[k++];
		elmdlv.sy[i][j]=itemp[k++];
		elmdlv.sz[i][j]=itemp[k++];
	    }
	    for(j=0;j<elmdlv.nrec[i];j++)
	    {
		elmdlv.gx[i][j]=itemp[k++];
		elmdlv.gy[i][j]=itemp[k++];
		elmdlv.gz[i][j]=itemp[k++];
	    }
	}
	for(i=0;i<elmdls.nrcmpnt;i++) elmdlv.rcmpnt_list[i] = itemp[k++];
	for(i=0;i<elmdls.shnloc;i++)
	{
	    elmdlv.shxsec[i] = itemp[k++];
	    elmdlv.shgnum[i] = itemp[k++];
	}

        elmdlv.oinfo =   (char *)calloc(elmdls.leftProc[0]+1,sizeof(char));
        elmdlv.wd =      (char *)calloc(elmdls.leftProc[1]+1,sizeof(char));
        strncpy(elmdlv.oinfo,  ctemp,elmdls.leftProc[0]);
        strncpy(elmdlv.wd,     ctemp+elmdls.leftProc[0],elmdls.leftProc[1]);
	if( elmdls.rightProc[0] > 0 )
	{
	    elmdlv.mfile =   (char *)calloc(elmdls.rightProc[0]+1,sizeof(char));
	    strncpy(elmdlv.mfile,  ctemp+elmdls.leftProc[0]+elmdls.leftProc[1],elmdls.rightProc[0]);
	}
	else elmdlv.mfile = NULL;
    }

    free(itemp);
    free(ctemp);
}

void comm_write_info(char *tbuf,char *fname,int type,int csize)
{
   int i;
   MPI_Status status;
   if( type == 0 )
   {
	if( (Gfp = echkfopen(fname,"wa",shwarn)) != NULL )
	{
	    fprintf(Gfp,tbuf);
	    fclose(Gfp);
	}
   }
   else
   {
	if( elmdls.rank[3] != 0 )
	    MPI_Send(tbuf,csize,MPI_CHAR,0,1,MPI_COMM_WORLD);
	else
	{
	    if( (Gfp = echkfopen(fname,"wa",shwarn)) != NULL )
	    {
		fprintf(Gfp,tbuf);
		for(i=1;i<elmdls.nProcs[3];i++)
		{
	    	    MPI_Recv(tbuf,csize,MPI_CHAR,i,1,MPI_COMM_WORLD,&status);
		    fprintf(Gfp,tbuf);
		}
		fclose(Gfp);
	    }
	}
   }
}

void comm_write_seismogram(void)
{
    int i,j,k,icmpnt,totalnrec,nrec,*trcv_loc,ncmpnt;
    float **rs_data,*tbuf_rsm,sum;
    MPI_Status status;

    ncmpnt = elmdls.nrcmpnt;

    totalnrec = 0;
    if( elmdls.outstep > 0 ) return;

    if( elmdlv.nrec[elmdls.ishot] > 0 )
    {
/*
        if( elmdlv.res_name ) sprintf(Gfname,"%s/%s.%d",elmdlv.wd,elmdlv.res_name,elmdls.rank[3]);
	else		 
*/
sprintf(Gfname,"%s/seismic_response%d.bin",elmdlv.wd,elmdls.rank[3]);

        if( (Gfp = echkfopen(Gfname,"w",shwarn)) != NULL )
        {
	    for(icmpnt=0;icmpnt<ncmpnt;icmpnt++)
		for(i=0;i<elmdlv.nrec[elmdls.ishot];i++)
			    fwrite(elmdlv.seismogrm[icmpnt][i],sizeof(float),elmdls.nt/elmdls.sample,Gfp);
	    fclose(Gfp);
	}
    }
/* For Gary
    MPI_Reduce(&elmdlv.nrec[elmdls.ishot], &totalnrec, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    if( elmdls.rank[3] == 0 && totalnrec > 0 )
    {
        rs_data = (float **)matrix(0,totalnrec*ncmpnt-1,0,elmdls.nt-1);
        for(i=0;i<elmdlv.nrec[elmdls.ishot];i++)
        {
	    for(icmpnt=0;icmpnt<ncmpnt;icmpnt++)
		memcpy( rs_data[ elmdlv.rec_loc[i]+totalnrec*icmpnt ], elmdlv.seismogrm[icmpnt][i], sizeof(float)*elmdls.nt);
        }

        for(i=1;i<elmdls.nProcs[3];i++)
        {
            MPI_Recv(&nrec,1,MPI_INT,i,1,MPI_COMM_WORLD,&status);
            if( nrec > 0 )
            {
                trcv_loc = (int *)calloc(nrec,sizeof(int));
                tbuf_rsm = (float *)calloc(nrec*elmdls.nt*ncmpnt,sizeof(float));
                MPI_Recv(trcv_loc,nrec,MPI_INT,i,1,MPI_COMM_WORLD,&status);
                MPI_Recv(tbuf_rsm,nrec*elmdls.nt*ncmpnt,MPI_FLOAT,i,1,MPI_COMM_WORLD,&status);
                for(j=0;j<nrec;j++)
                {
		    for(icmpnt=0;icmpnt<ncmpnt;icmpnt++)
			memcpy(rs_data[ trcv_loc[j]+totalnrec*icmpnt ], tbuf_rsm+(j*ncmpnt+icmpnt)*elmdls.nt, sizeof(float)*elmdls.nt);
                }
            }
        }

	if( elmdls.res_type == 2 )
	{
	    for(i=0;i<totalnrec*ncmpnt;i++)
	    {
		for(j=0,sum=0.;j<elmdls.nt;j++)
		{
		    sum += rs_data[i][j];
		    rs_data[i][j] = sum * elmdls.dt;
		}
	    }
	}

        if( elmdlv.res_name ) sprintf(Gfname,"%s/%s",elmdlv.wd,elmdlv.res_name);
	else		 sprintf(Gfname,"%s/seismic_response.bin",elmdlv.wd);

        if( (Gfp = echkfopen(Gfname,"w",shwarn)) != NULL )
        {
	    switch( elmdls.res_file_type )
	    {
		case 1:
			for(i=0;i<elmdls.nt;i++)
			{
			    fprintf(Gfp,"%g ",elmdls.dt*(i+1));
			    for(j=0;j<totalnrec*ncmpnt;j++) fprintf(Gfp,"%g ",rs_data[j][i]);
			    fprintf(Gfp,"\n");
			}
			break;
		case 2:
			for(i=0;i<totalnrec*ncmpnt;i++)
			    fwrite(rs_data[i],sizeof(float),elmdls.nt,Gfp);
			break;
	    }
            fclose(Gfp);
        }
        freemat(&rs_data, 0,totalnrec*ncmpnt-1,0,elmdls.nt);
    }
    else
    {
        nrec = elmdlv.nrec[elmdls.ishot];
        MPI_Send(&nrec,1,MPI_INT,0,1,MPI_COMM_WORLD);
        if( nrec > 0 )
        {
                tbuf_rsm = (float *)calloc(nrec*elmdls.nt*ncmpnt,sizeof(float));
                for(j=0;j<nrec;j++)
                {
		    for(icmpnt=0;icmpnt<ncmpnt;icmpnt++)
                    memcpy(tbuf_rsm+(j*ncmpnt+icmpnt)*elmdls.nt,elmdlv.seismogrm[icmpnt][j],sizeof(float)*elmdls.nt);
                }
                MPI_Send(elmdlv.rec_loc,nrec,MPI_INT,0,1,MPI_COMM_WORLD);
                MPI_Send(tbuf_rsm,nrec*elmdls.nt*ncmpnt,MPI_FLOAT,0,1,MPI_COMM_WORLD);
                free(tbuf_rsm);
        }
    }
*/
}

