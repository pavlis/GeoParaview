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

#include "read_setup.h"

int Read_setup(char *fname)
{
    int i,j,ix,iy,iz,n,nxrec,nyrec,nzrec,ntr;
    int *itmp,*tgx,*tgy,*tgz;
    char  **shdata,*tmp;

    n = Read_para_file(fname);
    /* Check_para(); */

    /* Get model parameters */
    Get_par_int("nx",&elmdls.nx);
    Get_par_int("ny",&elmdls.ny);
    Get_par_int("nz",&elmdls.nz);

    /* Get source/receiver locations */
    if( !Get_par_int("nshot",&elmdls.nshot) ) sherror("Number of shot is not defined\n");
    if( elmdls.nshot < 0 )
    {
	shwarn("Wrong value: Number of shot, Force to 1\n");
	elmdls.nshot = 1;
    }

    elmdlv.sx = (int **)calloc(elmdls.nshot,sizeof(int*));
    elmdlv.sy = (int **)calloc(elmdls.nshot,sizeof(int*));
    elmdlv.sz = (int **)calloc(elmdls.nshot,sizeof(int*));
    elmdlv.gx = (int **)calloc(elmdls.nshot,sizeof(int*));
    elmdlv.gy = (int **)calloc(elmdls.nshot,sizeof(int*));
    elmdlv.gz = (int **)calloc(elmdls.nshot,sizeof(int*));
        
    elmdlv.nsrc = (int *)calloc(elmdls.nshot,sizeof(int));
    elmdlv.nrec = (int *)calloc(elmdls.nshot,sizeof(int));
    elmdlv.tnrec =  (int *)malloc(elmdls.nshot*sizeof(int));

    if( !Get_par_int("receiver_type",&elmdls.rtype) ) elmdls.rtype = 0;

    for(i=0;i<elmdls.nshot;i++)
    {
	sprintf(Gfname,"sx%d",i+1);
	nxrec = Get_par_vint(Gfname,&(elmdlv.sx[i]));
	sprintf(Gfname,"sy%d",i+1);
	nyrec = Get_par_vint(Gfname,&(elmdlv.sy[i]));
	sprintf(Gfname,"sz%d",i+1);
	nzrec = Get_par_vint(Gfname,&(elmdlv.sz[i]));
	if( nxrec != nzrec || nxrec != nzrec || nxrec*nyrec*nzrec==0 ) 
	    sherror("%dth source array location error:sx:%d,sy:%d,sz:%d\n",i+1,nxrec,nyrec,nzrec);
	elmdlv.nsrc[i] = nxrec;

	if( elmdls.rtype == 0 )
	{
	    sprintf(Gfname,"gx%d",i+1);
	    nxrec = Get_par_vint(Gfname,&(elmdlv.gx[i]));
	    sprintf(Gfname,"gy%d",i+1);
	    nyrec = Get_par_vint(Gfname,&(elmdlv.gy[i]));
	    sprintf(Gfname,"gz%d",i+1);
	    nzrec = Get_par_vint(Gfname,&(elmdlv.gz[i]));

	    if( nxrec != nzrec || nxrec != nyrec || nxrec*nyrec*nzrec==0 ) 
		sherror("%dth receiver array location error:gx:%d,gy:%d,gz:%d\n",i+1,nxrec,nyrec,nzrec);

	    ntr = nzrec;
	}
	else
	{
	    sprintf(Gfname,"gx%d",i+1);
	    nxrec = Get_par_vint(Gfname,&tgx);
	    sprintf(Gfname,"gy%d",i+1);
	    nyrec = Get_par_vint(Gfname,&tgy);
	    sprintf(Gfname,"gz%d",i+1);
	    nzrec = Get_par_vint(Gfname,&tgz);

	    ntr = nxrec*nyrec*nzrec;

	    if( ntr <= 0 )
	    {
		sherror("%dth receiver array location error:gx:%d,gy:%d,gz:%d\n",i+1,nxrec,nyrec,nzrec);
		break;
	    }

	    elmdlv.gx[i] = (int *)calloc(ntr, sizeof(int));
	    elmdlv.gy[i] = (int *)calloc(ntr, sizeof(int));
	    elmdlv.gz[i] = (int *)calloc(ntr, sizeof(int));

	    for(iz=0,j=0; iz<nzrec; iz++)
	    {
		for(iy=0; iy<nyrec; iy++)
		{
		    for(ix=0; ix<nxrec; ix++)
		    {
			elmdlv.gx[i][j]   = tgx[ix];
			elmdlv.gy[i][j]   = tgy[iy];
			elmdlv.gz[i][j++] = tgz[iz];
		    }
		}
	    }
	    free(tgx); free(tgy); free(tgz);
	    tgx = tgy = tgz = NULL;
	}

	elmdlv.nrec[i] = ntr;
	elmdlv.tnrec[i] = ntr;
    }

    i=1;
    if( !Get_par_vchar("res_cmpnt",&tmp) ) elmdls.rcmpnt = RC_Z;
    else
    {
	i=0;
	elmdls.rcmpnt = 0;
	if( strchr(tmp,'x') ) elmdls.rcmpnt ^= RC_X,i++;
	if( strchr(tmp,'y') ) elmdls.rcmpnt ^= RC_Y,i++;
	if( strchr(tmp,'z') ) elmdls.rcmpnt ^= RC_Z,i++;
	if( strchr(tmp,'D') ) elmdls.rcmpnt ^= RC_D,i++;
	if( strchr(tmp,'C') ) 
	{
	    elmdls.rcmpnt ^= RC_CX;
	    elmdls.rcmpnt ^= RC_CY;
	    elmdls.rcmpnt ^= RC_CZ;
	    i+=3;
	}
	free(tmp);
    }
    if( (elmdls.rcmpnt & (RC_X | RC_Y | RC_D | RC_CZ | RC_CX | RC_CY | RC_CZ)) == 0 )
    {
	if( elmdls.rank[3] == 0 ) shwarn("\"res_cmpnt\" should be set by x, y, z, D, or C\nNow set to z");
	elmdls.rcmpnt ^= RC_Z;
    }
    /* 1: Velocity, 2: Displacement */
    Get_par_int("res_type",&elmdls.res_type);
    if( elmdls.res_type != 1 && elmdls.res_type != 2 )
    {
	shwarn("Wrong value: response type, Force to 1 (velocity)\n");
	elmdls.res_type = 1;
    }

    /* 1: ASCII, 2: Binary */
    Get_par_int("res_file_type",&elmdls.res_file_type);
    if( elmdls.res_file_type != 1 && elmdls.res_file_type != 2 )
    {
	shwarn("Wrong value: response file type, Force to 2 (Binary)\n");
	elmdls.res_file_type = 2;
    }
    Get_par_vchar("res_file_name",&elmdlv.res_name);

    elmdls.nrcmpnt = i;
    elmdlv.rcmpnt_list = (int *)malloc(sizeof(int)*i);
    i=0;
    if( elmdls.rcmpnt & RC_X ) elmdlv.rcmpnt_list[i] = RC_X,i++;
    if( elmdls.rcmpnt & RC_Y ) elmdlv.rcmpnt_list[i] = RC_Y,i++;
    if( elmdls.rcmpnt & RC_Z ) elmdlv.rcmpnt_list[i] = RC_Z,i++;
    if( elmdls.rcmpnt & RC_D ) elmdlv.rcmpnt_list[i] = RC_D,i++;
    if( elmdls.rcmpnt & RC_CX ) elmdlv.rcmpnt_list[i] = RC_CX,i++;
    if( elmdls.rcmpnt & RC_CY ) elmdlv.rcmpnt_list[i] = RC_CY,i++;
    if( elmdls.rcmpnt & RC_CZ ) elmdlv.rcmpnt_list[i] = RC_CZ,i++;

    Get_par_int("nt",&elmdls.nt);
    Get_par_float("dh",&elmdls.dh);
    Get_par_float("dt",&elmdls.dt);
    if( !Get_par_int("sampling",&elmdls.sample) || elmdls.sample<=0 ) elmdls.sample=1;
    if( !Get_par_int("outstep",&elmdls.outstep) || elmdls.outstep<0 ) elmdls.outstep=0;

    Get_par_int("model",&elmdls.model);
    if( elmdls.model != 1 && elmdls.model != 2 && elmdls.model != 8 && elmdls.model != 11 &&  elmdls.model != 12 )
    {
	    shwarn("Warning! Force to use a HOMOGENEOUS model, check model: %d\n",
			    elmdls.model);
	    elmdls.model = 1;
    }

    if( !Get_par_vchar("wdir",&elmdlv.wd) )
    {
	    elmdlv.wd = malloc(2*sizeof(char));
	    strcpy(elmdlv.wd,".");
    }
    else
    {
	    if( elmdlv.wd[strlen(elmdlv.wd)-1] == '/' )
	    {
		elmdlv.wd[strlen(elmdlv.wd)-1] = '\0';
	    }
    }
    /* Get Media parameter file name */
    Get_par_vchar("mfile",&elmdlv.mfile);
    if( !Get_par_vchar("oinfo",&elmdlv.oinfo) )
    {
	    elmdlv.oinfo = (char *)calloc(BUFSIZ,sizeof(char));
	    sprintf(elmdlv.oinfo,"out_info.txt");
    }

    Get_par_float("rho",&elmdls.rho);
    Get_par_float("Vp",&elmdls.Vp);
    Get_par_float("Vs",&elmdls.Vs);

    /* Get source parameters */
    Get_par_float("sfc",&elmdls.fc);
    Get_par_float("stc",&elmdls.tc);
    if( !Get_par_int("srctype",&elmdls.src_type) ) elmdls.src_type=2;
    if( !Get_par_int("sfctype",&elmdls.src_func_type) ) elmdls.src_func_type=1;

    /* Get Snapshot informations */
    if( !Get_par_int("interim_out",&elmdls.interimflag) ) elmdls.interimflag = 0;
    if( !Get_par_int("shot_start",&elmdls.shstart) || elmdls.shstart == -1 ) 
    elmdls.shstart = elmdls.tc/elmdls.dt;
    if( !Get_par_int("shot_step",&elmdls.shstep) ) elmdls.shstep = elmdls.nt/10;
    if( elmdls.shnloc = Get_par_vstr("shot_axis",&shdata) )
    {
	elmdlv.shxsec = (int *)calloc(elmdls.shnloc,sizeof(int));
	elmdlv.shgnum = (int *)calloc(elmdls.shnloc,sizeof(int));
	for(i=0;i<elmdls.shnloc;i++)
	{
	    switch(shdata[i][0])
	    {
		case 'x':
		case 'X':
		    elmdlv.shxsec[i] = 0;
		    break;
		case 'y':
		case 'Y':
		    elmdlv.shxsec[i] = 1;
		    break;
		case 'z':
		case 'Z':
		    elmdlv.shxsec[i] = 2;
		    break;
		default:
		    if( elmdls.rank[3] == 0 ) 
			    shwarn("\"shot_axis\" should be started by x, y or z\nNow set to y");
		    elmdlv.shxsec[i] = 1;
	    }
	    elmdlv.shgnum[i] = atoi(shdata[i]+1);
	    free(shdata[i]);
	}
	free(shdata);
	if( !Get_par_vchar("shot_cmpnt",&tmp) ) elmdls.shcmpnt = SN_Z;
	else
	{
	    if( strchr(tmp,'x') ) elmdls.shcmpnt ^= SN_X;
	    if( strchr(tmp,'y') ) elmdls.shcmpnt ^= SN_Y;
	    if( strchr(tmp,'z') ) elmdls.shcmpnt ^= SN_Z;
	    free(tmp);
	}
	if( (elmdls.shcmpnt & (SN_X | SN_Y | SN_Z)) == 0 )
	{
	    if( elmdls.rank[3] == 0 ) shwarn("\"shot_cmpnt\" should be set by x, y or z\nNow set to z\n");
	    elmdls.shcmpnt ^= SN_Z;
	}
    }

    itmp = NULL;
    elmdls.nOutPara = Get_par_vint("parameter_out",&itmp);
    for(i=0;i<elmdls.nOutPara;i++) 
    {
	elmdls.OutParaInfo[i] = itmp[i];
	if( elmdls.OutParaInfo[i] < 0 || elmdls.OutParaInfo[i] > 2 )
	{
	    shwarn("Wrong value: output parameter, Force to 0 (rho)\n");
	    elmdls.OutParaInfo[i] = 0;
	}
    }
    if( itmp != NULL ) free(itmp);
        
    Get_par_int("free_surface",&elmdls.freesurface);
    if( elmdls.freesurface != 0 && elmdls.freesurface != 1 )
    {
	shwarn("Wrong value: free surface type, Force to 1 (free surface)\n");
       	elmdls.freesurface = 1;
    }

    /* Get PML Boundary condition parameters */
    if( !Get_par_int("npml",&elmdls.npml) ) elmdls.npml=10;
    if( !Get_par_float("amax",&elmdls.am) ) elmdls.am=2.;
    if( !Get_par_float("wmax",&elmdls.wm) ) elmdls.wm=11.;
    if( !Get_par_float("npow",&elmdls.npow) ) elmdls.npow=2.;
    if( !Get_par_float("apow",&elmdls.apow) ) elmdls.apow=1.;

    Free_para();

    validate_input();

    return n;
}

void validate_input(void)
{
    int i,j;

    if( elmdls.nx * elmdls.ny * elmdls.nz * elmdls.nt == 0 ) 
	    sherror("Error! Check Grid sizes: nx, ny, nz, nt ==> %d, %d, %d, %d\n",
			    elmdls.nx, elmdls.ny, elmdls.nz, elmdls.nt);

    if( elmdls.dh * elmdls.dt == 0. ) 
	    sherror("Error! Check step sizes: dh, dt ==> %g, %g\n",
			    elmdls.dh, elmdls.dt);

    if( elmdls.fc * elmdls.tc == 0. ) 
	    sherror("Error! Check source time function: sfc, stc ==> %g, %g\n",
			    elmdls.fc, elmdls.tc);

    if( elmdls.src_type != 1 && elmdls.src_type != 2 && elmdls.src_type != 3 && elmdls.src_type != 4 )
    {
	    shwarn("Warning! Force srctype to 2: check value: %d\n",
			    elmdls.src_type);
	    elmdls.src_type = 2;
    }
    if( elmdls.src_func_type != 1 && elmdls.src_func_type != 2 )
    {
	    shwarn("Warning! Force sfctype to 1 (Ricker wavelet): check value: %d\n",
			    elmdls.src_func_type);
	    elmdls.src_func_type = 1;
    }

    for(i=0;i<elmdls.nshot;i++)
    {
	for(j=0; j<elmdlv.nsrc[i]; j++)
	{
	    if( !Is_inside3D( elmdlv.sx[i][j], elmdlv.sy[i][j], elmdlv.sz[i][j],
			1, 1, 1, elmdls.nx, elmdls.ny, elmdls.nz) )
	    {
		sherror("Error! Check %d th source location: %d, %d, %d\n",
			i,elmdlv.sx[i][j], elmdlv.sy[i][j], elmdlv.sz[i][j]);
	    }
	}
	for(j=0; j<elmdlv.nrec[i]; j++)
	{
	    if( !Is_inside3D( elmdlv.gx[i][j], elmdlv.gy[i][j], elmdlv.gz[i][j],
			1, 1, 1, elmdls.nx, elmdls.ny, elmdls.nz) )
	    {
		sherror("Error! Check %d th receiver location: %d, %d, %d\n",
			i,elmdlv.gx[i][j], elmdlv.gy[i][j], elmdlv.gz[i][j]);
	    }
	}
    }
	
    for(i=0;i<elmdls.shnloc;i++)
    {
	switch( elmdlv.shxsec[i] )
	{
	    case 0:
		if( elmdlv.shgnum[i] < 1 || elmdlv.shgnum[i] > elmdls.nx )
		    sherror("Error! Check %d th snap shot location: %d\n", i, elmdlv.shgnum[i]);
		break;
	    case 1:
		if( elmdlv.shgnum[i] < 1 || elmdlv.shgnum[i] > elmdls.ny )
		    sherror("Error! Check %d th snap shot location: %d\n", i, elmdlv.shgnum[i]);
		break;
	    case 2:
		if( elmdlv.shgnum[i] < 1 || elmdlv.shgnum[i] > elmdls.nz )
		    sherror("Error! Check %d th snap shot location: %d\n", i, elmdlv.shgnum[i]);
		break;
	}
    }
}
