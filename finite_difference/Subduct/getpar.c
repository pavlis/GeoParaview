/*###################################################################
	Parallel poroelastice wave propagation in 3D Version 1.1
	Last modified: April, 2005
  
	Author:	Dong-Hoon Sheen
		Ph.D. student at Seoul National University
	       	Research associate at Indiana University
	
	Reference: 
	Parallel implementation of a velocity-stress
	    staggered-grid finite-difference method 
	    for poroelastic wave propagation
	    Sheen, D.-H., Tuncay, K., Baag, C.-E. and Ortoleva, P. J.
###################################################################*/
/*
	Changes to v1.1
	Tail_cut: 
	    ADD 'if': if( i > 0 && i < strlen(str) ) 
*/

#include "getpar.h"

struct para{
	int Nvalue;
	char *Cname;
	char **Cvalue;
	struct para *next;
};
typedef struct para PARA;
PARA *tpara,*bpara;
int Npara;

int Is_digit(char str)
{
	if( str >= '0' && str <= '9' )	return TRUE;
	return 0;
}

int Is_char(char str)
{
	if((str>='a'&&str<='z')|| (str>='A'&&str<='Z'))	return TRUE;
	return 0;
}
int Is_blank(char str)
{
	if( str == ' ' || str == '\t' )	return TRUE;
	return 0;
}

int Where_char(char *str,char ch)
{
	register int i,n;
	n = strlen(str);
	for(i=0;i<n;i++) if( str[i] == ch ) return i;
	return -1;
}

void Tail_cut(char *str)
{
	int i,n;
	i = Where_char(str,'#');
	if( i > 0 && i < strlen(str) ) str[i] = '\0';
	n = strlen(str);
	for(i=n-1;i>=0;i--)
	{
		if( str[i] != ' ' && str[i] != '\t' && str[i] != '\n' )
			break;
	}
	if( i >= 0 && i < n-1 ) str[i+1] = '\0';
}
int Skip_blank(char *str)
{
	int i=0;
	while( str[i] != '\n' && str[i] != '\0' )
	{
		if( Is_blank(str[i]) ) i++;
		else if ( str[i] == '\n' || str[i] == '\0' ) return WARN;
		else return i;
	}
	return WARN;
}

int Save_para(char *str)
{
	register int i,j,n,nvalue,pindex;
	char *cname,**cvalue;

	n = strlen(str);

	/* check and save parameter's name */
	if( ((pindex=Where_char(str,'='))<=0) 
			|| (n<=pindex+1) ) return -1;
	str[pindex] = '\0';
	Tail_cut(str);
	i = strlen(str);
	cname = (char *)calloc(i+1,sizeof(char));
	strcpy(cname,str);
	/* null terminate string */

	/* move string position to check paramter's value */
	str += pindex+1;
	n -= pindex+1;
	nvalue = 1;
	for(i=0;i<n;i++) if(str[i]==',') nvalue++;
	cvalue =(char **)calloc(nvalue,sizeof(char *));
	for(i=0;i<nvalue-1;i++)
	{
		j = Skip_blank(str);
		if( j < 0 ) return -1;
		str += j;
		if( (pindex=Where_char(str,',')) <= 0 ) return -1;
		str[pindex] = '\0';
		Tail_cut(str);
		j = strlen(str);
		cvalue[i] = (char *)calloc(j+1,sizeof(char));
		strcpy(cvalue[i],str);
		str += pindex+1;
	}
	j = Skip_blank(str);
	if( j < 0 ) return -1;
	str += j;
	j = strlen(str);
	cvalue[i] = (char *)calloc(j+1,sizeof(char));
	strcpy(cvalue[i],str);

	/* top,bottom of parameter stack */
	if( tpara == NULL )
	{
		tpara = (PARA *)calloc(1,sizeof(PARA));
		bpara = tpara;
		tpara->Nvalue = nvalue;
		tpara->Cname  = cname;
		tpara->Cvalue = cvalue;
		tpara->next   = NULL;
	}
	else
	{
		bpara->next = (PARA *)calloc(1,sizeof(PARA));
		bpara = bpara->next;
		bpara->Nvalue = nvalue;
		bpara->Cname  = cname;
		bpara->Cvalue = cvalue;
		bpara->next   = NULL;
	}
	Npara++;
	return nvalue;
}

int Read_para_file(char *sname)
{
	int ipara,npara,t,iline;
	char *buf;
	FILE *fp;
	
	fp=echkfopen(sname,"r",sherror);
	if( fp == NULL ) return 0;
	
	tpara = NULL;
	bpara = NULL;

	npara = 0;
	Npara = 0;
	buf = (char *)calloc(BUFSIZ,sizeof(char));
	while( fgets(buf,BUFSIZ,fp) != NULL )
	{	
		if( (t = Skip_blank(buf)) != WARN )
			if( buf[t] != '#' ) npara++;
	}
	if( npara == 0 ) 
		sherror("File %s open error, check and try again\n",sname);

	ipara=0;
	iline=1;

	fseek(fp,0L,SEEK_SET);
	while( fgets(buf,BUFSIZ,fp) != NULL )
	{
		if( strlen(buf) > 1 )
		{
		   if( (t = Skip_blank(buf)) != WARN )
		   {
			if( buf[t] != '#' )
			{
		   		Tail_cut(buf);
				if( Save_para(buf+t) < 0 )
					printf("(%d)->%s\n",iline,buf);
				else ipara++;
			}
		    }
		}
		iline++;
	}

	fclose(fp);
	free(buf);
	return npara;
}

int Get_par_nchar(char *str,char *val,int pos)
{
	PARA *tmp;
	tmp = tpara;
	while(tmp != NULL)
	{
		if( strcmp(str,tmp->Cname) == 0 )
		{
			*val = tmp->Cvalue[pos][0];
			return 1;
		}
		tmp = tmp->next;
	}
	return 0;
}
int Get_par_char(char *str,char *val)
{
	return Get_par_nchar(str,val,0);
}
int Get_par_vchar(char *str,char **val)
{
	return Get_par_nstr(str,val,0);
}

int Get_par_nint(char *str,int *val,int pos)
{
	int n;
	int *tval;

	*val = 0;
	if( pos < 0 ) return 0;
	n = Get_par_vint(str,&tval);
	if( n > 0 ) 
	{
	    if( pos < n ) *val = tval[pos];
	    free(tval);
	    return 1;
	}
	return 0;
}
int Get_par_int(char *str,int *val)
{
	return Get_par_nint(str,val,0);
}
int Get_par_vint(char *str,int **val)
{
	int i,n,ilist,nlist,add;
	int start,end,step,tv;
	PARA *tmp;
	tmp = tpara;
	while(tmp != NULL)
	{
	    if( strcmp(str,tmp->Cname) == 0 )
	    {
		n = nlist = tmp->Nvalue;
		if( n <= 0 ) return 0;
		*val = (int *)calloc(nlist,sizeof(int));
		for(i=0,ilist=0;i<n;i++)
		{
		    if( tmp->Cvalue[i][0] == ':' )
		    {
			if( i==0 ) start=0;
			else if( i < n-1 ) start = (*val)[ilist-1];
			else			break;
		        sscanf(tmp->Cvalue[i]+1,"%d",&step);
		        sscanf(tmp->Cvalue[i+1],"%d",&end);
			if( step == 0 )		break;
			else if( step > 0 ) add = (end-start-1)/step;
			else 		    add = (end-start+1)/step;
			nlist += (add-1);
			if( add < 0 )		break;
			else if( add != 1 )
			    *val = (int *)realloc( *val, nlist*sizeof(int));
			if( step > 0 ) 
			    for(tv=start+step;tv<end;tv+=step) (*val)[ilist++] = tv;
			else
			    for(tv=start+step;tv>end;tv+=step) (*val)[ilist++] = tv;
		    }
		    else if( tmp->Cvalue[i][0] == ';' )
		    {
			if( i==0 ) start=0;
			else start = (*val)[ilist-1];
		        sscanf(tmp->Cvalue[i]+1,"%d",&add);
			if( add <= 0 ) break;
			nlist += (add-2);
			if( add != 2 )
			    *val = (int *)realloc( *val, nlist*sizeof(int));
			for(step=1;step<add;step++)
			    (*val)[ilist++] = start;
		    }
		    else 
			 sscanf(tmp->Cvalue[i],"%d",(*val)+(ilist++));
		}
		if( i != n )
		{
		    free( (*val) );
		    *val = NULL;
		    return 0;
		}
		return nlist;
	    }
	    tmp = tmp->next;
	}
	return 0;
}

int Get_par_nfloat(char *str,float *val,int pos)
{
	int n;
	float *tval;

	*val = 0;
	if( pos < 0 ) return 0;
	n = Get_par_vfloat(str,&tval);
	if( n > 0 ) 
	{
	    if( pos < n ) *val = tval[pos];
	    free(tval);
	    return 1;
	}
	return 0;
}
int Get_par_float(char *str,float *val)
{
	return Get_par_nfloat(str,val,0);
}
int Get_par_vfloat(char *str,float **val)
{
	int i,n,ilist,nlist,add;
	float start,end,step,tv;
	PARA *tmp;
	tmp = tpara;
	while(tmp != NULL)
	{
	    if( strcmp(str,tmp->Cname) == 0 )
	    {
		n = nlist = tmp->Nvalue;
		if( n <= 0 ) return 0;
		*val = (float *)calloc(nlist,sizeof(float));
		for(i=0,ilist=0;i<n;i++)
		{
		    if( tmp->Cvalue[i][0] == ':' )
		    {
			if( i==0 ) start=0.;
			else if( i < n-1 ) start = (*val)[ilist-1];
			else			break;
		        sscanf(tmp->Cvalue[i]+1,"%f",&step);
		        sscanf(tmp->Cvalue[i+1],"%f",&end);
			if( step == 0 )		break;
			else if( step > 0 ) add = (int)floor((end-start-step*1.e-5)/step);
			else 		    add = (int)floor((end-start+step*1.e-5)/step);
			nlist += (add-1);
			if( add < 0 )		break;
			else if( add != 1 )
			    *val = (float *)realloc( *val, nlist*sizeof(float));
			if( step > 0 ) 
			    for(tv=start+step;tv<end;tv+=step) (*val)[ilist++] = tv;
			else
			    for(tv=start+step;tv>end;tv+=step) (*val)[ilist++] = tv;
		    }
		    else if( tmp->Cvalue[i][0] == ';' )
		    {
			if( i==0 ) start=0.;
			else start = (*val)[ilist-1];
		        sscanf(tmp->Cvalue[i]+1,"%d",&add);
			if( add <= 0 ) break;
			nlist += (add-2);
			if( add != 2 )
			    *val = (float *)realloc( *val, nlist*sizeof(float));
			for(step=1;step<add;step++)
			    (*val)[ilist++] = start;
		    }
		    else 
			 sscanf(tmp->Cvalue[i],"%f",(*val)+(ilist++));
		}
		if( i != n )
		{
		    free( (*val) );
		    *val = NULL;
		    return 0;
		}
		return nlist;
	    }
	    tmp = tmp->next;
	}
	return 0;
}

int Get_par_ndouble(char *str,double *val,int pos)
{
	int n;
	double *tval;

	*val = 0;
	if( pos < 0 ) return 0;
	n = Get_par_vdouble(str,&tval);
	if( n > 0 ) 
	{
	    if( pos < n ) *val = tval[pos];
	    free(tval);
	    return 1;
	}
	return 0;
}
int Get_par_double(char *str,double *val)
{
	return Get_par_ndouble(str,val,0);
}
int Get_par_vdouble(char *str,double **val)
{
	int i,n,ilist,nlist,add;
	double start,end,step,tv;
	PARA *tmp;
	tmp = tpara;
	while(tmp != NULL)
	{
	    if( strcmp(str,tmp->Cname) == 0 )
	    {
		n = nlist = tmp->Nvalue;
		if( n <= 0 ) return 0;
		*val = (double *)calloc(nlist,sizeof(double));
		for(i=0,ilist=0;i<n;i++)
		{
		    if( tmp->Cvalue[i][0] == ':' )
		    {
			if( i==0 ) start=0.;
			else if( i < n-1 ) start = (*val)[ilist-1];
			else			break;
		        sscanf(tmp->Cvalue[i]+1,"%lf",&step);
		        sscanf(tmp->Cvalue[i+1],"%lf",&end);
			if( step == 0 )		break;
			else if( step > 0 ) add = (int)floor((end-start-step*1.e-5)/step);
			else 		    add = (int)floor((end-start+step*1.e-5)/step);
			nlist += (add-1);
			if( add < 0 )		break;
			else if( add != 1 )
			    *val = (double *)realloc( *val, nlist*sizeof(double));
			if( step > 0 ) 
			    for(tv=start+step;tv<end;tv+=step) (*val)[ilist++] = tv;
			else
			    for(tv=start+step;tv>end;tv+=step) (*val)[ilist++] = tv;
		    }
		    else if( tmp->Cvalue[i][0] == ';' )
		    {
			if( i==0 ) start=0.;
			else start = (*val)[ilist-1];
		        sscanf(tmp->Cvalue[i]+1,"%d",&add);
			if( add <= 0 ) break;
			nlist += (add-2);
			if( add != 2 )
			    *val = (double *)realloc( *val, nlist*sizeof(double));
			for(step=1;step<add;step++)
			    (*val)[ilist++] = start;
		    }
		    else 
			 sscanf(tmp->Cvalue[i],"%lf",(*val)+(ilist++));
		}
		if( i != n )
		{
		    free( (*val) );
		    *val = NULL;
		    return 0;
		}
		return nlist;
	    }
		tmp = tmp->next;
	}
	return 0;
}

int Get_par_nstr(char *str,char **val,int pos)
{
	PARA *tmp;
	tmp = tpara;
	while(tmp != NULL)
	{
		if( strcmp(str,tmp->Cname) == 0 )
		{
			if( pos > tmp->Nvalue ) return 0;
			*val = (char *)calloc(strlen(tmp->Cvalue[pos])+1,sizeof(char));
			strcpy(*val,tmp->Cvalue[pos]);
			return 1;
		}
		tmp = tmp->next;
	}
	return 0;
}
int Get_par_str(char *str,char **val)
{
	return Get_par_nstr(str,val,0);
}
int Get_par_vstr(char *str,char ***val)
{
	int i;
	PARA *tmp;
	tmp = tpara;
	while(tmp != NULL)
	{
		if( strcmp(str,tmp->Cname) == 0 )
		{
			*val = (char **)calloc(tmp->Nvalue,sizeof(char *));
			for(i=0;i<tmp->Nvalue;i++)
			{
				(*val)[i] = (char *)calloc(strlen(tmp->Cvalue[i])+1,sizeof(char));
				strcpy((*val)[i],tmp->Cvalue[i]);
			}
			return tmp->Nvalue;
		}
		tmp = tmp->next;
	}
	return 0;
}

void Check_para(void)
{
	register int i;
	PARA *tmp;
	tmp = tpara;
	while(tmp != NULL)
	{
		printf("========================\n");
		printf("%s\n",tmp->Cname);
		for(i=0;i<tmp->Nvalue;i++)
			printf("\t%s\n",tmp->Cvalue[i]);
		tmp = tmp->next;
	}
}

void Free_para(void)
{
	register int i;
	PARA *tmp,*ptmp;
	tmp = tpara;
	while(tmp != NULL)
	{
		free(tmp->Cname);
		for(i=0;i<tmp->Nvalue;i++)
			free(tmp->Cvalue[i]);
		free(tmp->Cvalue);
		ptmp = tmp;
		tmp = tmp->next;
		free(ptmp);
		ptmp = NULL;
	}
}

#ifdef GET_PAR_TEST

int main(int argc, char *argv[])
{
	int n,itest,*vitest;
	float ftest,*vftest;
	double dtest,etest,*vdtest;
	char *str,**vstest,*vctest;

	Read_para_file("input.txt");
	Check_para();
	n = Get_par_vfloat("rho",&vftest);
	Get_par_ndouble("int_test",&dtest,2);

	Free_para();
	printf("int_test:%g\n",dtest);
	printf("n:%d\n",n);
	for(itest=0;itest<n;itest++)
		printf("%3d:%g\n",itest,vftest[itest]);

	return 0;
}

#endif
