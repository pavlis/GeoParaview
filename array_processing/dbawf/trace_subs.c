#include <stdio.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <strings.h>

#include "db.h"

/* Blindly reads trace files for dbawf program.  Modified from original
Danny Harvey code in dbbeams.  Program is simple because it does
no checking for data gaps or other real data flaws.  This is because
this program is designed primarily for simulations, although one can
imagine using it for subracting beam traces.  In that case, fix the
gaps first. */

float *read_trace (fname, foff, datatype, nsamps)

char *fname;
int foff;
char *datatype;
int *nsamps;

{
	FILE *file;
	int n, size;
	unsigned char *buf;
	float *t;
	int i;
	int *lt;
	short *st;

	size = atoi(&datatype[strlen(datatype)-1]);
	file = fopen (fname, "r");
	if (file == NULL) {
		fprintf (stderr, "read_file: Unable to open '%s'\n",
							fname);
		return (NULL);
	}
	if (fseek(file, foff, 0) < 0) {
		fprintf (stderr, "read_file: fseek() error on '%s'\n",
							fname);
		fclose (file);
		return (NULL);
	}
	buf = (unsigned char *) malloc (size*(*nsamps));
	if (buf == NULL) {
		fprintf (stderr, "read_file: Malloc error.\n");
		fclose (file);
		return (NULL);
	}
	if ((n=fread (buf, size, *nsamps, file)) < *nsamps) {
		if (n < 1) {
			fprintf (stderr, "read_file: fread() error on '%s'\n",
							fname);
			fclose (file);
			free (buf);
			return (NULL);
		} else {
			fprintf (stderr, "read_file: Read %d samples but expected %d on '%s'\n",
							n, *nsamps, fname);
			*nsamps = n;
		}
	}
	fclose (file);
	/* convert integer formats to float*/
	if(!strncmp(datatype,"t4",2))
		t=(float *)buf;
	else if(!strncmp(datatype,"s4",2))
	{
		t = (float *)calloc(*nsamps,sizeof(float));
		for(i=0;i<*nsamps;++i) 
		{
			lt=(int *) (buf+i*4);
			t[i] = (float )*lt;
		}
		free(buf);
	}
	else if(!strncmp(datatype,"s2",2))
	{
		t = (float *)calloc(*nsamps,sizeof(float));
		for(i=0;i<*nsamps;++i)
		{
			st=(short *) (buf+i*2);
			t[i] = (float)*st;
		}
		free(buf);

	}
	else
	{
		fprintf(stderr,"read_file:  Unknown sample code %s for file %s\n",
			datatype,fname);
		fprintf(stderr,"trace skipped");
		return(NULL);
	}
	return (t);
}
/* Inverse of read_trace.  Note anything before foff is left intact. */
write_trace (fname, foff, datatype, nsamps,t)

char *fname;
int foff;
char *datatype;
int nsamps;  /* Note nsamps is not a pointer here like it was in read */
float *t;

{
	FILE *file;
	int n, size;
	int i;
	int *lt;
	short *st;

	file = fopen (fname, "r+");
	if (file == NULL) {
		fprintf (stderr, "write_file: Unable to open '%s'\n",
							fname);
		return (1);
	}
	if (fseek(file, foff, 0) < 0) {
		fprintf (stderr, "write_file: fseek() error on '%s'\n",
							fname);
		fclose (file);
		return (1);
	}
	/* convert back to original data representation*/
	if(!strncmp(datatype,"t4",2))
	{
		if((n=fwrite(t,4,nsamps,file)) < nsamps)
			fprintf(stderr,"write_file: fwrite() error on '%s'\n",
							fname);
		free(t);
	}	
	else if(!strncmp(datatype,"s4",2))
	{
		lt = (int *)calloc(nsamps,4);
		for(i=0;i<nsamps;++i) 
			lt[i] = (int) t[i];
		if((n=fwrite(lt,4,nsamps,file)) < nsamps)
			fprintf(stderr,"write_file: fwrite() error on '%s'\n",
							fname);
		free(t);
		free(lt);
	}
	else if(!strncmp(datatype,"s2",2))
	{
		st = (short *)calloc(nsamps,2);
		for(i=0;i<nsamps;++i) 
			st[i] = (short) t[i];
		if((n=fwrite(st,2,nsamps,file)) < nsamps)
			fprintf(stderr,"write_file: fwrite() error on '%s'\n",
							fname);
		free(t);
		free(st);
	}
	fclose (file);
}


/* Simple routine to scale trace t by factor "scale"*/

scale_trace(n,t,scale)
int n;
float *t;
float scale;
{
	int i;
	for(i=0;i<n;++i) t[i] *= scale;
}

/* The following routine adds a waveform in array "wf" to the signals
stored in "trace".  Arguments:
	nt - number of elements in trace array
	t - float array of signal into which wf is to be added.
	nwf - number of elements in wf array.
	weight - scale factor.  Each element of wf is multiplied by weight
		before being accumulated.  
	wf - waveform to sum in 
	lag - number of samples from start of trace to position wf.  i.e.
		the basic algorithm is that trace[lag+i] += weight*wf[i]
		Note that negative lag values are allowed, but wf will
		be trimmed at the front end this way.  Similarly, for
		large lag values either nothing happens or wf will
		be trimmed on the right.

The routine returns the number of elements of wf actually summed into the
trace array.  
*/
add_wf_shifted(nt,t,nwf,weight,wf,lag)
int nt;
float	*t;
int		nwf;
float			weight;
float				*wf;
int					lag;
{
	int i,j;
	int count;
	if(nt<nwf) nwf = nt;
	/* Large negative lag case */
	if((lag+nwf)<0) return(0);
	/* Negative lag case with trimming */	
	if(lag<0)
	{
		for(i=(-lag),j=0;i<nwf;++i,++j) 
			t[j] += weight*wf[i];
		return (j);
	}
	else
	{
		for(i=0,j=lag;(i<nwf) && (j<nt);++i,++j)
			t[j] += weight*wf[i];
		return(i);
	}
}
