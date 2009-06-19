#include <stdio.h>
#include <malloc.h>

int read_wf(wffile,w,datatype,foff)
char *wffile;
float **w;
char *datatype;
int foff;
{
	FILE *fp;
	int n;
	float tmp;
	long file_length;

	/* This is for a NULL file discription.  If the file name field is NULL,
	this is a signal to assume a nothing signal for that component */
	if(!strcmp(wffile,"NULL")) return(0);

	if((fp=fopen(wffile,"r"))==NULL)
	{
		fprintf(stderr,"Fatal (read_wf):  Cannot open waveform file %s\n",
			wffile);
		return(-1);
	}
	if(!strncmp(datatype,"as",2))
	{
	/*  Block for ascii data 
	first scan the input file to count number of samples to 
	be read.  This is preferable to constant reallocs */
		if(foff != 0) fprintf(stderr,
			"Warning:  foff value of %d in file %s ignored with ascii format specification\n",
			foff,wffile);
		for(n=0;;++n)
			if((fscanf(fp,"%g",&tmp))==EOF) break;
		*w = (float *) calloc(n,sizeof(float));
		rewind(fp);
		n=0;
		while(fscanf(fp,"%g",(*w+n))!=EOF) ++n;
		fclose(fp);
	}
	else if(  (!strcmp(datatype,"s2")) || (!strcmp(datatype,"s4"))
		|| (!strcmp(datatype,"t4"))  )
	{
		/* block for binary raw formats. Use crude trace read function
		under the assumption now none would use a glitched trace with
		this program.  However, because we aren't accessing this 
		trace data through a db, we don't know the number of samples
		beforehand.  To handle this, we have to look at the file
		length and calculate the number of samples */
		fseek(fp,0L,2);
		file_length = ftell(fp);
		if(!(strcmp(datatype,"s2")))
			n = (int) ((file_length-foff)/2);
		else
		{
			n = (int) ((file_length-foff)/4);
		/* we have to close fp because read_trace does an open
		internally */
			fclose(fp);
			*w = (float *) read_trace(wffile, foff, datatype, &n);
		}
	}
	else
	{
		fprintf(stderr,"Error (read_wf):  Unknown data type %s\n",
			datatype);
		fprintf(stderr,"No data read for waveform file %s\n",wffile);
		return (-2);
	}	 
	return(n);
}

 
int
zaccess (path, mode)
 
char *path;
int mode;
 
{
        char fname[1024];
 
        if (access(path, mode) == 0) return (0);
        sprintf (fname, "%s.Z", path);
        if (access(fname, mode) == 0) return (1);
        sprintf (fname, "%s.gz", path);
        if (access(fname, mode) == 0) return (2);
        return (-1);
}

