/* This library is designed to parse a pf stream for use mainly in data-driven, parallel
processing algorithms.  A pf stream is a stream of parameter file type blocks of 
metadata separated by a keyword defined in pfstream.h.  That is, each block of
ascii characters between the keyword sentinel is assumed to consist of keywords and
numbers as used in antelope parameter files.  

The format has two modes.  A single object mode will call pfcompile on each block of
text bounded by the input separator keyword.  In "ensemble" mode, the block of text 
is assumed to be in this structure:
ensemble &Arr{
---- constant parameters for this ensemble----
000001 &Arr {
--parameter for object 1----
}
000002 &Arr {
--parameters for object 2---
}
repeat for n objects

}

This file contains basic input and output functions for a pfstream.

Author:  GAry L. Pavlis
Written:  September 2002
*/
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include "stock.h"
#include "coords.h"
#include "pf.h"
#include "elog.h"
#include "pfstream.h"

/* This is a low level read routine that sucks in data from file fname 
until the  end of block line (symbol ENDPF_SENTINEL defined in pfstream.h)
is read.  It then closes the file to allow another process to access the
file under an assuption that fname is a fifo.  Finally, 
pfcompile is then called and the resulting pfobject is returned.

The algorithm is more complicated than one might think it needs to be because
the size of the input block of text is not predictable.  With the structure
of parameter file descriptions it is necessary to suck in the whole file
before trying to interpret it.  As a result the algorithm uses an initial
buffer size which grows as needed whenever the previous read gets close to 
a high water mark.  There may be a more elegant way to do this with some
alternative function, but I don't know about it.  

Another major complication is that an empirical study by the author found
that at least in Solaris the read has to be done with the unbuffered read
function.  This is necessary because stdio routines, which are buffered,
will extract data from a fifo beyond the ENDPF_SENTINEL causing dropped
data at the beginning of each block.  

This function depends upon two documented features in solaris unix that
could prove problematic in portability:
1. the man page says realloc will not alter the contents of memory allocated
and set previously.  I depend on this to not have to copy the previous buffer.
2.  The whole idea of working with the fifo here assumes something not clear
in the man pages, but which I found worked empirically.  That is once a 
process opens a fifo for read it gets all that comes through it until 
it closes the pipe.  Another process can successfully open the file but 
will block until the previous reader calls close.  

RETURNS:
normal return is a valid pointer to a pf.  A NULL pointer indicates 
a failure in pfcompile.  The function will die if memory allocation fails.  

Author:  Gary L. Pavlis
Written:  September 2002
*/
#define BUFSIZE0 10000
#define DELTA_BUF 5000
#define MAXLINE 512   /* realloc when this close to end of buffer = maximum line length*/
/* first this small unbuffered companion similar to fgets.  Reads up to maxline
characters or a new line into inline from open file with file descriptor fd. 

Returns:  
	normal = number of characters in returned line
	0 or -1 read errors.  0 may not be an error as it can mean the 
		writer closed the fifo.  
	maxline - input line did not encounter a newline

*/
int getline(char *line,int maxline, int fd)
{
	int i;
	int ierr;
	for(i=0;i<maxline;++i)
	{
		ierr=read(fd,line+i,1);
		if(ierr<=0) return(ierr);
		if(line[i]=='\n') 
		{
			line[i+1]='\0';
			return(i);
		}
	}
	return(maxline);
}
Pf *pfstream_read(char *fname)
{
	char *buffer;
	char line[MAXLINE];
	int current_buffer_size=BUFSIZE0;
	int fd;
	int space_in_buffer=BUFSIZE0;
	int high_water_mark=0;
	Pf *pf=NULL;
	int linecount=0,ncread,ierr;

	buffer=(char *)malloc(BUFSIZE0);
	if(buffer==NULL) elog_die(0,"pfstream_read:  cannot malloc input buffer\n");
	buffer[0]='\0';

	fd = open(fname,O_RDONLY | O_EXCL);

	while(ncread=getline(line,MAXLINE,fd)>0)
	{
		++linecount;
		if(ncread==MAXLINE)
			elog_complain(0,"pfstream_read encountered a long line of length %d that exceeded the internal buffer size of %d\nData may be dropped\n",
				ncread,MAXLINE);
		if(strstr(line,ENDPF_SENTINEL)!=NULL )break;

		strcat(buffer,line);
		high_water_mark=strlen(buffer);
		if(high_water_mark>(current_buffer_size-MAXLINE))
		{
			current_buffer_size+=DELTA_BUF;
			buffer=(char *)realloc(buffer,current_buffer_size);
			if(buffer==NULL)
			  elog_die(0,"pfstream_read:  realloc of %d byte buffer failed\n",
							current_buffer_size);
		}
	}
	close(fd);
	if(linecount<=0) return(NULL);
	ierr=pfcompile(buffer,&pf);
	if(ierr!=0) 
	{
		elog_complain(1,"pfstream_read:  pfcompile failed\n");
		pf=NULL;
	}
	free(buffer);
	return(pf);
}
/*  The next block of stuff are Pf_ensemble routines */


/* This is a C constructor for a Pf_ensemble object with hanging 
pointers.  It is called by the ensemble input routine below */

Pf_ensemble *create_Pf_ensemble(int nmembers,int ngroups)
{
	Pf_ensemble *pfe;
	int i;

	allot(Pf_ensemble *,pfe,1);
	/* I intentionally set these null to start */
	pfe->ensemble_keys=NULL;
	pfe->group_keys=NULL;
	if(nmembers>0)
	{
		allot(Pf **,pfe->pf,nmembers);
		pfe->nmembers=nmembers;
		/* explicit initialization a good idea here */
		for(i=0;i<nmembers;++i) pfe->pf[i]=NULL;
	}
	else
		pfe->nmembers=0;
	if(ngroups>0)
	{
		allot(int *,pfe->group_start,ngroups);
		allot(int *,pfe->group_end,ngroups);
		pfe->ngroups=ngroups;
	}
	else
	{
		pfe->group_start=NULL;
		pfe->group_end=NULL;
		pfe->ngroups=0;
	}
	return(pfe);
}
/* destructor */
void free_Pf_ensemble(Pf_ensemble *pfe)
{
	int i;
	if(pfe->nmembers > 0)
	{
		for(i=0;i<(pfe->nmembers);++i)
			if(pfe->pf[i] != NULL) pffree(pfe->pf[i]);
		free(pfe->pf);
	}
	if(pfe->ngroups>0)
	{
		if(pfe->group_start != NULL) free(pfe->group_start);
		if(pfe->group_end != NULL) free(pfe->group_end);
	}
	/* This is dangerous because reading these from a pf file
	normally does not duplicate the tbl entries read from
	the parent file.  This assumes these tbl's were created
	by a copy operation like that below */
	if(pfe->ensemble_keys!=NULL) freetbl(pfe->ensemble_keys,free);
	if(pfe->group_keys!=NULL) freetbl(pfe->group_keys,free);
	free(pfe);
}

		
/* This routine takes apart an ensemble pf description and returns a vector of
pf objects for the components of the ensemble.  (see above).  It should always
be called after a read function (pfstream_read above, pfread, or extracting pf's
from an orb) to take the parameter ensemble apart.  

The algorithm has some built in assumptions that the pieces of the ensemble
are tagged by 6 digit numbers that start at 000000.  6 digits are surely 
overkill, but given the way memory is going it may not be so silly.  The are
taken apart in the order Arr returns them in to build Tbl containing pointers
to pf objects of the components of the ensemble.  That is, the ith pf object
from this ensemble can be extracted from the output tbl, t, with a call like this:
pf=(Pf *)gettbl(t,i)

This structure, although complex, was done by design as it allows a complete
ensemble to be encapsulated in a single pf object.  This simplifies and 
abstracts the process of reading from a fifo or from an orb.  i.e. calls
two functions that do quite different things can yield the same result.
This could be overloaded in C++ if I wanted to go that route. 

Note that the input Pf object, pf, is only accessed to take apart the 
components of the "ensemble &Arr{...}" component of Pf.  If there are other parameters
in pf outside the range of this object the caller needs to handle these
separately.

Normal returns is a pointer to the Pf_ensemble object.  A NULL pointer
is returned if cracking the complex pf defined by the ensemble failed.

Author:  Gary L. Pavlis
Date:  September 2002
Modified:  October 2002
Rather than rewrite teh comments above note that a design change was
recognized that I needed a variable tag to name the ensemble.  Above
and at the top of this file I refer to the keyword "ensemble".  
There are cases, however, when multiple tags are needed.  The type
case is in converting a pfstream to a database rows.  If multiple
tables are needed to produce the required output at some point
in the running of a program this allows a way to dump tables with
different numbers of rows that may not even be related. 
So, primary thing is tag is a variable that needs to be passed downstream.
*/
Pf_ensemble *pfget_Pf_ensemble(Pf *pfin,char *tag)
{
	Pf_ensemble *pfe;
	Tbl *ttmp;
	int nmembers,ngroups;

	Tbl *t=newtbl(0);
	Pf *pferaw,*pf;
	Tbl *list_keys;
	int i;

	/* We duplicate the elements of this tbl so 
	we can release the space of the parent strings */
	ttmp=pfget_tbl(pfin,"group_keys");
	if(ttmp!=NULL)
	{
		pfe->group_keys=duptbl(ttmp,strdup);
		freetbl(ttmp,0);
	}
	/* group_records is a Tbl containing secondary grouping of 
	the ensemble.  If it is empty, we assume no grouping */
	ttmp = pfget_tbl(pfin,"group_records");
	if(ttmp==NULL)
		ngroups=0;
	else
		ngroups=maxtbl(ttmp);

	pferaw=NULL;

	/*This extracts the data enclosed by "ensemble &Arr {" to "}" */
	if(pfget(pfin,"ensemble",(void **)&pferaw) != PFARR) return(NULL);

	/* this loop takes apart the pieces of the ensemble tagged by unspecified keywords.*/
	/* output will be in the order determined by the tags (normally numbers like 000000, 000001, etc.*/
	list_keys=pfkeys(pferaw);
	/* We now know the size of the Pf_ensemble and can create one */
	nmembers=maxtbl(list_keys);
	pfe=create_Pf_ensemble(nmembers,ngroups);
	/* First we build the group array vectors */
	if(ngroups>0)
	{
		for(i=0;i<maxtbl(ttmp);++i)
		{
			int is,ie;
			char *line;
			line=gettbl(ttmp,i);
			sscanf(line,"%d%d",&is,&ie);
			pfe->group_start[i]=is;
			pfe->group_end[i]=ie;
		}
		freetbl(ttmp,0);
	}
	/* now the big work -- filling the Pf vector */
	for(i=0;i<maxtbl(list_keys);++i)
	{
		char *key;

		pf=NULL;
		key = gettbl(list_keys,i);
		if(pfget(pferaw,key,(void **)&pf) != PFARR)
			elog_complain(0,"syntax error parsing ensemble for key %s\nData for this member of this ensemble will be ignored\n", key);
		else
		{
			/* We want to duplicate this pf so we can 
			release the input string storage */
			pfe->pf[i]=pfdup(pf);
		}
	}
	freetbl(list_keys,0);
	/* This is going to leave some memory leaks of the pf's created
	during parsing the Arr heirarchy, but this seems unavoidable. */
	pffree(pfin);
	return(pfe);
}

/* Now we have the inverse functions.  The above were readers while the following 
are all writers */

/* First a function to open a fifo in a way that seems to guarantee exclusive
access.  fname is the name of the file/fifo to be openned for writing.  
(note it should work equivalently for ordinary files) 

Normal return is a valid file pointer for stream i/o.  A NULL pointer means
the file couldn't be openned.
*/

FILE *open_pfstream_output(char *fname)
{
	int fd;
	FILE *fp;

	fd = open(fname,O_WRONLY | O_EXCL);  /* O_EXCL seems to supply an implicit lock */
	if(fd<0)
	{
		elog_complain(0,"Cannot open %s\n",fname); 
		return(NULL);
	}
	/* necessary to get stream fp from file descriptor */
	fp=fdopen(fd,"w");
	if(fp==NULL) 
		return(NULL);
	else
		return(fp);
}
	

/* This function dumps pf to and output file associated with the stream fp.  
There are two ways this function might be called.  With a fifo I can be
useful to close fp after each call to this function as a way to block the
output into working chunks.  This is controlled by the 
boolean close_on_return variable.  If true (nonzero) fclose will be 
called before returning.  Otherwise the ENDPF_SENTINEL is written to 
output and fp is left open.   

If fp is NULL on entry the fp function above will be called on the file/stream
described logically by fname.  If fp is not NULL, fname is ignored.  
 
Normal completion returns 0.  -1 for errors that are posted on the 
elog.

Author:  Gary Pavlis
Written:  September 2002
*/

int pfwrite_stream(Pf *pf, char *fname, FILE *fp, int close_on_return)
{
	if(fp==NULL)
	{
		fp=open_pfstream_output(fname);
		/* assume open_pfstream_output posted an error in ths condition*/
		if(fp==NULL) return(-1);
	}
	pfout(fp,pf);
	if(close_on_return) fclose(fp);
	return(0);
}

/* this is the reciprocal of pfget_Pf_ensemble.  Basically it compiles
the pieces of a Pf_ensemble into a single pf object with a keyword
that identifies the block defined by the input variable tag.  

This is a large memory algorithm that works in a brute force way.  pf2string is
used to expand each of the components, the components are combined by string
manipulation, and the results are converted by to a pf by pfcompile.  

Arguments:
	ensemble_pf - (optional) set of parameters to dump at head
		of output stream that is build here.  That is these
		are global parameters.  This should be set to NULL
		before calling this function if no global parameters
		are required.
	pfe - Pf_ensemble that is to be compiled into single output pf
	tag - name used to tag the output Arr.

Author:  Gary L. Pavlis
Date:  September 2002
*/

Pf *build_ensemble(Pf *ensemble_pf,Pf_ensemble *pfe, char *tag)
{
	char **blocks;
	char *pfimage;
	int sum_block_sizes,pfimage_size;
	char *hdrblock;
	int hdrsize;
	int i;
	Pf *pfresult=NULL;
	int npf;  /* useful shorthand */

	npf=pfe->nmembers;
	allot(char **,blocks,npf);
	for(i=0,sum_block_sizes=0;i<npf;++i)
	{
		blocks[i]=pf2string(pfe->pf[i]);
		sum_block_sizes+=strlen(blocks[i]);
	}
	if(ensemble_pf==NULL)
		hdrblock=strdup("");
	else
		hdrblock=pf2string(ensemble_pf);
	hdrsize = strlen(hdrblock);
	/* this is a generous estimate of the size of the full pf image*/
	pfimage_size = sum_block_sizes + npf*20+hdrsize+100;
	pfimage = (char *) malloc(pfimage_size);
	if(pfimage==NULL) elog_die(0,"Cannot malloc %d bytes to build output ensemble pf image\n",
					pfimage_size);

	/*Now we use string functions to assemble this mess */
	strcpy(pfimage,hdrblock);
	strcat(pfimage,"\n");
	strcat(pfimage,tag);
	strcat(pfimage," &Arr{\n");
fprintf(stderr,"%s\n",pfimage);
	for(i=0;i<npf;++i)
	{
		char key[16];
		sprintf(key,"%0.5d ",i);
		strcat(pfimage,key);
		strcat(pfimage,blocks[i]);
	}
	strcat(pfimage,"}\n");

	for(i=0;i<npf;++i) free(blocks[i]);
	free(blocks);
	free(hdrblock);
	pfcompile(pfimage,&pfresult);
	free(pfimage);
	return(pfresult);
}
/*  This pair of get and put functions are required because the 
current version of antelope truncates doubles on a put to 6 significant
figures.  As a workaround times can be treated as date strings
and cracked with Dan Quinlan's time string functions.

Author:  Gary Pavlis
Written:  October 2002
*/
double pfget_time(Pf *pf,char *name)
{
	char *s;
	double time;
	s=pfget_string(pf,name);
	time=str2epoch(s);
	return(time);
}
void pfput_time(Pf *pf,char *name,double time)
{
	pfput_string(pf,name,strtime(time));
}
/* This collection of functions put a single attribute of a 
particular type to all members of a pfensemble.  If I did
this in C++ it would be a classic reason for function
overloading.  Here I just have to do this the hard way.
Author:  G Pavlis
Written:  October 2002
*/
void pfensemble_put_double(Pf_ensemble *pfe,char *name,double value)
{
	int i;

	for(i=0;i<(pfe->nmembers);++i)
		pfput_double(pfe->pf[i],name,value);
}
void pfensemble_put_time(Pf_ensemble *pfe,char *name,double value)
{
	int i;

	for(i=0;i<(pfe->nmembers);++i)
		pfput_time(pfe->pf[i],name,value);
}
void pfensemble_put_int(Pf_ensemble *pfe,char *name,int value)
{
	int i;

	for(i=0;i<(pfe->nmembers);++i)
		pfput_int(pfe->pf[i],name,value);
}
void pfensemble_put_string(Pf_ensemble *pfe,char *name,char *value)
{
	int i;

	for(i=0;i<(pfe->nmembers);++i)
		pfput_string(pfe->pf[i],name,value);
}
