#include <stdio.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sunmath.h>

#include "db.h"
#include "arrays.h"
#include "scv2.h"
/* These prototypes seem to be required */
void *my_malloc(char *,int);
void my_free(char *);
int SCV_trace_fixgap(Trace *,char*);
int
read_file (char *fname, int foff, char *datatype, int *nsamps, void **buf)

{
	FILE *file;
	int n, size;

	FILE *zopen();

	size = atoi(&datatype[strlen(datatype)-1]);
	file = zopen (fname, "r");
	if (file == NULL) {
		elog_notify(0, "read_file: Unable to open '%s'\n",
							fname);
		return (0);
	}
	if (fseek(file, foff, 0) < 0) {
		elog_notify(0, "read_file: fseek() error on '%s'\n",
							fname);
		fclose (file);
		return (0);
	}
	*buf = (void *) my_malloc ("read_file: *buf", size*(*nsamps));
	if (*buf == NULL) {
		elog_notify(0, "read_file: Malloc error.\n");
		fclose (file);
		return (0);
	}
	if ((n=fread (*buf, size, *nsamps, file)) < *nsamps) {
		if (n < 1) {
			elog_notify(0, "read_file: fread() error on '%s'\n",
							fname);
			fclose (file);
			my_free (*buf);
			return (0);
		} else {
			elog_notify(0, "read_file: Read %d samples but expected %d on '%s'\n",
							n, *nsamps, fname);
			*nsamps = n;
		}
	}
	fclose (file);
	return (1);
}

Trace *read_trace (Dbptr db, double tstart, double tend)

{
	char fname[1024];
	char dtype[8];
	char segtype[8];
	int foff, nsamp;
	void *data;
	Trace *trace;
	double time, dt, samprate;
	double calib, calper;
	int isamp, jsamp, size;

	if (dbextfile (db, "wfdisc", fname) < 1) {
		elog_notify(0, "read_trace: Unable to find input file '%s'\n",
						fname);
		return (NULL);
	}
	/* If this lookup failed, it needs to be a fatal error
		as something is seriously wrong */
	if( dbgetv (db, 0, "time", &time, "samprate", &samprate,
				"nsamp", &nsamp, "datatype", dtype, 
				"segtype", segtype, "foff", &foff, 
				"calib", &calib, "calper", &calper, 0)
				== dbINVALID)
		die(1,"read_trace error reading from database\n");


	/* isamp is number of samples from the start that is the time
	tstart.  The number to skip to reach this sample is isamp-1.
	jsamp is the index of the time tend */

	isamp = nint( (tstart-time)*samprate );
	jsamp = nint( (tend - time)*samprate);
	if (isamp < 0) isamp = 0;
	if (jsamp > nsamp) 
		nsamp -= isamp;
	else
		nsamp = nint( (tend - tstart)*samprate);
	if (nsamp < 0) nsamp = 0;
	size = atoi(&dtype[strlen(dtype)-1]);
	/* Note offset is isamp - because it specifies the number of
	samples to skip */
	foff += (isamp-1)*size;
	time += ( (double) isamp)/samprate;
	data = NULL;
	if (nsamp > 0) {
		if (!read_file (fname, foff, dtype, &nsamp, &data)) {
			elog_notify(0, "read_trace: Unable to read input file '%s'\n",
							fname);
			return (NULL);
		}
	}
	trace = (Trace *) my_malloc ("read_trace: trace", sizeof(Trace));
	if (trace == NULL) {
		elog_notify(0, "read_trace: Malloc error on Trace structure.\n");
		my_free (data);
		return (NULL);
	}
	trace->tstart = time;
	trace->dt = 1.0/samprate;
	trace->nsamps = nsamp;
	trace->calib = calib;
	trace->calper = calper;
	strcpy (trace->rawdata_format, dtype);
	strcpy (trace->rawdata_type, segtype);
	trace->data = NULL;
	trace->data_free = NULL;
	trace->data_malloc = 0;
	trace->raw_data = data;
	trace->rawdata_free = data;
	if (data) trace->rawdata_malloc = nsamp*size; else trace->rawdata_malloc = 0;
	trace->prev = NULL;
	trace->next = NULL;
	if (data) trace = (Trace *) SCV_trace_fixgaps (trace, "segment");
	return (trace);
}

