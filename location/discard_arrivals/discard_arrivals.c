#include <stdlib.h>
#include <math.h>
#include "stock.h"
#include "db.h"
#include "pf.h"

/* very crude program for randomly discarding data from 
a database produced by simcat.   The program loops through
the assoc table calling dbmark on a fraction of the arrivals
defined by a probability parameter "probability_of_deletion"".
Uses sun's random function as a simple random number generator.
Only other input parameter is the initial seed passed to
srandom */
void main (int argc, char **argv)
{
	char *dbin;
	Dbptr db;
	int nrec;
	int seed;
	long int test; 
	long int threshold;
	Pf *pf;
	double ptodel,dthresh;

	if(argc<2)
	{
		fprintf(stderr,"Usage:  discard_arrivals db\n");
		exit(-1);
	}
	dbin = argv[1];

	pfread("discard_arrivals",&pf);
	ptodel = pfget_double(pf,"probability_of_deletion");
	seed = pfget_int(pf,"random_number_seed");
	srandom((unsigned int)seed);
	dthresh = pow(2.0,31);
	dthresh -= 1.0;
	threshold = (long int )(dthresh*ptodel);

	dbopen(dbin,"r+",&db);
	db = dblookup(db,0,"assoc",0,0);
	dbquery(db,dbRECORD_COUNT,&nrec);

	for(db.record=0;db.record<nrec;++db.record)
	{
		test = random();
		if(test<=threshold)
			dbmark(db);
	}
	dbcrunch(db);
}
