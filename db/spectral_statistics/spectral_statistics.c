#include <malloc.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "stock.h"
#include "db.h"
#include "pf.h"
#include "elog.h"
#include "arrays.h"
#include "spectrum.h"


#define MAX_FREQ 10000
#define SN_CUTOFF 4.0  /* power signal-to-noise cutoff */
Dbptr copy_dbptr(Dbptr db)
{
	Dbptr dbo;
	dbo.database = db.database;
	dbo.table = db.table;
	dbo.field = db.field;
	dbo.record = db.record;
	return(dbo);
}
Dbptr build_working_views(Dbptr dbspec, Dbptr dbev,
	char *sta, char *chan, char *sname, char *phase, char *pstype)
{
	Tbl *joinkeys;
	Dbptr dbo;
	int nrecs;
	char string[256];

	joinkeys = strtbl("arid",0);
	sprintf(string,"( sta =~ /%s/ ) && ( chan =~ /%s/ ) && (phase =~ /%s/ ) && (sname =~ /%s/ ) && ( pstype =~ /%s/ )",
		sta,chan,phase,sname,pstype);

	dbo = dbsubset(dbspec,string,0);
	if(dbo.record == dbINVALID)
		die(0,"noise power spectra subset failed with condition\n%s\n",
			string);
	dbquery(dbo,dbRECORD_COUNT,&nrecs);
	if(nrecs<=0)return(dbo);
	dbo = dbjoin(dbo,dbev,&joinkeys,&joinkeys,0,0,0);

	return(dbo);
}

void usage()
{
	die(0, "Usage: sc_spec_statistics db [-sift expression -pf pffile -v]\n");
}
void main(argc, argv)
int argc;
char **argv;

{
	char *dbin;
	Dbptr db;  /* Base db pointer */
	Dbptr dbwork;  /* Working copy of dbnoise*/
	Dbptr dbevent;  /* event based view of event->origin->assoc */
	Tbl *dbptbl;  /* Tbl passed to dbprocess to build dbevent */
	char *sta,*chan;
	/* These three are major groupings - parent is original power spectrum other obvious*/
	Dbptr dbparent,dbratio,dbnoise;  
	char *parent_spectrum_sname,*parent_spectrum_phase,*parent_spectrum_pstype;
	char *noise_sname,*noise_phase,*noise_pstype;
	char *ratio_sname,*ratio_phase,*ratio_pstype;

	char string[1024];  /* subset string buffer */

	int nparent,nratio,nnoise; /* number of records for different subsets*/

	int i;

	FILE *fp;
	char *fname;

	Pf *pf;
	char *pffile=NULL;
	int verbose=0;
	/* sift expression used to subset event view */
	char *sift=NULL;
  
	/* These are used by dbmatches */
	static Hook *hook1=0, *hook2=0;
	Tbl *kmatch,*matches=NULL;

	/* These are two parallel main work arrays.  ns is an array defining
	the current number of spectral ratios in frequency bin [i].  This
	is necessary because the program uses a signal/noise cutoff that
	makes this count frequency dependent.  The second array is a series 
	of pointers to alloced float arrays for each frequency bin.  */

	int ns[MAX_FREQ];
	float *s[MAX_FREQ];
	/* Added to fix resample problem.  */
	Spectrum pattern;
	double df,freq0;



	if(argc < 2)
	{
		usage();
		exit(1);
	}
	elog_init(argc,argv);
	dbin = argv[1];
        for(i=2;i<argc;++i)
        {
                if(!strcmp(argv[i],"-pf"))
                {
                        ++i;
                        pffile=argv[i];
                }
		else if(!strcmp(argv[i],"-v"))
			verbose = 1;
		else if(!strcmp(argv[i],"-sift"))
		{
			++i;
			sift = argv[i];
		}
		else
			usage();
        }
        if(pffile==NULL) pffile=strdup("spectral_statistics");
        if(pfread(pffile,&pf)) die(0,"pfread error\n");
	if(verbose)fprintf(stdout,"Spectral Statistics Program \n");

	/* base parameters acquired from pf file we get right away */
	freq0=pfget_double(pf,"freq0_reference");
	df=pfget_double(pf,"df_reference");

	kmatch = strtbl("sta","chan","orid",0);

	/* Open input database */

	if(verbose)fprintf(stdout,"Opening database %s\n",dbin);
	if (dbopen(dbin, "r", &db) == dbINVALID) {
		clear_register(1);
		fprintf(stderr, "dbspectra: Unable to open database '%s'.\n", dbin);
		exit(1);
	}
	/* This view is joined to each psdisc subset below so we build
	 it first */
	dbptbl = strtbl("dbopen event",
		"dbjoin origin",
		"dbsubset orid==prefor",
		"dbjoin assoc",0);
	dbevent = dbprocess(db,dbptbl,0);
	if(sift!=NULL)dbevent = dbsubset(dbevent,sift,0);
	if(dbevent.record==dbINVALID) die(0,"Failure building working event->origin->assoc view\n");

	db = dblookup(db,0,"psdisc",0,0);
	if (db.record == dbINVALID)
		die(0,"Error in psdisc table lookup\n");
	dbquery(db,dbRECORD_COUNT,&nratio);
	if(verbose)fprintf(stdout,"Input db has a total of %d rows in psdisc table\n",
			nratio);
	sta=pfget_string(pf,"station");
	chan=pfget_string(pf,"channel");
	parent_spectrum_pstype=pfget_string(pf,"parent_spectrum_pstype");
	parent_spectrum_phase=pfget_string(pf,"parent_spectrum_phase");
	parent_spectrum_sname=pfget_string(pf,"parent_spectrum_sname");
	noise_phase=pfget_string(pf,"noise_phase");
	noise_sname=pfget_string(pf,"noise_sname");
	noise_pstype=pfget_string(pf,"noise_pstype");
	ratio_phase=pfget_string(pf,"ratio_phase");
	ratio_sname=pfget_string(pf,"ratio_sname");
	ratio_pstype=pfget_string(pf,"ratio_pstype");

	/* We now build three working views of the parent power spectrum, spectral ratio,
	and noise spectrum estimates */
	dbparent = build_working_views(db,dbevent,
		sta,chan,parent_spectrum_sname,parent_spectrum_phase,parent_spectrum_pstype);
	dbquery(dbparent,dbRECORD_COUNT,&nparent);
	if(nparent <= 0)
		die(0,"Working view of parent signal power spectra is empty\nCheck database and parameter file\n");
	dbratio = build_working_views(db,dbevent,
		sta,chan,ratio_sname,ratio_phase,ratio_pstype);
	dbquery(dbratio,dbRECORD_COUNT,&nratio);
	if(nratio <= 0)
		die(0,"Working view of power spectral ratios is empty\nCheck database and parameter file\n");
	dbnoise = build_working_views(db,dbevent,
		sta,chan,noise_sname,noise_phase,noise_pstype);
	dbquery(dbnoise,dbRECORD_COUNT,&nnoise);
	if(nnoise <= 0)
		die(0,"Working view of noise power power spectral estimates is empty\nCheck database and parameter file\n");


	fname = pfget_string(pf,"output_filename");
	if((fp=fopen(fname,"w")) == NULL) die(0,"Cannot open output file %s\n");
	if(verbose)fprintf(stdout,"Writing results to output file %s\n");

	/* Scan the data view to get the number of frequencies to use in work space */
	for(dbratio.record=0,pattern.nfreq=0;
		dbratio.record<nratio;++dbratio.record)
	{
		int nf_test;
		if( dbgetv(dbratio,0,"nfreq",&nf_test,
			0) == dbINVALID) die(0,"dbgetv error scanning input view of data subset\n");
		pattern.nfreq=MAX(pattern.nfreq,nf_test);
	}
	if(pattern.nfreq<=0) 
		die(0,"Number of frequencies found in scanning input = %d is invalid\n",pattern.nfreq);
	else if(verbose) fprintf(stdout,"Using %d frequencies in work space \n",pattern.nfreq);

	pattern.freq0 = freq0;
	pattern.df = df;
	for(i=0;i<pattern.nfreq;++i)
	{
		allot(float *,s[i],nratio);
	}
	for(i=0;i<pattern.nfreq;++i) ns[i] = 0;

	/* Now we enter main loop */
	for(dbratio.record=0;dbratio.record<nratio;++dbratio.record)
	{

		Spectrum sr,parent,noise; /* working current spectra in 
					C structure used by read_spectrum*/
		int nmatches,orid;


		dbgetv(dbratio,0,"orid",&orid,0);
		sr = read_spectrum(dbratio);
		if(sr.spec == NULL)die(0,"read_spectrum failed on record %d of ratio view\n",
					dbratio.record);
		/* This obscure line uses a standard safe floating point
		equals test */

		if( (pattern.nfreq != sr.nfreq)
			&& ( (df - sr.df)/(df) != 0.0 )
			&& ( (freq0 - sr.freq0)/(freq0) != 0.0 ) )
		{
			sr = resample_spectrum(sr,pattern);
		}
		/* We now find the associated noise spectrum using dbmatches, which 
		finds the set of rows that match this sta:chan:orid */
		dbmatches(dbratio,dbnoise,&kmatch,&kmatch,&hook1,&matches);
		nmatches = maxtbl(matches);
		if(nmatches<=0) 
		{
			elog_notify(0,"No matching noise record for station:chan %s:%s with orid %d\n",
				sr.sta,sr.chan,orid);
			continue;
		}
		else if(nmatches>1)
		{
			elog_notify(0,"Warning:  %d matches in noise list for sta:chan:orid=%s:%s:%d\nUsing first one found\n",
				nmatches,sr.sta,sr.chan,orid);
		}
		
		dbwork = copy_dbptr(dbnoise);
		dbwork.record=(int)gettbl(matches,0);
		noise = read_spectrum(dbwork);
		if(noise.spec == NULL)die(0,"read_spectrum failed on record %d of ratio view\n",
					dbnoise.record);

		if( (pattern.nfreq != noise.nfreq)
			&& ( (df - noise.df)/(df) != 0.0 )
			&& ( (freq0 - noise.freq0)/(freq0) != 0.0 ) )
		{
			noise = resample_spectrum(noise,sr);
		}
		freetbl(matches,0);

		/* Now the same for parent spectrum */
		dbmatches(dbratio,dbparent,&kmatch,&kmatch,&hook2,&matches);
		nmatches = maxtbl(matches);
		if(nmatches<=0) 
		{
			elog_notify(0,"No matching parent spectrum record for station:chan %s:%s with orid %d\n",
				sr.sta,sr.chan,orid);
			continue;
		}
		else if(nmatches>1)
		{
			elog_notify(0,"Warning:  %d matches in parent spectrum list for sta:chan:orid=%s:%s:%d\nUsing first one found\n",
				nmatches,sr.sta,sr.chan,orid);
		}
		
		dbwork = copy_dbptr(dbparent);
		dbwork.record=(int)gettbl(matches,0);
		parent = read_spectrum(dbwork);
		if(parent.spec == NULL)die(0,"read_spectrum failed on record %d of ratio view\n",
					dbparent.record);
		if( (pattern.nfreq != parent.nfreq)
			&& ( (df - parent.df)/(df) != 0.0 )
			&& ( (freq0 - parent.freq0)/(freq0) != 0.0 ) )
		{
			parent = resample_spectrum(parent,sr);
		}
		freetbl(matches,0);
		

		/* Convert parent.spec array values to signal-to-noise ratios */
		for(i=0;i<parent.nfreq;++i) parent.spec[i] /= noise.spec[i];

		/* Now accumulate spectral ratio for each frequency, but
		skipping frequencies where signal-to-noise is inadequate */
		for(i=0;i<parent.nfreq;++i)
		{
			if(parent.spec[i] > SN_CUTOFF)
			{
				*(s[ i ] + ns[i]) = sr.spec[i];
				++ns[i];
			}
		}
		free(sr.spec);
		free(parent.spec);
		free(noise.spec);
	}


	/* Now loop through s and calculate statistics */
	if(verbose)fprintf(stdout,"Freq  Median\n");
	for(i=0;i<pattern.nfreq;++i)
	{
		double freq;
		float srmean,srmedian,sr1_4,sr3_4,srmin,srmax;

		freq = freq0 + df*( (float) i);

		if(ns[i] >= 4)
		{
			calc_statistics(s[i],ns[i],&srmean,&srmedian,
				&sr1_4,&sr3_4,&srmin,&srmax);
			fprintf(fp,"%lg %g %g %g %g %g %g %d\n",
				freq,srmean,srmedian,sr1_4,sr3_4,srmin,srmax,ns[i]);
			if(verbose)fprintf(stdout,"%lg %g %d\n", freq,srmedian,ns[i]);
		}
		else
		{
			if(verbose)fprintf(stdout,"Insufficient data to calculate statistics for frequency %lg\n",freq);
		}
	}
	exit(0);
}
	
