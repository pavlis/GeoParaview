/* This program is a quick and dirty hack of the spectral_statistics program designed to 
produce data for plots of spectral amplitude ratios versus distance/travel time for
amplitude study of hindu kush events with Kopnichev.  Do not expect it to be elegant.
May 30, 2001
*/
#include <malloc.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "stock.h"
#include "db.h"
#include "pf.h"
#include "coords.h"
#include "elog.h"
#include "arrays.h"
#include "spectrum.h"
#include "tttaup.h"
#include "location.h"

Arr *MTSload_station_table(Dbptr);

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
	char *chan, char *sname, char *phase, char *pstype)
{
	Tbl *joinkeys;
	Dbptr dbo;
	int nrecs;
	char string[256];

	joinkeys = strtbl("arid",0);
	sprintf(string,"( chan =~ /%s/ ) && (phase =~ /%s/ ) && (sname =~ /%s/ ) && ( pstype =~ /%s/ )",
		chan,phase,sname,pstype);

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
	die(0, "Usage: sprt_band db [-sift expression -pf pffile -v]\n");
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

	int nfreq;
	double df,freq0;

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

        Arr *stalist;
	double fmin,fmax;  /*band to compute average results for */
	Station *reference,*current;
	char *refsta;
	double bavg;
	int count;
	double r0,r,d0,distance_ratio,hypod_ratio,t_ratio;
	double slat,slon,depth,distance,azimuth;


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
        if(pffile==NULL) pffile=strdup("sprt_band");
        if(pfread(pffile,&pf)) die(0,"pfread error\n");
	if(verbose)fprintf(stdout,"Spectral Statistics Program \n");

	/* base parameters acquired from pf file we get right away */
	freq0=pfget_double(pf,"freq0_reference");
	df=pfget_double(pf,"df_reference");
	fmin = pfget_double(pf,"fmin");
	fmax = pfget_double(pf,"fmax");

	kmatch = strtbl("sta","chan","orid",0);

	/* Open input database */

	if(verbose)fprintf(stdout,"Opening database %s\n",dbin);
	if (dbopen(dbin, "r", &db) == dbINVALID) {
		clear_register(1);
		fprintf(stderr, "dbspectra: Unable to open database '%s'.\n", dbin);
		exit(1);
	}
        stalist = MTSload_station_table(db);
	refsta = pfget_string(pf,"reference_station");
	reference=(Station *)getarr(stalist,refsta);
	if(reference==NULL) die(0,"Cannot find reference station %s in site table\n",
			refsta);
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

	if(!( (!strcmp(ratio_phase,"P")) || (!strcmp(ratio_phase,"S")) ))
		die(0,"ratio_phase parameter must be either P or S\nParameter file has %s\n",
			ratio_phase);

	/* We now build three working views of the parent power spectrum, spectral ratio,
	and noise spectrum estimates */
	dbparent = build_working_views(db,dbevent,
		chan,parent_spectrum_sname,parent_spectrum_phase,parent_spectrum_pstype);
	dbquery(dbparent,dbRECORD_COUNT,&nparent);
	if(nparent <= 0)
		die(0,"Working view of parent signal power spectra is empty\nCheck database and parameter file\n");
	dbratio = build_working_views(db,dbevent,
		chan,ratio_sname,ratio_phase,ratio_pstype);
	dbquery(dbratio,dbRECORD_COUNT,&nratio);
	if(nratio <= 0)
		die(0,"Working view of power spectral ratios is empty\nCheck database and parameter file\n");
	dbnoise = build_working_views(db,dbevent,
		chan,noise_sname,noise_phase,noise_pstype);
	dbquery(dbnoise,dbRECORD_COUNT,&nnoise);
	if(nnoise <= 0)
		die(0,"Working view of noise power power spectral estimates is empty\nCheck database and parameter file\n");


	fname = pfget_string(pf,"output_filename");
	if((fp=fopen(fname,"w")) == NULL) die(0,"Cannot open output file %s\n");
	if(verbose)fprintf(stdout,"Writing results to output file %s\n",fname);


	/* Now we enter main loop */
	for(dbratio.record=0;dbratio.record<nratio;++dbratio.record)
	{

		Spectrum sr,parent,noise; /* working current spectra in 
					C structure used by read_spectrum*/
		int nmatches,orid;


		dbgetv(dbratio,0,"orid",&orid,
			"origin.lat",&slat,
			"origin.lon",&slon,
			"origin.depth",&depth,0);
		sr = read_spectrum(dbratio);
		if(sr.spec == NULL)die(0,"read_spectrum failed on record %d of ratio view\n",
					dbratio.record);
		/* This obscure line uses a standard safe floating point
		equals test */
		nfreq = sr.nfreq;

/*
		if( (nfreq != sr.nfreq)
			&& ( (df - sr.df)/df != 0.0 )
			&& ( (freq0 - sr.freq0)/freq0 != 0.0 ) )
		{
			fprintf(stderr,"Warning:  parameter mismatch for row %d\n",
				dbratio.record);
			fprintf(stderr,"This entry will be skipped\n");
			free(sr.spec);
			continue;
		}
*/
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
				sr.sta,sr.chan,orid);
		}
		
		dbwork = copy_dbptr(dbnoise);
		dbwork.record=(int)gettbl(matches,0);
		noise = read_spectrum(dbwork);
		if(noise.spec == NULL)die(0,"read_spectrum failed on record %d of ratio view\n",
					dbnoise.record);

		if( (nfreq != noise.nfreq)
			&& ( (df - noise.df)/df != 0.0 )
			&& ( (freq0 - noise.freq0)/freq0 != 0.0 ) )
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
				sr.sta,sr.chan,orid);
		}
		
		dbwork = copy_dbptr(dbparent);
		dbwork.record=(int)gettbl(matches,0);
		parent = read_spectrum(dbwork);
		if(parent.spec == NULL)die(0,"read_spectrum failed on record %d of ratio view\n",
					dbparent.record);
		freetbl(matches,0);
 

		/* Convert parent.spec array values to signal-to-noise ratios */
		for(i=0;i<nfreq;++i) parent.spec[i] /= noise.spec[i];

		/* compute distance and time quantities */
		dist(rad(slat),rad(slon),rad(reference->lat),rad(reference->lon),
				&distance,&azimuth);
		d0 = deg2km(deg(distance));
		r0=hypot(depth,d0);
		current=(Station *) getarr(stalist,sr.sta);
		if(current==NULL)
		{
			elog_notify(0,"Station %s not in site table\nData skipped\n",
				sr.sta);
			free(sr.spec);
			free(parent.spec);
			free(noise.spec);
			continue;
		}
		dist(rad(slat),rad(slon),rad(current->lat),rad(current->lon),
			&distance,&azimuth);
                r=hypot(depth,deg2km(deg(distance)));
		distance_ratio = deg2km(deg(distance))/d0;
		hypod_ratio = r/r0;
		/* Very restrictive, but ok for a quick hack */
		if(ratio_phase[0]=='P')
			t_ratio = pphasetime(deg(distance),depth)
					/pphasetime(km2deg(d0),depth);
		else
			t_ratio = sphasetime(deg(distance),depth)
					/sphasetime(km2deg(d0),depth);

		/* Now accumulate spectral ratio for each frequency, but
		skipping frequencies where signal-to-noise is inadequate */
		for(i=0,bavg=0.0,count=0;i<nfreq;++i)
		{
			double freq;
			freq=freq0+df*((double)i);
			if( (freq<fmin) || (freq>fmax) ) continue;
			if(parent.spec[i] > SN_CUTOFF)
			{
				bavg += (double)(sr.spec[i]);
				++count;
			}
		}
		if(count>0)
		{
			bavg /= (double)count;
			fprintf(fp,"%s %lf %lf %lf %lf\n",
			  sr.sta,distance_ratio,hypod_ratio,t_ratio,bavg);
		}
		free(sr.spec);
		free(parent.spec);
		free(noise.spec);
	}

	exit(0);
}
	
