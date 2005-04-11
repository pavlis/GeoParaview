#include <stdio.h>
#include <math.h>

#include "stock.h"
#include "pf.h"
#include "db.h"

#define STATION_MODE 1
#define PLANE_WAVE_MODE 2
#define DEG_TO_RADIANS 0.01745329


/*  dbawf

Usage:
  dbawf dbin [-zero -pf pffile] 

where
	dbin - input db to process 
	when -zero appears, each trace is zeroed before accumulation
	-pf -pffile - use alternative parameter file to default dbawf.pf

This program takes its name from db Accumulate WaveForms.  It takes the 
traces associated with dbin and adds (accumulates) waveforms as 
specified in a parameter file.  The program
runs in one of two basic modes.  (1) It can take a base waveform that is
scaled and time shifted according to a specified plane wave azimuth
and slowness.  These waveforms are added to existing contents of all 
traces found in dbin.  (2)   In station mode specified waveform files
are scaled and added at specified time lags. 

Behaviour is controlled by parameter "mode".  A Tbl called
"process_list" contains lines of instructions for constructing
desired results.



When mode is "STATION", following lines should be a series
of lines in the following format:
When 

station channel  weight lag  filename format foff

where
	station, channel = station-channel code (must match with wfdisc) 
		specification
	weight = multiplier for waveform (i.e. waveform "filename" values
		are all multiplied by this constant)  Like a gain factor,
		but also can be used to flip polarities. 
	lag = time in seconds from start of trace to translate the base
		waveform to.
	filename = waveform trace input file name.  
	format = data file format.  Valid values are: as, t4, s2, and s4.
		The later three are standard CSS binary format descriptions.
		If "as" is set, the file is assumed to contain
		a sequence of ascii numbers defining the base waveform.
		In all cases it is blatantly assumed
		assumed that the sampling interval of this
		waveform equals the sample interval of the original trace.
		No checking is made for consistency in the wfdisc files either.
	foff - file offset ala CSS wfdisc files.  Ignored for ascii input.

Note that stations can be specified multiple times with different weights,
lags, and/or filenames.  Results are just continually accumulated.

When mode is PLANE_WAVE, following lines should be in this format:

channel slowness azimuth weight lag filename format foff

where
	channel is as above.  All stations with this channel will be 
		processed in this mode.  
	slowness and azimuth define how the plane wave time shifts will be
		calculated
	weight - is as above
	lag = time shift relative to start time of "array reference station" 
		for placing waveform to be accumulated.  At other stations
		lag = lag + plane_wave_timeshift
	filename, format, foff = as above

The mode switches (PLANE_WAVE or STATION) can appear multiple times
in the file.  Whenever they appear the operation mode is switched.

Author:  Gary Pavlis
Written:  June 1994

Note:  August 1995
A warning to potential users.  Keep in mind this program takes the
trace data you use and converts it back to the native format of the
original wfdisc table.  This means it is possible to get odd results
if, for example, the input wavelet has numbers of around 1 (in float)
and the trace data is stored as s2 or s4 (ints).  When the input is
converted back internally to the int format, the waveform can be 
truncated and loose precision.  
*/ 

float *read_trace();
int usage()
{
	fprintf(stderr,"Usage: \n  dbawf dbin [-zero -pf pffile]\n");
	fprintf(stderr,"WARNING:  traces pointed to by dbin will be altered\n");
	exit(-1);
}
/* this is a function only because I stole it from somewhere else. */
void read_wf_error(char *file,char *datatype,int error_code)
{
        fprintf(stderr,"Warning:  error reading input waveform file %s from read_wf\n",file);
        switch(error_code)
        {
        case (-1):
                fprintf(stderr,"Cannot open data file\n");
                break;
        case (-2):
                fprintf(stderr,"Unknown data type %s\n",datatype);
                break;
        default:
                fprintf(stderr,"Unrecognized error code %d",error_code);
        }
	exit(-1);
}


main (int argc, char **argv)

{
	int zero;  /* set to one if input traces are to be initialized to zero*/
	Dbptr dbin;  /* input database pointer */
	char *db_name;  /* database root name */
	/* These are values read from input (recipe) file */
	float weight, lag; 
	double slowness, azimuth;
	char station[8], channel[10];
	int mode;
	char *line;

	char wf_file[128];  /* Waveform input file name */
	FILE *wfp;  /* fp for wf_file */
	FILE *rfp;  /* recipe file pointer */
	double sta_x, sta_y;  /* station coordinates relative to origin */
	float total_lag;  /* computed total lag for current station */
	float *wf;  /* holds waveform to be accumulated */
	int nwf;  /* length of wf */
	int i,j;

	double ux, uy;
	char *pfname="dbawf";
	char *modestr;
	Tbl *ptbl;
	Pf *pf;

	zero = 0;
	rfp = stdin;
	if(argc < 2)
	{
		usage();
	}
	db_name = argv[1];
	for(i=2;i<argc;++i)
	{
		if(strcmp(argv[i],"-zero")==0) 
			zero = 1;
		else if(strcmp(argv[i],"-pf")==0)
		{
			++i;
			pfname=argv[i];
		}
		else
		{
			usage();
		}
			
	}
	elog_init(argc,argv);
	/* Open input database */
	
	if (dbopen (db_name, "r+", &dbin) == dbINVALID) {
                clear_register (1);
                fprintf (stderr, "dbawf: Unable to open database '%s'.\n", dbin);
                exit (1);
        }
	if(pfread(pfname,&pf))
		elog_die(1,"pfread failed\n");
	modestr=pfget_string(pf,"mode");
	if(modestr==NULL)
		elog_die(0,"Parameter file error: must specify mode=STATION or PLANE_WAVE\n");
	if(strncmp(modestr,"STATION",7)==0)
		mode = STATION_MODE;
	if(strncmp(modestr,"PLANE_WAVE",10)==0)
		mode = PLANE_WAVE_MODE;
	ptbl=pfget_tbl(pf,"process_list");
	if(ptbl==NULL)
		elog_die(1,"Parameter file error:  missing process_list Tbl\n");

	/* Now we loop though the process list line by line */
	for(i=0;i<maxtbl(ptbl);++i)
	{
		Dbptr dbs;  /* pointer to site table database */
		Dbptr dbwf;  /*pointer to wfdisc table subset */
		int nwfdisc; /* number of rows in dbwf */
		char fname[1024];  /* trace file name */

		int ns;
		int foff;
		char dataform[10];

		line = gettbl(ptbl,i);
		/* Probably not necessary as pfread strips these I think*/
		if(line[0]=='#') continue;
		switch(mode)
		{
		case STATION_MODE:
			sscanf(line,"%s%s%g%g%s%s%d",
				station,channel,&weight,&lag,wf_file,
				dataform,&foff);
			total_lag = lag;			
			break;
		case PLANE_WAVE_MODE:
			sscanf(line,"%s%lg%lg%g%g%s%s%d",
				channel,&slowness,&azimuth,&weight,
				&lag,wf_file,dataform,&foff);
			ux = slowness*sin(DEG_TO_RADIANS*azimuth);
			uy = slowness*cos(DEG_TO_RADIANS*azimuth);
			break;
		}
		/* bring in waveform to be accumulated */
		nwf = read_wf(wf_file,&wf,dataform,foff);
		if(nwf<0) read_wf_error(wf_file,dataform,nwf);
		/* for either mode we need to load the wfdisc table */
		dbwf = dblookup(dbin,0,"wfdisc",0,0);

		switch(mode)
		{
		case STATION_MODE:
		/*  Here we subset the wfdisc table to look only 
			at those records matching station and channel */
			sprintf(line,"(sta == \"%s\") && (chan == \"%s\")",
				station,channel);
			dbwf = dbsubset(dbwf,line,0);
			break;
		case PLANE_WAVE_MODE:
		/* Here we just select out the correct channel */
			sprintf(line,"(chan == \"%s\")",channel);
			dbwf = dbsubset(dbwf,line,0);
		}

		/* Now loop through the wfdisc file calculating lags for
		array mode */
		dbs = dblookup(dbin,0,"site",0,0);
		dbquery(dbs,dbRECORD_COUNT,&ns);
		dbquery(dbwf,dbRECORD_COUNT,&nwfdisc);
		for(dbwf.record=0;dbwf.record<nwfdisc;dbwf.record++)
		{
			double samprate,calib;
			int nsamp,foff;
			char dtype[8];
			float *trace;
			Dbptr db_this_station;
			int lag_in_samples;
			int nadd;
			char fname[1024];

			if(dbextfile(dbwf,"wfdisc",fname) < 1)
			{
				fprintf(stderr,"Unable to find trace file %s\n",
					fname);
				continue;
			}
			dbgetv(dbwf,0,"samprate",&samprate,
				"nsamp",&nsamp,"datatype",dtype,
				"calib",&calib,"foff",&foff,0);
			/* Because of the intended use of this routine,
			I do not worry about data gaps and that jazz here.*/
			trace = read_trace(fname,foff,dtype,&nsamp);
			if(trace == NULL)
				continue;  /* read_trace generates error
					messages so none are neede here*/
			if(zero && (i==0))
				for(j=0;j<nsamp;++j) trace[j]=0.0;
			else if(calib>0.0) 
				scale_trace(nsamp,trace,calib);
			if(mode == PLANE_WAVE_MODE)
			{
				dbgetv(dbwf,0,"sta",station,0);
				for(dbs.record=0;dbs.record<ns;++dbs.record)
				{
					char sta_test[8];
					dbgetv(dbs,0,"sta",sta_test,
					  "dnorth",&sta_y,"deast",&sta_x,0);
					if(!strcmp(station,sta_test)) break;
				}
/* the following section of code should be equivalent to the above, 
but for some unknown reason it does not work correctly */
/*
				sprintf(line,"sta =~ /%s/",station);
				dbs.record = 0;
				dbs=dbsubset(dbs,line,0);
				dbquery(dbs,dbRECORD_COUNT,&ns);
				dbs.record = 0;
				dbgetv(dbs,0,0,"dnorth",&sta_y,
					"deast",&sta_x,0);
*/
				total_lag = lag - (ux*sta_x+uy*sta_y);
			}
			lag_in_samples = nint(total_lag*samprate);
			fprintf(stdout,"%s->%s lag is %f sec = %d samples\n",
				station, channel, total_lag, lag_in_samples);
			if((nadd = add_wf_shifted(nsamp,trace,nwf,weight,
			wf,lag_in_samples)) != nwf) fprintf(stderr,
				"Warning:  only %d of %d samples of file %s added\n",
					nadd, nwf, fname);
			/* This rescales trace back to make calib still 
			valid */
			if(calib>0.0) scale_trace(nsamp,trace,1.0/calib);
			write_trace(fname,foff,dtype,nsamp,trace);
		}
	}
}
