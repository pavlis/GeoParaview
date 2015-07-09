#include "stock.h"
#include "coords.h"
#include "db.h"

/* This program reads PASSCAL pcf files to build a clock rating table
stored in a css3.0 extension table.  It uses a set of average REFTEK
clock accuracy figures supplied by Jim Fowler by email in fall 2000.

Usage: 
	reftek_dbtiming db [-rampout] < pcffile

where db is a css3.0 database with two required extension tables 
(see below) and pcffile is a PASSCAL pcf file.  Multiple pcf files
can be concatenated together and passed through this program in 
a single step.  

The -rampout flag controls what happens at the end of the time 
period defined by the pcf file.  If the -rampout flag appears a
series of entries are placed in the timing table (see below) that
reflect the decreasing clock accuracy 

Two extension tables are required:

timing
This table contains clock accuracy figures keyed by 
station, channel, time:endtime.  

reftek
This table contains a cross reference between reftek 
serial numbers and station names for a particular time:endtime
(Borrowed from iceworm)

The actual email used to design this:  
Jim Fowler, October 10, 2000:
I am not sure how you can fix the clocks unless you have a beginning time
and an end time that are correct.  The nominal accuracy of the Ref Tek
oscillator is .5 ppm or about 40 ms/day.  The test we ran at Stanford was
to set the clock in the instrument at the beginning of the period and then
let it free run for about a week.  We measured the total drift at the end
of the week and then applied a linear correction to the time.  By doing so
we were with in 10 msec at any time during the week.  If you have a good
measure of the drift over a long time, you can do something like this.  You
can not just remove the nominal drift from the unit as each unit is
different.


The clock quality flag in the seed data goes as such:

TL = time since last clock lock

Quality         Meaning
100%            TL <= 1 hr
80%             1 hr < TL <= 2 hr
60%             2 hr < TL <= 4 hr
40%             4 hr < TL <= 8 hr
20%             8 hr < TL <= 1 day
10%             1 day < TL <= 7 days
5%              7 days < TL
0%              NO LOCK EVER


Author:  GAry Pavlis
Date:  January 2001
*/
/* This is made a global variable out of laziness.  Data are 
static definitions used in main and one function below
0..5 as above, 6=7 to 30 days, 7=30 to 90 days. 
All use 0.5 ppm precision for clock.  Note that this
program will generate a LOT of records if the first
component of this structure array has high<1hr 
because the reftek nominally locks every hour.  
For this reason I fudged the first time period up 
by a few minutes. */

struct clockerr {
		double low;
		double high;
		double error;
		double seedqual;
};
struct clockerr clkerr_ranges[8] = {	0.0, 4000.0, 0.0018, 100.0,
	4000.0, 7200.0, 0.0036, 80.0,
	7200.0, 14400.0, 0.0072, 60.0,
	14400.0, 28800.0, 0.0144, 40.0,
	28800.0, 86400.0, 0.0432, 20.0,
	86400.0, 604800.0, 0.3024, 10.0,
	604800.0, 2592000.0, 1.296,5.0,
	2592000.0, 7776000.0, 3.888,5.0};

#define NCLKFIELDS 8
double error_nolock=20.0; /* value used for over largest time above
				Assumes about 1 year maximum PASSCAL exp*/

/*------------functions used by main ----------------*/

void usage(char *prog)
{
	elog_die(0,"Usage error:\n  %s db [-rampout] < pcffile\n",prog);
}
/* Small function to look up station name related to a given Reftek
DAS serial number.  Uses the stage tables.

Arguments:
	db - open database to search
	das - serial number of requested das
	time - epoch time to match
Author:  GAry Pavlis
DAte:  January 2001
*/
char *lookup_refteksta(Dbptr db, int das, double time)
{
	Dbptr dbm;
	static Tbl *matches,*kmatch;
	static Hook *hook=0;
	int nmatches;
	char sta[10],*retsta;
	char ssident[15];

	db = dblookup(db,0,"stage",0,0);
	sprintf(ssident,"%d",das);
	dbm = dblookup(db,0,"stage",0,0);
	dbm.record = dbSCRATCH;
	dbputv(dbm,0,"ssident",ssident,"time",time,0);
	kmatch=strtbl("serial","time:endtime",0);
	dbmatches(dbm,db,&kmatch,&kmatch,&hook,&matches);
	nmatches = maxtbl(matches);
	if(nmatches>1) 
		elog_complain(0,"lookup_refsta:  database inconsistency in stage table.\n  %d rows of stage match serial=%d at time %s\n",
			nmatches,das,strtime(time));
	else if(nmatches<=0) 
	{
		*retsta = NULL;
	}
	else
	{
		db.record = (int)gettbl(matches,0);
		dbgetv(db,0,"sta",sta,0);
		retsta = strdup(sta);			
	}
	freetbl(matches,free);
	freetbl(kmatch,free);
	return(retsta);
}
/* Finds the epoch time of the end date of station.  Returns a negative
epoch time if the station has no offdate.

Arguments: 
	db - database to access (site table is hit)
	sta - station to check

*/
double get_station_closing_time(Dbptr db, char *sta)
{
	Dbptr dbm;
	static Tbl *matches,*kmatch;
	static Hook *hook=0;
	Tbl *sortkeys;
	int nmatches;
	int offdate;
	double time;

	db = dblookup(db,0,"site",0,0);
	dbm = dblookup(db,0,"site",0,0);
	dbm.record = dbSCRATCH;
	dbputv(dbm,0,"sta",sta,0);
	kmatch=strtbl("sta",0);
	sortkeys=strtbl("sta","ondate::offdate",0);
	db = dbsort(db,sortkeys,0,0);
	dbmatches(dbm,db,&kmatch,&kmatch,&hook,&matches);
	nmatches = maxtbl(matches);
	db.record = (int)gettbl(matches,nmatches-1);
	dbgetv(db,0,"offdate",&offdate,0);
	if(offdate<0)
		time = -1.0;
	else
		time = epoch(offdate);
	freetbl(matches,free);
	freetbl(kmatch,free);
	freetbl(sortkeys,free);
	return(time);
}
/* This function could probably have been avoided with a different logic
in the major loop of main, but this works.  It implements the rampout
feature in a manner perhaps somewhat more complex than necessary.  It
actually looks at the final offdate of a station.  It uses that as
an ending time mark if it is defined.  If it is null it ramp out to 
one year after the last entry in the pcf file.

Arguments:
	das - serial number of this das
	dbtime - endtime of last db entry saved for this station
	lastlock - last lock time
	etime - final lock time in pcf file.
	db - database pointer 
	rampout - logical to turn rampout feature on or off.
Author:  GAry Pavlis
DAte:  January 2001
*/
void closeout_das(int das, double dbtime, double lastlock, double etime,
	Dbptr db, int rampout)
{
	char *sta;
	double station_close_time;
	int i;
	double tend;

	if(rampout)
	{
		sta = lookup_refteksta(db,das,lastlock);
		station_close_time = get_station_closing_time(db,sta);
		/* A negative return means no enddate was set.  This allows
		a full rampout */
		if(station_close_time<0.0) 
		   station_close_time=lastlock+clkerr_ranges[NCLKFIELDS].high;
	}
	else
		station_close_time = etime;
	/* This totally duplicates the main loop of main*/
	for(i=0;i<NCLKFIELDS;++i)
	{
		if((station_close_time-lastlock) > clkerr_ranges[i].high)
			tend = lastlock + clkerr_ranges[i].high;
		else
			tend = station_close_time;
		if(dbaddv(db,0,"sta",sta,
			"time",dbtime,
			"endtime",tend,
			"clkerr",clkerr_ranges[i].high,
			"seedtqual",clkerr_ranges[i].seedqual,
					0)<0);
		{
			elog_complain(0,"dbaddv error on timing table for station %s for time interval:\n%s to %s\n",
				sta,strtime(dbtime),strtime(tend));
		}
		dbtime=tend;
		if(tend>station_close_time) break;
	}
	if(i>=NCLKFIELDS)
	{
		if(dbaddv(db,0,"sta",sta,
			"time",dbtime,
			"endtime",station_close_time,
			"clkerr",error_nolock,
			"seedtqual",0.0,
				0)<0);
		{
			elog_complain(0,"dbaddv error on timing table for station %s for time interval:\n%s to %s\n",
				sta,
				strtime(dbtime),
				strtime(station_close_time));
		}
	}
	free(sta);
	return;
}

void main(int argc, char **argv)
{
	char *dbname;
	Dbptr db,dbt;
	int i;
	char linein[256];
	int npcf_records=0;
	
	char timestr[40];
	int das_sn;
	int last_das;   /* We need to track the DAS serial number to 
			allow multiple DAS entries in a pcf stream */
	double lastlock, lastdbentry;  /* times of last lock and endtime of
					last db record. */
	int rampout=0;
	char *sta;
	double time,tend;  

	elog_init(argc,argv);
        elog_notify (0, "$Revision: 1.1 $ $Date: 2001/02/01 17:01:53 $") ;
	if(argc<2) usage(argv[0]);

	dbname = argv[1];
	dbopen(dbname,"r+",&db);

	for(i=2;i<argc;++i)
	{
		if(!strcmp(argv[i],"-rampout"))
			rampout=1;
		else
			usage(argv[0]);
	}
	dbt = dblookup(db,0,"timing",0,0);
	if(dbt.record == dbINVALID)
		elog_die(0,"%s requires extension tables called timing\n",
			argv[0]);
	while(fgets(linein,256,stdin)!=NULL)
	{
		if(linein[0]=='#') continue;
		/* Redundant, but better safe than sorry */
		if(strstr(linein,"#!PCF") != NULL) continue;
		sscanf(linein,"%d %s",&das_sn,timestr);
		time = str2epoch(timestr);
		if(npcf_records==0) 
		{
			last_das = das_sn;
			lastlock = time;
			lastdbentry = time;
		}
		++npcf_records;
		if(das_sn == last_das)
		{
			if( time < lastlock ) die(0,"PCF file for station %d is corrupted\nTime jumps backward from %s to %s\n",
				das_sn,strtime(lastlock),strtime(time));
			if( (time-lastlock)<=clkerr_ranges[0].high)
			{
			/* When clock is locking normally we will skip
			through here until we miss a lock building a long
			time span */
				lastlock=time;
			}
			else
			{
				/* It is safer to use the lastlock
				time instead of the current time to avoid
				closeout errors at the end of experiment*/
				sta = lookup_refteksta(db,das_sn,lastlock);
				for(i=0;i<NCLKFIELDS;++i)
				{
					if((time-lastlock)
					   > clkerr_ranges[i].high)
						tend = lastlock
						     + clkerr_ranges[i].high;
					else
						tend = time;
					if(dbaddv(dbt,0,"sta",sta,
						"time",lastdbentry,
						"endtime",tend,
						"clkerr",clkerr_ranges[i].high,
						"seedtqual",
						  clkerr_ranges[i].seedqual,
							0)<0);
					{
						elog_complain(0,"dbaddv error on timing table for station %s for time interval:\n%s to %s\n",
							sta,
							strtime(lastdbentry),
							strtime(tend));
					}
					lastdbentry=tend;
					if((time-lastlock)
						<=clkerr_ranges[i].high)
					{
						lastlock=time;
						break;
					}
				}
				if(i>=NCLKFIELDS)
				{
					if(dbaddv(dbt,0,"sta",sta,
						"time",lastdbentry,
						"endtime",time,
						"clkerr",error_nolock,
						"seedtqual",0.0,
							0)<0);
					{
						elog_complain(0,"dbaddv error on timing table for station %s for time interval:\n%s to %s\n",
							sta,
							strtime(lastdbentry),
							strtime(time));
					}
					lastdbentry = time;
					lastlock=time;
				}
				free(sta);
			}
		}
		else
		{
			closeout_das(last_das,lastdbentry,lastlock,time,
				dbt,rampout);
			lastdbentry = time;
			lastlock = time;
			last_das = das_sn;
		}
	}
	closeout_das(last_das,lastdbentry,lastlock,time,dbt,rampout);
	exit(0);
}

