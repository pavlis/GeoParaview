/*Program to calculate spectral ratios from seismic 
array data.  Program does event oriented processing.  It assembles all
related spectra for each event.  From this it calculates the median 
array spectrum for that phase.  It then writes out this median spectrum,
calculates the spectral ratio for each station, and then stores the
spectral ratio data.  This is done by adding rows to the psdisc file and
creating a parallel set of disk files in a directory defined in 
the input parameter file.

Usage:

db_ratio dbin [-pf pffile -plot]

where
	dbin - input array database (unly psdisc is accessed)
	pffile - alternate parameter file to default db_ratio.pf
	-plot - turn on plotting (default is no plot)



Author:  Gary Pavlis
Written:  November 1994
Modified:  June 1996,  Changed schema to remove evid.  As a result the
comment above became invalid, and I removed the code that keyed on evid.
I discovered I misunderstood the schema and could not really join 
the origin table easily to the psdisc table by that method.  Instead
I now use arid in the same place.  With arid->assoc you can join origin.
Stupid, but then this schema was not designed for my needs.  
Modified:  May 2001
Major changes to allow use of this program with local and regional
network data.  Key changes;
1.  added parameter file interface (old was pure arg list)
2.  added distance amplitude correction term feature (see parameter file
and man page)
3.  The above required assimilation of a routine to read in a station
table.
4.  Got rid of all old style C function declarations
August 2001
REmoved plotting junk that just confused the code.
*/
#include <stdio.h>
#include <string.h>
#include "coords.h"
#include "stock.h"
#include "arrays.h"
#include "db.h"
#include "pf.h"
#include "location.h"
#include "db_ratio.h"

void usage()
{
	fprintf(stderr,"Usage:  \n db_ratio db [-plot -pf pffile]\n");
}

void free_spectrum_gather(Spectrum_gather g)
{
	int i;
	for(i=0;i<4;++i) free(g.chan[i]);
}

void free_spectrum(Spectrum s)
{
	free(s.spec);
}
void spectral_ratio(float *s1,float *s2,int ns)
{
	int i;
	for(i=0;i<ns;++i)
		s1[i] = s1[i]/s2[i];
}

int save_median(Dbptr db,Spectrum *s,int ns,float *m[4],char *threec[6])
{
	int i,k;
	int ii;
	int isave=-1;
	int ierr;

	/* Hunt for the first entry in s that matches each channel 
	name defined in the threec array.  Use this as the basic pattern
	and alter only parameters in s that are necessary for this
	entry.  */
	for(k=0;k<3;++k)
	{
		ii = -1;
		for(i=0;i<ns;++i)
		{
			if( (!strcmp(s[i].chan,threec[k])) && (s[i].index>=0) )
			{
				ii = i;
				isave = i;
				break;
			}
		}
		if(ii == -1) 
		{
			fprintf(stderr,"Warning (save_median):  error trying to save median for channel %s for event near time %s\n",
				threec[k],strtime(s[0].time));
			fprintf(stderr,"No match in loaded spectra.\n");
			continue;
		}
		/* Note before this code is called, dir is altered so we
		should stuff things in the right directly. All we should
		need to alter is sta, dfile, and pstype.  */
		strcpy(s[ii].sta,"MEDIAN");
		strcpy(s[ii].pstype,"psmedian");
		
		strcpy(s[ii].dfile,"");
		strcat(s[ii].dfile,epoch2str(s[ii].time,"%y%j%k%M%S"));
		strcat(s[ii].dfile,"_");
		strcat(s[ii].dfile,"MEDIAN");
		strcat(s[ii].dfile,"_");
		strcat(s[ii].dfile,threec[k]);
		strcat(s[ii].dfile,"_");
		strcat(s[ii].dfile,s[ii].sname);
		for(i=0;i<strlen(s[ii].dfile);++i) 
			if(s[ii].dfile[i] == ' ') s[ii].dfile[i] = '0';
		for(i=0;i<s[ii].nfreq;++i)
			s[ii].spec[i] = m[k][i];
		if((ierr=save_spectrum(db,s+ii)) != 0) 
			save_spectrum_error(s+ii,ierr);
	}
	/* Now we have to handle the 3c trace.  This time we just reuse
	the spec array that present ii points to.  We use the "isave"
	value set above as the last valid "ii".  If this value is
	negative, it means there were no matches, so we trap that
	case and return if that happens.  */
	if(isave < 0 )
	{
		fprintf(stderr,"No matching channel codes found to form 3C median output\nMedian 3C not save\n");
		return (4);
	}
	/* WARNING:  This only works with channel codes where the first
	character implies different types of channels (e.g. BHZ, EHZ, etc.)*/
	s[isave].chan[1] = '3';
	s[isave].chan[2] = 'C';
	s[isave].chan[3] = '\0';

	strcpy(s[isave].dfile,"");
	strcat(s[isave].dfile,epoch2str(s[isave].time,"%y%j%k%M%S"));
	strcat(s[isave].dfile,"_");
	strcat(s[isave].dfile,"MED");
	strcat(s[isave].dfile,"_");
	strcat(s[isave].dfile,s[isave].chan);
	strcat(s[isave].dfile,"_");
	strcat(s[isave].dfile,s[isave].sname);


	for(i=0;i<strlen(s[isave].dfile);++i) 
		if(s[isave].dfile[i] == ' ') s[isave].dfile[i] = '0';
	for(i=0;i<s[isave].nfreq;++i)
		s[isave].spec[i] = m[3][i];
	if((ierr=save_spectrum(db,s+isave)) != 0) 
		save_spectrum_error(s+isave,ierr);
	return(0);
}

#define MAX_CHANNELS 800 /* sets maximum number of spectra expected per gather*/

void main (argc,argv)
 
int argc;
char **argv;
 
{
	char *dbin;
	char *outdir=NULL;
	int station_ratio=0;  /* set to one as a flag to use a named 
				station to form ratios instead of median */
	char *sta_denom;  /*Set to argv value for -sta option */
	int i,j,k;
	int nspec;  /* total number of spectra in input db */
	int apply_distance_correction;
	double decay_constant;

	Dbptr db;
	int n_spectra;  /* count of number of spectra read this gather*/
	Tbl *process_list;
	Pf *pf;
	char *pffile=NULL;

	double time;
	int last_evid;
	int evid;
	char last_sname[8];
	char *threeC_channels[6];
	Tbl *chanlist_pf;

	Spectrum_gather gather;
	Spectrum threeCspec[MAX_CHANNELS];
	int n3cspec;

	struct ratio_indices {
		int chan1, chan2, chan3, threeC;
	}  ratio_indices;

	Spectrum spec[MAX_CHANNELS];
	int ierr;
	
	Station *current,*reference;


	float *frequencies;
	Arr *stalist;

	if(argc < 2)
	{
		usage();
		exit(1);
	}
	dbin = argv[1];
	for(i=2;i<argc;++i)
	{
		if(!strcmp(argv[i],"-pf"))
		{
			++i;
			pffile = argv[i];
		}
		else
		{
			fprintf(stderr,"Unrecognized argument %s", argv[i]);
			usage();
			exit(1);
		}
	}
	/* This reads stuff from the parameter file. Would be better form
	to make this a seperate function, but so what */
	if(pffile==NULL)pffile=strdup("db_ratio");
	if(pfread(pffile,&pf)) die(0,"Error in pf file %s\n",pffile);
	outdir = pfget_string(pf,"output_directory");
	chanlist_pf = pfget_tbl(pf,"channel_list");
	if(maxtbl(chanlist_pf)!=3)die(0,"channel_list must have exactly 3 lines\n in pf file\n");
	for(i=0;i<3;++i)
	{
		char *line,primary[10],auxiliary[10];
		line = (char *)gettbl(chanlist_pf,i);
		sscanf(line,"%s %s",primary,auxiliary);
		threeC_channels[i]=strdup(primary);
		threeC_channels[i+3]=strdup(auxiliary);
	}
		
	station_ratio = pfget_boolean(pf,"use_reference_station");
	sta_denom = pfget_string(pf,"reference_station");
	apply_distance_correction = pfget_boolean(pf,"apply_distance_correction");
	if(apply_distance_correction) 
		decay_constant = pfget_double(pf,"decay_constant");
	
	if(makedir(outdir) == -1)
	{
		die(0,"Cannot open spectral ratio output directory %s\n",
			outdir);
        }
        if (dbopen (dbin, "r+", &db) == dbINVALID) {
		die(0,"db_ratio:  Unable to open database '%s'.\n", dbin);
	}
	stalist = MTSload_station_table(db);
	
	reference = (Station *)getarr(stalist,sta_denom);
	if(reference==NULL) die(0,"Requested reference station %s not in site table\n",
					sta_denom);

	process_list = strtbl("dbopen event",
                "dbjoin origin",
                "dbsubset orid==prefor",
                "dbjoin assoc",
                "dbjoin arrival",
                "dbjoin psdisc arid",
		"dbsort evid phase sname sta psdisc.chan",0);
        db = dbprocess(db,process_list,0);

	dbquery(db,dbRECORD_COUNT,&nspec);
	if(nspec <= 0)
		die(0,"db_ratio:  Empty working db view\n");

/* A new event is defined by a change in evid */

	for(db.record=0; db.record<nspec;++db.record)
	{
		int k;  /* index for spec array */
		float *median[4];  /* Holds array medians for 3c and total */
		double slat,slon,sdepth,distance,azimuth,r0,r,ampfactor;
		char phase[10],last_phase[10];

		k=0;
		if( dbgetv(db,0,"evid",&evid,
			"origin.lat",&slat,
			"origin.lon",&slon,
			"origin.depth",&sdepth,
			"arrival.iphase",last_phase,
                	"sname",last_sname,0) ==dbINVALID)
		{
			fprintf(stderr,"Fatal error(db_ratio):  dbgetv failure reading working view\n");
			exit(1);
		}
		if(db.record==0) 
			last_evid=evid;
		if(evid!=last_evid)
		{
			dist(rad(slat),rad(slon),
			rad(reference->lat),rad(reference->lon),
			&distance,&azimuth);
			r0=hypot(sdepth,deg2km(deg(distance)));
			last_evid=evid;
		}
			
		do
		{
			spec[k] = read_spectrum(db);
			if(spec[k].spec == NULL)
			{
				fprintf(stderr,"Fatal error in read_spectrum\n");
				fprintf(stderr,"Problem reading record %d in sorted psdisc table\n",db.record);
				
				exit(1);
			}
			if(apply_distance_correction)
			{
				current=(Station *)getarr(stalist,spec[k].sta);
				if(current!=NULL)
				{
					dist(slat,slon,
					  current->lat,current->lon,&distance,&azimuth);
					r=hypot(sdepth,deg2km(deg(distance)));
					ampfactor = pow(r/r0,decay_constant);
					/* we could use the blas here, but why complicate
					the makefile */
					for(i=0;i<spec[k].nfreq;++i)
						spec[k].spec[i] *= ampfactor;				
				}
				else
					elog_complain(0,
					  "Warning:  station %s not found in site table\nDistance corrections cannot be applied\n",
					  spec[k].sta);
			}
			++db.record;
			++k;
			if(k>MAX_CHANNELS)
			{
				fprintf(stderr,"Fatal(db_ratio):  too many channels.\n");
				fprintf(stderr,"Static space is %d and program hit this wall\n",MAX_CHANNELS);
				exit(1);
			}
		}
		while( (!strcmp(spec[k-1].sname,last_sname) )
			&& (evid == last_evid)
			&& (!strcmp(spec[k-1].phase,last_phase))
			&& (db.record < nspec) );

		n_spectra = k - 1;
		if(db.record<nspec) db.record-=2;


						
		/* Now we have a "gather" of spectra for the same sname
		with all channels.  The routine called here does the nasty
		job of sorting out the right set of 3 components and 
		assembling the structure of 2-d arrays called a 
		"Spectrum_gather" in db_ratio.h  This routine does most
		of the messy work for this code.
		The conditional here on station_ration determines 
		handling of the pattern for mixed sample rates and
		window sizes.  That is, in station_ratio mode the
		reference station is used to determine the pattern
		of frequency binning.  Otherwise a global search
		of the group is used. */

		if(station_ratio)
			gather = data_gather(spec,n_spectra,sta_denom,
				threeC_channels,threeCspec,&n3cspec);
		else
			 gather = data_gather(spec,n_spectra,NULL,
				threeC_channels,threeCspec,&n3cspec);
		if(gather.n1gather<=1) continue;
		if(gather.nfreq < 0) 
		{
			fprintf(stderr,"Fatal error (db_ratio):  memory allocation failure in data_gather\n");
			exit(1);
		}

		/* This is an error log, and an useless comment */
		printf("%s %s %s %d %d %d %d\n", 
			epoch2str(spec[0].time,"%y%j:%k%M%S"),
			spec[0].phase, spec[0].sname,
			gather.nstation[0],
			gather.nstation[1],
			gather.nstation[2],
			gather.nstation[3]);

		/* If we are using a station for the spectral ratio 
		denominator, scan the spec array of structures to find
		the proper indices using a dumb linear search*/
		if(station_ratio) 
		{
			ratio_indices.chan1 = -1;
			ratio_indices.chan2 = -1;
			ratio_indices.chan3 = -1;
			ratio_indices.threeC = -1;
			for(k=0;k<n_spectra;++k)
			{
				if(!strcmp(spec[k].sta,sta_denom))
				{
					if(!strcmp(spec[k].chan,
					  threeC_channels[0])) ratio_indices.chan1 = k;
					if(!strcmp(spec[k].chan,
					  threeC_channels[1])) ratio_indices.chan2 = k;
					if(!strcmp(spec[k].chan,
					  threeC_channels[2])) ratio_indices.chan3 = k;
				}
			}
			/* Here we try to use alternative channel codes when the 
			primary channel codes are missing.  Same dumb linear search */
			if(ratio_indices.chan1 < 0)
			{
				for(k=0;k<n_spectra;++k)
					if( (!strcmp(spec[k].sta,sta_denom))
					  && (!strcmp(spec[k].chan,threeC_channels[3])) )
						ratio_indices.chan1 = k;
			}
			if(ratio_indices.chan2 < 0)
			{
				for(k=0;k<n_spectra;++k)
					if( (!strcmp(spec[k].sta,sta_denom))
					  && (!strcmp(spec[k].chan,threeC_channels[4])) )
						ratio_indices.chan2 = k;
			}
			if(ratio_indices.chan3 < 0)
			{
				for(k=0;k<n_spectra;++k)
					if( (!strcmp(spec[k].sta,sta_denom))
					  && (!strcmp(spec[k].chan,threeC_channels[5])) )
						ratio_indices.chan3 = k;
			}
			/* the same search is easier for the threeC spectra */
			for(k=0;k<n3cspec;++k)
			{
				if(!strcmp(threeCspec[k].sta,sta_denom))
				{
					ratio_indices.threeC = k;
					break;
				}
			}
		}


		/* now calculate array median at each frequency
		for each channel.  More messy pointer arithmetic,
		but the idea is each chan is stored in 2d C array
		with station being the last index.  i.e. each 
		frequency is stored in an adjacent element of the
		chan arrays.*/
		for(k=0;k<4;++k) if( (median[k]
		  = (float *) calloc(gather.nfreq,sizeof(float))) == NULL)
		{
			fprintf(stderr,"Fatal (db_ratio):  Cannot alloc median array of size %d\n", gather.n1gather);
			exit(1);
		}
		for(k=0;k<4;++k)
		{
			int jj;
			for(j=0;j<gather.nfreq;++j)
			{
				jj = j*gather.n1gather;
				median[k][j] = median_float((gather.chan[k])+jj,
						gather.nstation[k]);
			}
		}
		/* Now we have to convert spec to spectral ratios.  There are now
		two cases.  Station mode and median ratios.  In the first case
		the indexing is simpler, but we have to be careful to not apply the 
		ratios to the special spectrum that is used in the denominator for
		each channel.  In both cases, spec[].spec vectors are all altered.
		We do this in a stupid way passing through all the data
		three times to match corresponding channels.  it makes this slower,
		but much simpler and clearer. Note that for the denominator 
		trace is not saved. */

		if(station_ratio)
		{
		    if(ratio_indices.chan1 >= 0)
		    {
			for(j=0;j<n_spectra;++j)
			{
				if( (!strcmp(spec[j].chan,threeC_channels[0]))
				  || (!strcmp(spec[j].chan,threeC_channels[3])))
				{
				    if(strcmp(spec[j].sta,sta_denom))
				    {
					spectral_ratio(spec[j].spec,spec[ratio_indices.chan1].spec,
						spec[j].nfreq);
					strcpy(spec[j].dir,outdir);
					strcpy(spec[j].pstype,"1cpsr");
					if( (ierr=save_spectrum(db,spec+j)) != 0)
						save_spectrum_error(spec+j,ierr);
				    }
				}
			}
		    }
		    if(ratio_indices.chan2 >= 0)
		    {
			for(j=0;j<n_spectra;++j)
			{
				if( (!strcmp(spec[j].chan,threeC_channels[1]))
				  || (!strcmp(spec[j].chan,threeC_channels[4])))
				{
				    if(strcmp(spec[j].sta,sta_denom))
				    {
					spectral_ratio(spec[j].spec,spec[ratio_indices.chan2].spec,
						spec[j].nfreq);
					strcpy(spec[j].dir,outdir);
					strcpy(spec[j].pstype,"1cpsr");
					if( (ierr=save_spectrum(db,spec+j)) != 0)
						save_spectrum_error(spec+j,ierr);
				    }
				}
			}
		    }
		    if(ratio_indices.chan3 >= 0)
		    {
			for(j=0;j<n_spectra;++j)
			{
				if( (!strcmp(spec[j].chan,threeC_channels[2]))
				  || (!strcmp(spec[j].chan,threeC_channels[5])))
				{
				    if(strcmp(spec[j].sta,sta_denom))
				    {
					spectral_ratio(spec[j].spec,spec[ratio_indices.chan3].spec,
						spec[j].nfreq);
					strcpy(spec[j].dir,outdir);
					strcpy(spec[j].pstype,"1cpsr");
					if( (ierr=save_spectrum(db,spec+j)) != 0)
						save_spectrum_error(spec+j,ierr);
				    }
				}
			}
		    }
		}
		else
		{


		/*Block for median spectral ratio case.  It
		gets messy again because of different components we have
		to keep straight.  That is, we have to associate the correct
		median spectrum to form the spectral ratio.  We only 
		convert spectra marked as being used by data_gather through
		a nonnegative spec.index value.  Note the spectral_ratio 
		is destructive as it overwrites spec[j] values with 
		the spectral ratios. This block is for median ratios*/
		    for(j=0;j<n_spectra;++j)
		    {
			if(spec[j].index >= 0)
			{
				spectral_ratio(spec[j].spec,
					median[spec[j].index],spec[j].nfreq);
				strcpy(spec[j].dir,outdir);
				strcpy(spec[j].pstype,"1cpsr");
				/* for individual components we will leave the
				dfile names intact, but write ratios in a
				seperate directory */
				if( (ierr=save_spectrum(db,spec+j)) != 0)
					save_spectrum_error(spec+j,ierr);
			}
		    }
		}
		/* Now we have to do a slightly different thing for the total
		power spectra.  The variation is that we have to change the
		file name slightly to avoid collision with 1c values.  */
		for(j=0;j<gather.nstation[3];++j)
		{
			char *cptr;
			/*We first save the total power spectra.
			 We label this channels as first character of chan
			name + "3C".  Thus, for instance, EHZ -> E3C.  Not
			very good form because it works only with SEED style 
			channel codes */

			threeCspec[j].chan[1] = '3';
			threeCspec[j].chan[2] = 'C';
			threeCspec[j].chan[3] = '\0';
			strcpy(threeCspec[j].dir,outdir);
			strcpy(threeCspec[j].pstype,"sp3c");
			if( (ierr = save_spectrum(db,threeCspec+j)) != 0)
				save_spectrum_error(threeCspec+j,ierr);


			if(station_ratio)
			{
				if(ratio_indices.threeC < 0)  continue;
				spectral_ratio(threeCspec[j].spec,
					threeCspec[ratio_indices.threeC].spec,
					threeCspec[j].nfreq);
			}
			else
			{
			/* Here we don't need to check index because we 
			already figured out which stations don't have all
			three components alive */
				spectral_ratio(threeCspec[j].spec,median[3],threeCspec[j].nfreq);
			}
			strcpy(threeCspec[j].pstype,"3cpsr");

			if( (ierr = save_spectrum(db,threeCspec+j)) != 0)
				save_spectrum_error(threeCspec+j,ierr);
		}
		/* we want to save the median spectra too.  To do this 
		we call this function which recycles an appropriate 
		element of the spec array of structure.  This approach is
		a little weird, but it is, to my mind, the easiest way 
		to do this.  The function is defined above */
		save_median(db,spec,n_spectra,median,threeC_channels);

		for(j=0;j<n3cspec;++j) free_spectrum(threeCspec[j]);
		for(j=0;j<n_spectra;++j) free_spectrum(spec[j]);
		free_spectrum_gather(gather);
		for(k=0;k<4;++k) free(median[k]);	
	}
	exit(0);
}		
