#include <malloc.h>
#include <stdio.h>
#include <math.h>
#include "db.h"
#include "db_ratio.h"
#define FLOAT_TEST_RATIO 0.0001

Spectrum Spectrum_copy(Spectrum s)
{
	Spectrum sout;

	if( (sout.spec = (float *) calloc(s.nfreq,sizeof(float))) == NULL) return(sout);

	strcpy(sout.sta,s.sta);
	strcpy(sout.sname,s.sname);
	strcpy(sout.chan,s.chan);  /* Note db_ratio depends on this field
				being set this way to rename dfile. It
				then changes this for the db output */
	strcpy(sout.phase,s.phase);
	strcpy(sout.pstype,"3cps");
	strcpy(sout.dir,s.dir);
	strcpy(sout.dfile,s.dfile);
	strcpy(sout.datatype,"t4");
	strcpy(sout.auth,s.auth);

	sout.arid = s.arid;
	sout.time = s.time;
	sout.chanid = s.chanid;
	sout.jdate = s.jdate;
	sout.endtime = s.endtime;
	sout.freq0 = s.freq0;
	sout.nfreq = s.nfreq;
	sout.df = s.df;
	sout.rayleigh = s.rayleigh;
	sout.tbp = s.tbp;
	sout.sscale = 1.0;
	sout.index = 3;  /* Here this is always 3 for total power */
	return(sout);
}
	
/* This function does the messy job of assembling spectra from a "gather"
defined as a group of related spectra from one event (e.g. P spectra for
all components from event x).  It gets messy because of the desire I 
had to assemble the 3component traces and handle the most general case
with high and low gain channels.  I wanted to automatically use low 
gains if the high gains are clipped.  

Args:

spec is the input spectrum array of structures of length nspec.  These
are the data assembled as the output.

chanids - 6 element char array holding channel codes we key on.  That is,
	chanids[0] to 2 are normally high gain channel codes, while
	chanids[3] to 5 are low gain channel codes.
threeCspec is an output "Spectrum" array of structures containing 3component
	power spectra estimated from the sum of the three components.
nthreeCspec - is the number of elements of the output threeCspec array 

Main return is the "Spectrum_gather" structure.  See db_ratio.h.  Note
that n1gather defines the storage increment between each frequency
vector.  This is fixed for each channel.  
	chan[0] = chanid[0] or 2
	chan[1] = chanid[1] or 3
	chan[2] = chanid[2] or 4
	chan[3] = 3component total power spectrum.
Each of these is allocced of size sg.nfreq*sg.n1gather, but each is
filled to a level defined by the int nstation.  i.e. on return
nstation[2], for example, will contain the number of live components
found for chanid[2] and/or chanid[5].

Author:  Gary Pavlis
Written:  November 1994 
*/

Spectrum_gather data_gather(Spectrum *spec,
	int nspec,
	char *chanids[6],
	Spectrum threeCspec[],
	int *nthreeCspec)
{
	Spectrum_gather sg;
        int nstations;
        int i,j,k,l;
        /* test variables */
	char laststa[10];
        int nfreq;
        double df,freq0;
	int first_chan;

	int index[6];  /* channel index array */ 
	int threeC_index[3];  /* used for assembling 3component power spectrum*/

        /* Scan for largest value of nfreq in group and use this to set 
        size */
        nfreq = 0;
        for(i=0;i<nspec;++i) nfreq = MAX(spec[i].nfreq,nfreq);

	sg.nfreq = nfreq;

        /* Establish the number of actual stations in this group */
        nstations = 0;
        for(i=1;i<nspec;++i) if(strcmp(spec[i].sta,spec[i-1].sta))++nstations;

	/* This is the size of the leading dimension for each 2d array*/
	++nstations;
	sg.n1gather = nstations;

        /* Alloc space using the size just determined */
        for(i=0;i<4;++i)
        {
                sg.chan[i] = (float *) calloc(nfreq*nstations,sizeof(float));
                if(sg.chan[i] == NULL)
		{
			sg.nfreq = -1;
			return(sg);
		}
        }
        /* this won't always work, but it should unless the user really
        screws up and mixes grossly different data types with different
        sample rates.  We scan for the first entry with spec.nfreq == nfreq.
        We assume df and freq0 from this channel are representative.  This
        should always be true with dbspectra because it only trims windows
        down never up.  Problems are possible if one appends data from multiple
        runs, but you can't correct for all levels of stupidity. */
        for(i=0;i<nspec;++i)
        {
                if(spec[i].nfreq == nfreq)
                {
                        df = spec[i].df;
                        freq0 = spec[i].freq0;
                        break;
                }
        }

        strcpy(laststa,spec[0].sta);
	/* The following algorithm only works if spec array is ordered
	by station-channel.  We scan to get valid 3c samples when possible.
	and add to chan arrays only when something is found */

	for(j=0;j<6;++j)index[j] = -1;	
	for(j=0;j<4;++j) sg.nstation[j] = 0;
	for(i=0;i<nspec;++i) spec[i].index = -1;

        for(i=0,first_chan=0,*nthreeCspec=0; i<nspec; ++i)
        {
		/* Found earlier logic error that this was needed to
		not drop the last channel in the group */
		if((i+1)==nspec)
		{
			for(j=0;j<6;++j)
			{
				if(!strcmp(chanids[j],spec[i].chan))
				{
					index[j] = i;
					break;
				}
			}
		}	
		if(strcmp(laststa,spec[i].sta) || ( (i+1) == nspec) )
		{
		/* Land here when station name changes.  This triggers
		3-component assembly processing for this station.
		First scan for inconsistent parameters, and mark them.
		when they occur for any channel. */
			for(j=first_chan;j<i;++j)
			{
				if( (spec[j].nfreq != nfreq ) 
                		|| ( fabs( (df-spec[j].df)/df) > FLOAT_TEST_RATIO)
                		|| ( fabs( (freq0-spec[j].freq0)/freq0) > FLOAT_TEST_RATIO) )
				{
					for(k=0;k<6;++k)
					{ 
						if(index[k]==j) 
						{
							index[k]= (-1);
							break;
						}
					}

				}
			}
			/* Now copy appropriate spectra to sg.chan arrays.
			Auxiliary channels are copied if primary is absent
			(defined by index set to -1). The for loop makes this
			shorter, but somewhat more obscure */
			for(l=0;l<3;++l) threeC_index[l] = -1;
			for(k=0;k<3;++k)
			{
				int kk,kaux,iused;
				kaux = k + 3;

				if(index[k] != -1) 
					kk = k;
				else
				{
					if(index[kaux] != -1)
						kk = kaux;
					else
						continue;
				}
				threeC_index[k] = index[kk];

				iused = index[kk];
				for(j=0,l=sg.nstation[k];j<nfreq;
				   ++j,l+=sg.n1gather)
					sg.chan[k][l] = spec[iused].sscale
							* spec[iused].spec[j];
				spec[iused].index = kk%3;
				++sg.nstation[k];
			}
			/* Now we have to deal with the three component
			problem.  To assemble the 3c power get messy 
			because of high-low gain (primary-auxiliary) 
			channels.  Alg is:  if we have all 3cs, alloc new
			row for threeCspec structure, alloc spec array, compute
			3c spec, copy to sg */
			if( (threeC_index[0] >= 0) &&  (threeC_index[1] >= 0) 
							&& (threeC_index[2] >= 0))
			{
				/* This routine is above.  Note it allocs space
				and copies most, but not all parameters from the 
				pattern Spectrum passed through arg1 */

				threeCspec[*nthreeCspec] = Spectrum_copy(spec[threeC_index[0]]); 

				/* Beautifully obscure pointer code.  C you gotta 
				love it ! */

				for(l=0;l<nfreq;++l) threeCspec[*nthreeCspec].spec[l] = 0.0;
				for(k=0;k<3;++k)
					for(l=0;l<nfreq;++l)
						threeCspec[*nthreeCspec].spec[l] 
						  += (spec[threeC_index[k]].sscale
						     * spec[threeC_index[k]].spec[l]);
				for(j=0,l=sg.nstation[3];j<nfreq;++j,l+=sg.n1gather)
					sg.chan[3][l] = (threeCspec[*nthreeCspec].spec[j]);
				++(sg.nstation[3]);
				++(*nthreeCspec);
			}
			
			
			/* These are cleanup increments before closing
			this block */

			first_chan = i;
			for(j=0;j<6;++j)index[j] = -1;
		}
		/* Always do this because even after gather procesing
		for each station, we need to reset and handle the current
		row. */

		for(j=0;j<6;++j)
		{
			if(!strcmp(chanids[j],spec[i].chan))
			{
				index[j] = i;
				break;
			}
		}	
		strcpy(laststa,spec[i].sta);
	}
        return(sg);
}

