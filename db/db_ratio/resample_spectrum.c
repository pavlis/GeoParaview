/*This routine resamples the spectrum stored in "s" according to 
spectral parameters stored in pattern.  This routine always returns a 
revised copy of s with nfreq, freq0, and df set to pattern. s.spec is
also replaced by a pointer to a new float array containing the resampled
spectrum.  Note we also free the old s.spec before doing this.

Author:  Gary Pavlis
Written: Nov 1994
Revised:  October 2001
Major enhancement.  Earlier version was very limited.  This one
should handle all variants correctly.  freq0 and nfreq are 
adjusted as required to match overlapping area to be resampled.
That is, the pattern and s spectra can overlap in a variety of
ways.  The result is returned with the overlap area only set.
This varies nfreq and freq0.
*/

#include <malloc.h>
#include <sunmath.h>
#include "stock.h"
#include "db.h"
#include "db_ratio.h"

Spectrum resample_spectrum(s,pattern)
Spectrum s,pattern;
{
	float *work;  /* alloced array here to hold resampled spectrum*/
	int navg2;  /* half width to average around for decimation */
	int navg;  /* 2*navg2 + 1 */
	
	int i,j,is,k;
	double freq,freq0;

	/* Do nothing when df matches */
	if(s.df == pattern.df) return(s);

	work = (float *) calloc(pattern.nfreq,sizeof(float));
	if(work == NULL)
	{
		fprintf(stderr,"Fatal (resample_spectrum):  Cannot alloc %d floats\n",
			pattern.nfreq);
		exit(-1);
	}

	if(s.df > pattern.df)
	{
	/* This is block for interpolation.  Uses a 
	simple linear interpolation between points */
		if(pattern.freq0<s.freq0)
		{
			freq0=s.freq0;
			fprintf(stderr,"Warning(resample_spectrum:  First frequency mismatch:  pattern freq0=%lf is less than requested %lf\nUsing freq0=%lf\n",
				pattern.freq0,s.freq0,s.freq0);
		}
		else
		{
			freq0=pattern.freq0;
		}
		for(i=0,freq=freq0;i<pattern.nfreq;
				++i,freq+=pattern.df)
		{
			float grad,freq_grid;
			is=(int)((freq-s.freq0)/s.df);
			if((is+1)>=s.nfreq)break;
			grad = (s.spec[is+1]-s.spec[is])/s.df;
			freq_grid = s.freq0+((double)is)*s.df;
			work[i] = s.spec[is]
				+ grad*(freq-freq_grid);
		}
		s.freq0=freq0;
		s.nfreq=i;
	}
	else
	{
	/* this is block for decimation */
		if(pattern.freq0<s.freq0)
		{
			freq0=s.freq0;
			fprintf(stderr,"Warning(resample_spectrum:  First frequency mismatch:  pattern freq0=%lf is less than requested %lf\nUsing freq0=%lf\n",
				pattern.freq0,s.freq0,s.freq0);
		}
		else
		{
			freq0=pattern.freq0;
		}
		navg2 = nint ( (pattern.df)/(s.df));
		navg2 /= 2;
		navg = 2*navg2 + 1;
		if(navg>s.nfreq)
		{
			fprintf(stderr,"fatal (resample_spectrum): illegal parameter combination\n");
			fprintf(stderr,"averaging width %d exceeds number of frequencies = %d\n",
				navg,s.nfreq);
			exit(-1);
		}
		for(i=0,freq=freq0;i<pattern.nfreq;
				++i,freq+=pattern.df)
		{
			is = nint((freq-s.freq0)/s.df);
			if(is>=s.nfreq) break;
			if(is<navg2)
				is = 0;
			else if( (is+navg2) >= s.nfreq)
				is = s.nfreq - navg - 1;
			else
				is -= navg2;
			work[i] = 0.0;
			for(j=is,k=0;k<navg;++j,++k) work[i] += s.spec[j];
			work[i] /= navg; 
		}
		s.freq0=freq0;
		s.nfreq = i;
	}
	s.df = pattern.df;
	free(s.spec);
	s.spec = work;
	return(s);
}
