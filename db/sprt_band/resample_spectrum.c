/*This routine resamples the spectrum stored in "s" according to 
spectral parameters stored in pattern.  This routine always returns a 
revised copy of s with nfreq, freq0, and df set to pattern. s.spec is
also replaced by a pointer to a new float array containing the resampled
spectrum.  Note we also free the old s.spec before doing this.

Author:  Gary Pavlis
Written: Nov 1994
*/

#include <malloc.h>
#include "stock.h"
#include "db.h"
#include "spectrum.h"

Spectrum resample_spectrum(s,pattern)
Spectrum s,pattern;
{
	float *work;  /* alloced array here to hold resampled spectrum*/
	int navg2;  /* half width to average around for decimation */
	int navg;  /* 2*navg2 + 1 */
	
	int i,j,is,k;
	double freq;

	work = (float *) calloc(pattern.nfreq,sizeof(float));
	if(work == NULL)
	{
		fprintf(stderr,"Fatal (resample_spectrum):  Cannot alloc %d floats\n",
			pattern.nfreq);
		exit(-1);
	}

	if(s.df > pattern.df)
	{
		fprintf(stderr,"Upsampling not implemented\n");
		exit(-1);
	}
	else if(s.df == pattern.df)
		return(s);
	else
	{
	/* this is block for decimation */
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

		for(i=0;i<pattern.nfreq;++i)
		{
			freq = (pattern.df)*((double)i);
			is = nint(freq/s.df);
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
	}
	free(s.spec);
	s.nfreq = pattern.nfreq;
	s.df = pattern.df;
	s.spec = work;
	return(s);
}
	
			

