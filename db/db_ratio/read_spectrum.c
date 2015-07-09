#include <stdlib.h>
#include "db.h"
#include "db_ratio.h"
Spectrum read_spectrum(Dbptr db)
{
	Spectrum s;
	char *fname;
	int nchar_fname;
	FILE *fp;
	Dbvalue dbv;
	int foff;

	if(dbgetv(db,0,
		"sta",s.sta,
		"psdisc.chan",s.chan,
		"sname",s.sname,
		"psdisc.phase",s.phase,
		"pstype",s.pstype,
		"datatype",s.datatype,
		"dir",s.dir,
		"dfile",s.dfile,
		"foff",&foff,
		"arid",&s.arid,
		"chanid",&s.chanid,
		"jdate",&s.jdate,
		"time",&s.time,
		"endtime",&s.endtime,
		"freq0",&s.freq0,
		"nfreq",&s.nfreq,
		"df",&s.df,
		"rayleigh",&s.rayleigh,
		"tbp",&s.tbp,
		"auth",s.auth,
		"sscale",&s.sscale,0) == dbINVALID)
	{
		s.spec = NULL;
		return(s);
	}

	/* Now we need to read the spectral file */
	s.spec = (float *) calloc(s.nfreq,sizeof(float));
	if(s.spec == NULL) return(s);

	if(s.dir[0] == '/')
	{
		nchar_fname = strlen(s.dfile)+strlen(s.dir)+2;
		if((fname = (char *)malloc(nchar_fname)) == NULL) 
		{
			free(s.spec);
			s.spec = NULL;
			return(s); 
		}

		strcpy(fname,s.dir);
	}
	else
	{
		dbquery(db,dbTABLE_DIRNAME,&dbv);
		if((fname=(char *)malloc(strlen(dbv.t)+strlen(s.dir)
					+strlen(s.dfile)+4)) == NULL)
                {
                        free(s.spec);
                        s.spec = NULL;
                        return(s);   
                }
		strcpy(fname,dbv.t);
		strcat(fname,"/");
		strcat(fname,s.dir);
	}
	strcat(fname,"/");
	strcat(fname,s.dfile);
	
	if((fp = fopen(fname,"r")) == NULL)
	{
		free(s.spec);
		free(fname);
		s.spec = NULL;
		return(s);
	}
	if(foff > 0) fseek(fp, (long) foff, SEEK_SET);
	
	if(fread((void *)s.spec,sizeof(float),s.nfreq,fp) != s.nfreq)
        {
                free(s.spec);
                s.spec = NULL;
        }
	fclose(fp);
        free(fname);
        return(s);
}
