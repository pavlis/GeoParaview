#include <stdlib.h>
#include <string>
#include "stock.h"
#include "pf.h"
#include "pfstream.h"
#include "metadata.h"
#include "seispp.h"

using namespace std;

/*
 * Gets a Time_Series_Ensemble object from a pfstream.  gmdlist
 * is a list of names extracted from the input stream and placed into
 * the global metadata space for the ensemble.
 * This function is not very general but is a good illustration of
 * how we can load a generic object through a pfstream
 * Function returns NULL when the input is exhausted.  It throws 
 * a string exception containing a file name where i/o fails. 
 * It can also terminate internally with input errors on the pfstream.
 * This is not very elegant, but appropriate for planned usage.
 */
vector<Time_Series> *get_next_ensemble(Pfstream_handle *pfh, char *tag)
		throw(string)
{
	Pf *pfin;
	Pf_ensemble *pfe;
	int i;
	char *pfstr;
	string datasource;
	string dfile, dir;
	int foff;
	Tbl *samples_raw;
	int ns;

	//  This routine gets the data required to construct this
	//  from the input pfstream.  I returns everything 
	//  encapsulated in a single pf
	pfin = pfstream_get_next_ensemble(pfh);
	if(pfin==NULL) return(NULL);

	// We next parse the input pf encapsulated in a single
	// set of curly brackets into an ensemble
	pfe = pfget_Pf_ensemble(pfin,tag);
	if(pfe==NULL) elog_die(0,
		(char *)"get_next_ensemble:  input stream error\nParsing pf for name %s yielded nothing\n",
			tag);
	// For this type of ensemble groups are meaningless.  Just 
	// build a simple vector of the results and discard anything like
	// that if present.  We build a vector of 
	// Time_Series objects in sequence.  They are built and pushed 
	// at the bottom of this loop
	vector <Time_Series> *tse=new vector<Time_Series>;
	for(i=0;i<pfe->nmembers;++i)
	{
		double samprate;
		double t0;
		string stref;

		// Intentionally construct an empty Time_Series object
		Time_Series *ts=new Time_Series();
		
		pfstr = pf2string(pfe->pf[i]);
		try{
			ts->md.load_metadata((string)pfstr);
			free(pfstr);
		} 
		catch(int ierr)
		{
			elog_die(0,(char *)"get_next_ensemble:  pfcompile error on input\nInput stream has been corrupted\n");
		}
		try {
		// Name space here is frozen.  Not elegant, but a starter
		// We are extracting required parameters placed into the 
		// Time_Series object.  
			samprate = ts->md.get_double("samprate");
			ts->dt = 1.0/samprate;
			ts->t0 = ts->md.get_double("starttime");
			ts->ns = ts->md.get_int("nsamp");
			stref = ts->md.get_string("Time_Reference_Type");
			if(stref == "relative")
				ts->tref = relative;
			else
				ts->tref = absolute;
			//  This allows data to be transmitted inline
			//  through the pf or through external files
			datasource = ts->md.get_string("datasource");
			if(datasource == "external_files")
			{
				dir = ts->md.get_string("dir");
				dfile = ts->md.get_string("dfile");
				foff = ts->md.get_int("foff");
			}
			else
			{
				samples_raw = ts->md.get_list("time_series_samples");
				ns = maxtbl(samples_raw);
				for(int j=0;j<(ts->ns);++j)
				{
					char *ssr;
					ssr = (char *)gettbl(samples_raw,j);
					ts->s[j] = atof(ssr);
				}
			}
		} 
		catch(int ierr_mdparse)
		{
			elog_die(0,(char *)"get_next_ensemble:  error parsing required metadata\n");
		}
		// I do this here because of the try/catch block.
		// I want the function to throw an exception that
		// can be caught be the caller an handled in the 
		// event of i/o errors.  The above should only cause
		// problems with usage errors.
		if(datasource =="external_files")
		{
			FILE *fp;
			string fname=dir+"/"+dfile;
			if((fp=fopen(fname.c_str(),"r")) == NULL) throw(fname);
			if (foff>0)fseek(fp,(long)foff,SEEK_SET);
			ts->s = new double[ts->ns];
			if(fread((void *)(ts[i].s),sizeof(double),ts[i].ns,fp)
					!= (ts[i].ns) ) throw(fname);
			fclose(fp);
		}
		tse->push_back(*ts);
	}
	return(tse);
}
