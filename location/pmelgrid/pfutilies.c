/* This function takes a vector Pf *'s from is to ie and 
 * builds an output list of Arrival objects.  This should be
 * thought of like a database get routine that works from row is
 * to ie extracting information from the arrival table associated
 * with a given event with an assuption that the table is sorted
 * so that is to ie is associated with a common event.  
 * 
 * phases and stations are associative arrays for phase handles
 * and station object (defined in libgenloc).  This routine
 * is assumed to be called at a high level and little error
 * checking is done.  It always returns a valid Tbl handle, but
 * that Tbl could be empty.  Any problems in cracking a single
 * Pf item will lead to skipping that entry with a diagnostic
 * issued to elog.
 *
 * Author:  Gary L. Pavlis
 * Written:  October 2002
 */
Tbl *pfextract_arrivals(Pf **pfall, int is, int ie, 
	Arr *phases,Arr *stations)
{
	Arrival *a;
	Tbl *tout;
	int i;
	Pf *pf;

	tout=newtbl(0);
	for(i=is;i<=ie;++i)
	{
		char *sta,*phase_name;
		/* a convenient shorthand */
		pf=pfall[i];
		allot(Arrival *,a,1);
		sta = pfget_string(pf,"sta");
		if(sta==NULL)
		{
			register_error(0,"No sta parameter for ensemble member %d\nSkipped\n",
				i);
			free(a);
			continue;
		}
                a->sta = (Station *) getarr(stations,sta);
                if(a->sta == NULL)
                {
                        register_error(0,"Warning (read_arrivals):  Can't find coordin
ates for station %s\n%s phase arrival for this station skipped\n",
                                sta, phase_name);
                        free(a);
                        continue;
                }
		phase_name=pfget_string(pf,"phase");
                a->phase = (Phase_handle *) getarr(phases,phase_name);
                if(a->phase == NULL)
                {
                        register_error(0,"Warning (read_arrivals):  No phase handle fo
r phase name %s\nArrival for station %s skipped\n",
                                phase_name,sta);
                        free(a);
                        continue;
                }
		a->arid = pfget_int(pf,"arid");
		a->time = str2epoch(pfget_string(pf,"arrival.time"));
		a->deltat = pfget_double(pf,"deltim");
                /* Here set set deltat to default when this is required */
                if( (a->deltat) <= 0.0 ) a->deltat = (double)a->phase->deltat0;
                pushtbl(tout,a);
	}
	return(tout);
}
int pfstream_save_results(char *streamout,
		char *tag,
		int nevents,
		int *evid,
		Hypocenter *h,
		Tbl **ta,
		double sswrodf,
		)
	
