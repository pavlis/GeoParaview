/* This function builds an input database view from deconvolution
output.  The view is created as:
1.  the join:
wfprocess->psdecon->site
where wfprocess is the core processed waveform table used 
as a foundation for this and a sequence of related programs. 
This is done in a generalized way through the dbpwfjoin function.
2.  the view is subsetted to a unique "gridname" and "method"
3.  view is sorted by evid:ix1:ix2:sta:chan and grouped by
evid:ix1:ix2 

The function returns the output of 3.

Arguments:
	db - input database name
	gridname, method - the view is subsetting to have
		unique gridname and method values.  Gridname is
		the name of the pseudostation grid that defines
		the geometry.  method is the migration method.

Returns a pointer to the grouped table described above.  If
there are any errors in the processing the db.record field
is set to dbINVALID.

Author:  Gary Pavlis
Written:  August 2000
*/
#define DECON_ALGORITHM_FIELD ""cwdecon""  /* pneumonic (ha ha ) for
				converted wave deconvolution*/
Dbptr create_decon_input_view(Dbptr db, char *gridname, char *method)
{
	Dbptr dbva,dbvb;
	char sset[100];
	Tbl *sorkeys,*grpkeys;

	dbva = dbpwfjoin(db,DECON_ALGORITHM_FIELD);

	sprintf(sset,"(gridname=~/%s/)&&(method=~/%s)",
			gridname,method);
