#include <stdio.h>
#include "stock.h"
#include "pwmig.h"
/* Small companion function to one immediately below.  It times shifts
trace by shift samples inserting 0s at front or end depending on sign
of shift.   n is the length of trace.
*/
void trace_sample_shift(int n, float *trace,int shift)
{
	int i;
	/* Note we intentionally do nothing if shift=0 */
	if(shift>0)
	{
		for(i=n-1;i>=shift:--i)
			trace[i] = trace[i-shift];
		for(i=0;i<shift;++i) trace[i]=0.0;
	}
	else if(shift<0)
	{
		for(i=0;i<n+shift;++i)
			trace[i] = trace[i-shift];
		for(i=n+shift;i<n;++i) trace[i] = 0.0;
	}
}
/* this function reads a gather of raw traces from a evid grouped 
complicated view used within this program.  Barring further changes
the view is formed the the join:
wfprocess->psdecon->site
where wfprocess is the core processed waveform table used 
as a foundation for this and a sequence of related programs.  
Several assumptions about this view are made:
1.  the view is subsetted to a unique "gridname" and "method"
2.  view is sorted by evid:ix1:ix2:sta:chan and grouped by
evid:ix1:ix2 
3.  the db pointer passed to the routine is the output of 
dbgroup applied to 2.

This processing means, in words, that the range of rows passed
are for one event at one point in a grid of pseudostations.  

The function returns a pointer to a gather object containing all
the deconvolved seismograms from point ix1, ix2 and evid.  

The function sorts out channels in a not very robust fashion under
an assumption this function is only called after data has passed
a lot of steps.  Specifically, three things must be satisfied or
this routine will most likely return nothing or (worse) garbage:
1.  Every station must have all three components and only three 
components defined
2.  The channel codes must be defined so the sort yields an order
in standard coordinates (E,N,Z for example do this as will 
X1, X2, X3)  
3.  All traces must have the same sample rate and trace length in
the gather.  

If t0 in the input table is a nonzero value, the data are shifted
appropriately leaving zeros at the start or end of the trace.
t0 is a way to apply post decon statics and may never be used, but
was added in anticipation of it being a useful entity in future 
versions.  

Author:  G Pavlis
Written:  June 2000
Modification:  Aug 2000
Changed gather structure to carry along a base slowness vector 
that defines the moveout of the input data.  Original was wrong
because I hadn't recognized at the time the imaging condition 
depends upon the slowness vector relative to this base slowness
vector.  This added pf as a parameter.  From pf we use the
array name and access the arrival table to get a slowness vector
estimate.  
*/
#define KMPERDEG 111.31996
PWMIGraw_trace_gather *dbread_decon_data(Dbptr dbgrp, Pf *pf)
{
	Dbptr db_bundle;
	Dbptr dbs;
	int is, ie;
	int evid,nsta;
	PWMIGraw_trace_gather *g;
	int nsamp,nsamptest,samprate,sampratetest;
	int i,ista;
	char sta[10];
	char chan[10];
	char gridname[10];
	int ix1,ix2;
	double lat,lon;
	double dnorth,deast,t0;
	char datatype[4];
	double foff;
	char filename[128];
	FILE *fp;
	double t0;
	int it0;
	char ssexpr[50],ssexpr2[50];
	Tbl *proctbl;
	int narvrecs;
	double slo,azimuth;
	

	/* Get the bundle pointer from the group db pointer */
	if(dbgetv(dbgrp,0,"evid",&evid,"bundle",&db_bundle,0) == dbINVALID)
	{
		elog_complain(0,"PWMIGtrace_gather: dbgetv error reading row %d of evid group view\nData for this event will be skipped\n",
			dbgrp.record);
		return(NULL);
	}
	dbget_range(db_bundle,&is,&ie);

	allot(PWMIGraw_trace_gather *,g,1);
	nsta = (ie - is + 1)/3;

	/* Check that the number of traces in this group is exactly 
	divisible by three.  If it isn't, something is seriously 
	wrong and we will not attempt to recover but issue a complain
	and return a NULL */
	if( (nsta*3) != (ie-is+1) )
	{
		free(g);
		elog_complain(0,"Missing data for evid %d\nEvery station must have all three components operational\nNumber of traces for this evid = %d which is not divisible by 3\n",
			evid,ie-is+1);
		return(NULL);
	}
	/* Before we fill the rest of g we need to get the slowness vector
	from the arrival table.  This is subject to many potential problems
	in database setup, so we do it before we do the main work later.
	*/
	array_name = pfget_string(pf,"array_name");
	sprintf(ssexpr,"dbsubset evid==%d",evid);
	sprintf(ssexpr2,"dbsubset sta=~/%s/",array_name);
	proctbl = strtbl("dbopen event",
			ssexpr,
			"dbjoin origin",
			"dbjoin assoc",
			"dbjoin arrival",
			ssexpr2,0);
	dbs = dbprocess(dbgrp,proctbl,0);
	freetbl(proctbl,0);
	dbquery(dbs,dbRECORD_COUNT,&arvrecs);
	if(arvrecs<=0)
	{
		elog_complain(0,"Lookup for array name %s in arrival table failed looking up slowness vector for evid %d\nCannot process this event\n",
			array_name,evid);
		free(g);
		return(NULL);
	}
	if(arvrecs>1) 
		elog_notify(0,"WARNING:  %d matching slowness vectors found in arrival table for evid %d with sta/array name = %s\nThere should only be one.  Running dbverify is recommended\n",
			arvrecs,evid,array_name);
	dbs.record = 0;
	dbgetv(dbs,0,"slo",&slo,"azimuth",&azimuth,0);
	/* slo is in units of s/deg but internally we use s/km so we
	have to convert */
	slo /= KMPERDEG;
	azimuth = rad(azimuth);  
	g->ux = slo*sin(azimuth);
	g->uy = slo*cos(azimuth);

	g->nsta = nsta;
	g->evid = evid;

	allot(char **,g->sta,nsta);
	allot(double *,g->lat,nsta);
	allot(double *,g->lon,nsta);
	allot(float **,g->x,nsta);
	allot(float **,g->y,nsta);
	allot(float **,g->z,nsta);

	/* these data are presumed to be in standard coordinates
	so we set the transformation matrix to an identity */
	for(i=0;i<9;++i)g->U[i] = 0.0;
	g->U[0] = 1.0;
	g->U[3] = 1.0;
	g->U[6] = 1.0;

	/* We read through the data being inflexible in what we 
	accept.*/

	for(db_bundle.record=is,i=0,ista=0;db_bundle.record<ie;
			++db_bundle.record,++i)
	{
		dbgetv(db_bundle,0,
			"sta",sta,
			"chan",chan,
			"nsamp",&nsamp,
			"samprate",&samprate,
			"gridname",(g->gridname),
			"ix1",&ix1,
			"ix2",&ix2,
			"lat",lat+i,
			"lon",lon+i,
			"deast",&deast,
			"dnorth",&dnorth,
			"t0",&t0,
			"datatype",datatype,
			"foff",foff,0);
		if(i==0)
		{
			sampratetest=samprate;
			g->si = 1.0/samprate;
			nsamptest = nsamp;
			g->nsamp = nsamp;
			g->ix1 = ix1;
			g->ix2 = ix2;
		}
		/* Probably should trap ix1 and ix2 mismatch also, but
		if the database is correctly sorted this should not
		happen.  Potential problem only if this routine is
		recycled in a different context*/
		if((samprate != sampratetest) 
			|| (nsamp != nsamptest) )
		{
			elog_complain(0,"Trace mismatch reading data for evid %d in row %d\nSample rate = %lf expected %lf and/or nsamp = %d expected %d\nData for this event will be skipped\n",
				evid,i,samprate,sampratetest,nsamp,nsamptest);
			/* This is a trick to be sure we don't try to 
			free null pointers */
			g->nsta = ista;
			destroy_PWMIGraw_trace_gather(g);
			return(NULL);
		}
		if(dbextfile(db_bundle,0,filename) != 0)
		{
			elog_complain(0,"File %s does not exist\nProcessing aborted for evid %d\n",
				filename,evid);
			g->nsta = ista;
			destroy_PWMIGraw_trace_gather(g);
			return(NULL);
		}
		if(fp=fopen(filename,"r") == NULL)
		{
		{
			elog_complain(0,"Cannot open file %s for reading\nProcessing aborted for evid %d\n",
				filename,evid);
			g->nsta = ista;
			destroy_PWMIGraw_trace_gather(g);
			return(NULL);
		}
		fseek(fp,foff,SEEK_SET);
		it0 = nint(t0*samprate);
		ichan = i%3;
		switch(ichan)
		{
		case (0):
			g->sta[ista] = strdup(sta);
			/* Hate this backward logic, but necessary */
			if( !((chan[0]=='E') || (chan[0]=='X')) )
				elog_notify(0,"Warning(dbread_decon_data):  potential channel mismatch on record %d = %s\nChannel inserted as EW component\n",
					db_bundle.record,chan);
			allot(float *,g->x[ista],nsamp);
			if(fread(g->x[ista],sizeof(float),nsamp) != nsamp)
				elog_notify(0,"Warning(dbread_decon_data): read error on file %s for sta:chan=%s:%s\nGather may be corrupted\n",
					filename,sta,chan);
			if(it0) trace_sample_shift(nsamp,g->x[ista],it0);
			break;
		case (1):
			/* Near identical logic for NS component */
			if( !((chan[0]=='N') || (chan[0]=='Y')) )
				elog_notify(0,"Warning(dbread_decon_data):  potential channel mismatch on record %d = %s\nChannel inserted as a NS component\n",
					db_bundle.record,chan);
			allot(float *,g->y[ista],nsamp);
			if(fread(g->y[ista],sizeof(float),nsamp) != nsamp)
				elog_notify(0,"Warning(dbread_decon_data): read error on file %s for sta:chan=%s:%s\nGather may be corrupted\n",
					filename,sta,chan);
			if(it0) trace_sample_shift(nsamp,g->y[ista],it0);
			break;
		case (2):
			/* and for z */
			if( !((chan[0]=='Z') || (chan[0]=='U')) )
				elog_notify(0,"Warning(dbread_decon_data):  potential channel mismatch on record %d = %s\nChannel inserted as a Vertical component\n",
					db_bundle.record,chan);
			allot(float *,g->z[ista],nsamp);
			if(fread(g->z[ista],sizeof(float),nsamp) != nsamp)
				elog_notify(0,"Warning(dbread_decon_data): read error on file %s for sta:chan=%s:%s\nGather may be corrupted\n",
					filename,sta,chan);			
			if(it0) trace_sample_shift(nsamp,g->z[ista],it0);
			++ista;
		}
		fclose(fp);
	}
	return(g);
}

void destroy_PWMIGraw_trace_gather(PWMIGraw_trace_gather *g)
{
	int i;
	for(i=0;i<g->nsta;++i)
	{
		free(x[i]);
		free(y[i]);
		free(z[i]);
		free(sta[i]);
	}
	free(g->sta);
	free(g->lat);
	free(g->lon);
	free(g->x);
	free(g->y);
	free(g->z);
	free(g);
}
/* This short function simultaneously applies a vector of moveout 
corrections and computes a weighted stack of a gather of seismograms.
Sybolically, it returns

<- sum i s[i][0:nsamp-1] w[i]

Arguments:
	s - matrix of seismograms to be stacked with s[i] being a pointer
		to the data for the ith station
	nsta - number of rows in s = number of stations
	nsamp - number of columns in s = trace length 
		(Note the returned trace is of length nsamp)
	si - sample interval in seconds
	w - parallel vector to s of weights to be applied to each trace.
	moveout - vector of moveouts (in units of seconds ) to apply to 
		each trace.  + shifts lead to zeros at the front of the
		trace and - shifts lead to zeros at the end of the trace.

The function uses the BLAS function saxpy and returns a pointer to 
a vector of length nsamp containing the stack.  Note the stack will
be irregular at the front and end over a range that depends on the
size of the elements of the moveout vector. 

NOTE: the returned vector is alloced internally here so the calling
function must take care to avoid a memory leak.

Author:  G Pavlis
Written:  June 2000
*/
float *weighted_stack_with_moveout(float **s, int nsta, int nsamp,
	double si, double *w, double *moveout)
{
	int i;
	float *stack;
	int im;
	int ntostack;

	allot(float *,stack, nsamp);
	for(i=0;i<nsamp;++i) stack[i]=0.0;
	for(i=0;i<nsta;++i)
	{
		im = nint(moveout[i]/si);
		/* This is done for efficiency with w[i]=0.0 a
		likely possibility in some cases */
		if( (abs(im)>nsamp) || (w[i] == 0.0) ) continue;
		if(im>=0)
		{
			ntostack = nsamp - im;
			saxpy(ntostack,(float)w[i],s[i],1,stack+im,1);
		}
		else
		{
			ntostack = nsamp + im;
			saxpy(ntostack,(float)w[i],s[i]-im,1,stack,1);
		}
	}
	return(stack);
}	
/* Apply a front end mute to a single trace.  

Arguments:
	s = seismogram to apply mutes to 
	i1 and i2 define the mute zone in samples (see below)

The mute zone is defined in samples.  The trace is zeroed from s[0]
THROUGH s[i1].  A linear ramp of weights is applied between i1 and i2 
such that s[i1] will be zero and s[i2] will be unaltered.  The 
routine does not check if i1 and i2 are valid and seg faults are 
likely if s[0:i2]  is outside the bounds of the vector s.  It
does return immediately if i1 is negative.  If i2 < i1 no taper
is performed.

Author:  G Pavlis
Written: June 2000
*/
void apply_fend_mute(float *s,int i1, int i2)
{
	int i;
	float w_per_sample,weight;
	if(i1<0) return;
	for(i=0;i<i1;++i) s[i] = 0.0;
	if(i2<=i1) 
	{
		s[i1]=0.0;
		return;
	}
	w_per_sample = 1.0/((float)(i2-i1));
	for(i=i1;i<i2;++i)
	{
		weight = ((float)(i-i1))*w_per_sample;
		s[i]*=weight;
	}
}

/* Higher level version of above.  It applies the front end mute function
defined above to a full gather of traces stored in the input matrices
x, y, and z.  This is not fancy as fixed mutes are applied across the
whole gather.  It is appropriate for p to s conversion imaging when
applied to deconvolved traces, but little else as normal front end 
mutes used, for example, in reflection processing are trace dependent.  

The parameters that define the mute zone are passed through a parameter
file object.

Arguments:
	nsta - number of stations
	nsamp - number of samples in each seismogram of the gather.
	x,y,z - matrices of traces to be muted.  These are nsta by
		nsamp C matrices that hold the 3 component traces
		to be muted.
	si - sample interval in seconds.
	pf - parameter file object defining mute zone parameters

Author:  G Pavlis
Written:  June 2000
*/
void mute_gather(int nsta, int nsamp, 
	float **x, float **y, float **z,
	double si, Pf *pf)
{
	int i;
	double t1, t2, dt;
	int i1, i2;

	t1 = pfget_double(pf,"mute_zero_length");
	dt = pfget_double(pf,"mute_taper_length");
	t2 = t1 + dt;
	i1 = nint(t1/si);
	i2 = nint(t2/si);
	if(i1>=nsamp)
	{
		elog_complain(0,"mute_gather:  requested mute zone length is %d samples while the input traces are only %d long\nMutes not applied\n",
			i1, nsamp);
		return;
	}
	if(i2>=nsamp)
	{
		elog_complain(0,"mute_gather:  requested end of tapered area of front end mute zone of %d samples is longer than the input trace length %d\nReset to end of trace\n",
			i2 = nsamp-1;
	}
	for(i=0;i<nsta;++i)
	{
		apply_fend_mute(x[i],i1,i2); 
		apply_fend_mute(y[i],i1,i2); 
		apply_fend_mute(z[i],i1,i2); 
	}
}
/* Companion to below.  Duplicates contents of input gather, gi, to
output gather, go*/
void copy_pwmigraw_gather(PWMIGraw_trace_gather *gi,PWMIGraw_trace_gather *go)
{
	int i;
	go->nsta = gi->nsta;
	go->nsamp = gi->nsamp;
	go->ix1 = gi->ix1;
	go->ix2 = gi->ix2;
	go->si = gi->si;
	go->ux = gi->ux;
	go->uy = gi->uy;
	for(i=0;i<9;++i)go->U[i] = gi->U[i];
	strcpy(go->gridname,gi->gridname);
	allot(char **,go->sta,nsta);
	allot(double *,go->lat,nsta);
	allot(double *,go->lon,nsta);
	allot(float **,go->x,nsta);
	allot(float **,go->y,nsta);
	allot(float **,go->z,nsta);
	for(i=0;i<nsta;++i) go->sta[i] = strdup(gi->sta[i]);
	for(i=0;i<nsta;++i) go->lat[i] = gi->lat[i];
	for(i=0;i<nsta;++i) go->lon[i] = gi->lon[i];
	for(i=0;i<nsta;++i) scopy(gi->nsamp,gi->x[i],1,go->x[i],1);
	for(i=0;i<nsta;++i) scopy(gi->nsamp,gi->y[i],1,go->y[i],1);
	for(i=0;i<nsta;++i) scopy(gi->nsamp,gi->z[i],1,go->z[i],1);
}
/* This function takes an input gather g and returns a copy of
it but with the traces rotated to ray coordinated defined by
the spherical coordinate input variable s.  A few details:
1.  The angles in s are defined from geographical coordinates
with x1 + east, x2 + north, x3 + up 
2.  angles are in radians
3.  The transformation matrix used internally is in fortran order

Author:  G Pavlis
Written: August 2000
*/
PWMIGraw_trace_gather *rotate_pwgather(PWMIGraw_trace_gather *g,
		Spherical_Coordinate s)
{
	PWMIGraw_trace_gather *newgath;
	double U[9];   /*3x3 coordinate transformation matrix in fortran form*/
	float *work;
	int i;

	allot(PWMIGraw_trace_gather *,newgath,1);
	allot(float *,work,gi->nsamp);
	copy_pwmigraw_gather(g,newgath);

	ray_coordinate_trans(s,U);
	for(i=0;i<9;++i) newgath->U[i] = U[i];
	/* This would have been simpler had I chosen to make the gather
	structure a matrix, but I hope this is clearer */
	for(i=0;i<g->nsta;++i)
	{
		scopy(g->nsamp,g->x[i],1,work,1);
		sscal(g->nsamp,(float)U[0],work,1);
		saxpy(g->nsamp,(float)U[1],g->y[i],1,work,1);
		saxpy(g->nsamp,(float)U[2],g->z[i],1,work,1);
		scopy(g->nsamp,work,1,newgath->x[i],1);

		scopy(g->nsamp,g->x[i],1,work,1);
		sscal(g->nsamp,(float)U[3],work,1);
		saxpy(g->nsamp,(float)U[4],g->y[i],1,work,1);
		saxpy(g->nsamp,(float)U[5],g->z[i],1,work,1);
		scopy(g->nsamp,work,1,newgath->y[i],1);

		scopy(g->nsamp,g->x[i],1,work,1);
		sscal(g->nsamp,(float)U[6],work,1);
		saxpy(g->nsamp,(float)U[7],g->y[i],1,work,1);
		saxpy(g->nsamp,(float)U[8],g->z[i],1,work,1);
		scopy(g->nsamp,work,1,newgath->z[i],1);
	}

	free(work);
}
