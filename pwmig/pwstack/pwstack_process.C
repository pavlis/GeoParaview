#include <sunmath.h>
#include <sunperf.h>
#include "coords.h"
#include "pwstack.h"

void dzero(int n, double *d, int inc)
{
	int i, ii;
	for(i=0,ii=0;i<n;++i,ii+=inc) d[ii]=0.0;
}
void vscal(int n, double *w, int incw, double *y, int incy)
{
	int i, iw, iy;
	for(i=0,iw=0,iy=0;i<n;++i,iw+=incw,iy+=incy)
		y[iy] = w[iw]*y[iy];
}
void vadd(int n, double *x,int  incx, double *y,int  incy)
{
	int i, ix, iy;
	for(i=0,ix=0,iy=0;i<n;++i,ix+=incx,iy+=incy)
		y[iy] = x[ix] + y[iy];
}
/*
Main processing function for this program.  Takes an input 
data ensemble and produce a complete suite of plane-wave stacks
defined by the Rectangualr_Slowness_Grid object.  
This function uses and enhancement from the original Neal and Pavlis
and Poppeliers and Pavlis papers.  It allows the aperture to 
be tim variable.  This is generalized, but the normal expectation
is that the apeture would grow wider with time to compensate 
somewhat for diffraction.

Arguments
	indata - input raw data ensemble (see object definition
		the enscapsulates all this requires)
	ugrid - object to define slowness grid for stacking
	mute - mute applied to data before stacking
	stackmute - mute applied to data after stacking
		(This stack is aligned relative to latest
		mute time of raw data in stack.)
	lat0, lon0 - pseudostation grid point (IN RADIANS)
	ux0, uy0 - slowness vector of input data (stacking
		is relative to this vector)
	tstart:tend - define time period for output stack
		The routine pretty much assumes relative timing
		so this is normally time wrt the start of 
		each trace in the ensemble.
	aperture - defines variable aperture stack weighting
		(see above)
	mdlist - defines metadata to copy from the raw data
		ensemble to each stack output
	pfsho - output stream handle.  Data are pushed to 
		this routine for output.
*/
int pwstack_ensemble(Three_Component_Ensemble& indata,
	Rectangular_Slowness_Grid& ugrid,
	Top_Mute& mute,
	Top_Mute& stackmute,
	double lat0,
	double lon0,
	double ux0,
	double uy0,
	double tstart,
	double tend,
	Depth_Dependent_Aperture& aperture,
	list<Metadata_typedef>& mdlist,
	Pfstream_handle *pfsho)
{
	int i,j;
	vector<Three_Component_Seismogram>::iterator iv,ov;
	int nsta = indata.nsta;   // used so much having this alias is useful
	// This computes the output gather size.  It assumes all data have
	// a common sample rate and we can just grab the first one in the
	// list.  It also quietly assumes a relative time base 
	// so all times are computed relative to the start of
	// each trace.  Caller should guarantee this.
	//
	double dt=indata.tcse[0].dt;
	int nsin = indata.tcse[0].ns;


	Three_Component_Stack *stackout;

	int istart = nint(tstart/dt);
	int iend = nint(tend/dt);
	int nsout = iend-istart+1;
	int ismute = nint(mute.t1/dt);
	int ismute_this, ie_this;

	if(istart>=indata.tcse[0].ns)
		elog_die(0,(char *)"Irreconcilable window request:  Requested stack time window = %lf to %lf\nThis is outside range of input data\n",
			tstart,tend);

	/* Apply front end mutes to all traces */
	apply_top_mute_3c_ensemble(indata,mute);

	/* We need dnorth, deast vectors to compute moveout sensibly
	for this program.  Since we use them repeatedly we will
	extract them once from the gather.*/
	double *dnorth=new double[nsta];
	double *deast=new double[nsta];
	double *elev = new double[nsta]; // need elevations to compute geometric statics
	for(i=0,iv=indata.tcse.begin();iv!=indata.tcse.end();++iv,++i)
	{
		double lat,lon;
		int ierr;
		try{
			lat = (*iv).x[0].md.get_double("lat");
			lon = (*iv).x[0].md.get_double("lon");
			// assume metadata store these in degrees so we have to convert
			lat = rad(lat);
			lon = rad(lon);
			geographic_to_dne(lat0, lon0, lat, lon, dnorth+i, deast+i);
			elev[i]=(*iv).x[0].md.get_double("elev");
		} catch (Metadata_error merr)
		{
			throw merr;
		}
	}
	//
	// We want to make sure that all stations are in standard coordinates
	//
	for(iv=indata.tcse.begin();iv!=indata.tcse.end();++iv)
	{
		(*iv).rotate_to_standard();
	}
	// 
	// the weights become are a nsta by nsamp matrix to allow 
	// variable length apertures.  This algorithm assumes the
	// aperture object input is an nsta by 
	// Should use a matrix class for this, but it is easy enough
	// for this simple case to just code inline.

	double **weights = new double*[nsta];
	for(i=0;i<nsta;++i) weights[i]=new double[nsout];
	int stack_count=0;  //maximum number of traces not zeroed
	double *work = new double[nsta];
	for(i=0;i<nsout;++i)
	{
		int nused;
		
		nused=compute_pseudostation_weights(nsta, dnorth,deast,
		    aperture.get_aperture(tstart+dt*(double)i),
		    aperture.get_cutoff(tstart+dt*(double)i),work);
		for(j=0;j<nsta;++i) weights[j][i]=work[j];
		stack_count=max(nused,stack_count);
	}
	///
	// ERROR RETURN THAT NEEDS TO BE HANDLED GRACEFULLY
	// I don't throw an exception because this should not be viewed as
	// an exception.  It is a case that has to be handled gracefully.
	// it will happen often at the edges of arrays
	//
	if(stack_count<=0) return(-1);
	double *moveout=new double[nsta];
	double *stack[3];
	for(i=0;i<3;++i)stack[i]=new double[nsout];
	double *stack_weight=new double[nsout];
	double *twork=new double[nsout];

	double avg_elev,sum_wgt;  // computed now as weighted sum of station elevations
	for(i=0,avg_elev=0.0,sum_wgt=0.0;i<nsout;++i)
	{
		avg_elev += weights[i][0]*elev[i];
		sum_wgt += weights[i][0];
	}
	avg_elev /= sum_wgt;

	//
	// Loop over slowness grid range storing results in new output ensemble
	//

	for(i=0;i<ugrid.nux;++i)
	{
	    for(j=0;j<ugrid.nuy;++j)
	    {
		double ux,uy,dux,duy;
		int iend_this;
		int ismin;

		for(i=0;i<3;++i) dzero(nsout,stack[i],1);
		dzero(nsout,stack_weight,1);

		ux = ugrid.ux(i);
		uy = ugrid.uy(j);	
	
		/* The input gather is assumed prealigned with the slowness
		vector defined by raygath->ux,uy.  We use relative moveouts
		from this base moveout for the actual stacks */
		dux = ux - ux0;
		duy = uy - uy0;
		// moveout computed here and used below assumes the
		// data are aligned on the P arrival
		compute_moveout(nsta,deast,dnorth,dux,duy,moveout);
		for(i=0,iv=indata.tcse.begin(),ismute_this=0,iend_this=0,ismin=nsout;
			iv<indata.tcse.end();++iv,++i)
		{
			int is0, is, ismin, ns_to_copy; 
			is0 = istart + moveout[i];
			// need this to define latest front mute time
			ismute_this = max(ismute+nint(moveout[i]),ismute_this);
			// This computes starting position for the input trace  
			if(is0<0)
			{ 
				is=0;
				ismin = 0;
				ns_to_copy = min(nsin+is0,nsout);
				iend_this = max(iend_this,ns_to_copy-1);
			}
			else
			{
				is = is0;
				ismin = min(is,ismin);
				ns_to_copy = min(nsin-is0,nsout-is0);
				iend_this = max(iend_this,ns_to_copy+is0-1);
			}

			for(j=0;j<3;++j)
			{
				dzero(nsout,twork,1);
				dzero(nsout,work,1);
				// ugly expression below is an is offset from
				// start of vector component x[j] of 3c trace
				dcopy(ns_to_copy,((*iv).x[j].s)+is,1,
					twork,1);
				vscal(iend_this,weights[i],1,twork,1);
				vadd(iend_this,twork,1,stack[j],1);
				// Need to accumulate the sum of weights to normalize
				if(j==0) vadd(iend_this,weights[i],1,stack_weight,1);
			}
		}
		// normalize the stack.
		// not so trivial because of machine zero issue
		// Use a hard limit based on ismax and iend_this and then 
		// assume mutes kill possible transients in mute zone
		for(j=0;j<3;++j)
			for(i=ismin;i<=iend_this;++i)
			{
				if(stack_weight[i]>0)
					stack[j][i]/=stack_weight[i];
			}
		// this syntax might be wrong.  Need memory for a special 
		// 3c trace object inherited from Three_Component_Seismogram
		// this object add a set of metadata needed to define it
		stackout = new Three_Component_Stack(nsout);
		for(j=0;j<3;++j) dcopy(nsout,stack[j],1,stackout->x[j].s,1);
		copy_selected_metadata(indata.tcse[0].x[0].md,stackout->md,mdlist);
		//MISSING
		// somewhere here we need to load metadata for the stack object.
		// What needs to be loaded and how requires some thought.
		// putting this aside to work on it in background
		stackout->md.put_metadata("ux",ux);
		stackout->md.put_metadata("uy",uy);
		stackout->md.put_metadata("dux",dux);
		stackout->md.put_metadata("duy",duy);
		// may want to output a static here, but it is probably better to 
		// just keep a good estimate of elevation and deal with this in the
		// migration algorithm. 
		stackout->md.put_metadata("elev",avg_elev);
		for(j=0;j<3;++j) apply_top_mute(stackout->x[j],stackmute);

		// place holder for now.  Dependent on final mpi design
		// concept is that each stack is pushed to output stream as
		// individual 3component trace objects.
		save_stack(stackout,pfsho);
	    }
	}
	for(i=0;i<nsta;++i) delete [] weights[i];
	delete [] weights;
	delete [] dnorth;
	delete [] deast;
	delete [] moveout;
	delete [] stack_weight;
	delete [] work;
	delete [] twork;
	for(i=0;i<nsta;++i) delete[]stack[i];
	delete [] elev;
	return(0);
}
	
