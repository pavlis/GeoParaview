#include "pwstack.h"

void dzero(int n, double *d, inc)
{
	int i, ii;
	for(i=0,ii=0;i<n;++i,ii+=inc) d[ii]=0.0;
}
void vscal(int n, w, incw, y, incy)
{
	int i, iw, iy;
	for(i=0,iw=0,iy=0;i<n;++i,iw+=incw.iy+=incy)
		y[iy] = w[iw]*y[iy];
}
void vadd(int n, x, incx, y, incy)
{
	int i, ix, iy;
	for(i=0,ix=0,iy=0;i<n;++i,ix+=incx,iy+=incy)
		y[iy] = x[ix] + y[iy];
}
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
	list<Metadata_typedef *>mdlist&,
	Pfstream_handle *pfsho)
{
	int i,j;
	vector::iterator iv,ov;
	int nsta = indata.nsta;   // used so much having this alias is useful
	// This computes the output gather size.  It assumes all data have
	// a common sample rate and we can just grab the first one in the
	// list.  It also quietly assumes a relative time base 
	// so all times are computed relative to the start of
	// each trace.  Caller should guarantee this.
	//
	double dt=indata.tcsd[0].x[0].dt;
	int nsin = indata.tcsd[0].nsin;


	Three_Component_Stack *stackout;

	int istart = nint(tstart/dt);
	int iend = nint(tend/dt);
	int nsout = iend-istart+1;
	int ismute = nint(t1/dt);
	int ismute_this, ie_this;

	if(istart>=indata.tcsd[0].nsamp)
		elog_die(0,(char *)"Irreconcilable window request:  Requested stack time window = %lf to %lf\nThis is outside range of input data\n",
			tstart,tend);

	/* Apply front end mutes to all traces */
	apply_top_mute_3c_ensemble(indata,mute);

	/* We need dnorth, deast vectors to compute moveout sensibly
	for this program.  Since we use them repeatedly we will
	extract them once from the gather.  Assume all have the
	same dnorth, deast loaded for the same station so use x[0] */
	double dnorth[nsta];
	double deast[nsta];
	for(iv=indata.begin();iv<indata.end;++iv)
	{
		double lat,lon;
		int ierr;
		try{
			lat=indata.tcse[iv].x[0].md.get_double("lat");
			lon=indata.tcse[iv].x[0].md.get_double("lon");
			geographic_to_dne(lat0, lon0, lat, lon, dnorth+iv, deast+iv);
		} catch (int ierr)
		{
			elog_die(0,"Data lack required parameters lat and lon\n");
		}
	}
	//
	// We want to make sure that all stations are in standard coordinates
	//
	for(iv=indata.begin();iv<indata.end;++iv)
	{
		indata.tcse[iv].rotate_to_standard();
	}
	// 
	// the weights become are a nsta by nsamp matrix to allow 
	// variable length apertures.  This algorithm assumes the
	// aperture object input is an nsta by 

	double *weights[nsta];
	for(i=0;i<nsta;++i) 
		new double *weights[i]=double[nsout];
	int stack_count=0;  //maximum number of traces not zeroed
	for(i=0;i<nsout;++i)
	{
		double work[nsta],nused;
		
		if(i>=aperture.npoints)
		{
			// backfill with the last value if the vector is short 
			// we do this silently as it is should be viewed as a 
			// trivial problem we handle automatically
			compute_pseudostation_weights(nsta,dnorth,deast,
			    aperture.get_aperture(tstart+dt*(double)i),
			    aperture.get_cutoff(tstart+dt*(double)i),
			    work,&nused);
		}
		else
		{
			compute_pseudostation_weights(nsta,dnorth,deast,
			    aperture.aperture[i],aperture.cutoff[i],
			    work,&nused);
		}
		for(j=0;j<nsta;++i) weights[j][i]=work[j];
		stack_count=max(nused,stack_count);
	}
	///
	// ERROR RETURN THAT NEEDS TO BE HANDLED GRACEFULLY
	// I don't throw an exception because this should not be viewed as
	// an exception.  It is a case that has to be handled gracefully.
	// it will happen often at the edges of arrays
	//
	if(stack_count<=0) return(-1)


	//
	// Loop over slowness grid range storing results in new output ensemble
	//
	int nslow = ugrid.nux*ugrid.nuy;

	for(i=0;i<ugrid.nux;++i)
	for(j=0;j<ugrid.nuy,++j)
	{
		double moveout[nsta];
		double stack[3][nsout];
		double stack_weight[nsout];
		double twork[nsout];

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
		compute_pmmoveout(nsta,deast,dnorth,dux,duy,moveout);
		for(iv=indata.begin(),ismute_this=0,iend_this=0,ismin=nsout;
			iv<indata.end;++iv)
		{
			int is0, is, ismin, ns_to_copy; 
			is0 = istart + moveout[iv];
			// need this to define latest front mute time
			ismute_this = max(ismute+moveout[iv],ismuth_this);
			// This computes starting position for the input trace  
			if(is0<0)
			{ 
				is=0;
				ismin = 0;
				ns_to_copy = min(ns+is0,nsout);
				iend_this = max(iend_this,ns_to_copy-1);
			}
			else
			{
				is = is0;
				ismin = min(is,ismin);
				ns_to_copy = min(ns-is0,nsout-is0);
				iend_this = max(iend_this,ns_to_copy+is0-1);
			}

			for(j=0;j<3;++j)
			{
				dzero(nsout,twork,1);
				dzero(nsout,wwork,1);
				dcopy(ns_to_copy,&(indata.tcse[iv].x[j][is]),1,twork,1);
				vscal(ns_to_copy,weights[iv],1,twork,1);
				vadd(ns_to_copy,twork,1,stack[j],1);
				// Need to accumulate the sum of weights to normalize
				if(j==0) vadd(nthis,weights[iv],1,stack_weight,1);
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
		for(j=0;j<3;++j) dcopy(nsout,stack[j],1,stackout.x[j],1)
		copy_selected_metadata(tsce[0][0].md,stackout.md,mdlist);
		//MISSING
		// somewhere here we need to load metadata for the stack object.
		// What needs to be loaded and how requires some thought.
		// putting this aside to work on it in background
		for(j=0;j<3;++j) apply_top_mute(stackout->x[j],stackmute);

		// place holder for now.  Dependent on final mpi design
		// concept is that each stack is pushed to output stream as
		// individual 3component trace objects.
		save_stack(stackout,pfsho);
	}
	for(i=0;i<nsta;++i) delete weights[i];
	return(0);
}
	
