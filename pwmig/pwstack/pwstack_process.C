#include <math.h>
#include <cstdio>
#include <vector>
#include <sstream>
#include "perf.h"
#include "coords.h"
#include "seispp.h"
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
string virtual_station_name(int ix1, int ix2)
{
	char name[8];
	sprintf(name,"%03.3d%03.3d",ix1,ix2);
	return(string(name));
}
// make station name for pwstack using 3 integer ids
string MakeDfileName(int evid, int x1, int x2)
{
	string dfile;
	ostringstream sbuf(dfile);
	sbuf<<"pwstack_"<<evid<<"_"<<x1<<"_"<<x2;
	return(sbuf.str());
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
	mdlcopy - defines metadata to copy from the raw data
		ensemble to each stack output
	mdlout - list of metadata to be saved to output database 
		(fixed for a given schema)
	am - Attribute_Map object defining mapping of internal to
		external namespace (invariant between calls)
	dir and dfile - define file name where output written (dir/dfile)
	dbh - Database_Handle object for output

Normal return is 0.  Returns nonzero if stack was null.
(This is not an exception as it will happen outside the boundaries 
of an array and this case needs to be handled cleanly.

Throws a Metadata_error exception object if there are problems
parsing required metadata from any trace.  Current caller will
abort the program on this condition, but evolution might want
to produce a handler.
Change Jan 15,2004
Used to pass lat0, lon0 by arg list, now passed through the
ensemble metadata.
*/
int pwstack_ensemble(Three_Component_Ensemble& indata,
	Rectangular_Slowness_Grid& ugrid,
	Top_Mute& mute,
	Top_Mute& stackmute,
	double tstart,
	double tend,
	Depth_Dependent_Aperture& aperture,
	Metadata_list& mdlcopy,
	Metadata_list& mdlout,
	Attribute_Map& am,
	string dir,
	Database_Handle& dbh) 
{
	// lat0 and lon0 are location of target pseudostation grid point
	// elev0 is elevation of datum to use for geometric statics
	//
	double lat0,lon0,elev0;
	// Incident wavefield slowness vector components
	double ux0,uy0;
	string gridname;
	int ix1, ix2;
	int evid;
	string dfile;
	try {
		lat0=indata.get_double("lat0");
		lon0=indata.get_double("lon0");
		elev0=indata.get_double("elev0");
		ux0=indata.get_double("ux0");
		uy0=indata.get_double("uy0");
		ix1=indata.get_int("ix1");
		ix2=indata.get_int("ix2");
		evid=indata.get_int("evid");
		gridname=indata.get_string("gridname");
	} catch (Metadata_get_error& mderr)
	{
		// The above are all set by main so this is the 
		// correct error message.  Could perhaps be dropped
		// after debugging finished.
		mderr.log_error();
		cerr << "Coding error requiring a bug fix"<<endl;
		exit(-1);
	}
	// This routine builds a sta name from the index positions.
	// Made a function to allow consistency between programs.
	//
	string sta=virtual_station_name(ix1,ix2);
	int i,j,k;
	vector<Three_Component_Seismogram>::iterator iv,ov;
	int nsta = indata.tcse.size();
	// This computes the output gather size.  It assumes all data have
	// a common sample rate and we can just grab the first one in the
	// list.  It also quietly assumes a relative time base 
	// so all times are computed relative to the start of
	// each trace.  Caller should guarantee this.
	//
	double dt=indata.tcse[0].dt;
	int nsin = indata.tcse[0].ns;


	Three_Component_Seismogram *stackout;

	int istart = SEISPP::nint(tstart/dt);
	int iend = SEISPP::nint(tend/dt);
	int nsout = iend-istart+1;
	int ismute = SEISPP::nint(mute.t1/dt);
	int ismute_this, ie_this;
	Datascope_Handle& dshandle=dynamic_cast<Datascope_Handle&>(dbh);
	Datascope_Handle dbstack(dshandle);
	dbstack.lookup("pwstack");
	Datascope_Handle dbevl(dshandle);
	dbevl.lookup("evlink");
	Datascope_Handle dbwf(dshandle);
	dbwf.lookup("wfprocess");

	if(istart>=indata.tcse[0].ns)
		elog_die(0,(char *)"Irreconcilable window request:  Requested stack time window = %lf to %lf\nThis is outside range of input data\n",
			tstart,tend);

	/* Apply front end mutes to all traces */
	apply_top_mute(indata,mute);

	/* We need dnorth, deast vectors to compute moveout sensibly
	for this program.  Since we use them repeatedly we will
	extract them once from the gather.*/
	vector <double> dnorth;
	vector <double> deast;
	vector <double> elev;
	dnorth.resize(nsta);   deast.resize(nsta);  elev.resize(nsta);
	for(i=0,iv=indata.tcse.begin();iv!=indata.tcse.end();++iv,++i)
	{
		double lat,lon;
		int ierr;
		try{
			lat = (*iv).get_double("site.lat");
			lon = (*iv).get_double("site.lon");
			// assume metadata store these in degrees so we have to convert
			lat = rad(lat);
			lon = rad(lon);
			geographic_to_dne(lat0, lon0, lat, lon, &(dnorth[i]),&(deast[i]));
			elev[i]=(*iv).get_double("site.elev");
		} catch (Metadata_error& merr)
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

	dmatrix weights(nsta,nsout);
	int stack_count=0;  //maximum number of traces not zeroed
	vector <double> work(nsta);
	for(i=0;i<nsout;++i)
	{
		int nused;
		
		nused=compute_pseudostation_weights(nsta, &(dnorth[0]),&(deast[0]),
		    aperture.get_aperture(tstart+dt*(double)i),
		    aperture.get_cutoff(tstart+dt*(double)i),&(work[0]));
		for(j=0;j<nsta;++j) weights(j,i)=work[j];
		stack_count=max(nused,stack_count);
	}
	///
	// ERROR RETURN THAT NEEDS TO BE HANDLED GRACEFULLY
	// I don't throw an exception because this should not be viewed as
	// an exception.  It is a case that has to be handled gracefully.
	// it will happen often at the edges of arrays
	//
	if(stack_count<=0) return(-1);
	vector <double>moveout(nsta);
	dmatrix stack(3,nsout);
	vector<double>stack_weight(nsout);
	vector<double>twork(nsout);

	double avg_elev,sum_wgt;  // computed now as weighted sum of station elevations
	for(i=0,avg_elev=0.0,sum_wgt=0.0;i<nsta;++i)
	{
		avg_elev += weights(i,0)*elev[i];
		sum_wgt += weights(i,0);
	}
	avg_elev /= sum_wgt;

	//
	// Loop over slowness grid range storing results in new output ensemble
	//

	for(int iu=0;iu<ugrid.nux;++iu)
	{
	    for(int ju=0;ju<ugrid.nuy;++ju)
	    {
		double ux,uy,dux,duy;
		int iend_this;
		int ismin;
		char buffer[128];
		string table_name("wfprocess");

		stack.zero();

		dzero(nsout,&(stack_weight[0]),1);

		ux = ugrid.ux(iu);
		uy = ugrid.uy(ju);	
	
		/* The input gather is assumed prealigned with the slowness
		vector defined by raygath->ux,uy.  We use relative moveouts
		from this base moveout for the actual stacks */
		dux = ux - ux0;
		duy = uy - uy0;
		// moveout computed here and used below assumes the
		// data are aligned on the P arrival
		compute_pwmoveout(nsta,&(deast[0]),&(dnorth[0]),dux,duy,&(moveout[0]));
		for(i=0,iv=indata.tcse.begin(),ismute_this=0,iend_this=0,ismin=nsout;
			iv!=indata.tcse.end();++iv,++i)
		{
			int is0,ie0, is, ns_to_copy; 
			int j0;
			double lag;
			lag = tstart - (iv->t0) + moveout[i];
			is0=SEISPP::nint(lag);
			if(is0>=0)
			{
				j0=0;
				is=is0;
				ismin = min(is,ismin);
				ie0=is0+nsout-1;
				if(ie0<=nsin)
				{
					ns_to_copy=nsout;
					iend_this=nsout;
				}
				else
				{
					ns_to_copy=nsin-is;
					iend_this=max(ns_to_copy,iend_this);
				}
			}
			else
			{
				j0=-is0;
				is=0;
				ismin=0;
				ie0=nsout+is0-1;  // a subtraction by is0 since is0 is negative
				ns_to_copy=nsout+is0;
				if(ns_to_copy>nsin) 
				{
					ns_to_copy=nsin+is0;
					iend_this=max(nsin,iend_this);
				}
				else
				{
					iend_this=nsout;
				}
			}

			for(j=0;j<3;++j)
			{
				dzero(nsout,&(twork[0]),1);
				// ugly expression below depends on data being
				// stored in FORTRAN order in the dmatrix object
				dcopy(ns_to_copy,(*iv).u.get_address(j,is),3,
					&(twork[j0]),1);
				vscal(nsout,weights.get_address(i,0),nsta,&(twork[0]),1);
				vscal(nsout,weights.get_address(i,0),nsta,&(twork[0]),1);
				vadd(nsout,&(twork[0]),1,stack.get_address(j,0),3);
				// Need to accumulate the sum of weights to normalize
				if(j==0) vadd(nsout,weights.get_address(i,0),nsta,&(stack_weight[0]),1);
			}
		}
		// normalize the stack.
		// not so trivial because of machine zero issue
		// Use a hard limit based on ismax and iend_this and then 
		// assume mutes kill possible transients in mute zone
		for(j=0;j<3;++j)
			for(i=ismin;i<iend_this;++i)
			{
				if(stack_weight[i]>0)
					stack(j,i)/=stack_weight[i];
			}
		// Create the output stack as a 3c trace object and copy
		// metadata from the input into the output object.
		stackout = new Three_Component_Seismogram(nsout);
		stackout->dt=dt;
		stackout->t0=tstart;
		stackout->live=true;
		stackout->u=stack;
		copy_selected_metadata(dynamic_cast<Metadata&>(indata),
			dynamic_cast<Metadata&>(*stackout),mdlcopy);
		stackout->put_metadata("ux",ux);
		stackout->put_metadata("uy",uy);
		stackout->put_metadata("dux",dux);
		stackout->put_metadata("duy",duy);
		// These could be copied with the copy_selected_metadata 
		// call above, but it is better to do this explicitly 
		// as these are completely internal parameters
		//
		stackout->put_metadata("lat0",lat0);
		stackout->put_metadata("lon0",lon0);
		stackout->put_metadata("elev0",elev0);
		stackout->put_metadata("ux0",ux0);
		stackout->put_metadata("uy0",uy0);
		// may want to output a static here, but it is probably better to 
		// just keep a good estimate of elevation and deal with this in the
		// migration algorithm. 
		stackout->put_metadata("elev",avg_elev);
		// These are needed in out waveform output file and 
		// with this method they need to be pushed to the Metadata
		// area of the stackout object.
		stackout->put_metadata("wfprocess.time",stackout->t0);
		stackout->put_metadata("wfprocess.endtime",stackout->endtime());
		stackout->put_metadata("wfprocess.timetype","r"); // always relative
		stackout->put_metadata("wfprocess.dir",dir);
		dfile=MakeDfileName(evid, ix1, ix2);
		stackout->put_metadata("wfprocess.dfile",dfile);
		stackout->put_metadata("wfprocess.datatype","3c");
		stackout->put_metadata("samprate",1.0/dt);
		stackout->put_metadata("nsamp",nsout);
		stackout->put_metadata("wfprocess.algorithm","pwstack");

		apply_top_mute(*stackout,stackmute);


		try{
			int pwfid;
			int rec;
			rec=dbsave(*stackout,dshandle.db,table_name,mdlout,am);
			// negative rec means nothings saved.  
			// skip silently because this happens when 
			// marked dead
			if(rec>=0)
			{
				dbwf.db.record=rec;
				dbgetv(dbwf.db,0,"pwfid",&pwfid,0);
				// now fill in the special table for this program
				dbstack.append();
				dbstack.put("gridname",gridname);
				dbstack.put("ix1",ix1);
				dbstack.put("ix2",ix2);
				dbstack.put("ux",ux);
				dbstack.put("uy",uy);
				dbstack.put("ux0",ux0);
				dbstack.put("uy0",uy0);
				dbstack.put("pwfid",pwfid);
				dbstack.put("sta",sta);
				// Perhaps should put this out as a static, but we
				// save the elevation instead and compute it in 
				// the migration algorithm where a velocity available.
				// Here we'd have to get that information to this point
				dbstack.put("elev",avg_elev);
				//
				// Now for event link table -- awful database table, but 
				// see now way around it if we are to retain a css back
				// compatibility
				//
				dbevl.append();
				dbevl.put("evid",evid);
				dbevl.put("pwfid",pwfid);
			}

		}
		catch(seispp_error& err)
		{
			err.log_error();
			cerr << "Write failure abort:  cannot continue"<<endl;
			exit(-1);
		}
	    }
	}
	return(0);
}
	
