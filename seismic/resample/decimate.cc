/* This group of routines implement the Fir_Decimator_Filter object */
#include "decimate.h"
// Constructor for Decimation_Definitions object using a pf
Decimation_Definitions::Decimation_Definitions(Pf *pf)
{
	Tbl *dlist;
	char *line;
	int i;
	char *white=" \t";
	Decimator *d;
	Decimator& dr=*d;
	char *s;
	char *decfiles;

	dlist=pfget_tbl(pf,"decimators");
	for(i=0;i<maxtbl(dlist);++i)
	{
		d = new Decimator();
		line=gettbl(dlist,i);
		s = strtok(line,white);
		dr.exact=atof(s);
		s = strtok(NULL,white);
		dr.start = atof(s);
		s = strtok(NULL,white);
		dr.end = atof(s);
		s = strtok(NULL,white);
		do
		{
			decfile=strdup(s);
			pushtbl(dr.declist,decfile);
		} while((s=strtok(NULL,white))!=NULL);
		decmap[dr.exact]=dr;  //THIS IS MOST CERTAINLY WRONG.,  NEED TO CHECK c++ MANUAL
		
	}
}
		



/* This is the primary application function of the group of functions 
contained in this file.  It takes an input time series object,
finds the correct decimation definition for the given trace, and 
if a decimator is needed it is applied.  The function returns a 
new Time_Series object (it may be he same as the original if no
decimation or resampling was required) that is decimate/resampled
onto the target sample rate defined by the Decimation_Definition
object.  

*/

Time_Series  Decimate_Time_Series(Time_Series& indata,Decimation_Definitions& dec)
{
	double samprate;
	Decimator d;
	int decresult;

	samprate=1.0/indata.dt;

	d = dec.decmap[samprate];
	// Return immediately if no decimation or interpolation is
	// required.  Use a safe test for samprate because it is 
	// a computed quantity.  Can't just use ==
	if((d.decfac==1.0) && (fabs(d.exact-samprate)<=DBL_EPSILON))
		return(indata);
	// Apply decimation defined by d.
	// This uses an existing routine in libmultiwavelet that uses
	// floats which is why the float arrays are needed.
	float *intmp=new float[indata.ns];
	float *otmp;
	double odt,ot0;
	int i;
	for(i=0;i<indata.ns;++i) intmp[i]=(float)indata.s[i];

	decresult = decimate_trace(d.declist,intmp,indata.ns,
			indata.dt,indata.t0,&otmp,&ons,&odt,&ot0);

	samprate = 1.0/odt;
	Time_Series outdata(ons);
	outdata.live=true;
	outdata.md=indata.md;
	outdata.dt = odt;
	outdata.t0=ot0;
	outdata.ns = ons;
	outdata.tref = indata.tref;
	if(fabs(samprate-d.exact)<=DBL_EPSILON)
	{
		// come here when resampling is not needed
		for(i=0;i<ons;++i)outdata.s[i]=otmp[i];
	}
	else
	{
		double te_computed;
		double *stmp = new double[ons];
		for(i=0;i<ons;++i)outdata.s[i]=otmp[i];
		// compute the size of the resampled vector 
		// If sample rate is far off won't match
		outdata.dt = 1.0/d.exact;
		te_computed = ot0+((double)(ons-1))*odt;
		outdata.ns = nint((te_computed-t0)/outdata.dt)-1;
		if(outdata.ns!=ons)
		{
			delete [] outdata.s;
			outdata.s = new double[outdata.ns];
		}
		linear_scalar_regular_to_regular(ons,ot0,odt,otmp,
				outdata.ns,outdata.t0,outdata.dt,outdata.s);
	}
	free(otmp);
	delete [] intmp;
	return(outdata);
}
		
