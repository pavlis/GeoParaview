#ifndef _SEISPP_H_
#define _SEISPP_H_
#include <vector>
#include "metadata.h"
#include "pfstream.h"

// This is used to define gaps and/or time variable weights
//
class Time_Window
{
public:
	bool gap;  // set true if this window is a gap
	double tstart, tend;  // time period of window (relative or abs)
	double weight(double t)
	{
		if(gap) return(0.0);
		return(w0+(t-tstart)*wgrad);
	}
private:
	double w0;  // weight at t0
	double wgrad;  // dw/dt 
}

enum Time_Reference_Type {absolute,relative};
class Time_Series
{
public:
	bool live;
	Metadata md;
	double dt,t0;
	int ns;
	Time_Reference_Type tref;
	double *s;
	Time_Series() {lie=0;s=NULL;dt=0.0; t0=0.0; ns=0; tref=absolute;};
	Time_Series(int nsin) {live=0; dt=0.0; t0=0.0; ns=nsin;
		tref=absolute; s = new double[nsin];};
	Time_Series(const Time_Series&);
	Time_Series(const Time_Series *);
	~Time_Series() {  if(s!=NULL) delete [] s; };
	Time_Series& operator=(const Time_Series&);
};

// 
// Spherical coordinate angles are convenient for trace rotation
// definitions
//
typedef struct Spherical_Coordinate_{
        double radius;
	// in radians 
        double theta;
        double phi;
} Spherical_Coordinate;
class Three_Component_Seismogram
{
public:
	Metadata md;
	double dt,t0;
	int ns;
	Time_Reference_Type tref;
	bool components_are_orthogonal;  
	bool components_are_cardinal;  // true if x1=e, x2=n, x3=up
	// This is the transformation matrix applied relative to standard
	double tmatrix[3][3]; 
	//Using a Time_Series object adds overhead and memory use 
	//by duplicating metadata, but adds generality
	Time_Series x[3];
	Three_Component_Seismogram();
	Three_Component_Seismogram(int nsamp);
	Three_Component_Seismogram(const Three_Component_Seismogram&);
	~Three_Component_Seismogram(){}; // intentionally null because default ok
	void rotate_to_standard();
	// This overloaded pair do the same thing for a vector
	// specified as a unit vector nu or as spherical coordinate angles
	void rotate(Spherical_Coordinate);
	void rotate(double nu[3]);
	// This applies a general transform with a 3x3 matrix.  
	// User should set components_are_orthogonal true if they
	// are after this transformation as the default assumes no
	// Note this is routine does NOT return to standard before
	// applying the transformation so this is accumulative.
	void apply_transformation_matrix(double a[3][3]);
};
class Time_Series_Ensemble
{
public:  
	// These are nominal sizes.  Individual members should not be assumed
	// necessarily to be the same
	//
	int nts, nsamp;
	Metadata md;
	vector <Time_Series> tse;

	Time_Series_Ensemble();
	Time_Series_Ensemble(int ntsin, int nsampin);
	Time_Series_Ensemble(Pf_ensemble *pfe);
	// This function is bad form, but I see no alernative
	// It's purpose is to initialize the namespace for the 
	// global metadata list.  It should not be advertised
	void add_to_ensemble_mdlist(char *s){pushtbl(mdlist,s);};
private:
	// mdlist contains a list of names to be pulled from input pf 
	// and copied to global Metadata area for the ensemble
	// Note no typing is needed as this is a string copy
	Tbl *mdlist;
};

class Three_Component_Ensemble
{
public:
	int nsta, nsamp;
	Metadata md;
	vector <Three_Component_Seismogram> tcse;

	Three_Component_Ensemble();
	Three_Component_Ensemble(int nsta, int nsamp);
	Three_Component_Ensemble(Pf_ensemble *pfe);
	void add_to_ensemble_mdlist(char *s){pushtbl(mdlist,s);};
private:
	Tbl *mdlist;
};
//
//  Mute definitions
//
class Top_Mute
{
public:
	double t0e, t1;  // t0e is end of to 0 zone, t1 time when weight goes to 1
	 Time_Reference_Type reftype;   // from seispp is enum as absolute or relative
	Top_Mute(){t0e=1.0; t1=2.0; reftype=relative;};
	Top_Mute(Pf *pf);
};
//
// We use an abstract base class and derived classes to allow variations
// in errors thrown.  This uses methods described in books by Stroustrup
//
class seispp_error
{
public:
	string message;
	seispp_error(){message="seispp library error\n";};
	seispp_error(const string mess){message=mess;};
	virtual void log_error(){cerr << "seispp error: "<<message<<endl;};
};

class SAC_data_error : public seispp_error
{
public:
	SAC_data_error(const string mess){message=mess;};
	virtual void log_error()
	{
		cerr<<"Error processing SAC format time series"<<endl;
		cerr<<"Error message = "<<message;
	}
};
// Warning the following names collide with location.h
class Slowness_vector
{
public:
	double ux,uy;   // base vector stored as components in s/km units
	double mag(){return(hypot(ux,uy));};
	double azimuth(){
		double phi;
		phi=M_PI_2-atan2(uy,ux);
		if(phi>M_PI)
			return(phi-2.0*M_PI);
		else
			return(phi);
	};
	double baz(){
		double phi;
		phi = M_PI_2-atan2(-uy,-ux);
		if(phi>M_PI)
			return(phi-2.0*M_PI);
		else
			return(phi);
	};
};
class Rectangular_Slowness_Grid
{
public:
	string name;
	double uxlow, uylow;
	double dux, duy;
	int nux, nuy;
	Rectangular_Slowness_Grid(Pf *);
	double ux(int i) {return(uxlow+i*dux);};
	double uy(int i) {return(uylow+i*duy);};
};
// An assignment operator is not necessary for this object as it
//stands since it is constructed of simple types.
class Hypocenter
{
public:
	double lat,lon,z;
	double time;
	Hypocenter(){method=strdup("tttaup"); model=strdup("iasp91");};  // default
	double distance(double lat0, double lon0);
	double esaz(double lat0, double lon0);
	double seaz(double lat0, double lon0);
	double ptime(double lat0, double lon0, double elev);
	Slowness_vector pslow(double lat0, double lon0, double elev);
	double phasetime(double lat0, double lon0, double elev, string phase);
	Slowness_vector phaseslow(double lat0, double lon0, double elev, string phase);
	void tt_setup(string meth, string mod); // change default method:model
private:
	char *method;
	char *model;
};
//
//Helpers
//
void apply_top_mute(Time_Series &ts,Top_Mute& mute);
void apply_top_mute(Time_Series_Ensemble& t, Top_Mute& mute);
void apply_top_mute(Three_Component_Ensemble &t3c, Top_Mute& mute);


#endif
