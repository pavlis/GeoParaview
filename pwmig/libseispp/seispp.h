#ifndef _SEISPP_H_
#define _SEISPP_H_
#include <vector>
#include "metadata.h"


enum Time_Reference_Type {absolute,relative};
class Time_Series
{
public:
	Metadata md;
	double dt,t0;
	int ns;
	Time_Reference_Type tref;
	double *s;
	Time_Series() {s=NULL;dt=0.0; t0=0.0; ns=0; tref=absolute;};
	Time_Series(int nsin) {dt=0.0; t0=0.0; ns=nsin;
		tref=absolute; s = new double(nsin);};
	~Time_Series() {  if(s!=NULL) delete s; };
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
	~Three_Component_Seismogram();
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
	void apply_transformation_matrix(a[3][3]);
};
class Time_Series_Ensemble
{
public:  
	Metadata md;
	vector <Time_Series *> tse;

	Time_Series_Ensemble();
	Time_Series_Ensemble(int nts, int nsamp);
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
	Metadata md;
	vector <Three_Component_Seismogram> tcse;

	Three_Component_Ensemble();
	Three_Component_Ensemble(int nsta, int nsamp);
	Three_Component_Ensemble(Pf_ensemble *pfe);
	void add_to_ensemble_mdlist(char *s){pushtbl(mdlist,s);};
private:
	Tbl *mdlist;
};

#endif
