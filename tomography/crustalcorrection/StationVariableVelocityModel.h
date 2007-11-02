#include <map>
#include <list>
#include <iostream>
#include "VelocityModel_1d.h"
#include "dbpp.h"
namespace SEISPP {
/*! \brief Provides a way to associate variable velocity models with keywords (normally station names).

Some processing methods benefit from allowing variable velocity models at individual seismic
stations.  A good example for which this was originally written is crustal corrections for
body wave tomography.  The object provides a convenient interface to utilize variable velocity
models of this form.  To use this object first construct it.  There is only one constructor
from a database, but one could readily extend this using some file-based methodology.  I have
chosen to not have that in the interface as it would only encourage bad behaviour.  
\author Gary L. Pavlis
*/
class StationVariableVelocityModel
{
public:
	/*! Construct from a database. 

	\param dbh database containing velocity models.  Abstracted to be totally general.
		through passing of a base class.   Implementation will normally pass a 
		handle to a child of this abstract base class.
	\param def list of strings contain pairs of two keywords.  First token in the string
		is assumed to be the station name and the second is a keyword defining a 
		particular velocity model to be attached to that station name. 
	\param default_vname is the default velocity model name.
	\param type is the type of model requested.  Should be either "P" or "S".
		Anything else will cause and exception to be thrown.
	\exception throws a SeisppError exception if anything goes wrong.  The most common
		such error is referencing a velocity model not actually in the database.
	*/
	StationVariableVelocityModel(DatabaseHandle& dbh,list<string> vdef,
		string default_vname, string type);
	/*! Standard copy constructor. */
	StationVariableVelocityModel(const StationVariableVelocityModel& parent);
	/*! Returns default velocity model.*/
	VelocityModel_1d getmod();
	/*! \brief Return a model by name.
	
	This is the primary method of this object.  Returns a velocity model keyed by the
	name given.  If the requested model is not know the default model will be returned.
	Note this references SEISPP::SeisppVerbose.  When true missing model errors will
	cause messages to be sent to stderr.  Otherwise the default will be silenetly returned.
	*/
	VelocityModel_1d getmod(string name);
private:
	map<string,VelocityModel_1d> models;
	string default_vmodel_name;  
};
} /* End namespace declaration */
