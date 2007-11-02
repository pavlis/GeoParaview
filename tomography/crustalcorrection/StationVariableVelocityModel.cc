#include <stdio.h>
#include "seispp.h"
#include "StationVAriableVelocityModel.h"
namespace SEISPP {
using namespace std;
using namespace SEISPP;

StationVariableVelocityModel::StationVariableVelocityModel(DatabaseHandle& dbhraw,
	list<string> vdefs, string default_vname, string type)
{
	default_vmodel_name=default_vname;
	list<string>::iterator inputline;
	char sta[30],modelname[40];
	string property;
	DatascopeHandle dbh=dynamic_cast<DatascopeHandle&>(dbhraw);
	if(type=="P")
		property=string("Pvelocity");
	else if(type=="S")
		property=string("Svelocity");
	else
		throw SeisppError(string("StationVariableVelocityModel database")
		 + string(" constructor.  Illegal model type ")
		+ type + string(" requested") );
	/* This creates a new entry for the default model */
	string defline("default ");
	defline=defline+default_vname;
	vdefs.push_back(defline);
	for(inputline=vdefs.begin();inputline!=vdefs.end();++inputline)
	{
		sscanf(inputline->c_str(),"%s%s",sta,modelname);
		try {
			models[string(sta)] = VelocityModel_1d(dbh.db,
				string(modelname),property);
		} catch (VelocityModel_1d_Dberror vdbe)
		{
			throw SeisppError(string("StationVariableVelocityModel database constructor:  ")	
				+ string("model ")+string(modelname)
				+ string(" not in database") );
		}
	}
}
StationVariableVelocityModel::StationVariableVelocityModel(const StationVariableVelocityModel& old)
{
	models=old.models;
	default_vmodel_name=old.default_vmodel_name;
}
VelocityModel_1d StationVariableVelocityModel::getmod(string mod)
{
	map<string,VelocityModel_1d>::iterator modptr;
	modptr=models.find(mod);
	if(modptr==models.end())
	{
		if(SEISPP_verbose)
			cerr << "Requested model "<<mod
				<< " not found.  Using default model with name="
				<< default_vmodel_name<<endl;
		modptr=models.find(default_vmodel_name);
	}
	return(modptr->second);
}
VelocityModel_1d StationVariableVelocityModel::getmod()
{
	return(models[default_vmodel_name]);
}
}  // End SEISPP namespace declaration
