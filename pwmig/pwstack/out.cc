/* This file contains output functions used by pwstack.  */
#include "seispp.h"
#include "pwstack.h"
void save_stack(Three_Component_Seismogram& stack,
	Database_Handle& dbhg,
		Metadata_list& mdl,
			Attribute_Map& am)
{
    try {
	Datascope_Handle& dbh_base=dynamic_cast<Datascope_Handle&>(dbhg);
	// These will point at each table.  It is a bit inefficient
	// to do this for each output trace, but I'm assuming 
	// until proven otherwise that the overhead is small
	Datascope_Handle dbwf(dbh_base);
	Datascope_Handle dbsite(dbh_base);
	Datascope_Handle dblink(dbh_base);
	Datascope_Handle dbpws(dbh_base);
	dbwf.lookup("wfprocess");
	dbsite.lookup("vsite");
	dblink.lookup("vsitelink");
	dbpws.lookup("pwstack");

	int pwfid,vsid;
	pwfid=dbnextid(dbwf.db,"pwfid");
	vsid = dbnextid(dbsite.db,"vsid");
	// Save the waveform data first as it is most likely to 
	// generate an error.
	stack.put_metadata("pwfid",pwfid);
	dbsave(stack,dbwf.db,"wfprocess",mdl,am);

	double lat,lon,elev;
	lat=stack.get_double("lat");
	lon=stack.get_double("lon");
	elev=stack.get_double("elev");
	dbsite.append();
	dbsite.put("vsid",vsid);
	dbsite.put("dlat",lat);
	dbsite.put("dlon",lon);
	dbsite.put("elev",elev);

	dblink.append();
	dblink.put("pwfid",pwfid);
	dblink.put("vsid",vsid);

	double ux,uy;
	ux=stack.get_double("ux");
	uy=stack.get_double("uy");
	dbpws.append();
	dbpws.put("ux",ux);
	dbpws.put("uy",uy);
	dbpws.put("pwfid",pwfid);
	dbpws.put("vsid",vsid);

    } catch (... ) { throw; };
}
