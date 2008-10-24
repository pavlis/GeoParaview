#include <libgen.h>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <algorithm>
#include "iterator"
#include "stock.h"
#include "pf.h"
#include "seispp.h"
#include "Hypocenter.h"
#include "AttributeMap.h"
#include "Metadata.h"
#include "TimeSeries.h"
#include "dmatrix.h"
#include "ThreeComponentSeismogram.h"
#include "EventCatalog.h"
using namespace std;
using namespace SEISPP;
/* HERE TO END OF NEXT PROCEDURE:  Eventually will go to libseispp*/
#include <exception>
/*! \brief Write a database tuple from a Metadata object.

Sometimes it is useful to write a set of attributes tagged to a set of
data objects that are subclasses of Metadata.  Examples in SEISPP are
TimeSeries and ThreeComponentSeismogram objects.  This method saves
only a specified list of attributes and not associated data such
as seismic data stored within a TimeSeries object.  A common
reason for this is indexing an existing data file where the
data already exist and all that is needed are appropriae entries
for wfdisc or wfprocess.  Another is that Metadata provide simply
a convenient way to store attributes for one row a database in
a convenient package.

Be aware this routine will always add a new row to the db.
Errors can render this row incomplete as it is filled out
sequentially based on entrie in the MetadataList.

Note this procedure is very similar to an internal function
in readwrite.cc called save_metadata_for_object.  The later
may eventually be depricated.

\param md Metadata containing attributes to be written to db.
\param dbh DatabaseHandle to which contents of md are to be written.
This is a past as the base class for a generic db, but
the current implementation only works for a DatascopeHandle.
Current version immediately casts dbh to a DatascopeHandle.
\param mdl Internal to external namespace mapper that defines what
attributes will be written to the db.  The function will
attempt to update ALL attributes listed.
\return number record number of added tuple
\exception SeisppError is thrown for some problems.  In virtually
all cases this means the new row was probably created but
is probably incomplete.
*/
int dbsave(Metadata& md, DatabaseHandle& dbh, MetadataList& mdl,
AttributeMap& am)
throw(SeisppError)
{
    MetadataList::iterator mdli;
    map<string,AttributeProperties>::iterator ami,amie=am.attributes.end();
    map<string,AttributeProperties> aliasmap;
    string cval;
    const string base_message("dbsave:  ");

    /* We try to cast from the virtual DatabaseHandle to
    the a DatascopeHandle */
    int rec;
    DatascopeHandle dsdbh;
    try
    {
        dsdbh=dynamic_cast<DatascopeHandle&>(dbh);
        rec=dsdbh.append();
    } catch (exception& e)
    {
        string message;
        message=base_message
            + "dynamic cast failed.  Probable coding error. stdlib error:\n"
            + e.what();
        throw SeisppError(message);
    }
    catch (SeisppError serr){throw serr;};

    Dbptr db=dsdbh.db;
    Dbvalue dbval;
    if(dbquery(db,dbTABLE_NAME,&dbval)<0)
    {
        throw SeisppError(base_message
            + "dbquery failed asking for dbTABLE_NAME");
    }
    char *sbuf=dbval.t;
    string table(sbuf);

    for(mdli=mdl.begin();mdli!=mdl.end();++mdli)
    {
        double dval;
        int ival;
        string mdkey;
        if(am.is_alias((*mdli).tag))
        {
            aliasmap=am.aliases((*mdli).tag);
            ami=aliasmap.find(table);
            if(ami==aliasmap.end())
            {
                dbmark(db);
                throw SeisppError(base_message
                    + string("Alias name=")
                    + (*mdli).tag
                    + string(" is not associated with table=")
                    + table
                    + string("\nVerify output specification against schema") );
            }
            mdkey=(*mdli).tag;                    // in this case the alias is the key
        }
        else
        {
            mdkey=(*mdli).tag;
            ami = am.attributes.find(mdkey);
            if(ami==amie)
            {
                dbmark(db);
                throw SeisppError(
                    string("Required attribute ")
                    +(*mdli).tag
                    +string(" cannot be mapped to output namespace"));
            }
            if( (ami->second.db_table_name) != table)
            {
                dbmark(db);
                throw SeisppError(
                    string("dbsave (database table mismatch): attribute ")
                    + ami->second.db_attribute_name
                    + string(" is tagged with table name ")
                    + ami->second.db_table_name
                    + string("expected to find ")
                    + table);
            }
            /* In this case the key we use the name from the Attribute map as the key */
            mdkey=ami->second.internal_name;
        }
        try
        {
            switch(ami->second.mdt)
            {
                case MDint:
                    if(ami->second.is_key)
                    {
                        ival = dbnextid(db,
                            const_cast<char *>
                            (ami->second.db_attribute_name.c_str()) );
                        if(ival<0)throw SeisppError(
                                string("dbsave:  ")
                                + ami->second.db_attribute_name
                                + string(" is defined as integer key for table ")
                                + ami->second.db_table_name
                                + string(" but dbnextid failed") );

                    }
                    else
                        ival = md.get_int(mdkey);
                    dbputv(db,0,ami->second.db_attribute_name.c_str(),
                        ival,0);
                    // In this case we need to push this back to metadata
                    // so it can be used downstream
                    md.put(ami->second.db_attribute_name,ival);
                    break;
                case MDreal:
                    dval = md.get_double(mdkey);
                    dbputv(db,0,ami->second.db_attribute_name.c_str(),
                        dval,0);
                    break;
                case MDstring:
                    cval = md.get_string(mdkey);
                    dbputv(db,0,ami->second.db_attribute_name.c_str(),
                        cval.c_str(),0);
                    break;
                case MDboolean:
                    // treat booleans as ints for external representation
                    // This isn't really necessary as Antelope
                    // doesn't support boolean attributes
                    if(md.get_bool(mdkey))
                        ival = 1;
                    else
                        ival = 0;
                    dbputv(db,0,ami->second.db_attribute_name.c_str(),
                        ival,0);
                    break;

                case MDinvalid:
                    cerr << "dbsave: database attribute "
                        << ami->second.db_attribute_name
                        << " was marked as invalid\n"
                        << "Data for this attribute not saved"
                        << endl;
                    break;

                default:
                    cerr << "dbsave: database attribute "
                        << ami->second.db_attribute_name
                        << " has unknown data type\n"
                        << "Data for this attribute not saved"
                        << endl;
            }

        }
        catch (MetadataGetError& mderr)
        {
            mderr.log_error();
            throw SeisppError(
                string("dbsave object failure from problem in metadata components"));
        }
    }
    return(rec);
}


/*********END SECTION FOR SEISPP ***********/
/* This is used for receiver cooordinate metadata*/
const string rlatkey("rlat"), rlonkey("rlon"),relevkey("relev");

class SHHdrAttribute
{
    public:
        SHHdrAttribute(){key="";};
        SHHdrAttribute(string line);
        string key;
        string shtype;
        int qfnum;
        string shname;
        string table;
        string outname;
        MDtype mdt;
};
SHHdrAttribute::SHHdrAttribute(string line)
{
    istringstream ss(line);
    ss >> key;
    ss >> shtype;
    ss >> qfnum;
    ss >> shname;
    ss >> table;
    if(table=="NULL") table="";
    ss >> outname;
    string mdtname;
    ss >> mdtname;
    if(mdtname=="int")
        mdt=MDint;
    else if(mdtname=="real")
        mdt=MDreal;
    else if(mdtname=="string")
        mdt=MDstring;
    else
        mdt=MDinvalid;
}


typedef struct map<string,SHHdrAttribute> SHHdrMap;
SHHdrMap load_shhdrmap(Pf *pf)
{
    SHHdrMap result;
    Tbl *t=pfget_tbl(pf,"SHHdrAttributes");
    if(t==NULL)
    {
        cerr << "Parameter file is missing required Tbl SHHdrAttributes"
            <<endl;
        exit(-1);
    }
    int i;
    for(i=0;i<maxtbl(t);++i)
    {
        char *s;
        s=(char *) gettbl(t,i);
        SHHdrAttribute sattr(s);
        result[sattr.key]=sattr;
    }
    return(result);
}


/* This uses the istream getline method to parse the SH header and return
a vector of strings with the newlines removed and arranged as a vector of
strings.  That is, the return result will contain one single string with
header data for each member of the vector of strings.  e.g. if we asked
for result[2] this would contain a single std string for entries tagged
with "02|" in the original header file.

Note that both the newline and the leading index (00|, 01|, etc.) are s
stripped to make a string that can be parsed with simpler rules. */

vector<string> parse_hdr(istream& in)
{
    /* buffer for getline method */
    char buf[128];                                /* Probably could be 80, but better safe than sorry.*/
    vector<string> result;
    /* This string is used to build up each component of the output */
    string outline;
    int linecount(0);
    string linestr;
    int nseis,nlastseis;
    /* Always compare the first line to the magic string here.  It seems
    all Q files have this magic string as the first line.  We'll print
    a warning if there is a mismatch.*/
    if(!in.getline(buf,128))
        throw SeisppError("Empty input stream.  EOF encountered immediately.");
    const string Qmagic("43981");
    string sbuf=string(buf);
    sbuf=sbuf.substr(0,5);
    if(Qmagic!=sbuf)
        cerr << "Warning:  Q file magic string not found on first line"
            << " of header file.  "<<endl
            << "Read this: "<< sbuf<<endl
            << "Expected to find: "<<Qmagic<<endl;

    /* getline returns *this or a NULL on EOF, so the books says this should work*/
    while(in.getline(buf,128))                    // default unix newline character since SH is unix always
    {
        string wholeline(buf);
        if(linecount==0)
        {
            nseis=1;
            nlastseis=1;
        }
        /* This strips leading "xx|" used for trace numbers*/
        string thispiece=wholeline.substr(3);
        /* this turns trace number string to an int */
        linestr=wholeline.substr(0,2);
        nseis=atoi(linestr.c_str());
        if(nseis!=nlastseis)
        {
            result.push_back(outline);
            nlastseis=nseis;
            outline=thispiece;
        }
        else
        {
            outline += thispiece;
        }
        ++linecount;
    }
    /* Have to push the last entry before we return */
    result.push_back(outline);
    return result;
}


/* Procedure for function immediately below. Seismic Handler uses an
odd string for storing time fields.  This procedure takes that string apart
and converts it to an epoch time.  The algorith depends on the fact that str2epoch is
pretty smart and seems to work fine by converting underscore ("_") characters to blanks.*/
double shtime_to_epoch(string s)
{
    char *buffer;
    buffer=strdup(s.c_str());
    int i;
    // tacitly assume s was correctly copied and buffer has same length
    for(i=0;i<s.length();++i)
    {
        if(buffer[i]=='_') buffer[i]=' ';
        if(buffer[i]=='-') buffer[i]=' ';
    }
    double result=str2epoch(buffer);
    free(buffer);
    return(result);
}


/* Takes a string produced from the parse_hdr procedure above and converts it
to a Metadata object using hdrmap.  A number of things had to be hard coded into
this that are specific odd properties of SH.  That is mainly the separator and
the oddity of storing time as a string*/
Metadata shstr2md(SHHdrMap& hdrmap, string hdrstr)
{
    Metadata result;
    SHHdrMap::iterator hmptr,hmeptr;
    hmeptr=hdrmap.end();
    string shkey,hdrval;
    size_t found;                                 // index position returned by find method of std::string
    size_t current(0);
    int ival;
    double dval;
    int endstr=hdrstr.size();
    do
    {
        found=hdrstr.find(":",current);
        if(found==string::npos)
        {
            throw SeisppError("The following line is scrambled and cannot be parsed\n"
                + hdrstr);
        }
        shkey=hdrstr.substr(current,found-1);
        shkey.assign(hdrstr,current,found-current);
        current=found+1;
        found=hdrstr.find("~ ",current);
        hdrval.assign(hdrstr,current,found-current);
        current=found+2;
        /* Now we have to take these two strings and convert them to something
        we can post into the result. */
        hmptr=hdrmap.find(shkey);
        if(hmptr==hmeptr)
        {
            cerr << "shstr2md:  I do not know how to handle SH key = "<<shkey
                << endl << "Check your parameter file and the SH documentation"
                <<endl;
        }
        else
        {
            switch(hmptr->second.mdt)
            {
                case MDint:
                    ival=atoi(hdrval.c_str());
                    result.put(hmptr->second.outname,ival);
                    if(hmptr->second.table.length()>0)
                        result.put(hmptr->second.table
                            +"."+hmptr->second.outname,ival);
                    break;
                case MDreal:
                    if(hmptr->second.shtype=="time")
                        dval=shtime_to_epoch(hdrval.c_str());
                    else
                        dval=atof(hdrval.c_str());
                    result.put(hmptr->second.outname,dval);
                    if(hmptr->second.table.length()>0)
                        result.put(hmptr->second.table
                            +"."+hmptr->second.outname,dval);
                    break;
                case MDstring:
                default:
                    result.put(hmptr->second.outname,hdrval);
            };

        }
    }
    while((found!=string::npos) && (current<endstr));
    return(result);
}


/* Names and some attribute concepts clash in any format and this is no
exception.  This takes QHD file contents and uses them to generate required
Metadata for downstream procesing.  A SeisppError can be thrown if
a required attribute is missing.  datatype,dir, and dfile have to be passed
in as they are not a concept used in Seismic Handler but is a core
idea in Antelope */
void set_required(Metadata& md,string dir, string dfile, string datatype)
{
    int err=1;
    try
    {
        int nsamp=md.get_int("nsamp");
        double dt=md.get_double("dt");
        double samprate=1.0/dt;
        md.put("samprate",samprate);
        double time,endtime;
        time=md.get_double("time");
        endtime=time+dt*static_cast<double>(nsamp-1);
        md.put("endtime",endtime);
        string comp,chan1,chan2;
        comp=md.get_string("comp");
        /* pchan is used in wfprocess.  We use comp as an equivalent
        concept here */
        md.put("pchan",comp);
        chan1=md.get_string("chan1");
        chan2=md.get_string("chan2");
        string chan=chan1+chan2+comp;
        md.put("chan",chan);
        md.put("dir",dir);
        md.put("dfile",dfile);
        md.put("pchan",comp);
        md.put("datatype",datatype);
        /* As far as I can tell SH only saves data in absolute
        time frames */
        md.put("timetype","a");
        md.put("wfprocess.algorithm","SH");
    } catch(MetadataGetError mderr)
    {
        mderr.log_error();
        cerr << "Problems getting required attributes in this Metadata needed to build data objects"<<endl;
        cerr << md;
        throw err;
    }
}


void set_hang_vang(Metadata& mx,double ema,double seaz)
{
    /* LQT in SH seems to be defined as follows:  x1=L = up and away,
    x2=Q up and in seaz direction (backazimuth), and x3=T=LXQ (cross product)*/
    double vang,hang,az;
    az=seaz+M_PI;
    if(az>2.0*M_PI) az -= 2.0*M_PI;
    /* Need to convert all these angles to degrees for external storage*/
    az=deg(az);
    ema=deg(ema);
    seaz=deg(seaz);
    /* intentionally do not use a try catch here */
    string comp=mx.get_string("comp");
    if(comp=="L")
    {
        vang=ema;
        hang=az;
    }
    else if(comp=="Q")
    {
        vang=90.0-ema;
        hang=seaz+90.0;
    }
    else if(comp=="T")
    {
        vang=90.0;
        hang=az+90.0;
        if(hang>360.0) hang -= 360.0;
    }
    else if(comp=="Z")
    {
        hang=0.0;
        vang=0.0;
    }
    else if(comp=="N")
    {
        hang=0.0;
        vang=90.0;
    }
    else if(comp=="E")
    {
        hang=90.0;
        vang=90.0;
    }
    else
    {
        hang=0.0;
        vang=0.0;
    }
    mx.put("hang",hang);
    mx.put("vang",vang);
}

void set_receiver_coordinates(Metadata& md, SeismicArray& receivers)
{
	try {
            string sta=md.get_string("sta");
            map<string,SeismicStationLocation>::iterator rptr;
            rptr=receivers.array.find(sta);
            if(rptr==receivers.array.end())
            {
                cerr << "Warning (set_orientation):  station coordinates for "<<sta
                    << " were not found.  "<<endl
		    << "Coordinates are not set in trace header.  Expect downstream problems"
                    <<endl;
            }
            else
            {
                double rlat,rlon,relev;
                rlat=rptr->second.lat;
                rlon=rptr->second.lon;
                relev=rptr->second.elev;
		md.put(rlatkey,rlat);
		md.put(rlonkey,rlon);
		md.put(relevkey,relev);
            }
		
	}catch(MetadataGetError mderr)
	{
		mderr.log_error();
		cerr << "This should not happen in this code so will exit."
			<< "Contact the author:  pavlis@indiana.edu"<<endl;
		exit(-1);
	}

}

/* This internal procedure ultimately sets the hang and vang metadata attributes
that are needed to define orientation on output channels.  The overall behaviour
is strongly controlled by two parameters:  ema_key and and threecmode.
If ema_key is an empty string (i.e. "") this triggers a computed estimate of
the emergence angle for an incident P wave from the passed hypocenter and using
surface velocity v0.  If it is ema_key is not empty the procedure assumes it
can fetch ema from the passed Metadata objects.  If that fails it will ALWAYS
complain with a message to stderr and fall back to the computed ema method.
(Note if processing a large amount of daa this will generate at error message
for every seismogram that doesn't really may ema_key set. )

If threecmode is true mx1,mx2, and mx3 are assumed to be headers from a
three-component seismogram with comp set either to E,N, and Z (order is
not important) or L,Q, and T (order not important).

arguments:
mx1,mx2,mx3 are Metadata objects created by parsing the SH header.  All
three are hit if threecmode is true.  Only mx1 is used when false.
h - hypocenter of the event corresponding to these data.
array - receiver array that is associated with these data.  The procedure
fetches sta and tries to find that sta in the STL map internal to
array.  Note if the lookup fails ema and back azimuth (seaz) are
set to zero and an error is posted.  If you get this message the db will
need to be fixed.
ema_key - see above
v0 - see above
threecmode - see above
*/
void set_orientation(Metadata& mx1, Metadata& mx2, Metadata& mx3,
Hypocenter& h, SeismicArray& receivers,string ema_key,double v0,
bool threecmode)
{
    try
    {
        double ema,seaz;
        /*Use a zero length ema_key to signal it should be computed. */
        if(ema_key.size()>0)
        {
            try
            {
                ema=mx1.get_double(ema_key);
                /* have to convert to radians if this worked */
                ema=rad(ema);
                mx1.put("orientation_type",string("measured"));
                if(threecmode)
                {
                    mx2.put("orientation_type",string("measured"));
                    mx3.put("orientation_type",string("measured"));
                }
            } catch (MetadataGetError mderr)
            {
                cerr << "Warning (set_orientation):  ema_key parameter is "<<ema_key<<endl
                    <<"This key was not set in the converted SH header.  "
                    <<"Reverting to computed ema method."<<endl;
                ema_key=string("");
            }
        }
        /* Note this is not an else clause.  ema_key can be changed in the exception
        handler above */
        if(ema_key.size()==0)
        {
            mx1.put("orientation_type",string("computed"));
            if(threecmode)
            {
                mx2.put("orientation_type",string("computed"));
                mx3.put("orientation_type",string("computed"));
            }
            string sta=mx1.get_string("sta");
            map<string,SeismicStationLocation>::iterator rptr;
            rptr=receivers.array.find(sta);
            if(rptr==receivers.array.end())
            {
                cerr << "Warning (set_orientation):  station coordinates for "<<sta
                    << " were not found.  Setting ema an seaz for this stations to 0"
                    <<endl;
                ema=0.0;
                seaz=0.0;
            }
            else
            {
		/* These could be retrieved from Metadata but
		this is probably actually faster */
                double rlat,rlon,relev;
                rlat=rptr->second.lat;
                rlon=rptr->second.lon;
                relev=rptr->second.elev;
                seaz=h.seaz(rlat,rlon);
                SlownessVector up=h.pslow(rlat,rlon,relev);
                double u=up.mag();
                double uprod=u*v0;
                if(uprod>1.0)
                    ema=M_PI_2;
                else
                    ema=asin(u*v0);
            }
        }
        set_hang_vang(mx1,ema,seaz);
        /* simply skip these if no running 3c*/
        if(threecmode)
        {
            set_hang_vang(mx2,ema,seaz);
            set_hang_vang(mx3,ema,seaz);
        }
    }
    catch (SeisppError serr)
    {
        throw serr;
    };
}


/* Get Hypocenter information and update EventCatalog if needed.*/
Hypocenter handle_event_data(Metadata& mx1, EventCatalog& evcat, Dbptr db,
string ttmethod,string ttmodel)
{

    double slat,slon,sz,otime;
    /*WARNING:  this is very error prone since it depends
    on name mapping in the parameter file.  Users should
    be warned to not mess with the name mapping Tbl */
    try
    {
        slat=mx1.get_double("lat");
        slon=mx1.get_double("lon");
        sz=mx1.get_double("depth");
        otime=mx1.get_double("otime");
    }
    catch (MetadataGetError mderr)
    {
        cerr << "Problems parsing header data for event data"<<endl
            << "Required location information is missing or incomplete"
            <<endl<< "Details:"<<endl;
        mderr.log_error();
        cerr << "Fatal eroror.  Exiting"<<endl;
        exit(-1);
    }
    Hypocenter h(rad(slat),rad(slon),sz,otime,ttmethod,ttmodel);
    /* This is a truly evil way to accomplish this, so a comment
    is absolutely required given how evil it is.  To simplify
    database updating we handle this here by posting a boolean
    "is_new_evid".  When true the downstream writer will save
    event and origin table entries for event so marked in
    its Metadata.  I suppose a spin doctor would put this a
    different way. This isn't evil it's clever and innovative.*/
    int evid,orid;
    if(evcat.find(h))
    {
        Metadata maux=evcat.current_aux();
        try
        {
            evid=maux.get_int("evid");
            orid=maux.get_int("orid");
            mx1.put("evid",evid);
            mx1.put("orid",orid);
        }
        catch(...)
        {
            throw;
        };
    }
    else
    {
        evid=dbnextid(db,"evid");
        orid=dbnextid(db,"orid");
        mx1.put("evid",evid);
        mx1.put("orid",orid);
        Metadata md;
        md.put("evid",evid);
        md.put("orid",orid);
        md.put("is_new_evid",true);
        try
        {
            double mb=mx1.get_double("mb");
            md.put("mb",mb);
        }
        catch (...)
        {
        };
        evcat.add(h,md);
    }
    return(h);
}
double RFskewtime(vector<TimeSeries>& components, string skewcomp, double dtmax,
	Hypocenter& h, bool normalize)
{
    /* First find the skewcomp channel */
    int i,iskew(-1);
    for(i=0;i<components.size();++i)
    {
    	string comp=components[i].get_string("comp");
	if(comp==skewcomp)
	{
		iskew=i;
		break;
	}
    }
    if(iskew<0) throw SeisppError(string("RRskewtime: required component")
    			+ skewcomp + " not found in 3c set");
    vector<double>::iterator maxamp;
    maxamp=max_element(components[iskew].s.begin(),components[iskew].s.end());
    /* Normalize using a trick.  Change nothing, just alter calib */
    if(normalize)
    {
	double calib;
    	double gain=(*maxamp);  
	for(i=0;i<components.size();++i)
	{
	    try {
		calib=components[i].get_double("calib");
	     } catch (MetadataGetError mderr)
	    {
		calib=1.0;
	    }
	    //calib*=gain;
	    calib=gain;
	    components[i].put("calib",calib);
	}
    }
    else
    {
	for(i=0;i<components.size();++i) components[i].put("calib",1.0);
    }
    //1difference_type idist;  // STL distance function returns this object
    int ip=distance(components[iskew].s.begin(),maxamp);
    double ptime=components[iskew].time(ip);
    double predP;
    double lat,lon,elev;
    try {
    	lat=components[iskew].get_double(rlatkey);
	lon=components[iskew].get_double(rlonkey);
	elev=components[iskew].get_double(relevkey);
	predP=h.ptime(lat,lon,elev) + h.time;
    } catch (MetadataGetError mderr)
    {
    	cerr << "RFskewtime:  error fetching receiver coordinates"<<endl;
	mderr.log_error();
	cerr << "Will return 0 skew time"<<endl;
	return(0.0);
    }
    double skewtime=predP-ptime;
    for(i=0;i<components.size();++i)
    {
    	components[i].t0 += skewtime;
	/* safer to reset these.  Required actually to skew
	scalar traces since we only reference the metatdata */
	components[i].put("time",components[i].t0);
	components[i].put("endtime",components[i].endtime());
	/* Need to reset this even if it is set correctly in the header.
	No guarantee this calculator matches that used in SH */
	components[i].put("Ptime",predP);
    }
    return(skewtime);
}


/* Procedure to save all the attributes or one scalar trace to both wfdisc and wfprocess.
Can fail for a lot fo reason leading to throwing a SeisppError object.

arguments
    d - data to save (writes only db attributes writes no seismic data hence needs
        only the Metadata object)
    dbhwd - wfdisc handle
    dbhwp - wfprocess handle
    dbhevl - evlink handle
    dbhscl - sclink handle
    dbhor - orient handle
mdlwd - MetadataList to drive wfdisc save
mdlwp - MetadataList to drive wfprocess save
*/
void save_scalar_data(Metadata& d,
DatascopeHandle& dbhwd,
DatascopeHandle& dbhwp,
DatascopeHandle& dbhevl,
DatascopeHandle& dbhscl,
DatascopeHandle& dbhor,
MetadataList& mdlwd,
MetadataList& mdlwp,
AttributeMap& am)
{
    try
    {
        /* First save the trace data because if that fails the secondary
        tables are completely junk*/
        int rec=dbsave(d,dynamic_cast<DatabaseHandle&>(dbhwd),
            mdlwd,am);
        /* Seem to need to do this.  Not sure why */
        dbhwd.db.record=rec;
        /* We get everything we need before we try to save wfprocess
        and the set of link tables connect to it.  If this fails we
        will produce less debris */
        int evid=d.get_int("evid");
        string sta=d.get_string("sta");
        string chan=d.get_string("chan");
        string pchan=d.get_string("pchan");
        double hang=d.get_double("hang");
        double vang=d.get_double("vang");
        /* This is used to mark orientation as measured or computed */
        string otype=d.get_string("orientation_type");
        rec=dbsave(d,dynamic_cast<DatabaseHandle&>(dbhwp),mdlwp,am);
        dbhwp.db.record=rec;
        int pwfid=dbhwp.get_int("pwfid");
        dbhevl.append();
        dbhevl.put("pwfid",pwfid);
        dbhevl.put("evid",evid);
        dbhscl.append();
        dbhscl.put("pwfid",pwfid);
        dbhscl.put("sta",sta);
        dbhscl.put("chan",chan);
        dbhor.append();
        dbhor.put("pwfid",pwfid);
        dbhor.put("pchan",pchan);
        dbhor.put("hang",hang);
        dbhor.put("vang",vang);
        dbhor.put("meastype",otype);
	/* this is an oddity done as a local workaround.
	If found in a release version, consider removing this */
	try {
	  string comment=d.get_string("comment");
	  dbhor.put("auth",comment);
	} catch(...){}; // intentional do nothing
    }
    catch (SeisppError serr){ throw serr;}
    catch (MetadataGetError mderr)
    {
        cerr << "save_scalar_data:  MetadataGet error:";
        mderr.log_error();
        throw SeisppError(string("Database save failure.  Output db is definitely incomplete."));
    }
}
/* Painfully parallel routine for 3c data.  Key difference is does
not map to wfdisc so there is no argument for wfdisc handle.  Similarly
there is no need for the orient table. */
void save_vector_data(ThreeComponentSeismogram& d,
DatascopeHandle& dbhwp,
DatascopeHandle& dbhevl,
DatascopeHandle& dbhscl,
MetadataList& mdlwp,
AttributeMap& am)
{
    try
    {
        string sdtype;
        if(IntelByteOrder())
                sdtype=string("c3");
        else
                sdtype=string("3c");
	d.put("datatype",sdtype);
        int rec=dbsave(d,dbhwp.db,
            string("wfprocess"),mdlwp,am);
        dbhwp.db.record=rec;
        /* We get everything we need before we try to save wfprocess
        and the set of link tables connect to it.  If this fails we
        will produce less debris */
        int evid=d.get_int("evid");
        string sta=d.get_string("sta");
        string chan=sdtype;
        string pchan=sdtype;
        int pwfid=dbhwp.get_int("pwfid");
        dbhevl.append();
        dbhevl.put("pwfid",pwfid);
        dbhevl.put("evid",evid);
        dbhscl.append();
        dbhscl.put("pwfid",pwfid);
        dbhscl.put("sta",sta);
        dbhscl.put("chan",chan);
    }
    catch (SeisppError serr){ throw serr;}
    catch (MetadataGetError mderr)
    {
        cerr << "save_scalar_data:  MetadataGet error:";
        mderr.log_error();
        throw SeisppError(string("Database save failure.  Output db is definitely incomplete."));
    }
}


void save_catalog(DatascopeHandle& dbh, EventCatalog& evcat)
{
    evcat.rewind();
    int i;
    int nrec=evcat.size();
    DatascopeHandle dbhev(dbh);
    dbhev.lookup("event");
    DatascopeHandle dbhorigin(dbh);
    dbhorigin.lookup("origin");

    for(i=0;i<nrec;++i,++evcat)
    {
        Metadata md=evcat.current_aux();
        double mb(-10.0);
        /* put this is a separate try block in
        case it isn't set */
        try
        {
            mb=md.get_double("mb");
        }
        catch(...)
        {
        };

        try
        {
            bool is_new=md.get_bool("is_new_evid");

            if(is_new)
            {
                int evid,orid;
                evid=md.get_int("evid");
                orid=md.get_int("orid");
                Hypocenter h=evcat.current();
                dbaddv(dbhev.db,0,"evid",evid,
                    "prefor",orid,0);
                int iret=dbaddv(dbhorigin.db,0,
                    "lat",deg(h.lat),
                    "lon",deg(h.lon),
                    "time",h.time,
                    "depth",h.z,
		    "evid",evid,
                    "orid",orid,0);
                if(iret==dbINVALID)
                {
                    cerr << "Warning:  dbaddv failure "
                        << "appending new origin rows"
                        <<endl;
                }
                else if(mb>0.0)
                {
                    dbhorigin.db.record=iret;
                    dbputv(dbhorigin.db,0,"mb",mb,0);
                }
            }
        }                                         // handle this silently
        catch (MetadataGetError)
        {
        };
    }
}
void apply_calib(TimeSeries& d)
{
	try {
		double calib;
		calib=d.get_double("calib");
		if(calib!=1.0 & calib>0.0)
		{
			int i;
			for(i=0;i<d.ns;++i) d.s[i]*=calib;
		}
	} catch (MetadataGetError mderr)
	{
		cerr << "Warning:  calib not initialized.  apply_calib does nothing"<<endl;
	}
}
void save_arrivals(Metadata& d,DatascopeHandle& dbh, string Pcomp, string Scomp)
{
	string sta;
	int arid;
	double time;
	int jdate;
	string iphase;
	int ierr;
	try {
		sta=d.get_string("sta");
	} catch (MetadataGetError mderr)
	{
		throw SeisppError(string("save_arrivals:  sta is not defined for data passed\n")
			+"Cannot save arrivals");
	}
	try {
		time=d.get_double("Ptime");
		jdate=yearday(time);
		iphase="P";
		ierr=dbaddv(dbh.db,0,"sta",sta.c_str(),
			"time",time,
			"jdate",jdate,
			"chan",Pcomp.c_str(),
			"iphase",iphase.c_str(),0);
		if(ierr==dbINVALID) cerr << "save_arrivals(WARNING):   "
			<< "dbaddv failed writing P arrival record for station="

			<< sta<<endl;
	} catch(...){};  //Intentionally do nothing here.  Ptime not always defined
	try {
		time=d.get_double("Stime");
		jdate=yearday(time);
		iphase="S";
		ierr=dbaddv(dbh.db,0,"sta",sta.c_str(),
			"time",time,
			"jdate",jdate,
			"chan",Scomp.c_str(),
			"iphase",iphase.c_str(),0);
		if(ierr==dbINVALID) cerr << "save_arrivals(WARNING):   "
			<< "dbaddv failed writing S arrival record for station="

			<< sta<<endl;
	} catch(...){};  //Intentionally do nothing here.  Ptime not always defined
		
}

void usage()
{
    cerr << "sh2db db [-3c -e "
        << "-RFmode"
        << "-pf pffile -v] < file_list"  <<endl
        << "Default all flags false"<<endl;
    exit(-1);
}


bool SEISPP::SEISPP_verbose(false);
int main(int argc, char **argv)
{
    const string qdfext("QBN"), qhdrext("QHD"),threecext("3C");
    const string cont_mess("Trying next input line\n");
    if(argc<2) usage();
    string dbname(argv[1]);
    int i;
    bool threecmode(false);
    bool saveorigin(false);
    bool RFmode(false);
    string pffile("sh2db");
    bool verbose(false);
    for(i=2;i<argc;++i)
    {
        string sarg(argv[i]);
        if(sarg=="-3c")
            threecmode=true;
        else if(sarg=="-e")
        {
            saveorigin=true;
        }
        else if(sarg=="-RFmode")
        {
            RFmode=true;
        }
        else if(sarg=="-pf")
        {
            ++i;
            if(i>=argc) usage();
            pffile=string(argv[i]);
        }
        else if(sarg=="-v")
        {
            SEISPP::SEISPP_verbose=true;
            verbose=true;
        }
        else
            usage();
    }
    if(RFmode)
    	if(!threecmode)
	{
	    cerr << "Warning:  -3c is a required option of -RFmode option"<<endl
	    	<< "Turning 3c mode on and continuing"<<endl;
	    threecmode=true;
        }
    Pf *pf;
    if(pfread(const_cast<char *>(pffile.c_str()),&pf))
    {
        cerr << "pfread failed for pffile="<<pffile<<endl;
        exit(-1);
    }
    char *stmp;
    stmp=pfget_string(pf,"datatype");
    if(stmp==NULL)
    {
        cerr << "Missing required parameter datatype"<<endl;
        exit(-1);
    }
    string datatype(stmp);
    if(datatype!="u4" && datatype!="t4")
    {
        cerr << "Illegal value for datatype parameter = "<<datatype<<endl;
        cerr << "SH Q files are either u4 or t4"<<endl;
    }
    SHHdrMap shhdrmap;
    shhdrmap=load_shhdrmap(pf);
    try
    {
        DatascopeHandle dbh;
        dbh=DatascopeHandle(dbname,false);
        dbh.lookup("wfdisc");
        DatascopeHandle dbhwfp(dbh);
        dbhwfp.lookup("wfprocess");
        DatascopeHandle dbhevl(dbh);
        dbhevl.lookup("evlink");
        DatascopeHandle dbhscl(dbh);
        dbhscl.lookup("sclink");
        DatascopeHandle dbhorient(dbh);
        dbhorient.lookup("orient");
	DatascopeHandle dbharr(dbh);
	dbharr.lookup("arrival");

        /* indent tear here.  Eventually need to fix this with bcpp*/
        char *schema=pfget_string(pf,"schema");
        // We default this silently if it isn't defined.
        if(schema==NULL) schema=strdup("css3.0");
        AttributeMap am;
        am=AttributeMap(string(schema));
        string tag("wfdisc_attribute_list");
        MetadataList mdlwd=pfget_mdlist(pf,tag);
        if(mdlwd.size()<=0)
        {
            cerr << "Parameter file is missing Tbl with tag = "<<tag<<endl;
            cerr << "Check parameter file including setting of waveform_output_table"
                << endl;
            exit(-1);
        }
        /* seriously repetitious, but oh well */
        tag="wfprocess_attribute_list";
        MetadataList mdlwp=pfget_mdlist(pf,tag);
        if(mdlwp.size()<=0)
        {
            cerr << "Parameter file is missing Tbl with tag = "<<tag<<endl;
            cerr << "Check parameter file including setting of waveform_output_table"
                << endl;
            exit(-1);
        }
        tag="eventdb_attribute_list";
        MetadataList mdlev=pfget_mdlist(pf,tag);
        if(mdlev.size()<=0)
        {
            cerr << "Parameter file is missing Tbl with tag = "<<tag<<endl;
            cerr << "Check parameter file including setting of waveform_output_table"
                << endl;
            exit(-1);
        }

        /*This key is used to try to extract inclination angle
        If set to "none" or "NULL" program never attempts to
        fetch this angle and assumes it has to compute it
        Note this MUST match a Metadata key in the SHHdrAttributes
        Tbl in this pf*/
        string ema_key=pfget_string(pf,"ema_key");
        bool use_header_ema;
        if(ema_key=="none" || ema_key=="NULL")
        {
            use_header_ema=false;
            ema_key="";
        }
        else
            use_header_ema=true;
        double v0=pfget_double(pf,"v0");
        /* Set up the handle for saving the event data */
        DatascopeHandle evdbh(dbh);
        evdbh.lookup("event");
        evdbh.natural_join("origin");
        evdbh.subset(string("orid==prefor"));
        EventCatalog evcat(mdlev);
        if(evdbh.number_tuples()>0 )
        {
            cout << "Loading existing catalog with "<<evdbh.number_tuples()
                << " events"<<endl;
            evcat=EventCatalog(evdbh,mdlev,am);

        }
        /* We need this because we cannot guarantee SH stores receiver coordinates
        in a Q file.  We require a db to have a shot at this. */
        char *netstr=pfget_string(pf,"netname");
        char *tsstr=pfget_string(pf,"initial_timestamp");
        double initial_timestamp=str2epoch(tsstr);
        SeismicArray receivers(dynamic_cast<DatabaseHandle&> (dbh),initial_timestamp,
            string(netstr));
        char *s=pfget_string(pf,"ttmethod");
        string ttmethod(s);
        s=pfget_string(pf,"ttmodel");
        string ttmodel(s);
	/* These are used only in rfmode */
	bool time_skewing(false);
	string time_skew_component("L");
	double maximum_time_skew;
	time_skewing=pfget_boolean(pf,"enable_time_skewing");
	time_skew_component=string(pfget_string(pf,"time_zero_component"));
	bool normalize(false);
	normalize=pfget_boolean(pf,"peak_normalization");
	maximum_time_skew=pfget_double(pf,"maximum_time_skew");
	bool save_relative_time_data;
	save_relative_time_data=pfget_boolean(pf,"save_relative_time_data");
	TimeWindow atime_window;
	if(RFmode && time_skewing && save_relative_time_data)
	{
		double ts,te;
		ts=pfget_double(pf,"relative_time_window_start");
		te=pfget_double(pf,"relative_time_window_end");
		atime_window=TimeWindow(ts,te);
	}
	/* these determine which component used to post P and S  arrivals */
	string Pcomp,Scomp;
	Pcomp=pfget_string(pf,"P_arrival_component");
	Scomp=pfget_string(pf,"S_arrival_component");

        /* In the loop below we have to split dir and dfile of the input files.  The
        interface to the plain C function requires a pointer and warnings to not
        free said pointer.*/
        char *dir;
        char *dfilefull;
        char *dfilebase;                          // contains dfilefull stripping ".QBN" at end of name
        string qfdfile,hdrdfile;
        char inbuf[1024],tstr[1024];
        while(cin.getline(inbuf,1024))
        {
            cout << "Processing file="<<inbuf<<endl;
	    /* a copy is required.  man page says so and found out the 
	    hard way why */
	    strcpy(tstr,inbuf);
            dir=dirname(tstr);
            dfilefull=basename(inbuf);
            if(dir==NULL || dfilefull==NULL)
            {
                cerr << "Failure splitting the following input line into dir and dfile:"
                    <<endl
                    << inbuf<<endl;
                cerr << cont_mess;
                continue;
            }
            int lenfull=strlen(dfilefull);
            if(lenfull<5)
            {
                cerr << "Failure splitting the following line to dir and dfile:"
                    <<endl
                    <<"Parsed dfile name="<<dfilefull<< " is not a legal filename"<<endl;
                cerr << cont_mess;
                continue;
            }
            //dfilebase=new char[lenfull-3];  //-3 to allow for \0
            //strncpy(dfilebase,dfilefull,lenfull-4);
            string dfilebase;
            dfilebase.assign(string(dfilefull),0,lenfull-4);
            qfdfile=string(dir)+"/"+string(dfilebase)+"."+qdfext;
            hdrdfile=string(dir)+"/"+string(dfilebase)+"."+qhdrext;
            string dfbsave(dfilebase); // need this below, but better to use a string than char *
            ifstream headerstream;
            try
            {
                headerstream.open(hdrdfile.c_str(),ios::in);
            }
            catch (ios::failure& var)
            {
                string mess;
                mess=var.what();
                cerr << "Open failure on Q file pair="
                    <<hdrdfile<<" and " <<qfdfile<<endl;
                cerr << "System error message:  "
                    << mess<<endl;
                cerr << cont_mess;
                continue;
            }
            vector<string> hdrstrvec = parse_hdr(headerstream);
            headerstream.close();
            int ntraces=hdrstrvec.size();
            if(ntraces<=0)
            {
                cerr << "Failure in procedure parse_hdr.  Return vector is empty"
                    <<endl<< "Check file="<<hdrdfile<<endl;
                cerr << cont_mess;
                continue;
            }
            if(threecmode)
            {
                /* Dogmatically refuse to process a file when ntrace is not a
                multiple of 3 */
                if(ntraces%3 == 0)
                {
                    /* The test of ntraces-2 is safer, but somewhat redundant
                    with the mod 3 conditional */
                    int foff;
                    for(i=0,foff=0;i<(ntraces-2);i+=3)
                    {
                        Metadata mx1=shstr2md(shhdrmap,hdrstrvec[i]);
                        Metadata mx2=shstr2md(shhdrmap,hdrstrvec[i+1]);
                        Metadata mx3=shstr2md(shhdrmap,hdrstrvec[i+2]);
                        string qdfbase=string(dfbsave)+"."+qdfext;
			string threecdfile=string(dfbsave)+"."+threecext;
                        try
                        {
                            set_required(mx1,string(dir),qdfbase,datatype);
                            set_required(mx2,string(dir),qdfbase,datatype);
                            set_required(mx3,string(dir),qdfbase,datatype);
                        } catch (int ierr)
                        {
                            cerr << cont_mess;
                            continue;
                        }
                        vector<TimeSeries> components;
                        /* computer foff needed to read actual data into memory in this mode */
                        int nsamp=mx1.get_int("nsamp");
                        mx1.put("foff",foff);
                        foff+=4*nsamp;
                        /* This may have an off by one error.  It also assumes Q files are 32 bit */
                        mx2.put("foff",foff);
                        foff+=4*nsamp;
                        mx3.put("foff",foff);
                        foff+=4*nsamp;
                        Hypocenter h=handle_event_data(mx1,evcat,dbh.db,ttmethod,ttmodel);
			/* We can assume the are set here */
			int evid,orid;
			evid=mx1.get_int("evid");
			orid=mx1.get_int("orid");
			mx2.put("evid",evid);
			mx2.put("orid",orid);
			mx3.put("evid",evid);
			mx3.put("orid",orid);

                        /* Refresh array geometry if necessary.  No try/catch
                        here as we can assume set_required set these.  */
                        TimeWindow tracerange(mx1.get_double("time"),
                            mx1.get_double("endtime"));
                        try
                        {
                            if(!receivers.GeometryIsValid(tracerange))
                            {
                                receivers=SeismicArray(
                                    dynamic_cast<DatabaseHandle&> (dbh),
                                    tracerange.start,string(netstr));
                            }
			    set_receiver_coordinates(mx1,receivers);
			    set_receiver_coordinates(mx2,receivers);
			    set_receiver_coordinates(mx3,receivers);
                            /* This procedure sets hang and vang
                            or the three components making some
                            rather restrictive assumptions.  See procedure
                            code above or details */
                            set_orientation(mx1,mx2,mx3,h,receivers,ema_key,v0,true);

                            TimeSeries sx1,sx2,sx3;
                            sx1=TimeSeries(mx1,true);
                            components.push_back(sx1);
                            sx2=TimeSeries(mx2,true);
                            components.push_back(sx2);
                            sx3=TimeSeries(mx3,true);
                            components.push_back(sx3);
			    if(RFmode && time_skewing) 
			    {
			    	 double skew;
			    	 skew=RFskewtime(components,time_skew_component,maximum_time_skew,
				 	h,normalize);
				 if(verbose) cout << "station " << sx1.get_string("sta") 
				 		<< " for event "<< sx1.get_int("evid")
						<< " skewed time by "<<skew<<" seconds"<<endl;
			    }
			    else if(RFmode)
			    {
				// only do this if skewing off and RFmode on
				for(int ic=0;ic<3;++ic)
					components[ic].put("calib",1.0);
			    }
			 
                            /* Save these scalar components to wfdisc and
                            wfprocess.  Use just one catch block assuming
                            if one fails they all fail.  Not optimal, but
                            cleaner. */
                            try
                            {
				int icmp;
				for(icmp=0;icmp<3;++icmp)
                                	save_scalar_data(
					 dynamic_cast<Metadata&>(components[icmp]),
                                    	 dbh,dbhwfp,dbhevl,
                                    	 dbhscl,dbhorient,mdlwd,mdlwp,am);
                            } catch (SeisppError serr)
                            {
                                cerr << "Error saving db attributes from"
                                    << " Q file="<< inbuf <<endl
                                    << "SeisppError message:";
                                serr.log_error();
                            }
			    save_arrivals(dynamic_cast<Metadata&>(components[0]),
				dbharr,Pcomp,Scomp);
			    if(normalize)
			    {
				int icmp;
				for(icmp=0;icmp<3;++icmp)
					apply_calib(components[icmp]);
			    }
                            ThreeComponentSeismogram d3c(components);
			    /* This gets a different name, but is placed in same dir */
			    d3c.put("dfile",threecdfile);
			    save_vector_data(d3c,dbhwfp,dbhevl,dbhscl,mdlwp,am);
                            if(verbose)
                                cout << dynamic_cast<Metadata&>(d3c);
			    if(save_relative_time_data)
			    {
				auto_ptr<ThreeComponentSeismogram> d3cr=ArrivalTimeReference(
						d3c,string("Ptime"),atime_window);
				/* This must be reset to tag this as relative time */
				d3cr->put("timetype","r");
				save_vector_data(*d3cr,dbhwfp,dbhevl,dbhscl,mdlwp,am);
			    }
                        } catch (SeisppError serr)
                        {
                            cerr << "Error creating ThreeComponentSeismogram object."<<endl;
                            serr.log_error();
                            cerr << cont_mess;
                            continue;
                        }
                    }
                }
                else
                {
                    cerr << hdrdfile << " and associated Q file="<<qfdfile
                        << " were not processed."<<endl
                        << "Running in three component mode requires number of trace to be "
                        << "a multiple of 3.  Actual ntraces="<<ntraces<<endl;
                }

            }
            else
            {
                /* When not in 3c mode this algorithm is simpler.  We just work through the
                list of header strings, create TimeSeries objects, and write them to the db */
                try
                {
                    string qdfbase=string(dfbsave)+"."+qdfext;
                    int foff;
                    for(i=0,foff=0;i<ntraces;++i)
                    {
                        Metadata mx=shstr2md(shhdrmap,hdrstrvec[i]);
                        set_required(mx,string(dir),qdfbase,datatype);
                        Hypocenter h=handle_event_data(mx,evcat,dbh.db,ttmethod,ttmodel);
                        /* repetitious with above, but judged better
                        to repeat than add overhead of a new call */
                        TimeWindow tracerange(mx.get_double("time"),
                            mx.get_double("endtime"));
                        if(!receivers.GeometryIsValid(tracerange))
                        {
                            receivers=SeismicArray(
                                dynamic_cast<DatabaseHandle&> (dbh),
                                tracerange.start,string(netstr));
                        }
			set_receiver_coordinates(mx,receivers);
                        /* This works because the false arg causes this to
                        not touch arg 2 and 3. */
                        set_orientation(mx,mx,mx,h,receivers,ema_key,v0,false);
                        /* no try block here as we can
                        normally assume this attribute is set*/
                        int nsamp=mx.get_int("nsamp");
                        mx.put("foff",foff);
                        foff+=4*nsamp;
                        save_scalar_data(mx,dbh,dbhwfp,dbhevl,
                            dbhscl,dbhorient,mdlwd,mdlwp,am);
			save_arrivals(mx,dbharr,Pcomp,Scomp);
                    }
                } catch(SeisppError serr)
                {
                    serr.log_error();
                    cerr << cont_mess;
                    continue;
                }
            }
        }
        if(saveorigin)  save_catalog(evdbh,evcat);

    }
    catch(SeisppError serr)
    {
        serr.log_error();
        usage();
    }
}
