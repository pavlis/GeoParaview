/* The following are some basic assumptions that permeate this implementation.
1.  data or the db are cachd one tuple at a time.  
2.  In constructors current_tuple is initialized to the first tuple in group and/or the file.
Note this is a guaranteed to be true only if the db is read only.  In write mode this is 
not guaranteed and in fact you should assume the opposite.  If writing the tuple buffer should
be assumed to be initialized to defaults.
3.  The default constructor initializes current_tuple to NULL.  Destructor checks this before releasing 
4.  The FILE pointing at the database table is always open.  It is closed only when the last copy of the 
handle is destroyed.  The file handle is managed by a counter in the private area of this object.
5.  When a tuple is created it is initialzed with a default defined by all constructors.
6.  In read mode current_tuple is loaded in one of two ways. First, as noted above all constructors
must load the first tuple with which it is assoiated (not necessarily row 0 since secondary 
grouping is allowed) on creation.  Second, next loads the tuple.
7.  In write mode current_tuple is assumed to be filled out as needed and likely in multiple places
in a program.  The contents are pushed to the output file in one of two ways.  First, the most common
will be to call the append method.  This will always flush and clear current_tuple.  Second, 
the destructor has to handle the case of dealing with last tuple created.  This is handled
through a private boolean save_me that we have to manage to make sure this works correctly. 
This creates an inefficiency for which I've been unable to figure out an alternative.  
All put methods have to set this attribute true on each and every put call.
8.  All methods that read or write current_tuple must never trust that the file is properly
positioned.  They must ALWAY call fseek before calling fread or fwrite.  Writers need to always
put a lock on the file to avoid collisions in a multiprocessor environment.
9.  fwrite calls are enclosed in flock calls to avoid collisions from multiple processors that
we assume may be writing to the same file simultaneously.  We do this carefully to make sure
write exceptions do not make a multiprocessor run freeze. i.e. we must avoid a fwrite error
exception to cause a lock to be not released.  
*/
/* Regularize method for creating index file name from base (db) name */
string HPCDB_index_file_name(string base)
{
	const string dbindexfext("idx");
	return(base+"."+dbindexfext);
}
string HPCDB_queue_file_name(string base);
{
	const string qfext("q");
	return(base+"."+qfext);
}
/* All constructors use this private method one way or another */
void HPCDBHandle::load_index()
{
	const string base_error("HPCDBHandle::load_index:  ");
	/* Note dbname is assumed to have been set already.  It is a private
	attribute.  Maintenance issue as this evolves. */
	string fname=HPCDB_index_file_name(dbname);
	FILE *fpidx=fopen(fname.c_str(),"r");
	if(fpidx==NULL) throw SeisppError(base_error
		+ "fopen failed on index file = "
		+ fname);
	/* Compute the size of the vector needed by positioning to eof, calling ftell,
	and computing size from number of bytes.  BE WARNED this depends on a simple
	long int structure for the index. */
	fseek(fpidx,SEEK_END,0);
	long fsize=ftell(fpidx);
	int ntuples;
	/* if the index is empty, assume this is to be treated as a simple, ungrouped
	table.  In this case, we can initialize the file_offsets and nrecs vectors
	internally and don't even reference the index file */
	if(fsize<=0)
	{
		this_is_a_group_handle=false;
		/* This is a bit intrisically dangerous, but I'll do it anyway.  
		We assume fp and the HeaderMap are already defined.  We can then
		compute the index directly by size of the data file and use this to
		build a simple index.  The main danger here is that the entire object
		MUST consistently assume fp has to ALWAY be positioned with fseek 
		before any read operation. i.e. no method should ever assume fp
		is properly positioned.  This is ok because of the planned use of 
		this object by the author, but is a very very important maintenance
		issue as evolution could create mysterious errors.  */
		long dfsize;
		fseek(fp,SEEK_END,0);
		dfsize=ftell(fp);
		int recsize=attributes.size();
		ntuples=dfsize/recsize;
		file_offsets.reserve(ntuples);
		long int foff;
		for(int i=0,foff=0;<i<ntuples;++i,foff+=recsize) file_offsets.push_back(foff);
	}
	else
	{
		rewind(fpidx);
		/* The file format has two control values at the beginning */
		if(fread(&nuples,sizeof(int),1,fpidx)!=1)
			throw SeisppError(base_error
				+ "read error on first access of file="
					+fname);
		/* Intentionally don't check for read error here.  If this doesn't work the
		loop below will abort anyway. */
		fread(&this_is_a_group_handle,sizeof(bool),1,fpidx);
		file_offsets.reserve(ntuples);
		if(this_is_a_group_handle) nrecs.reserve(ntuples);
		/* The index file is binary pairs of foff and number of records per group
		parameters.  We load this only if this is known to be a grouped table */
		long foff;
		int grpsize;
		while(fread(&foff,sizeof(long),1,fpidx)==1) 
		{
			file_offsets.push_back(foff);
			if(fread(&grpsize,sizeof(int),1,fpidx)!=1)
				throw SeisppError("load_index:  read error on file="
					+fname);
			if(this_is_a_group_handle) nrecs.push_back(grpsize);
		}
		if(file_offsets.size()!=ntuples)
		{
			stringstream ss;
			ss << base_error << "Index file size mismatch"<<endl
				<<"Expected "<<ntuples<<" grouping definitions.  Found only "
				<< file_offsets.size()<<endl;
			throw SeisppError(ss.str());
		}
	}
	fclose(fpidx);
}
/* This writes the index file.  This is the source of the file read by
by the load_index method immediately above. */
void HPCDBHandle::save_index()
{
	const string base_error("HPCDBHandle::save_index:  ");
	/* Silently return if the index is empty */
	int ntuples=file_offsets.size();
	if(ntuples<=0) return;
	/* Similarly return immediately if this is not defined as a group pointer. 
	This is redundant, but better safe than sorry */
	if(!this_is_a_group_handle)return;
	/* As for the load method we assume dbname is already set. */
	string fname=HPCDB_index_file_name(dbname);
	FILE *fpidx=fopen(fname.c_str(),"w");
	if(fpidx==NULL) throw SeisppError(base_error
		+ " Cannot create index file="
		+ fname);
	fwrite(&ntuples,sizeof(int),1,fpidx)!=1)
		throw SeisppError(base_error
			+ "write failed on first attempt to file="
			+ fname);
	fwrite(&this_is_a_group_handle,sizeof(bool),1,fpidx);
	/* now write the data for this index */
	long int foff;
	int nr;
	int c1,c2;
	for(int i=0;i<ntuples;++i)
	{
		/* Copying to these two temporary variables is a big inefficient
		but this method is only called once per run so am doing it this
		way to avoid the confusing syntax that is necessary otherwise */
		foff=file_offsets[i];
		nr=nrecs[i];
		c1=fwrite(&foff,sizeof(long int),1,fpidx);
		c2=fwrite(&nr,sizeof(int),1,fpidx);
		if((c1!=1) || (c2!1) )
			throw SeisppErrro(base_error
				+ "fwrite error while saving index data");
	}
}
/* This routine builds the index.  IMPORTANT:  dbAll is a negative number defined in
the include file.  As in Datascope this is used a flag to say use the full table,
Note that if nrecs==dbAll the foff field is ignored and immediately set to 0  so the
entire table will be scanned. */
void HPCDBHandle::build_index(long foff, int nr)
{
	const string base_error("HPCDBHandle::build_index:  ");
	if(foff<0) throw SeisppError(base_error+"illegal negative file offset");
	/* For efficiency and because this is a private method this routine assumes the FILE fp 
	is already set and the HeaderMap attributes has already been loaded successfully.  
	Normal practice would be to check this, but in this context we can trust this is true.
	As noted elsewhere in this file, however, this is a maintenance issue so watch out. */
	
	int tuplesize=attributes.size();
	long dfsize;
	if(nrecs==dbAll)
	{
		foff=0;
		fseek(fp,SEEK_END,0);
		dfsize=ftell(fp);
		/* Faster to just compute number of tuples because here we used fixed record lengths.
		Alternative would require calling ftell or break on eof.  This is cleaner but less general.*/
		int ntuples=dfsize/tuplesize;
		/* Sanity check*/
		if(ntuples<=0) throw SeisppError (base_error+"no data in database.\n"
			+ string("Check name passed to HPCDBHandle constructor."));
		/* handle 1 tuple case as a weird special case.  Will not create an index file,
		which is a potential oddity but will do it this way unless it proves to be problematic*/
		if(ntuples==1)
		{
			nrecs.push_back(1);
			file_offsets.push_back(0);
			return;
		}
		rewind(fp);
	}
	else
	{
		ntuples=nrecs;
		fseek(fp,SEEK_SET,foff);
	}
	/* done for efficiency - avoids reallocs */
	nrecs.reserve(ntuples);
	file_offsets.reserve(ntuples);
	/*We assume current_tuple has already been created to hold data from one tuple.
	For the algorithm used here, however, we need to to keep this record plus the previous one.  
	The approach is basically to define a new group whenever one of the group keys
	changes value from the previous */
	unsigned char *previous_tuple=new unsigned char[tuplesize];
	int i;
	if(fread(previous_tuple,sizeof(unsigned char),tuplesize,fp)!=tuplesize) )
		throw SeisppError(base_error + "read error on first tuple in database\n"
			+ string("Check file permissions."));
	list<string>::iterator  gpkeyiter,gpkeyiterend;
	grpkeyiterend=group_keys.end();
	int nrecs_this_group(1);
	long int current_file_offset,start_this_group;
	current_file_offset=foff;
	start_this_group=foff;
	for(i=1;i<ntuples;++i)
	{
		AttributeType rdtype;
		float fval1,fval2;
		double dval1,dval2;
		strin sval1,sval2;
		short int sival1,sival2;
		long int ival1,ival2;
		long long lival1,lival2;
		bool bval1,bval2;
		bool value_changed;
		if(fread(current_tuple,sizeof(unsigned char),tuplesize,fp)!=tuplesize)(
			throw SeisppError (base_error
				+ "read error while reading database file\n");
		for(value_changed=false,
		    grpkeyiter=group_keys.begin();grpkeyiter!=grpkeyiterend;++grpkeyiter)
		{
			rdtype=attributes.dtype();
			switch(rdtype)
			{
			case REAL32:
				fval1=attributes.get<float>(*grpkeyiter,previous_tuple);
				fval2=attributes.get<float>(*grpkeyiter,current_tuple);
				if(fval1!=fval2) value_changed=true;
				break;
			case REAL64:
				dval1=attributes.get<double>(*grpkeyiter,previous_tuple);
				dval2=attributes.get<double>(*grpkeyiter,current_tuple);
				if(dval1!=dval2) value_changed=true;
				break;
			case INT16:
				sival1=attributes.get<short int>(*grpkeyiter,previous_tuple);
				sival2=attributes.get<short int>(*grpkeyiter,current_tuple);
				if(sival1!=sival2) value_changed=true;
				break;
			case INT32:
				ival1=attributes.get<int>(*grpkeyiter,previous_tuple);
				ival2=attributes.get<int>(*grpkeyiter,current_tuple);
				if(ival1!=ival2) value_changed=true;
				break;
			case INT64:
				lival1=attributes.get<long int>(*grpkeyiter,previous_tuple);
				lival2=attributes.get<long int>(*grpkeyiter,current_tuple);
				if(lival1!=lival2) value_changed=true;
				break;
			case STRING:
				sval1=attributes.get_string(*grpkeyiter,previous_tuple);
				sval2=attributes.get_string(*grpkeyiter,current_tuple);
				if(sval1!=sval2) value_changed=true;
				break;
			case BOOL:
				bval1=attributes.get_bool(*grpkeyiter,previous_tuple);
				bval2=attributes.get_bool(*grpkeyiter,current_tuple);
				if(bval1!=bval2) value_changed=true;
				break;
			case HDRINVALID:
			default:
				throw SeisppError(base_error
					+ "group_key="
					+ *grpkeyiter
					+ "is undefined\nCheck db definition file");

			}
			if(value_changed)break;
		}
		current_file_offset+=tuplesize;
		if(value_changed)
		{
			nrecs.push_back(nrecs_this_group);
			file_offsets.push_back(start_this_group);
			nrecs_this_group=0;
			start_this_group=current_file_offset;
		}
		else
		{
			++nrecs_this_group;
		}
		memcpy((void *)previous_tuple,(void *)current_tuple,tuplesize);

	}
}
void HPCDBHandle::load_default(Pf *pf)
{
	/* Note this method assumes the attributes object was created for the current
	handle before this method was called.  */
	size_t tuplesize=attributes.size();
	default_tuple=new unsigned char[tuplesize];
	Tbl *t=pfget_tbl(pf,"Attribute_defaults");
	if(maxtbl(t)<1) throw SeisppError(string("load_default:  ")
		+ "missing required Tbl Attribute_defaults in pf");
	int i;
	for(i=0;i<maxtbl(t);++i)
	{
		char *line;
		string name;
		int ival;
		short int sival;
		long int lival;
		float fal;
		double dval;
		string sval;
		bool bval;
		unsigned char ucval;
		AttributeType rdtype;
		line=(char *)gettbl(t,i);
		stringstream ss(line);
		ss >> name;
		rdtype=attributes.dtype(name);
		switch(rdtype)
		{
		case INT16:
			ss >> sival;
			attributes.put<short>(name,sival,default_tuple);
			break;
		case INT32:
			ss >> ival;
			attributes.put<int>(name,ival,default_tuple);
			break;
		case INT64:
			ss >> lival;
			attributes.put<long>(name,lival,default_tuple);
			break;
		case REAL32:
			ss >> fval;
			attributes.put<float>(name,fval,default_tuple);
			break;
		case REAL64:
			ss >> dval;
			attributes.put<double>(name,dval,default_tuple);
			break;
		case STRING:
			ss >> sval;
			attributes.put_string(name,sval,default_tuple);
			break;
		case BYTE:
			ss >> ucval;
			attributes.put<unsigned char>(name,ucval,default_tuple);
			break;
		case BOOL:
			ss >> bval;
			attributes.put_bool(name,bval,default_tuple);
			break;
		case HDRINVALID:
		default:
			cerr << "load_default(WARNING):  "
			  << "Error in HPCDB pf default definition.\nAttribute="
			  << name
			  << " is not used in this table.  Review pf definitions."
			  <<endl;
			
		}
	}
	freetbl(t,0);
}
HPCDBhandle::HPCDBHandle(string dbname, list<string> grp_keys) : DatabaseHandle(), 
		HPCDBHandle(dbname,true)
{
	/* Make it very clear that when created this way the handle produced must
	be read only */
	read_only=true;
	save_me=false;
	const string base_error("HPCDBHandle(Initialization constructor):  ");
	/* Note this constructor uses common code in a related HPCDBHandle constructor
	called before we get here.  This is assumed to do these things:
	1.  open the db dable (sets fp)
	2.  creates the HeaderMap object.
	It will also somewhat inefficiently create the queue and load the index
	(if present).  This was considered a small cost because this constructor
	is required to be called once and only once on a particular dataset.  */
	try {
		group_keys=grp_keys;
		if(group_keys.size()==0)
			this_is_a_group_handle=false;
		else
			this_is_a_group_handle=true;
		/* Note this assumes earlier constructor initializes the index
		to the trivial form for a nongrouped pointer and we only need
		to build a new one when a new set of group keys is given */
		if(this_is_a_group_handle)
		{
			this->build_index((long int) 0, dbAll);
			/* The above just creates the index file.  For consistency 
			we need to always load it. */
			load_index(dbname);
		}
		current_record=0;
		current_tuple=NULL;
		grpptr=NULL;
		string qfname=HPCDB_queue_file_name(dbname);
		qfp=fopen(qfname.c_str(),"r");
	} catch (SeisppError serr)
	{
		fclose(fp);
		throw serr;
	}
}
HPCDBhandle::HPCDBHandle(string dbname, bool ro) : DatabaseHandle()
{
	const string base_error("HPCDBHandle(standard constructor):  ");
	read_only=ro;
	/* dbname is a base name.  This is the file name containing the db data */
	string dbdfile=dbname+"."+dbdatafext;
	fp=fopen(dbdfile.c_str(),"r");
	if(fp==NULL)
	   throw SeisppError(base_error + "Cannot open db file="+dbdfile);
	/* In read mode check to make sure there is data in this file */
	if(read_only)
	{
		fseek(*fp,0L,SEEK_END);
		if(ftell(fp)==0)
		{
		    fclose(fp);
		    throw SeisppError(base_error+"database file="+dbdfile
		    		+ "is empty");
		}
		rewind(fp);
	}
	save_me=false;
	Pf *pf;
	string pfname=dbname+dbdefext;
	if(pfread(const_cast<char *>(pfname.c_str()),&pf)
		throw SeisppErro(base_error+"pfread failed for file="+pfname);
	try {
		attributes=HeaderMap(pf,dbpftag);
		load_default(pf);
		int tuplesize=attributes.size();
		current_tuple=new unsigned char[tuplesize];
		/* Note this routine sets this_is_a_group_handle.   This is potentially
		error prone, but this can only be established by attempting to access
		an existing index.*/
		load_index(dbname);
		current_record=0;
		grpptr=NULL;
		dbqueue=HPProcessingQueue(dbname);
	} catch (SeisppError serr)
	{
		fclose(fp);
		throw serr;
	}
	/* Watch this one carefully.  Here we create the handle_count and initialize
	it to 1.  */
	handle_count=new int;
	*handle_count=1;
	/* Hopefully clearing this won't create problems */
	pffree(pf);
}
/* Secondary grouping constructor */
HPCDBHandle::HPCDBHandle(HPCDBHandle& parent, list<string> grpkeys) 
	: HeaderMap(parent.attributes), HPCProcessingQueue(parent.dbqueue)
{
	const string base_error("HPCDBHandle secondary grouping constructor:  ");
	/* First we copy attributes in common with the more generic copy
	constructor. i.e these are not special for case of secondary grouping*/
	this_is_group_pointer=parent.this_is_group_pointer;
	read_only=parent.read_only;
	save_me=parent.save_me;
	int tuplesize=attributes.size();
	memcpy((void *)default_tuple,(void *)parent.default_tuple,tuplesize);
	fp=parent.fp;
	/* these are ones that have to be handled specially here.  NEEDED:  copy `
	when finished. 
	current_tuple blah blah
	*/
	try {
		/* Now we have to build a new index from this group using these new keys.  */
		group_keys=grpkeys;
		this->build_index(parent.file_offset[current_record],nrecs[current_record]);
		/* seek to the first record and load it to initialize */
		if(fseek(fp,file_offset[0],SEED_SET))
		{
			stringstream ss;
			ss << base_error << "fseek failed to offset = "
				<< file_offset[0];
			throw SeisppError(ss.str());
		}
		if(fread(current_tuple,sizeof(unsigned char),tuplesize,fp)!=tuplesize)
		{
			stringstream ss;
			ss << base_error << "fread failed loading current tuple at file offset="
				<, file_offset[0];
			throw SeisppErro(ss.str());
		}
	} catch (...) {throw;};
	current_record=0;
	/* WARNING:  significant maintenance issue.  In the present form the primary 
	purpose of the handle_count attribute is to manage the the open file handle.
	The index attributes can and will be different when this constructor is called
	because secondary grouping produces a rather different handle than the one
	pointing at the full table.  If this approach to the implementation changes
	this WILL require a change. */
	/* Don't think this is needed because the copy constructor increments this and
	constructor adds onto the copy not creating yet another handle */
	//++(*handle_count);
}
HPCDBHandle::~HPCDBHandle()
{
	--(*handle_count);
	if(handle_count<=0)
	{
		if(save_me)
		{
		    try {
			write_current();
		    } catch(SeisppError serr)
		    {
		    	cerr << "Warning:  problems writing last record in closing"
				<<endl
				<< "write_current method posted this message:"<<endl
			serr.log_error();
			cerr << "output db table may be corrupted"<<endl;
		    }
		}
		fclose(fp);
		delete handle_count;
		if(!read_only) fclose(qfp);
	}
	if(current_tuple!=NULL) delete [] current_tuple;
	delete [] default_tuple;
}
HPCDBHandle::HPCDBHandle(const HPCDBHandle& parent)
	: HeaderMap(dbh.attributes), HPCProcessingQueue(dbh.dbqueue)
{
	this_is_group_pointer=parent.this_is_group_pointer;
	read_only=parent.read_only;
	save_me=parent.save_me;
	current_record=parent.current_record;
	group_keys=parent.group_keys;
	int tuplesize=attributes.size();
	if(parent.current_tuple!=NULL) 
	{
		current_tuple=new unsigned char[tuplesize];
		memcpy((void *)current_tuple,(void *)parent.current_tuple,tuplesize);
	}
	memcpy((void *)default_tuple,(void *)parent.default_tuple,tuplesize);
	file_offsets=parent.file_offsets;
	nrecs=parent.nrecs;
	fp=parent.fp;
	++(*handle_count);
}
int HPCDBHandle::number_tuples()
{
	return(file_offsets.size());
}
double HPCDBHandle::get_double(string name)
{
	AttributeType rdtype=attributes.dtype(name);
	double result;
	float fval;
	switch(rdtype)
	{
	case REAL32:
		fval=this->attributes.get<float>(name,current_tuple);
		result=static_cast<double>(vfal);
		break;
	case REAL64:
		result=this->attributes.get<double>(name,current_tuple);
		break;
	default:
		throw SeisppError(string("HPCDBHandle::get_double: failed trying to fetch attribute=")
			+ name);
	}
	return(result);
}
string HPCDBHandle::get_string(string name)
{
	AttributeType rdtype=attributes.dtype(name);
	string result;
	switch(rdtype)
	{
	case STRING:
		result=this->attributes.get_string(name,current_tuple);
		break;
	default:
		throw SeisppError(string("HPCDBHandle::get_string: failed trying to fetch attribute=")
			+ name);
	}
	return(result);
}
int HPCDBHandle::get_int(int name)
{
	AttributeType rdtype=attributes.dtype(name);
	int result;
	short int sival;
	long int lival;
	/* This assumes the LP64 model for int in a 64 bit environment */
	switch(rdtype)
	{
	case INT16:
		sival=this->attributes.get<short int>(name,current_tuple);
		result=static_cast<int>(sival);
		break;
	case INT64:
		lival=this->attributes.get<long>(name,current_tuple);
		result=static_cast<double>(lival);
		break;
	case INT32:
		result=this->attributes.get<int>(name,current_tuple);
		break;
	default:
		throw SeisppError(string("HPCDBHandle::get_double: failed trying to fetch attribute=")
			+ name);
	}
	return(result);
}
bool HPCDBHandle::get_bool(string name)
{
	AttributeType rdtype=attributes.dtype(name);
	bool result;
	switch(rdtype)
	{
	case BOOL:
		result=this->attributes.get_bool(name,current_tuple);
		break;
	default:
		throw SeisppError(string("HPCDBHandle::get_bool: failed trying to fetch attribute=")
			+ name);
	}
	return(result);
}
void HPCDBHandle::put(string name,const char *s)
{
	this->attributes.put_string(name,string(s),current_tuple);
}
void HPCDBHandle::put(string name,char *s)
{
	this->attributes.put_string(name,string(s),current_tuple);
}
void HPCDBHandle::put(string name,string s)
{
	this->attributes.put_string(name,s,current_tuple);
}
void HPCDBHandle::put(string name,int value)
{
	this->attribures.put<int>(name,value,current_tuple);
}
void HPCDBHandle::put(string name,float value)
{
	this->attribures.put<float>(name,value,current_tuple);
}
void HPCDBHandle::put(string name,double value)
{
	this->attribures.put<double>(name,value,current_tuple);
}
void HPCDBHandle::put(string name,bool value)
{
	this->attribures.put_bool(name,value,current_tuple);
}
void HPCDBHandle::put(Metadata& md, MetadataList& mdl, AttributeMap& am)
{
    const string base_error("HPCDBHandle::put:  ");
    MetadataList::iterator mdli;
    map<string,AttributeProperties>::iterator ami,amie=am.attributes.end();
    map<string,AttributeProperties> aliasmap;
    /* This loop is essentially identical to a similar loop in dbsave for Metadata to 
    a DatascopeHandle */
    for(mdli=mdl.begin();mdli!=mdl.end();++mdli)
    {
        double dval;
        int ival;
	string sval;
	bool bval;
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
                throw SeisppError(base_error + "attribute "
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
                    ival = md.get_int(mdkey);
		    attributes.put<int>(second.db_attribute_name,ival,current_tuple);
                    break;
                case MDreal:
                    dval = md.get_double(mdkey);
		    attributes.put<double>(second.db_attribute_name,dval,current_tuple);
                    break;
                case MDstring:
                    sval = md.get_string(mdkey);
		    attributes.put_string(second.db_attribute_name,sval,current_tuple);
                    break;
                case MDboolean:
                    bval = md.get_bool(mdkey);
		    attributes.put_bool(second.db_attribute_name,bval,current_tuple);
                    break;

                case MDinvalid:
                    cerr << base_error << "database attribute "
                        << ami->second.db_attribute_name
                        << " was marked as invalid\n"
                        << "Data for this attribute not saved"
                        << endl;
                    break;

                default:
                    cerr << base_error << "database attribute "
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
}
	
int HPCDBHandle::next()
{
	const string base_error("HPCDBHandle::next() method:  ");
	if(read_only)
	{
		if(this_is_a_group_handle)
		{
		}
		else
		{
		}
	}
	/* note if save_me is false nothing gets written in write mode.
	This is a safety value to avoid null records. */
	else if(save_me)
	{
		try {
			this->append();
		} catch SeisppError (serr)
		{
			throw serr;
		};
	}
}
int HPCDBHandle::append()
{
	const string base_error("HPCDBHandle::append:  ");
	if(read_only)
		throw SeisppError(base_error+"attempting to write to a read only db");
	try {
		if(save_me) this->write_current();
		size_t tuplesize=attributes.size();
		save_me=false;
		memcpy(current_tuple,default_tuple,tuplesize);
	}
}
HPCDBHandle& HPCDBHandle::current()
{
}
HPCDBHandle& HPCDBHandle::operator=(const HPCDBHandle& parent)
{
}
void HPCDBHandle::operator++()
{
}
void HPCDBHandle::write_current()
{
	const string base_error("HPCDBHandle::write_current:  ");
	size_t nbytes=attributes.size();
	int fd=fileno(fp);
	if(lockf(fd,F_LOCK,(off_t)0)
		throw SeisppError(base_error+"Could not lock db table file");
	if(fseek(fp,0L,SEEK_END))
	{
		lockf(fd,F_ULOCK,(off_t)0);
		throw SeisppError(base_error+"fseek to eof failed");
	}
	if(fwrite((void *)current_tuple,sizeof(unsigned char),nbytes,fp)
		!= nbytes)
	{
		lockf(fd,F_ULOCK,(off_t)0);
		throw SeisppErro(base_error+"fwrite to output data file failed");
	}
	lockf(fd,F_ULOCK,(off_t)0);
}
