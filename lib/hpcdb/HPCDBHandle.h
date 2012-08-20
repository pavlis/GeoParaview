#ifndef _HPCDB_H_
#define _HPCDB_H_

#include <string>
#include <list>
#include <vector>
#include "DatabaseHandle.h"
#include "Metadata.h"
#include "AttributeMap.h"
#include "HeaderMap.h"
#include "ProcessingQueue.h"
namespace SEISPP
{

typedef struct {
	int begin;
	int end;
} HPCDBBundle;
using namespace std;
using namespace SEISPP;
/*! Used to define a full table.  Similar to datascope constant with the same name. */
const int dbAll(-999);
/*! Constant extension used or hpcdb data table. */
const string dbdatafext("db");
/*! Constant extension used or hpcdb data table. */
const string dbindexfext("idx");
/*! File containing queue data to mark state of processing of each dataset */
const string dbqueuefext("Q");
/*! File defining table layouts.  Currently an Antelope parameter file. */
const string dbdefext("pf");
/*! tag on pf used to build HeaderMap */
const string dbpftab("RecordLayout");
/*! \brief Database handle for a simple database for high performance computing applcations.

There are a class of data driven processing algorithms that are well suited to massively
parallel computing.  They are of a class sometimes called "embarassingly parallel" because
the scaling to multiple processors is close to linear.  The class of algorithms of note
here are best thought of as table driven.  If the problem can be reduced to reading 
one or more rows of a table to drive reasonably balanced processing of many similar
data sets this approach makes sense.  The underlying concept is a simple table.  
The only allowed database concept is "group by".  A method exists in the handle to 
group by specific attribute keys.  The idea behind this is ensemble processing.
That is it is assumed the algorithm can be driven by single instances of a common data object or
ensembles (groups) of the same thing with some common attributes.  A classic example
in reflection processing are shot gathers (group by shot id), CDP gathers (group by
CDP bin), and common offset (group by offset bins).  One might want to generalize
the group concept through generic templates that supply compare operators like those
use in sorting in the STL, but I viewed that as unnecessarily complicated due to 
the total flexibility this object provides in defining what is in a table.
That is, in the worst case if you want to group by some odd compare method do 
this externally and define an auxiliary integer key for internal use to characterize
the grouping.  

The definition of the table that drives the process is generic.  The description 
of layout of the table is specified by a file derived automatically from the 
database name.  This is a bit like Antelope, but without the concept of multiple
tables.  That is, if we have a database with a name mydb, the actual data for 
the table is assumed to be in a file called mydb.dat.  The dictionary to define
attributes in mydb will be assumed to be present in a file called mydb.pf.  
This name is what it is because the current implementation uses an Antelope
parameter file (pf) to define the table layout.  It is expected that this might
evolve someday to an xml format that would be transparent.  That is, if the 
dictionary were xml we would expect mydb.xml instead of mydb.pf.

This object works in a parallel computing environment.  The current implementation
requires every processor to see a common file system AND all be a common architecture.
The interface, however, is totally generic and should not depend on this choice. 
Upgrades to full grid computing on multiple platforms would require huge modifications
but I aimed to produce an interface that was extensible to such an environment. 
Undoubtedly this is not likely to work, but I've given it my best shot.  Call it
a prototype interface. */
class HPCDBHandle : public DatabaseHandle : public ProcessingQueue
{
public:
	/*! Switch to define if this is or is not a group handle.

	This object can define a grouped table in reading (not allowed for
	writing).  This public boolean is used to test for this condition.*/
	bool this_is_a_group_handle;
	/*! \brief Create a new handle for initial processing setup.

	This is a one step creation method that both builds a handle
	and initializes the db for processing.  It should be called
	once and only once to initialize processing on a new dataset.

	This constructor looks for 2 files linked to the 
	input parameter dbname.  dbname.db is a binary file
	that contains the processing table data that contains
	the data that define this "database" and  dbname.pf is
	a parameter file that provides the dictionary of mapping
	of the table data to named attributes.  This constructor
	will ALWAYS build a third file with the name dbname.idx.  
	This file is an index that defines grouping of the table
	using the attributes passed through the group_keys argument.
	If group_keys is empty the the assumption is that data will be 
	processed one row at a time and the index reduces to 
	one entry per row of the database.  Otherwise the size of
	the index is smaller than the number of rows in the 
	database table.  

	Note the pf file contains the dictionary that describes the
	table layout.  Internally this uses the HeaderMap object.
	See documentation for the HeaderMap to see how the pf should
	be built.  The HeaderMap section of the pf file are encased
	in the (frozen) tag "RecordLayout".  This constructor will
	modify the pf file to add a string entry "group_keys" that 
	is a blanks seperated list of the group_keys argument passed
	to this constructor.  
	
	If the above is successful this constructor initializes an
	HPCProcessingQueue object.  In the current implementation this
	creates a local file with the name dbname.Q.  This file contains
	two things used in subsequent processing. First, when grouping
	is used it is a cache of the computed bundle structure.  Second,
	it contains a list of processing status indicators used to mark
	a component of the data in progress, completed, or left to process.*/ 
	HPCDBHandle(string dbname,list<string> group_keys);
	/*! \Create a secondary group bundled handle.

	This constructor is intended to provide a way to produce
	a higher level grouping of an existing bundled handle.
	It will abort if applied to an existing handle that is not
	already a group (bundled) handle.  This is necessary because
	in this interface the lowest level group is always under
	control of the system used to manage multiple processor 
	negotiation for the next block of data.  In contrast, the
	result of this constructor is NOT under control of this 
	system and can only be incremented one group record at
	a time.  A type example in seismology for this approach
	is a three-component ensemble.  To reference all channels
	requires a double grouping:  one by station and the other
	by whatever grouping method is used to define the ensemble.  

	\param dbh is the parent handle used to build this secondary
		grouped handle.  This must be an already grouped 
		handle or this constructor will throw an exception.
	\param group_keys defines the keys used for grouping.
	*/
	HPCDBHandle(HPCDBHandle& dbh, list<string> group_keys);
	/*! \brief Constructor to open partially processed dataset. 

	This constructor should be used by most processes to create this 
	handle.  It assumes the initialization constructor was called 
	earlier and the three files that define this db structure have
	already been created.  If any are lacking this constructor will
	throw a SeisppError exception.  

	This constructor abstracts this basic idea.  When successfully created in
	a multiprocessor environment the handle can be assumed to be position to
	the next dataset defined by the handle.  The negotiations for what defines
	next in a multiprocessor environment is considered to be implemenation
	dependent.  

	Finally, the read_only boolean has a very critical role here.  This
	handle behaves totally different in read and write mode.  The above
	describes read_only mode.  When set false read_only really implies write_only.
	The perspective of the author in creating this interface is that in the
	current HPC environment characterized by clusters and grid computing 
	database updates should not be allowed.  If that is needed, extend this
	object by an inheritance mechanism.  The goal of this interface is to 
	maximize performance on a class of algorithms that can be reduced to 
	a simple data-driven algorithm that could be reduced to this:
	while(read_data_ok()); do_something(); write_result(); endwhile;.
	read_only mode should be used in the read_data_ok to read next dataset.
	read_only false should be the mode for the write_result abstraction.
	*/
	HPCDBHandle(string dname, bool read_only);

	HPCDBHandle(const HPCDBHandle& parent);
	int number_tuples();
	double get_double(string name);
	int get_int(string name);
	string get_string(string name);
	string get_filename(string name);
	bool get_bool(string name);
	template <class T> T get(string name);
	template <class T> void put(string name,T value);
	void put(string name, const char *s);
	void put(string name, char *s);
	void put(string name, string s);
	void put(string name, int value);
	void put(string name, float value);
	void put(string name, double value);
	void put(string name, bool value);
	/*! \brief List driven push of multiple attributes. */
	void put(Metadata& md, MetadataList& mdl, AttributeMap& am);
	/*! \brief Position the handle to the next record. 

	This method positions the handle to the next record in the simple db table it
	references.  This can be a group or a single row.  It is more complicated than
	simply incrementing the pointer because in  multiprocessor environment the next
	row is more than likely already being handled by another processor.  
	When this method is successful the current record (our group depending on context)
	pointer is set to a valid value = integer index returned.  If it fails the 
	position should be considered undefined and the handle should be considered
	no longer usable.

	In read mode involking this method causes the next available db row to be loaded.
	In write mode this causes two events.  First, the contents of the 
	current record are flushed to output (method is not specified.
	It could be through a file with locking or through a dbwriter using mpi).  
	Then a new record is created with defaults loaded.  (like dbaddnull
	in Datascope).  Under construction.  Several ways this could 
	create an exception.

	\return record number if there are data remaining to be processed.  Returns
	a -1 if the processing queue linked to this object has no more data 
	*/
	int next();
	/* Implemented to mesh with virtual DatabaseHandle. Should throw an
	exception if used in read_only mode. In write mode is should be 
	an alias for next()*/
	int append();
	/*! \brief Get handle to current group.  
	
	This is a core method.  When the object defines a bundled group the current
	record counter is group number, not record number.  This drops into the 
	heirarchy and returns a handle that can be addressed at the record level.
	
	\exception SeisppError is thrown if the parent is not a handle to a grouped
	table or if the handle is empty (an example of this is one attempts to 
	access data beyond the end of the processing table)*/
	HPCDBHandle& current();
	/* These are implementations of virtual methods in ProcessingQueue*/
	void mark(ProcessingStatus ps);
	ProcessingStatus status();
	bool has_data();
	HPCDBHandle& operator=(const HPCDBHandle& parent);
	/*! \brief Position the handle to the next tuple. 

	This operator is effectively an alias for the next() method.  It can
	be uses if one does not need the integer returned by the next method.
	*/
	void operator++();  
private:
	/*! database name ala Datascope (root with no extension) */
	string dbname;  
	/*! Set true if read_only mode.  */
	bool read_only;
	/*! Flags if current_tuple has data that needs saving.

	In write mode this variable is necessary mainly to handle proper saving of the
	last tuple.  The destructor has to be told if the current_tuple needs to be 
	pushed to output.  All put methods set this attribute true.
	*/
	bool save_me;
	HeaderMap attributes;
	int current_record;
	list<string> group_keys;
	/* This is always a pointer to the data in the current record (one row) */
	unsigned char *current_tuple;
	/* Each copy of the handle keeps a copy of a null record created from the
	definition in the pf.  This holds that data.*/
	unsigned char *default_tuple;
	/* This contains a list of offsets defining the starting positions 
	(file offsets) of each row in a group.  This can index the entire
	table or a more limited range of rows for a secondary grouping. */
	vector<long int> file_offsets;
	/* This is a parallel vector to file_offsets containing number of
	records per group.  Purists might prefer a class, but this allows
	this to be empty for ungrouped tables. */
	vector<int> nrecs;
	/* Open file holding the database table. */
	FILE *fp;
	/* This int ptr is used to manage copy count.  When it's contents drop to 
	zero the database table is closed and the content of this ptr is released.
	Otherwise the destructor decrements this. */
	int *handle_count;
	/* queue file pointer */
	FILE *qfp;
	/* This private method loads the index file.  Necessary because the
	index was intentionally made a private attribute to make the 
	interface independent of this implemenation detail.*/
	void load_index(string dbn);
	/* related private method to build a group index and write it to disk. */
	void build_index_file();
	/* This method creates and loads default values in default_tuple*/
	void load_default(Pf *pf);
	/* saves the current tuple */
	void write_current();
};

}
#endif
