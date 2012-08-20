#ifndef _PROCESSINGQUEUE_H_
#define _PROCESSINGQUEUE_H_

#include "dbpp.h"
namespace SEISPP
{
/*! Used in ProcessingQueue objects to mark state of each record to be processed.

Any ProcessingQueue object can mark a record as one of these items.
As the names sugested, TODO means processing has not yet been completed,
FINISHED means processing is completed and was successful, SKIPPED means
the records was skipped for some reason, and PROCESSING is used to 
mark records locked and being processed by this or another process.
*/
enum ProcessingStatus {TODO, FINISHED, SKIPPED, PROCESSING};
/*! Virtual base class for objects used to implement the concept of a processing queue.

*/
class ProcessingQueue
{
public:
	virtual void mark(ProcessingStatus ps)=0;
	virtual ProcessingStatus status()=0;
	virtual bool has_data()=0;
	/*! Position queue to next record.

	A common concept is the need to position the queue to the
	next record in the database.  This defines this concept
	in a generic way.  The assumption of the interface is that
	when this method is called a position pointer in the 
	(generic) DatabaseHandle is incremented to the next record
	the queue says is available for processing.  The implementation
	should then mark that tuple as being processed (PROCESSING)
	and return the record number of that tuple.  This is 
	essential in a multiprocessor enviroment because the 
	next tuple is rarely +1 from the current one.  

	Note an implementation may find it necessary to overload
	this method, but it should retain this basic concept if
	it becomes a subclass of this virtual base.  The main 
	reason for overloading would be that one might choose
	to provide a copy of the handle in a private area of
	the object or use a handle that is not a subclass of
	the virtual class DatabaseHandle. 

	\return A concrete implementation should return a negative
	number if there is no more data to process.  Otherwise if
	the implementation finds this useful it should return a record
	number.  
	*/
	virtual int next(DatabaseHandle& dbh)=0;
};
}
#endif
