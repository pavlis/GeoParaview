#ifndef _HPC_PROCESSINGQUEUE_H_
#define _HPC_PROCESSINGQUEUE_H_

#include <stdio.h>
#include <string>
#include <vector>
#include "ProcessingQueue.h"
#include "HeaderMap.h"
namespace SEISPP
{
using namespace std;
using namespace SEISPP;
class HPCProcessingQueue : public ProcessingQueue
{
public: 
	/*! \brief Construct in an initialization mode.

	This should be viewed as an initialization constructor.  It will
	build and initialize the queue.  It should be called once and 
	only once for processing of a particular database.  Normal mode
	is to open a database, build a group index, and initialize the 
	queue system (currently a file).  There are, however, two special
	cases that are important to recognize.  First, if database is empty
	the constructor immediately returns doing nothing.  This is done 
	to allow using this object seamlessly on a new table intended 
	completely for output.  In this state the object is initialized
	internally to a state to prevent incorrect downstream usage.
	The second special case is if the group_keys list is empty.
	In that situation it is assumed the queue is to be constructed
	without grouping (effectively one tuple per group). */
	HPCProcessingQueue(string dbname,list<string> group_keys,HeaderMap& hmin);
	//HPCDBBundle next(HPCDBHandle& dbh);
	int next(DatabaseHandle& dbhbase);
	HPCDBBundle current();
	void mark(ProcessingStatus ps);
	ProcessingStatus status();
	bool has_data();
private:
	FILE *fp;
	/* This boolean is set false when the table is empty on 
	creation.  This is a common situation in write mode. */
	bool is_valid; 
};
}
#endif
