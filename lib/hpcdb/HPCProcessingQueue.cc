#include "HPCProcessingQueue.h"
#include "HPCDBHandle.h"  
namespace SEISPP
{
using namespace std;
using namespace SEISPP;

HPCProcessingQueue::HPCProcessingQueue(string dbname,
	list<string> group_keys,
		HeaderMap& hmin)
{
	const string base_error("HPCProcessingQueue initializing constructor:  ");
	/* This constructor always initializes.  Open the queue file in creation mode.*/
	string qfilename=dbname+"."+dbqueuefext;
	fp=fopen(qfilename.c_str(),"w");
	if(fp==NULL)
	    throw SeisppError(base_error+"Cannot open db queue file="+qfilename);
	/* Open the index file and read it.  If it doesn't exist throw an 
	exception */
	string indexfile;
	indexfile=dbname+"."+dbindexfext;
	FILE *fpidx=fopen(indexfile.c_str(),"r");
	if(fpidx==NULL)
		throw SeisppError(base_error+"Cannot open db index file="
			+ qfilename);
}
}
