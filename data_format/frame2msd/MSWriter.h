#include <stdio.h>
#include <stdlib.h>
#include <typeinfo>
#include <string>
#include <stddef.h>
#include "libmseed.h"
using namespace std;
/* \brief Error class thrown for errors in this package.

Simple child of stdlib abstract exception class.  Handler should catch a
generic exception and call the what() method to print the message.
*/
class MSWexception : public exception
{
public:
    MSWexception(string mess){message=mess;};
    MSWexception(const char *s){message=string(s);};
    virtual const char* what() const throw()
   {
    return message.c_str();
   };
   /* This is necessary baggage to mesh with std::exception declaration
      for this virtual method.  Not required on all compilers */
   ~MSWexception() throw ()
   {};
private:
    string message;
};
static void record_handler(char *record, int reclen, void *prvdata)
{
    size_t retcode;
    FILE *fp=reinterpret_cast<FILE *>(prvdata);
    retcode=fwrite(record,reclen,1,fp);
    if(retcode<=0) throw MSWexception(string("MSWriter::record handler:  ")
            + "fwrite error on miniseed output file");
};

/* \brief Miniseed writing interface.

   There are multiple implementations of miniseed writers. This object
   is designed to abstract the interface so one can, for example, use
   antelope or the libmseed library from IRIS with the same interface. 

   Note this interface INTENTIONALLY does NOT contain a copy constructor.
   The reason is the writer is expected to normally hold an open 
   file handle that will be closed when the desctructor is called.  
   A copy doesn't make much sense in that context and is best avoided
   to avoid chaos.

   This implementation is as general a would be desired for a completely
   generic implementation.  Several potential variables in miniseed
   are frozen.  These include:  (1) byteorder is always big endian
   (chosen because libmseed documentation warns little endian is 
   potentially problematic and readers all handle this automatically); 
   and (2) encoding is frozen for each allowed sample data type.  
   For the later only four data types are accepted:  (1) 16 bit
   ints (declared int32_t or short); (2) 32 bit ints (declared int32_t
   or int); (3) float (32 bit real); or (4) double (64bit real).  
   These are written (respectively) with libmseed encoding:  DE_STEIM2, 
   DE_STEIM2, DE_FLOAT32, and DE_FLOAT64.  

   Finally a very important limitation of this implementation is that 
   because of the internal structure of libmseed this object opens 
   it's output file on creation and does not close that file until
   the destructor is called.  The best way to understand this is to
   just always view an MSWriter as a specialized file handle.  In the
   same way you take care to close files when no longer needed, destroy
   MSWriter objects when they are no longer needed.
   */
class MSWriter
{
public:
    /*! \brief Simple constructor.

       This is only constructor.  Creates a handle
       that allows writing data.  The file is openned for writing and
       remains open until destruction.  

       \param fname  is the file name to which the miniseed data should
          be written.
       \param reclen is the SEED logical record size to use for this file. 
            (defaults to 4096 - highly recommended) 
       \exception MSWexception is thrown if fname open fails.
          */
    MSWriter(string fname, int reclen=4096);
    /*! \brief Standard destructor. 
     
     The destructor here is nontrivial as it closes the output file.*/
    ~MSWriter();  
    /*! \brief Writes a vector of data holding a unique net:sta:chan:loc
       to output miniseed volume.

       This is the main method of this object.  It writes a vector of
       data of any allowed type to the output miniseed file.

       \param d is the data vector assumed to be of specified type
       \param dt is sample interval in seconds
       \param nsamp is number of samples for d
       \param t0 is start time of first sample in posix standard epoch time
       \param net is the network code for this data channel
       \param sta is the station code for this data channel
       \param chan is the channel code to be used for this data channel
       \param loc is the seed location code for this channe
         (an empty string is used when no loc code should be set)

       \return number of miniseed records written

       \exception throws MSWexception object for any problem. 
         Handler should display the what() method output for explanation.
         */
    template <class T> int write(T *d, double dt, int nsamp, double t0,
            string net, string sta, string chan, string loc);
    /*! \brief Simplified writer for multiple data segments.

       Sometimes you can assume you are working on the same data channel
       in a single loop.  This is a convenience routine that writes the
       next block of data using the same net,sta,chan,loc as the previous
       call.  

       \param d is the data vector assumed to be of specified type
       \param dt is sample interval in seconds
       \param nsamp is number of samples for d
       \param t0 is start time of first sample in posix standard epoch time

       \return number of miniseed records written

       \exception throws MSWexception object for any problem. 
         Handler should display the what() method output for explanation.

       */
    template <class T>int write(T *d, double dt, int nsamp, double t0);
    /*! \brief Return current file size. */
    long current_filesize();
private:
    FILE *MSWfp;
    int32_t record_length;
    string net,sta,chan,loc;
};
/* Implementation of generic write method which assumes net:sta:chan:loc 
   are the same as a previous call.  Used when appending. */
template <class T> int MSWriter::write(T *d, double dt, int nsamp, double t0)
{
    const string base_error("MSWriter::write:  ");
    /* Require net, sta, and chan to be set.  Allow loc to be null */
    if(net.length()<0 || sta.length()<0 || chan.length()<0) 
        throw MSWexception(base_error + "net sta and chan must all be set");
    /* This section closely follows mst_pack example */
    int32_t npackets;
    int32_t precords;
    MSTrace *mst;
    mst=mst_init(NULL);
    strcpy(mst->network,net.c_str());
    strcpy(mst->station,sta.c_str());
    strcpy(mst->channel,chan.c_str());
    strcpy(mst->network,net.c_str());
    /* We could probably use a long int here, but this alias is used in
       libmseed for epoch time in a microsecond format. Our input is
       an epoch time stored as a double.  We have to convert it. */
    hptime_t mststart;
    hptime_t mstbase,mstremainder;
    mstbase=static_cast<hptime_t>(t0);
    double tr=t0-static_cast<double>(mstbase);
    mstremainder=static_cast<hptime_t>(1000000*tr);;
    mstbase *= 1000000;
    mststart=mstbase+mstremainder;
    mst->starttime = mststart;
    mst->samprate = dt;
    mst->numsamples = nsamp;
    mst->samplecnt  = nsamp;
    /* Incredibly the libmseed routine calls free on the datasamples
       pointer.  Hence we have to copy d and let the pack routine
       free it.  Also has to use plain C allocs.*/
    T *d2;
    if((d2=(T*)calloc(nsamp,sizeof(T)))==NULL) throw MSWexception(base_error
            +"  calloc failed for internal buffer");
    for(int j=0;j<nsamp;++j) d2[j]=d[j];
    mst->datasamples = static_cast<void *>(d2);
    /* static seems necessary because this is in a template.  
    50 is longer than needed but what is recommended in libmseed doc */
    static char srcname[50];  
    mst_srcname(mst,srcname,0);
    /* Not the best style for C++, but necessary evil here.  This 
       is branch for acceptable types.  Throw an exception if no
       match.  Note all integer formats are automatically steim
       compressed while floats and doubles are just written. */
    if(typeid(T)==typeid(int32_t) || (typeid(T)==typeid(int)))
    {
        mst->sampletype = 'i';
        try{
            precords=mst_pack(mst,&record_handler,static_cast<void *>(MSWfp),
                    record_length,DE_STEIM2,1,&npackets,1,0,NULL);
        }
        catch (...) {
            mst_free(&mst);
            throw;
        }
    }
    else if(typeid(T)==typeid(int16_t) || (typeid(T)==typeid(short int)))
    {
        mst->sampletype = 'i';
        /* libmseed does not support int16_t directly. We will promote
           16s to 32s without any space issue because we'll write the
           output as steim compressed packets */
        int32_t *d32;
        d32=static_cast<int32_t *>(calloc(nsamp,sizeof(int32_t)));
        if(d32==NULL) throw MSWexception(base_error
                            +"  calloc failed for internal buffer");
        for(int i=0;i<nsamp;++i) d32[i]=static_cast<int32_t>(d[i]);
        mst->datasamples=static_cast<void *>(d32);
        try{
            precords=mst_pack(mst,&record_handler,
                    static_cast<void *>(MSWfp),record_length,DE_STEIM2,
                1,&npackets,1,0,NULL);
        }catch (...) {
            mst_free(&mst);
            throw;
        }
    }
    else if(typeid(T)==typeid(long))
    {
        mst_free(&mst);
        throw MSWexception(base_error
                + "SEED does not allow long int sample data");
    }
    else if(typeid(T)==typeid(float))
    {
        mst->sampletype = 'f';
        try{
            precords=mst_pack(mst,&record_handler,
                static_cast<void *>(MSWfp),record_length,DE_FLOAT32,
                1,&npackets,1,0,NULL);
        }catch (...){
            mst_free(&mst);
            throw;
        }
    }
    else if(typeid(T)==typeid(double))
    {
        mst->sampletype = 'd';
        try{
            precords=mst_pack(mst,&record_handler,
                static_cast<void *>(MSWfp),record_length,DE_FLOAT64,
                1,&npackets,1,0,NULL);
        }catch (...) {
            mst_free(&mst);
            throw;
        }
    }
    else
    {
        mst_free(&mst);
        throw MSWexception(base_error
                +"Coding error.  Unsupported sample type.");
    }
    mst_free(&mst);
    if(precords<0) throw MSWexception(base_error 
            + "mst_pack returned -1 error code");
    return(npackets);
}
template <class T> int MSWriter::write(T *d, double dt, int nsamp, double t0,
	string n, string s, string c, string l)
{
    net=n;
    sta=s;
    chan=c;
    loc=l;
    int npackets;
    npackets=this->write<T>(d,dt,nsamp,t0);
    return(npackets);
}
