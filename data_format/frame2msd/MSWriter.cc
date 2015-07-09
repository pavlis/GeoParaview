#include "MSWriter.h"
MSWriter::MSWriter(string fname, int reclen)
    : net(),sta(),chan(),loc()
{
    const string base_error("MWWriter constructor:  ");
    record_length=static_cast<int32_t>(reclen);
    MSWfp=fopen(fname.c_str(),"a");
    if(MSWfp==NULL) throw MSWexception(base_error
    			+"fopen failed on file name="
    			+fname);
}
MSWriter::~MSWriter()
{
    fclose(MSWfp);
}
long MSWriter::current_filesize()
{
    if(fseek(MSWfp,0L,SEEK_END))
            throw MSWexception("MSWriter::current_filesize:  fseek failure on output file");

    return(ftell(MSWfp));
}
