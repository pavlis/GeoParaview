/* This little test function was written to test revised error handlers in
   the SEISPP library.  All now inherit the what method from std::exception
   and some of the higher order ones now inherit SeisppError.  The goals is 
   to allow all errors to be caught with a SeisppError& and which will print
   the proper log_error result. 
*/
#include <string>
#include "SeisppError.h"
#include "Metadata.h"
#include "VelocityModel_1d.h"
#include "seispp.h"
using namespace std;
using namespace SEISPP;
bool SEISPP::SEISPP_verbose(true);
double foo[250000];
int main(int argc, char **argv)
{
    try {
        throw SeisppError(string("I am a SeisppError object"));
    } catch(SeisppError& serr)
    {
        cout << "Testing SeisppError object catch"<<endl;
        cout << "Calling log_error method"<<endl;
        serr.log_error();
        cout << "Same with what() method"<<endl;
        cout << serr.what()<<endl;
    }
    cout <<"////////////////////////////////////////////////////"<<endl;
    try {
        Dbptr db={1,2,3,4};
        throw SeisppDberror(string("I am a SeisppDberror"),db);
    } catch(SeisppError& serr)
    {
        cout << "Testing SeisppDberror object catch"<<endl;
        cout << "Calling log_error method"<<endl;
        serr.log_error();
        cout << "Same with what() method"<<endl;
        cout << serr.what()<<endl;
    }
    cout <<"////////////////////////////////////////////////////"<<endl;
    try {
        throw MetadataError(string("I am a MetadataError"));
    } catch(SeisppError& serr)
    {
        cout << "Testing Metadata object catch"<<endl;
        cout << "Calling log_error method"<<endl;
        serr.log_error();
        cout << "Same with what() method"<<endl;
        cout << serr.what()<<endl;
    }
    cout <<"////////////////////////////////////////////////////"<<endl;
    try {
        string mdtype("test"),name("foo");
        throw MetadataGetError(mdtype,name,string("I am a MetadataGetError"));
    } catch(SeisppError& serr)
    {
        cout << "Testing MetadataGetError object catch"<<endl;
        cout << "Calling log_error method"<<endl;
        serr.log_error();
        cout << "Same with what() method"<<endl;
        cout << serr.what()<<endl;
    }
    cout <<"////////////////////////////////////////////////////"<<endl;
    try {
        int ecode(-1);
        throw MetadataParseError(ecode,string("I am a MetadataParseError"));
    } catch(SeisppError& serr)
    {
        cout << "Testing MetadataParseError object catch"<<endl;
        cout << "Calling log_error method"<<endl;
        serr.log_error();
        cout << "Same with what() method"<<endl;
        cout << serr.what()<<endl;
    }
    cout <<"////////////////////////////////////////////////////"<<endl;
    try {
        throw VelocityModel_1d_Error(string("I am a VelocityModel_1d_Error error"));
    } catch(SeisppError& serr)
    {
        cout << "Testing VelocityModel_1d_Error object catch"<<endl;
        cout << "Calling log_error method"<<endl;
        serr.log_error();
        cout << "Same with what() method"<<endl;
        cout << serr.what()<<endl;
    }
    cout <<"////////////////////////////////////////////////////"<<endl;
    try {
        string modname("testmod");
        throw VelocityModel_1d_Dberror(modname,string("I am a MetadataParseError"));
    } catch(SeisppError& serr)
    {
        cout << "Testing VelocityModel_1d_Dberror object catch"<<endl;
        cout << "Calling log_error method"<<endl;
        serr.log_error();
        cout << "Same with what() method"<<endl;
        cout << serr.what()<<endl;
    }
    cout <<"////////////////////////////////////////////////////"<<endl;
    try {
        string ioerr("test i/o detail message");
        throw VelocityModel_1d_IOerror(string("I am a MetadataParseError"),ioerr);
    } catch(SeisppError& serr)
    {
        cout << "Testing VelocityModel_1d_IOerror object catch"<<endl;
        cout << "Calling log_error method"<<endl;
        serr.log_error();
        cout << "Same with what() method"<<endl;
        cout << serr.what()<<endl;
    }
}

