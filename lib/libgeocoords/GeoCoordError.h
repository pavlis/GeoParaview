#ifndef _GEOCOORDERROR_H_
#define _GEOCOORDERROR_H_
#include <stdlib.h>
#include <string>
using namespace std;
/* \brief Error class thrown for errors in this package.

Simple child of stdlib abstract exception class.  Handler should catch a
generic exception and call the what() method to print the message.
*/
class GeoCoordError : public exception
{
public:
    GeoCoordError(string mess){message="GeooCoordError:  "+mess;};
    GeoCoordError(const char *s){message="GeooCoordError:  "+string(s);};
    virtual const char* what() const throw()
   {
    return message.c_str();
   };
   /* This is necessary baggage to mesh with std::exception declaration
      for this virtual method.  Not required on all compilers */
   ~GeoCoordError() throw ()
   {};
protected:
    string message;
};
#endif
