#include "initBftPrintf.h"
FILE* listing;
int initBftPrintf
(
 const char     *const format,
       va_list         arg_ptr
)
{
  //  assert(listing != NULL);
  return vfprintf(listing, format, arg_ptr);
}
