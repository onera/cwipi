#ifndef __INITBFTPRINTF__
#define __INITBFTPRINTF__

#include <stdio.h>
#include <stdarg.h>
extern FILE* listing;

int initBftPrintf (const char     *const format, va_list arg_ptr);

#endif /* __BAR_COORDS_H__ */
