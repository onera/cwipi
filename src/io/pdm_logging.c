/*============================================================================
 * Base user-definable PDM_printf() wrapper or replacement.
 *============================================================================*/

/*
  This file is part of the CWIPI library.

  Copyright (C) 2017 ONERA

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

/*-----------------------------------------------------------------------------*/

#include "pdm_config.h"

/*
 * Standard C library headers
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <time.h>


/*
 * Optional library and BFT headers
 */

#include "pdm_logging.h"
#include "pdm_mpi.h"

/*-----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*-----------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------
 * Local type definitions
 *-----------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------
 * Local function prototypes
 *-----------------------------------------------------------------------------*/

static FILE* logging_file = NULL;

static
void
free_logging_file
(
)
{
  if(logging_file != NULL)
    fclose(logging_file);
}

static struct {
  void *udata;
  log_lock_fn lock;
  FILE *fp;
  int level;
  int quiet;
} L;


// static const char *level_names[] = {
//   "TRACE", "DEBUG", "INFO", "WARN", "ERROR", "FATAL"
// };

#ifdef LOG_USE_COLOR
static const char *level_colors[] = {
  "\x1b[94m", "\x1b[36m", "\x1b[32m", "\x1b[33m", "\x1b[31m", "\x1b[35m"
};
#endif


static void lock(void)   {
  if (L.lock) {
    L.lock(L.udata, 1);
  }
}


static void unlock(void) {
  if (L.lock) {
    L.lock(L.udata, 0);
  }
}


void log_set_udata(void *udata) {
  L.udata = udata;
}


void log_set_lock(log_lock_fn fn) {
  L.lock = fn;
}


void log_set_fp(FILE *fp) {
  L.fp = fp;
}


void log_set_level(int level) {
  L.level = level;
}


void log_set_quiet(int enable) {
  L.quiet = enable ? 1 : 0;
}


void log_log(int level, const char *file, int line, const char *fmt, ...) {
  if (level < L.level) {
    return;
  }

  if(logging_file == NULL){
    char filename[50];
    int i_rank;
    PDM_MPI_Comm_rank (PDM_MPI_COMM_WORLD, &i_rank);
    sprintf(filename, "paradigm_%d.log", i_rank);
    logging_file = fopen(filename, "w");
    atexit(free_logging_file);
  }

  /* Acquire lock */
  lock();

  /* Log to file */
  if (L.fp) {
    va_list args;
    va_start(args, fmt);
    // fprintf(L.fp, "-- %s -- | ", level_names[level]);
    vfprintf(L.fp, fmt, args);
    va_end(args);
    fflush(L.fp);
  }

  va_list args;
  va_start(args, fmt);
  // fprintf(logging_file, "-- %s -- | ", level_names[level]);
  vfprintf(logging_file, fmt, args);
  va_end(args);
  fflush(logging_file);

  /* Release lock */
  unlock();
}

// void log_log(int level, const char *file, int line, const char *fmt, ...) {
//   if (level < L.level) {
//     return;
//   }

//   if(logging_file == NULL){
//     char filename[50];
//     int i_rank;
//     PDM_MPI_Comm_rank (PDM_MPI_COMM_WORLD, &i_rank);
//     sprintf(filename, "paradigm_%d.log", i_rank);
//     logging_file = fopen(filename, "w");
//     atexit(free_logging_file);
//   }

//   /* Acquire lock */
//   lock();

//   /* Get current time */
//   time_t t = time(NULL);
//   struct tm *lt = localtime(&t);

//   /* Log to stderr */
//   if (!L.quiet) {
//     va_list args;
//     char buf[16];
//     buf[strftime(buf, sizeof(buf), "%H:%M:%S", lt)] = '\0';
// #ifdef LOG_USE_COLOR
//     fprintf(
//       stderr, "%s %s%-5s\x1b[0m \x1b[90m%s:%d:\x1b[0m ",
//       buf, level_colors[level], level_names[level], file, line);
// #else
//     fprintf(stderr, "%s %-5s %s:%d: ", buf, level_names[level], file, line);
// #endif
//     va_start(args, fmt);
//     vfprintf(stderr, fmt, args);
//     va_end(args);
//     fprintf(stderr, "\n");
//     fflush(stderr);
//   }

//   /* Log to file */
//   if (L.fp) {
//     va_list args;
//     // char buf[32];
//     // buf[strftime(buf, sizeof(buf), "%Y-%m-%d %H:%M:%S", lt)] = '\0';
//     // fprintf(L.fp, "%s %-5s %s:%d: ", buf, level_names[level], file, line);
//     va_start(args, fmt);
//     vfprintf(L.fp, fmt, args);
//     va_end(args);
//     // fprintf(L.fp, "\n");
//     fflush(L.fp);
//   }

//   /* Release lock */
//   unlock();
// }



void
PDM_log_trace_array_int
(
 int* array,
 int  larray,
 const char* header
)
{
  log_trace(header);
  for(int i = 0; i < larray; ++i){
    log_trace("%d ", array[i]);
  }
  log_trace("\n");
}

void
PDM_log_trace_array_long
(
 PDM_g_num_t* array,
 int          larray,
 const char*  header
)
{
  log_trace(header);
  for(int i = 0; i < larray; ++i){
    log_trace("%d ", array[i]);
  }
  log_trace("\n");
}


#ifdef __cplusplus
}
#endif /* __cplusplus */
