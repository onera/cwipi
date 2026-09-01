/*============================================================================
 * Logging adapted from "pdm_logging"
 *============================================================================*/

/*
  This file is part of the CWIPI library.

  Copyright (C) 2022-2023  ONERA

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library. If not, see <http://www.gnu.org/licenses/>.
*/

/*-----------------------------------------------------------------------------*/

/*
 * Standard C library headers
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <stdbool.h>
#include <string.h>

#include "cwp_logging.h"
#include "mpi.h"
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
free_logging_file(void)
{
  if (logging_file != NULL)
    fclose(logging_file);
}

static const char* _filename_without_path(const char* path) {
  const char* slash     = strrchr(path, '/');
  const char* backslash = strrchr(path, '\\');
  if (slash && backslash) {
    return (slash > backslash) ? slash + 1 : backslash + 1;
  } else if (slash) {
    return slash + 1;
  } else if (backslash) {
    return backslash + 1;
  }
  return path;
}

static struct {
  int level;
  int verbose_console;
  int verbose_logfile;
} L = {CWP_LOG_FATAL, 0, 0};

static void _cwp_log_va(int level, bool add_line_break, const char* func, const char *file, int line, const char *fmt, va_list args){

  // Supress unused warning
  (void)(level);
  (void)(func);
  (void)(file);
  (void)(line);

  va_list args_copy;
  va_copy(args_copy, args);

  if (L.verbose_logfile == 1) {
    if (logging_file == NULL) {
      char filename[50];
      int i_rank;
      MPI_Comm_rank (MPI_COMM_WORLD, &i_rank);
      sprintf(filename, "cwp_%d.log", i_rank);
      logging_file = fopen(filename, "w");
      atexit(free_logging_file);
    }

    vfprintf(logging_file, fmt, args);
    if (add_line_break) fprintf(logging_file, "\n");
    fflush(logging_file);
  }

  if (L.verbose_console == 1) {
    vfprintf(stdout, fmt, args_copy);
    if (add_line_break) fprintf(stdout, "\n");
    fflush(stdout);
  }

  va_end(args_copy);
}

/*-----------------------------------------------------------------------------
 * Public function prototypes
 *-----------------------------------------------------------------------------*/

void CWP_Log_level_set(int level) {
  L.level = level;
}

void CWP_Log_console_enable(int enable) {
  L.verbose_console = enable ? 1 : 0;
}

void CWP_Log_logfile_enable(int enable) {
  L.verbose_logfile = enable ? 1 : 0;
}

void _cwp_log(int level, const char* func, const char *file, int line, const char *fmt, ...) {
  if ((L.verbose_console == 0 && L.verbose_logfile == 0) || level < L.level) {
    return;
  }

  bool add_line_break = false;
  size_t len = strlen(fmt);
  if (len == 0){
    add_line_break = true;
  }
  else if (strcmp(&fmt[len-1], "\n") != 0) {
    add_line_break = true;
  }

  va_list args;
  va_start(args, fmt);
  _cwp_log_va(level, add_line_break, func, file, line, fmt, args);
  va_end(args);
}

void _cwp_log_assert(bool condition, const char* func, const char *file, int line, const char *fmt, ...) {
  if (condition) return;
  _cwp_log(CWP_LOG_FATAL, func, file, line, "Assertion error in \"%s\" (file %s, line %d)", func, _filename_without_path(file), line);
  va_list args;
  va_start(args, fmt);
  _cwp_log_va(CWP_LOG_FATAL, true, func, file, line, fmt, args);
  va_end(args);
  fflush(stdout);
  fflush(stderr);
  exit(1);
}

void _cwp_log_fail(const char* func, const char *file, int line, const char *fmt, ...) {
  _cwp_log(CWP_LOG_FATAL, func, file, line, "Fatal error in \"%s\" (file %s, line %d)", func, _filename_without_path(file), line);
  va_list args;
  va_start(args, fmt);
  _cwp_log_va(CWP_LOG_FATAL, true, func, file, line, fmt, args);
  va_end(args);
  fflush(stdout);
  fflush(stderr);
  exit(1);
}

void
CWP_Log_array_int
(
 const int   level,
 const int*  array,
 const int   larray,
 const char* header
)
{
  CWP_Log(level, header);
  CWP_Log(level, "(%i) -> ", larray);
  for(int i = 0; i < larray; ++i){
    CWP_Log(level, "%d ", array[i]);
  }
  CWP_Log(level, "\n");
}

void
CWP_Log_array_double
(
 const int     level,
 const double* array,
 const int     larray,
 const char*   header
)
{
  CWP_Log(level, header);
  CWP_Log(level, "(%i) -> ", larray);
  for(int i = 0; i < larray; ++i){
    CWP_Log(level, "%20.16e ", array[i]);
  }
  CWP_Log(level, "\n");
}

void
CWP_Log_array_size_t
(
 const int     level,
 const size_t *array,
 const int     larray,
 const char   *header
)
{
  CWP_Log(level, header);
  for(int i = 0; i < larray; ++i){
    CWP_Log(level, "%lu ", array[i]);
  }
  CWP_Log(level, "\n");
}

#ifdef __cplusplus
}
#endif /* __cplusplus */
