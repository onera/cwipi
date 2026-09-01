/*
 * \file
 */

#ifndef __CWP_LOGGING_H__
#define __CWP_LOGGING_H__

/*
  This file is part of the CWIPI library.

  Copyright (C) 2021-2023  ONERA

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
/* Standard C library headers */
#include <stdio.h>
#include <stdbool.h>
/*-----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Public types
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

typedef enum {

  CWP_LOG_TRACE,  /*!< Trace       log */
  CWP_LOG_DEBUG,  /*!< Debug       log */
  CWP_LOG_INFO,   /*!< Info        log */
  CWP_LOG_WARN,   /*!< Warning     log */
  CWP_LOG_ERROR,  /*!< Error       log */
  CWP_LOG_FATAL   /*!< Fatal error log */

} CWP_Log_t;

void _cwp_log(int level, const char* func, const char *file, int line, const char *fmt, ...);
void _cwp_log_assert(bool condition, const char* func, const char *file, int line, const char *fmt, ...);
void _cwp_log_fail(const char* func, const char *file, int line, const char *fmt, ...);

#define CWP_Log(level, ...) _cwp_log(level, __func__, __FILE__, __LINE__, __VA_ARGS__)
#define CWP_Log_assert(condition, ...) _cwp_log_assert(condition, __func__, __FILE__, __LINE__, __VA_ARGS__)
#define CWP_Log_fail(...) _cwp_log_fail(__func__, __FILE__, __LINE__, __VA_ARGS__)

/**
 *
 * \brief Set logging level on local rank
 *
 * \param [in] level Logging level
 *
 */
void CWP_Log_level_set(int level);

/**
 *
 * \brief Enable/disable console logging on local rank
 *
 * \param [in] enable Enable console logging on local rank (0 = quiet [default], 1 = verbose)
 *
 */
void CWP_Log_console_enable(int enable);

/**
 *
 * \brief Enable/disable logfile logging on local rank
 *
 * \param [in] enable Enable logfile logging on local rank (0 = quiet [default], 1 = verbose)
 *
 */
void CWP_Log_logfile_enable(int enable);

/**
 *
 * \brief Pretty print of array in log
 *
 * \param [in]    level  Logging level
 * \param [inout] array  Array to print
 * \param [in]    larray Array length
 * \param [inout] header First line of log
 *
 */
void
CWP_Log_array_int
(
 const int   level,
 const int*  array,
 const int   larray,
 const char* header
);

/**
 *
 * \brief Pretty print of array in log
 *
 * \param [in]    level  Logging level
 * \param [inout] array  Array to print
 * \param [in]    larray Array length
 * \param [inout] header First line of log
 *
 */
void
CWP_Log_array_double
(
 const int     level,
 const double* array,
 const int     larray,
 const char*   header
);

/**
 *
 * \brief Pretty print of array in log
 *
 * \param [in]    level  Logging level
 * \param [inout] array  Array to print
 * \param [in]    larray Array length
 * \param [inout] header First line of log
 *
 */
void
CWP_Log_array_size_t
(
 const int     level,
 const size_t *array,
 const int     larray,
 const char   *header
);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __CWP_LOGGING_H__ */
