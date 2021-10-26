#ifndef __PDM_LOGGING_H__
#define __PDM_LOGGING_H__

/*-----------------------------------------------------------------------------*/

/* Standard C library headers */
#include <stdio.h>
#include <stdarg.h>

/* BFT library headers */

#include "pdm.h"
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

typedef void (*log_lock_fn)(void *udata, int lock);

enum { LOG_TRACE, LOG_DEBUG, LOG_INFO, LOG_WARN, LOG_ERROR, LOG_FATAL };

#define log_trace(...) log_log(LOG_TRACE, __FILE__, __LINE__, __VA_ARGS__)
#define log_debug(...) log_log(LOG_DEBUG, __FILE__, __LINE__, __VA_ARGS__)
#define log_info(...)  log_log(LOG_INFO,  __FILE__, __LINE__, __VA_ARGS__)
#define log_warn(...)  log_log(LOG_WARN,  __FILE__, __LINE__, __VA_ARGS__)
#define log_error(...) log_log(LOG_ERROR, __FILE__, __LINE__, __VA_ARGS__)
#define log_fatal(...) log_log(LOG_FATAL, __FILE__, __LINE__, __VA_ARGS__)

void log_set_udata(void *udata);
void log_set_lock(log_lock_fn fn);
void log_set_fp(FILE *fp);
void log_set_level(int level);
void log_set_quiet(int enable);

void log_log(int level, const char *file, int line, const char *fmt, ...);

/**
 *
 * \brief Pretty print of array in trace_log
 *
 * \param [inout] array        Array to print
 * \param [in]    lArray       Array length
 * \param [inout] header       First line of log
 *
 */
void
PDM_log_trace_array_int
(
 const int*  array,
 const int   larray,
 const char* header
);

/**
 *
 * \brief Pretty print of array in trace_log
 *
 * \param [inout] array        Array to print
 * \param [in]    lArray       Array length
 * \param [inout] header       First line of log
 *
 */
void
PDM_log_trace_array_long
(
 const PDM_g_num_t* array,
 const int          larray,
 const char*        header
);

/**
 *
 * \brief Pretty print of array in trace_log
 *
 * \param [inout] array        Array to print
 * \param [in]    lArray       Array length
 * \param [inout] header       First line of log
 *
 */
void
PDM_log_trace_array_double
(
 const double*      array,
 const int          larray,
 const char*        header
);


/**
 *
 * \brief Pretty print of array in trace_log
 *
 * \param [inout] array        Array to print
 * \param [in]    lArray       Array length
 * \param [inout] header       First line of log
 *
 */
void
PDM_log_trace_array_size_t
(
 const size_t *array,
 const int     larray,
 const char   *header
);

/**
 *
 * \brief Pretty print of array in trace_log
 *
 * \param [inout] array        Array to print
 * \param [in]    lArray       Array length
 * \param [inout] header       First line of log
 *
 */
void
PDM_log_trace_connectivity_long
(
 const int         *array_idx,
 const PDM_g_num_t *array,
 const int          larray,
 const char*        header
);

void
PDM_log_trace_connectivity_int
(
 const int         *array_idx,
 const int         *array,
 const int          larray,
 const char*        header
);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_PRINTF_H__ */
