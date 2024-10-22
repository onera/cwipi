/*============================================================================
 * Base memory usage information (System and Library dependent)
 *============================================================================*/

/*
  This file is part of the "Base Functions and Types" library, intended to
  simplify and enhance portability, memory and I/O use for scientific codes.

  Copyright (C) 2004-2009  EDF

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

#include "config_priv.h"
#include "bftc_config_defs.h"

/* On Solaris, procfs may not be compiled in a largefile environment,
 * so we redefine macros before including any system header file. */

#if defined(bftc_OS_Solaris) && defined(HAVE_UNISTD_H) \
   && defined(HAVE_SYS_PROCFS_H) && !defined(__cplusplus)
#define _STRUCTURED_PROC 1
#undef _FILE_OFFSET_BITS
#define _FILE_OFFSET_BITS 32
#endif

/*
 * Standard C library headers
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined (bftc_OS_Linux) && defined(HAVE_SYS_STAT_H) \
                           && defined(HAVE_SYS_TYPES_H)
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#elif defined(bftc_OS_OSF1) && defined(_OSF_SOURCE) && defined(HAVE_UNISTD_H)
#include <fcntl.h>
#include <sys/types.h>
#include <sys/signal.h>
#include <sys/fault.h>
#include <sys/syscall.h>
#include <sys/procfs.h>
#include <unistd.h>

#elif defined(bftc_OS_Solaris) && defined(HAVE_UNISTD_H) \
   && defined(HAVE_SYS_PROCFS_H) && !defined(__cplusplus)
#include <sys/types.h>
#include <sys/procfs.h>
#include <unistd.h>

#elif (defined(bftc_OS_IRIX64) || defined(bftc_OS_UXPV))
#if defined(HAVE_UNISTD_H) && defined(HAVE_SYS_TYPES_H) \
                           && defined(HAVE_SYS_STAT_H)
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#endif

#elif defined (bftc_OS_AIX) && defined(HAVE_GETRUSAGE)
#include <sys/times.h>
#include <sys/resource.h>

#elif defined(HAVE_GETRUSAGE)
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#endif

#if defined(HAVE_UNISTD_H) && defined(HAVE_SBRK) 
#if defined(__blrts__) || defined(__bgp_)
#define USE_SBRK 1
#elif defined (bftc_OS_Linux)
#define __USE_MISC 1
#endif
#include <unistd.h>
#endif

#if defined(HAVE_MALLOC_HOOKS)
#include <malloc.h>
#endif

#if defined(HAVE_STDDEF_H)
#include <stddef.h>
#endif

/*
 * Optional library and BFT headers
 */

#include "bftc_mem_usage.h"

/*-----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*-------------------------------------------------------------------------------
 * Local type definitions
 *-----------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------
 * Local static variable definitions
 *-----------------------------------------------------------------------------*/

static int  _bftc_mem_usage_global_initialized = 0;

static size_t _bftc_mem_usage_global_max_pr = 0;

#if defined(USE_SBRK)
static void  *_bftc_mem_usage_global_init_sbrk = NULL;
#endif

#if defined (bftc_OS_Linux) && defined(HAVE_SYS_STAT_H) \
                           && defined(HAVE_SYS_TYPES_H)
static int    _bftc_mem_usage_proc_file_init = 0;
static int    _bftc_mem_usage_proc_file_peak = 0;
static int    _bftc_mem_usage_proc_file_fd   = -1;
static long   _bftc_mem_usage_proc_file_pos  = 0;
#endif

#if defined(HAVE_MALLOC_HOOKS)
static __malloc_ptr_t
(* _bftc_mem_usage_old_malloc_hook)   (size_t,
                                      const __malloc_ptr_t);
static __malloc_ptr_t
(* _bftc_mem_usage_old_realloc_hook)  (void *,
                                      size_t,
                                      const __malloc_ptr_t);
static void
(* _bftc_mem_usage_old_free_hook)     (void *,
                                      const __malloc_ptr_t);

static int _bftc_mem_usage_global_use_hooks = 0;
static int _bftc_mem_usage_global_up = 0;

#endif /* (HAVE_MALLOC_HOOKS) */

/*-----------------------------------------------------------------------------
 * Local function definitions
 *-----------------------------------------------------------------------------*/

#if defined(HAVE_MALLOC_HOOKS)

/*
 * Memory size update to be called by malloc_hook functions.
 */

static void
_bftc_mem_usage_size_update(void)
{
  (void) bftc_mem_usage_pr_size();
}

/*
 * Test malloc_hook function.
 *
 * This function does not allocate memory. When it is called, it sets
 *  the _bftc_mem_usage_global_use_hooks global counter to 1, indicating
 * the malloc hooks are effective and may be called. This should be the
 * usual case when linking with the glibc, except when we prelink with
 * some specific allocation library, such as is the case when using
 * Electric Fence.
 *
 * returns:
 *   1 (NULL preferred, 1 should avoid TotalView warning).
 */

static __malloc_ptr_t
_bftc_mem_usage_malloc_hook_test(size_t size,
                                const __malloc_ptr_t _ptr)
{
  _bftc_mem_usage_global_use_hooks = 1;
  return (__malloc_ptr_t)1;
}

/*
 * Memory counting malloc_hook function.
 *
 * This function calls the regular realloc function, but tries to keep
 * track of whether allocated memory is increasing or decreasing;
 * just before allocated memory decreases, memory counting functions
 * are called unless deactivated through bftc_mem_usage_set_options().
 *
 * returns:
 *   Pointer to allocated memory.
 */

static __malloc_ptr_t
_bftc_mem_usage_malloc_hook(size_t size,
                           const __malloc_ptr_t _ptr)
{
  void *result;

  __malloc_hook = _bftc_mem_usage_old_malloc_hook;

  if (_bftc_mem_usage_global_up == 0)
    _bftc_mem_usage_global_up = 1;

  result = malloc(size);

  __malloc_hook = _bftc_mem_usage_malloc_hook;

  return result;
}

/*
 * Memory counting realloc_hook function.
 *
 * This function calls the regular realloc function, but tries to keep
 * track of whether allocated memory is increasing or decreasing;
 * just before allocated memory decreases, memory counting functions
 * are called unless deactivated through bftc_mem_usage_set_options().
 */

static __malloc_ptr_t
_bftc_mem_usage_realloc_hook(void *ptr,
                            size_t size,
                            const __malloc_ptr_t _ptr)
{
  void *result;

  /* Protect __malloc_hook as well as __realloc hook, in case
     realloc() uses malloc(). If we do not reset __malloc_hook
     before exiting here, the __malloc_hook may be unset after
     a realloc() on some systems */

  __realloc_hook = _bftc_mem_usage_old_realloc_hook;
  __malloc_hook  = _bftc_mem_usage_old_malloc_hook;

  /* Memory usage may be going down, so update */
  if (_bftc_mem_usage_global_up == 1)
    _bftc_mem_usage_size_update();

  result = realloc(ptr, size);

  /* Memory usage may be going up; set flag just in case */
  if (_bftc_mem_usage_global_up == 0)
    _bftc_mem_usage_global_up = 1;

  /* Reset hooks */
  __realloc_hook = _bftc_mem_usage_realloc_hook;
  __malloc_hook  = _bftc_mem_usage_malloc_hook;

  return result;
}

/*
 * Memory counting free_hook function.
 *
 * This function calls the regular free function, but tries to keep
 * track of whether allocated memory is increasing or decreasing;
 * just before allocated memory decreases, memory counting functions
 * are called unless deactivated through bftc_mem_usage_set_options().
 */

static void
_bftc_mem_usage_free_hook(void *ptr,
                         const __malloc_ptr_t _ptr)
{
  __free_hook = _bftc_mem_usage_old_free_hook;

  if (_bftc_mem_usage_global_up == 1) {
    _bftc_mem_usage_global_up = 0;
    _bftc_mem_usage_size_update();
  }

  free(ptr);

  __free_hook = _bftc_mem_usage_free_hook;
}

/*
 * Set this library's memory counting malloc hooks if possible.
 */

static void
_bftc_mem_usage_set_hooks(void)
{
  /* Test if hooks may really be used (i.e. if there
     is no prelinking with some other allocation library) */

  if (_bftc_mem_usage_global_use_hooks == 0) {

    static __malloc_ptr_t
      (* old_malloc_hook) (size_t, const __malloc_ptr_t);
    void *ptr_test;

    old_malloc_hook = __malloc_hook;
    __malloc_hook = _bftc_mem_usage_malloc_hook_test;
    ptr_test = malloc(128);
    /* We have really allocated ptr_test if hooks are not active */
    if (ptr_test != NULL && _bftc_mem_usage_global_use_hooks == 0)
      free(ptr_test);
    __malloc_hook = old_malloc_hook;

  }

  /* Set memory counting hooks */

  if (_bftc_mem_usage_global_use_hooks != 0) {

    if (__malloc_hook != _bftc_mem_usage_malloc_hook) {
      _bftc_mem_usage_old_malloc_hook = __malloc_hook;
      __malloc_hook = _bftc_mem_usage_malloc_hook;
    }
    if (__realloc_hook != _bftc_mem_usage_realloc_hook) {
      _bftc_mem_usage_old_realloc_hook = __realloc_hook;
      __realloc_hook = _bftc_mem_usage_realloc_hook;
    }
    if (__free_hook != _bftc_mem_usage_free_hook) {
      _bftc_mem_usage_old_free_hook = __free_hook;
      __free_hook = _bftc_mem_usage_free_hook;
    }

  }
}

/*
 * Unset this library's memory counting malloc hooks if possible.
 */

static void
_bftc_mem_usage_unset_hooks(void)
{
  if (_bftc_mem_usage_global_use_hooks != 0) {
    __malloc_hook  = _bftc_mem_usage_old_malloc_hook;
    __realloc_hook = _bftc_mem_usage_old_realloc_hook;
    __free_hook    = _bftc_mem_usage_old_free_hook;
    _bftc_mem_usage_global_use_hooks = 0;
  }
}

#endif /* defined(HAVE_MALLOC_HOOKS) */

#if defined (bftc_OS_Linux) && defined(HAVE_SYS_STAT_H) \
                           && defined(HAVE_SYS_TYPES_H)

/*!
 * \brief Initialize current process memory use count depending on system.
 */

static void
_bftc_mem_usage_pr_size_init(void)
{
  char  buf[512]; /* should be large enough for "/proc/%lu/status"
                     then beginning of file content */
  size_t  r_size, i;
  _Bool   status_has_peak = false;
  const pid_t  pid = getpid();

  /*
    Under Linux with procfs, one line of the pseudo-file "/proc/pid/status"
    (where pid is the process number) is of the following form:
    VmSize:     xxxx kB
    This line may be the 10th for a 2.4 kernel, or the 12th to 13th for
    a 2.6.x kernel. On more recent 2.6.x kernels, another line (the 12th)
    is of the form:
    VmPeak:     xxxx kB
    When VmSize is available but not VmPeak, we use low-level file I/O on
    this file (using read() instead of fread()) to avoid update issues due
    to buffering when keeping the file open and re-reading the area
    containing the "VmSize" information.
    When VmPeak is also available, memory use tracking is not necessary, so
    the number of calls to bftc_mem_usage_pr_size() is presumably much lower
    (being based solely on the user code, not on bftc_mem_usage's tracking),
    so we prefer to open and close "/proc/pid/status" on each
    call to bftc_mem_usage_pr_size().
  */

  if (_bftc_mem_usage_proc_file_init != 0)
    return;

  sprintf(buf, "/proc/%lu/status", (unsigned long) pid);

  _bftc_mem_usage_proc_file_fd = open(buf, O_RDONLY);

  if (_bftc_mem_usage_proc_file_fd != -1) {

    r_size = read(_bftc_mem_usage_proc_file_fd, buf, 512);

    if (r_size > 32) { /* Leave a margin for "VmPeak" or "VmSize:" line */
      r_size -= 32;
      for (i = 0; i < r_size; i++) {
        if (buf[i] == 'V' && strncmp(buf+i, "VmPeak:", 7) == 0) {
          status_has_peak = true;
          break;
        }
      }
      for (i = 0; i < r_size; i++) {
        if (buf[i] == 'V' && strncmp(buf+i, "VmSize:", 7) == 0)
          break;
      }
      /* If VmSize was found, proc file may be used */
      if (i < r_size) {
        if (status_has_peak == true) {
          _bftc_mem_usage_proc_file_peak = 1;
        }
        else {
          /* set search position, leaving small margin ahead of position
             in case file content length before "VmSize:" changes
             (typically by 1 or 2 chars) */
          i = (i > 16) ? i -16 : 0;
          _bftc_mem_usage_proc_file_pos = i;
        }
        _bftc_mem_usage_proc_file_init = 1;
      }
    }

    /* Close file if unusable or if VmPeak is available */
    if (   _bftc_mem_usage_proc_file_init == -1
        || _bftc_mem_usage_proc_file_peak == 1) {
      (void)close(_bftc_mem_usage_proc_file_fd);
      _bftc_mem_usage_proc_file_fd = -1;
    }
  }

  /* If initialization failed for some reason (proc file unavailable or does
     or does not contain the required fields), mark method as unusable */
  if (_bftc_mem_usage_proc_file_init == 0)
    _bftc_mem_usage_proc_file_init = -1;

#if defined(HAVE_MALLOC_HOOKS)
  /* Unset tracking hooks if unneeded */
  _bftc_mem_usage_unset_hooks();
#endif
}

/*!
 * \brief Finalize current process memory use count depending on system.
 */

static void
_bftc_mem_usage_pr_size_end(void)
{
  if (_bftc_mem_usage_proc_file_init != 1)
    return;

  if (_bftc_mem_usage_proc_file_fd != -1) {
    (void)close(_bftc_mem_usage_proc_file_fd);
    _bftc_mem_usage_proc_file_init = 0;
    _bftc_mem_usage_proc_file_fd = -1;
  }
}

#else  /* defined (bftc_OS_Linux) && ... */

#define _bftc_mem_usage_pr_size_init()
#define _bftc_mem_usage_pr_size_end()

#endif /* defined (bftc_OS_Linux) && ... */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*!
 * \brief Initialize memory usage count depending on system.
 *
 * This functions checks if it has already been called, so
 * it is safe to call more than once (though it is not
 * thread-safe). Only the first call is effective.
 */

void
bftc_mem_usage_init(void)
{
  if (_bftc_mem_usage_global_initialized != 0)
    return;

#if defined(USE_SBRK)

  /*
    We use sbrk() to know the size of the heap. This is not of any use
    to guess at allocated memory when some part of the memory may
    be allocated with mmap(), such as with glibc on Linux.
  */

  _bftc_mem_usage_global_init_sbrk = (void *) sbrk(0);

#endif /* (USE_SBRK) */

#if defined(HAVE_MALLOC_HOOKS)

  _bftc_mem_usage_set_hooks();

#endif

  _bftc_mem_usage_global_initialized = 1;
}

/*!
 * \brief End memory usage count depending on system.
 */

void
bftc_mem_usage_end(void)
{
#if defined(HAVE_MALLOC_HOOKS)

  _bftc_mem_usage_unset_hooks();

#endif

  _bftc_mem_usage_pr_size_end();
}

/*!
 * \brief Indicates if bftc_mem_usage_...() functions are initialized.
 *
 * \returns 1 if bftc_mem_usage_init has been called, 0 otherwise.
 */

int
bftc_mem_usage_initialized(void)
{
  return _bftc_mem_usage_global_initialized;
}

/*!
 * \brief Return current process memory use (in kB) depending on system.
 *
 * If the information is not available (depending on availability of
 * non-portable function calls), 0 is returned.
 */

#if defined (bftc_OS_Linux) && defined(HAVE_SYS_STAT_H) \
                           && defined(HAVE_SYS_TYPES_H)

size_t
bftc_mem_usage_pr_size(void)
{
  size_t sys_mem_usage = 0;

  /*
    Under Linux with procfs, one line of the pseudo-file "/proc/pid/status"
    (where pid is the process number) is of the following form:
    VmSize:     xxxx kB
    With more recent kernels, we also have a line of the form:
    VmPeak:     xxxx kB
  */

  {
    if (_bftc_mem_usage_proc_file_init == 0)
      _bftc_mem_usage_pr_size_init();

    if (_bftc_mem_usage_proc_file_init == 1) {

      /* If "VmPeak:" unavailable, find "VmSize:" in already opened file */

      if (_bftc_mem_usage_proc_file_peak == 0) {

        char  buf[64]; /* should be large enough for "VmSize:" area */
        size_t  r_size, i;
        unsigned long  val;

        if (lseek(_bftc_mem_usage_proc_file_fd,
                  _bftc_mem_usage_proc_file_pos,
                  SEEK_SET) > -1) {

          r_size = read(_bftc_mem_usage_proc_file_fd, buf, 64);

          if (r_size > 32) { /* Leave a margin for "VmSize:" line */
            buf[r_size - 1] = '\0';
            r_size -= 32;
            for (i = 0; i < r_size; i++)
              if (buf[i] == 'V' && strncmp(buf+i, "VmSize:", 7) == 0)
                break;
            sscanf (buf + i + 7, "%lu", &val);
            sys_mem_usage = (size_t) val;
          }

        }

      }

      /* If "VmPeak:" available, open and close file */

      else {

        char  buf[81]; /* should be large enough for "/proc/%lu/status" */
        const pid_t  pid = getpid();

        FILE *fp;
        unsigned long val;
        char *s;

        sprintf(buf, "/proc/%lu/status", (unsigned long) pid);
        fp = fopen(buf, "r");

        if (fp != NULL) {

          int fields_read = 0;

          while (fields_read < 2) {
            s = fgets(buf, 80, fp);
            if (s == NULL)
              break;
            if (strncmp(s, "VmSize:", 7) == 0) {
              sscanf (s + 7, "%lu", &val);
              sys_mem_usage = (size_t) val;
              fields_read += 1;
            }
            else if (strncmp(s, "VmPeak:", 7) == 0) {
              sscanf (s + 7, "%lu", &val);
              if ((size_t) val > _bftc_mem_usage_global_max_pr)
                _bftc_mem_usage_global_max_pr = (size_t) val;
              fields_read += 1;
            }
          }

          fclose(fp);

        }

      } /* End of condition on "VmPeak:" availability */

    }

#if !defined(HAVE_MALLOC_HOOKS)
    _bftc_mem_usage_pr_size_end();
#endif

  }

  if (sys_mem_usage > _bftc_mem_usage_global_max_pr)
    _bftc_mem_usage_global_max_pr = sys_mem_usage;

  return sys_mem_usage;
}

#elif defined (bftc_OS_OSF1) && defined(_OSF_SOURCE) && defined(HAVE_UNISTD_H)

size_t
bftc_mem_usage_pr_size(void)
{
  size_t sys_mem_usage = 0;

  /* On Compaq Tru64 Unix */
  {
    char        buf[81];  /* should be large enough for "/proc/%lu/status" */
    int         procfile;
    prpsinfo_t  p;

    const  pid_t  pid = getpid();

    sprintf (buf, "/proc/%05lu", (unsigned long) pid);

    procfile = open(buf, O_RDONLY);

    if (procfile != -1) {

      if (ioctl(procfile, PIOCPSINFO, &p) != -1)
        sys_mem_usage  = (p.pr_size * getpagesize()) / 1024;

      close(procfile);

    }

  }

  if (sys_mem_usage > _bftc_mem_usage_global_max_pr)
    _bftc_mem_usage_global_max_pr = sys_mem_usage;

  return sys_mem_usage;
}

#elif defined(bftc_OS_Solaris) && defined(HAVE_UNISTD_H) \
   && defined(HAVE_SYS_PROCFS_H) && !defined(__cplusplus)

size_t
bftc_mem_usage_pr_size(void)
{
  size_t sys_mem_usage = 0;

  {
    /* We have binary pseudo-files /proc/pid/status and /proc/pid/psinfo */

    char   buf[81];     /* should be large enough for "/proc/%lu/status" */
    const  unsigned long  pid = getpid ();

    FILE     *fp;
    int       val;
    char     *s ;
    size_t    ret;
    psinfo_t  pid_info;

    sprintf (buf, "/proc/%lu/psinfo", pid);

    fp = fopen (buf, "r");
    if (fp != NULL) {
      ret = fread(&pid_info, sizeof(pid_info), 1, fp);
      if (ret == 1)
        sys_mem_usage = pid_info.pr_size;
      fclose (fp);
    }

  }

  if (sys_mem_usage > _bftc_mem_usage_global_max_pr)
    _bftc_mem_usage_global_max_pr = sys_mem_usage;

  return sys_mem_usage;
}

#elif (defined(bftc_OS_IRIX64) || defined(bftc_OS_UXPV))

size_t
bftc_mem_usage_pr_size(void)
{
  size_t sys_mem_usage = 0;

#if defined(HAVE_UNISTD_H) && defined(HAVE_SYS_TYPES_H) \
 && defined(HAVE_SYS_STAT_H)
  /* On SGI IRIX and Fujitsu VPP 5000, what follows should work */

  {
    char   buf[81];     /* should be large enough for "/proc/%lu/status" */
    const  pid_t  pid = getpid ();

    struct stat file_stat;

    sprintf (buf, "/proc/%05lu", (unsigned long) pid);

    if (stat (buf, &file_stat) != -1)
      sys_mem_usage = file_stat.st_size / 1024;

  }
#endif /* HAVE_UNISTD_H and SYS_TYPES_H and SYS_STAT_H */

  if (sys_mem_usage > _bftc_mem_usage_global_max_pr)
    _bftc_mem_usage_global_max_pr = sys_mem_usage;

  return sys_mem_usage;
}

#elif defined(USE_SBRK)

size_t
bftc_mem_usage_pr_size(void)
{
  size_t alloc_size = 0;

  if (_bftc_mem_usage_global_initialized) {
    void    *end_addr;

    end_addr = (void *) sbrk(0);

#if defined(HAVE_PTRDIFF_T)
    alloc_size = (size_t)(  (ptrdiff_t)end_addr
                          - (ptrdiff_t)_bftc_mem_usage_global_init_sbrk) / 1024;
#else
    alloc_size = (end_addr - _bftc_mem_usage_global_init_sbrk) / 1024;
#endif

  }

  if (alloc_size > _bftc_mem_usage_global_max_pr)
    _bftc_mem_usage_global_max_pr = alloc_size;

  return alloc_size;
}

#elif defined(HAVE_GETRUSAGE)

size_t
bftc_mem_usage_pr_size(void)
{
  size_t sys_mem_usage = 0;
  struct rusage usage;

  getrusage(RUSAGE_SELF, &usage);

  sys_mem_usage = usage.ru_maxrss / 1024;

  return sys_mem_usage;
}

#else /* Default case */

size_t
bftc_mem_usage_pr_size(void)
{
  return 0;
}

#endif /* bftc_OS_Linux, bftc_OS_OSF1, ... */

/*
 * \brief Return maximum process memory use (in kB) depending on OS.
 *
 * The returned value is the maximum returned by bftc_mem_usage_pr_size()
 * during the program's lifetime. With memory allocations which return
 * memory to the system (such as the GNU glibc on Linux systems),
 * this value will be correct only if allocation is tracked. This should
 * be the case if malloc hooks are used with the glibc allocation
 * functions (BFT library's default configuration/installation option),
 * but may give results lower than the true maximum in other cases.
 */

size_t
bftc_mem_usage_max_pr_size(void)
{
  (void) bftc_mem_usage_pr_size();

  return _bftc_mem_usage_global_max_pr;
}

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
