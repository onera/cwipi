#ifndef __FVM_LOCATOR_H__
#define __FVM_LOCATOR_H__

/*============================================================================
 * Locate points in a nodal representation associated with a mesh
 *============================================================================*/

/*
  This file is part of the "Finite Volume Mesh" library, intended to provide
  finite volume mesh and associated fields I/O and manipulation services.

  Copyright (C) 2005-2008  EDF

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

/*----------------------------------------------------------------------------*/

#include "fvm_config.h"

#if defined(FVM_HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "fvm_defs.h"
#include "fvm_nodal.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Function pointer type for user definable logging/profiling type functions
 *----------------------------------------------------------------------------*/

typedef int
(fvm_locator_log_t) (int         event,
                     int         data,
                     const char *string);

/*----------------------------------------------------------------------------
 * Structure defining a locator
 *----------------------------------------------------------------------------*/

typedef struct _fvm_locator_t fvm_locator_t;

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Creation of a locator structure.
 *
 * Note that depending on the choice of ranks of the associated communicator,
 * distant ranks may in fact be truly distant or not. If n_ranks = 1 and
 * start_rank is equal to the current rank in the communicator, the locator
 * will work only locally.
 *
 * parameters:
 *   tolerance  <-- addition to local extents of each element:
 *                  extent = base_extent * (1 + tolerance)
 *   comm       <-- associated MPI communicator
 *   n_ranks    <-- number of MPI ranks associated with distant location
 *   start_rank <-- first MPI rank associated with distant location
 *
 * returns:
 *   pointer to locator
 *----------------------------------------------------------------------------*/

#if defined(FVM_HAVE_MPI)

fvm_locator_t *
fvm_locator_create(double    tolerance,
                   MPI_Comm  comm,
                   int       n_ranks,
                   int       start_rank);

#else

fvm_locator_t *
fvm_locator_create(double    tolerance);

#endif

/*----------------------------------------------------------------------------
 * Destruction of a locator structure.
 *
 * parameters:
 *   this_locator <-> locator to destroy
 *
 * returns:
 *   NULL pointer
 *----------------------------------------------------------------------------*/

fvm_locator_t *
fvm_locator_destroy(fvm_locator_t  * this_locator);

/*----------------------------------------------------------------------------
 * Prepare locator for use with a given nodal mesh representation.
 *
 * parameters:
 *   this_locator      <-> pointer to locator structure
 *   this_nodal        <-- pointer to mesh representation structure
 *   locate_on_parents <-- location relative to parent element numbers if
 *                         1, id of element + 1 in concatenated sections
 *                         of same element dimension if 0
 *   dim               <-- space dimension of points to locate
 *   n_points          <-- number of points to locate
 *   point_list        <-- optional indirection array to point_coords
 *                         (1 to n_points numbering)
 *   point_coords      <-- coordinates of points to locate
 *                         (dimension: dim * n_points)
 *----------------------------------------------------------------------------*/

void
fvm_locator_set_nodal(fvm_locator_t       *this_locator,
                      const fvm_nodal_t   *this_nodal,
                      int                  locate_on_parents,
                      int                  dim,
                      fvm_lnum_t           n_points,
                      const fvm_lnum_t     point_list[],
                      const fvm_coord_t    point_coords[]);

/*----------------------------------------------------------------------------
 * Return number of distant points after locator initialization.
 *
 * parameters:
 *   this_locator <-- pointer to locator structure
 *
 * returns:
 *   number of distant points.
 *----------------------------------------------------------------------------*/

fvm_lnum_t
fvm_locator_get_n_dist_points(const fvm_locator_t  *this_locator);

/*----------------------------------------------------------------------------
 * Return an array of local element numbers containing (or nearest to)
 * each distant point after locator initialization.
 *
 * parameters:
 *   this_locator <-- pointer to locator structure
 *
 * returns:
 *   local element numbers associated with distant points (1 to n numbering).
 *----------------------------------------------------------------------------*/

const fvm_lnum_t *
fvm_locator_get_dist_locations(const fvm_locator_t  *this_locator);

/*----------------------------------------------------------------------------
 * Return an array of coordinates of each distant point after
 * locator initialization.
 *
 * parameters:
 *   this_locator <-- pointer to locator structure
 *
 * returns:
 *   coordinate array associated with distant points (interlaced).
 *----------------------------------------------------------------------------*/

const fvm_coord_t *
fvm_locator_get_dist_coords(const fvm_locator_t  *this_locator);

/*----------------------------------------------------------------------------
 * Return number of points located after locator initialization.
 *
 * parameters:
 *   this_locator <-- pointer to locator structure
 *
 * returns:
 *   number of points located.
 *----------------------------------------------------------------------------*/

fvm_lnum_t
fvm_locator_get_n_interior(const fvm_locator_t  *this_locator);

/*----------------------------------------------------------------------------
 * Return list of points located after locator initialization.
 * This list defines a subset of the point set used at initialization.
 *
 * parameters:
 *   this_locator <-- pointer to locator structure
 *
 * returns:
 *   list of points located (1 to n numbering).
 *----------------------------------------------------------------------------*/

const fvm_lnum_t *
fvm_locator_get_interior_list(const fvm_locator_t  *this_locator);

/*----------------------------------------------------------------------------
 * Return number of points not located after locator initialization.
 *
 * parameters:
 *   this_locator <-- pointer to locator structure
 *
 * returns:
 *   number of points not located.
 *----------------------------------------------------------------------------*/

fvm_lnum_t
fvm_locator_get_n_exterior(const fvm_locator_t  *this_locator);

/*----------------------------------------------------------------------------
 * Return list of points not located after locator initialization.
 * This list defines a subset of the point set used at initialization.
 *
 * parameters:
 *   this_locator <-- pointer to locator structure
 *
 * returns:
 *   list of points not located (1 to n numbering).
 *----------------------------------------------------------------------------*/

const fvm_lnum_t *
fvm_locator_get_exterior_list(const fvm_locator_t  *this_locator);

/*----------------------------------------------------------------------------
 * Discard list of points not located after locator initialization.
 * This list defines a subset of the point set used at initialization.
 *
 * parameters:
 *   this_locator <-- pointer to locator structure
 *
 * returns:
 *   list of points not located (1 to n numbering).
 *----------------------------------------------------------------------------*/

void
fvm_locator_discard_exterior(fvm_locator_t  *this_locator);

#if defined(FVM_HAVE_MPI)

/*----------------------------------------------------------------------------
 * Non blocking send dist_var or (local_var if reverse == true)
 *
 * parameters:
 *   this_locator  <-- pointer to locator structure
 *   var           <-- variable defined on distant points (distant_var)
 *                     or variable defined on local points (local_var)
 *                     size: n_dist_points*stride
 *   local_list    <-- optional indirection list (1 to n) for local_var
 *   type_size     <-- sizeof (float or double) variable type
 *   stride        <-- dimension (1 for scalar, 3 for interlaced vector)
 *   reverse       <-- if nonzero, exchange is reversed
 *   tag           <-- tag for MPI_issend
 *   request       <-> communication request
 *----------------------------------------------------------------------------*/

void
fvm_locator_issend_point_var(fvm_locator_t     *this_locator,
                             void              *var,
                             const fvm_lnum_t  *local_list,
                             size_t             type_size,
                             size_t             stride,
                             int                reverse,
                             int                tag,
                             int               *request);


/*----------------------------------------------------------------------------
 * Wait for fvm_locator_iissed
 *
 * parameters:
 *   this_locator  <-- pointer to locator structure
 *   request       <-> communication request
 *----------------------------------------------------------------------------*/

void
fvm_locator_issend_wait(fvm_locator_t     *this_locator,
                        int               request);

/*----------------------------------------------------------------------------
 * Non blocking receive dist_var or (local_var if reverse == true)
 *
 * parameters:
 *   this_locator  <-- pointer to locator structure
 *   var           --> variable defined on distant points (distant_var)
 *                     or variable defined on local points (local_var)
 *                     size: n_dist_points*stride
 *   local_list    <-- optional indirection list (1 to n) for local_var
 *   type_size     <-- sizeof (float or double) variable type
 *   stride        <-- dimension (1 for scalar, 3 for interlaced vector)
 *   reverse       <-- if nonzero, exchange is reversed
 *   tag           <-- tag for MPI_irecv
 *   request       <-> communication request
 *----------------------------------------------------------------------------*/

void
fvm_locator_irecv_point_var(fvm_locator_t     *this_locator,
                            void              *var,
                            const fvm_lnum_t  *local_list,
                            size_t             type_size,
                            size_t             stride,
                            int                reverse,
                            int                tag,
                            int               *request);

/*----------------------------------------------------------------------------
 * Wait for fvm_locator_irecv
 *
 * parameters:
 *   this_locator  <-- pointer to locator structure
 *   request       <-> communication request
 *----------------------------------------------------------------------------*/

void
fvm_locator_irecv_wait(fvm_locator_t     *this_locator,
                       int                request);

#endif

/*----------------------------------------------------------------------------
 * Distribute variable defined on distant points to processes owning
 * the original points (i.e. distant processes).
 *
 * The exchange is symmetric if both variables are defined, receive
 * only if distant_var is NULL, or send only if local_var is NULL.
 *
 * The caller should have defined the values of distant_var[] for the
 * distant points, whose coordinates are given by
 * fvm_locator_get_dist_coords(), and which are located in the elements
 * whose numbers are given by fvm_locator_get_dist_locations().
 *
 * The local_var[] is defined at the located points (those whose
 * numbers are returned by fvm_locator_get_interior_list().
 *
 * parameters:
 *   this_locator  <-- pointer to locator structure
 *   distant_var   <-> variable defined on distant points (ready to send)
 *                     size: n_dist_points*stride
 *   local_var     <-> variable defined on located local points (received)
 *                     size: n_interior*stride
 *   local_list    <-- optional indirection list (1 to n) for local_var
 *   type_size     <-- sizeof (float or double) variable type
 *   stride        <-- dimension (1 for scalar, 3 for interlaced vector)
 *   reverse       <-- if nonzero, exchange is reversed
 *                     (receive values associated with distant points
 *                     from the processes owning the original points)
 *----------------------------------------------------------------------------*/

void
fvm_locator_exchange_point_var(fvm_locator_t     *this_locator,
                               void              *distant_var,
                               void              *local_var,
                               const fvm_lnum_t  *local_list,
                               size_t             type_size,
                               size_t             stride,
                               int                reverse);

/*----------------------------------------------------------------------------
 * Return timing information.
 *
 * In parallel mode, this includes communication time.
 *
 * parameters:
 *   this_locator      <-- pointer to locator structure
 *   location_wtime    --> Location Wall-clock time (or NULL)
 *   location_cpu_time --> Location CPU time (or NULL)
 *   exchange_wtime    --> Variable exchange Wall-clock time (or NULL)
 *   exchange_cpu_time --> Variable exchange CPU time (or NULL)
 *   issend_wtime      --> Variable exchange Wall-clock time (or NULL)
 *   issend_cpu_time   --> Variable exchange CPU time (or NULL)
 *   irecv_wtime       --> Variable exchange Wall-clock time (or NULL)
 *   irecv_cpu_time    --> Variable exchange CPU time (or NULL)
 *----------------------------------------------------------------------------*/

void
fvm_locator_get_times(const fvm_locator_t  *this_locator,
                      double               *location_wtime,
                      double               *location_cpu_time,
                      double               *exchange_wtime,
                      double               *exchange_cpu_time,
                      double               *issend_wtime,
                      double               *issend_cpu_time,
                      double               *irecv_wtime,
                      double               *irecv_cpu_time);


/*----------------------------------------------------------------------------
 * Return communication timing information.
 *
 * In serial mode, returned times are always zero..
 *
 * parameters:
 *   this_locator      <-- pointer to locator structure
 *   location_wtime    --> Location Wall-clock time (or NULL)
 *   location_cpu_time --> Location CPU time (or NULL)
 *   exchange_wtime    --> Variable exchange Wall-clock time (or NULL)
 *   exchange_cpu_time --> Variable exchange CPU time (or NULL)
 *   issend_wtime      --> Variable exchange Wall-clock time (or NULL)
 *   issend_cpu_time   --> Variable exchange CPU time (or NULL)
 *   irecv_wtime       --> Variable exchange Wall-clock time (or NULL)
 *   irecv_cpu_time    --> Variable exchange CPU time (or NULL)
 *----------------------------------------------------------------------------*/

void
fvm_locator_get_comm_times(const fvm_locator_t  *this_locator,
                           double               *location_wtime,
                           double               *location_cpu_time,
                           double               *exchange_wtime,
                           double               *exchange_cpu_time,
                           double               *issend_wtime,
                           double               *issend_cpu_time,
                           double               *irecv_wtime,
                           double               *irecv_cpu_time);
 

/*----------------------------------------------------------------------------
 * Dump printout of a locator structure.
 *
 * parameters:
 *   this_locator  <-- pointer to structure that should be dumped
 *----------------------------------------------------------------------------*/

void
fvm_locator_dump(const fvm_locator_t  *this_locator);

/*----------------------------------------------------------------------------
 * Register communication logging functions for locator instrumentation.
 *
 * By default, locators are not instrumented.

 * Functions using MPE may be defined and used, but other similar systems
 * may be used.
 *
 * parameters:
 *   fct           <-- pointer to logging function
 *   start_p_comm  <-- point to point communication start event number
 *   end_p_comm    <-- point to point communication end event number
 *   start_g_comm  <-- global communication start event number
 *   end_g_comm    <-- global communication end event number
 *----------------------------------------------------------------------------*/

#if defined(FVM_HAVE_MPI)

void
fvm_locator_set_comm_log(fvm_locator_log_t  *log_function,
                         int                 start_p_comm,
                         int                 end_p_comm,
                         int                 start_g_comm,
                         int                 end_g_comm);

#endif /* defined(FVM_HAVE_MPI) */

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __FVM_LOCATOR_H__ */
