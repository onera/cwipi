#ifndef __FVM_CONVERT_ARRAY_H__
#define __FVM_CONVERT_ARRAY_H__

/*============================================================================
 * Functions related to the transformation of data arrays for import
 * or export of meshes and fields.
 *
 * All "reasonable" combinations of datatypes are handled here.
 * (templates would be useful here).
 *============================================================================*/

/*
  This file is part of the "Finite Volume Mesh" library, intended to provide
  finite volume mesh and associated fields I/O and manipulation services.

  Copyright (C) 2011       ONERA

  Copyright (C) 2004-2006  EDF

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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "fvm_defs.h"
#include "fvm_nodal.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
#ifdef FVM_CPPCALLER

namespace fvm {
#else
extern "C" {
#endif
#if 0
}} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Convert an array representation of one type to that of another type, with
 * possible indirection, interlacing, de-interlacing, or change of data
 * dimension (i.e. projection or filling extra dimension with zeroes).
 *
 * Floating point (real or double) source and destination arrays may be
 * multidimensional (interlaced or not), but with and integer type
 * for source or destination, only 1-D arrays are allowed (no use for
 * integer "vector fields" being currently required or apparent).
 *
 * Integer type destination arrays may be converted to floating point
 * (for output formats with no integer datatype, such as EnSight),
 * but floating point values may not be converted to integer values
 * (no use for this operation being currently apparent).
 *
 * parameters:
 *   src_dim          <-- dimension of source data
 *   src_dim_shift    <-- source data dimension shift (start index)
 *   dest_dim         <-- destination data dimension (1 if non interlaced)
 *   src_idx_start    <-- start index in source data
 *   src_idx_end      <-- past-the-end index in source data
 *   src_interlace    <-- indicates if source data is interlaced
 *   src_datatype     <-- source data type (float, double, or int)
 *   dest_datatype    <-- destination data type (float, double, or int)
 *   n_parent_lists   <-- number of parent lists (if parent_num != NULL)
 *   parent_num_shift <-- parent number to value array index shifts;
 *                        size: n_parent_lists
 *   parent_num       <-- if n_parent_lists > 0, parent entity numbers
 *   src_data         <-- array of source arrays (at least one, with one per
 *                        source dimension if non interlaced, times one per
 *                        parent list if multiple parent lists, with
 *                        x_parent_1, y_parent_1, ..., x_parent_2, ...) order
 *   dest_data        --> destination buffer
 *----------------------------------------------------------------------------*/

void
fvm_convert_array(const int                     src_dim,
                  const int                     src_dim_shift,
                  const int                     dest_dim,
                  const fvm_lnum_t              src_idx_start,
                  const fvm_lnum_t              src_idx_end,
                  const fvm_interlace_t         src_interlace,
                  const fvm_datatype_t          src_datatype,
                  const fvm_datatype_t          dest_datatype,
                  const int                     n_parent_lists,
                  const fvm_lnum_t              parent_num_shift[],
                  const fvm_lnum_t              parent_num[],
                  const void             *const src_data[],
                  void                   *const dest_data);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __FVM_CONVERT_ARRAY_H__ */
