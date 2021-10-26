#ifndef __PDM_POLY_VOL_GEN_H__
#define __PDM_POLY_VOL_GEN_H__
/*
  This file is part of the CWIPI library.

  Copyright (C) 2011  ONERA

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

#include "pdm.h"
#include "pdm_mpi.h"

/*----------------------------------------------------------------------------
 * Macro for handling of different symbol names (underscored or not,
 * lowercase or uppercase) between C and Fortran, for link resolution.
 *----------------------------------------------------------------------------*/

//#define PROCF(x, y) x
#if !defined (__hpux) && !defined (_AIX)
#define PROCF(x, y) x##_
#else
#define PROCF(x, y) x
#endif

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

void
PDM_poly_vol_gen
(
 PDM_MPI_Comm  comm,
 double        xmin,
 double        ymin,
 double        zmin,
 double        lengthx,
 double        lengthy,
 double        lengthz,
 PDM_g_num_t   nx,
 PDM_g_num_t   ny,
 PDM_g_num_t   nz,
 int           randomize,
 int           random_seed,
 PDM_g_num_t  *ng_cell,
 PDM_g_num_t  *ng_face,
 PDM_g_num_t  *ng_vtx,
 int          *n_face_group,
 int          *dn_cell,
 int          *dn_face,
 int          *dn_vtx,
 int         **dcell_face_idx,
 PDM_g_num_t **dcell_face,
 PDM_g_num_t **dface_cell,
 int         **dface_vtx_idx,
 PDM_g_num_t **dface_vtx,
 double      **dvtx_coord,
 int         **dface_group_idx,
 PDM_g_num_t **dface_group
 );

#ifdef __cplusplus
}
#endif
#endif
