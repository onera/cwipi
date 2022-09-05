#ifndef __PDM_READER_GAMMA_H__
#define __PDM_READER_GAMMA_H__
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
#include "pdm_dmesh_nodal.h"

/*----------------------------------------------------------------------------
 * Macro for handling of different symbol names (underscored or not,
 * lowercase or uppercase) between C and Fortran, for link resolution.
 *----------------------------------------------------------------------------*/


#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */



/**
 *
 * \brief Create a dmesh nodal from a file in ASCII GAMMA mesh format
 *
 * \param[in]  comm                MPI communicator
 * \param[in]  filename            Filename
 * \param[in]  fix_orientation_2d  Ensure positive area for 2d faces
 * \param[in]  fix_orientation_3d  Ensure positive volume for 3d cells
 *
 * \return Pointer to PDM_dmesh_nodal object
 *
 */

PDM_dmesh_nodal_t *
PDM_reader_gamma_dmesh_nodal
(
 PDM_MPI_Comm   comm,
 const char    *filename,
 int            fix_orientation_2d,
 int            fix_orientation_3d
 );


#ifdef __cplusplus
}
#endif
#endif