#ifndef __PDM_EXTRACT_PART_H__
#define __PDM_EXTRACT_PART_H__

/*
  This file is part of the ParaDiGM library.

  Copyright (C) 2017       ONERA

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

#include "pdm.h"
#include "pdm_mpi.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

typedef struct _pdm_extract_part_t PDM_extract_part_t;


/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/**
 *
 * \brief Build an iso surface struct
 *
 *
 * \return     Identifier
 */
PDM_extract_part_t*
PDM_extract_part_create
(
 const int                    dim,
 const int                    n_part_in,
 const int                    n_part_out,
       PDM_split_dual_t       split_dual_method,
       PDM_ownership_t        ownership,
       PDM_MPI_Comm           comm
);


void
PDM_extract_part_compute
(
  PDM_extract_part_t        *extrp
);


void
PDM_extract_part_selected_lnum_set
(
  PDM_extract_part_t       *extrp,
  int                       i_part,
  int                       n_extract,
  int                      *extract_lnum
);

void
PDM_extract_part_part_set
(
  PDM_extract_part_t        *extrp,
  int                       i_part,
  int                       n_cell,
  int                       n_face,
  int                       n_edge,
  int                       n_vtx,
  int                      *cell_face_idx,
  int                      *cell_face,
  int                      *face_edge_idx,
  int                      *face_edge,
  int                      *edge_vtx,
  int                      *face_vtx_idx,
  int                      *face_vtx,
  PDM_g_num_t              *cell_ln_to_gn,
  PDM_g_num_t              *face_ln_to_gn,
  PDM_g_num_t              *edge_ln_to_gn,
  PDM_g_num_t              *vtx_ln_to_gn,
  double                   *vtx_coord
);



void
PDM_extract_part_n_entity_get
(
 PDM_extract_part_t       *extrp,
 PDM_mesh_entities_t       entity_type,
 int                     **pn_entity
);


void
PDM_extract_part_connectivity_get
(
 PDM_extract_part_t        *extrp,
 PDM_connectivity_type_t    connectivity_type,
 int                     ***connect,
 int                     ***connect_idx,
 PDM_ownership_t           ownership
);


void
PDM_extract_part_ln_to_gn_get
(
 PDM_extract_part_t        *extrp,
 PDM_mesh_entities_t       entity_type,
 PDM_g_num_t            ***pentity_ln_to_gn,
 PDM_ownership_t           ownership
);


void
PDM_extract_part_parent_ln_to_gn_get
(
 PDM_extract_part_t        *extrp,
 PDM_mesh_entities_t       entity_type,
 PDM_g_num_t            ***pentity_parent_ln_to_gn,
 PDM_ownership_t           ownership
);

void
PDM_extract_part_vtx_coord_get
(
 PDM_extract_part_t         *extrp,
 double                   ***pvtx_coord,
 PDM_ownership_t           ownership
);

void
PDM_extract_part_free
(
  PDM_extract_part_t  *extrp
);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_EXTRACT_PART_H__ */
