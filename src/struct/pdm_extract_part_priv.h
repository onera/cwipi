#ifndef __PDM_EXTRACT_PART_PRIV_H__
#define __PDM_EXTRACT_PART_PRIV_H__

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

/**
 * \struct _pdm_extract_part_t
 * \brief  Define a partition mesh. Arrays are shared
 *
 */

struct _pdm_extract_part_t
{
  int                    dim;
  int                    n_part_in;
  int                    n_part_out;

  PDM_split_dual_t       split_dual_method;
  PDM_ownership_t        ownership;
  PDM_MPI_Comm           comm;

  /* Partitioned view - To do with extract for selected gnum + part_to_part */
  int                 *n_cell;
  int                 *n_face;
  int                 *n_edge;
  int                 *n_vtx;
  int                **pcell_face;
  int                **pcell_face_idx;
  int                **pface_edge;
  int                **pface_edge_idx;
  int                **pedge_vtx;
  int                **pface_vtx;
  int                **pface_vtx_idx;

  PDM_g_num_t        **cell_ln_to_gn;
  PDM_g_num_t        **face_ln_to_gn;
  PDM_g_num_t        **edge_ln_to_gn;
  PDM_g_num_t        **vtx_ln_to_gn;

  double             **pvtx_coord;

  /* Which cell or face is selected */
  int                 *n_extract;
  int                **extract_lnum;

  /* Extracted part (Intermediate distributed result) */
  int                   dn_equi_cell;
  int                   dn_equi_face;
  int                   dn_equi_edge;
  int                   dn_equi_vtx;
  int                  *dequi_cell_face;
  int                  *dequi_cell_face_idx;
  int                  *dequi_face_edge;
  int                  *dequi_face_edge_idx;
  int                  *dequi_edge_vtx;
  int                  *dequi_edge_vtx_idx;
  int                  *dequi_face_vtx;
  int                  *dequi_face_vtx_idx;

  PDM_g_num_t          *dequi_parent_cell_ln_to_gn;
  PDM_g_num_t          *dequi_parent_face_ln_to_gn;
  PDM_g_num_t          *dequi_parent_edge_ln_to_gn;
  PDM_g_num_t          *dequi_parent_vtx_ln_to_gn;

  double               *dequi_vtx_coord;

  PDM_part_to_block_t  *ptb_equi_cell;
  PDM_part_to_block_t  *ptb_equi_face;
  PDM_part_to_block_t  *ptb_equi_edge;
  PDM_part_to_block_t  *ptb_equi_vtx;


  /* Extrated part */
  double             **pextract_vtx_coord;

  PDM_bool_t         *is_owner_connectivity;
  PDM_bool_t         *is_owner_ln_to_gn;
  PDM_bool_t         *is_owner_parent_ln_to_gn;
  PDM_bool_t          is_owner_vtx_coord;

  /* Only for mapping and clear API */
  int                *pextract_n_entity              [PDM_MESH_ENTITY_MAX];
  int               **pextract_connectivity          [PDM_CONNECTIVITY_TYPE_MAX];
  int               **pextract_connectivity_idx      [PDM_CONNECTIVITY_TYPE_MAX];
  PDM_g_num_t       **pextract_entity_ln_to_gn       [PDM_MESH_ENTITY_MAX];
  PDM_g_num_t       **pextract_entity_parent_ln_to_gn[PDM_MESH_ENTITY_MAX];


};


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_EXTRACT_PART_PRIV_H__ */
