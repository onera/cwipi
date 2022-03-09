#ifndef __PDM_PART_MESH_H__
#define __PDM_PART_MESH_H__

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

typedef struct _pdm_part_mesh_t PDM_part_mesh_t;

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/**
 *
 * \brief Build a distributed mesh structure
 *
 * \param [in]   dn_cell             Number of distributed cells
 * \param [in]   dn_face             Number of distributed faces
 * \param [in]   dn_vtx              Number of distributed vertices
 * \param [in]   dn_bnd              Number of boundaries
 * \param [in]   n_join              Number of interfaces with other zones
 *
 * \return     Identifier
 */

PDM_part_mesh_t*
PDM_part_mesh_create
(
 const int             n_part,
       PDM_MPI_Comm    comm
);

void
PDM_part_mesh_free
(
 PDM_part_mesh_t        *dmesh
);


void
PDM_part_mesh_n_entity_set
(
 PDM_part_mesh_t          *pmesh,
 PDM_mesh_entities_t       entity_type,
 int                      *pn_entity
);


void
PDM_part_mesh_n_entity_get
(
 PDM_part_mesh_t          *pmesh,
 PDM_mesh_entities_t       entity_type,
 int                     **pn_entity
);

void
PDM_part_mesh_connectivity_set
(
 PDM_part_mesh_t          *pmesh,
 PDM_connectivity_type_t   connectivity_type,
 int                     **connect,
 int                     **connect_idx,
 PDM_ownership_t           ownership
);

void
PDM_part_mesh_connectivity_get
(
 PDM_part_mesh_t           *pmesh,
 PDM_connectivity_type_t    connectivity_type,
 int                     ***connect,
 int                     ***connect_idx,
 PDM_ownership_t           ownership
);


void
PDM_part_mesh_bound_set
(
 PDM_part_mesh_t          *pmesh,
 PDM_bound_type_t          bound_type,
 int                       n_bound,
 int                     **connect,
 int                     **connect_idx,
 PDM_ownership_t           ownership
);

void
PDM_part_mesh_bound_get
(
 PDM_part_mesh_t           *pmesh,
 PDM_bound_type_t           bound_type,
 int                       *n_bound,
 int                     ***connect,
 int                     ***connect_idx,
 PDM_ownership_t           ownership
);


void
PDM_part_mesh_bound_ln_to_gn_set
(
 PDM_part_mesh_t          *pmesh,
 PDM_bound_type_t          bound_type,
 PDM_g_num_t             **bound_ln_to_gn,
 PDM_ownership_t           ownership
);


void
PDM_part_mesh_bound_ln_to_gn_get
(
 PDM_part_mesh_t           *pmesh,
 PDM_bound_type_t           bound_type,
 PDM_g_num_t            ***bound_ln_to_gn,
 PDM_ownership_t           ownership
);


/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_PART_MESH_H__ */
