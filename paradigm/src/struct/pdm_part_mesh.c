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

/*============================================================================
 * Interface structure to represent a distributed mesh
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_part_mesh.h"
#include "pdm_part_mesh_priv.h"
#include "pdm_priv.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_dmesh.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Local structure definitions
 *============================================================================*/

/*============================================================================
 * Global variable
 *============================================================================*/

/*=============================================================================
 * Private function definitions
 *============================================================================*/

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief Build a distributed mesh structure
 *
 * \param [in]   dn_cell             Number of distributed cells
 * \param [in]   dn_face             Number of distributed faces
 * \param [in]   dn_vtx              Number of distributed vertices
 * \param [in]   n_bnd               Number of boundaries
 * \param [in]   n_join              Number of interfaces with other zones
 *
 * \return     Identifier
 */

PDM_part_mesh_t*
PDM_part_mesh_create
(
 const int             n_part,
       PDM_MPI_Comm    comm
)
{
  PDM_part_mesh_t *pmesh = (PDM_part_mesh_t *) malloc(sizeof(PDM_part_mesh_t));

  pmesh->n_part = n_part;
  pmesh->comm   = comm;

  pmesh->pconnectivity          = (int         *** ) malloc( PDM_CONNECTIVITY_TYPE_MAX * sizeof(int         **) );
  pmesh->pconnectivity_idx      = (int         *** ) malloc( PDM_CONNECTIVITY_TYPE_MAX * sizeof(int         **) );
  pmesh->pentity_ln_to_gn       = (PDM_g_num_t *** ) malloc( PDM_MESH_ENTITY_MAX       * sizeof(PDM_g_num_t **) );

  pmesh->pn_entity              = (int          ** ) malloc( PDM_MESH_ENTITY_MAX       * sizeof(int          *) );

  pmesh->is_owner_connectivity  = malloc( PDM_CONNECTIVITY_TYPE_MAX * sizeof(PDM_bool_t   ) );
  pmesh->is_owner_ln_to_gn      = malloc( PDM_MESH_ENTITY_MAX       * sizeof(PDM_bool_t   ) );

  for(int i = 0; i < PDM_CONNECTIVITY_TYPE_MAX; ++i) {
    pmesh->is_owner_connectivity[i] = PDM_FALSE;
    pmesh->pconnectivity        [i] = NULL;
    pmesh->pconnectivity_idx    [i] = NULL;
  }

  for(int i = 0; i < PDM_MESH_ENTITY_MAX; ++i) {
    pmesh->is_owner_ln_to_gn[i] = PDM_FALSE;
    pmesh->pentity_ln_to_gn [i] = NULL;
    pmesh->pn_entity        [i] = NULL;
  }

  pmesh->pbound                  = malloc( PDM_BOUND_TYPE_MAX * sizeof(PDM_g_num_t *) );
  pmesh->pbound_idx              = malloc( PDM_BOUND_TYPE_MAX * sizeof(int         *) );
  pmesh->is_owner_bound          = malloc( PDM_BOUND_TYPE_MAX * sizeof(PDM_bool_t   ) );
  pmesh->is_owner_bound_ln_to_gn = malloc( PDM_BOUND_TYPE_MAX * sizeof(PDM_bool_t   ) );

  for(int i = 0; i < PDM_BOUND_TYPE_MAX; ++i ) {
    pmesh->n_group_bnd[i] = 0;
  }

  for(int i = 0; i < PDM_BOUND_TYPE_MAX; ++i) {
    pmesh->is_owner_bound         [i] = PDM_FALSE;
    pmesh->is_owner_bound_ln_to_gn[i] = PDM_FALSE;
    pmesh->pbound                 [i] = NULL;
    pmesh->pbound_idx             [i] = NULL;
  }

  pmesh->vtx_coords = malloc(pmesh->n_part * sizeof(double));
  for(int i_part = 0; i_part < n_part; ++i_part) {
    pmesh->vtx_coords[i_part] = NULL;
  }
  pmesh->is_owner_vtx_coord = PDM_FALSE;

  return pmesh;
}


void
PDM_part_mesh_n_entity_set
(
 PDM_part_mesh_t          *pmesh,
 int                       i_part,
 PDM_mesh_entities_t       entity_type,
 int                       pn_entity
)
{
  if(pmesh->pn_entity[entity_type] == NULL) {
    pmesh->pn_entity[entity_type] = malloc(pmesh->n_part * sizeof(int));
    for(int i = 0; i < pmesh->n_part; ++i) {
      pmesh->pn_entity[entity_type][i] = 0;
    }
  }
  pmesh->pn_entity[entity_type][i_part] = pn_entity;
}


int
PDM_part_mesh_n_entity_get
(
 PDM_part_mesh_t          *pmesh,
 int                       i_part,
 PDM_mesh_entities_t       entity_type
)
{
  if(pmesh->pn_entity[entity_type] != NULL) {
    return pmesh->pn_entity[entity_type][i_part];
  } else {
    return 0;
  }

}


void
PDM_part_mesh_entity_ln_to_gn_set
(
 PDM_part_mesh_t          *pmesh,
 int                       i_part,
 PDM_mesh_entities_t       entity_type,
 PDM_g_num_t              *pentity_ln_to_gn,
 PDM_ownership_t           ownership
)
{
  if(pmesh->pentity_ln_to_gn [entity_type] == NULL) {
    pmesh->pentity_ln_to_gn [entity_type] = malloc(pmesh->n_part * sizeof(int *));
    for(int i = 0; i < pmesh->n_part; ++i) {
      pmesh->pentity_ln_to_gn [entity_type][i] = NULL;
    }
  }
  pmesh->pentity_ln_to_gn [entity_type][i_part] = pentity_ln_to_gn;
  if(ownership == PDM_OWNERSHIP_USER || ownership == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE) {
    pmesh->is_owner_ln_to_gn[entity_type] = PDM_FALSE;
  } else {
    pmesh->is_owner_ln_to_gn[entity_type] = PDM_TRUE;
  }
}


void
PDM_part_mesh_entity_ln_to_gn_get
(
 PDM_part_mesh_t          *pmesh,
 int                       i_part,
 PDM_mesh_entities_t       entity_type,
 PDM_g_num_t             **pentity_ln_to_gn,
 PDM_ownership_t           ownership
)
{
  if(pmesh->pentity_ln_to_gn[entity_type] != NULL) {
    *pentity_ln_to_gn = pmesh->pentity_ln_to_gn[entity_type][i_part];
  } else {
    *pentity_ln_to_gn = NULL;
  }
  if(ownership == PDM_OWNERSHIP_USER || ownership == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE) {
    pmesh->is_owner_ln_to_gn[entity_type] = PDM_FALSE;
  } else {
    pmesh->is_owner_ln_to_gn[entity_type] = PDM_TRUE;
  }
}



void
PDM_part_mesh_connectivity_set
(
 PDM_part_mesh_t          *pmesh,
 int                       i_part,
 PDM_connectivity_type_t   connectivity_type,
 int                      *connect,
 int                      *connect_idx,
 PDM_ownership_t           ownership
)
{
  if(pmesh->pconnectivity [connectivity_type] == NULL) {
    pmesh->pconnectivity    [connectivity_type] = malloc(pmesh->n_part * sizeof(int *));
    pmesh->pconnectivity_idx[connectivity_type] = malloc(pmesh->n_part * sizeof(int *));
    for(int i = 0; i < pmesh->n_part; ++i) {
      pmesh->pconnectivity    [connectivity_type][i] = NULL;
      pmesh->pconnectivity_idx[connectivity_type][i] = NULL;
    }
  }
  pmesh->pconnectivity        [connectivity_type][i_part] = connect;
  pmesh->pconnectivity_idx    [connectivity_type][i_part] = connect_idx;

  if(ownership == PDM_OWNERSHIP_USER || ownership == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE) {
    pmesh->is_owner_connectivity[connectivity_type] = PDM_FALSE;
  } else {
    pmesh->is_owner_connectivity[connectivity_type] = PDM_TRUE;
  }
}

void
PDM_part_mesh_vtx_coord_set
(
 PDM_part_mesh_t   *pmesh,
 int                i_part,
 double            *vtx_coord,
 PDM_ownership_t    ownership
)
{
  pmesh->vtx_coords[i_part] = vtx_coord;
  if(ownership == PDM_OWNERSHIP_USER || ownership == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE) {
    pmesh->is_owner_vtx_coord = PDM_FALSE;
  } else {
    pmesh->is_owner_vtx_coord = PDM_TRUE;
  }
}

void
PDM_part_mesh_vtx_coord_get
(
 PDM_part_mesh_t   *pmesh,
 int                i_part,
 double           **vtx_coord,
 PDM_ownership_t    ownership
)
{
  if(pmesh->vtx_coords[i_part] != NULL) {
    *vtx_coord = pmesh->vtx_coords[i_part];
  } else {
    *vtx_coord = NULL;
  }
  if(ownership == PDM_OWNERSHIP_USER || ownership == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE) {
    pmesh->is_owner_vtx_coord = PDM_FALSE;
  } else {
    pmesh->is_owner_vtx_coord = PDM_TRUE;
  }
}

void
PDM_part_mesh_connectivity_get
(
 PDM_part_mesh_t           *pmesh,
 int                        i_part,
 PDM_connectivity_type_t    connectivity_type,
 int                      **connect,
 int                      **connect_idx,
 PDM_ownership_t           ownership
)
{
  if(pmesh->pconnectivity    [connectivity_type][i_part] != NULL) {
    *connect     = pmesh->pconnectivity    [connectivity_type][i_part];
    *connect_idx = pmesh->pconnectivity_idx[connectivity_type][i_part];
  } else {
    *connect     = NULL;
    *connect_idx = NULL;
  }

  if(ownership == PDM_OWNERSHIP_USER || ownership == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE) {
    pmesh->is_owner_connectivity[connectivity_type] = PDM_FALSE;
  } else {
    pmesh->is_owner_connectivity[connectivity_type] = PDM_TRUE;
  }
}


void
PDM_part_mesh_bound_set
(
 PDM_part_mesh_t          *pmesh,
 PDM_bound_type_t          bound_type,
 int                       n_bound,
 int                     **connect,
 int                     **connect_idx,
 PDM_ownership_t           ownership
)
{
  pmesh->n_group_bnd   [bound_type] = n_bound;
  pmesh->pbound        [bound_type] = connect;
  pmesh->pbound_idx    [bound_type] = connect_idx;

  if(ownership == PDM_OWNERSHIP_USER || ownership == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE) {
    pmesh->is_owner_bound[bound_type] = PDM_FALSE;
  } else {
    pmesh->is_owner_bound[bound_type] = PDM_TRUE;
  }
}


void
PDM_part_mesh_bound_get
(
 PDM_part_mesh_t           *pmesh,
 PDM_bound_type_t           bound_type,
 int                       *n_bound,
 int                     ***connect,
 int                     ***connect_idx,
 PDM_ownership_t           ownership
)
{
  *n_bound     = pmesh->n_group_bnd[bound_type];
  *connect     = pmesh->pbound     [bound_type];
  *connect_idx = pmesh->pbound_idx [bound_type];

  if(ownership == PDM_OWNERSHIP_USER || ownership == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE) {
    pmesh->is_owner_bound[bound_type] = PDM_FALSE;
  } else {
    pmesh->is_owner_bound[bound_type] = PDM_TRUE;
  }
}



void
PDM_part_mesh_bound_ln_to_gn_set
(
 PDM_part_mesh_t          *pmesh,
 PDM_bound_type_t          bound_type,
 PDM_g_num_t             **bound_ln_to_gn,
 PDM_ownership_t           ownership
)
{
  pmesh->pbound_ln_to_gn        [bound_type] = bound_ln_to_gn;
  if(ownership == PDM_OWNERSHIP_USER || ownership == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE) {
    pmesh->is_owner_bound_ln_to_gn[bound_type] = PDM_FALSE;
  } else {
    pmesh->is_owner_bound_ln_to_gn[bound_type] = PDM_TRUE;

  }
}


void
PDM_part_mesh_bound_ln_to_gn_get
(
 PDM_part_mesh_t           *pmesh,
 PDM_bound_type_t           bound_type,
 PDM_g_num_t            ***bound_ln_to_gn,
 PDM_ownership_t           ownership
)
{
  *bound_ln_to_gn = pmesh->pbound_ln_to_gn[bound_type];

  if(ownership == PDM_OWNERSHIP_USER || ownership == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE) {
    pmesh->is_owner_bound_ln_to_gn[bound_type] = PDM_FALSE;
  } else {
    pmesh->is_owner_bound_ln_to_gn[bound_type] = PDM_TRUE;
  }
}




void
PDM_part_mesh_free
(
 PDM_part_mesh_t        *pmesh
)
{

  /* Free connectivity */
  for(int i = 0; i < PDM_CONNECTIVITY_TYPE_MAX; ++i) {
    if(pmesh->is_owner_connectivity[i] == PDM_TRUE) {
      for(int i_part = 0; i_part < pmesh->n_part; ++i_part) {
        if(pmesh->pconnectivity[i][i_part] != NULL) {
          free(pmesh->pconnectivity[i][i_part]);
        }
        if(pmesh->pconnectivity_idx[i][i_part] != NULL) {
          free(pmesh->pconnectivity_idx[i][i_part]);
        }
      }
    }

    if(pmesh->pconnectivity[i] != NULL) {
      free(pmesh->pconnectivity[i]);
      pmesh->pconnectivity[i] = NULL;
    }

    if(pmesh->pconnectivity_idx[i] != NULL) {
      free(pmesh->pconnectivity_idx[i]);
      pmesh->pconnectivity_idx[i] = NULL;
    }
  }

  /* Free ln_to_gn */
  for(int i = 0; i < PDM_MESH_ENTITY_MAX; ++i) {
    if(pmesh->is_owner_ln_to_gn[i] == PDM_TRUE) {
      for(int i_part = 0; i_part < pmesh->n_part; ++i_part) {
        if(pmesh->pentity_ln_to_gn[i][i_part] != NULL) {
          free(pmesh->pentity_ln_to_gn[i][i_part]);
        }
      }
    }

    if(pmesh->pentity_ln_to_gn[i] != NULL) {
      free(pmesh->pentity_ln_to_gn[i]);
      pmesh->pentity_ln_to_gn[i] = NULL;
    }
  }

  /* Free vtx__coord */
  if(pmesh->is_owner_vtx_coord == PDM_TRUE) {
    for(int i_part = 0; i_part < pmesh->n_part; ++i_part) {
      free(pmesh->vtx_coords[i_part]);
    }
  }
  free(pmesh->vtx_coords);

  /* Free group */
  for(int i = 0; i < PDM_BOUND_TYPE_MAX; ++i) {
    if(pmesh->is_owner_bound[i] == PDM_TRUE) {
      for(int i_part = 0; i_part < pmesh->n_part; ++i_part) {
        if(pmesh->pbound[i][i_part] != NULL) {
          free(pmesh->pbound[i][i_part]);
        }
        if(pmesh->pbound_idx[i][i_part] != NULL) {
          free(pmesh->pbound_idx[i][i_part]);
        }
      }

      if(pmesh->pbound[i] != NULL) {
        free(pmesh->pbound[i]);
        pmesh->pbound[i] = NULL;
      }

      if(pmesh->pbound_idx[i] != NULL) {
        free(pmesh->pbound_idx[i]);
        pmesh->pbound_idx[i] = NULL;
      }
    }

  }

  /* Free bound_ln_to_gn */
  for(int i = 0; i < PDM_BOUND_TYPE_MAX; ++i) {
    if(pmesh->is_owner_bound_ln_to_gn[i] == PDM_TRUE) {
      for(int i_part = 0; i_part < pmesh->n_part; ++i_part) {
        if(pmesh->pbound_ln_to_gn[i][i_part] != NULL) {
          free(pmesh->pbound_ln_to_gn[i][i_part]);
        }
      }

      if(pmesh->pbound_ln_to_gn[i] != NULL) {
        free(pmesh->pbound_ln_to_gn[i]);
        pmesh->pbound_ln_to_gn[i] = NULL;
      }
    }
  }


  for(int i = 0; i < PDM_MESH_ENTITY_MAX; ++i) {
    if(pmesh->pn_entity[i] !=NULL){
      free(pmesh->pn_entity[i]);
    }
  }

  free(pmesh->pn_entity);
  free(pmesh->pconnectivity);
  free(pmesh->pconnectivity_idx);
  free(pmesh->pentity_ln_to_gn);
  free(pmesh->is_owner_connectivity);
  free(pmesh->is_owner_ln_to_gn    );
  free(pmesh->pbound                 );
  free(pmesh->pbound_idx             );
  free(pmesh->is_owner_bound         );
  free(pmesh->is_owner_bound_ln_to_gn);


  free(pmesh);
}

