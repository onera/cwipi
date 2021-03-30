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
#include "pdm_dmesh_priv.h"
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

PDM_dmesh_t*
PDM_dmesh_create
(
       PDM_ownership_t owner,
 const int             dn_cell,
 const int             dn_face,
 const int             dn_edge,
 const int             dn_vtx,
 const int             n_bnd,
 const int             n_join,
       PDM_MPI_Comm    comm
)
{
  PDM_dmesh_t *dmesh = (PDM_dmesh_t *) malloc(sizeof(PDM_dmesh_t));

  dmesh->comm              = comm;
  dmesh->owner             = owner;
  dmesh->results_is_getted = malloc( PDM_CONNECTIVITY_TYPE_MAX * sizeof(PDM_bool_t) );
  dmesh->dn_cell           = dn_cell;
  dmesh->dn_face           = dn_face;
  dmesh->dn_edge           = dn_edge;
  dmesh->dn_vtx            = dn_vtx;
  dmesh->n_bnd             = n_bnd;
  dmesh->n_join            = n_join;

  dmesh->_dface_cell       = NULL;
  dmesh->_dface_vtx_idx    = NULL;
  dmesh->_dface_vtx        = NULL;
  dmesh->_dvtx_coord       = NULL;

  dmesh->_dedge_vtx_idx    = NULL;
  dmesh->_dedge_vtx        = NULL;

  dmesh->_dedge_face_idx   = NULL;
  dmesh->_dedge_face       = NULL;

  dmesh->_dedge_bound_idx  = NULL;
  dmesh->_dedge_bound      = NULL;

  dmesh->_dface_bound_idx  = NULL;
  dmesh->_dface_bound      = NULL;
  dmesh->_joins_glob_id    = NULL;
  dmesh->_dface_join_idx   = NULL;
  dmesh->_dface_join       = NULL;

  dmesh->cell_distrib       = NULL;
  dmesh->face_distrib       = NULL;
  dmesh->edge_distrib       = NULL;
  dmesh->vtx_distrib        = NULL;

  dmesh->dconnectivity          = malloc( PDM_CONNECTIVITY_TYPE_MAX * sizeof(PDM_g_num_t *) );
  dmesh->dconnectivity_idx      = malloc( PDM_CONNECTIVITY_TYPE_MAX * sizeof(int         *) );
  dmesh->is_owner_connectivity  = malloc( PDM_CONNECTIVITY_TYPE_MAX * sizeof(PDM_bool_t   ) );

  for(int i = 0; i < PDM_CONNECTIVITY_TYPE_MAX; ++i) {
    dmesh->is_owner_connectivity[i] = PDM_FALSE;
    dmesh->dconnectivity        [i] = NULL;
    dmesh->dconnectivity_idx    [i] = NULL;
  }

  dmesh->dbound          = malloc( PDM_BOUND_TYPE_MAX * sizeof(PDM_g_num_t *) );
  dmesh->dbound_idx      = malloc( PDM_BOUND_TYPE_MAX * sizeof(int         *) );
  dmesh->is_owner_bound  = malloc( PDM_BOUND_TYPE_MAX * sizeof(PDM_bool_t   ) );

  for(int i = 0; i < PDM_BOUND_TYPE_MAX; ++i) {
    dmesh->is_owner_bound[i] = PDM_FALSE;
    dmesh->dbound        [i] = NULL;
    dmesh->dbound_idx    [i] = NULL;
  }

  return dmesh;
}

/**
 *
 * \brief Set the arrays into the distributed mesh structure
 *
 * \param [in]   id                 id of the dmesh to be set
 * \param [in]   dvtx_coord          Coordinates of  vertices (size = 3 * dn_vtx)
 * \param [in]   dface_vtx_idx        Face-vertex connectivity index of
 *                                    faces (size = dn_face + 1)
 * \param [in]   dface_vtx           Face-vertex connectivity of faces
 *                                    (size = dface_vtx_idx[dn_face])
 * \param [in]   dface_cell          Face-cell connectivity of faces (size =
 *                                    2 * dn_face). If iface is a boundary face,
 *                                    dface_cell[2*iface + 1] = 0
 * \param [in]   dface_bound_idx      Index of faces list of each boundary
 *                                    (size = n_bnd + 1)
 * \param [in]   dface_bound         Faces list of each boundary
 *                                    (size = dface_bound_idx[n_bnd])
 * \param [in]   joins_glob_id       Global id of each join (size = n_join)
 * \param [in]   dface_join_idx       Index of faces list of each join
 *                                    (size = n_join + 1)
 * \param [in]   dface_join          Faces list of each join
 *                                    (size = dface_join_idx[n_join])
 */

void
PDM_dmesh_set
(
       PDM_dmesh_t  *dmesh,
 const double       *dvtx_coord,
 const int          *dface_vtx_idx,
 const PDM_g_num_t  *dface_vtx,
 const PDM_g_num_t  *dface_cell,
 const int          *dface_bound_idx,
 const PDM_g_num_t  *dface_bound,
 const int          *joins_glob_id,
 const int          *dface_join_idx,
 const PDM_g_num_t  *dface_join
)
{
  dmesh->_dvtx_coord      = (double      *) dvtx_coord;
  dmesh->_dface_vtx_idx   = (int         *) dface_vtx_idx;
  dmesh->_dface_vtx       = (PDM_g_num_t *) dface_vtx;
  dmesh->_dface_cell      = (PDM_g_num_t *) dface_cell;
  dmesh->_dface_bound_idx = (int         *) dface_bound_idx;
  dmesh->_dface_bound     = (PDM_g_num_t *) dface_bound;
  dmesh->_joins_glob_id   = (int         *) joins_glob_id;
  dmesh->_dface_join_idx  = (int         *) dface_join_idx;
  dmesh->_dface_join      = (PDM_g_num_t *) dface_join;

}

/**
 *
 * \brief Get the dimensions of the distributed mesh
 *
 * \param [in]    id                id of the dmesh requested
 * \param [out]   dn_cell            Number of distributed cells
 * \param [out]   dn_face            Number of distributed faces
 * \param [out]   dn_vtx             Number of distributed vertices
 * \param [out]   n_bnd              Number of boundaries
 * \param [out]   n_join             Number of interfaces with other zones
 */

void
PDM_dmesh_dims_get
(
 PDM_dmesh_t *dmesh,
 int         *dn_cell,
 int         *dn_face,
 int         *dn_edge,
 int         *dn_vtx,
 int         *n_bnd,
 int         *n_join
)
{
  *dn_cell = dmesh->dn_cell;
  *dn_face = dmesh->dn_face;
  *dn_edge = dmesh->dn_edge;
  *dn_vtx  = dmesh->dn_vtx;
  *n_bnd   = dmesh->n_bnd;
  *n_join  = dmesh->n_join;
}

/**
 *
 * \brief Get the data (arrays) of the distributed mesh
 *
 * \param [in]    id                 id of the requested dmesh
 * \param [out]   dvtx_coord          Coordinates of  vertices
 * \param [out]   dface_vtx_idx        Face-vertex connectivity indices
 * \param [out]   dface_vtx           Face-vertex connectivity
 * \param [out]   dface_cell          Face-cell connectivity of faces
 * \param [out]   dface_bound_idx      Indices of faces list of each boundary
 * \param [out]   dface_bound         Faces list of each boundary
 * \param [out]   joins_glob_id       Global Id  of each join
 * \param [out]   dface_join_idx       Indices of faces list of each join
 * \param [out]   dface_join          Faces list of each join
 */

void
PDM_dmesh_data_get
(
       PDM_dmesh_t   *dmesh,
 const double       **dvtx_coord,
 const int          **dface_vtx_idx,
 const PDM_g_num_t  **dface_vtx,
 const PDM_g_num_t  **dface_cell,
 const int          **dface_bound_idx,
 const PDM_g_num_t  **dface_bound,
 const int          **joins_glob_id,
 const int          **dface_join_idx,
 const PDM_g_num_t  **dface_join
)
{
  *dvtx_coord      = dmesh->_dvtx_coord;
  *dface_vtx_idx   = dmesh->_dface_vtx_idx;
  *dface_vtx       = dmesh->_dface_vtx;
  *dface_cell      = dmesh->_dface_cell;
  *dface_bound_idx = dmesh->_dface_bound_idx;
  *dface_bound     = dmesh->_dface_bound;
  *joins_glob_id   = dmesh->_joins_glob_id;
  *dface_join_idx  = dmesh->_dface_join_idx;
  *dface_join      = dmesh->_dface_join;
}


int
PDM_dmesh_connectivity_get
(
 PDM_dmesh_t              *dmesh,
 PDM_connectivity_type_t   connectivity_type,
 PDM_g_num_t             **connect,
 int                     **connect_idx,
 PDM_ownership_t           ownership
)
{
  PDM_UNUSED(ownership);
  assert(dmesh != NULL);

  // assert(dmesh->dconnectivity[connectivity_type] != NULL);

  *connect     = dmesh->dconnectivity    [connectivity_type];
  *connect_idx = dmesh->dconnectivity_idx[connectivity_type];

  if(ownership == PDM_OWNERSHIP_USER || ownership == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE) {
    dmesh->is_owner_connectivity[connectivity_type] = PDM_FALSE;
  }

  int dn_entity = -1;
  if( connectivity_type == PDM_CONNECTIVITY_TYPE_CELL_ELMT ||
      connectivity_type == PDM_CONNECTIVITY_TYPE_CELL_CELL ||
      connectivity_type == PDM_CONNECTIVITY_TYPE_CELL_FACE ||
      connectivity_type == PDM_CONNECTIVITY_TYPE_CELL_EDGE ||
      connectivity_type == PDM_CONNECTIVITY_TYPE_CELL_VTX)
  {
    dn_entity = dmesh->dn_cell;
  } else if( connectivity_type == PDM_CONNECTIVITY_TYPE_FACE_ELMT ||
             connectivity_type == PDM_CONNECTIVITY_TYPE_FACE_CELL ||
             connectivity_type == PDM_CONNECTIVITY_TYPE_FACE_FACE ||
             connectivity_type == PDM_CONNECTIVITY_TYPE_FACE_EDGE ||
             connectivity_type == PDM_CONNECTIVITY_TYPE_FACE_VTX )
  {
    dn_entity = dmesh->dn_face;
  } else if( connectivity_type == PDM_CONNECTIVITY_TYPE_EDGE_ELMT ||
             connectivity_type == PDM_CONNECTIVITY_TYPE_EDGE_CELL ||
             connectivity_type == PDM_CONNECTIVITY_TYPE_EDGE_FACE ||
             connectivity_type == PDM_CONNECTIVITY_TYPE_EDGE_EDGE ||
             connectivity_type == PDM_CONNECTIVITY_TYPE_EDGE_VTX )
  {
    dn_entity = dmesh->dn_edge;
  } else if( connectivity_type == PDM_CONNECTIVITY_TYPE_VTX_ELMT ||
             connectivity_type == PDM_CONNECTIVITY_TYPE_VTX_CELL ||
             connectivity_type == PDM_CONNECTIVITY_TYPE_VTX_FACE ||
             connectivity_type == PDM_CONNECTIVITY_TYPE_VTX_EDGE ||
             connectivity_type == PDM_CONNECTIVITY_TYPE_VTX_VTX )
  {
    dn_entity = dmesh->dn_vtx;
  } else if( connectivity_type == PDM_CONNECTIVITY_TYPE_ELMT_CELL ||
             connectivity_type == PDM_CONNECTIVITY_TYPE_ELMT_FACE ||
             connectivity_type == PDM_CONNECTIVITY_TYPE_ELMT_EDGE ||
             connectivity_type == PDM_CONNECTIVITY_TYPE_ELMT_VTX )
  {
    dn_entity = -1;
  }

  return dn_entity;
}


int
PDM_dmesh_bound_get
(
 PDM_dmesh_t       *dmesh,
 PDM_bound_type_t   bound_type,
 PDM_g_num_t      **connect,
 int              **connect_idx,
 PDM_ownership_t    ownership
)
{
  PDM_UNUSED(ownership);
  assert(dmesh != NULL);

  // assert(dmesh->dbound[bound_type] != NULL);

  *connect     = dmesh->dbound    [bound_type];
  *connect_idx = dmesh->dbound_idx[bound_type];

  return dmesh->n_bnd;
}


int
PDM_dmesh_distrib_get
(
 PDM_dmesh_t              *dmesh,
 PDM_mesh_entities_t       entity_type,
 PDM_g_num_t             **distrib
)
{
  switch (entity_type) {
   case PDM_MESH_ENTITY_CELL:
     *distrib = dmesh->cell_distrib;
     break;
   case PDM_MESH_ENTITY_FACE:
     *distrib = dmesh->face_distrib;
     break;
   case PDM_MESH_ENTITY_EDGE:
     *distrib = dmesh->edge_distrib;
     break;
   case PDM_MESH_ENTITY_VERTEX:
     *distrib = dmesh->vtx_distrib;
     break;
   }
   int n_rank;
   PDM_MPI_Comm_size(dmesh->comm, &n_rank);
   return n_rank;
}


/**
 *
 * \brief Free
 *
 * \param [in]   id           Identifier
 *
 */

void
PDM_dmesh_free
(
 PDM_dmesh_t         *dmesh
)
{
  dmesh->dn_cell           = 0;
  dmesh->dn_face           = 0;
  dmesh->dn_edge           = 0;
  dmesh->dn_vtx            = 0;
  dmesh->n_bnd             = 0;
  dmesh->n_join            = 0;

  dmesh->_dface_cell       = NULL;
  dmesh->_dface_vtx_idx    = NULL;
  dmesh->_dface_vtx        = NULL;
  dmesh->_dvtx_coord       = NULL;
  dmesh->_dface_bound_idx  = NULL;
  dmesh->_dface_bound      = NULL;

  dmesh->_dedge_bound_idx  = NULL;
  dmesh->_dedge_bound      = NULL;

  dmesh->_joins_glob_id    = NULL;
  dmesh->_dface_join_idx   = NULL;
  dmesh->_dface_join       = NULL;

  // On doit gérer les cas ou la structure est partagé en python et auquel cas
  // On est owner des resultats et il faut free le reste
  // Donc il faut un is_getted + is_owner pour s'en sortir

  if(( dmesh->owner == PDM_OWNERSHIP_KEEP ) ||
     ( dmesh->owner == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE)){
    for(int i = 0; i < PDM_CONNECTIVITY_TYPE_MAX; ++i) {

      if(dmesh->is_owner_connectivity[i] == PDM_TRUE) {

        if(dmesh->dconnectivity[i] != NULL){
          free(dmesh->dconnectivity[i]);
        }
        if(dmesh->dconnectivity_idx[i] != NULL){
          free(dmesh->dconnectivity_idx[i]);
        }
        dmesh->dconnectivity    [i] = NULL;
        dmesh->dconnectivity_idx[i] = NULL;

      }
    }

    for(int i = 0; i < PDM_BOUND_TYPE_MAX; ++i) {

      if(dmesh->is_owner_bound[i] == PDM_TRUE) {

        printf(" dmesh_free :: %i \n", i);
        assert(dmesh->dbound[i] != NULL);
        free(dmesh->dbound[i]);
        if(dmesh->dbound_idx[i] != NULL){
          free(dmesh->dbound_idx[i]);
        }
        dmesh->dbound    [i] = NULL;
        dmesh->dbound_idx[i] = NULL;

      }
    }
  }

  free(dmesh->results_is_getted    );
  free(dmesh->dconnectivity        );
  free(dmesh->dconnectivity_idx    );
  free(dmesh->is_owner_connectivity);

  free(dmesh->dbound        );
  free(dmesh->dbound_idx    );
  free(dmesh->is_owner_bound);

  /* This result is never getted so we can free them */
  if(dmesh->cell_distrib != NULL) {
    free(dmesh->cell_distrib);
    dmesh->cell_distrib = NULL;
  }

  if(dmesh->face_distrib != NULL) {
    free(dmesh->face_distrib);
    dmesh->face_distrib = NULL;
  }

  if(dmesh->edge_distrib != NULL) {
    free(dmesh->edge_distrib);
    dmesh->edge_distrib = NULL;
  }

  if(dmesh->vtx_distrib != NULL) {
    free(dmesh->vtx_distrib);
    dmesh->vtx_distrib = NULL;
  }

  free (dmesh);
}


#ifdef __cplusplus
}
#endif /* __cplusplus */
