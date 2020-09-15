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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"
#include "pdm_handles.h"
#include "pdm_printf.h"
#include "pdm_error.h"
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


/**
 * \struct _pdm_dmesh_t
 * \brief  Define a distributed mesh. Arrays are shared: this structure does not
 *         holds the data.
 *
 */

typedef struct
{
  int               dn_cell;          /*!< Number of distributed cells         */
  int               dn_face;          /*!< Number of distributed faces         */
  int               dn_vtx;           /*!< Number of distributed vertices      */
  int               n_bnd;            /*!< Number of boundaries                */
  int               n_join;           /*!< Number of interfaces with other zone*/
  const PDM_g_num_t *_dface_cell;     /*!< Face-cell connectivity of distributed
                                        faces (size = 2 * dn_face)
                                        if iface is a boundary face,
                                        _dface_cell[2*iface + 1] = 0           */
  const int         *_dface_vtx_idx;   /*!< Face-vertex connectivity index of
                                        distributed faces (size = dn_face + 1) */
  const PDM_g_num_t *_dface_vtx;      /*!< Face-vertex connectivity of
                                        distributed faces (size = dface_vtx_idx[
                                        dn_face])                              */
  const double      *_dvtx_coord;     /*!< Coordinates of ditributed vertices
                                        (size = 3 * dn_vtx)                    */
  const int         *_dface_bound_idx; /*!< Index of distributed faces list of
                                        each boundary (size = n_bnd + 1)
                                        or NULL                               */
  const PDM_g_num_t *_dface_bound;    /*!< Distributed faces list of each
                                       boundary (size = dface_bound_idx[n_bnd])
                                        or NULL                               */
  const int         *_joins_glob_id;  /*!< Global id of each joi (size=n_join)
                                           or NULL. Same data for all procs   */
  const int         *_dface_join_idx;  /*!< Index of distributed faces list of
                                        each join (size = n_join + 1)
                                        or NULL                               */
  const PDM_g_num_t *_dface_join;     /*!< Distributed faces list of each
                                       join (size = dface_join_idx[n_join])
                                        or NULL                               */
} _pdm_dmesh_t;

/*============================================================================
 * Global variable
 *============================================================================*/

static PDM_Handles_t *_dmeshes   = NULL;

/*=============================================================================
 * Private function definitions
 *============================================================================*/

/**
 *
 * \brief Return ppart object from it identifier
 *
 * \param [in]   ppart_id        ppart identifier
 *
 */

static _pdm_dmesh_t *
_get_from_id
(
 int  id
)
{

  _pdm_dmesh_t *dmesh = (_pdm_dmesh_t *) PDM_Handles_get (_dmeshes, id);

  if (dmesh == NULL) {
    PDM_error(__FILE__, __LINE__, 0, "PDM_dmesh error : Bad identifier\n");
  }

  return dmesh;
}


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

int
PDM_dmesh_create
(
 const int          dn_cell,
 const int          dn_face,
 const int          dn_vtx,
 const int          n_bnd,
 const int          n_join
)
{

  /*
   * Search a ppart free id
   */
  if (_dmeshes == NULL) {
    _dmeshes = PDM_Handles_create (4);
  }

  _pdm_dmesh_t *dmesh = (_pdm_dmesh_t *) malloc(sizeof(_pdm_dmesh_t));
  int id = PDM_Handles_store (_dmeshes, dmesh);

  dmesh->dn_cell          = dn_cell;
  dmesh->dn_face          = dn_face;
  dmesh->dn_vtx           = dn_vtx;
  dmesh->n_bnd            = n_bnd;
  dmesh->n_join           = n_join;
  dmesh->_dface_cell      = NULL;
  dmesh->_dface_vtx_idx   = NULL;
  dmesh->_dface_vtx       = NULL;
  dmesh->_dvtx_coord      = NULL;
  dmesh->_dface_bound_idx = NULL;
  dmesh->_dface_bound     = NULL;
  dmesh->_joins_glob_id   = NULL;
  dmesh->_dface_join_idx  = NULL;
  dmesh->_dface_join      = NULL;

  return id;

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
 const int           id,
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
  _pdm_dmesh_t *dmesh = _get_from_id (id);

  dmesh->_dvtx_coord      = dvtx_coord;
  dmesh->_dface_vtx_idx   = dface_vtx_idx;
  dmesh->_dface_vtx       = dface_vtx;
  dmesh->_dface_cell      = dface_cell;
  dmesh->_dface_bound_idx = dface_bound_idx;
  dmesh->_dface_bound     = dface_bound;
  dmesh->_joins_glob_id   = joins_glob_id;
  dmesh->_dface_join_idx  = dface_join_idx;
  dmesh->_dface_join      = dface_join;

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
 const int   id,
 int        *dn_cell,
 int        *dn_face,
 int        *dn_vtx,
 int        *n_bnd,
 int        *n_join
)
{
  _pdm_dmesh_t *dmesh = _get_from_id (id);
  *dn_cell = dmesh->dn_cell;
  *dn_face = dmesh->dn_face;
  *dn_vtx  = dmesh->dn_vtx;
  *n_bnd  = dmesh->n_bnd;
  *n_join = dmesh->n_join;
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
 const int          id,
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
  _pdm_dmesh_t *dmesh = _get_from_id (id);

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
 const int id
)
{
  _pdm_dmesh_t *dmesh = _get_from_id (id);

  dmesh->dn_cell           = 0;
  dmesh->dn_face           = 0;
  dmesh->dn_vtx            = 0;
  dmesh->n_bnd             = 0;
  dmesh->n_join            = 0;
  dmesh->_dface_cell       = NULL;
  dmesh->_dface_vtx_idx    = NULL;
  dmesh->_dface_vtx        = NULL;
  dmesh->_dvtx_coord       = NULL;
  dmesh->_dface_bound_idx  = NULL;
  dmesh->_dface_bound      = NULL;
  dmesh->_joins_glob_id    = NULL;
  dmesh->_dface_join_idx   = NULL;
  dmesh->_dface_join       = NULL;

  free (dmesh);

  PDM_Handles_handle_free (_dmeshes, id, PDM_FALSE);

  const int n_dmesh = PDM_Handles_n_get (_dmeshes);

  if (n_dmesh == 0) {
    _dmeshes = PDM_Handles_free (_dmeshes);
  }

}


/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
