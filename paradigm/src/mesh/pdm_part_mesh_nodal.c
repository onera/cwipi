/*----------------------------------------------------------------------------
 *  System headers
 *----------------------------------------------------------------------------*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_mpi.h"
#include "pdm_part_mesh_nodal.h"
#include "pdm_part_mesh_nodal_priv.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_gnum.h"
#include "pdm_geom_elem.h"
#include "pdm_array.h"

#include "pdm_writer.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Fortran function header
 *============================================================================*/

/*============================================================================
 * Local macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Maximum number of blocks depending of block type
 *----------------------------------------------------------------------------*/

/*============================================================================
 * Global variable
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

static
void
_vtx_free
(
 PDM_Mesh_nodal_vtx_t *vtx
)
{
  if (vtx != NULL) {
    if (vtx->parent != NULL) {
      _vtx_free (vtx->parent);
    }

    if (vtx->coords != NULL) {
      free (vtx->coords);
      vtx->coords = NULL;
    }

    free (vtx);
  }
}

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/**
 * \brief Create a Mesh nodal structure
 *
 * \param [in]   n_part   Number of partition on the current process
 * \param [in]   comm     MPI communicator
 *
 * \return       New mesh nodal handle
 *
 */

PDM_part_mesh_nodal_t*
PDM_part_mesh_nodal_create
(
 const int          mesh_dimension,
 const int          n_part,
 const PDM_MPI_Comm comm
)
{
  PDM_part_mesh_nodal_t *pmn = (PDM_part_mesh_nodal_t *) malloc (sizeof(PDM_part_mesh_nodal_t));

  pmn->comm           = comm;
  pmn->mesh_dimension = mesh_dimension;
  pmn->n_part         = n_part;

  pmn->vtx      = malloc(n_part * sizeof(PDM_Mesh_nodal_vtx_t *));
  for (int i = 0; i < n_part; i++) {
    pmn->vtx[i] = malloc(sizeof(PDM_Mesh_nodal_vtx_t));
    pmn->vtx[i]->_coords    = NULL;
    pmn->vtx[i]->_numabs    = NULL;
    pmn->vtx[i]->_numparent = NULL;
    pmn->vtx[i]->n_vtx      = 0;
    pmn->vtx[i]->parent     = NULL;
    pmn->vtx[i]->coords     = NULL;
    pmn->vtx[i]->owner      = PDM_OWNERSHIP_USER;
  }

  pmn->n_vol    = (int *) malloc( n_part * sizeof(int));
  pmn->n_surf   = (int *) malloc( n_part * sizeof(int));
  pmn->n_ridge  = (int *) malloc( n_part * sizeof(int));
  pmn->n_corner = (int *) malloc( n_part * sizeof(int));

  pmn->volumic  = NULL;
  pmn->surfacic = NULL;
  pmn->ridge    = NULL;
  pmn->corner   = NULL;

  pmn->is_owner_volumic  = PDM_OWNERSHIP_USER;
  pmn->is_owner_surfacic = PDM_OWNERSHIP_USER;
  pmn->is_owner_ridge    = PDM_OWNERSHIP_USER;
  pmn->is_owner_corner   = PDM_OWNERSHIP_USER;

  return pmn;
}


void
PDM_part_mesh_nodal_coord_set
(
       PDM_part_mesh_nodal_t *pmn,
 const int                    id_part,
 const int                    n_vtx,
 const PDM_real_t            *coords,
 const PDM_g_num_t           *numabs,
       PDM_ownership_t        owner
)
{

  if (id_part >= pmn->n_part) {
    PDM_error (__FILE__, __LINE__, 0, "Bad part identifier\n");
  }

  PDM_Mesh_nodal_vtx_t *vtx = pmn->vtx[id_part];

  if ((vtx->_coords != NULL) ||
      (vtx->_numabs != NULL)) {
    PDM_error(__FILE__, __LINE__, 0, "these partition vertices are already defined\n");
  }

  /* Mapping memoire */
  vtx->n_vtx   = n_vtx;
  vtx->_coords = coords;
  vtx->_numabs = numabs;
  vtx->owner   = owner;

}

void
PDM_part_mesh_nodal_add_part_mesh_nodal_elmts
(
 PDM_part_mesh_nodal_t       *pmn,
 PDM_part_mesh_nodal_elmts_t *pmne,
 PDM_ownership_t              owner
)
{
  assert(pmn->n_part == pmne->n_part);
  assert(pmn->mesh_dimension >= pmne->mesh_dimension);
  if(pmne->mesh_dimension == 3) {
    pmn->volumic          = pmne;
    pmn->is_owner_volumic = owner;
  } else if(pmne->mesh_dimension == 2){
    pmn->surfacic          = pmne;
    pmn->is_owner_surfacic = owner;
  } else if(pmne->mesh_dimension == 1){
    pmn->ridge          = pmne;
    pmn->is_owner_ridge = owner;
  } else if(pmne->mesh_dimension == 0){
    pmn->corner          = pmne;
    pmn->is_owner_corner = owner;
  } else {
    PDM_error (__FILE__, __LINE__, 0, "PDM_Mesh_nodal_add_dmesh_nodal_elmts bad mesh_dimension\n");
  }
}

void
PDM_part_mesh_nodal_free
(
 PDM_part_mesh_nodal_t* pmn
)
{

  if(pmn->is_owner_volumic == PDM_OWNERSHIP_KEEP) {
    PDM_part_mesh_nodal_elmts_free(pmn->volumic);
  }
  if(pmn->is_owner_surfacic == PDM_OWNERSHIP_KEEP) {
    PDM_part_mesh_nodal_elmts_free(pmn->surfacic);
  }
  if(pmn->is_owner_ridge == PDM_OWNERSHIP_KEEP) {
    PDM_part_mesh_nodal_elmts_free(pmn->ridge);
  }
  if(pmn->is_owner_corner == PDM_OWNERSHIP_KEEP) {
    PDM_part_mesh_nodal_elmts_free(pmn->corner);
  }

  if (pmn->vtx != NULL) {
    for (int i = 0; i < pmn->n_part; i++) {
      if(pmn->vtx[i]->owner == PDM_OWNERSHIP_KEEP){
        _vtx_free (pmn->vtx[i]);
      }
    }

    free(pmn->vtx);
    pmn->vtx = NULL;
  }

  free(pmn->n_vol   );
  free(pmn->n_surf  );
  free(pmn->n_ridge );
  free(pmn->n_corner);

  free(pmn);
}
