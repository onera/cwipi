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

static
PDM_part_mesh_nodal_elmts_t*
_get_from_geometry_kind
(
 PDM_part_mesh_nodal_t *pmn,
 PDM_geometry_kind_t    geom_kind
)
{
  PDM_part_mesh_nodal_elmts_t* pmne = NULL;
  if(geom_kind == PDM_GEOMETRY_KIND_VOLUMIC){
    assert(pmn->mesh_dimension == 3);
    pmne = pmn->volumic;
  } else if( geom_kind == PDM_GEOMETRY_KIND_SURFACIC){
    assert(pmn->mesh_dimension >= 2);
    pmne = pmn->surfacic;
  } else if( geom_kind == PDM_GEOMETRY_KIND_RIDGE){
    assert(pmn->mesh_dimension >= 1);
    pmne = pmn->ridge;
  } else if( geom_kind == PDM_GEOMETRY_KIND_CORNER){
    pmne = pmn->corner;
  } else {
    PDM_error(__FILE__, __LINE__, 0, "Bad geom_kind in _get_from_geometry_kind \n");
  }
  return pmne;
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


int
PDM_part_mesh_nodal_n_part_get
(
       PDM_part_mesh_nodal_t *pmn
)
{
  return pmn->n_part;
}


int
PDM_part_mesh_nodal_n_vtx_get
(
       PDM_part_mesh_nodal_t *pmn,
 const int                    id_part
)
{
  if (id_part >= pmn->n_part) {
    PDM_error (__FILE__, __LINE__, 0, "Bad part identifier\n");
  }

  PDM_Mesh_nodal_vtx_t *vtx = pmn->vtx[id_part];
  return vtx->n_vtx;
}

double*
PDM_part_mesh_nodal_vtx_coord_get
(
       PDM_part_mesh_nodal_t *pmn,
 const int                    id_part
)
{
  if (id_part >= pmn->n_part) {
    PDM_error (__FILE__, __LINE__, 0, "Bad part identifier\n");
  }

  PDM_Mesh_nodal_vtx_t *vtx = pmn->vtx[id_part];
  return (double *) vtx->_coords;
}

PDM_g_num_t*
PDM_part_mesh_nodal_vtx_g_num_get
(
       PDM_part_mesh_nodal_t *pmn,
 const int                    id_part
)
{
  if (id_part >= pmn->n_part) {
    PDM_error (__FILE__, __LINE__, 0, "Bad part identifier\n");
  }

  PDM_Mesh_nodal_vtx_t *vtx = pmn->vtx[id_part];
  return (PDM_g_num_t*) vtx->_numabs;
}


int
PDM_part_mesh_nodal_n_section_get
(
 PDM_part_mesh_nodal_t *pmn,
 PDM_geometry_kind_t    geom_kind
)
{
  if (pmn == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad part mesh nodal identifier\n");
  }
  PDM_part_mesh_nodal_elmts_t* pmne = _get_from_geometry_kind(pmn, geom_kind);
  if(pmne){
    return pmne->n_section;
  } else {
    return 0;
  }
}


int *
PDM_part_mesh_nodal_sections_id_get
(
 PDM_part_mesh_nodal_t *pmn,
 PDM_geometry_kind_t    geom_kind
)
{
  if (pmn == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad part mesh nodal identifier\n");
  }
  PDM_part_mesh_nodal_elmts_t* pmne = _get_from_geometry_kind(pmn, geom_kind);
  return pmne->sections_id;
}

PDM_Mesh_nodal_elt_t
PDM_part_mesh_nodal_section_elt_type_get
(
        PDM_part_mesh_nodal_t *pmn,
        PDM_geometry_kind_t    geom_kind,
  const int                    id_section
)
{
  PDM_part_mesh_nodal_elmts_t* dmne = _get_from_geometry_kind(pmn, geom_kind);
  if (id_section <= PDM_BLOCK_ID_BLOCK_POLY2D) { // std
    int _id_section = id_section - PDM_BLOCK_ID_BLOCK_STD;
    return dmne->sections_std[_id_section]->t_elt;
  }
  assert(0); // only useful for std elements
}



PDM_Mesh_nodal_elt_t
PDM_part_mesh_nodal_section_type_get
(
      PDM_part_mesh_nodal_t *pmn,
      PDM_geometry_kind_t    geom_kind,
const int                    id_section
)
{
  if (pmn == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }
  PDM_part_mesh_nodal_elmts_t* pmne = _get_from_geometry_kind(pmn, geom_kind);

  PDM_Mesh_nodal_elt_t t_elt = PDM_MESH_NODAL_POLY_3D;

  if (id_section < PDM_BLOCK_ID_BLOCK_POLY2D) {

    t_elt = PDM_MESH_NODAL_POLY_3D;
    int _id_section = id_section - PDM_BLOCK_ID_BLOCK_STD;
    PDM_Mesh_nodal_block_std_t *section = pmne->sections_std[_id_section];

    if (section == NULL) {
      PDM_error (__FILE__, __LINE__, 0, "Bad section identifier\n");
    }

    t_elt = section->t_elt;
  }

  else if (id_section < PDM_BLOCK_ID_BLOCK_POLY3D) {

    t_elt = PDM_MESH_NODAL_POLY_2D;

  }

  else {

    t_elt = PDM_MESH_NODAL_POLY_3D;

  }

  return t_elt;
}


int
PDM_part_mesh_nodal_section_add
(
      PDM_part_mesh_nodal_t *pmn,
      PDM_geometry_kind_t    geom_kind,
const PDM_Mesh_nodal_elt_t   t_elt
)
{
  if( _get_from_geometry_kind(pmn, geom_kind) == NULL) {
    if(geom_kind == PDM_GEOMETRY_KIND_VOLUMIC) {
      pmn->volumic = PDM_part_mesh_nodal_elmts_create(pmn->mesh_dimension, pmn->n_part, pmn->comm);
    } else if( geom_kind == PDM_GEOMETRY_KIND_SURFACIC) {
      pmn->surfacic = PDM_part_mesh_nodal_elmts_create(pmn->mesh_dimension, pmn->n_part, pmn->comm);
    } else if( geom_kind == PDM_GEOMETRY_KIND_RIDGE) {
      pmn->ridge = PDM_part_mesh_nodal_elmts_create(pmn->mesh_dimension, pmn->n_part, pmn->comm);
    } else if( geom_kind == PDM_GEOMETRY_KIND_CORNER) {
      pmn->corner = PDM_part_mesh_nodal_elmts_create(pmn->mesh_dimension, pmn->n_part, pmn->comm);
    }
  }

  PDM_part_mesh_nodal_elmts_t* pmne = _get_from_geometry_kind(pmn, geom_kind);
  assert(pmne != NULL);
  return PDM_part_mesh_nodal_elmts_add(pmne, t_elt);
}


void
PDM_part_mesh_nodal_section_std_set
(
      PDM_part_mesh_nodal_t *pmn,
      PDM_geometry_kind_t    geom_kind,
const int                    id_block,
const int                    id_part,
const int                    n_elt,
const int                   *connec,
const PDM_g_num_t           *numabs,
const int                   *parent_num,
      PDM_ownership_t        owner
)
{
  PDM_part_mesh_nodal_elmts_t* pmne = _get_from_geometry_kind(pmn, geom_kind);
  assert(pmne != NULL);
  PDM_part_mesh_nodal_elmts_std_set(pmne, id_block, id_part, n_elt, connec, numabs, parent_num, owner);
}

int
PDM_part_mesh_nodal_block_n_elt_get
(
      PDM_part_mesh_nodal_t  *pmn,
      PDM_geometry_kind_t     geom_kind,
const int                     id_block,
const int                     id_part
)
{
  PDM_part_mesh_nodal_elmts_t* pmne = _get_from_geometry_kind(pmn, geom_kind);
  assert(pmne != NULL);
  return PDM_part_mesh_nodal_elmts_block_n_elt_get(pmne, id_block, id_part);
}



void
PDM_part_mesh_nodal_block_std_get
(
      PDM_part_mesh_nodal_t  *pmn,
      PDM_geometry_kind_t     geom_kind,
const int                     id_block,
const int                     id_part,
      int                   **connec,
      PDM_g_num_t           **numabs,
      int                   **parent_num
)
{
  PDM_part_mesh_nodal_elmts_t* pmne = _get_from_geometry_kind(pmn, geom_kind);
  assert(pmne != NULL);
  PDM_part_mesh_nodal_elmts_block_std_get(pmne, id_block, id_part, connec, numabs, parent_num);
}


PDM_Mesh_nodal_elt_t
PDM_part_mesh_nodal_block_type_get
(
      PDM_part_mesh_nodal_t  *pmn,
      PDM_geometry_kind_t     geom_kind,
const int                     id_block
)
{
  PDM_part_mesh_nodal_elmts_t* pmne = _get_from_geometry_kind(pmn, geom_kind);
  assert(pmne != NULL);
  PDM_part_mesh_nodal_elmts_block_type_get(pmne, id_block);
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



// Faire extraction de maillage depuis mesh_nodal
// Donc par exemple rendre un nouveau maillage (avec numero absolu des rigdes à partir du mesh_nodal )
// PDM_extract_from_indices()
