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
_block_std_free_partial
(
 PDM_Mesh_nodal_block_std_t *_block_std
)
{

  if (_block_std == NULL) {
    return;
  }

  if (_block_std->_connec != NULL) {
    if (_block_std->owner == PDM_OWNERSHIP_KEEP) {
      for (int i = 0; i < _block_std->n_part; i++) {
        if (_block_std->_connec[i] != NULL)
          free(_block_std->_connec[i]);
        _block_std->_connec[i] = NULL;
      }
    }
    free(_block_std->_connec);
    _block_std->_connec = NULL;
  }

  if (_block_std->_numabs != NULL) {

    if (_block_std->owner == PDM_OWNERSHIP_KEEP) {
      for (int i = 0; i < _block_std->n_part; i++) {
        if (_block_std->_numabs[i] != NULL)
          free(_block_std->_numabs[i]);
        _block_std->_numabs[i] = NULL;
      }
    }
    free(_block_std->_numabs);
    _block_std->_numabs = NULL;
  }

}


static
PDM_Mesh_nodal_block_std_t *
_block_std_free
(
 PDM_Mesh_nodal_block_std_t *_block_std
)
{

  if (_block_std == NULL) {
    return NULL;
  }

  _block_std_free_partial(_block_std);

  if (_block_std->n_elt != NULL) {
    free(_block_std->n_elt);
    _block_std->n_elt = NULL;
  }

  if (_block_std->numabs_int != NULL) {
    for (int j = 0; j < _block_std->n_part; j++) {
      if (_block_std->numabs_int[j] != NULL) {
        free(_block_std->numabs_int[j]);
      }
    }
    free(_block_std->numabs_int);
    _block_std->numabs_int = NULL;
  }

  if (_block_std->cell_centers != NULL) {
    for (int j = 0; j < _block_std->n_part; j++) {
      if (_block_std->cell_centers[j] != NULL) {
        free(_block_std->cell_centers[j]);
      }
    }
    free(_block_std->cell_centers);
    _block_std->cell_centers = NULL;
  }

  if (_block_std->_parent_num != NULL) {
    if (_block_std->owner == PDM_OWNERSHIP_KEEP) {
      for (int i = 0; i < _block_std->n_part; i++) {
        if (_block_std->_parent_num[i] != NULL)
          free(_block_std->_parent_num[i]);
        _block_std->_parent_num[i] = NULL;
      }
    }
    free(_block_std->_parent_num);
    _block_std->_parent_num = NULL;
  }

  if (_block_std->_parent_entity_g_num != NULL) {
    if (_block_std->owner == PDM_OWNERSHIP_KEEP) {
      for (int i = 0; i < _block_std->n_part; i++) {
        if (_block_std->_parent_entity_g_num[i] != NULL)
          free(_block_std->_parent_entity_g_num[i]);
        _block_std->_parent_entity_g_num[i] = NULL;
      }
    }
    free(_block_std->_parent_entity_g_num);
    _block_std->_parent_entity_g_num = NULL;
  }

  free(_block_std);
  return NULL;
}

/**
 *
 * \brief Update blocks identifier list
 *
 * \param [inout]  mesh        Mesh
 */

static void
_update_elmt_sections_id
(
 PDM_part_mesh_nodal_elmts_t *pmne
)
{
  int n_section = 0;

  if (pmne->sections_std != NULL) {
    n_section += pmne->n_section_std;
  }

  if (pmne->sections_poly2d != NULL) {
    n_section += pmne->n_section_poly2d;
  }

  if (pmne->sections_poly3d != NULL) {
    n_section += pmne->n_section_poly3d;
  }

  if (pmne->n_section < n_section) {
    pmne->sections_id = (int *) realloc(pmne->sections_id, sizeof(int) * n_section);
  }

  int k = 0;
  if (pmne->sections_std != NULL) {
    for (int i = 0; i < pmne->n_section_std; i++) {
      pmne->sections_id[k++] = i + PDM_BLOCK_ID_BLOCK_STD;
    }
  }

  if (pmne->sections_poly2d != NULL) {
    for (int i = 0; i < pmne->n_section_poly2d; i++) {
      pmne->sections_id[k++] = i + PDM_BLOCK_ID_BLOCK_POLY2D;
    }
  }

  if (pmne->sections_poly3d != NULL) {
    for (int i = 0; i < pmne->n_section_poly3d; i++) {
      pmne->sections_id[k++] = i + PDM_BLOCK_ID_BLOCK_POLY3D;
    }
  }

  pmne->n_section = n_section;
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

PDM_part_mesh_nodal_elmts_t*
PDM_part_mesh_nodal_elmts_create
(
 const int          mesh_dimension,
 const int          n_part,
 const PDM_MPI_Comm comm
)
{
  PDM_part_mesh_nodal_elmts_t *pmne = (PDM_part_mesh_nodal_elmts_t *) malloc (sizeof(PDM_part_mesh_nodal_elmts_t));

  pmne->comm             = comm;
  pmne->mesh_dimension   = mesh_dimension;
  pmne->n_part           = n_part;

  pmne->n_elmts          = PDM_array_zeros_int(n_part);

  pmne->n_section        = 0;
  pmne->n_section_std    = 0;
  pmne->n_section_poly2d = 0;
  pmne->n_section_poly3d = 0;

  pmne->sections_id      = NULL;
  pmne->sections_std     = NULL;
  pmne->sections_poly2d  = NULL;
  pmne->sections_poly3d  = NULL;

  return pmne;
}


int
PDM_part_mesh_nodal_elmts_add
(
      PDM_part_mesh_nodal_elmts_t *pmne,
const PDM_Mesh_nodal_elt_t         t_elt
)
{
  if (pmne == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad pmne nodal identifier\n");
  }

  if(t_elt == PDM_MESH_NODAL_POINT) {
    if(pmne->mesh_dimension != 0){
      PDM_error (__FILE__, __LINE__, 0, "Bad mesh_dimension in PDM_DMesh_nodal_elmts_section_add = expected = %i and given = %i \n", pmne->mesh_dimension, 0);
    }
  } else if(t_elt == PDM_MESH_NODAL_BAR2) {
    if(pmne->mesh_dimension != 1){
      PDM_error (__FILE__, __LINE__, 0, "Bad mesh_dimension in PDM_DMesh_nodal_elmts_section_add = expected = %i and given = %i \n", pmne->mesh_dimension, 1);
    }
  } else if(t_elt == PDM_MESH_NODAL_TRIA3 || t_elt == PDM_MESH_NODAL_QUAD4 || t_elt == PDM_MESH_NODAL_POLY_2D) {
    if(pmne->mesh_dimension != 2){
      PDM_error (__FILE__, __LINE__, 0, "Bad mesh_dimension in PDM_DMesh_nodal_elmts_section_add = expected = %i and given = %i \n", pmne->mesh_dimension, 2);
    }
  } else {
    if(pmne->mesh_dimension != 3){
      PDM_error (__FILE__, __LINE__, 0, "Bad mesh_dimension in PDM_DMesh_nodal_elmts_section_add = expected = %i and given = %i \n", pmne->mesh_dimension, 3);
    }
  }

  int id_block = -1;

  switch (t_elt) {

  case PDM_MESH_NODAL_POINT    :
  case PDM_MESH_NODAL_BAR2     :
  case PDM_MESH_NODAL_TRIA3    :
  case PDM_MESH_NODAL_QUAD4    :
  case PDM_MESH_NODAL_TETRA4   :
  case PDM_MESH_NODAL_PYRAMID5 :
  case PDM_MESH_NODAL_PRISM6   :
  case PDM_MESH_NODAL_HEXA8    :
    {
      /* Mise a jour du tableau de stockage */

      pmne->n_section_std++;

      pmne->sections_std = realloc(pmne->sections_std, pmne->n_section_std * sizeof(PDM_Mesh_nodal_block_std_t *));

      id_block = pmne->n_section_std-1;

      /* Intialisation du bloc */
      pmne->sections_std[id_block] = malloc( sizeof(PDM_Mesh_nodal_block_std_t) );
      pmne->sections_std[id_block]->t_elt        = t_elt;
      pmne->sections_std[id_block]->n_part       = pmne->n_part;

      pmne->sections_std[id_block]->n_elt                 = (int  *) malloc(sizeof(int  ) * pmne->sections_std[id_block]->n_part);
      pmne->sections_std[id_block]->_connec               = (int **) malloc(sizeof(int *) * pmne->sections_std[id_block]->n_part);
      pmne->sections_std[id_block]->_numabs               = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t *) * pmne->sections_std[id_block]->n_part);
      pmne->sections_std[id_block]->numabs_int            = NULL;
      pmne->sections_std[id_block]->_parent_num           = NULL;
      pmne->sections_std[id_block]->_parent_entity_g_num  = NULL;
      pmne->sections_std[id_block]->cell_centers          = NULL;
      pmne->sections_std[id_block]->owner                 = PDM_OWNERSHIP_KEEP;
      pmne->sections_std[id_block]->order                 = 1;
      pmne->sections_std[id_block]->ho_ordering           = NULL;

      for (int i = 0; i < pmne->sections_std[id_block]->n_part; i++) {
        pmne->sections_std[id_block]->n_elt    [i] = 0;
        pmne->sections_std[id_block]->_connec  [i] = NULL;
        pmne->sections_std[id_block]->_numabs  [i] = NULL;
      }

      id_block += PDM_BLOCK_ID_BLOCK_STD;
      if (id_block >= PDM_BLOCK_ID_BLOCK_POLY2D) {
        PDM_error(__FILE__, __LINE__, 0, "The number of standard blocks must be less than %d\n",
                  PDM_BLOCK_ID_BLOCK_POLY2D);
        abort();
      }
    }

    break;

  case PDM_MESH_NODAL_POLY_2D  :
    {
      /* Mise a jour du tableau de stockage */

      pmne->n_section_poly2d++;

      pmne->sections_poly2d = realloc(pmne->sections_poly2d, pmne->n_section_poly2d * sizeof(PDM_Mesh_nodal_block_poly2d_t *));

      id_block = pmne->n_section_poly2d-1;

      /* Intialisation du bloc */
      pmne->sections_poly2d[id_block] = malloc( sizeof(PDM_Mesh_nodal_block_poly2d_t) );
      pmne->sections_poly2d[id_block]->n_part            = pmne->n_part;

      pmne->sections_poly2d[id_block]->n_elt                 = (int * ) malloc(sizeof(int  ) * pmne->sections_poly2d[id_block]->n_part);
      pmne->sections_poly2d[id_block]->_connec_idx           = (int **) malloc(sizeof(int *) * pmne->sections_poly2d[id_block]->n_part);
      pmne->sections_poly2d[id_block]->_connec               = (int **) malloc(sizeof(int *) * pmne->sections_poly2d[id_block]->n_part);
      pmne->sections_poly2d[id_block]->_numabs               = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t *) * pmne->sections_poly2d[id_block]->n_part);
      pmne->sections_poly2d[id_block]->numabs_int            = NULL;
      pmne->sections_poly2d[id_block]->cell_centers          = NULL;
      pmne->sections_poly2d[id_block]->_parent_num           = NULL;
      pmne->sections_poly2d[id_block]->_parent_entity_g_num  = NULL;
      pmne->sections_poly2d[id_block]->owner                 = PDM_OWNERSHIP_KEEP;

      for (int i = 0; i < pmne->sections_poly2d[id_block]->n_part; i++) {
        pmne->sections_poly2d[id_block]->n_elt      [i] = 0;
        pmne->sections_poly2d[id_block]->_connec_idx[i] = NULL;
        pmne->sections_poly2d[id_block]->_connec    [i] = NULL;
        pmne->sections_poly2d[id_block]->_numabs    [i] = NULL;
      }

      id_block += PDM_BLOCK_ID_BLOCK_POLY2D;
      if (id_block >= PDM_BLOCK_ID_BLOCK_POLY3D) {
        PDM_error(__FILE__, __LINE__, 0, "The number of polygon blocks must be less than %d\n",
                  PDM_BLOCK_ID_BLOCK_POLY3D - PDM_BLOCK_ID_BLOCK_POLY2D);
      }
    }

    break;

  case PDM_MESH_NODAL_POLY_3D  :
    {
      pmne->n_section_poly3d++;

      pmne->sections_poly3d = realloc(pmne->sections_poly3d, pmne->n_section_poly3d * sizeof(PDM_Mesh_nodal_block_poly3d_t *));

      id_block = pmne->n_section_poly3d-1;

      /* Intialisation du bloc */

      pmne->sections_poly3d[id_block] = malloc( sizeof(PDM_Mesh_nodal_block_poly3d_t) );
      pmne->sections_poly3d[id_block]->n_part       = pmne->n_part;

      pmne->sections_poly3d[id_block]->n_elt                = (int * ) malloc(sizeof(int  ) * pmne->sections_poly3d[id_block]->n_part);
      pmne->sections_poly3d[id_block]->n_face               = (int * ) malloc(sizeof(int  ) * pmne->sections_poly3d[id_block]->n_part);
      pmne->sections_poly3d[id_block]->_facvtx_idx          = (int **) malloc(sizeof(int *) * pmne->sections_poly3d[id_block]->n_part);
      pmne->sections_poly3d[id_block]->_facvtx              = (int **) malloc(sizeof(int *) * pmne->sections_poly3d[id_block]->n_part);
      pmne->sections_poly3d[id_block]->_cellfac_idx         = (int **) malloc(sizeof(int *) * pmne->sections_poly3d[id_block]->n_part);
      pmne->sections_poly3d[id_block]->_cellfac             = (int **) malloc(sizeof(int *) * pmne->sections_poly3d[id_block]->n_part);
      pmne->sections_poly3d[id_block]->_cellvtx_idx         = (int **) malloc(sizeof(int *) * pmne->sections_poly3d[id_block]->n_part);
      pmne->sections_poly3d[id_block]->_cellvtx             = (int **) malloc(sizeof(int *) * pmne->sections_poly3d[id_block]->n_part);
      pmne->sections_poly3d[id_block]->_numabs              = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t *) * pmne->sections_poly3d[id_block]->n_part);
      pmne->sections_poly3d[id_block]->numabs_int           = NULL;
      pmne->sections_poly3d[id_block]->cell_centers         = NULL;
      pmne->sections_poly3d[id_block]->_parent_num          = NULL;
      pmne->sections_poly3d[id_block]->_parent_entity_g_num = NULL;

      pmne->sections_poly3d[id_block]->owner        = PDM_OWNERSHIP_KEEP;

      for (int i = 0; i < pmne->sections_poly3d[id_block]->n_part; i++) {
        pmne->sections_poly3d[id_block]->n_elt       [i] = 0;
        pmne->sections_poly3d[id_block]->n_face      [i] = 0;
        pmne->sections_poly3d[id_block]->_facvtx_idx [i] = NULL;
        pmne->sections_poly3d[id_block]->_facvtx     [i] = NULL;
        pmne->sections_poly3d[id_block]->_cellfac_idx[i] = NULL;
        pmne->sections_poly3d[id_block]->_cellfac    [i] = NULL;
        pmne->sections_poly3d[id_block]->_cellvtx_idx[i] = NULL;
        pmne->sections_poly3d[id_block]->_cellvtx    [i] = NULL;
        pmne->sections_poly3d[id_block]->_numabs     [i] = NULL;
      }

      id_block += PDM_BLOCK_ID_BLOCK_POLY3D;

    }

    break;

  default :
    PDM_error(__FILE__, __LINE__, 0, "Unknown element type\n");
    break;

  }

  _update_elmt_sections_id (pmne);
  return id_block ;

}



int
PDM_part_mesh_nodal_elmts_ho_add
(
      PDM_part_mesh_nodal_elmts_t *pmne,
const PDM_Mesh_nodal_elt_t         t_elt,
const int                          order,
const char                        *ho_ordering
)
{
  if (pmne == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad pmne nodal identifier\n");
  }

  if(t_elt == PDM_MESH_NODAL_POINT) {
    if(pmne->mesh_dimension != 0){
      PDM_error (__FILE__, __LINE__, 0, "Bad mesh_dimension in PDM_DMesh_nodal_elmts_section_add = expected = %i and given = %i \n", pmne->mesh_dimension, 0);
    }
  } else if(t_elt == PDM_MESH_NODAL_BAR2 || t_elt == PDM_MESH_NODAL_BARHO) {
    if(pmne->mesh_dimension != 1){
      PDM_error (__FILE__, __LINE__, 0, "Bad mesh_dimension in PDM_DMesh_nodal_elmts_section_add = expected = %i and given = %i \n", pmne->mesh_dimension, 1);
    }
  } else if(t_elt == PDM_MESH_NODAL_TRIA3 || t_elt == PDM_MESH_NODAL_QUAD4 || t_elt == PDM_MESH_NODAL_POLY_2D ||
            t_elt == PDM_MESH_NODAL_TRIAHO || t_elt == PDM_MESH_NODAL_QUADHO) {
    if(pmne->mesh_dimension != 2){
      PDM_error (__FILE__, __LINE__, 0, "Bad mesh_dimension in PDM_DMesh_nodal_elmts_section_add = expected = %i and given = %i \n", pmne->mesh_dimension, 2);
    }
  } else {
    if(pmne->mesh_dimension != 3){
      PDM_error (__FILE__, __LINE__, 0, "Bad mesh_dimension in PDM_DMesh_nodal_elmts_section_add = expected = %i and given = %i \n", pmne->mesh_dimension, 3);
    }
  }

  int id_block = -1;

  switch (t_elt) {

  case PDM_MESH_NODAL_BAR2      :
  case PDM_MESH_NODAL_BARHO     :
  case PDM_MESH_NODAL_TRIA3     :
  case PDM_MESH_NODAL_TRIAHO    :
  case PDM_MESH_NODAL_QUAD4     :
  case PDM_MESH_NODAL_QUADHO    :
  case PDM_MESH_NODAL_TETRA4    :
  case PDM_MESH_NODAL_TETRAHO   :
  case PDM_MESH_NODAL_PYRAMID5  :
  case PDM_MESH_NODAL_PYRAMIDHO :
  case PDM_MESH_NODAL_PRISM6    :
  case PDM_MESH_NODAL_PRISMHO   :
  case PDM_MESH_NODAL_HEXA8     :
  case PDM_MESH_NODAL_HEXAHO    :
    {
      /* Mise a jour du tableau de stockage */

      pmne->n_section_std++;

      pmne->sections_std = realloc(pmne->sections_std, pmne->n_section_std * sizeof(PDM_Mesh_nodal_block_std_t *));

      id_block = pmne->n_section_std-1;

      /* Intialisation du bloc */
      pmne->sections_std[id_block] = malloc( sizeof(PDM_Mesh_nodal_block_std_t) );
      pmne->sections_std[id_block]->t_elt        = t_elt;
      pmne->sections_std[id_block]->n_part       = pmne->n_part;

      pmne->sections_std[id_block]->n_elt                 = (int  *) malloc(sizeof(int  ) * pmne->sections_std[id_block]->n_part);
      pmne->sections_std[id_block]->_connec               = (int **) malloc(sizeof(int *) * pmne->sections_std[id_block]->n_part);
      pmne->sections_std[id_block]->_numabs               = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t *) * pmne->sections_std[id_block]->n_part);
      pmne->sections_std[id_block]->numabs_int            = NULL;
      pmne->sections_std[id_block]->_parent_num           = NULL;
      pmne->sections_std[id_block]->_parent_entity_g_num  = NULL;
      pmne->sections_std[id_block]->cell_centers          = NULL;
      pmne->sections_std[id_block]->owner                 = PDM_OWNERSHIP_KEEP;
      pmne->sections_std[id_block]->order                 = order;
      pmne->sections_std[id_block]->ho_ordering           = ho_ordering;

      for (int i = 0; i < pmne->sections_std[id_block]->n_part; i++) {
        pmne->sections_std[id_block]->n_elt    [i] = 0;
        pmne->sections_std[id_block]->_connec  [i] = NULL;
        pmne->sections_std[id_block]->_numabs  [i] = NULL;
      }

      id_block += PDM_BLOCK_ID_BLOCK_STD;
      if (id_block >= PDM_BLOCK_ID_BLOCK_POLY2D) {
        PDM_error(__FILE__, __LINE__, 0, "The number of standard blocks must be less than %d\n",
                  PDM_BLOCK_ID_BLOCK_POLY2D);
        abort();
      }
    }

    break;
  default :
    PDM_error(__FILE__, __LINE__, 0, "Unknown element type\n");
    break;

  }

  _update_elmt_sections_id (pmne);
  return id_block ;

}



void
PDM_part_mesh_nodal_elmts_std_set
(
      PDM_part_mesh_nodal_elmts_t *pmne,
const int                          id_block,
const int                          id_part,
const int                          n_elt,
const int                         *connec,
const PDM_g_num_t                 *numabs,
const int                         *parent_num,
const PDM_g_num_t                 *parent_entity_g_num,
      PDM_ownership_t              owner
)
{

  if (pmne == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }

  int _id_block = id_block - PDM_BLOCK_ID_BLOCK_STD;

  PDM_Mesh_nodal_block_std_t *block = pmne->sections_std[_id_block];

  if (block == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad standard block identifier\n");
  }

  if (id_part >= block->n_part) {
    PDM_error(__FILE__, __LINE__, 0, "Partition identifier too big\n");
  }

  /* Mapping */
  block->n_elt  [id_part] += n_elt;
  block->_connec[id_part]  = (PDM_l_num_t *) connec;
  block->_numabs[id_part]  = (PDM_g_num_t *) numabs;
  block->owner             = owner;

  if (parent_num != NULL) {
    if (block->_parent_num == NULL) {
      block->_parent_num = malloc (sizeof(PDM_l_num_t *) * block->n_part);
      for (int i = 0; i < block->n_part; i++) {
        block->_parent_num[i] = NULL;
      }
    }
    block->_parent_num[id_part] = (PDM_l_num_t *) parent_num;
  }

  if (parent_entity_g_num != NULL) {
    if (block->_parent_entity_g_num == NULL) {
      block->_parent_entity_g_num = malloc (sizeof(PDM_l_num_t *) * block->n_part);
      for (int i = 0; i < block->n_part; i++) {
        block->_parent_entity_g_num[i] = NULL;
      }
    }
    block->_parent_entity_g_num[id_part] = (PDM_g_num_t *) parent_entity_g_num;
  }
}

void
PDM_part_mesh_nodal_elmts_block_std_get
(
      PDM_part_mesh_nodal_elmts_t  *pmne,
const int                           id_block,
const int                           id_part,
      int                         **connec,
      PDM_g_num_t                 **numabs,
      int                         **parent_num,
      PDM_g_num_t                 **parent_entity_g_num
)
{
  if (pmne == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad pmne nodal identifier\n");
  }

  int _id_block = id_block - PDM_BLOCK_ID_BLOCK_STD;

  PDM_Mesh_nodal_block_std_t *block = pmne->sections_std[_id_block];

  if (block == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad standard block identifier\n");
  }

  if (id_part >= block->n_part) {
    PDM_error(__FILE__, __LINE__, 0, "Partition identifier too big\n");
  }

  *connec              = block->_connec             [id_part];
  *numabs              = block->_numabs             [id_part];
  *parent_num          = block->_parent_num         [id_part];
  *parent_entity_g_num = NULL;
  if(block->_parent_entity_g_num != NULL) {
    *parent_entity_g_num = block->_parent_entity_g_num[id_part];
  }
}


void
PDM_part_mesh_nodal_elmts_block_std_ho_get
(
      PDM_part_mesh_nodal_elmts_t  *pmne,
const int                           id_block,
const int                           id_part,
      int                         **connec,
      PDM_g_num_t                 **numabs,
      int                         **parent_num,
      PDM_g_num_t                 **parent_entity_g_num,
      int                          *order,
const char                        **ho_ordering
)
{
  if (pmne == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad pmne nodal identifier\n");
  }

  int _id_block = id_block - PDM_BLOCK_ID_BLOCK_STD;

  PDM_Mesh_nodal_block_std_t *block = pmne->sections_std[_id_block];

  if (block == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad standard block identifier\n");
  }

  if (id_part >= block->n_part) {
    PDM_error(__FILE__, __LINE__, 0, "Partition identifier too big\n");
  }

  *connec              = block->_connec             [id_part];
  *numabs              = block->_numabs             [id_part];
  *parent_num          = block->_parent_num         [id_part];
  *parent_entity_g_num = NULL;
  if(block->_parent_entity_g_num != NULL) {
    *parent_entity_g_num = block->_parent_entity_g_num[id_part];
  }
  *order       = block->order;
  *ho_ordering = block->ho_ordering;


}

PDM_Mesh_nodal_elt_t
PDM_part_mesh_nodal_elmts_block_type_get
(
      PDM_part_mesh_nodal_elmts_t *pmne,
const int                          id_block
)
{

  PDM_Mesh_nodal_elt_t t_elt;

  if (pmne == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }

  if (id_block < PDM_BLOCK_ID_BLOCK_POLY2D) {

    t_elt = PDM_MESH_NODAL_POLY_3D;
    const PDM_Mesh_nodal_block_std_t *block = pmne->sections_std[id_block];

    if (block == NULL) {
      PDM_error (__FILE__, __LINE__, 0, "Bad block identifier\n");
    }

    t_elt = block->t_elt;
  }

  else if (id_block < PDM_BLOCK_ID_BLOCK_POLY3D) {

    t_elt = PDM_MESH_NODAL_POLY_2D;

  }

  else {

    t_elt = PDM_MESH_NODAL_POLY_3D;

  }

  return t_elt;

}

int
PDM_part_mesh_nodal_elmts_block_n_elt_get
(
      PDM_part_mesh_nodal_elmts_t *pmne,
const int                          id_block,
const int                          id_part
)
{

  if (pmne == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }

  int _id_block;

  if (id_block >= PDM_BLOCK_ID_BLOCK_POLY3D) {

    _id_block = id_block - PDM_BLOCK_ID_BLOCK_POLY3D;

    PDM_Mesh_nodal_block_poly3d_t *block = pmne->sections_poly3d[_id_block];

    if (block == NULL) {
      PDM_error (__FILE__, __LINE__, 0, "Bad standard block identifier\n");
    }

    if (id_part >= block->n_part) {
      PDM_error(__FILE__, __LINE__, 0, "Partition identifier too big\n");
    }

    return block->n_elt[id_part];
  }

  else if (id_block >= PDM_BLOCK_ID_BLOCK_POLY2D) {

    _id_block = id_block - PDM_BLOCK_ID_BLOCK_POLY2D;

    PDM_Mesh_nodal_block_poly2d_t *block = pmne->sections_poly2d[_id_block];

    if (block == NULL) {
      PDM_error (__FILE__, __LINE__, 0, "Bad standard block identifier\n");
    }

    if (id_part >= block->n_part) {
      PDM_error(__FILE__, __LINE__, 0, "Partition identifier too big\n");
    }

    return block->n_elt[id_part];
  }

  else {

    _id_block = id_block - PDM_BLOCK_ID_BLOCK_STD;

    PDM_Mesh_nodal_block_std_t *block = pmne->sections_std[_id_block];

    if (block == NULL) {
      PDM_error (__FILE__, __LINE__, 0, "Bad standard block identifier\n");
    }

    if (id_part >= block->n_part) {
      PDM_error(__FILE__, __LINE__, 0, "Partition identifier too big\n");
    }

    return block->n_elt[id_part];
  }

}


int
PDM_part_mesh_nodal_elmts_n_section_get
(
  PDM_part_mesh_nodal_elmts_t *pmne
)
{
  return pmne->n_section;
}


int *
PDM_part_mesh_nodal_elmts_sections_id_get
(
  PDM_part_mesh_nodal_elmts_t *pmne
)
{
  return pmne->sections_id;
}


void
PDM_part_mesh_nodal_elmts_free
(
 PDM_part_mesh_nodal_elmts_t* pmne
)
{

  if(pmne->n_elmts != NULL) {
    free(pmne->n_elmts);
  }

  /* free standard blocks */
  if (pmne->sections_std != NULL) {
    for (int i = 0; i < pmne->n_section_std; i++) {
      _block_std_free(pmne->sections_std[i]);
    }
    free(pmne->sections_std);
  }

  assert(pmne->n_section_poly2d == 0);
  assert(pmne->n_section_poly3d == 0);

  /* Free polygon blocks */
  // if (pmne->sections_poly2d != NULL) {
  //   for (int i = 0; i < pmne->n_section_poly2d; i++) {
  //     _block_poly2d_free(pmne->sections_poly2d[i]);
  //   }
  //   free(pmne->sections_poly2d);
  // }

  /* Free polyhedron blocks */
  // if (pmne->sections_poly3d != NULL) {
  //   for (int i = 0; i < pmne->n_section_poly3d; i++) {
  //     _block_poly3d_free(pmne->sections_poly3d[i]);
  //   }
  //   free(pmne->sections_poly3d);
  // }

  if(pmne->sections_id != NULL) {
    free(pmne->sections_id);
  }

  free(pmne);
}
