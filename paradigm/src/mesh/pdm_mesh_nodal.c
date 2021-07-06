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
#include "pdm_mesh_nodal.h"
#include "pdm_mesh_nodal_priv.h"
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
PDM_Mesh_nodal_vtx_t *
_vtx_free
(
 PDM_Mesh_nodal_vtx_t *vtx
 )
{

  if (vtx != NULL) {
    if (vtx->parent != NULL) {
      vtx->parent =_vtx_free (vtx->parent);
    }

    if (vtx->coords != NULL) {
      free (vtx->coords);
      vtx->coords = NULL;
    }

    free (vtx);
  }
  return NULL;
}


/**
 *
 * \brief Initialize a mesh
 *
 * \param [inout]  mesh        Mesh
 * \param [in]     n_part      Number of partitions
 */

static void
_mesh_init
(
 PDM_Mesh_nodal_t   *mesh,
 const int           n_part,
 const PDM_MPI_Comm  comm
 )
{

  mesh->n_part                   = n_part;
  mesh->numabs                   = NULL;
  mesh->vtx                      = malloc(n_part * sizeof(PDM_Mesh_nodal_vtx_t *));
  mesh->n_cell                   = malloc(n_part * sizeof(int));
  for (int i = 0; i < n_part; i++) {
    mesh->vtx[i] = malloc(sizeof(PDM_Mesh_nodal_vtx_t));
    mesh->vtx[i]->_coords    = NULL;
    mesh->vtx[i]->_numabs    = NULL;
    mesh->vtx[i]->_numparent = NULL;
    mesh->vtx[i]->n_vtx      = 0;
    mesh->vtx[i]->parent     = NULL;
    mesh->vtx[i]->coords     = NULL;
    mesh->n_cell[i]          = 0;
  }

  mesh->n_block_std              = 0;
  mesh->n_block_poly2d           = 0;
  mesh->n_block_poly3d           = 0;

  mesh->blocks_std               = NULL;
  mesh->blocks_poly2d            = NULL;
  mesh->blocks_poly3d            = NULL;

  mesh->pdm_mpi_comm             = comm;
  mesh->prepa_blocks             = NULL;
  mesh->num_cell_parent_to_local = NULL;
  mesh->blocks_id                = NULL;
  mesh->n_blocks                 = 0;
  mesh->is_vtx_def_from_parent   = 0;
}

/**
 *
 * \brief Update blocks identifier list
 *
 * \param [inout]  mesh        Mesh
 */

static void
_update_blocks_id
(
 PDM_Mesh_nodal_t *mesh
)
{
  int n_blocks = 0;

  if (mesh->blocks_std != NULL) {
    n_blocks += mesh->n_block_std;
  }

  if (mesh->blocks_poly2d != NULL) {
    n_blocks += mesh->n_block_poly2d;
  }

  if (mesh->blocks_poly3d != NULL) {
    n_blocks += mesh->n_block_poly3d;
  }

  if (mesh->n_blocks < n_blocks) {
    mesh->blocks_id = (int *) realloc(mesh->blocks_id, sizeof(int) * n_blocks);
  }

  int k = 0;
  if (mesh->blocks_std != NULL) {
    for (int i = 0; i < mesh->n_block_std; i++) {
      mesh->blocks_id[k++] = i + PDM_BLOCK_ID_BLOCK_STD;
    }
  }

  if (mesh->blocks_poly2d != NULL) {
    for (int i = 0; i < mesh->n_block_poly2d; i++) {
      mesh->blocks_id[k++] = i + PDM_BLOCK_ID_BLOCK_POLY2D;
    }
  }

  if (mesh->blocks_poly3d != NULL) {
    for (int i = 0; i < mesh->n_block_poly3d; i++) {
      mesh->blocks_id[k++] = i + PDM_BLOCK_ID_BLOCK_POLY3D;
    }
  }

  mesh->n_blocks = n_blocks;

}


/**
 *
 * \brief Free partially a standard block
 *
 * \param [inout]  _bloc_std    Standard block
 *
 */

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
    if (_block_std->st_free_data == PDM_TRUE) {
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
    if (_block_std->st_free_data == PDM_TRUE) {
      for (int i = 0; i < _block_std->n_part; i++) {
        if (_block_std->_numabs[i] != NULL)
          free(_block_std->_numabs[i]);
        _block_std->_numabs[i] = NULL;
      }
    }
    free(_block_std->_numabs);
    _block_std->_numabs = NULL;
  }

  if (_block_std->_num_part != NULL) {
    if (_block_std->st_free_data == PDM_TRUE) {
      for (int i = 0; i < _block_std->n_part; i++) {
        if (_block_std->_num_part[i] != NULL)
          free(_block_std->_num_part[i]);
        _block_std->_num_part[i] = NULL;
      }
    }
    free(_block_std->_num_part);
    _block_std->_num_part = NULL;
  }

}


/**
 *
 * \brief Free a standard block
 *
 * \param [inout]  _bloc_std    Standard block
 *
 * \return         Null
 *
 */

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
    if (_block_std->st_free_data == PDM_TRUE) {
      for (int i = 0; i < _block_std->n_part; i++) {
        if (_block_std->_parent_num[i] != NULL)
          free(_block_std->_parent_num[i]);
        _block_std->_parent_num[i] = NULL;
      }
    }
    free(_block_std->_parent_num);
    _block_std->_parent_num = NULL;
  }

  free(_block_std);
  return NULL;
}


/**
 *
 * \brief Free partially a polygon block
 *
 * \param [inout]  _block_poly2d   polygon block
 *
 */

static
void
_block_poly2d_free_partial
(
 PDM_Mesh_nodal_block_poly2d_t *_block_poly2d
 )
{

  if (_block_poly2d->_connec_idx != NULL) {
    if (_block_poly2d->st_free_data == PDM_TRUE) {
      for (int i = 0; i < _block_poly2d->n_part; i++) {
        if (_block_poly2d->_connec_idx[i] != NULL)
          free(_block_poly2d->_connec_idx[i]);
        _block_poly2d->_connec_idx[i] = NULL;
      }
    }
    free(_block_poly2d->_connec_idx);
    _block_poly2d->_connec_idx = NULL;
  }

  if (_block_poly2d->_connec != NULL) {
    if (_block_poly2d->st_free_data == PDM_TRUE) {
      for (int i = 0; i < _block_poly2d->n_part; i++) {
        if (_block_poly2d->_connec[i] != NULL)
          free(_block_poly2d->_connec[i]);
        _block_poly2d->_connec[i] = NULL;
      }
    }
    free(_block_poly2d->_connec);
    _block_poly2d->_connec = NULL;
  }

  if (_block_poly2d->_num_part != NULL) {
    if (_block_poly2d->st_free_data == PDM_TRUE) {
      for (int i = 0; i < _block_poly2d->n_part; i++) {
        if (_block_poly2d->_num_part[i] != NULL)
          free(_block_poly2d->_num_part[i]);
        _block_poly2d->_num_part[i] = NULL;
      }
    }
    free(_block_poly2d->_num_part);
    _block_poly2d->_num_part = NULL;
  }

  if (_block_poly2d->_numabs != NULL) {
    if (_block_poly2d->st_free_data == PDM_TRUE) {
      for (int i = 0; i < _block_poly2d->n_part; i++) {
        if (_block_poly2d->_numabs[i] != NULL)
          free(_block_poly2d->_numabs[i]);
        _block_poly2d->_numabs[i] = NULL;
      }
    }
    free(_block_poly2d->_numabs);
    _block_poly2d->_numabs = NULL;
  }

}


/**
 *
 * \brief Free a polygon block
 *
 * \param [inout]  _bloc_poly2d    Polygon block
 *
 * \return         Null
 *
 */

static
PDM_Mesh_nodal_block_poly2d_t *
_block_poly2d_free
(
 PDM_Mesh_nodal_block_poly2d_t *_block_poly2d
 )
{
  _block_poly2d_free_partial(_block_poly2d);

  if (_block_poly2d->n_elt != NULL) {
    free(_block_poly2d->n_elt);
    _block_poly2d->n_elt = NULL;
  }

  if (_block_poly2d->numabs_int != NULL) {
    for (int j = 0; j < _block_poly2d->n_part; j++) {
      if (_block_poly2d->numabs_int[j] != NULL) {
        free(_block_poly2d->numabs_int[j]);
      }
    }
    free(_block_poly2d->numabs_int);
    _block_poly2d->numabs_int = NULL;
  }

  if (_block_poly2d->cell_centers != NULL) {
    for (int j = 0; j < _block_poly2d->n_part; j++) {
      if (_block_poly2d->cell_centers[j] != NULL) {
        free(_block_poly2d->cell_centers[j]);
      }
    }
    free(_block_poly2d->cell_centers);
    _block_poly2d->cell_centers = NULL;
  }

  if (_block_poly2d->_parent_num != NULL) {
    if (_block_poly2d->st_free_data == PDM_TRUE) {
      for (int i = 0; i < _block_poly2d->n_part; i++) {
        if (_block_poly2d->_parent_num[i] != NULL)
          free(_block_poly2d->_parent_num[i]);
        _block_poly2d->_parent_num[i] = NULL;
      }
    }
    free(_block_poly2d->_parent_num);
    _block_poly2d->_parent_num = NULL;
  }

  free(_block_poly2d);

  return NULL;
}


/**
 *
 * \brief Free partially a polyhedron block
 *
 * \param [inout]  _block_poly3d   polyhedron block
 *
 */

static
void
_block_poly3d_free_partial
(
 PDM_Mesh_nodal_block_poly3d_t *_block_poly3d
 )
{

  if (_block_poly3d->_facvtx_idx != NULL) {
    if (_block_poly3d->st_free_data == PDM_TRUE) {
      for (int i = 0; i < _block_poly3d->n_part; i++) {
        if (_block_poly3d->_facvtx_idx[i] != NULL)
          free(_block_poly3d->_facvtx_idx[i]);
        _block_poly3d->_facvtx_idx[i] = NULL;
      }
    }
    free(_block_poly3d->_facvtx_idx);
    _block_poly3d->_facvtx_idx = NULL;
  }

  if (_block_poly3d->_facvtx != NULL) {
    if (_block_poly3d->st_free_data == PDM_TRUE) {
      for (int i = 0; i < _block_poly3d->n_part; i++) {
        if (_block_poly3d->_facvtx[i] != NULL)
          free(_block_poly3d->_facvtx[i]);
        _block_poly3d->_facvtx[i] = NULL;
      }
    }
    free(_block_poly3d->_facvtx);
    _block_poly3d->_facvtx = NULL;
  }

  if (_block_poly3d->_cellfac_idx != NULL) {
    if (_block_poly3d->st_free_data == PDM_TRUE) {
      for (int i = 0; i < _block_poly3d->n_part; i++) {
        if (_block_poly3d->_cellfac_idx[i] != NULL)
          free(_block_poly3d->_cellfac_idx[i]);
        _block_poly3d->_cellfac_idx[i] = NULL;
      }
    }
    free(_block_poly3d->_cellfac_idx);
    _block_poly3d->_cellfac_idx = NULL;
  }

  if (_block_poly3d->_cellfac != NULL) {
    if (_block_poly3d->st_free_data == PDM_TRUE) {
      for (int i = 0; i < _block_poly3d->n_part; i++) {
        if (_block_poly3d->_cellfac[i] != NULL)
          free(_block_poly3d->_cellfac[i]);
        _block_poly3d->_cellfac[i] = NULL;
      }
    }
    free(_block_poly3d->_cellfac);
    _block_poly3d->_cellfac = NULL;
  }

  if (_block_poly3d->_cellvtx_idx != NULL) {
    if (_block_poly3d->st_free_data == PDM_TRUE) {
      for (int i = 0; i < _block_poly3d->n_part; i++) {
        if (_block_poly3d->_cellvtx_idx[i] != NULL)
          free(_block_poly3d->_cellvtx_idx[i]);
        _block_poly3d->_cellvtx_idx[i] = NULL;
      }
    }
    free(_block_poly3d->_cellvtx_idx);
    _block_poly3d->_cellvtx_idx = NULL;
  }

  if (_block_poly3d->_cellvtx != NULL) {
    if (_block_poly3d->st_free_data == PDM_TRUE) {
      for (int i = 0; i < _block_poly3d->n_part; i++) {
        if (_block_poly3d->_cellvtx[i] != NULL)
          free(_block_poly3d->_cellvtx[i]);
        _block_poly3d->_cellvtx[i] = NULL;
      }
    }
    free(_block_poly3d->_cellvtx);
    _block_poly3d->_cellvtx = NULL;
  }

  if (_block_poly3d->_numabs != NULL) {
    if (_block_poly3d->st_free_data == PDM_TRUE) {
      for (int i = 0; i < _block_poly3d->n_part; i++) {
        if (_block_poly3d->_numabs[i] != NULL)
          free(_block_poly3d->_numabs[i]);
        _block_poly3d->_numabs[i] = NULL;
      }
    }
    free(_block_poly3d->_numabs);
    _block_poly3d->_numabs = NULL;
  }
}


/**
 *
 * \brief Free a polyhedron block
 *
 * \param [inout]  _block_poly3d    Polyhedron block
 *
 * \return         Null
 *
 */

static
PDM_Mesh_nodal_block_poly3d_t *
_block_poly3d_free
(
 PDM_Mesh_nodal_block_poly3d_t *_block_poly3d
 )
{
  _block_poly3d_free_partial(_block_poly3d);

  if (_block_poly3d->n_elt != NULL) {
    free(_block_poly3d->n_elt);
    _block_poly3d->n_elt = NULL;
  }

  if (_block_poly3d->n_face!= NULL) {
    free(_block_poly3d->n_face);
    _block_poly3d->n_face= NULL;
  }

  if (_block_poly3d->numabs_int != NULL) {
    for (int j = 0; j < _block_poly3d->n_part; j++) {
      if (_block_poly3d->numabs_int[j] != NULL) {
        free(_block_poly3d->numabs_int[j]);
      }
    }
    free(_block_poly3d->numabs_int);
    _block_poly3d->numabs_int = NULL;
  }

  if (_block_poly3d->cell_centers != NULL) {
    for (int j = 0; j < _block_poly3d->n_part; j++) {
      if (_block_poly3d->cell_centers[j] != NULL) {
        free(_block_poly3d->cell_centers[j]);
      }
    }
    free(_block_poly3d->cell_centers);
    _block_poly3d->cell_centers = NULL;
  }

  if (_block_poly3d->_parent_num != NULL) {
    if (_block_poly3d->st_free_data == PDM_TRUE) {
      for (int i = 0; i < _block_poly3d->n_part; i++) {
        if (_block_poly3d->_parent_num[i] != NULL)
          free(_block_poly3d->_parent_num[i]);
        _block_poly3d->_parent_num[i] = NULL;
      }
    }
    free(_block_poly3d->_parent_num);
    _block_poly3d->_parent_num = NULL;
  }

  free(_block_poly3d);
  return NULL;
}


/**
 *
 * \brief  Cross product
 *
 * \param[in]     a    First vector
 * \param[in]     b    Second vector
 * \param[inout]  c    \ref a X \ref b vector
 *
 */

static inline void
_p_cross
(
 const double a[3],
 const double b[3],
 double c[3]
 )
{
  c[0] = a[1] * b[2] - b[1] * a[2];
  c[1] = b[0] * a[2] - a[0] * b[2];
  c[2] = a[0] * b[1] - b[0] * a[1];
}


/**
 *
 * Dot product
 *
 * \param[in]     a    First vector
 * \param[in]     b    Second Vector
 *
 * \return    Dot product
 *
 */

static inline double
_p_dot
(
 const double a[3],
 const double b[3]
 )
{
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}


/**
 *
 * \brief Build tetrahedron nodal connectivity from faces connectivity
 *
 *   \param[in] vtx         Vertices coordinates
 *   \param[in] tria_vtx    Faces connectivity
 *   \param[out] tetra_vtx   Tetrahedron connectivity
 *
 */

static void
_connec_tetra
(
 PDM_Mesh_nodal_vtx_t *vtx,
 PDM_l_num_t *tria_vtx,
 PDM_l_num_t tetra_vtx[]
 )
{

  /* Initialization */

  tetra_vtx[0] = tria_vtx[0];
  tetra_vtx[1] = tria_vtx[1];
  tetra_vtx[2] = tria_vtx[2];

  for (int i = 3; i < 11; i++) {
    if ((tria_vtx[i] != tetra_vtx[0]) &&
        (tria_vtx[i] != tetra_vtx[1]) &&
        (tria_vtx[i] != tetra_vtx[2]))
      tetra_vtx[3] = tria_vtx[i];
  }

  /* Orientation */

  const PDM_real_t *_coords = vtx->_coords;
  double v1[3];
  double v2[3];
  double v3[3];
  double n[3];

  for (int i = 0; i < 3; i++) {
    v1[i] = _coords[3*(tetra_vtx[1] - 1) + i] - _coords[3*(tetra_vtx[0] - 1) + i];
    v2[i] = _coords[3*(tetra_vtx[2] - 1) + i] - _coords[3*(tetra_vtx[0] - 1) + i];
    v3[i] = _coords[3*(tetra_vtx[3] - 1) + i] - _coords[3*(tetra_vtx[0] - 1) + i];
  }

  _p_cross(v1, v2, n);
  double orient = _p_dot(v3, n);

  if (orient < 0) {
    tetra_vtx[0] = tria_vtx[2];
    tetra_vtx[1] = tria_vtx[1];
    tetra_vtx[2] = tria_vtx[0];
  }
}

/**
 *
 * \brief Build prism nodal connectivity from faces connectivity
 *
 *   \param[in] vtx         Vertices coordinates
 *   \param[in] tria_vtx    Faces connectivity
 *   \param[in] quad_vtx    Faces connectivity
 *   \param[out] prism_vtx   Prism connectivity
 *
 */

static void
_connec_prism
(
 PDM_Mesh_nodal_vtx_t *vtx,
 PDM_l_num_t *tria_vtx,
 PDM_l_num_t *quad_vtx,
 PDM_l_num_t prism_vtx[]
 )
{

  /* Initialisation */

  for (int i = 0; i < 6; i++)
    prism_vtx[i] = tria_vtx[i];

  /* Orientation des faces */

  const PDM_real_t *_coords = vtx->_coords;

  double c[6];
  double n[6];

  for (int i = 0; i < 2; i++) {
    for (int k = 0; k < 3; k++)
      c[3*i+k] = 0.;
    for (int j = 0; j < 3; j++) {
      int isom = prism_vtx[3*i+j] - 1;
      for (int k = 0; k < 3; k++)
        c[3*i+k] += _coords[3*isom+k];
    }
    for (int k = 0; k < 3; k++)
      c[3*i+k] *= 1.0/3.0;

    for (int k = 0; k < 3; k++)
      n[3*i+k] = 0.;

    double v1[3];
    double v2[3];
    int isom3 = prism_vtx[3*i+2] - 1 ;
    int isom2 = prism_vtx[3*i+1] - 1;
    int isom1 = prism_vtx[3*i] - 1;

    for (int k = 0; k < 3; k++) {
      v1[k] = _coords[3*isom2+k] - _coords[3*isom1+k];
      v2[k] = _coords[3*isom3+k] - _coords[3*isom1+k];
    }
    _p_cross(v1, v2, n + 3*i);
  }

  double cc[3];
  for (int k = 0; k < 3; k++)
    cc[k] = c[3+k] - c[k];

  double orientation = _p_dot(cc, n);
  double orientation2 = _p_dot(cc, n+3);

  if (orientation < 0) {
    int tmp = prism_vtx[1];
    prism_vtx[1] = prism_vtx[2];
    prism_vtx[2] = tmp;
  }

  if (orientation2 < 0) {
    int tmp = prism_vtx[4];
    prism_vtx[4] = prism_vtx[5];
    prism_vtx[5] = tmp;
  }

  /* Permutation circulaire */

  int id1 = -1;
  for (int j = 0; j < 12; j++) {
    if (quad_vtx[j] == prism_vtx[0]) {
      id1 = j;
      break;
    }
  }

  int id2 = (id1 / 4) * 4 + (id1 + 1) % 4;
  if ((quad_vtx[id2] == prism_vtx[1]) ||
      (quad_vtx[id2] == prism_vtx[2]))
    id2 =  (id1 / 4) * 4 + (id1 + 3) % 4;

  int id_deb;
  for (int j = 0; j < 3; j++) {
    if (quad_vtx[id2] == prism_vtx[3+j]) {
      id_deb = j;
      break;
    }
  }

  int tmp[3];
  for (int j = 0; j < 3; j++)
    tmp[j] = prism_vtx[3+j];

  for (int j = 0; j < 3; j++) {
    int idx = (id_deb + j) % 3;
    prism_vtx[3+j] = tmp[idx];
  }

}


/**
 *
 * \brief Build pyramid nodal connectivity from faces connectivity
 *
 *   \param[in] vtx         Vertices coordinates
 *   \param[in] tria_vtx    Faces connectivity
 *   \param[in] quad_vtx    Faces connectivity
 *   \param[out] pyramid_vtx Pyramid connectivity
 *
 */

static void
_connec_pyramid
(
 PDM_Mesh_nodal_vtx_t *vtx,
 PDM_l_num_t *tria_vtx,
 PDM_l_num_t *quad_vtx,
 PDM_l_num_t pyramid_vtx[]
 )
{

  /* Initialisation */

  pyramid_vtx[0] = quad_vtx[0];
  pyramid_vtx[1] = quad_vtx[1];
  pyramid_vtx[2] = quad_vtx[2];
  pyramid_vtx[3] = quad_vtx[3];

  for (int i = 0; i < 9; i++) {
    if ((tria_vtx[i] != pyramid_vtx[0]) &&
        (tria_vtx[i] != pyramid_vtx[1]) &&
        (tria_vtx[i] != pyramid_vtx[2]) &&
        (tria_vtx[i] != pyramid_vtx[3])) {
      pyramid_vtx[4] = tria_vtx[i];
      break;
    }
  }

  /* Orientation */

  const PDM_real_t *_coords = vtx->_coords;

  double c[3];
  double n[3];

  for (int k = 0; k < 3; k++)
    c[k] = 0.;
  for (int j = 0; j < 4; j++) {
    int isom = pyramid_vtx[j] - 1;
    for (int k = 0; k < 3; k++)
      c[k] += _coords[3*isom+k];
  }
  for (int k = 0; k < 3; k++)
    c[k] *= 0.25;

  for (int k = 0; k < 3; k++)
    n[k] = 0.;

  for (int j = 0; j < 4; j++) {
    int isom = pyramid_vtx[j] - 1;
    int suiv = (j+1) % 4;
    int isom_suiv = pyramid_vtx[suiv] - 1;

    double v1[3];
    double v2[3];
    for (int k = 0; k < 3; k++) {
      v1[k] = _coords[3*isom+k] -  c[k];
      v2[k] = _coords[3*isom_suiv+k] -  c[k];
    }

    _p_cross(v1, v2, n);

  }

  double cc[3];
  for (int k = 0; k < 3; k++)
    cc[k] = _coords[3*(pyramid_vtx[3] - 1) + k] - c[k];

  /* Inversion eventuelle des sens de rotation des faces*/

  double orientation = _p_dot(cc, n);

  if (orientation < 0) {
    int tmp = pyramid_vtx[0];
    pyramid_vtx[0] = pyramid_vtx[3];
    pyramid_vtx[3] = tmp;
    tmp = pyramid_vtx[1];
    pyramid_vtx[1] = pyramid_vtx[2];
    pyramid_vtx[2] = tmp;
  }

}


/**
 *
 * \brief Build hexahedron nodal connectivity from faces connectivity
 *
 *   \param[in] vtx         Vertices coordinates
 *   \param[in] quad_vtx    Faces connectivity
 *   \param[out] hexa_vtx    Hexahedron connectivity
 *
 */

static void
_connec_hexa
(
 PDM_Mesh_nodal_vtx_t *vtx,
 PDM_l_num_t *quad_vtx,
 PDM_l_num_t hexa_vtx[]
)
{

  /* Initialization */

  hexa_vtx[0] = quad_vtx[0];
  hexa_vtx[1] = quad_vtx[1];
  hexa_vtx[2] = quad_vtx[2];
  hexa_vtx[3] = quad_vtx[3];

  PDM_l_num_t face_contact[4];

  for (int i = 1; i < 6; i++) {
    int cpt = 0;
    for (int j = 0; j < 4; j++) {
      PDM_l_num_t som_courant = quad_vtx[4*i+j];
      if ((som_courant != hexa_vtx[0]) &&
          (som_courant != hexa_vtx[1]) &&
          (som_courant != hexa_vtx[2]) &&
          (som_courant != hexa_vtx[3]))
        cpt += 1;
    }
    if (cpt == 4) {
      hexa_vtx[4] = quad_vtx[4*i];
      hexa_vtx[5] = quad_vtx[4*i+1];
      hexa_vtx[6] = quad_vtx[4*i+2];
      hexa_vtx[7] = quad_vtx[4*i+3];
    }
    if (cpt == 2) {
      face_contact[0] = quad_vtx[4*i];
      face_contact[1] = quad_vtx[4*i+1];
      face_contact[2] = quad_vtx[4*i+2];
      face_contact[3] = quad_vtx[4*i+3];
    }
  }

  /* Calcul des centres et normales de la base et de la face opposee */

  const PDM_real_t *_coords = vtx->_coords;

  double c[6];
  double n[6];

  for (int i = 0; i < 2; i++) {
    for (int k = 0; k < 3; k++)
      c[3*i+k] = 0.;
    for (int j = 0; j < 4; j++) {
      int isom = hexa_vtx[4*i+j] - 1;
      for (int k = 0; k < 3; k++)
        c[3*i+k] += _coords[3*isom+k];
    }
    for (int k = 0; k < 3; k++)
      c[3*i+k] *= 0.25;

    for (int k = 0; k < 3; k++)
      n[3*i+k] = 0.;

    for (int j = 0; j < 4; j++) {
      int isom = hexa_vtx[4*i+j] - 1;
      int suiv = (j+1) % 4;
      int isom_suiv = hexa_vtx[4*i+suiv] - 1;

      double v1[3];
      double v2[3];
      for (int k = 0; k < 3; k++) {
        v1[k] = _coords[3*isom+k] -  c[3*i+k];
        v2[k] = _coords[3*isom_suiv+k] -  c[3*i+k];
      }

      _p_cross(v1, v2, n + 3*i);

    }

  }

  double cc[3];
  for (int k = 0; k < 3; k++)
    cc[k] = c[3+k] - c[k];

  /* Inversion eventuelle des sens de rotation des faces*/

  double orientation = _p_dot(cc, n);
  double orientation2 = _p_dot(cc, n+3);

  if (orientation < 0) {
    int tmp = hexa_vtx[0];
    hexa_vtx[0] = hexa_vtx[3];
    hexa_vtx[3] = tmp;
    tmp = hexa_vtx[1];
    hexa_vtx[1] = hexa_vtx[2];
    hexa_vtx[2] = tmp;
  }

  if (orientation2 < 0) {
    int tmp = hexa_vtx[4];
    hexa_vtx[4] = hexa_vtx[7];
    hexa_vtx[7] = tmp;
    tmp = hexa_vtx[5];
    hexa_vtx[5] = hexa_vtx[6];
    hexa_vtx[6] = tmp;
  }

  /* Permutation circulaire eventuelle de la face sup */

  int id1 = -1;
  int k1 = -1;
  for (int k = 0; k < 4; k++) {
    for (int j = 0; j < 4; j++) {
      if (face_contact[j] == hexa_vtx[k]) {
        id1 = j;
        k1 = k;
        break;
      }
      if (id1 != -1)
        break;
    }
  }

  if (k1 == -1) {
    PDM_printf("Error connect_hexa : %d %d %d %d %d %d %d %d\n",
               hexa_vtx[0],
               hexa_vtx[1],
               hexa_vtx[2],
               hexa_vtx[3],
               hexa_vtx[4],
               hexa_vtx[5],
               hexa_vtx[6],
               hexa_vtx[7]);

    for (int i10 = 0; i10 < 4; i10++) {
      PDM_printf("   face %d : %d %d %d %d\n", i10+1, quad_vtx[4*i10],
                 quad_vtx[4*i10+1],
                 quad_vtx[4*i10+2],
                 quad_vtx[4*i10+3]);
    }
    abort();

  }

  int id2 = (id1 + 1) % 4;
  int k2 = (k1 + 1) % 4;
  int k3 = (k1 + 3) % 4;

  if ((face_contact[id2] == hexa_vtx[k2]) ||
      (face_contact[id2] == hexa_vtx[k3]))
    id2 = (id1 + 3) % 4;

  int id_deb;
  for (int j = 0; j < 4; j++) {
    if (face_contact[id2] == hexa_vtx[4+j]) {
      id_deb = (j - k1);
      if (id_deb < 0)
        id_deb += 4;
      id_deb = id_deb % 4;
      break;
    }
  }

  int tmp[4];
  for (int j = 0; j < 4; j++)
    tmp[j] = hexa_vtx[4+j];

  for (int j = 0; j < 4; j++) {
    int idx = (id_deb + j) % 4;
    hexa_vtx[4+j] = tmp[idx];
  }
}


/**
 * \brief Get element type
 *
 *   \param[in] n_face_cell     Number of faces in the current cell
 *   \param[in]  face_cell      Face to cell connectivity
 *   \param[in]  face_cell_idx  Face to vertex connectivity index
 *   \param[in]  face_cell_n    Number of vertices for each face
 *   \param[in]  face_vtx       Face to vertex connectivity
 *   \param[out] tria_vtx       Quadrangles connectivity
 *   \param[out] quad_vtx       Quadrangles connectivity
 *
 *  \return   Cell type
 */

inline static
PDM_Mesh_nodal_elt_t
_type_cell_3D
(
 const int             n_face_cell,
 const PDM_l_num_t    *cell_face,
 const PDM_l_num_t    *face_vtx_idx,
 const PDM_l_num_t    *face_vtx_nb,
 const PDM_l_num_t    *face_vtx,
 PDM_l_num_t           tria_vtx[],
 PDM_l_num_t           quad_vtx[]
)
{
  int adjust = 0;
  if (face_vtx_idx[0] == 1) {
    adjust = 1;
  }


  PDM_l_num_t  n_trias = 0;
  PDM_l_num_t  n_quads = 0;

  if (n_face_cell > 6) {
    return PDM_MESH_NODAL_POLY_3D;
  }

  for (int i = 0; i < n_face_cell; i++) {

    const int face_id = PDM_ABS(cell_face[i]) - 1;
    const int n_som_face = face_vtx_nb[face_id];
    PDM_l_num_t idx = face_vtx_idx[face_id] - adjust;

    if (n_som_face == 3) {
      PDM_l_num_t *cell_som_tria_courant = tria_vtx + 3*n_trias;
      for (int j = idx; j < idx + n_som_face; j++) {
        cell_som_tria_courant[j-idx] = face_vtx[j];
      }
      n_trias += 1;
    }
    else if (n_som_face == 4) {
      PDM_l_num_t *cell_som_quad_courant = quad_vtx + 4*n_quads;
      for (int j = idx; j < idx + n_som_face; j++) {
        cell_som_quad_courant[j-idx] = face_vtx[j];
      }
      n_quads += 1;
    }
    else
      return PDM_MESH_NODAL_POLY_3D;

  }

  PDM_Mesh_nodal_elt_t cell_type;

  if ((n_quads == 0) && (n_trias == 4))
    cell_type = PDM_MESH_NODAL_TETRA4;
  else if (n_quads == 6)
    cell_type = PDM_MESH_NODAL_HEXA8;
  else if ((n_quads == 1) && (n_trias == 4))
    cell_type = PDM_MESH_NODAL_PYRAMID5;
  else if ((n_quads == 3) && (n_trias == 2)) {
    int trias[6];
    n_trias = 0;
    for (int i = 0; i < n_face_cell; i++) {

      const int face_id = PDM_ABS(cell_face[i]) - 1;
      const int ideb = face_vtx_idx[face_id] - adjust;

      const int n_som_face = face_vtx_nb[face_id];

      if (n_som_face == 3) {
        for (int j = 0; j < 3; j++) {
          trias[3*n_trias+j] = face_vtx[ideb+j];
        }
        n_trias += 1;
      }
      if (n_trias >= 2)
        break;
    }

    cell_type = PDM_MESH_NODAL_PRISM6;
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        if (trias[i] == trias[3+j]) {
          cell_type = PDM_MESH_NODAL_POLY_3D;
          break;
        }
      }
      if (cell_type == PDM_MESH_NODAL_POLY_3D)
        break;
    }
  }

  else {
    cell_type = PDM_MESH_NODAL_POLY_3D;
  }

  return cell_type;

}


/*=============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief Free vtx structure
 *
 * \param[inout]  vtx    Vertices
 *
 * \return        NULL
 *
 */

int
PDM_Mesh_nodal_is_2D_element
(
  PDM_Mesh_nodal_elt_t type
)
{
  return type==PDM_MESH_NODAL_TRIA3
      || type==PDM_MESH_NODAL_QUAD4
      || type==PDM_MESH_NODAL_POLY_2D;
}

int
PDM_Mesh_nodal_is_3D_element
(
  PDM_Mesh_nodal_elt_t type
)
{
  return type==PDM_MESH_NODAL_TETRA4
      || type==PDM_MESH_NODAL_PYRAMID5
      || type==PDM_MESH_NODAL_PRISM6
      || type==PDM_MESH_NODAL_HEXA8
      || type==PDM_MESH_NODAL_POLY_3D;
}


/**
 * \brief Get the number of vertices of an element type
 *
 * \param [in]   type     Element type
 * \param [in]   comm     Element order
 *
 * \return       Number of vertices
 *
 */

int
PDM_Mesh_nodal_n_vtx_elt_get
(
  PDM_Mesh_nodal_elt_t type,
  const int order
)
{
  if (order != 1) {
     PDM_error (__FILE__, __LINE__, 0, "Not implemented yet for order > 1\n");
  }

  switch (type) {
  case PDM_MESH_NODAL_POINT :
    return 1;
    break;
  case PDM_MESH_NODAL_BAR2 :
    return 2;
    break;
  case PDM_MESH_NODAL_TRIA3 :
    return 3;
    break;
  case PDM_MESH_NODAL_QUAD4 :
    return 4;
    break;
  case PDM_MESH_NODAL_TETRA4 :
    return 4;
    break;
  case PDM_MESH_NODAL_PYRAMID5 :
    return 5;
    break;
  case PDM_MESH_NODAL_PRISM6 :
    return 6;
    break;
  case PDM_MESH_NODAL_HEXA8 :
    return 8;
    break;
  default :
    PDM_error (__FILE__, __LINE__, 0, "Unknown for order Poly2D and Polyu3D\n");
  }
  return -1;
}



/**
 * \brief Create a Mesh nodal structure
 *
 * \param [in]   n_part   Number of partition on the current process
 * \param [in]   comm     MPI communicator
 *
 * \return       New mesh nodal handle
 *
 */

PDM_Mesh_nodal_t*
PDM_Mesh_nodal_create
(
 const int     n_part,
 const PDM_MPI_Comm comm
)
{
  PDM_Mesh_nodal_t *mesh = (PDM_Mesh_nodal_t *) malloc (sizeof(PDM_Mesh_nodal_t));

  _mesh_init (mesh, n_part, comm);

  return mesh;
}


/**
 * \brief  Return number of partitions
 *
 * \param [in]  idx            Nodal mesh handle
 *
 * \return  Number of partitions
 *
 */

int
PDM_Mesh_nodal_n_part_get
(
PDM_Mesh_nodal_t *mesh
)
{
  return mesh->n_part;
}


/**
 * \brief Free partially a nodal mesh structure
 *
 * \param [in]  idx   Nodal mesh handle
 *
 * \return      NULL
 *
 */

void
PDM_Mesh_nodal_partial_free
(
PDM_Mesh_nodal_t *mesh
)
{

  if (mesh->blocks_std != NULL) {
    for (int i = 0; i < mesh->n_block_std; i++) {
      _block_std_free_partial(mesh->blocks_std[i]);
    }
  }

  if (mesh->blocks_poly2d != NULL) {
    for (int i = 0; i < mesh->n_block_poly2d; i++) {
      _block_poly2d_free_partial(mesh->blocks_poly2d[i]);
    }
  }

  if (mesh->blocks_poly3d != NULL) {
    for (int i = 0; i < mesh->n_block_poly3d; i++) {
      _block_poly3d_free_partial(mesh->blocks_poly3d[i]);
    }
  }
}


/**
 * \brief Get the cell global numbering taking into account parent_num
 *
 * \param [in]  idx   Nodal mesh handle
 *
 * \return      NULL
 *
 */

PDM_g_num_t *
PDM_Mesh_nodal_g_num_get_from_part
(
PDM_Mesh_nodal_t *mesh,
const int i_part
)
{


  if (mesh->numabs == NULL) {
    mesh->numabs = malloc (sizeof(PDM_g_num_t*)*mesh->n_part);
    int is_not_parent_num = (PDM_Mesh_nodal_block_parent_num_get(mesh, mesh->blocks_id[0], 0) == NULL);
    for (int i = 0; i < mesh->n_part; i++) {
      for (int i1 = 0; i1 < mesh->n_blocks; i1++) {
        assert (is_not_parent_num == (PDM_Mesh_nodal_block_parent_num_get(mesh, mesh->blocks_id[i1], i) == NULL));
      }
    }

    if (is_not_parent_num) {

      for (int i = 0; i < mesh->n_part; i++) {
        int k = 0;
        mesh->numabs[i] = malloc (sizeof(PDM_g_num_t)*mesh->n_cell[i]);
        for (int i1 = 0; i1 < mesh->n_block_std; i1++) {
          for (int i2 = 0; i2 < mesh->blocks_std[i1]->n_elt[i]; i2++) {
            mesh->numabs[i][k++] = mesh->blocks_std[i1]->_numabs[i][i2];
          }
        }
        for (int i1 = 0; i1 < mesh->n_block_poly2d; i1++) {
          for (int i2 = 0; i2 < mesh->blocks_poly2d[i1]->n_elt[i]; i2++) {
            mesh->numabs[i][k++] = mesh->blocks_poly2d[i1]->_numabs[i][i2];
          }
        }
        for (int i1 = 0; i1 < mesh->n_block_poly3d; i1++) {
          for (int i2 = 0; i2 < mesh->blocks_poly3d[i1]->n_elt[i]; i2++) {
            mesh->numabs[i][k++] = mesh->blocks_poly3d[i1]->_numabs[i][i2];
          }
        }
      }
    }

    else {
      for (int i = 0; i < mesh->n_part; i++) {
        mesh->numabs[i] = malloc (sizeof(PDM_g_num_t)*mesh->n_cell[i]);
        for (int i1 = 0; i1 < mesh->n_block_std; i1++) {
          for (int i2 = 0; i2 < mesh->blocks_std[i1]->n_elt[i]; i2++) {
            mesh->numabs[i][mesh->blocks_std[i1]->_parent_num[i][i2]] = mesh->blocks_std[i1]->_numabs[i][i2];
          }
        }
        for (int i1 = 0; i1 < mesh->n_block_poly2d; i1++) {
          for (int i2 = 0; i2 < mesh->blocks_poly2d[i1]->n_elt[i]; i2++) {
            mesh->numabs[i][mesh->blocks_poly2d[i1]->_parent_num[i][i2]] = mesh->blocks_poly2d[i1]->_numabs[i][i2];
          }
        }
        for (int i1 = 0; i1 < mesh->n_block_poly3d; i1++) {
          for (int i2 = 0; i2 < mesh->blocks_poly3d[i1]->n_elt[i]; i2++) {
            mesh->numabs[i][mesh->blocks_poly3d[i1]->_parent_num[i][i2]] = mesh->blocks_poly3d[i1]->_numabs[i][i2];
          }
        }
      }
    }
  }
  return mesh->numabs[i_part];
}

/**
 * \brief Free a nodal mesh structure
 *
 * \param [in]  idx   Nodal mesh handle
 *
 * \return      NULL
 *
 */

void
PDM_Mesh_nodal_free
(
PDM_Mesh_nodal_t *mesh
)
{

  PDM_Mesh_nodal_partial_free (mesh);

  if (mesh != NULL) {

    if (mesh->numabs != NULL) {
      for (int i = 0; i < mesh->n_part; i++) {
        free (mesh->numabs[i]);
      }
      free (mesh->numabs);
    }

    if (mesh->blocks_id != NULL) {
      free (mesh->blocks_id);
    }

    mesh->blocks_id = NULL;

    /* Free vertices */

    if (mesh->vtx != NULL) {
      for (int i = 0; i < mesh->n_part; i++) {
        mesh->vtx[i] = _vtx_free (mesh->vtx[i]);
      }

      free(mesh->vtx);
      mesh->vtx = NULL;
    }

    /* free standard blocks */
    if (mesh->blocks_std != NULL) {
      for (int i = 0; i < mesh->n_block_std; i++) {
        _block_std_free(mesh->blocks_std[i]);
      }
      free(mesh->blocks_std);
    }

    /* Free polygon blocks */
    if (mesh->blocks_poly2d != NULL) {
      for (int i = 0; i < mesh->n_block_poly2d; i++) {
        _block_poly2d_free(mesh->blocks_poly2d[i]);
      }
      free(mesh->blocks_poly2d);
    }

    /* Free polyhedron blocks */
    if (mesh->blocks_poly3d != NULL) {
      for (int i = 0; i < mesh->n_block_poly3d; i++) {
        _block_poly3d_free(mesh->blocks_poly3d[i]);
      }
      free(mesh->blocks_poly3d);
    }

    /* Free structure */
    if (mesh->num_cell_parent_to_local != NULL) {
      for (int i_part = 0; i_part < mesh->n_part; i_part++) {
        if (mesh->num_cell_parent_to_local[i_part] != NULL)
          free(mesh->num_cell_parent_to_local[i_part]);
      }
      free(mesh->num_cell_parent_to_local);
      mesh->num_cell_parent_to_local = NULL;
    }

    free(mesh->n_cell);
    mesh->n_cell = NULL;

    if (mesh->blocks_id != NULL) {
      free (mesh->blocks_id);
    }

    free(mesh);

  }
}


/**
 * \brief Define partition vertices
 *
 * \param [in]  idx       Nodal mesh handle
 * \param [in]  id_part   Partition identifier
 * \param [in]  n_vtx     Number of vertices
 * \param [in]  coords    Interlaced coordinates (size = 3 * \ref n_vtx)
 * \param [in]  numabs    Global numbering
 *
 */

void
PDM_Mesh_nodal_coord_set
(
       PDM_Mesh_nodal_t *mesh,
 const int               id_part,
 const int               n_vtx,
 const PDM_real_t       *coords,
 const PDM_g_num_t      *numabs
)
{

  if (id_part >= mesh->n_part) {
    PDM_error (__FILE__, __LINE__, 0, "Bad part identifier\n");
  }

  PDM_Mesh_nodal_vtx_t *vtx = mesh->vtx[id_part];

  if ((vtx->_coords != NULL) ||
      (vtx->_numabs != NULL)) {
    PDM_error(__FILE__, __LINE__, 0, "these partition vertices are already defined\n");
  }

  /* Mapping memoire */

  vtx->n_vtx   = n_vtx;
  vtx->_coords = coords;
  vtx->_numabs = numabs;

}


/**
 * \brief  Return number of vertices
 *
 * \param [in]  idx       Nodal mesh handle
 * \param [in]  id_part   Partition identifier
 *
 * \return  Number of vertices
 *
 */

int
PDM_Mesh_nodal_n_vertices_get
(
       PDM_Mesh_nodal_t *mesh,
 const int               id_part
)
{
  if (id_part >= mesh->n_part) {
    PDM_error (__FILE__, __LINE__, 0, "Bad part identifier %d %d\n", id_part, mesh->n_part);
  }

  PDM_Mesh_nodal_vtx_t *vtx = mesh->vtx[id_part];

  return vtx->n_vtx;
}


/**
 * \brief  Return parent num of vertices
 *
 * \param [in]  mesh           Nodal mesh
 *
 * \return  Parent of vertices
 *
 */

const int *
PDM_Mesh_nodal_vertices_parent_get
(
       PDM_Mesh_nodal_t *mesh,
 const int               id_part
 )
{
  if (id_part >= mesh->n_part) {
    PDM_error (__FILE__, __LINE__, 0, "Bad part identifier\n");
  }

  PDM_Mesh_nodal_vtx_t *vtx = mesh->vtx[id_part];

  return vtx->_numparent;
}


/**
 * \brief  Return cell centers
 *
 * \param [in]  idx       Nodal mesh handle
 * \param [in]  id_part   Partition identifier
 * \param [in]  id_block   Block identifier
 *
 * \return  Return cell centers
 *
 */

const double *
PDM_Mesh_cell_centers_get
(
       PDM_Mesh_nodal_t *mesh,
 const int               id_block,
 const int               id_part
)
{


  double* elt_centers;
  int _id_block;
  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }

  if (id_block < PDM_BLOCK_ID_BLOCK_POLY2D) {
    _id_block = id_block;

    const PDM_Mesh_nodal_block_std_t *block = mesh->blocks_std[_id_block];

    if (block == NULL) {
      PDM_error (__FILE__, __LINE__, 0, "Bad block identifier\n");
    }

    elt_centers = block->cell_centers[id_part] ;
  }

  else if (id_block < PDM_BLOCK_ID_BLOCK_POLY3D) {

    _id_block = id_block - PDM_BLOCK_ID_BLOCK_POLY2D;

    const PDM_Mesh_nodal_block_poly2d_t *block = mesh->blocks_poly2d[_id_block];

    if (block == NULL) {
      PDM_error (__FILE__, __LINE__, 0, "Bad block identifier\n");
    }
    elt_centers = block->cell_centers[id_part] ;
  }

  else {
    _id_block = id_block - PDM_BLOCK_ID_BLOCK_POLY3D;

    const PDM_Mesh_nodal_block_poly3d_t *block = mesh->blocks_poly3d[_id_block];

    elt_centers = block->cell_centers[id_part] ;

  }

  return elt_centers;

}


/**
 * \brief  Return coordinates of vertices
 *
 * \param [in]  mesh           Nodal mesh
 *
 * \return  Coordinates of vertices
 *
 */

const double *
PDM_Mesh_nodal_vertices_get
(
       PDM_Mesh_nodal_t *mesh,
 const int               id_part
 )
{
  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }

  if (id_part >= mesh->n_part) {
    PDM_error (__FILE__, __LINE__, 0, "Bad part identifier\n");
  }

  PDM_Mesh_nodal_vtx_t *vtx = mesh->vtx[id_part];

  return vtx->_coords;
}


/**
 * \brief  Return global numbering of vertices
 *
 * \param [in]  idx       Nodal mesh handle
 * \param [in]  id_part   Partition identifier
 *
 * \return  Global numbering of vertices
 *
 */

const PDM_g_num_t *
PDM_Mesh_nodal_vertices_g_num_get
(
       PDM_Mesh_nodal_t *mesh,
 const int               id_part
 )
{
  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }

  if (id_part >= mesh->n_part) {
    PDM_error (__FILE__, __LINE__, 0, "Bad part identifier\n");
  }

  PDM_Mesh_nodal_vtx_t *vtx = mesh->vtx[id_part];

  return vtx->_numabs;
}


/**
 * \brief Extract vertices from parent vertices
 *
 * \param [in]  mesh           Nodal mesh
 *
 * \return true if the vertices are defined from parents
 */

int
PDM_Mesh_nodal_is_set_coord_from_parent
(
 PDM_Mesh_nodal_t *mesh
)
{

  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }

  return mesh->is_vtx_def_from_parent;

}

/**
 * \brief Extract vertices from parent vertices
 *
 * \param [in]  mesh           Nodal mesh
 * \param [in]  id_part        Partition identifier
 * \param [in]  n_vtx          Number of vertices
 * \param [in]  n_vtx_parent   Number of parent vertices
 * \param [in]  numabs         Global numbering (size = \ref n_vtx)
 * \param [in]  num_parent     Numbering in the parent numbering (size = \ref n_vtx)
 * \param [in]  coords_parent  Parent interlaced coordinates (size = 3 * \ref n_vtx_parent)
 * \param [in]  numabs_parent  Parent global numbering (size = \ref n_vtx_parent)
 *
 */

void
PDM_Mesh_nodal_coord_from_parent_set
(
       PDM_Mesh_nodal_t *mesh,
 const int               id_part,
 const int               n_vtx,
 const int               n_vtx_parent,
 const PDM_g_num_t      *numabs,
 const int              *num_parent,
 const PDM_real_t       *coords_parent,
 const PDM_g_num_t      *numabs_parent
)
{

  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }

  if (id_part >= mesh->n_part) {
    PDM_error (__FILE__, __LINE__, 0, "Bad part identifier\n");
  }

  PDM_Mesh_nodal_vtx_t *vtx = mesh->vtx[id_part];

  if ((vtx->_coords != NULL) ||
      (vtx->_numabs != NULL)) {
    PDM_error(__FILE__, __LINE__, 0, "Vertices are already defined\n");
  }

  vtx->parent = (PDM_Mesh_nodal_vtx_t *) malloc (sizeof (PDM_Mesh_nodal_vtx_t));
  PDM_Mesh_nodal_vtx_t *_parent = vtx->parent;
  _parent->parent = NULL;
  _parent->n_vtx = n_vtx_parent;
  _parent->coords = NULL;
  _parent->_coords = coords_parent ;
  _parent->_numabs = numabs_parent;
  _parent->_numparent = NULL;

  vtx->n_vtx      = n_vtx;
  vtx->coords     = malloc (sizeof(double) * 3 * n_vtx);
  vtx->_coords    = vtx->coords;
  vtx->_numabs    = numabs;
  vtx->_numparent = num_parent;

  for (int i = 0; i < n_vtx; i++) {
    int i_parent = num_parent[i] - 1;
    for (int j = 0; j < 3; j++) {
      vtx->coords[3*i+j] = _parent->_coords[3*i_parent+j];
    }
  }
  mesh->is_vtx_def_from_parent   = 1;

}


/**
 * \brief  Return number of blocks
 *
 * \param [in]  idx            Nodal mesh handle
 *
 * \return  Number of blocks
 *
 */

int
PDM_Mesh_nodal_n_blocks_get
(
 PDM_Mesh_nodal_t *mesh
)
{
  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }

  return mesh->n_blocks;
}



/**
 * \brief  Return blocks identifier
 *
 * \param [in]  idx            Nodal mesh handle
 *
 * \return  Blocks identifier
 *
 */

int *
PDM_Mesh_nodal_blocks_id_get
(
 PDM_Mesh_nodal_t *mesh
)
{
  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }

  return mesh->blocks_id;

}


/**
 * \brief  Return type of block
 *
 * \param [in]  idx            Nodal mesh handle
 * \param [in]  id_block   Block identifier
 *
 * \return  Type of block
 *
 */

PDM_Mesh_nodal_elt_t
PDM_Mesh_nodal_block_type_get
(
      PDM_Mesh_nodal_t *mesh,
const int               id_block
)
{

  PDM_Mesh_nodal_elt_t t_elt;

  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }

  if (id_block < PDM_BLOCK_ID_BLOCK_POLY2D) {

    t_elt = PDM_MESH_NODAL_POLY_3D;
    const PDM_Mesh_nodal_block_std_t *block = mesh->blocks_std[id_block];

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


/**
 * \brief  Add a new block to the current mesh
 *
 * \param [in]  idx            Nodal mesh handle
 * \param [in]  st_free_data   Status of Release of the memory
 *                             when the block is destroyed
 * \param [in]  id_block       Block identifier
 *
 * \return Block identifier
 *
 */

int
PDM_Mesh_nodal_block_add
(
      PDM_Mesh_nodal_t     *mesh,
      PDM_bool_t            st_free_data,
const PDM_Mesh_nodal_elt_t  t_elt
)
{

  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
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

      mesh->n_block_std++;

      mesh->blocks_std = realloc(mesh->blocks_std, mesh->n_block_std * sizeof(PDM_Mesh_nodal_block_std_t *));

      id_block = mesh->n_block_std-1;

      /* Intialisation du bloc */
      mesh->blocks_std[id_block] = malloc( sizeof(PDM_Mesh_nodal_block_std_t) );
      mesh->blocks_std[id_block]->t_elt        = t_elt;
      mesh->blocks_std[id_block]->st_free_data = st_free_data;
      mesh->blocks_std[id_block]->n_part       = mesh->n_part;

      mesh->blocks_std[id_block]->n_elt        = (PDM_l_num_t  *) malloc(sizeof(PDM_l_num_t  ) * mesh->blocks_std[id_block]->n_part);
      mesh->blocks_std[id_block]->_connec      = (PDM_l_num_t **) malloc(sizeof(PDM_l_num_t *) * mesh->blocks_std[id_block]->n_part);
      mesh->blocks_std[id_block]->_num_part    = (PDM_l_num_t **) malloc(sizeof(PDM_l_num_t *) * mesh->blocks_std[id_block]->n_part);
      mesh->blocks_std[id_block]->_numabs      = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t *) * mesh->blocks_std[id_block]->n_part);
      mesh->blocks_std[id_block]->numabs_int   = NULL;
      mesh->blocks_std[id_block]->_parent_num  = NULL;
      mesh->blocks_std[id_block]->cell_centers = NULL;

      for (int i = 0; i < mesh->blocks_std[id_block]->n_part; i++) {
        mesh->blocks_std[id_block]->n_elt[i]     = 0;
        mesh->blocks_std[id_block]->_connec[i]   = NULL;
        mesh->blocks_std[id_block]->_num_part[i] = NULL;
        mesh->blocks_std[id_block]->_numabs[i]   = NULL;
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

      mesh->n_block_poly2d++;

      mesh->blocks_poly2d = realloc(mesh->blocks_poly2d, mesh->n_block_poly2d * sizeof(PDM_Mesh_nodal_block_poly2d_t *));

      id_block = mesh->n_block_poly2d-1;

      /* Intialisation du bloc */
      mesh->blocks_poly2d[id_block] = malloc( sizeof(PDM_Mesh_nodal_block_poly2d_t) );
      mesh->blocks_poly2d[id_block]->st_free_data = st_free_data;
      mesh->blocks_poly2d[id_block]->n_part            = mesh->n_part;

      mesh->blocks_poly2d[id_block]->n_elt       = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t ) * mesh->blocks_poly2d[id_block]->n_part);
      mesh->blocks_poly2d[id_block]->_connec_idx = (PDM_l_num_t **) malloc(sizeof(PDM_l_num_t *) * mesh->blocks_poly2d[id_block]->n_part);
      mesh->blocks_poly2d[id_block]->_connec     = (PDM_l_num_t **) malloc(sizeof(PDM_l_num_t *) * mesh->blocks_poly2d[id_block]->n_part);
      mesh->blocks_poly2d[id_block]->_num_part   = (PDM_l_num_t **) malloc(sizeof(PDM_l_num_t *) * mesh->blocks_poly2d[id_block]->n_part);
      mesh->blocks_poly2d[id_block]->_numabs     = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t *) * mesh->blocks_poly2d[id_block]->n_part);
      mesh->blocks_poly2d[id_block]->numabs_int = NULL;
      mesh->blocks_poly2d[id_block]->cell_centers = NULL;
      mesh->blocks_poly2d[id_block]->_parent_num = NULL;

      for (int i = 0; i < mesh->blocks_poly2d[id_block]->n_part; i++) {
        mesh->blocks_poly2d[id_block]->n_elt[i]       = 0;
        mesh->blocks_poly2d[id_block]->_connec_idx[i] = NULL;
        mesh->blocks_poly2d[id_block]->_connec[i]     = NULL;
        mesh->blocks_poly2d[id_block]->_num_part[i]   = NULL;
        mesh->blocks_poly2d[id_block]->_numabs[i]     = NULL;
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
      mesh->n_block_poly3d++;

      mesh->blocks_poly3d = realloc(mesh->blocks_poly3d, mesh->n_block_poly3d * sizeof(PDM_Mesh_nodal_block_poly3d_t *));

      id_block = mesh->n_block_poly3d-1;

      /* Intialisation du bloc */

      mesh->blocks_poly3d[id_block] = malloc( sizeof(PDM_Mesh_nodal_block_poly3d_t) );
      mesh->blocks_poly3d[id_block]->n_part       = mesh->n_part;
      mesh->blocks_poly3d[id_block]->st_free_data = st_free_data;

      mesh->blocks_poly3d[id_block]->n_elt        = (PDM_l_num_t * ) malloc(sizeof(PDM_l_num_t  ) * mesh->blocks_poly3d[id_block]->n_part);
      mesh->blocks_poly3d[id_block]->n_face       = (PDM_l_num_t * ) malloc(sizeof(PDM_l_num_t  ) * mesh->blocks_poly3d[id_block]->n_part);
      mesh->blocks_poly3d[id_block]->_facvtx_idx  = (PDM_l_num_t **) malloc(sizeof(PDM_l_num_t *) * mesh->blocks_poly3d[id_block]->n_part);
      mesh->blocks_poly3d[id_block]->_facvtx      = (PDM_l_num_t **) malloc(sizeof(PDM_l_num_t *) * mesh->blocks_poly3d[id_block]->n_part);
      mesh->blocks_poly3d[id_block]->_cellfac_idx = (PDM_l_num_t **) malloc(sizeof(PDM_l_num_t *) * mesh->blocks_poly3d[id_block]->n_part);
      mesh->blocks_poly3d[id_block]->_cellfac     = (PDM_l_num_t **) malloc(sizeof(PDM_l_num_t *) * mesh->blocks_poly3d[id_block]->n_part);
      mesh->blocks_poly3d[id_block]->_cellvtx_idx = (PDM_l_num_t **) malloc(sizeof(PDM_l_num_t *) * mesh->blocks_poly3d[id_block]->n_part);
      mesh->blocks_poly3d[id_block]->_cellvtx     = (PDM_l_num_t **) malloc(sizeof(PDM_l_num_t *) * mesh->blocks_poly3d[id_block]->n_part);
      mesh->blocks_poly3d[id_block]->_numabs      = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t *) * mesh->blocks_poly3d[id_block]->n_part);
      mesh->blocks_poly3d[id_block]->numabs_int   = NULL;
      mesh->blocks_poly3d[id_block]->cell_centers = NULL;
      mesh->blocks_poly3d[id_block]->_parent_num  = NULL;

      for (int i = 0; i < mesh->blocks_poly3d[id_block]->n_part; i++) {
        mesh->blocks_poly3d[id_block]->n_elt[i]        = 0;
        mesh->blocks_poly3d[id_block]->n_face[i]       = 0;
        mesh->blocks_poly3d[id_block]->_facvtx_idx[i]  = NULL;
        mesh->blocks_poly3d[id_block]->_facvtx[i]      = NULL;
        mesh->blocks_poly3d[id_block]->_cellfac_idx[i] = NULL;
        mesh->blocks_poly3d[id_block]->_cellfac[i]     = NULL;
        mesh->blocks_poly3d[id_block]->_cellvtx_idx[i] = NULL;
        mesh->blocks_poly3d[id_block]->_cellvtx[i]     = NULL;
        mesh->blocks_poly3d[id_block]->_numabs[i]      = NULL;
      }

      id_block += PDM_BLOCK_ID_BLOCK_POLY3D;

    }

    break;

  default :
    PDM_error(__FILE__, __LINE__, 0, "Unknown element type\n");
    break;

  }

  _update_blocks_id (mesh);
  return id_block ;

}


/**
 * \brief Define a standard block
 *
 *  - PDM_MESH_NODAL_POINT :
 *
 *   1 x
 *
 *  - PDM_MESH_NODAL_BAR2 :
 *
 *   1 x-------x 2
 *
 *  - PDM_MESH_NODAL_TRIA3 :
 *
 *   1 x-------x 3
 *      \     /
 *       \   /
 *        \ /
 *         x 2
 *
 *  - PDM_MESH_NODAL_QUAD4 :
 *
 *      4 x-------x 3
 *       /       /
 *      /       /
 *   1 x-------x2
 *
 *   - PDM_MESH_NODAL_TETRA4 :
 *
 *         x 4
 *        /|\
 *       / | \
 *      /  |  \
 *   1 x- -|- -x 3
 *      \  |  /
 *       \ | /
 *        \|/
 *         x 2
 *
 *   - PDM_MESH_NODAL_PYRAMID5 :
 *
 *          5 x
 *           /|\
 *          //| \
 *         // |  \
 *      4 x/--|---x 3
 *       //   |  /
 *      //    | /
 *   1 x-------x 2
 *
 *  - PDM_MESH_NODAL_PRSIM6 :
 *
 *   4 x-------x 6
 *     |\     /|
 *     | \   / |
 *   1 x- \-/ -x 3
 *      \ 5x  /
 *       \ | /
 *        \|/
 *         x 2
 *
 *  - PDM_MESH_NODAL_HEXA8 :
 *
 *      8 x-------x 7
 *       /|      /|
 *      / |     / |
 *   5 x-------x6 |
 *     | 4x----|--x 3
 *     | /     | /
 *     |/      |/
 *   1 x-------x 2
 *
 * \param [in]  idx            Nodal mesh handle
 * \param [in]  id_block       Block identifier
 * \param [in]  id_part        Partition identifier
 * \param [in]  n_elt          Number of elements
 * \param [in]  connect        Connectivity
 * \param [in]  numabs         Global numbering
 * \param [in]  parent_num     Parent numbering or NULL
 *
 */

void
PDM_Mesh_nodal_block_std_set
(
      PDM_Mesh_nodal_t *mesh,
const int               id_block,
const int               id_part,
const int               n_elt,
const PDM_l_num_t      *connec,
const PDM_g_num_t      *numabs,
const PDM_l_num_t      *parent_num
)
{

  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }

  int _id_block = id_block - PDM_BLOCK_ID_BLOCK_STD;

  PDM_Mesh_nodal_block_std_t *block = mesh->blocks_std[_id_block];

  if (block == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad standard block identifier\n");
  }

  if (id_part >= block->n_part) {
    PDM_error(__FILE__, __LINE__, 0, "Partition identifier too big\n");
  }

  /* Mapping */

  mesh->n_cell[id_part] += -block->n_elt[id_part];
  mesh->n_cell[id_part] += n_elt;
  block->n_elt[id_part] = n_elt;
  block->_connec[id_part] = (PDM_l_num_t *) connec;
  block->_numabs[id_part] = (PDM_g_num_t *) numabs;

  if (parent_num != NULL) {
    if (block->_parent_num == NULL) {
      block->_parent_num = malloc (sizeof(PDM_l_num_t *) * block->n_part);
      for (int i = 0; i < block->n_part; i++) {
        block->_parent_num[i] = NULL;
      }
    }
    block->_parent_num[id_part] = (PDM_l_num_t *) parent_num;
  }

  /* for (int i = 0; i < n_elt; i++) { */
  /*   mesh->n_elt_abs = PDM_MAX (mesh->n_elt_abs, numabs[i]); */
  /* } */

}


/**
 * \brief Return standard block description
 *
 *  - PDM_MESH_NODAL_POINT :
 *
 *   1 x
 *
 *  - PDM_MESH_NODAL_BAR2 :
 *
 *   1 x-------x 2
 *
 *  - PDM_MESH_NODAL_TRIA3 :
 *
 *   1 x-------x 3
 *      \     /
 *       \   /
 *        \ /
 *         x 2
 *
 *  - PDM_MESH_NODAL_QUAD4 :
 *
 *      4 x-------x 3
 *       /       /
 *      /       /
 *   1 x-------x2
 *
 *   - PDM_MESH_NODAL_TETRA4 :
 *
 *         x 4
 *        /|\
 *       / | \
 *      /  |  \
 *   1 x- -|- -x 3
 *      \  |  /
 *       \ | /
 *        \|/
 *         x 2
 *
 *   - PDM_MESH_NODAL_PYRAMID5 :
 *
 *          5 x
 *           /|\
 *          //| \
 *         // |  \
 *      4 x/--|---x 3
 *       //   |  /
 *      //    | /
 *   1 x-------x 2
 *
 *  - PDM_MESH_NODAL_PRSIM6 :
 *
 *   4 x-------x 6
 *     |\     /|
 *     | \   / |
 *   1 x- \-/ -x 3
 *      \ 5x  /
 *       \ | /
 *        \|/
 *         x 2
 *
 *  - PDM_MESH_NODAL_HEXA8 :
 *
 *      8 x-------x 7
 *       /|      /|
 *      / |     / |
 *   5 x-------x6 |
 *     | 4x----|--x 3
 *     | /     | /
 *     |/      |/
 *   1 x-------x 2
 *
 * \param [in]  idx            Nodal mesh handle
 * \param [in]  id_block       Block identifier
 * \param [in]  id_part        Partition identifier
 * \param [out]  n_elt          Number of elements
 * \param [out]  connect        Connectivity
 * \param [out]  numabs         Global numbering
 * \param [out] numabs_block   Global numbering in the block or NULL (if not computed)
 * \param [out] parent_num     Parent numbering or NULL
 *
 */

void
PDM_Mesh_nodal_block_std_get
(
      PDM_Mesh_nodal_t  *mesh,
const int                id_block,
const int                id_part,
      PDM_l_num_t      **connec
)
{

  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }

  int _id_block = id_block - PDM_BLOCK_ID_BLOCK_STD;

  PDM_Mesh_nodal_block_std_t *block = mesh->blocks_std[_id_block];

  if (block == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad standard block identifier\n");
  }

  if (id_part >= block->n_part) {
    PDM_error(__FILE__, __LINE__, 0, "Partition identifier too big\n");
  }

  *connec = block->_connec[id_part];
}


/**
 * \brief Get number of block elements
 *
 * \param [in]  idx            Nodal mesh handle
 * \param [in]  id_block       Block identifier
 * \param [in]  id_part        Partition identifier
 *
 * \return      Number of elements
 *
 */

int
PDM_Mesh_nodal_block_n_elt_get
(
      PDM_Mesh_nodal_t *mesh,
const int               id_block,
const int               id_part
)
{

  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }

  int _id_block;

  if (id_block >= PDM_BLOCK_ID_BLOCK_POLY3D) {

    _id_block = id_block - PDM_BLOCK_ID_BLOCK_POLY3D;

    PDM_Mesh_nodal_block_poly3d_t *block = mesh->blocks_poly3d[_id_block];

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

    PDM_Mesh_nodal_block_poly2d_t *block = mesh->blocks_poly2d[_id_block];

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

    PDM_Mesh_nodal_block_std_t *block = mesh->blocks_std[_id_block];

    if (block == NULL) {
      PDM_error (__FILE__, __LINE__, 0, "Bad standard block identifier\n");
    }

    if (id_part >= block->n_part) {
      PDM_error(__FILE__, __LINE__, 0, "Partition identifier too big\n");
    }

    return block->n_elt[id_part];
  }

}


/**
 * \brief Get global element numbering of block elements inside the block
 *
 * \param [in]  idx            Nodal mesh handle
 * \param [in]  id_block       Block identifier
 * \param [in]  id_part        Partition identifier
 *
 * \return      Return global numbering of block elements inside the block
 *
 */

PDM_g_num_t *
PDM_Mesh_nodal_block_g_num_get
(
      PDM_Mesh_nodal_t *mesh,
const int               id_block,
const int               id_part
)
{

  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }

  int _id_block;

  if (id_block >= PDM_BLOCK_ID_BLOCK_POLY3D) {

    _id_block = id_block - PDM_BLOCK_ID_BLOCK_POLY3D;

    PDM_Mesh_nodal_block_poly3d_t *block = mesh->blocks_poly3d[_id_block];

    if (block == NULL) {
      PDM_error (__FILE__, __LINE__, 0, "Bad standard block identifier\n");
    }

    if (id_part >= block->n_part) {
      PDM_error(__FILE__, __LINE__, 0, "Partition identifier too big\n");
    }

    return block->numabs_int[id_part];
  }

  else if (id_block >= PDM_BLOCK_ID_BLOCK_POLY2D) {

    _id_block = id_block - PDM_BLOCK_ID_BLOCK_POLY2D;

    PDM_Mesh_nodal_block_poly2d_t *block = mesh->blocks_poly2d[_id_block];

    if (block == NULL) {
      PDM_error (__FILE__, __LINE__, 0, "Bad standard block identifier\n");
    }

    if (id_part >= block->n_part) {
      PDM_error(__FILE__, __LINE__, 0, "Partition identifier too big\n");
    }

    return block->numabs_int[id_part];
  }

  else {

    _id_block = id_block - PDM_BLOCK_ID_BLOCK_STD;

    PDM_Mesh_nodal_block_std_t *block = mesh->blocks_std[_id_block];

    if (block == NULL) {
      PDM_error (__FILE__, __LINE__, 0, "Bad standard block identifier\n");
    }

    if (id_part >= block->n_part) {
      PDM_error(__FILE__, __LINE__, 0, "Partition identifier too big\n");
    }
    if (block->numabs_int == NULL) printf("!!! id_block = %d\n", id_block);
    return block->numabs_int[id_part];
  }
}


/**
 * \brief Get global element numbering of block elements
 *
 * \param [in]  idx            Nodal mesh handle
 * \param [in]  id_block       Block identifier
 * \param [in]  id_part        Partition identifier
 *
 * \return      Return global element numbering of block elements
 *
 */

PDM_g_num_t *
PDM_Mesh_nodal_g_num_get
(
      PDM_Mesh_nodal_t *mesh,
const int               id_block,
const int               id_part
)
{

  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }

  int _id_block;

  PDM_g_num_t *_numabs = NULL;

  if (id_block >= PDM_BLOCK_ID_BLOCK_POLY3D) {

    _id_block = id_block - PDM_BLOCK_ID_BLOCK_POLY3D;

    PDM_Mesh_nodal_block_poly3d_t *block = mesh->blocks_poly3d[_id_block];

    if (block == NULL) {
      PDM_error (__FILE__, __LINE__, 0, "Bad standard block identifier\n");
    }

    if (id_part >= block->n_part) {
      PDM_error(__FILE__, __LINE__, 0, "Partition identifier too big\n");
    }

    if (block->_numabs != NULL) {
      _numabs = block->_numabs[id_part];
    }
  }

  else if (id_block >= PDM_BLOCK_ID_BLOCK_POLY2D) {

    _id_block = id_block - PDM_BLOCK_ID_BLOCK_POLY2D;

    PDM_Mesh_nodal_block_poly2d_t *block = mesh->blocks_poly2d[_id_block];

    if (block == NULL) {
      PDM_error (__FILE__, __LINE__, 0, "Bad standard block identifier\n");
    }

    if (id_part >= block->n_part) {
      PDM_error(__FILE__, __LINE__, 0, "Partition identifier too big\n");
    }

    if (block->_numabs != NULL) {
      _numabs = block->_numabs[id_part];
    }
  }

  else {

    _id_block = id_block - PDM_BLOCK_ID_BLOCK_STD;

    PDM_Mesh_nodal_block_std_t *block = mesh->blocks_std[_id_block];

    if (block == NULL) {
      PDM_error (__FILE__, __LINE__, 0, "Bad standard block identifier\n");
    }

    if (id_part >= block->n_part) {
      PDM_error(__FILE__, __LINE__, 0, "Partition identifier too big\n");
    }

    if (block->_numabs != NULL) {
      _numabs = block->_numabs[id_part];
    }
  }

  return _numabs;
}


/**
 * \brief Get parent numbering of block elements
 *
 * \param [in]  idx            Nodal mesh handle
 * \param [in]  id_block       Block identifier
 * \param [in]  id_part        Partition identifier
 *
 * \return      Return parent numbering of block elements
 *
 */

int *
PDM_Mesh_nodal_block_parent_num_get
(
      PDM_Mesh_nodal_t *mesh,
const int               id_block,
const int               id_part
)
{

  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }

  int _id_block;

  int *_parent_num = NULL;

  if (id_block >= PDM_BLOCK_ID_BLOCK_POLY3D) {

    _id_block = id_block - PDM_BLOCK_ID_BLOCK_POLY3D;

    PDM_Mesh_nodal_block_poly3d_t *block = mesh->blocks_poly3d[_id_block];

    if (block == NULL) {
      PDM_error (__FILE__, __LINE__, 0, "Bad standard block identifier\n");
    }

    if (id_part >= block->n_part) {
      PDM_error(__FILE__, __LINE__, 0, "Partition identifier too big\n");
    }

    if (block->_parent_num != NULL) {
      _parent_num = block->_parent_num[id_part];
    }
  }

  else if (id_block >= PDM_BLOCK_ID_BLOCK_POLY2D) {

    _id_block = id_block - PDM_BLOCK_ID_BLOCK_POLY2D;

    PDM_Mesh_nodal_block_poly2d_t *block = mesh->blocks_poly2d[_id_block];

    if (block == NULL) {
      PDM_error (__FILE__, __LINE__, 0, "Bad standard block identifier\n");
    }

    if (id_part >= block->n_part) {
      PDM_error(__FILE__, __LINE__, 0, "Partition identifier too big\n");
    }

    if (block->_parent_num != NULL) {
      _parent_num = block->_parent_num[id_part];
    }
  }

  else {

    _id_block = id_block - PDM_BLOCK_ID_BLOCK_STD;

    PDM_Mesh_nodal_block_std_t *block = mesh->blocks_std[_id_block];

    if (block == NULL) {
      PDM_error (__FILE__, __LINE__, 0, "Bad standard block identifier\n");
    }

    if (id_part >= block->n_part) {
      PDM_error(__FILE__, __LINE__, 0, "Partition identifier too big\n");
    }

    if (block->_parent_num != NULL) {
      _parent_num = block->_parent_num[id_part];
    }
  }

  return _parent_num;
}


/**
 * \brief Define a polygon block
 *
 * \param [in]  idx            Nodal mesh handle
 * \param [in]  id_block       Block identifier
 * \param [in]  id_part        Partition identifier
 * \param [in]  n_elt          Number of elements
 * \param [in]  connect_idx    Connectivity index (size = \ref n_elt + 1)
 * \param [in]  connect        Connectivity (size = \ref connect_idx[\ref n_elt])
 * \param [in]  numabs         Global numbering
 * \param [in]  parent_num     Parent numbering or NULL
 *
 */

void
PDM_Mesh_nodal_block_poly2d_set
(
      PDM_Mesh_nodal_t *mesh,
const int               id_block,
const int               id_part,
const PDM_l_num_t       n_elt,
const PDM_l_num_t      *connec_idx,
const PDM_l_num_t      *connec,
const PDM_g_num_t      *numabs,
const PDM_l_num_t      *parent_num
)
{

  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }

  int _id_block = id_block - PDM_BLOCK_ID_BLOCK_POLY2D;

  PDM_Mesh_nodal_block_poly2d_t *block = mesh->blocks_poly2d[_id_block];

  if (block == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad standard block identifier\n");
  }

  if (id_part >= block->n_part) {
    PDM_error(__FILE__, __LINE__, 0, "Partition identifier too big\n");
  }

  /* Mapping */

  mesh->n_cell[id_part]      += -block->n_elt[id_part];
  mesh->n_cell[id_part]      += n_elt;
  block->n_elt[id_part]       = n_elt;
  block->_connec_idx[id_part] = (PDM_l_num_t *) connec_idx;
  block->_connec[id_part]     = (PDM_l_num_t *) connec;
  block->_numabs[id_part]     = (PDM_g_num_t *) numabs;

  /* for (int i = 0; i < n_elt; i++) { */
  /*   mesh->n_elt_abs = PDM_MAX(mesh->n_elt_abs, numabs[i]); */
  /* } */

  if (parent_num != NULL) {
    if (block->_parent_num == NULL) {
      block->_parent_num = malloc (sizeof(PDM_l_num_t *) * block->n_part);
      for (int i = 0; i < block->n_part; i++) {
        block->_parent_num[i] = NULL;
      }
    }
    block->_parent_num[id_part] = (PDM_l_num_t *) parent_num;
  }

}



/**
 * \brief Return a polygon block description
 *
 * \param [in]  idx            Nodal mesh handle
 * \param [in]  id_block       Block identifier
 * \param [in]  id_part        Partition identifier
 * \param [out] connect_idx    Connectivity index (size = \ref n_elt + 1)
 * \param [out] connect        Connectivity (size = \ref connect_idx[\ref n_elt])
 *
 */

void
PDM_Mesh_nodal_block_poly2d_get
(
      PDM_Mesh_nodal_t  *mesh,
 const int               id_block,
 const int               id_part,
       PDM_l_num_t     **connec_idx,
       PDM_l_num_t     **connec
)
{

  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }

  int _id_block = id_block - PDM_BLOCK_ID_BLOCK_POLY2D;

  PDM_Mesh_nodal_block_poly2d_t *block = mesh->blocks_poly2d[_id_block];

  if (block == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad standard block identifier\n");
  }

  if (id_part >= block->n_part) {
    PDM_error(__FILE__, __LINE__, 0, "Partition identifier too big\n");
  }

  *connec_idx = block->_connec_idx[id_part];
  *connec     = block->_connec[id_part];

}



static int
_binary_search
(
 const PDM_l_num_t  elem,
 const PDM_l_num_t  array[],
 const PDM_l_num_t  n,
 PDM_bool_t        *in_array
 )
{
  int l = 0;
  int r = n;

  *in_array = PDM_FALSE;

  if (n < 1)
    return 0;

  while (l + 1 < r) {
    int m = l + (r - l)/2;

    if (elem < array[m])
      r = m;
    else
      l = m;
  }



  if (array[l] == elem) {
    *in_array = PDM_TRUE;
    return l;

  }

  else if (array[l] < elem)
    return l + 1;

  else
    return l;
}


static void _compute_cell_vtx_connectivity
(
 const PDM_l_num_t   n_cell,
 const PDM_l_num_t   n_face,
 const PDM_l_num_t  *face_vtx_idx,
 const PDM_l_num_t  *face_vtx,
 const PDM_l_num_t  *cell_face_idx,
 const PDM_l_num_t  *cell_face,
 PDM_l_num_t       **cell_vtx_idx,
 PDM_l_num_t       **cell_vtx
 )
{
  PDM_UNUSED(n_face);

  const int DEBUG = 0;

  *cell_vtx_idx = malloc (sizeof(int) * (n_cell + 1));
  PDM_l_num_t *_cell_vtx_idx = *cell_vtx_idx;

  _cell_vtx_idx[0] = 0;

  size_t s_cell_vtx = 10 * n_cell;
  *cell_vtx = malloc (sizeof(PDM_l_num_t) * s_cell_vtx);

  PDM_bool_t already_in_cell;
  int pos, i;
  PDM_l_num_t icell, iface, ivtx, id_face, id_vtx;

  /* Loop on cells */
  PDM_l_num_t n_vtx_cell;
  for (icell = 0; icell < n_cell; icell++) {

    PDM_l_num_t *_cell_vtx = *cell_vtx + _cell_vtx_idx[icell];
    n_vtx_cell = 0;

    /* Loop on current cell's faces */
    for (iface = cell_face_idx[icell]; iface < cell_face_idx[icell+1]; iface++) {
      id_face = PDM_ABS (cell_face[iface]) - 1;

      /* Loop on current face's vertices */
      for (ivtx = face_vtx_idx[id_face]; ivtx < face_vtx_idx[id_face+1]; ivtx++) {
        id_vtx = face_vtx[ivtx];

        pos = _binary_search (id_vtx,
                              _cell_vtx,
                              n_vtx_cell,
                             &already_in_cell);

        if (already_in_cell == PDM_TRUE) {
          continue;
        }

        if (n_vtx_cell + _cell_vtx_idx[icell] >= (int) s_cell_vtx) {
          s_cell_vtx = PDM_MAX ((int) (2*s_cell_vtx), n_vtx_cell + _cell_vtx_idx[icell]);
          *cell_vtx = realloc (*cell_vtx, sizeof(PDM_l_num_t) * s_cell_vtx);
          _cell_vtx = *cell_vtx + _cell_vtx_idx[icell];
        }

        for (i = n_vtx_cell; i > pos; i--) {
          _cell_vtx[i] = _cell_vtx[i-1];
        }
        _cell_vtx[pos] = id_vtx;
        n_vtx_cell++;

      } // End of loop on current face's vertices

    } // End of loop on current cell's faces

    _cell_vtx_idx[icell+1] = _cell_vtx_idx[icell] + n_vtx_cell;

    if (DEBUG) {
      printf("cell #%d vtx =", icell);
      for (int j = _cell_vtx_idx[icell]; j < _cell_vtx_idx[icell+1]; j++) {
        printf(" %d", (*cell_vtx)[j]);
      }
      printf("\n");
    }

  } // End of loop on cells

  *cell_vtx = realloc (*cell_vtx, sizeof(PDM_l_num_t) * _cell_vtx_idx[n_cell]);
}



/**
 * \brief Define a polyhedra block
 *
 * \param [in]  idx            Nodal mesh handle
 * \param [in]  id_block       Block identifier
 * \param [in]  id_part        Partition identifier
 * \param [in]  n_elt          Number of polyhedra
 * \param [in]  n_face         Number of faces used to describe polyhedra
 * \param [in]  facvtx_idx     Index of face vertex connectivity
 * \param [in]  facvtx         Face vertex connectivity
 * \param [in]  cellfac_idx    Index of cell face connectivity
 * \param [in]  cellfac        Cell face connectivity
 * \param [in]  numabs         Global numbering
 * \param [in]  parent_num     Parent numbering or NULL
 *
 */

void
PDM_Mesh_nodal_block_poly3d_set
(
      PDM_Mesh_nodal_t *mesh,
const int               id_block,
const int               id_part,
const PDM_l_num_t       n_elt,
const PDM_l_num_t       n_face,
const PDM_l_num_t      *facvtx_idx,
const PDM_l_num_t      *facvtx,
const PDM_l_num_t      *cellfac_idx,
const PDM_l_num_t      *cellfac,
const PDM_g_num_t      *numabs,
const PDM_l_num_t      *parent_num
)
{

  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }

  int _id_block = id_block - PDM_BLOCK_ID_BLOCK_POLY3D;


  PDM_Mesh_nodal_block_poly3d_t *block = mesh->blocks_poly3d[_id_block];

  if (block == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad standard block identifier\n");
  }

  if (id_part >= block->n_part) {
    PDM_error(__FILE__, __LINE__, 0, "Partition identifier too big\n");
  }

  mesh->n_cell[id_part]       += -block->n_elt[id_part];
  mesh->n_cell[id_part]       += n_elt;
  block->n_elt[id_part]        = n_elt;
  block->n_face[id_part]       = n_face;
  block->_facvtx_idx[id_part]  = (PDM_l_num_t *) facvtx_idx;
  block->_facvtx[id_part]      = (PDM_l_num_t *) facvtx;
  block->_cellfac_idx[id_part] = (PDM_l_num_t *) cellfac_idx;
  block->_cellfac[id_part]     = (PDM_l_num_t *) cellfac;
  block->_numabs[id_part]      = (PDM_g_num_t *) numabs;


  /* Compute cell-vertex connectivity */
  _compute_cell_vtx_connectivity (n_elt,
                                  n_face,
                                  facvtx_idx,
                                  facvtx,
                                  cellfac_idx,
                                  cellfac,
                                  &(block->_cellvtx_idx[id_part]),
                                  &(block->_cellvtx[id_part]));

  /* for (int i = 0; i < n_elt; i++) { */
  /*   mesh->n_elt_abs = PDM_MAX (mesh->n_elt_abs, numabs[i]); */
  /* } */

  if (parent_num != NULL) {
    if (block->_parent_num == NULL) {
      block->_parent_num = malloc (sizeof(PDM_l_num_t *) * block->n_part);
      for (int i = 0; i < block->n_part; i++) {
        block->_parent_num[i] = NULL;
      }
    }
    block->_parent_num[id_part] = (PDM_l_num_t *) parent_num;
  }

}


/**
 * \brief Define a polyhedra block
 *
 * \param [in]  idx            Nodal mesh handle
 * \param [in]  id_block       Block identifier
 * \param [in]  id_part        Partition identifier
 * \param [out]  n_face         Number of faces used to describe polyhedra
 * \param [out]  facvtx_idx     Index of face vertex connectivity
 * \param [out]  facvtx         Face vertex connectivity
 * \param [out]  cellfac_idx    Index of cell face connectivity
 * \param [out]  cellfac        Cell face connectivity
 *
 */

void
PDM_Mesh_nodal_block_poly3d_get
(
      PDM_Mesh_nodal_t  *mesh,
const int                id_block,
const int                id_part,
      PDM_l_num_t       *n_face,
      PDM_l_num_t      **facvtx_idx,
      PDM_l_num_t      **facvtx,
      PDM_l_num_t      **cellfac_idx,
      PDM_l_num_t      **cellfac
)
{

  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }

  int _id_block = id_block - PDM_BLOCK_ID_BLOCK_POLY3D;


  PDM_Mesh_nodal_block_poly3d_t *block = mesh->blocks_poly3d[_id_block];

  if (block == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad standard block identifier\n");
  }

  if (id_part >= block->n_part) {
    PDM_error(__FILE__, __LINE__, 0, "Partition identifier too big\n");
  }

  *n_face      = block->n_face[id_part];
  *facvtx_idx  = block->_facvtx_idx[id_part];
  *facvtx      = block->_facvtx[id_part];
  *cellfac_idx = block->_cellfac_idx[id_part];
  *cellfac     = block->_cellfac[id_part];

}


/**
 * \brief Get the cell-vertex connectivity of a polyhedra block
 *
 * \param [in]  idx            Nodal mesh handle
 * \param [in]  id_block       Block identifier
 * \param [in]  id_part        Partition identifier
 * \param [out] cellvtx_idx    Index of cell vertex connectivity
 * \param [out] cellvtx        Cell vertex connectivity
 *
 */

void
PDM_Mesh_nodal_block_poly3d_cell_vtx_connect_get
(
      PDM_Mesh_nodal_t  *mesh,
 const int               id_block,
 const int               id_part,
 PDM_l_num_t           **cellvtx_idx,
 PDM_l_num_t           **cellvtx
)
{

  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }

  int _id_block = id_block - PDM_BLOCK_ID_BLOCK_POLY3D;


  PDM_Mesh_nodal_block_poly3d_t *block = mesh->blocks_poly3d[_id_block];

  if (block == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad standard block identifier\n");
  }

  if (id_part >= block->n_part) {
    PDM_error(__FILE__, __LINE__, 0, "Partition identifier too big\n");
  }

  *cellvtx_idx = block->_cellvtx_idx[id_part];
  *cellvtx     = block->_cellvtx[id_part];

}


/**
 * \brief  Add some 3D cells from cell face conectivity.
 *
 * For each cell, this function searchs the type of the cell (tetrahedra, hexahedra, ...)
 * and stores it in the corresponding block. \ref ind_num gives the indirection
 * between old and new numbering.
 *
 * \param [in]  idx            Nodal mesh handle
 * \param [in]  id_part        Partition identifier
 * \param [in]  n_elt          Number of polyhedra
 * \param [in]  n_face         Number of faces used to describe polyhedra
 * \param [in]  face_vtx_idx   Index of face vertex connectivity
 * \param [in]  face_vtx_nb    Number of vertices for each face
 * \param [in]  face_vtx       Face vertex connectivity
 * \param [in]  cell_face_idx  Index of cell face connectivity
 * \param [in]  cell_face_nb   Number of faces for each cell
 * \param [in]  cell_face      Cell face connectivity
 * \param [in]  numabs         Global numbering
 * \param [out] ind_num        old numbering to new numbering
 *
 */

void
PDM_Mesh_nodal_cell3d_cellface_add
(
      PDM_Mesh_nodal_t *mesh,
const int               id_part,
const int               n_cell,
const int               n_face,
const PDM_l_num_t      *face_vtx_idx,
const PDM_l_num_t      *face_vtx_nb,
const PDM_l_num_t      *face_vtx,
const PDM_l_num_t      *cell_face_idx,
const PDM_l_num_t      *cell_face_nb,
const PDM_l_num_t      *cell_face,
const PDM_g_num_t      *numabs
)
{

  int adjust = 0;
  if (n_cell > 0) {
    if (cell_face_idx[0] == 1) {
      adjust = 1;
    }
  }

  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }

  if (id_part >= mesh->n_part) {
    PDM_error(__FILE__, __LINE__, 0, "Partition identifier too big\n");
  }

  int n_part = 0;

  if (mesh->num_cell_parent_to_local == NULL) {
    mesh->num_cell_parent_to_local = (PDM_l_num_t **) malloc(sizeof(PDM_l_num_t *) * mesh->n_part);
    for (int i_part = 0; i_part < mesh->n_part; i_part++) {
      mesh->num_cell_parent_to_local[i_part] = NULL;
    }
  }

  mesh->num_cell_parent_to_local[id_part] = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t) * n_cell);
  for (int i = 0; i < n_cell; i++) {
    mesh->num_cell_parent_to_local[id_part][i] = 0;
  }

  if (mesh->prepa_blocks == NULL) {
    mesh->prepa_blocks = (PDM_Mesh_nodal_prepa_blocks_t *) malloc(sizeof(PDM_Mesh_nodal_prepa_blocks_t));
    mesh->prepa_blocks->t_add = 1;
    mesh->prepa_blocks->n_tria_proc    = 0;  /* Nb de triangles par proc */
    mesh->prepa_blocks->n_quad_proc    = 0;  /* Nb de quads par proc     */
    mesh->prepa_blocks->n_poly2d_proc  = 0;  /* Nb de poly2d par proc    */
    mesh->prepa_blocks->n_tetra_proc   = 0;  /* Nb de tetra par proc     */
    mesh->prepa_blocks->n_hexa_proc    = 0;  /* Nb d'hexa par proc       */
    mesh->prepa_blocks->n_prism_proc   = 0;  /* Nb de prisme par proc    */
    mesh->prepa_blocks->n_pyramid_proc = 0;  /* Nb de pyramide par proc  */
    mesh->prepa_blocks->n_poly3d_proc  = 0;  /* Nb de poly3d par proc    */

    mesh->prepa_blocks->n_cell        = (PDM_l_num_t  *) malloc(sizeof(PDM_l_num_t  ) * mesh->n_part);
    mesh->prepa_blocks->n_face        = (PDM_l_num_t  *) malloc(sizeof(PDM_l_num_t  ) * mesh->n_part);
    mesh->prepa_blocks->n_tetra       = (PDM_l_num_t  *) malloc(sizeof(PDM_l_num_t  ) * mesh->n_part);
    mesh->prepa_blocks->n_hexa        = (PDM_l_num_t  *) malloc(sizeof(PDM_l_num_t  ) * mesh->n_part);
    mesh->prepa_blocks->n_prism       = (PDM_l_num_t  *) malloc(sizeof(PDM_l_num_t  ) * mesh->n_part);
    mesh->prepa_blocks->n_pyramid     = (PDM_l_num_t  *) malloc(sizeof(PDM_l_num_t  ) * mesh->n_part);
    mesh->prepa_blocks->n_poly3d      = (PDM_l_num_t  *) malloc(sizeof(PDM_l_num_t  ) * mesh->n_part);
    mesh->prepa_blocks->face_vtx_idx  = (PDM_l_num_t **) malloc(sizeof(PDM_l_num_t *) * mesh->n_part);
    mesh->prepa_blocks->face_vtx_nb   = (PDM_l_num_t **) malloc(sizeof(PDM_l_num_t *) * mesh->n_part);
    mesh->prepa_blocks->face_vtx      = (PDM_l_num_t **) malloc(sizeof(PDM_l_num_t *) * mesh->n_part);
    mesh->prepa_blocks->cell_face_idx = (PDM_l_num_t **) malloc(sizeof(PDM_l_num_t *) * mesh->n_part);
    mesh->prepa_blocks->cell_face_nb  = (PDM_l_num_t **) malloc(sizeof(PDM_l_num_t *) * mesh->n_part);
    mesh->prepa_blocks->cell_face     = (PDM_l_num_t **) malloc(sizeof(PDM_l_num_t *) * mesh->n_part);
    mesh->prepa_blocks->add_etat      = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t)*mesh->n_part);
    mesh->prepa_blocks->numabs = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t *)*mesh->n_part);
    for (int i = 0; i < mesh->n_part; i++) {
      mesh->prepa_blocks->add_etat[i] = 0;
    }
  }

  if (mesh->prepa_blocks->t_add != 1) {
    PDM_error(__FILE__, __LINE__, 0, "Erreur Cs_geom_cell3d_cell_face_add : Un autre type d'ajout est en cours\n");
    abort();
  }

  /* Determination du type de chaque element */

  PDM_l_num_t cell_som_tria[18]; /* 6 triangles max in _type_cell_3D */
  PDM_l_num_t cell_som_quad[24]; /* 6 quadrangles max in _type_cell_3D */
  PDM_l_num_t n_tetra   = 0;
  PDM_l_num_t n_hexa    = 0;
  PDM_l_num_t n_prism   = 0;
  PDM_l_num_t n_pyramid = 0;
  PDM_l_num_t n_poly3d  = 0;

  if (0 == 1) {
    printf("1 cell_face %d %d : \n",id_part, n_cell);
    for (int i = 0; i < n_cell; i++) {
      for (int j = cell_face_idx[i] -adjust; j < cell_face_idx[i] -adjust + cell_face_nb[i]; j++) {
        printf(" %d", cell_face[j]);
      }
      printf("\n");
    }

    printf("1 face_vtx %d %d: \n", id_part, n_face);
    for (int i = 0; i < n_face; i++) {
      for (int j = face_vtx_idx[i] -adjust; j < face_vtx_idx[i] -adjust + face_vtx_nb[i] ; j++) {
        printf(" %d", face_vtx[j]);
      }
      printf("\n");
    }
  }

  for (int i = 0; i < n_cell; i++) {

    PDM_Mesh_nodal_elt_t cell_type = _type_cell_3D(cell_face_nb[i],
                                                   cell_face + cell_face_idx[i] - adjust,
                                                   face_vtx_idx,
                                                   face_vtx_nb,
                                                   face_vtx,
                                                   cell_som_tria,
                                                   cell_som_quad);
    switch(cell_type) {
    case PDM_MESH_NODAL_TETRA4 :
      n_tetra += 1;
      break;
    case PDM_MESH_NODAL_PYRAMID5 :
      n_pyramid += 1;
      break;
    case PDM_MESH_NODAL_PRISM6 :
      n_prism += 1;
      break;
    case PDM_MESH_NODAL_HEXA8 :
      n_hexa += 1;
      break;
    case PDM_MESH_NODAL_POLY_3D :
      n_poly3d += 1;
      break;
    default :
      break;
    }
  }

  mesh->prepa_blocks->n_tetra_proc          += n_tetra;
  mesh->prepa_blocks->n_hexa_proc           += n_hexa;
  mesh->prepa_blocks->n_prism_proc          += n_prism;
  mesh->prepa_blocks->n_pyramid_proc        += n_pyramid;
  mesh->prepa_blocks->n_poly3d_proc         += n_poly3d;
  mesh->prepa_blocks->n_tetra[id_part]       = n_tetra;
  mesh->prepa_blocks->n_hexa[id_part]        = n_hexa;
  mesh->prepa_blocks->n_prism[id_part]       = n_prism;
  mesh->prepa_blocks->n_pyramid[id_part]     = n_pyramid;
  mesh->prepa_blocks->n_poly3d[id_part]      = n_poly3d;
  mesh->prepa_blocks->face_vtx_idx[id_part]  = (PDM_l_num_t *) face_vtx_idx;
  mesh->prepa_blocks->face_vtx_nb[id_part]   = (PDM_l_num_t *) face_vtx_nb;
  mesh->prepa_blocks->face_vtx[id_part]      = (PDM_l_num_t *) face_vtx;
  mesh->prepa_blocks->cell_face_idx[id_part] = (PDM_l_num_t *) cell_face_idx;
  mesh->prepa_blocks->cell_face_nb[id_part]  = (PDM_l_num_t *) cell_face_nb;
  mesh->prepa_blocks->cell_face[id_part]     = (PDM_l_num_t *) cell_face;
  mesh->prepa_blocks->numabs[id_part]        = (PDM_g_num_t *) numabs;
  mesh->prepa_blocks->add_etat[id_part]      = 1;
  mesh->prepa_blocks->n_face[id_part]        = n_face;
  mesh->prepa_blocks->n_cell[id_part]        = n_cell;

  /* Creation des blocs si toutes les parts sont remplies */

  for (int i = 0; i < mesh->n_part; i++) {
    if (mesh->prepa_blocks->add_etat[i] == 1)
      n_part += 1;
  }

  if (mesh->n_part == n_part) {

    /* Creation des blocs */

    PDM_l_num_t elts[5];
    PDM_l_num_t som_elts[5];

    elts[0] = mesh->prepa_blocks->n_tetra_proc   > 0;
    elts[1] = mesh->prepa_blocks->n_hexa_proc    > 0;
    elts[2] = mesh->prepa_blocks->n_prism_proc   > 0;
    elts[3] = mesh->prepa_blocks->n_pyramid_proc > 0;
    elts[4] = mesh->prepa_blocks->n_poly3d_proc  > 0;

    PDM_MPI_Allreduce(elts, som_elts, 5, PDM_MPI_INT, PDM_MPI_SUM, mesh->pdm_mpi_comm);

    int id_bloc_tetra4   = -1;
    int id_bloc_hexa8    = -1;
    int id_bloc_prism6   = -1;
    int id_bloc_pyramid5 = -1;
    int id_bloc_poly_3d  = -1;

    if (som_elts[0] > 0)
      id_bloc_tetra4 = PDM_Mesh_nodal_block_add(mesh,
                                                PDM_TRUE,
                                                PDM_MESH_NODAL_TETRA4);

    if (som_elts[1] > 0)
      id_bloc_hexa8 = PDM_Mesh_nodal_block_add(mesh,
                                               PDM_TRUE,
                                               PDM_MESH_NODAL_HEXA8);

    if (som_elts[2] > 0)
      id_bloc_prism6 = PDM_Mesh_nodal_block_add(mesh,
                                                PDM_TRUE,
                                                PDM_MESH_NODAL_PRISM6);

    if (som_elts[3] > 0)
      id_bloc_pyramid5 = PDM_Mesh_nodal_block_add(mesh,
                                                  PDM_TRUE,
                                                  PDM_MESH_NODAL_PYRAMID5);

    if (som_elts[4] > 0)
      id_bloc_poly_3d = PDM_Mesh_nodal_block_add(mesh,
                                                 PDM_TRUE,
                                                 PDM_MESH_NODAL_POLY_3D);

    /* Determination de la connectivite de chaque element */


    for (int i_part = 0; i_part < mesh->n_part; i_part++) {

      PDM_l_num_t n_cell_courant = mesh->prepa_blocks->n_cell[i_part];
      PDM_l_num_t *num_cell_parent_to_local_courant = mesh->num_cell_parent_to_local[i_part];
      PDM_l_num_t *face_som_idx_courant = mesh->prepa_blocks->face_vtx_idx[i_part];
      PDM_l_num_t *face_som_nb_courant = mesh->prepa_blocks->face_vtx_nb[i_part];
      PDM_l_num_t *face_som_courant = mesh->prepa_blocks->face_vtx[i_part];
      PDM_l_num_t *cell_face_idx_courant = mesh->prepa_blocks->cell_face_idx[i_part];
      PDM_l_num_t *cell_face_nb_courant = mesh->prepa_blocks->cell_face_nb[i_part];
      PDM_l_num_t *cell_face_courant = mesh->prepa_blocks->cell_face[i_part];
      PDM_g_num_t *numabs_courant = mesh->prepa_blocks->numabs[i_part];
      PDM_l_num_t n_face_part   = mesh->prepa_blocks->n_face[i_part];

      PDM_l_num_t n_tetra_part   = mesh->prepa_blocks->n_tetra  [i_part];
      PDM_l_num_t n_hexa_part    = mesh->prepa_blocks->n_hexa   [i_part];
      PDM_l_num_t n_prism_part   = mesh->prepa_blocks->n_prism  [i_part];
      PDM_l_num_t n_pyramid_part = mesh->prepa_blocks->n_pyramid[i_part];
      PDM_l_num_t n_poly3d_part  = mesh->prepa_blocks->n_poly3d [i_part];

      PDM_l_num_t *connec_tetra   = NULL;
      PDM_l_num_t *connec_hexa    = NULL;
      PDM_l_num_t *connec_prism   = NULL;
      PDM_l_num_t *connec_pyramid = NULL;

      PDM_g_num_t *numabs_tetra   = NULL;
      PDM_g_num_t *numabs_hexa    = NULL;
      PDM_g_num_t *numabs_prism   = NULL;
      PDM_g_num_t *numabs_pyramid = NULL;
      PDM_g_num_t *numabs_poly3d  = NULL;

      PDM_l_num_t *num_parent_tetra   = NULL;
      PDM_l_num_t *num_parent_hexa    = NULL;
      PDM_l_num_t *num_parent_prism   = NULL;
      PDM_l_num_t *num_parent_pyramid = NULL;
      PDM_l_num_t *num_parent_poly3d  = NULL;

      adjust = 0;
      if (n_cell_courant > 0) {
        if (cell_face_idx_courant[0] == 1) {
          adjust = 1;
        }
      }

      if (0 == 1) {
        printf("2 cell_face %d %d: \n",i_part, n_cell_courant);
        for (int i = 0; i < n_cell_courant; i++) {
          for (int j = cell_face_idx_courant[i] -adjust; j < cell_face_idx_courant[i]  -adjust+ cell_face_nb_courant[i]; j++) {
            printf(" %d", cell_face_courant[j]);
          }
          printf("\n");
        }

        printf("2 face_vtx %d %d: \n", i_part, n_face_part);
        for (int i = 0; i < n_face_part; i++) {
          for (int j = face_som_idx_courant[i] -adjust; j < face_som_idx_courant[i] -adjust + face_som_nb_courant[i] ; j++) {
            printf(" %d", face_som_courant[j]);
          }
          printf("\n");
        }
      }

//      if (n_tetra_part > 0) {
      if (som_elts[0] > 0) {
        connec_tetra = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t) * 4 *n_tetra_part);
        numabs_tetra = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * n_tetra_part);
        num_parent_tetra = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t) * n_tetra_part);
      }

//      if (n_hexa_part > 0) {
      if (som_elts[1] > 0) {
        connec_hexa = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t) * 8 * n_hexa_part);
        numabs_hexa = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * n_hexa_part);
        num_parent_hexa = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t) * n_hexa_part);
      }

//      if (n_prism_part > 0) {
      if (som_elts[2] > 0) {
        connec_prism = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t) * 6 * n_prism_part);
        numabs_prism = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * n_prism_part);
        num_parent_prism = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t) * n_prism_part);
      }

//      if (n_pyramid_part > 0) {
      if (som_elts[3] > 0) {
        connec_pyramid = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t) * 5 * n_pyramid_part);
        numabs_pyramid = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * n_pyramid_part);
        num_parent_pyramid = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t) * n_pyramid_part);
      }

//      if (n_poly3d_part > 0) {
      if (som_elts[4] > 0) {
        numabs_poly3d = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * n_poly3d_part);
        num_parent_poly3d = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t) * n_poly3d_part);
      }

      PDM_l_num_t *num_parent_tetra_courant = num_parent_tetra;
      PDM_l_num_t *num_parent_hexa_courant = num_parent_hexa;
      PDM_l_num_t *num_parent_prism_courant = num_parent_prism;
      PDM_l_num_t *num_parent_pyramid_courant = num_parent_pyramid;
      PDM_l_num_t *num_parent_poly3d_courant = num_parent_poly3d;

      PDM_l_num_t *connec_tetra_courant = connec_tetra;
      PDM_l_num_t *connec_hexa_courant = connec_hexa;
      PDM_l_num_t *connec_prism_courant = connec_prism;
      PDM_l_num_t *connec_pyramid_courant = connec_pyramid;

      PDM_g_num_t *numabs_tetra_courant = numabs_tetra;
      PDM_g_num_t *numabs_hexa_courant = numabs_hexa;
      PDM_g_num_t *numabs_prism_courant = numabs_prism;
      PDM_g_num_t *numabs_pyramid_courant = numabs_pyramid;
      PDM_g_num_t *numabs_poly3d_courant = numabs_poly3d;

      PDM_l_num_t *tag_face_poly3d = NULL;
      PDM_l_num_t  n_face_poly = 0;
      PDM_l_num_t *facsom_poly_idx = NULL;
      PDM_l_num_t *facsom_poly = NULL;
      PDM_l_num_t *cellfac_poly_idx = NULL;
      PDM_l_num_t *cellfac_poly = NULL;
      PDM_l_num_t l_cellfac_poly = 0;

      if (n_poly3d_part > 0) {
        tag_face_poly3d = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t) * n_face_part);
        for (int i = 0; i < n_face_part; i++) {
          tag_face_poly3d[i] = -1;
        }
        cellfac_poly_idx = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t) * (n_poly3d_part + 1));
        cellfac_poly_idx[0] = 0;
      }

      PDM_l_num_t idx_tetra = 0;
      PDM_l_num_t idx_hexa = n_tetra_part;
      PDM_l_num_t idx_prism = idx_hexa + n_hexa_part;
      PDM_l_num_t idx_pyramid = idx_prism + n_prism_part;
      PDM_l_num_t idx_poly3d = idx_pyramid + n_pyramid_part;

      n_poly3d_part = 0;
      for (int i = 0; i < n_cell_courant; i++) {
        num_cell_parent_to_local_courant[i] = 0;
        PDM_Mesh_nodal_elt_t cell_type = _type_cell_3D(cell_face_nb_courant[i],
                                                       cell_face_courant + cell_face_idx_courant[i] - adjust,
                                                       face_som_idx_courant,
                                                       face_som_nb_courant,
                                                       face_som_courant,
                                                       cell_som_tria,
                                                       cell_som_quad);

        switch(cell_type) {
        case PDM_MESH_NODAL_TETRA4 :
          _connec_tetra(mesh->vtx[i_part],
                        cell_som_tria,
                        connec_tetra_courant);
          *numabs_tetra_courant = numabs_courant[i];
          numabs_tetra_courant += 1;
          connec_tetra_courant += 4;
          *num_parent_tetra_courant = i;
          num_parent_tetra_courant += 1;
          num_cell_parent_to_local_courant[i] = idx_tetra++;
          break;
        case PDM_MESH_NODAL_HEXA8 :
          _connec_hexa(mesh->vtx[i_part],
                       cell_som_quad,
                       connec_hexa_courant);
          *numabs_hexa_courant = numabs_courant[i];
          numabs_hexa_courant += 1;
          connec_hexa_courant += 8;
          *num_parent_hexa_courant = i;
          num_parent_hexa_courant += 1;
          num_cell_parent_to_local_courant[i] = idx_hexa++;
          break;
        case PDM_MESH_NODAL_PRISM6 :
          _connec_prism(mesh->vtx[i_part],
                        cell_som_tria,
                        cell_som_quad,
                        connec_prism_courant);
          *numabs_prism_courant = numabs_courant[i];
          numabs_prism_courant += 1;
          connec_prism_courant += 6;
          *num_parent_prism_courant = i;
          num_parent_prism_courant += 1;
          num_cell_parent_to_local_courant[i] = idx_prism++;
          break;
        case PDM_MESH_NODAL_PYRAMID5 :
          _connec_pyramid(mesh->vtx[i_part],
                          cell_som_tria,
                          cell_som_quad,
                          connec_pyramid_courant);
          *numabs_pyramid_courant = numabs_courant[i];
          numabs_pyramid_courant += 1;
          connec_pyramid_courant += 5;
          *num_parent_pyramid_courant = i;
          num_parent_pyramid_courant += 1;
          num_cell_parent_to_local_courant[i] = idx_pyramid++;
          break;
        case PDM_MESH_NODAL_POLY_3D :
          {
            PDM_l_num_t *cell_face_cell = cell_face_courant + cell_face_idx_courant[i] - adjust;
            for (int j = 0; j < cell_face_nb_courant[i]; j++) {
              tag_face_poly3d[PDM_ABS(cell_face_cell[j]) - 1] = 0;
            }
            *numabs_poly3d_courant = numabs_courant[i];
            numabs_poly3d_courant += 1;
            l_cellfac_poly += cell_face_nb_courant[i];
            cellfac_poly_idx[n_poly3d_part+1] = l_cellfac_poly;
            n_poly3d_part += 1;
            *num_parent_poly3d_courant = i;
            num_parent_poly3d_courant += 1;
            num_cell_parent_to_local_courant[i] = idx_poly3d++;
            break;
          }
        default :
          break;
        }
      }

      if (n_poly3d_part > 0) {
        cellfac_poly = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t) * l_cellfac_poly);

        /* Stockage des faces du bloc */

        n_face_poly = 0;
        PDM_l_num_t l_facsom_poly = 0;
        for (int i = 0; i < n_face_part; i++) {
          if (tag_face_poly3d[i] == 0) {
            tag_face_poly3d[i] = n_face_poly++;
            l_facsom_poly += face_som_nb_courant[i];
          }
        }

        facsom_poly_idx = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t) * (n_face_poly + 1));
        facsom_poly = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t) * l_facsom_poly);

        facsom_poly_idx[0] = 0;
        PDM_l_num_t idx_facsom_poly = 0;
        PDM_l_num_t idx_facsom = 0;
        for (int i = 0; i < n_face_part; i++) {
          if (tag_face_poly3d[i] >= 0) {
            PDM_l_num_t ideb = face_som_idx_courant[i] - adjust;
            PDM_l_num_t ifin = ideb + face_som_nb_courant[i];
            facsom_poly_idx[idx_facsom+1] = facsom_poly_idx[idx_facsom] + face_som_nb_courant[i];
            idx_facsom += 1;
            for (int j = ideb; j < ifin; j++) {
              facsom_poly[idx_facsom_poly++] = face_som_courant[j];
            }
          }
        }

        /* Remplissage de la structure cellfac_poly */

        l_cellfac_poly = 0;
        for (int i = 0; i < n_cell_courant; i++) {
          PDM_Mesh_nodal_elt_t cell_type = _type_cell_3D(cell_face_nb_courant[i],
                                                         cell_face_courant + cell_face_idx_courant[i] - adjust,
                                                         face_som_idx_courant,
                                                         face_som_nb_courant,
                                                         face_som_courant,
                                                         cell_som_tria,
                                                         cell_som_quad);

          switch(cell_type) {

          case PDM_MESH_NODAL_POLY_3D :
            {
              PDM_l_num_t *cell_face_cell = cell_face_courant + cell_face_idx_courant[i] - adjust;
              for (int j = 0; j < cell_face_nb_courant[i]; j++) {
                cellfac_poly[l_cellfac_poly++] = tag_face_poly3d[PDM_ABS(cell_face_cell[j]) - 1] + 1;

                if (cell_face_cell[j] < 0) {
                  cellfac_poly[l_cellfac_poly-1] = -cellfac_poly[l_cellfac_poly-1];
                }

              }
              break;
            }
          default:
            break;
          }
        }
        free(tag_face_poly3d);
      }

      if (som_elts[0] > 0)
        PDM_Mesh_nodal_block_std_set(mesh,
                                     id_bloc_tetra4,
                                     i_part,
                                     n_tetra_part,
                                     connec_tetra,
                                     numabs_tetra,
                                     num_parent_tetra);

      if (som_elts[1] > 0)
        PDM_Mesh_nodal_block_std_set(mesh,
                                     id_bloc_hexa8,
                                     i_part,
                                     n_hexa_part,
                                     connec_hexa,
                                     numabs_hexa,
                                     num_parent_hexa);

      if (som_elts[2] > 0)
        PDM_Mesh_nodal_block_std_set(mesh,
                                     id_bloc_prism6,
                                     i_part,
                                     n_prism_part,
                                     connec_prism,
                                     numabs_prism,
                                     num_parent_prism);

      if (som_elts[3] > 0)
        PDM_Mesh_nodal_block_std_set(mesh,
                                     id_bloc_pyramid5,
                                     i_part,
                                     n_pyramid_part,
                                     connec_pyramid,
                                     numabs_pyramid,
                                     num_parent_pyramid);

      if (som_elts[4] > 0)
        PDM_Mesh_nodal_block_poly3d_set(mesh,
                                        id_bloc_poly_3d,
                                        i_part,
                                        n_poly3d_part,
                                        n_face_poly,
                                        facsom_poly_idx,
                                        facsom_poly,
                                        cellfac_poly_idx,
                                        cellfac_poly,
                                        numabs_poly3d,
                                        num_parent_poly3d);
    }

    if (mesh->prepa_blocks != NULL) {
      free(mesh->prepa_blocks->n_cell);
      free(mesh->prepa_blocks->n_face);
      free(mesh->prepa_blocks->n_tetra);
      free(mesh->prepa_blocks->n_hexa);
      free(mesh->prepa_blocks->n_prism);
      free(mesh->prepa_blocks->n_pyramid);
      free(mesh->prepa_blocks->n_poly3d);
      free(mesh->prepa_blocks->face_vtx_idx);
      free(mesh->prepa_blocks->face_vtx_nb);
      free(mesh->prepa_blocks->face_vtx);
      free(mesh->prepa_blocks->cell_face_idx);
      free(mesh->prepa_blocks->cell_face_nb);
      free(mesh->prepa_blocks->cell_face);
      free(mesh->prepa_blocks->add_etat);
      free(mesh->prepa_blocks->numabs);
      free(mesh->prepa_blocks);
      mesh->prepa_blocks = NULL;
    }
  }
}



/**
 * \brief  Add some 2D cells from cell edge conectivity.
 *
 * For each cell, this function searchs the type of the cell (tetrahedra, hexahedra, ...)
 * and stores it in the corresponding block. \ref ind_num gives the indirection
 * between old and new numbering.
 *
 * \param [in]  idx            Nodal mesh handle
 * \param [in]  id_part        Partition identifier
 * \param [in]  n_elt          Number of polyhedra
 * \param [in]  n_edge         Number of edges used to describe polyhedra
 * \param [in]  edge_vtx_idx   Index of edge vertex connectivity
 * \param [in]  edge_vtx_nb    Number of vertices for each edge
 * \param [in]  edge_vtx       Edge vertex connectivity
 * \param [in]  cell_edge_idx  Index of cell edge connectivity
 * \param [in]  cell_edge_nb   Number of edges for each cell
 * \param [in]  cell_edge      Cell edge connectivity
 * \param [in]  numabs         Global numbering
 *
 */

void
PDM_Mesh_nodal_cell2d_celledge_add
(
      PDM_Mesh_nodal_t *mesh,
const int               id_part,
const int               n_cell,
const int               n_edge,
const PDM_l_num_t      *edge_vtx_idx,
const PDM_l_num_t      *edge_vtx_nb,
const PDM_l_num_t      *edge_vtx,
const PDM_l_num_t      *cell_edge_idx,
const PDM_l_num_t      *cell_edge_nb,
const PDM_l_num_t      *cell_edge,
const PDM_g_num_t      *numabs
)
{
  int adjust = 0;
  if (n_cell > 0) {
    if (edge_vtx_idx[0] == 1) {
      adjust = 1;
    }
  }

  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }

  if (id_part >= mesh->n_part) {
    PDM_error(__FILE__, __LINE__, 0, "Partition identifier too big\n");
  }

  int n_part = 0;

  if (mesh->num_cell_parent_to_local == NULL) {
    mesh->num_cell_parent_to_local = (PDM_l_num_t **) malloc(sizeof(PDM_l_num_t *) * mesh->n_part);
    for (int i_part = 0; i_part < mesh->n_part; i_part++) {
      mesh->num_cell_parent_to_local[i_part] = NULL;
    }
  }

  mesh->num_cell_parent_to_local[id_part] = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t) * n_cell);
  for (int i = 0; i < n_cell; i++) {
    mesh->num_cell_parent_to_local[id_part][i] = 0;
  }

  if (mesh->prepa_blocks == NULL) {
    mesh->prepa_blocks = (PDM_Mesh_nodal_prepa_blocks_t *) malloc(sizeof(PDM_Mesh_nodal_prepa_blocks_t));
    mesh->prepa_blocks->t_add = 2;
    mesh->prepa_blocks->n_tria_proc = 0;    /* Nb de triangles par proc */
    mesh->prepa_blocks->n_quad_proc = 0;    /* Nb de quads par proc */
    mesh->prepa_blocks->n_poly2d_proc = 0;  /* Nb de poly2d par proc */
    mesh->prepa_blocks->n_cell = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t)*mesh->n_part);
    mesh->prepa_blocks->n_face = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t)*mesh->n_part);
    mesh->prepa_blocks->n_tria = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t)*mesh->n_part);
    mesh->prepa_blocks->n_quad = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t)*mesh->n_part);
    mesh->prepa_blocks->n_poly2d = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t)*mesh->n_part);
    mesh->prepa_blocks->l_connec_poly2d = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t)*mesh->n_part);
    mesh->prepa_blocks->face_vtx_idx = (PDM_l_num_t **) malloc(sizeof(PDM_l_num_t *)*mesh->n_part);
    mesh->prepa_blocks->face_vtx_nb = (PDM_l_num_t **) malloc(sizeof(PDM_l_num_t *)*mesh->n_part);
    mesh->prepa_blocks->face_vtx = (PDM_l_num_t **) malloc(sizeof(PDM_l_num_t *)*mesh->n_part);
    mesh->prepa_blocks->cell_face_idx = (PDM_l_num_t **) malloc(sizeof(PDM_l_num_t *)*mesh->n_part);
    mesh->prepa_blocks->cell_face_nb = (PDM_l_num_t **) malloc(sizeof(PDM_l_num_t *)*mesh->n_part);
    mesh->prepa_blocks->cell_face = (PDM_l_num_t **) malloc(sizeof(PDM_l_num_t *)*mesh->n_part);
    mesh->prepa_blocks->add_etat  = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t)*mesh->n_part);
    mesh->prepa_blocks->numabs = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t *)*mesh->n_part);
    for (int i = 0; i < mesh->n_part; i++) {
      mesh->prepa_blocks->add_etat[i] = 0;
    }
  }

  if (mesh->prepa_blocks->t_add != 2) {
    PDM_error(__FILE__, __LINE__, 0, "Erreur Cs_geom_cell2d_cell_face_add : Un autre type d'ajout est en cours\n");
    abort();
  }

  PDM_l_num_t n_tria    = 0;
  PDM_l_num_t n_quad    = 0;
  PDM_l_num_t n_poly2d  = 0;
  PDM_l_num_t l_connec_poly2d = 0;

  for (int i = 0; i < n_cell; i++) {

    PDM_l_num_t n_face_cell = cell_edge_nb[i];
    if (n_face_cell == 3)
      n_tria += 1;
    else if (n_face_cell == 4)
      n_quad += 1;
    else {
      n_poly2d  += 1;
      l_connec_poly2d += cell_edge_nb[i];
    }
  }

  mesh->prepa_blocks->n_tria_proc           += n_tria;
  mesh->prepa_blocks->n_quad_proc           += n_quad;
  mesh->prepa_blocks->n_poly2d_proc         += n_poly2d;
  mesh->prepa_blocks->add_etat[id_part]      = 1;
  mesh->prepa_blocks->n_cell[id_part]        = n_cell;
  mesh->prepa_blocks->n_tria[id_part]        = n_tria;
  mesh->prepa_blocks->n_quad[id_part]        = n_quad;
  mesh->prepa_blocks->n_poly2d[id_part]      = n_poly2d;
  mesh->prepa_blocks->l_connec_poly2d[id_part] = l_connec_poly2d;
  mesh->prepa_blocks->face_vtx_idx[id_part]  = (PDM_l_num_t *) edge_vtx_idx;
  mesh->prepa_blocks->face_vtx_nb[id_part]   = (PDM_l_num_t *) edge_vtx_nb;
  mesh->prepa_blocks->face_vtx[id_part]      = (PDM_l_num_t *) edge_vtx;
  mesh->prepa_blocks->cell_face_idx[id_part] = (PDM_l_num_t *) cell_edge_idx;
  mesh->prepa_blocks->cell_face_nb[id_part]  = (PDM_l_num_t *) cell_edge_nb;
  mesh->prepa_blocks->cell_face[id_part]     = (PDM_l_num_t *) cell_edge;
  mesh->prepa_blocks->numabs[id_part]        = (PDM_g_num_t *) numabs;
  mesh->prepa_blocks->add_etat[id_part]      = 1;
  mesh->prepa_blocks->n_face[id_part]        = n_edge;

  /* Creation des blocs si toutes les parts sont remplies */

  for (int i = 0; i < mesh->n_part; i++) {
    if (mesh->prepa_blocks->add_etat[i] == 1)
      n_part += 1;
  }

  if (mesh->n_part == n_part) {

    /* Creation des blocs */

    PDM_l_num_t elts[3];
    PDM_l_num_t som_elts[3];

    elts[0] = mesh->prepa_blocks->n_tria_proc > 0;
    elts[1] = mesh->prepa_blocks->n_quad_proc > 0;
    elts[2] = mesh->prepa_blocks->n_poly2d_proc > 0;

    PDM_MPI_Allreduce(elts, som_elts, 3, PDM_MPI_INT, PDM_MPI_SUM, mesh->pdm_mpi_comm);

    int id_bloc_tria3 = -1;
    int id_bloc_quad4 = -1;
    int id_bloc_poly_2d = -1;

    if (som_elts[0] > 0)
      id_bloc_tria3 = PDM_Mesh_nodal_block_add (mesh,
                                                PDM_TRUE,
                                                PDM_MESH_NODAL_TRIA3);

    if (som_elts[1] > 0)
      id_bloc_quad4 = PDM_Mesh_nodal_block_add (mesh,
                                                PDM_TRUE,
                                                PDM_MESH_NODAL_QUAD4);

    if (som_elts[2] > 0)
      id_bloc_poly_2d = PDM_Mesh_nodal_block_add (mesh,
                                                  PDM_TRUE,
                                                  PDM_MESH_NODAL_POLY_2D);

    /* Determination de la connectivite de chaque element */

    for (int i_part = 0; i_part < mesh->n_part; i_part++) {

      PDM_l_num_t n_cell_courant = mesh->prepa_blocks->n_cell[i_part];
      PDM_l_num_t *num_cell_parent_to_local_courant = mesh->num_cell_parent_to_local[i_part];
      PDM_l_num_t *face_som_courant = mesh->prepa_blocks->face_vtx[i_part];
      PDM_l_num_t *cell_face_idx_courant = mesh->prepa_blocks->cell_face_idx[i_part];
      PDM_l_num_t *cell_face_nb_courant = mesh->prepa_blocks->cell_face_nb[i_part];
      PDM_l_num_t *cell_face_courant = mesh->prepa_blocks->cell_face[i_part];
      PDM_g_num_t *numabs_courant = mesh->prepa_blocks->numabs[i_part];

      adjust = 0;
      if (n_cell_courant > 0) {
        if (cell_face_idx_courant[0] == 1) {
          adjust = 1;
        }
      }

      n_tria   = mesh->prepa_blocks->n_tria[i_part];
      n_quad    = mesh->prepa_blocks->n_quad[i_part];
      n_poly2d  = mesh->prepa_blocks->n_poly2d[i_part];
      l_connec_poly2d = mesh->prepa_blocks->l_connec_poly2d[i_part];

      PDM_l_num_t *connec_tria = NULL;
      PDM_l_num_t *connec_quad = NULL;
      PDM_l_num_t *connec_poly2d = NULL;
      PDM_l_num_t *connec_poly2d_idx = NULL;

      PDM_g_num_t *numabs_tria = NULL;
      PDM_g_num_t *numabs_quad = NULL;
      PDM_g_num_t *numabs_poly2d = NULL;

      PDM_l_num_t *num_parent_tria = NULL;
      PDM_l_num_t *num_parent_quad = NULL;
      PDM_l_num_t *num_parent_poly2d = NULL;

      if (som_elts[0] > 0) {
        connec_tria = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t) * 3 *n_tria);
        numabs_tria = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * n_tria);
        num_parent_tria = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t) * n_tria);
      }

      if (som_elts[1] > 0) {
        connec_quad = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t) * 4 * n_quad);
        numabs_quad = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * n_quad);
        num_parent_quad = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t) * n_quad);
      }

      if (som_elts[2] > 0) {
        connec_poly2d_idx = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t) * (n_poly2d + 1));
        connec_poly2d_idx[0] = 0;
        connec_poly2d = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t) * l_connec_poly2d);
        numabs_poly2d = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * n_poly2d);
        num_parent_poly2d = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t) * n_poly2d);
      }


      PDM_l_num_t *connec_tria_courant = connec_tria;
      PDM_l_num_t *connec_quad_courant = connec_quad;
      PDM_l_num_t *connec_poly2d_idx_courant = connec_poly2d_idx + 1;
      PDM_l_num_t *connec_poly2d_courant = connec_poly2d;

      PDM_g_num_t *numabs_tria_courant = numabs_tria;
      PDM_g_num_t *numabs_quad_courant = numabs_quad;
      PDM_g_num_t *numabs_poly2d_courant = numabs_poly2d;

      PDM_l_num_t *num_parent_tria_courant = num_parent_tria;
      PDM_l_num_t *num_parent_quad_courant = num_parent_quad;
      PDM_l_num_t *num_parent_poly2d_courant = num_parent_poly2d;

      /* Construction de la connectivit� sommet-> arrete */

      PDM_l_num_t *connec_som_are =
        (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t) * 2 * mesh->vtx[i_part]->n_vtx);

      PDM_l_num_t idx_tria   = 0;
      PDM_l_num_t idx_quad   = n_tria;
      PDM_l_num_t idx_poly2d = idx_quad + n_quad;

      for (int j = 0; j < 2 * mesh->vtx[i_part]->n_vtx; j++) {
        connec_som_are[j] = -1;
      }

      for (int i = 0; i < n_cell_courant; i++) {

        PDM_l_num_t ideb = cell_face_idx_courant[i] - adjust;
        PDM_l_num_t n_face_cell = cell_face_nb_courant[i];
        PDM_l_num_t ifin = ideb + n_face_cell;

        for (int j = ideb; j < ifin; j++) {
          PDM_l_num_t ifac = PDM_ABS(cell_face_courant[j]) - 1;
          PDM_l_num_t isom1 = face_som_courant[2*ifac] - 1;
          PDM_l_num_t isom2 = face_som_courant[2*ifac+1] - 1;

          if (connec_som_are[2*isom1] == -1)
            connec_som_are[2*isom1] = ifac;
          else
            connec_som_are[2*isom1+1] = ifac;

          if (connec_som_are[2*isom2] == -1)
            connec_som_are[2*isom2] = ifac;
          else
            connec_som_are[2*isom2+1] = ifac;
        }

        PDM_l_num_t *connec_courant;
        if (n_face_cell == 3) {
          *num_parent_tria_courant = i;
          num_parent_tria_courant += 1;
          num_cell_parent_to_local_courant[i] = idx_tria++;
          *numabs_tria_courant = numabs_courant[i];
          numabs_tria_courant += 1;
          connec_courant = connec_tria_courant;
          connec_tria_courant += n_face_cell;
        }
        else if (n_face_cell == 4) {
          *num_parent_quad_courant = i;
          num_parent_quad_courant += 1;
          num_cell_parent_to_local_courant[i] = idx_quad++;;
          *numabs_quad_courant = numabs_courant[i];
          numabs_quad_courant += 1;
          connec_courant = connec_quad_courant;
          connec_quad_courant += n_face_cell;
        }
        else {
          *num_parent_poly2d_courant = i;
          num_parent_poly2d_courant += 1;
          num_cell_parent_to_local_courant[i] = idx_poly2d++;
          *numabs_poly2d_courant = numabs_courant[i];
          numabs_poly2d_courant += 1;
          connec_courant = connec_poly2d_courant;
          *connec_poly2d_idx_courant = *(connec_poly2d_idx_courant - 1) +  n_face_cell;
          connec_poly2d_idx_courant += 1;
          connec_poly2d_courant += n_face_cell;
        }

        /* Remplissage de la connectivite */

        PDM_l_num_t idx_som = 0;
        PDM_l_num_t face_courant = PDM_ABS(cell_face_courant[ideb]) - 1;
        PDM_l_num_t isom1 = face_som_courant[2*face_courant] - 1;
        PDM_l_num_t isom_suiv = face_som_courant[2*face_courant + 1] - 1;
        connec_courant[idx_som++] = isom1 + 1;

        while (isom1 != isom_suiv) {
          assert(idx_som <= n_face_cell);
          connec_courant[idx_som++] = isom_suiv + 1;

          /* Face suivante */

          PDM_l_num_t face_suiv = connec_som_are[2*isom_suiv];
          if (face_suiv == face_courant)
            face_suiv = connec_som_are[2*isom_suiv + 1];
          face_courant = face_suiv;

          /* Sommet suivant */

          PDM_l_num_t isom_tmp = face_som_courant[2*face_courant] - 1;
          if (isom_tmp == isom_suiv)
            isom_tmp = face_som_courant[2*face_courant + 1] - 1;
          isom_suiv = isom_tmp;
        }

        for (int j= 0; j < n_face_cell; j++) {
          connec_som_are[2*(connec_courant[j] -1)] = - 1;
          connec_som_are[2*(connec_courant[j] -1) + 1] = - 1;
        }
      }

      free(connec_som_are);

      if (som_elts[0] > 0)
        PDM_Mesh_nodal_block_std_set(mesh,
                                     id_bloc_tria3,
                                     i_part,
                                     n_tria,
                                     connec_tria,
                                     numabs_tria,
                                     num_parent_tria);

      if (som_elts[1] > 0)
        PDM_Mesh_nodal_block_std_set(mesh,
                                     id_bloc_quad4,
                                     i_part,
                                     n_quad,
                                     connec_quad,
                                     numabs_quad,
                                     num_parent_quad);

      if (som_elts[2] > 0)
        PDM_Mesh_nodal_block_poly2d_set(mesh,
                                        id_bloc_poly_2d,
                                        i_part,
                                        n_poly2d,
                                        connec_poly2d_idx,
                                        connec_poly2d,
                                        numabs_poly2d,
                                        num_parent_poly2d);
    }
    if (mesh->prepa_blocks != NULL) {
      free(mesh->prepa_blocks->n_cell);
      free(mesh->prepa_blocks->n_face);
      free(mesh->prepa_blocks->n_tria);
      free(mesh->prepa_blocks->n_quad);
      free(mesh->prepa_blocks->n_poly2d);
      free(mesh->prepa_blocks->l_connec_poly2d);
      free(mesh->prepa_blocks->face_vtx_idx);
      free(mesh->prepa_blocks->face_vtx_nb);
      free(mesh->prepa_blocks->face_vtx);
      free(mesh->prepa_blocks->cell_face_idx);
      free(mesh->prepa_blocks->cell_face_nb);
      free(mesh->prepa_blocks->cell_face);
      free(mesh->prepa_blocks->add_etat);
      free(mesh->prepa_blocks->numabs);
      free(mesh->prepa_blocks);
      mesh->prepa_blocks = NULL;
    }
  }
}


/**
 * \brief  Add some 2D cells from cell vertex connectivity.
 *
 * For each cell, this function searchs the type of the cell (tetrahedra, hexahedra, ...)
 * and stores it in the corresponding block. \ref ind_num gives the indirection
 * between old and new numbering.
 *
 * \param [in]  idx            Nodal mesh handle
 * \param [in]  id_part        Partition identifier
 * \param [in]  n_face         Number of polygon
 * \param [in]  face_vtx_idx   Index of edge vertex connectivity
 * \param [in]  face_vtx_nb    Number of vertices for each edge
 * \param [in]  face_vtx       Edge vertex connectivity
 *
 */

void
PDM_Mesh_nodal_faces_facevtx_add
(
      PDM_Mesh_nodal_t *mesh,
const int               id_part,
const int               n_face,
const PDM_l_num_t      *face_vtx_idx,
const PDM_l_num_t      *face_vtx_nb,
const PDM_l_num_t      *face_vtx,
const PDM_g_num_t      *numabs
)
{

  int adjust = 0;
  if (n_face > 0) {
    if (face_vtx_idx[0] == 1) {
      adjust = 1;
    }
  }


  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }

  if (id_part >= mesh->n_part) {
    PDM_error(__FILE__, __LINE__, 0, "Partition identifier too big\n");
  }

  int n_part = 0;

  if (mesh->num_cell_parent_to_local == NULL) {
    mesh->num_cell_parent_to_local = (PDM_l_num_t **) malloc(sizeof(PDM_l_num_t *) * mesh->n_part);
    for (int i_part = 0; i_part < mesh->n_part; i_part++) {
      mesh->num_cell_parent_to_local[i_part] = NULL;
    }
  }

  mesh->num_cell_parent_to_local[id_part] = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t) * n_face);
  for (int i = 0; i < n_face; i++) {
    mesh->num_cell_parent_to_local[id_part][i] = 0;
  }

  if (mesh->prepa_blocks == NULL) {
    mesh->prepa_blocks = (PDM_Mesh_nodal_prepa_blocks_t *) malloc(sizeof(PDM_Mesh_nodal_prepa_blocks_t));
    mesh->prepa_blocks->t_add = 3;
    mesh->prepa_blocks->n_tria_proc = 0;    /* Nb de triangles par proc */
    mesh->prepa_blocks->n_quad_proc = 0;    /* Nb de quads par proc */
    mesh->prepa_blocks->n_poly2d_proc = 0;  /* Nb de poly2d par proc */
    mesh->prepa_blocks->n_face = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t)*mesh->n_part);
    mesh->prepa_blocks->n_tria = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t)*mesh->n_part);
    mesh->prepa_blocks->n_quad = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t)*mesh->n_part);
    mesh->prepa_blocks->n_poly2d = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t)*mesh->n_part);
    mesh->prepa_blocks->l_connec_poly2d = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t)*mesh->n_part);
    mesh->prepa_blocks->face_vtx_idx = (PDM_l_num_t **) malloc(sizeof(PDM_l_num_t *)*mesh->n_part);
    mesh->prepa_blocks->face_vtx_nb = (PDM_l_num_t **) malloc(sizeof(PDM_l_num_t *)*mesh->n_part);
    mesh->prepa_blocks->face_vtx = (PDM_l_num_t **) malloc(sizeof(PDM_l_num_t *)*mesh->n_part);
    mesh->prepa_blocks->add_etat  = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t)*mesh->n_part);
    mesh->prepa_blocks->numabs = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t *)*mesh->n_part);
    for (int i = 0; i < mesh->n_part; i++) {
      mesh->prepa_blocks->add_etat[i] = 0;
    }
  }

  if (mesh->prepa_blocks->t_add != 3) {
    PDM_error(__FILE__, __LINE__, 0, "Erreur Cs_geom_cell2d_cell_face_add : Un autre type d'ajout est en cours\n");
    abort();
  }

  PDM_l_num_t n_tria    = 0;
  PDM_l_num_t n_quad    = 0;
  PDM_l_num_t n_poly2d  = 0;
  PDM_l_num_t l_connec_poly2d  = 0;

  for (int i = 0; i < n_face; i++) {

    PDM_l_num_t n_som_face = face_vtx_nb[i];
    if (n_som_face == 3)
      n_tria += 1;
    else if (n_som_face == 4)
      n_quad += 1;
    else {
      n_poly2d  += 1;
      l_connec_poly2d += n_som_face;
    }
  }

  mesh->prepa_blocks->n_tria_proc           += n_tria;
  mesh->prepa_blocks->n_quad_proc           += n_quad;
  mesh->prepa_blocks->n_poly2d_proc         += n_poly2d;
  mesh->prepa_blocks->add_etat[id_part]      = 1;
  mesh->prepa_blocks->n_tria[id_part]        = n_tria;
  mesh->prepa_blocks->n_quad[id_part]        = n_quad;
  mesh->prepa_blocks->n_poly2d[id_part]      = n_poly2d;
  mesh->prepa_blocks->l_connec_poly2d[id_part] = l_connec_poly2d;
  mesh->prepa_blocks->face_vtx_idx[id_part]  = (PDM_l_num_t *) face_vtx_idx;
  mesh->prepa_blocks->face_vtx_nb[id_part]   = (PDM_l_num_t *) face_vtx_nb;
  mesh->prepa_blocks->face_vtx[id_part]      = (PDM_l_num_t *) face_vtx;
  mesh->prepa_blocks->numabs[id_part]        = (PDM_g_num_t *) numabs;
  mesh->prepa_blocks->add_etat[id_part]      = 1;
  mesh->prepa_blocks->n_face[id_part]        = n_face;

  /* Creation des blocs si toutes les parts sont remplies */

  for (int i = 0; i < mesh->n_part; i++) {
    if (mesh->prepa_blocks->add_etat[i] == 1)
      n_part += 1;
  }

  if (mesh->n_part == n_part) {

    /* Creation des blocs */

    PDM_l_num_t elts[3];
    PDM_l_num_t som_elts[3];

    elts[0] = mesh->prepa_blocks->n_tria_proc > 0;
    elts[1] = mesh->prepa_blocks->n_quad_proc > 0;
    elts[2] = mesh->prepa_blocks->n_poly2d_proc > 0;

    PDM_MPI_Allreduce(elts, som_elts, 3, PDM_MPI_INT, PDM_MPI_SUM, mesh->pdm_mpi_comm);

    int id_bloc_tria3 = -1;
    int id_bloc_quad4 = -1;
    int id_bloc_poly_2d = -1;

    if (som_elts[0] > 0) {
      id_bloc_tria3 = PDM_Mesh_nodal_block_add (mesh,
                                                PDM_TRUE,
                                                PDM_MESH_NODAL_TRIA3);
    }

    if (som_elts[1] > 0) {
      id_bloc_quad4 = PDM_Mesh_nodal_block_add (mesh,
                                                PDM_TRUE,
                                                PDM_MESH_NODAL_QUAD4);
    }

    if (som_elts[2] > 0) {
      id_bloc_poly_2d = PDM_Mesh_nodal_block_add (mesh,
                                                  PDM_TRUE,
                                                  PDM_MESH_NODAL_POLY_2D);
    }

    /* Determination de la connectivite de chaque element */

    for (int i_part = 0; i_part < mesh->n_part; i_part++) {

      PDM_l_num_t n_face_courant = mesh->prepa_blocks->n_face[i_part];
      PDM_l_num_t *num_cell_parent_to_local_courant = mesh->num_cell_parent_to_local[i_part];
      PDM_l_num_t *face_som_idx_courant = mesh->prepa_blocks->face_vtx_idx[i_part];
      PDM_l_num_t *face_som_nb_courant = mesh->prepa_blocks->face_vtx_nb[i_part];
      PDM_l_num_t *face_som_courant = mesh->prepa_blocks->face_vtx[i_part];
      PDM_g_num_t *numabs_courant = mesh->prepa_blocks->numabs[i_part];

      adjust = 0;
      if (n_face_courant > 0) {
        if (face_som_idx_courant[0] == 1) {
          adjust = 1;
        }
      }

      n_tria   = mesh->prepa_blocks->n_tria[i_part];
      n_quad    = mesh->prepa_blocks->n_quad[i_part];
      n_poly2d  = mesh->prepa_blocks->n_poly2d[i_part];
      l_connec_poly2d  = mesh->prepa_blocks->l_connec_poly2d[i_part];

      PDM_l_num_t *connec_tria = NULL;
      PDM_l_num_t *connec_quad = NULL;
      PDM_l_num_t *connec_poly2d = NULL;
      PDM_l_num_t *connec_poly2d_idx = NULL;

      PDM_g_num_t *numabs_tria = NULL;
      PDM_g_num_t *numabs_quad = NULL;
      PDM_g_num_t *numabs_poly2d = NULL;

      PDM_l_num_t *num_parent_tria = NULL;
      PDM_l_num_t *num_parent_quad = NULL;
      PDM_l_num_t *num_parent_poly2d = NULL;


      if (som_elts[0] > 0) {
        connec_tria = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t) * 3 *n_tria);
        numabs_tria = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * n_tria);
        num_parent_tria = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t) * n_tria);
      }

      if (som_elts[1] > 0) {
        connec_quad = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t) * 4 * n_quad);
        numabs_quad = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * n_quad);
        num_parent_quad = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t) * n_quad);
      }

      if (som_elts[2] > 0) {
        connec_poly2d_idx = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t) * (n_poly2d + 1));
        connec_poly2d_idx[0] = 0;
        connec_poly2d = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t) * l_connec_poly2d);
        numabs_poly2d = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * n_poly2d);
        num_parent_poly2d = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t) * n_poly2d);
      }

      PDM_l_num_t *connec_tria_courant = connec_tria;
      PDM_l_num_t *connec_quad_courant = connec_quad;
      PDM_l_num_t *connec_poly2d_idx_courant = connec_poly2d_idx + 1;
      PDM_l_num_t *connec_poly2d_courant = connec_poly2d;

      PDM_l_num_t *num_parent_tria_courant = num_parent_tria;
      PDM_l_num_t *num_parent_quad_courant = num_parent_quad;
      PDM_l_num_t *num_parent_poly2d_courant = num_parent_poly2d;

      PDM_g_num_t *numabs_tria_courant = numabs_tria;
      PDM_g_num_t *numabs_quad_courant = numabs_quad;
      PDM_g_num_t *numabs_poly2d_courant = numabs_poly2d;

      PDM_l_num_t idx_tria   = 0;
      PDM_l_num_t idx_quad   = n_tria;
      PDM_l_num_t idx_poly2d = idx_quad + n_quad;

      PDM_l_num_t n_som_face = 0;

      for (int i = 0; i < n_face_courant; i++) {
        n_som_face = face_som_nb_courant[i];
        PDM_l_num_t idx_som_face = face_som_idx_courant[i] - adjust;
        PDM_l_num_t *connec_courant;

        if (n_som_face == 3) {
          *num_parent_tria_courant = i;
          num_parent_tria_courant += 1;
          num_cell_parent_to_local_courant[i] = idx_tria++;
          *numabs_tria_courant = numabs_courant[i];
          numabs_tria_courant += 1;
          connec_courant = connec_tria_courant;
          connec_tria_courant += n_som_face;
        }
        else if (n_som_face == 4) {
          *num_parent_quad_courant = i;
          num_parent_quad_courant += 1;
          num_cell_parent_to_local_courant[i] = idx_quad++;;
          *numabs_quad_courant = numabs_courant[i];
          numabs_quad_courant += 1;
          connec_courant = connec_quad_courant;
          connec_quad_courant += n_som_face;
        }
        else {
          *num_parent_poly2d_courant = i;
          num_parent_poly2d_courant += 1;
          num_cell_parent_to_local_courant[i] = idx_poly2d++;
          *numabs_poly2d_courant = numabs_courant[i];
          numabs_poly2d_courant += 1;
          *connec_poly2d_idx_courant = *(connec_poly2d_idx_courant - 1) + n_som_face;
          connec_poly2d_idx_courant += 1;
          connec_courant = connec_poly2d_courant;
          connec_poly2d_courant += n_som_face;
        }

        /* Remplissage de la connectivite */

        for (int j = 0; j < n_som_face; j++)
          connec_courant[j] = face_som_courant[idx_som_face++];
      }

      if (som_elts[0] > 0)
        PDM_Mesh_nodal_block_std_set(mesh,
                                     id_bloc_tria3,
                                     i_part,
                                     n_tria,
                                     connec_tria,
                                     numabs_tria,
                                     num_parent_tria);

      if (som_elts[1] > 0)
        PDM_Mesh_nodal_block_std_set(mesh,
                                     id_bloc_quad4,
                                     i_part,
                                     n_quad,
                                     connec_quad,
                                     numabs_quad,
                                     num_parent_quad);

      if (som_elts[2] > 0)
        PDM_Mesh_nodal_block_poly2d_set(mesh,
                                        id_bloc_poly_2d,
                                        i_part,
                                        n_poly2d,
                                        connec_poly2d_idx,
                                        connec_poly2d,
                                        numabs_poly2d,
                                        num_parent_poly2d);
    }
    if (mesh->prepa_blocks != NULL) {
      free(mesh->prepa_blocks->n_face);
      free(mesh->prepa_blocks->n_tria);
      free(mesh->prepa_blocks->n_quad);
      free(mesh->prepa_blocks->n_poly2d);
      free(mesh->prepa_blocks->l_connec_poly2d);
      free(mesh->prepa_blocks->face_vtx_idx);
      free(mesh->prepa_blocks->face_vtx_nb);
      free(mesh->prepa_blocks->face_vtx);
      free(mesh->prepa_blocks->add_etat);
      free(mesh->prepa_blocks->numabs);
      free(mesh->prepa_blocks);
      mesh->prepa_blocks = NULL;
    }
  }
}


/**
 * \brief  Compute a global numbering in a block
 *
 * \param [in]  idx            Nodal mesh handle
 * \param [in]  id_block       Block identifier
 *
 */

void
PDM_Mesh_nodal_g_num_in_block_compute
(
      PDM_Mesh_nodal_t *mesh,
const int               id_block
)
{

  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }

  PDM_gen_gnum_t *gnum_gen = PDM_gnum_create (3, mesh->n_part,
                                              PDM_FALSE,
                                              1e-3,
                                              mesh->pdm_mpi_comm,
                                              PDM_OWNERSHIP_USER); /* The result is getted and you are owner */

  if (id_block >= PDM_BLOCK_ID_BLOCK_POLY3D) {

    int _id_block = id_block - PDM_BLOCK_ID_BLOCK_POLY3D;

    PDM_Mesh_nodal_block_poly3d_t *block = mesh->blocks_poly3d[_id_block];

    if (block == NULL) {
      PDM_error (__FILE__, __LINE__, 0, "Bad standard block identifier\n");
    }

    if (block->numabs_int == NULL) {
      block->numabs_int = (PDM_g_num_t **) malloc (sizeof(PDM_g_num_t *) * mesh->n_part);
      for (int i = 0; i < block->n_part; i++) {
        block->numabs_int[i] = NULL;
      }
    }
    else {
      return;
    }

    for (int i = 0; i < mesh->n_part; i++) {
      PDM_gnum_set_from_parents (gnum_gen, i, block->n_elt[i], block->_numabs[i]);
    }

  }

  else if (id_block >= PDM_BLOCK_ID_BLOCK_POLY2D) {

    int _id_block = id_block - PDM_BLOCK_ID_BLOCK_POLY2D;


    PDM_Mesh_nodal_block_poly2d_t *block = mesh->blocks_poly2d[_id_block];

    if (block == NULL) {
      PDM_error (__FILE__, __LINE__, 0, "Bad standard block identifier\n");
    }

    if (block->numabs_int == NULL) {
      block->numabs_int = (PDM_g_num_t **) malloc (sizeof(PDM_g_num_t *) * mesh->n_part);
      for (int i = 0; i < block->n_part; i++) {
        block->numabs_int[i] = NULL;
      }
    }
    else {
      return;
    }

    for (int i = 0; i < mesh->n_part; i++) {
      PDM_gnum_set_from_parents (gnum_gen, i, block->n_elt[i], block->_numabs[i]);
    }

  }

  else {

    int _id_block = id_block;

    PDM_Mesh_nodal_block_std_t *block = mesh->blocks_std[_id_block];

    if (block == NULL) {
      PDM_error (__FILE__, __LINE__, 0, "Bad standard block identifier\n");
    }

    if (block->numabs_int == NULL) {
      block->numabs_int = (PDM_g_num_t **) malloc (sizeof(PDM_g_num_t *) * mesh->n_part);
      for (int i = 0; i < block->n_part; i++) {
        block->numabs_int[i] = NULL;
      }
    }
    else {
      return;
    }

    for (int i = 0; i < mesh->n_part; i++) {
      PDM_gnum_set_from_parents (gnum_gen, i, block->n_elt[i], block->_numabs[i]);
    }

  }

  PDM_gnum_compute (gnum_gen);

  if (id_block >= PDM_BLOCK_ID_BLOCK_POLY3D) {
    int _id_block = id_block - PDM_BLOCK_ID_BLOCK_POLY3D;

    PDM_Mesh_nodal_block_poly3d_t *block = mesh->blocks_poly3d[_id_block];

    for (int i = 0; i < mesh->n_part; i++) {
      block->numabs_int[i] = (PDM_g_num_t *) PDM_gnum_get (gnum_gen, i);
    }
  }

  else if (id_block >= PDM_BLOCK_ID_BLOCK_POLY2D) {
    int _id_block = id_block - PDM_BLOCK_ID_BLOCK_POLY2D;

    PDM_Mesh_nodal_block_poly2d_t *block = mesh->blocks_poly2d[_id_block];

    for (int i = 0; i < mesh->n_part; i++) {
      block->numabs_int[i] = (PDM_g_num_t *) PDM_gnum_get (gnum_gen, i);
    }
  }

  else {

    int _id_block = id_block;

    PDM_Mesh_nodal_block_std_t *block = mesh->blocks_std[_id_block];

    for (int i = 0; i < mesh->n_part; i++) {
      block->numabs_int[i] = (PDM_g_num_t *) PDM_gnum_get (gnum_gen, i);
    }
  }

  PDM_gnum_free (gnum_gen);

}


/**
 * \brief Compute cell centers of a block
 *
 * \param [in]  idx            Nodal mesh handle
 * \param [in]  id_block       Block identifier
 * \param [in]  id_part        Partition identifier
 *
 */

void
PDM_Mesh_nodal_cell_centers_compute
(
      PDM_Mesh_nodal_t *mesh,
const int               id_block,
const int               i_part
)
{

  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }

  double**  characteristicLength = (double**)malloc(sizeof(double*)*mesh->n_part);
  int**     isDegenerated        = (int**)malloc(sizeof(int*)*mesh->n_part);

  if (id_block >= PDM_BLOCK_ID_BLOCK_POLY3D) {

    int _id_block = id_block - PDM_BLOCK_ID_BLOCK_POLY3D;


    PDM_Mesh_nodal_block_poly3d_t *block = mesh->blocks_poly3d[_id_block];

    if (block == NULL) {
      PDM_error (__FILE__, __LINE__, 0, "Bad standard block identifier\n");
    }

    double**  volume  = (double**)malloc(sizeof(double*)*mesh->n_part);

    double* coords = (double *) PDM_Mesh_nodal_vertices_get(mesh,i_part);
    volume[i_part] = (double*)malloc(sizeof(double)*block->n_elt[i_part]);

    if (block->cell_centers == NULL) {
      block->cell_centers = (double **) malloc (sizeof(double *) * mesh->n_part);
      for (int i = 0; i < block->n_part; i++) {
        block->cell_centers[i] = NULL;
      }
    }
    else {
      return;
    }

    if (block->cell_centers[i_part] == NULL) {
      block->cell_centers[i_part] = (double*) malloc(sizeof(double)*3*block->n_elt[i_part]);
    }

    characteristicLength[i_part] = (double*)malloc(sizeof(double)*block->n_elt[i_part]);
    isDegenerated[i_part] = (int*)malloc(sizeof(int)*block->n_elt[i_part]);

    PDM_geom_elem_polyhedra_properties(0,
                                       block->n_elt[i_part],
                                       block->n_face[i_part],
                                       block->_facvtx_idx[i_part],
                                       block->_facvtx[i_part],
                                       block->_cellfac_idx[i_part],
                                       block->_cellfac[i_part],
                                       mesh ->vtx[i_part]->n_vtx,
                                       coords,
                                       volume[i_part],
                                       block->cell_centers[i_part],
                                       characteristicLength[i_part],
                                       isDegenerated[i_part] );

    free(volume[i_part]);
    free(characteristicLength[i_part]);
    free(isDegenerated[i_part]);

    free(volume);
    free(characteristicLength);
    free(isDegenerated);

  }

  else if (id_block >= PDM_BLOCK_ID_BLOCK_POLY2D) {
    int _id_block = id_block - PDM_BLOCK_ID_BLOCK_POLY2D;


    PDM_Mesh_nodal_block_poly2d_t *block = mesh->blocks_poly2d[_id_block];

    if (block == NULL) {
      PDM_error (__FILE__, __LINE__, 0, "Bad standard block identifier\n");
    }

    double**  surface_vector  = (double**)malloc(sizeof(double*)*mesh->n_part);
    double* coords = (double *) PDM_Mesh_nodal_vertices_get(mesh,i_part);
    surface_vector[i_part] = (double*)malloc(sizeof(double)*3*block->n_elt[i_part]);

    if (block->cell_centers == NULL) {
      block->cell_centers = (double **) malloc (sizeof(double *) * mesh->n_part);
      for (int i = 0; i < block->n_part; i++) {
        block->cell_centers[i] = NULL;
      }
    }
    else {
      return;
    }

    if (block->cell_centers[i_part]==NULL) {
      block->cell_centers[i_part] = (double*)malloc(sizeof(double)*3*block->n_elt[i_part]);
    }

    characteristicLength[i_part] = (double*)malloc(sizeof(double)*block->n_elt[i_part]);
    isDegenerated[i_part] = (int*)malloc(sizeof(int)*block->n_elt[i_part]);

    PDM_geom_elem_polygon_properties(block->n_elt[i_part],
                                     block->_connec_idx[i_part],
                                     block->_connec[i_part],
                                     coords,
                                     surface_vector[i_part],
                                     block->cell_centers[i_part],
                                     characteristicLength[i_part],
                                     isDegenerated[i_part] );

    free(surface_vector[i_part]);
    free(characteristicLength[i_part]);
    free(isDegenerated[i_part]);

    free(surface_vector);
    free(characteristicLength);
    free(isDegenerated);

  }
  else {

    int _id_block = id_block;

    PDM_Mesh_nodal_block_std_t *block = mesh->blocks_std[_id_block];

    if (block == NULL) {
      PDM_error (__FILE__, __LINE__, 0, "Bad standard block identifier\n");
    }

    double**  surface_vector  = (double**)malloc(sizeof(double*)*mesh->n_part);
    double**  volume         = (double**)malloc(sizeof(double*)*mesh->n_part);
    double**  length         = (double**)malloc(sizeof(double*)*mesh->n_part);

    double* coords = (double *) PDM_Mesh_nodal_vertices_get(mesh,i_part);

    if (block->cell_centers == NULL) {
      block->cell_centers = (double **) malloc (sizeof(double *) * mesh->n_part);
      for (int i = 0; i < block->n_part; i++) {
        block->cell_centers[i] = NULL;
      }
    }
    else {
      return;
    }

    if (block->cell_centers[i_part] == NULL) {
      block->cell_centers[i_part] = (double*)malloc(sizeof(double)*3*block->n_elt[i_part]);
    }

    characteristicLength[i_part] = (double*)malloc(sizeof(double)*block->n_elt[i_part]);
    isDegenerated[i_part] = (int*)malloc(sizeof(int)*block->n_elt[i_part]);

    switch (block->t_elt) {
    case PDM_MESH_NODAL_POINT:
      memcpy(block->cell_centers[i_part],
             mesh->vtx[i_part]->coords,
             sizeof(double)*3*(block->n_elt[i_part]) );
      break;

    case PDM_MESH_NODAL_BAR2:
      length[i_part] = (double*)malloc(sizeof(double)*block->n_elt[i_part]);
      PDM_geom_elem_edges_properties(block->n_elt[i_part],
                                     block->_connec[i_part],
                                     coords,
                                     length[i_part],
                                     block->cell_centers[i_part],
                                     characteristicLength[i_part],
                                     isDegenerated[i_part]);
      break;

    case PDM_MESH_NODAL_TRIA3:
      surface_vector[i_part] = (double*)malloc(sizeof(double)*3*block->n_elt[i_part]);
      PDM_geom_elem_tria_properties(block->n_elt[i_part],
                                    block->_connec[i_part],
                                    coords,
                                    surface_vector[i_part],
                                    block->cell_centers[i_part],
                                    characteristicLength[i_part],
                                    isDegenerated[i_part]);


      break;

    case PDM_MESH_NODAL_QUAD4:
      surface_vector[i_part] = (double*)malloc(sizeof(double)*3*block->n_elt[i_part]);
      PDM_geom_elem_quad_properties(block->n_elt[i_part],
                                    block->_connec[i_part],
                                    coords,
                                    surface_vector[i_part],
                                    block->cell_centers[i_part],
                                    characteristicLength[i_part],
                                    isDegenerated[i_part]);
      break;

    case  PDM_MESH_NODAL_TETRA4:
      volume[i_part] = (double*)malloc(sizeof(double)*block->n_elt[i_part]);
      PDM_geom_elem_tetra_properties(block->n_elt[i_part],
                                     block->_connec[i_part],
                                     coords,
                                     volume[i_part],
                                     block->cell_centers[i_part],
                                     characteristicLength[i_part],
                                     isDegenerated[i_part]);
      break;

    case PDM_MESH_NODAL_PYRAMID5:
      volume[i_part] = (double*)malloc(sizeof(double)*block->n_elt[i_part]);
      PDM_geom_elem_pyramid_properties(block->n_elt[i_part],
                                       block->_connec[i_part],
                                       mesh->vtx[i_part]->n_vtx,
                                       coords,
                                       volume[i_part],
                                       block->cell_centers[i_part],
                                       characteristicLength[i_part],
                                       isDegenerated[i_part]);
      break;

    case PDM_MESH_NODAL_PRISM6:
      volume[i_part] = (double*)malloc(sizeof(double)*block->n_elt[i_part]);
      PDM_geom_elem_prism_properties(block->n_elt[i_part],
                                     block->_connec[i_part],
                                     mesh->vtx[i_part]->n_vtx,
                                     coords,
                                     volume[i_part],
                                     block->cell_centers[i_part],
                                     characteristicLength[i_part],
                                     isDegenerated[i_part]);
      break;
    case PDM_MESH_NODAL_HEXA8:
      volume[i_part] = (double*)malloc(sizeof(double)*block->n_elt[i_part]);
      PDM_geom_elem_hexa_properties(block->n_elt[i_part],
                                    block->_connec[i_part],
                                    mesh->vtx[i_part]->n_vtx,
                                    coords,
                                    volume[i_part],
                                    block->cell_centers[i_part],
                                    characteristicLength[i_part],
                                    isDegenerated[i_part]);
      break;
    case PDM_MESH_NODAL_POLY_2D:
    case PDM_MESH_NODAL_POLY_3D:
      break;
    default:
      break;
    }//end switch t_elt

    if ((block->t_elt ==  PDM_MESH_NODAL_TRIA3) ||
        (block->t_elt == PDM_MESH_NODAL_QUAD4)) {
      free(surface_vector[i_part]);
    }
    else if (block->t_elt == PDM_MESH_NODAL_BAR2) {
      free(length[i_part]);
    }
    else if ((block->t_elt == PDM_MESH_NODAL_TETRA4) ||
             (block->t_elt == PDM_MESH_NODAL_PYRAMID5) ||
             (block->t_elt == PDM_MESH_NODAL_PRISM6) ||
             (block->t_elt == PDM_MESH_NODAL_HEXA8)) {
      free(volume[i_part]);
    }

    free(characteristicLength[i_part]);
    free(isDegenerated[i_part]);

    free(characteristicLength);
    free(isDegenerated);
    free(surface_vector);
    free(length);
    free(volume);

  } // if id_block

}

/**
 * \brief Get global inside numbering of block elements
 *
 * \param [in]  idx            Nodal mesh handle
 * \param [in]  id_block       Block identifier
 * \param [in]  id_part        Partition identifier
 *
 * \return      Return global inside numbering of block elements
 *
 */

PDM_g_num_t *
PDM_Mesh_nodal_block_inside_g_num_get
(
      PDM_Mesh_nodal_t *mesh,
const int               id_block,
const int               id_part
)
{

  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }

  int _id_block;

  PDM_g_num_t *_numabs_int = NULL;

  if (id_block >= PDM_BLOCK_ID_BLOCK_POLY3D) {

    _id_block = id_block - PDM_BLOCK_ID_BLOCK_POLY3D;

    PDM_Mesh_nodal_block_poly3d_t *block = mesh->blocks_poly3d[_id_block];

    if (block == NULL) {
      PDM_error (__FILE__, __LINE__, 0, "Bad standard block identifier\n");
    }

    if (id_part >= block->n_part) {
      PDM_error(__FILE__, __LINE__, 0, "Partition identifier too big\n");
    }

    if (block->numabs_int != NULL) {
      _numabs_int = block->numabs_int[id_part];
    }
  }

  else if (id_block >= PDM_BLOCK_ID_BLOCK_POLY2D) {

    _id_block = id_block - PDM_BLOCK_ID_BLOCK_POLY2D;

    PDM_Mesh_nodal_block_poly2d_t *block = mesh->blocks_poly2d[_id_block];

    if (block == NULL) {
      PDM_error (__FILE__, __LINE__, 0, "Bad standard block identifier\n");
    }

    if (id_part >= block->n_part) {
      PDM_error(__FILE__, __LINE__, 0, "Partition identifier too big\n");
    }

    if (block->numabs_int != NULL) {
      _numabs_int = block->numabs_int[id_part];
    }
  }

  else {

    _id_block = id_block - PDM_BLOCK_ID_BLOCK_STD;

    PDM_Mesh_nodal_block_std_t *block = mesh->blocks_std[_id_block];

    if (block == NULL) {
      PDM_error (__FILE__, __LINE__, 0, "Bad standard block identifier\n");
    }

    if (id_part >= block->n_part) {
      PDM_error(__FILE__, __LINE__, 0, "Partition identifier too big\n");
    }

    if (block->numabs_int != NULL) {
      _numabs_int = block->numabs_int[id_part];
    }
  }

  return _numabs_int;
}

/**
 * \brief  Return parent cell number to local number
 *
 * \param [in]  idx       Nodal mesh handle
 * \param [in]  id_part   Partition identifier
 *
 * \return  Parent cell number to local number
 *
 */

int *
PDM_Mesh_nodal_num_cell_parent_to_local_get
(
      PDM_Mesh_nodal_t *mesh,
const int               id_part
)
{

  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }

  if (id_part >= mesh->n_part) {
    PDM_error(__FILE__, __LINE__, 0, "Partition identifier too big\n");
  }

  if (mesh->num_cell_parent_to_local != NULL)
    return mesh->num_cell_parent_to_local[id_part];
  else
    return NULL;

}


/**
 * \brief  Return number elements of a partition
 *
 * \param [in]  idx       Nodal mesh handle
 * \param [in]  id_part   Partition identifier
 *
 * \return  Return number elements of a partition
 *
 */

int
PDM_Mesh_nodal_n_cell_get
(
      PDM_Mesh_nodal_t *mesh,
const int               id_part
)
{

  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }

  if (id_part >= mesh->n_part) {
    PDM_error(__FILE__, __LINE__, 0, "Partition identifier too big\n");
  }

  return mesh->n_cell[id_part];

}


/**
 * \brief  Return parent  absolute number
 *
 * \param [in]  mesh           Nodal mesh
 *
 * \return  Parent of vertices
 *
 */

const PDM_g_num_t *
PDM_Mesh_nodal_vertices_g_num_parent_get
(
       PDM_Mesh_nodal_t *mesh,
 const int               id_part
 )
{
  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }

  if (id_part >= mesh->n_part) {
    PDM_error (__FILE__, __LINE__, 0, "Bad part identifier\n");
  }

  PDM_Mesh_nodal_vtx_t *vtx = mesh->vtx[id_part];

  assert(vtx->parent != NULL);

  return vtx->parent->_numabs;
}

/**
 * \brief Reset a nodal mesh structure
 *
 * \param [in]  idx   Nodal mesh handle
 *
 * \return      NULL
 *
 */

void
PDM_Mesh_nodal_reset
(
 PDM_Mesh_nodal_t *mesh
)
{
  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }

  if (mesh->blocks_std != NULL) {
    for (int i = 0; i < mesh->n_block_std; i++) {
      _block_std_free(mesh->blocks_std[i]);
    }
    free(mesh->blocks_std);
  }

  if (mesh->blocks_poly2d != NULL) {
    for (int i = 0; i < mesh->n_block_poly2d; i++) {
      _block_poly2d_free(mesh->blocks_poly2d[i]);
    }
    free(mesh->blocks_poly2d);
  }

  if (mesh->blocks_poly3d != NULL) {
    for (int i = 0; i < mesh->n_block_poly3d; i++) {
      _block_poly3d_free(mesh->blocks_poly3d[i]);
    }
    free(mesh->blocks_poly3d);
  }

  mesh->n_block_std              = 0;
  mesh->n_block_poly2d           = 0;
  mesh->n_block_poly3d           = 0;

  mesh->blocks_std               = NULL;
  mesh->blocks_poly2d            = NULL;
  mesh->blocks_poly3d            = NULL;
  if (mesh->blocks_id != NULL) {
    free (mesh->blocks_id);
  }
  mesh->blocks_id                = NULL;
  mesh->n_blocks                 = 0;
  mesh->prepa_blocks             = NULL;
  mesh->is_vtx_def_from_parent   = 0;

  if (mesh->num_cell_parent_to_local != NULL) {
    for (int i_part = 0; i_part < mesh->n_part; i_part++) {
      if (mesh->num_cell_parent_to_local[i_part] != NULL)
        free(mesh->num_cell_parent_to_local[i_part]);
    }
    free(mesh->num_cell_parent_to_local);
    mesh->num_cell_parent_to_local = NULL;
  }

  for (int i = 0; i < mesh->n_part; i++) {
    mesh->n_cell[i] = 0;
  }

  if (mesh->vtx != NULL) {
    for (int i = 0; i < mesh->n_part; i++) {
      mesh->vtx[i]->_coords = NULL;
      mesh->vtx[i]->_numabs = NULL;
      mesh->vtx[i]->_numparent = NULL;
      mesh->vtx[i]->n_vtx   = 0;
      if (mesh->vtx[i]->parent != NULL) {
        mesh->vtx[i]->parent =_vtx_free (mesh->vtx[i]->parent);
        mesh->vtx[i]->parent = NULL;
      }
      if (mesh->vtx[i]->coords != NULL) {
        free (mesh->vtx[i]->coords);
        mesh->vtx[i]->coords = NULL;
      }
    }
  }
}






/**
 * \brief Returns the number of vertices in an element
 *
 * \param [in]  type   Element type
 * \param [in]  order  Element order
 *
 * \return      Number of vertices in element
 *
 */

int
PDM_Mesh_nodal_n_vertices_element
(
 const PDM_Mesh_nodal_elt_t type,
 const int                  order
 )
{
  int n_vtx = 0;
  int _order = order;
  if (order == -1) {
    _order = 1;
  }

  switch(type) {
  case PDM_MESH_NODAL_POINT:          /* Point */
    n_vtx = 1;
    break;
  case PDM_MESH_NODAL_BAR2:           /* Edge */
    n_vtx = (_order+1);
    break;
  case PDM_MESH_NODAL_TRIA3:          /* Triangle */
    n_vtx = (_order+1)*(_order+2)/2;
    break;
  case PDM_MESH_NODAL_QUAD4:          /* Quadrangle */
    n_vtx = (_order+1)*(_order+1);
    break;
  case PDM_MESH_NODAL_POLY_2D:        /* Simple Polygon */
    n_vtx = -1;
    break;
  case PDM_MESH_NODAL_TETRA4:         /* Tetrahedron */
    n_vtx = (_order+1)*(_order+2)*(_order+3)/6;
    break;
  case PDM_MESH_NODAL_PYRAMID5:       /* Pyramid */
    n_vtx = (_order+1)*(_order+2)*(2*_order+3)/6;
    break;
  case PDM_MESH_NODAL_PRISM6:         /* Prism (pentahedron) */
    n_vtx = (_order+1)*(_order+1)*(_order+2)/2;
    break;
  case PDM_MESH_NODAL_HEXA8:          /* Hexahedron (brick) */
    n_vtx = (_order+1)*(_order+1)*(_order+1);
    break;
  case PDM_MESH_NODAL_POLY_3D:        /* Simple Polyhedron (convex or quasi-convex) */
    n_vtx = -1;
    break;
  default:
    n_vtx = -1;
  }

  return n_vtx;
}

void
PDM_Mesh_nodal_compute_cell_extents
(
       PDM_Mesh_nodal_t *mesh,
 const int               id_block,
 const int               id_part,
 const double            tolerance,
       double           *extents
 )
{
  const double eps_extents = 1.e-7;

  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }

  PDM_l_num_t *cell_vtx     = NULL;
  PDM_l_num_t *cell_vtx_idx = NULL;
  PDM_l_num_t  n_elt, n_vtx_elt = 0;

  int _id_block;

  /* Polyhedra */
  if (id_block >= PDM_BLOCK_ID_BLOCK_POLY3D) {

    _id_block = id_block - PDM_BLOCK_ID_BLOCK_POLY3D;

    PDM_Mesh_nodal_block_poly3d_t *block = mesh->blocks_poly3d[_id_block];

    if (block == NULL) {
      PDM_error (__FILE__, __LINE__, 0, "Bad standard block identifier\n");
    }

    n_elt = block->n_elt[id_part];
    cell_vtx_idx = block->_cellvtx_idx[id_part];
    cell_vtx     = block->_cellvtx[id_part];
  }

  /* Polygons */
  else if (id_block >= PDM_BLOCK_ID_BLOCK_POLY2D) {

    _id_block = id_block - PDM_BLOCK_ID_BLOCK_POLY2D;

    PDM_Mesh_nodal_block_poly2d_t *block = mesh->blocks_poly2d[_id_block];

    if (block == NULL) {
      PDM_error (__FILE__, __LINE__, 0, "Bad standard block identifier\n");
    }

    n_elt = block->n_elt[id_part];
    cell_vtx_idx = block->_connec_idx[id_part];
    cell_vtx     = block->_connec[id_part];
  }

  /* Standard elements */
  else {

    _id_block = id_block;

    PDM_Mesh_nodal_block_std_t *block = mesh->blocks_std[_id_block];

    if (block == NULL) {
      PDM_error (__FILE__, __LINE__, 0, "Bad standard block identifier\n");
    }

    n_elt = block->n_elt[id_part];
    cell_vtx = block->_connec[id_part];

    const int order = 1;//
    n_vtx_elt = PDM_Mesh_nodal_n_vertices_element (block->t_elt, order);
  }

  /* Get vertices coordinates of current part */
  double *coords = (double *) PDM_Mesh_nodal_vertices_get (mesh, id_part);


  /* Loop on elements */
  int idx = 0;
  for (PDM_l_num_t ielt = 0; ielt < n_elt; ielt++) {

    double *_extents = extents + 6 * ielt;

    for (int idim = 0; idim < 3; idim++) {
      _extents[idim]   =  HUGE_VAL;
      _extents[3+idim] = -HUGE_VAL;
    }


    if (id_block >= PDM_BLOCK_ID_BLOCK_POLY2D) {
      idx = cell_vtx_idx[ielt];
      n_vtx_elt = cell_vtx_idx[ielt+1] - cell_vtx_idx[ielt];
    }

    for (int ivtx = 0; ivtx < n_vtx_elt; ivtx++) {
      PDM_l_num_t id_vtx = cell_vtx[idx++] - 1;

      for (int idim = 0; idim < 3; idim++) {
        double x = coords[3*id_vtx + idim];

        if (x < _extents[idim]) {
          _extents[idim] = x;
        }
        if (x > _extents[3+idim]) {
          _extents[3+idim] = x;
        }
      }
    }


    /* Expand bounding box */
    double delta = 0.;
    for (int idim = 0; idim < 3; idim++) {
      double x = _extents[3+idim] - _extents[idim];

      if (delta < x) {
        delta = x;
      }
    }

    if (delta > eps_extents) {
      delta *= tolerance;
    } else {
      delta = eps_extents;
    }

    for (int idim = 0; idim < 3; idim++) {
      _extents[idim]   -= delta;
      _extents[3+idim] += delta;
    }
  } // End of loop on elements

}




PDM_l_num_t
PDM_Mesh_nodal_poly3d_cell_vtx_get
(
 const PDM_l_num_t   icell,
 const PDM_l_num_t   face_vtx_idx[],
 const PDM_l_num_t   face_vtx[],
 const PDM_l_num_t   cell_face_idx[],
 const PDM_l_num_t   cell_face[],
 PDM_l_num_t       **cell_vtx
 )
{
  PDM_l_num_t n_vtx_cell = 0;

  int n_face_cell = cell_face_idx[icell+1] - cell_face_idx[icell];
  int s_cell_vtx = 6 * n_face_cell;

  *cell_vtx = malloc (sizeof(PDM_l_num_t) * s_cell_vtx);
  PDM_l_num_t *_cell_vtx = *cell_vtx;

  /* Loop on faces */
  PDM_bool_t already_in_cell;
  for (PDM_l_num_t iface = cell_face_idx[icell]; iface < cell_face_idx[icell+1]; iface++) {
    int id_face = PDM_ABS (cell_face[iface]) - 1;

    for (int ivtx = face_vtx_idx[id_face]; ivtx < face_vtx_idx[id_face+1]; ivtx++) {
      PDM_l_num_t id_vtx = face_vtx[ivtx];// - 1;

      int pos =_binary_search (id_vtx,
                               _cell_vtx,
                               n_vtx_cell,
                               &already_in_cell);

      if (already_in_cell == PDM_TRUE) {
        continue;
      }

      if (pos >= s_cell_vtx) {
        s_cell_vtx *= 2;
        *cell_vtx = realloc (*cell_vtx, sizeof(PDM_l_num_t) * s_cell_vtx);
        _cell_vtx = *cell_vtx;
      }

      for (int j = n_vtx_cell; j > pos; j--) {
        _cell_vtx[j] = _cell_vtx[j-1];
      }
      _cell_vtx[pos] = id_vtx;
      n_vtx_cell++;

    }

  }

  if (s_cell_vtx > n_vtx_cell) {
    *cell_vtx = realloc (*cell_vtx, sizeof(PDM_l_num_t) * n_vtx_cell);
  }

  return n_vtx_cell;
}


/**
 * \brief Create a new Mesh nodal from elements selected in a parent Mesh nodal
 *
 * \param [in]   parent_mesh       Parent Mesh nodal structure
 * \param [in]   n_select_elt      Number of selected element for each partition of each nodal block
 * \param [in]   select_elt_l_num  Local numbers of selected elements (for each partition of each nodal block)
 *
 * \return       New mesh nodal
 *
 */

PDM_Mesh_nodal_t *
PDM_Mesh_nodal_extract_selection
(
 PDM_Mesh_nodal_t  *parent_mesh,
 const int        **n_select_elt,
 const int       ***select_elt_l_num
 )
{
  PDM_Mesh_nodal_t *child_mesh = PDM_Mesh_nodal_create (parent_mesh->n_part,
                                                        parent_mesh->pdm_mpi_comm);

  PDM_UNUSED (n_select_elt);
  PDM_UNUSED (select_elt_l_num);
#if 0
  int  *n_select_vtx = PDM_array_zeros_int (parent_mesh->n_part);
  int **select_vtx_l_num      = malloc (sizeof(int *) * parent_mesh->n_part);
  int **parent_to_child_l_num = malloc (sizeof(int *) * parent_mesh->n_part);

  /* Tag vertices incident to selected elements */
  PDM_l_num_t *cell_vtx     = NULL;
  PDM_l_num_t *cell_vtx_idx = NULL;
  for (int ipart = 0; ipart < parent_mesh->n_part; ipart++) {

    int part_n_vtx = PDM_Mesh_nodal_n_vertices_get (parent_mesh, ipart);
g
    parent_to_child_l_num[ipart] = PDM_array_zeros_int (part_n_vtx);

    for (int iblock = 0; iblock < parent_mesh->n_blocks; iblock++) {
      int id_block = parent_mesh->blocks_id[iblock];

      /* Polyhedra */
      if (id_block >= PDM_BLOCK_ID_BLOCK_POLY3D) {
        id_block -= PDM_BLOCK_ID_BLOCK_POLY3D;
        PDM_Mesh_nodal_block_poly3d_t *block = parent_mesh->blocks_poly3d[id_block];

        if (block == NULL) {
          PDM_error (__FILE__, __LINE__, 0, "Bad poly3d block identifier\n");
        }

        cell_vtx_idx = block->_cellvtx_idx[ipart];
        cell_vtx     = block->_cellvtx[ipart];

        for (int i = 0; i < n_select_elt[iblock][ipart]; i++) {
          int ielt = select_elt_l_num[iblock][ipart][i];
          for (int j = cell_vtx_idx[ielt]; j < cell_vtx_idx[ielt+1]; j++) {
            int ivtx = cell_vtx[j] - 1;
            if (parent_to_child_l_num[ipart][ivtx] <= 0) {
              parent_to_child_l_num[ipart][ivtx] = ++n_select_vtx[ipart];
            }
          }
        }
      }

      /* Polygons */
      else if (id_block >= PDM_BLOCK_ID_BLOCK_POLY2D) {
        id_block -= PDM_BLOCK_ID_BLOCK_POLY2D;
        PDM_Mesh_nodal_block_poly2d_t *block = parent_mesh->blocks_poly2d[id_block];

        if (block == NULL) {
          PDM_error (__FILE__, __LINE__, 0, "Bad poly2d block identifier\n");
        }

        cell_vtx_idx = block->_connec_idx[ipart];
        cell_vtx     = block->_connec[ipart];

        for (int i = 0; i < n_select_elt[iblock][ipart]; i++) {
          int ielt = select_elt_l_num[iblock][ipart][i];
          for (int j = cell_vtx_idx[ielt]; j < cell_vtx_idx[ielt+1]; j++) {
            int ivtx = cell_vtx[j] - 1;
            if (parent_to_child_l_num[ipart][ivtx] <= 0) {
              parent_to_child_l_num[ipart][ivtx] = ++n_select_vtx[ipart];
            }
          }
        }
      }

      /* Standard elements */
      else {
        PDM_Mesh_nodal_block_std_t *block = parent_mesh->blocks_std[id_block];

        if (block == NULL) {
          PDM_error (__FILE__, __LINE__, 0, "Bad standard block identifier\n");
        }

        cell_vtx = block->_connec[ipart];

        const int order = 1;//
        int n_vtx_elt = PDM_Mesh_nodal_n_vertices_element (block->t_elt, order);

        for (int i = 0; i < n_select_elt[iblock][ipart]; i++) {
          int ielt = select_elt_l_num[iblock][ipart][i];
          for (int j = 0; j < n_vtx_elt; j++) {
            int ivtx = cell_vtx[ielt*n_vtx_elt + j] - 1;
            if (parent_to_child_l_num[ipart][ivtx] <= 0) {
              parent_to_child_l_num[ipart][ivtx] = ++n_select_vtx[ipart];
            }
          }
        }
      }


      select_vtx_l_num[ipart] = malloc (sizeof(int) * n_select_vtx[ipart]);
      int idx = 0;
      for (int i = 0; i < part_n_vtx; i++) {
        if (parent_to_child_l_num[ipart][i] > 0) {
          select_vtx_l_num[ipart][idx++] = i;// +1?
        }
      }

      const double      *part_vtx_coord = PDM_Mesh_nodal_vertices_get (parent_mesh, ipart);
      const PDM_g_num_t *part_vtx_g_num = PDM_Mesh_nodal_vertices_g_num_get (parent_mesh, ipart);

      /*PDM_Mesh_nodal_coord_from_parent_set (child_mesh,
                                            ipart,
                                            n_select_vtx[ipart],
                                            part_n_vtx,
                                            numabs,
                                            select_vtx_l_num[ipart],
                                            part_vtx_coord,
                                            part_vtx_g_num);*/

    } // End of loop on nodal blocks
  } // End of loop on parts

#endif
  return child_mesh;
}














void
PDM_Mesh_nodal_write
(
 const char       *filename,
 PDM_Mesh_nodal_t *mesh
 )
{
  int i_rank;
  PDM_MPI_Comm_rank (mesh->pdm_mpi_comm, &i_rank);

  int id_cs = PDM_writer_create("Ensight",
                                PDM_WRITER_FMT_ASCII,
                                PDM_WRITER_TOPO_CONSTANTE,
                                PDM_WRITER_OFF,
                                "test_3d_ens",
                                filename,
                                PDM_MPI_COMM_WORLD,
                                PDM_IO_ACCES_MPI_SIMPLE,
                                1.,
                                NULL);

  /* Creation de la geometrie */
  int id_geom = PDM_writer_geom_create_from_mesh_nodal (id_cs,
                                                        "test3d_geom",
                                                        PDM_WRITER_OFF,
                                                        PDM_WRITER_OFF,
                                                        mesh);

  /* Creation des variables */
  int id_var_elt_part = PDM_writer_var_create(id_cs,
                                              PDM_WRITER_OFF,
                                              PDM_WRITER_VAR_SCALAIRE,
                                              PDM_WRITER_VAR_ELEMENTS,
                                              "num_part");

  int id_var_elt_gnum = PDM_writer_var_create(id_cs,
                                              PDM_WRITER_OFF,
                                              PDM_WRITER_VAR_SCALAIRE,
                                              PDM_WRITER_VAR_ELEMENTS,
                                              "elt_gnum");

  PDM_real_t **val_elt_part = (PDM_real_t **) malloc(sizeof(PDM_real_t *) * mesh->n_part);
  PDM_real_t **val_elt_gnum = (PDM_real_t **) malloc(sizeof(PDM_real_t *) * mesh->n_part);
  for (int ipart = 0; ipart < mesh->n_part; ipart++) {

    int n_cell = PDM_Mesh_nodal_n_cell_get (mesh, ipart);


    val_elt_part[ipart] = (PDM_real_t *) malloc(sizeof(PDM_real_t) * n_cell);
    val_elt_gnum[ipart] = (PDM_real_t *) malloc(sizeof(PDM_real_t) * n_cell);

    for (int i = 0; i < n_cell; i++) {
      val_elt_part[ipart][i] = 1 + i_rank*mesh->n_part + ipart;
    }
  }



  int ielt = 0;
  for (int iblock = 0; iblock < mesh->n_blocks; iblock++) {
    for (int ipart = 0; ipart < mesh->n_part; ipart++) {
      int n_elt = PDM_Mesh_nodal_block_n_elt_get (mesh,
                                                  mesh->blocks_id[iblock],
                                                  ipart);
      PDM_g_num_t *elt_gnum = PDM_Mesh_nodal_g_num_get (mesh,
                                                        mesh->blocks_id[iblock],
                                                        ipart);
      for (int i = 0; i < n_elt; i++) {
        val_elt_gnum[ipart][ielt++] = (PDM_real_t) elt_gnum[i];
      }
    }
  }



  /* Ecriture des variables */
  PDM_writer_step_beg(id_cs, 0.);

  for (int ipart = 0; ipart < mesh->n_part; ipart++) {
    PDM_writer_var_set(id_cs,
                       id_var_elt_part,
                       id_geom,
                       ipart,
                       val_elt_part[ipart]);

    PDM_writer_var_set(id_cs,
                       id_var_elt_gnum,
                       id_geom,
                       ipart,
                       val_elt_gnum[ipart]);
  }

  PDM_writer_var_write(id_cs,
                       id_var_elt_part);
  PDM_writer_var_write(id_cs,
                       id_var_elt_gnum);

  PDM_writer_var_free(id_cs,
                      id_var_elt_part);
  PDM_writer_var_free(id_cs,
                      id_var_elt_gnum);


  PDM_writer_step_end(id_cs);


  PDM_writer_geom_write(id_cs,
                        id_geom);

  for (int ipart = 0; ipart < mesh->n_part; ipart++) {
    free(val_elt_part[ipart]);
    free(val_elt_gnum[ipart]);
  }
  free(val_elt_part);
  free(val_elt_gnum);

  // Free...
}










#ifdef __cplusplus
}
#endif /* __cplusplus */
