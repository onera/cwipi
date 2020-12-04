
/*----------------------------------------------------------------------------
 *  System headers
 *----------------------------------------------------------------------------*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_error.h"
#include "pdm_dmesh_nodal_priv.h"
#include "pdm_dmesh_nodal_elements_utils.h"

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

/*============================================================================
 * Global variable
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief PDM_section_size_elt_faces_get
 *
 * \param [in]     mesh               Current mesh
 * \param [in]     id_section         Section identifier
 * \param [inout]  elt_face_vtx_idx   Index of element faces connectivity (preallocated)
 * \param [inout]  elt_face_vtx       Element faces connectivity (preallocated)
 *
 */
int
PDM_section_size_elt_faces_get
(
  PDM_dmesh_nodal_t *dmesh_nodal,
  int               *s_elt_face_vtx_idx,
  int               *s_elt_face_vtx,
  int               *s_elt_face_cell
)
{

  // You should do it for each level of the mesh !!!!
  // abort();
  printf("PDM_section_size_elt_faces_get WARNING NOT WORKING \n");

  int _s_elt_face_vtx_idx = 0;
  int _s_elt_face_vtx     = 0;

  for (int i = 0; i < dmesh_nodal->n_section_std; i++) {
    int n_face_elt     = PDM_n_face_elt_per_elmt    (dmesh_nodal->sections_std[i]->t_elt);
    int n_sum_vtx_face = PDM_n_sum_vtx_face_per_elmt(dmesh_nodal->sections_std[i]->t_elt);

    _s_elt_face_vtx_idx += dmesh_nodal->sections_std[i]->n_elt * n_face_elt;
    _s_elt_face_vtx     += dmesh_nodal->sections_std[i]->n_elt * n_sum_vtx_face;
  }

  for (int i = 0; i < dmesh_nodal->n_section_poly3d; i++) {
    int _n_face = dmesh_nodal->sections_poly3d[i]->n_face;
    _s_elt_face_vtx_idx += _n_face;
    _s_elt_face_vtx     += dmesh_nodal->sections_poly3d[i]->_face_vtx[dmesh_nodal->sections_poly3d[i]->_face_vtx_idx[_n_face]];
  }

  for (int i = 0; i < dmesh_nodal->n_section_poly2d; i++) {
    _s_elt_face_vtx_idx +=     dmesh_nodal->sections_poly2d[i]->_connec_idx[dmesh_nodal->sections_poly2d[i]->n_elt];
    _s_elt_face_vtx     += 2 * dmesh_nodal->sections_poly2d[i]->_connec_idx[dmesh_nodal->sections_poly2d[i]->n_elt];
  }

  *s_elt_face_cell    = _s_elt_face_vtx_idx;
  *s_elt_face_vtx_idx = _s_elt_face_vtx_idx + 1;
  *s_elt_face_vtx     = _s_elt_face_vtx     + 1;

  return *s_elt_face_vtx - 1;
}


/**
 *
 * \brief PDM_section_size_elt_edges_get
 *
 * \param [in]     mesh               Current mesh
 * \param [in]     id_section         Section identifier
 * \param [inout]  elt_edge_vtx_idx   Index of element faces connectivity (preallocated)
 * \param [inout]  elt_edge_vtx       Element faces connectivity (preallocated)
 *
 */
int
PDM_section_size_elt_edges_get
(
  PDM_dmesh_nodal_t *dmesh_nodal,
  int               *s_elt_edge_vtx_idx,
  int               *s_elt_edge_vtx,
  int               *s_elt_edge_cell
)
{
  printf("PDM_section_size_elt_edges_get WARNING NOT WORKING \n");

  int _s_elt_edge_vtx_idx = 0;
  int _s_elt_edge_vtx     = 0;

  for (int i = 0; i < dmesh_nodal->n_section_std; i++) {
    int n_edge_elt     = PDM_n_nedge_elt_per_elmt   (dmesh_nodal->sections_std[i]->t_elt);
    int n_sum_vtx_edge = PDM_n_sum_vtx_edge_per_elmt(dmesh_nodal->sections_std[i]->t_elt);

    _s_elt_edge_vtx_idx += dmesh_nodal->sections_std[i]->n_elt * n_edge_elt;
    _s_elt_edge_vtx     += dmesh_nodal->sections_std[i]->n_elt * n_sum_vtx_edge;
  }

  assert(dmesh_nodal->n_section_poly3d == 0); // Not implemented to test
  for (int i = 0; i < dmesh_nodal->n_section_poly3d; i++) {
    int _n_face = dmesh_nodal->sections_poly3d[i]->n_face;
    _s_elt_edge_vtx_idx +=     dmesh_nodal->sections_poly3d[i]->_face_vtx_idx[_n_face];
    _s_elt_edge_vtx     += 2 * dmesh_nodal->sections_poly3d[i]->_face_vtx_idx[_n_face];
  }

  assert(dmesh_nodal->n_section_poly2d == 0); // Not implemented
  for (int i = 0; i < dmesh_nodal->n_section_poly2d; i++) {
    _s_elt_edge_vtx_idx +=     dmesh_nodal->sections_poly2d[i]->_connec_idx[dmesh_nodal->sections_poly2d[i]->n_elt];
    _s_elt_edge_vtx     += 2 * dmesh_nodal->sections_poly2d[i]->_connec_idx[dmesh_nodal->sections_poly2d[i]->n_elt];
  }

  *s_elt_edge_cell    = _s_elt_edge_vtx_idx;
  *s_elt_edge_vtx_idx = _s_elt_edge_vtx_idx + 1;
  *s_elt_edge_vtx     = _s_elt_edge_vtx     + 1;

  return *s_elt_edge_vtx - 1;
}



/**
 * \brief Return for standard elements the number of face that build this element
 *
 */
int
PDM_n_face_elt_per_elmt
(
  PDM_Mesh_nodal_elt_t t_elt
)
{
  int n_face_elt = -1;
  switch (t_elt) {
   case PDM_MESH_NODAL_TRIA3:
     n_face_elt = 1;
     break;
   case PDM_MESH_NODAL_QUAD4:
     n_face_elt = 1;
     break;
   case PDM_MESH_NODAL_TETRA4:
     n_face_elt = 4;
     break;
   case PDM_MESH_NODAL_PYRAMID5:
     n_face_elt = 5;
     break;
   case PDM_MESH_NODAL_PRISM6:
     n_face_elt = 5;
     break;
   case PDM_MESH_NODAL_HEXA8:
     n_face_elt = 6;
     break;
   default:
     n_face_elt = -1;
     PDM_error(__FILE__, __LINE__, 0, "Error PDM_n_face_elt_per_elmt : Element type is not taking int account\n");
  }
  return n_face_elt;
}

/**
 * \brief Return for standard elements the number of edge that build this element
 *
 */
int
PDM_n_nedge_elt_per_elmt
(
  PDM_Mesh_nodal_elt_t t_elt
)
{
  int n_nedge_elt = -1;
  switch (t_elt) {
   case PDM_MESH_NODAL_POINT:
     n_nedge_elt = 0;
     break;
   case PDM_MESH_NODAL_BAR2:
     n_nedge_elt = 1;
     break;
   case PDM_MESH_NODAL_TRIA3:
     n_nedge_elt = 3;
     break;
   case PDM_MESH_NODAL_QUAD4:
     n_nedge_elt = 4;
     break;
   case PDM_MESH_NODAL_TETRA4:
     n_nedge_elt = 6;
     break;
   case PDM_MESH_NODAL_PYRAMID5:
     n_nedge_elt = 8;
     break;
   case PDM_MESH_NODAL_PRISM6:
     n_nedge_elt = 9;
     break;
   case PDM_MESH_NODAL_HEXA8:
     n_nedge_elt = 12;
     break;
   default:
     n_nedge_elt = -1;
     PDM_error(__FILE__, __LINE__, 0, "Error PDM_n_nedge_elt_per_elmt : Element type is not taking int account\n");
  }
  return n_nedge_elt;
}

/**
 * \brief Return for standard elements the total number of face vtx connectivity that build this element
 *
 */
int
PDM_n_sum_vtx_face_per_elmt
(
  PDM_Mesh_nodal_elt_t t_elt
)
{
  int n_sum_vtx_face = -1;
  switch (t_elt) {
   case PDM_MESH_NODAL_TRIA3:
     n_sum_vtx_face = 3;
     break;
   case PDM_MESH_NODAL_QUAD4:
     n_sum_vtx_face = 4;
     break;
   case PDM_MESH_NODAL_TETRA4:
     n_sum_vtx_face = 12;
     break;
   case PDM_MESH_NODAL_PYRAMID5:
     n_sum_vtx_face = 16;
     break;
   case PDM_MESH_NODAL_PRISM6:
     n_sum_vtx_face = 18;
     break;
   case PDM_MESH_NODAL_HEXA8:
     n_sum_vtx_face = 24;
     break;
   default:
     n_sum_vtx_face = -1;
     PDM_error(__FILE__, __LINE__, 0, "Error PDM_n_sum_vtx_face_per_elmt : Element type is not taking int account\n");
  }
  return n_sum_vtx_face;
}


/**
 * \brief Return for standard elements the total number of edge vtx connectivity that build this element
 *
 */
int
PDM_n_sum_vtx_edge_per_elmt
(
  PDM_Mesh_nodal_elt_t t_elt
)
{
  int n_sum_vtx_edge = -1;
  switch (t_elt) {
   case PDM_MESH_NODAL_POINT:
     n_sum_vtx_edge = 0;
     break;
   case PDM_MESH_NODAL_BAR2:
     n_sum_vtx_edge = 2;
     break;
   case PDM_MESH_NODAL_TRIA3:
     n_sum_vtx_edge = 6;
     break;
   case PDM_MESH_NODAL_QUAD4:
     n_sum_vtx_edge = 8;
     break;
   case PDM_MESH_NODAL_TETRA4:
     n_sum_vtx_edge = 12;
     break;
   case PDM_MESH_NODAL_PYRAMID5:
     n_sum_vtx_edge = 16;
     break;
   case PDM_MESH_NODAL_PRISM6:
     n_sum_vtx_edge = 18;
     break;
   case PDM_MESH_NODAL_HEXA8:
     n_sum_vtx_edge = 24;
     break;
   default:
     n_sum_vtx_edge = -1;
     PDM_error(__FILE__, __LINE__, 0, "Error PDM_n_sum_vtx_edge_per_elmt : Element type is not taking int account\n");
  }
  return n_sum_vtx_edge;
}

/**
*
* \brief Decompose quad cell_vtx connectivity to a flatten view of edges
*/
void
PDM_bar_decomposes_edges
(
       int          n_elt,
       int         *n_elt_current,
       int         *n_edge_current,
       PDM_g_num_t  beg_gnum_elt_current,
       PDM_g_num_t  beg_gnum_face_current,
 const PDM_g_num_t *connectivity_elmt_vtx,
       int         *elmt_edge_vtx_idx,
       PDM_g_num_t *elmt_edge_vtx,
       PDM_g_num_t *elmt_edge_cell,
       int         *elmt_cell_edge_idx,
       PDM_g_num_t *elmt_cell_edge
)
{
  PDM_UNUSED(elmt_cell_edge_idx);
  PDM_UNUSED(elmt_cell_edge);
  PDM_UNUSED(beg_gnum_face_current);

  const int n_edge_elt        = 1;
  const int n_sum_vtx_edge    = 2;
  const int n_sum_vtx_elt     = 2;

  int _n_edge_current = *n_edge_current;
  int         *_current_elmt_edge_vtx_idx = elmt_edge_vtx_idx + _n_edge_current;
  PDM_g_num_t *_current_elmt_edge_vtx     = elmt_edge_vtx + elmt_edge_vtx_idx[_n_edge_current];
  PDM_g_num_t *_current_elmt_edge_cell    = elmt_edge_cell + _n_edge_current;

  /*
   * For each element we flaten all connectivities in one array
   */
  for (int ielt = 0; ielt < n_elt; ielt++) {

    for (int i_edge = 0; i_edge < n_edge_elt; i_edge++) {
      _current_elmt_edge_vtx_idx[ielt * n_edge_elt + i_edge + 1] = _current_elmt_edge_vtx_idx[ielt * n_edge_elt + i_edge] + 2;
      _current_elmt_edge_cell   [ielt * n_edge_elt + i_edge    ] =  beg_gnum_elt_current + ielt + 1;
    }

    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 0]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt    ];
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 1]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 1];

  }

  *n_elt_current  += n_elt;
  *n_edge_current += n_elt * n_edge_elt;

}

/**
*
* \brief Decompose quad cell_vtx connectivity to a flatten view of edges
*/
void
PDM_tri_decomposes_edges
(
       int          n_elt,
       int         *n_elt_current,
       int         *n_edge_current,
       PDM_g_num_t  beg_gnum_elt_current,
       PDM_g_num_t  beg_gnum_face_current,
 const PDM_g_num_t *connectivity_elmt_vtx,
       int         *elmt_edge_vtx_idx,
       PDM_g_num_t *elmt_edge_vtx,
       PDM_g_num_t *elmt_edge_cell,
       int         *elmt_cell_edge_idx,
       PDM_g_num_t *elmt_cell_edge
)
{
  PDM_UNUSED(elmt_cell_edge_idx);
  PDM_UNUSED(elmt_cell_edge);
  PDM_UNUSED(beg_gnum_face_current);

  const int n_edge_elt        = 3;
  const int n_sum_vtx_edge    = 6;
  const int n_sum_vtx_elt     = 3;

  int _n_edge_current = *n_edge_current;
  int         *_current_elmt_edge_vtx_idx = elmt_edge_vtx_idx + _n_edge_current;
  PDM_g_num_t *_current_elmt_edge_vtx     = elmt_edge_vtx + elmt_edge_vtx_idx[_n_edge_current];
  PDM_g_num_t *_current_elmt_edge_cell    = elmt_edge_cell + _n_edge_current;

  /*
   * For each element we flaten all connectivities in one array
   */
  for (int ielt = 0; ielt < n_elt; ielt++) {

    for (int i_edge = 0; i_edge < n_edge_elt; i_edge++) {
      _current_elmt_edge_vtx_idx[ielt * n_edge_elt + i_edge + 1] = _current_elmt_edge_vtx_idx[ielt * n_edge_elt + i_edge] + 2;
      _current_elmt_edge_cell   [ielt * n_edge_elt + i_edge    ] =  beg_gnum_elt_current + ielt + 1;
    }

    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 0]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt    ];
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 1]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 1];

    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 2]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 1];
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 3]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 2];

    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 4]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 2];
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 5]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt    ];

  }

  *n_elt_current  += n_elt;
  *n_edge_current += n_elt * n_edge_elt;

}

/**
*
* \brief Decompose quad cell_vtx connectivity to a flatten view of edges
*/
void
PDM_quad_decomposes_edges
(
       int          n_elt,
       int         *n_elt_current,
       int         *n_edge_current,
       PDM_g_num_t  beg_gnum_elt_current,
       PDM_g_num_t  beg_gnum_face_current,
 const PDM_g_num_t *connectivity_elmt_vtx,
       int         *elmt_edge_vtx_idx,
       PDM_g_num_t *elmt_edge_vtx,
       PDM_g_num_t *elmt_edge_cell,
       int         *elmt_cell_edge_idx,
       PDM_g_num_t *elmt_cell_edge
)
{
  PDM_UNUSED(elmt_cell_edge_idx);
  PDM_UNUSED(elmt_cell_edge);
  PDM_UNUSED(beg_gnum_face_current);

  const int n_edge_elt        = 4;
  const int n_sum_vtx_edge    = 8;
  const int n_sum_vtx_elt     = 4;

  int _n_edge_current = *n_edge_current;
  int         *_current_elmt_edge_vtx_idx = elmt_edge_vtx_idx + _n_edge_current;
  PDM_g_num_t *_current_elmt_edge_vtx     = elmt_edge_vtx + elmt_edge_vtx_idx[_n_edge_current];
  PDM_g_num_t *_current_elmt_edge_cell    = elmt_edge_cell + _n_edge_current;

  /*
   * For each element we flaten all connectivities in one array
   */
  for (int ielt = 0; ielt < n_elt; ielt++) {

    for (int i_edge = 0; i_edge < n_edge_elt; i_edge++) {
      _current_elmt_edge_vtx_idx[ielt * n_edge_elt + i_edge + 1] = _current_elmt_edge_vtx_idx[ielt * n_edge_elt + i_edge] + 2;
      _current_elmt_edge_cell   [ielt * n_edge_elt + i_edge    ] =  beg_gnum_elt_current + ielt + 1;
    }

    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 0]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt    ];
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 1]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 1];

    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 2]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 1];
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 3]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 2];

    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 4]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 2];
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 5]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 3];

    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 6]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 3];
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 7]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt    ];
  }

  *n_elt_current  += n_elt;
  *n_edge_current += n_elt * n_edge_elt;

}

/**
*
* \brief Decompose quad cell_vtx connectivity to a flatten view of edges
*/
void
PDM_tri_decomposes_faces
(
       int          n_elt,
       int         *n_elt_current,
       int         *n_face_current,
       PDM_g_num_t  beg_gnum_elt_current,
       PDM_g_num_t  beg_gnum_face_current,
 const PDM_g_num_t *connectivity_elmt_vtx,
       int         *elmt_face_vtx_idx,
       PDM_g_num_t *elmt_face_vtx,
       PDM_g_num_t *elmt_face_cell,
       int         *elmt_cell_face_idx,
       PDM_g_num_t *elmt_cell_face
)
{
  PDM_UNUSED(elmt_cell_face_idx);
  PDM_UNUSED(elmt_cell_face);
  PDM_UNUSED(beg_gnum_face_current);

  const int n_face_elt        = 1;
  const int n_sum_vtx_face    = 3;
  const int n_sum_vtx_elt     = 3;

  int _n_face_current = *n_face_current;
  int         *_current_elmt_face_vtx_idx = elmt_face_vtx_idx + _n_face_current;
  PDM_g_num_t *_current_elmt_face_vtx     = elmt_face_vtx + elmt_face_vtx_idx[_n_face_current];
  PDM_g_num_t *_current_elmt_face_cell    = elmt_face_cell + _n_face_current;

  /*
   * For each element we flaten all connectivities in one array
   */
  for (int ielt = 0; ielt < n_elt; ielt++) {

    for (int i_face = 0; i_face < n_face_elt; i_face++) {
      _current_elmt_face_vtx_idx[ielt * n_face_elt + i_face + 1] = _current_elmt_face_vtx_idx[ielt * n_face_elt + i_face] + 3;
      _current_elmt_face_cell   [ielt * n_face_elt + i_face    ] = beg_gnum_elt_current + ielt + 1;
    }

    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 0]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt    ];
    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 1]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 1];
    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 2]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 2];

  }

  *n_elt_current  += n_elt;
  *n_face_current += n_elt * n_face_elt;
}

/**
*
* \brief Decompose quad cell_vtx connectivity to a flatten view of edges
*/
void
PDM_quad_decomposes_faces
(
       int          n_elt,
       int         *n_elt_current,
       int         *n_face_current,
       PDM_g_num_t  beg_gnum_elt_current,
       PDM_g_num_t  beg_gnum_face_current,
 const PDM_g_num_t *connectivity_elmt_vtx,
       int         *elmt_face_vtx_idx,
       PDM_g_num_t *elmt_face_vtx,
       PDM_g_num_t *elmt_face_cell,
       int         *elmt_cell_face_idx,
       PDM_g_num_t *elmt_cell_face
)
{
  PDM_UNUSED(elmt_cell_face_idx);
  PDM_UNUSED(elmt_cell_face);
  PDM_UNUSED(beg_gnum_face_current);

  const int n_face_elt        = 1;
  const int n_sum_vtx_face    = 4;
  const int n_sum_vtx_elt     = 4;

  int _n_face_current = *n_face_current;
  int         *_current_elmt_face_vtx_idx = elmt_face_vtx_idx + _n_face_current;
  PDM_g_num_t *_current_elmt_face_vtx     = elmt_face_vtx + elmt_face_vtx_idx[_n_face_current];
  PDM_g_num_t *_current_elmt_face_cell    = elmt_face_cell + _n_face_current;

  /*
   * For each element we flaten all connectivities in one array
   */
  for (int ielt = 0; ielt < n_elt; ielt++) {

    for (int i_face = 0; i_face < n_face_elt; i_face++) {
      _current_elmt_face_vtx_idx[ielt * n_face_elt + i_face + 1] = _current_elmt_face_vtx_idx[ielt * n_face_elt + i_face] + 4;
      _current_elmt_face_cell   [ielt * n_face_elt + i_face    ] = beg_gnum_elt_current + ielt + 1;
    }

    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 0]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt    ];
    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 1]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 1];
    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 2]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 2];
    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 3]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 3];

  }

  *n_elt_current  += n_elt;
  *n_face_current += n_elt * n_face_elt;
}

/**
*
* \brief Decompose quad cell_vtx connectivity to a flatten view of edges
*/
void
PDM_poly2d_decomposes_faces
(
       int          n_elt,
       int         *n_elt_current,
       int         *n_face_current,
       PDM_g_num_t  beg_gnum_elt_current,
       PDM_g_num_t  beg_gnum_face_current,
 const PDM_g_num_t *connectivity_elmt_vtx,
 const PDM_g_num_t *connectivity_elmt_vtx_idx,
       int         *elmt_face_vtx_idx,
       PDM_g_num_t *elmt_face_vtx,
       PDM_g_num_t *elmt_face_cell,
       int         *elmt_cell_face_idx,
       PDM_g_num_t *elmt_cell_face
)
{
  PDM_UNUSED(elmt_cell_face_idx);
  PDM_UNUSED(elmt_cell_face);
  PDM_UNUSED(beg_gnum_face_current);

  int _n_face_current = *n_face_current;
  int         *_current_elmt_face_vtx_idx = elmt_face_vtx_idx + _n_face_current;
  PDM_g_num_t *_current_elmt_face_vtx     = elmt_face_vtx + elmt_face_vtx_idx[_n_face_current];
  PDM_g_num_t *_current_elmt_face_cell    = elmt_face_cell + _n_face_current;

  /*
   * For each element we flaten all connectivities in one array
   */
  int idx = 0;
  for (int ielt = 0; ielt < n_elt; ielt++) {
    int beg = connectivity_elmt_vtx_idx[ielt];
    int n_vtx_on_face = connectivity_elmt_vtx_idx[ielt+1] - beg;
    _current_elmt_face_vtx_idx[idx + 1] = _current_elmt_face_vtx_idx[idx] + n_vtx_on_face;
    _current_elmt_face_cell   [idx    ] = beg_gnum_elt_current + ielt + 1;
    for(int ivtx = 0; ivtx < n_vtx_on_face; ++ivtx ) {
       _current_elmt_face_vtx[idx++] = connectivity_elmt_vtx[beg+ivtx];
    }
  }

  *n_elt_current  += n_elt;
  *n_face_current += n_elt;
}

/**
*
* \brief Decompose quad cell_vtx connectivity to a flatten view of edges
*/
void
PDM_poly2d_decomposes_edges
(
       int          n_elt,
       int         *n_elt_current,
       int         *n_edge_current,
       PDM_g_num_t  beg_gnum_elt_current,
       PDM_g_num_t  beg_gnum_edge_current,
 const PDM_g_num_t *connectivity_elmt_vtx,
 const PDM_g_num_t *connectivity_elmt_vtx_idx,
       int         *elmt_edge_vtx_idx,
       PDM_g_num_t *elmt_edge_vtx,
       PDM_g_num_t *elmt_edge_cell,
       int         *elmt_cell_edge_idx,
       PDM_g_num_t *elmt_cell_edge
)
{
  PDM_UNUSED(elmt_cell_edge_idx);
  PDM_UNUSED(elmt_cell_edge);
  PDM_UNUSED(beg_gnum_edge_current);

  int _n_edge_current = *n_edge_current;
  int         *_current_elmt_edge_vtx_idx = elmt_edge_vtx_idx + _n_edge_current;
  PDM_g_num_t *_current_elmt_edge_vtx     = elmt_edge_vtx + elmt_edge_vtx_idx[_n_edge_current];
  PDM_g_num_t *_current_elmt_edge_cell    = elmt_edge_cell + _n_edge_current;

  /*
   * For each element we flaten all connectivities in one array
   */
  int idx = 0;
  for (int ielt = 0; ielt < n_elt; ielt++) {
    // Reminder for poly2d -> Number of vertex = Number of edge
    int n_edge_elt = connectivity_elmt_vtx_idx[ielt+1] - connectivity_elmt_vtx_idx[ielt];
    *n_edge_current += n_edge_elt;

    int idx2 = connectivity_elmt_vtx_idx[ielt];
    for (int i_edge = 0; i_edge < n_edge_elt; i_edge++) {
      _current_elmt_edge_vtx_idx[idx + 1] = _current_elmt_edge_vtx_idx[idx] + 2;
      _current_elmt_edge_cell   [idx    ] = beg_gnum_elt_current + ielt + 1;

      int inext = (i_edge + 1) % n_edge_elt;
      _current_elmt_edge_vtx[2 * idx    ]  = connectivity_elmt_vtx[idx2 + i_edge];
      _current_elmt_edge_vtx[2 * idx + 1]  = connectivity_elmt_vtx[idx2 + inext ];

      idx += 1;
    }
  }

  *n_elt_current  += n_elt;
}



/**
*
* \brief Decompose quad cell_vtx connectivity to a flatten view of edges
*/
void
PDM_poly3d_decomposes_faces
(
       int          n_elt,
       int         *n_elt_current,
       int         *n_face_current,
       PDM_g_num_t  beg_gnum_elt_current,
       PDM_g_num_t  beg_gnum_face_current,
 const PDM_g_num_t *connectivity_elmt_vtx,
 const PDM_g_num_t *connectivity_elmt_vtx_idx,
       int         *elmt_face_vtx_idx,
       PDM_g_num_t *elmt_face_vtx,
       PDM_g_num_t *elmt_face_cell,
       int         *elmt_cell_face_idx,
       PDM_g_num_t *elmt_cell_face
)
{
  PDM_UNUSED(n_elt);
  PDM_UNUSED(n_elt_current);
  PDM_UNUSED(n_face_current);
  PDM_UNUSED(beg_gnum_elt_current);
  PDM_UNUSED(beg_gnum_face_current);
  PDM_UNUSED(connectivity_elmt_vtx);
  PDM_UNUSED(connectivity_elmt_vtx_idx);
  PDM_UNUSED(elmt_face_vtx_idx);
  PDM_UNUSED(elmt_face_vtx);
  PDM_UNUSED(elmt_face_cell);
  PDM_UNUSED(elmt_cell_face_idx);
  PDM_UNUSED(elmt_cell_face);

  abort();
}

/**
*
* \brief Decompose quad cell_vtx connectivity to a flatten view of edges
*/
void
PDM_poly3d_decomposes_edges
(
       int          n_elt,
       int         *n_elt_current,
       int         *n_edge_current,
       PDM_g_num_t  beg_gnum_elt_current,
       PDM_g_num_t  beg_gnum_edge_current,
 const PDM_g_num_t *connectivity_elmt_vtx,
 const PDM_g_num_t *connectivity_elmt_vtx_idx,
       int         *elmt_edge_vtx_idx,
       PDM_g_num_t *elmt_edge_vtx,
       PDM_g_num_t *elmt_edge_cell,
       int         *elmt_cell_edge_idx,
       PDM_g_num_t *elmt_cell_edge
)
{
  PDM_UNUSED(n_elt);
  PDM_UNUSED(n_elt_current);
  PDM_UNUSED(n_edge_current);
  PDM_UNUSED(beg_gnum_elt_current);
  PDM_UNUSED(beg_gnum_edge_current);
  PDM_UNUSED(connectivity_elmt_vtx);
  PDM_UNUSED(connectivity_elmt_vtx_idx);
  PDM_UNUSED(elmt_edge_vtx_idx);
  PDM_UNUSED(elmt_edge_vtx);
  PDM_UNUSED(elmt_edge_cell);
  PDM_UNUSED(elmt_cell_edge_idx);
  PDM_UNUSED(elmt_cell_edge);
  abort();
}

/**
*
* \brief Decompose tetra cell_vtx connectivity to a flatten view of faces
*/
void
PDM_tetra_decomposes_faces
(
       int          n_elt,
       int         *n_elt_current,
       int         *n_face_current,
       PDM_g_num_t  beg_gnum_elt_current,
       PDM_g_num_t  beg_gnum_face_current,
 const PDM_g_num_t *connectivity_elmt_vtx,
       int         *elmt_face_vtx_idx,
       PDM_g_num_t *elmt_face_vtx,
       PDM_g_num_t *elmt_face_cell,
       int         *elmt_cell_face_idx,
       PDM_g_num_t *elmt_cell_face
)
{
  PDM_UNUSED(elmt_cell_face_idx);
  PDM_UNUSED(elmt_cell_face);
  PDM_UNUSED(beg_gnum_elt_current);
  PDM_UNUSED(beg_gnum_face_current);

  const int n_face_elt        = 4;
  const int n_sum_vtx_face    = 12;
  const int n_sum_vtx_elt     = 4;

  int _n_face_current = *n_face_current;
  int         *_current_elmt_face_vtx_idx = elmt_face_vtx_idx + _n_face_current;
  PDM_g_num_t *_current_elmt_face_vtx     = elmt_face_vtx + elmt_face_vtx_idx[_n_face_current];
  PDM_g_num_t *_current_elmt_face_cell    = elmt_face_cell + _n_face_current;

  /*
   * For each element we flaten all connectivities in one array
   */
  for (int ielt = 0; ielt < n_elt; ielt++) {

    for (int i_face = 0; i_face < n_face_elt; i_face++) {
      _current_elmt_face_vtx_idx[ielt * n_face_elt + i_face + 1] = _current_elmt_face_vtx_idx[ielt * n_face_elt + i_face] + 3;
      _current_elmt_face_cell   [ielt * n_face_elt + i_face    ] = beg_gnum_elt_current + ielt + 1;
    }

    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 0]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt    ];
    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 1]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 1];
    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 2]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 2];

    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 3]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt    ];
    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 4]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 3];
    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 5]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 1];

    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 6]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt    ];
    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 7]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 2];
    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 8]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 3];

    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 9]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 1];
    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 10] = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 3];
    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 11] = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 2];

  }

  *n_elt_current  += n_elt;
  *n_face_current += n_elt * n_face_elt;

}

/**
*
* \brief Decompose tetra cell_vtx connectivity to a flatten view of faces
*/
void
PDM_tetra_decomposes_edges
(
       int          n_elt,
       int         *n_elt_current,
       int         *n_edge_current,
       PDM_g_num_t  beg_gnum_elt_current,
       PDM_g_num_t  beg_gnum_edge_current,
 const PDM_g_num_t *connectivity_elmt_vtx,
       int         *elmt_edge_vtx_idx,
       PDM_g_num_t *elmt_edge_vtx,
       PDM_g_num_t *elmt_edge_cell,
       int         *elmt_cell_edge_idx,
       PDM_g_num_t *elmt_cell_edge
)
{
  PDM_UNUSED(elmt_cell_edge_idx);
  PDM_UNUSED(elmt_cell_edge);
  PDM_UNUSED(beg_gnum_elt_current);
  PDM_UNUSED(beg_gnum_edge_current);

  const int n_edge_elt        = 6;
  const int n_sum_vtx_edge    = 12;
  const int n_sum_vtx_elt     = 4;

  int _n_edge_current = *n_edge_current;
  int         *_current_elmt_edge_vtx_idx = elmt_edge_vtx_idx + _n_edge_current;
  PDM_g_num_t *_current_elmt_edge_vtx     = elmt_edge_vtx + elmt_edge_vtx_idx[_n_edge_current];
  PDM_g_num_t *_current_elmt_edge_cell    = elmt_edge_cell + _n_edge_current;

  /*
   * For each element we flaten all connectivities in one array
   */
  for (int ielt = 0; ielt < n_elt; ielt++) {

    for (int i_edge = 0; i_edge < n_edge_elt; i_edge++) {
      _current_elmt_edge_vtx_idx[ielt * n_edge_elt + i_edge + 1] = _current_elmt_edge_vtx_idx[ielt * n_edge_elt + i_edge] + 2;
      _current_elmt_edge_cell   [ielt * n_edge_elt + i_edge    ] = beg_gnum_elt_current + ielt + 1;
    }

    // E1 = N1 N2
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 0]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt    ];
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 1]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 1];

    // E2 = N2 N3
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 2]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 1];
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 3]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 2];

    // E3 = N3 N1
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 4]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 2];
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 5]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt    ];

    // E4 = N1 N4
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 6]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt   ];
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 7]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 3];

    // E5 = N2 N4
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 8]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 1];
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 9]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 3];

    // E6 = N3 N4
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 10]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 2];
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 11]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 3];

  }

  *n_elt_current  += n_elt;
  *n_edge_current += n_elt * n_edge_elt;

}

/**
*
* \brief Decompose pyra cell_vtx connectivity to a flatten view of faces
*/
void
PDM_pyra_decomposes_faces
(
       int          n_elt,
       int         *n_elt_current,
       int         *n_face_current,
       PDM_g_num_t  beg_gnum_elt_current,
       PDM_g_num_t  beg_gnum_face_current,
 const PDM_g_num_t *connectivity_elmt_vtx,
       int         *elmt_face_vtx_idx,
       PDM_g_num_t *elmt_face_vtx,
       PDM_g_num_t *elmt_face_cell,
       int         *elmt_cell_face_idx,
       PDM_g_num_t *elmt_cell_face
)
{
  PDM_UNUSED(elmt_cell_face_idx);
  PDM_UNUSED(elmt_cell_face);
  PDM_UNUSED(beg_gnum_elt_current);
  PDM_UNUSED(beg_gnum_face_current);

  const int n_face_elt        = 5;
  const int n_sum_vtx_face    = 1*4 + 4*3;
  const int n_sum_vtx_elt     = 5;

  int _n_face_current = *n_face_current;
  int         *_current_elmt_face_vtx_idx = elmt_face_vtx_idx + _n_face_current;
  PDM_g_num_t *_current_elmt_face_vtx     = elmt_face_vtx + elmt_face_vtx_idx[_n_face_current];
  PDM_g_num_t *_current_elmt_face_cell    = elmt_face_cell + _n_face_current;

  /*
   * For each element we flaten all connectivities in one array
   */
  for (int ielt = 0; ielt < n_elt; ielt++) {

    for (int i_face = 0; i_face < n_face_elt; i_face++) {
      // _current_elmt_face_cell[ielt * n_face_elt + i_face    ] = *n_elt_current + ielt + 1;
      _current_elmt_face_cell[ielt * n_face_elt + i_face    ] = beg_gnum_elt_current + ielt + 1;
    }

    _current_elmt_face_vtx_idx[ielt * n_face_elt + 1]  = _current_elmt_face_vtx_idx[ielt * n_face_elt    ] + 4;
    _current_elmt_face_vtx_idx[ielt * n_face_elt + 2]  = _current_elmt_face_vtx_idx[ielt * n_face_elt + 1] + 3;
    _current_elmt_face_vtx_idx[ielt * n_face_elt + 3]  = _current_elmt_face_vtx_idx[ielt * n_face_elt + 2] + 3;
    _current_elmt_face_vtx_idx[ielt * n_face_elt + 4]  = _current_elmt_face_vtx_idx[ielt * n_face_elt + 3] + 3;
    _current_elmt_face_vtx_idx[ielt * n_face_elt + 5]  = _current_elmt_face_vtx_idx[ielt * n_face_elt + 4] + 3;

    // printf("Pyra::_current_elmt_face_vtx_idx[%i] = %i \n", _n_face_current + ielt * n_face_elt + 1, _current_elmt_face_vtx_idx[ielt * n_face_elt    ] + 4);
    // printf("Pyra::_current_elmt_face_vtx_idx[%i] = %i \n", _n_face_current + ielt * n_face_elt + 2, _current_elmt_face_vtx_idx[ielt * n_face_elt + 1] + 4);
    // printf("Pyra::_current_elmt_face_vtx_idx[%i] = %i \n", _n_face_current + ielt * n_face_elt + 3, _current_elmt_face_vtx_idx[ielt * n_face_elt + 2] + 4);
    // printf("Pyra::_current_elmt_face_vtx_idx[%i] = %i \n", _n_face_current + ielt * n_face_elt + 4, _current_elmt_face_vtx_idx[ielt * n_face_elt + 3] + 4);
    // printf("Pyra::_current_elmt_face_vtx_idx[%i] = %i \n", _n_face_current + ielt * n_face_elt + 5, _current_elmt_face_vtx_idx[ielt * n_face_elt + 4] + 4);

    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 0]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 3];
    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 1]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 2];
    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 2]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 1];
    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 3]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt    ];

    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 4]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 4];
    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 5]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt    ];
    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 6]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 1];

    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 7]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 4];
    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 8]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 1];
    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 9]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 2];

    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 10] = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 4];
    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 11] = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 2];
    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 12] = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 3];

    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 13] = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 4];
    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 14] = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 3];
    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 15] = connectivity_elmt_vtx[n_sum_vtx_elt * ielt    ];
  }

  *n_elt_current  += n_elt;
  *n_face_current += n_elt * n_face_elt;

}


/**
*
* \brief Decompose pyra cell_vtx connectivity to a flatten view of faces
*/
void
PDM_pyra_decomposes_edges
(
       int          n_elt,
       int         *n_elt_current,
       int         *n_edge_current,
       PDM_g_num_t  beg_gnum_elt_current,
       PDM_g_num_t  beg_gnum_edge_current,
 const PDM_g_num_t *connectivity_elmt_vtx,
       int         *elmt_edge_vtx_idx,
       PDM_g_num_t *elmt_edge_vtx,
       PDM_g_num_t *elmt_edge_cell,
       int         *elmt_cell_edge_idx,
       PDM_g_num_t *elmt_cell_edge
)
{
  PDM_UNUSED(elmt_cell_edge_idx);
  PDM_UNUSED(elmt_cell_edge);
  PDM_UNUSED(beg_gnum_elt_current);
  PDM_UNUSED(beg_gnum_edge_current);

  const int n_edge_elt        = 8;
  const int n_sum_vtx_edge    = 16;
  const int n_sum_vtx_elt     = 5;

  int _n_edge_current = *n_edge_current;
  int         *_current_elmt_edge_vtx_idx = elmt_edge_vtx_idx + _n_edge_current;
  PDM_g_num_t *_current_elmt_edge_vtx     = elmt_edge_vtx + elmt_edge_vtx_idx[_n_edge_current];
  PDM_g_num_t *_current_elmt_edge_cell    = elmt_edge_cell + _n_edge_current;

  /*
   * For each element we flaten all connectivities in one array
   */
  for (int ielt = 0; ielt < n_elt; ielt++) {

    for (int i_edge = 0; i_edge < n_edge_elt; i_edge++) {
      _current_elmt_edge_vtx_idx[ielt * n_edge_elt + i_edge + 1] = _current_elmt_edge_vtx_idx[ielt * n_edge_elt + i_edge] + 2;
      _current_elmt_edge_cell   [ielt * n_edge_elt + i_edge    ] = beg_gnum_elt_current + ielt + 1;
    }

    // E1 = N1 N2
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 0]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt    ];
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 1]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 1];

    // E2 = N2 N3
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 2]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 1];
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 3]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 2];

    // E3 = N3 N4
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 4]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 2];
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 5]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 3];

    // E4 = N4 N1
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 6]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 3];
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 7]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt    ];

    // E5 = N1 N5
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 8]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt    ];
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 9]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 4];

    // E6 = N2 N5
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 10] = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 1];
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 11] = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 4];

    // E7 = N3 N5
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 12] = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 2];
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 13] = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 4];

    // E8 = N4 N5
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 14] = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 3];
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 15] = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 4];

  }

  *n_elt_current  += n_elt;
  *n_edge_current += n_elt * n_edge_elt;

}

/**
*
* \brief Decompose tetra cell_vtx connectivity to a flatten view of faces
*/
void
PDM_prism_decomposes_faces
(
       int          n_elt,
       int         *n_elt_current,
       int         *n_face_current,
       PDM_g_num_t  beg_gnum_elt_current,
       PDM_g_num_t  beg_gnum_face_current,
 const PDM_g_num_t *connectivity_elmt_vtx,
       int         *elmt_face_vtx_idx,
       PDM_g_num_t *elmt_face_vtx,
       PDM_g_num_t *elmt_face_cell,
       int         *elmt_cell_face_idx,
       PDM_g_num_t *elmt_cell_face
)
{
  PDM_UNUSED(elmt_cell_face_idx);
  PDM_UNUSED(elmt_cell_face);
  PDM_UNUSED(beg_gnum_elt_current);
  PDM_UNUSED(beg_gnum_face_current);

  const int n_face_elt        = 5;
  const int n_sum_vtx_face    = 3*4 + 2*3;
  const int n_sum_vtx_elt     = 6;

  int _n_face_current = *n_face_current;
  int         *_current_elmt_face_vtx_idx = elmt_face_vtx_idx + _n_face_current;
  PDM_g_num_t *_current_elmt_face_vtx     = elmt_face_vtx + elmt_face_vtx_idx[_n_face_current];
  PDM_g_num_t *_current_elmt_face_cell    = elmt_face_cell + _n_face_current;

  /*
   * For each element we flaten all connectivities in one array
   */
  for (int ielt = 0; ielt < n_elt; ielt++) {

    for (int i_face = 0; i_face < n_face_elt; i_face++) {
      // _current_elmt_face_cell[ielt * n_face_elt + i_face    ] = *n_elt_current + ielt + 1;
      _current_elmt_face_cell[ielt * n_face_elt + i_face    ] = beg_gnum_elt_current + ielt + 1;
    }

    // elmt_cell_face_idx[ielt + ]
    // for (int i_face = 0; i_face < n_face_elt; i_face++) {
    //   _current_elmt_face_cell[ielt * n_face_elt + i_face    ] = *n_elt_current + ielt + 1;
    // }

    _current_elmt_face_vtx_idx[ielt * n_face_elt + 1]  = _current_elmt_face_vtx_idx[ielt * n_face_elt    ] + 3;
    _current_elmt_face_vtx_idx[ielt * n_face_elt + 2]  = _current_elmt_face_vtx_idx[ielt * n_face_elt + 1] + 3;
    _current_elmt_face_vtx_idx[ielt * n_face_elt + 3]  = _current_elmt_face_vtx_idx[ielt * n_face_elt + 2] + 4;
    _current_elmt_face_vtx_idx[ielt * n_face_elt + 4]  = _current_elmt_face_vtx_idx[ielt * n_face_elt + 3] + 4;
    _current_elmt_face_vtx_idx[ielt * n_face_elt + 5]  = _current_elmt_face_vtx_idx[ielt * n_face_elt + 4] + 4;

    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 0]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 2];
    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 1]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 1];
    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 2]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt    ];

    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 3]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 4];
    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 4]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 5];
    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 5]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 3];

    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 6]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 5];
    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 7]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 4];
    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 8]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 1];
    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 9]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 2];

    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 10] = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 4];
    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 11] = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 3];
    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 12] = connectivity_elmt_vtx[n_sum_vtx_elt * ielt    ];
    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 13] = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 1];

    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 14] = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 3];
    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 15] = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 5];
    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 16] = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 2];
    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 17] = connectivity_elmt_vtx[n_sum_vtx_elt * ielt    ];

  }

  *n_elt_current  += n_elt;
  *n_face_current += n_elt * n_face_elt;

}




/**
*
* \brief Decompose tetra cell_vtx connectivity to a flatten view of faces
*/
void
PDM_prism_decomposes_edges
(
       int          n_elt,
       int         *n_elt_current,
       int         *n_edge_current,
       PDM_g_num_t  beg_gnum_elt_current,
       PDM_g_num_t  beg_gnum_edge_current,
 const PDM_g_num_t *connectivity_elmt_vtx,
       int         *elmt_edge_vtx_idx,
       PDM_g_num_t *elmt_edge_vtx,
       PDM_g_num_t *elmt_edge_cell,
       int         *elmt_cell_edge_idx,
       PDM_g_num_t *elmt_cell_edge
)
{
  PDM_UNUSED(elmt_cell_edge_idx);
  PDM_UNUSED(elmt_cell_edge);
  PDM_UNUSED(beg_gnum_elt_current);
  PDM_UNUSED(beg_gnum_edge_current);

  const int n_edge_elt        = 9;
  const int n_sum_vtx_edge    = 18;
  const int n_sum_vtx_elt     = 6;

  int _n_edge_current = *n_edge_current;
  int         *_current_elmt_edge_vtx_idx = elmt_edge_vtx_idx + _n_edge_current;
  PDM_g_num_t *_current_elmt_edge_vtx     = elmt_edge_vtx + elmt_edge_vtx_idx[_n_edge_current];
  PDM_g_num_t *_current_elmt_edge_cell    = elmt_edge_cell + _n_edge_current;

  /*
   * For each element we flaten all connectivities in one array
   */
  for (int ielt = 0; ielt < n_elt; ielt++) {

    for (int i_edge = 0; i_edge < n_edge_elt; i_edge++) {
      _current_elmt_edge_vtx_idx[ielt * n_edge_elt + 1     ] = _current_elmt_edge_vtx_idx[ielt * n_edge_elt    ] + 2;
      _current_elmt_edge_cell   [ielt * n_edge_elt + i_edge] = beg_gnum_elt_current + ielt + 1;
    }

    // E1 = N1 N2
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 0]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt    ];
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 1]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 1];

    // E2 = N2 N3
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 2]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 1];
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 3]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 2];

    // E3 = N3 N1
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 4]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 2];
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 5]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt    ];

    // E4 = N1 N4
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 6]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt    ];
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 7]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 3];

    // E5 = N2 N5
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 8]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 1];
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 9]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 4];

    // E6 = N3 N6
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 10] = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 2];
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 11] = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 5];

    // E7 = N4 N5
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 12] = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 3];
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 13] = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 4];

    // E8 = N5 N6
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 14] = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 4];
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 15] = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 5];

    // E9 = N6 N4
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 16] = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 5];
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 17] = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 3];

  }

  *n_elt_current  += n_elt;
  *n_edge_current += n_elt * n_edge_elt;

}

/**
*
* \brief Decompose hexa cell_vtx connectivity to a flatten view of faces
*/
void
PDM_hexa_decomposes_faces
(
       int          n_elt,
       int         *n_elt_current,
       int         *n_face_current,
       PDM_g_num_t  beg_gnum_elt_current,
       PDM_g_num_t  beg_gnum_face_current,
 const PDM_g_num_t *connectivity_elmt_vtx,
       int         *elmt_face_vtx_idx,
       PDM_g_num_t *elmt_face_vtx,
       PDM_g_num_t *elmt_face_cell,
       int         *elmt_cell_face_idx,
       PDM_g_num_t *elmt_cell_face
)
{
  PDM_UNUSED(elmt_cell_face_idx);
  PDM_UNUSED(elmt_cell_face);
  PDM_UNUSED(beg_gnum_elt_current);
  PDM_UNUSED(beg_gnum_face_current);

  const int n_face_elt        = 6;
  const int n_sum_vtx_face    = 24;
  const int n_sum_vtx_elt     = 8;

  int _n_face_current = *n_face_current;
  int         *_current_elmt_face_vtx_idx = elmt_face_vtx_idx + _n_face_current;
  PDM_g_num_t *_current_elmt_face_vtx     = elmt_face_vtx + elmt_face_vtx_idx[_n_face_current];
  PDM_g_num_t *_current_elmt_face_cell    = elmt_face_cell + _n_face_current;


 printf("_n_face_current:: %i\n", _n_face_current);
 printf("elt_face_vtx_idx[%i]:: %i \n", _n_face_current, elmt_face_vtx_idx[_n_face_current]);
  /*
   * For each element we flaten all connectivities in one array
   */
  for (int ielt = 0; ielt < n_elt; ielt++) {

    /* Store the face_cell */
    for (int i_face = 0; i_face < n_face_elt; i_face++) {
      _current_elmt_face_vtx_idx[ielt * n_face_elt + i_face + 1] = _current_elmt_face_vtx_idx[ielt * n_face_elt + i_face] + 4;
      _current_elmt_face_cell   [ielt * n_face_elt + i_face    ] = beg_gnum_elt_current + ielt + 1;
    }

    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 0]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 3];
    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 1]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 2];
    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 2]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 1];
    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 3]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt    ];

    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 4]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 6];
    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 5]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 7];
    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 6]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 4];
    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 7]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 5];

    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 8]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 4];
    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 9]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 7];
    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 10] = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 3];
    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 11] = connectivity_elmt_vtx[n_sum_vtx_elt * ielt    ];

    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 12] = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 7];
    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 13] = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 6];
    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 14] = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 2];
    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 15] = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 3];

    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 16] = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 2];
    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 17] = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 6];
    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 18] = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 5];
    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 19] = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 1];

    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 20] = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 1];
    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 21] = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 5];
    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 22] = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 4];
    _current_elmt_face_vtx[n_sum_vtx_face * ielt + 23] = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 0];

  }

 *n_elt_current  += n_elt;
 *n_face_current += n_elt * n_face_elt;

}


/**
*
* \brief Decompose hexa cell_vtx connectivity to a flatten view of faces
*/
void
PDM_hexa_decomposes_edges
(
       int          n_elt,
       int         *n_elt_current,
       int         *n_edge_current,
       PDM_g_num_t  beg_gnum_elt_current,
       PDM_g_num_t  beg_gnum_edge_current,
 const PDM_g_num_t *connectivity_elmt_vtx,
       int         *elmt_edge_vtx_idx,
       PDM_g_num_t *elmt_edge_vtx,
       PDM_g_num_t *elmt_edge_cell,
       int         *elmt_cell_edge_idx,
       PDM_g_num_t *elmt_cell_edge
)
{
  PDM_UNUSED(elmt_cell_edge_idx);
  PDM_UNUSED(elmt_cell_edge);
  PDM_UNUSED(beg_gnum_elt_current);
  PDM_UNUSED(beg_gnum_edge_current);

  const int n_edge_elt        = 12;
  const int n_sum_vtx_edge    = 24; // 2 vtx * 12 edge
  const int n_sum_vtx_elt     = 8;

  int _n_edge_current = *n_edge_current;
  int         *_current_elmt_edge_vtx_idx = elmt_edge_vtx_idx + _n_edge_current;
  PDM_g_num_t *_current_elmt_edge_vtx     = elmt_edge_vtx + elmt_edge_vtx_idx[_n_edge_current];
  PDM_g_num_t *_current_elmt_edge_cell    = elmt_edge_cell + _n_edge_current;


  /*
   * For each element we flaten all connectivities in one array
   */
  for (int ielt = 0; ielt < n_elt; ielt++) {

    /* Store the edge_cell */
    for (int i_edge = 0; i_edge < n_edge_elt; i_edge++) {
      _current_elmt_edge_vtx_idx[ielt * n_edge_elt + i_edge + 1] = _current_elmt_edge_vtx_idx[ielt * n_edge_elt + i_edge] + 2;
      _current_elmt_edge_cell   [ielt * n_edge_elt + i_edge    ] = beg_gnum_elt_current + ielt + 1;
    }

    // E1 = N1 N2
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 0]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 0];
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 1]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 1];

    // E2 = N2 N3
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 2]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 1];
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 3]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 2];

    // E3 = N3 N4
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 4]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 2];
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 5]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 3];

    // E4 = N4 N1
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 6]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 3];
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 7]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt    ];

    // E5 = N1 N5
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 8]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt    ];
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 9]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 4];

    // E6 = N2 N6
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 10]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 1];
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 11]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 5];

    // E7 = N3 N7
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 12]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 2];
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 13]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 6];

    // E8 = N4 N8
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 14]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 3];
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 15]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 7];

    // E9 = N5 N6
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 16]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 4];
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 17]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 5];

    // E10 = N6 N7
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 18]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 5];
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 19]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 6];

    // E11 = N7 N8
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 20]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 6];
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 21]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 7];

    // E12 = N8 N5
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 22]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 7];
    _current_elmt_edge_vtx[n_sum_vtx_edge * ielt + 23]  = connectivity_elmt_vtx[n_sum_vtx_elt * ielt + 4];

  }

 *n_elt_current  += n_elt;
 *n_edge_current += n_elt * n_edge_elt;

}

/**
*
* \brief PDM_sections_decompose_faces
*
* \param [in]     mesh               Current mesh
* \param [inout]  elt_face_vtx_idx   Index of element faces connectivity (preallocated)
* \param [inout]  elt_face_vtx       Element faces connectivity (preallocated)
* \param [inout]  elmt_face_cell     Element faces connectivity (preallocated or NULL )
* \param [inout]  elmt_cell_face     Element faces connectivity (preallocated or NULL )
*
*/
void
PDM_sections_decompose_faces
(
  PDM_dmesh_nodal_t *dmesh_nodal,
  int               *elmt_face_vtx_idx,
  PDM_g_num_t       *elmt_face_vtx,
  PDM_g_num_t       *elmt_face_cell,
  int               *elmt_cell_face_idx,
  PDM_g_num_t       *elmt_cell_face
)
{

  // A faire : local_num_in_parent_element

  int n_elt_current  = 0;
  int n_face_current = 0;

  // The order of all following matter : following the global numebring
  for (int i_section = 0; i_section < dmesh_nodal->n_section_std; i_section++) {

    PDM_g_num_t beg_elmt_gnum = dmesh_nodal->sections_std[i_section]->distrib[dmesh_nodal->i_rank] + dmesh_nodal->section_distribution[i_section];
    PDM_g_num_t beg_face_gnum = 0; // Useless in this context
    switch (dmesh_nodal->sections_std[i_section]->t_elt) {
     case PDM_MESH_NODAL_POINT:
       abort();
       break;
     case PDM_MESH_NODAL_BAR2:
       abort();
       break;
     case PDM_MESH_NODAL_TRIA3:
       PDM_tri_decomposes_faces(dmesh_nodal->sections_std[i_section]->n_elt,
                                &n_elt_current,
                                &n_face_current,
                                beg_elmt_gnum,
                                beg_face_gnum,
                                dmesh_nodal->sections_std[i_section]->_connec,
                                elmt_face_vtx_idx,
                                elmt_face_vtx,
                                elmt_face_cell,
                                elmt_cell_face_idx,
                                elmt_cell_face);
       break;
     case PDM_MESH_NODAL_QUAD4:
       PDM_quad_decomposes_faces(dmesh_nodal->sections_std[i_section]->n_elt,
                                  &n_elt_current,
                                  &n_face_current,
                                  beg_elmt_gnum,
                                  beg_face_gnum,
                                  dmesh_nodal->sections_std[i_section]->_connec,
                                  elmt_face_vtx_idx,
                                  elmt_face_vtx,
                                  elmt_face_cell,
                                  elmt_cell_face_idx,
                                  elmt_cell_face);
       break;
     case PDM_MESH_NODAL_TETRA4:
       PDM_tetra_decomposes_faces(dmesh_nodal->sections_std[i_section]->n_elt,
                                  &n_elt_current,
                                  &n_face_current,
                                  beg_elmt_gnum,
                                  beg_face_gnum,
                                  dmesh_nodal->sections_std[i_section]->_connec,
                                  elmt_face_vtx_idx,
                                  elmt_face_vtx,
                                  elmt_face_cell,
                                  elmt_cell_face_idx,
                                  elmt_cell_face);
       break;
     case PDM_MESH_NODAL_PYRAMID5:
       PDM_pyra_decomposes_faces(dmesh_nodal->sections_std[i_section]->n_elt,
                                 &n_elt_current,
                                 &n_face_current,
                                 beg_elmt_gnum,
                                 beg_face_gnum,
                                 dmesh_nodal->sections_std[i_section]->_connec,
                                 elmt_face_vtx_idx,
                                 elmt_face_vtx,
                                 elmt_face_cell,
                                 elmt_cell_face_idx,
                                 elmt_cell_face);
       break;
     case PDM_MESH_NODAL_PRISM6:
       PDM_prism_decomposes_faces(dmesh_nodal->sections_std[i_section]->n_elt,
                                  &n_elt_current,
                                  &n_face_current,
                                  beg_elmt_gnum,
                                  beg_face_gnum,
                                  dmesh_nodal->sections_std[i_section]->_connec,
                                  elmt_face_vtx_idx,
                                  elmt_face_vtx,
                                  elmt_face_cell,
                                  elmt_cell_face_idx,
                                  elmt_cell_face);
       break;
     case PDM_MESH_NODAL_HEXA8:
       PDM_hexa_decomposes_faces(dmesh_nodal->sections_std[i_section]->n_elt,
                                 &n_elt_current,
                                 &n_face_current,
                                 beg_elmt_gnum,
                                 beg_face_gnum,
                                 dmesh_nodal->sections_std[i_section]->_connec,
                                 elmt_face_vtx_idx,
                                 elmt_face_vtx,
                                 elmt_face_cell,
                                 elmt_cell_face_idx,
                                 elmt_cell_face);
       break;
     default:
       PDM_error(__FILE__, __LINE__, 0, "Error PDM_sections_decompose_faces : Element type is not taking int account\n");
    }
  }

  // The order of all following matter : following the global numebring
  // for (int i_section = 0; i_section < dmesh_nodal->n_section_std_l1; i_section++) {

  //   PDM_g_num_t beg_elmt_gnum = dmesh_nodal->sections_std[i_section]->distrib[dmesh_nodal->i_rank] + dmesh_nodal->section_distribution_l1[i_section];
  //   PDM_g_num_t beg_face_gnum = 0; // Useless in this context
  //   // printf("section_std_l1 --> beg_elmt_gnum : "PDM_FMT_G_NUM" \n", beg_elmt_gnum);
  //   switch (dmesh_nodal->sections_std[i_section]->t_elt) {

  //    default:
  //      PDM_error(__FILE__, __LINE__, 0, "Error PDM_sections_decompose_faces : Element type is not taking int account\n");
  //   }
  // }

  // // The order of all following matter : following the global numebring
  // for (int i_section = 0; i_section < dmesh_nodal->n_section_std_l2; i_section++) {
  //   // PDM_g_num_t beg_elmt_gnum = dmesh_nodal->sections_std_l2[i_section]->distrib[dmesh_nodal->i_rank] + dmesh_nodal->section_distribution[i_section];
  //   // PDM_g_num_t beg_face_gnum = 0; // Useless in this context
  //   // printf("section_std_l2 --> beg_elmt_gnum : "PDM_FMT_G_NUM" \n", beg_elmt_gnum);
  //   switch (dmesh_nodal->sections_std_l2[i_section]->t_elt) {
  //    case PDM_MESH_NODAL_POINT:
  //      abort();
  //      break;
  //    case PDM_MESH_NODAL_BAR2:
  //      abort();
  //      break;     default:
  //      PDM_error(__FILE__, __LINE__, 0, "Error PDM_sections_decompose_faces : Element type is not taking int account\n");
  //   }
  // }

  /* Not implemented */
  if (dmesh_nodal->n_section_poly2d != 0) {
    // int n_section_poly2d = PDM_Handles_n_get  (dmesh_nodal->sections_poly2d);
    // const int *list_ind   = PDM_Handles_idx_get(dmesh_nodal->sections_poly2d);
    abort();
  }

  if (dmesh_nodal->n_section_poly3d != 0) {
    // int n_section_poly3d = PDM_Handles_n_get  (dmesh_nodal->sections_poly3d);
    // const int *list_ind   = PDM_Handles_idx_get(dmesh_nodal->sections_poly3d);
    abort();
  }
}

/**
*
* \brief PDM_sections_decompose_faces
*
* \param [in]     mesh               Current mesh
* \param [inout]  elmt_edge_vtx_idx  Index of element faces connectivity (preallocated)
* \param [inout]  elmt_edge_vtx      Element faces connectivity (preallocated)
* \param [inout]  elmt_edge_cell     Element faces connectivity (preallocated or NULL )
* \param [inout]  elmt_cell_edge     Element faces connectivity (preallocated or NULL )
*
*/
void
PDM_sections_decompose_edges
(
  PDM_dmesh_nodal_t *dmesh_nodal,
  int               *elmt_edge_vtx_idx,
  PDM_g_num_t       *elmt_edge_vtx,
  PDM_g_num_t       *elmt_edge_cell,
  int               *elmt_cell_edge_idx,
  PDM_g_num_t       *elmt_cell_edge
)
{
  PDM_UNUSED(dmesh_nodal);
  PDM_UNUSED(elmt_edge_vtx_idx);
  PDM_UNUSED(elmt_edge_vtx);
  PDM_UNUSED(elmt_edge_cell);
  PDM_UNUSED(elmt_cell_edge_idx);
  PDM_UNUSED(elmt_cell_edge);

  int n_elt_current  = 0;
  int n_dedge_current = 0;

  for (int i_section = 0; i_section < dmesh_nodal->n_section_std; i_section++) {

    PDM_g_num_t beg_elmt_gnum = dmesh_nodal->sections_std[i_section]->distrib[dmesh_nodal->i_rank] + dmesh_nodal->section_distribution[i_section];
    PDM_g_num_t beg_edge_gnum = 0; // Useless in this context

    switch (dmesh_nodal->sections_std[i_section]->t_elt) {
     case PDM_MESH_NODAL_POINT:
       abort();
       break;
     case PDM_MESH_NODAL_BAR2:
       PDM_bar_decomposes_edges(dmesh_nodal->sections_std[i_section]->n_elt,
                                &n_elt_current,
                                &n_dedge_current,
                                beg_elmt_gnum,
                                beg_edge_gnum,
                                dmesh_nodal->sections_std[i_section]->_connec,
                                elmt_edge_vtx_idx,
                                elmt_edge_vtx,
                                elmt_edge_cell,
                                elmt_cell_edge_idx,
                                elmt_cell_edge);
       break;
     case PDM_MESH_NODAL_TRIA3:
       PDM_tri_decomposes_edges(dmesh_nodal->sections_std[i_section]->n_elt,
                                 &n_elt_current,
                                 &n_dedge_current,
                                 beg_elmt_gnum,
                                 beg_edge_gnum,
                                 dmesh_nodal->sections_std[i_section]->_connec,
                                 elmt_edge_vtx_idx,
                                 elmt_edge_vtx,
                                 elmt_edge_cell,
                                 elmt_cell_edge_idx,
                                 elmt_cell_edge);
       break;
     case PDM_MESH_NODAL_QUAD4:
       PDM_quad_decomposes_edges(dmesh_nodal->sections_std[i_section]->n_elt,
                                 &n_elt_current,
                                 &n_dedge_current,
                                 beg_elmt_gnum,
                                 beg_edge_gnum,
                                 dmesh_nodal->sections_std[i_section]->_connec,
                                 elmt_edge_vtx_idx,
                                 elmt_edge_vtx,
                                 elmt_edge_cell,
                                 elmt_cell_edge_idx,
                                 elmt_cell_edge);
       break;
     case PDM_MESH_NODAL_TETRA4:
       PDM_tetra_decomposes_edges(dmesh_nodal->sections_std[i_section]->n_elt,
                                  &n_elt_current,
                                  &n_dedge_current,
                                  beg_elmt_gnum,
                                  beg_edge_gnum,
                                  dmesh_nodal->sections_std[i_section]->_connec,
                                  elmt_edge_vtx_idx,
                                  elmt_edge_vtx,
                                  elmt_edge_cell,
                                  elmt_cell_edge_idx,
                                  elmt_cell_edge);
       break;
     case PDM_MESH_NODAL_PYRAMID5:
       PDM_pyra_decomposes_edges(dmesh_nodal->sections_std[i_section]->n_elt,
                                 &n_elt_current,
                                 &n_dedge_current,
                                 beg_elmt_gnum,
                                 beg_edge_gnum,
                                 dmesh_nodal->sections_std[i_section]->_connec,
                                 elmt_edge_vtx_idx,
                                 elmt_edge_vtx,
                                 elmt_edge_cell,
                                 elmt_cell_edge_idx,
                                 elmt_cell_edge);
       break;
     case PDM_MESH_NODAL_PRISM6:
       PDM_prism_decomposes_edges(dmesh_nodal->sections_std[i_section]->n_elt,
                                  &n_elt_current,
                                  &n_dedge_current,
                                  beg_elmt_gnum,
                                  beg_edge_gnum,
                                  dmesh_nodal->sections_std[i_section]->_connec,
                                  elmt_edge_vtx_idx,
                                  elmt_edge_vtx,
                                  elmt_edge_cell,
                                  elmt_cell_edge_idx,
                                  elmt_cell_edge);
       break;
     case PDM_MESH_NODAL_HEXA8:
       PDM_hexa_decomposes_edges(dmesh_nodal->sections_std[i_section]->n_elt,
                                 &n_elt_current,
                                 &n_dedge_current,
                                 beg_elmt_gnum,
                                 beg_edge_gnum,
                                 dmesh_nodal->sections_std[i_section]->_connec,
                                 elmt_edge_vtx_idx,
                                 elmt_edge_vtx,
                                 elmt_edge_cell,
                                 elmt_cell_edge_idx,
                                 elmt_cell_edge);
       break;
     default:
       PDM_error(__FILE__, __LINE__, 0, "Error PDM_sections_decompose_edges : Element type is not taking int account\n");
    }
  }

  // The order of all following matter : following the global numebring
  // for (int i_section = 0; i_section < dmesh_nodal->n_section_std_l1; i_section++) {

  //   PDM_g_num_t beg_elmt_gnum = dmesh_nodal->sections_std[i_section]->distrib[dmesh_nodal->i_rank] + dmesh_nodal->section_distribution_l1[i_section];
  //   PDM_g_num_t beg_edge_gnum = 0; // Useless in this context

  //   // printf("section_std_l1 --> beg_elmt_gnum : "PDM_FMT_G_NUM" \n", beg_elmt_gnum);
  //   switch (dmesh_nodal->sections_std[i_section]->t_elt) {
  //    case PDM_MESH_NODAL_TRIA3:
  //      PDM_tri_decomposes_edges(dmesh_nodal->sections_std[i_section]->n_elt,
  //                               &n_elt_current,
  //                               &n_dedge_current,
  //                               beg_elmt_gnum,
  //                               beg_edge_gnum,
  //                               dmesh_nodal->sections_std[i_section]->_connec,
  //                               elmt_edge_vtx_idx,
  //                               elmt_edge_vtx,
  //                               elmt_edge_cell,
  //                               elmt_cell_edge_idx,
  //                               elmt_cell_edge);
  //      break;
  //    case PDM_MESH_NODAL_QUAD4:
  //      PDM_quad_decomposes_edges(dmesh_nodal->sections_std[i_section]->n_elt,
  //                                &n_elt_current,
  //                                &n_dedge_current,
  //                                beg_elmt_gnum,
  //                                beg_edge_gnum,
  //                                dmesh_nodal->sections_std[i_section]->_connec,
  //                                elmt_edge_vtx_idx,
  //                                elmt_edge_vtx,
  //                                elmt_edge_cell,
  //                                elmt_cell_edge_idx,
  //                                elmt_cell_edge);
  //      break;
  //    default:
  //      PDM_error(__FILE__, __LINE__, 0, "Error PDM_sections_decompose_faces : Element type is not taking int account\n");
  //   }
  // }


  /* Not implemented */
  if (dmesh_nodal->n_section_poly2d != 0) {
    // int n_section_poly2d = PDM_Handles_n_get  (dmesh_nodal->sections_poly2d);
    // const int *list_ind   = PDM_Handles_idx_get(dmesh_nodal->sections_poly2d);
    abort();
  }

  if (dmesh_nodal->n_section_poly3d != 0) {
    // int n_section_poly3d = PDM_Handles_n_get  (dmesh_nodal->sections_poly3d);
    // const int *list_ind   = PDM_Handles_idx_get(dmesh_nodal->sections_poly3d);
    abort();
  }

}

