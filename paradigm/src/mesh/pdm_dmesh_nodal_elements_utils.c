
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

/**
*
* \brief _section_elt_faces_get
*
* \param [in]     mesh               Current mesh
* \param [in]     id_section         Section identifier
* \param [inout]  elt_face_vtx_idx   Index of element faces connectivity (preallocated)
* \param [inout]  elt_face_vtx       Element faces connectivity (preallocated)
*
*/

// static int
// _section_size_elt_faces_get
// (
//   PDM_DMesh_nodal_t *mesh,
//   int               *s_elt_face_vtx_idx,
//   int               *s_elt_face_vtx,
//   int               *s_elt_face_cell
// )
// {

//  int _s_elt_face_vtx_idx = 0;
//  int _s_elt_face_vtx = 0;

//  int n_sections_std  = PDM_Handles_n_get (mesh->sections_std);
//  const int *list_ind = PDM_Handles_idx_get (mesh->sections_std);

//  for (int i = 0; i < n_sections_std; i++) {
//    PDM_DMesh_nodal_section_std_t *section =
//      (PDM_DMesh_nodal_section_std_t *) PDM_Handles_get (mesh->sections_std, list_ind[i]);
//    int n_face_elt = 0;
//    int n_sum_vtx_face = 0;

//    switch (section->t_elt) {
//    case PDM_MESH_NODAL_TRIA3:
//      n_face_elt = 3;
//      n_sum_vtx_face = 6;
//      break;
//    case PDM_MESH_NODAL_TETRA4:
//      n_face_elt = 4;
//      n_sum_vtx_face = 12;
//      break;
//    case PDM_MESH_NODAL_QUAD4:
//      n_face_elt = 4;
//      n_sum_vtx_face = 8;
//      break;
//    case PDM_MESH_NODAL_HEXA8:
//      n_face_elt = 6;
//      n_sum_vtx_face = 24;
//      break;
//    case PDM_MESH_NODAL_PYRAMID5:
//      n_face_elt = 5;
//      n_sum_vtx_face = 16;
//      break;
//    case PDM_MESH_NODAL_PRISM6:
//      n_face_elt = 5;
//      n_sum_vtx_face = 18;
//      break;
//    default:
//      PDM_error(__FILE__, __LINE__, 0, "Error _section_size_elt_faces_get : Element type is not taking int account\n");
//      exit(EXIT_FAILURE);
//    }

//    _s_elt_face_vtx_idx += section->n_elt * n_face_elt;
//    _s_elt_face_vtx += section->n_elt * n_sum_vtx_face;
//  }

//  int n_sections_poly2d = PDM_Handles_n_get (mesh->sections_poly2d);
//  list_ind = PDM_Handles_idx_get (mesh->sections_poly2d);

//  for (int i = 0; i < n_sections_poly2d; i++) {
//    PDM_DMesh_nodal_section_poly2d_t *section =
//      (PDM_DMesh_nodal_section_poly2d_t *) PDM_Handles_get (mesh->sections_poly2d, list_ind[i]);
//    _s_elt_face_vtx_idx += section->_connec_idx[section->n_elt];
//    _s_elt_face_vtx += 2 * section->_connec_idx[section->n_elt];
//  }

//  int n_sections_poly3d = PDM_Handles_n_get (mesh->sections_poly3d);
//  list_ind = PDM_Handles_idx_get (mesh->sections_poly3d);

//  for (int i = 0; i < n_sections_poly3d; i++) {
//    PDM_DMesh_nodal_section_poly3d_t *section =
//      (PDM_DMesh_nodal_section_poly3d_t *) PDM_Handles_get (mesh->sections_poly3d, list_ind[i]);
//    _s_elt_face_vtx_idx += section->n_face;
//    _s_elt_face_vtx += section->_facvtx[section->_facvtx_idx[section->n_face]];
//  }

//  *s_elt_face_cell    =  _s_elt_face_vtx_idx;
//  *s_elt_face_vtx_idx = _s_elt_face_vtx_idx + 1;
//  *s_elt_face_vtx     = _s_elt_face_vtx + 1;


//  return *s_elt_face_vtx - 1;

// }

// /**
// *
// * \brief _section_elt_faces_get
// *
// * \param [in]     mesh               Current mesh
// * \param [in]     id_section         Section identifier
// * \param [inout]  elt_face_vtx_idx   Index of element faces connectivity (preallocated)
// * \param [inout]  elt_face_vtx       Element faces connectivity (preallocated)
// *
// */

// static void
// _section_elt_faces_add
// (
//       PDM_DMesh_nodal_t *mesh,
// const int                id_section,
//       int               *n_elt_current,
//       int               *n_face_current,
//       int               *elt_face_vtx_idx,
//       PDM_g_num_t       *elt_face_vtx,
//       PDM_g_num_t       *elt_face_cell
// )
// {
//  int _n_face_current = *n_face_current;

//  int         *_current_elmt_face_vtx_idx = elt_face_vtx_idx + _n_face_current;
//  PDM_g_num_t *_current_elmt_face_vtx     = elt_face_vtx + elt_face_vtx_idx[_n_face_current];
//  PDM_g_num_t *_current_elmt_face_cell    = elt_face_cell + _n_face_current;

//  if (id_section < PDM_BLOCK_ID_BLOCK_POLY2D) {

//    PDM_DMesh_nodal_section_std_t *section =
//            (PDM_DMesh_nodal_section_std_t *) PDM_Handles_get (mesh->sections_std, id_section);

//    if (section == NULL) {
//      PDM_error (__FILE__, __LINE__, 0, "Bad section identifier\n");
//    }

//    switch (section->t_elt) {

//    case PDM_MESH_NODAL_TRIA3:
//      {
//        const int n_face_elt        = 3;
//        const int n_sum_vtx_face    = 6;
//        const int n_sum_vtx_elt     = 3;

//        for (int ielt = 0; ielt < section->n_elt; ielt++) {
//          for (int i_face = 0; i_face < n_face_elt; i_face++) {
//            _current_elmt_face_vtx_idx[ielt * n_face_elt + i_face + 1] =
//              _current_elmt_face_vtx_idx[ielt * n_face_elt + i_face] + 3;
//            _current_elmt_face_cell[ielt * n_face_elt + i_face    ] = *n_elt_current + ielt + 1;
//          }

//          _current_elmt_face_vtx[n_sum_vtx_face * ielt + 0] = section->_connec[n_sum_vtx_elt * ielt    ];
//          _current_elmt_face_vtx[n_sum_vtx_face * ielt + 1] = section->_connec[n_sum_vtx_elt * ielt + 1];

//          _current_elmt_face_vtx[n_sum_vtx_face * ielt + 2] = section->_connec[n_sum_vtx_elt * ielt + 1];
//          _current_elmt_face_vtx[n_sum_vtx_face * ielt + 3] = section->_connec[n_sum_vtx_elt * ielt + 2];

//          _current_elmt_face_vtx[n_sum_vtx_face * ielt + 4] = section->_connec[n_sum_vtx_elt * ielt + 2];
//          _current_elmt_face_vtx[n_sum_vtx_face * ielt + 5] = section->_connec[n_sum_vtx_elt * ielt + 3];
//        }

//        *n_elt_current += section->n_elt;
//        *n_face_current += section->n_elt * n_face_elt;
//      }
//      break;

//    case PDM_MESH_NODAL_TETRA4:
//      {
//        const int n_face_elt        = 4;
//        const int n_sum_vtx_face    = 12;
//        const int n_sum_vtx_elt     = 4;

//        for (int ielt = 0; ielt < section->n_elt; ielt++) {
//          for (int i_face = 0; i_face < n_face_elt; i_face++) {
//            _current_elmt_face_vtx_idx[ielt * n_face_elt + i_face + 1] =
//              _current_elmt_face_vtx_idx[ielt * n_face_elt + i_face] + 3;
//            _current_elmt_face_cell[ielt * n_face_elt + i_face    ] = *n_elt_current + ielt + 1;
//          }

//          _current_elmt_face_vtx[n_sum_vtx_face * ielt + 0]  = section->_connec[n_sum_vtx_elt * ielt    ];
//          _current_elmt_face_vtx[n_sum_vtx_face * ielt + 1]  = section->_connec[n_sum_vtx_elt * ielt + 1];
//          _current_elmt_face_vtx[n_sum_vtx_face * ielt + 2]  = section->_connec[n_sum_vtx_elt * ielt + 2];

//          _current_elmt_face_vtx[n_sum_vtx_face * ielt + 3]  = section->_connec[n_sum_vtx_elt * ielt    ];
//          _current_elmt_face_vtx[n_sum_vtx_face * ielt + 4]  = section->_connec[n_sum_vtx_elt * ielt + 3];
//          _current_elmt_face_vtx[n_sum_vtx_face * ielt + 5]  = section->_connec[n_sum_vtx_elt * ielt + 1];

//          _current_elmt_face_vtx[n_sum_vtx_face * ielt + 6]  = section->_connec[n_sum_vtx_elt * ielt    ];
//          _current_elmt_face_vtx[n_sum_vtx_face * ielt + 7]  = section->_connec[n_sum_vtx_elt * ielt + 2];
//          _current_elmt_face_vtx[n_sum_vtx_face * ielt + 8]  = section->_connec[n_sum_vtx_elt * ielt + 3];

//          _current_elmt_face_vtx[n_sum_vtx_face * ielt + 9]  = section->_connec[n_sum_vtx_elt * ielt + 1];
//          _current_elmt_face_vtx[n_sum_vtx_face * ielt + 10] = section->_connec[n_sum_vtx_elt * ielt + 3];
//          _current_elmt_face_vtx[n_sum_vtx_face * ielt + 11] = section->_connec[n_sum_vtx_elt * ielt + 2];

//        }

//        *n_elt_current += section->n_elt;
//        *n_face_current += section->n_elt * n_face_elt;

//    }
//      break;
//    case PDM_MESH_NODAL_QUAD4:
//      {
//        const int n_face_elt        = 4;
//        const int n_sum_vtx_face    = 8;
//        const int n_sum_vtx_elt     = 4;

//        for (int ielt = 0; ielt < section->n_elt; ielt++) {
//          for (int i_face = 0; i_face < n_face_elt; i_face++) {
//            _current_elmt_face_vtx_idx[ielt * n_face_elt + i_face + 1] =
//              _current_elmt_face_vtx_idx[ielt * n_face_elt + i_face] + 3;
//            _current_elmt_face_cell[ielt * n_face_elt + i_face    ] = *n_elt_current + ielt + 1;
//          }

//          _current_elmt_face_vtx[n_sum_vtx_face * ielt + 0] = section->_connec[n_sum_vtx_elt * ielt    ];
//          _current_elmt_face_vtx[n_sum_vtx_face * ielt + 1] = section->_connec[n_sum_vtx_elt * ielt + 1];

//          _current_elmt_face_vtx[n_sum_vtx_face * ielt + 2] = section->_connec[n_sum_vtx_elt * ielt + 1];
//          _current_elmt_face_vtx[n_sum_vtx_face * ielt + 3] = section->_connec[n_sum_vtx_elt * ielt + 2];

//          _current_elmt_face_vtx[n_sum_vtx_face * ielt + 4] = section->_connec[n_sum_vtx_elt * ielt + 2];
//          _current_elmt_face_vtx[n_sum_vtx_face * ielt + 5] = section->_connec[n_sum_vtx_elt * ielt + 3];
//        }

//        *n_elt_current  += section->n_elt;
//        *n_face_current += section->n_elt * n_face_elt;

//      }
//      break;

//    case PDM_MESH_NODAL_HEXA8:
//      {
//        const int n_face_elt        = 6;
//        const int n_sum_vtx_face    = 24;
//        const int n_sum_vtx_elt     = 8;

//        for (int ielt = 0; ielt < section->n_elt; ielt++) {
//          for (int i_face = 0; i_face < n_face_elt; i_face++) {
//            _current_elmt_face_vtx_idx[ielt * n_face_elt + i_face + 1] =
//              _current_elmt_face_vtx_idx[ielt * n_face_elt + i_face] + 4;
//            _current_elmt_face_cell[ielt * n_face_elt + i_face    ] = *n_elt_current + ielt + 1;
//          }

//          _current_elmt_face_vtx[n_sum_vtx_face * ielt + 0]  = section->_connec[n_sum_vtx_elt * ielt + 3];
//          _current_elmt_face_vtx[n_sum_vtx_face * ielt + 1]  = section->_connec[n_sum_vtx_elt * ielt + 2];
//          _current_elmt_face_vtx[n_sum_vtx_face * ielt + 2]  = section->_connec[n_sum_vtx_elt * ielt + 1];
//          _current_elmt_face_vtx[n_sum_vtx_face * ielt + 3]  = section->_connec[n_sum_vtx_elt * ielt    ];

//          _current_elmt_face_vtx[n_sum_vtx_face * ielt + 4]  = section->_connec[n_sum_vtx_elt * ielt + 6];
//          _current_elmt_face_vtx[n_sum_vtx_face * ielt + 5]  = section->_connec[n_sum_vtx_elt * ielt + 7];
//          _current_elmt_face_vtx[n_sum_vtx_face * ielt + 6]  = section->_connec[n_sum_vtx_elt * ielt + 4];
//          _current_elmt_face_vtx[n_sum_vtx_face * ielt + 7]  = section->_connec[n_sum_vtx_elt * ielt + 5];

//          _current_elmt_face_vtx[n_sum_vtx_face * ielt + 8]  = section->_connec[n_sum_vtx_elt * ielt + 4];
//          _current_elmt_face_vtx[n_sum_vtx_face * ielt + 9]  = section->_connec[n_sum_vtx_elt * ielt + 7];
//          _current_elmt_face_vtx[n_sum_vtx_face * ielt + 10] = section->_connec[n_sum_vtx_elt * ielt + 3];
//          _current_elmt_face_vtx[n_sum_vtx_face * ielt + 11] = section->_connec[n_sum_vtx_elt * ielt    ];

//          _current_elmt_face_vtx[n_sum_vtx_face * ielt + 12] = section->_connec[n_sum_vtx_elt * ielt + 7];
//          _current_elmt_face_vtx[n_sum_vtx_face * ielt + 13] = section->_connec[n_sum_vtx_elt * ielt + 6];
//          _current_elmt_face_vtx[n_sum_vtx_face * ielt + 14] = section->_connec[n_sum_vtx_elt * ielt + 2];
//          _current_elmt_face_vtx[n_sum_vtx_face * ielt + 15] = section->_connec[n_sum_vtx_elt * ielt + 3];

//          _current_elmt_face_vtx[n_sum_vtx_face * ielt + 16] = section->_connec[n_sum_vtx_elt * ielt + 2];
//          _current_elmt_face_vtx[n_sum_vtx_face * ielt + 17] = section->_connec[n_sum_vtx_elt * ielt + 6];
//          _current_elmt_face_vtx[n_sum_vtx_face * ielt + 18] = section->_connec[n_sum_vtx_elt * ielt + 5];
//          _current_elmt_face_vtx[n_sum_vtx_face * ielt + 19] = section->_connec[n_sum_vtx_elt * ielt + 1];

//          _current_elmt_face_vtx[n_sum_vtx_face * ielt + 20] = section->_connec[n_sum_vtx_elt * ielt + 1];
//          _current_elmt_face_vtx[n_sum_vtx_face * ielt + 21] = section->_connec[n_sum_vtx_elt * ielt + 5];
//          _current_elmt_face_vtx[n_sum_vtx_face * ielt + 22] = section->_connec[n_sum_vtx_elt * ielt + 4];
//          _current_elmt_face_vtx[n_sum_vtx_face * ielt + 23] = section->_connec[n_sum_vtx_elt * ielt + 0];
//        }

//        *n_elt_current  += section->n_elt;
//        *n_face_current += section->n_elt * n_face_elt;

//      }
//      break;

//    case PDM_MESH_NODAL_PYRAMID5:
//      {
//        const int n_face_elt        = 5;
//        const int n_sum_vtx_face    = 1*4 + 4*3;
//        const int n_sum_vtx_elt     = 5;

//        elt_face_vtx_idx[0] = 0;

//        for (int ielt = 0; ielt < section->n_elt; ielt++) {

//          for (int i_face = 0; i_face < n_face_elt; i_face++) {
//            _current_elmt_face_cell[ielt * n_face_elt + i_face    ] = *n_elt_current + ielt + 1;
//          }

//          _current_elmt_face_vtx_idx[ielt * n_face_elt + 1]  = elt_face_vtx_idx[ielt * n_face_elt    ] + 4;
//          _current_elmt_face_vtx_idx[ielt * n_face_elt + 2]  = elt_face_vtx_idx[ielt * n_face_elt + 1] + 3;
//          _current_elmt_face_vtx_idx[ielt * n_face_elt + 3]  = elt_face_vtx_idx[ielt * n_face_elt + 2] + 3;
//          _current_elmt_face_vtx_idx[ielt * n_face_elt + 4]  = elt_face_vtx_idx[ielt * n_face_elt + 3] + 3;
//          _current_elmt_face_vtx_idx[ielt * n_face_elt + 5]  = elt_face_vtx_idx[ielt * n_face_elt + 4] + 3;

//          _current_elmt_face_vtx[n_sum_vtx_face * ielt + 0]  = section->_connec[n_sum_vtx_elt * ielt + 3];
//          _current_elmt_face_vtx[n_sum_vtx_face * ielt + 1]  = section->_connec[n_sum_vtx_elt * ielt + 2];
//          _current_elmt_face_vtx[n_sum_vtx_face * ielt + 2]  = section->_connec[n_sum_vtx_elt * ielt + 1];
//          _current_elmt_face_vtx[n_sum_vtx_face * ielt + 3]  = section->_connec[n_sum_vtx_elt * ielt    ];

//          _current_elmt_face_vtx[n_sum_vtx_face * ielt + 4]  = section->_connec[n_sum_vtx_elt * ielt + 4];
//          _current_elmt_face_vtx[n_sum_vtx_face * ielt + 5]  = section->_connec[n_sum_vtx_elt * ielt    ];
//          _current_elmt_face_vtx[n_sum_vtx_face * ielt + 6]  = section->_connec[n_sum_vtx_elt * ielt + 1];

//          _current_elmt_face_vtx[n_sum_vtx_face * ielt + 7]  = section->_connec[n_sum_vtx_elt * ielt + 4];
//          _current_elmt_face_vtx[n_sum_vtx_face * ielt + 8]  = section->_connec[n_sum_vtx_elt * ielt + 1];
//          _current_elmt_face_vtx[n_sum_vtx_face * ielt + 9]  = section->_connec[n_sum_vtx_elt * ielt + 2];

//          _current_elmt_face_vtx[n_sum_vtx_face * ielt + 10] = section->_connec[n_sum_vtx_elt * ielt + 4];
//          _current_elmt_face_vtx[n_sum_vtx_face * ielt + 11] = section->_connec[n_sum_vtx_elt * ielt + 2];
//          _current_elmt_face_vtx[n_sum_vtx_face * ielt + 12] = section->_connec[n_sum_vtx_elt * ielt + 3];

//          _current_elmt_face_vtx[n_sum_vtx_face * ielt + 13] = section->_connec[n_sum_vtx_elt * ielt + 4];
//          _current_elmt_face_vtx[n_sum_vtx_face * ielt + 14] = section->_connec[n_sum_vtx_elt * ielt + 3];
//          _current_elmt_face_vtx[n_sum_vtx_face * ielt + 15] = section->_connec[n_sum_vtx_elt * ielt    ];

//        }

//        *n_elt_current  += section->n_elt;
//        *n_face_current += section->n_elt * n_face_elt;

//      }
//      break;

//    case PDM_MESH_NODAL_PRISM6:
//      {
//        const int n_face_elt        = 5;
//        const int n_sum_vtx_face    = 3*4 + 2*3;
//        const int n_sum_vtx_elt     = 6;

//        for (int ielt = 0; ielt < section->n_elt; ielt++) {

//          for (int i_face = 0; i_face < n_face_elt; i_face++) {
//            _current_elmt_face_cell[ielt * n_face_elt + i_face    ] = *n_elt_current + ielt + 1;
//          }

//          _current_elmt_face_vtx_idx[ielt * n_face_elt + 1]  = elt_face_vtx_idx[ielt * n_face_elt    ] + 3;
//          _current_elmt_face_vtx_idx[ielt * n_face_elt + 2]  = elt_face_vtx_idx[ielt * n_face_elt + 1] + 3;
//          _current_elmt_face_vtx_idx[ielt * n_face_elt + 3]  = elt_face_vtx_idx[ielt * n_face_elt + 2] + 4;
//          _current_elmt_face_vtx_idx[ielt * n_face_elt + 4]  = elt_face_vtx_idx[ielt * n_face_elt + 3] + 4;
//          _current_elmt_face_vtx_idx[ielt * n_face_elt + 5]  = elt_face_vtx_idx[ielt * n_face_elt + 4] + 4;

//          _current_elmt_face_vtx[n_sum_vtx_face * ielt + 0]  = section->_connec[n_sum_vtx_elt * ielt + 2];
//          _current_elmt_face_vtx[n_sum_vtx_face * ielt + 1]  = section->_connec[n_sum_vtx_elt * ielt + 1];
//          _current_elmt_face_vtx[n_sum_vtx_face * ielt + 2]  = section->_connec[n_sum_vtx_elt * ielt    ];

//          _current_elmt_face_vtx[n_sum_vtx_face * ielt + 3]  = section->_connec[n_sum_vtx_elt * ielt + 4];
//          _current_elmt_face_vtx[n_sum_vtx_face * ielt + 4]  = section->_connec[n_sum_vtx_elt * ielt + 5];
//          _current_elmt_face_vtx[n_sum_vtx_face * ielt + 5]  = section->_connec[n_sum_vtx_elt * ielt + 3];

//          _current_elmt_face_vtx[n_sum_vtx_face * ielt + 6]  = section->_connec[n_sum_vtx_elt * ielt + 5];
//          _current_elmt_face_vtx[n_sum_vtx_face * ielt + 7]  = section->_connec[n_sum_vtx_elt * ielt + 4];
//          _current_elmt_face_vtx[n_sum_vtx_face * ielt + 8]  = section->_connec[n_sum_vtx_elt * ielt + 1];
//          _current_elmt_face_vtx[n_sum_vtx_face * ielt + 9]  = section->_connec[n_sum_vtx_elt * ielt + 2];

//          _current_elmt_face_vtx[n_sum_vtx_face * ielt + 10] = section->_connec[n_sum_vtx_elt * ielt + 4];
//          _current_elmt_face_vtx[n_sum_vtx_face * ielt + 11] = section->_connec[n_sum_vtx_elt * ielt + 3];
//          _current_elmt_face_vtx[n_sum_vtx_face * ielt + 12] = section->_connec[n_sum_vtx_elt * ielt    ];
//          _current_elmt_face_vtx[n_sum_vtx_face * ielt + 13] = section->_connec[n_sum_vtx_elt * ielt + 1];

//          _current_elmt_face_vtx[n_sum_vtx_face * ielt + 14] = section->_connec[n_sum_vtx_elt * ielt + 3];
//          _current_elmt_face_vtx[n_sum_vtx_face * ielt + 15] = section->_connec[n_sum_vtx_elt * ielt + 5];
//          _current_elmt_face_vtx[n_sum_vtx_face * ielt + 16] = section->_connec[n_sum_vtx_elt * ielt + 2];
//          _current_elmt_face_vtx[n_sum_vtx_face * ielt + 17] = section->_connec[n_sum_vtx_elt * ielt    ];

//        }

//        *n_elt_current  += section->n_elt;
//        *n_face_current += section->n_elt * n_face_elt;

//      }
//      break;

//      default:
//        break;
//    }
//  }

//  else if (id_section < PDM_BLOCK_ID_BLOCK_POLY3D) {

//    int _id = id_section - PDM_BLOCK_ID_BLOCK_POLY2D;
//    PDM_DMesh_nodal_section_poly2d_t *section =
//            (PDM_DMesh_nodal_section_poly2d_t *) PDM_Handles_get (mesh->sections_poly2d, _id);

//    if (section == NULL) {
//      PDM_error (__FILE__, __LINE__, 0, "Bad section identifier\n");
//    }

//    int idx = 0;

//    for (int ielt = 0; ielt < section->n_elt; ielt++) {
//      int n_face_elt = section->_connec_idx[section->n_elt];
//      *n_face_current += n_face_elt;
//      int idx2 = section->_connec_idx[ielt];
//      for (int i_face = 0; i_face < n_face_elt; i_face++) {
//        _current_elmt_face_vtx_idx[idx + 1]  = _current_elmt_face_vtx_idx[idx] + 2;
//        int inext = (i_face + 1) % n_face_elt;
//        _current_elmt_face_vtx[2 * idx    ]  = section->_connec[idx2 + i_face];
//        _current_elmt_face_vtx[2 * idx + 1]  = section->_connec[idx2 + inext];
//        _current_elmt_face_cell[idx   ]  = *n_elt_current + ielt + 1;
//        idx += 1;
//      }
//    }

//    *n_elt_current += section->n_elt;

//  }

//  else {

//    PDM_error (__FILE__, __LINE__, 0, "PDM_BLOCK_ID_BLOCK_POLY3D : Not implemented yet\n");

//    //TODO: Compliqué car il faut faire des échanges Block_to_part pour recuperer les
//    // definitions des faces
//    // Il faut redupliquer les faces et les stocker comme pour les autres types de
//    // section
//    // Attention : Dans le cas d'un maillage avec une seule section poly3d, il ne faut
//    // rien faire.et ne pas passer dans cette fonction

//  }

// }


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
  PDM_DMesh_nodal_t *mesh,
  int               *s_elt_face_vtx_idx,
  int               *s_elt_face_vtx,
  int               *s_elt_face_cell
)
{

 int _s_elt_face_vtx_idx = 0;
 int _s_elt_face_vtx     = 0;

 int n_sections_std  = PDM_Handles_n_get (mesh->sections_std);
 const int *list_ind = PDM_Handles_idx_get (mesh->sections_std);

 for (int i = 0; i < n_sections_std; i++) {
   PDM_DMesh_nodal_section_std_t *section =
     (PDM_DMesh_nodal_section_std_t *) PDM_Handles_get (mesh->sections_std, list_ind[i]);
   int n_face_elt     = PDM_n_face_elt_per_elmt(section->t_elt);
   int n_sum_vtx_face = PDM_n_sum_vtx_face_per_elmt(section->t_elt);

   _s_elt_face_vtx_idx += section->n_elt * n_face_elt;
   _s_elt_face_vtx     += section->n_elt * n_sum_vtx_face;
 }

 int n_sections_poly2d = PDM_Handles_n_get (mesh->sections_poly2d);
 list_ind = PDM_Handles_idx_get (mesh->sections_poly2d);

 for (int i = 0; i < n_sections_poly2d; i++) {
   PDM_DMesh_nodal_section_poly2d_t *section =
     (PDM_DMesh_nodal_section_poly2d_t *) PDM_Handles_get (mesh->sections_poly2d, list_ind[i]);
   _s_elt_face_vtx_idx +=     section->_connec_idx[section->n_elt];
   _s_elt_face_vtx     += 2 * section->_connec_idx[section->n_elt];
 }

 int n_sections_poly3d = PDM_Handles_n_get (mesh->sections_poly3d);
 list_ind = PDM_Handles_idx_get (mesh->sections_poly3d);

 for (int i = 0; i < n_sections_poly3d; i++) {
   PDM_DMesh_nodal_section_poly3d_t *section =
     (PDM_DMesh_nodal_section_poly3d_t *) PDM_Handles_get (mesh->sections_poly3d, list_ind[i]);
   _s_elt_face_vtx_idx += section->n_face;
   _s_elt_face_vtx     += section->_facvtx[section->_facvtx_idx[section->n_face]];
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
  PDM_DMesh_nodal_t *mesh,
  int               *s_elt_edge_vtx_idx,
  int               *s_elt_edge_vtx,
  int               *s_elt_edge_cell
)
{

 int _s_elt_edge_vtx_idx = 0;
 int _s_elt_edge_vtx     = 0;

 int n_sections_std  = PDM_Handles_n_get (mesh->sections_std);
 const int *list_ind = PDM_Handles_idx_get (mesh->sections_std);

 for (int i = 0; i < n_sections_std; i++) {
   PDM_DMesh_nodal_section_std_t *section =
     (PDM_DMesh_nodal_section_std_t *) PDM_Handles_get (mesh->sections_std, list_ind[i]);
   int n_edge_elt     = PDM_n_nedge_elt_per_elmt(section->t_elt);
   int n_sum_vtx_edge = PDM_n_sum_vtx_edge_per_elmt(section->t_elt);

   _s_elt_edge_vtx_idx += section->n_elt * n_edge_elt;
   _s_elt_edge_vtx     += section->n_elt * n_sum_vtx_edge;
 }

 int n_sections_poly2d = PDM_Handles_n_get (mesh->sections_poly2d);
 // list_ind = PDM_Handles_idx_get (mesh->sections_poly2d);
 assert(n_sections_poly2d == 0); // Not implemented

 // for (int i = 0; i < n_sections_poly2d; i++) {
 //   PDM_DMesh_nodal_section_poly2d_t *section =
 //     (PDM_DMesh_nodal_section_poly2d_t *) PDM_Handles_get (mesh->sections_poly2d, list_ind[i]);
 //   _s_elt_edge_vtx_idx +=     section->_connec_idx[section->n_elt];
 //   _s_elt_edge_vtx     += 2 * section->_connec_idx[section->n_elt];
 // }

 int n_sections_poly3d = PDM_Handles_n_get (mesh->sections_poly3d);
 // list_ind = PDM_Handles_idx_get (mesh->sections_poly3d);
 assert(n_sections_poly3d == 0); // Not implemented

 // for (int i = 0; i < n_sections_poly3d; i++) {
 //   PDM_DMesh_nodal_section_poly3d_t *section =
 //     (PDM_DMesh_nodal_section_poly3d_t *) PDM_Handles_get (mesh->sections_poly3d, list_ind[i]);
 //   _s_elt_edge_vtx_idx += section->n_face;
 //   _s_elt_edge_vtx     += section->_facvtx[section->_facvtx_idx[section->n_face]];
 // }

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
     n_face_elt = 1;  // TODO Eric have initialy 3
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
      // _current_elmt_edge_cell   [ielt * n_edge_elt + i_edge    ] = *n_elt_current + beg_gnum_elt_current + ielt + 1;
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
      // _current_elmt_edge_cell   [ielt * n_edge_elt + i_edge    ] = *n_elt_current + beg_gnum_elt_current + ielt + 1;
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

    _current_elmt_face_vtx_idx[ielt * n_face_elt + 1]  = elmt_face_vtx_idx[ielt * n_face_elt    ] + 4;
    _current_elmt_face_vtx_idx[ielt * n_face_elt + 2]  = elmt_face_vtx_idx[ielt * n_face_elt + 1] + 3;
    _current_elmt_face_vtx_idx[ielt * n_face_elt + 3]  = elmt_face_vtx_idx[ielt * n_face_elt + 2] + 3;
    _current_elmt_face_vtx_idx[ielt * n_face_elt + 4]  = elmt_face_vtx_idx[ielt * n_face_elt + 3] + 3;
    _current_elmt_face_vtx_idx[ielt * n_face_elt + 5]  = elmt_face_vtx_idx[ielt * n_face_elt + 4] + 3;

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
      // _current_elmt_face_cell   [ielt * n_face_elt + i_face    ] = *n_elt_current + ielt + 1;
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
  PDM_DMesh_nodal_t *mesh,
  int               *elmt_face_vtx_idx,
  PDM_g_num_t       *elmt_face_vtx,
  PDM_g_num_t       *elmt_face_cell,
  int               *elmt_cell_face_idx,
  PDM_g_num_t       *elmt_cell_face
)
{
  PDM_UNUSED(mesh);
  PDM_UNUSED(elmt_face_vtx_idx);
  PDM_UNUSED(elmt_face_vtx);
  PDM_UNUSED(elmt_face_cell);
  PDM_UNUSED(elmt_cell_face);

  int n_sections_std  = PDM_Handles_n_get  (mesh->sections_std);
  const int *list_ind = PDM_Handles_idx_get(mesh->sections_std);

  int n_elt_current  = 0;
  int n_face_current = 0;

  for (int i_section = 0; i_section < n_sections_std; i_section++) {
    PDM_DMesh_nodal_section_std_t *section = (PDM_DMesh_nodal_section_std_t *) PDM_Handles_get (mesh->sections_std, list_ind[i_section]);

    PDM_g_num_t beg_elmt_gnum = section->distrib[mesh->i_proc] + mesh->section_distribution[i_section];
    PDM_g_num_t beg_face_gnum = 0; // Useless in this context
    switch (section->t_elt) {
     case PDM_MESH_NODAL_POINT:
       abort();
       break;
     case PDM_MESH_NODAL_BAR2:
       abort();
       break;
     case PDM_MESH_NODAL_TRIA3:
       abort();
       break;
     case PDM_MESH_NODAL_QUAD4:
       abort();
       break;
     case PDM_MESH_NODAL_TETRA4:
       PDM_tetra_decomposes_faces(section->n_elt,
                                  &n_elt_current,
                                  &n_face_current,
                                  beg_elmt_gnum,
                                  beg_face_gnum,
                                  section->_connec,
                                  elmt_face_vtx_idx,
                                  elmt_face_vtx,
                                  elmt_face_cell,
                                  elmt_cell_face_idx,
                                  elmt_cell_face);
       break;
     case PDM_MESH_NODAL_PYRAMID5:
       PDM_pyra_decomposes_faces(section->n_elt,
                                 &n_elt_current,
                                 &n_face_current,
                                 beg_elmt_gnum,
                                 beg_face_gnum,
                                 section->_connec,
                                 elmt_face_vtx_idx,
                                 elmt_face_vtx,
                                 elmt_face_cell,
                                 elmt_cell_face_idx,
                                 elmt_cell_face);
       break;
     case PDM_MESH_NODAL_PRISM6:
       PDM_prism_decomposes_faces(section->n_elt,
                                  &n_elt_current,
                                  &n_face_current,
                                  beg_elmt_gnum,
                                  beg_face_gnum,
                                  section->_connec,
                                  elmt_face_vtx_idx,
                                  elmt_face_vtx,
                                  elmt_face_cell,
                                  elmt_cell_face_idx,
                                  elmt_cell_face);
       break;
     case PDM_MESH_NODAL_HEXA8:
       PDM_hexa_decomposes_faces(section->n_elt,
                                 &n_elt_current,
                                 &n_face_current,
                                 beg_elmt_gnum,
                                 beg_face_gnum,
                                 section->_connec,
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

  /* Not implemented */
  if (mesh->sections_poly2d != NULL) {
    // int n_sections_poly2d = PDM_Handles_n_get  (mesh->sections_poly2d);
    // const int *list_ind   = PDM_Handles_idx_get(mesh->sections_poly2d);
    abort();
  }

  if (mesh->sections_poly3d != NULL) {
    // int n_sections_poly3d = PDM_Handles_n_get  (mesh->sections_poly3d);
    // const int *list_ind   = PDM_Handles_idx_get(mesh->sections_poly3d);
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
  PDM_DMesh_nodal_t *mesh,
  int               *elmt_edge_vtx_idx,
  PDM_g_num_t       *elmt_edge_vtx,
  PDM_g_num_t       *elmt_edge_cell,
  int               *elmt_cell_edge_idx,
  PDM_g_num_t       *elmt_cell_edge
)
{
  PDM_UNUSED(mesh);
  PDM_UNUSED(elmt_edge_vtx_idx);
  PDM_UNUSED(elmt_edge_vtx);
  PDM_UNUSED(elmt_edge_cell);
  PDM_UNUSED(elmt_cell_edge_idx);
  PDM_UNUSED(elmt_cell_edge);

  int n_sections_std  = PDM_Handles_n_get  (mesh->sections_std);
  const int *list_ind = PDM_Handles_idx_get(mesh->sections_std);

  int n_elt_current  = 0;
  int n_dedge_current = 0;

  for (int i_section = 0; i_section < n_sections_std; i_section++) {
    PDM_DMesh_nodal_section_std_t *section = (PDM_DMesh_nodal_section_std_t *) PDM_Handles_get (mesh->sections_std, list_ind[i_section]);

    PDM_g_num_t beg_elmt_gnum = section->distrib[mesh->i_proc] + mesh->section_distribution[i_section];
    PDM_g_num_t beg_edge_gnum = 0; // Useless in this context

    switch (section->t_elt) {
     case PDM_MESH_NODAL_POINT:
       abort();
       break;
     case PDM_MESH_NODAL_BAR2:
       abort();
       break;
     case PDM_MESH_NODAL_TRIA3:
       PDM_tri_decomposes_edges(section->n_elt,
                                 &n_elt_current,
                                 &n_dedge_current,
                                 beg_elmt_gnum,
                                 beg_edge_gnum,
                                 section->_connec,
                                 elmt_edge_vtx_idx,
                                 elmt_edge_vtx,
                                 elmt_edge_cell,
                                 elmt_cell_edge_idx,
                                 elmt_cell_edge);
       break;
     case PDM_MESH_NODAL_QUAD4:
       PDM_quad_decomposes_edges(section->n_elt,
                                 &n_elt_current,
                                 &n_dedge_current,
                                 beg_elmt_gnum,
                                 beg_edge_gnum,
                                 section->_connec,
                                 elmt_edge_vtx_idx,
                                 elmt_edge_vtx,
                                 elmt_edge_cell,
                                 elmt_cell_edge_idx,
                                 elmt_cell_edge);
       break;
     case PDM_MESH_NODAL_TETRA4:
       // PDM_tetra_decomposes_edges(section->n_elt,
       //                            &n_elt_current,
       //                            &n_dedge_current,
       //                            section->_connec,
       //                            elmt_edge_vtx_idx,
       //                            elmt_edge_vtx,
       //                            elmt_edge_cell,
       //                            elmt_cell_edge);
       abort();
       break;
     case PDM_MESH_NODAL_PYRAMID5:
       // PDM_pyra_decomposes_edges(section->n_elt,
       //                           &n_elt_current,
       //                           &n_dedge_current,
       //                           section->_connec,
       //                           elmt_edge_vtx_idx,
       //                           elmt_edge_vtx,
       //                           elmt_edge_cell,
       //                           elmt_cell_edge);
       abort();
       break;
     case PDM_MESH_NODAL_PRISM6:
       // PDM_prism_decomposes_edges(section->n_elt,
       //                            &n_elt_current,
       //                            &n_dedge_current,
       //                            section->_connec,
       //                            elmt_edge_vtx_idx,
       //                            elmt_edge_vtx,
       //                            elmt_edge_cell,
       //                            elmt_cell_edge);
       abort();
       break;
     case PDM_MESH_NODAL_HEXA8:
       // PDM_hexa_decomposes_edges(section->n_elt,
       //                           &n_elt_current,
       //                           &n_dedge_current,
       //                           section->_connec,
       //                           elmt_edge_vtx_idx,
       //                           elmt_edge_vtx,
       //                           elmt_edge_cell,
       //                           elmt_cell_edge);
       abort();
       break;
     default:
       PDM_error(__FILE__, __LINE__, 0, "Error PDM_sections_decompose_edges : Element type is not taking int account\n");
    }
  }

  /* Not implemented */
  if (mesh->sections_poly2d != NULL) {
    // int n_sections_poly2d = PDM_Handles_n_get  (mesh->sections_poly2d);
    // const int *list_ind   = PDM_Handles_idx_get(mesh->sections_poly2d);
    abort();
  }

  if (mesh->sections_poly3d != NULL) {
    // int n_sections_poly3d = PDM_Handles_n_get  (mesh->sections_poly3d);
    // const int *list_ind   = PDM_Handles_idx_get(mesh->sections_poly3d);
    abort();
  }

}



// void
// PDM_hexa_section_decompose_elemt_to_face
// (
//       PDM_g_num_t  n_elmt,
// const PDM_g_num_t *elmt_vtx,
//       int         *elmt_face_vtx_idx,
//       PDM_g_num_t *elmt_face_vtx,
//       PDM_g_num_t *elmt_face_cell,
//       PDM_g_num_t *elmt_cell_face
// )
// {
//   PDM_UNUSED(n_elmt);
//   PDM_UNUSED(elmt_vtx);
//   PDM_UNUSED(elmt_face_vtx_idx);
//   PDM_UNUSED(elmt_face_vtx);
//   PDM_UNUSED(elmt_face_cell);
//   PDM_UNUSED(elmt_cell_face);

//   const int n_face_elt     = 6;
//   const int n_sum_vtx_face = 24;
//   const int n_sum_vtx_elt  = 8;



// }
