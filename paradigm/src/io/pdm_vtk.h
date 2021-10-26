#ifndef __PDM_VTK_H__
#define __PDM_VTK_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"
#include "pdm_io.h"
#include "pdm_error.h"
#include "pdm_mesh_nodal.h"

/*=============================================================================
 * Macro definition
 *============================================================================*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Types definition
 *============================================================================*/

/*=============================================================================
 * Global variables
 *============================================================================*/

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function interfaces
 *============================================================================*/

void
PDM_vtk_write_boxes
(
 const char        *filename,
 const int          n_box,
 const double      *box_extents,
 const PDM_g_num_t *box_g_num
);

void
PDM_vtk_write_circles
(
 const char        *filename,
 const int          n_circles,
 const double      *center,
 const double      *radius,
 const PDM_g_num_t *g_num,
 const int         *color,
 const int          resolution
);


void
PDM_vtk_write_polydata
(
 const char        *filename,
 const int          n_vtx,
 const double       vtx_coord[],
 const PDM_g_num_t  vtx_g_num[],
 const int          n_face,
 const int          face_vtx_idx[],
 const int          face_vtx[],
 const PDM_g_num_t  face_g_num[],
 const int          face_color[]
 );


void
PDM_vtk_write_point_cloud
(
 const char        *filename,
 const int          n_vtx,
 const double       vtx_coord[],
 const PDM_g_num_t  vtx_g_num[],
 const int          color[]
 );


void
PDM_vtk_write_lines
(
 const char        *filename,
 const int          n_line,
 const double      *coord,
 const PDM_g_num_t *g_num,
 const int         *color
 );



void
PDM_vtk_write_std_elements
(
 const char                 *filename,
 const int                   n_vtx,
 const double                vtx_coord[],
 const PDM_g_num_t           vtx_g_num[],
 const PDM_Mesh_nodal_elt_t  elt_type,
 const int                   n_elt,
 const int                   elt_vtx[],
 const PDM_g_num_t           elt_g_num[],
 const int                   n_elt_ifield,
 const char                 *elt_ifield_name[],
 const int                  *elt_ifield[]
 );


void
PDM_vtk_write_std_elements_double
(
 const char                 *filename,
 const int                   n_vtx,
 const double                vtx_coord[],
 const PDM_g_num_t           vtx_g_num[],
 const PDM_Mesh_nodal_elt_t  elt_type,
 const int                   n_elt,
 const int                   elt_vtx[],
 const PDM_g_num_t           elt_g_num[],
 const int                   n_elt_field,
 const char                 *elt_field_name[],
 const double               *elt_field[]
 );

void
PDM_vtk_write_std_elements_ho
(
 const char                 *filename,
 const int                   order,
 const int                   n_vtx,
 const double                vtx_coord[],
 const PDM_g_num_t           vtx_g_num[],
 const PDM_Mesh_nodal_elt_t  elt_type,
 const int                   n_elt,
 const int                   elt_vtx[],
 const PDM_g_num_t           elt_g_num[],
 const int                   n_elt_field,
 const char                 *elt_field_name[],
 const double               *elt_field[]
 );

void
PDM_vtk_write_ellipses
(
 const char        *filename,
 const int          n_ellipse,
 const double      *center,
 const double      *axes,
 const double      *radii,
 const PDM_g_num_t *g_num,
 const int         *color,
 const int          resolution
 );


int *
PDM_vtk_lagrange_to_ijk
(
 const PDM_Mesh_nodal_elt_t elt_type,
 const int                  order
 );

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_VTK_H__ */
