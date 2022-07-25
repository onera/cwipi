#ifndef __PDM_LAGRANGE_TO_BEZIER_H__
#define __PDM_LAGRANGE_TO_BEZIER_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mesh_nodal.h"
#include "pdm_mpi.h"
#include "pdm_io.h"

/*=============================================================================
 * Macro definition
 *============================================================================*/


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


/* Get Bezier coordinates from Lagrange coordinates for a bar (copied from pdm_t_dcube_nodal_gen.c) */
void
PDM_lagrange_to_bezier_bar
(
 const int  order,
 double    *lag,
 double    *bez
);


/* Get Bezier coordinates from Lagrange coordinates for a triangle (copied from pdm_t_dcube_nodal_gen.c) */
void
PDM_lagrange_to_bezier_tria
(
 const int  order,
 double    *lag,
 double    *bez
);


/* Get bezier bounding box (copied from pdm_t_dcube_nodal_gen.c) */
void
PDM_bezier_bounding_boxes
(
 const PDM_Mesh_nodal_elt_t  t_elt,
 const int                   order,
 const int                   n_nodes,
 int                         n_elt,
 double                     *lagrange_coord,
 double                    **extents
);

#endif /* __PDM_LAGRANGE_TO_BEZIER_H__ */
