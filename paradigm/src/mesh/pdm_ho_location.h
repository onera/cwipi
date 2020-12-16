#ifndef __PDM_HO_LOCATION_H__
#define __PDM_HO_LOCATION_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mesh_nodal.h"

/*=============================================================================
 * Macro definitions
 *============================================================================*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Type
 *============================================================================*/

/*----------------------------------------------------------------------------
 *
 * Callback to define location in a high order element
 *
 * parameters:
 *   order             <-- element order
 *   n_nodes           <-- number of nodes of the element
 *   nodes_coords      <-- nodes coordinates
 *   point_coords      <-- point to locate coordinates
 *   projected_coords  --> projected point coordinates (if point is outside)
 *   projected_uvw     --> parametric coordinates of the projected point
 *
 * return:
 *   distance to the cell (distance <= 0 if point is inside)
 *
 *----------------------------------------------------------------------------*/

typedef double (*PDM_ho_location_fct_t)
(const int     entities_dim,
 const int     order,
 const int     n_nodes,
 const double *nodes_coords,
 const double *point_coords,
 double       *projected_coords,
 double       *uvw);

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 *
 * Point location in a high order cell
 *
 * parameters:
 *   type             <-- element type
 *   order            <-- element order
 *   n_nodes          <-- number of nodes
 *   nodes_coords     <-- nodes coordinates (size = 3 * n_nodes)
 *   point_coords     <-- point to locate coordinates (size = 3)
 *   projected_coords --> projected point coordinates (size = 3)
 *   uvw              --> parametric coordinates of the projected point on the element
 *
 * return:
 *   distance to the cell
 *
 *----------------------------------------------------------------------------*/

double
PDM_ho_location
(
 const PDM_Mesh_nodal_elt_t  type,
 const int                   order,
 const int                   n_nodes,
 const double               *nodes_coords,
 const double               *point_coords,
 double                     *projected_coords,
 double                     *uvw
 );


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_HO_LOCATION_H__ */
