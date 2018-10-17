#ifndef __FVMC_HO_H__
#define __FVMC_HO_H__

/*============================================================================
 * Functions about high order meshes
 *============================================================================*/

/*
  This file is part of the "Finite Volume Mesh" library, intended to provide
  finite volume mesh and associated fields I/O and manipulation services.

  Copyright (C) 2018       ONERA

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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "fvmc_config_defs.h"
#include "fvmc_defs.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
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

typedef double (*fvmc_ho_location_fct_t)
(const int order,
 const int n_nodes,
 const double *nodes_coords,
 const double *point_coords,
 double *projected_coords,
 double *uvw);

/*----------------------------------------------------------------------------
 * 
 * Callback to define the basis functions of an high order 
 * element
 * 
 * parameters:
 *   order             <-- element order
 *   n_nodes           <-- number of nodes of the element
 *   n_pts             <-- number of points 
 *   uvw               <-- Parametric coordinates of points
 *   projected_uvw     --> Interpolation weights associated to uvw coordinates
 * 
 *----------------------------------------------------------------------------*/

typedef void (*fvmc_ho_basis_fct_t)
(const int order,
 const int n_nodes,
 const int n_pts,
 const double *uvw,
 double *weights);

/*----------------------------------------------------------------------------
 * 
 * Callback to define parametric coordinates of the element nodes
 * 
 * parameters:
 *   order             <-- element order
 *   n_nodes           <-- number of nodes of the element
 *   xsi_uvv           --> Parametric coordinates of a the element nodes
 *                         (size = dim of element * n_nodes)
 * 
 *----------------------------------------------------------------------------*/

typedef void (*fvmc_ho_xsi_fct_t)
(const int order,
 const int n_nodes,
 double *xsi_coords);

/*----------------------------------------------------------------------------
 * Function pointer to define an high order interpolation
 *
 * parameters:
 *   type              <-- element type
 *   order             <-- element order
 *   n_nodes           <-- number of nodes
 *   local_to_user     <-- local to user ordering (for type)
 *   nodes_coords      <-- nodes coordinates
 *   point_coords      <-- point coordinates 
 *   distance          <-- distance to the element
 *   point_proj_coords  <-- projected point coordinates
 *   weight             <-- weights
 *   stride_field      <-- field stride
 *   source_field      <-- source field (user ordering) 
 *   target_field      --> target field (defined to point_coords)
 *
 *----------------------------------------------------------------------------*/

typedef void (*fvmc_ho_interp_fct_t)
(const int order,
 const int n_nodes,
 const int *ho_vertex_num,
 const int *local_to_user,
 const double *nodes_coords,
 const double *point_coords,
 const float *distance,
 const double *point_proj_coords,
 const double *weight,
 const int stride_field,
 const double *src_field,
 double *target_field);

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * 
 * Set elementary functions
 * 
 * parameters:
 *   location_tetra    <-- Location in a tetrahedron
 *   location_prism    <-- Location in a prism
 *   location_pyramid  <-- Location in a pyramid
 *   location_hexa     <-- Location in a hexaedron
 *   location_tria     <-- Location on a triangle
 *   location_quad     <-- Location on a quandragle
 *   location_edge     <-- Location on a edge
 *   interp_tetra       <-- Interpolation in a tetrahedron
 *   interp_prism       <-- Interpolation in a prism
 *   interp_pyramid     <-- Interpolation in a pyramid
 *   interp_hexa        <-- Interpolation in a hexaedron
 *   interp_tria        <-- Interpolation on a triangle
 *   interp_quad        <-- Interpolation on a quandragle
 *   interp_edge        <-- Interpolation on a edge
 *
 *----------------------------------------------------------------------------*/

void
fvmc_ho_user_elt_set (fvmc_element_t elt_type,
                      fvmc_ho_basis_fct_t elt_basis,
                      fvmc_ho_xsi_fct_t xsi_coords,
                      fvmc_ho_location_fct_t location_in_elt);


/*----------------------------------------------------------------------------
 * 
 * Unset elementary functions
 * 
 *----------------------------------------------------------------------------*/

void
fvmc_ho_user_elt_unset (fvmc_element_t elt_type);

void
fvmc_ho_user_elts_unset (void);

/*----------------------------------------------------------------------------
 * 
 * high order basis
 * 
 * parameters:
 *   type            <-- element type
 *   order           <-- order
 *   n_nodes         <-- number of nodes
 *   n_pts           <-- number of points 
 *   uvw             <-- uvw (size = elt_dim * n_pts)
 *   weights         --> weights (size = n_nodes * n_pts)
 *
 *----------------------------------------------------------------------------*/

void
fvmc_ho_basis
(
const fvmc_element_t type,
const int order,
const int n_nodes,
const int n_pts,
const double *uvw,
      double *weights 
);


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
fvmc_ho_location
(
 const fvmc_element_t type,
 const int order,
 const int n_nodes,
 const double *nodes_coords,
 const double *point_coords,
 double *projected_coords,
 double *uvw
 );

/*----------------------------------------------------------------------------
 * 
 * Free static variables
 * 
 *----------------------------------------------------------------------------*/

void
fvmc_ho_free
(
 void
 );

/*----------------------------------------------------------------------------*/


void 
fvmc_ho_interp_on_cell_2d (const fvmc_element_t type,
                           const int order,
                           const int n_node,
                           const int *ho_vertex_num,
                           const int *local_to_user,
                           const double *vertex_coords,
                           const double *point_coords,
                           const float *distance,
                           const double *point_proj_coords,
                           const double *weight,
                           const int stride_field,
                           const double *src_field,
                           double *target_field);




#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __FVMC_HO_H__ */
