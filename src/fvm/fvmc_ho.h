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
 * Function pointer to locate in a high order cell 3d
 * 
 * parameters:
 *   type             <-- element type
 *   order            <-- element order
 *   n_nodes          <-- number of nodes
 *   nodes_coords     <-- nodes coordinates
 *   point_coords     <-- point to locate coordinates
 *   projected_coords --> projected point coordinates (if point is outside) 
 *   uvw              --> parametric coordinates of the projected points in the element 
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


typedef void (*fvmc_ho_basis_fct_t)
(const int order,
 const int n_nodes,
 const double *uvw,
 double *weights);

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
 * high order basis
 * 
 * parameters:
 *   type            <-- element type
 *   order           <-- order
 *   uvw             <-- uvw
 *   weights         --> weights (size = n_nodes)
 *
 *----------------------------------------------------------------------------*/

void
fvmc_ho_basis_pn
(
const fvmc_element_t type,
const int order,
const double *uvw,
      double *weights 
);


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

void  
fvmc_ho_user_elementary_functions_set (fvmc_ho_location_fct_t location_tetra,
                                       fvmc_ho_location_fct_t location_prism,
                                       fvmc_ho_location_fct_t location_pyramid,
                                       fvmc_ho_location_fct_t location_hexa,
                                       fvmc_ho_location_fct_t location_tria,
                                       fvmc_ho_location_fct_t location_quad,
                                       fvmc_ho_location_fct_t location_edge,
                                       fvmc_ho_interp_fct_t interp_tetra,
                                       fvmc_ho_interp_fct_t interp_prism,
                                       fvmc_ho_interp_fct_t interp_pyramid,
                                       fvmc_ho_interp_fct_t interp_hexa,
                                       fvmc_ho_interp_fct_t interp_tria,
                                       fvmc_ho_interp_fct_t interp_quad,
                                       fvmc_ho_interp_fct_t interp_edge);


/*----------------------------------------------------------------------------
 * 
 * Unset elementary functions
 * 
 *----------------------------------------------------------------------------*/

void
fvmc_ho_user_elt_unset (fvmc_element_t elt_type);

void
fvmc_ho_user_elementary_functions_unset (void);

/*----------------------------------------------------------------------------
 * 
 * Point location in a high order cell 3d
 * 
 * parameters:
 *   type             <-- element type
 *   order            <-- element order
 *   n_nodes          <-- number of nodes
 *   nodes_coords     <-- nodes coordinates (size = 3 * n_nodes)
 *   point_coords     <-- point to locate coordinates (size = 3)
 *   projected_coords --> projected point coordinates if outside (size = 3)
 *   uvw              --> parametric coordinates (point if inside the element
 *                                                projected point if outside)
 *
 * return: 
 *   distance to the cell (distance <= 0 if point is inside)
 *
 *----------------------------------------------------------------------------*/

double 
fvmc_ho_location_in_cell_3d (const fvmc_element_t type,
                             const int order,
                             const int n_nodes,
                             const double *nodes_coords,
                             const double *point_coords,
                             double *projected_coords,
                             double* uvw);


/*----------------------------------------------------------------------------
 * 
 * Point location on a high order cell 2d
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
fvmc_ho_location_on_cell_2d (const fvmc_element_t type,
                             const int order,
                             const int n_nodes,
                             const double *nodes_coords,
                             const double *point_coords,
                             double *projected_coords,
                             double *uvw);


/*----------------------------------------------------------------------------
 * 
 * Point location on a high order cell 1d
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
fvmc_ho_location_on_cell_1d (const fvmc_element_t type,
                             const int order,
                             const int n_nodes,
                             const double *nodes_coords,
                             const double *point_coords,
                             double *projected_coords,
                             double *uvw);

/*----------------------------------------------------------------------------
 * 
 * Interpolate field to the target point_coords in a 3D cell
 * 
 * parameters:
 *   type              <-- element type
 *   order             <-- element order
 *   n_node            <-- number of nodes
 *   ho_vertex_num     <-- high order vertex num (internal ordering)
 *   local_to_user     <-- local to user ordering (for type)
 *   vertex_coords     <-- vertex coordinates
 *   point_coords      <-- point inside cell (or on boundary) 
 *   distance          <-- distance to the element
 *   point_proj_coords  <-- projected point coordinates
 *   weight             <-- weights
 *   stride_field      <-- field stride
 *   source_field      <-- source field (user ordering) 
 *   target_field      --> target field (defined to point_coords)
 *
 *----------------------------------------------------------------------------*/

void 
fvmc_ho_interp_in_cell_3d (const fvmc_element_t type,
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

/*----------------------------------------------------------------------------
 * 
 * Interpolate field to the target point_coords on a 2D cell
 * 
 * parameters:
 *   type              <-- element type
 *   order             <-- element order
 *   n_node            <-- number of nodes
 *   ho_vertex_num     <-- high order vertex num (internal ordering)
 *   local_to_user     <-- local to user ordering (for type)
 *   vertex_coords     <-- vertex coordinates
 *   point_coords      <-- point inside cell (or on boundary) 
 *   distance          <-- distance to the element
 *   point_proj_coords  <-- projected point coordinates
 *   weight             <-- weights
 *   stride_field      <-- field stride
 *   source_field      <-- source field (user ordering) 
 *   target_field      --> target field (defined to point_coords)
 *
 *----------------------------------------------------------------------------*/

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


/*----------------------------------------------------------------------------
 * 
 * Interpolate field to the target point_coords on a 2D cell
 * 
 * parameters:
 *   type              <-- element type
 *   order             <-- element order
 *   n_node            <-- number of nodes
 *   ho_vertex_num     <-- high order vertex num (internal ordering)
 *   local_to_user     <-- local to user ordering (for type)
 *   vertex_coords     <-- vertex coordinates
 *   point_coords      <-- point inside cell (or on boundary)
 *   distance          <-- distance to the element
 *   point_proj_coords  <-- projected point coordinates
 *   weight             <-- barycenter's coordinates
 *   stride_field      <-- field stride
 *   source_field      <-- source field (user ordering) 
 *   target_field      --> target field (defined to point_coords)
 *
 *----------------------------------------------------------------------------*/

void 
fvmc_ho_interp_on_cell_1d (const fvmc_element_t type,
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


/*----------------------------------------------------------------------------*/



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

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __FVMC_HO_H__ */
