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

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <math.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

#include <bftc_error.h>
#include <bftc_mem.h>
#include <bftc_printf.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "fvmc_config_defs.h"
#include "fvmc_defs.h"
#include "fvmc_nodal.h"
#include "fvmc_ho.h"
#include "fvmc_nodal_priv.h"
#include "fvmc_triangulate.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "fvmc_ho.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Local macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*============================================================================
 * Public function definitions
 *============================================================================*/


/*----------------------------------------------------------------------------
 * 
 * Point location in a high order cell 3d
 * 
 * parameters:
 *   type             <-- element type
 *   order            <-- element order
 *   ho_vertex_num    <-- high order vertex num (internal ordering)
 *   vertex_coords    <-- vertex coordinates
 *   point_coords     <-- point to locate coordinates
 *   projected_coords --> projected point coordinates if outside (or NULL)
 * 
 * return: 
 *   distance to the cell (distance <= 0 if point is inside)
 *
 *----------------------------------------------------------------------------*/

double 
fvmc_ho_location_in_cell_3d (const fvmc_element_t type,
                             const int order,
                             const int *ho_vertex_num,
                             const double *vertex_coords,
                             const double *point_coords,
                             double *projected_coords)
{
  double dist = 0.;
  bftc_error(__FILE__, __LINE__, 0,
             _("fvmc_ho_location_in_cell_3d : Not implement yet\n"));
  return dist;
}

/*----------------------------------------------------------------------------
 * 
 * Point location on a high order cell 2d
 * 
 * parameters:
 *   type             <-- element type
 *   order            <-- element order
 *   ho_vertex_num    <-- high order vertex num (internal ordering)
 *   vertex_coords    <-- vertex coordinates
 *   point_coords     <-- point to locate coordinates
 *   projected_coords --> projected point coordinates (or NULL)
 * 
 * return: 
 *   distance to the cell
 *
 *----------------------------------------------------------------------------*/

double 
fvmc_ho_location_on_cell_2d (const fvmc_element_t type,
                             const int order,
                             const int *ho_vertex_num,
                             const double *vertex_coords,
                             const double *point_coords,
                             double *projected_coords)
{
  double dist = 0.;
  bftc_error(__FILE__, __LINE__, 0,
             _("fvmc_ho_location_on_cell_2d : Not implement yet\n"));
  return dist;
}

/*----------------------------------------------------------------------------
 * 
 * Point location on a high order cell 1d
 * 
 * parameters:
 *   type             <-- element type
 *   order            <-- element order
 *   ho_vertex_num    <-- high order vertex num (internal ordering)
 *   vertex_coords    <-- vertex coordinates
 *   point_coords     <-- point to locate coordinates
 *   projected_coords --> projected point coordinates (or NULL)
 * 
 * return: 
 *   distance to the cell
 *
 *----------------------------------------------------------------------------*/

double 
fvmc_ho_location_on_cell_1d (const fvmc_element_t type,
                             const int order,
                             const int *ho_vertex_num,
                             const double *vertex_coords,
                             const double *point_coords,
                             double *projected_coords)
{
  double dist = 0.;
  bftc_error(__FILE__, __LINE__, 0,
             _("fvmc_ho_location_on_cell_1d : Not implement yet\n"));
  return dist;
}

/*----------------------------------------------------------------------------
 * 
 * Compute shape in a high order cell 3d
 * 
 * parameters:
 *   type             <-- element type
 *   order            <-- element order
 *   ho_vertex_num    <-- high order vertex num (internal ordering)
 *   vertex_coords    <-- vertex coordinates
 *   point_coords     <-- point inside cell  
 *   shape            --> barycenter's coordinates
 * 
 *----------------------------------------------------------------------------*/

void 
fvmc_ho_shape_in_cell_3d (const fvmc_element_t type,
                          const int order,
                          const int *ho_vertex_num,
                          const double *vertex_coords,
                          const double *point_coords,
                          double *shape)
{
  bftc_error(__FILE__, __LINE__, 0,
             _("fvmc_ho_shape_in_cell_3d : Not implement yet\n"));
}

/*----------------------------------------------------------------------------
 * 
 * Compute shape on a high order cell 2d
 * 
 * parameters:
 *   type             <-- element type
 *   order            <-- element order
 *   ho_vertex_num    <-- high order vertex num (internal ordering)
 *   vertex_coords    <-- vertex coordinates
 *   point_coords     <-- point on cell  
 *   shape            --> barycenter's coordinates
 * 
 *----------------------------------------------------------------------------*/

void 
fvmc_ho_shape_on_cell_2d (const fvmc_element_t type,
                          const int order,
                          const int *ho_vertex_num,
                          const double *vertex_coords,
                          const double *point_coords,
                          double *shape)
{
  bftc_error(__FILE__, __LINE__, 0,
             _("fvmc_ho_shape_on_cell_2d : Not implement yet\n"));
}

/*----------------------------------------------------------------------------
 * 
 * Compute shape on a high order cell 1d
 * 
 * parameters:
 *   type             <-- element type
 *   order            <-- element order
 *   ho_vertex_num    <-- high order vertex num (internal ordering)
 *   vertex_coords    <-- vertex coordinates
 *   point_coords     <-- point on cell  
 *   shape            --> barycenter's coordinates
 * 
 *----------------------------------------------------------------------------*/

void 
fvmc_ho_shape_on_cell_1d (const fvmc_element_t type,
                          const int order,
                          const int *ho_vertex_num,
                          const double *vertex_coords,
                          const double *point_coords,
                          double *shape)
{
  bftc_error(__FILE__, __LINE__, 0,
             _("fvmc_ho_shape_on_cell_1d : Not implement yet\n"));
}


/*----------------------------------------------------------------------------
 * 
 * Interpolate field to the target point_coords in a 3D cell
 * 
 * parameters:
 *   type              <-- element type
 *   order             <-- element order
 *   ho_vertex_num     <-- high order vertex num (internal ordering)
 *   local_to_user     <-- local to user ordering (for type)
 *   vertex_coords     <-- vertex coordinates
 *   point_coords      <-- point inside cell (or on boundary) 
 *   shape             <-- barycenter's coordinates
 *   stride_field      <-- field stride
 *   source_field      <-- source field (user ordering) 
 *   target_field      --> target field (defined to point_coords)
 *
 *----------------------------------------------------------------------------*/

void 
fvmc_ho_interp_in_cell_3d (const fvmc_element_t type,
                           const int order,
                           const int *ho_vertex_num,
                           const int *local_to_user,
                           const double *vertex_coords,
                           const double *point_coords,
                           const double *shape,
                           const int stride_field,
                           const double *src_field,
                           double *target_field)
                           
{
  bftc_error(__FILE__, __LINE__, 0,
             _("fvmc_ho_interp_in_cell_3d : Not implement yet\n"));
}


/*----------------------------------------------------------------------------
 * 
 * Interpolate field to the target point_coords on a 2D cell
 * 
 * parameters:
 *   type              <-- element type
 *   order             <-- element order
 *   ho_vertex_num     <-- high order vertex num (internal ordering)
 *   local_to_user     <-- local to user ordering (for type)
 *   vertex_coords     <-- vertex coordinates
 *   point_coords      <-- point inside cell (or on boundary) 
 *   shape             <-- barycenter's coordinates
 *   stride_field      <-- field stride
 *   source_field      <-- source field (user ordering) 
 *   target_field      --> target field (defined to point_coords)
 *
 *----------------------------------------------------------------------------*/

void 
fvmc_ho_interp_on_cell_2d (const fvmc_element_t type,
                           const int order,
                           const int *ho_vertex_num,
                           const int *local_to_user,
                           const double *vertex_coords,
                           const double *point_coords,
                           const double *shape,
                           const int stride_field,
                           const double *src_field,
                           double *target_field)
{
  bftc_error(__FILE__, __LINE__, 0,
             _("fvmc_ho_interp_on_cell_2d : Not implement yet\n"));
}


/*----------------------------------------------------------------------------
 * 
 * Interpolate field to the target point_coords on a 2D cell
 * 
 * parameters:
 *   type              <-- element type
 *   order             <-- element order
 *   ho_vertex_num     <-- high order vertex num (internal ordering)
 *   local_to_user     <-- local to user ordering (for type)
 *   vertex_coords     <-- vertex coordinates
 *   point_coords      <-- point inside cell (or on boundary) 
 *   shape             <-- barycenter's coordinates
 *   stride_field      <-- field stride
 *   source_field      <-- source field (user ordering) 
 *   target_field      --> target field (defined to point_coords)
 *
 *----------------------------------------------------------------------------*/

void 
fvmc_ho_interp_on_cell_1d (const fvmc_element_t type,
                           const int order,
                           const int *ho_vertex_num,
                           const int *local_to_user,
                           const double *vertex_coords,
                           const double *point_coords,
                           const double *shape,
                           const int stride_field,
                           const double *src_field,
                           double *target_field)
{
  bftc_error(__FILE__, __LINE__, 0,
             _("fvmc_ho_shape_on_cell_1d : Not implement yet\n"));
}
