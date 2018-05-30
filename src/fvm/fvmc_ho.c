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

typedef struct _fvmc_ho_user_fcts_t {
  fvmc_ho_location_fct_t location_tetra;
  fvmc_ho_location_fct_t location_prism;
  fvmc_ho_location_fct_t location_pyramid;
  fvmc_ho_location_fct_t location_hexa;
  fvmc_ho_location_fct_t location_tria;
  fvmc_ho_location_fct_t location_quad;
  fvmc_ho_location_fct_t location_edge;
  
  fvmc_ho_shape_fct_t shape_tetra;
  fvmc_ho_shape_fct_t shape_prism;
  fvmc_ho_shape_fct_t shape_pyramid;
  fvmc_ho_shape_fct_t shape_hexa;
  fvmc_ho_shape_fct_t shape_tria;
  fvmc_ho_shape_fct_t shape_quad;
  fvmc_ho_shape_fct_t shape_edge;


  fvmc_ho_interp_fct_t interp_tetra;
  fvmc_ho_interp_fct_t interp_prism;
  fvmc_ho_interp_fct_t interp_pyramid;
  fvmc_ho_interp_fct_t interp_hexa;
  fvmc_ho_interp_fct_t interp_tria;
  fvmc_ho_interp_fct_t interp_quad;
  fvmc_ho_interp_fct_t interp_edge;

} fvmc_ho_user_fcts_t;

/*============================================================================
 * Static global variables
 *============================================================================*/

static fvmc_ho_user_fcts_t *_user_fcts = NULL;

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * 
 * Default point location in a high order cell 3d
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

static double 
_default_location_in_cell_3d (const fvmc_element_t type,
                             const int order,
                             const int *ho_vertex_num,
                             const double *vertex_coords,
                             const double *point_coords,
                             double *projected_coords)
{
  double dist = 0.;
  bftc_error(__FILE__, __LINE__, 0,
             _("_default_location_in_cell_3d : Not implement yet\n"));
  return dist;
}

/*----------------------------------------------------------------------------
 * 
 * Default point location on a high order cell 2d
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

static double 
_default_location_on_cell_2d (const fvmc_element_t type,
                             const int order,
                             const int *ho_vertex_num,
                             const double *vertex_coords,
                             const double *point_coords,
                             double *projected_coords)
{
  double dist = 0.;
  bftc_error(__FILE__, __LINE__, 0,
             _("_default_location_on_cell_2d : Not implement yet\n"));
  return dist;
}

/*----------------------------------------------------------------------------
 * 
 * Default point location on a high order cell 1d
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

static double 
_default_location_on_cell_1d (const fvmc_element_t type,
                             const int order,
                             const int *ho_vertex_num,
                             const double *vertex_coords,
                             const double *point_coords,
                             double *projected_coords)
{
  double dist = 0.;
  bftc_error(__FILE__, __LINE__, 0,
             _("_default_location_on_cell_1d : Not implement yet\n"));
  return dist;
}

/*----------------------------------------------------------------------------
 * 
 * Default shape computataion in a high order cell 3d
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
static _default_shape_in_cell_3d (const fvmc_element_t type,
                          const int order,
                          const int *ho_vertex_num,
                          const double *vertex_coords,
                          const double *point_coords,
                          double *shape)
{
  bftc_error(__FILE__, __LINE__, 0,
             _("_default_shape_in_cell_3d : Not implement yet\n"));
}


/*----------------------------------------------------------------------------
 * 
 * Default shape computataion on a high order cell 2d
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
static _default_shape_on_cell_2d (const fvmc_element_t type,
                          const int order,
                          const int *ho_vertex_num,
                          const double *vertex_coords,
                          const double *point_coords,
                          double *shape)
{
  bftc_error(__FILE__, __LINE__, 0,
             _("_default_shape_on_cell_2d : Not implement yet\n"));
}

/*----------------------------------------------------------------------------
 * 
 * Default shape computataion on a high order cell 1d
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
static _default_shape_on_cell_1d (const fvmc_element_t type,
                          const int order,
                          const int *ho_vertex_num,
                          const double *vertex_coords,
                          const double *point_coords,
                          double *shape)
{
  bftc_error(__FILE__, __LINE__, 0,
             _("_default_shape_on_cell_1d : Not implement yet\n"));
}


/*----------------------------------------------------------------------------
 * 
 * Default field interpolation field to the target point_coords in a 3D cell
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
static _default_interp_in_cell_3d (const fvmc_element_t type,
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
             _("_default_interp_in_cell_3d : Not implement yet\n"));
}


/*----------------------------------------------------------------------------
 * 
 * Default field interpolation field to the target point_coords on a 2D cell
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

static void 
_default_interp_on_cell_2d (const fvmc_element_t type,
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
             _("_default_interp_on_cell_2d : Not implement yet\n"));
}


/*----------------------------------------------------------------------------
 * 
 * Default field interpolation field to the target point_coords on a 1D cell
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

static void 
_default_interp_on_cell_1d (const fvmc_element_t type,
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
             _("_default_shape_on_cell_1d : Not implement yet\n"));
}

/*============================================================================
 * Public function definitions
 *============================================================================*/


/*----------------------------------------------------------------------------
 * 
 * Unset elementary functions
 * 
 *----------------------------------------------------------------------------*/

void
fvmc_ho_user_elementary_functions_unset (void)
{
  if (_user_fcts != NULL) {
    free (_user_fcts);
    _user_fcts = NULL;
  }
}


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
 *   shape_tetra       <-- Shape computation in a tetrahedron
 *   shape_prism       <-- Shape computation in a prism
 *   shape_pyramid     <-- Shape computation in a pyramid
 *   shape_hexa        <-- Shape computation in a hexaedron
 *   shape_tria        <-- Shape computation on a triangle
 *   shape_quad        <-- Shape computation on a quandragle
 *   shape_edge        <-- Shape computation on a edge
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
fvmc_ho_user_elementary_functions_set (fvmc_ho_location_fct_t location_tetra,
                                       fvmc_ho_location_fct_t location_prism,
                                       fvmc_ho_location_fct_t location_pyramid,
                                       fvmc_ho_location_fct_t location_hexa,
                                       fvmc_ho_location_fct_t location_tria,
                                       fvmc_ho_location_fct_t location_quad,
                                       fvmc_ho_location_fct_t location_edge,
                                       fvmc_ho_shape_fct_t shape_tetra,
                                       fvmc_ho_shape_fct_t shape_prism,
                                       fvmc_ho_shape_fct_t shape_pyramid,
                                       fvmc_ho_shape_fct_t shape_hexa,
                                       fvmc_ho_shape_fct_t shape_tria,
                                       fvmc_ho_shape_fct_t shape_quad,
                                       fvmc_ho_shape_fct_t shape_edge,
                                       fvmc_ho_interp_fct_t interp_tetra,
                                       fvmc_ho_interp_fct_t interp_prism,
                                       fvmc_ho_interp_fct_t interp_pyramid,
                                       fvmc_ho_interp_fct_t interp_hexa,
                                       fvmc_ho_interp_fct_t interp_tria,
                                       fvmc_ho_interp_fct_t interp_quad,
                                       fvmc_ho_interp_fct_t interp_edge)
{

  if (_user_fcts == NULL) {
    _user_fcts = (fvmc_ho_user_fcts_t *) malloc (sizeof(fvmc_ho_user_fcts_t));
  }

  _user_fcts->location_tetra   = location_tetra;
  _user_fcts->location_prism   = location_prism;
  _user_fcts->location_pyramid = location_pyramid;
  _user_fcts->location_hexa    = location_hexa;
  _user_fcts->location_tria    = location_tria;
  _user_fcts->location_quad    = location_quad;
  _user_fcts->location_edge    = location_edge;

  _user_fcts->shape_tetra   = shape_tetra;
  _user_fcts->shape_prism   = shape_prism;
  _user_fcts->shape_pyramid = shape_pyramid;
  _user_fcts->shape_hexa    = shape_hexa;
  _user_fcts->shape_tria    = shape_tria;
  _user_fcts->shape_quad    = shape_quad;
  _user_fcts->shape_edge    = shape_edge;

  _user_fcts->interp_tetra   = interp_tetra;
  _user_fcts->interp_prism   = interp_prism;
  _user_fcts->interp_pyramid = interp_pyramid;
  _user_fcts->interp_hexa    = interp_hexa;
  _user_fcts->interp_tria    = interp_tria;
  _user_fcts->interp_quad    = interp_quad;
  _user_fcts->interp_edge    = interp_edge;
}


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

  if (_user_fcts != NULL) {

    switch (type) {
      
    case FVMC_CELL_TETRA:
      return (_user_fcts->location_tetra) (order,
                                           ho_vertex_num,
                                           vertex_coords,
                                           point_coords,
                                           projected_coords);
      break;
      
    case FVMC_CELL_PRISM:
      return (_user_fcts->location_prism) (order,
                                           ho_vertex_num,
                                           vertex_coords,
                                           point_coords,
                                           projected_coords);
      break;
      
    case FVMC_CELL_PYRAM:
      return (_user_fcts->location_pyramid) (order,
                                             ho_vertex_num,
                                             vertex_coords,
                                             point_coords,
                                             projected_coords);
      break;

    case FVMC_CELL_HEXA:
      return (_user_fcts->location_hexa) (order,
                                          ho_vertex_num,
                                          vertex_coords,
                                          point_coords,
                                          projected_coords);
      break;

    default:

      bftc_error(__FILE__, __LINE__, 0,
                 _("fvmc_ho_location_in_cell_3d : Not a high order 3D element type\n"));
    } 
  }

  else {

    return _default_location_in_cell_3d (type,
                                         order,
                                         ho_vertex_num,
                                         vertex_coords,
                                         point_coords,
                                         projected_coords);
    
  }
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

  if (_user_fcts != NULL) {
    switch (type) {
      
    case FVMC_FACE_TRIA:
      return (_user_fcts->location_tria) (order,
                                          ho_vertex_num,
                                          vertex_coords,
                                          point_coords,
                                          projected_coords);
      break;
      
    case FVMC_FACE_QUAD:
      return (_user_fcts->location_quad) (order,
                                          ho_vertex_num,
                                          vertex_coords,
                                          point_coords,
                                          projected_coords);
      break;
      
    default:

      bftc_error(__FILE__, __LINE__, 0,
                 _("fvmc_ho_location_on_cell_2d : Not a high order 2D element type\n"));
    }
  }

  else {

    return _default_location_on_cell_2d (type,
                                         order,
                                         ho_vertex_num,
                                         vertex_coords,
                                         point_coords,
                                         projected_coords);
    
  }
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

  if (_user_fcts != NULL) {
    switch (type) {
      
    case FVMC_EDGE:
      return (_user_fcts->location_edge) (order,
                                          ho_vertex_num,
                                          vertex_coords,
                                          point_coords,
                                          projected_coords);
      break;
      
    default:

      bftc_error(__FILE__, __LINE__, 0,
                 _("fvmc_ho_location_on_cell_1d : Not a high order 1D element type\n"));
    }
  }

  else {

    return _default_location_on_cell_1d (type,
                                         order,
                                         ho_vertex_num,
                                         vertex_coords,
                                         point_coords,
                                         projected_coords);

  }
  
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
  if (_user_fcts != NULL) {

    switch (type) {
      
    case FVMC_CELL_TETRA:
      (_user_fcts->shape_tetra) (order,
                                 ho_vertex_num,
                                 vertex_coords,
                                 point_coords,
                                 shape);
      break;
      
    case FVMC_CELL_PRISM:
      (_user_fcts->shape_prism) (order,
                                 ho_vertex_num,
                                 vertex_coords,
                                 point_coords,
                                 shape);
      break;
      
    case FVMC_CELL_PYRAM:
      (_user_fcts->shape_pyramid) (order,
                                   ho_vertex_num,
                                   vertex_coords,
                                   point_coords,
                                   shape);
      break;

    case FVMC_CELL_HEXA:
      (_user_fcts->shape_hexa) (order,
                                ho_vertex_num,
                                vertex_coords,
                                point_coords,
                                shape);
      break;

    default:

      bftc_error(__FILE__, __LINE__, 0,
                 _("fvmc_ho_shape_in_cell_3d : Not a high order 3D element type\n"));
    } 

  }

  else {

    _default_shape_in_cell_3d (type,
                               order,
                               ho_vertex_num,
                               vertex_coords,
                               point_coords,
                               shape); 

  }
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
  if (_user_fcts != NULL) {
    switch (type) {
      
    case FVMC_FACE_TRIA:
      (_user_fcts->shape_tria) (order,
                                ho_vertex_num,
                                vertex_coords,
                                point_coords,
                                shape);
      break;
      
    case FVMC_FACE_QUAD:
      (_user_fcts->shape_quad) (order,
                                ho_vertex_num,
                                vertex_coords,
                                point_coords,
                                shape);
      break;
      
    default:

      bftc_error(__FILE__, __LINE__, 0,
                 _("fvmc_ho_shape_on_cell_2d : Not a high order 2D element type\n"));
    }
  }

  else {

    _default_shape_on_cell_2d (type,
                               order,
                               ho_vertex_num,
                               vertex_coords,
                               point_coords,
                               shape); 

  }
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
  if (_user_fcts != NULL) {
    switch (type) {
      
    case FVMC_EDGE:
      (_user_fcts->shape_edge) (order,
                                ho_vertex_num,
                                vertex_coords,
                                point_coords,
                                shape);
      break;
      
    default:

      bftc_error(__FILE__, __LINE__, 0,
                 _("fvmc_ho_shape_on_cell_1d : Not a high order 1D element type\n"));
    }

  }

  else {

    _default_shape_on_cell_1d (type,
                               order,
                               ho_vertex_num,
                               vertex_coords,
                               point_coords,
                               shape); 

  }
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
  if (_user_fcts != NULL) {

    switch (type) {
      
    case FVMC_CELL_TETRA:
      (_user_fcts->interp_tetra) (order,
                                  ho_vertex_num,
                                  local_to_user, 
                                  vertex_coords,
                                  point_coords,
                                  shape,
                                  stride_field,
                                  src_field,
                                  target_field);

      break;
      
    case FVMC_CELL_PRISM:

      (_user_fcts->interp_prism) (order,
                                  ho_vertex_num,
                                  local_to_user, 
                                  vertex_coords,
                                  point_coords,
                                  shape,
                                  stride_field,
                                  src_field,
                                  target_field);

      break;
      
    case FVMC_CELL_PYRAM:
      (_user_fcts->interp_pyramid ) (order,
                                     ho_vertex_num,
                                     local_to_user, 
                                     vertex_coords,
                                     point_coords,
                                     shape,
                                     stride_field,
                                     src_field,
                                     target_field);

      break;

    case FVMC_CELL_HEXA:
      (_user_fcts->interp_hexa) (order,
                                 ho_vertex_num,
                                 local_to_user, 
                                 vertex_coords,
                                 point_coords,
                                 shape,
                                 stride_field,
                                 src_field,
                                 target_field);
      break;

    default:

      bftc_error(__FILE__, __LINE__, 0,
                 _("fvmc_ho_interp_in_cell_3d : Not a high order 3D element type\n"));
    } 

  }

  else {

    _default_interp_in_cell_3d (type,
                                order,
                                ho_vertex_num,
                                local_to_user,
                                vertex_coords,
                                point_coords,
                                shape,
                                stride_field,
                                src_field,
                                target_field);

  }
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
  if (_user_fcts != NULL) {
    switch (type) {
      
    case FVMC_FACE_TRIA:
      (_user_fcts->interp_tria) (order,
                                 ho_vertex_num,
                                 local_to_user, 
                                 vertex_coords,
                                 point_coords,
                                 shape,
                                 stride_field,
                                 src_field,
                                 target_field);
      break;
      
    case FVMC_FACE_QUAD:
      (_user_fcts->interp_quad) (order,
                                 ho_vertex_num,
                                 local_to_user, 
                                 vertex_coords,
                                 point_coords,
                                 shape,
                                 stride_field,
                                 src_field,
                                 target_field);
      break;
      
    default:

      bftc_error(__FILE__, __LINE__, 0,
                 _("fvmc_ho_interp_on_cell_2d : Not a high order 2D element type\n"));
    }

  }

  else {

    _default_interp_on_cell_2d (type,
                                order,
                                ho_vertex_num,
                                local_to_user,
                                vertex_coords,
                                point_coords,
                                shape,
                                stride_field,
                                src_field,
                                target_field);
  }
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
  if (_user_fcts != NULL) {
    switch (type) {
      
    case FVMC_EDGE:
      (_user_fcts->interp_edge) (order,
                                 ho_vertex_num,
                                 local_to_user, 
                                 vertex_coords,
                                 point_coords,
                                 shape,
                                 stride_field,
                                 src_field,
                                 target_field);
      break;
      
    default:

      bftc_error(__FILE__, __LINE__, 0,
                 _("fvmc_ho_interp_on_cell_1d : Not a high order 1D element type\n"));
    }

  }

  else {

    _default_interp_on_cell_1d (type,
                                order,
                                ho_vertex_num,
                                local_to_user,
                                vertex_coords,
                                point_coords,
                                shape,
                                stride_field,
                                src_field,
                                target_field);

  }
}
