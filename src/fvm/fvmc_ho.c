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

#include "bftc_error.h"
#include "bftc_mem.h"
#include "bftc_printf.h"
#include "fvmc_point_location.h"

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
 *   n_node           <-- number of nodes
 *   ho_vertex_num    <-- high order vertex num (internal ordering)
 *   vertex_coords    <-- vertex coordinates
 *   point_coords     <-- point to locate coordinates
 *   projected_coords --> projected point coordinates if outside (or NULL)
 *   weights          --> interpolation weights in the element
 * 
 * return: 
 *   distance to the cell (distance <= 0 if point is inside)
 *
 *----------------------------------------------------------------------------*/

static double 
_default_location_in_cell_3d
(
 const fvmc_element_t type,
 const int order,
 const int n_node,
 const int *ho_vertex_num,
 const double *vertex_coords,
 const double *point_coords,
 double *projected_coords,
 double* weights
)
{
  type;
  order;
  n_node;
  ho_vertex_num;
  vertex_coords;
  point_coords;
  projected_coords;
  weights;
  
  double dist = 0.;
  bftc_error(__FILE__, __LINE__, 0,
             _("_default_location_in_cell_3d : Not implemented yet\n"));
  return dist;
}


/*----------------------------------------------------------------------------
 * 
 * Get uv coordinates of ho triangle nodes
 * 
 * parameters:
 *   order            <-- element order
 *   umin             <-- u min
 *   umax             <-- u max
 *   vmin             <-- v min
 *   vmax             <-- v max
 *   uv               --> uv (size = 2*n_node)
 *
 *----------------------------------------------------------------------------*/

static void
_uv_ho_tria_nodes
(
 const int order,
 const double umin,
 const double umax,
 const double vmin,
 const double vmax,
 double *uv
)
{
  int k = 0;
  int _order = order+1;

  double ustep = (umax - umin) / order;
  double vstep = (vmax - vmin) / order;
  
  for (int j = 0; j < _order; j++) {
    double v = vmin + j * vstep; 
    for (int i = 0; i < _order - j; i++) {
      double u = umin + i * ustep;

      uv[k++] = u;
      uv[k++] = v;
      
    }
  }
}


/*----------------------------------------------------------------------------
 * 
 * Get uv coordinates of ho quadrangle nodes
 * 
 * parameters:
 *   order            <-- element order
 *   umin             <-- u min
 *   umax             <-- u max
 *   vmin             <-- v min
 *   vmax             <-- v max
 *   uv               --> uv (size = 2*n_node)
 *
 *----------------------------------------------------------------------------*/

static void
_default_uv_ho_quad_nodes
(
 const int order,
 const double umin,
 const double umax,
 const double vmin,
 const double vmax,
 double *uv
)
{
  int k = 0;
  int _order = order+1;

  double ustep = (umax - umin) / order;
  double vstep = (vmax - vmin) / order;
  
  for (int j = 0; j < _order; j++) {
    double v = vmin + j * vstep; 
    for (int i = 0; i < _order; i++) {
      double u = umin + i * ustep;

      uv[k++] = u;
      uv[k++] = v;
      
    }
  }
}


/*----------------------------------------------------------------------------
 * 
 * Triangle Pn basis
 * 
 * parameters:
 *   order           <-- order
 *   u               <-- u
 *   v               <-- v
 *   weights         --> weights (size = n_nodes)
 *
 *
 *----------------------------------------------------------------------------*/

static void
_base_tria_pn
(
 const int order,
 const double u,
 const double v,
 double *weights
)
{

  if (order == 1) {

    weights[0] = 1. - u - v;
    weights[1] = u;
    weights[2] = v;

  }

  else if (order == 2) {

    double w  = 1. - u - v; 
    double u2 = 2. * u;
    double v2 = 2. * v;
    double w2 = 2. * w;
    
    weights[0] = w * (-1. + w2);  /* (i,j,k)=(0,0,2) */
    weights[1] = u2 * w2;         /* (i,j,k)=(1,0,1) */
    weights[2] = u * (-1. + u2);  /* (i,j,k)=(2,0,0) */
    weights[3] = u2 * v2;         /* (i,j,k)=(1,1,0) */
    weights[4] = v * (-1. + v2);  /* (i,j,k)=(0,2,0) */
    weights[5] = v2 * w2;         /* (i,j,k)=(0,1,1) */
  
  }

  else {

    bftc_error(__FILE__, __LINE__, 0,
               _("_base_tria_pn not yet implemented for order > 2\n"));

  }
}



/*----------------------------------------------------------------------------
 * 
 * Default point location on a high order triangle
 * 
 * parameters:
 *   order            <-- element order
 *   n_node           <-- number of nodes
 *   ho_vertex_num    <-- high order vertex num (internal ordering)
 *   vertex_coords    <-- vertex coordinates
 *   point_coords     <-- point to locate coordinates
 *   projected_coords --> projected point coordinates (or NULL)
 *   weights          --> interpolation weights in the element
 * 
 * return: 
 *   distance to the cell
 *
 *----------------------------------------------------------------------------*/

static double 
_default_location_on_tria_2d
(
 const int order,
 const int n_node,
 const int *ho_vertex_num,
 const double *vertex_coords,
 const double *point_coords,
 double *projected_coords,
 double *weights
)
{
  int _order = order;
  const double *_vertex_coords = vertex_coords;

  printf (" \n\n === _default_location_on_tria_2d === \n\n");
  
  if (1 == 1) {
    printf("order : %d\n", order);
    printf("n_node : %d\n", n_node);
    printf("ho_vertex_node : ");
    for (int i = 0; i < n_node; i++) {
      printf (" %d", ho_vertex_num[i]);
    }
    printf("\n");
    printf("vertex_coords : \n");
    for (int i = 0; i < n_node; i++) {
      int j = ho_vertex_num[i] - 1;
      printf ("%12.5e %12.5e %12.5e",
              vertex_coords[3*j],
              vertex_coords[3*j+1],
              vertex_coords[3*j+2]);
      printf("\n");
    }
    printf("\n");
    printf("point_coords : %12.5e %12.5e %12.5e\n",
           point_coords[0],
           point_coords[1],
           point_coords[2]);
  }
  
  /* Build sub-triangles */

  double uvP1[3];
  double weightsP1[3];
  double closest_pointP1[3];

  int selected_triaP1[3];
  
  double min_dist2 = HUGE_VAL;

  int ibeg = 0;
  int iend = _order;

  int proj_pt_in_selected_tria = 0;
  int itria = 0;
  
  for (int j = 0; j < _order; j++) {
    int k1 = 0;
    for (int i = ibeg; i < iend - 1; i++) {

      int _vtx1 = ho_vertex_num[i] - 1;
      int _vtx2 = ho_vertex_num[i + 1] - 1;
      int _vtx3 = ho_vertex_num[iend + 1 + k1] - 1;
      int _vtx4 = ho_vertex_num[iend + 2 + k1] - 1;

      double x1 = _vertex_coords[3*_vtx1];
      double y1 = _vertex_coords[3*_vtx1 + 1];
      double z1 = _vertex_coords[3*_vtx1 + 2];

      double x2 = _vertex_coords[3*_vtx2];
      double y2 = _vertex_coords[3*_vtx2 + 1];
      double z2 = _vertex_coords[3*_vtx2 + 2];

      double x3 = _vertex_coords[3*_vtx3];
      double y3 = _vertex_coords[3*_vtx3 + 1];
      double z3 = _vertex_coords[3*_vtx3 + 2];

      double x4 = _vertex_coords[3*_vtx4];
      double y4 = _vertex_coords[3*_vtx4 + 1];
      double z4 = _vertex_coords[3*_vtx4 + 2];
        
      double __vertex_coords[9] = {x1, y1, z1,
                                   x2, y2, z2,
                                   x3, y3, z3};
      double _uvP1[3];
      double _weightsP1[3];
      double _closest_pointP1[3];
      double _uvClosestPointP1[3];
      double _weightsClosestPointP1[3];
      double _dist2;

      printf ("   * itria : %d\n", ++itria);

      int proj_in_tria = fvmc_triangle_evaluate_Position ((double *)point_coords,
                                                          __vertex_coords,
                                                          _closest_pointP1,
                                                          _uvClosestPointP1,
                                                          _uvP1,
                                                          &_dist2,
                                                          _weightsClosestPointP1,
                                                          _weightsP1);

      printf ("     * _vertex_coords \n");
      for (int i1 = 0; i1 < 3; i1++) {
        printf ("        %12.5e %12.5e %12.5e",
                __vertex_coords[3*i1],
                __vertex_coords[3*i1+1],
                __vertex_coords[3*i1+2]);
        printf("\n");
      }

      printf("     * _closest_pointP1 : %12.5e %12.5e %12.5e\n",
             _closest_pointP1[0],
             _closest_pointP1[1],
             _closest_pointP1[2]);

      printf("     * _uvClosestPointP1 : %12.5e %12.5e\n",
             _uvClosestPointP1[0],
             _uvClosestPointP1[1]);
      
      printf("     * _uvP1 : %12.5e %12.5e\n",
             _uvP1[0],
             _uvP1[1]);

      printf("     * _dist2 : %12.5e\n", _dist2);

      
      printf("     * _weightsClosestPointP1 :");
      for (int i1 = 0; i1 < n_node; i1++) {
        printf (" %12.5e ", _weightsClosestPointP1[i1]);
      }
      printf("\n");

      printf("     * _weightsP1 :");
      for (int i1 = 0; i1 < n_node; i1++) {
        printf (" %12.5e ", _weightsP1[i1]);
      }
      printf("\n");

      if (_dist2 <= min_dist2) {
        min_dist2 = _dist2;
        for (int i1 = 0; i1 < 3; i1++) {
          uvP1[i1] = _uvClosestPointP1[i1];
          weightsP1[i1] = _weightsClosestPointP1[i1];
          closest_pointP1[i1] = _closest_pointP1[i1];
        }
        selected_triaP1[0] = _vtx1;
        selected_triaP1[1] = _vtx2;
        selected_triaP1[2] = _vtx3;
        proj_pt_in_selected_tria = proj_in_tria;
      }

      __vertex_coords[0] = x2;
      __vertex_coords[1] = y2;
      __vertex_coords[2] = z2;
      __vertex_coords[3] = x4;
      __vertex_coords[4] = y4;
      __vertex_coords[5] = z4;
      __vertex_coords[6] = x3;
      __vertex_coords[7] = y3;
      __vertex_coords[8] = z3;

      printf ("   * itria : %d\n", ++itria);

      proj_in_tria = fvmc_triangle_evaluate_Position ((double *) point_coords,
                                                       __vertex_coords,
                                                       _closest_pointP1,
                                                       _uvClosestPointP1,
                                                       _uvP1,
                                                       &_dist2,
                                                       _weightsClosestPointP1,
                                                       _weightsP1);
      printf ("     * _vertex_coords \n");
      for (int i1 = 0; i1 < 3; i1++) {
        printf ("        %12.5e %12.5e %12.5e",
                __vertex_coords[3*i1],
                __vertex_coords[3*i1+1],
                __vertex_coords[3*i1+2]);
        printf("\n");
      }

      printf("     * _closest_pointP1 : %12.5e %12.5e %12.5e\n",
             _closest_pointP1[0],
             _closest_pointP1[1],
             _closest_pointP1[2]);

      printf("     * _uvClosestPointP1 : %12.5e %12.5e\n",
             _uvClosestPointP1[0],
             _uvClosestPointP1[1]);
      
      printf("     * _uvP1 : %12.5e %12.5e\n",
             _uvP1[0],
             _uvP1[1]);

      printf("     * _dist2 : %12.5e\n", _dist2);

      
      printf("     * _weightsClosestPointP1 :");
      for (int i1 = 0; i1 < n_node; i1++) {
        printf (" %12.5e ", _weightsClosestPointP1[i1]);
      }
      printf("\n");

      printf("     * _weightsP1 :");
      for (int i1 = 0; i1 < n_node; i1++) {
        printf (" %12.5e ", _weightsP1[i1]);
      }
      printf("\n");

      if (_dist2 <= min_dist2) {
        min_dist2 = _dist2;
        for (int i1 = 0; i1 < 3; i1++) {
          uvP1[i1] = _uvClosestPointP1[i1];
          weightsP1[i1] = _weightsClosestPointP1[i1];
          closest_pointP1[i1] = _closest_pointP1[i1];
        }
        selected_triaP1[0] = _vtx2;
        selected_triaP1[1] = _vtx4;
        selected_triaP1[2] = _vtx3;
        proj_pt_in_selected_tria = proj_in_tria;
      }

      k1++;
    }

    int _vtx1 = ho_vertex_num[iend - 1] - 1;
    int _vtx2 = ho_vertex_num[iend - 1 + 1] - 1;
    int _vtx3 = ho_vertex_num[iend + 1 + k1]- 1;
      
    double x1 = vertex_coords[3*_vtx1];
    double y1 = vertex_coords[3*_vtx1 + 1];
    double z1 = vertex_coords[3*_vtx1 + 2];
      
    double x2 = vertex_coords[3*_vtx2];
    double y2 = vertex_coords[3*_vtx2 + 1];
    double z2 = vertex_coords[3*_vtx2 + 2];
      
    double x3 = vertex_coords[3*_vtx3];
    double y3 = vertex_coords[3*_vtx3 + 1];
    double z3 = vertex_coords[3*_vtx3 + 2];
      
    double __vertex_coords[9] = {x1, y1, z1,
                                x2, y2, z2,
                                x3, y3, z3};
      
    double _uvP1[3];
    double _weightsP1[3];
    double _closest_pointP1[3];
    double _uvClosestPointP1[3];
    double _weightsClosestPointP1[3];
   
    double _dist2;
      
    printf ("   * itria : %d\n", ++itria);

    int proj_in_tria =fvmc_triangle_evaluate_Position ((double *)point_coords,
                                                       __vertex_coords,
                                                       _closest_pointP1,
                                                       _uvClosestPointP1,
                                                       _uvP1,
                                                       &_dist2,
                                                       _weightsClosestPointP1,
                                                       _weightsP1);
      
      printf ("     * _vertex_coords \n");
      for (int i1 = 0; i1 < 3; i1++) {
        printf ("        %12.5e %12.5e %12.5e",
                __vertex_coords[3*i1],
                __vertex_coords[3*i1+1],
                __vertex_coords[3*i1+2]);
        printf("\n");
      }

      printf("     * _closest_pointP1 : %12.5e %12.5e %12.5e\n",
             _closest_pointP1[0],
             _closest_pointP1[1],
             _closest_pointP1[2]);

      printf("     * _uvClosestPointP1 : %12.5e %12.5e\n",
             _uvClosestPointP1[0],
             _uvClosestPointP1[1]);
      
      printf("     * _uvP1 : %12.5e %12.5e\n",
             _uvP1[0],
             _uvP1[1]);

      printf("     * _dist2 : %12.5e\n", _dist2);

      
      printf("     * _weightsClosestPointP1 :");
      for (int i1 = 0; i1 < n_node; i1++) {
        printf (" %12.5e ", _weightsClosestPointP1[i1]);
      }
      printf("\n");

      printf("     * _weightsP1 :");
      for (int i1 = 0; i1 < n_node; i1++) {
        printf (" %12.5e ", _weightsP1[i1]);
      }
      printf("\n");

      if (_dist2 <= min_dist2) {
      min_dist2 = _dist2;
      for (int i1 = 0; i1 < 3; i1++) {
        uvP1[i1] = _uvClosestPointP1[i1];
        weightsP1[i1] = _weightsClosestPointP1[i1];
        closest_pointP1[i1] = _closest_pointP1[i1];
      }
      selected_triaP1[0] = _vtx1;
      selected_triaP1[1] = _vtx2;
      selected_triaP1[2] = _vtx3;
      proj_pt_in_selected_tria = proj_in_tria;
    }

    ibeg = iend + 1;
    iend += _order - j;
  }

  /* uv proj P1 -> uv proj P2 */ 

  //_base_tria_pn (1    , uvP1[0], uvP1[1], weightsP1);

  double *uvNodes   = malloc (sizeof(double) * 2 * n_node);

  double uvP1inP2[2];
  
  _uv_ho_tria_nodes (order, 0., 1., 0, 1., uvNodes);

  for (int i = 0; i < 2; i++) {
    uvP1inP2[i] = weightsP1[0] * uvNodes[i] + weightsP1[1] * uvNodes[2+i];
  }
  
  free (uvNodes);
  
  _base_tria_pn (order   , uvP1inP2[0], uvP1inP2[1], weights);

  double _projected_coords[3];
  for (int j = 0; j < 3; j++) {
    _projected_coords[j] = 0;
  }

  for (int i = 0; i < n_node; i++) {

    const double *node_coords = vertex_coords + 3 * (ho_vertex_num[i] - 1);
    
    for (int j = 0; j < 3; j++) {
      _projected_coords[j] += weights[i] * node_coords[j]; 
    }
  }

  double dist2 = 0;

  for (int j = 0; j < 3; j++) {
    double comp = _projected_coords[j] - point_coords[j];
    dist2 += comp * comp; 
  }

  if (projected_coords != NULL) {
    for (int j = 0; j < 3; j++) {
      projected_coords[j] = _projected_coords[j];
    }
  }

  if (1 == 1) {
    printf(" --- resultats ---\n\n");
    printf("weights : ");
    for (int i = 0; i < n_node; i++) {
      printf (" %12.5e", weights[i]);
    }
    printf("\n");
    printf("min_dist2 : %12.5e\n", min_dist2);
    printf("projected_coords : %12.5e %12.5e %12.5e\n",
           projected_coords[0],
           projected_coords[1],
           projected_coords[2]);
  }
   
  
  printf ("\n\n === _default_location_on_tria_2d === \n\n");
  exit(1);
  return dist2;

 
}


/*----------------------------------------------------------------------------
 * 
 * Default point location on a high order quadrangle
 * 
 * parameters:
 *   order            <-- element order
 *   n_node           <-- number of nodes
 *   ho_vertex_num    <-- high order vertex num (internal ordering)
 *   vertex_coords    <-- vertex coordinates
 *   point_coords     <-- point to locate coordinates
 *   projected_coords --> projected point coordinates (or NULL)
 *   weights          --> interpolation weights in the element
 * 
 * return: 
 *   distance to the cell
 *
 *----------------------------------------------------------------------------*/

static double 
_default_location_on_quad_2d
(
 const int order,
 const int n_node,
 const int *ho_vertex_num,
 const double *vertex_coords,
 const double *point_coords,
 double *projected_coords,
 double* weights
)
{
  bftc_error(__FILE__, __LINE__, 0,
             _("_default_location_on_quad_2d : Not implemented yet\n"));

  int _order = order;

  double dist2 = HUGE_VAL;

  /* Build sub-triangles */

  double uvP1[3];
  double weightsP1[3];
  double closest_pointP1[3];

  int selected_triaP1[3];
  
  double min_dist2 = HUGE_VAL;

  int ibeg = 0;
  int iend = _order;

  for (int j = 0; j < _order; j++) {
    int k1 = 0;
      
    for (int i = ibeg; i < iend; i++) {

      int _vtx1 = ho_vertex_num[i] - 1;
      int _vtx2 = ho_vertex_num[i + 1] - 1;
      int _vtx3 = ho_vertex_num[iend + 1 + k1] - 1;
      int _vtx4 = ho_vertex_num[iend + 2 + k1] - 1;

      double x1 = vertex_coords[3*_vtx1];
      double y1 = vertex_coords[3*_vtx1 + 1];
      double z1 = vertex_coords[3*_vtx1 + 2];

      double x2 = vertex_coords[3*_vtx2];
      double y2 = vertex_coords[3*_vtx2 + 1];
      double z2 = vertex_coords[3*_vtx2 + 2];

      double x3 = vertex_coords[3*_vtx3];
      double y3 = vertex_coords[3*_vtx3 + 1];
      double z3 = vertex_coords[3*_vtx3 + 2];

      double x4 = vertex_coords[3*_vtx4];
      double y4 = vertex_coords[3*_vtx4 + 1];
      double z4 = vertex_coords[3*_vtx4 + 2];
        
      double _vertex_coords[9] = {x1, y1, z1,
                                  x2, y2, z2,
                                  x3, y3, z3};
      double _uvP1[3];
      double _weightsP1[3];
      double _closest_pointP1[3];
      double _uvClosestPointP1[3];
      double _weightsClosestPointP1[3];
      double _dist2;

      fvmc_triangle_evaluate_Position ((double *) point_coords,
                                       _vertex_coords,
                                       _closest_pointP1,
                                       _uvClosestPointP1,
                                       _uvP1,
                                       &_dist2,
                                       _weightsClosestPointP1,
                                       _weightsP1);

      if (_dist2 <= min_dist2) {
        min_dist2 = _dist2;
        for (int i1 = 0; i1 < 3; i1++) {
          uvP1[i1] = _uvClosestPointP1[i1];
          weightsP1[i1] = _weightsClosestPointP1[i1];
          closest_pointP1[i1] = _closest_pointP1[i1];
        }
        selected_triaP1[0] = _vtx1;
        selected_triaP1[1] = _vtx2;
        selected_triaP1[2] = _vtx3;
      }

      _vertex_coords[0] = x2;
      _vertex_coords[1] = y2;
      _vertex_coords[2] = z2;
      _vertex_coords[3] = x4;
      _vertex_coords[4] = y4;
      _vertex_coords[5] = z4;
      _vertex_coords[6] = x3;
      _vertex_coords[7] = y3;
      _vertex_coords[8] = z3;

      fvmc_triangle_evaluate_Position ((double *) point_coords,
                                       _vertex_coords,
                                       _closest_pointP1,
                                       _uvClosestPointP1,
                                       _uvP1,
                                       &_dist2,
                                       _weightsClosestPointP1,
                                       _weightsP1);

      if (_dist2 <= min_dist2) {
        min_dist2 = _dist2;
        for (int i1 = 0; i1 < 3; i1++) {
          uvP1[i1] = _uvClosestPointP1[i1];
          weightsP1[i1] = _weightsClosestPointP1[i1];
          closest_pointP1[i1] = _closest_pointP1[i1];
        }
        selected_triaP1[0] = _vtx1;
        selected_triaP1[1] = _vtx2;
        selected_triaP1[2] = _vtx3;
      }
        
      k1++;
    }

    ibeg = iend + 1;
    iend = ibeg + _order ;

  }

  /* uv P1 -> uv P2 */

  /* call baseTrianP1(u=uv1(1),v=uv1(2),ai=ai(1:3)) */
  
  /*         uv1_(1:2,iTria)= uvTriangle(1:2,tria1(1,iTria))*ai(1) & */
  /*         &               +uvTriangle(1:2,tria1(2,iTria))*ai(2) & */
  /*         &               +uvTriangle(1:2,tria1(3,iTria))*ai(3) */
  /*         call baseTrianP2(u=uv1_(1,iTria),v=uv1_(2,iTria),ai=ai(1:6)) */

  /* calcul des coordonnees et faire la nomrme L2 (voir base geom avec christophe) */

  return dist2;
  
}



/*----------------------------------------------------------------------------
 * 
 * Default point location on a high order cell 2d
 * 
 * parameters:
 *   type             <-- element type
 *   order            <-- element order
 *   n_node           <-- number of nodes
 *   ho_vertex_num    <-- high order vertex num (internal ordering)
 *   vertex_coords    <-- vertex coordinates
 *   point_coords     <-- point to locate coordinates
 *   projected_coords --> projected point coordinates (or NULL)
 *   weights          --> interpolation weights in the element
 * 
 * return: 
 *   distance to the cell
 *
 *----------------------------------------------------------------------------*/

static double 
_default_location_on_cell_2d
(
 const fvmc_element_t type,
 const int order,
 const int n_node,
 const int *ho_vertex_num,
 const double *vertex_coords,
 const double *point_coords,
 double *projected_coords,
 double* weights
)
{
  fvmc_element_t _type = type;
  double dist2;
  
  switch (_type) {

  case FVMC_FACE_TRIA:

    dist2 = _default_location_on_tria_2d (order,
                                          n_node,
                                          ho_vertex_num,
                                          vertex_coords,
                                          point_coords,
                                          projected_coords,
                                          weights);

    break;

  case FVMC_FACE_QUAD: 

    dist2 = _default_location_on_quad_2d (order,
                                          n_node,
                                          ho_vertex_num,
                                          vertex_coords,
                                          point_coords,
                                          projected_coords,
                                          weights);
    break;

  }
 
  return dist2;

}

/*----------------------------------------------------------------------------
 * 
 * Default point location on a high order cell 1d
 * 
 * parameters:
 *   type             <-- element type
 *   order            <-- element order
 *   n_node           <-- number of nodes
 *   ho_vertex_num    <-- high order vertex num (internal ordering)
 *   vertex_coords    <-- vertex coordinates
 *   point_coords     <-- point to locate coordinates
 *   projected_coords --> projected point coordinates (or NULL)
 *   weights          --> interpolation weights in the element
 * 
 * return: 
 *   distance to the cell
 *
 *----------------------------------------------------------------------------*/

static double 
_default_location_on_cell_1d
(
 const fvmc_element_t type,
 const int order,
 const int n_node,
 const int *ho_vertex_num,
 const double *vertex_coords,
 const double *point_coords,
 double *projected_coords,
 double* weights
)
{

  type;
  order;
  n_node;
  ho_vertex_num;
  vertex_coords;
  point_coords;
  projected_coords;
    
  double dist = 0.;
  bftc_error(__FILE__, __LINE__, 0,
             _("_default_location_on_cell_1d : Not implemented yet\n"));
  return dist;
}

/*----------------------------------------------------------------------------
 * 
 * Default field interpolation field to the target point_coords in a 3D cell
 * 
 * parameters:
 *   type              <-- element type
 *   order             <-- element order
 *   n_node            <-- number of nodes
 *   ho_vertex_num     <-- high order vertex num (internal ordering)
 *   local_to_user     <-- local to user ordering (for type)
 *   vertex_coords     <-- vertex coordinates
 *   point_coords      <-- point coordinates
 *   distance          <-- distance to the element
 *   point_proj_coords  <-- projected point coordinates
 *   weight             <-- weights
 *   stride_field      <-- field stride
 *   source_field      <-- source field (user ordering) 
 *   target_field      --> target field (defined to point_coords)
 *
 *----------------------------------------------------------------------------*/

static void 
_default_interp_in_cell_3d
(
 const fvmc_element_t type,
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
 double *target_field
 )                          
{
  for (int j = 0; j < stride_field; j++) {
    target_field[j] += 0.;
  }

  for (int i = 0; i < n_node; i++) {

    int i_node = ho_vertex_num[i] - 1;
    
    for (int j = 0; j < stride_field; j++) {
      target_field[j] += weight[i] * src_field[stride_field * i_node + j];
    }
  }
}


/*----------------------------------------------------------------------------
 * 
 * Default field interpolation field to the target point_coords on a 2D cell
 * 
 * parameters:
 *   type              <-- element type
 *   order             <-- element order
 *   n_node            <-- number of nodes
 *   ho_vertex_num     <-- high order vertex num (internal ordering)
 *   local_to_user     <-- local to user ordering (for type)
 *   vertex_coords     <-- vertex coordinates
 *   point_coords      <-- point coordinates
 *   distance          <-- distance to the element
 *   point_proj_coords  <-- projected point coordinates
 *   weight             <-- weights
 *   stride_field      <-- field stride
 *   source_field      <-- source field (user ordering) 
 *   target_field      --> target field (defined to point_coords)
 *
 *----------------------------------------------------------------------------*/

static void 
_default_interp_on_cell_2d
(
 const fvmc_element_t type,
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
 double *target_field
)
{
  for (int j = 0; j < stride_field; j++) {
    target_field[j] += 0.;
  }

  for (int i = 0; i < n_node; i++) {

    int i_node = ho_vertex_num[i] - 1;
    
    for (int j = 0; j < stride_field; j++) {
      target_field[j] += weight[i] * src_field[stride_field * i_node + j];
    }
  }
}


/*----------------------------------------------------------------------------
 * 
 * Default field interpolation field to the target point_coords on a 1D cell
 * 
 * parameters:
 *   type              <-- element type
 *   order             <-- element order
 *   n_node            <-- number of nodes
 *   ho_vertex_num     <-- high order vertex num (internal ordering)
 *   local_to_user     <-- local to user ordering (for type)
 *   vertex_coords     <-- vertex coordinates
 *   point_coords      <-- point coordinates
 *   distance          <-- distance to the element
 *   point_proj_coords  <-- projected point coordinates
 *   weight             <-- weights
 *   stride_field      <-- field stride
 *   source_field      <-- source field (user ordering) 
 *   target_field      --> target field (defined to point_coords)
 *
 *----------------------------------------------------------------------------*/

static void 
_default_interp_on_cell_1d
(
 const fvmc_element_t type,
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
 double *target_field
)
{
  for (int j = 0; j < stride_field; j++) {
    target_field[j] += 0.;
  }

  for (int i = 0; i < n_node; i++) {

    int i_node = ho_vertex_num[i] - 1;
    
    for (int j = 0; j < stride_field; j++) {
      target_field[j] += weight[i] * src_field[stride_field * i_node + j];
    }
  }
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
 *   n_node            <-- number of nodes
 *   ho_vertex_num    <-- high order vertex num (internal ordering)
 *   vertex_coords    <-- vertex coordinates
 *   point_coords     <-- point to locate coordinates
 *   projected_coords --> projected point coordinates if outside (or NULL)
 *   weights          --> interpolation weights in the element
 * 
 * return: 
 *   distance to the cell (distance <= 0 if point is inside)
 *
 *----------------------------------------------------------------------------*/

double 
fvmc_ho_location_in_cell_3d
(
 const fvmc_element_t type,
 const int order,
 const int n_node,
 const int *ho_vertex_num,
 const double *vertex_coords,
 const double *point_coords,
 double *projected_coords,
 double *weights
)
{

  if (_user_fcts != NULL) {

    switch (type) {
      
    case FVMC_CELL_TETRA:
      return (_user_fcts->location_tetra) (order,
                                           n_node,
                                           ho_vertex_num,
                                           vertex_coords,
                                           point_coords,
                                           projected_coords,
                                           weights);
      break;
      
    case FVMC_CELL_PRISM:
      return (_user_fcts->location_prism) (order,
                                           n_node,
                                           ho_vertex_num,
                                           vertex_coords,
                                           point_coords,
                                           projected_coords,
                                           weights);
      break;
      
    case FVMC_CELL_PYRAM:
      return (_user_fcts->location_pyramid) (order,
                                             n_node,
                                             ho_vertex_num,
                                             vertex_coords,
                                             point_coords,
                                             projected_coords,
                                             weights);
      break;

    case FVMC_CELL_HEXA:
      return (_user_fcts->location_hexa) (order,
                                          n_node,
                                          ho_vertex_num,
                                          vertex_coords,
                                          point_coords,
                                          projected_coords,
                                          weights);
      break;

    default:

      bftc_error(__FILE__, __LINE__, 0,
                 _("fvmc_ho_location_in_cell_3d : Not a high order 3D element type\n"));
    } 
  }

  else {

    return _default_location_in_cell_3d (type,
                                         order,
                                         n_node,
                                         ho_vertex_num,
                                         vertex_coords,
                                         point_coords,
                                         projected_coords,
                                         weights);
    
  }

  return HUGE_VAL;
}

/*----------------------------------------------------------------------------
 * 
 * Point location on a high order cell 2d
 * 
 * parameters:
 *   type             <-- element type
 *   order            <-- element order
 *   n_node            <-- number of nodes
 *   ho_vertex_num    <-- high order vertex num (internal ordering)
 *   vertex_coords    <-- vertex coordinates
 *   point_coords     <-- point to locate coordinates
 *   projected_coords --> projected point coordinates (or NULL)
 *   weights          --> interpolation weights in the element
 * 
 * return: 
 *   distance to the cell
 *
 *----------------------------------------------------------------------------*/

double 
fvmc_ho_location_on_cell_2d
(
 const fvmc_element_t type,
 const int order,
 const int n_node,
 const int *ho_vertex_num,
 const double *vertex_coords,
 const double *point_coords,
 double *projected_coords,
 double *weights
)
{

  double dist2 = HUGE_VAL;
  
  if (_user_fcts != NULL) {
    switch (type) {
      
    case FVMC_FACE_TRIA:
      dist2 = (_user_fcts->location_tria) (order,
                                           n_node,
                                           ho_vertex_num,
                                           vertex_coords,
                                           point_coords,
                                           projected_coords,
                                           weights);
      break;
      
    case FVMC_FACE_QUAD:
      dist2 =  (_user_fcts->location_quad) (order,
                                            n_node,
                                            ho_vertex_num,
                                            vertex_coords,
                                            point_coords,
                                            projected_coords,
                                            weights);
      break;
      
    default:

      bftc_error(__FILE__, __LINE__, 0,
                 _("fvmc_ho_location_on_cell_2d : Not a high order 2D element type\n"));
    }
  }

  else {

    dist2 = _default_location_on_cell_2d (type,
                                          order,
                                          n_node,
                                          ho_vertex_num,
                                          vertex_coords,
                                          point_coords,
                                          projected_coords,
                                          weights);
    
  }
  return dist2;
}

/*----------------------------------------------------------------------------
 * 
 * Point location on a high order cell 1d
 * 
 * parameters:
 *   type             <-- element type
 *   order            <-- element order
 *   n_node            <-- number of nodes
 *   ho_vertex_num    <-- high order vertex num (internal ordering)
 *   vertex_coords    <-- vertex coordinates
 *   point_coords     <-- point to locate coordinates
 *   projected_coords --> projected point coordinates (or NULL)
 *   weights          --> interpolation weights in the element
 * 
 * return: 
 *   distance to the cell
 *
 *----------------------------------------------------------------------------*/

double 
fvmc_ho_location_on_cell_1d
(
 const fvmc_element_t type,
 const int order,
 const int n_node,
 const int *ho_vertex_num,
 const double *vertex_coords,
 const double *point_coords,
 double *projected_coords,
 double *weights
)
{

  if (_user_fcts != NULL) {
    switch (type) {
      
    case FVMC_EDGE:
      return (_user_fcts->location_edge) (order,
                                          n_node,
                                          ho_vertex_num,
                                          vertex_coords,
                                          point_coords,
                                          projected_coords,
                                          weights);
      break;
      
    default:

      bftc_error(__FILE__, __LINE__, 0,
                 _("fvmc_ho_location_on_cell_1d : Not a high order 1D element type\n"));
    }
  }

  else {

    return _default_location_on_cell_1d (type,
                                         order,
                                         n_node,
                                         ho_vertex_num,
                                         vertex_coords,
                                         point_coords,
                                         projected_coords,
                                         weights);

  }

  return HUGE_VAL;
}


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
 *   point_coords      <-- point coordinates 
 *   distance          <-- distance to the element
 *   point_proj_coords  <-- projected point coordinates
 *   weight             <-- weights
 *   stride_field      <-- field stride
 *   source_field      <-- source field (user ordering) 
 *   target_field      --> target field (defined to point_coords)
 *
 *----------------------------------------------------------------------------*/

void 
fvmc_ho_interp_in_cell_3d
(
 const fvmc_element_t type,
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
 double *target_field
 )                          
{
  if (_user_fcts != NULL) {

    switch (type) {
      
    case FVMC_CELL_TETRA:
      (_user_fcts->interp_tetra) (order,
                                  n_node,
                                  ho_vertex_num,
                                  local_to_user, 
                                  vertex_coords,
                                  point_coords,
                                  distance,
                                  point_proj_coords,
                                  weight,
                                  stride_field,
                                  src_field,
                                  target_field);

      break;
      
    case FVMC_CELL_PRISM:

      (_user_fcts->interp_prism) (order,
                                  n_node,
                                  ho_vertex_num,
                                  local_to_user, 
                                  vertex_coords,
                                  point_coords,
                                  distance,
                                  point_proj_coords,
                                  weight,
                                  stride_field,
                                  src_field,
                                  target_field);

      break;
      
    case FVMC_CELL_PYRAM:
      (_user_fcts->interp_pyramid ) (order,
                                     n_node,
                                     ho_vertex_num,
                                     local_to_user, 
                                     vertex_coords,
                                     point_coords,
                                     distance,
                                     point_proj_coords,
                                     weight,
                                     stride_field,
                                     src_field,
                                     target_field);

      break;

    case FVMC_CELL_HEXA:
      (_user_fcts->interp_hexa) (order,
                                 n_node,
                                 ho_vertex_num,
                                 local_to_user, 
                                 vertex_coords,
                                 point_coords,
                                 distance,
                                 point_proj_coords,
                                 weight,
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
                                n_node,
                                ho_vertex_num,
                                local_to_user,
                                vertex_coords,
                                point_coords,
                                distance,
                                point_proj_coords,
                                weight,
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
 *   n_node            <-- number of nodes
 *   ho_vertex_num     <-- high order vertex num (internal ordering)
 *   local_to_user     <-- local to user ordering (for type)
 *   vertex_coords     <-- vertex coordinates
 *   point_coords      <-- point coordinates 
 *   distance          <-- distance to the element
 *   point_proj_coords  <-- projected point coordinates
 *   weight             <-- weights
 *   stride_field      <-- field stride
 *   source_field      <-- source field (user ordering) 
 *   target_field      --> target field (defined to point_coords)
 *
 *----------------------------------------------------------------------------*/

void 
fvmc_ho_interp_on_cell_2d
(
 const fvmc_element_t type,
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
 double *target_field
)
{
  if (_user_fcts != NULL) {
    switch (type) {
      
    case FVMC_FACE_TRIA:
      (_user_fcts->interp_tria) (order,
                                 n_node,
                                 ho_vertex_num,
                                 local_to_user, 
                                 vertex_coords,
                                 point_coords,
                                 distance,
                                 point_proj_coords,
                                 weight,
                                 stride_field,
                                 src_field,
                                 target_field);
      break;
      
    case FVMC_FACE_QUAD:
      (_user_fcts->interp_quad) (order,
                                 n_node,
                                 ho_vertex_num,
                                 local_to_user, 
                                 vertex_coords,
                                 point_coords,
                                 distance,
                                 point_proj_coords,
                                 weight,
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
                                n_node,
                                ho_vertex_num,
                                local_to_user,
                                vertex_coords,
                                point_coords,
                                distance,
                                point_proj_coords,
                                weight,
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
 *   n_node            <-- number of nodes
 *   ho_vertex_num     <-- high order vertex num (internal ordering)
 *   local_to_user     <-- local to user ordering (for type)
 *   vertex_coords     <-- vertex coordinates
 *   point_coords      <-- point coordinates 
 *   distance          <-- distance to the element
 *   point_proj_coords  <-- projected point coordinates
 *   weight             <-- weights
 *   stride_field      <-- field stride
 *   source_field      <-- source field (user ordering) 
 *   target_field      --> target field (defined to point_coords)
 *
 *----------------------------------------------------------------------------*/

void 
fvmc_ho_interp_on_cell_1d
(
 const fvmc_element_t type,
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
 double *target_field
 )
{
  if (_user_fcts != NULL) {
    switch (type) {
      
    case FVMC_EDGE:
      (_user_fcts->interp_edge) (order,
                                 n_node,
                                 ho_vertex_num,
                                 local_to_user, 
                                 vertex_coords,
                                 point_coords,
                                 distance,
                                 point_proj_coords,
                                 weight,
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
                                n_node,
                                ho_vertex_num,
                                local_to_user,
                                vertex_coords,
                                point_coords,
                                distance,
                                point_proj_coords,
                                weight,
                                stride_field,
                                src_field,
                                target_field);

  }
}
