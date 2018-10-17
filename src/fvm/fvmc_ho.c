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
#include "cwipi_config.h"

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

#if defined (HAVE_SPACE_BASIS) 
#include "spacebasis.h"
#endif

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

#define S_HEAP 50

/*============================================================================
 * Type definitions
 *============================================================================*/

typedef struct _fvmc_ho_user_elt_t {
  fvmc_ho_basis_fct_t elt_basis;
  fvmc_ho_xsi_fct_t   xsi_coords;
  fvmc_ho_location_fct_t location_in_elt;
} fvmc_ho_user_elt_t;

/*----------------------------------------------------------------------------
 * Sorted heap for sub-triangle storage
 *----------------------------------------------------------------------------*/

typedef struct {

  int idx;

  double vtx_tria[3*3*S_HEAP];
  double uvInPn_tria[2*3*S_HEAP];

  double closest_pt[3*S_HEAP]; 
  double closest_pt_uvP1[2*S_HEAP];
  double closest_pt_uvInPn[2*S_HEAP];
  double dist2[S_HEAP];

  int child[S_HEAP];

  int free_idx[S_HEAP];
  int sorted_idx[S_HEAP];  

} _heap_t;


/*----------------------------------------------------------------------------
 * Function pointer to define a initial tesselation and push it in the heap
 *
 * parameters:
 *   heap              <-> Heap
 *   order             <-- element order
 *   n_nodes           <-- number of nodes
 *   local_to_user     <-- local to user ordering (for type)
 *   nodes_coords      <-- nodes coordinates
 *   point_coords      <-- point coordinates 
 *
 *----------------------------------------------------------------------------*/

typedef void (*_heap_fill_init_sub_tria_t)
(
 _heap_t  *heap,
 const int order,
 const int n_nodes,
 const double *nodes_coords,
 const double *point_coords
);


/*----------------------------------------------------------------------------
 * Function pointer to define a basis for a 2D element
 *
 * parameters:
 *   order           <-- order
 *   n_pts           <-- number of points
 *   uv              <-- uv coordinates of points
 *   weights         --> weights (size = n_nodes)
 *
 *----------------------------------------------------------------------------*/

typedef void (*_basis_generic_2D_t)
(
 const int order,
 const int n_pts,
 const double *uv,
 double *weights
);


/*============================================================================
 * Static global variables
 *============================================================================*/

static fvmc_ho_user_elt_t *_user_edge = NULL;

static fvmc_ho_user_elt_t *_user_tria = NULL;
static fvmc_ho_user_elt_t *_user_quad = NULL;

static fvmc_ho_user_elt_t *_user_tetra = NULL;
static fvmc_ho_user_elt_t *_user_hexa = NULL;
static fvmc_ho_user_elt_t *_user_prism = NULL;
static fvmc_ho_user_elt_t *_user_pyra = NULL;

static int                 _idebug = 0;

static int                 _n_ijk_tria_space = 0;
static int               **_ijk_tria_space = NULL;

static int                 _n_ijk_quad_space = 0;
static int               **_ijk_quad_space = NULL;

static int                 _n_ijk_1D_space = 0;
static int               **_ijk_1D_space = NULL;

static int isWarningPrinted1 = 0;
static int isWarningPrinted2 = 0;

/*============================================================================
 * Private function definitions
 *============================================================================*/

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
_uv_ho_quad_nodes
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
 *   n_pts           <-- number of points
 *   uv              <-- u (size = 2 * n_pts)
 *   weights         --> weights (size = n_nodes * n_pts)
 *
 *
 *----------------------------------------------------------------------------*/

static void
_basis_tria_pn
(
 const int order,
 const int n_pts,
 const double *uv,
 double *weights
)
{

  if (order == 1) {

    if (n_pts != 1) {

      double *u  = malloc (sizeof(double) *n_pts);
      double *v  = malloc (sizeof(double) *n_pts);
      
      for (int i = 0; i < n_pts; i++) {
        u[i]  = uv[2*i];
        v[i]  = uv[2*i+1];
      }
      
      for (int i = 0; i < n_pts; i++) {
        weights[3*i]   = 1. - u[i] - v[i];
        weights[3*i+1] = u[i];
        weights[3*i+2] = v[i];
      }
      
      free (u);
      free (v);
    }

    else {

      double u = uv[0];
      double v = uv[1];
      
      weights[0] = 1. - u - v;
      weights[1] = u;
      weights[2] = v;

    }
  }
  
  else if (order == 2) {

    if (n_pts != 1) {
      double *u  = malloc (sizeof(double) *n_pts);
      double *v  = malloc (sizeof(double) *n_pts);
      double *w  = malloc (sizeof(double) *n_pts);
      double *u2 = malloc (sizeof(double) *n_pts);
      double *v2 = malloc (sizeof(double) *n_pts);
      double *w2 = malloc (sizeof(double) *n_pts);
    

      for (int i = 0; i < n_pts; i++) {
        u[i]  = uv[2*i];
        v[i]  = uv[2*i+1];
        w[i]  = 1. - u[i] - v[i];
        u2[i] = 2. * u[i];
        v2[i] = 2. * v[i];
        w2[i] = 2. * w[i];
      }
      
      for (int i = 0; i < n_pts; i++) {
        
        weights[6*i+0] = w[i] * (-1. + w2[i]);  /* (i,j,k)=(0,0,2) */
        weights[6*i+1] = u2[i] * w2[i];         /* (i,j,k)=(1,0,1) */
        weights[6*i+2] = u[i] * (-1. + u2[i]);  /* (i,j,k)=(2,0,0) */
        weights[6*i+3] = v2[i] * w2[i];         /* (i,j,k)=(0,1,1) */
        weights[6*i+4] = u2[i] * v2[i];         /* (i,j,k)=(1,1,0) */
        weights[6*i+5] = v[i] * (-1. + v2[i]);  /* (i,j,k)=(0,2,0) */
      }
    
      free (u);
      free (v);
      free (w);
      free (u2);
      free (v2);
      free (w2);
    }
    else {

      double u = uv[0];
      double v = uv[1];
      
      double w  = 1. - u - v;
      double u2 = 2. * u;
      double v2 = 2. * v;
      double w2 = 2. * w;
      
      weights[0] = w * (-1. + w2);  /* (i,j,k)=(0,0,2) */
      weights[1] = u2 * w2;         /* (i,j,k)=(1,0,1) */
      weights[2] = u * (-1. + u2);  /* (i,j,k)=(2,0,0) */
      weights[3] = v2 * w2;         /* (i,j,k)=(0,1,1) */
      weights[4] = u2 * v2;         /* (i,j,k)=(1,1,0) */
      weights[5] = v * (-1. + v2);  /* (i,j,k)=(0,2,0) */

    }
  }

  /* else if (order == 3) { */

  /*   SNB_setT3Basis_P3 (n_pts, uv, weights); */

  /* } */
  
  else {

#if defined (HAVE_SPACE_BASIS)
    
    if (_ijk_tria_space == NULL) {
      _n_ijk_tria_space = FVMC_MAX (order-1, 10);
      _ijk_tria_space = malloc (sizeof(int *) * _n_ijk_tria_space);
      for (int i = 0; i < _n_ijk_tria_space; i++) {
        _ijk_tria_space[i] = NULL;
      }
    }

    if (order > _n_ijk_tria_space) {
      _ijk_tria_space = realloc (_ijk_tria_space, sizeof(int *) * _n_ijk_tria_space);
      for (int i = _n_ijk_tria_space; i < order; i++) {
        _ijk_tria_space[i] = NULL;
      }
    }
    
    int *__ijk_tria_space = _ijk_tria_space[order-1];

    if (__ijk_tria_space == NULL) {
      const int n_nodes = (order+1)*(order+2)/2;
      __ijk_tria_space = malloc (sizeof(int) * 2 * n_nodes);
      int k = 0;
      for (int j = 0; j < order+1; j++) {
        for (int i = 0; i < order+1-j; i++) {
          __ijk_tria_space[2*k]   = i;
          __ijk_tria_space[2*k+1] = j;
          k += 1;
        }
      }
    }

    SNB_setT3BasisEqui_uv (order, n_pts, __ijk_tria_space, uv, weights);

#else
    bftc_error(__FILE__, __LINE__, 0,
               _("_basis_tria_pn not yet implemented for order > 2 without space basis \n"));
#endif

  }
}




/*----------------------------------------------------------------------------
 * 
 * Quadrangle Qn basis
 * 
 * parameters:
 *   order           <-- order
 *   n_pts           <-- number of points
 *   uv              <-- u (size = 2 * n_pts)
 *   weights         --> weights (size = n_nodes * n_pts)
 *
 *
 *----------------------------------------------------------------------------*/


static void
_basis_quad_qn
(
 const int order,
 const int n_pts,
 const double *uv,
 double *weights
)
{
  
  if (order == 1) {

    double *u = malloc (sizeof(double) *n_pts);
    double *v = malloc (sizeof(double) *n_pts);
    double *u1 = malloc (sizeof(double) *n_pts);
    double *v1 = malloc (sizeof(double) *n_pts);

    for (int i = 0; i < n_pts; i++) {
      u[i] = uv[2*i];
      v[i] = uv[2*i+1];
      u1[i] = (1 - u[i]);
      v1[i] = (1 - v[i]);
    }
    
    for (int i = 0; i < n_pts; i++) {
      weights[4*i+0] = u1[i] * v1[i];
      weights[4*i+1] = u[i] * v1[i];
      weights[4*i+2] = u1[i] * v[i];
      weights[4*i+3] = u[i] * v[i];
    }

    free(u1);
    free(v1);
    free(u);
    free(v);

  }

  else if (order == 2) {
    
    double *u = malloc (sizeof(double) *n_pts);
    double *v = malloc (sizeof(double) *n_pts);

    double *uM = malloc (sizeof(double) *n_pts);
    double *uP = malloc (sizeof(double) *n_pts);
    double *u0 = malloc (sizeof(double) *n_pts);

    double *au1 = malloc (sizeof(double) *n_pts);
    double *au2 = malloc (sizeof(double) *n_pts);
    double *au3 = malloc (sizeof(double) *n_pts);
    
    double *vM = malloc (sizeof(double) *n_pts);
    double *vP = malloc (sizeof(double) *n_pts);
    double *v0 = malloc (sizeof(double) *n_pts);

    double *av1 = malloc (sizeof(double) *n_pts);
    double *av2 = malloc (sizeof(double) *n_pts);
    double *av3 = malloc (sizeof(double) *n_pts);

    
    for (int i = 0; i < n_pts; i++) {
      u[i] = uv[2*i];
      v[i] = uv[2*i+1];

      uM[i] = 2*(1-u[i]);
      uP[i] = 2*u[i];
      u0[i] = u[i]-0.5;
      
      au1[i] = -uM[i] * u0[i]; 
      au2[i] =  uM[i] * uP[i];
      au3[i] =  u0[i] * uP[i];
    
      vM[i] = 2*(1-v[i]);
      vP[i] = 2*v[i];
      v0[i] = v[i]-0.5;

      av1[i] = -vM[i] * v0[i]; 
      av2[i] =  vM[i] * vP[i];
      av3[i] =  v0[i] * vP[i];
    }
    
    for (int i = 0; i < n_pts; i++) {
      weights[9*i+0]=au1[i]*av1[i];
      weights[9*i+1]=au2[i]*av1[i];
      weights[9*i+2]=au3[i]*av1[i];
      weights[9*i+3]=au1[i]*av2[i];
      weights[9*i+4]=au2[i]*av2[i];
      weights[9*i+5]=au3[i]*av2[i];
      weights[9*i+6]=au1[i]*av3[i];
      weights[9*i+7]=au2[i]*av3[i];
      weights[9*i+8]=au3[i]*av3[i];
    }
    
    free(u);
    free(v);

    free(uM);
    free(uP);
    free(u0);

    free(au1);
    free(au2);
    free(au3);
    
    free(vM);
    free(vP);
    free(v0);

    free(av1);
    free(av2);
    free(av3);

  }

  else {

#if defined (HAVE_SPACE_BASIS)
    
    if (_ijk_quad_space == NULL) {
      _n_ijk_quad_space = FVMC_MAX (order-1, 10);
      _ijk_quad_space = malloc (sizeof(int *) * _n_ijk_quad_space);
      for (int i = 0; i < _n_ijk_quad_space; i++) {
        _ijk_quad_space[i] = NULL;
      }
    }

    if (order > _n_ijk_quad_space) {
      _ijk_quad_space = realloc (_ijk_quad_space, sizeof(int *) * _n_ijk_quad_space);
      for (int i = _n_ijk_quad_space; i < order; i++) {
        _ijk_quad_space[i] = NULL;
      }
    }
    
    int *__ijk_quad_space = _ijk_quad_space[order-1];

    if (__ijk_quad_space == NULL) {
      const int n_nodes = (order+1)*(order+1);
      __ijk_quad_space = malloc (sizeof(int) * 2 * n_nodes);
      int k = 0;
      for (int j = 0; j < order+1; j++) {
        for (int i = 0; i < order+1; i++) {
          __ijk_quad_space[k++] = i;
          __ijk_quad_space[k++] = j;
        }
      }
    }

    double *u = malloc (sizeof(double) *n_pts);
    double *v = malloc (sizeof(double) *n_pts);

    for (int i = 0; i < n_pts; i++) {
      u[i] = 2 * uv[2*i]   - 1;
      v[i] = 2 * uv[2*i+1] - 1;
    }
    
    SNB_setQ4BasisEqui_uv (order, n_pts, __ijk_quad_space,
                           u, v,
                           weights);

    free (u);
    free (v);
    
#else
    bftc_error(__FILE__, __LINE__, 0,
               _("_basis_quad_pn not yet implemented for order > 2 without space basis \n"));
#endif
  }
}
  
/*----------------------------------------------------------------------------
 * 
 * Init heap
 * 
 * parameters:
 *   heap             <-- heap to initialize
 *
 *----------------------------------------------------------------------------*/

static void
_heap_init
(
_heap_t *heap
)
{
  heap->idx = -1;

  for (int i = 0; i < S_HEAP; i++) {
    heap->free_idx[i] = i;
  }
}


/*----------------------------------------------------------------------------
 * 
 * Get top of the heap
 * 
 * parameters:
 *   heap             <-> heap to initialize
 *   order            <-> element order
 *   n_node           <-> number of nodes
 *   ho_vertex_num    <-> high order vertex num (internal ordering)
 *   vertex_coords    <-> vertex coordinates
 *   point_coords     <-> point to locate coordinates
 *
 *----------------------------------------------------------------------------*/

static int
_heap_top_get
(
 _heap_t *heap,
 double *vtx_tria_current,
 double *uvInPn_tria_current,
 double *closest_pt_current,
 double *closest_pt_uvP1_current,
 double *closest_pt_uvInPn_current,
 double *dist2_current,
 int *child
)
{
  if (heap->idx < 0) {
    return 1;
  }

  int idx = heap->sorted_idx[heap->idx];
  heap->free_idx[heap->idx--]= idx;
  
  double *_vtx_tria_current = heap->vtx_tria + 9 * idx;
  for (int i = 0; i < 9; i++) {
    vtx_tria_current[i] = _vtx_tria_current[i];
  }

  double *_uvInPn_tria_current = heap->uvInPn_tria + 6 *idx;
  for (int i = 0; i < 6; i++) {
    uvInPn_tria_current[i] = _uvInPn_tria_current[i];
  }

  double *_closest_pt_current = heap->closest_pt + 3 *idx;
  for (int i = 0; i < 3; i++) {
    closest_pt_current[i] = _closest_pt_current[i];
  }

  double *_closest_pt_uvP1_current = heap->closest_pt_uvP1 + 2 *idx;
  double *_closest_pt_uvInPn_current = heap->closest_pt_uvInPn + 2 *idx;
  for (int i = 0; i < 2; i++) {
    closest_pt_uvP1_current[i] = _closest_pt_uvP1_current[i];
    closest_pt_uvInPn_current[i] = _closest_pt_uvInPn_current[i];
  }

  *child = heap->child[idx];
  *dist2_current = heap->dist2[idx];

  return 0;
}

/*----------------------------------------------------------------------------
 * 
 * Add sub-triangles of a pn-triangle in the heap
 * 
 * parameters:
 *   heap             <-- heap to initialize
 *   order            <-- element order
 *   n_node           <-- number of nodes
 *   ho_vertex_num    <-- high order vertex num (internal ordering)
 *   vertex_coords    <-- vertex coordinates
 *   point_coords     <-- point to locate coordinates
 *
 *----------------------------------------------------------------------------*/

static void
_heap_insert
(
 _heap_t *heap,
 double *vtx_tria,
 double *uvInPn_tria,
 double *closest_pt,
 double *closest_pt_uvP1,
 double *closest_pt_uvInPn,
 double dist2,
 int child
)
{
  /* Look for index (dicothomy) */

  if (0 == 1) {
    printf ("distances in heap deb :");
    for (int i = 0; i < heap->idx + 1; i++) {
      int _idx2 = heap->sorted_idx[i];
      printf (" %12.5e", heap->dist2[_idx2]);

    }
    printf ("\n");
  }
  
  int curr_idx = heap->idx;
  int *sorted_idx = heap->sorted_idx;
  double *sorted_dist2 = heap->dist2;

  int beg = 0;
  int end = curr_idx;

  while (beg <= end) {
    double dist2_beg = sorted_dist2[sorted_idx[beg]]; 
    double dist2_end = sorted_dist2[sorted_idx[end]];
    
    if (dist2 >= dist2_beg) {
      end = beg - 1;
    }
    
    else if (dist2 <= dist2_end) {
      beg = end + 1;
    }
    
    else {
      
      const int middle = (end + beg) / 2;
      if (beg == middle) {
        beg = beg + 1;
        end = beg - 1;
      }
      else {
        const double dist2_middle = sorted_dist2[sorted_idx[middle]];
        if (dist2 > dist2_middle) {
          end = middle;
        }
        else {
          beg = middle;
        }
      }
    }
  }


  /* If the heap is full remove the most distant */

  if (curr_idx >= (S_HEAP - 1)) {
    if (beg == 0) {
      return;
    }
    else {
      const int free_idx = sorted_idx[0];

      for (int i = 1; i < S_HEAP; i++) {
        sorted_idx[i-1] = sorted_idx[i]; 
      }

      heap->free_idx[heap->idx] = free_idx;
      heap->idx--;
      
      beg = beg - 1;
      end = end - 1;
      
    }
  }

  /* Add the element to the heap */

  heap->idx++;
  assert (heap->free_idx[heap->idx] != -1);

  for (int j = heap->idx; j > beg; j--) {
    sorted_idx[j] = sorted_idx[j-1];
  }

  sorted_idx[beg] = heap->free_idx[heap->idx];
  
  heap->free_idx[heap->idx] = -1;

  int _idx = sorted_idx[beg];

  for (int j = 0; j < 9; j++) {
    heap->vtx_tria[9*_idx+j] = vtx_tria[j];
  }

  for (int j = 0; j < 6; j++) {
    heap->uvInPn_tria[6*_idx+j] = uvInPn_tria[j]; 
  }

  for (int j = 0; j < 3; j++) {
    heap->closest_pt[3*_idx+j]  = closest_pt[j]; 
  }

  for (int j = 0; j < 2; j++) {
    heap->closest_pt_uvP1[2*_idx+j] = closest_pt_uvP1[j];
    heap->closest_pt_uvInPn[2*_idx+j] = closest_pt_uvInPn[j];
  }

  heap->dist2[_idx] = dist2;
  heap->child[_idx] = child;
  

  if (0 == 1) {
    printf ("distances in heap fin :");
    for (int i = 0; i < heap->idx + 1; i++) {
      int _idx2 = sorted_idx[i];
      printf (" %12.5e", heap->dist2[_idx2]);
    }
    printf ("\n");
    
    for (int i = 0; i < heap->idx + 1; i++) {
      int _idx2 = sorted_idx[i];
      printf (" %d", _idx2);
    }
    printf ("\n\n");

  }
  
}


/*----------------------------------------------------------------------------
 * 
 * Add sub-triangles of a pn-triangle in the heap
 * 
 * parameters:
 *   heap             <-- heap to initialize
 *   order            <-- element order
 *   n_node           <-- number of nodes
 *   ho_vertex_num    <-- high order vertex num (internal ordering)
 *   vertex_coords    <-- vertex coordinates
 *   point_coords     <-- point to locate coordinates
 *
 *----------------------------------------------------------------------------*/

static void
_heap_fill_pn_sub_tria 
(
 _heap_t *heap,
 const int order,
 const int n_nodes,
 const double *nodes_coords,
 const double *point_coords
)
{
  int ibeg = 0;
  int iend = order;

  double *uvNodes   = malloc (sizeof(double) * 2 * n_nodes);

  _uv_ho_tria_nodes (order, 0., 1., 0, 1., uvNodes);
  
  int child = 0;
  for (int j = 0; j < order; j++) {
    int k1 = 0;
    for (int i = ibeg; i < iend - 1; i++) {

      int idx1 = i;
      int idx2 = i+1;
      int idx3 = iend + 1 + k1;
      int idx4 = iend + 2 + k1;
      
      double x1 = nodes_coords[3*idx1];
      double y1 = nodes_coords[3*idx1 + 1];
      double z1 = nodes_coords[3*idx1 + 2];

      double x2 = nodes_coords[3*idx2];
      double y2 = nodes_coords[3*idx2 + 1];
      double z2 = nodes_coords[3*idx2 + 2];

      double x3 = nodes_coords[3*idx3];
      double y3 = nodes_coords[3*idx3 + 1];
      double z3 = nodes_coords[3*idx3 + 2];

      double x4 = nodes_coords[3*idx4];
      double y4 = nodes_coords[3*idx4 + 1];
      double z4 = nodes_coords[3*idx4 + 2];
        
      double __vertex_coords[9] = {x1, y1, z1,
                                   x2, y2, z2,
                                   x3, y3, z3};
      double _closest_pointP1[3];
      double _uvClosestPointP1[2];
      double _uvClosestPointPn[2];
      double _weightsClosestPointP1[3];
      double _dist2;

      int isDegenerated = fvmc_triangle_evaluate_Position ((double *)point_coords,
                                                           __vertex_coords,
                                                           _closest_pointP1,
                                                           _uvClosestPointP1,
                                                           &_dist2,
                                                           _weightsClosestPointP1);

      double _uvPn_sub_tria[6];
      
      _uvPn_sub_tria[0] = uvNodes[2*idx1];
      _uvPn_sub_tria[1] = uvNodes[2*idx1+1];
      _uvPn_sub_tria[2] = uvNodes[2*idx2];
      _uvPn_sub_tria[3] = uvNodes[2*idx2+1];
      _uvPn_sub_tria[4] = uvNodes[2*idx3];
      _uvPn_sub_tria[5] = uvNodes[2*idx3+1];
      
      if (isDegenerated != -1) {
        
        for (int j1 = 0; j1 < 2; j1++) {
          _uvClosestPointPn[j1] = 0;
        }
        for (int j1 = 0; j1 < 2; j1++) {
          for (int k = 0; k < 3; k++) {
            _uvClosestPointPn[j1] += _weightsClosestPointP1[k] * _uvPn_sub_tria[2*k + j1];
          }
        }
        
        if (1 == 0) {
          printf("_uvClosestPointP1 : %12.5e %12.5e\n",
                 _uvClosestPointP1[0], _uvClosestPointP1[1]);
          printf("__vertex_coords + uvpn 1 : %12.5e %12.5e %12.5e // %12.5e %12.5e\n",
                 __vertex_coords[0], __vertex_coords[1], __vertex_coords[2],
                 _uvPn_sub_tria[0], _uvPn_sub_tria[1]);
          printf("__vertex_coords + uvpn 2 : %12.5e %12.5e %12.5e // %12.5e %12.5e\n",
                 __vertex_coords[3], __vertex_coords[4], __vertex_coords[5],
                 _uvPn_sub_tria[2], _uvPn_sub_tria[3]);
          printf("__vertex_coords + uvpn 3 : %12.5e %12.5e %12.5e // %12.5e %12.5e\n",
                 __vertex_coords[6], __vertex_coords[7], __vertex_coords[8],
                 _uvPn_sub_tria[4], _uvPn_sub_tria[5]);
        }
        
        _heap_insert (heap,
                      __vertex_coords,
                      _uvPn_sub_tria,
                      _closest_pointP1,
                      _uvClosestPointP1,
                      _uvClosestPointPn,
                      _dist2, child++);
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

      isDegenerated = fvmc_triangle_evaluate_Position ((double *) point_coords,
                                                       __vertex_coords,
                                                       _closest_pointP1,
                                                       _uvClosestPointP1,
                                                       &_dist2,
                                                       _weightsClosestPointP1);

      _uvPn_sub_tria[0] = uvNodes[2*idx2];
      _uvPn_sub_tria[1] = uvNodes[2*idx2+1];
      _uvPn_sub_tria[2] = uvNodes[2*idx4];
      _uvPn_sub_tria[3] = uvNodes[2*idx4+1];
      _uvPn_sub_tria[4] = uvNodes[2*idx3];
      _uvPn_sub_tria[5] = uvNodes[2*idx3+1];

      if (isDegenerated != -1) {
      
        for (int j1 = 0; j1 < 2; j1++) {
          _uvClosestPointPn[j1] = 0;
        }
        for (int j1 = 0; j1 < 2; j1++) {
          for (int k = 0; k < 3; k++) {
            _uvClosestPointPn[j1] += _weightsClosestPointP1[k] * _uvPn_sub_tria[2*k + j1];
          }
        }
      
        if (1 == 0) {
          printf("_uvClosestPointP1 : %12.5e %12.5e\n",
                 _uvClosestPointP1[0], _uvClosestPointP1[1]);
          printf("__vertex_coords 1 + uvpn  : %12.5e %12.5e %12.5e // %12.5e %12.5e\n",
                 __vertex_coords[0], __vertex_coords[1], __vertex_coords[2],
                 _uvPn_sub_tria[0], _uvPn_sub_tria[1]);
          printf("__vertex_coords 2 + uvpn : %12.5e %12.5e %12.5e // %12.5e %12.5e\n",
                 __vertex_coords[3], __vertex_coords[4], __vertex_coords[5],
                 _uvPn_sub_tria[2], _uvPn_sub_tria[3]);
          printf("__vertex_coords 3 + uvpn: %12.5e %12.5e %12.5e // %12.5e %12.5e\n",
                 __vertex_coords[6], __vertex_coords[7], __vertex_coords[8],
                 _uvPn_sub_tria[4], _uvPn_sub_tria[5]);
        }
      
        _heap_insert (heap,
                      __vertex_coords,
                      _uvPn_sub_tria,
                      _closest_pointP1,
                      _uvClosestPointP1,
                      _uvClosestPointPn,
                      _dist2, child++);
      }

      k1++;
    }

    int idx1 = iend - 1;
    int idx2 = iend - 1 + 1;
    int idx3 = iend + 1 + k1;

    double x1 = nodes_coords[3*idx1];
    double y1 = nodes_coords[3*idx1 + 1];
    double z1 = nodes_coords[3*idx1 + 2];
      
    double x2 = nodes_coords[3*idx2];
    double y2 = nodes_coords[3*idx2 + 1];
    double z2 = nodes_coords[3*idx2 + 2];
      
    double x3 = nodes_coords[3*idx3];
    double y3 = nodes_coords[3*idx3 + 1];
    double z3 = nodes_coords[3*idx3 + 2];
      
    double __vertex_coords[9] = {x1, y1, z1,
                                 x2, y2, z2,
                                 x3, y3, z3};
      
    double _closest_pointP1[3];
    double _uvClosestPointP1[2];
    double _uvClosestPointPn[2];
    
    double _weightsClosestPointP1[3];
   
    double _dist2;
    
    int isDegenerated = fvmc_triangle_evaluate_Position ((double *)point_coords,
                                                         __vertex_coords,
                                                         _closest_pointP1,
                                                         _uvClosestPointP1,
                                                         &_dist2,
                                                         _weightsClosestPointP1);
      
    double _uvPn_sub_tria[6];
    
    _uvPn_sub_tria[0] = uvNodes[2*idx1];
    _uvPn_sub_tria[1] = uvNodes[2*idx1+1];
    _uvPn_sub_tria[2] = uvNodes[2*idx2];
    _uvPn_sub_tria[3] = uvNodes[2*idx2+1];
    _uvPn_sub_tria[4] = uvNodes[2*idx3];
    _uvPn_sub_tria[5] = uvNodes[2*idx3+1];

    if (isDegenerated != -1) {
    
      for (int j1 = 0; j1 < 2; j1++) {
        _uvClosestPointPn[j1] = 0;
      }
      for (int j1 = 0; j1 < 2; j1++) {
        for (int k = 0; k < 3; k++) {
          _uvClosestPointPn[j1] += _weightsClosestPointP1[k] * _uvPn_sub_tria[2*k + j1];
        }
      }
      
      if (1 == 0) {
        printf("_uvClosestPointP1 : %12.5e %12.5e\n",
               _uvClosestPointP1[0], _uvClosestPointP1[1]);
        printf("__vertex_coords 1 + uvpn: %12.5e %12.5e %12.5e // %12.5e %12.5e\n",
               __vertex_coords[0], __vertex_coords[1], __vertex_coords[2],
               _uvPn_sub_tria[0], _uvPn_sub_tria[1]);
        printf("__vertex_coords 2 + uvpn: %12.5e %12.5e %12.5e // %12.5e %12.5e\n",
               __vertex_coords[3], __vertex_coords[4], __vertex_coords[5],
               _uvPn_sub_tria[2], _uvPn_sub_tria[3]);
        printf("__vertex_coords 3 + uvpn: %12.5e %12.5e %12.5e // %12.5e %12.5e\n",
               __vertex_coords[6], __vertex_coords[7], __vertex_coords[8],
               _uvPn_sub_tria[4], _uvPn_sub_tria[5]);
      }
    
      _heap_insert (heap,
                    __vertex_coords,
                    _uvPn_sub_tria,
                    _closest_pointP1,
                    _uvClosestPointP1,
                    _uvClosestPointPn,
                    _dist2, child++);
    }
      
    ibeg = iend + 1;
    iend += order - j;
  }

  free (uvNodes);
}


/*----------------------------------------------------------------------------
 * 
 * Add sub-triangles of a qn-quadrangle in the heap
 * 
 * parameters:
 *   heap             <-- heap to initialize
 *   order            <-- element order
 *   n_node           <-- number of nodes
 *   ho_vertex_num    <-- high order vertex num (internal ordering)
 *   vertex_coords    <-- vertex coordinates
 *   point_coords     <-- point to locate coordinates
 *
 *----------------------------------------------------------------------------*/

static void
_heap_fill_qn_sub_tria 
(
 _heap_t *heap,
 const int order,
 const int n_nodes,
 const double *nodes_coords,
 const double *point_coords
)
{
  double *uvNodes   = malloc (sizeof(double) * 2 * n_nodes);

  _uv_ho_quad_nodes (order, 0., 1., 0, 1., uvNodes);
  
  int child = 0;
  int step = order + 1;
  
  for (int j = 0; j < order; j++) {
    int j1 = j+1;
    
    for (int i = 0; i < order; i++) {

      int i1 = i+1;
      
      int idx1 = j*step  + i;
      int idx2 = j*step  + i1;
      int idx3 = j1*step + i;
      int idx4 = j1*step + i1;
      
      double x1 = nodes_coords[3*idx1];
      double y1 = nodes_coords[3*idx1 + 1];
      double z1 = nodes_coords[3*idx1 + 2];

      double x2 = nodes_coords[3*idx2];
      double y2 = nodes_coords[3*idx2 + 1];
      double z2 = nodes_coords[3*idx2 + 2];

      double x3 = nodes_coords[3*idx3];
      double y3 = nodes_coords[3*idx3 + 1];
      double z3 = nodes_coords[3*idx3 + 2];

      double x4 = nodes_coords[3*idx4];
      double y4 = nodes_coords[3*idx4 + 1];
      double z4 = nodes_coords[3*idx4 + 2];
        
      double __vertex_coords[9] = {x1, y1, z1,
                                   x2, y2, z2,
                                   x3, y3, z3};
      double _closest_pointP1[3];
      double _uvClosestPointP1[2];
      double _uvClosestPointPn[2];
      double _weightsClosestPointP1[3];
      double _dist2;

      int isDegenerated = fvmc_triangle_evaluate_Position ((double *)point_coords,
                                                           __vertex_coords,
                                                           _closest_pointP1,
                                                           _uvClosestPointP1,
                                                           &_dist2,
                                                           _weightsClosestPointP1);

      double _uvPn_sub_tria[6];
      
      _uvPn_sub_tria[0] = uvNodes[2*idx1];
      _uvPn_sub_tria[1] = uvNodes[2*idx1+1];
      _uvPn_sub_tria[2] = uvNodes[2*idx2];
      _uvPn_sub_tria[3] = uvNodes[2*idx2+1];
      _uvPn_sub_tria[4] = uvNodes[2*idx3];
      _uvPn_sub_tria[5] = uvNodes[2*idx3+1];
      
      if (isDegenerated != -1) {
        
        for (int j2 = 0; j2 < 2; j2++) {
          _uvClosestPointPn[j2] = 0;
        }
        for (int j2 = 0; j2 < 2; j2++) {
          for (int k = 0; k < 3; k++) {
            _uvClosestPointPn[j2] += _weightsClosestPointP1[k] * _uvPn_sub_tria[2*k + j2];
          }
        }
        
        if (1 == 0) {
          printf("_uvClosestPointP1 : %12.5e %12.5e\n",
                 _uvClosestPointP1[0],
                 _uvClosestPointP1[1]);
          printf("__vertex_coords + uvpn 1 : %12.5e %12.5e %12.5e // %12.5e %12.5e\n",
                 __vertex_coords[0],
                 __vertex_coords[1],
                 __vertex_coords[2],
                 _uvPn_sub_tria[0],
                 _uvPn_sub_tria[1]);
          printf("__vertex_coords + uvpn 2 : %12.5e %12.5e %12.5e // %12.5e %12.5e\n",
                 __vertex_coords[3],
                 __vertex_coords[4],
                 __vertex_coords[5],
                 _uvPn_sub_tria[2],
                 _uvPn_sub_tria[3]);
          printf("__vertex_coords + uvpn 3 : %12.5e %12.5e %12.5e // %12.5e %12.5e\n",
                 __vertex_coords[6],
                 __vertex_coords[7],
                 __vertex_coords[8],
                 _uvPn_sub_tria[4],
                 _uvPn_sub_tria[5]);
        }
        
        _heap_insert (heap,
                      __vertex_coords,
                      _uvPn_sub_tria,
                      _closest_pointP1,
                      _uvClosestPointP1,
                      _uvClosestPointPn,
                      _dist2, child++);
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

      isDegenerated = fvmc_triangle_evaluate_Position ((double *) point_coords,
                                                       __vertex_coords,
                                                       _closest_pointP1,
                                                       _uvClosestPointP1,
                                                       &_dist2,
                                                       _weightsClosestPointP1);

      _uvPn_sub_tria[0] = uvNodes[2*idx2];
      _uvPn_sub_tria[1] = uvNodes[2*idx2+1];
      _uvPn_sub_tria[2] = uvNodes[2*idx4];
      _uvPn_sub_tria[3] = uvNodes[2*idx4+1];
      _uvPn_sub_tria[4] = uvNodes[2*idx3];
      _uvPn_sub_tria[5] = uvNodes[2*idx3+1];

      if (isDegenerated != -1) {
      
        for (int j2 = 0; j2 < 2; j2++) {
          _uvClosestPointPn[j2] = 0;
        }
        for (int j2 = 0; j2 < 2; j2++) {
          for (int k = 0; k < 3; k++) {
            _uvClosestPointPn[j2] +=
              _weightsClosestPointP1[k] * _uvPn_sub_tria[2*k + j2];
          }
        }
      
        if (1 == 0) {
          printf("_uvClosestPointP1 : %12.5e %12.5e\n",
                 _uvClosestPointP1[0],
                 _uvClosestPointP1[1]);
          printf("__vertex_coords 1 + uvpn  : %12.5e %12.5e %12.5e // %12.5e %12.5e\n",
                 __vertex_coords[0],
                 __vertex_coords[1],
                 __vertex_coords[2],
                 _uvPn_sub_tria[0],
                 _uvPn_sub_tria[1]);
          printf("__vertex_coords 2 + uvpn : %12.5e %12.5e %12.5e // %12.5e %12.5e\n",
                 __vertex_coords[3],
                 __vertex_coords[4],
                 __vertex_coords[5],
                 _uvPn_sub_tria[2],
                 _uvPn_sub_tria[3]);
          printf("__vertex_coords 3 + uvpn: %12.5e %12.5e %12.5e // %12.5e %12.5e\n",
                 __vertex_coords[6],
                 __vertex_coords[7],
                 __vertex_coords[8],
                 _uvPn_sub_tria[4],
                 _uvPn_sub_tria[5]);
        }
      
        _heap_insert (heap,
                      __vertex_coords,
                      _uvPn_sub_tria,
                      _closest_pointP1,
                      _uvClosestPointP1,
                      _uvClosestPointPn,
                      _dist2, child++);
      }
    }
  }

  free (uvNodes);
}

/*-----------------------------------------------------------------------------
 * 
 * Add children to a heap
 * 
 * parameters:
 *   heap             <-> heap
 *   order            <-- element order
 *   n_nodes          <-- number of nodes
 *   nodes_coords     <-- nodes coordinates
 *   point_coords     <-- point to locate coordinates
 *   weightsPn        <-> work array
 *   vtx_tria_current <-- current triangle
 *   uvPn_tria_current<-- uv of current triangle vertices in the on element
 *   _basis_generic   <-- generic basis
 * 
 *----------------------------------------------------------------------------*/

static void
_insert_subtria
(
 _heap_t *heap,
 const int order,
 const int n_nodes,
 const double nodes_coords[],
 const double point_coords[],
 double weightsPn[],
 double vtx_tria_current[],
 double uvPn_tria_current[],
 _basis_generic_2D_t _basis_generic
 )
{
  double _vtx_tria_children[18];
  double _uvPn_tria_children[12];

  const int idx_sub_tria[12] = {0, 3, 5,
                                3, 4, 5,
                                3, 1, 4,
                                5, 4, 2};
  
  /* Compute middle vertices */
    
  for (int i = 0; i < 9; i++) {
    _vtx_tria_children[i] = vtx_tria_current[i];
  }
  
  for (int i = 0; i < 6; i++) {
    _uvPn_tria_children[i] = uvPn_tria_current[i];
  }
  
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 2; j++) {
      _uvPn_tria_children[6+2*i+j] =
        (uvPn_tria_current[2*i+j] + uvPn_tria_current[2*((i+1)%3)+j])/2;
    }
    
    (_basis_generic) (order   ,
                      1,
                      _uvPn_tria_children + 6 + 2*i,
                      weightsPn);
    
    for (int j = 0; j < 3; j++) {
      _vtx_tria_children[9+3*i+j] = 0;
    }
    for (int k = 0; k < n_nodes; k++) {
      const double *_node_coords = nodes_coords + 3 * k;
      for (int j = 0; j < 3; j++) {
        _vtx_tria_children[9+3*i+j] += weightsPn[k] * _node_coords[j];
      }
    }
  }
  
  int child = 0;
  for (int i = 0; i < 4; i++) {

    double _vtx_tria_child[9];
    double _uvPn_tria_child[6];
    
    for (int j = 0; j < 3; j++) {
      int _j = idx_sub_tria[3 * i + j];
      for (int k = 0; k < 3; k++) {
        _vtx_tria_child[3*j + k] = _vtx_tria_children[3*_j + k];
        }
      for (int k = 0; k < 2; k++) {
        _uvPn_tria_child[2*j + k] = _uvPn_tria_children[2*_j + k];
      }
    }
    
    double _closest_pt_child[3];
    double _closest_pt_uvP1_child[2];
    double _closest_pt_uvPn_child[2];
    double _dist2_child = 0;
    double _closest_pt_weights_child[3];
    
    int isDegenerated = fvmc_triangle_evaluate_Position ((double *)point_coords,
                                                         _vtx_tria_child,
                                                         _closest_pt_child,
                                                         _closest_pt_uvP1_child,
                                                         &_dist2_child,
                                                         _closest_pt_weights_child);
    
    
    if (isDegenerated == -1) {
      continue;
    }
      
    for (int j = 0; j < 2; j++) {
      _closest_pt_uvPn_child[j] = 0;
    }
    for (int j = 0; j < 2; j++) {
      for (int k = 0; k < 3; k++) {
        _closest_pt_uvPn_child[j] +=
          _closest_pt_weights_child[k] * _uvPn_tria_child[2*k + j];
      }
    }
    
    if (1 == 0) {
      printf("isDegenrated : %d\n", isDegenerated);
      printf("_closest_pt_weights_child : %12.5e %12.5e %12.5e\n"
             ,_closest_pt_weights_child[0]
             ,_closest_pt_weights_child[1]
             ,_closest_pt_weights_child[2]);
      
      printf("_closest_pt_child : %12.5e %12.5e %12.5e\n"
             ,_closest_pt_child[0]
             ,_closest_pt_child[1]
             ,_closest_pt_child[2]);
      
      printf("_closest_pt_uvP1_child : %12.5e %12.5e\n"
             ,_closest_pt_uvP1_child[0]
             ,_closest_pt_uvP1_child[1]);
      
      
    }
    
    _heap_insert (heap,
                  _vtx_tria_child,
                  _uvPn_tria_child,
                  _closest_pt_child,
                  _closest_pt_uvP1_child,
                  _closest_pt_uvPn_child,
                  _dist2_child, child++);
    
  }

}  


/*-----------------------------------------------------------------------------
 * 
 * compute distance from closest triangle subdivision
 * 
 * parameters:
 *   heap             <-> heap
 *   order            <-- element order
 *   n_nodes           <-- number of nodes
 *   n_it_max         <-- maximum of iteration to compute distance
 *   err_max          <-- maximum error of the projected point
 *   err_proj         --> projected error 
 *   nodes_coords     <-- nodes coordinates
 *   point_coords     <-- point to locate coordinates
 *   weightsPn        <-> work array
 *   projected_coords --> current triangle
 *   uvw              --> uvw
 *   n_it             --> number of iterations
 *   err_proj         --> error of the projected point
 *   uncertain_result --> 1 if the result is uncertain
 *   _basis_generic   <-- generic basis
 * 
 *----------------------------------------------------------------------------*/

static double
_compute_dist2_from_closest_tria_subdivision
(
 _heap_t *heap,
 const int order,
 const int n_nodes,
 const int n_it_max,
 const double err_max,
 const double nodes_coords[],
 const double point_coords[],
 double weightsPn[],
 double projected_coords[],
 double uvw[],
 int    *n_it,
 double *err_proj,
 int *uncertain_result,
 _basis_generic_2D_t _basis_generic
 )
{
  *uncertain_result = 0;
  *n_it = 0;
  *err_proj = HUGE_VAL;
  double dist2 = HUGE_VAL;

  double dist2_min_min = HUGE_VAL;
  double dist2_pre = HUGE_VAL;
  int distance_extension = 0;

  while (1) {

    double _vtx_tria_current[9];
    double _uvPn_tria_current[6];
    
    double _closest_pt_current[3];
    double _closest_pt_uvP1_current[2];
    double _closest_pt_uvPn_current[2];
    double _dist2_current;
    
    /* Get closest triangle stored in the heap */

    int _child;
    int isEmpty = _heap_top_get (heap,
                                 _vtx_tria_current,
                                 _uvPn_tria_current,
                                 _closest_pt_current,
                                 _closest_pt_uvP1_current,
                                 _closest_pt_uvPn_current,
                                 &_dist2_current, &_child);

      
    if (isEmpty) {
      bftc_error(__FILE__, __LINE__, 0,
                 _("Heap is empty %s\n"));
      abort();
    }

    if ((distance_extension == 0) && (_dist2_current > dist2_pre)) {
      distance_extension = 1;
    }

    else if (distance_extension == 1) {
      if (_dist2_current <= dist2_min_min) {
        distance_extension = 0;
      }
    }

    dist2_min_min = FVMC_MIN (dist2_min_min, _dist2_current);
    dist2_pre = _dist2_current;

    /* Compute projected from current P1 triangle */
    
    double _projected_coords_from_p1[3];

    for (int j = 0; j < 3; j++) {
      _projected_coords_from_p1[j] = 0;
    }

    double weightsP1[3];
    
    _basis_tria_pn (1,
                    1,
                    _closest_pt_uvP1_current,
                    weightsP1);

    if (0 == 1) {
      printf("\n\n ========= get heap =========\n");
    
      printf ("uv Pn : %22.15e %22.15e\n",
              _closest_pt_uvPn_current[0], _closest_pt_uvPn_current[1]);
      printf ("uv P1 : %22.15e %22.15e\n",
              _closest_pt_uvP1_current[0], _closest_pt_uvP1_current[1]);
      printf ("_dist2 child : %22.15e %d\n", _dist2_current, _child);
    
      printf("Weights : %12.5e %12.5e %12.5e\n",
             weightsP1[0], weightsP1[1], weightsP1[2]);
      printf("vtx_tria_current 1 : %12.5e %12.5e %12.5e\n",
             _vtx_tria_current[0], _vtx_tria_current[1], _vtx_tria_current[2]);
      printf("vtx_tria_current 2 : %12.5e %12.5e %12.5e\n",
             _vtx_tria_current[3], _vtx_tria_current[4], _vtx_tria_current[5]);
      printf("vtx_tria_current 3 : %12.5e %12.5e %12.5e\n",
             _vtx_tria_current[6], _vtx_tria_current[7], _vtx_tria_current[8]);
    }

    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        _projected_coords_from_p1[k] += weightsP1[j] * _vtx_tria_current[3*j+k]; 
      }
    }

    /* Compute projected from current Pn triangle */

    double _projected_coords_from_pn[3];

    for (int j = 0; j < 3; j++) {
      _projected_coords_from_pn[j] = 0;
    }

    (_basis_generic) (order,
                      1,
                      _closest_pt_uvPn_current,
                      weightsPn);

    for (int j = 0; j < n_nodes; j++) {
      const double *_node_coords = nodes_coords + 3 * j;
      for (int k = 0; k < 3; k++) {
        _projected_coords_from_pn[k] += weightsPn[j] * _node_coords[k]; 
      }
    }

    /* Compute distance between two projected */    

    *err_proj = 0;
    for (int i = 0; i < 3; i++) {
      double val = _projected_coords_from_pn[i] - _projected_coords_from_p1[i];
      *err_proj += val * val;
    }

    /* Break if error is ok */

    if (*err_proj <= err_max || (*n_it)++ >= n_it_max) {
       
      for (int j = 0; j < 3; j++) {
        projected_coords[j] = _projected_coords_from_pn[j];
      }

      dist2 = 0;
      for (int j = 0; j < 3; j++) {
        double comp = projected_coords[j] - point_coords[j];
        dist2 += comp * comp; 
      }

      uvw[0] = _closest_pt_uvPn_current[0];
      uvw[1] = _closest_pt_uvPn_current[1];
        
      break;
    }
    
    /* 
     * Insert sub-triangles in the heap 
     */

    _insert_subtria (heap,
                     order,
                     n_nodes,
                     nodes_coords,
                     point_coords,
                     weightsPn,
                     _vtx_tria_current,
                     _uvPn_tria_current,
                     _basis_generic);
    
  }

  if (*n_it >= n_it_max) {

    if (isWarningPrinted2 == 0) {
      isWarningPrinted2 = 1;

      bftc_printf("warning _compute_dist2_from_closest_tria_subdivision : "
                  "compute of projected point is not converged (error = %17.10e > error max %17.10e\n",
                *err_proj, err_max);
    }
  }

  if (distance_extension) {
    *uncertain_result = 1;
  }
  
  return dist2;
}



/*-----------------------------------------------------------------------------
 * 
 * compute distance from closest triangle subdivision
 * 
 * parameters:
 *   heap1            <-> heap
 *   heap2            <-> work heap
 *   order            <-- element order
 *   n_nodes           <-- number of nodes
 *   n_it_max         <-- maximum of iteration to compute distance
 *   err_max          <-- maximum error of the projected point
 *   err_proj         --> projected error 
 *   ho_vertex_num    <-- high order vertex num (internal ordering)
 *   nodes_coords     <-- nodes coordinates
 *   point_coords     <-- point to locate coordinates
 *   weightsPn        <-> work array
 *   projected_coords --> current triangle
 *   uvw              --> uvw
 *   n_it             --> number of iterations
 *   err_proj         --> error of the projected point
 *   _basis_generic   <-- generic basis
 * 
 *----------------------------------------------------------------------------*/

static double
_compute_dist2_from_uniform_tria_subdivision
(
 _heap_t *heap1,
 _heap_t *heap2,
 const int order,
 const int n_nodes,
 const int n_it_max,
 const double err_max,
 const double nodes_coords[],
 const double point_coords[],
 double weightsPn[],
 double projected_coords[],
 double uvw[],
 int    *n_it,
 double *err_proj,
 _basis_generic_2D_t _basis_generic
 )
{
  *n_it = 0;
  *err_proj = HUGE_VAL;
  double dist2 = HUGE_VAL;

  double dist2_min_min = HUGE_VAL;

  _heap_t *heap      = heap1;
  _heap_t *next_heap = heap2;

  while (1) {

    double _vtx_tria_current[9];
    double _uvPn_tria_current[6];
    
    double _closest_pt_current[3];
    double _closest_pt_uvP1_current[2];
    double _closest_pt_uvPn_current[2];
    double _dist2_current;
    
    /* Get closest triangle stored in the heap */

    int _child;
    int isEmpty = _heap_top_get (heap,
                                 _vtx_tria_current,
                                 _uvPn_tria_current,
                                 _closest_pt_current,
                                 _closest_pt_uvP1_current,
                                 _closest_pt_uvPn_current,
                                 &_dist2_current, &_child);

      
    if (isEmpty) {
      bftc_error(__FILE__, __LINE__, 0,
                 _("Heap is empty %s\n"));
      abort();
    }

    dist2_min_min = FVMC_MIN (dist2_min_min, _dist2_current);

    /* Compute projected from current P1 triangle */
    
    double _projected_coords_from_p1[3];

    for (int j = 0; j < 3; j++) {
      _projected_coords_from_p1[j] = 0;
    }

    double weightsP1[3];
    
    _basis_tria_pn (1,
                    1,
                    _closest_pt_uvP1_current,
                    weightsP1);

    if (0 == 1) {
      printf("\n\n ========= get heap =========\n");
    
      printf ("uv Pn : %22.15e %22.15e\n", _closest_pt_uvPn_current[0],
              _closest_pt_uvPn_current[1]);
      printf ("uv P1 : %22.15e %22.15e\n", _closest_pt_uvP1_current[0],
              _closest_pt_uvP1_current[1]);
      printf ("_dist2 child : %22.15e %d\n", _dist2_current, _child);
    
      printf("Weights : %12.5e %12.5e %12.5e\n",
             weightsP1[0], weightsP1[1], weightsP1[2]);
      printf("vtx_tria_current 1 : %12.5e %12.5e %12.5e\n",
             _vtx_tria_current[0], _vtx_tria_current[1], _vtx_tria_current[2]);
      printf("vtx_tria_current 2 : %12.5e %12.5e %12.5e\n",
             _vtx_tria_current[3], _vtx_tria_current[4], _vtx_tria_current[5]);
      printf("vtx_tria_current 3 : %12.5e %12.5e %12.5e\n",
             _vtx_tria_current[6], _vtx_tria_current[7], _vtx_tria_current[8]);
    }

    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        _projected_coords_from_p1[k] += weightsP1[j] * _vtx_tria_current[3*j+k]; 
      }
    }

    /* Compute projected from current Pn triangle */

    double _projected_coords_from_pn[3];

    for (int j = 0; j < 3; j++) {
      _projected_coords_from_pn[j] = 0;
    }

    (_basis_generic) (order,
                      1,
                     _closest_pt_uvPn_current,
                     weightsPn);
    
    for (int j = 0; j < n_nodes; j++) {
      const double *_node_coords = nodes_coords + 3 * j;
      for (int k = 0; k < 3; k++) {
        _projected_coords_from_pn[k] += weightsPn[j] * _node_coords[k]; 
      }
    }

    /* Compute distance between two projected */    

    *err_proj = 0;
    for (int i = 0; i < 3; i++) {
      double val = _projected_coords_from_pn[i] - _projected_coords_from_p1[i];
      *err_proj += val * val;
    }

    /* Break if error is ok */

    if (*err_proj <= err_max || (*n_it)++ >= n_it_max) {
       
      for (int j = 0; j < 3; j++) {
        projected_coords[j] = _projected_coords_from_pn[j];
      }

      dist2 = 0;
      for (int j = 0; j < 3; j++) {
        double comp = projected_coords[j] - point_coords[j];
        dist2 += comp * comp; 
      }

      uvw[0] = _closest_pt_uvPn_current[0];
      uvw[1] = _closest_pt_uvPn_current[1];
        
      break;
    }
    
    /* 
     * Insert sub-triangles in the next heap
     */
    
    _heap_init (next_heap);

    _insert_subtria (next_heap,
                     order,
                     n_nodes,
                     nodes_coords,
                     point_coords,
                     weightsPn,
                     _vtx_tria_current,
                     _uvPn_tria_current,
                     _basis_generic);

    double _vtx_tria_current2[9];
    double _uvPn_tria_current2[6];
    
    double _closest_pt_current2[3];
    double _closest_pt_uvP1_current2[2];
    double _dist2_current2;
    int _child_current2;
    
    while ( !_heap_top_get (heap,
                            _vtx_tria_current2,
                            _uvPn_tria_current2,
                            _closest_pt_current2,
                            _closest_pt_uvP1_current2,
                            _closest_pt_uvPn_current,
                            &_dist2_current2, &_child_current2)) {


      _insert_subtria (next_heap,
                       order,
                       n_nodes,
                       nodes_coords,
                       point_coords,
                       weightsPn,
                       _vtx_tria_current2,
                       _uvPn_tria_current2,
                       _basis_generic);

    }

    _heap_t *heap_tmp = heap;
    heap = next_heap;
    next_heap = heap_tmp;
    
  }

  if (*n_it >= n_it_max) {

    if (isWarningPrinted1 == 0) {
          isWarningPrinted1 = 1;

      bftc_printf("warning _compute_dist2_from_uniform_tria_subdivision : "
                  "compute of projected point is not converged (error = %17.10e > error max %17.10e\n",
                  *err_proj, err_max);
    }
  }

  return dist2;
}


/*-----------------------------------------------------------------------------
 * 
 * Default point location on a high order triangle
 * 
 * parameters:
 *   order            <-- element order
 *   n_nodes          <-- number of nodes
 *   nodes_coords     <-- nodes coordinates
 *   point_coords     <-- point to locate coordinates
 *   projected_coords --> projected point coordinates (or NULL)
 *   uvw              --> parametric coordinates in the element
 * 
 * return: 
 *   distance to the cell
 *
 *----------------------------------------------------------------------------*/

static double
_default_location_generic_2d
(
 const int order,
 const double char_size,
 const int n_nodes,
 const double *nodes_coords,
 const double *point_coords,
 double *projected_coords,
 double *uvw,
 _heap_fill_init_sub_tria_t fill_init_fct,
 _basis_generic_2D_t basis_generic
)
{

  if (1 == 0) {
    printf ("\n\n***********\n");
  
    printf ("n_node %d\n", n_nodes);
  }
 
  //  const int n_it_max = 100;
  const int n_it_max = 100;
  double err_max = FVMC_MAX (char_size * 1e-6, 1e-15);
  
  double dist2 = HUGE_VAL;

  _heap_t heap;
  _heap_t heap2;

  double *weightsPn = malloc(sizeof(double) * n_nodes);

  /* Initialize heap */
  
  _heap_init (&heap);
  
  /* Build initial sub-triangles and store them in the heap */

  (fill_init_fct) (&heap,
                   order,
                   n_nodes,
                   nodes_coords,
                   point_coords);

  /* 
   *  While error > error_max
   *    - Get closest triangle in the heap
   *    - Cut it in sub-triangles
   *    - Store them in the heap
   */


  const int method = 0;

  int n_it;
  double err_proj = HUGE_VAL;
  int uncertain_result = 0;;

  if (_idebug == 1) {
    printf ("=====================debut=========================\n");
  }
  
  if (method == 0) {
    dist2 = _compute_dist2_from_closest_tria_subdivision (&heap,
                                                          order,
                                                          n_nodes,
                                                          n_it_max,
                                                          err_max,
                                                          nodes_coords,
                                                          point_coords,
                                                          weightsPn,
                                                          projected_coords,
                                                          uvw,
                                                          &n_it,
                                                          &err_proj,
                                                          &uncertain_result,
                                                          basis_generic);

    if (1 == 0) {
      printf("\nCalcul distance triangle premier essai :\n");
      printf("          point_coords : %22.15e %22.15e %22.15e         \n", point_coords[0]    , point_coords[1]    , point_coords[2]);
      printf("  project point_coords : %22.15e %22.15e %22.15e - %22.15e\n", projected_coords[0], projected_coords[1], projected_coords[2], dist2);
      printf("  iteration number, error^2 : %d %22.15e\n", n_it, err_proj);
    }
    
    if (uncertain_result) {
      
      /* Initialize heap */
  
      _heap_init (&heap);
      _heap_init (&heap2);
  
      /* Build initial sub-triangles and store them in the heap */
  
      (fill_init_fct) (&heap,
                       order,
                       n_nodes,
                       nodes_coords,
                       point_coords);


      dist2 = _compute_dist2_from_uniform_tria_subdivision (&heap,
                                                            &heap2,
                                                            order,
                                                            n_nodes,
                                                            n_it_max,
                                                            err_max,
                                                            nodes_coords,
                                                            point_coords,
                                                            weightsPn,
                                                            projected_coords,
                                                            uvw,
                                                            &n_it,
                                                            &err_proj,
                                                            basis_generic);

      if (1 == 0) {
        printf("\nCalcul distance triangle deuxieme essai :\n");
        printf("          point_coords : %22.15e %22.15e %22.15e         \n", point_coords[0]    , point_coords[1]    , point_coords[2]);
        printf("  project point_coords : %22.15e %22.15e %22.15e - %22.15e\n", projected_coords[0], projected_coords[1], projected_coords[2], dist2);
        printf("  iteration number, error^2 : %d %22.15e\n", n_it, err_proj);
      }
      
    }
      if (_idebug == 1) {

        printf ("====================fin============================\n\n\n");
      }
  }
  
  else {

    _heap_init (&heap2);
    dist2 = _compute_dist2_from_uniform_tria_subdivision (&heap,
                                                          &heap2,
                                                          order,
                                                          n_nodes,
                                                          n_it_max,
                                                          err_max,
                                                          nodes_coords,
                                                          point_coords,
                                                          weightsPn,
                                                          projected_coords,
                                                          uvw,
                                                          &n_it,
                                                          &err_proj,
                                                          basis_generic);

  }

  free (weightsPn);

  //  printf("\nCalcul distance triangle :\n");
  // printf("          point_coords : %22.15e %22.15e %22.15e         \n", point_coords[0]    , point_coords[1]    , point_coords[2]);
  //printf("  project point_coords : %22.15e %22.15e %22.15e - %22.15e\n", projected_coords[0], projected_coords[1], projected_coords[2], dist2);
  //printf("  iteration number, error^2 : %d %22.15e\n", n_it, err_proj);

  return dist2;
  
}


/*----------------------------------------------------------------------------
 * 
 * Compute the radius of a triangle inscribed circle
 * 
 * parameters:
 *   coords           <-- coordinates of vertices
 * 
 * return: 
 *   radius
 *
 *----------------------------------------------------------------------------*/

static double 
_radius_inscribed_circle
(
 const double *coords
)
{
  double a = sqrt ((coords[3*1    ] - coords[3*0    ]) * (coords[3*1    ] - coords[3*0    ]) +
                   (coords[3*1 + 1] - coords[3*0 + 1]) * (coords[3*1 + 1] - coords[3*0 + 1]) +
                   (coords[3*1 + 2] - coords[3*0 + 2]) * (coords[3*1 + 2] - coords[3*0 + 2]));

  double b = sqrt ((coords[3*2    ] - coords[3*1    ]) * (coords[3*2    ] - coords[3*1    ]) +
                   (coords[3*2 + 1] - coords[3*1 + 1]) * (coords[3*2 + 1] - coords[3*1 + 1]) +
                   (coords[3*2 + 2] - coords[3*1 + 2]) * (coords[3*2 + 2] - coords[3*1 + 2]));

  double c = sqrt ((coords[3*0    ] - coords[3*2    ]) * (coords[3*0    ] - coords[3*2    ]) +
                   (coords[3*0 + 1] - coords[3*2 + 1]) * (coords[3*0 + 1] - coords[3*2 + 1]) +
                   (coords[3*0 + 2] - coords[3*2 + 2]) * (coords[3*0 + 2] - coords[3*2 + 2]));

  double p = a + b + c;
  double S = sqrt (p*(p-a)*(p-b)*(p-c));
  
  return S/p;
}


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

static double 
_default_location
(
 const fvmc_element_t type,
 const int order,
 const int n_nodes,
 const double *nodes_coords,
 const double *point_coords,
 double *projected_coords,
 double* uvw
)
{
  fvmc_element_t _type = type;
  double dist2 = HUGE_VAL;

  
  switch (_type) {

  case FVMC_FACE_TRIA: {

    int v1 = 0;
    int v2 = order;
    int v3 = (order + 2) * (order + 1) / 2 - 1;
      
    double p1_coords[9] = {nodes_coords[3*v1], nodes_coords[3*v1+1], nodes_coords[3*v1+2],
                           nodes_coords[3*v2], nodes_coords[3*v2+1], nodes_coords[3*v2+2],
                           nodes_coords[3*v3], nodes_coords[3*v3+1], nodes_coords[3*v3+2]};
    
    double char_size = _radius_inscribed_circle (p1_coords);
    
    dist2 = _default_location_generic_2d (order,
                                          char_size,
                                          n_nodes,
                                          nodes_coords,
                                          point_coords,
                                          projected_coords,
                                          uvw,
                                          _heap_fill_pn_sub_tria,
                                          _basis_tria_pn);

    break;

  }

  case FVMC_FACE_QUAD: {

    const int order1 = order + 1;
    
    const int quadrangle_vertices[4] = {1,
                                        order1,
                                        order1 * order + 1,
                                        order1 * order1};
    int triangle_vertices[6];
    
    int n_sub_tria =  fvmc_triangulate_quadrangle(3,
                                                  nodes_coords,
                                                  NULL,
                                                  quadrangle_vertices,
                                                  triangle_vertices);

    double char_size = HUGE_VAL;
    
    for (int i = 0; i < n_sub_tria; i++) {

      int v1 = triangle_vertices[3*i    ] - 1;
      int v2 = triangle_vertices[3*i + 1] - 1;
      int v3 = triangle_vertices[3*i + 2] - 1;

      double p1_coords[9] = {nodes_coords[3*v1], nodes_coords[3*v1+1], nodes_coords[3*v1+2],
                             nodes_coords[3*v2], nodes_coords[3*v2+1], nodes_coords[3*v2+2],
                             nodes_coords[3*v3], nodes_coords[3*v3+1], nodes_coords[3*v3+2]};
    
      double _char_size = _radius_inscribed_circle (p1_coords);
     
      char_size = FVMC_MAX (char_size, _char_size);
      
    }
    
    dist2 = _default_location_generic_2d (order,
                                          char_size,
                                          n_nodes,
                                          nodes_coords,
                                          point_coords,
                                          projected_coords,
                                          uvw,
                                          _heap_fill_qn_sub_tria,
                                          _basis_quad_qn);
    break;
  }
    
  default:
    bftc_error(__FILE__, __LINE__, 0,
               _("_default_location : Element not implemented yet\n"));

  }
  
  return dist2;

}


/*----------------------------------------------------------------------------
 * 
 * high order basis
 * 
 * parameters:
 *   type            <-- element type
 *
 * return:
 *
 *----------------------------------------------------------------------------*/

static fvmc_ho_user_elt_t *
_get_user_elt (fvmc_element_t elt_type)
{

  fvmc_ho_user_elt_t *user_elt = NULL;
  
  switch(elt_type) {

  case FVMC_EDGE:
    user_elt = _user_edge;
    break;
    
  case FVMC_FACE_TRIA:
    user_elt = _user_tria;
    break;
    
  case FVMC_FACE_QUAD:
    user_elt = _user_quad;
    break;
    
  case FVMC_CELL_TETRA:
    user_elt = _user_tetra;
    break;
    
  case FVMC_CELL_PYRAM:
    user_elt = _user_pyra;
    break;
    
  case FVMC_CELL_PRISM:
    user_elt = _user_prism;
    break;
    
  case FVMC_CELL_HEXA:
    user_elt = _user_hexa;
    break;

  default:
    bftc_error(__FILE__, __LINE__, 0,
               _("fvmc_ho_user_elt_unset : Unvailable element type\n"));
  }

  return user_elt;
}

/*----------------------------------------------------------------------------
 * 
 * default high order basis
 * 
 * parameters:
 *   type            <-- element type
 *   order           <-- order
 *   n_pts           <-- number of points 
 *   uvw             <-- uvw (size = elt_dim * n_pts)
 *   weights         --> weights (size = n_nodes * n_pts)
 *
 *----------------------------------------------------------------------------*/

static void
_default_elt_basis
(
const fvmc_element_t type,
const int order,
const int n_pts,
const double *uvw,
      double *weights 
)
{
  switch (type) {

  case FVMC_FACE_TRIA:
    _basis_tria_pn (order, n_pts, uvw, weights);
    break;

  case FVMC_FACE_QUAD:
    _basis_quad_qn (order, n_pts, uvw, weights);
    break;

  case FVMC_EDGE:
  case FVMC_CELL_TETRA:
  case FVMC_CELL_PRISM:
  case FVMC_CELL_PYRAM:
  case FVMC_CELL_HEXA:
  default: 
    bftc_error(__FILE__, __LINE__, 0,
               _("_default_elts_basis : '%d' element type not yet implemented\n"),
               type);
  }
}

/*============================================================================
 * Public function definitions
 *============================================================================*/


/*----------------------------------------------------------------------------
 * 
 * Unset a user element
 * 
 *----------------------------------------------------------------------------*/

void
fvmc_ho_user_elt_unset (fvmc_element_t elt_type)
{

  fvmc_ho_user_elt_t *_user_elt = _get_user_elt (elt_type);
  
  if (_user_elt != NULL) {
    free (_user_elt);
    _user_elt = NULL;
  }
  
}


void
fvmc_ho_user_elts_unset (void)
{
  fvmc_ho_user_elt_unset (FVMC_EDGE);
  fvmc_ho_user_elt_unset (FVMC_FACE_TRIA);
  fvmc_ho_user_elt_unset (FVMC_FACE_QUAD);
  fvmc_ho_user_elt_unset (FVMC_CELL_TETRA);
  fvmc_ho_user_elt_unset (FVMC_CELL_HEXA);
  fvmc_ho_user_elt_unset (FVMC_CELL_PRISM);
  fvmc_ho_user_elt_unset (FVMC_CELL_PYRAM);
}

/*----------------------------------------------------------------------------
 * 
 * Unset a user element
 * 
 *----------------------------------------------------------------------------*/


void
fvmc_ho_user_elt_set (fvmc_element_t elt_type,
                      fvmc_ho_basis_fct_t elt_basis,
                      fvmc_ho_xsi_fct_t xsi_coords,
                      fvmc_ho_location_fct_t location_in_elt)
{
  fvmc_ho_user_elt_t *user_elt = _get_user_elt (elt_type);
  
  if (user_elt == NULL) {
    user_elt = (fvmc_ho_user_elt_t *) malloc (sizeof(fvmc_ho_user_elt_t));
  }

  user_elt->elt_basis = elt_basis;
  user_elt->xsi_coords = xsi_coords;
  user_elt->location_in_elt = location_in_elt;

}

/*----------------------------------------------------------------------------
 * 
 * high order basis
 * 
 * parameters:
 *   type            <-- element type
 *   order           <-- order
 *   n_nodes         <-- number of nodes
 *   n_pts           <-- number of points 
 *   uvw             <-- uvw
 *   uvw             <-- uvw (size = n_pts)
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
)
{
  fvmc_ho_user_elt_t *user_elt = _get_user_elt (type);
  
  if (user_elt != NULL) {
    if (user_elt->elt_basis != NULL) {
      (user_elt->elt_basis) (order,
                             n_nodes,
                             n_pts,
                             uvw,
                             weights);
    }
    else {
      _default_elt_basis (type,
                          order,
                          n_pts,
                          uvw,
                          weights);
    }
  }

  else {
    _default_elt_basis (type,
                        order,
                        n_pts,
                        uvw,
                        weights);

  }
}


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
)
{
  
  fvmc_ho_user_elt_t *user_elt = _get_user_elt (type);
  
  if (user_elt != NULL) {
    if (user_elt->location_in_elt != NULL) {
      return (user_elt->location_in_elt ) (order,
                                           n_nodes,
                                           nodes_coords,
                                           point_coords,
                                           projected_coords,
                                           uvw);
    }
    else {
      return _default_location (type,
                                order,
                                n_nodes,
                                nodes_coords,
                                point_coords,
                                projected_coords,
                                uvw);
    }
  }

  else {

    return _default_location (type,
                              order,
                              n_nodes,
                              nodes_coords,
                              point_coords,
                              projected_coords,
                              uvw);
  }

  return HUGE_VAL;
}


/*----------------------------------------------------------------------------
 * 
 * Free static variables
 * 
 *----------------------------------------------------------------------------*/


void
fvmc_ho_free
(
 void
)
{
#if defined (HAVE_SPACE_BASIS) 
  if (_ijk_tria_space != NULL) {

    for (int i = 0; i < _n_ijk_tria_space; i++) {
      if (_ijk_tria_space[i] != NULL) {
        free(_ijk_tria_space[i]);
        _ijk_tria_space[i] = NULL;
      }
    }
    free (_ijk_tria_space);
    _ijk_tria_space = NULL;
  
  }

  if (_ijk_quad_space != NULL) {

    for (int i = 0; i < _n_ijk_quad_space; i++) {
      if (_ijk_quad_space[i] != NULL) {
        free(_ijk_quad_space[i]);
        _ijk_quad_space[i] = NULL;
      }
    }
    free (_ijk_quad_space);
    _ijk_quad_space = NULL;
    
  }
  
  if (_ijk_1D_space != NULL) {

    for (int i = 0; i < _n_ijk_1D_space; i++) {
      if (_ijk_1D_space[i] != NULL) {
        free(_ijk_1D_space[i]);
        _ijk_1D_space[i] = NULL;
      }
    }
    
    free (_ijk_1D_space);
    _ijk_1D_space = NULL;
  }

#endif
  
}

#ifdef __cplusplus
}
#endif /* __cplusplus */
