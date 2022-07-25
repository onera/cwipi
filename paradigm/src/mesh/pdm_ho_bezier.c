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
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm_error.h"
#include "pdm_printf.h"
#include "pdm_logging.h"
#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_mpi.h"
#include "pdm_vtk.h"
#include "pdm_ho_location.h"
#include "pdm_mesh_nodal.h"
#include "pdm_triangle.h"
#include "pdm_ho_bezier_basis.h"
#include "pdm_ho_bezier.h"

/*=============================================================================
 * Macro definition
 *============================================================================*/

#define ij2idx(i, j, n) ((i) + (j)*((n) + 1) - ((j)-1)*(j)/2)

/*============================================================================
 * Private function definitions
 *============================================================================*/

static int
_newton_tria
(
 const int     order,
 const double *target,
       double *u,
       double *v,
       double *p,
       double *dp_du,
       double *dp_dv,
       double *xyz
)
{
  const int vb = 1;
  if (vb) {
    log_trace(">>> Newton\n");
  }

  const int    it_max  = 10;
  const double tol_res = 1e-6;
  const double tol_uv2 = 1e-12;

  double _u = *u;
  double _v = *v;

  double dxyz_du[3], dxyz_dv[3];

  int converged = 0;

  for (int it = 0; it < it_max; it++) {
    PDM_ho_bezier_de_casteljau_tria(3,
                                    order,
                                    _u,
                                    _v,
                                    p,
                                    xyz,
                                    NULL, NULL, NULL);

    double vec[3] = {
      target[0] - xyz[0],
      target[1] - xyz[1],
      target[2] - xyz[2]
    };

    // Jacobian
    PDM_ho_bezier_de_casteljau_tria(3,
                                    order-1,
                                    _u,
                                    _v,
                                    dp_du,
                                    dxyz_du,
                                    NULL, NULL, NULL);

    PDM_ho_bezier_de_casteljau_tria(3,
                                    order-1,
                                    _u,
                                    _v,
                                    dp_dv,
                                    dxyz_dv,
                                    NULL, NULL, NULL);

    double a00 = PDM_DOT_PRODUCT(dxyz_du, dxyz_du);
    double a01 = PDM_DOT_PRODUCT(dxyz_du, dxyz_dv);
    double a11 = PDM_DOT_PRODUCT(dxyz_dv, dxyz_dv);

    double b0 = PDM_DOT_PRODUCT(dxyz_du, vec);
    double b1 = PDM_DOT_PRODUCT(dxyz_dv, vec);

    double res = PDM_MAX(PDM_ABS(b0)/a00, PDM_ABS(b1)/a11);

    if (vb) {
      log_trace("  it %d, u = %f, v = %f, w = %f, res = %e\n",
                it, _u, _v, 1-_u-_v, res);
    }

    if (res < tol_res) {
      converged = 1;
      if (vb) {
        log_trace("  converged (res << 1)\n");
      }
      break;
    }

    double det = a00*a11 - a01*a01;

    // TODO: check if Jacobian is singular

    double idet = 1. / det;

    double du = (b0*a11 - b1*a01)*idet;
    double dv = (b1*a00 - b0*a01)*idet;

    _u += du;
    _v += dv;
    if (vb) {
      log_trace("  du = %e, dv = %e (--> %f %f %f)\n",
                du, dv, _u, _v, 1-_u-_v);
    }

    if (_u >= 0 && _u <= 1 && _v >= 0 && _v <= 1 && _u + _v <= 1) {
      if (du*du + dv*dv < tol_uv2) {
        converged = 1; // ?
        if (vb) {
          log_trace("  converged (duv << 1)\n");
        }
        break;
      }
    }
    else {
      // Outside triangle
      converged = 1; // ?
      if (vb) {
        log_trace("  outside domain\n");
      }
      break;
    }

  } // End of Newton iteration

  if (converged) {
    *u = _u;
    *v = _v;
  }

  return converged;
}


/*=============================================================================
 * Global variables
 *============================================================================*/

/*=============================================================================
 * Public function interfaces
 *============================================================================*/

// TODO: implement de_casteljau_bar


/**
 *
 * \brief De Casteljau algorithm for Bézier triangles
 *
 * Evaluates the point (u,v) on the Bézier triangle
 * and (optionally) builds the control points of the three
 * subtriangles sharing the evaluated point as a common vertex.
 *
 * \param[in]   dim     Dimension
 * \param[in]   order   Order
 * \param[in]   u       Parametric coordinate u
 * \param[in]   v       Parametric coordinate v
 * \param[in]   b       Bézier control points
 * \param[out]  val     Evaluated point
 * \param[out]  atr     Control points of the 1st subtriangle
 * \param[out]  ars     Control points of the 2nd subtriangle
 * \param[out]  ast     Control points of the 3rd subtriangle
 *
 */

void
PDM_ho_bezier_de_casteljau_tria
(
 const int     dim,
 const int     order,
 const double  u,
 const double  v,
 double       *b,
 double       *val,
 double       *atr,
 double       *ars,
 double       *ast
 )
{
  const int n = (order+1)*(order+2)/2;

  double w = 1. - u - v;

  double p0[n*dim];
  double p1[(n-order-1)*dim];
  double *p[2] = {p0, p1};

  // initialize
  if (b != NULL) {
    memcpy(p[0], b, sizeof(double) * n * dim);
  }
  else {
    int idx = 0;
    for (int j = 0; j <= order; j++) {
      for (int i = 0; i <= order - j; i++) {
        for (int k = 0; k < dim; k++) {
          if (k == idx) {
            p[0][dim*idx + k] = 1.;
          } else {
            p[0][dim*idx + k] = 0.;
          }
        }

        idx++;
      }
    }
  }


  // subdivision
  if (atr != NULL || ars != NULL || ast != NULL) {
    int idx = 0;
    for (int j = 0; j <= order; j++) {
      for (int i = 0; i <= order-j; i++) {
        if (atr != NULL && j == 0) {
          int idx_atr = ij2idx(order-i, i, order);
          memcpy(atr + dim*idx_atr, b + dim*idx, sizeof(double)*dim);
        }
        if (ars != NULL) {
          memcpy(ars + dim*idx, b + dim*idx, sizeof(double)*dim);
        }
        if (ast != NULL && i == 0) {
          int idx_ast = ij2idx(j, order-j, order);
          memcpy(ast + dim*idx_ast, b + dim*idx, sizeof(double)*dim);
        }

        idx++;
      }
    }
  }



  for (int l = 1; l <= order; l++) {

    int idx = 0;
    for (int j = 0; j <= order-l; j++) {
      for (int i = 0; i <= order-l-j; i++) {
        int idxu = ij2idx(i+1, j,   order-l+1);
        int idxv = ij2idx(i,   j+1, order-l+1);
        int idxw = ij2idx(i,   j,   order-l+1);

        for (int k = 0; k < dim; k++) {
          p[1][dim*idx + k] =
          u*p[0][dim*idxu + k] +
          v*p[0][dim*idxv + k] +
          w*p[0][dim*idxw + k];
        }


        // subdivision
        if (atr != NULL && j == 0) {
          int idx_atr = ij2idx(order-l-i, i, order);
          memcpy(atr + dim*idx_atr, p[1] + dim*idx, sizeof(double)*dim);
        }
        if (ars != NULL) {
          int idx_ars = ij2idx(i, j, order);
          memcpy(ars + dim*idx_ars, p[1] + dim*idx, sizeof(double)*dim);
        }
        if (ast != NULL && i == 0) {
          int idx_ast = ij2idx(j, order-l-j, order);
          memcpy(ast + dim*idx_ast, p[1] + dim*idx, sizeof(double)*dim);
        }


        idx++;
      }
    }

    if (l < order) {
      // swap pointers
      double *tmp = p[1];
      p[1] = p[0];
      p[0] = tmp;
    }
  }

  if (val != NULL) {
    memcpy(val, p[1], sizeof(double) * dim);
  }
}


/**
 *
 * \brief Build control points for partial derivatives of a Bézier triangle
 *
 * \param[in]   dim     Dimension
 * \param[in]   order   Order
 * \param[in]   b       Bézier control points
 * \param[out]  bu      Bézier control points of 1st partial derivative
 * \param[out]  bv      Bézier control points of 2nd partial derivative
 *
 */

void
PDM_ho_bezier_triangle_derivatives
(
 const int     dim,
 const int     order,
       double *b,
       double *bu,
       double *bv
 )
{
  int idx = 0;
  for (int j = 0; j < order; j++) {
    for (int i = 0; i < order-j; i++) {

      int idxij  = ij2idx(i, j, order);
      int idxi1j = idxij + 1;            //ij2idx(i+1, j, order);
      int idxij1 = idxij + order - j + 1;//ij2idx(i, j+1, order);

      for (int k = 0; k < dim; k++) {
        bu[dim*idx + k] = order * (b[dim*idxi1j + k] - b[dim*idxij + k]);
        bv[dim*idx + k] = order * (b[dim*idxij1 + k] - b[dim*idxij + k]);
      }

      idx++;
    }
  }
}


/**
 *
 * \brief Point location in a high-order Bézier triangle
 *
 * \param[in]   order             Order
 * \param[in]   n_node            Number of nodes
 * \param[in]   node_coord        Coordinates of the Bézier control points (size = 3 * \ref n_node)
 * \param[in]   point_coord       Coordinates of the point to locate (size = 3)
 * \param[out]  projected_coords  Coordinates of the projection on the Bézier triangle (size = 3)
 * \param[out]  uvw               Parametric coordinates of the projection on the Bézier triangle
 *
 */

double
PDM_ho_bezier_tria_location
(
 const int     order,
 const int     n_node,
       double *node_coord,
       double *point_coord,
       double *projected_coord,
       double *uvw
 )
{
  double P1_coord[9];
  memcpy(P1_coord,     node_coord,                sizeof(double) * 3);
  memcpy(P1_coord + 3, node_coord + 3*order,      sizeof(double) * 3);
  memcpy(P1_coord + 6, node_coord + 3*(n_node-1), sizeof(double) * 3);

  double distance, weight[3];
  PDM_triangle_evaluate_position(point_coord,
                                 P1_coord,
                                 projected_coord,
                                 &distance,
                                 weight);
  uvw[0] = weight[2];
  uvw[1] = weight[0];
  uvw[2] = weight[1];

  const int n = order*(order+1)/2;
  double db_du[3*n], db_dv[3*n];
  PDM_ho_bezier_triangle_derivatives(3,
                                     order,
                                     node_coord,
                                     db_du,
                                     db_dv);

  int converged = _newton_tria(order,
                               point_coord,
                               &uvw[0],
                               &uvw[1],
                               node_coord,
                               db_du,
                               db_dv,
                               projected_coord);
  PDM_UNUSED(converged);

  uvw[2] = 1 - uvw[0] - uvw[1];

  distance = 0;
  for (int i = 0; i < 3; i++) {
    double d = point_coord[i] - projected_coord[i];
    distance += d*d;
  }

  return distance;
}





#undef ij2idx
