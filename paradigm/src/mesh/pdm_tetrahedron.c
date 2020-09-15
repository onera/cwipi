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
 * PDM library headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_error.h"
#include "pdm_printf.h"

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm_triangle.h"
#include "pdm_geom_elem.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_tetrahedron.h"

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

/**
 * \brief Evaluates the position in a tetrahedron
 *
 * \param [in]  x                         Point coordinates to evaluate position
 * \param [in]  vtx_coord                 Tetrahedron vertices coordinates
 * \param [out] closest_point             Closest Point in Tetrahedron or NULL
 * \param [out] closest_point_dist2       Square of the distance
 * \param [out] closest_point_weights     Vertices weights or NULL
 *
 * \return      -1 if the tetrahedron is degenerate, 0 else
 *
 */

int PDM_tetrahedron_evaluate_position
(
 const double  x[3],
 const double  vtx_coord[12],
 double        closest_point[3],
 double       *closest_point_dist2,
 double        closest_point_weights[4]
 )
{
  int i, j, k;

  double weights_local[4];
  double *_weights = weights_local;
  if (closest_point_weights != NULL) {
    _weights = closest_point_weights;
  }

  double cp_local[3];
  double *_closest_point = cp_local;
  if (closest_point != NULL) {
    _closest_point = closest_point;
  }

  double v[3][3];
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      v[i][j] = vtx_coord[3*(i+1) + j] - vtx_coord[j];
    }
  }

  double vol6 =
    v[0][0] * (v[1][1]*v[2][2] - v[1][2]*v[2][1]) +
    v[0][1] * (v[1][2]*v[2][0] - v[1][0]*v[2][2]) +
    v[0][2] * (v[1][0]*v[2][1] - v[1][1]*v[2][0]);

  if (fabs(vol6) < 1e-16){
    return -1;
  }

  double r[3][3];
  PDM_CROSS_PRODUCT (r[0], v[2], v[1]);
  PDM_CROSS_PRODUCT (r[1], v[0], v[2]);
  PDM_CROSS_PRODUCT (r[2], v[1], v[0]);

  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      r[i][j] /= vol6;
    }
  }

  double xv0[3] = {vtx_coord[0] - x[0],
                   vtx_coord[1] - x[1],
                   vtx_coord[2] - x[2]};

  double uvw[3], s = 1.;
  for (i = 0; i < 3; i++) {
    uvw[i] = PDM_DOT_PRODUCT (xv0, r[i]);
    s -= uvw[i];
  }

  /* Point inside tetrahedron */
  if (s      >= 0. &&
      uvw[0] >= 0. &&
      uvw[1] >= 0. &&
      uvw[2] >= 0.) {

    _weights[0] = s;
    for (j = 0; j < 3; j++) {
      _closest_point[j] = x[j];
      _weights[j+1] = uvw[j];
    }

    *closest_point_dist2 = 0.;
  }

  /* Point inside tetrahedron */
  else {

    *closest_point_dist2 = HUGE_VAL;

    int orientation = vol6 > 0.;
    const int tetra_vtx[4] = {0, 1, 2, 3};
    int tri_vtx_idx[5], tri_vtx[12];
    PDM_geom_elem_tetra_faces (1,
                               orientation,
                               tetra_vtx,
                               tri_vtx_idx,
                               tri_vtx);

    int *_tri_vtx = NULL;
    double tri_coord[9];
    double tri_weights[3];
    double tri_closest_point[3];
    double tri_dist2;

    for (i = 0; i < 3; i++) {

      if (uvw[i] < 0.) {

        _tri_vtx = tri_vtx + 3*(2-i);

        for (j = 0; j < 3; j++) {
          for (k = 0; k < 3; k++) {
            tri_coord[3*j + k] = vtx_coord[3*_tri_vtx[j] + k];
          }
        }

        PDM_triangle_evaluate_position (x,
                                        tri_coord,
                                        tri_closest_point,
                                        &tri_dist2,
                                        tri_weights);

        if (*closest_point_dist2 > tri_dist2) {
          *closest_point_dist2 = tri_dist2;

          for (j = 0; j < 3; j++) {
            _closest_point[j] = x[j];
          }


          _weights[i+1] = 0.;
          _weights[_tri_vtx[0]] = tri_weights[1];
          _weights[_tri_vtx[1]] = tri_weights[2];
          _weights[_tri_vtx[2]] = tri_weights[0];
        }

      }
    }


    if (s < 0.) {

      _tri_vtx = tri_vtx + 9;

      for (j = 0; j < 3; j++) {
        for (k = 0; k < 3; k++) {
          tri_coord[3*j + k] = vtx_coord[3*_tri_vtx[j] + k];
        }
      }

      PDM_triangle_evaluate_position (x,
                                      tri_coord,
                                      tri_closest_point,
                                      &tri_dist2,
                                      tri_weights);

      if (*closest_point_dist2 > tri_dist2) {
        *closest_point_dist2 = tri_dist2;

        for (j = 0; j < 3; j++) {
          _closest_point[j] = x[j];
        }

        _weights[0] = 0.;
        _weights[_tri_vtx[0]] = tri_weights[1];
        _weights[_tri_vtx[1]] = tri_weights[2];
        _weights[_tri_vtx[2]] = tri_weights[0];
      }

    }
  }

  return 0;
}


#ifdef __cplusplus
}
#endif /* __cplusplus */

