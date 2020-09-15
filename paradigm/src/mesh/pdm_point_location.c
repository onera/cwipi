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
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_plane.h"
#include "pdm_line.h"
#include "pdm_triangle.h"
#include "pdm_polygon.h"
#include "pdm_triangulate.h"
#include "pdm_mesh_nodal.h"
#include "pdm_mean_values.h"
#include "pdm_geom_elem.h"
#include "pdm_binary_search.h"
#include "pdm_ho_location.h"

#include "pdm_point_location.h"

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
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

static double _epsilon_denom = 1.e-30;       /* Minimum denominator */

/*=============================================================================
 * Private function definition
 *============================================================================*/

static inline double
_determinant_3x3
(
 const double a[3],
 const double b[3],
 const double c[3]
 )
{
  return a[0] * (b[1]*c[2] - b[2]*c[1])
    +    a[1] * (b[2]*c[0] - b[0]*c[2])
    +    a[2] * (b[0]*c[1] - b[1]*c[0]);
}

/*---------------------------------------------------------------------------
 * Solve the equation "A.x = b" with Cramer's rule.
 *
 * parameters:
 *   A[3][3] <-- equation matrix
 *   b[3]    <-- b equation right hand side
 *   x[3]    <-> equation solution (unchanged if matrix is singular)
 *
 * returns:
 *   PDM_FALSE if matrix is singular, PDM_TRUE otherwise
 *----------------------------------------------------------------------------*/

static int
_solve_3x3(double  A[3][3],
           double  b[3],
           double  x[3])
{
  double det, det_inv, x0, x1, x2;

  det = A[0][0]*(A[1][1]*A[2][2] - A[2][1]*A[1][2])
    -   A[1][0]*(A[0][1]*A[2][2] - A[2][1]*A[0][2])
    +   A[2][0]*(A[0][1]*A[1][2] - A[1][1]*A[0][2]);

  if (PDM_ABS(det) < _epsilon_denom) {
    printf("_solve_3x3: det = %e\n", det);
    return PDM_FALSE;
  }

  else {
    det_inv = 1./det;
  }

  /* Use local variables to ensure no aliasing */

  x0 = (  b[0]*(A[1][1]*A[2][2] - A[2][1]*A[1][2])
          - b[1]*(A[0][1]*A[2][2] - A[2][1]*A[0][2])
          + b[2]*(A[0][1]*A[1][2] - A[1][1]*A[0][2])) * det_inv;

  x1 = (  A[0][0]*(b[1]*A[2][2] - b[2]*A[1][2])
          - A[1][0]*(b[0]*A[2][2] - b[2]*A[0][2])
          + A[2][0]*(b[0]*A[1][2] - b[1]*A[0][2])) * det_inv;

  x2 = (  A[0][0]*(A[1][1]*b[2] - A[2][1]*b[1])
          - A[1][0]*(A[0][1]*b[2] - A[2][1]*b[0])
          + A[2][0]*(A[0][1]*b[1] - A[1][1]*b[0])) * det_inv;

  /* Copy local variables to output */
  x[0] = x0;
  x[1] = x1;
  x[2] = x2;

  return PDM_TRUE;
}


/*----------------------------------------------------------------------------
 * Compute 3d shape functions and their derivatives given element
 * parametric coordinates.
 *
 * This function is adapted from the CGNS interpolation tool.
 *
 * parameters:
 *   elt_type    <-- type of element
 *   uvw[]       <-- parametric coordinates
 *   shapef[]    --> barycenter's coordinates
 *   deriv [][]  --> derivative of shape function
 *----------------------------------------------------------------------------*/

static void
_compute_shapef_3d
(
 const PDM_Mesh_nodal_elt_t elt_type,
 const double               uvw[3],
 double                     shapef[8],
 double                     deriv[8][3]
 )
{
  switch (elt_type) {


  case PDM_MESH_NODAL_PYRAMID5: {

    double u = uvw[0];
    double v = uvw[1];
    double w = uvw[2];

    double w1 = 1. - w;
    if (fabs(w1) > 1.e-6) {
      w1 = 1. / w1;
    }

    shapef[0] = (1. - u - w) * (1. - v - w) * w1;
    shapef[1] =            u * (1. - v - w) * w1;
    shapef[2] =            u *            v * w1;
    shapef[3] = (1. - u - w) *            v * w1;
    shapef[4] = w;

    if (deriv != NULL) {
      deriv[0][0] = (v + w - 1.) * w1;
      deriv[0][1] = (u + w - 1.) * w1;
      deriv[0][2] = shapef[0]*w1 + deriv[0][0] + deriv[0][1];

      deriv[1][0] = (1. - v - w) * w1;
      deriv[1][1] = -u * w1;
      deriv[1][2] = shapef[1]*w1 + deriv[1][1];

      deriv[2][0] = v * w1;
      deriv[2][1] = u * w1;
      deriv[2][2] = shapef[2]*w1;

      deriv[3][0] = -v * w1;
      deriv[3][1] = (1. - u - w) * w1;
      deriv[3][2] = shapef[3]*w1 + deriv[3][0];

      deriv[4][0] = 0.;
      deriv[4][1] = 0.;
      deriv[4][2] = 1.;
    }

    break;
  }



  case PDM_MESH_NODAL_PRISM6: {

    shapef[0] = (1.0 - uvw[0] - uvw[1]) * (1.0 - uvw[2]);
    shapef[1] = uvw[0] * (1.0 - uvw[2]);
    shapef[2] = uvw[1] * (1.0 - uvw[2]);
    shapef[3] = (1.0 - uvw[0] - uvw[1]) * uvw[2];
    shapef[4] = uvw[0] * uvw[2];
    shapef[5] = uvw[1] * uvw[2];

    if (deriv != NULL) {
      deriv[0][0] = -(1.0 - uvw[2]);
      deriv[0][1] = -(1.0 - uvw[2]);
      deriv[0][2] = -(1.0 - uvw[0] - uvw[1]);
      deriv[1][0] =  (1.0 - uvw[2]);
      deriv[1][1] =  0.0;
      deriv[1][2] = -uvw[0];
      deriv[2][0] =  0.0;
      deriv[2][1] =  (1.0 - uvw[2]);
      deriv[2][2] = -uvw[1];
      deriv[3][0] = -uvw[2];
      deriv[3][1] = -uvw[2];
      deriv[3][2] =  (1.0 - uvw[0] - uvw[1]);
      deriv[4][0] =  uvw[2];
      deriv[4][1] =  0.0;
      deriv[4][2] =  uvw[0];
      deriv[5][0] =  0.0;
      deriv[5][1] =  uvw[2];
      deriv[5][2] =  uvw[1];
    }

    break;
  }



  case PDM_MESH_NODAL_HEXA8: {

    shapef[0] = (1.0 - uvw[0]) * (1.0 - uvw[1]) * (1.0 - uvw[2]);
    shapef[1] =        uvw[0]  * (1.0 - uvw[1]) * (1.0 - uvw[2]);
    shapef[2] =        uvw[0]  *        uvw[1]  * (1.0 - uvw[2]);
    shapef[3] = (1.0 - uvw[0]) *        uvw[1]  * (1.0 - uvw[2]);
    shapef[4] = (1.0 - uvw[0]) * (1.0 - uvw[1]) *        uvw[2];
    shapef[5] =        uvw[0]  * (1.0 - uvw[1]) *        uvw[2];
    shapef[6] =        uvw[0]  *        uvw[1]  *        uvw[2];
    shapef[7] = (1.0 - uvw[0]) *        uvw[1]  *        uvw[2];

    if (deriv != NULL) {
      deriv[0][0] = -(1.0 - uvw[1]) * (1.0 - uvw[2]);
      deriv[0][1] = -(1.0 - uvw[0]) * (1.0 - uvw[2]);
      deriv[0][2] = -(1.0 - uvw[0]) * (1.0 - uvw[1]);
      deriv[1][0] =  (1.0 - uvw[1]) * (1.0 - uvw[2]);
      deriv[1][1] = -uvw[0] * (1.0 - uvw[2]);
      deriv[1][2] = -uvw[0] * (1.0 - uvw[1]);
      deriv[2][0] =  uvw[1] * (1.0 - uvw[2]);
      deriv[2][1] =  uvw[0] * (1.0 - uvw[2]);
      deriv[2][2] = -uvw[0] * uvw[1];
      deriv[3][0] = -uvw[1] * (1.0 - uvw[2]);
      deriv[3][1] =  (1.0 - uvw[0]) * (1.0 - uvw[2]);
      deriv[3][2] = -(1.0 - uvw[0]) * uvw[1];
      deriv[4][0] = -(1.0 - uvw[1]) * uvw[2];
      deriv[4][1] = -(1.0 - uvw[0]) * uvw[2];
      deriv[4][2] =  (1.0 - uvw[0]) * (1.0 - uvw[1]);
      deriv[5][0] =  (1.0 - uvw[1]) * uvw[2];
      deriv[5][1] = -uvw[0] * uvw[2];
      deriv[5][2] =  uvw[0] * (1.0 - uvw[1]);
      deriv[6][0] =  uvw[1] * uvw[2];
      deriv[6][1] =  uvw[0] * uvw[2];
      deriv[6][2] =  uvw[0] * uvw[1];
      deriv[7][0] = -uvw[1] * uvw[2];
      deriv[7][1] =  (1.0 - uvw[0]) * uvw[2];
      deriv[7][2] =  (1.0 - uvw[0]) * uvw[1];
    }

    break;
  }



  default:
    PDM_error (__FILE__, __LINE__, 0, "Wrong element type\n");

  }

}


/*----------------------------------------------------------------------------
 * Compute hexahedron, pyramid, or prism parametric coordinates for a
 * given point.
 *
 * This function is adapted from the CGNS interpolation tool.
 *
 * parameters:
 *   elt_type            <-- type of element
 *   point_coords        <-- point coordinates
 *   vertex_coords[]     <-- pointer to element vertex coordinates
 *   tolerance           <-- location tolerance factor
 *   uvw[]               --> parametric coordinates of point in element
 *----------------------------------------------------------------------------*/

static PDM_bool_t
_compute_uvw
(
 const PDM_Mesh_nodal_elt_t elt_type,
 const double               point_coords[3],
 const double               vertex_coords[],
 const double               tolerance,
 double                     uvw[3]
 )
{
  int i, j, n_elt_vertices, iter;
  const int max_iter = 30;
  const double tolerance2 = tolerance * tolerance;
  double dist;
  double a[3][3], b[3], x[3], shapef[8], dw[8][3];

  /* Get number of vertices */
  const int order = 1;
  n_elt_vertices = PDM_Mesh_nodal_n_vertices_element (elt_type, order);

  assert (elt_type == PDM_MESH_NODAL_PYRAMID5 ||
          elt_type == PDM_MESH_NODAL_PRISM6   ||
          elt_type == PDM_MESH_NODAL_HEXA8);

  /* Use Newton-method to determine parametric coordinates and shape function */
  for (i = 0; i < 3; i++) {
    uvw[i] = 0.5;
  }

  for (iter = 0; iter < max_iter; iter++) {

    _compute_shapef_3d (elt_type, uvw, shapef, dw);

    b[0] = - point_coords[0];
    b[1] = - point_coords[1];
    b[2] = - point_coords[2];

    for (i = 0; i < 3; i++) {
      for (j = 0; j < 3; j++) {
        a[i][j] = 0.0;
      }
    }

    for (i = 0; i < n_elt_vertices; i++) {

      b[0] += (shapef[i] * vertex_coords[3*i]);
      b[1] += (shapef[i] * vertex_coords[3*i+1]);
      b[2] += (shapef[i] * vertex_coords[3*i+2]);

      for (j = 0; j < 3; j++) {
        a[0][j] -= (dw[i][j] * vertex_coords[3*i]);
        a[1][j] -= (dw[i][j] * vertex_coords[3*i+1]);
        a[2][j] -= (dw[i][j] * vertex_coords[3*i+2]);
      }

    }

    if (_solve_3x3(a, b, x) == PDM_FALSE) {
      printf("_compute_uvw: singular matrix\n");
      return PDM_FALSE;
    }

    dist = 0.0;

    for (i = 0; i < 3; i++) {
      dist += x[i] * x[i];
      uvw[i] += x[i];
    }

    if (dist <= tolerance2) {
      return PDM_TRUE;
    }

  }

  return PDM_FALSE;
}


/*----------------------------------------------------------------------------
 * Locate points on an edge.
 *
 * parameters:
 *   edge_vtx           <-- ids of edge vertices
 *   parent_vertex_num  <-- pointer to parent vertex numbers (or NULL)
 *   vtx_coord          <-- pointer to vertex coordinates
 *   n_pts              <-- number of points to locate
 *   pts_coord          <-- coordinates of points to locate (size: dim * n_pts)
 *   distance           --> distance from point to edge (size: n_pts)
 *   bar_coord          --> barycentric coordinates of closest points (size: 2*n_pts)
 *----------------------------------------------------------------------------*/
static void
_locate_on_edge
(
 const double       vtx_coord[],
 const int          n_pts,
 const double       pts_coord[],
 float              distance[],
 double             bar_coord[]
 )
{
  int idim, ipt;

  /* Calculate edge vector and length */
  double u[3], uu = 0.;
  for (idim = 0; idim < 3; idim++) {
    u[idim] = vtx_coord[3 + idim] - vtx_coord[idim];
    uu += u[idim] * u[idim];
  }

  if (uu < _epsilon_denom){
    PDM_printf("warning _locate_on_edge : Reduce _epsilon_denom criteria : %12.5e < %12.5e\n", uu, _epsilon_denom);
    PDM_printf_flush();
    return;
  }

  double inv_uu = 1. / uu;

  /* Loop on points to locate on edge */
  double v[3], uv = 0., t;
  for (ipt = 0; ipt < n_pts; ipt++) {

    const double *_pt = pts_coord + 3 * ipt;
    double       *_bc = bar_coord + 2 * ipt;

    /* Calculate linear coordinates of projection of point on edge axis */
    for (idim = 0; idim < 3; idim++) {
      v[idim] = _pt[idim] - vtx_coord[idim];
      uv += u[idim] * v[idim];
    }

    t = uv * inv_uu;


    /* Set v to be the vector from the point to the closest point on
       the segment (if t < 0, v is already that vector) */
    if (t >= 1.) {
      for (idim = 0; idim < 3; idim++) {
        v[idim] = _pt[idim] - vtx_coord[3 + idim];
      }
    }

    else if (t > 0.) {
      for (idim = 0; idim < 3; idim++) {
        v[idim] -= t * u[idim];
      }
    }

    /* Distance between point to locate and its projection */
    double dist2 = 0.;
    for (idim = 0; idim < 3; idim++) {
      dist2 += v[idim] * v[idim];
    }
    distance[ipt] = (float) dist2;

    /* Barycentric coordinates */
    _bc[0] = 1. - t;
    _bc[1] =      t;

  } // End of loop on points

}

/*----------------------------------------------------------------------------
 * Locate points in a given set of triangles.
 *
 * This function is called for sets of triangles belonging to the subdivision
 * of a given 3d face. Barycentric coordinates are used to locate the
 * projection of points.
 *
 * parameters:
 *   n_tri               <-- number of triangles
 *   tri_vtx             <-- triangles connectivity (size: n_tri * 3)
 *   parent_vertex_num   <-- pointer to parent vertex numbers (or NULL)
 *   vtx_coord           <-- pointer to vertex coordinates (size: 9)
 *   n_pts               <-- number of points to locate
 *   pts_coord           <-- point coordinates (size: 3 * n_pts)
 *   location            <-> lnum of element containing or closest to each point (size: n_pts)
 *   distance            <-> distance from point to element (size: n_pts)
 *   bar_coord           <-> barcyentric coordinates of closest points (size: n_pts * 3)
 *----------------------------------------------------------------------------*/

static void
_locate_on_triangles
(
 const int          n_tri,
 const PDM_l_num_t  tri_vtx[],
 const double       vtx_coord[],
 const int          n_pts,
 const double       pts_coord[],
 int                location[],
 float              distance[],
 double             bar_coord[]
 )
{
  int itri, ipt, ivtx, idim;

  /* Initialize distance of points to locate */
  for (ipt = 0; ipt < n_pts; ipt++) {
    distance[ipt] = HUGE_VAL;
  }

  if (location != NULL) {
    for (ipt = 0; ipt < n_pts; ipt++) {
      location[ipt] = -1;
    }
  }

  /* const int _order = 1;
     const int n_vtx_tri = (_order+1)*(_order+2)/2;*/
  const int n_vtx_tri = 3;

  double tri_coord[9];

  PDM_l_num_t id[3];

  double weights[3];

  /* Loop on triangles */
  for (itri = 0; itri < n_tri; itri++) {

    /* vertex index of current triangle */
    for (ivtx = 0; ivtx < 3; ivtx++) {
      id[ivtx] = tri_vtx[itri*n_vtx_tri + ivtx] - 1;
    }

    /* Calculate triangle-constant values for barycentric coordinates */
    for (ivtx = 0; ivtx < 3; ivtx++) {
      for (idim = 0; idim < 3; idim++) {
        tri_coord[3*ivtx + idim] = vtx_coord[3*id[ivtx] + idim];
      }
    }

    /* Loop on points to locate */
    for (ipt = 0; ipt < n_pts; ipt++) {

      const double *_pt = pts_coord + 3 * ipt;
      double       *_bc = bar_coord + 3 * ipt;

      double closest_point[3];
      double dist2;

      PDM_triangle_status_t stat = PDM_triangle_evaluate_position (_pt,
                                                                   tri_coord,
                                                                   closest_point,
                                                                   &dist2,
                                                                   weights);
      if (stat == PDM_TRIANGLE_DEGENERATED) {
        continue;
      }

      if (dist2 < distance[ipt]) {
        if (bar_coord != NULL) {
          _bc[0] = weights[1];
          _bc[1] = weights[2];
          _bc[2] = weights[0];
        }

        distance[ipt] = (float) dist2;

        if (location != NULL) {
          location[ipt] = itri;
        }
      }

    } // End of loop on points

  } // End of loop on triangles
}


/*----------------------------------------------------------------------------
 * Locate points in a given quadrangle.
 *
 * Barycentric coordinates are used to locate the projection of points.
 *
 * parameters:
 *   quad_vtx            <-- quadrangle connectivity (size: 4)
 *   parent_vertex_num   <-- pointer to parent vertex numbers (or NULL)
 *   vtx_coord           <-- pointer to vertex coordinates (size: dim * 4)
 *   n_pts               <-- number of points to locate
 *   pts_coord           <-- point coordinates (size: dim * n_pts)
 *   distance            <-> distance from point to element (size: n_pts)
 *   bar_coord           <-> barcyentric coordinates of closest points (size: n_pts * 4)
 *----------------------------------------------------------------------------*/
static void
_locate_on_quadrangle
(
 const double       quad_coord[],
 const int          n_pts,
 const double       pts_coord[],
 float              distance[],
 double             bar_coord[]
 )
{
  int ipt, ivtx, idim;


  PDM_mean_value_coordinates_polygon_3d (4,
                                         quad_coord,
                                         n_pts,
                                         pts_coord,
                                         bar_coord);


  for (ipt = 0; ipt < n_pts; ipt++) {

    const double *_pt = pts_coord + 3 * ipt;
    double *_bc = bar_coord + 4 * ipt;

    double v_cp_p[3];
    for (idim = 0; idim < 3; idim++) {
      v_cp_p[idim] = _pt[idim];
    }

    for (ivtx = 0; ivtx < 4; ivtx++) {
      for (idim = 0; idim < 3; idim++) {
        v_cp_p[idim] -= _bc[ivtx] * quad_coord[3*ivtx + idim];
      }
    }

    double dist2 = 0.;
    for (idim = 0; idim < 3; idim++) {
      dist2 += v_cp_p[idim] * v_cp_p[idim];
    }

    distance[ipt] = (float) dist2;
  }
}

/*----------------------------------------------------------------------------
 * Locate points in a tetrahedron.
 *
 * parameters:
 *   tetra_coords        <-- tetrahedra vertex coordinates
 *   n_pts               <-- number of points to locate
 *   point_coords        <-- point coordinates
 *   distance            <-> distance from point to element (size: n_pts)
 *   bar_coord           <-> barcyentric coordinates of closest points (size: n_pts * 4)
 *----------------------------------------------------------------------------*/
static void
_locate_in_tetrahedron
(
 const double        tetra_coord[12],
 const int           n_pts,
 const double        pts_coord[],
 float              *distance,
 double             *bar_coord
 )
{
  int ivtx, idim, ipt, i, j, k;

  int n_pts_out = 0;
  int *pts_out = malloc (sizeof(int) * n_pts);

  double v[3][3];
  for (ivtx = 0; ivtx < 3; ivtx++) {
    for (idim = 0; idim < 3; idim++) {
      v[ivtx][idim] = tetra_coord[3*(ivtx+1) + idim] - tetra_coord[idim];
    }
  }

  double vol6 = v[0][0] * (v[1][1]*v[2][2] - v[1][2]*v[2][1]) +
    v[0][1] * (v[1][2]*v[2][0] - v[1][0]*v[2][2]) +
    v[0][2] * (v[1][0]*v[2][1] - v[1][1]*v[2][0]);


  if (fabs(vol6) < _epsilon_denom){
    PDM_printf("warning _locate_in_tetrahedron : Reduce _epsilon_denom criteria : %12.5e < %12.5e\n", vol6, _epsilon_denom);
    PDM_printf_flush();
    return;
  }

  int orientation = vol6 > 0.;

  double r[3][3];
  for (i = 0; i < 3; i++) {
    j = (i + 1) % 3;
    k = (i + 2) % 3;

    PDM_CROSS_PRODUCT (r[i], v[k], v[j]);

    for (idim = 0; idim < 3; idim++) {
      r[i][idim] /= vol6;
    }
  }

  /*
   *  Locate points inside and identify points outside
   */
  double vp0[3];
  for (ipt = 0; ipt < n_pts; ipt++) {
    const double *_pt = pts_coord + 3 * ipt;
    double       *_bc = bar_coord + 4 * ipt;

    for (idim = 0; idim < 3; idim++) {
      vp0[idim] = tetra_coord[idim] - _pt[idim];
    }

    /* Compute barycentric coordinates of current point in tetrahedron */
    _bc[0] = 1.;
    for (i = 0; i < 3; i++) {
      _bc[i+1] = PDM_DOT_PRODUCT (vp0, r[i]);
      _bc[0] -= _bc[i+1];
    }

    /* Compute 'distance' from current point to tetrahedron */
    double min_bc = HUGE_VAL;
    for (ivtx = 0; ivtx < 4; ivtx++) {
      min_bc = PDM_MIN (min_bc, _bc[ivtx]);
    }

    distance[ipt] = (float) (min_bc * min_bc);
      if (min_bc > 0.) {
        distance[ipt] = -distance[ipt];
      }

    /* Point outside tetrahedron */
    if (distance[ipt] > 0.) {
      pts_out[n_pts_out++] = ipt;
    }

  } // End loop on points


  if (n_pts_out == 0) {
    free (pts_out);
    return;
  }


  /*
   *  Locate points outside (closest points on boundary)
   */
  double *pts_out_coord = malloc (sizeof(double) * n_pts_out * 3);
  for (ipt = 0; ipt < n_pts_out; ipt++) {
    int id_pt = pts_out[ipt];
    for (idim = 0; idim < 3; idim++) {
      pts_out_coord[3*ipt + idim] = pts_coord[3*id_pt + idim];
    }
  }

  /* Get tetrahedron's faces */
  const int cell_vtx[4] = {1, 2, 3, 4};
  int face_vtx_idx[5], face_vtx[12];
  PDM_geom_elem_tetra_faces (1,
                             orientation,
                             cell_vtx,
                             face_vtx_idx,
                             face_vtx);

  int *id_face = malloc (sizeof(int) * n_pts_out);
  double *bar_coord_face = malloc (sizeof(double) * n_pts_out * 3);
  float *distance_face = malloc (sizeof(float) * n_pts_out);
  _locate_on_triangles (4,
                        face_vtx,
                        tetra_coord,
                        n_pts_out,
                        pts_out_coord,
                        id_face,
                        distance_face,
                        bar_coord_face);

  for (ipt = 0; ipt < n_pts_out; ipt++) {
    int id_pt = pts_out[ipt];
    double *_bc = bar_coord + 4 * id_pt;

    for (ivtx = 0; ivtx < 4; ivtx++) {
      _bc[ivtx] = 0.;
    }

    for (ivtx = 0; ivtx < 3; ivtx++) {
      int id_vtx = face_vtx[face_vtx_idx[id_face[ipt]] + ivtx] - 1;
      _bc[id_vtx] = bar_coord_face[3*ipt + ivtx];
    }

    distance[id_pt] = distance_face[ipt];
  }

  free (pts_out);
  free (pts_out_coord);
  free (id_face);
  free (bar_coord_face);
  free (distance_face);
}


/*----------------------------------------------------------------------------
 * Locate points in a standard cell (i.e. tetrahedron, pyramid, prism or hexahedron).
 *
 * parameters:
 *   elt_type            <-- type of cell (tetrahedron, pyramid, prism or hexahedron)
 *   cell_vtx            <-- cell-vertex connectivity
 *   parent_vertex_num   <-- pointer to parent vertex numbers (or NULL)
 *   vtx_coord           <-- pointer to vertex coordinates
 *   n_pts               <-- number of points to locate
 *   pts_coord           <-- point coordinates (size: dim * n_pts)
 *   tolerance           <-- tolerance (used to check coincidence between points and vertices)
 *   distance            <-> distance from point to element (size: n_pts)
 *   bar_coord           <-> barcyentric coordinates of closest points (size: n_pts * 4)
 *----------------------------------------------------------------------------*/
static void
_locate_in_cell_3d
(
 const PDM_Mesh_nodal_elt_t  elt_type,
 /*const PDM_l_num_t           cell_vtx[],
   const PDM_l_num_t          *parent_vertex_num,*/
 const double                cell_coord[],
 const int                   n_pts,
 const double                pts_coord[],
 const double                tolerance,
 float                      *distance,
 double                     *bar_coord
 )
{
  double uvw[3];
  double eps_vtx2 = 1.e-6 * tolerance;
  eps_vtx2 *= eps_vtx2;

  const int order = 1;
  const int n_vtx = PDM_Mesh_nodal_n_vertices_element (elt_type, order);

  int *pts_out = malloc (sizeof(int) * n_pts);
  int n_pts_out = 0;



  /* Tetrahedron */
  if (elt_type == PDM_MESH_NODAL_TETRA4) {
    /* Shape functions may be computed directly */
    _locate_in_tetrahedron (cell_coord,
                            n_pts,
                            pts_coord,
                            distance,
                            bar_coord);
    return;
  }

  double *_cell_coord = NULL;
  if (elt_type == PDM_MESH_NODAL_PRISM6) {
    _cell_coord = (double *) cell_coord;
  }
  else {
    _cell_coord = malloc (sizeof(double) * n_vtx * 3);
    _cell_coord[ 0] = cell_coord[ 0];
    _cell_coord[ 1] = cell_coord[ 1];
    _cell_coord[ 2] = cell_coord[ 2];

    _cell_coord[ 3] = cell_coord[ 6];
    _cell_coord[ 4] = cell_coord[ 7];
    _cell_coord[ 5] = cell_coord[ 8];

    _cell_coord[ 6] = cell_coord[ 3];
    _cell_coord[ 7] = cell_coord[ 4];
    _cell_coord[ 8] = cell_coord[ 5];

    _cell_coord[ 9] = cell_coord[ 9];
    _cell_coord[10] = cell_coord[10];
    _cell_coord[11] = cell_coord[11];

    _cell_coord[12] = cell_coord[12];
    _cell_coord[13] = cell_coord[13];
    _cell_coord[14] = cell_coord[14];

    if (elt_type == PDM_MESH_NODAL_HEXA8) {
      _cell_coord[15] = cell_coord[18];
      _cell_coord[16] = cell_coord[19];
      _cell_coord[17] = cell_coord[20];

      _cell_coord[18] = cell_coord[15];
      _cell_coord[19] = cell_coord[16];
      _cell_coord[20] = cell_coord[17];

      _cell_coord[21] = cell_coord[21];
      _cell_coord[22] = cell_coord[22];
      _cell_coord[23] = cell_coord[23];
    }
  }


  /* Other cell types, shape functions must be computed iteratively */
  for (int ipt = 0; ipt < n_pts; ipt++) {

    const double *_pt = pts_coord + 3 * ipt;
    double *_bc = bar_coord + n_vtx * ipt;

    /* Check vertices (To avoid singularity with pyramids) */
    PDM_bool_t on_vtx = PDM_FALSE;
    for (int ivtx = 0; ivtx < n_vtx; ivtx++) {

      double dist2 = 0.;
      for (int idim = 0; idim < 3; idim++) {
        double delta = cell_coord[3*ivtx + idim] - _pt[idim];
        dist2 += delta * delta;
      }

      if (dist2 < eps_vtx2) {
        distance[ipt] = 0.;

        for (int i = 0; i < n_vtx; i++) {
          _bc[i] = 0.;
        }
        _bc[ivtx] = 1.;

        on_vtx = PDM_TRUE;
        break;
      }

    } // End loop on vertices

    if (on_vtx == PDM_TRUE) {
      continue;
    }

    /* Compute parametric coordinates */
    PDM_bool_t stat_uvw = _compute_uvw (elt_type,
                                        _pt,
                                        cell_coord,
                                        tolerance,
                                        uvw);
    if (stat_uvw == PDM_TRUE) {
      _compute_shapef_3d (elt_type, uvw, _bc, NULL);

      double min_bc = HUGE_VAL;
      for (int ivtx = 0; ivtx < n_vtx; ivtx++) {
        min_bc = PDM_MIN (min_bc, _bc[ivtx]);
      }

      distance[ipt] = (float) (min_bc * min_bc);
      if (min_bc > 0.) {
        distance[ipt] = -distance[ipt];
      }
    }

    /* Failed to compute parametric coordinates */
    else {
      /*
       * Use hierarchical subdivision
       */
      double proj_coord[3];
      double dist2 = PDM_ho_location (elt_type,
                                      order,
                                      n_vtx,
                                      _cell_coord,
                                      _pt,
                                      proj_coord,
                                      uvw);
      /* Point inside */
      if (dist2 < 1.e-12) {
        _compute_shapef_3d (elt_type, uvw, _bc, NULL);

        double min_bc = HUGE_VAL;
        for (int ivtx = 0; ivtx < n_vtx; ivtx++) {
          min_bc = PDM_MIN (min_bc, _bc[ivtx]);
        }

        distance[ipt] = (float) (-min_bc * min_bc);
      }

      /* Point outside */
      else {
        distance[ipt] = (float) dist2;
      }
    }

    /* Point outside cell */
    if (distance[ipt] > 0.) {
      pts_out[n_pts_out++] = ipt;
      distance[ipt] = HUGE_VAL;
    }

  } // End loop on points

  if (elt_type != PDM_MESH_NODAL_PRISM6) {
    free (_cell_coord);
  }


  /* Locate points outside cell (closest point on boundary) */
  if (n_pts_out > 0) {

    int    *closest_face  = malloc (sizeof(int)    * n_pts_out);
    double *closest_point = malloc (sizeof(double) * n_pts_out * 3);

    int n_face = -1;
    int face_vtx_idx[7];
    int face_vtx[24];
    int _cell_vtx[8] = {1, 2, 3, 4, 5, 6, 7, 8};
    int orientation = 0;
    int n_vtx_face;

    switch (elt_type) {
    case PDM_MESH_NODAL_PYRAMID5:
      n_face = 5;

      PDM_geom_elem_pyramid_faces (1,
                                   orientation,
                                   _cell_vtx,
                                   face_vtx_idx,
                                   face_vtx);
      break;

    case PDM_MESH_NODAL_PRISM6:
      n_face = 5;

      PDM_geom_elem_prism_faces (1,
                                 orientation,
                                 _cell_vtx,
                                 face_vtx_idx,
                                 face_vtx);
      break;

    case PDM_MESH_NODAL_HEXA8:
      n_face = 6;

      PDM_geom_elem_hexa_faces (1,
                                orientation,
                                _cell_vtx,
                                face_vtx_idx,
                                face_vtx);
      break;

    default:
      PDM_error (__FILE__, __LINE__, 0, "Wrong standard element type\n");
      break;

    }

    PDM_l_num_t n_tri;
    double tri_coord[9];
    PDM_l_num_t _tri_vtx[3], tri_vtx[6];
    PDM_triangulate_state_t *state = PDM_triangulate_state_create (4);

    /* Find closest face/point for each point outside the cell */
    for (int iface = 0; iface < n_face; iface++) {

      const PDM_l_num_t *_face_vtx = face_vtx + face_vtx_idx[iface];
      n_vtx_face = face_vtx_idx[iface+1] - face_vtx_idx[iface];

      /* Triangular face */
      if (n_vtx_face == 3) {
        n_tri = 1;
        for (int ivtx = 0; ivtx < 3; ivtx++) {
          tri_vtx[ivtx] = _face_vtx[ivtx];
        }
      }

      /* Quadrilateral face */
      else {
        n_tri = PDM_triangulate_quadrangle (3,
                                            (double *) cell_coord,
                                            NULL,
                                            _face_vtx,
                                            tri_vtx);
      }

      /* Loop on face triangles */
      for (int itri = 0; itri < n_tri; itri++) {

        for (int ivtx = 0; ivtx < 3; ivtx++) {
          _tri_vtx[ivtx] = tri_vtx[3*itri + ivtx] - 1;
        }

        for (int ivtx = 0; ivtx < 3; ivtx++) {
          for (int idim = 0; idim < 3; idim++) {
            tri_coord[3*ivtx + idim] = cell_coord[3*_tri_vtx[ivtx] + idim];
          }
        }

        /* Loop on points */
        for (int ipt = 0; ipt < n_pts_out; ipt++) {

          int _ipt = pts_out[ipt];
          const double *_pt = pts_coord + 3 * _ipt;
          double *_cp = closest_point + 3 * ipt;

          double min_dist2, closest[3];
          PDM_triangle_status_t error = PDM_triangle_evaluate_position (_pt,
                                                                        tri_coord,
                                                                        closest,
                                                                        &min_dist2,
                                                                        NULL);

          if (error == PDM_TRIANGLE_DEGENERATED) {
            continue;
          }

          if (distance[_ipt] > min_dist2) {
            distance[_ipt] = (float) min_dist2;

            closest_face[ipt] = iface;
            for (int idim = 0; idim < 3; idim++) {
              _cp[idim] = closest[idim];
            }
          }

        } // End loop on points

      } // End loop on triangles

    } // End loop on faces

    state = PDM_triangulate_state_destroy(state);

    /* Loctate closest points on closest faces */
    double face_coord[12];
    double bar_coord_face[4];
    for (int ipt = 0; ipt < n_pts_out; ipt++) {

      int _ipt = pts_out[ipt];
      double *_cp = closest_point + 3 * ipt;
      double *_bc = bar_coord + n_vtx * _ipt;

      for (int ivtx = 0; ivtx < n_vtx; ivtx++) {
        _bc[ivtx] = 0.;
      }

      int iface = closest_face[ipt];
      n_vtx_face = face_vtx_idx[iface+1] - face_vtx_idx[iface];

      for (int ivtx = 0; ivtx < n_vtx_face; ivtx++) {
        int _ivtx = face_vtx[face_vtx_idx[iface] + ivtx] - 1;

        for (int idim = 0; idim < 3; idim++) {
          face_coord[3*ivtx + idim] = cell_coord[3*_ivtx + idim];
        }
      }

      PDM_mean_value_coordinates_polygon_3d (n_vtx_face,
                                             face_coord,
                                             1,
                                             _cp,
                                             bar_coord_face);

      for (int ivtx = 0; ivtx < n_vtx_face; ivtx++) {
        int _ivtx = face_vtx[face_vtx_idx[iface] + ivtx] - 1;
        _bc[_ivtx] = bar_coord_face[ivtx];
      }

      const double *_pt = pts_coord + 3 * _ipt;
      double v_p_cp[3] = {_pt[0], _pt[1], _pt[2]};
      for (int ivtx = 0; ivtx < n_vtx_face; ivtx++) {
        int _ivtx = face_vtx[face_vtx_idx[iface] + ivtx] - 1;
        for (int idim = 0; idim < 3; idim++) {
          v_p_cp[idim] -= bar_coord_face[ivtx] * cell_coord[3*_ivtx + idim];
        }
      }

      double min_dist2 = PDM_DOT_PRODUCT (v_p_cp, v_p_cp);
      distance[_ipt] = (float) min_dist2;
    }

    free (closest_face);
    free (closest_point);
  }

  free (pts_out);

}



/*----------------------------------------------------------------------------
 * Locate points in a given 3d polygon.
 *
 * Barycentric coordinates are used to locate the projection of points.
 *
 * parameters:
 *   n_vtx               <-- number of vertices
 *   vtx_coord           <-- pointer to vertex coordinates (size: 3 * n_vtx)
 *   n_pts               <-- number of points to locate
 *   pts_coord           <-- point coordinates (size: 3 * n_pts)
 *   distance            <-> distance from point to element (size: n_pts)
 *   bar_coord           <-> barcyentric coordinates of closest points (size: n_pts * n_vtx)
 *----------------------------------------------------------------------------*/
static void
_locate_in_polygon
(
 const PDM_l_num_t n_vtx,
 const double      vtx_coord[],
 const int         n_pts,
 const double      pts_coord[],
 float             distance[],
 double            bar_coord[]
 )
{
  /* Compute mean value coordinates of closest points on polygon */
  PDM_mean_value_coordinates_polygon_3d (n_vtx,
                                         vtx_coord,
                                         n_pts,
                                         pts_coord,
                                         bar_coord);

  /* Compute distances */
  for (int ipt = 0; ipt < n_pts; ipt++) {
    const double *_pt = pts_coord + ipt * 3;
    double       *_bc = bar_coord + ipt * n_vtx;

    double v_cp_p[3] = {_pt[0], _pt[1], _pt[2]};
    for (int ivtx = 0; ivtx < n_vtx; ivtx++) {
      for (int idim = 0; idim < 3; idim++) {
        v_cp_p[idim] -= _bc[ivtx] * vtx_coord[3*ivtx + idim];
      }
    }
    double dist2 = PDM_DOT_PRODUCT (v_cp_p, v_cp_p);
    distance[ipt] = (float) dist2;
  }
}


/*----------------------------------------------------------------------------
 * Locate points in a given polyhedron.
 *
 * parameters:
 *   n_vtx               <-- number of vertices
 *   vtx_coord           <-- pointer to vertex coordinates (size: 3 * n_vtx)
 *   n_face              <-- number of faces
 *   face_vtx_idx        <-- index of face-vertex connectivity (size: n_face + 1)
 *   face_vtx            <-- face-vertex connectivity (size: face_vtx_idx[n_face])
 *   face_orientation    <-- orientation of faces in current polyhedron (size: n_face)
 *   n_pts               <-- number of points to locate
 *   pts_coord           <-- point coordinates (size: 3 * n_pts)
 *   char_length         <-- characteristic length of polyhedron
 *   distance            <-> distance from point to element (size: n_pts)
 *   bar_coord           <-> barcyentric coordinates of closest points (size: n_pts * n_vtx)
 *----------------------------------------------------------------------------*/
static int
_locate_in_polyhedron
(
 const PDM_l_num_t n_vtx,
 const double      vtx_coord[],
 const PDM_l_num_t n_face,
 const PDM_l_num_t face_vtx_idx[],
 const PDM_l_num_t face_vtx[],
 const int         face_orientation[],
 const int         n_pts,
 const double      pts_coord[],
 float             distance[],
 double            bar_coord[]
 )
{
  const double four_PI = 4. * PDM_PI;

  const double eps_solid_angle = 1e-5;
  const float eps_distance2 = 1e-10;//char_length...
  const double eps2 = 1e-20;

  /*
   * Identify points inside/outside polyhedron
   */
  double *solid_angle = malloc (sizeof(double) * n_pts);
  for (int ipt = 0; ipt < n_pts; ipt++) {
    solid_angle[ipt] = 0.;
    distance[ipt] = HUGE_VAL;
  }

  /* Count max nb of vertices per face */
  PDM_l_num_t n_vtx_face, n_vtx_face_max = 0;
  for (int iface = 0; iface < n_face; iface++) {
    n_vtx_face = face_vtx_idx[iface+1] - face_vtx_idx[iface];
    if (n_vtx_face > n_vtx_face_max) {
      n_vtx_face_max = n_vtx_face;
    }
  }

  PDM_l_num_t n_tri;
  double tri_coord[9];
  PDM_l_num_t _tri_vtx[3];
  PDM_l_num_t *tri_vtx = malloc (sizeof(PDM_l_num_t) * (n_vtx_face_max - 2)*3);
  PDM_triangulate_state_t *state = PDM_triangulate_state_create (n_vtx_face_max);

  int *closest_face = malloc (sizeof(int) * n_pts);
  double *closest_point = malloc (sizeof(double) * n_pts * 3);

  /* Loop on faces */
  for (int iface = 0; iface < n_face; iface++) {

    const PDM_l_num_t *_face_vtx = face_vtx + face_vtx_idx[iface];
    n_vtx_face = face_vtx_idx[iface+1] - face_vtx_idx[iface];

    /* Triangular face */
    if (n_vtx_face == 3) {
      n_tri = 1;
      for (int ivtx = 0; ivtx < 3; ivtx++) {
        tri_vtx[ivtx] = _face_vtx[ivtx];
      }
    }

    /* Quadrilateral face */
    else if (n_vtx_face == 4) {
      n_tri = PDM_triangulate_quadrangle (3,
                                          vtx_coord,
                                          NULL,
                                          _face_vtx,
                                          tri_vtx);
    }

    /* Polygonal face */
    else {
      n_tri = PDM_triangulate_polygon(3,
                                      n_vtx_face,
                                      vtx_coord,
                                      NULL,
                                      _face_vtx,
                                      PDM_TRIANGULATE_MESH_DEF,
                                      tri_vtx,
                                      state);
    }

    /* Loop on face triangles so as to loop on tetrahedra
       built by joining face triangles and pseudo-center */
    for (int itri = 0; itri < n_tri; itri++) {

      for (int ivtx = 0; ivtx < 3; ivtx++) {
        _tri_vtx[ivtx] = tri_vtx[3*itri + ivtx] - 1;
      }

      if (face_orientation[iface] < 0) {
        PDM_l_num_t tmp = _tri_vtx[0];
        _tri_vtx[0] = _tri_vtx[2];
        _tri_vtx[2] = tmp;
      }

      for (int ivtx = 0; ivtx < 3; ivtx++) {
        for (int idim = 0; idim < 3; idim++) {
          tri_coord[3*ivtx + idim] = vtx_coord[3*_tri_vtx[ivtx] + idim];
        }
      }

      /* Loop on points (compute solid angle and min_dist) */
      for (int ipt = 0; ipt < n_pts; ipt++) {

        if (distance[ipt] < eps2) {
          continue;
        }

        const double *_pt = pts_coord + 3*ipt;
        double *_cp = closest_point + 3*ipt;

        double min_dist2, closest[3];
        PDM_triangle_status_t error = PDM_triangle_evaluate_position (_pt,
                                                                      tri_coord,
                                                                      closest,
                                                                      &min_dist2,
                                                                      NULL);

        if (error == PDM_TRIANGLE_DEGENERATED) {
          continue;
        }

        if (distance[ipt] > min_dist2) {
          distance[ipt] = (float) min_dist2;

          closest_face[ipt] = iface;
          for (int idim = 0; idim < 3; idim++) {
            _cp[idim] = closest[idim];
          }

        }

        /* Solid angle */
        if (distance[ipt] < eps2) {
          continue;
        }

        double v[3][3], lv[3];
        double denom = 1.;
        for (int i = 0; i < 3; i++) {
          for (int idim = 0; idim < 3; idim++) {
            v[i][idim] = tri_coord[3*i + idim] - _pt[idim];
          }
          lv[i] = PDM_MODULE (v[i]);
          denom *= lv[i];
        }
        double det_abc = _determinant_3x3 (v[0], v[1], v[2]);

        denom += lv[0] * PDM_DOT_PRODUCT (v[1], v[2]);
        denom += lv[1] * PDM_DOT_PRODUCT (v[2], v[0]);
        denom += lv[2] * PDM_DOT_PRODUCT (v[0], v[1]);

        double half_angle = atan2 (det_abc, denom);

        if ((half_angle < 0.) && (det_abc > 0)) {
          half_angle = 2.*PDM_PI - half_angle;
        }
        else if ((half_angle > 0.) && (det_abc < 0)){
          half_angle = -2.*PDM_PI + half_angle;
        }

        if (0) {
          printf("face %d, tri %d, pt %d : solid angle = %f * PI\n",
                 iface, itri, ipt, 2. * half_angle / PDM_PI);
        }
        solid_angle[ipt] += 2. * half_angle;

      } // End of loop on points

    } // End of loop on triangles

  } // End of loop on faces

  state = PDM_triangulate_state_destroy (state);
  free (tri_vtx);


  /*
   * Locate points (compute mean value coordinates)
   */
  double *face_coord = malloc (sizeof(double) * n_vtx_face_max * 3);
  double *bar_coord_face = malloc (sizeof(double) * n_vtx_face_max);

  int stat = 1;

  for (int ipt = 0; ipt < n_pts; ipt++) {

    if (solid_angle[ipt] > eps_solid_angle &&
        solid_angle[ipt] < four_PI - eps_solid_angle &&
        distance[ipt] > eps_distance2) {
      /* Non-closed polyhedron */
      const double *_pt = pts_coord + 3 * ipt;
      printf("!! pt %d (%f, %f, %f) solid_angle/PI = %g, dist = %g\n",
             ipt, _pt[0], _pt[1], _pt[2], solid_angle[ipt]/PDM_PI, distance[ipt]);
      stat = 0;
      continue;
    }

    const double *_pt = pts_coord + 3 * ipt;
    double       *_bc = bar_coord + n_vtx * ipt;

    /* Point outside polyhedron or very close from its boundary */
    if (solid_angle[ipt] < eps_solid_angle || distance[ipt] < eps_distance2) {
      /* Compute mean value coords in closest face */
      for (int ivtx = 0; ivtx < n_vtx; ivtx++) {
        _bc[ivtx] = 0.;
      }

      int iface = closest_face[ipt];


      n_vtx_face = face_vtx_idx[iface+1] - face_vtx_idx[iface];

      for (int ivtx = 0; ivtx < n_vtx_face; ivtx++) {
        int _ivtx = face_vtx[face_vtx_idx[iface] + ivtx] - 1;

        for (int idim = 0; idim < 3; idim++) {
          face_coord[3*ivtx + idim] = vtx_coord[3*_ivtx + idim];
        }
      }

      PDM_mean_value_coordinates_polygon_3d (n_vtx_face,
                                             face_coord,
                                             1,
                                             closest_point + 3*ipt,
                                             bar_coord_face);

      for (int ivtx = 0; ivtx < n_vtx_face; ivtx++) {
        int _ivtx = face_vtx[face_vtx_idx[iface] + ivtx] - 1;
        _bc[_ivtx] = bar_coord_face[ivtx];
      }

    }

    /* Point strictly inside polyhedron */
    else {

      PDM_mean_values_polyhedron (n_vtx,
                                  vtx_coord,
                                  n_face,
                                  face_vtx_idx,
                                  face_vtx,
                                  face_orientation,
                                  _pt,
                                  _bc);

      distance[ipt] = -distance[ipt];
    }


  } // End of loop on points

  free (solid_angle);
  free (closest_face);
  free (closest_point);
  free (bar_coord_face);
  free (face_coord);

  return stat;
}

/*============================================================================
 * Public function definitions
 *============================================================================*/
void
PDM_point_location_nodal
(
 const int           type_idx[],
 const int           elt_vtx_idx[],
 const double        elt_vtx_coord[],
 const PDM_l_num_t   poly3d_face_idx[],
 const PDM_l_num_t   face_vtx_idx[],
 const PDM_l_num_t   face_vtx[],
 const int           face_orientation[],
 const int           pts_idx[],
 const double        pts_coord[],
 const double        tolerance,
 float             **distance,
 double            **projected_coord,
 int               **bar_coord_idx,
 double            **bar_coord
 )
{
  /*
   * TO DO: add support for high-order elements
   */
  //const int order = 1;


  int ielt = 0;
  int ipt = 0;


  /*
   * Allocate arrays
   */
  int n_elt = type_idx[PDM_MESH_NODAL_N_ELEMENT_TYPES];
  int n_pts = pts_idx[n_elt];
  *distance = malloc (sizeof(float) * n_pts);
  *projected_coord = malloc (sizeof(double) * n_pts * 3);
  *bar_coord_idx = malloc (sizeof(int) * (n_pts+1));
  (*bar_coord_idx)[0] = 0;

  for (ielt = 0; ielt < n_elt; ielt++) {
    int n_vtx = elt_vtx_idx[ielt+1] - elt_vtx_idx[ielt];

    for (ipt = pts_idx[ielt]; ipt < pts_idx[ielt+1]; ipt++) {
      (*bar_coord_idx)[ipt+1] = (*bar_coord_idx)[ipt] + n_vtx;
    }
  }

  *bar_coord = malloc (sizeof(double) * ((*bar_coord_idx)[n_pts]));


  /********************************
   * Perform elementary locations *
   ********************************/

  /*
   * Single nodes
   */
  for (ielt = type_idx[0]; ielt < type_idx[PDM_MESH_NODAL_BAR2]; ielt++) {
    const double *vtx_coord = elt_vtx_coord + 3*elt_vtx_idx[ielt];

    for (ipt = pts_idx[ielt]; ipt < pts_idx[ielt+1]; ipt++) {
      double dist2 = 0.;

      for (int idim = 0; idim < 3; idim++) {
        double delta = pts_coord[3*ipt + idim] - vtx_coord[idim];
        dist2 += delta * delta;
      }
      (*distance)[ipt] = (float) dist2;
      (*bar_coord)[ipt] = 1.;
    }
  }


  /*
   * Edges
   */
  for (ielt = type_idx[PDM_MESH_NODAL_BAR2]; ielt < type_idx[PDM_MESH_NODAL_TRIA3]; ielt++) {

    _locate_on_edge (elt_vtx_coord + elt_vtx_idx[ielt] * 3,
                     pts_idx[ielt+1] - pts_idx[ielt],
                     pts_coord + pts_idx[ielt] * 3,
                     *distance + pts_idx[ielt],
                     *bar_coord + (*bar_coord_idx)[pts_idx[ielt]]);
  }


  /*
   * Triangles
   */
  const int _connec[3] = {1, 2, 3};
  for (ielt = type_idx[PDM_MESH_NODAL_TRIA3]; ielt < type_idx[PDM_MESH_NODAL_QUAD4]; ielt++) {

    _locate_on_triangles (1,
                          _connec,
                          elt_vtx_coord + elt_vtx_idx[ielt] * 3,
                          pts_idx[ielt+1] - pts_idx[ielt],
                          pts_coord + pts_idx[ielt] * 3,
                          NULL,
                          *distance + pts_idx[ielt],
                          *bar_coord + (*bar_coord_idx)[pts_idx[ielt]]);
  }


  /*
   * Quadrangles
   */
  for (ielt = type_idx[PDM_MESH_NODAL_QUAD4]; ielt < type_idx[PDM_MESH_NODAL_POLY_2D]; ielt++) {

    _locate_on_quadrangle (elt_vtx_coord + elt_vtx_idx[ielt] * 3,
                           pts_idx[ielt+1] - pts_idx[ielt],
                           pts_coord + pts_idx[ielt] * 3,
                           *distance + pts_idx[ielt],
                           *bar_coord + (*bar_coord_idx)[pts_idx[ielt]]);
  }


  /*
   * Polygons
   */
  for (ielt = type_idx[PDM_MESH_NODAL_POLY_2D]; ielt < type_idx[PDM_MESH_NODAL_TETRA4]; ielt++) {

    _locate_in_polygon (elt_vtx_idx[ielt+1] - elt_vtx_idx[ielt],
                        elt_vtx_coord + elt_vtx_idx[ielt] * 3,
                        pts_idx[ielt+1] - pts_idx[ielt],
                        pts_coord + pts_idx[ielt] * 3,
                        *distance + pts_idx[ielt],
                        *bar_coord + (*bar_coord_idx)[pts_idx[ielt]]);
  }



  /*
   * Tetrahedra
   */
  for (ielt = type_idx[PDM_MESH_NODAL_TETRA4]; ielt < type_idx[PDM_MESH_NODAL_PYRAMID5]; ielt++) {

    _locate_in_tetrahedron (elt_vtx_coord + elt_vtx_idx[ielt] * 3,
                            pts_idx[ielt+1] - pts_idx[ielt],
                            pts_coord + pts_idx[ielt] * 3,
                            *distance + pts_idx[ielt],
                            *bar_coord + (*bar_coord_idx)[pts_idx[ielt]]);
  }


  /*
   * Other standard cells
   */
  for (PDM_Mesh_nodal_elt_t type = PDM_MESH_NODAL_PYRAMID5;
       type < PDM_MESH_NODAL_POLY_3D;
       type++) {

    for (ielt = type_idx[type]; ielt < type_idx[type+1]; ielt++) {

      _locate_in_cell_3d (type,
                          elt_vtx_coord + elt_vtx_idx[ielt] * 3,
                          pts_idx[ielt+1] - pts_idx[ielt],
                          pts_coord + pts_idx[ielt] * 3,
                          tolerance,
                          *distance + pts_idx[ielt],
                          *bar_coord + (*bar_coord_idx)[pts_idx[ielt]]);
    }
  }


  /*
   * Polyhedra
   */
  int ipoly = 0;
  for (ielt = type_idx[PDM_MESH_NODAL_POLY_3D]; ielt < n_elt; ielt++) {

    int n_vtx = elt_vtx_idx[ielt+1] - elt_vtx_idx[ielt];
    _locate_in_polyhedron (n_vtx,
                           elt_vtx_coord + elt_vtx_idx[ielt] * 3,
                           poly3d_face_idx[ipoly+1] - poly3d_face_idx[ipoly],
                           face_vtx_idx + poly3d_face_idx[ipoly],
                           face_vtx,
                           face_orientation + poly3d_face_idx[ipoly],
                           pts_idx[ielt+1] - pts_idx[ielt],
                           pts_coord + pts_idx[ielt] * 3,
                           //                           poly3d_char_length[ipoly],
                           *distance + pts_idx[ielt],
                           *bar_coord + (*bar_coord_idx)[pts_idx[ielt]]);

    ipoly++;
  }

}

#ifdef __cplusplus
}
#endif /* __cplusplus */
