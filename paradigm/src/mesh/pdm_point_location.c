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
#include "pdm_logging.h"
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
#include "pdm_ho_basis.h"
#include "pdm_ho_ordering.h"
#include "pdm_part_mesh_nodal_elmts.h"
#include "pdm_array.h"

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
  double u = uvw[0];
  double v = uvw[1];
  double w = uvw[2];
  double u1 = 1. - u;
  double v1 = 1. - v;
  double w1 = 1. - w;

  switch (elt_type) {

  case PDM_MESH_NODAL_PYRAMID5: {

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

    shapef[0] = (1. - u - v) * w1;
    shapef[1] = u * w1;
    shapef[2] = v * w1;
    shapef[3] = (1. - u - v) * w;
    shapef[4] = u * w;
    shapef[5] = v * w;

    if (deriv != NULL) {
      deriv[0][0] = -w1;
      deriv[0][1] = -w1;
      deriv[0][2] = -(1. - u - v);
      deriv[1][0] =  w1;
      deriv[1][1] =  0.;
      deriv[1][2] = -u;
      deriv[2][0] =  0.;
      deriv[2][1] =  w1;
      deriv[2][2] = -v;
      deriv[3][0] = -w;
      deriv[3][1] = -w;
      deriv[3][2] =  (1. - u - v);
      deriv[4][0] =  w;
      deriv[4][1] =  0.;
      deriv[4][2] =  u;
      deriv[5][0] =  0.;
      deriv[5][1] =  w;
      deriv[5][2] =  v;
    }

    break;
  }



  case PDM_MESH_NODAL_HEXA8: {

    shapef[0] = u1 * v1 * w1;
    shapef[1] = u  * v1 * w1;
    shapef[2] = u  * v  * w1;
    shapef[3] = u1 * v  * w1;
    shapef[4] = u1 * v1 * w;
    shapef[5] = u  * v1 * w;
    shapef[6] = u  * v  * w;
    shapef[7] = u1 * v  * w;

    if (deriv != NULL) {
      deriv[0][0] = -v1 * w1;
      deriv[0][1] = -u1 * w1;
      deriv[0][2] = -u1 * v1;
      deriv[1][0] =  v1 * w1;
      deriv[1][1] = -u  * w1;
      deriv[1][2] = -u  * v1;
      deriv[2][0] =  v  * w1;
      deriv[2][1] =  u  * w1;
      deriv[2][2] = -u  * v;
      deriv[3][0] = -v  * w1;
      deriv[3][1] =  u1 * w1;
      deriv[3][2] = -u1 * v;
      deriv[4][0] = -v1 * w;
      deriv[4][1] = -u1 * w;
      deriv[4][2] =  u1 * v1;
      deriv[5][0] =  v1 * w;
      deriv[5][1] = -u  * w;
      deriv[5][2] =  u  * v1;
      deriv[6][0] =  v  * w;
      deriv[6][1] =  u  * w;
      deriv[6][2] =  u  * v;
      deriv[7][0] = -v  * w;
      deriv[7][1] =  u1 * w;
      deriv[7][2] =  u1 * v;
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
 *   proj_coord         --> coordinates of closest points (size: 3*n_pts)
 *   bar_coord          --> barycentric coordinates of closest points (size: 2*n_pts)
 *----------------------------------------------------------------------------*/
static void
_locate_on_edge
(
 const double       vtx_coord[],
 const int          n_pts,
 const double       pts_coord[],
 double             distance[],
 double             proj_coord[],
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
    double       *_cp = proj_coord + 3 * ipt;

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
        _cp[idim] = vtx_coord[3 + idim];
      }
    }

    else if (t > 0.) {
      for (idim = 0; idim < 3; idim++) {
        v[idim] -= t * u[idim];
        _cp[idim] = vtx_coord[idim] + v[idim];
      }
    }

    else {
      for (idim = 0; idim < 3; idim++) {
        _cp[idim] = vtx_coord[idim];
      }
    }

    /* Distance between point to locate and its projection */
    distance[ipt] = 0.;
    for (idim = 0; idim < 3; idim++) {
      distance[ipt] += v[idim] * v[idim];
    }

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
 *   proj_coord          --> coordinates of closest points (size: 3*n_pts)
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
 double             distance[],
 double             proj_coord[],
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
      double       *_cp = proj_coord + 3 * ipt;

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

        if (proj_coord != NULL) {
          _cp[0] = closest_point[0];
          _cp[1] = closest_point[1];
          _cp[2] = closest_point[2];
        }

        distance[ipt] = dist2;

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
 *   proj_coord          --> coordinates of closest points (size: 3*n_pts)
 *   bar_coord           <-> barcyentric coordinates of closest points (size: n_pts * 4)
 *----------------------------------------------------------------------------*/
static void
_locate_on_quadrangle
(
 const double       quad_coord[],
 const int          n_pts,
 const double       pts_coord[],
 double             distance[],
 double             proj_coord[],
 double             bar_coord[]
 )
{
  int ipt, ivtx, idim;


  PDM_mean_values_polygon_3d (4,
                              quad_coord,
                              n_pts,
                              pts_coord,
                              bar_coord);


  for (ipt = 0; ipt < n_pts; ipt++) {

    const double *_pt = pts_coord + 3 * ipt;
    double       *_bc = bar_coord + 4 * ipt;
    double       *_cp = proj_coord + 3 * ipt;

    /*double v_cp_p[3];
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

    distance[ipt] = dist2;*/
    for (idim = 0; idim < 3; idim++) {
      _cp[idim] = 0.;
    }

    for (ivtx = 0; ivtx < 4; ivtx++) {
      for (idim = 0; idim < 3; idim++) {
        _cp[idim] += _bc[ivtx] * quad_coord[3*ivtx + idim];
      }
    }

    double v_cp_p[3] = {_pt[0] - _cp[0], _pt[1] - _cp[1], _pt[2] - _cp[2]};
    distance[ipt] = PDM_DOT_PRODUCT (v_cp_p, v_cp_p);
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
 *   proj_coord          --> coordinates of closest points (size: 3*n_pts)
 *   bar_coord           <-> barcyentric coordinates of closest points (size: n_pts * 4)
 *----------------------------------------------------------------------------*/
static void
_locate_in_tetrahedron
(
 const double        tetra_coord[12],
 const int           n_pts,
 const double        pts_coord[],
 double             *distance,
 double             *proj_coord,
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
    double       *_cp = proj_coord + 3 * ipt;

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

    distance[ipt] = min_bc * min_bc;
    if (min_bc > 0.) {
      distance[ipt] = -distance[ipt];
    }

    /* Point outside tetrahedron */
    if (distance[ipt] > 0.) {
      pts_out[n_pts_out++] = ipt;
    }

    /* Point inside tetrahedron */
    else {
      for (i = 0; i < 3; i++) {
        _cp[i] = _pt[i];
      }
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
  double *closest_point_face = malloc (sizeof(double) * n_pts_out * 3);
  double *distance_face = malloc (sizeof(double) * n_pts_out);
  _locate_on_triangles (4,
                        face_vtx,
                        tetra_coord,
                        n_pts_out,
                        pts_out_coord,
                        id_face,
                        distance_face,
                        closest_point_face,
                        bar_coord_face);

  for (ipt = 0; ipt < n_pts_out; ipt++) {
    int id_pt = pts_out[ipt];
    double *_bc = bar_coord + 4 * id_pt;
    double *_cp = proj_coord + 3 * id_pt;

    for (ivtx = 0; ivtx < 4; ivtx++) {
      _bc[ivtx] = 0.;
    }

    for (ivtx = 0; ivtx < 3; ivtx++) {
      int id_vtx = face_vtx[face_vtx_idx[id_face[ipt]] + ivtx] - 1;
      _bc[id_vtx] = bar_coord_face[3*ipt + ivtx];
    }

    for (i = 0; i < 3; i++) {
      _cp[i] = closest_point_face[3*ipt + i];
    }

    distance[id_pt] = distance_face[ipt];
  }

  free (pts_out);
  free (pts_out_coord);
  free (id_face);
  free (bar_coord_face);
  free (distance_face);
  free (closest_point_face);
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
 *   proj_coord          --> coordinates of closest points (size: 3*n_pts)
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
 double                     *distance,
 double                     *proj_coord,
 double                     *bar_coord,
 double                     *uvw
 )
{
  double _uvw[3];
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
                            proj_coord,
                            bar_coord);
    // TODO: compute uvw...
    if (uvw != NULL) {
      for (int ipt = 0; ipt < n_pts; ipt++) {
        uvw[3*ipt  ] = -1.;
        uvw[3*ipt+1] = -1.;
        uvw[3*ipt+2] = -1.;
      }
    }
    return;
  }

  double *_cell_coord = NULL;
  if (elt_type == PDM_MESH_NODAL_PRISM6) {
    _cell_coord = (double *) cell_coord;
  }
  else {

    //FIXME: Allocation a chaque passage. Prevoir buffer en argument alloue a 8 * 3

    _cell_coord = malloc (sizeof(double) * n_vtx * 3);
    _cell_coord[ 0] = cell_coord[ 0];
    _cell_coord[ 1] = cell_coord[ 1];
    _cell_coord[ 2] = cell_coord[ 2];

    _cell_coord[ 3] = cell_coord[ 3];
    _cell_coord[ 4] = cell_coord[ 4];
    _cell_coord[ 5] = cell_coord[ 5];

    _cell_coord[ 6] = cell_coord[ 9];
    _cell_coord[ 7] = cell_coord[10];
    _cell_coord[ 8] = cell_coord[11];

    _cell_coord[ 9] = cell_coord[ 6];
    _cell_coord[10] = cell_coord[ 7];
    _cell_coord[11] = cell_coord[ 8];


    _cell_coord[12] = cell_coord[12];
    _cell_coord[13] = cell_coord[13];
    _cell_coord[14] = cell_coord[14];

    if (elt_type == PDM_MESH_NODAL_HEXA8) {
      _cell_coord[15] = cell_coord[15];
      _cell_coord[16] = cell_coord[16];
      _cell_coord[17] = cell_coord[17];

      _cell_coord[18] = cell_coord[21];
      _cell_coord[19] = cell_coord[22];
      _cell_coord[20] = cell_coord[23];

      _cell_coord[21] = cell_coord[18];
      _cell_coord[22] = cell_coord[19];
      _cell_coord[23] = cell_coord[20];
    }
  }


  /* Other cell types, shape functions must be computed iteratively */
  for (int ipt = 0; ipt < n_pts; ipt++) {

    const double *_pt = pts_coord + 3 * ipt;
    double *_bc = bar_coord + n_vtx * ipt;
    double *_cp = proj_coord + 3 * ipt;

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

        for (int i = 0; i < 3; i++) {
          _cp[i] = cell_coord[3*ivtx + i];
        }

        on_vtx = PDM_TRUE;
        break;
      }

    } // End loop on vertices

    if (on_vtx == PDM_TRUE) {
      continue;
    }

    /* Compute parametric coordinates : Try Newton */
    PDM_bool_t stat_uvw = _compute_uvw (elt_type,
                                        _pt,
                                        cell_coord,
                                        tolerance,
                                        _uvw);
    if (stat_uvw == PDM_TRUE) {
      _compute_shapef_3d (elt_type, _uvw, _bc, NULL);

      double min_bc = HUGE_VAL;
      for (int ivtx = 0; ivtx < n_vtx; ivtx++) {
        min_bc = PDM_MIN (min_bc, _bc[ivtx]);
      }

      distance[ipt] = min_bc * min_bc;
      if (min_bc > 0.) {
        distance[ipt] = -distance[ipt];
      }
    }

    /* Failed to compute parametric coordinates */
    else {
      // Retrieve ho type
      PDM_Mesh_nodal_elt_t t_elt_ho;
      switch (elt_type) {
        case PDM_MESH_NODAL_BAR2:
          t_elt_ho = PDM_MESH_NODAL_BARHO;
          break;

        case PDM_MESH_NODAL_TRIA3:
          t_elt_ho = PDM_MESH_NODAL_TRIAHO;
          break;

        case PDM_MESH_NODAL_QUAD4:
          t_elt_ho = PDM_MESH_NODAL_QUADHO;
          break;

        case PDM_MESH_NODAL_TETRA4:
          t_elt_ho = PDM_MESH_NODAL_TETRAHO;
          break;

        case PDM_MESH_NODAL_PYRAMID5:
          t_elt_ho = PDM_MESH_NODAL_PYRAMIDHO;
          break;

        case PDM_MESH_NODAL_PRISM6:
          t_elt_ho = PDM_MESH_NODAL_PRISMHO;
          break;

        case PDM_MESH_NODAL_HEXA8:
          t_elt_ho = PDM_MESH_NODAL_HEXAHO;
          break;

        default:
          t_elt_ho = elt_type;
      }
      /*
       * Use hierarchical subdivision
       */
      double _proj_coord[3];
      double dist2 = PDM_ho_location (t_elt_ho,
                                      order,
                                      n_vtx,
                                      _cell_coord,
                                      _pt,
                                      _proj_coord,
                                      _uvw);
      /* Point inside */
      if (dist2 < 1.e-12) {
        _compute_shapef_3d (elt_type, _uvw, _bc, NULL);

        double min_bc = HUGE_VAL;
        for (int ivtx = 0; ivtx < n_vtx; ivtx++) {
          min_bc = PDM_MIN (min_bc, _bc[ivtx]);
        }

        distance[ipt] = -min_bc * min_bc;
      }

      /* Point outside */
      else {
        distance[ipt] = dist2;
      }
    }

    /* Point outside cell */
    if (distance[ipt] > 0.) {
      pts_out[n_pts_out++] = ipt;
      distance[ipt] = HUGE_VAL;
    }

     /* Point inside cell */
    else {
      for (int i = 0; i < 3; i++) {
        _cp[i] = _pt[i];
      }
    }

    if (uvw != NULL) {
      memcpy(uvw + 3*ipt, _uvw, sizeof(double) * 3);
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
            distance[_ipt] = min_dist2;

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

      PDM_mean_values_polygon_3d (n_vtx_face,
                                  face_coord,
                                  1,
                                  _cp,
                                  bar_coord_face);

      for (int ivtx = 0; ivtx < n_vtx_face; ivtx++) {
        int _ivtx = face_vtx[face_vtx_idx[iface] + ivtx] - 1;
        _bc[_ivtx] = bar_coord_face[ivtx];
      }

      const double *_pt = pts_coord + 3 * _ipt;
      double *_pp = proj_coord + 3 * _ipt;
      for (int idim = 0; idim < 3; idim++) {
        _pp[idim] = 0.;
      }

      for (int ivtx = 0; ivtx < n_vtx_face; ivtx++) {
        int _ivtx = face_vtx[face_vtx_idx[iface] + ivtx] - 1;
        for (int idim = 0; idim < 3; idim++) {
          _pp[idim] += bar_coord_face[ivtx] * cell_coord[3*_ivtx + idim];
        }
      }

      double v_p_cp[3] = {_pt[0] - _pp[0], _pt[1] - _pp[1], _pt[2] - _pp[2]};
      distance[_ipt] = PDM_DOT_PRODUCT (v_p_cp, v_p_cp);
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
 *   proj_coord          --> coordinates of closest points (size: 3*n_pts)
 *   bar_coord           <-> barcyentric coordinates of closest points (size: n_pts * n_vtx)
 *----------------------------------------------------------------------------*/
static void
_locate_in_polygon
(
 const PDM_l_num_t n_vtx,
 const double      vtx_coord[],
 const int         n_pts,
 const double      pts_coord[],
 double            distance[],
 double            proj_coord[],
 double            bar_coord[]
 )
{
  /* Compute mean value coordinates of closest points on polygon */
  PDM_mean_values_polygon_3d (n_vtx,
                              vtx_coord,
                              n_pts,
                              pts_coord,
                              bar_coord);

  /* Compute distances */
  for (int ipt = 0; ipt < n_pts; ipt++) {
    const double *_pt = pts_coord  + ipt * 3;
    double       *_bc = bar_coord  + ipt * n_vtx;
    double       *_cp = proj_coord + ipt * 3;

    for (int idim = 0; idim < 3; idim++) {
      _cp[idim] = 0.;
    }

    for (int ivtx = 0; ivtx < n_vtx; ivtx++) {
      for (int idim = 0; idim < 3; idim++) {
        _cp[idim] += _bc[ivtx] * vtx_coord[3*ivtx + idim];
      }
    }

    double v_cp_p[3] = {_pt[0] - _cp[0], _pt[1] - _cp[1], _pt[2] - _cp[2]};
    distance[ipt] = PDM_DOT_PRODUCT (v_cp_p, v_cp_p);
    int dbg = 0;
    // int dbg = (_pt[0] > -0.424127 && _pt[0] < -0.424125 &&
    //            _pt[1] > -0.029731 && _pt[1] < -0.029729);
    if (dbg) {
      log_trace("distance = %e\n", distance[ipt]);
    }

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
 *   distance            <-> distance from point to element (size: n_pts)
 *   proj_coord          --> coordinates of closest points (size: 3*n_pts)
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
 double            distance[],
 double            proj_coord[],
 double            bar_coord[]
 )
{
  const double four_PI = 4. * PDM_PI;

  const double eps_solid_angle = 1e-5;
  const double eps_distance2 = 1e-10;
  const double eps2 = 1e-20;

  int dbg = 0;
  // for (int i = 0; i < n_pts; i++) {
  //   const double *p = pts_coord + 3*i;
  //   if (p[0] > 0.36666 && p[0] < 0.36667 &&
  //       p[1] > 0.49999 && p[1] < 0.50001 &&
  //       p[2] > 0.83332 && p[2] < 0.83334) {
  //     dbg = 1;
  //     break;
  //   }
  // }

  // Visu for debugging
  if (0) {
    log_trace("# vtk DataFile Version 2.0\n");
    log_trace("mesh\n");
    log_trace("ASCII\n");
    log_trace("DATASET POLYDATA\n");

    log_trace("POINTS %d double\n", n_vtx);
    for (int i = 0; i < n_vtx; i++) {
      for (int j = 0; j < 3; j++) {
        log_trace("%.20lf ", vtx_coord[3*i+j]);
      }
      log_trace("\n");
    }


    log_trace("POLYGONS %d %d\n", n_face, n_face + face_vtx_idx[n_face] - face_vtx_idx[0]);
    for (int i = 0; i < n_face; i++) {
      log_trace("%d", face_vtx_idx[i+1] - face_vtx_idx[i]);
      if (face_orientation[i] < 0) {
        for (int j = face_vtx_idx[i+1]-1; j >= face_vtx_idx[i]; j--) {
          log_trace(" %d", face_vtx[j] - 1);
        }
      } else {
        for (int j = face_vtx_idx[i]; j < face_vtx_idx[i+1]; j++) {
          log_trace(" %d", face_vtx[j] - 1);
        }
      }
      log_trace("\n");
    }
    log_trace("-------\n");
  }

  /*
   * Identify points inside/outside polyhedron
   */
  double *solid_angle = malloc (sizeof(double) * n_pts);
  for (int ipt = 0; ipt < n_pts; ipt++) {
    solid_angle[ipt] = 0.;
    distance[ipt] = HUGE_VAL;
  }

  /* Count max nb of vertices per face */
  PDM_l_num_t n_vtx_face, n_vtx_face_max = 3;
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
          distance[ipt] = min_dist2;

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

    const double *p = pts_coord + 3 * ipt;
    // dbg = (p[0] > 0.36666 && p[0] < 0.36667 &&
    //        p[1] > 0.49999 && p[1] < 0.50001 &&
    //        p[2] > 0.83332 && p[2] < 0.83334);
    if (dbg) {
      log_trace("p %f %f %f : solid_angle = %f*PI, dist = %e\n",
                p[0], p[1], p[2], solid_angle[ipt]/PDM_PI, distance[ipt]);
    }

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
    double       *_cp = proj_coord + 3 * ipt;

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

      PDM_mean_values_polygon_3d (n_vtx_face,
                                  face_coord,
                                  1,
                                  closest_point + 3*ipt,
                                  bar_coord_face);

      for (int idim = 0; idim < 3; idim++) {
        _cp[idim] = 0.;
      }

      for (int ivtx = 0; ivtx < n_vtx_face; ivtx++) {
        int _ivtx = face_vtx[face_vtx_idx[iface] + ivtx] - 1;
        _bc[_ivtx] = bar_coord_face[ivtx];
        for (int idim = 0; idim < 3; idim++) {
          _cp[idim] += bar_coord_face[ivtx] * face_coord[3*ivtx + idim];
        }
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

      for (int idim = 0; idim < 3; idim++) {
        _cp[idim] = _pt[idim];
      }

      // check
      if (1) {
        double pp[3] = {-_pt[0], -_pt[1], -_pt[2]};
        for (int ivtx = 0; ivtx < n_vtx; ivtx++) {
          for (int idim = 0; idim < 3; idim++) {
            pp[idim] += _bc[ivtx] * vtx_coord[3*ivtx + idim];
          }
        }

        double err = PDM_MODULE(pp);
        if (err > 1e-6) {
          log_trace("!!! _pt %20.16f %20.16f %20.16f\n",
                    _pt[0], _pt[1], _pt[2]);
          log_trace("!!! error = %e\n", err);
        }
      }

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





static void
_compute_mean_value_coord_polyhedron_point_inside
(
 const int                      n_vtx,
 const double                   vtx_coord[],
 const int                      n_face,
 const int                      face_vtx_idx[],
 const int                      face_vtx[],
 const int                      face_orientation[],
 const double                   pt_coord[],
       double                   weight[],
       PDM_triangulate_state_t *tri_state,
       int                     *tri_vtx
 )
 {
  double eps_in_triangle = 1e-12;

  /* Initialize weights */
  for (int i = 0; i < n_vtx; i++) {
    weight[i] = 0.;
  }

  /*
   * Compute distances d between p and the vertices
   * and unit vectors u from point to vertex
   */
  double d[n_vtx];
  double u[n_vtx][3];

  for (int ivtx = 0; ivtx < n_vtx; ivtx++) {
    d[ivtx] = 0.;
    for (int j = 0; j < 3; j++) {
      u[ivtx][j] = vtx_coord[3*ivtx+j] - pt_coord[j];
      d[ivtx] += u[ivtx][j]*u[ivtx][j];
    }

    d[ivtx] = sqrt(d[ivtx]);
    double invd = 1/d[ivtx];
    for (int j = 0; j < 3; j++) {
      u[ivtx][j] *= invd;
    }
  } // End of loop on vertices

  /* Loop on faces */
  for (int iface = 0; iface < n_face; iface++) {

    const int *_face_vtx = face_vtx + face_vtx_idx[iface];
    int face_vtx_n = face_vtx_idx[iface+1] - face_vtx_idx[iface];

    /* Triangulate face if necessary */
    int n_tri;
    if (face_vtx_n == 3) {
      /* Triangular face */
      n_tri = 1;
      memcpy(tri_vtx, _face_vtx, sizeof(int) * 3);
    }
    else if (face_vtx_n == 4) {
      /* Quadrilateral face */
      n_tri = PDM_triangulate_quadrangle(3,
                                         vtx_coord,
                                         NULL,
                                         _face_vtx,
                                         tri_vtx);
    }
    else {
      /* Polygonal face */
      n_tri = PDM_triangulate_polygon(3,
                                      face_vtx_n,
                                      vtx_coord,
                                      NULL,
                                      _face_vtx,
                                      PDM_TRIANGULATE_MESH_DEF,
                                      tri_vtx,
                                      tri_state);
    }


    /* Loop on triangles */
    double *tri_u[3];
    // double tri_eps2[3];
    for (int itri = 0; itri < n_tri; itri++) {

      int __tri_vtx[3];
      if (face_orientation[iface] > 0) {
        for (int ivtx = 0; ivtx < 3; ivtx++) {
          __tri_vtx[ivtx] = tri_vtx[3*itri + ivtx] - 1;
        }
      }
      else {
        for (int ivtx = 0; ivtx < 3; ivtx++) {
          __tri_vtx[ivtx] = tri_vtx[3*itri + 2-ivtx] - 1;
        }
      }

      for (int ivtx = 0; ivtx < 3; ivtx++) {
        tri_u   [ivtx] = u[__tri_vtx[ivtx]];
        // tri_eps2[ivtx] = vtx_eps2[__tri_vtx[ivtx]];
      }

      double l[3];
      double theta[3];
      double sin_theta[3];
      double h = 0.;
      for (int ivtx = 0; ivtx < 3; ivtx++) {
        int ip = (ivtx+1)%3;
        int im = (ivtx+2)%3;
        l[ivtx] = 0.;
        for (int j = 0; j < 3; j++) {
          double delta = tri_u[ip][j] - tri_u[im][j];
          l[ivtx] += delta*delta;
        }

        // assert(l[ivtx] > tri_eps2[ip] && l[ivtx] > tri_eps2[im]);

        l[ivtx] = sqrt(l[ivtx]);

        theta[ivtx] = asin(PDM_MIN(1., 0.5*l[ivtx]));
        h += theta[ivtx];
        theta[ivtx] *= 2;

        sin_theta[ivtx] = sin(theta[ivtx]);
      }

      /* Check if p is coplanar with current triangle to avoid division by zero */
      int coplanar = (PDM_PI - h < eps_in_triangle);

      assert(!coplanar);

      double sdet = PDM_SIGN(_determinant_3x3(tri_u[0],
                                              tri_u[1],
                                              tri_u[2]));

      double c[3];
      double s[3];
      double sin_h = sin(h);
      int coplanar_outside = 0;
      for (int ivtx = 0; ivtx < 3; ivtx++) {
        int ip = (ivtx+1)%3;
        int im = (ivtx+2)%3;
        c[ivtx] = -1 + (2 * sin_h * sin(h - theta[ivtx])) /
        (sin_theta[ip] * sin_theta[im]);

        s[ivtx] = sdet * sqrt(PDM_MAX(0., 1 - c[ivtx]*c[ivtx]));

        if (PDM_ABS(s[ivtx]) < eps_in_triangle) {
          coplanar_outside = 1;
          break;
        }
      }

      if (coplanar_outside) {
        /* p lies outside current triangle, on the same plane */
        continue;
      }

      for (int i = 0; i < 3; i++) {
        int ip = (i+1)%3;
        int im = (i+2)%3;
        weight[__tri_vtx[i]] += (theta[i] - c[ip]*theta[im] - c[im]*theta[ip]) /
        (d[__tri_vtx[i]] * sin_theta[ip] * s[im]);
      }

    } // End of loop on triangles

  } // End of loop on faces


  /* Normalize weights */
  double sum_w = 0.;
  for (int ivtx = 0; ivtx < n_vtx; ivtx++) {
    sum_w += weight[ivtx];
  }
  double isum_w = 1./sum_w;
  for (int ivtx = 0; ivtx < n_vtx; ivtx++) {
    weight[ivtx] *= isum_w;
  }
 }




static int
_locate_in_polyhedron2
(
 const int                      n_vtx,
 const double                   vtx_coord[],
 const int                      n_face,
 const int                      face_vtx_idx[],
 const int                      face_vtx[],
 const int                      face_orientation[],
 const int                      n_pts,
 const double                   pts_coord[],
       double                   distance[],
       double                   proj_coord[],
       double                   weight[],
       PDM_triangulate_state_t *tri_state,
       int                     *tri_vtx
 )
{
  int dbg = 0;

  double eps_on_face2 = 1.e-24;

  /* Prepare for triangulation */
  PDM_triangulate_state_t *_tri_state = tri_state;
  int                     *_tri_vtx   = tri_vtx;
  if (_tri_state == NULL || _tri_vtx == NULL) {
    int max_face_vtx_n = 0;
    for (int iface = 0; iface < n_face; iface++) {
      int face_vtx_n = face_vtx_idx[iface+1] - face_vtx_idx[iface];
      max_face_vtx_n = PDM_MAX(max_face_vtx_n, face_vtx_n);
    }

    if (_tri_state == NULL && max_face_vtx_n > 4) {
      _tri_state = PDM_triangulate_state_create(max_face_vtx_n);
    }

    if (_tri_vtx == NULL) {
      _tri_vtx = malloc(sizeof(int) * (max_face_vtx_n - 2) * 3);
    }
  }


  /*
   * First, identify points inside/outside polyhedron
   *  - compute sum of solid angles
   *  - compute distance from closest subtriangle
   */
  double solid_angle[n_pts];
  for (int ipt = 0; ipt < n_pts; ipt++) {
    solid_angle [ipt] = 0.;
    distance    [ipt] = HUGE_VAL;
  }

  /* Loop on faces */
  for (int iface = 0; iface < n_face; iface++) {

    const int *_face_vtx = face_vtx + face_vtx_idx[iface];
    int face_vtx_n = face_vtx_idx[iface+1] - face_vtx_idx[iface];

    /* Triangulate face if necessary */
    int n_tri;
    if (face_vtx_n == 3) {
      /* Triangular face */
      n_tri = 1;
      memcpy(_tri_vtx, _face_vtx, sizeof(int) * 3);
    }
    else if (face_vtx_n == 4) {
      /* Quadrilateral face */
      n_tri = PDM_triangulate_quadrangle(3,
                                         vtx_coord,
                                         NULL,
                                         _face_vtx,
                                         _tri_vtx);
    }
    else {
      /* Polygonal face */
      n_tri = PDM_triangulate_polygon(3,
                                      face_vtx_n,
                                      vtx_coord,
                                      NULL,
                                      _face_vtx,
                                      PDM_TRIANGULATE_MESH_DEF,
                                      _tri_vtx,
                                      _tri_state);
    }

    /* Loop on subtriangles */
    for (int itri = 0; itri < n_tri; itri++) {

      double tri_coord[9];
      if (face_orientation[iface] > 0) {
        for (int idx_vtx = 0; idx_vtx < 3; idx_vtx++) {
          int vtx_id = _tri_vtx[3*itri + idx_vtx] - 1;
          memcpy(tri_coord + 3*idx_vtx, vtx_coord + 3*vtx_id, sizeof(double) * 3);
        }
      } else {
        for (int idx_vtx = 0; idx_vtx < 3; idx_vtx++) {
          int vtx_id = _tri_vtx[3*itri + 2-idx_vtx] - 1;
          memcpy(tri_coord + 3*idx_vtx, vtx_coord + 3*vtx_id, sizeof(double) * 3);
        }
      }

      /* Loop on points */
      for (int ipt = 0; ipt < n_pts; ipt++) {

        if (distance[ipt] < eps_on_face2) {
          continue;
        }

        const double *pt_coord = pts_coord + 3*ipt;

        /* Distance */
        double tri_dist2, tri_closest_point[3], tri_weight[3];
        PDM_triangle_status_t stat = PDM_triangle_evaluate_position(pt_coord,
                                                                    tri_coord,
                                                                    tri_closest_point,
                                                                    &tri_dist2,
                                                                    tri_weight);

        if (stat == PDM_TRIANGLE_DEGENERATED) {
          /* Raise error? */
          continue;
        }

        if (tri_dist2 < distance[ipt]) {
          distance    [ipt] = tri_dist2;

          memcpy(proj_coord + 3*ipt, tri_closest_point, sizeof(double) * 3);

          double *w = weight + ipt*n_vtx;
          for (int i = 0; i < n_vtx; i++) {
            w[i] = 0;
          }

          for (int i = 0; i < 3; i++) {
            int vtx_id = _tri_vtx[3*itri + i] - 1;
            w[vtx_id] = tri_weight[i];
          }

          if (distance[ipt] < eps_on_face2) {
            continue;
          }
        }

        /* Solid angle */
        double v[3][3], lv[3];
        double denom = 1.;
        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            v[i][j] = tri_coord[3*i + j] - pt_coord[j];
          }
          lv[i] = PDM_MODULE(v[i]);
          denom *= lv[i];
        }
        for (int i = 0; i < 3; i++) {
          denom += lv[i] * PDM_DOT_PRODUCT(v[(i+1)%3], v[(i+2)%3]);
        }

        double numer = _determinant_3x3(v[0], v[1], v[2]);

        double half_angle = atan2(numer, denom);
        if (half_angle < 0. && numer > 0) {
          half_angle = 2.*PDM_PI - half_angle;
        }
        else if (half_angle > 0. && numer < 0){
          half_angle = -2.*PDM_PI + half_angle;
        }

        solid_angle[ipt] += 2*half_angle;

      } // End of loop on points

    } // End of loop on subtriangles

  } // End of loop on faces


  /*
   * Locate points (compute mean value coordinates)
   */
  double eps_solid_angle = 1.e-6;
  double threshold_inside = 4*PDM_PI*(1 - eps_solid_angle);

  for (int ipt = 0; ipt < n_pts; ipt++) {

    const double *p  = pts_coord  + ipt*3;
    double       *w  = weight     + ipt*n_vtx;
    double       *pc = proj_coord + ipt*3;

    if (dbg) {
      log_trace("point %d (%f %f %f), solid_angle = %f*PI, dist2 = %e\n",
                ipt, p[0], p[1], p[2], solid_angle[ipt]/PDM_PI, distance[ipt]);
    }


    if (solid_angle[ipt] > threshold_inside) {
      /* Point strictly inside polyhedron */
      distance[ipt] = -distance[ipt];

      // Compute mean value coordinates
      _compute_mean_value_coord_polyhedron_point_inside(n_vtx,
                                                        vtx_coord,
                                                        n_face,
                                                        face_vtx_idx,
                                                        face_vtx,
                                                        face_orientation,
                                                        p,
                                                        w,
                                                        _tri_state,
                                                        _tri_vtx);

      for (int j = 0; j < 3; j++) {
        pc[j] = 0.;
      }

      for (int ivtx = 0; ivtx < n_vtx; ivtx++) {
        const double *vc = vtx_coord + 3*ivtx;
        for (int j = 0; j < 3; j++) {
          pc[j] += w[ivtx] * vc[j];
        }
      }

      if (dbg) {
        log_trace("  interior\n");
        double e[3] = {pc[0] - p[0], pc[1] - p[1], pc[2] - p[2]};
        double err = PDM_DOT_PRODUCT(e, e);
        // PDM_log_trace_array_double(w, n_vtx, "  mean_value_coord : ");
        log_trace("  proj_coord = %f %f %f, dist2 = %e\n", pc[0], pc[1], pc[2], err);
      }

    }

    else if (distance[ipt] > eps_on_face2) {
      if (dbg) {
        /* Point strictly outside polyhedron */
        log_trace("  exterior\n");
        // PDM_log_trace_array_double(w, n_vtx, "  mean_value_coord : ");
        log_trace("  proj_coord = %f %f %f, dist2 = %e\n", pc[0], pc[1], pc[2], distance[ipt]);
      }
    }

    else {
      if (dbg) {
        /* Point on a polyhedron face */
        log_trace("  on face\n");
        // PDM_log_trace_array_double(w, n_vtx, "  mean_value_coord : ");
        log_trace("  proj_coord = %f %f %f, dist2 = %e\n", pc[0], pc[1], pc[2], distance[ipt]);
      }
    }


  } // End of loop on points



  if (_tri_state != tri_state && _tri_state != NULL) {
    PDM_triangulate_state_destroy(_tri_state);
  }

  if (_tri_vtx != tri_vtx) {
    free(_tri_vtx);
  }


  return 0;
}





/*============================================================================
 * Public function definitions
 *============================================================================*/

/**
 * \brief Locate a set points inside a set of elements
 *
 * Elements are ordered by type (points, lines, triangles, quadrangles,
 * polygons, tetrahedra, pyramids, prisms, hexahedra, polyhedra).
 *
 * \param [in]   type_idx           Index for the element types (size = 11)
 * \param [in]   elt_vtx_idx        Index of the element-vertex connectivity
 * \param [in]   elt_vtx_coord      Coordinates of the elements' vertices
 * \param [in]   poly3d_face_idx    Index of the element-face connectivity (only for polyhedra)
 * \param [in]   face_vtx_idx       Index for the face-vertex connectivity
 * \param [in]   face_vtx           Face-vertex connectivity
 * \param [in]   face_orientation   Orientation of the faces
 * \param [in]   pts_idx            Index of points (size = n_elt + 1)
 * \param [in]   pts_coord          Coordinates of the points to locate
 * \param [in]   tolerance          Geometric tolerance
 * \param [out]  distance           Distance from points to elements (< 0 if inside, > 0 if outside)
 * \param [out]  projected_coord    Coordinates of the projection of the points on the elements
 * \param [out]  bar_coord_idx      Index for the mean-value coordinates of the projections
 * \param [out]  bar_coord          Mean-value coordinates of the projections
 *
 */

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
 double            **distance,
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
  *distance = malloc (sizeof(double) * n_pts);
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
      (*distance)[ipt] = dist2;
      (*bar_coord)[ipt] = 1.;
    }
  }


  /*
   * Edges
   */
  for (ielt = type_idx[PDM_MESH_NODAL_BAR2]; ielt < type_idx[PDM_MESH_NODAL_TRIA3]; ielt++) {

    if (pts_idx[ielt+1] == pts_idx[ielt]) continue;

    _locate_on_edge (elt_vtx_coord + elt_vtx_idx[ielt] * 3,
                     pts_idx[ielt+1] - pts_idx[ielt],
                     pts_coord + pts_idx[ielt] * 3,
                     *distance + pts_idx[ielt],
                     *projected_coord + pts_idx[ielt] * 3,
                     *bar_coord + (*bar_coord_idx)[pts_idx[ielt]]);
  }


  /*
   * Triangles
   */
  const int _connec[3] = {1, 2, 3};
  for (ielt = type_idx[PDM_MESH_NODAL_TRIA3]; ielt < type_idx[PDM_MESH_NODAL_QUAD4]; ielt++) {

    if (pts_idx[ielt+1] == pts_idx[ielt]) continue;

    _locate_on_triangles (1,
                          _connec,
                          elt_vtx_coord + elt_vtx_idx[ielt] * 3,
                          pts_idx[ielt+1] - pts_idx[ielt],
                          pts_coord + pts_idx[ielt] * 3,
                          NULL,
                          *distance + pts_idx[ielt],
                          *projected_coord + pts_idx[ielt] * 3,
                          *bar_coord + (*bar_coord_idx)[pts_idx[ielt]]);
  }


  /*
   * Quadrangles
   */
  for (ielt = type_idx[PDM_MESH_NODAL_QUAD4]; ielt < type_idx[PDM_MESH_NODAL_POLY_2D]; ielt++) {

    if (pts_idx[ielt+1] == pts_idx[ielt]) continue;

    _locate_on_quadrangle (elt_vtx_coord + elt_vtx_idx[ielt] * 3,
                           pts_idx[ielt+1] - pts_idx[ielt],
                           pts_coord + pts_idx[ielt] * 3,
                           *distance + pts_idx[ielt],
                           *projected_coord + pts_idx[ielt] * 3,
                           *bar_coord + (*bar_coord_idx)[pts_idx[ielt]]);
  }


  /*
   * Polygons
   */
  for (ielt = type_idx[PDM_MESH_NODAL_POLY_2D]; ielt < type_idx[PDM_MESH_NODAL_TETRA4]; ielt++) {

    if (pts_idx[ielt+1] == pts_idx[ielt]) continue;

    _locate_in_polygon (elt_vtx_idx[ielt+1] - elt_vtx_idx[ielt],
                        elt_vtx_coord + elt_vtx_idx[ielt] * 3,
                        pts_idx[ielt+1] - pts_idx[ielt],
                        pts_coord + pts_idx[ielt] * 3,
                        *distance + pts_idx[ielt],
                        *projected_coord + pts_idx[ielt] * 3,
                        *bar_coord + (*bar_coord_idx)[pts_idx[ielt]]);
  }



  /*
   * Tetrahedra
   */
  for (ielt = type_idx[PDM_MESH_NODAL_TETRA4]; ielt < type_idx[PDM_MESH_NODAL_PYRAMID5]; ielt++) {

if (pts_idx[ielt+1] == pts_idx[ielt]) continue;

    _locate_in_tetrahedron (elt_vtx_coord + elt_vtx_idx[ielt] * 3,
                            pts_idx[ielt+1] - pts_idx[ielt],
                            pts_coord + pts_idx[ielt] * 3,
                            *distance + pts_idx[ielt],
                            *projected_coord + pts_idx[ielt] * 3,
                            *bar_coord + (*bar_coord_idx)[pts_idx[ielt]]);
  }


  /*
   * Other standard cells
   */
  for (PDM_Mesh_nodal_elt_t type = PDM_MESH_NODAL_PYRAMID5;
       type < PDM_MESH_NODAL_POLY_3D;
       type++) {

    for (ielt = type_idx[type]; ielt < type_idx[type+1]; ielt++) {

      if (pts_idx[ielt+1] == pts_idx[ielt]) continue;

      _locate_in_cell_3d (type,
                          elt_vtx_coord + elt_vtx_idx[ielt] * 3,
                          pts_idx[ielt+1] - pts_idx[ielt],
                          pts_coord + pts_idx[ielt] * 3,
                          tolerance,
                          *distance + pts_idx[ielt],
                          *projected_coord + pts_idx[ielt] * 3,
                          *bar_coord + (*bar_coord_idx)[pts_idx[ielt]],
                          NULL);
    }
  }


  /*
   * Polyhedra
   */
  int ipoly = 0;
  for (ielt = type_idx[PDM_MESH_NODAL_POLY_3D]; ielt < n_elt; ielt++) {

    if (pts_idx[ielt+1] == pts_idx[ielt]) continue;

    int n_vtx = elt_vtx_idx[ielt+1] - elt_vtx_idx[ielt];
    _locate_in_polyhedron (n_vtx,
                           elt_vtx_coord + elt_vtx_idx[ielt] * 3,
                           poly3d_face_idx[ipoly+1] - poly3d_face_idx[ipoly],
                           face_vtx_idx + poly3d_face_idx[ipoly],
                           face_vtx,
                           face_orientation + poly3d_face_idx[ipoly],
                           pts_idx[ielt+1] - pts_idx[ielt],
                           pts_coord + pts_idx[ielt] * 3,
                           *distance + pts_idx[ielt],
                           *projected_coord + pts_idx[ielt] * 3,
                           *bar_coord + (*bar_coord_idx)[pts_idx[ielt]]);

    ipoly++;
  }

}


/* Awkward to pass n_part and pvtx_coord aside from pmne... */
void
PDM_point_location_nodal2
(
 PDM_part_mesh_nodal_elmts_t   *pmne,
 const int                      n_part,
 const double                 **pvtx_coord,
 const int                    **pts_idx,
 const double                 **pts_coord,
 const double                   tolerance,
 double                      ***distance,
 double                      ***projected_coord,
 int                         ***bar_coord_idx,
 double                      ***bar_coord,
 double                      ***uvw
 )
{
  int n_section = PDM_part_mesh_nodal_elmts_n_section_get(pmne);

  int *sections_id = PDM_part_mesh_nodal_elmts_sections_id_get(pmne);

  *distance        = malloc(sizeof(double *) * n_part);
  *projected_coord = malloc(sizeof(double *) * n_part);
  *bar_coord_idx   = malloc(sizeof(int    *) * n_part);
  *bar_coord       = malloc(sizeof(double *) * n_part);
  *uvw             = malloc(sizeof(double *) * n_part);

  /* First loop to allocate */
  for (int ipart = 0; ipart < n_part; ipart++) {

    int n_cell = 0;
    for (int isection = 0; isection < n_section; isection++) {
      int n_elt = PDM_part_mesh_nodal_elmts_block_n_elt_get(pmne,
                                                            sections_id[isection],
                                                            ipart);

      n_cell += n_elt;
    }

    int n_pts = pts_idx[ipart][n_cell];

    (*distance       )[ipart] = malloc (sizeof(double) * n_pts);
    (*projected_coord)[ipart] = malloc (sizeof(double) * n_pts * 3);
    (*uvw            )[ipart] = malloc (sizeof(double) * n_pts * 3);
    (*bar_coord_idx  )[ipart] = PDM_array_zeros_int(n_pts+1);

    int *_bar_coord_idx = (*bar_coord_idx)[ipart];

    for (int isection = 0; isection < n_section; isection++) {

      int id_section = sections_id[isection];

      int n_elt = PDM_part_mesh_nodal_elmts_block_n_elt_get(pmne,
                                                            id_section,
                                                            ipart);

      PDM_Mesh_nodal_elt_t t_elt = PDM_part_mesh_nodal_elmts_block_type_get(pmne,
                                                                            id_section);

      int *parent_num = PDM_part_mesh_nodal_elmts_parent_num_get(pmne,
                                                                 id_section,
                                                                 ipart);

      if (t_elt == PDM_MESH_NODAL_POLY_2D) {
        /* Polygonal section */
        int *connec_idx;
        int *connec;
        PDM_part_mesh_nodal_elmts_block_poly2d_get(pmne,
                                                   id_section,
                                                   ipart,
                                                   &connec_idx,
                                                   &connec);

        for (int ielt = 0; ielt < n_elt; ielt++) {
          int icell = ielt;
          if (parent_num != NULL) {
            icell = parent_num[ielt];
          }

          int n_vtx = connec_idx[ielt+1] - connec_idx[ielt];

          for (int idx_pt = pts_idx[ipart][icell]; idx_pt < pts_idx[ipart][icell+1]; idx_pt++) {
            _bar_coord_idx[idx_pt+1] = n_vtx;
          }
        }
      }

      else if (t_elt == PDM_MESH_NODAL_POLY_3D) {
        /* Polyhedral section */
        int *connec_idx;
        int *connec;
        PDM_part_mesh_nodal_elmts_block_poly3d_cell_vtx_connect_get(pmne,
                                                                    id_section,
                                                                    ipart,
                                                                    &connec_idx,
                                                                    &connec);

        for (int ielt = 0; ielt < n_elt; ielt++) {
          int icell = ielt;
          if (parent_num != NULL) {
            icell = parent_num[ielt];
          }

          int n_vtx = connec_idx[ielt+1] - connec_idx[ielt];
          for (int idx_pt = pts_idx[ipart][icell]; idx_pt < pts_idx[ipart][icell+1]; idx_pt++) {
            _bar_coord_idx[idx_pt+1] = n_vtx;
          }
        }
      }

      else {
        /* Standard section */
              int         *connec              = NULL;
              PDM_g_num_t *numabs              = NULL;
              int         *_parent_num         = NULL;
              PDM_g_num_t *parent_entity_g_num = NULL;
              int          order               = 0;
        const char        *ho_ordering         = NULL;
        PDM_part_mesh_nodal_elmts_block_std_ho_get(pmne,
                                                   id_section,
                                                   ipart,
                                                   &connec,
                                                   &numabs,
                                                   &_parent_num,
                                                   &parent_entity_g_num,
                                                   &order,
                                                   &ho_ordering);

        int n_vtx = PDM_Mesh_nodal_n_vtx_elt_get(t_elt,
                                                 order);

        for (int ielt = 0; ielt < n_elt; ielt++) {
          int icell = ielt;
          if (parent_num != NULL) {
            icell = parent_num[ielt];
          }

          for (int idx_pt = pts_idx[ipart][icell]; idx_pt < pts_idx[ipart][icell+1]; idx_pt++) {
            _bar_coord_idx[idx_pt+1] = n_vtx;
          }
        }
      }

    } // End of loop on current part's sections

    for (int i = 0; i < n_pts; i++) {
      _bar_coord_idx[i+1] += _bar_coord_idx[i];
    }
    (*bar_coord)[ipart] = malloc(sizeof(double) * _bar_coord_idx[n_pts]);

  } // End of loop on parts



  /* Second loop to perform location */
  for (int ipart = 0; ipart < n_part; ipart++) {

    double *_distance        = (*distance)       [ipart];
    double *_projected_coord = (*projected_coord)[ipart];
    int    *_bar_coord_idx   = (*bar_coord_idx)  [ipart];
    double *_bar_coord       = (*bar_coord)      [ipart];
    double *_uvw             = (*uvw)            [ipart];

    const double *vtx_coord = pvtx_coord[ipart];

    for (int isection = 0; isection < n_section; isection++) {

      int id_section = sections_id[isection];

      int n_elt = PDM_part_mesh_nodal_elmts_block_n_elt_get(pmne,
                                                            id_section,
                                                            ipart);

      PDM_Mesh_nodal_elt_t t_elt = PDM_part_mesh_nodal_elmts_block_type_get(pmne,
                                                                            id_section);

      int *parent_num = PDM_part_mesh_nodal_elmts_parent_num_get(pmne,
                                                                 id_section,
                                                                 ipart);

      if (t_elt == PDM_MESH_NODAL_POLY_2D) {
        /* Polygonal section */
        int *connec_idx;
        int *connec;
        PDM_part_mesh_nodal_elmts_block_poly2d_get(pmne,
                                                   id_section,
                                                   ipart,
                                                   &connec_idx,
                                                   &connec);

        int n_vtx_max = 0;
        for (int ielt = 0; ielt < n_elt; ielt++) {
          n_vtx_max = PDM_MAX(n_vtx_max, connec_idx[ielt+1] - connec_idx[ielt]);
        }

        double *poly_coord = malloc(sizeof(double) * n_vtx_max * 3);

        for (int ielt = 0; ielt < n_elt; ielt++) {
          int icell = ielt;
          if (parent_num != NULL) {
            icell = parent_num[ielt];
          }

          int n_vtx = connec_idx[ielt+1] - connec_idx[ielt];
          int *_connec = connec + connec_idx[ielt];
          for (int i = 0; i < n_vtx; i++) {
            int vtx_id = _connec[i] - 1;
            memcpy(poly_coord + 3*i, vtx_coord + 3*vtx_id, sizeof(double) * 3);
          }

          int idx_pt = pts_idx[ipart][icell];
          int _n_pts = pts_idx[ipart][icell+1] - idx_pt;

          _locate_in_polygon(n_vtx,
                             poly_coord,
                             _n_pts,
            (const double *) pts_coord[ipart] + idx_pt * 3,
                             _distance        + idx_pt,
                             _projected_coord + idx_pt * 3,
                             _bar_coord       + _bar_coord_idx[idx_pt]);

          for (int ipt = idx_pt; ipt < pts_idx[ipart][icell+1]; ipt++) {
            _uvw[3*ipt] = _uvw[3*ipt+1] = _uvw[3*ipt+2] = -1.;
          }
        } // End of loop on polygons
        free(poly_coord);

      } // end if PDM_MESH_NODAL_POLY_2D


      else if (t_elt == PDM_MESH_NODAL_POLY_3D) {
        /* Polyhedral section */
        int *cell_vtx_idx;
        int *cell_vtx;
        PDM_part_mesh_nodal_elmts_block_poly3d_cell_vtx_connect_get(pmne,
                                                                    id_section,
                                                                    ipart,
                                                                    &cell_vtx_idx,
                                                                    &cell_vtx);

        int *cell_face_idx;
        int *cell_face;
        int  n_face;
        int *face_vtx_idx;
        int *face_vtx;
        PDM_g_num_t *face_ln_to_gn            = NULL;
        int         *_parent_num              = NULL;
        PDM_g_num_t *numabs                   = NULL;
        PDM_g_num_t *parent_entitity_ln_to_gn = NULL;
        PDM_part_mesh_nodal_elmts_block_poly3d_get(pmne,
                                                   id_section,
                                                   ipart,
                                                   &n_face,
                                                   &face_ln_to_gn,
                                                   &face_vtx_idx,
                                                   &face_vtx,
                                                   &numabs,
                                                   &cell_face_idx,
                                                   &cell_face,
                                                   &_parent_num,
                                                   &parent_entitity_ln_to_gn);

        int pn_vtx = 0;
        for (int i = 0; i < cell_vtx_idx[n_elt]; i++) {
          pn_vtx = PDM_MAX(pn_vtx, cell_vtx[i]);
        }
        int *_vtx_id = PDM_array_zeros_int(pn_vtx);

        int n_vtx_max      = 0;
        int n_face_max     = 0;
        int s_face_vtx_max = 0;
        for (int ielt = 0; ielt < n_elt; ielt++) {
          n_vtx_max  = PDM_MAX(n_vtx_max,  cell_vtx_idx [ielt+1] - cell_vtx_idx [ielt]);
          n_face_max = PDM_MAX(n_face_max, cell_face_idx[ielt+1] - cell_face_idx[ielt]);

          int s_face_vtx = 0;
          for (int idx = cell_face_idx[ielt]; idx < cell_face_idx[ielt+1]; idx++) {
            int face_id = PDM_ABS(cell_face[idx]) - 1;
            s_face_vtx += face_vtx_idx[face_id+1] - face_vtx_idx[face_id];
          }

          s_face_vtx_max = PDM_MAX(s_face_vtx_max, s_face_vtx);
        }

        int *_face_orientation = malloc(sizeof(int) * n_face_max);
        int *_face_vtx_idx = malloc(sizeof(int) * (n_face_max + 1));
        _face_vtx_idx[0] = 0;

        int *_face_vtx = malloc(sizeof(int) * s_face_vtx_max);
        double *_vtx_coord = malloc(sizeof(double) * n_vtx_max * 3);

        PDM_triangulate_state_t *_tri_state = NULL;
        _tri_state = PDM_triangulate_state_create(s_face_vtx_max);
        int *_tri_vtx = malloc(sizeof(int) * (s_face_vtx_max - 2) * 3);


        for (int ielt = 0; ielt < n_elt; ielt++) {
          int icell = ielt;
          if (parent_num != NULL) {
            icell = parent_num[ielt];
          }

          int _n_vtx  = cell_vtx_idx [ielt+1] - cell_vtx_idx [ielt];
          int _n_face = cell_face_idx[ielt+1] - cell_face_idx[ielt];

          for (int ivtx = 0; ivtx < _n_vtx; ivtx++) {
            int vtx_id = cell_vtx[cell_vtx_idx[ielt] + ivtx] - 1;
            _vtx_id[vtx_id] = ivtx+1;
            memcpy(_vtx_coord + 3*ivtx, pvtx_coord[ipart] + 3*vtx_id, sizeof(double) * 3);
          }

          for (int iface = 0; iface < _n_face; iface++) {
            int idx_face = cell_face_idx[ielt] + iface;
            int face_id = PDM_ABS(cell_face[idx_face]) - 1;
            _face_orientation[iface] = PDM_SIGN(cell_face[idx_face]);
            _face_vtx_idx[iface+1] = _face_vtx_idx[iface];
            for (int idx_vtx = face_vtx_idx[face_id]; idx_vtx < face_vtx_idx[face_id+1]; idx_vtx++) {
              int vtx_id = face_vtx[idx_vtx] - 1;
              assert(_vtx_id[vtx_id] > 0);
              _face_vtx[_face_vtx_idx[iface+1]++] = _vtx_id[vtx_id];
            }
          }

          int idx_pt = pts_idx[ipart][icell];
          int _n_pts = pts_idx[ipart][icell+1] - idx_pt;

          // _locate_in_polyhedron(_n_vtx,
          //      (const double *) _vtx_coord,
          //                       _n_face,
          //                       _face_vtx_idx,
          //      (const int    *) _face_vtx,
          //      (const int    *) _face_orientation,
          //                       _n_pts,
          //      (const double *) pts_coord[ipart]  + idx_pt * 3,
          //                       _distance         + idx_pt,
          //                       _projected_coord  + idx_pt * 3,
          //                       _bar_coord        + _bar_coord_idx[idx_pt]);
          _locate_in_polyhedron2(_n_vtx,
                (const double *) _vtx_coord,
                                 _n_face,
                                 _face_vtx_idx,
                (const int    *) _face_vtx,
                (const int    *) _face_orientation,
                                 _n_pts,
                (const double *) pts_coord[ipart]  + idx_pt * 3,
                                 _distance         + idx_pt,
                                 _projected_coord  + idx_pt * 3,
                                 _bar_coord        + _bar_coord_idx[idx_pt],
                                 _tri_state,
                                 _tri_vtx);

          for (int ipt = idx_pt; ipt < pts_idx[ipart][icell+1]; ipt++) {
            _uvw[3*ipt] = _uvw[3*ipt+1] = _uvw[3*ipt+2] = -1.;
          }

          // Reset local vtx_id
          for (int idx = cell_vtx_idx[ielt]; idx < cell_vtx_idx[ielt+1]; idx++) {
            int vtx_id = cell_vtx[idx] - 1;
            _vtx_id[vtx_id] = 0;
          }
        } // End of loop on polyhedra

        free(_vtx_id);
        free(_face_orientation);
        free(_face_vtx_idx);
        free(_face_vtx);
        free(_vtx_coord);
        free(_tri_vtx);
        PDM_triangulate_state_destroy(_tri_state);

      } // end if PDM_MESH_NODAL_POLY_3D


      else {
        /* Standard section */
              int         *connec              = NULL;
              PDM_g_num_t *numabs              = NULL;
              int         *_parent_num         = NULL;
              PDM_g_num_t *parent_entity_g_num = NULL;
              int          order               = 0;
        const char        *ho_ordering         = NULL;
        PDM_part_mesh_nodal_elmts_block_std_ho_get(pmne,
                                                   id_section,
                                                   ipart,
                                                   &connec,
                                                   &numabs,
                                                   &_parent_num,
                                                   &parent_entity_g_num,
                                                   &order,
                                                   &ho_ordering);


        switch (t_elt) {

          case PDM_MESH_NODAL_POINT: {
            for (int ielt = 0; ielt < n_elt; ielt++) {
              int icell = ielt;
              if (parent_num != NULL) {
                icell = parent_num[ielt];
              }

              const double *node_coord = pvtx_coord[ipart] + 3*(connec[ielt] - 1);
              for (int idx_pt = pts_idx[ipart][icell]; idx_pt < pts_idx[ipart][icell+1]; idx_pt++) {
                double dist2 = 0.;
                for (int j = 0; j < 3; j++) {
                  double dj = pts_coord[ipart][3*idx_pt+j] - node_coord[j];
                  dist2 += dj * dj;
                }
                _distance [idx_pt] = dist2;
                _bar_coord[_bar_coord_idx[idx_pt]] = 1.;
                _uvw      [3*idx_pt  ] = 0.;
                _uvw      [3*idx_pt+1] = -1.;
                _uvw      [3*idx_pt+2] = -1.;
              }
            }
            break;
          } // end case PDM_MESH_NODAL_POINT


          case PDM_MESH_NODAL_BAR2: {
            double edge_coord[6];
            for (int ielt = 0; ielt < n_elt; ielt++) {
              int icell = ielt;
              if (parent_num != NULL) {
                icell = parent_num[ielt];
              }

              int *_connec = connec + 2*ielt;
              for (int i = 0; i < 2; i++) {
                int vtx_id = _connec[i] - 1;
                memcpy(edge_coord + 3*i, vtx_coord + 3*vtx_id, sizeof(double) * 3);
              }

              int idx_pt = pts_idx[ipart][icell];

              _locate_on_edge(edge_coord,
                              pts_idx[ipart][icell+1] - idx_pt,
             (const double *) pts_coord[ipart] + idx_pt * 3,
                              _distance        + idx_pt,
                              _projected_coord + idx_pt * 3,
                              _bar_coord       + _bar_coord_idx[idx_pt]);

              for (int ipt = idx_pt; ipt < pts_idx[ipart][icell+1]; ipt++) {
                _uvw[3*ipt  ] = _bar_coord[_bar_coord_idx[ipt]+1];
                _uvw[3*ipt+1] = -1.;
                _uvw[3*ipt+2] = -1.;
              }

            } // End of loop on edges
            break;
          } // end case PDM_MESH_NODAL_BAR2


          case PDM_MESH_NODAL_TRIA3: {
            for (int ielt = 0; ielt < n_elt; ielt++) {
              int icell = ielt;
              if (parent_num != NULL) {
                icell = parent_num[ielt];
              }

              int *_connec = connec + 3*ielt;
              int idx_pt = pts_idx[ipart][icell];

              _locate_on_triangles(1,
                                   _connec,
                                   pvtx_coord[ipart],
                                   pts_idx[ipart][icell+1] - idx_pt,
                  (const double *) pts_coord[ipart] + idx_pt * 3,
                                   NULL,
                                   _distance        + idx_pt,
                                   _projected_coord + idx_pt * 3,
                                   _bar_coord       + _bar_coord_idx[idx_pt]);

              for (int ipt = idx_pt; ipt < pts_idx[ipart][icell+1]; ipt++) {
                _uvw[3*ipt  ] = _bar_coord[_bar_coord_idx[ipt]+1];
                _uvw[3*ipt+1] = _bar_coord[_bar_coord_idx[ipt]+2];
                _uvw[3*ipt+2] = -1.;
              }

            } // End of loop on triangles

            break;
          } // end case PDM_MESH_NODAL_TRIA3


          case PDM_MESH_NODAL_QUAD4: {
            double quad_coord[12];
            for (int ielt = 0; ielt < n_elt; ielt++) {
              int icell = ielt;
              if (parent_num != NULL) {
                icell = parent_num[ielt];
              }

              int *_connec = connec + 4*ielt;
              for (int i = 0; i < 4; i++) {
                int vtx_id = _connec[i] - 1;
                memcpy(quad_coord + 3*i, vtx_coord + 3*vtx_id, sizeof(double) * 3);
              }

              int idx_pt = pts_idx[ipart][icell];

              _locate_on_quadrangle(quad_coord,
                                    pts_idx[ipart][icell+1] - idx_pt,
                   (const double *) pts_coord[ipart] + idx_pt * 3,
                                    _distance        + idx_pt,
                                    _projected_coord + idx_pt * 3,
                                    _bar_coord       + _bar_coord_idx[idx_pt]);

              // TODO: compute uvw
              for (int ipt = idx_pt; ipt < pts_idx[ipart][icell+1]; ipt++) {
                _uvw[3*ipt  ] = -1;
                _uvw[3*ipt+1] = -1;
                _uvw[3*ipt+2] = -1.;
              }

            } // End of loop on quadrangles
            break;
          } // end case PDM_MESH_NODAL_QUAD4

          case PDM_MESH_NODAL_TETRA4: {
            double tetra_coord[12];
            for (int ielt = 0; ielt < n_elt; ielt++) {
              int icell = ielt;
              if (parent_num != NULL) {
                icell = parent_num[ielt];
              }

              int *_connec = connec + 4*ielt;
              for (int i = 0; i < 4; i++) {
                int vtx_id = _connec[i] - 1;
                memcpy(tetra_coord + 3*i, vtx_coord + 3*vtx_id, sizeof(double) * 3);
              }

              int idx_pt = pts_idx[ipart][icell];

              _locate_in_tetrahedron(tetra_coord,
                                     pts_idx[ipart][icell+1] - idx_pt,
                    (const double *) pts_coord[ipart] + idx_pt * 3,
                                     _distance        + idx_pt,
                                     _projected_coord + idx_pt * 3,
                                     _bar_coord       + _bar_coord_idx[idx_pt]);

              for (int ipt = idx_pt; ipt < pts_idx[ipart][icell+1]; ipt++) {
                _uvw[3*ipt  ] = _bar_coord[_bar_coord_idx[ipt]+1];
                _uvw[3*ipt+1] = _bar_coord[_bar_coord_idx[ipt]+2];
                _uvw[3*ipt+2] = _bar_coord[_bar_coord_idx[ipt]+3];
              }

            } // End of loop on tetrahedra
            break;
          } // end case PDM_MESH_NODAL_TETRA4


          case PDM_MESH_NODAL_PYRAMID5:
          case PDM_MESH_NODAL_PRISM6  :
          case PDM_MESH_NODAL_HEXA8   : {
            double elt_coord[24];
            int n_vtx = PDM_Mesh_nodal_n_vtx_elt_get(t_elt,
                                                     order);
            for (int ielt = 0; ielt < n_elt; ielt++) {
              int icell = ielt;
              if (parent_num != NULL) {
                icell = parent_num[ielt];
              }

              int *_connec = connec + n_vtx*ielt;
              for (int i = 0; i < n_vtx; i++) {
                int vtx_id = _connec[i] - 1;
                memcpy(elt_coord + 3*i, vtx_coord + 3*vtx_id, sizeof(double) * 3);
              }

              int idx_pt = pts_idx[ipart][icell];
              // for (int ipt = idx_pt; ipt < pts_idx[ipart][icell+1]; ipt++) {
              //   _uvw[3*ipt] = _uvw[3*ipt+1] = _uvw[3*ipt+2] = -1.;
              // }

              _locate_in_cell_3d(t_elt,
                                 elt_coord,
                                 pts_idx[ipart][icell+1] - idx_pt,
                (const double *) pts_coord[ipart] + idx_pt * 3,
                                 tolerance,
                                 _distance        + idx_pt,
                                 _projected_coord + idx_pt * 3,
                                 _bar_coord       + _bar_coord_idx[idx_pt],
                                 _uvw             + idx_pt * 3);

            } // End of loop on pyramids/prisms/hexahedra
            break;
          } // end case PDM_MESH_NODAL_PYRAMID5/PRISM6/HEXA8

          case PDM_MESH_NODAL_BARHO    :
          case PDM_MESH_NODAL_TRIAHO   :
          case PDM_MESH_NODAL_QUADHO   :
          case PDM_MESH_NODAL_TETRAHO  :
          case PDM_MESH_NODAL_PYRAMIDHO:
          case PDM_MESH_NODAL_PRISMHO  :
          case PDM_MESH_NODAL_HEXAHO   : {

            int n_node = PDM_Mesh_nodal_n_vtx_elt_get(t_elt,
                                                      order);
            double *elt_coord = malloc(sizeof(double) * n_node * 3);

            int *ijk_to_user = NULL;
            if (ho_ordering != NULL) {
              ijk_to_user = PDM_ho_ordering_ijk_to_user_get(ho_ordering,
                                                            t_elt,
                                                            order);
              assert(ijk_to_user != NULL);
            }


            for (int ielt = 0; ielt < n_elt; ielt++) {
              int icell = ielt;
              if (parent_num != NULL) {
                icell = parent_num[ielt];
              }

              int *_connec = connec + n_node*ielt;
              for (int idx_ijk = 0; idx_ijk < n_node; idx_ijk++) {
                int idx_user = idx_ijk;
                if (ijk_to_user != NULL) {
                  idx_user = ijk_to_user[idx_ijk];
                }
                int node_id = _connec[idx_user] - 1;
                memcpy(elt_coord + 3*idx_ijk, vtx_coord + 3*node_id, sizeof(double) * 3);
              }

              for (int idx_pt = pts_idx[ipart][icell]; idx_pt < pts_idx[ipart][icell+1]; idx_pt++) {
                _distance[idx_pt] = PDM_ho_location(t_elt,
                                                    order,
                                                    n_node,
                                                    elt_coord,
                                   (const double *) pts_coord[ipart] + idx_pt * 3,
                                                    _projected_coord + idx_pt * 3,
                                                    _uvw             + idx_pt * 3);

              } // End of loop on current elt's points

              int idx_pt0 = pts_idx[ipart][icell];
              PDM_ho_basis(t_elt,
                           order,
                           n_node,
                           pts_idx[ipart][icell+1] - idx_pt0,
                           (const double *) _uvw + idx_pt0 * 3,
                           _bar_coord + _bar_coord_idx[idx_pt0]);

            } // End of loop on elt
            free(elt_coord);
            break;
          } // end case HO LAGRANGE


          // HO_BEZIER...

          default: {
            PDM_error(__FILE__, __LINE__, 0,
                      "PDM_point_location_nodal : element type %d not supported yet\n",
                      (int) t_elt);
          }

        } // End switch on std elt type


      } // End if standard section

    } // End of loop on current part's sections

  } // End of loop on parts

}


/**
 * \brief Compute hexahedron, pyramid, or prism parametric coordinates for a
 * given point.
 *
 * This function is adapted from the CGNS interpolation tool.
 *
 * \param [in]   elt_type       Type of element
 * \param [in]   point_coords   Point coordinates
 * \param [in]   vertex_coords  Pointer to element vertex coordinates
 * \param [in]   tolerance      Location tolerance factor
 * \param [out]  uvw            Parametric coordinates of point in element
 *
 */

PDM_bool_t
PDM_point_location_compute_uvw
(
 const PDM_Mesh_nodal_elt_t elt_type,
 const double               point_coords[3],
 const double               vertex_coords[],
 const double               tolerance,
 double                     uvw[3]
 )
{
  return _compute_uvw (elt_type,
                       point_coords,
                       vertex_coords,
                       tolerance,
                       uvw);
}

#ifdef __cplusplus
}
#endif /* __cplusplus */
