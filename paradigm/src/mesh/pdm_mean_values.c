/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_printf.h"
#include "pdm_plane.h"
#include "pdm_line.h"
#include "pdm_polygon.h"
#include "pdm_geom_elem.h"
#include "pdm_triangulate.h"

#include "pdm_mean_values.h"

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
 * Type
 *============================================================================*/

/*=============================================================================
 * Private function definition
 *============================================================================*/

/*=============================================================================
 * Public function definition
 *============================================================================*/

/**
 * \brief Compute mean values of a list of points in a polygon
 *
 * \param [in]    n_pts            Number of points to locate
 * \param [in]    pts              xyz-Coordinates of points locate
 * \param [in]    n_vtx            Number of polygon vertices
 * \param [in]    poly_vtx         Polygon connectivity
 * \param [in]    mesh_vtx_coords  Coordinates of mesh vertices
 * \param [inout] mesh_vtx_coords  Mean value coordinates of points to locate
 *
 *
 */

void
PDM_mean_value_coordinates_polygon_2d
(
 const int    n_vtx,
 const double vtx_coord[],
 const int    n_pts,
 const double pts_coord[],
 double       mean_value_coord[]
)
{
  /*
    See "Mean value coordinates for arbitrary planar polygons",
    Kai Hormann, and Michael S. Floater. (2006)
  */
  const double eps_base = 1e-12;

  int ipt, ivtx, jvtx, idim;

  /* Compute polygon bounds */
  double bounds[4] = {DBL_MAX, -DBL_MAX,
                      DBL_MAX, -DBL_MAX};
  double char_length = eps_base;
  for (ivtx = 0; ivtx < n_vtx; ivtx++) {
    for (idim = 0; idim < 2; idim++) {
      bounds[2*idim]   = PDM_MIN (bounds[2*idim],   vtx_coord[2*ivtx+idim]);
      bounds[2*idim+1] = PDM_MAX (bounds[2*idim+1], vtx_coord[2*ivtx+idim]);

      char_length = PDM_MAX (char_length, bounds[2*idim+1] - bounds[2*idim]);
    }
  }

  const double eps = eps_base * char_length;
  const double eps2 = eps * eps;

  double *s = malloc (sizeof(double) * n_vtx * 2);
  double *r = malloc (sizeof(double) * n_vtx);
  double *A = malloc (sizeof(double) * n_vtx);
  double *D = malloc (sizeof(double) * n_vtx);

  /* Loop on points */
  for (ipt = 0; ipt < n_pts; ipt++) {
    const double *_pt = pts_coord + 2 * ipt;
    double       *_bc = mean_value_coord + n_vtx * ipt;

    /* If current point is outside the polygon,
       consider its projection on the polygon boundary */
        if (PDM_polygon_point_in_2d (_pt,
                                     n_vtx,
                                     vtx_coord,
                                     //                             char_length,
                                     bounds) != PDM_POLYGON_INSIDE) {
      double dist2_min = DBL_MAX;
      int i_min;
      double t_min;

      for (ivtx = 0; ivtx < n_vtx; ivtx++) {
        jvtx = (ivtx + 1) % n_vtx;
        double t;
        double closest[2];
        double dist2 = PDM_line_distance_2d (_pt,
                                             vtx_coord + 2*ivtx,
                                             vtx_coord + 2*jvtx,
                                             &t,
                                             closest);

        if (dist2 < dist2_min) {
          dist2_min = dist2;
          i_min = ivtx;

          if (t < 0.) {
            t_min = 0.;
          } else if (t > 1.) {
            t_min = 1.;
          } else  {
            t_min = t;
          }
        }
      }

      for (int i = 0; i < n_vtx; i++) {
        _bc[i] = 0.;
      }
      _bc[i_min]           = 1.0 - t_min;
      _bc[(i_min+1)%n_vtx] = t_min;
      continue;
    }

    PDM_bool_t special_case = PDM_FALSE;
    for (ivtx = 0; ivtx < n_vtx; ivtx++) {
      double *vec = s + 2*ivtx;
      for (idim = 0; idim < 2; idim++) {
        vec[idim] = vtx_coord[2*ivtx + idim] - _pt[idim];
      }
      r[ivtx] = PDM_DOT_PRODUCT_2D (vec, vec);

      if (r[ivtx] < eps2) {
        /* Point coincident with vertex */
        for (int i = 0; i < n_vtx; i++) {
          _bc[i] = 0.;
        }
        _bc[ivtx] = 1.;

        special_case = PDM_TRUE;
        break;
      }

      else {
        r[ivtx] = sqrt (r[ivtx]);
      }
    }


    if (special_case == PDM_TRUE) {
      continue;
    }

    for (ivtx = 0; ivtx < n_vtx; ivtx++) {
      jvtx = (ivtx + 1) % n_vtx;
      double *veci = s + 2*ivtx;
      double *vecj = s + 2*jvtx;

      A[ivtx] = veci[0] * vecj[1] - veci[1] * vecj[0];
      D[ivtx] = PDM_DOT_PRODUCT_2D (veci, vecj);

      /* Point on edge */
      if (fabs(A[ivtx]) < eps) {
        for (int i = 0; i < n_vtx; i++) {
          _bc[i] = 0.;
        }

        _bc[ivtx] = r[jvtx] / (r[ivtx] + r[jvtx]);
        _bc[jvtx] = r[ivtx] / (r[ivtx] + r[jvtx]);

        special_case = PDM_TRUE;
        break;
      }
    }

    if (special_case == PDM_TRUE) {
      continue;
    }

    /* General case (point strictly inside polygon) */
    double sum_w = 0.;
    for (int i = 0; i < n_vtx; i++) {
      int ip = (i + 1) % n_vtx;
      int im = (i - 1 + n_vtx) % n_vtx;

      _bc[i] = (r[ip] - D[i]/r[i]) / A[i] + (r[im] - D[im]/r[i]) / A[im];

      sum_w += _bc[i];
    }

    if (fabs(sum_w) > eps_base) {
      sum_w = 1. / sum_w;
      for (int i = 0; i < n_vtx; i++) {
        _bc[i] *= sum_w;
      }
    }

    else {
      printf("!!! sum_w = %g\n", sum_w);
    }
  }

  free (s);
  free (r);
  free (A);
  free (D);
}




void
PDM_mean_value_coordinates_polygon_3d
(
 const int    n_vtx,
 const double vtx_coord[],
 const int    n_pts,
 const double pts_coord[],
 double       mean_value_coord[]
)
{

  double *vtx_uv = malloc (sizeof(double) * n_vtx * 2);
  double *pts_uv = malloc (sizeof(double) * n_pts * 2);

  PDM_polygon_3d_to_2d (n_vtx,
                        vtx_coord,
                        vtx_uv,
                        n_pts,
                        pts_coord,
                        pts_uv,
                        NULL);

  PDM_mean_value_coordinates_polygon_2d (n_vtx,
                                         vtx_uv,
                                         n_pts,
                                         pts_uv,
                                         mean_value_coord);

  free (vtx_uv);
  free (pts_uv);
}





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


/**
 * See "Mean value coordinates for closed triangular meshes", T. Ju et al. (2005)
 *
 **/
void
PDM_mean_value_coordinates_polyhedron
(
 const int         n_vtx,
 const double      vtx_coord[],
 const PDM_l_num_t n_face,
 const PDM_l_num_t face_vtx_idx[],
 const PDM_l_num_t face_vtx[],
 const int         face_orientation[],
 const double      pt_coord[],
 double            mean_value_coord[]
 )
{
  const int DEBUG = 0;//(pt_coord[0] < 0);
  if (DEBUG) {
    printf("\n\n-- PDM_mean_value_coordinates_polyhedron3 --\n");
    printf("pt_coord = %f %f %f\n", pt_coord[0], pt_coord[1], pt_coord[2]);
  }
  const int LOCATE_ON_TRIANGLES = 1;

  const double eps = 1.e-9;
  const double eps2 = eps * eps;

  double *u = malloc (sizeof(double) * n_vtx * 3);
  double *d = malloc (sizeof(double) * n_vtx);

  for (int ivtx = 0; ivtx < n_vtx; ivtx++) {
    mean_value_coord[ivtx] = 0.;
  }

  /*
   *  Offset vertices and project on unit sphere centered at point to locate
   */
  for (int ivtx = 0; ivtx < n_vtx; ivtx++) {
    double *_u = u + 3*ivtx;
    for (int idim = 0; idim < 3; idim++) {
      _u[idim] = vtx_coord[3*ivtx + idim] - pt_coord[idim];
    }
    double uu = PDM_DOT_PRODUCT (_u, _u);

    /* Point coincident with vertex */
    if (uu < eps2) {
      mean_value_coord[ivtx] = 1.;

      free (u);
      free (d);
      return;
    }

    d[ivtx] = sqrt(uu);

    for (int idim = 0; idim < 3; idim++) {
      _u[idim] = _u[idim] / d[ivtx];
    }
  } // End of loop on vertices



  /*
   *  Prepare face triangulation
   */
  /* Count max nb of vertices per face */
  PDM_l_num_t n_vtx_face, n_vtx_face_max = 0;
  for (int iface = 0; iface < n_face; iface++) {
    n_vtx_face = face_vtx_idx[iface+1] - face_vtx_idx[iface];
    if (n_vtx_face > n_vtx_face_max) {
      n_vtx_face_max = n_vtx_face;
    }
  }

  /*
   *  Loop on faces
   */
  PDM_l_num_t n_tri;
  PDM_l_num_t _tri_vtx[3];
  PDM_l_num_t *tri_vtx = malloc (sizeof(PDM_l_num_t) * (n_vtx_face_max - 2)*3);
  PDM_triangulate_state_t *state = PDM_triangulate_state_create (n_vtx_face_max);

  for (int iface = 0; iface < n_face; iface++) {

    const PDM_l_num_t *_face_vtx = face_vtx + face_vtx_idx[iface];

    /*
     *  Triangulate current face
     */
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


    /* Loop on triangles */
    for (int itri = 0; itri < n_tri; itri++) {
      for (int ivtx = 0; ivtx < 3; ivtx++) {
        _tri_vtx[ivtx] = tri_vtx[3*itri + ivtx] - 1;
      }

      if (face_orientation[iface] < 0) {
        PDM_l_num_t tmp = _tri_vtx[0];
        _tri_vtx[0] = _tri_vtx[2];
        _tri_vtx[2] = tmp;
      }

      double l[3] = {0., 0., 0.}, theta[3], sint[3], h = 0.;
      for (int i = 0; i < 3; i++) {
        int ip = (i+1)%3;
        int im = (i+2)%3;
        for (int idim = 0; idim < 3; idim++) {
          double delta = u[3*_tri_vtx[ip] + idim] - u[3*_tri_vtx[im] + idim];
          l[i] += delta * delta;
        }
        l[i] = sqrt(l[i]);

        theta[i] = asin(0.5 * l[i]);
        h += theta[i];
        theta[i] *= 2.;
        sint[i] = sin(theta[i]);
      }

      if (PDM_PI - h < eps) {
        /*
         *  point lies on current tirangle, use 2D barycentric coordinates
         */

        /* Triangular face */
        if (n_tri == 1 || LOCATE_ON_TRIANGLES) {
          for (int ivtx = 0; ivtx < n_vtx; ivtx++) {
            mean_value_coord[ivtx] = 0.;
          }

          double sum_w = 0.;
          for (int i = 0; i < 3; i++) {
            int ivtx = _tri_vtx[i];
            int ip = _tri_vtx[(i+1)%3];
            int im = _tri_vtx[(i+2)%3];
            double w = sint[i] * d[ip] * d[im];
            sum_w += w;
            mean_value_coord[ivtx] = w;
          }

          for (int i = 0; i < 3; i++) {
            mean_value_coord[_tri_vtx[i]] /= sum_w;
          }
        }

        /* Polygonal face */
        else {
          n_vtx_face = face_vtx_idx[iface+1] - face_vtx_idx[iface];
          double *face_coord = malloc (sizeof(double) * n_vtx_face * 3);
          double *mean_value_coord_face = malloc (sizeof(double) * n_vtx_face);

          for (int j = 0; j < n_vtx_face; j++) {
            int id_vtx = face_vtx[face_vtx_idx[iface] + j] - 1;
            for (int idim = 0; idim < 3; idim++) {
              face_coord[3*j + idim] = vtx_coord[3*id_vtx + idim];
            }
          }

          PDM_mean_value_coordinates_polygon_3d (n_vtx_face,
                                                 face_coord,
                                                 1,
                                                 pt_coord,
                                                 mean_value_coord_face);

          for (int j = 0; j < n_vtx; j++) {
            mean_value_coord[j] = 0.;
          }

          for (int j = 0; j < n_vtx_face; j++) {
            int id_vtx = face_vtx[face_vtx_idx[iface] + j] - 1;
            mean_value_coord[id_vtx] = mean_value_coord_face[j];
          }

          free (face_coord);
          free (mean_value_coord_face);
        }

        free (u);
        free (d);
        free (tri_vtx);

        if (DEBUG) {
          printf("point located on face %d (PDM_PI - h = %g)\n", iface, PDM_PI - h);
        }
        return;
      }

      double c[3], s[3];
      double det = _determinant_3x3 ((u + 3*_tri_vtx[0]),
                                     (u + 3*_tri_vtx[1]),
                                     (u + 3*_tri_vtx[2]));
      PDM_bool_t ignore_triangle = PDM_FALSE;
      for (int i = 0; i < 3; i++) {
        c[i] = 2. * sin(h) * sin(h - theta[i]) / (sint[(i+1)%3] * sint[(i+2)%3]) - 1.;
        s[i] = sqrt(1. - c[i]*c[i]);

        if (s[i] < eps) {
          /* point lies outside current triangle on the same plane, ignore current triangle */
          ignore_triangle = PDM_TRUE;
          break;
        }

        if (det < 0.) {
          s[i] = -s[i];
        }
      }

      if (ignore_triangle == PDM_TRUE) {
        if (DEBUG) {
          printf("ignore triangle %d of face %d\n", itri, iface);
        }
        continue;
      }

      for (int i = 0; i < 3; i++) {
        int ip = (i+1)%3;
        int im = (i+2)%3;
        double w = (theta[i] - c[ip]*theta[im] - c[im]*theta[ip]) / (sint[ip] * s[im]);

        mean_value_coord[_tri_vtx[i]] += w;
      }

    } // End of loop on triangles

  } // End of loop on faces

  state = PDM_triangulate_state_destroy (state);

  /* Normalize */
  double sum = 0.;
  for (int ivtx = 0; ivtx < n_vtx; ivtx++) {
    mean_value_coord[ivtx] /= d[ivtx];
    sum += mean_value_coord[ivtx];
  }

  if (fabs(sum) > 1.e-15) {
    sum = 1. / sum;
    for (int ivtx = 0; ivtx < n_vtx; ivtx++) {
      mean_value_coord[ivtx] *= sum;
    }
  }

  free (u);
  free (d);
  free (tri_vtx);
}




static inline double _distance2(const double a[3], const double b[3]) {
  return (a[0] - b[0]) * (a[0] - b[0])
    +    (a[1] - b[1]) * (a[1] - b[1])
    +    (a[2] - b[2]) * (a[2] - b[2]);
}


void
PDM_mean_values_polyhedron
(
 const int         n_vtx,
 const double      vtx_coord[],
 const PDM_l_num_t n_face,
 const PDM_l_num_t face_vtx_idx[],
 const PDM_l_num_t face_vtx[],
 const int         face_orientation[],
 const double      pt_coord[],
 double            weights[]
 )
{
  const int DEBUG = 0;
  /*printf("\n\n\nPDM_mean_values_polyhedron for pt %f %f %f\n",
    pt_coord[0], pt_coord[1], pt_coord[2]);*/
  // Begin by initializing weights.
  int ivtx, jvtx, idim, i, im;
  double mag, l, l2;
  double *ui, *uip;

  const double eps = 1.e-8;
  const double eps2 = eps * eps;

  for (i = 0; i < n_vtx; i++) {
    weights[i] = 0.;
  }

  // create local array for storing point-to-vertex vectors and distances
  double *inv_dist = malloc (sizeof(double) * n_vtx);
  double *u = malloc (sizeof(double) * n_vtx * 3);

  for (i = 0; i < n_vtx; i++) {
    // point-to-vertex vector
    ui = u + 3 * i;
    for (idim = 0; idim < 3; idim++) {
      ui[idim] = vtx_coord[3*i + idim] - pt_coord[idim];
    }

    // distance
    l2 = PDM_DOT_PRODUCT (ui, ui);

    // handle special case when the point is really close to a vertex
    if (l2 < eps2) {
      weights[i] = 1.0;
      free (inv_dist);
      free (u);
      return;
    }


    // project onto unit sphere
    inv_dist[i] = 1. / sqrt(l2);
    for (idim = 0; idim < 3; idim++) {
      ui[idim] *= inv_dist[i];
    }
  }


  // Now loop over all triangle to compute weights
  double *tan_half_alpha = malloc (sizeof(double) * n_vtx);
  double *theta = malloc (sizeof(double) * n_vtx);

  for (int iface = 0; iface < n_face; iface++) {
    int n_vtx_face = face_vtx_idx[iface+1] - face_vtx_idx[iface];
    const PDM_l_num_t *_face_vtx = face_vtx + face_vtx_idx[iface];

    assert (n_vtx_face < n_vtx);

    // unit vector v.
    double v[3] = {0., 0., 0.};
    double angle;
    double temp[3];
    for (i = 0; i < n_vtx_face; i++) {
      if (face_orientation[iface] < 0) {
        ui  = u + 3*(_face_vtx[i] - 1);
        uip = u + 3*(_face_vtx[(i+n_vtx_face-1)%n_vtx_face] - 1);
      } else {
        ui  = u + 3*(_face_vtx[i] - 1);
        uip = u + 3*(_face_vtx[(i+1)%n_vtx_face] - 1);
      }

      PDM_CROSS_PRODUCT (temp, ui, uip);
      mag = PDM_MODULE (temp);
      if (mag < eps) {
        printf("!!! face %d, mag(u[%d] x u[%d]) = %f\n",
               iface, _face_vtx[i] - 1, _face_vtx[(i+1)%n_vtx_face] - 1, mag);
      }
      temp[0] /= mag;
      temp[1] /= mag;
      temp[2] /= mag;

      l = sqrt (_distance2 (ui, uip));
      angle = asin(0.5 * l);

      v[0] += angle * temp[0];
      v[1] += angle * temp[1];
      v[2] += angle * temp[2];
    }

    const double mag_v = PDM_MODULE (v);
    if (mag_v < eps) {
      printf("!!! face %d, mag(v) = %f\n", iface, mag_v);
    }
    v[0] /= mag_v;
    v[1] /= mag_v;
    v[2] /= mag_v;


    // angles between edges
    double n0[3], n1[3];
    for (i = 0; i < n_vtx_face; i++) {
      if (face_orientation[iface] < 0) {
        ui  = u + 3*(_face_vtx[i] - 1);
        uip = u + 3*(_face_vtx[(i+n_vtx_face-1)%n_vtx_face] - 1);
      } else {
        ui  = u + 3*(_face_vtx[i] - 1);
        uip = u + 3*(_face_vtx[(i+1)%n_vtx_face] - 1);
      }

      // alpha
      PDM_CROSS_PRODUCT (n0, ui, v);
      mag = PDM_MODULE (n0);
      if (mag < eps) {
        printf("!!! face %d, i = %d, mag(n0) = %f\n", iface, i, mag);
      }
      n0[0] /= mag;
      n0[1] /= mag;
      n0[2] /= mag;

      PDM_CROSS_PRODUCT (n1, uip, v);
      mag = PDM_MODULE (n1);
      if (mag < eps) {
        printf("!!! face %d, i = %d, mag(n1) = %f\n", iface, i, mag);
      }
      n1[0] /= mag;
      n1[1] /= mag;
      n1[2] /= mag;

      l2 = _distance2 (n0, n1);
      tan_half_alpha[i] = 0.5 * sqrt(l2 / (1.0 - 0.25 * l2));
      if (PDM_DOT_PRODUCT (temp, v) < 0) {
        tan_half_alpha[i] = -tan_half_alpha[i];
      }

      // theta
      l = sqrt (_distance2 (ui, v));
      theta[i] = 2.0 * asin(0.5 * l);
    }


    PDM_bool_t outlier = PDM_FALSE;
    for (i = 0; i < n_vtx_face; i++) {
      if (theta[i] < eps) {
        outlier = PDM_TRUE;
        ivtx = _face_vtx[i] - 1;
        weights[ivtx] += mag_v * inv_dist[ivtx];
        break;
      }
    }

    if (outlier == PDM_TRUE) {
      if (DEBUG) {
        printf("Outlier, iface = %d\n", iface);
      }
      continue;
    }

    double sum = 0.0;
    for (i = 0; i < n_vtx_face; i++) {
      if (face_orientation[iface] < 0) {
        im = (i + 1)%n_vtx_face;
      } else {
        im = (i + n_vtx_face - 1)%n_vtx_face;
      }
      sum += (tan_half_alpha[i] + tan_half_alpha[im]) / tan(theta[i]);
    }


    // the special case when x lies on the polygon, handle it using 2D mvc.
    // in the 2D case, alpha = theta
    if (fabs(sum) < eps) {
      for (ivtx = 0; ivtx < n_vtx; ivtx++) {
        weights[ivtx] = 0.;
      }

      if (DEBUG) {
        printf("Point lies on polygon\n");
      }

      // recompute theta, the theta computed previously are not robust
      for (i = 0; i < n_vtx_face; i++) {

        ivtx = _face_vtx[i] - 1;

        if (face_orientation[iface] < 0) {
          jvtx = _face_vtx[(i + n_vtx_face - 1)%n_vtx_face] - 1;
        } else {
          jvtx = _face_vtx[(i + 1)%n_vtx_face] - 1;
        }

        ui  = u + 3 * ivtx;
        uip = u + 3 * jvtx;

        l2 = _distance2 (ui, uip);
        tan_half_alpha[i] = 0.5 * sqrt(l2 / (1.0 - 0.25 * l2));
      }


      double sum_w = 0.0;
      for (i = 0; i < n_vtx_face; i++) {
        ivtx = _face_vtx[i] - 1;
        if (face_orientation[iface] < 0) {
          im = (i + 1)%n_vtx_face;
        } else {
          im = (i + n_vtx_face - 1)%n_vtx_face;
        }

        weights[ivtx] = (tan_half_alpha[im] + tan_half_alpha[i]) * inv_dist[ivtx];
        sum_w += weights[ivtx];
      }


      if (sum_w > eps) {
        for (i = 0; i < n_vtx_face; i++) {
          weights[_face_vtx[i] - 1] /= sum_w;
        }
      }

      free (inv_dist);
      free (u);
      free (tan_half_alpha);
      free (theta);
      return;
    }


    // weights
    for (i = 0; i < n_vtx_face; i++) {
      ivtx = _face_vtx[i] - 1;
      if (face_orientation[iface] < 0) {
        im = (i + 1)%n_vtx_face;
      } else {
        im = (i + n_vtx_face - 1)%n_vtx_face;
      }

      weights[ivtx] += inv_dist[ivtx] * mag_v * (tan_half_alpha[i] + tan_half_alpha[im]) /
        (sum * sin(theta[i]));
    }
  }

  // normalize weights
  double sum_w = 0.0;
  for (ivtx = 0; ivtx < n_vtx; ivtx++) {
    sum_w += weights[ivtx];
  }

  if (fabs(sum_w) > eps) {
    sum_w = 1. / sum_w;
    for (ivtx = 0; ivtx < n_vtx; ivtx++) {
      weights[ivtx] *= sum_w;
    }
  }

  free (inv_dist);
  free (u);
  free (tan_half_alpha);
  free (theta);
}


void
PDM_mean_values_polygon
(
 const int         n_vtx,
 const PDM_l_num_t poly_vtx[],
 const double      vtx_coord[],
 const double      pt_coord[],
 double            weights[]
)
{
  int ivtx, jvtx, i, ip, j, idim;
  const double eps = 1e-12;
  const double eps2 = eps * eps;

  double *dist = malloc (sizeof(double) * n_vtx);
  double *u = malloc (sizeof(double) * n_vtx * 3);
  double barycenter[3] = {0., 0., 0.};
  for (i = 0; i < n_vtx; i++) {
    // point-to-vertex vector
    double *_u = u + 3*i;

    if (poly_vtx == NULL) {
      ivtx = i;
    } else {
      ivtx = poly_vtx[i] - 1;
    }

    for (idim = 0; idim < 3; idim++) {
      _u[idim] = vtx_coord[3*ivtx + idim] - pt_coord[idim];
      barycenter[idim] += vtx_coord[3*ivtx + idim];
    }

    // distance
    double dist2 = PDM_DOT_PRODUCT (_u, _u);

    // handle special case when the point is really close to a vertex
    if (dist2 < eps2) {
      for (j = 0; j < n_vtx; j++) {
        weights[j] = 0.;
      }

      weights[i] = 1.0;
      free (dist);
      free (u);
      return;
    }


    // project onto unit sphere
    dist[i] = sqrt(dist2);
    for (idim = 0; idim < 3; idim++) {
      _u[idim] /= dist[i];
    }
  }

  for (idim = 0; idim < 3; idim++) {
    barycenter[idim] /= (float) n_vtx;
  }

  double poly_normal[3] = {0., 0., 0.};
  double tmp[3];
  double *normal = malloc (sizeof(double) * 3 * n_vtx);
  for (i = 0; i < n_vtx; i++) {
    ip = (i+1)%n_vtx;

    const double *ui  = u + 3*i;
    const double *uip = u + 3*ip;
    double *ni = normal + 3*i;

    PDM_CROSS_PRODUCT (ni, ui, uip);

    if (poly_vtx == NULL) {
      ivtx = i;
      jvtx = ip;
    } else {
      ivtx = poly_vtx[i] - 1;
      jvtx = poly_vtx[ip] - 1;
    }
    double vi[3] = {vtx_coord[3*ivtx]     - barycenter[0],
                    vtx_coord[3*ivtx + 1] - barycenter[1],
                    vtx_coord[3*ivtx + 2] - barycenter[2]};
    double vip[3] = {vtx_coord[3*jvtx]     - barycenter[0],
                     vtx_coord[3*jvtx + 1] - barycenter[1],
                     vtx_coord[3*jvtx + 2] - barycenter[2]};

    PDM_CROSS_PRODUCT (tmp, vi, vip);
    poly_normal[0] += tmp[0];
    poly_normal[1] += tmp[1];
    poly_normal[2] += tmp[2];
  }
  //printf("%f %f %f\n", poly_normal[0], poly_normal[1], poly_normal[2]);

  double *tan_half_theta = malloc (sizeof(double) * n_vtx);
  double l2;
  for (i = 0; i < n_vtx; i++) {
    const double *ui  = u + 3*i;
    const double *uip = u + 3*((i+1)%n_vtx);

    l2 = _distance2 (ui, uip);
    double denom = 1.0 - 0.25 * l2;
    if (denom < eps) {
      //tan_half_theta[i] = 0.;
      /* l = 2 <=> theta = PI <=> point on edge (i, ip) */
      ip = (i+1)%n_vtx;
      for (j = 0; j < n_vtx; j++) {
        weights[j] = 0.;
      }

      weights[i] = dist[ip] / (dist[i] + dist[ip]);
      weights[ip] = dist[i] / (dist[i] + dist[ip]);
      free (dist);
      free (u);
      free (normal);
      return;

    } else {
      tan_half_theta[i] = 0.5 * sqrt(l2 / denom);
      if (0) {//PDM_DOT_PRODUCT (poly_normal, (normal + 3*i)) < 0) {
        tan_half_theta[i] = -tan_half_theta[i];
      }
    }
  }
  free (normal);

  double sum_w = 0.0;
  for (i = 0; i < n_vtx; i++) {
    weights[i] = (tan_half_theta[(i + n_vtx - 1)%n_vtx] + tan_half_theta[i]) / dist[i];
    sum_w += weights[i];
  }

  if (sum_w > eps) {
    sum_w = 1. / sum_w;
    for (i = 0; i < n_vtx; i++) {
      weights[i] *= sum_w;
    }
  }

  free (dist);
  free (u);
  free (tan_half_theta);
}




#ifdef __cplusplus

}
#endif /* __cplusplus */
