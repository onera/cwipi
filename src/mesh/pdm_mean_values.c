/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <float.h>
#include <math.h>
#include <stdlib.h>

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

#include "fvmc_defs.h"
#include "fvmc_triangulate.h"
#include "fvmc_point_location.h"
#include "fvmc_ho_basis.h"

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


void
computePolygonMeanValues(const int           n_dist_points,
                         const fvmc_lnum_t  *dist_locations,
                         const fvmc_coord_t *dist_coords,
                         const int          *meshConnectivityIndex,
                         const int          *meshConnectivity,
                         const double       *meshVertexCoords,
                         /* const std::vector <int>& nDistBarCoords, */
                         /* std::vector <double>& distBarCoords */
                         const int nDistBarCoords[],
                         double    *distBarCoords)

{
  /* Boucle sur les points distants */

  fvmc_coord_t coo_point_dist[3];

  /* Tableaux locaux */

  const double eps_base = 1e-10;
  double *coo_som_fac = NULL;
  double *s           = NULL;
  double *dist        = NULL;
  double *aire        = NULL;
  double *proScal     = NULL;
  int size = 0;

  for (int ipoint =  0; ipoint < n_dist_points; ipoint++ ) {

    double *_distBarCoords = distBarCoords + nDistBarCoords[ipoint];

    /* Initialisation - Copie locale */

    int isOnEdge = 0;
    int isVertex = 0;
    int ielt = dist_locations[ipoint] - 1;

    int nbr_som_fac =  meshConnectivityIndex[ielt+1] -
                       meshConnectivityIndex[ielt];
    coo_point_dist[0] = dist_coords[3*ipoint];
    coo_point_dist[1] = dist_coords[3*ipoint + 1];
    coo_point_dist[2] = dist_coords[3*ipoint + 2];

    if (ipoint == 0) {
      size = nbr_som_fac;

      coo_som_fac = malloc (sizeof(double) * size * 3);
      s           = malloc (sizeof(double) * size * 3);
      dist        = malloc (sizeof(double) * size);
      aire        = malloc (sizeof(double) * size);
      proScal     = malloc (sizeof(double) * size);
    }
    else {
      if (size < nbr_som_fac) {
        size = nbr_som_fac;

        coo_som_fac = realloc (coo_som_fac, sizeof(double) * size * 3);
        s           = realloc (s,           sizeof(double) * size * 3);
        dist        = realloc (dist,        sizeof(double) * size);
        aire        = realloc (aire,        sizeof(double) * size);
        proScal     = realloc (proScal,     sizeof(double) * size);
      }
    }

    for (int isom = 0; isom < nbr_som_fac; isom++) {
      coo_som_fac[3*isom]   =
        meshVertexCoords[3*(meshConnectivity[meshConnectivityIndex[ielt]+isom]-1)];

      coo_som_fac[3*isom+1] =
        meshVertexCoords[3*(meshConnectivity[meshConnectivityIndex[ielt]+isom]-1)+1];

      coo_som_fac[3*isom+2] =
        meshVertexCoords[3*(meshConnectivity[meshConnectivityIndex[ielt]+isom]-1)+2];
    }

    /* Projection sur un plan moyen */

    double bary[3];
    PDM_polygon_compute_barycenter (nbr_som_fac, &(coo_som_fac[0]), bary);

    double n[3] = {0, 0, 1};
    PDM_plane_normal (nbr_som_fac, &(coo_som_fac[0]), n);

    PDM_plane_projection2 (coo_point_dist, bary, n, coo_point_dist);

    for (int isom = 0; isom < nbr_som_fac; isom++) {

      double *pt1 = &(coo_som_fac[0]) + 3 *isom;
      PDM_plane_projection2 (pt1, bary, n, pt1);

    }

    double bounds[6] = {DBL_MAX, -DBL_MAX,
                        DBL_MAX, -DBL_MAX,
                        DBL_MAX, -DBL_MAX};

    for (int isom = 0; isom < nbr_som_fac; isom++) {
      bounds[0] = PDM_MIN (bounds[0], coo_som_fac[3*isom]);
      bounds[1] = PDM_MAX (bounds[1], coo_som_fac[3*isom]);

      bounds[2] = PDM_MIN (bounds[2], coo_som_fac[3*isom + 1]);
      bounds[3] = PDM_MAX (bounds[3], coo_som_fac[3*isom + 1]);

      bounds[4] = PDM_MIN (bounds[4], coo_som_fac[3*isom + 2]);
      bounds[5] = PDM_MAX (bounds[5], coo_som_fac[3*isom + 2]);
    }


    /* Verification que le point est dans l'element */

    if (PDM_polygon_point_in (coo_point_dist,
                              nbr_som_fac,
                              &(coo_som_fac[0]),
                              bounds,
                              n) != 1) {

      double closestPoint[3];
      double dist_min = DBL_MAX;

      for (int k = 0; k < nbr_som_fac; k++) {
        double *p1 = &(coo_som_fac[3 * k]);
        double *p2 = &(coo_som_fac[3 * ((k+1) % nbr_som_fac)]);
        double closest[3];
        double t;

        double dist2 = PDM_line_distance (coo_point_dist,
                                          p1,
                                          p2,
                                          &t,
                                          closest);

        if (dist2 < dist_min) {
          dist_min = dist2;
          closestPoint[0] = closest[0];
          closestPoint[1] = closest[1];
          closestPoint[2] = closest[2];
        }
      }

      coo_point_dist[0] = closestPoint[0];
      coo_point_dist[1] = closestPoint[1];
      coo_point_dist[2] = closestPoint[2];

    }

    /* Calcul des coordonnnees barycentriques */

    double min_dist = DBL_MAX;
    for (int isom = 0; isom < nbr_som_fac; isom++) {

      int inext = (isom + 1) % nbr_som_fac;
      double *vect = &s[0] + 3*isom;
      double l_edge;
      vect[0] = coo_som_fac[3*inext]   - coo_som_fac[3*isom];
      vect[1] = coo_som_fac[3*inext+1] - coo_som_fac[3*isom+1];
      vect[2] = coo_som_fac[3*inext+2] - coo_som_fac[3*isom+2];
      l_edge  = PDM_MODULE (vect);

      min_dist = PDM_MIN (l_edge, min_dist);
    }

    double eps = PDM_MAX (min_dist * eps_base, 1.e-30);

    for (int isom = 0; isom < nbr_som_fac; isom++) {

      double *vect = &s[0] + 3*isom;
      vect[0] = coo_som_fac[3*isom]   - coo_point_dist[0];
      vect[1] = coo_som_fac[3*isom+1] - coo_point_dist[1];
      vect[2] = coo_som_fac[3*isom+2] - coo_point_dist[2];
      dist[isom] = PDM_MODULE (vect);

    }

    int currentVertex;
    for (int isom = 0; isom < nbr_som_fac; isom++) {
      int inext = (isom + 1) % nbr_som_fac;
      double *vect1 = &s[0] + 3 * isom;
      double *vect2 = &s[0] + 3 * inext;
      double pvect[3];

      proScal[isom] = PDM_DOT_PRODUCT (vect1, vect2);
      PDM_CROSS_PRODUCT(pvect, vect1, vect2);

      double sign = PDM_DOT_PRODUCT (pvect, n);
      aire[isom] = PDM_MODULE(pvect);

      if (sign < 0) {
        aire[isom] = -aire[isom];
      }

      if (dist[isom] <= eps) {

        isVertex = 1;
        currentVertex = isom;
        break;
      }

      else if ((fabs(aire[isom]) <= eps)  && (proScal[isom] < 0)) {

        isOnEdge = 1;
        currentVertex = isom;
        break;

      }

    }

    /* Le point distant est un sommet */

    if (isVertex) {
      for (int isom = 0; isom < nbr_som_fac; isom++)
        _distBarCoords[isom] = 0.;
      _distBarCoords[currentVertex] = 1.;
    }

    /* Le point distant est sur arete */

    else if (isOnEdge) {

      for (int isom = 0; isom < nbr_som_fac; isom++)
        _distBarCoords[isom] = 0.;

      int nextPoint = (currentVertex + 1) % nbr_som_fac;

      _distBarCoords[currentVertex] =
        dist[nextPoint]     / (dist[nextPoint]+dist[currentVertex]);
      _distBarCoords[nextPoint]     =
        dist[currentVertex] / (dist[nextPoint]+dist[currentVertex]);

    }

    /* Cas general */

    else {

      double sigma = 0;
      for (int isom = 0; isom < nbr_som_fac; isom++) {
        double coef = 0.;
        int previousVertex = (isom - 1 + nbr_som_fac) % nbr_som_fac;
        int nextVertex = (isom + 1) % nbr_som_fac;

        if (fabs(aire[previousVertex]) > eps)
          coef += (dist[previousVertex] - proScal[previousVertex]/dist[isom]) / aire[previousVertex];
        if (fabs(aire[isom]) > eps)
          coef += (dist[nextVertex]     - proScal[isom]/dist[isom])           / aire[isom];
        sigma += coef;
        _distBarCoords[isom] = coef;

      }

      if (fabs(sigma) >= eps ) {
        for (int isom = 0; isom < nbr_som_fac; isom++) {
          _distBarCoords[isom] /= sigma;
        }
      }

      else {

        double abs_sigma = fabs(sigma);
        printf("Warning : Mise à NAN %f %f\n", abs_sigma,  eps);
        for (int isom = 0; isom < nbr_som_fac; isom++) {
          _distBarCoords[isom] = NAN;
        }
      }

      /* Check Result */

      for (int isom = 0; isom <  nbr_som_fac; isom++) {
        if ( _distBarCoords[isom] != _distBarCoords[isom] ||
             _distBarCoords[isom] < 0. ||
             _distBarCoords[isom] > 1. ) {

          double dist_min = DBL_MAX;
          int k_min = 0;
          double t_min;

          for (int k = 0; k < nbr_som_fac; k++) {
            _distBarCoords[k] = 0.0;
          }

          for (int k = 0; k < nbr_som_fac; k++) {
            double *p1 = &(coo_som_fac[3 * k]);
            double *p2 = &(coo_som_fac[3 * ((k+1) % nbr_som_fac)]);
            double closest[3];
            double t;

            double dist2 = PDM_line_distance (coo_point_dist,
                                              p1,
                                              p2,
                                              &t,
                                              closest);
            if (dist2 < dist_min) {
              t_min = t;
              k_min = k;
            }
          }

          _distBarCoords[k_min] = 1 - t_min;
          _distBarCoords[(k_min + 1) % nbr_som_fac] = t_min;

          break;

        }

      }

    }

    if (0 == 1) {
      if ((n_dist_points == 1) && (dist_locations[0] == 1)) {

        PDM_printf("coord %i %i :", ipoint+1, ielt+1);
        PDM_printf(" %12.5e %12.5e %12.5e", dist_coords[3*ipoint],
                    dist_coords[3*ipoint+1],
                    dist_coords[3*ipoint+2] );
        PDM_printf("\n");

        PDM_printf("coo b %i :", ipoint+1);
        for (int isom = 0; isom < nbr_som_fac; isom++) {
          PDM_printf(" %f", _distBarCoords[isom]);
        }
        PDM_printf("\n");
      }
    }
  }

  free (coo_som_fac);
  free (s);
  free (aire);
  free (dist);
  free (proScal);
}


//-->> duplicates from pdm_geom_elem.c --> change scope?
static const double GEOM_EPS_VOL  = 1e-9; /*!< Constant value used to compute geomtric epsilon for volume */

static const double GEOM_EPS_SURF = 1e-9; /*!< Constant value used to compute geomtric epsilon for surface */
//<<--

void
compute3DMeanValuesPoly(const double point_coords[],
                        const int    n_poly_faces,
                        const int    n_poly_vertex,
                        const int    faceDirection[],
                        const int    faceToVertexIdx[],
                        const int    faceToVertex[],
                        const double vertex_coords[],
                        const double characteristicLength,
                        const float  distElt,
                        double       distBarCoords[])
{

  //
  // Polyhedron
  //


  // Mise a jour des tableaux locaux

  /* double *coo_som_face = malloc (sizeof(double) * 3 * n_poly_vertex); */
  double *dist = malloc (sizeof(double) * n_poly_vertex);
  double *s = malloc (sizeof(double) * 3 * n_poly_vertex);
  double angle[3]; //angle[  v(i) v v(i+1) ]
  double normale[9]; //normale

  double sigma = 0;

  int isOnFace = 0;

  double eps_loc = PDM_geom_elem_geometric_epsilon (characteristicLength, GEOM_EPS_VOL);

  /**** Inialisation du tableau des coordonnees temporaires a 0 ****/

  for (int isom = 0; isom < n_poly_vertex; isom++)
    distBarCoords[isom] = 0.;

  for (int isom = 0; isom < n_poly_vertex; isom++) {
    s[3 * isom    ] = vertex_coords[3 * isom    ] - point_coords[0];
    s[3 * isom + 1] = vertex_coords[3 * isom + 1] - point_coords[1];
    s[3 * isom + 2] = vertex_coords[3 * isom + 2] - point_coords[2];
  }

  if (distElt > 1.) {

    //
    // Search clostest face
    //

    int n_vtx_max = 0;
    for (int i = 0; i < n_poly_faces; i++) {
      const int n_vtx = faceToVertexIdx[i+1] - faceToVertexIdx[i];
      if (n_vtx > n_vtx_max) {
        n_vtx_max = n_vtx;
      }
    }

    double * poly_vertex = (double *) malloc (3 * sizeof(double) * n_vtx_max);

    double dist_min = DBL_MAX;
    int iface = -1;

    for (int i = 0; i < n_poly_faces; i++) {
      const int n_vtx = faceToVertexIdx[i+1] - faceToVertexIdx[i];
      int k = 0;

      for (int j = faceToVertexIdx[i]; j < faceToVertexIdx[i+1]; j++) {
        int i_vertex = faceToVertex[j] - 1;
        poly_vertex[k++] = vertex_coords[3 * i_vertex];
        poly_vertex[k++] = vertex_coords[3 * i_vertex + 1];
        poly_vertex[k++] = vertex_coords[3 * i_vertex + 2];

      }

      double closest[3];
      double dist_face;

      PDM_polygon_evaluate_position (point_coords,
                                     n_vtx, poly_vertex, closest,
                                     &dist_face);

      if (dist_face < dist_min) {
        dist_min = dist_face;
        iface = i;
      }

    }

    free (poly_vertex);

    //
    // Compute mean value for this face if point is on face
    //

    if (iface >= 0) {

      const int n_face_vertex = faceToVertexIdx[iface + 1]
                              - faceToVertexIdx[iface];

      //
      // Copy

      int distBarCoordsFaceIdx[2];
      distBarCoordsFaceIdx[0] = 0;
      distBarCoordsFaceIdx[1] = n_face_vertex;
      double *distBarCoordsFace = malloc (sizeof(double) * n_face_vertex);
      int face_location = iface + 1;

      computePolygonMeanValues(1,
                               &face_location,
                               point_coords,
                               faceToVertexIdx,
                               faceToVertex,
                               vertex_coords,
                               distBarCoordsFaceIdx,
                               distBarCoordsFace);

      for (int j = 0; j < n_poly_vertex; j++)
        distBarCoords[j] = 0.;

      for (int j = 0; j < n_face_vertex; j++) {
        int vertex = faceToVertex[faceToVertexIdx[iface]+j] - 1;
        distBarCoords[vertex] = distBarCoordsFace[j];
      }

      /* distBarCoordsFaceIdx.clear(); */
      free (distBarCoordsFace);

    }

    else {

      //
      // If point is not in a face, get the closest vertex
      //

      double normS = sqrt(s[3*0]*s[3*0] + s[3*0+1]*s[3*0+1] + s[3*0+2]*s[3*0+2]);
      int closestVertex = 0;
      for (int isom = 1; isom < n_poly_vertex; isom++) {

        double nextNormS = sqrt (s[3*isom]*s[3*isom]
                                 + s[3*isom+1]*s[3*isom+1]
                                 + s[3*isom+2]*s[3*isom+2]);
        if (nextNormS < normS) {
          closestVertex = isom;
          normS = nextNormS;
        }
      }

      distBarCoords[closestVertex] = 1;

    }

  }

  else if (distElt > 0 ) {

    //
    // Check if point is on a face
    //

    const int   n_points  = 1;
    fvmc_lnum_t point_ids = 0;

    //
    // Search clostest face
    //

    fvmc_lnum_t face_location = -1;
    float face_distance = -1.;

    fvmc_point_dist_closest_polygon(3,
                                    n_poly_faces,
                                    faceToVertexIdx,
                                    faceToVertex,
                                    vertex_coords,
                                    n_points,
                                    &point_ids,
                                    point_coords,
                                    &face_location,
                                    &face_distance);

    //
    // Compute mean value for this face if point is on face
    //

    if (face_location > 0 && face_distance < eps_loc) {

      isOnFace = 1;

      const int face_location_idx = face_location - 1;
      const int n_face_vertex = faceToVertexIdx[face_location_idx+1]
                              - faceToVertexIdx[face_location_idx];

      //
      // Copy

      int distBarCoordsFaceIdx[2];
      distBarCoordsFaceIdx[0] = 0;
      distBarCoordsFaceIdx[1] = n_face_vertex;
      double *distBarCoordsFace = malloc (sizeof(double) * n_face_vertex);

      computePolygonMeanValues(1,
                               &face_location,
                               point_coords,
                               faceToVertexIdx,
                               faceToVertex,
                               vertex_coords,
                               distBarCoordsFaceIdx,
                               distBarCoordsFace);

      for (int j = 0; j < n_poly_vertex; j++)
        distBarCoords[j] = 0.;

      for (int j = 0; j < n_face_vertex; j++) {
        int vertex = faceToVertex[faceToVertexIdx[face_location_idx]+j] - 1;

        distBarCoords[vertex] = distBarCoordsFace[j];

        /* distBarCoordsFaceIdx.clear(); */
        free (distBarCoordsFace);

      }

    }

    //
    // General alogorithm for point in polyhedron
    //

    if (!isOnFace) {

      for (int isom = 0; isom < n_poly_vertex; isom++) {

        dist[isom] = sqrt(s[3*isom    ] * s[3*isom    ]
                          + s[3*isom + 1] * s[3*isom + 1]
                          + s[3*isom + 2] * s[3*isom + 2]);

        s[3*isom]     /= dist[isom];
        s[3*isom + 1] /= dist[isom];
        s[3*isom + 2] /= dist[isom];

      }

      //
      // Second loop on faces to commpute barycentric coordinates
      //


      int s_triangle_vertices = 9;
      int *triangle_vertices = malloc (sizeof(int) * s_triangle_vertices); //Nombre de sommets apres decoupage en triangle

      for (int iface = 0; iface < n_poly_faces; iface++) {

        const int n_vertex_fac = faceToVertexIdx[iface + 1]
                               - faceToVertexIdx[iface];

        const int ind_fac_som = faceToVertexIdx[iface];

        fvmc_triangulate_state_t *fvmc_triangulate = fvmc_triangulate_state_create(n_vertex_fac);

        int triangle_vertice_size = (n_vertex_fac-2) * 3;

        if (s_triangle_vertices < triangle_vertice_size) {
          s_triangle_vertices = triangle_vertice_size;
          triangle_vertices = realloc (triangle_vertices, sizeof(int) * s_triangle_vertices);
        }

        //
        // Face triangulation
        //

        int n_triangles;

        if (n_vertex_fac == 4) {

          n_triangles = fvmc_triangulate_quadrangle(3,
                                                    vertex_coords,
                                                    NULL,
                                                    faceToVertex + ind_fac_som,
                                                    &(triangle_vertices[0]));

        }

        else if (n_vertex_fac > 4) {

          n_triangles = fvmc_triangulate_polygon(3,
                                                 n_vertex_fac,
                                                 vertex_coords,
                                                 NULL,
                                                 faceToVertex + ind_fac_som,
                                                 FVMC_TRIANGULATE_MESH_DEF,
                                                 &(triangle_vertices[0]),
                                                 fvmc_triangulate);
        }

        else {
          n_triangles = 1;
          for (int i = 0; i < 3; i++) {
            triangle_vertices[i] = faceToVertex[ind_fac_som + i];
          }
        }

        //
        // Loop on triangles
        //

        for (int itri = 0; itri < n_triangles; itri++) {

          //
          // Check triangle surface
          //

          const int i = triangle_vertices[3*itri    ] - 1;
          int j, k;
          if (faceDirection[iface] < 0) {
            j = triangle_vertices[3*itri + 1] - 1;
            k = triangle_vertices[3*itri + 2] - 1;
          }
          else {
            j = triangle_vertices[3*itri + 2] - 1;
            k = triangle_vertices[3*itri + 1] - 1;
          }

          const double coo_ijx = vertex_coords[3*j]   - vertex_coords[3*i];
          const double coo_ijy = vertex_coords[3*j+1] - vertex_coords[3*i+1];
          const double coo_ijz = vertex_coords[3*j+2] - vertex_coords[3*i+2];
          const double coo_ikx = vertex_coords[3*k]   - vertex_coords[3*i];
          const double coo_iky = vertex_coords[3*k+1] - vertex_coords[3*i+1];
          const double coo_ikz = vertex_coords[3*k+2] - vertex_coords[3*i+2];

          const double areaTri_ijk = 0.5 * sqrt((coo_ijy * coo_ikz - coo_ijz * coo_iky)
                                                * (coo_ijy * coo_ikz - coo_ijz * coo_iky)
                                                + (coo_ijz * coo_ikx - coo_ijx * coo_ikz)
                                                * (coo_ijz * coo_ikx - coo_ijx * coo_ikz)
                                                + (coo_ijx * coo_iky - coo_ijy * coo_ikx)
                                                * (coo_ijx * coo_iky - coo_ijy * coo_ikx));

          double eps_face = PDM_geom_elem_geometric_epsilon (characteristicLength, GEOM_EPS_SURF);

          if (fabs(areaTri_ijk) > eps_face) {

            for (int isom = 0; isom < 3; isom++) {

              int isuiv;
              int iprec;
              double prod_scal;
              double mod;

              if (faceDirection[iface] < 0) {
                iprec = triangle_vertices[3*itri + (isom + 2) % 3] - 1;
                isuiv = triangle_vertices[3*itri + (isom + 1) % 3] - 1;
              }
              else {
                iprec = triangle_vertices[3*itri + (isom + 1) % 3] - 1;
                isuiv = triangle_vertices[3*itri + (isom + 2) % 3] - 1;
              }

              prod_scal = s[3*iprec    ] * s[3*isuiv    ]
                        + s[3*iprec + 1] * s[3*isuiv + 1]
                        + s[3*iprec + 2] * s[3*isuiv + 2];

              angle[isom] = acos(prod_scal); //s

              normale[3 * isom    ] =  s[3*iprec + 1] * s[3*isuiv + 2]
                                     - s[3*iprec + 2] * s[3*isuiv + 1];
              normale[3 * isom + 1] =  s[3*iprec + 2] * s[3*isuiv    ]
                                     - s[3*iprec    ] * s[3*isuiv + 2];
              normale[3 * isom + 2] =  s[3*iprec    ] * s[3*isuiv + 1]
                                     - s[3*iprec + 1] * s[3*isuiv    ];

              /// verifier norm

              mod = sqrt(normale[3*isom    ] * normale[3*isom    ]
                         + normale[3*isom + 1] * normale[3*isom + 1]
                         + normale[3*isom + 2] * normale[3*isom + 2]);

              if (mod <  eps_face) {
                normale[3*isom    ] = 0.;
                normale[3*isom + 1] = 0.;
                normale[3*isom + 2] = 0.;
              }

              else {

                normale[3*isom    ] /= mod;
                normale[3*isom + 1] /= mod;
                normale[3*isom + 2] /= mod;
              }
            }

            for (int isom = 0; isom < 3; isom++) {

              double ps_nij_njk; //a ameliorer
              double ps_nki_njk; //a ameliorer
              double ps_ei_njk;  //a ameliorer

              const int iprec = (isom + 2) % 3;
              const int isuiv = (isom + 1) % 3;

              ps_nij_njk = normale[3 * isom    ] * normale[3 * isuiv    ]
                + normale[3 * isom + 1] * normale[3 * isuiv + 1]
                + normale[3 * isom + 2] * normale[3 * isuiv + 2];

              ps_nki_njk = normale[3 * isom    ] * normale[3 * iprec    ]
                + normale[3 * isom + 1] * normale[3 * iprec + 1]
                + normale[3 * isom + 2] * normale[3 * iprec + 2];

              // ps_ei_njk --> sur la face


              const int ivertex_tri = triangle_vertices[3*itri + isom] - 1;
              ps_ei_njk = s[3*ivertex_tri    ] * normale[3*isom]
                        + s[3*ivertex_tri + 1] * normale[3*isom + 1]
                        + s[3*ivertex_tri + 2] * normale[3*isom + 2];

              // vérifier ps_ei_njk

              if (fabs(ps_ei_njk) >  eps_face) {
                distBarCoords[ivertex_tri] +=
                  (angle[isom] + angle[isuiv] * ps_nij_njk + angle[iprec] * ps_nki_njk)
                  / (2 * ps_ei_njk);
               }

            } // Loop on vertices

          } // Good triangle

        } // Loop on triangles

        fvmc_triangulate_state_destroy(fvmc_triangulate);

      } // Loop on faces

      for (int isom = 0; isom < n_poly_vertex; isom++) {

        distBarCoords[isom] /= dist[isom];
        sigma += distBarCoords[isom];

      }

      for (int isom = 0; isom < n_poly_vertex; isom++)
        distBarCoords[isom] = distBarCoords[isom] / sigma;

    } // End of general algorithm (if (!isonface))

    //
    // Output results
    //

    if (0 == 1) {

      double test[3];

      for (int i = 0; i < 3; i++)
        test[i] = 0;

      for (int isom = 0; isom < n_poly_vertex; isom++){

        test[0] += distBarCoords[isom] * vertex_coords[3*isom];
        test[1] += distBarCoords[isom] * vertex_coords[3*isom + 1];
        test[2] += distBarCoords[isom] * vertex_coords[3*isom + 2];

      }

      PDM_printf("point distant | verification \n");

      double dd = 0;
      for (int i = 0; i < 3; i++) {
        PDM_printf("  %f       |    %f \n",point_coords[i],test[i]);
        dd += (point_coords[i] - test[i]) * (point_coords[i] - test[i]);
      }

      if (sqrt(dd) > 1e-3)
        PDM_printf(" !!!! Erreur sur les coordonnees baryc directionf: %12.5e %i !!!!\n",sqrt(dd), isOnFace);
      else
        PDM_printf(" ++++ ok                                         : %12.5e %i ++++\n",sqrt(dd), isOnFace);

      PDM_printf("coord :");
      PDM_printf(" %12.5e %12.5e %12.5e", point_coords[0],
                  point_coords[1],
                  point_coords[2] );
      PDM_printf("\n");

      PDM_printf("coo b :");
      for (int isom = 0; isom < n_poly_vertex; isom++)
        PDM_printf(" %f", distBarCoords[isom]);

      PDM_printf("\n");

    }

  }

}



/*=============================================================================
 * Public function definition
 *============================================================================*/



#ifdef __cplusplus
}
#endif /* __cplusplus */
