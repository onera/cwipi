#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "pdm.h"
#include "pdm_config.h"
#include "pdm_mpi.h"
#include "pdm_partitioning_algorithm.h"
#include "pdm_para_graph_dual.h"
#include "pdm_dmesh_nodal_to_dmesh.h"
#include "pdm_dmesh_nodal_elements_utils.h"
#include "pdm_dconnectivity_transform.h"
#include "pdm_part_connectivity_transform.h"
#include "pdm_part_mesh_nodal_elmts.h"
#include "pdm_dcube_nodal_gen.h"
#include "pdm_multipart.h"
#include "pdm_printf.h"
#include "pdm_sort.h"
#include "pdm_distrib.h"
#include "pdm_error.h"
#include "pdm_vtk.h"
#include "pdm_logging.h"
#include "pdm_priv.h"
#include "pdm_plugin.h"
#include "pdm_part_connectivity_transform.h"
#include "pdm_polygon.h"
#include "pdm_array.h"

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/**
 *
 * \brief  Usage
 *
 */

static void
_usage(int exit_code)
{
  PDM_printf
    ("\n"
     "  Usage: \n\n"
     "  -n      <level>  Number of vertices on the cube side.\n\n"
     "  -l      <level>  Cube length.\n\n"
     "  -n_part <level>  Number of partitions par process.\n\n"
     "  -post            Ensight outputs (only if n_part == 1). \n\n"
     "  -parmetis        Call ParMETIS.\n\n"
     "  -pt-scocth       Call PT-Scotch.\n\n"
     "  -h               This message.\n\n");

  exit(exit_code);
}


/**
 *
 * \brief  Read arguments from the command line
 *
 * \param [in]      argc     Number of arguments
 * \param [in]      argv     Arguments
 * \param [inout]   n_vtx_seg  Number of vertices on the cube side
 * \param [inout]   length   Cube length
 * \param [inout]   n_part   Number of partitions par process
 * \param [inout]   post     Ensight outputs status
 * \param [inout]   part_method Partitioner (1 ParMETIS, 2 Pt-Scotch)
 *
 */

static void
_read_args(int                    argc,
           char                 **argv,
           PDM_g_num_t           *n_vtx_seg,
           double                *length,
           int                   *n_part,
           int                   *post,
           int                   *part_method,
           PDM_Mesh_nodal_elt_t  *elt_type)
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {

    if (strcmp(argv[i], "-h") == 0)
      _usage(EXIT_SUCCESS);

    else if (strcmp(argv[i], "-n") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        long _n_vtx_seg = atol(argv[i]);
        *n_vtx_seg = (PDM_g_num_t) _n_vtx_seg;
      }
    }
    else if (strcmp(argv[i], "-l") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *length = atof(argv[i]);
    }
    else if (strcmp(argv[i], "-n_part") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *n_part = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-post") == 0) {
      *post = 1;
    }
    else if (strcmp(argv[i], "-pt-scotch") == 0) {
      *part_method = 2;
    }
    else if (strcmp(argv[i], "-parmetis") == 0) {
      *part_method = 1;
    }
    else if (strcmp(argv[i], "-t") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *elt_type = (PDM_Mesh_nodal_elt_t) atoi(argv[i]);
    }
    else
      _usage(EXIT_FAILURE);
    i++;
  }
}


static void
_compute_cell_vtx_hexa
(
 int  n_cell,
 int *cell_face,
 int *face_vtx_idx,
 int *face_vtx,
 int *cell_vtx
 )
{
  int dbg = 0;


  for (int icell = 0; icell < n_cell; icell++) {

    int *cf = cell_face + 6*icell;
    int *cv = cell_vtx  + 8*icell;
    for (int i = 0; i < 8; i++) {
      cv[i] = 0;
    }

    if (dbg) {
      log_trace("\nCell %d\n", icell);
      for (int i = 0; i < 6; i++) {
        int iface = PDM_ABS(cf[i]) - 1;
        int *fv = face_vtx + face_vtx_idx[iface];
        log_trace(" Face %6d : %6d %6d %6d %6d\n", cf[i], fv[0], fv[1], fv[2], fv[3]);
      }
    }


    // first face
    int iface = cf[0];

    if (iface < 0) {
      iface = -iface - 1;
      int *fv = face_vtx + face_vtx_idx[iface];
      assert(face_vtx_idx[iface+1] - face_vtx_idx[iface] == 4);
      for (int i = 0; i < 4; i++) {
        cv[i] = fv[i];
      }
    }

    else {
      iface = iface - 1;
      int *fv = face_vtx + face_vtx_idx[iface];
      assert(face_vtx_idx[iface+1] - face_vtx_idx[iface] == 4);
      for (int i = 0; i < 4; i++) {
        cv[i] = fv[3-i];
      }
    }

    if (dbg) {
      log_trace("first 4 vtx : %d %d %d %d\n", cv[0], cv[1], cv[2], cv[3]);
    }

    // opposite face
    int count = 0;

    for (int idx_face = 1; idx_face < 6; idx_face++) {
      int jface = PDM_ABS(cf[idx_face]) - 1;
      int sgn   = PDM_SIGN(cf[idx_face]);

      if (dbg) {
        log_trace("  face %6d\n", cf[idx_face]);
      }

      int *fv = face_vtx + face_vtx_idx[jface];
      assert(face_vtx_idx[jface+1] - face_vtx_idx[jface] == 4);

      for (int i = 0; i < 4; i++) {
        int ivtx1, ivtx2;
        if (sgn > 0) {
          ivtx1 = fv[i];
          ivtx2 = fv[(i+1)%4];
        } else {
          ivtx2 = fv[i];
          ivtx1 = fv[(i+1)%4];
        }

        if (dbg) {
          log_trace("    edge %6d %6d\n", ivtx1, ivtx2);
        }

        if (ivtx1 == cv[0] && ivtx2 == cv[1]) {
          if (sgn < 0) {
            cv[4] = fv[(i+2)%4];
            cv[5] = fv[(i+3)%4];
          } else {
            cv[4] = fv[(i+3)%4];
            cv[5] = fv[(i+2)%4];
          }
          count++;
        }

        else if (ivtx1 == cv[2] && ivtx2 == cv[3]) {
          if (sgn < 0) {
            cv[6] = fv[(i+2)%4];
            cv[7] = fv[(i+3)%4];
          } else {
            cv[6] = fv[(i+3)%4];
            cv[7] = fv[(i+2)%4];
          }
          count++;
        }
      }

      if (count == 2) {
        break;
      }

    }

    if (dbg) {
      log_trace("count = %d\n", count);
      log_trace("cv = %d %d %d %d %d %d %d %d\n",
                cv[0], cv[1], cv[2], cv[3],
                cv[4], cv[5], cv[6], cv[7]);
    }
    assert(count == 2);
  }
}


static
void
_setup_edge_upwind_and_downwind
(
 int     n_cell,
 int     n_face,
 int     n_edge,
 int     n_vtx,
 int    *cell_face_idx,
 int    *cell_face,
 int    *face_vtx_idx,
 int    *face_vtx,
 int    *vtx_cell_idx,
 int    *vtx_cell,
 int    *edge_vtx,
 double *vtx_coord
)
{
  PDM_Mesh_nodal_elt_t elt_type;
  int vtk_type_cell;
  int n_face_cell = cell_face_idx[1] - cell_face_idx[0];
  if (n_face_cell == 6) {
    elt_type = PDM_MESH_NODAL_HEXA8;
    vtk_type_cell = 12;
  } else if (n_face_cell == 4) {
    elt_type = PDM_MESH_NODAL_TETRA4;
    vtk_type_cell = 10;
  } else {
    printf("Wrong type");
    abort();
  }

  if (1) {
    PDM_vtk_write_polydata("check_faces.vtk",
                           n_vtx,
                           vtx_coord,
                           NULL,
                           n_face,
                           face_vtx_idx,
                           face_vtx,
                           NULL,
                           NULL);

    PDM_vtk_write_std_elements("check_edges.vtk",
                               n_vtx,
                               vtx_coord,
                               NULL,
                               PDM_MESH_NODAL_BAR2,
                               n_edge,
                               edge_vtx,
                               NULL,
                               0,
                               NULL,
                               NULL);
  }

  int *upwind_face   = PDM_array_const_int(n_edge, -1);
  int *downwind_face = PDM_array_const_int(n_edge, -1);
  int *upwind_cell   = PDM_array_const_int(n_edge, -1);
  int *downwind_cell = PDM_array_const_int(n_edge, -1);

  double *upwind_point   = (double *) malloc(sizeof(double) * n_edge * 3);
  double *downwind_point = (double *) malloc(sizeof(double) * n_edge * 3);

  /*
   *  Compute face centers and normals
   */
  double *face_center = (double *) malloc(sizeof(double) * n_face * 3);
  double *face_normal = (double *) malloc(sizeof(double) * n_face * 3);
  int max_n_vtx_on_face = 0;

  for (int iface = 0; iface < n_face; iface++) {

    int    *fv = face_vtx + face_vtx_idx[iface];
    double *fc = face_center + 3*iface;
    double *fn = face_normal + 3*iface;

    for (int i = 0; i < 3; i++) {
      fc[i] = fn[i] = 0;
    }

    int n_vtx_on_face = face_vtx_idx[iface+1] - face_vtx_idx[iface];
    max_n_vtx_on_face = PDM_MAX(max_n_vtx_on_face,
                                n_vtx_on_face);

    for (int i = 0; i < n_vtx_on_face; i++) {
      int ivtx = fv[i] - 1;
      for (int j = 0; j < 3; j++) {
        fc[j] += vtx_coord[3*ivtx+j];
      }
    }

    double pond = 1. / (double) n_vtx_on_face;
    for (int i = 0; i < 3; i++) {
      fc[i] *= pond;
    }


    for (int i = 0; i < n_vtx_on_face; i++) {
      int ivtx1 = fv[i] - 1;
      int ivtx2 = fv[(i+1)%n_vtx_on_face] - 1;

      double vec1[3], vec2[3];
      for (int j = 0; j < 3; j++) {
        vec1[j] = vtx_coord[3*ivtx1+j] - fc[j];
        vec2[j] = vtx_coord[3*ivtx2+j] - fc[j];
      }

      double fni[3];
      PDM_CROSS_PRODUCT(fni, vec1, vec2);

      for (int j = 0; j < 3; j++) {
        fn[j] += 0.5 * fni[j];
      }
    }
  }

  PDM_log_trace_connectivity_int(cell_face_idx, cell_face, n_cell, "cell_face : ");


  int *is_visited_face = PDM_array_zeros_int(n_face);
  int *visited_faces   = (int *) malloc(sizeof(int) * n_face);
  int n_visited_face = 0;

  double *poly_coord = (double *) malloc(sizeof(double) * max_n_vtx_on_face * 3);

  const double epsilon = 1.e-12;
  for (int iedge = 0; iedge < n_edge; iedge++) {

    log_trace("iedge = %d/%d\n", iedge, n_edge);

    int ivtx1 = edge_vtx[2*iedge  ] - 1;
    int ivtx2 = edge_vtx[2*iedge+1] - 1;

    // resest visited faces
    for (int i = 0; i < n_visited_face; i++) {
      is_visited_face[visited_faces[i]] = 0;
    }
    n_visited_face = 0;

    double edge_vec[3] = {
      vtx_coord[3*ivtx2  ] - vtx_coord[3*ivtx1  ],
      vtx_coord[3*ivtx2+1] - vtx_coord[3*ivtx1+1],
      vtx_coord[3*ivtx2+2] - vtx_coord[3*ivtx1+2],
    };

    int found[2] = {0};
    int sgn = 1;
    for (int jvtx = 0; jvtx < 2; jvtx++) {

      int ivtx = edge_vtx[2*iedge + jvtx] - 1;
      log_trace("    ivtx = %d/%d\n", ivtx, n_vtx);

      for (int idx_cell = vtx_cell_idx[ivtx]; idx_cell < vtx_cell_idx[ivtx+1]; idx_cell++) {

        int icell = PDM_ABS(vtx_cell[idx_cell]) - 1;
        log_trace("    icell = %d/%d\n", icell, n_cell);

        for (int idx_face = cell_face_idx[icell]; idx_face < cell_face_idx[icell+1]; idx_face++) {

          int iface = PDM_ABS(cell_face[idx_face]) - 1;
          log_trace("      idx_face = %d, iface = %d/%d\n", idx_face, iface, n_face);

          if (is_visited_face[iface]) {
            continue;
          }

          is_visited_face[iface] = 1;
          visited_faces[n_visited_face++] = iface;

          // eliminate faces that contain current vtx
          int has_current_vtx = 0;
          for (int idx_vtx = face_vtx_idx[iface]; idx_vtx < face_vtx_idx[iface+1]; idx_vtx++) {
            if (face_vtx[idx_vtx] - 1 == ivtx) {
              has_current_vtx = 1;
              break;
            }
          }

          if (has_current_vtx) {
            continue;
          }

          // intersect edge's line with face's plane
          int stat = -1;
          double intersection[3];
          double vec[3] = {
            face_center[3*iface  ] - vtx_coord[3*ivtx  ],
            face_center[3*iface+1] - vtx_coord[3*ivtx+1],
            face_center[3*iface+2] - vtx_coord[3*ivtx+2]
          };

          double denom = PDM_DOT_PRODUCT(edge_vec, face_normal + 3*iface);
          double numer = PDM_DOT_PRODUCT(vec, face_normal + 3*iface);

          if (PDM_ABS(denom) < epsilon) {   // Epsilon is here to avoid division by 0
            if (PDM_ABS(numer) < epsilon) { // Le edge contenu dans le plan de la face
              stat = 2;
              memcpy(intersection, vtx_coord+3*ivtx, sizeof(double) * 3);
            } else {// Edge parallèle mais pas dans le même plan
              stat = 0;
            }
          } else {  // Ca intersecte nickel
            double t = sgn * numer/denom;

            if (t < 0) {
              stat = 1;
              for (int j = 0; j < 3; j++) {
                intersection[j] = vtx_coord[3*ivtx+j] + sgn*t*edge_vec[j];
              }
            }
          }

          if (stat <= 0) {
            continue;
          }

          // We found an intersection point, now check if it is inside the face
          int *fv = face_vtx + face_vtx_idx[iface];
          int n_vtx_on_face = face_vtx_idx[iface+1] - face_vtx_idx[iface];

          double poly_bound[6] = {
            HUGE_VAL, -HUGE_VAL,
            HUGE_VAL, -HUGE_VAL,
            HUGE_VAL, -HUGE_VAL
          };
          for (int i = 0; i < n_vtx_on_face; i++) {
            int kvtx = fv[i] - 1;
            double *vc = vtx_coord + 3*kvtx;
            for (int j = 0; j < 3; j++) {
              poly_coord[3*i+j] = vc[j];
              poly_bound[2*j  ] = PDM_MIN(poly_bound[2*j  ], vc[j]);
              poly_bound[2*j+1] = PDM_MAX(poly_bound[2*j+1], vc[j]);
            }
          }

          PDM_polygon_status_t in_poly = PDM_polygon_point_in_new(intersection,
                                                                  n_vtx_on_face,
                                                                  poly_coord,
                                                                  poly_bound,
                                                                  face_normal + 3*iface);
          if (in_poly == PDM_POLYGON_INSIDE) {

            found[jvtx] = 1;
            if (jvtx == 0) {
              upwind_cell[iedge] = icell;
              upwind_face[iedge] = iface;
              memcpy(upwind_point + 3*iedge,
                     intersection,
                     sizeof(double) * 3);
            } else {
              downwind_cell[iedge] = icell;
              downwind_face[iedge] = iface;
              memcpy(downwind_point + 3*iedge,
                     intersection,
                     sizeof(double) * 3);
            }

          }

          if (found[jvtx]) {
            break;
          }

        } // end of loop on current cell's faces

        if (found[jvtx]) {
          break;
        }

      } // end of loop on cells incident to current vtx

      sgn = -sgn;
      // assert(found[jvtx]);
      log_trace("found[%d] = %d\n", jvtx, found[jvtx]);

    } // end of loop on current edge's vtx


  } // end of loop on edges

  free(face_center);
  free(face_normal);
  free(is_visited_face);
  free(visited_faces  );
  free(poly_coord);


  for (int iedge = 0; iedge < n_edge; iedge++) {
    log_trace("edge %d (%d %d), up : c %d, f %d, p (%f %f %f), down : c %d, f %d, p (%f %f %f)\n",
              iedge,
              edge_vtx[2*iedge]-1, edge_vtx[2*iedge+1]-1,
              upwind_cell[iedge], upwind_face[iedge],
              upwind_point[3*iedge], upwind_point[3*iedge+1], upwind_point[3*iedge+2],
              downwind_cell[iedge], downwind_face[iedge],
              downwind_point[3*iedge], downwind_point[3*iedge+1], downwind_point[3*iedge+2]);
  }


  // const char* field_name[] = {"upwind face", "downwind face", "upwind cell", "downwind cell", 0 };

  // const int *field[2] = {upwind_face, downwind_face};

  // PDM_vtk_write_std_elements("edges_up_down.vtk",
  //                            n_vtx,
  //                            vtx_coord,
  //                            NULL,
  //                            PDM_MESH_NODAL_BAR2,
  //                            n_edge,
  //                            edge_vtx,
  //                            NULL,
  //                            0,//2,
  //                            NULL,//field_name,
  //                            NULL);//field);

  double *line_coord = (double *) malloc(sizeof(double ) * n_edge * 6);
  for (int i = 0; i < n_edge; i++) {
    memcpy(line_coord + 6*i,     upwind_point   + 3*i, sizeof(double) * 3);
    memcpy(line_coord + 6*i + 3, downwind_point + 3*i, sizeof(double) * 3);
  }

  PDM_vtk_write_lines("edges_up_down.vtk",
                      n_edge,
                      line_coord,
                      NULL,
                      NULL);

  // PDM_vtk_write_point_cloud("upwind_points.vtk",
  //                           n_edge,
  //                           upwind_point,
  //                           NULL,
  //                           NULL);

  // PDM_vtk_write_point_cloud("downwind_points.vtk",
  //                           n_edge,
  //                           downwind_point,
  //                           NULL,
  //                           NULL);


  const char* field_name[] = {"iface", "icell", 0};

  int *field[2];

  int *connec = (int *) malloc(sizeof(int) * n_edge);
  for (int i = 0; i < n_edge; i++) {
    connec[i] = i+1;
  }

  field[0] = upwind_face;
  field[1] = upwind_cell;
  PDM_vtk_write_std_elements("upwind_points.vtk",
                             n_edge,
                             upwind_point,
                             NULL,
                             PDM_MESH_NODAL_POINT,
                             n_edge,
                             connec,
                             NULL,
                             2,
                             field_name,
              (const int **) field);

  field[0] = downwind_face;
  field[1] = downwind_cell;
  PDM_vtk_write_std_elements("downwind_points.vtk",
                             n_edge,
                             downwind_point,
                             NULL,
                             PDM_MESH_NODAL_POINT,
                             n_edge,
                             connec,
                             NULL,
                             2,
                             field_name,
              (const int **) field);

  free(connec);




  // -->> !! Only works with tetra
  int n_vtx_cell = PDM_Mesh_nodal_n_vertices_element(elt_type, 1);
  int *cell_vtx = (int *) malloc(sizeof(int) * n_vtx_cell * n_cell);

  if (elt_type == PDM_MESH_NODAL_TETRA4) {
    for (int i = 0; i < n_cell; i++) {
      int iface = cell_face[4*i];

      if (iface < 0) {
        iface = -iface - 1;
        cell_vtx[4*i    ] = face_vtx[3*iface    ];
        cell_vtx[4*i + 1] = face_vtx[3*iface + 1];
        cell_vtx[4*i + 2] = face_vtx[3*iface + 2];
      } else {
        iface = iface - 1;
        cell_vtx[4*i    ] = face_vtx[3*iface + 2];
        cell_vtx[4*i + 1] = face_vtx[3*iface + 1];
        cell_vtx[4*i + 2] = face_vtx[3*iface    ];
      }

      iface = PDM_ABS(cell_face[4*i+1]) - 1;
      for (int j = 0; j < 3; j++) {
        int ivtx = face_vtx[3*iface + j];
        if (ivtx != cell_vtx[4*i    ] &&
            ivtx != cell_vtx[4*i + 1] &&
            ivtx != cell_vtx[4*i + 2] ) {
          cell_vtx[4*i + 3] = ivtx;
        }
      }
    }
  }

  else if (elt_type == PDM_MESH_NODAL_HEXA8) {
    _compute_cell_vtx_hexa(n_cell,
                           cell_face,
                           face_vtx_idx,
                           face_vtx,
                           cell_vtx);
  }
  // <<--

  int *cell_vtx_idx = PDM_array_new_idx_from_const_stride_int(n_vtx_cell, n_cell);
  PDM_log_trace_connectivity_int(cell_vtx_idx, cell_vtx, n_cell, "cell_vtx : ");

  int n_elt = 0;
  int l_elt_vtx = 0;
  for (int i = 0; i < n_edge; i++) {
    n_elt += 1; // edge
    l_elt_vtx += 2;

    if (upwind_cell[i] >= 0) {
      n_elt += 3; // pt, face, cell
      l_elt_vtx += 1;
      l_elt_vtx += face_vtx_idx[upwind_face[i]+1] - face_vtx_idx[upwind_face[i]];
      l_elt_vtx += cell_vtx_idx[upwind_cell[i]+1] - cell_vtx_idx[upwind_cell[i]];
    }

    if (downwind_cell[i] >= 0) {
      n_elt += 3; // pt, face, cell
      l_elt_vtx += 1;
      l_elt_vtx += face_vtx_idx[downwind_face[i]+1] - face_vtx_idx[downwind_face[i]];
      l_elt_vtx += cell_vtx_idx[downwind_cell[i]+1] - cell_vtx_idx[downwind_cell[i]];
    }
  }


  printf("n_vtx = %d\n", n_vtx);
  FILE *f = fopen("visu_V4.vtk", "w");

  fprintf(f, "# vtk DataFile Version 2.0\n");
  fprintf(f, "circles\n");
  fprintf(f, "ASCII\n");
  fprintf(f, "DATASET UNSTRUCTURED_GRID\n");

  fprintf(f, "POINTS %d double\n", n_vtx + 2*n_edge);
  for (int i = 0; i < n_vtx; i++) {
    for (int j = 0; j < 3; j++) {
      fprintf(f, "%.20lf ", vtx_coord[3*i+j]);
    }
    fprintf(f, "\n");
  }
  for (int iedge = 0; iedge < n_edge; iedge++) {
    for (int j = 0; j < 3; j++) {
      fprintf(f, "%.20lf ", upwind_point[3*iedge+j]);
    }
    fprintf(f, "\n");
  }
  for (int iedge = 0; iedge < n_edge; iedge++) {
    for (int j = 0; j < 3; j++) {
      fprintf(f, "%.20lf ", downwind_point[3*iedge+j]);
    }
    fprintf(f, "\n");
  }

  fprintf(f, "CELLS %d %d\n", n_elt, n_elt + l_elt_vtx);
  for (int iedge = 0; iedge < n_edge; iedge++) {
    fprintf(f, "2 %d %d\n", edge_vtx[2*iedge]-1, edge_vtx[2*iedge+1]-1);

    if (upwind_cell[iedge] >= 0) {
      fprintf(f, "1 %d\n", n_vtx + iedge);

      int iface = upwind_face[iedge];
      fprintf(f, "%d ", face_vtx_idx[iface+1] - face_vtx_idx[iface]);
      for (int j = face_vtx_idx[iface]; j < face_vtx_idx[iface+1]; j++) {
        fprintf(f, "%d ", face_vtx[j]-1);
      }
      fprintf(f, "\n");

      int icell = upwind_cell[iedge];
      fprintf(f, "%d ", cell_vtx_idx[icell+1] - cell_vtx_idx[icell]);
      for (int j = cell_vtx_idx[icell]; j < cell_vtx_idx[icell+1]; j++) {
        fprintf(f, "%d ", cell_vtx[j]-1);
      }
      fprintf(f, "\n");

    }

    if (downwind_cell[iedge] >= 0) {
      fprintf(f, "1 %d\n", n_vtx + n_edge + iedge);

      int iface = downwind_face[iedge];
      fprintf(f, "%d ", face_vtx_idx[iface+1] - face_vtx_idx[iface]);
      for (int j = face_vtx_idx[iface]; j < face_vtx_idx[iface+1]; j++) {
        fprintf(f, "%d ", face_vtx[j]-1);
      }
      fprintf(f, "\n");

      int icell = downwind_cell[iedge];
      fprintf(f, "%d ", cell_vtx_idx[icell+1] - cell_vtx_idx[icell]);
      for (int j = cell_vtx_idx[icell]; j < cell_vtx_idx[icell+1]; j++) {
        fprintf(f, "%d ", cell_vtx[j]-1);
      }
      fprintf(f, "\n");

    }
  }


  fprintf(f, "CELL_TYPES %d\n", n_elt);
  for (int iedge = 0; iedge < n_edge; iedge++) {
    fprintf(f, "3\n"); // line

    if (upwind_cell[iedge] >= 0) {
      fprintf(f, "1\n"); // point
      fprintf(f, "7\n"); // polygon
      fprintf(f, "%d\n", vtk_type_cell); // tetra
    }

    if (downwind_cell[iedge] >= 0) {
      fprintf(f, "1\n"); // point
      fprintf(f, "7\n"); // polygon
      fprintf(f, "%d\n", vtk_type_cell); // tetra
    }

  }

  fprintf(f, "CELL_DATA %d\n", n_elt);
  // fprintf(f, "SCALARS i_edge int 1\n");
  // fprintf(f, "LOOKUP_TABLE default\n");
  // for (int iedge = 0; iedge < n_edge; iedge++) {
  //   fprintf(f, "%d\n", iedge);

  //   if (upwind_cell[iedge] >= 0) {
  //     fprintf(f, "%d\n", iedge);
  //     fprintf(f, "%d\n", iedge);
  //     fprintf(f, "%d\n", iedge);
  //   }

  //   if (downwind_cell[iedge] >= 0) {
  //     fprintf(f, "%d\n", iedge);
  //     fprintf(f, "%d\n", iedge);
  //     fprintf(f, "%d\n", iedge);
  //   }

  // }
  fprintf(f, "FIELD fields 2\n");
  fprintf(f, "i_edge 1 %d int\n", n_elt);
  for (int iedge = 0; iedge < n_edge; iedge++) {
    fprintf(f, "%d\n", iedge);
    if (upwind_cell[iedge] >= 0) {
      fprintf(f, "%d\n", iedge);
      fprintf(f, "%d\n", iedge);
      fprintf(f, "%d\n", iedge);
    }
    if (downwind_cell[iedge] >= 0) {
      fprintf(f, "%d\n", iedge);
      fprintf(f, "%d\n", iedge);
      fprintf(f, "%d\n", iedge);
    }
  }

  fprintf(f, "up_down 1 %d int\n", n_elt);
  for (int iedge = 0; iedge < n_edge; iedge++) {
    fprintf(f, "0\n");
    if (upwind_cell[iedge] >= 0) {
      fprintf(f, "1\n");
      fprintf(f, "1\n");
      fprintf(f, "1\n");
    }
    if (downwind_cell[iedge] >= 0) {
      fprintf(f, "-1\n");
      fprintf(f, "-1\n");
      fprintf(f, "-1\n");
    }
  }

  fclose(f);

  free(cell_vtx);
  free(cell_vtx_idx);


  free(line_coord  );
  free(upwind_face  );
  free(downwind_face);
  free(upwind_cell  );
  free(downwind_cell);
  free(upwind_point  );
  free(downwind_point);
}


static void
_compute_face_vtx
(
 const int   n_face,
 int        *pface_edge_idx,
 int        *pface_edge,
 int        *pedge_vtx,
 int       **pface_vtx
 )
{
  *pface_vtx = (int *) malloc(sizeof(int) * pface_edge_idx[n_face]);

  for (int i = 0; i < n_face; i++) {

    int *_pface_vtx = *pface_vtx + pface_edge_idx[i];

    int cur_vtx, next_vtx;
    int cur_edge = pface_edge[pface_edge_idx[i]];
    if (cur_edge < 0) {
      cur_edge = -cur_edge - 1;
      cur_vtx  = pedge_vtx[2*cur_edge+1];
      next_vtx = pedge_vtx[2*cur_edge  ];
    } else {
      cur_edge = cur_edge - 1;
      cur_vtx  = pedge_vtx[2*cur_edge  ];
      next_vtx = pedge_vtx[2*cur_edge+1];
    }

    for (int ivtx = 0; ivtx < pface_edge_idx[i+1] - pface_edge_idx[i]; ivtx++) {
      _pface_vtx[ivtx] = cur_vtx;

      for (int iedg = pface_edge_idx[i]; iedg < pface_edge_idx[i+1]; iedg++) {
        cur_edge = pface_edge[iedg];
        int vtx1, vtx2;
        if (cur_edge < 0) {
          cur_edge = -cur_edge - 1;
          vtx1 = pedge_vtx[2*cur_edge+1];
          vtx2 = pedge_vtx[2*cur_edge  ];
        } else {
          cur_edge = cur_edge - 1;
          vtx1 = pedge_vtx[2*cur_edge  ];
          vtx2 = pedge_vtx[2*cur_edge+1];
        }

        if (vtx1 == next_vtx) {
          cur_vtx  = next_vtx;
          next_vtx = vtx2;
          break;
        }
      }
    }

  }
}

/**
 *
 * \brief  Main
 *
 */

int main(int argc, char *argv[])
{

  /*
   *  Set default values
   */

  PDM_g_num_t        n_vtx_seg = 10;
  double             length    = 1.;
  int                n_part    = 1;
  int                post      = 0;
#ifdef PDM_HAVE_PARMETIS
  PDM_split_dual_t part_method  = PDM_SPLIT_DUAL_WITH_PARMETIS;
#else
#ifdef PDM_HAVE_PTSCOTCH
  PDM_split_dual_t part_method  = PDM_SPLIT_DUAL_WITH_PTSCOTCH;
#endif
#endif
  PDM_Mesh_nodal_elt_t elt_type  = PDM_MESH_NODAL_TETRA4;
  //  5 -> tetra
  //  6 -> pyramid
  //  7 -> prism
  //  8 -> hexa
  //  9 -> poly3d

  /*
   *  Read args
   */
  _read_args(argc,
             argv,
             &n_vtx_seg,
             &length,
             &n_part,
             &post,
     (int *) &part_method,
             &elt_type);
  /*
   *  Init
   */

  int i_rank;
  int n_rank;

  PDM_MPI_Init(&argc, &argv);
  PDM_MPI_Comm_rank(PDM_MPI_COMM_WORLD, &i_rank);
  PDM_MPI_Comm_size(PDM_MPI_COMM_WORLD, &n_rank);

  /*
   *  Create distributed cube
   */

  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;

  PDM_dcube_nodal_t* dcube = PDM_dcube_nodal_gen_create (comm,
                                                         n_vtx_seg,
                                                         n_vtx_seg,
                                                         n_vtx_seg,
                                                         length,
                                                         0.,
                                                         0.,
                                                         0.,
                                                         elt_type,
                                                         1,
                                                         PDM_OWNERSHIP_KEEP);
  PDM_dcube_nodal_gen_build (dcube);


  PDM_dmesh_nodal_t* dmn = PDM_dcube_nodal_gen_dmesh_nodal_get(dcube);
  PDM_dmesh_nodal_generate_distribution(dmn);

  if (1) {
    PDM_g_num_t *vtx_distrib = PDM_dmesh_nodal_vtx_distrib_get(dmn);
    double      *dvtx_coord  = PDM_DMesh_nodal_vtx_get(dmn);
    int dn_vtx = vtx_distrib[i_rank+1] - vtx_distrib[i_rank];

    double noise = 0.2 / (double) (n_vtx_seg - 1);
    for (int i = 0; i < 3*dn_vtx; i++) {
      dvtx_coord[i] += noise * (rand() / (double) RAND_MAX - 0.5);
    }
  }

  if(1 == 1) {
    PDM_dmesh_nodal_dump_vtk(dmn, PDM_GEOMETRY_KIND_VOLUMIC , "out_volumic");
    PDM_dmesh_nodal_dump_vtk(dmn, PDM_GEOMETRY_KIND_SURFACIC, "out_surfacic");
  }

  /*
   * Partitionnement
   */
  int n_zone = 1;
  int n_part_zones = n_part;
  PDM_multipart_t *mpart_id = PDM_multipart_create(n_zone,
               &n_part_zones,
               PDM_FALSE,
               part_method,
               PDM_PART_SIZE_HOMOGENEOUS,
               NULL,
               comm,
               PDM_OWNERSHIP_KEEP);

  PDM_multipart_set_reordering_options(mpart_id, -1, "PDM_PART_RENUM_CELL_NONE",
                                                     NULL,
                                                     "PDM_PART_RENUM_FACE_NONE");

  PDM_multipart_register_dmesh_nodal(mpart_id, 0, dmn);
  PDM_multipart_run_ppart(mpart_id);

  for (int i_zone = 0; i_zone < n_zone; i_zone++){
    for (int i_part = 0; i_part < n_part; i_part++){
      int  n_proc, tn_part;
      int  n_cell, n_face, n_vtx, n_bounds, n_joins, n_part_joins;
      int  scell_face, sface_vtx, sface_bound, sface_join;
      int  n_section;
      int *n_elt;

      PDM_multipart_part_dim_get(mpart_id,
                                 i_zone,
                                 i_part,
                                 &n_section,
                                 &n_elt,
                                 &n_cell,
                                 &n_face,
                                 &n_part_joins,
                                 &n_vtx,
                                 &n_proc,
                                 &tn_part,
                                 &scell_face,
                                 &sface_vtx,
                                 &sface_bound,
                                 &n_bounds,
                                 &sface_join,
                                 &n_joins);
      double       *vtx;
      int          *cell_face_idx, *cell_face, *face_cell, *face_vtx_idx, *face_vtx;
      int          *face_bound_idx, *face_bound, *face_join_idx, *face_join;
      int          *face_part_bound_proc_idx, *face_part_bound_part_idx, *face_part_bound;
      PDM_g_num_t  *cell_ln_to_gn, *face_ln_to_gn, *vtx_ln_to_gn, *face_bound_ln_to_gn, *face_join_ln_to_gn;
      int          *cell_tag, *face_tag, *vtx_tag;
      int         **elt_vtx_idx;
      int         **elt_vtx;
      PDM_g_num_t **elt_section_ln_to_gn;

      PDM_multipart_part_val_get(mpart_id, i_zone, i_part, &elt_vtx_idx, &elt_vtx, &elt_section_ln_to_gn,
                                 &cell_tag, &cell_face_idx, &cell_face, &cell_ln_to_gn,
                                 &face_tag, &face_cell, &face_vtx_idx, &face_vtx, &face_ln_to_gn,
                                 &face_part_bound_proc_idx, &face_part_bound_part_idx, &face_part_bound,
                                 &vtx_tag, &vtx, &vtx_ln_to_gn, &face_bound_idx, &face_bound,
                                 &face_bound_ln_to_gn, &face_join_idx, &face_join, &face_join_ln_to_gn);

      // PDM_log_trace_array_long(cell_ln_to_gn, n_cell, "cell_ln_to_gn : ");

      // PDM_log_trace_connectivity_int(cell_face_idx, cell_face, n_cell, "cell_face : ");

      int *vtx_part_bound_proc_idx = NULL;
      int *vtx_part_bound_part_idx = NULL;
      int *vtx_part_bound          = NULL;
      PDM_multipart_part_graph_comm_vtx_data_get(mpart_id,
                                                 i_zone,
                                                 i_part,
                                                 &vtx_part_bound_proc_idx,
                                                 &vtx_part_bound_part_idx,
                                                 &vtx_part_bound);

      int *face_edge     = NULL;
      int *face_edge_idx = NULL;
      int n_face2 = PDM_multipart_part_connectivity_get(mpart_id,
                                                        i_zone,
                                                        i_part,
                                                        PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                                        &face_edge,
                                                        &face_edge_idx,
                                                        PDM_OWNERSHIP_KEEP);
      assert(n_face == n_face2);
      int *edge_vtx     = NULL;
      int *edge_vtx_idx = NULL;
      int n_edge = PDM_multipart_part_connectivity_get(mpart_id,
                                                       i_zone,
                                                       i_part,
                                                       PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                                       &edge_vtx,
                                                       &edge_vtx_idx,
                                                       PDM_OWNERSHIP_KEEP);
      assert(edge_vtx_idx == NULL);
      PDM_g_num_t* edge_ln_to_gn = NULL;
      int n_edge2 = PDM_multipart_part_ln_to_gn_get(mpart_id,
                                                    i_zone,
                                                    i_part,
                                                    PDM_MESH_ENTITY_EDGE,
                                                    &edge_ln_to_gn,
                                                    PDM_OWNERSHIP_KEEP);
      assert(n_edge2 == n_edge);

      /*
       *  Compute additionnal connectivity
       */
      assert(face_vtx_idx == NULL);
      assert(face_vtx     == NULL);

      // int* tmp_edge_vtx_idx = malloc( (n_edge+1) * sizeof(int));
      // tmp_edge_vtx_idx[0] = 0;
      // for(int i = 0; i < n_edge; ++i) {
      //   tmp_edge_vtx_idx[i+1] = tmp_edge_vtx_idx[i] + 2;
      // }

      // PDM_combine_connectivity(n_face, face_edge_idx, face_edge, tmp_edge_vtx_idx, edge_vtx, &face_vtx_idx, &face_vtx);
      // free(tmp_edge_vtx_idx);
      _compute_face_vtx(n_face,
                        face_edge_idx,
                        face_edge,
                        edge_vtx,
                        &face_vtx);
      face_vtx_idx = (int *) malloc(sizeof(int) * (n_face + 1));
      memcpy(face_vtx_idx, face_edge_idx, sizeof(int) * (n_face + 1));

      // PDM_log_trace_connectivity_int(face_vtx_idx, face_vtx, n_face, "face_vtx : ");


      int* cell_vtx_idx = NULL;
      int* cell_vtx     = NULL;
      PDM_combine_connectivity(n_cell, cell_face_idx, cell_face, face_vtx_idx, face_vtx, &cell_vtx_idx, &cell_vtx);

      int* vtx_cell_idx = NULL;
      int* vtx_cell     = NULL;
      PDM_connectivity_transpose(n_cell, n_vtx, cell_vtx_idx, cell_vtx, &vtx_cell_idx, &vtx_cell);

      // if (elt_type == PDM_MESH_NODAL_HEXA8) {
      //   int *_cell_vtx = (int *) malloc(sizeof(int) * n_cell * 8);
      //   _compute_cell_vtx_hexa(n_cell,
      //                          cell_face,
      //                          face_vtx_idx,
      //                          face_vtx,
      //                          &_cell_vtx);

      //   PDM_vtk_write_std_elements("check_hexa.vtk",
      //                              n_vtx,
      //                              vtx,
      //                              NULL,
      //                              PDM_MESH_NODAL_HEXA8,
      //                              n_cell,
      //                              _cell_vtx,
      //                              NULL,
      //                              0,
      //                              NULL,
      //                              NULL);
      //   free(_cell_vtx);

      //   abort();
      // }



      _setup_edge_upwind_and_downwind(n_cell,
                                      n_face,
                                      n_edge,
                                      n_vtx,
                                      cell_face_idx,
                                      cell_face,
                                      face_vtx_idx,
                                      face_vtx,
                                      vtx_cell_idx,
                                      vtx_cell,
                                      edge_vtx,
                                      vtx);



      free(vtx_cell_idx);
      free(vtx_cell);
      free(cell_vtx_idx);
      free(cell_vtx);
      free(face_vtx_idx);
      free(face_vtx);

    }
  }


  PDM_multipart_free(mpart_id);


  PDM_dcube_nodal_gen_free(dcube);


  if (i_rank == 0) {
    printf("-- End\n");
    fflush(stdout);
  }

  PDM_MPI_Finalize();

  return 0;

}
