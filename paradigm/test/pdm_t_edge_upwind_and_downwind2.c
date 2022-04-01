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
_read_args
(
 int    argc,
 char **argv,
 char **filename
 )
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {

    if (strcmp(argv[i], "-h") == 0)
      _usage(EXIT_SUCCESS);

    else if (strcmp(argv[i], "-f") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *filename = argv[i];
      }
    }

    else
      _usage(EXIT_FAILURE);
    i++;
  }
}



static void
_read_mesh
(
 const char  *filename,
 int         *n_cell,
 int         *n_face,
 int         *n_edge,
 int         *n_vtx,
 int        **cell_face_idx,
 int        **cell_face,
 int        **face_vtx_idx,
 int        **face_vtx,
 int        **vtx_cell_idx,
 int        **vtx_cell,
 int        **edge_vtx,
 double     **vtx_coord,
 int        **cell_vtx_idx,
 int        **cell_vtx
 )
{
  FILE *f = fopen(filename, "r");
  assert(f != NULL);


  fscanf(f, "%d %d %d %d\n", n_vtx, n_edge, n_face, n_cell);

  *vtx_coord = malloc (sizeof(double) * (*n_vtx) * 3);
  for (int i = 0; i < *n_vtx; i++) {
    fscanf(f, "%lf %lf %lf\n",
           *vtx_coord + 3*i, *vtx_coord + 3*i+1, *vtx_coord + 3*i+2);
  }


  *cell_face_idx = malloc (sizeof(int) * (*n_cell + 1));
  for (int i = 0; i <= *n_cell; i++) {
    fscanf(f, "%d", *cell_face_idx + i);
  }

  *cell_face = malloc (sizeof(int) * (*cell_face_idx)[*n_cell]);
  for (int i = 0; i < (*cell_face_idx)[*n_cell]; i++) {
    fscanf(f, "%d", *cell_face + i);
  }


  *face_vtx_idx = malloc (sizeof(int) * (*n_face + 1));
  for (int i = 0; i <= *n_face; i++) {
    fscanf(f, "%d", *face_vtx_idx + i);
  }

  *face_vtx = malloc (sizeof(int) * (*face_vtx_idx)[*n_face]);
  for (int i = 0; i < (*face_vtx_idx)[*n_face]; i++) {
    fscanf(f, "%d", *face_vtx + i);
  }


  *edge_vtx = malloc (sizeof(int) * 2*(*n_edge));
  for (int i = 0; i < 2*(*n_edge); i++) {
    fscanf(f, "%d", *edge_vtx + i);
  }


  *vtx_cell_idx = malloc (sizeof(int) * (*n_vtx + 1));
  for (int i = 0; i <= *n_vtx; i++) {
    fscanf(f, "%d", *vtx_cell_idx + i);
  }

  *vtx_cell = malloc (sizeof(int) * (*vtx_cell_idx)[*n_vtx]);
  for (int i = 0; i < (*vtx_cell_idx)[*n_vtx]; i++) {
    fscanf(f, "%d", *vtx_cell + i);
  }


  *cell_vtx_idx = malloc (sizeof(int) * (*n_cell + 1));
  for (int i = 0; i <= *n_cell; i++) {
    fscanf(f, "%d", *cell_vtx_idx + i);
  }

  *cell_vtx = malloc (sizeof(int) * (*cell_vtx_idx)[*n_cell]);
  for (int i = 0; i < (*cell_vtx_idx)[*n_cell]; i++) {
    fscanf(f, "%d", *cell_vtx + i);
  }

  fclose(f);
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

          PDM_polygon_status_t in_poly = PDM_polygon_point_in(intersection,
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


  free(upwind_face  );
  free(downwind_face);
  free(upwind_cell  );
  free(downwind_cell);
  free(upwind_point  );
  free(downwind_point);
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

  char *filemesh = NULL;

  /*
   *  Read args
   */
  _read_args(argc,
             argv,
             &filemesh);

  assert(filemesh != NULL);

  /*
   *  Init
   */

  int i_rank;
  int n_rank;


  PDM_MPI_Init(&argc, &argv);
  PDM_MPI_Comm_rank(PDM_MPI_COMM_WORLD, &i_rank);
  PDM_MPI_Comm_size(PDM_MPI_COMM_WORLD, &n_rank);

  assert(n_rank == 1);

  int     n_cell = 0;
  int     n_face = 0;
  int     n_edge = 0;
  int     n_vtx  = 0;
  int    *cell_face_idx =  NULL;
  int    *cell_face     =  NULL;
  int    *face_vtx_idx  =  NULL;
  int    *face_vtx      =  NULL;
  int    *vtx_cell_idx  =  NULL;
  int    *vtx_cell      =  NULL;
  int    *edge_vtx      =  NULL;
  double *vtx_coord     =  NULL;
  int    *cell_vtx_idx  =  NULL;
  int    *cell_vtx      =  NULL;

  _read_mesh (filemesh,
              &n_cell,
              &n_face,
              &n_edge,
              &n_vtx,
              &cell_face_idx,
              &cell_face,
              &face_vtx_idx,
              &face_vtx,
              &vtx_cell_idx,
              &vtx_cell,
              &edge_vtx,
              &vtx_coord,
              &cell_vtx_idx,
              &cell_vtx);


  PDM_vtk_write_polydata("check_faces.vtk",
                         n_vtx,
                         vtx_coord,
                         NULL,
                         n_face,
                         face_vtx_idx,
                         face_vtx,
                         NULL,
                         NULL);

  PDM_vtk_write_std_elements("check_cells.vtk",
                             n_vtx,
                             vtx_coord,
                             NULL,
                             PDM_MESH_NODAL_TETRA4,
                             n_cell,
                             cell_vtx,
                             NULL,
                             0,
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
                                  vtx_coord);



  free(cell_face_idx);
  free(cell_face    );
  free(face_vtx_idx );
  free(face_vtx     );
  free(vtx_cell_idx );
  free(vtx_cell     );
  free(edge_vtx     );
  free(vtx_coord    );
  free(cell_vtx_idx );
  free(cell_vtx     );

  if (i_rank == 0) {
    printf("-- End\n");
    fflush(stdout);
  }

  PDM_MPI_Finalize();

  return 0;

}
