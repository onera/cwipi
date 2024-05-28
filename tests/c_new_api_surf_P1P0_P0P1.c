/*
  This file is part of the CWIPI library.

  Copyright (C) 2021-2023  ONERA

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library. If not, see <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "cwp.h"
#include "cwp_priv.h"
#include "grid_mesh.h"

#define ABS(a)   ((a) <  0  ? -(a) : (a))
#define MAX(a,b) ((a) > (b) ?  (a) : (b))

/*----------------------------------------------------------------------
 *
 * Display usage
 *
 * parameters:
 *   exit code           <-- Exit code
 *---------------------------------------------------------------------*/

static void
_display_usage(int exit_code) {
  printf("\n"
         "  Usage: \n\n"
         "  -n     <level>  Number of vertices in band width.\n\n"
         "  -rand  <level>  Random level ( > 0 and < 0.4) \n\n"
         "  -h             this message.\n\n");
  exit(exit_code);
}


/*----------------------------------------------------------------------
 *
 * Read args from the command line
 *
 * parameters:
 *   nVertex             <-- Number of vertices in bandwidth
 *   randLevel           <-- Random level
 *---------------------------------------------------------------------*/

static void
_read_args(int argc, char **argv, int *nVertex, double *randLevel) {
  int i = 1;

  while (i < argc) {
    if (strcmp(argv[i], "-h") == 0) {_display_usage(EXIT_SUCCESS);}
    else if (strcmp(argv[i], "-n") == 0) {
      i++;
      if (i >= argc) {_display_usage(EXIT_FAILURE);}
      else {*nVertex = atoi(argv[i]);}
    }
    else if (strcmp(argv[i], "-rand") == 0) {
      i++;
      if (i >= argc) {_display_usage(EXIT_FAILURE);}
      else {*randLevel = atof(argv[i]);}
    }
    i++;
  }
}


static void
_gen_mesh
(
 const MPI_Comm            comm,
 const int                 n_part,
 const CWPT_split_dual_t   part_method,
 const CWP_g_num_t         n_vtx_seg,
 const int                 randomize,
 const int                 random_seed,
 int                     **pn_face,
 int                     **pn_vtx,
 int                    ***pface_vtx_idx,
 int                    ***pface_vtx,
 double                 ***pvtx_coord,
 CWP_g_num_t            ***pface_ln_to_gn,
 CWP_g_num_t            ***pvtx_ln_to_gn
 )
{
  CWP_UNUSED(randomize);
  CWP_UNUSED(random_seed);

  int          *pn_edge        = NULL;
  int         **pface_edge     = NULL;
  int         **pedge_vtx      = NULL;
  CWP_g_num_t **pedge_ln_to_gn = NULL;
  CWPT_generate_mesh_rectangle_ngon(comm,
                                    (CWPT_Mesh_nodal_elt_t) CWP_BLOCK_FACE_QUAD4,
                                    0.,
                                    0.,
                                    0.,
                                    1.,
                                    1.,
                                    n_vtx_seg,
                                    n_vtx_seg,
                                    n_part,
                                    part_method,
                                    0,
                                    pn_vtx,
                                    &pn_edge,
                                    pn_face,
                                    pvtx_coord,
                                    &pedge_vtx,
                                    pface_vtx_idx,
                                    &pface_edge,
                                    pface_vtx,
                                    pvtx_ln_to_gn,
                                    &pedge_ln_to_gn,
                                    pface_ln_to_gn);

  for (int i = 0; i < n_part; i++) {
    free(pface_edge    [i]);
    free(pedge_vtx     [i]);
    free(pedge_ln_to_gn[i]);
  }
  free(pn_edge);
  free(pface_edge    );
  free(pedge_vtx     );
  free(pedge_ln_to_gn);
}


/*----------------------------------------------------------------------
 *
 * Main : surface coupling test : P1P0_P0P1
 *
 *---------------------------------------------------------------------*/

int
main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);

  int rank;
  int commWorldSize;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &commWorldSize);

  // Read args from command line
  int nVertexSeg = 10;
  double randLevel = 0.4;

  _read_args(argc, argv, &nVertexSeg, &randLevel);

  // Initialization
  int n_code = 1;
  const char **code_name = malloc(sizeof(char *) * n_code);
  const char **coupled_code_name = malloc(sizeof(char *) * n_code);
  CWP_Status_t is_active_rank = CWP_STATUS_ON;
  MPI_Comm *intra_comm = malloc(sizeof(MPI_Comm) * n_code);
  MPI_Comm *connectableLocalComm = malloc(sizeof(MPI_Comm) * n_code);
  int *connectableLocalCommSize = malloc(sizeof(int) * n_code);

  if (rank % 2 == 0) {
    printf("%d - Working for code1\n", rank);
    code_name[0] = "code1";
    coupled_code_name[0] = "code2";
  }

  else {
    printf("%d - Working for code2\n", rank);
    code_name[0] = "code2";
    coupled_code_name[0] = "code1";
  }

  CWP_Init(MPI_COMM_WORLD,
           n_code,
           (const char **) code_name,
           is_active_rank,
           intra_comm);

  printf("%d - Create coupling\n", rank);
  const char *cpl_name = "c_new_api_surf_P1P0_P0P1";
  int nb_part = 1;
  // CWP_Spatial_interp_t spatial_interp = CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE;
  // CWP_Spatial_interp_t spatial_interp = CWP_SPATIAL_INTERP_FROM_NEAREST_SOURCES_LEAST_SQUARES;
  CWP_Spatial_interp_t spatial_interp = CWP_SPATIAL_INTERP_FROM_INTERSECTION;
  CWP_Cpl_create(code_name[0],                                          // Code name
                 cpl_name,                                              // Coupling id
                 coupled_code_name[0],                                  // Coupled application id
                 CWP_INTERFACE_SURFACE,
                 CWP_COMM_PAR_WITH_PART,                                // Coupling type
                 spatial_interp,
                 nb_part,                                               // Number of partitions
                 CWP_DYNAMIC_MESH_STATIC,                               // Mesh type
                 CWP_TIME_EXCH_USER_CONTROLLED);                        // Postprocessing frequency

  printf("%d - Set visu\n", rank);
  CWP_Visu_set(code_name[0],            // Code name
               cpl_name,                // Coupling id
               1,                       // Postprocessing frequency
               CWP_VISU_FORMAT_ENSIGHT, // Postprocessing format
               "text");                 // Postprocessing option

  // Mesh definition
  printf("%d - Create mesh\n", rank);

  int nVertex;
  double **coords = NULL;
  int nElts;
  int **eltsConnecPointer = NULL;
  int **eltsConnec = NULL;

  // Domain bounds
  const double xmin = -10;
  const double xmax = 10;
  const double ymin = -10;
  const double ymax = 10;

  nVertex = nVertexSeg * nVertexSeg;
  nElts = (nVertexSeg - 1) * (nVertexSeg - 1);

  coords            = (double **) malloc(sizeof(double *) * n_code);
  eltsConnecPointer = (int    **) malloc(sizeof(int    *) * n_code);
  eltsConnec        = (int    **) malloc(sizeof(int    *) * n_code);


  randLevel = 0.4;
  int randSeed = 0;

  connectableLocalComm[0] = CWP_Connectable_comm_get((char *) code_name[0]);
  MPI_Comm_size(connectableLocalComm[0], &connectableLocalCommSize[0]);

  srand(0);//time(NULL));


  if (strcmp(code_name[0], "code2") == 0) {
    randLevel = 0.2;
    randSeed  = 1;
  }

  if (0) {
    int          *pn_face        = NULL;
    int          *pn_vtx         = NULL;
    int         **pface_vtx_idx  = NULL;
    int         **pface_vtx      = NULL;
    double      **pvtx_coord     = NULL;
    CWP_g_num_t **pface_ln_to_gn = NULL;
    CWP_g_num_t **pvtx_ln_to_gn  = NULL;
    _gen_mesh(intra_comm[0],
              nb_part,
              CWPT_SPLIT_DUAL_WITH_HILBERT,
              nVertexSeg,
              1,        // Warning : unused
              randSeed, // Warning : unused
              &pn_face,
              &pn_vtx,
              &pface_vtx_idx,
              &pface_vtx,
              &pvtx_coord,
              &pface_ln_to_gn,
              &pvtx_ln_to_gn);
    nVertex              = pn_vtx       [0];
    nElts                = pn_face      [0];
    coords           [0] = pvtx_coord   [0];
    eltsConnecPointer[0] = pface_vtx_idx[0];
    eltsConnec       [0] = pface_vtx    [0];
    free(pvtx_ln_to_gn [0]);
    free(pface_ln_to_gn[0]);
    free(pvtx_ln_to_gn );
    free(pface_ln_to_gn);
    free(pn_face       );
    free(pn_vtx        );
    free(pface_vtx_idx );
    free(pface_vtx     );
    free(pvtx_coord    );
  }
  else {
    coords[0]            = (double *) malloc(sizeof(double) * 3 * nVertex);
    eltsConnecPointer[0] = (int    *) malloc(sizeof(int   ) * (nElts + 1));
    eltsConnec[0]        = (int    *) malloc(sizeof(int   ) * 4 * nElts);

    grid_mesh(xmin,
              xmax,
              ymin,
              ymax,
              randLevel,
              nVertexSeg,
              (int) sqrt(connectableLocalCommSize[0]),
              coords[0],
              eltsConnecPointer[0],
              eltsConnec[0],
              connectableLocalComm[0]);
  }

  printf("%d - Number of vertex   : %d\n", rank, nVertex);
  printf("%d - Number of elements : %d\n", rank, nElts);

  CWP_Mesh_interf_vtx_set(code_name[0], cpl_name, 0, nVertex, coords[0], NULL);

  int block_id = CWP_Mesh_interf_block_add(code_name[0], cpl_name, CWP_BLOCK_FACE_POLY);

  CWP_Mesh_interf_f_poly_block_set(code_name[0],
                                   cpl_name,
                                   0,
                                   block_id,
                                   nElts,
                                   eltsConnecPointer[0],
                                   eltsConnec[0],
                                   NULL);

  CWP_Mesh_interf_finalize(code_name[0], cpl_name);

  printf("%d - Exchange Code1 <-> Code2\n", rank);

  double **sendValues = (double **) malloc(sizeof(double *) * n_code);
  double **recvValues = (double **) malloc(sizeof(double *) * n_code);

  if (strcmp(code_name[0], "code1") == 0) {
    sendValues[0] = (double *) malloc(sizeof(double) * nVertex);
    recvValues[0] = (double *) malloc(sizeof(double) * nElts);
    for (int i = 0 ; i < nVertex ; i++) {
      sendValues[0][i] = coords[0][3 * i];
    }
  }
  else {
    sendValues[0] = (double *) malloc(sizeof(double) * nElts);
    recvValues[0] = (double *) malloc(sizeof(double) * nVertex);
    for (int i = 0 ; i < nElts ; i++) {
      sendValues[0][i] = (double) rand() / (double) RAND_MAX;//rank;
    }
  }

  // Exchange
  // field1: code1 -> code2
  // field2: code2 -> code1
  const char *field_name1 = "cooX";
  const char *field_name2 = "rank";

  CWP_Status_t visu_status = CWP_STATUS_ON;
  printf("%d - Defining fields\n", rank);
  if (strcmp(code_name[0], "code1") == 0) {
    if (spatial_interp != CWP_SPATIAL_INTERP_FROM_INTERSECTION) { // TODO: remove if when intersection implemented for node-based fields
      CWP_Field_create(code_name[0],
                       cpl_name,
                       field_name1,
                       CWP_DOUBLE,
                       CWP_FIELD_STORAGE_INTERLEAVED,
                       1,
                       CWP_DOF_LOCATION_NODE,
                       CWP_FIELD_EXCH_SEND,
                       visu_status);
    }
    CWP_Field_create(code_name[0],
                     cpl_name,
                     field_name2,
                     CWP_DOUBLE,
                     CWP_FIELD_STORAGE_INTERLEAVED,
                     1,
                     CWP_DOF_LOCATION_CELL_CENTER,
                     CWP_FIELD_EXCH_RECV,
                     visu_status);

    CWP_Time_step_beg(code_name[0],
                      0.0);

    if (spatial_interp != CWP_SPATIAL_INTERP_FROM_INTERSECTION) {
      CWP_Field_data_set(code_name[0], cpl_name, field_name1, 0, CWP_FIELD_MAP_SOURCE, sendValues[0]);
    }
    CWP_Field_data_set(code_name[0], cpl_name, field_name2, 0, CWP_FIELD_MAP_TARGET, recvValues[0]);
  }
  
  else if (strcmp(code_name[0], "code2") == 0) {
    if (spatial_interp != CWP_SPATIAL_INTERP_FROM_INTERSECTION) {
      CWP_Field_create(code_name[0],
                       cpl_name,
                       field_name1,
                       CWP_DOUBLE,
                       CWP_FIELD_STORAGE_INTERLEAVED,
                       1,
                       CWP_DOF_LOCATION_NODE,
                       CWP_FIELD_EXCH_RECV,
                       visu_status);
    }
    CWP_Field_create(code_name[0],
                     cpl_name,
                     field_name2,
                     CWP_DOUBLE,
                     CWP_FIELD_STORAGE_INTERLEAVED,
                     1,
                     CWP_DOF_LOCATION_CELL_CENTER,
                     CWP_FIELD_EXCH_SEND,
                     visu_status);

    CWP_Time_step_beg(code_name[0],
                      0.0);

    CWP_Field_data_set(code_name[0], cpl_name, field_name2, 0, CWP_FIELD_MAP_SOURCE, sendValues[0]);
    if (spatial_interp != CWP_SPATIAL_INTERP_FROM_INTERSECTION) {
      CWP_Field_data_set(code_name[0], cpl_name, field_name1, 0, CWP_FIELD_MAP_TARGET, recvValues[0]);
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);
  printf("%d - Before compute\n", rank);
  // CWP_next_recv_time_set(code_name[0],
  //                        cpl_name,
  //                        0.);
  CWP_Spatial_interp_weights_compute(code_name[0], cpl_name);

  int n_uncomputed = 0;

  if (strcmp(code_name[0], "code1") == 0) {
    n_uncomputed = CWP_N_uncomputed_tgts_get (code_name[0], cpl_name, field_name2, 0);
  }

  else if (strcmp(code_name[0], "code2") == 0) {
    if (spatial_interp != CWP_SPATIAL_INTERP_FROM_INTERSECTION) {
      n_uncomputed = CWP_N_uncomputed_tgts_get (code_name[0], cpl_name, field_name1, 0);
    }
  }

  printf("%d - After compute %d\n", rank, n_uncomputed);

  if (strcmp(code_name[0], "code1") == 0) {
    if (spatial_interp != CWP_SPATIAL_INTERP_FROM_INTERSECTION) {
      CWP_Field_issend(code_name[0], cpl_name, field_name1);
    }
    CWP_Field_irecv(code_name[0], cpl_name, field_name2);
  }
  else if (strcmp(code_name[0], "code2") == 0) {
    if (spatial_interp != CWP_SPATIAL_INTERP_FROM_INTERSECTION) {
      CWP_Field_irecv(code_name[0], cpl_name, field_name1);
    }
    CWP_Field_issend(code_name[0], cpl_name, field_name2);
  }

  if (strcmp(code_name[0], "code1") == 0) {
    if (spatial_interp != CWP_SPATIAL_INTERP_FROM_INTERSECTION) {
      CWP_Field_wait_issend(code_name[0], cpl_name, field_name1);
    }
    CWP_Field_wait_irecv(code_name[0], cpl_name, field_name2);
  }
  else if (strcmp(code_name[0], "code2") == 0) {
    if (spatial_interp != CWP_SPATIAL_INTERP_FROM_INTERSECTION) {
      CWP_Field_wait_irecv(code_name[0], cpl_name, field_name1);
    }
    CWP_Field_wait_issend(code_name[0], cpl_name, field_name2);
  }

  int check = EXIT_SUCCESS;
  if (spatial_interp == CWP_SPATIAL_INTERP_FROM_INTERSECTION) {
    /* Check mass conservation */
    double *field_value = recvValues[0];
    int     icode       = -1;
    if (strcmp(code_name[0], "code1") == 0) {
      field_value = recvValues[0];
      icode = 0;
    }
    else if (strcmp(code_name[0], "code2") == 0) {
      field_value = sendValues[0];
      icode = 1;
    }

    // Compute polygon area and mass
    int    n_vtx_face = 0;
    double area       = 0.;
    double face_center[3];
    double edge_center[3];
    double v1v2[3];
    double vectFECenter[3];
    double l_integral[2] = {0., 0.};
    for (int i_face = 0; i_face < nElts; i_face++) {

      area = 0;
      // Compute face center
      for (int j = 0; j < 3; j++) {
        face_center[j]  = 0.;
      }

      n_vtx_face = eltsConnecPointer[0][i_face+1]-eltsConnecPointer[0][i_face];
      for (int i_vtx = 0; i_vtx < n_vtx_face; i_vtx++) {
        int vtx_id = eltsConnec[0][eltsConnecPointer[0][i_face] + i_vtx]-1;

        for (int j = 0; j < 3; j++) {
          face_center[j] += coords[0][3*vtx_id + j];
        }
      }

      for (int j = 0; j < 3; j++) {
        face_center[j] /= n_vtx_face;
      }

      // Compute area
      for (int i_vtx = 0; i_vtx < n_vtx_face; i_vtx++) {
        int vtx_id = eltsConnec[0][eltsConnecPointer[0][i_face] + i_vtx]-1;
        int vtx_jd = eltsConnec[0][eltsConnecPointer[0][i_face] + (i_vtx+1)%n_vtx_face]-1;

        for (int j = 0; j < 3; j++) {
          edge_center[j] = 0.5 * (coords[0][3*vtx_jd + j] + coords[0][3*vtx_id + j]);
          v1v2[j]        = coords[0][3*vtx_jd + j] - coords[0][3*vtx_id + j];
        }

        for (int j = 0; j < 3; j++) {
          vectFECenter[j] = edge_center[j] - face_center[j];
        }

        double surface_vectorTria[3];
        surface_vectorTria[0] = vectFECenter[1] * v1v2[2] - v1v2[1] * vectFECenter[2];
        surface_vectorTria[1] = v1v2[0] * vectFECenter[2] - vectFECenter[0] * v1v2[2];
        surface_vectorTria[2] = vectFECenter[0] * v1v2[1] - v1v2[0] * vectFECenter[1];

        for (int i = 0; i < 3; i++)
          surface_vectorTria[i] *= 0.5;

        const double areaTri = sqrt(surface_vectorTria[0] * surface_vectorTria[0] + surface_vectorTria[1] * surface_vectorTria[1] + surface_vectorTria[2] * surface_vectorTria[2]);

        area += areaTri;
        // area += areaTri; // TO DO: twice??
      }

      l_integral[icode] += field_value[i_face] * area;
    }

    double g_integral[2];
    MPI_Allreduce(l_integral, g_integral, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    double rel_error = ABS(g_integral[0] - g_integral[1])/ABS(g_integral[1]);

    if (rank == 0) {
      printf("g_integral = %20.16e / %20.16e, error = %f%%\n",
             g_integral[0], g_integral[1], 100.*rel_error);
    }

    if (rel_error > 1e-12) {
      check = EXIT_FAILURE;
    }
  }

  CWP_Time_step_end(code_name[0]);

  printf("%d - Delete mesh\n", rank);
  CWP_Mesh_interf_del(code_name[0], cpl_name);

  printf("%d - Delete coupling\n", rank);
  CWP_Cpl_del(code_name[0], cpl_name);

  // Freeing memory
  free(coords[0]           );
  free(eltsConnecPointer[0]);
  free(eltsConnec[0]       );

  free(coords);
  free(eltsConnecPointer);
  free(eltsConnec       );

  free(sendValues[0]);
  free(recvValues[0]);
  free(sendValues);
  free(recvValues);

  free(code_name);
  free(coupled_code_name);
  free(intra_comm);
  free(connectableLocalComm);
  free(connectableLocalCommSize);


  // Finalize
  CWP_Finalize();
  MPI_Finalize();
  return check;
}
