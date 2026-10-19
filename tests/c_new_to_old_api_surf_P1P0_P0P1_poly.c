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
#include <assert.h>
#include <math.h>
#include <time.h>

#include "cwp.h"
#include "cwp_priv.h"
#include "pdm_timer.h"

#include "cwipi.h"


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
         "  -nx     <val>    Global number of vertices in the side of domaine.\n\n"
         "  -new             Use new CWIPI API.\n\n"
         "  -old             Use old CWIPI API.\n\n"
         "  -h               this message.\n\n");

  exit(exit_code);
}


/*----------------------------------------------------------------------
 *
 * Read args from the command line
 *
 * parameters:
 *   nx             --> Global number of vertices in the side of domaine
 *   new            --> New CWIPI API is used
 *   old            --> Old CWIPI API is used
 *
 *---------------------------------------------------------------------*/

static void
_read_args(int argc, char **argv, CWP_g_num_t *nx, int *new, int *old,
           int *n_compute) {
  int i = 1;
  int isNew = 0;
  int isOld = 0;

  // Parse and check command line
  while (i < argc) {
    if (strcmp(argv[i], "-h") == 0) {
      _display_usage(EXIT_SUCCESS);
    }
    else if (strcmp(argv[i], "-nx") == 0) {
      i++;
      if (i >= argc + 1) {
        _display_usage(EXIT_FAILURE);
      }
      else {
        *nx = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-old") == 0) {
      if (isNew != 0) {
        printf("Error : new CWIPI is already selected\n");
        exit(1);
      }
      isOld = 1;
      if (i >= argc + 1) {
        _display_usage(EXIT_FAILURE);
      }
      else {
        *old = 1;
        *new = 0;
      }
    }
    else if (strcmp(argv[i], "-nc") == 0) {
      i++;
      if (i >= argc) {
        _display_usage(EXIT_FAILURE);
      }
      else {
        *n_compute = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-new") == 0) {
      if (isOld != 0) {
        printf("Error : old CWIPI is already selected\n");
        exit(1);
      }
      isNew = 1;
      if (i >= argc + 1) {
        _display_usage(EXIT_FAILURE);
      }
      else {
        *new = 1;
        *old = 0;
      }
    }
    i++;
  }
}


static void
_gen_mesh
(
  const MPI_Comm       comm,
  const int            n_part,
  const CWP_g_num_t    n_vtx_seg,
  int                **pn_face,
  int                **pn_vtx,
  int               ***pface_vtx_idx,
  int               ***pface_vtx,
  double            ***pvtx_coord,
  CWP_g_num_t       ***pface_ln_to_gn,
  CWP_g_num_t       ***pvtx_ln_to_gn
)
{
  int i_rank;
  MPI_Comm_rank(comm, &i_rank);

  int          *pn_edge        = NULL;
  int         **pedge_vtx      = NULL;
  int         **pface_edge     = NULL;
  CWP_g_num_t **pedge_ln_to_gn = NULL;

  CWPT_generate_mesh_rectangle_ngon(comm,
                                    CWPT_MESH_NODAL_POLY_2D,
                                    0.,
                                    0.,
                                    0.,
                                    1.,
                                    1.,
                                    n_vtx_seg,
                                    n_vtx_seg,
                                    n_part,
                                    CWPT_SPLIT_DUAL_WITH_HILBERT,
                                    0.,
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

  for (int i_part = 0; i_part < n_part; i_part++) {
    free(pedge_vtx     [i_part]);
    free(pface_edge    [i_part]);
    free(pedge_ln_to_gn[i_part]);
  }
  free(pn_edge       );
  free(pedge_vtx     );
  free(pface_edge    );
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
  CWP_g_num_t nx = 10;
  int new = 0;
  int old = 1;
  const int nPart = 1;
  int n_compute = 10;

  _read_args(argc, argv, &nx, &new, &old, &n_compute);

  // Init + create coupling
  const char *codeName;
  int codeId;
  const char *codeCoupledName;

  assert (commWorldSize >= 2);

  if (rank < commWorldSize / 2) {
    codeName = "code1";
    codeId = 1;
    codeCoupledName = "code2";
  }
  else {
    codeName = "code2";
    codeId = 2;
    codeCoupledName = "code1";
  }

  assert ((old == 1) || (new == 1));

  MPI_Comm localComm = MPI_COMM_NULL;

  const int nb_part = 1;

  const char *cpl_name;

  if (old) {
    cwipi_init(MPI_COMM_WORLD, codeName, &localComm);

    cwipi_solver_type_t solver_type;

    if (codeId == 1) {
      solver_type = CWIPI_SOLVER_CELL_VERTEX;
    }
    else {
      solver_type = CWIPI_SOLVER_CELL_CENTER;
    }

    cpl_name = "old_cpl";
    cwipi_create_coupling(cpl_name,                                  // Coupling id
                          CWIPI_COUPLING_PARALLEL_WITH_PARTITIONING, // Coupling type
                          codeCoupledName,                           // Coupled application id
                          2,                                         // Geometric entities dimension
                          0.01,                                      // Geometric tolerance
                          CWIPI_STATIC_MESH,                         // Mesh type
                          solver_type,                               // Solver type
                          1,                                         // Postprocessing frequency
                          "EnSight Gold",                            // Postprocessing format
                          "text");                                   // Postprocessing option
  }
  else {
    const int n_code = 1;
    const CWP_Status_t is_coupled_rank = CWP_STATUS_ON;

    CWP_Init(MPI_COMM_WORLD,
             n_code,
             (const char **) &(codeName),
             is_coupled_rank,
             &localComm);

    cpl_name = "new_cpl";
    CWP_Cpl_create(codeName,                                              // Code name
                   cpl_name,                                              // Coupling id
                   codeCoupledName,                                       // Coupled application id
                   CWP_INTERFACE_SURFACE,
                   CWP_COMM_PAR_WITH_PART,                                // Coupling type
                   CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE, // Solver type
                   nb_part,                                               // Partition number
                   CWP_DYNAMIC_MESH_STATIC,                               // Mesh type
                   CWP_TIME_EXCH_USER_CONTROLLED);                        // Postprocessing frequency

    CWP_Visu_set(codeName,                // Code name
                 cpl_name,                // Coupling id
                 1,                       // Postprocessing frequency
                 CWP_VISU_FORMAT_ENSIGHT, // Postprocessing format
                 "text");                 // Postprocessing option
  }

  int          *pn_face        = NULL;
  int          *pn_vtx         = NULL;
  int         **pface_vtx_idx  = NULL;
  int         **pface_vtx      = NULL;
  double      **pvtx_coord     = NULL;
  CWP_g_num_t **pface_ln_to_gn = NULL;
  CWP_g_num_t **pvtx_ln_to_gn  = NULL;

  _gen_mesh(localComm,
            nPart,
            nx,
            &pn_face,
            &pn_vtx,
            &pface_vtx_idx,
            &pface_vtx,
            &pvtx_coord,
            &pface_ln_to_gn,
            &pvtx_ln_to_gn);

  PDM_timer_t *timer = PDM_timer_create(MPI_COMM_WORLD);

  if (old) {
    cwipi_define_mesh(cpl_name, pn_vtx[0], pn_face[0], pvtx_coord[0], pface_vtx_idx[0], pface_vtx[0]);

    MPI_Barrier(MPI_COMM_WORLD);

    for (int i = 0 ; i < n_compute ; i++) {

      PDM_timer_start(timer, "update_location", 1 /* Force synchro */);
      cwipi_update_location(cpl_name);
      PDM_timer_end(timer, "update_location", 1 /* Force synchro */);

      PDM_timer_start(timer, "locate", 1 /* Force synchro */);
      cwipi_locate(cpl_name);
      PDM_timer_end(timer, "locate", 1 /* Force synchro */);

    }
  }
  else {
    CWP_Mesh_interf_vtx_set(codeName,
                            cpl_name,
                            0,
                            pn_vtx[0],
                            pvtx_coord[0],
                            pvtx_ln_to_gn[0]);

    int block_id = CWP_Mesh_interf_block_add(codeName,
                                             cpl_name,
                                             CWP_BLOCK_FACE_POLY);

    CWP_Mesh_interf_f_poly_block_set(codeName,
                                     cpl_name,
                                     0,
                                     block_id,
                                     pn_face[0],
                                     pface_vtx_idx[0],
                                     pface_vtx[0],
                                     pface_ln_to_gn[0]);

    CWP_Mesh_interf_finalize(codeName, cpl_name);
  }

  // Fields exchange
  printf("        Exchange Code1 <-> Code2 %i\n", rank);

  double *sendValues;
  double *recvValues;

  if (codeId == 1) {
    sendValues = (double *) malloc(sizeof(double) * pn_vtx [0]);
    recvValues = (double *) malloc(sizeof(double) * pn_face[0]);
    for (int i = 0 ; i < pn_vtx [0] ; i++) {
      sendValues[i] = pvtx_coord[0][3 * i];
    }
  }

  else {
    sendValues = (double *) malloc(sizeof(double) * pn_face[0]);
    recvValues = (double *) malloc(sizeof(double) * pn_vtx [0]);
    for (int i = 0 ; i < pn_face[0] ; i++) {
      sendValues[i] = rank;
    }
  }

  // Exchange
  int n_int = 1;
  const char *fieldName1 = "cooX";

  if (old) {
    int sRequest, rRequest;
    int tag = 1;

    PDM_timer_start(timer, "timer1", 1 /* Force synchro */);

    for (int i = 0 ; i < n_int ; i++) {
      PDM_timer_start(timer, "timer2", 1 /* Force synchro */);

      if (codeId == 1) {
        cwipi_issend(cpl_name, "ech1", tag, 1, 1, 0.1, fieldName1, sendValues, &sRequest);
        cwipi_wait_issend(cpl_name, sRequest);
      }

      else {
        cwipi_irecv(cpl_name, "ech1", tag, 1, 1, 0.1, fieldName1, recvValues, &rRequest);
        cwipi_wait_irecv(cpl_name, rRequest);
      }

      PDM_timer_end(timer, "timer2", 1 /* Force synchro */);
    }

    PDM_timer_end(timer, "timer1", 1 /* Force synchro */);

  }
  else {
    CWP_Status_t visu_status = CWP_STATUS_ON;

    if (strcmp(codeName, "code1") == 0) {
      CWP_Field_create(codeName,
                       cpl_name,
                       fieldName1,
                       CWP_DOUBLE,
                       CWP_FIELD_STORAGE_INTERLEAVED,
                       1,
                       CWP_DOF_LOCATION_NODE,
                       CWP_FIELD_EXCH_SEND,
                       visu_status);

      CWP_Time_step_beg(codeName,
                        0.0);

      CWP_Field_data_set(codeName, cpl_name, fieldName1, 0, CWP_FIELD_MAP_SOURCE, sendValues);
    }
    else {
      CWP_Field_create(codeName,
                       cpl_name,
                       fieldName1,
                       CWP_DOUBLE,
                       CWP_FIELD_STORAGE_INTERLEAVED,
                       1,
                       CWP_DOF_LOCATION_NODE,
                       CWP_FIELD_EXCH_RECV,
                       visu_status);

      CWP_Time_step_beg(codeName,
                        0.0);

      CWP_Field_data_set(codeName, cpl_name, fieldName1, 0, CWP_FIELD_MAP_TARGET, recvValues);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    PDM_timer_start(timer, "timer1", 1 /* Force synchro */);

    for (int i = 0 ; i < n_compute ; i++) {

      PDM_timer_start(timer, "timer2", 1 /* Force synchro */);
      CWP_Spatial_interp_weights_compute(codeName, cpl_name);
      PDM_timer_end(timer, "timer2", 1 /* Force synchro */);
    }

    PDM_timer_end(timer, "timer1", 1 /* Force synchro */);


    for (int i = 0 ; i < n_int ; i++) {
      PDM_timer_start(timer, "timer2", 1 /* Force synchro */);

      if (strcmp(codeName, "code1") == 0) {
        CWP_Field_issend(codeName, cpl_name, fieldName1);
        CWP_Field_wait_issend(codeName, cpl_name, fieldName1);
      }
      else {
        CWP_Field_irecv(codeName, cpl_name, fieldName1);
        CWP_Field_wait_irecv(codeName, cpl_name, fieldName1);
      }
      PDM_timer_end(timer, "timer2", 1 /* Force synchro */);

    }

  }

  if (old) {
    cwipi_delete_coupling(cpl_name);
  }
  else {
    CWP_Time_step_end(codeName);
    CWP_Mesh_interf_del(codeName, cpl_name);
    CWP_Cpl_del(codeName, cpl_name);
  }
  PDM_timer_free(timer);

  // Free
  for (int i_part = 0; i_part < nPart; i_part++) {
    free(pface_vtx_idx [i_part]);
    free(pface_vtx     [i_part]);
    free(pvtx_coord    [i_part]);
    free(pface_ln_to_gn[i_part]);
    free(pvtx_ln_to_gn [i_part]);
  }
  free(pn_face       );
  free(pn_vtx        );
  free(pface_vtx_idx );
  free(pface_vtx     );
  free(pvtx_coord    );
  free(pface_ln_to_gn);
  free(pvtx_ln_to_gn );
  free(sendValues    );
  free(recvValues    );

  // Finalize
  if (old) {
    cwipi_finalize();
  }
  else {
    CWP_Finalize();
  }
  MPI_Finalize();

  return EXIT_SUCCESS;
}
