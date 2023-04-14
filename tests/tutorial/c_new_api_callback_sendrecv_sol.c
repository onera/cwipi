/*
  This file is part of the CWIPI library.

  Copyright (C) 2023  ONERA

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
#include <assert.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include "cwp.h"
#include "cwp_priv.h"

#include "pdm_logging.h"
#include "pdm_generate_mesh.h"

/*----------------------------------------------------------------------
 *
 * Main : advanced test : Callback (code 1)
 *
 *---------------------------------------------------------------------*/

int
main
(
 int   argc,
 char *argv[]
 )
{
  // Initialize MPI :
  MPI_Init(&argc, &argv);
  int i_rank;
  int n_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &i_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &n_rank);

  // Initialize CWIPI :
  int n_code = 1;

  const char  **code_name      = malloc(sizeof(char *) * n_code);
  CWP_Status_t *is_active_rank = malloc(sizeof(CWP_Status_t) * n_code);
  double       *time_init      = malloc(sizeof(double) * n_code);
  MPI_Comm     *intra_comm     = malloc(sizeof(MPI_Comm) * n_code);

  code_name[0]      = "code 1";
  is_active_rank[0] = CWP_STATUS_ON;
  time_init[0]      = 0.;

  CWP_Init(MPI_COMM_WORLD,
           n_code,
           (const char **) code_name,
           is_active_rank,
           time_init,
           intra_comm);

  // Create the coupling with code 2 :
  int n_part = 1;
  const char  *coupling12_name     = "coupling_12";
  const char **coupled_code_name = malloc(sizeof(char *) * n_code);
  coupled_code_name[0] = "code 2";
  CWP_Cpl_create(code_name[0],
                 coupling12_name,
                 coupled_code_name[0],
                 CWP_INTERFACE_SURFACE,
                 CWP_COMM_PAR_WITH_PART,
                 CWP_SPATIAL_INTERP_FROM_CLOSEST_SOURCES_LEAST_SQUARES,
                 n_part,
                 CWP_DYNAMIC_MESH_STATIC,
                 CWP_TIME_EXCH_USER_CONTROLLED);

  // Create mesh :
  int     n_vtx = 0;
  int     n_elt = 0;
  double *coords      = NULL;
  int    *elt_vtx_idx = NULL;
  int    *elt_vtx     = NULL;
  PDM_generate_mesh_rectangle_simplified(PDM_MPI_mpi_2_pdm_mpi_comm((void *) &intra_comm[0]),
                                         10,
                                         &n_vtx,
                                         &n_elt,
                                         &coords,
                                         &elt_vtx_idx,
                                         &elt_vtx);

  // Set mesh :
  CWP_Mesh_interf_vtx_set(code_name[0],
                          coupling12_name,
                          0,
                          n_vtx,
                          coords,
                          NULL);

  int block_id = CWP_Mesh_interf_block_add(code_name[0],
                                           coupling12_name,
                                           CWP_BLOCK_FACE_POLY);

  CWP_Mesh_interf_f_poly_block_set(code_name[0],
                                   coupling12_name,
                                   0,
                                   block_id,
                                   n_elt,
                                   elt_vtx_idx,
                                   elt_vtx,
                                   NULL);

  CWP_Mesh_interf_finalize(code_name[0],
                           coupling12_name);

  // Create and set field :
  CWP_Field_create(code_name[0],
                   coupling12_name,
                   send_field_name,
                   CWP_DOUBLE,
                   CWP_FIELD_STORAGE_INTERLACED,
                   n_components,
                   CWP_DOF_LOCATION_NODE,
                   CWP_FIELD_EXCH_SEND,
                   CWP_STATUS_ON);

  CWP_Field_data_set(code_name[0],
                     coupling12_name,
                     send_field_name,
                     0,
                     CWP_FIELD_MAP_SOURCE,
                     send_field_data);

  CWP_Field_data_set(code_name[0],
                     coupling12_name,
                     send_field_name,
                     0,
                     CWP_FIELD_MAP_SOURCE,
                     send_field_data);

  // Exchange

  // Delete the coupling :
  CWP_Cpl_del(code_name[0],
              coupling12_name);

  // free

  // Finalize CWIPI :
  CWP_Finalize();

  // Finalize MPI :
  MPI_Finalize();

  return EXIT_SUCCESS;

}
