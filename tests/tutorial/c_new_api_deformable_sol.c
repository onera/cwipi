/*
  This file is part of the CWIPI library.

  Copyright (C) 2011  ONERA

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
#include <time.h>
#include <math.h>

#include "cwp.h"
#include "pdm_error.h"
#include "pdm_priv.h"
#include "pdm_generate_mesh.h"

/*----------------------------------------------------------------------
 *
 * User interpolation function
 *
 *---------------------------------------------------------------------*/

static void
_user_interpolation_function
(
 const char           *local_code_name,
 const char           *cpl_id,
 const char           *field_id,
 int                   i_part,
 CWP_Spatial_interp_t  spatial_interp_algorithm,
 CWP_Field_storage_t   storage,
 double               *buffer_in,
 double               *buffer_out
)
{
  PDM_UNUSED(spatial_interp_algorithm);
  PDM_UNUSED(storage);
  // Get the location of the fields degrees of freedom :
  CWP_Dof_location_t location = CWP_Field_target_dof_location_get(local_code_name,
                                                                  cpl_id,
                                                                  field_id);

  if (location == CWP_DOF_LOCATION_USER) {

    // Get the number of components of the field :
    int n_components = CWP_Field_n_component_get(local_code_name,
                                                 cpl_id,
                                                 field_id);

    // TO DO: is this what Florent Duchaine wanted to do ?
    int           n_elt_src       = 0;
    int          *src_to_tgt_idx  = NULL;
    CWP_Interp_src_data_get(local_code_name,
                            cpl_id,
                            field_id,
                            i_part,
                            &n_elt_src,
                            &src_to_tgt_idx);

    int ival = 0;
    for (int i = 0; i < n_elt_src; i++) {
      for (int j = src_to_tgt_idx[i]; j < src_to_tgt_idx[i+1]; j++) {
        for (int k1 = 0; k1 < n_components; k1++) {
          buffer_out[ival++] = buffer_in[i*n_components + k1];
        }
      }
  }

  } else {
    PDM_error(__FILE__, __LINE__, 0, "Error user interpolation not implemented for dof location %d\n", location);
  }
}

/*----------------------------------------------------------------------
 *
 * Main : advanced test : Deformable (code 2)
 *
 *---------------------------------------------------------------------*/

int
main(int argc, char *argv[]) {

  // Initialize MPI :
  MPI_Init(&argc, &argv);
  int i_rank;
  int n_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &i_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &n_rank);

  // Check running on correct number of MPI ranks :
  int n_partition = 0;
  while(2 * pow(n_partition, 2) < n_rank) n_partition++;

  const int two = 2;
  int n2 = two * (int) pow(n_partition, two);

  if (n2 != n_rank) {
    if (i_rank == 0)
      printf("      Not executed : only available if the number of processus in the form of '2 * n^2' \n");
    exit(1);
    return EXIT_SUCCESS;
  }

  // Initialize CWIPI :
  int n_code = 1;

  const char  **code_name      = malloc(sizeof(char *) * n_code);
  CWP_Status_t *is_active_rank = malloc(sizeof(CWP_Status_t) * n_code);
  double       *time_init      = malloc(sizeof(double) * n_code);
  MPI_Comm     *intra_comm     = malloc(sizeof(MPI_Comm) * n_code);

  code_name[0]      = "code2";
  is_active_rank[0] = CWP_STATUS_ON;
  time_init[0]      = 0.;

  CWP_Init(MPI_COMM_WORLD,
           n_code,
           (const char **) code_name,
           is_active_rank,
           time_init,
           intra_comm);

  printf("C - CWP_Init : OK\n");
  fflush(stdout);

  // Create the coupling :
  // CWP_DYNAMIC_MESH_DEFORMABLE allows us to take into account the modifications
  // to the mesh over the coupling steps.
  int n_part = 1;
  const char  *coupling_name     = "code1_code2";
  const char **coupled_code_name = malloc(sizeof(char *) * n_code);
  coupled_code_name[0] = "code1";
  CWP_Cpl_create(code_name[0],
                 coupling_name,
                 coupled_code_name[0],
                 CWP_INTERFACE_SURFACE,
                 CWP_COMM_PAR_WITH_PART,
                 CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE,
                 n_part,
                 CWP_DYNAMIC_MESH_DEFORMABLE,
                 CWP_TIME_EXCH_USER_CONTROLLED);

  printf("C - CWP_Cpl_create : OK\n");
  fflush(stdout);

  // Set coupling visualisation:
  CWP_Visu_set(code_name[0],
               coupling_name,
               1,
               CWP_VISU_FORMAT_ENSIGHT,
               "text");

  printf("C - CWP_Visu_set : OK\n");
  fflush(stdout);

  // Create mesh :
  int     n_vtx = 0;
  int     n_elt = 0;
  double *coords     = NULL;
  int    *elt_vtx_idx = NULL;
  int    *elt_vtx     = NULL;
  PDM_generate_mesh_rectangle_simplified(PDM_MPI_mpi_2_pdm_mpi_comm((void *) &intra_comm[0]),
                                         &n_vtx,
                                         &n_elt,
                                         &coords,
                                         &elt_vtx_idx,
                                         &elt_vtx);

  printf("C - grid_mesh : OK\n");
  fflush(stdout);

  // Interations :
  const char *send_field_name = "chinchilla";
  const char *recv_field_name = "girafe";
  int         n_components    = 1;
  double     *send_field_data = malloc(sizeof(double) * n_vtx);
  double     *recv_field_data = malloc(sizeof(double) * n_vtx);

  const int    itdeb =  1;
  const int    itend =  50;
  const double freq  =  0.20;
  const double ampl  =  0.012;
  const double phi   =  0.1;
  double       ttime = 0.0;
  double       dt    = 0.1;

  double omega = 2.0*acos(-1.0)*freq;

  for (int it = itdeb; it <= itend; it ++) {

    ttime = (it-itdeb)*dt;

    for (int i = 0; i < n_vtx; i++) {
      coords[3 * i + 2]  = ampl * (coords[3 * i]*coords[3 * i]+coords[1 + 3 * i]*coords[1 + 3 * i])*cos(omega*ttime+phi);
      send_field_data[i] = coords[3 * i + 2];
    }


    if (it == itdeb) {

      // Set the mesh vertices coordinates :
      CWP_Mesh_interf_vtx_set(code_name[0],
                              coupling_name,
                              0,
                              n_vtx,
                              coords,
                              NULL);

      printf("C - CWP_Mesh_interf_vtx_set : OK\n");
      fflush(stdout);

      // Set the mesh polygons connectivity :
      int block_id = CWP_Mesh_interf_block_add(code_name[0],
                                               coupling_name,
                                               CWP_BLOCK_FACE_POLY);

      printf("C - CWP_Mesh_interf_block_add : OK\n");
      fflush(stdout);

      CWP_Mesh_interf_f_poly_block_set(code_name[0],
                                       coupling_name,
                                       0,
                                       block_id,
                                       n_elt,
                                       elt_vtx_idx,
                                       elt_vtx,
                                       NULL);

      printf("C - CWP_Mesh_interf_f_poly_block_set : OK\n");
      fflush(stdout);

      // Finalize mesh :
      CWP_Mesh_interf_finalize(code_name[0],
                               coupling_name);

      printf("C - CWP_Mesh_interf_finalize : OK\n");
      fflush(stdout);

      // Create field :
      CWP_Field_create(code_name[0],
                       coupling_name,
                       send_field_name,
                       CWP_DOUBLE,
                       CWP_FIELD_STORAGE_INTERLACED,
                       n_components,
                       CWP_DOF_LOCATION_NODE,
                       CWP_FIELD_EXCH_SEND,
                       CWP_STATUS_ON);

      printf("C - CWP_Field_create send : OK\n");
      fflush(stdout);

      CWP_Field_create(code_name[0],
                       coupling_name,
                       recv_field_name,
                       CWP_DOUBLE,
                       CWP_FIELD_STORAGE_INTERLACED,
                       n_components,
                       CWP_DOF_LOCATION_NODE,
                       CWP_FIELD_EXCH_RECV,
                       CWP_STATUS_ON);

      printf("C - CWP_Field_create recv : OK\n");
      fflush(stdout);

      // Set interpolation property :
      CWP_Spatial_interp_property_set(code_name[0],
                                      coupling_name,
                                      "tolerance",
                                      "double",
                                      "0.1");

      printf("C - CWP_Spatial_interp_property_set : OK\n");
      fflush(stdout);

      // Set user interpolation function :
      // CWP_Interp_function_set(code_name[0],
      //                         coupling_name,
      //                         recv_field_name,
      //                         _user_interpolation_function);

      // printf("C - CWP_Interp_function_set : OK\n");
      // fflush(stdout);

    } else {
      // Update mesh :
      // Nothing to do since CWIPI stores the pointers of the arrays passed. Thus if the data
      // in the pointer is changed, it is automatically in the CWIPI code.
    }

    // Compute interpolation weights :
    CWP_Spatial_interp_weights_compute(code_name[0],
                                       coupling_name);

    printf("C - CWP_Spatial_interp_weights_compute : OK\n");
    fflush(stdout);

    // Set and exchange the field values :
    CWP_Field_data_set(code_name[0],
                       coupling_name,
                       send_field_name,
                       0,
                       CWP_FIELD_MAP_SOURCE,
                       send_field_data);

    printf("C - CWP_Field_data_set send : OK\n");
    fflush(stdout);

    CWP_Field_issend(code_name[0],
                     coupling_name,
                     send_field_name);

    printf("C - CWP_Field_issend : OK\n");
    fflush(stdout);

    CWP_Field_data_set(code_name[0],
                       coupling_name,
                       recv_field_name,
                       0,
                       CWP_FIELD_MAP_TARGET,
                       recv_field_data);

    printf("C - CWP_Field_data_set recv : OK\n");
    fflush(stdout);

    CWP_Field_irecv(code_name[0],
                    coupling_name,
                    recv_field_name);

    printf("C - CWP_Field_irecv : OK\n");
    fflush(stdout);

    CWP_Field_wait_issend(code_name[0],
                          coupling_name,
                          send_field_name);

    printf("C - CWP_Field_wait_issend : OK\n");
    fflush(stdout);

    CWP_Field_wait_irecv(code_name[0],
                         coupling_name,
                         recv_field_name);

    printf("C - CWP_Field_wait_irecv : OK\n");
    fflush(stdout);

  // Check interpolation :
  int n_uncomputed_tgts = CWP_N_uncomputed_tgts_get(code_name[0],
                                                    coupling_name,
                                                    recv_field_name,
                                                    0);

  printf("C - CWP_N_uncomputed_tgts_get : OK\n");
  fflush(stdout);

  int *uncomputed_tgts = NULL;
  if (n_uncomputed_tgts != 0) {
    uncomputed_tgts = CWP_Uncomputed_tgts_get(code_name[0],
                                              coupling_name,
                                              recv_field_name,
                                              0);

    printf("C - CWP_Uncomputed_tgts_get : OK\n");
    fflush(stdout);
  }

  PDM_UNUSED(n_uncomputed_tgts);
  PDM_UNUSED(uncomputed_tgts);
  } // end interations

  // Delete field :
  CWP_Field_del(code_name[0],
                coupling_name,
                send_field_name);

  printf("C - CWP_Field_del send : OK\n");
  fflush(stdout);

  CWP_Field_del(code_name[0],
                coupling_name,
                recv_field_name);

  printf("C - CWP_Field_del recv : OK\n");
  fflush(stdout);

  // Delete Mesh :
  CWP_Mesh_interf_del(code_name[0],
                      coupling_name);

  printf("C - CWP_Mesh_interf_del : OK\n");
  fflush(stdout);

  // Delete the coupling :
  CWP_Cpl_del(code_name[0],
              coupling_name);

  printf("C - CWP_Cpl_del : OK\n");
  fflush(stdout);

  // free
  free(code_name);
  free(is_active_rank);
  free(time_init);
  free(coupled_code_name);
  free(coords);
  free(elt_vtx_idx);
  free(elt_vtx);
  free(send_field_data);
  free(recv_field_data);

  // Finalize CWIPI :
  CWP_Finalize();

  printf("C - CWP_Finalize : OK\n");
  fflush(stdout);

  // Finalize MPI :
  MPI_Finalize();

  return EXIT_SUCCESS;
}
