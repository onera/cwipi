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

#include "grid_mesh.h"

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
  // Get the location of the fields degrees of freedom :
  CWP_Dof_location_t location = CWP_Field_target_dof_location_get(local_code_name,
                                                                  cpl_id,
                                                                  field_id);

  if (location == CWP_DOF_LOCATION_CELL_CENTER) {

    // Get the number of components of the field :
    int n_components = CWP_Field_n_component_get(local_code_name,
                                                 cpl_id,
                                                 field_id);

    // TO DO: pourquoi on définit dans le code cell_vertex la fonction d'interpolation
    //        qui va être utilisé côté code cell_center alors?

    // for (int i = 0; i < ; i++) {
    //   ((double *) buffer_out)[i] = ((double *) buffer_in)[closest_src_pt[i]-1];
    // }

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

  // Set coupling visualisation:
  CWP_Visu_set(code_name[0],
               coupling_name,
               1,
               CWP_VISU_FORMAT_ENSIGHT,
               "text");

  // Set mesh data :
  int i_intra_rank;
  MPI_Comm_rank(intra_comm[0], &i_intra_rank);
  srand(i_intra_rank+time(0));

  const int    itdeb     =  1;
  const int    itend     =  50;
  const double freq      =  0.20;
  const double ampl      =  0.012;
  const double phi       =  0.1;
  const double xmin      = -10;
  const double xmax      =  10;
  const double ymin      = -10;
  const double ymax      =  10;
  const int    n_vtx_seg =  10;
  const double randLevel =  0.4;

  int n_vtx = n_vtx_seg * n_vtx_seg;
  int n_elt = (n_vtx_seg - 1) * (n_vtx_seg - 1);

  double *coords     = malloc(sizeof(double) * 3 * n_vtx);
  int    *connec_idx = malloc(sizeof(int)    * (n_elt + 1));
  int    *connec     = malloc(sizeof(int)    * 4 * n_elt);

  grid_mesh(xmin,
            xmax,
            ymin,
            ymax,
            randLevel,
            n_vtx_seg,
            n_partition,
            coords,
            connec_idx,
            connec,
            intra_comm[0]);

  // Interations :
  const char *field_name      = "a super fancy field";
  int         n_components    = 1;
  double     *send_field_data = malloc(sizeof(double) * n_vtx);
  double     *recv_field_data = malloc(sizeof(double) * n_vtx);

  double ttime = 0.0;
  double dt = 0.1;

  double omega = 2.0*acos(-1.0)*freq;

  for (int it = itdeb; it <= itend; it ++) {

    ttime = (it-itdeb)*dt;

    for (int i = 0; i < n_vtx; i++) {
      coords[2 + 3 * i] = ampl * (coords[3 * i]*coords[3 * i]+coords[1 + 3 * i]*coords[1 + 3 * i])*cos(omega*ttime+phi);
      send_field_data[i] = coords[2 + 3 * i];
    }


    if (it == itdeb) {

      // Set the mesh vertices coordinates :
      CWP_Mesh_interf_vtx_set(code_name[0],
                              coupling_name,
                              0,
                              n_vtx,
                              coords,
                              NULL);

      // Set the mesh polygons connectiviy :
      int block_id = CWP_Mesh_interf_block_add(code_name[0],
                                               coupling_name,
                                               CWP_BLOCK_FACE_POLY);

      CWP_Mesh_interf_f_poly_block_set(code_name[0],
                                       coupling_name,
                                       0,
                                       block_id,
                                       n_elt,
                                       connec_idx,
                                       connec,
                                       NULL);

      // Finalize mesh :
      CWP_Mesh_interf_finalize(code_name[0],
                               coupling_name);

      // Create field :
      CWP_Field_create(code_name[0],
                       coupling_name,
                       field_name,
                       CWP_DOUBLE,
                       CWP_FIELD_STORAGE_INTERLACED,
                       n_components,
                       CWP_DOF_LOCATION_NODE,
                       CWP_FIELD_EXCH_SEND,
                       CWP_STATUS_ON);

      CWP_Field_create(code_name[0],
                       coupling_name,
                       field_name,
                       CWP_DOUBLE,
                       CWP_FIELD_STORAGE_INTERLACED,
                       n_components,
                       CWP_DOF_LOCATION_NODE,
                       CWP_FIELD_EXCH_RECV,
                       CWP_STATUS_ON);

      // Set interpolation property :
      CWP_Spatial_interp_property_set(code_name[0],
                                  coupling_name,
                                  "tolerance",
                                  "double",
                                  "0.1");

    } else {
      // Update mesh :
      // Nothing to do since CWIPI stores the pointers of the arrays passed. Thus if the data
      // in the pointer is changed, it is automatically in the CWIPI code.
    }

    // Compute interpolation weights :
    CWP_Spatial_interp_weights_compute(code_name[0],
                                       coupling_name);

    // Set and exchange the field values :
    CWP_Field_data_set(code_name[0],
                       coupling_name,
                       field_name,
                       0,
                       CWP_FIELD_MAP_SOURCE,
                       send_field_data);

    CWP_Field_issend(code_name[0],
                     coupling_name,
                     field_name);

    CWP_Field_data_set(code_name[0],
                       coupling_name,
                       field_name,
                       0,
                       CWP_FIELD_MAP_TARGET,
                       recv_field_data);

    CWP_Field_irecv(code_name[0],
                    coupling_name,
                    field_name);

    CWP_Field_wait_issend(code_name[0],
                          coupling_name,
                          field_name);

    CWP_Field_wait_irecv(code_name[0],
                         coupling_name,
                         field_name);

  // Check interpolation :
  int n_uncomputed_tgts = CWP_N_uncomputed_tgts_get(code_name[0],
                                                    coupling_name,
                                                    field_name,
                                                    0);
  int *uncomputed_tgts = CWP_Uncomputed_tgts_get(code_name[0],
                                                 coupling_name,
                                                 field_name,
                                                 0);

  } // end interations

  // Delete field :
  CWP_Field_del(code_name[0],
                coupling_name,
                field_name);

  // Delete Mesh :
  CWP_Mesh_interf_del(code_name[0],
                      coupling_name);

  // Delete the coupling :
  CWP_Cpl_del(code_name[0],
              coupling_name);

  // free
  free(code_name);
  free(is_active_rank);
  free(time_init);
  free(coupled_code_name);
  free(coords);
  free(connec_idx);
  free(connec);
  free(send_field_data);
  free(recv_field_data);

  // Finalize CWIPI :
  CWP_Finalize();

  // Finalize MPI :
  MPI_Finalize();

  return EXIT_SUCCESS;
}
