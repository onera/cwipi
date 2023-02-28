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
#include <string.h>
#include <time.h>

#include "cwipi.h"
#include "cwp.h"

#include "pdm.h"
#include "pdm_multipart.h"
#include "pdm_logging.h"
#include "pdm_error.h"

#include "pdm_dmesh_nodal_priv.h"
#include "pdm_block_to_part.h"
#include "pdm_mesh_nodal.h"
#include "pdm_part_connectivity_transform.h"
#include "pdm_distrib.h"


/*----------------------------------------------------------------------
 *
 * Display usage
 *
 * parameters:
 *   exit code           <-- Exit code
 *---------------------------------------------------------------------*/

static void
_usage(int exit_code) {
  printf("\n"
         "  Usage: \n\n"
         "  -f1             first filename.\n\n"
         "  -f2             second filename.\n\n"
         "  -v              verbose.\n\n"
         "  -h              this message.\n\n");

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
_read_args
(
  int            argc,
  char         **argv,
  int           *verbose,
  char         **filename1,
  char         **filename2
)
{
  int i = 1;

  // Parse and check command line
  while (i < argc) {
    if (strcmp(argv[i], "-h") == 0) {
      _usage(EXIT_SUCCESS);
    }
    else if (strcmp(argv[i], "-v") == 0) {
      *verbose = 1;
    }
    else if (strcmp(argv[i], "-f1") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *filename1 = argv[i];
      }
    }
    else if (strcmp(argv[i], "-f2") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *filename2 = argv[i];
      }
    }
    else
      _usage(EXIT_FAILURE);
    i++;
  }
}

/*----------------------------------------------------------------------
 *
 * Read VTK file
 *
 *---------------------------------------------------------------------*/

static void
_read_vtk
(
 const char        *filename,
       int          *n_vtx,
       int          *n_face,
       PDM_g_num_t **face_vtx,
       double      **vtx_coord,
       double      **vtx_field
)
{
  char line[999];


  FILE *f = fopen(filename, "r");

  if (f == NULL) {
    PDM_error(__FILE__, __LINE__, 0, "Could not read file %s\n", filename);
  }

  // PASS 1
  while (1) {

    int stat = fscanf(f, "%s", line);

    if (stat == EOF) {
      // End of file
      break;
    }

    if (strstr(line, "POINTS") != NULL) {
      // Get dimension
      fscanf(f, "%d", &(*n_vtx));
      char *data_type = malloc(sizeof(char) * 10);
      fscanf(f, "%s", data_type);
      free(data_type);

      // Malloc
      (*vtx_coord) = malloc(sizeof(double) * 3 * ((*n_vtx)));

      // Get coordinates
      for (int i = 0; i < (*n_vtx); i++) {
        for (int j = 0; j < 3; j++) {
          fscanf(f, "%le", &(*vtx_coord)[3*i + j]);
        }
      }
    } // end if POINTS

    if (strstr(line, "CELLS") != NULL) {
      // Get dimension
      int n_cell = 0;
      fscanf(f, "%d", &n_cell);
      int size = 0;
      fscanf(f, "%d", &size);

      // Malloc
      (*face_vtx) = malloc(sizeof(PDM_g_num_t) * size);

      // Get connectivity
      int n_vtx_cell = 0;
      for (int i = 0; i < n_cell; i++) {
        fscanf(f, "%d", &n_vtx_cell);
        if (n_vtx_cell == 3) {
          for (int j = 0; j < 3; j++) {
            int id = -1;
           // fscanf(f, PDM_FMT_G_NUM, &(*face_vtx)[3*(*n_face)+j]);
            fscanf(f, "%d", &id);
            (*face_vtx)[3*(*n_face)+j] = (PDM_g_num_t) id + 1;
          }
          (*n_face)++;
        } else {
          PDM_g_num_t other_cell_vtx = 0;
          for (int j = 0; j < n_vtx_cell; j++) {
            fscanf(f, PDM_FMT_G_NUM, &other_cell_vtx);
          }
        }
      }

      (*face_vtx) = realloc((*face_vtx), 3 * (*n_face) * sizeof(PDM_g_num_t));
    } // end if CELLS

    if (strstr(line, "LOOKUP_TABLE") != NULL) {
      char *table = malloc(sizeof(char) * 10);
      fscanf(f, "%s", table);
      free(table);

      // Malloc
      (*vtx_field) = malloc(sizeof(double) * (*n_vtx));

      // Get connectivity
      for (int i = 0; i < (*n_vtx); i++) {
        fscanf(f, "%lf", &(*vtx_field)[i]);
      }
    } // end if POINT_DATA

  } // while not EOF
}


static PDM_dmesh_nodal_t *
_set_dmesh_nodal
(
 const PDM_MPI_Comm  comm,
       double       *dvtx_coord,
       PDM_g_num_t  *dface_vtx,
       PDM_g_num_t  *distrib_vtx,
       PDM_g_num_t  *distrib_face
 )
{
  int i_rank, n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);



  /*
   *  Create dmesh nodal
   */
  PDM_g_num_t gn_vtx  = distrib_vtx [n_rank];
  PDM_g_num_t gn_face = distrib_face[n_rank];
  int dn_vtx  = distrib_vtx [i_rank+1] - distrib_vtx [i_rank];
  int dn_face = distrib_face[i_rank+1] - distrib_face[i_rank];

  PDM_dmesh_nodal_t *dmn = PDM_DMesh_nodal_create(comm,
                                                  2,
                                                  gn_vtx,
                                                  0,
                                                  gn_face,
                                                  0);

  PDM_DMesh_nodal_coord_set(dmn,
                            dn_vtx,
                            dvtx_coord,
                            PDM_OWNERSHIP_KEEP);

  dmn->surfacic->n_g_elmts = gn_face;
  int id_section = PDM_DMesh_nodal_elmts_section_add(dmn->surfacic,
                                                     PDM_MESH_NODAL_TRIA3);
  PDM_DMesh_nodal_elmts_section_std_set(dmn->surfacic,
                                        id_section,
                                        dn_face,
                                        dface_vtx,
                                        PDM_OWNERSHIP_KEEP);

  int n_group = 1;
  int *dgroup_elt_idx = (int *) malloc(sizeof(int) * (n_group + 1));
  dgroup_elt_idx[0] = 0;
  dgroup_elt_idx[1] = dn_face;

  PDM_g_num_t *dgroup_elt = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * dgroup_elt_idx[n_group]);
  for (int i = 0; i < dn_face; i++) {
    dgroup_elt[i] = distrib_face[i_rank] + i + 1;
  }
  PDM_DMesh_nodal_elmts_group_set(dmn->surfacic,
                                  n_group,
                                  dgroup_elt_idx,
                                  dgroup_elt,
                                  PDM_OWNERSHIP_KEEP);

  PDM_dmesh_nodal_generate_distribution(dmn);

  return dmn;
}


static void
_gen_part_data
(
 const char               *filename,
 const PDM_MPI_Comm        comm,
 const int                 n_part,
 const PDM_split_dual_t    part_method,
       int               **pn_face,
       int               **pn_vtx,
       int              ***pface_vtx_idx,
       int              ***pface_vtx,
       double           ***pvtx_coord,
       PDM_g_num_t      ***pface_ln_to_gn,
       PDM_g_num_t      ***pvtx_ln_to_gn,
       double           ***pvtx_field

 )
{
  int i_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);

  int            dn_vtx     = 0;
  int            dn_face    = 0;
  PDM_g_num_t  *dface_vtx  = NULL;
  double       *dvtx_coord = NULL;
  double       *dvtx_field = NULL;

  if (i_rank == 0) {
    _read_vtk(filename,
              &dn_vtx,
              &dn_face,
              &dface_vtx,
              &dvtx_coord,
              &dvtx_field);
  }

  PDM_log_trace_array_long(dface_vtx, 3*dn_face, "dface_vtx : ");

  PDM_g_num_t *distrib_vtx  = PDM_compute_entity_distribution(comm, dn_vtx);
  PDM_g_num_t *distrib_face = PDM_compute_entity_distribution(comm, dn_face);

  // DMesh nodal
  PDM_dmesh_nodal_t *dmn = _set_dmesh_nodal(comm,
                                            dvtx_coord,
                                            dface_vtx,
                                            distrib_vtx,
                                            distrib_face);

  int n_zone = 1;
  int n_part_zones = n_part;
  PDM_multipart_t *mpart = PDM_multipart_create(n_zone,
                                                &n_part_zones,
                                                PDM_FALSE,
                                                part_method,
                                                PDM_PART_SIZE_HOMOGENEOUS,
                                                NULL,
                                                comm,
                                                PDM_OWNERSHIP_KEEP);

  PDM_multipart_set_reordering_options(mpart,
                                       -1,
                                       "PDM_PART_RENUM_CELL_NONE",
                                       NULL,
                                       "PDM_PART_RENUM_FACE_NONE");

  PDM_multipart_register_dmesh_nodal(mpart, 0, dmn);

  PDM_multipart_run_ppart(mpart);


  *pn_vtx         = malloc(sizeof(int          ) * n_part);
  *pvtx_coord     = malloc(sizeof(double      *) * n_part);
  *pvtx_ln_to_gn  = malloc(sizeof(PDM_g_num_t *) * n_part);
  *pn_face        = malloc(sizeof(int          ) * n_part);
  *pface_vtx_idx  = malloc(sizeof(int         *) * n_part);
  *pface_vtx      = malloc(sizeof(int         *) * n_part);
  *pface_ln_to_gn = malloc(sizeof(PDM_g_num_t *) * n_part);

  for (int ipart = 0; ipart < n_part; ipart++) {

    /* Vertices */
    PDM_g_num_t *_vtx_ln_to_gn;
    (*pn_vtx)[ipart] = PDM_multipart_part_ln_to_gn_get(mpart,
                                                       0,
                                                       ipart,
                                                       PDM_MESH_ENTITY_VERTEX,
                                                       &_vtx_ln_to_gn,
                                                       PDM_OWNERSHIP_USER);
    (*pvtx_ln_to_gn)[ipart] = malloc(sizeof(PDM_g_num_t) * (*pn_vtx)[ipart]);
    memcpy((*pvtx_ln_to_gn)[ipart], _vtx_ln_to_gn, sizeof(PDM_g_num_t) * (*pn_vtx)[ipart]);

    double *_vtx_coord;
    PDM_multipart_part_vtx_coord_get(mpart,
                                     0,
                                     ipart,
                                     &_vtx_coord,
                                     PDM_OWNERSHIP_USER);
    (*pvtx_coord)[ipart] = malloc(sizeof(double) * (*pn_vtx)[ipart] * 3);
    memcpy((*pvtx_coord)[ipart], _vtx_coord, sizeof(double) * (*pn_vtx)[ipart] * 3);


    /* Faces */
    PDM_g_num_t *_face_ln_to_gn;
    (*pn_face)[ipart] = PDM_multipart_part_ln_to_gn_get(mpart,
                                                       0,
                                                       ipart,
                                                       PDM_MESH_ENTITY_FACE,
                                                       &_face_ln_to_gn,
                                                       PDM_OWNERSHIP_USER);
    (*pface_ln_to_gn)[ipart] = malloc(sizeof(PDM_g_num_t) * (*pn_face)[ipart]);
    memcpy((*pface_ln_to_gn)[ipart], _face_ln_to_gn, sizeof(PDM_g_num_t) * (*pn_face)[ipart]);

    int *_face_vtx;
    int *_face_vtx_idx;
    PDM_multipart_part_connectivity_get(mpart,
                                        0,
                                        ipart,
                                        PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                        &_face_vtx,
                                        &_face_vtx_idx,
                                        PDM_OWNERSHIP_USER);

    if (_face_vtx != NULL) {
      (*pface_vtx_idx)[ipart] = malloc(sizeof(int) * ((*pn_face)[ipart]+1));
      memcpy((*pface_vtx_idx)[ipart], _face_vtx_idx, sizeof(int) * ((*pn_face)[ipart]+1));

      (*pface_vtx)[ipart] = malloc(sizeof(int) * _face_vtx_idx[(*pn_face)[ipart]]);
      memcpy((*pface_vtx)[ipart], _face_vtx, sizeof(int) * _face_vtx_idx[(*pn_face)[ipart]]);
    }

    else {
      int *_face_edge;
      int *_face_edge_idx;
      PDM_multipart_part_connectivity_get(mpart,
                                          0,
                                          ipart,
                                          PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                          &_face_edge,
                                          &_face_edge_idx,
                                          PDM_OWNERSHIP_KEEP);

      int *_edge_vtx;
      int *_edge_vtx_idx;
      PDM_multipart_part_connectivity_get(mpart,
                                          0,
                                          ipart,
                                          PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                          &_edge_vtx,
                                          &_edge_vtx_idx,
                                          PDM_OWNERSHIP_KEEP);

      (*pface_vtx_idx)[ipart] = malloc(sizeof(int) * ((*pn_face)[ipart]+1));
      memcpy((*pface_vtx_idx)[ipart], _face_edge_idx, sizeof(int) * ((*pn_face)[ipart]+1));


      PDM_compute_face_vtx_from_face_and_edge((*pn_face)[ipart],
                                              _face_edge_idx,
                                              _face_edge,
                                              _edge_vtx,
                                              &(*pface_vtx)[ipart]);
    }

  }

  PDM_DMesh_nodal_free(dmn);
  PDM_multipart_free(mpart);



  /* Exchange field to partitions */
  PDM_block_to_part_t *btp = PDM_block_to_part_create(distrib_vtx,
                               (const PDM_g_num_t **) *pvtx_ln_to_gn,
                                                      *pn_vtx,
                                                      n_part,
                                                      comm);

  int one = 1;

  PDM_block_to_part_exch(btp,
                         sizeof(double),
                         PDM_STRIDE_CST_INTERLACED,
                         &one,
              (void   *) dvtx_field,
                         NULL,
              (void ***) pvtx_field);

  PDM_block_to_part_free(btp);

  free(distrib_vtx);
  free(distrib_face);
  if (dvtx_field != NULL) {
    free(dvtx_field);
  }
}


/*----------------------------------------------------------------------
 *
 * Main : preCICE wind turbine blade test case
 *
 *---------------------------------------------------------------------*/

int
main(int argc, char *argv[]) {

  int verbose = 0;
  char *filename1 = NULL;
  char *filename2 = NULL;

  PDM_split_dual_t part_method = PDM_SPLIT_DUAL_WITH_HILBERT;
  int disjoint_comm = 1;

  _read_args (argc,
              argv,
              &verbose,
              &filename1,
              &filename2);

  // Initialize MPI
  MPI_Init(&argc, &argv);

  MPI_Comm comm = MPI_COMM_WORLD;

  int i_rank;
  int n_rank;
  MPI_Comm_rank(comm, &i_rank);
  MPI_Comm_size(comm, &n_rank);



  const char *all_code_names[2] = {"code1", "code2"};
  const char *all_file_names[2] = {filename1, filename2};
  int all_n_part[2] = {1, 1};
  int has_code[2] = {0, 0};
  if (disjoint_comm) {
    has_code[0] = i_rank < n_rank/2;
    has_code[1] = !has_code[0];
  }
  else {
    has_code[0] = 1;
    has_code[1] = 1;
  }

  int n_code = has_code[0] + has_code[1];

  int           *code_id           = malloc(sizeof(int         ) * n_code);
  int           *n_part            = malloc(sizeof(int         ) * n_code);
  const char   **code_name         = malloc(sizeof(char       *) * n_code);
  const char   **coupled_code_name = malloc(sizeof(char       *) * n_code);
  CWP_Status_t  *is_active_rank    = malloc(sizeof(CWP_Status_t) * n_code);
  double        *time_init         = malloc(sizeof(double      ) * n_code);
  MPI_Comm      *intra_comm        = malloc(sizeof(MPI_Comm    ) * n_code);
  const char   **file_name         = malloc(sizeof(char       *) * n_code);

  n_code = 0;
  for (int icode = 0; icode < 2; icode++) {
    if (has_code[icode]) {
      code_id          [n_code] = icode+1;
      code_name        [n_code] = all_code_names[icode];
      coupled_code_name[n_code] = all_code_names[(icode+1)%2];
      n_part           [n_code] = all_n_part    [icode];
      is_active_rank   [n_code] = CWP_STATUS_ON;
      time_init        [n_code] = 0.;
      file_name        [n_code] = all_file_names[icode];

      if (verbose) {
        log_trace("Running %s, coupled with %s, n_part = %d\n",
                  code_name[n_code], coupled_code_name[n_code], n_part[n_code]);
      }
      n_code++;
    }
  }


  // Set up
  CWP_Init(comm,
           n_code,
           (const char **) code_name,
           is_active_rank,
           time_init,
           intra_comm);

  MPI_Barrier(comm);
  if (i_rank == 0) {
    printf("CWIPI Init OK\n");
    fflush(stdout);
  }

  /* Create coupling */
  const char *cpl_name = "c_new_api_wind_turbine_blade";


  CWP_Spatial_interp_t spatial_interp = CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE;
  // CWP_Spatial_interp_t spatial_interp = CWP_SPATIAL_INTERP_FROM_CLOSEST_POINT_LEAST_SQUARES;
  // CWP_Spatial_interp_t spatial_interp = CWP_SPATIAL_INTERP_FROM_INTERSECTION;
  // CWP_Spatial_interp_t spatial_interp = CWP_SPATIAL_INTERP_FROM_IDENTITY;

  for (int icode = 0; icode < n_code; icode++) {
    CWP_Cpl_create(code_name[icode],
                   cpl_name,
                   coupled_code_name[icode],
                   CWP_INTERFACE_SURFACE,
                   CWP_COMM_PAR_WITH_PART,
                   spatial_interp,
                   n_part[icode],
                   CWP_DYNAMIC_MESH_STATIC,
                   CWP_TIME_EXCH_USER_CONTROLLED);

  }

  for (int icode = 0; icode < n_code; icode++) {
    CWP_Visu_set(code_name[icode],        // Code name
                 cpl_name,                // Coupling id
                 1,                       // Postprocessing frequency
                 CWP_VISU_FORMAT_ENSIGHT, // Postprocessing format
                 "text");                 // Postprocessing option
  }

  MPI_Barrier(comm);
  if (i_rank == 0) {
    printf("Create coupling OK\n");
    fflush(stdout);
  }

  // Mesh and field
  int          **pn_face        = malloc(sizeof(int          *) * n_code);
  int          **pn_vtx         = malloc(sizeof(int          *) * n_code);
  int         ***pface_vtx_idx  = malloc(sizeof(int         **) * n_code);
  int         ***pface_vtx      = malloc(sizeof(int         **) * n_code);
  double      ***pvtx_coord     = malloc(sizeof(double      **) * n_code);
  PDM_g_num_t ***pface_ln_to_gn = malloc(sizeof(PDM_g_num_t **) * n_code);
  PDM_g_num_t ***pvtx_ln_to_gn  = malloc(sizeof(PDM_g_num_t **) * n_code);
  double      ***pvtx_field     = malloc(sizeof(double      **) * n_code);

  for (int icode = 0; icode < n_code; icode++) {
    PDM_MPI_Comm mesh_comm = PDM_MPI_mpi_2_pdm_mpi_comm((void *) &intra_comm[icode]);

    _gen_part_data(file_name[icode],
                   mesh_comm,
                   n_part[icode],
                   part_method,
                   &pn_face       [icode],
                   &pn_vtx        [icode],
                   &pface_vtx_idx [icode],
                   &pface_vtx     [icode],
                   &pvtx_coord    [icode],
                   &pface_ln_to_gn[icode],
                   &pvtx_ln_to_gn [icode],
                   &pvtx_field    [icode]);

    int block_id = CWP_Mesh_interf_block_add(code_name[icode],
                                             cpl_name,
                                             CWP_BLOCK_FACE_POLY);

    for (int ipart = 0; ipart < n_part[icode]; ipart++) {
      CWP_Mesh_interf_vtx_set(code_name[icode],
                              cpl_name,
                              ipart,
                              pn_vtx       [icode][ipart],
                              pvtx_coord   [icode][ipart],
                              pvtx_ln_to_gn[icode][ipart]);

      CWP_Mesh_interf_f_poly_block_set(code_name[icode],
                                       cpl_name,
                                       ipart,
                                       block_id,
                                       pn_face       [icode][ipart],
                                       pface_vtx_idx [icode][ipart],
                                       pface_vtx     [icode][ipart],
                                       pface_ln_to_gn[icode][ipart]);
    }

    CWP_Mesh_interf_finalize(code_name[icode], cpl_name);
  }

   MPI_Barrier(comm);
  if (i_rank == 0) {
    printf("Set mesh OK\n");
    fflush(stdout);
  }


  // Field
  CWP_Status_t visu_status = CWP_STATUS_ON;
  const char *field_name = "field";

  int stride = 1;
  double ***recv_val  = malloc(sizeof(double **) * n_code);

  for (int icode = 0; icode < n_code; icode++) {

    recv_val[icode] = malloc(sizeof(double *) * n_part[icode]);
    for (int ipart = 0; ipart < n_part[icode]; ipart++) {
      recv_val[icode][ipart] = malloc(sizeof(double) * pn_vtx[icode][ipart]);
    }

    CWP_Field_create(code_name[icode],
                     cpl_name,
                     field_name,
                     CWP_DOUBLE,
                     CWP_FIELD_STORAGE_INTERLACED,
                     stride,
                     CWP_DOF_LOCATION_NODE,
                     CWP_FIELD_EXCH_SENDRECV,
                     visu_status);

    for (int ipart = 0; ipart < n_part[icode]; ipart++) {
      CWP_Field_data_set(code_name[icode],
                         cpl_name,
                         field_name,
                         ipart,
                         CWP_FIELD_MAP_SOURCE,
                         pvtx_field[icode][ipart]);

      CWP_Field_data_set(code_name[icode],
                         cpl_name,
                         field_name,
                         ipart,
                         CWP_FIELD_MAP_TARGET,
                         recv_val[icode][ipart]);
    }
  }

  // Exchange
  for (int icode = 0; icode < n_code; icode++) {
    CWP_Spatial_interp_property_set(code_name[icode],
                                    cpl_name,
                                    "tolerance",
                                    "double",
                                    "1e-1");
    CWP_Spatial_interp_weights_compute(code_name[icode], cpl_name);
  }

  for (int icode = 0; icode < n_code; icode++) {
      CWP_Field_issend(code_name[icode], cpl_name, field_name);
      CWP_Field_irecv (code_name[icode], cpl_name, field_name);
  }

  for (int icode = 0; icode < n_code; icode++) {
      CWP_Field_wait_issend(code_name[icode], cpl_name, field_name);
      CWP_Field_wait_irecv (code_name[icode], cpl_name, field_name);
  }

  MPI_Barrier(comm);
  if (i_rank == 0) {
    printf("Exchange fields OK\n");
    fflush(stdout);
  }


  // TO DO : check interpolation error
  //...


  /* Free memory */
  for (int icode = 0; icode < n_code; icode++) {
    CWP_Mesh_interf_del(code_name[icode], cpl_name);

    CWP_Cpl_del(code_name[icode], cpl_name);
  }

  for (int icode = 0; icode < n_code; icode++) {
    for (int ipart = 0; ipart < n_part[icode]; ipart++) {
      free(pface_vtx_idx [icode][ipart]);
      free(pface_vtx     [icode][ipart]);
      free(pvtx_coord    [icode][ipart]);
      free(pface_ln_to_gn[icode][ipart]);
      free(pvtx_ln_to_gn [icode][ipart]);
      free(pvtx_field    [icode][ipart]);
      free(recv_val      [icode][ipart]);
    }
    free(pn_face       [icode]);
    free(pn_vtx        [icode]);
    free(pface_vtx_idx [icode]);
    free(pface_vtx     [icode]);
    free(pvtx_coord    [icode]);
    free(pface_ln_to_gn[icode]);
    free(pvtx_ln_to_gn [icode]);
    free(pvtx_field    [icode]);
    free(recv_val      [icode]);
  }
  free(pn_face       );
  free(pn_vtx        );
  free(pface_vtx_idx );
  free(pface_vtx     );
  free(pvtx_coord    );
  free(pface_ln_to_gn);
  free(pvtx_ln_to_gn );
  free(pvtx_field    );
  free(recv_val      );

  free(code_id);
  free(n_part);
  free(coupled_code_name);
  free(code_name);
  free(is_active_rank);
  free(intra_comm);
  free(time_init);

  // Finalize CWIPI
  CWP_Finalize();

  // Finalize MPI
  MPI_Finalize();

  return EXIT_SUCCESS;
}
