#include "doctest/doctest.h"
#include "doctest/extensions/doctest_mpi.h"
#include "cwp.h"

#include <string.h>
#include <vector>

static double _max_error = 1e-6;
static std::vector<const char *> _forward  = {"code1", "code2"};
static std::vector<const char *> _backward = {"code2", "code1"};

#define CHECK_EPS(array, array_expected, size) \
  for(size_t i = 0; i < size; ++i){ CHECK(std::abs(array[i]-array_expected[i]) < _max_error); }

static void
_intersection_user_interpolation
(
  const char   *local_code_name, // Local code name
  const char   *cpl_id,          // Coupling identifier
  const char   *field_id,        // Field identifier
        int     i_part,          // Partition
        double *buffer_in,       // Input buffer
        double *buffer_out       // Output buffer
)
{
  // GENERAL GEOMETRIC DATA

  // Get the number of components of the field
  int n_components = CWP_Field_n_components_get(local_code_name,
                                                cpl_id,
                                                field_id);

  // Get the field data storage type
  CWP_Field_storage_t storage = CWP_Field_storage_get(local_code_name,
                                                      cpl_id,
                                                      field_id);

  // Get data from the target of the exchange
  int n_elt_tgt      = 0;
  int n_computed_tgt = 0;

  int *computed_tgt   = nullptr;
  int *tgt_to_src_idx = nullptr;

  CWP_Field_tgt_data_properties_get(local_code_name,
                                    cpl_id,
                                    field_id,
                                    i_part,
                                    &n_elt_tgt,
                                    &n_computed_tgt,
                                    &computed_tgt,
                                    &tgt_to_src_idx);

  // INTERSECTION SPATIAL INTERPOLATION ALGORITHM

  // Get the intersection volumes
  double *volumes = nullptr;

  CWP_Field_intersection_volumes_get(local_code_name,
                                     cpl_id,
                                     field_id,
                                     i_part,
                                     &volumes);

  // Fill buffer out
  for (int i = 0; i < n_components * n_elt_tgt; i++) {
    buffer_out[i] = 0.;
  }

  for (int i_ref = 0; i_ref < n_computed_tgt; i_ref++) {
    int i = computed_tgt[i_ref] - 1;
    double sum_w = 0;

    for (int k = tgt_to_src_idx[i_ref]; k < tgt_to_src_idx[i_ref + 1]; k++) {
      double w = volumes[k];
      sum_w += w;

      if (storage == CWP_FIELD_STORAGE_INTERLEAVED) {
        for (int j = 0; j < n_components; j++) {
          buffer_out[n_elt_tgt * j + i] += w * buffer_in[tgt_to_src_idx[n_computed_tgt] * j + k];
        }
      }
      else {
        for (int j = 0; j < n_components; j++) {
          buffer_out[n_components * i + j] += w * buffer_in[n_components * k + j];
        }
      }
    }

    sum_w = 1. / sum_w;
    if (storage == CWP_FIELD_STORAGE_INTERLEAVED) {
      for (int j = 0; j < n_components; j++) {
        buffer_out[n_elt_tgt * j + i] *= sum_w;
      }
    }
    else {
      for (int j = 0; j < n_components; j++) {
        buffer_out[n_components * i + j] *= sum_w;
      }
    }
  }
}

static void
_nearest_user_interpolation
(
  const char   *local_code_name, // Local code name
  const char   *cpl_id,          // Coupling identifier
  const char   *field_id,        // Field identifier
        int     i_part,          // Partition
        double *buffer_in,       // Input buffer
        double *buffer_out       // Output buffer
)
{
  // GENERAL GEOMETRIC DATA

  // Get the number of components of the field
  int n_components = CWP_Field_n_components_get(local_code_name,
                                                cpl_id,
                                                field_id);

  // Get the field data storage type
  CWP_Field_storage_t storage = CWP_Field_storage_get(local_code_name,
                                                      cpl_id,
                                                      field_id);

  // Get data from the target of the exchange
  int n_elt_tgt      = 0;
  int n_computed_tgt = 0;

  int *computed_tgt   = nullptr;
  int *tgt_to_src_idx = nullptr;

  CWP_Field_tgt_data_properties_get(local_code_name,
                                    cpl_id,
                                    field_id,
                                    i_part,
                                    &n_elt_tgt,
                                    &n_computed_tgt,
                                    &computed_tgt,
                                    &tgt_to_src_idx);

  // NEAREST NEIGHBORS SPATIAL INTERPOLATION ALGORITHM

  // Get squared distances to nearest source points
  double *distances2 = nullptr;

  CWP_Field_nearest_neighbors_distances_get(local_code_name,
                                            cpl_id,
                                            field_id,
                                            i_part,
                                            &distances2);

  // Fill buffer out
  for (int i = 0; i < n_components * n_elt_tgt; i++) {
    buffer_out[i] = 0.;
  }

  for (int i_ref = 0; i_ref < n_computed_tgt; i_ref++) {
    int i = computed_tgt[i_ref] - 1;
    double sum_w = 0;

    for (int k = tgt_to_src_idx[i_ref]; k < tgt_to_src_idx[i_ref + 1]; k++) {
      double w = 1. / std::max(1e-30, distances2[k]);
      sum_w += w;

      if (storage == CWP_FIELD_STORAGE_INTERLEAVED) {
        for (int j = 0; j < n_components; j++) {
          buffer_out[n_elt_tgt * j + i] += w * buffer_in[tgt_to_src_idx[n_computed_tgt] * j + k];
        }
      }
      else {
        for (int j = 0; j < n_components; j++) {
          buffer_out[n_components * i + j] += w * buffer_in[n_components * k + j];
        }
      }
    }

    sum_w = 1. / sum_w;
    if (storage == CWP_FIELD_STORAGE_INTERLEAVED) {
      for (int j = 0; j < n_components; j++) {
        buffer_out[n_elt_tgt * j + i] *= sum_w;
      }
    }
    else {
      for (int j = 0; j < n_components; j++) {
        buffer_out[n_components * i + j] *= sum_w;
      }
    }
  }
}

static void
_location_user_interpolation
(
  const char   *local_code_name, // Local code name
  const char   *cpl_id,          // Coupling identifier
  const char   *field_id,        // Field identifier
        int     i_part,          // Partition
        double *buffer_in,       // Input buffer
        double *buffer_out       // Output buffer
)
{
  // GENERAL GEOMETRIC DATA

  // Get the number of components of the field
  int n_components = CWP_Field_n_components_get(local_code_name,
                                                cpl_id,
                                                field_id);

  // Get the field data storage type
  CWP_Field_storage_t storage = CWP_Field_storage_get(local_code_name,
                                                      cpl_id,
                                                      field_id);

  // Get the location of the degrees of freedom
  CWP_Dof_location_t dof_location = CWP_Field_dof_location_get(local_code_name,
                                                               cpl_id,
                                                               field_id);

  // Get data from the target of the exchange
  int n_elt_tgt      = 0;
  int n_computed_tgt = 0;

  int *computed_tgt   = nullptr;
  int *tgt_to_src_idx = nullptr;

  CWP_Field_tgt_data_properties_get(local_code_name,
                                    cpl_id,
                                    field_id,
                                    i_part,
                                    &n_elt_tgt,
                                    &n_computed_tgt,
                                    &computed_tgt,
                                    &tgt_to_src_idx);

  // LOCATION SPATIAL INTERPOLATION ALGORITHM

  // Get the spatial interpolation internal cell->vertex connectivity
  int *cell_vtx_idx = nullptr;
  int *cell_vtx     = nullptr;

  CWP_Field_location_internal_cell_vtx_get(local_code_name,
                                           cpl_id,
                                           field_id,
                                           i_part,
                                           &cell_vtx_idx,
                                           &cell_vtx);

  // Get the weights computed by the location spatial interpolation algorithm
  double *weights = nullptr;

  CWP_Field_location_weights_get(local_code_name,
                                 cpl_id,
                                 field_id,
                                 i_part,
                                 &weights);

  // Fill buffer out
  for (int i = 0; i < n_components * n_elt_tgt; i++) {
    buffer_out[i] = 0.;
  }

  if (dof_location == CWP_DOF_LOCATION_CELL_CENTER) {
    for (int i_ref = 0; i_ref < n_computed_tgt; i_ref++) {
      int i = computed_tgt[i_ref] - 1;

      for (int k = tgt_to_src_idx[i_ref]; k < tgt_to_src_idx[i_ref + 1]; k++) {
        if (storage == CWP_FIELD_STORAGE_INTERLEAVED) {
          for (int j = 0; j < n_components; j++) {
            buffer_out[n_elt_tgt * j + i] = buffer_in[tgt_to_src_idx[n_computed_tgt] * j + k];
          }
        }
        else {
          for (int j = 0; j < n_components; j++) {
            buffer_out[n_components * i + j] = buffer_in[n_components * k + j];
          }
        }
      }
    }
  }
  else {
    for (int i_ref = 0; i_ref < n_computed_tgt; i_ref++) {
      int i = computed_tgt[i_ref] - 1;

      for (int k = tgt_to_src_idx[i_ref]; k < tgt_to_src_idx[i_ref + 1]; k++) {
        for (int l = cell_vtx_idx[k]; l < cell_vtx_idx[k + 1]; l++) {
          int i_vtx = cell_vtx[l] - 1;
          if (storage == CWP_FIELD_STORAGE_INTERLEAVED) {
            for (int j = 0; j < n_components; j++) {
              buffer_out[n_elt_tgt * j + i] += weights[l] * buffer_in[tgt_to_src_idx[n_computed_tgt] * j + i_vtx];
            }
          }
          else {
            for (int j = 0; j < n_components; j++) {
              buffer_out[n_components * i + j] += weights[l] * buffer_in[n_components * i_vtx + j];
            }
          }
        }
      }
    }
  }
}

static void
_gen_mesh
(
  const MPI_Comm             comm,           // MPI communicator
  const int                  n_part,         // Number of partitions
  const CWPT_split_dual_t    part_method,    // Partitioning method
  const CWP_g_num_t          n_vtx_seg,      // Number of vertices on square side
        int                **pn_face,        // Number of faces
        int                **pn_vtx,         // Number of vertices
        int               ***pface_vtx_idx,  // Face vertex connectivity index
        int               ***pface_vtx,      // Face vertex connectivity
        double            ***pvtx_coord,     // Vertex coordinates
        CWP_g_num_t       ***pface_ln_to_gn, // Face global numbering
        CWP_g_num_t       ***pvtx_ln_to_gn   // Vertex global numbering
)
{
  int          *pn_edge        = nullptr;
  int         **pface_edge     = nullptr;
  int         **pedge_vtx      = nullptr;
  CWP_g_num_t **pedge_ln_to_gn = nullptr;
  CWPT_generate_mesh_rectangle_ngon(comm,
                                    static_cast<CWPT_Mesh_nodal_elt_t>(CWP_BLOCK_FACE_POLY),
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

static double
_exchange_field
(
  const char                **order,          // Order [send code, recv code]
        int                   n_code,         // Number of codes
  const char                **code_names,     // Code names
  const char                 *cpl_name,       // Coupling name
  const char                 *field_name,     // Field name
        CWP_Dof_location_t    loc_send,       // DOF location for send code
        int                   n_part,         // Number of partitions
        int                 **pn_vtx,         // Number of vertices
        int                 **pn_face,        // Number of faces
        CWP_g_num_t        ***pvtx_ln_to_gn,  // Vertex global numbering
        CWP_g_num_t        ***pface_ln_to_gn, // Face global numbering
        double             ***send_val,       // Send buffers
        double             ***recv_val,       // Receive buffers
        bool                  interleaved     // Interleaved storage
)
{
  // Send source values
  for (int i_code = 0; i_code < n_code; i_code++) {
    if (strcmp(code_names[i_code], order[0]) == 0) {
      for (int i_part = 0; i_part < n_part; i_part++) {
        if (loc_send == CWP_DOF_LOCATION_CELL_CENTER) {
          for (int i = 0; i < pn_face[i_code][i_part]; i++) {
            if (interleaved) {
              send_val[i_code][i_part][i]                           = static_cast<double>(pface_ln_to_gn[i_code][i_part][i]);
              send_val[i_code][i_part][i + pn_face[i_code][i_part]] = static_cast<double>(-2 * pface_ln_to_gn[i_code][i_part][i]);
            }
            else {
              send_val[i_code][i_part][2 * i]     = static_cast<double>(pface_ln_to_gn[i_code][i_part][i]);
              send_val[i_code][i_part][2 * i + 1] = static_cast<double>(-2 * pface_ln_to_gn[i_code][i_part][i]);
            }
          }
        }
        else {
          for (int i = 0; i < pn_vtx[i_code][i_part]; i++) {
            if (interleaved) {
              send_val[i_code][i_part][i]                          = static_cast<double>(pvtx_ln_to_gn[i_code][i_part][i]);
              send_val[i_code][i_part][i + pn_vtx[i_code][i_part]] = static_cast<double>(-2 * pvtx_ln_to_gn[i_code][i_part][i]);
            }
            else {
              send_val[i_code][i_part][2 * i]     = static_cast<double>(pvtx_ln_to_gn[i_code][i_part][i]);
              send_val[i_code][i_part][2 * i + 1] = static_cast<double>(-2 * pvtx_ln_to_gn[i_code][i_part][i]);
            }
          }
        }
      }
    }
  }

  for (int i_code = 0; i_code < n_code; i_code++) {
    if (strcmp(code_names[i_code], order[0]) == 0) {

      CWP_Field_issend(code_names[i_code], cpl_name, field_name);
    }
    if (strcmp(code_names[i_code], order[1]) == 0) {
      CWP_Field_irecv (code_names[i_code], cpl_name, field_name);
    }

    if (strcmp(code_names[i_code], order[0]) == 0) {
      CWP_Field_wait_issend(code_names[i_code], cpl_name, field_name);
    }
    if (strcmp(code_names[i_code], order[1]) == 0) {
      CWP_Field_wait_irecv (code_names[i_code], cpl_name, field_name);
    }
  }

  // Accumulate received values on target
  double res = 0.;
  for (int i_code = 0; i_code < n_code; i_code++) {
    if (strcmp(code_names[i_code], order[1]) == 0) {
      for (int i_part = 0; i_part < n_part; i_part++) {
        int n_loc = CWP_N_computed_tgts_get(code_names[i_code], cpl_name, field_name, i_part);
        for (int i = 0; i < n_loc; i++) {
          if (interleaved) {
            res += recv_val[i_code][i_part][i];
          }
          else {
            res += recv_val[i_code][i_part][2 * i];
          }
        }
      }
    }
  }

  return res;
}

static double
_get_value_cpl
(
  const char                **order,         // Order [send code, recv code]
        int                   n_code,        // Number of codes
  const char                **code_names,    // Code names
  const char                 *cpl_name,      // Coupling name
  const char                 *field_name1,   // Field name for exchange 1 (interlaced)
  const char                 *field_name2,   // Field name for exchange 2 (interleaved)
        CWP_Dof_location_t    loc_send,      // DOF location for send code
        CWP_Dof_location_t    loc_recv,      // DOF location for recv code
        int                   n_part,        // Number of partitions
        int                 **pn_vtx,        // Number of vertices
        int                 **pn_face,       // Number of faces
        CWP_g_num_t        ***pvtx_ln_to_gn, // Vertex global numbering
        CWP_g_num_t        ***pface_ln_to_gn // Face global numbering
)
{
  // Allocate buffers and set field data
  double ***send_val = (double ***) malloc(sizeof(double **) * n_code);
  double ***recv_val = (double ***) malloc(sizeof(double **) * n_code);
  for (int i_code = 0; i_code < n_code; i_code++) {
    if (strcmp(code_names[i_code], order[0]) == 0) {
      send_val[i_code] = (double **) malloc(sizeof(double *) * n_part);
      for (int i = 0; i < n_part; i++) {
        if (loc_send == CWP_DOF_LOCATION_CELL_CENTER) {
          send_val[i_code][i] = (double *) malloc(sizeof(double) * 2 * pn_face[i_code][i]);
        }
        else {
          send_val[i_code][i] = (double *) malloc(sizeof(double) * 2 * pn_vtx[i_code][i]);
        }

        CWP_Field_data_set(code_names[i_code],
                           cpl_name,
                           field_name1,
                           i,
                           CWP_FIELD_MAP_SOURCE,
                           send_val[i_code][i]);
        CWP_Field_data_set(code_names[i_code],
                           cpl_name,
                           field_name2,
                           i,
                           CWP_FIELD_MAP_SOURCE,
                           send_val[i_code][i]);
      }
    }

    if (strcmp(code_names[i_code], order[1]) == 0) {
      recv_val[i_code] = (double **) malloc(sizeof(double *) * n_part);
      for (int i = 0; i < n_part; i++) {
        if (loc_recv == CWP_DOF_LOCATION_CELL_CENTER) {
          recv_val[i_code][i] = (double *) malloc(sizeof(double) * 2 * pn_face[i_code][i]);
        }
        else {
          recv_val[i_code][i] = (double *) malloc(sizeof(double) * 2 * pn_vtx[i_code][i]);
        }

        CWP_Field_data_set(code_names[i_code],
                           cpl_name,
                           field_name1,
                           i,
                           CWP_FIELD_MAP_TARGET,
                           recv_val[i_code][i]);
        CWP_Field_data_set(code_names[i_code],
                           cpl_name,
                           field_name2,
                           i,
                           CWP_FIELD_MAP_TARGET,
                           recv_val[i_code][i]);
      }
    }
  }

  double res_interlaced = _exchange_field(order,
                                          n_code,
                                          code_names,
                                          cpl_name,
                                          field_name1,
                                          loc_send,
                                          n_part,
                                          pn_vtx,
                                          pn_face,
                                          pvtx_ln_to_gn,
                                          pface_ln_to_gn,
                                          send_val,
                                          recv_val,
                                          false);

  double res_interleaved = _exchange_field(order,
                                           n_code,
                                           code_names,
                                           cpl_name,
                                           field_name2,
                                           loc_send,
                                           n_part,
                                           pn_vtx,
                                           pn_face,
                                           pvtx_ln_to_gn,
                                           pface_ln_to_gn,
                                           send_val,
                                           recv_val,
                                           true);

  // Free buffers
  for (int i_code = 0; i_code < n_code; i_code++) {
    if (strcmp(code_names[i_code], order[0]) == 0) {
      for (int i_part = 0; i_part < n_part; i_part++) {
        free(send_val[i_code][i_part]);
      }
      free(send_val[i_code]);
    }
    if (strcmp(code_names[i_code], order[1]) == 0) {
      for (int i_part = 0; i_part < n_part; i_part++) {
        free(recv_val[i_code][i_part]);
      }
      free(recv_val[i_code]);
    }
  }
  free(send_val);
  free(recv_val);

  // Compare local interleaved and interlaced values
  CHECK(res_interlaced == doctest::Approx(res_interleaved).epsilon(_max_error));

  return res_interlaced;
}

static std::vector<double>
_get_values
(
  MPI_Comm             comm,           // MPI communicator
  CWP_Comm_t           part1,          // Type of partitioning support for code 1
  CWP_Comm_t           part2,          // Type of partitioning support for code 2
  bool                 share_ranks,    // Share ranks for codes 1 and 2
  CWP_Spatial_interp_t spatial_interp, // Spatial interpolation
  bool                 user_interp     // User interpolation
)
{
  std::vector<double> res_loc = {0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                 0., 0., 0., 0., 0., 0., 0., 0., 0.};
  std::vector<double> res(res_loc.size());

  int i_rank;
  MPI_Comm_rank(comm, &i_rank);
  int n_rank;
  MPI_Comm_size(comm, &n_rank);
  assert(n_rank <= 2);

  std::vector<const char *> code_names;
  CWP_Status_t is_active_rank = CWP_STATUS_ON;

  // 1 rank = 2 codes on current rank
  if (n_rank == 1) {
    code_names = {"code1", "code2"};
  }
  // 2 ranks = dispatch codes to share ranks if share_ranks == true depending on partitioning support
  else {
    if (!share_ranks && i_rank == 1) {
      code_names.push_back("code1");
    }
    else if (share_ranks) {
      code_names.push_back("code1");
    }
    if (!share_ranks && i_rank == 0) {
      code_names.push_back("code2");
    }
    else if (share_ranks) {
      code_names.push_back("code2");
    }
  }
  int n_code = static_cast<int>(code_names.size());

  MPI_Comm *intra_comm = (MPI_Comm *) malloc(sizeof(MPI_Comm) * n_code);
  CWP_Init(comm,
           n_code,
           code_names.data(),
           is_active_rank,
           intra_comm);

  bool f_cell_cell = false;
  bool f_cell_vtx  = false;
  bool f_cell_usr  = false;
  bool f_vtx_cell  = false;
  bool f_vtx_vtx   = false;
  bool f_vtx_usr   = false;
  bool f_usr_cell  = false;
  bool f_usr_vtx   = false;
  bool f_usr_usr   = false;
  bool b_cell_cell = false;
  bool b_cell_vtx  = false;
  bool b_cell_usr  = false;
  bool b_vtx_cell  = false;
  bool b_vtx_vtx   = false;
  bool b_vtx_usr   = false;
  bool b_usr_cell  = false;
  bool b_usr_vtx   = false;
  bool b_usr_usr   = false;
  switch (spatial_interp) {
  case CWP_SPATIAL_INTERP_FROM_IDENTITY:
    f_cell_cell = true;
    f_vtx_vtx   = true;
    f_usr_usr   = true;
    b_cell_cell = true;
    b_vtx_vtx   = true;
    b_usr_usr   = true;
    break;
  case CWP_SPATIAL_INTERP_FROM_INTERSECTION:
    f_cell_cell = true;
    b_cell_cell = true;
    break;
  case CWP_SPATIAL_INTERP_FROM_NEAREST_SOURCES_LEAST_SQUARES:
  case CWP_SPATIAL_INTERP_FROM_NEAREST_TARGETS_LEAST_SQUARES:
    f_cell_cell = true;
    f_cell_vtx  = true;
    f_cell_usr  = true;
    f_vtx_cell  = true;
    f_vtx_vtx   = true;
    f_vtx_usr   = true;
    f_usr_cell  = true;
    f_usr_vtx   = true;
    f_usr_usr   = true;
    b_cell_cell = true;
    b_cell_vtx  = true;
    b_cell_usr  = true;
    b_vtx_cell  = true;
    b_vtx_vtx   = true;
    b_vtx_usr   = true;
    b_usr_cell  = true;
    b_usr_vtx   = true;
    b_usr_usr   = true;
    break;
  case CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE:
  case CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_BOXTREE:
  case CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_LOCATE_ALL_TGT:
    f_cell_cell = true;
    f_cell_vtx  = true;
    f_cell_usr  = true;
    f_vtx_cell  = true;
    f_vtx_vtx   = true;
    f_vtx_usr   = true;
    b_cell_cell = true;
    b_cell_vtx  = true;
    b_vtx_cell  = true;
    b_vtx_vtx   = true;
    b_usr_cell  = true;
    b_usr_vtx   = true;
    break;
  }

  CWP_Interp_function_t user_interpolation = nullptr;
  if (user_interp) {
    switch (spatial_interp) {
    case CWP_SPATIAL_INTERP_FROM_INTERSECTION:
      user_interpolation = _intersection_user_interpolation;
      break;
    case CWP_SPATIAL_INTERP_FROM_NEAREST_TARGETS_LEAST_SQUARES:
      user_interpolation = _nearest_user_interpolation;
      break;
    case CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE:
      user_interpolation = _location_user_interpolation;
      break;
    default:
      break;
    }
  }

  // Create coupling
  int         n_part   = 2;
  int         n_vtx    = 9;
  const char *cpl_name = "spatialInterpTest";

  for (int i_code = 0; i_code < n_code; i_code++) {
    const char *coupled_code_name = "code2";
    CWP_Comm_t part = part1;
    if (strcmp(code_names[i_code], "code2") == 0) {
      coupled_code_name = "code1";
      part = part2;
    }
    CWP_Cpl_create(code_names[i_code],
                   cpl_name,
                   coupled_code_name,
                   CWP_INTERFACE_SURFACE,
                   part,
                   spatial_interp,
                   n_part,
                   CWP_DYNAMIC_MESH_DEFORMABLE,
                   CWP_TIME_EXCH_USER_CONTROLLED);

  }

  // Define mesh
  int          **pn_face        = (int          **) malloc(sizeof(int          *) * n_code);
  int          **pn_vtx         = (int          **) malloc(sizeof(int          *) * n_code);
  int         ***pface_vtx_idx  = (int         ***) malloc(sizeof(int         **) * n_code);
  int         ***pface_vtx      = (int         ***) malloc(sizeof(int         **) * n_code);
  double      ***pvtx_coord     = (double      ***) malloc(sizeof(double      **) * n_code);
  CWP_g_num_t ***pface_ln_to_gn = (CWP_g_num_t ***) malloc(sizeof(CWP_g_num_t **) * n_code);
  CWP_g_num_t ***pvtx_ln_to_gn  = (CWP_g_num_t ***) malloc(sizeof(CWP_g_num_t **) * n_code);

  for (int i_code = 0; i_code < n_code; i_code++) {
    MPI_Comm mesh_comm;
    if ((part1 == CWP_COMM_PAR_WITHOUT_PART && strcmp(code_names[i_code], "code1") == 0) ||
        (part2 == CWP_COMM_PAR_WITHOUT_PART && strcmp(code_names[i_code], "code2") == 0)) {
      int color;
      MPI_Comm_rank(intra_comm[i_code], &color);
      MPI_Comm_split(intra_comm[i_code], color, color, &mesh_comm);
    }
    else {
      MPI_Comm_dup(intra_comm[i_code], &mesh_comm);
    }
    _gen_mesh(mesh_comm,
              n_part,
              CWPT_SPLIT_DUAL_WITH_HILBERT,
              n_vtx,
              &pn_face       [i_code],
              &pn_vtx        [i_code],
              &pface_vtx_idx [i_code],
              &pface_vtx     [i_code],
              &pvtx_coord    [i_code],
              &pface_ln_to_gn[i_code],
              &pvtx_ln_to_gn [i_code]);
    MPI_Comm_free(&mesh_comm);

    if (strcmp(code_names[i_code], "code1") == 0) {
      for (int i = 0; i < n_part; i++) {
        for (int i_vtx = 0; i_vtx < pn_vtx[i_code][i]; i_vtx++) {
          pvtx_coord[i_code][i][3*i_vtx+1] *= 1.5;
        }
      }
    }

    if (strcmp(code_names[i_code], "code2") == 0) {
      for (int i = 0; i < n_part; i++) {
        for (int i_vtx = 0; i_vtx < pn_vtx[i_code][i]; i_vtx++) {
          pvtx_coord[i_code][i][3*i_vtx] *= 1.5;
        }
      }
    }

    for (int i = 0; i < n_part; i++) {
      CWP_Mesh_interf_vtx_set(code_names[i_code],
                              cpl_name,
                              i,
                              pn_vtx       [i_code][i],
                              pvtx_coord   [i_code][i],
                              pvtx_ln_to_gn[i_code][i]);

      CWP_Mesh_interf_from_facevtx_set(code_names[i_code],
                                       cpl_name,
                                       i,
                                       pn_face       [i_code][i],
                                       pface_vtx_idx [i_code][i],
                                       pface_vtx     [i_code][i],
                                       pface_ln_to_gn[i_code][i]);

      CWP_User_tgt_pts_set(code_names[i_code],
                           cpl_name,
                           i,
                           pn_vtx       [i_code][i],
                           pvtx_coord   [i_code][i],
                           pvtx_ln_to_gn[i_code][i]);
    }

    CWP_Mesh_interf_finalize(code_names[i_code], cpl_name);
  }

  // Define fields
  std::vector<const char *> f_fields_interlaced;
  std::vector<const char *> f_fields_interleaved;
  std::vector<const char *> b_fields_interlaced;
  std::vector<const char *> b_fields_interleaved;
  std::vector<CWP_Dof_location_t> f_dof_location1;
  std::vector<CWP_Dof_location_t> f_dof_location2;
  std::vector<CWP_Dof_location_t> b_dof_location1;
  std::vector<CWP_Dof_location_t> b_dof_location2;
  if (f_cell_cell) {
    f_fields_interlaced .push_back("f_fieldinterlacedcellcell");
    f_fields_interleaved.push_back("f_fieldinterleavedcellcell");
    f_dof_location1     .push_back(CWP_DOF_LOCATION_CELL_CENTER);
    f_dof_location2     .push_back(CWP_DOF_LOCATION_CELL_CENTER);
  }
  if (b_cell_cell) {
    b_fields_interlaced .push_back("b_fieldinterlacedcellcell");
    b_fields_interleaved.push_back("b_fieldinterleavedcellcell");
    b_dof_location1     .push_back(CWP_DOF_LOCATION_CELL_CENTER);
    b_dof_location2     .push_back(CWP_DOF_LOCATION_CELL_CENTER);
  }
  if (f_cell_vtx) {
    f_fields_interlaced .push_back("f_fieldinterlacedcellvtx");
    f_fields_interleaved.push_back("f_fieldinterleavedcellvtx");
    f_dof_location1     .push_back(CWP_DOF_LOCATION_CELL_CENTER);
    f_dof_location2     .push_back(CWP_DOF_LOCATION_NODE);
  }
  if (b_cell_vtx) {
    b_fields_interlaced .push_back("b_fieldinterlacedcellvtx");
    b_fields_interleaved.push_back("b_fieldinterleavedcellvtx");
    b_dof_location1     .push_back(CWP_DOF_LOCATION_CELL_CENTER);
    b_dof_location2     .push_back(CWP_DOF_LOCATION_NODE);
  }
  if (f_cell_usr) {
    f_fields_interlaced .push_back("f_fieldinterlacedcellusr");
    f_fields_interleaved.push_back("f_fieldinterleavedcellusr");
    f_dof_location1     .push_back(CWP_DOF_LOCATION_CELL_CENTER);
    f_dof_location2     .push_back(CWP_DOF_LOCATION_USER);
  }
  if (b_cell_usr) {
    b_fields_interlaced .push_back("b_fieldinterlacedcellusr");
    b_fields_interleaved.push_back("b_fieldinterleavedcellusr");
    b_dof_location1     .push_back(CWP_DOF_LOCATION_CELL_CENTER);
    b_dof_location2     .push_back(CWP_DOF_LOCATION_USER);
  }
  if (f_vtx_cell) {
    f_fields_interlaced .push_back("f_fieldinterlacedvtxcell");
    f_fields_interleaved.push_back("f_fieldinterleavedvtxcell");
    f_dof_location1     .push_back(CWP_DOF_LOCATION_NODE);
    f_dof_location2     .push_back(CWP_DOF_LOCATION_CELL_CENTER);
  }
  if (b_vtx_cell) {
    b_fields_interlaced .push_back("b_fieldinterlacedvtxcell");
    b_fields_interleaved.push_back("b_fieldinterleavedvtxcell");
    b_dof_location1     .push_back(CWP_DOF_LOCATION_NODE);
    b_dof_location2     .push_back(CWP_DOF_LOCATION_CELL_CENTER);
  }
  if (f_vtx_vtx) {
    f_fields_interlaced .push_back("f_fieldinterlacedvtxvtx");
    f_fields_interleaved.push_back("f_fieldinterleavedvtxvtx");
    f_dof_location1     .push_back(CWP_DOF_LOCATION_NODE);
    f_dof_location2     .push_back(CWP_DOF_LOCATION_NODE);
  }
  if (b_vtx_vtx) {
    b_fields_interlaced .push_back("b_fieldinterlacedvtxvtx");
    b_fields_interleaved.push_back("b_fieldinterleavedvtxvtx");
    b_dof_location1     .push_back(CWP_DOF_LOCATION_NODE);
    b_dof_location2     .push_back(CWP_DOF_LOCATION_NODE);
  }
  if (f_vtx_usr) {
    f_fields_interlaced .push_back("f_fieldinterlacedvtxusr");
    f_fields_interleaved.push_back("f_fieldinterleavedvtxusr");
    f_dof_location1     .push_back(CWP_DOF_LOCATION_NODE);
    f_dof_location2     .push_back(CWP_DOF_LOCATION_USER);
  }
  if (b_vtx_usr) {
    b_fields_interlaced .push_back("b_fieldinterlacedvtxusr");
    b_fields_interleaved.push_back("b_fieldinterleavedvtxusr");
    b_dof_location1     .push_back(CWP_DOF_LOCATION_NODE);
    b_dof_location2     .push_back(CWP_DOF_LOCATION_USER);
  }
  if (f_usr_cell) {
    f_fields_interlaced .push_back("f_fieldinterlacedusrcell");
    f_fields_interleaved.push_back("f_fieldinterleavedusrcell");
    f_dof_location1     .push_back(CWP_DOF_LOCATION_USER);
    f_dof_location2     .push_back(CWP_DOF_LOCATION_CELL_CENTER);
  }
  if (b_usr_cell) {
    b_fields_interlaced .push_back("b_fieldinterlacedusrcell");
    b_fields_interleaved.push_back("b_fieldinterleavedusrcell");
    b_dof_location1     .push_back(CWP_DOF_LOCATION_USER);
    b_dof_location2     .push_back(CWP_DOF_LOCATION_CELL_CENTER);
  }
  if (f_usr_vtx) {
    f_fields_interlaced .push_back("f_fieldinterlacedusrvtx");
    f_fields_interleaved.push_back("f_fieldinterleavedusrvtx");
    f_dof_location1     .push_back(CWP_DOF_LOCATION_USER);
    f_dof_location2     .push_back(CWP_DOF_LOCATION_NODE);
  }
  if (b_usr_vtx) {
    b_fields_interlaced .push_back("b_fieldinterlacedusrvtx");
    b_fields_interleaved.push_back("b_fieldinterleavedusrvtx");
    b_dof_location1     .push_back(CWP_DOF_LOCATION_USER);
    b_dof_location2     .push_back(CWP_DOF_LOCATION_NODE);
  }
  if (f_usr_usr) {
    f_fields_interlaced .push_back("f_fieldinterlacedusrusr");
    f_fields_interleaved.push_back("f_fieldinterleavedusrusr");
    f_dof_location1     .push_back(CWP_DOF_LOCATION_USER);
    f_dof_location2     .push_back(CWP_DOF_LOCATION_USER);
  }
  if (b_usr_usr) {
    b_fields_interlaced .push_back("b_fieldinterlacedusrusr");
    b_fields_interleaved.push_back("b_fieldinterleavedusrusr");
    b_dof_location1     .push_back(CWP_DOF_LOCATION_USER);
    b_dof_location2     .push_back(CWP_DOF_LOCATION_USER);
  }

  for (int i_code = 0; i_code < n_code; i_code++) {
    if (strcmp(code_names[i_code], "code1") == 0) {
      for (size_t i_field = 0; i_field < f_fields_interlaced.size(); i_field++) {
        CWP_Field_create(code_names[i_code],
                         cpl_name,
                         f_fields_interlaced[i_field],
                         CWP_DOUBLE,
                         CWP_FIELD_STORAGE_INTERLACED,
                         2,
                         f_dof_location1[i_field],
                         CWP_FIELD_EXCH_SEND,
                         CWP_STATUS_ON);
        if (part1 == CWP_COMM_PAR_WITHOUT_PART) {
          CWP_Involved_srcs_bcast_enable(code_names[i_code],
                                         cpl_name,
                                         f_fields_interlaced[i_field]);
        }
      }
      for (size_t i_field = 0; i_field < f_fields_interleaved.size(); i_field++) {
        CWP_Field_create(code_names[i_code],
                         cpl_name,
                         f_fields_interleaved[i_field],
                         CWP_DOUBLE,
                         CWP_FIELD_STORAGE_INTERLEAVED,
                         2,
                         f_dof_location1[i_field],
                         CWP_FIELD_EXCH_SEND,
                         CWP_STATUS_ON);
        if (part1 == CWP_COMM_PAR_WITHOUT_PART) {
          CWP_Involved_srcs_bcast_enable(code_names[i_code],
                                         cpl_name,
                                         f_fields_interleaved[i_field]);
        }
      }
      for (size_t i_field = 0; i_field < b_fields_interlaced.size(); i_field++) {
        CWP_Field_create(code_names[i_code],
                         cpl_name,
                         b_fields_interlaced[i_field],
                         CWP_DOUBLE,
                         CWP_FIELD_STORAGE_INTERLACED,
                         2,
                         b_dof_location1[i_field],
                         CWP_FIELD_EXCH_RECV,
                         CWP_STATUS_ON);
        if (user_interp) {
          CWP_Field_interp_function_set(code_names[i_code],
                                        cpl_name,
                                        b_fields_interlaced[i_field],
                                        user_interpolation);
        }
        if (part1 == CWP_COMM_PAR_WITHOUT_PART) {
          CWP_Computed_tgts_bcast_enable(code_names[i_code],
                                         cpl_name,
                                         b_fields_interlaced[i_field]);
        }
      }
      for (size_t i_field = 0; i_field < b_fields_interleaved.size(); i_field++) {
        CWP_Field_create(code_names[i_code],
                         cpl_name,
                         b_fields_interleaved[i_field],
                         CWP_DOUBLE,
                         CWP_FIELD_STORAGE_INTERLEAVED,
                         2,
                         b_dof_location1[i_field],
                         CWP_FIELD_EXCH_RECV,
                         CWP_STATUS_ON);
        if (user_interp) {
          CWP_Field_interp_function_set(code_names[i_code],
                                        cpl_name,
                                        b_fields_interleaved[i_field],
                                        user_interpolation);
        }
        if (part1 == CWP_COMM_PAR_WITHOUT_PART) {
          CWP_Computed_tgts_bcast_enable(code_names[i_code],
                                         cpl_name,
                                         b_fields_interleaved[i_field]);
        }
      }
    }

    if (strcmp(code_names[i_code], "code2") == 0) {
      for (size_t i_field = 0; i_field < f_fields_interlaced.size(); i_field++) {
        CWP_Field_create(code_names[i_code],
                         cpl_name,
                         f_fields_interlaced[i_field],
                         CWP_DOUBLE,
                         CWP_FIELD_STORAGE_INTERLACED,
                         2,
                         f_dof_location2[i_field],
                         CWP_FIELD_EXCH_RECV,
                         CWP_STATUS_ON);
        if (user_interp) {
          CWP_Field_interp_function_set(code_names[i_code],
                                        cpl_name,
                                        f_fields_interlaced[i_field],
                                        user_interpolation);
        }
        if (part2 == CWP_COMM_PAR_WITHOUT_PART) {
          CWP_Computed_tgts_bcast_enable(code_names[i_code],
                                         cpl_name,
                                         f_fields_interlaced[i_field]);
        }
      }
      for (size_t i_field = 0; i_field < f_fields_interleaved.size(); i_field++) {
        CWP_Field_create(code_names[i_code],
                         cpl_name,
                         f_fields_interleaved[i_field],
                         CWP_DOUBLE,
                         CWP_FIELD_STORAGE_INTERLEAVED,
                         2,
                         f_dof_location2[i_field],
                         CWP_FIELD_EXCH_RECV,
                         CWP_STATUS_ON);
        if (user_interp) {
          CWP_Field_interp_function_set(code_names[i_code],
                                        cpl_name,
                                        f_fields_interleaved[i_field],
                                        user_interpolation);
        }
        if (part2 == CWP_COMM_PAR_WITHOUT_PART) {
          CWP_Computed_tgts_bcast_enable(code_names[i_code],
                                         cpl_name,
                                         f_fields_interleaved[i_field]);
        }
      }
      for (size_t i_field = 0; i_field < b_fields_interlaced.size(); i_field++) {
        CWP_Field_create(code_names[i_code],
                         cpl_name,
                         b_fields_interlaced[i_field],
                         CWP_DOUBLE,
                         CWP_FIELD_STORAGE_INTERLACED,
                         2,
                         b_dof_location2[i_field],
                         CWP_FIELD_EXCH_SEND,
                         CWP_STATUS_ON);
        if (part2 == CWP_COMM_PAR_WITHOUT_PART) {
          CWP_Involved_srcs_bcast_enable(code_names[i_code],
                                         cpl_name,
                                         b_fields_interlaced[i_field]);
        }
      }
      for (size_t i_field = 0; i_field < b_fields_interleaved.size(); i_field++) {
        CWP_Field_create(code_names[i_code],
                         cpl_name,
                         b_fields_interleaved[i_field],
                         CWP_DOUBLE,
                         CWP_FIELD_STORAGE_INTERLEAVED,
                         2,
                         b_dof_location2[i_field],
                         CWP_FIELD_EXCH_SEND,
                         CWP_STATUS_ON);
        if (part2 == CWP_COMM_PAR_WITHOUT_PART) {
          CWP_Involved_srcs_bcast_enable(code_names[i_code],
                                         cpl_name,
                                         b_fields_interleaved[i_field]);
        }
      }
    }

  }

  for (int i_step = 0; i_step < 2; i_step++) {

    // Begin time step
    for (int i_code = 0; i_code < n_code; i_code++) {
      CWP_Time_step_beg(code_names[i_code], static_cast<double>(i_step));
    }

    // Compute interpolation weights
    for (int i_code = 0; i_code < n_code; i_code++) {
      CWP_Spatial_interp_weights_compute(code_names[i_code], cpl_name);
    }

    if (f_cell_cell) {
      res_loc[ 0] = _get_value_cpl(_forward.data(),
                                   n_code,
                                   code_names.data(),
                                   cpl_name,
                                   "f_fieldinterlacedcellcell",
                                   "f_fieldinterleavedcellcell",
                                   CWP_DOF_LOCATION_CELL_CENTER,
                                   CWP_DOF_LOCATION_CELL_CENTER,
                                   n_part,
                                   pn_vtx,
                                   pn_face,
                                   pvtx_ln_to_gn,
                                   pface_ln_to_gn);
    }
    if (f_cell_vtx) {
      res_loc[ 1] = _get_value_cpl(_forward.data(),
                                   n_code,
                                   code_names.data(),
                                   cpl_name,
                                   "f_fieldinterlacedcellvtx",
                                   "f_fieldinterleavedcellvtx",
                                   CWP_DOF_LOCATION_CELL_CENTER,
                                   CWP_DOF_LOCATION_NODE,
                                   n_part,
                                   pn_vtx,
                                   pn_face,
                                   pvtx_ln_to_gn,
                                   pface_ln_to_gn);
    }
    if (f_cell_usr) {
      res_loc[ 2] = _get_value_cpl(_forward.data(),
                                   n_code,
                                   code_names.data(),
                                   cpl_name,
                                   "f_fieldinterlacedcellusr",
                                   "f_fieldinterleavedcellusr",
                                   CWP_DOF_LOCATION_CELL_CENTER,
                                   CWP_DOF_LOCATION_USER,
                                   n_part,
                                   pn_vtx,
                                   pn_face,
                                   pvtx_ln_to_gn,
                                   pface_ln_to_gn);
    }
    if (f_vtx_cell) {
      res_loc[ 3] = _get_value_cpl(_forward.data(),
                                   n_code,
                                   code_names.data(),
                                   cpl_name,
                                   "f_fieldinterlacedvtxcell",
                                   "f_fieldinterleavedvtxcell",
                                   CWP_DOF_LOCATION_NODE,
                                   CWP_DOF_LOCATION_CELL_CENTER,
                                   n_part,
                                   pn_vtx,
                                   pn_face,
                                   pvtx_ln_to_gn,
                                   pface_ln_to_gn);
    }
    if (f_vtx_vtx) {
      res_loc[ 4] = _get_value_cpl(_forward.data(),
                                   n_code,
                                   code_names.data(),
                                   cpl_name,
                                   "f_fieldinterlacedvtxvtx",
                                   "f_fieldinterleavedvtxvtx",
                                   CWP_DOF_LOCATION_NODE,
                                   CWP_DOF_LOCATION_NODE,
                                   n_part,
                                   pn_vtx,
                                   pn_face,
                                   pvtx_ln_to_gn,
                                   pface_ln_to_gn);
    }
    if (f_vtx_usr) {
      res_loc[ 5] = _get_value_cpl(_forward.data(),
                                   n_code,
                                   code_names.data(),
                                   cpl_name,
                                   "f_fieldinterlacedvtxusr",
                                   "f_fieldinterleavedvtxusr",
                                   CWP_DOF_LOCATION_NODE,
                                   CWP_DOF_LOCATION_USER,
                                   n_part,
                                   pn_vtx,
                                   pn_face,
                                   pvtx_ln_to_gn,
                                   pface_ln_to_gn);
    }
    if (f_usr_cell) {
      res_loc[ 6] = _get_value_cpl(_forward.data(),
                                   n_code,
                                   code_names.data(),
                                   cpl_name,
                                   "f_fieldinterlacedusrcell",
                                   "f_fieldinterleavedusrcell",
                                   CWP_DOF_LOCATION_USER,
                                   CWP_DOF_LOCATION_CELL_CENTER,
                                   n_part,
                                   pn_vtx,
                                   pn_face,
                                   pvtx_ln_to_gn,
                                   pface_ln_to_gn);
    }
    if (f_usr_vtx) {
      res_loc[ 7] = _get_value_cpl(_forward.data(),
                                   n_code,
                                   code_names.data(),
                                   cpl_name,
                                   "f_fieldinterlacedusrvtx",
                                   "f_fieldinterleavedusrvtx",
                                   CWP_DOF_LOCATION_USER,
                                   CWP_DOF_LOCATION_NODE,
                                   n_part,
                                   pn_vtx,
                                   pn_face,
                                   pvtx_ln_to_gn,
                                   pface_ln_to_gn);
    }
    if (f_usr_usr) {
      res_loc[ 8] = _get_value_cpl(_forward.data(),
                                   n_code,
                                   code_names.data(),
                                   cpl_name,
                                   "f_fieldinterlacedusrusr",
                                   "f_fieldinterleavedusrusr",
                                   CWP_DOF_LOCATION_USER,
                                   CWP_DOF_LOCATION_USER,
                                   n_part,
                                   pn_vtx,
                                   pn_face,
                                   pvtx_ln_to_gn,
                                   pface_ln_to_gn);
    }
    if (b_cell_cell) {
      res_loc[ 9] = _get_value_cpl(_backward.data(),
                                   n_code,
                                   code_names.data(),
                                   cpl_name,
                                   "b_fieldinterlacedcellcell",
                                   "b_fieldinterleavedcellcell",
                                   CWP_DOF_LOCATION_CELL_CENTER,
                                   CWP_DOF_LOCATION_CELL_CENTER,
                                   n_part,
                                   pn_vtx,
                                   pn_face,
                                   pvtx_ln_to_gn,
                                   pface_ln_to_gn);
    }
    if (b_cell_vtx) {
      res_loc[10] = _get_value_cpl(_backward.data(),
                                   n_code,
                                   code_names.data(),
                                   cpl_name,
                                   "b_fieldinterlacedcellvtx",
                                   "b_fieldinterleavedcellvtx",
                                   CWP_DOF_LOCATION_NODE,
                                   CWP_DOF_LOCATION_CELL_CENTER,
                                   n_part,
                                   pn_vtx,
                                   pn_face,
                                   pvtx_ln_to_gn,
                                   pface_ln_to_gn);
    }
    if (b_cell_usr) {
      res_loc[11] = _get_value_cpl(_backward.data(),
                                   n_code,
                                   code_names.data(),
                                   cpl_name,
                                   "b_fieldinterlacedcellusr",
                                   "b_fieldinterleavedcellusr",
                                   CWP_DOF_LOCATION_USER,
                                   CWP_DOF_LOCATION_CELL_CENTER,
                                   n_part,
                                   pn_vtx,
                                   pn_face,
                                   pvtx_ln_to_gn,
                                   pface_ln_to_gn);
    }
    if (b_vtx_cell) {
      res_loc[12] = _get_value_cpl(_backward.data(),
                                   n_code,
                                   code_names.data(),
                                   cpl_name,
                                   "b_fieldinterlacedvtxcell",
                                   "b_fieldinterleavedvtxcell",
                                   CWP_DOF_LOCATION_CELL_CENTER,
                                   CWP_DOF_LOCATION_NODE,
                                   n_part,
                                   pn_vtx,
                                   pn_face,
                                   pvtx_ln_to_gn,
                                   pface_ln_to_gn);
    }
    if (b_vtx_vtx) {
      res_loc[13] = _get_value_cpl(_backward.data(),
                                   n_code,
                                   code_names.data(),
                                   cpl_name,
                                   "b_fieldinterlacedvtxvtx",
                                   "b_fieldinterleavedvtxvtx",
                                   CWP_DOF_LOCATION_NODE,
                                   CWP_DOF_LOCATION_NODE,
                                   n_part,
                                   pn_vtx,
                                   pn_face,
                                   pvtx_ln_to_gn,
                                   pface_ln_to_gn);
    }
    if (b_vtx_usr) {
      res_loc[14] = _get_value_cpl(_backward.data(),
                                   n_code,
                                   code_names.data(),
                                   cpl_name,
                                   "b_fieldinterlacedvtxusr",
                                   "b_fieldinterleavedvtxusr",
                                   CWP_DOF_LOCATION_USER,
                                   CWP_DOF_LOCATION_NODE,
                                   n_part,
                                   pn_vtx,
                                   pn_face,
                                   pvtx_ln_to_gn,
                                   pface_ln_to_gn);
    }
    if (b_usr_cell) {
      res_loc[15] = _get_value_cpl(_backward.data(),
                                   n_code,
                                   code_names.data(),
                                   cpl_name,
                                   "b_fieldinterlacedusrcell",
                                   "b_fieldinterleavedusrcell",
                                   CWP_DOF_LOCATION_CELL_CENTER,
                                   CWP_DOF_LOCATION_USER,
                                   n_part,
                                   pn_vtx,
                                   pn_face,
                                   pvtx_ln_to_gn,
                                   pface_ln_to_gn);
    }
    if (b_usr_vtx) {
      res_loc[16] = _get_value_cpl(_backward.data(),
                                   n_code,
                                   code_names.data(),
                                   cpl_name,
                                   "b_fieldinterlacedusrvtx",
                                   "b_fieldinterleavedusrvtx",
                                   CWP_DOF_LOCATION_NODE,
                                   CWP_DOF_LOCATION_USER,
                                   n_part,
                                   pn_vtx,
                                   pn_face,
                                   pvtx_ln_to_gn,
                                   pface_ln_to_gn);
    }
    if (b_usr_usr) {
      res_loc[17] = _get_value_cpl(_backward.data(),
                                   n_code,
                                   code_names.data(),
                                   cpl_name,
                                   "b_fieldinterlacedusrusr",
                                   "b_fieldinterleavedusrusr",
                                   CWP_DOF_LOCATION_USER,
                                   CWP_DOF_LOCATION_USER,
                                   n_part,
                                   pn_vtx,
                                   pn_face,
                                   pvtx_ln_to_gn,
                                   pface_ln_to_gn);
    }
    MPI_Allreduce(res_loc.data(), res.data(), res.size(), MPI_DOUBLE, MPI_SUM, comm);

    // End time step
    for (int i_code = 0; i_code < n_code; i_code++) {
      CWP_Time_step_end(code_names[i_code]);
    }

  }

  // Free all
  for (int i_code = 0; i_code < n_code; i_code++) {
    if (strcmp(code_names[i_code], "code1") == 0) {
      for (size_t i_field = 0; i_field < f_fields_interlaced.size(); i_field++) {
        CWP_Field_del(code_names[i_code],
                      cpl_name,
                      f_fields_interlaced[i_field]);
      }
      for (size_t i_field = 0; i_field < f_fields_interleaved.size(); i_field++) {
        CWP_Field_del(code_names[i_code],
                      cpl_name,
                      f_fields_interleaved[i_field]);
      }
      for (size_t i_field = 0; i_field < b_fields_interlaced.size(); i_field++) {
        CWP_Field_del(code_names[i_code],
                      cpl_name,
                      b_fields_interlaced[i_field]);
      }
      for (size_t i_field = 0; i_field < b_fields_interleaved.size(); i_field++) {
        CWP_Field_del(code_names[i_code],
                      cpl_name,
                      b_fields_interleaved[i_field]);
      }
    }

    if (strcmp(code_names[i_code], "code2") == 0) {
      for (size_t i_field = 0; i_field < f_fields_interlaced.size(); i_field++) {
        CWP_Field_del(code_names[i_code],
                      cpl_name,
                      f_fields_interlaced[i_field]);
      }
      for (size_t i_field = 0; i_field < f_fields_interleaved.size(); i_field++) {
        CWP_Field_del(code_names[i_code],
                      cpl_name,
                      f_fields_interleaved[i_field]);
      }
      for (size_t i_field = 0; i_field < b_fields_interlaced.size(); i_field++) {
        CWP_Field_del(code_names[i_code],
                      cpl_name,
                      b_fields_interlaced[i_field]);
      }
      for (size_t i_field = 0; i_field < b_fields_interleaved.size(); i_field++) {
        CWP_Field_del(code_names[i_code],
                      cpl_name,
                      b_fields_interleaved[i_field]);
      }
    }

    CWP_Mesh_interf_del(code_names[i_code], cpl_name);
    CWP_Cpl_del(code_names[i_code], cpl_name);

  }

  for (int i_code = 0; i_code < n_code; i_code++) {
    for (int i_part = 0; i_part < n_part; i_part++) {
      free(pface_vtx_idx [i_code][i_part]);
      free(pface_vtx     [i_code][i_part]);
      free(pvtx_coord    [i_code][i_part]);
      free(pface_ln_to_gn[i_code][i_part]);
      free(pvtx_ln_to_gn [i_code][i_part]);
    }
    free(pn_face       [i_code]);
    free(pn_vtx        [i_code]);
    free(pface_vtx_idx [i_code]);
    free(pface_vtx     [i_code]);
    free(pvtx_coord    [i_code]);
    free(pface_ln_to_gn[i_code]);
    free(pvtx_ln_to_gn [i_code]);
  }
  free(pn_face);
  free(pn_vtx);
  free(pface_vtx_idx);
  free(pface_vtx);
  free(pvtx_coord);
  free(pface_ln_to_gn);
  free(pvtx_ln_to_gn);

  free(intra_comm);

  CWP_Finalize();

  return res;
}

MPI_TEST_CASE("[spatialInterp] - 1p - Identity #1", 1) {
  CWP_Spatial_interp_t spatial_interp = CWP_SPATIAL_INTERP_FROM_IDENTITY;

  std::vector<double> res_expected;
  std::vector<double> res;

  res_expected = {10585., 0., 0., 0., 47770., 0., 0., 0., 47770., 10585., 0., 0., 0., 47770., 0., 0., 0., 47770.};
  res = _get_values(test_comm,
                    CWP_COMM_PAR_WITH_PART,
                    CWP_COMM_PAR_WITH_PART,
                    false,
                    spatial_interp,
                    false);
  CHECK_EPS(res.data(), res_expected.data(), res.size());

  res_expected = {10585., 0., 0., 0., 47770., 0., 0., 0., 47770., 10585., 0., 0., 0., 47770., 0., 0., 0., 47770.};
  res = _get_values(test_comm,
                    CWP_COMM_PAR_WITHOUT_PART,
                    CWP_COMM_PAR_WITHOUT_PART,
                    false,
                    spatial_interp,
                    false);
  CHECK_EPS(res.data(), res_expected.data(), res.size());

  res_expected = {10585., 0., 0., 0., 47770., 0., 0., 0., 47770., 10585., 0., 0., 0., 47770., 0., 0., 0., 47770.};
  res = _get_values(test_comm,
                    CWP_COMM_PAR_WITH_PART,
                    CWP_COMM_PAR_WITHOUT_PART,
                    false,
                    spatial_interp,
                    false);
  CHECK_EPS(res.data(), res_expected.data(), res.size());

  res_expected = {10585., 0., 0., 0., 47770., 0., 0., 0., 47770., 10585., 0., 0., 0., 47770., 0., 0., 0., 47770.};
  res = _get_values(test_comm,
                    CWP_COMM_PAR_WITHOUT_PART,
                    CWP_COMM_PAR_WITH_PART,
                    false,
                    spatial_interp,
                    false);
  CHECK_EPS(res.data(), res_expected.data(), res.size());

}

MPI_TEST_CASE("[spatialInterp] - 2p - Identity #2", 2) {
  CWP_Spatial_interp_t spatial_interp = CWP_SPATIAL_INTERP_FROM_IDENTITY;

  std::vector<double> res_expected;
  std::vector<double> res;

  res_expected = {10585., 0., 0., 0., 47770., 0., 0., 0., 47770., 10585., 0., 0., 0., 47770., 0., 0., 0., 47770.};
  res = _get_values(test_comm,
                    CWP_COMM_PAR_WITH_PART,
                    CWP_COMM_PAR_WITH_PART,
                    false,
                    spatial_interp,
                    false);
  CHECK_EPS(res.data(), res_expected.data(), res.size());

  res_expected = {10585., 0., 0., 0., 47770., 0., 0., 0., 47770., 10585., 0., 0., 0., 47770., 0., 0., 0., 47770.};
  res = _get_values(test_comm,
                    CWP_COMM_PAR_WITH_PART,
                    CWP_COMM_PAR_WITHOUT_PART,
                    false,
                    spatial_interp,
                    false);
  CHECK_EPS(res.data(), res_expected.data(), res.size());

  res_expected = {10585., 0., 0., 0., 47770., 0., 0., 0., 47770., 10585., 0., 0., 0., 47770., 0., 0., 0., 47770.};
  res = _get_values(test_comm,
                    CWP_COMM_PAR_WITHOUT_PART,
                    CWP_COMM_PAR_WITHOUT_PART,
                    false,
                    spatial_interp,
                    false);
  CHECK_EPS(res.data(), res_expected.data(), res.size());

  res_expected = {10585., 0., 0., 0., 47770., 0., 0., 0., 47770., 10585., 0., 0., 0., 47770., 0., 0., 0., 47770.};
  res = _get_values(test_comm,
                    CWP_COMM_PAR_WITHOUT_PART,
                    CWP_COMM_PAR_WITH_PART,
                    false,
                    spatial_interp,
                    false);
  CHECK_EPS(res.data(), res_expected.data(), res.size());

  res_expected = {10585., 0., 0., 0., 51194., 0., 0., 0., 51194., 10585., 0., 0., 0., 51194., 0., 0., 0., 51194.};
  res = _get_values(test_comm,
                    CWP_COMM_PAR_WITH_PART,
                    CWP_COMM_PAR_WITH_PART,
                    true,
                    spatial_interp,
                    false);
  CHECK_EPS(res.data(), res_expected.data(), res.size());

  res_expected = {21170., 0., 0., 0., 95540., 0., 0., 0., 95540., 10585., 0., 0., 0., 51194., 0., 0., 0., 51194.};
  res = _get_values(test_comm,
                    CWP_COMM_PAR_WITH_PART,
                    CWP_COMM_PAR_WITHOUT_PART,
                    true,
                    spatial_interp,
                    false);
  CHECK_EPS(res.data(), res_expected.data(), res.size());

  res_expected = {10585., 0., 0., 0., 51194., 0., 0., 0., 51194., 21170., 0., 0., 0., 95540., 0., 0., 0., 95540.};
  res = _get_values(test_comm,
                    CWP_COMM_PAR_WITHOUT_PART,
                    CWP_COMM_PAR_WITH_PART,
                    true,
                    spatial_interp,
                    false);
  CHECK_EPS(res.data(), res_expected.data(), res.size());

}

MPI_TEST_CASE("[spatialInter] - 1p - Intersection #1", 1) {
  CWP_Spatial_interp_t spatial_interp = CWP_SPATIAL_INTERP_FROM_INTERSECTION;

  std::vector<double> res_expected;
  std::vector<double> res;

  res_expected = {6206.880769, 0., 0., 0., 0., 0., 0., 0., 0., 7768.53956, 0., 0., 0., 0., 0., 0., 0., 0.};
  res = _get_values(test_comm,
                    CWP_COMM_PAR_WITH_PART,
                    CWP_COMM_PAR_WITH_PART,
                    false,
                    spatial_interp,
                    false);
  CHECK_EPS(res.data(), res_expected.data(), res.size());

  res_expected = {6206.880769, 0., 0., 0., 0., 0., 0., 0., 0., 7768.53956, 0., 0., 0., 0., 0., 0., 0., 0.};
  res = _get_values(test_comm,
                    CWP_COMM_PAR_WITHOUT_PART,
                    CWP_COMM_PAR_WITHOUT_PART,
                    false,
                    spatial_interp,
                    false);
  CHECK_EPS(res.data(), res_expected.data(), res.size());

  res_expected = {6206.880769, 0., 0., 0., 0., 0., 0., 0., 0., 7768.53956, 0., 0., 0., 0., 0., 0., 0., 0.};
  res = _get_values(test_comm,
                    CWP_COMM_PAR_WITH_PART,
                    CWP_COMM_PAR_WITHOUT_PART,
                    false,
                    spatial_interp,
                    false);
  CHECK_EPS(res.data(), res_expected.data(), res.size());

  res_expected = {6206.880769, 0., 0., 0., 0., 0., 0., 0., 0., 7768.53956, 0., 0., 0., 0., 0., 0., 0., 0.};
  res = _get_values(test_comm,
                    CWP_COMM_PAR_WITHOUT_PART,
                    CWP_COMM_PAR_WITH_PART,
                    false,
                    spatial_interp,
                    false);
  CHECK_EPS(res.data(), res_expected.data(), res.size());

}

MPI_TEST_CASE("[spatialInterp] - 2p - Intersection #2", 2) {
  CWP_Spatial_interp_t spatial_interp = CWP_SPATIAL_INTERP_FROM_INTERSECTION;

  std::vector<double> res_expected;
  std::vector<double> res;

  res_expected = {6206.880769, 0., 0., 0., 0., 0., 0., 0., 0., 7768.53956, 0., 0., 0., 0., 0., 0., 0., 0.};
  res = _get_values(test_comm,
                    CWP_COMM_PAR_WITH_PART,
                    CWP_COMM_PAR_WITH_PART,
                    false,
                    spatial_interp,
                    false);
  CHECK_EPS(res.data(), res_expected.data(), res.size());

  res_expected = {6206.880769, 0., 0., 0., 0., 0., 0., 0., 0., 7768.53956, 0., 0., 0., 0., 0., 0., 0., 0.};
  res = _get_values(test_comm,
                    CWP_COMM_PAR_WITH_PART,
                    CWP_COMM_PAR_WITHOUT_PART,
                    false,
                    spatial_interp,
                    false);
  CHECK_EPS(res.data(), res_expected.data(), res.size());

  res_expected = {6206.880769, 0., 0., 0., 0., 0., 0., 0., 0., 7768.53956, 0., 0., 0., 0., 0., 0., 0., 0.};
  res = _get_values(test_comm,
                    CWP_COMM_PAR_WITHOUT_PART,
                    CWP_COMM_PAR_WITHOUT_PART,
                    false,
                    spatial_interp,
                    false);
  CHECK_EPS(res.data(), res_expected.data(), res.size());

  res_expected = {6206.880769, 0., 0., 0., 0., 0., 0., 0., 0., 7768.53956, 0., 0., 0., 0., 0., 0., 0., 0.};
  res = _get_values(test_comm,
                    CWP_COMM_PAR_WITHOUT_PART,
                    CWP_COMM_PAR_WITH_PART,
                    false,
                    spatial_interp,
                    false);
  CHECK_EPS(res.data(), res_expected.data(), res.size());

  res_expected = {7164.30989, 0., 0., 0., 0., 0., 0., 0., 0., 7664.778022, 0., 0., 0., 0., 0., 0., 0., 0.};
  res = _get_values(test_comm,
                    CWP_COMM_PAR_WITH_PART,
                    CWP_COMM_PAR_WITH_PART,
                    true,
                    spatial_interp,
                    false);
  CHECK_EPS(res.data(), res_expected.data(), res.size());

  res_expected = {12413.761538, 0., 0., 0., 0., 0., 0., 0., 0., 7664.778022, 0., 0., 0., 0., 0., 0., 0., 0.};
  res = _get_values(test_comm,
                    CWP_COMM_PAR_WITH_PART,
                    CWP_COMM_PAR_WITHOUT_PART,
                    true,
                    spatial_interp,
                    false);
  CHECK_EPS(res.data(), res_expected.data(), res.size());

  res_expected = {7164.30989, 0., 0., 0., 0., 0., 0., 0., 0., 15537.079121, 0., 0., 0., 0., 0., 0., 0., 0.};
  res = _get_values(test_comm,
                    CWP_COMM_PAR_WITHOUT_PART,
                    CWP_COMM_PAR_WITH_PART,
                    true,
                    spatial_interp,
                    false);
  CHECK_EPS(res.data(), res_expected.data(), res.size());

}

MPI_TEST_CASE("[spatialInter] - 1p - LS Source #1", 1) {
  CWP_Spatial_interp_t spatial_interp = CWP_SPATIAL_INTERP_FROM_NEAREST_SOURCES_LEAST_SQUARES;

  std::vector<double> res_expected;
  std::vector<double> res;

  res_expected = {-15291.364656, -32281.27028, -32281.27028, 24777.532217, -631026.13617, -631026.13617, 24777.532217, -631026.13617, -631026.13617, -27985.490068, 21353.406116, 21353.406116, -52265.229202, 210588.336342, 210588.336342, -52265.229202, 210588.336342, 210588.336342};
  res = _get_values(test_comm,
                    CWP_COMM_PAR_WITH_PART,
                    CWP_COMM_PAR_WITH_PART,
                    false,
                    spatial_interp,
                    false);
  CHECK_EPS(res.data(), res_expected.data(), res.size());

  res_expected = {-15291.364656, -32281.27028, -32281.27028, 24777.532217, -631026.13617, -631026.13617, 24777.532217, -631026.13617, -631026.13617, -27985.490068, 21353.406116, 21353.406116, -52265.229202, 210588.336342, 210588.336342, -52265.229202, 210588.336342, 210588.336342};
  res = _get_values(test_comm,
                    CWP_COMM_PAR_WITHOUT_PART,
                    CWP_COMM_PAR_WITHOUT_PART,
                    false,
                    spatial_interp,
                    false);
  CHECK_EPS(res.data(), res_expected.data(), res.size());

  res_expected = {-15291.364656, -32281.27028, -32281.27028, 24777.532217, -631026.13617, -631026.13617, 24777.532217, -631026.13617, -631026.13617, -27985.490068, 21353.406116, 21353.406116, -52265.229202, 210588.336342, 210588.336342, -52265.229202, 210588.336342, 210588.336342};
  res = _get_values(test_comm,
                    CWP_COMM_PAR_WITH_PART,
                    CWP_COMM_PAR_WITHOUT_PART,
                    false,
                    spatial_interp,
                    false);
  CHECK_EPS(res.data(), res_expected.data(), res.size());

  res_expected = {-15291.364656, -32281.27028, -32281.27028, 24777.532217, -631026.13617, -631026.13617, 24777.532217, -631026.13617, -631026.13617, -27985.490068, 21353.406116, 21353.406116, -52265.229202, 210588.336342, 210588.336342, -52265.229202, 210588.336342, 210588.336342};
  res = _get_values(test_comm,
                    CWP_COMM_PAR_WITHOUT_PART,
                    CWP_COMM_PAR_WITH_PART,
                    false,
                    spatial_interp,
                    false);
  CHECK_EPS(res.data(), res_expected.data(), res.size());

}

MPI_TEST_CASE("[spatialInterp] - 2p - LS Source #2", 2) {
  CWP_Spatial_interp_t spatial_interp = CWP_SPATIAL_INTERP_FROM_NEAREST_SOURCES_LEAST_SQUARES;

  std::vector<double> res_expected;
  std::vector<double> res;

  res_expected = {-15291.364656, -32281.27028, -32281.27028, 24777.532217, -631026.13617, -631026.13617, 24777.532217, -631026.13617, -631026.13617, -27985.490068, 22509.631163, 22509.631163, -52265.229202, 210588.336342, 210588.336342, -52265.229202, 210588.336342, 210588.336342};
  res = _get_values(test_comm,
                    CWP_COMM_PAR_WITH_PART,
                    CWP_COMM_PAR_WITH_PART,
                    false,
                    spatial_interp,
                    false);
  CHECK_EPS(res.data(), res_expected.data(), res.size());

  res_expected = {-15291.364656, -32281.27028, -32281.27028, 24777.532217, -631026.13617, -631026.13617, 24777.532217, -631026.13617, -631026.13617, -27985.490068, 22509.631163, 22509.631163, -52265.229202, 210588.336342, 210588.336342, -52265.229202, 210588.336342, 210588.336342};
  res = _get_values(test_comm,
                    CWP_COMM_PAR_WITH_PART,
                    CWP_COMM_PAR_WITHOUT_PART,
                    false,
                    spatial_interp,
                    false);
  CHECK_EPS(res.data(), res_expected.data(), res.size());

  res_expected = {-15291.364656, -32281.27028, -32281.27028, 24777.532217, -631026.13617, -631026.13617, 24777.532217, -631026.13617, -631026.13617, -27985.490068, 22509.631163, 22509.631163, -52265.229202, 210588.336342, 210588.336342, -52265.229202, 210588.336342, 210588.336342};
  res = _get_values(test_comm,
                    CWP_COMM_PAR_WITHOUT_PART,
                    CWP_COMM_PAR_WITHOUT_PART,
                    false,
                    spatial_interp,
                    false);
  CHECK_EPS(res.data(), res_expected.data(), res.size());

  res_expected = {-15291.364656, -32281.27028, -32281.27028, 24777.532217, -631026.13617, -631026.13617, 24777.532217, -631026.13617, -631026.13617, -27985.490068, 22509.631163, 22509.631163, -52265.229202, 210588.336342, 210588.336342, -52265.229202, 210588.336342, 210588.336342};
  res = _get_values(test_comm,
                    CWP_COMM_PAR_WITHOUT_PART,
                    CWP_COMM_PAR_WITH_PART,
                    false,
                    spatial_interp,
                    false);
  CHECK_EPS(res.data(), res_expected.data(), res.size());

  res_expected = {-15291.364656, -31120.694313, -31120.694313, 24777.532217, -628635.796759, -628635.796759, 24777.532217, -628635.796759, -628635.796759, -27985.490068, 22509.631163, 22509.631163, -57113.597334, 200599.29343, 200599.29343, -57113.597334, 200599.29343, 200599.29343};
  res = _get_values(test_comm,
                    CWP_COMM_PAR_WITH_PART,
                    CWP_COMM_PAR_WITH_PART,
                    true,
                    spatial_interp,
                    false);
  CHECK_EPS(res.data(), res_expected.data(), res.size());

  res_expected = {-30582.729312, -64562.54056, -64562.54056, 49555.064434, -1262052.27234, -1262052.27234, 49555.064434, -1262052.27234, -1262052.27234, -27985.490068, 22509.631163, 24051.754111, -57113.597334, 200599.29343, 200599.54343, -57113.597334, 200599.29343, 200599.54343};
  res = _get_values(test_comm,
                    CWP_COMM_PAR_WITH_PART,
                    CWP_COMM_PAR_WITHOUT_PART,
                    true,
                    spatial_interp,
                    false);
  CHECK_EPS(res.data(), res_expected.data(), res.size());

  res_expected = {-15291.364656, -31120.694313, -31120.694313, 24777.532217, -628635.796759, -628635.796759, 25251.836315, -628639.296759, -628639.296759, -55970.980137, 45019.262326, 45019.262326, -104530.458403, 421176.672683, 421176.672683, -104530.458403, 421176.672683, 421176.672683};
  res = _get_values(test_comm,
                    CWP_COMM_PAR_WITHOUT_PART,
                    CWP_COMM_PAR_WITH_PART,
                    true,
                    spatial_interp,
                    false);
  CHECK_EPS(res.data(), res_expected.data(), res.size());

}

MPI_TEST_CASE("[spatialInter] - 1p - LS Target #1", 1) {
  CWP_Spatial_interp_t spatial_interp = CWP_SPATIAL_INTERP_FROM_NEAREST_TARGETS_LEAST_SQUARES;

  std::vector<double> res_expected;
  std::vector<double> res;

  res_expected = {4790.247942, 10648.084076, 10648.084076, 8455.259911, 15122.233416, 15122.233416, 8455.259911, 15122.233416, 15122.233416, 5929.407298, 12571.427454, 12571.427454, 17743.650767, 36475.093495, 36475.093495, 17743.650767, 36475.093495, 36475.093495};
  res = _get_values(test_comm,
                    CWP_COMM_PAR_WITH_PART,
                    CWP_COMM_PAR_WITH_PART,
                    false,
                    spatial_interp,
                    false);
  CHECK_EPS(res.data(), res_expected.data(), res.size());

  res_expected = {4790.247942, 10648.084076, 10648.084076, 8455.259911, 15122.233416, 15122.233416, 8455.259911, 15122.233416, 15122.233416, 5929.407298, 12571.427454, 12571.427454, 17743.650767, 36475.093495, 36475.093495, 17743.650767, 36475.093495, 36475.093495};
  res = _get_values(test_comm,
                    CWP_COMM_PAR_WITHOUT_PART,
                    CWP_COMM_PAR_WITHOUT_PART,
                    false,
                    spatial_interp,
                    false);
  CHECK_EPS(res.data(), res_expected.data(), res.size());

  res_expected = {4790.247942, 10648.084076, 10648.084076, 8455.259911, 15122.233416, 15122.233416, 8455.259911, 15122.233416, 15122.233416, 5929.407298, 12571.427454, 12571.427454, 17743.650767, 36475.093495, 36475.093495, 17743.650767, 36475.093495, 36475.093495};
  res = _get_values(test_comm,
                    CWP_COMM_PAR_WITH_PART,
                    CWP_COMM_PAR_WITHOUT_PART,
                    false,
                    spatial_interp,
                    false);
  CHECK_EPS(res.data(), res_expected.data(), res.size());

  res_expected = {4790.247942, 10648.084076, 10648.084076, 8455.259911, 15122.233416, 15122.233416, 8455.259911, 15122.233416, 15122.233416, 5929.407298, 12571.427454, 12571.427454, 17743.650767, 36475.093495, 36475.093495, 17743.650767, 36475.093495, 36475.093495};
  res = _get_values(test_comm,
                    CWP_COMM_PAR_WITHOUT_PART,
                    CWP_COMM_PAR_WITH_PART,
                    false,
                    spatial_interp,
                    false);
  CHECK_EPS(res.data(), res_expected.data(), res.size());

}

MPI_TEST_CASE("[spatialInterp] - 2p - LS Target #2", 2) {
  CWP_Spatial_interp_t spatial_interp = CWP_SPATIAL_INTERP_FROM_NEAREST_TARGETS_LEAST_SQUARES;

  std::vector<double> res_expected;
  std::vector<double> res;

  res_expected = {4790.247942, 10648.084076, 10648.084076, 8455.259911, 15122.233416, 15122.233416, 8455.259911, 15122.233416, 15122.233416, 5929.407298, 12571.47694, 12571.47694, 17743.650767, 36475.093495, 36475.093495, 17743.650767, 36475.093495, 36475.093495};
  res = _get_values(test_comm,
                    CWP_COMM_PAR_WITH_PART,
                    CWP_COMM_PAR_WITH_PART,
                    false,
                    spatial_interp,
                    false);
  CHECK_EPS(res.data(), res_expected.data(), res.size());

  res_expected = {4790.247942, 10648.084076, 10648.084076, 8455.259911, 15122.233416, 15122.233416, 8455.259911, 15122.233416, 15122.233416, 5929.407298, 12571.47694, 12571.47694, 17743.650767, 36475.093495, 36475.093495, 17743.650767, 36475.093495, 36475.093495};
  res = _get_values(test_comm,
                    CWP_COMM_PAR_WITH_PART,
                    CWP_COMM_PAR_WITHOUT_PART,
                    false,
                    spatial_interp,
                    false);
  CHECK_EPS(res.data(), res_expected.data(), res.size());

  res_expected = {4790.247942, 10648.084076, 10648.084076, 8455.259911, 15122.233416, 15122.233416, 8455.259911, 15122.233416, 15122.233416, 5929.407298, 12571.47694, 12571.47694, 17743.650767, 36475.093495, 36475.093495, 17743.650767, 36475.093495, 36475.093495};
  res = _get_values(test_comm,
                    CWP_COMM_PAR_WITHOUT_PART,
                    CWP_COMM_PAR_WITHOUT_PART,
                    false,
                    spatial_interp,
                    false);
  CHECK_EPS(res.data(), res_expected.data(), res.size());

  res_expected = {4790.247942, 10648.084076, 10648.084076, 8455.259911, 15122.233416, 15122.233416, 8455.259911, 15122.233416, 15122.233416, 5929.407298, 12571.47694, 12571.47694, 17743.650767, 36475.093495, 36475.093495, 17743.650767, 36475.093495, 36475.093495};
  res = _get_values(test_comm,
                    CWP_COMM_PAR_WITHOUT_PART,
                    CWP_COMM_PAR_WITH_PART,
                    false,
                    spatial_interp,
                    false);
  CHECK_EPS(res.data(), res_expected.data(), res.size());

  res_expected = {5615.997122, 14725.499347, 14725.499347, 9485.458443, 21618.607358, 21618.607358, 9485.458443, 21618.607358, 21618.607358, 5838.181807, 12705.446687, 12705.446687, 19160.664586, 38876.768729, 38876.768729, 19160.664586, 38876.768729, 38876.768729};
  res = _get_values(test_comm,
                    CWP_COMM_PAR_WITH_PART,
                    CWP_COMM_PAR_WITH_PART,
                    true,
                    spatial_interp,
                    false);
  CHECK_EPS(res.data(), res_expected.data(), res.size());

  res_expected = {9580.495883, 21296.168153, 21296.168153, 16910.519822, 30244.466832, 30244.466832, 16910.519822, 30244.466832, 30244.466832, 5838.181807, 12705.446687, 12705.446687, 19160.664586, 38876.768729, 38876.768729, 19160.664586, 38876.768729, 38876.768729};
  res = _get_values(test_comm,
                    CWP_COMM_PAR_WITH_PART,
                    CWP_COMM_PAR_WITHOUT_PART,
                    true,
                    spatial_interp,
                    false);
  CHECK_EPS(res.data(), res_expected.data(), res.size());

  res_expected = {5615.997122, 14725.499347, 14725.499347, 9485.458443, 21618.607358, 21618.607358, 9485.458443, 21618.607358, 21618.607358, 11858.814595, 25142.953881, 25142.953881, 35487.301533, 72950.186991, 72950.186991, 35494.277768, 72950.186991, 72950.186991};
  res = _get_values(test_comm,
                    CWP_COMM_PAR_WITHOUT_PART,
                    CWP_COMM_PAR_WITH_PART,
                    true,
                    spatial_interp,
                    false);
  CHECK_EPS(res.data(), res_expected.data(), res.size());

}

MPI_TEST_CASE("[spatialInter] - 1p - Location Octree #1", 1) {
  CWP_Spatial_interp_t spatial_interp = CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE;

  std::vector<double> res_expected;
  std::vector<double> res;

  res_expected = {7969., 18517., 18517., 9375.770798, 22829.288685, 22829.288685, 0., 0., 0., 8599., 12723.407406, 0., 20965., 35057.307515, 0., 20965., 35057.307515, 0.};
  res = _get_values(test_comm,
                    CWP_COMM_PAR_WITH_PART,
                    CWP_COMM_PAR_WITH_PART,
                    false,
                    spatial_interp,
                    false);
  CHECK_EPS(res.data(), res_expected.data(), res.size());

  res_expected = {7969., 18517., 18517., 9375.770798, 22829.288685, 22829.288685, 0., 0., 0., 8599., 12723.407406, 0., 20965., 35057.307515, 0., 20965., 35057.307515, 0.};
  res = _get_values(test_comm,
                    CWP_COMM_PAR_WITHOUT_PART,
                    CWP_COMM_PAR_WITHOUT_PART,
                    false,
                    spatial_interp,
                    false);
  CHECK_EPS(res.data(), res_expected.data(), res.size());

  res_expected = {7969., 18517., 18517., 9375.770798, 22829.288685, 22829.288685, 0., 0., 0., 8599., 12723.407406, 0., 20965., 35057.307515, 0., 20965., 35057.307515, 0.};
  res = _get_values(test_comm,
                    CWP_COMM_PAR_WITH_PART,
                    CWP_COMM_PAR_WITHOUT_PART,
                    false,
                    spatial_interp,
                    false);
  CHECK_EPS(res.data(), res_expected.data(), res.size());

  res_expected = {7969., 18517., 18517., 9375.770798, 22829.288685, 22829.288685, 0., 0., 0., 8599., 12723.407406, 0., 20965., 35057.307515, 0., 20965., 35057.307515, 0.};
  res = _get_values(test_comm,
                    CWP_COMM_PAR_WITHOUT_PART,
                    CWP_COMM_PAR_WITH_PART,
                    false,
                    spatial_interp,
                    false);
  CHECK_EPS(res.data(), res_expected.data(), res.size());

}

MPI_TEST_CASE("[spatialInterp] - 2p - Location Octree #2", 2) {
  CWP_Spatial_interp_t spatial_interp = CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE;

  std::vector<double> res_expected;
  std::vector<double> res;

  res_expected = {7969., 18307., 18307., 9375.770798, 22829.288685, 22829.288685, 0., 0., 0., 8608., 12723.407406, 0., 20901., 35057.307515, 0., 20901., 35057.307515, 0.};
  res = _get_values(test_comm,
                    CWP_COMM_PAR_WITH_PART,
                    CWP_COMM_PAR_WITH_PART,
                    false,
                    spatial_interp,
                    false);
  CHECK_EPS(res.data(), res_expected.data(), res.size());

  res_expected = {7969., 18307., 18307., 9375.770798, 22829.288685, 22829.288685, 0., 0., 0., 8608., 12723.407406, 0., 20901., 35057.307515, 0., 20901., 35057.307515, 0.};
  res = _get_values(test_comm,
                    CWP_COMM_PAR_WITH_PART,
                    CWP_COMM_PAR_WITHOUT_PART,
                    false,
                    spatial_interp,
                    false);
  CHECK_EPS(res.data(), res_expected.data(), res.size());

  res_expected = {7969., 18307., 18307., 9375.770798, 22829.288685, 22829.288685, 0., 0., 0., 8608., 12723.407406, 0., 20901., 35057.307515, 0., 20901., 35057.307515, 0.};
  res = _get_values(test_comm,
                    CWP_COMM_PAR_WITHOUT_PART,
                    CWP_COMM_PAR_WITHOUT_PART,
                    false,
                    spatial_interp,
                    false);
  CHECK_EPS(res.data(), res_expected.data(), res.size());

  res_expected = {7969., 18307., 18307., 9375.770798, 22829.288685, 22829.288685, 0., 0., 0., 8608., 12723.407406, 0., 20901., 35057.307515, 0., 20901., 35057.307515, 0.};
  res = _get_values(test_comm,
                    CWP_COMM_PAR_WITHOUT_PART,
                    CWP_COMM_PAR_WITH_PART,
                    false,
                    spatial_interp,
                    false);
  CHECK_EPS(res.data(), res_expected.data(), res.size());

  res_expected = {7969., 20012., 20012., 9375.770798, 24973.127738, 24973.127738, 0., 0., 0., 8608., 12723.407406, 0., 22664., 37180.15072, 0., 22664., 37180.15072, 0.};
  res = _get_values(test_comm,
                    CWP_COMM_PAR_WITH_PART,
                    CWP_COMM_PAR_WITH_PART,
                    true,
                    spatial_interp,
                    false);
  CHECK_EPS(res.data(), res_expected.data(), res.size());

  res_expected = {15938., 36614., 36614., 18751.541595, 45658.577369, 45658.577369, 0., 0., 0., 8608., 12723.407406, 0., 22664., 37180.15072, 0., 22664., 37180.15072, 0.};
  res = _get_values(test_comm,
                    CWP_COMM_PAR_WITH_PART,
                    CWP_COMM_PAR_WITHOUT_PART,
                    true,
                    spatial_interp,
                    false);
  CHECK_EPS(res.data(), res_expected.data(), res.size());

  res_expected = {7969., 20012., 20012., 9375.770798, 24973.127738, 24973.127738, 0., 0., 0., 17216., 25446.814812, 0., 41802., 70114.61503, 0., 41802., 70114.61503, 0.};
  res = _get_values(test_comm,
                    CWP_COMM_PAR_WITHOUT_PART,
                    CWP_COMM_PAR_WITH_PART,
                    true,
                    spatial_interp,
                    false);
  CHECK_EPS(res.data(), res_expected.data(), res.size());

}

MPI_TEST_CASE("[spatialInter] - 1p - Location Boxtree #1", 1) {
  CWP_Spatial_interp_t spatial_interp = CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_BOXTREE;

  std::vector<double> res_expected;
  std::vector<double> res;

  res_expected = {7969., 18464., 18464., 9375.770798, 22829.288685, 22829.288685, 0., 0., 0., 8651., 12723.407406, 0., 21086., 35057.307515, 0., 21086., 35057.307515, 0.};
  res = _get_values(test_comm,
                    CWP_COMM_PAR_WITH_PART,
                    CWP_COMM_PAR_WITH_PART,
                    false,
                    spatial_interp,
                    false);
  CHECK_EPS(res.data(), res_expected.data(), res.size());

  res_expected = {7969., 18464., 18464., 9375.770798, 22829.288685, 22829.288685, 0., 0., 0., 8651., 12723.407406, 0., 21086., 35057.307515, 0., 21086., 35057.307515, 0.};
  res = _get_values(test_comm,
                    CWP_COMM_PAR_WITHOUT_PART,
                    CWP_COMM_PAR_WITHOUT_PART,
                    false,
                    spatial_interp,
                    false);
  CHECK_EPS(res.data(), res_expected.data(), res.size());

  res_expected = {7969., 18464., 18464., 9375.770798, 22829.288685, 22829.288685, 0., 0., 0., 8651., 12723.407406, 0., 21086., 35057.307515, 0., 21086., 35057.307515, 0.};
  res = _get_values(test_comm,
                    CWP_COMM_PAR_WITH_PART,
                    CWP_COMM_PAR_WITHOUT_PART,
                    false,
                    spatial_interp,
                    false);
  CHECK_EPS(res.data(), res_expected.data(), res.size());

  res_expected = {7969., 18464., 18464., 9375.770798, 22829.288685, 22829.288685, 0., 0., 0., 8651., 12723.407406, 0., 21086., 35057.307515, 0., 21086., 35057.307515, 0.};
  res = _get_values(test_comm,
                    CWP_COMM_PAR_WITHOUT_PART,
                    CWP_COMM_PAR_WITH_PART,
                    false,
                    spatial_interp,
                    false);
  CHECK_EPS(res.data(), res_expected.data(), res.size());

}

MPI_TEST_CASE("[spatialInterp] - 2p - Location Boxtree #2", 2) {
  CWP_Spatial_interp_t spatial_interp = CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_BOXTREE;

  std::vector<double> res_expected;
  std::vector<double> res;

  res_expected = {7969., 18254., 18254., 9375.770798, 22829.288685, 22829.288685, 0., 0., 0., 8653., 12723.407406, 0., 21089., 35057.307515, 0., 21089., 35057.307515, 0.};
  res = _get_values(test_comm,
                    CWP_COMM_PAR_WITH_PART,
                    CWP_COMM_PAR_WITH_PART,
                    false,
                    spatial_interp,
                    false);
  CHECK_EPS(res.data(), res_expected.data(), res.size());

  res_expected = {7969., 18254., 18254., 9375.770798, 22829.288685, 22829.288685, 0., 0., 0., 8653., 12723.407406, 0., 21089., 35057.307515, 0., 21089., 35057.307515, 0.};
  res = _get_values(test_comm,
                    CWP_COMM_PAR_WITH_PART,
                    CWP_COMM_PAR_WITHOUT_PART,
                    false,
                    spatial_interp,
                    false);
  CHECK_EPS(res.data(), res_expected.data(), res.size());

  res_expected = {7969., 18254., 18254., 9375.770798, 22829.288685, 22829.288685, 0., 0., 0., 8653., 12723.407406, 0., 21089., 35057.307515, 0., 21089., 35057.307515, 0.};
  res = _get_values(test_comm,
                    CWP_COMM_PAR_WITHOUT_PART,
                    CWP_COMM_PAR_WITHOUT_PART,
                    false,
                    spatial_interp,
                    false);
  CHECK_EPS(res.data(), res_expected.data(), res.size());

  res_expected = {7969., 18254., 18254., 9375.770798, 22829.288685, 22829.288685, 0., 0., 0., 8653., 12723.407406, 0., 21089., 35057.307515, 0., 21089., 35057.307515, 0.};
  res = _get_values(test_comm,
                    CWP_COMM_PAR_WITHOUT_PART,
                    CWP_COMM_PAR_WITH_PART,
                    false,
                    spatial_interp,
                    false);
  CHECK_EPS(res.data(), res_expected.data(), res.size());

  res_expected = {7969., 19959., 19959., 9375.770798, 24973.127738, 24973.127738, 0., 0., 0., 8653., 12723.407406, 0., 22857., 37180.15072, 0., 22857., 37180.15072, 0.};
  res = _get_values(test_comm,
                    CWP_COMM_PAR_WITH_PART,
                    CWP_COMM_PAR_WITH_PART,
                    true,
                    spatial_interp,
                    false);
  CHECK_EPS(res.data(), res_expected.data(), res.size());

  res_expected = {15938., 36508., 36508., 18751.541595, 45658.577369, 45658.577369, 0., 0., 0., 8653., 12723.407406, 0., 22857., 37180.15072, 0., 22857., 37180.15072, 0.};
  res = _get_values(test_comm,
                    CWP_COMM_PAR_WITH_PART,
                    CWP_COMM_PAR_WITHOUT_PART,
                    true,
                    spatial_interp,
                    false);
  CHECK_EPS(res.data(), res_expected.data(), res.size());

  res_expected = {7969., 19959., 19959., 9375.770798, 24973.127738, 24973.127738, 0., 0., 0., 17306., 25446.814812, 0., 42178., 70114.61503, 0., 42178., 70114.61503, 0.};
  res = _get_values(test_comm,
                    CWP_COMM_PAR_WITHOUT_PART,
                    CWP_COMM_PAR_WITH_PART,
                    true,
                    spatial_interp,
                    false);
  CHECK_EPS(res.data(), res_expected.data(), res.size());

}

MPI_TEST_CASE("[spatialInter] - 1p - Location All Tgt #1", 1) {
  CWP_Spatial_interp_t spatial_interp = CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_LOCATE_ALL_TGT;

  std::vector<double> res_expected;
  std::vector<double> res;

  res_expected = {12967., 28575., 28575., 14512.766524, 33346.878428, 33346.878428, 0., 0., 0., 12430., 27465.343303, 0., 29420., 62459.43572, 0., 29420., 62459.43572, 0.};
  res = _get_values(test_comm,
                    CWP_COMM_PAR_WITH_PART,
                    CWP_COMM_PAR_WITH_PART,
                    false,
                    spatial_interp,
                    false);
  CHECK_EPS(res.data(), res_expected.data(), res.size());

  res_expected = {12967., 28575., 28575., 14512.766524, 33346.878428, 33346.878428, 0., 0., 0., 12430., 27465.343303, 0., 29420., 62459.43572, 0., 29420., 62459.43572, 0.};
  res = _get_values(test_comm,
                    CWP_COMM_PAR_WITHOUT_PART,
                    CWP_COMM_PAR_WITHOUT_PART,
                    false,
                    spatial_interp,
                    false);
  CHECK_EPS(res.data(), res_expected.data(), res.size());

  res_expected = {12967., 28575., 28575., 14512.766524, 33346.878428, 33346.878428, 0., 0., 0., 12430., 27465.343303, 0., 29420., 62459.43572, 0., 29420., 62459.43572, 0.};
  res = _get_values(test_comm,
                    CWP_COMM_PAR_WITH_PART,
                    CWP_COMM_PAR_WITHOUT_PART,
                    false,
                    spatial_interp,
                    false);
  CHECK_EPS(res.data(), res_expected.data(), res.size());

  res_expected = {12967., 28575., 28575., 14512.766524, 33346.878428, 33346.878428, 0., 0., 0., 12430., 27465.343303, 0., 29420., 62459.43572, 0., 29420., 62459.43572, 0.};
  res = _get_values(test_comm,
                    CWP_COMM_PAR_WITHOUT_PART,
                    CWP_COMM_PAR_WITH_PART,
                    false,
                    spatial_interp,
                    false);
  CHECK_EPS(res.data(), res_expected.data(), res.size());

}

MPI_TEST_CASE("[spatialInterp] - 2p - Location All Tgt #2", 2) {
  CWP_Spatial_interp_t spatial_interp = CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_LOCATE_ALL_TGT;

  std::vector<double> res_expected;
  std::vector<double> res;

  res_expected = {12911., 28418., 28418., 14512.766524, 33346.878428, 33346.878428, 0., 0., 0., 12430., 27465.343303, 0., 29186., 62459.43572, 0., 29186., 62459.43572, 0.};
  res = _get_values(test_comm,
                    CWP_COMM_PAR_WITH_PART,
                    CWP_COMM_PAR_WITH_PART,
                    false,
                    spatial_interp,
                    false);
  CHECK_EPS(res.data(), res_expected.data(), res.size());

  res_expected = {12911., 28418., 28418., 14512.766524, 33346.878428, 33346.878428, 0., 0., 0., 12430., 27465.343303, 0., 29186., 62459.43572, 0., 29186., 62459.43572, 0.};
  res = _get_values(test_comm,
                    CWP_COMM_PAR_WITH_PART,
                    CWP_COMM_PAR_WITHOUT_PART,
                    false,
                    spatial_interp,
                    false);
  CHECK_EPS(res.data(), res_expected.data(), res.size());

  res_expected = {12911., 28418., 28418., 14512.766524, 33346.878428, 33346.878428, 0., 0., 0., 12430., 27465.343303, 0., 29186., 62459.43572, 0., 29186., 62459.43572, 0.};
  res = _get_values(test_comm,
                    CWP_COMM_PAR_WITHOUT_PART,
                    CWP_COMM_PAR_WITHOUT_PART,
                    false,
                    spatial_interp,
                    false);
  CHECK_EPS(res.data(), res_expected.data(), res.size());

  res_expected = {12911., 28418., 28418., 14512.766524, 33346.878428, 33346.878428, 0., 0., 0., 12430., 27465.343303, 0., 29186., 62459.43572, 0., 29186., 62459.43572, 0.};
  res = _get_values(test_comm,
                    CWP_COMM_PAR_WITHOUT_PART,
                    CWP_COMM_PAR_WITH_PART,
                    false,
                    spatial_interp,
                    false);
  CHECK_EPS(res.data(), res_expected.data(), res.size());

  res_expected = {12911., 30486., 30486., 14512.766524, 35756.294405, 35756.294405, 0., 0., 0., 12430., 27465.343303, 0., 31756., 66969.073797, 0., 31756., 66969.073797, 0.};
  res = _get_values(test_comm,
                    CWP_COMM_PAR_WITH_PART,
                    CWP_COMM_PAR_WITH_PART,
                    true,
                    spatial_interp,
                    false);
  CHECK_EPS(res.data(), res_expected.data(), res.size());

  //res_expected = {25822., 56836., 56836., 29025.533048, 66693.756856, 66693.756856, 0., 0., 0., 12430., 27465.343303, 0., 31756., 66969.073797, 0., 31756., 66969.073797, 0.};
  //res = _get_values(test_comm,
  //                  CWP_COMM_PAR_WITH_PART,
  //                  CWP_COMM_PAR_WITHOUT_PART,
  //                  true,
  //                  spatial_interp,
  //                  false);
  //CHECK_EPS(res.data(), res_expected.data(), res.size());

  //res_expected = {12911., 30486., 30486., 14512.766524, 35756.294405, 35756.294405, 0., 0., 0., 24860., 54930.686607, 0., 58372., 124918.87144, 0., 58372., 124918.87144, 0.};
  //res = _get_values(test_comm,
  //                  CWP_COMM_PAR_WITHOUT_PART,
  //                  CWP_COMM_PAR_WITH_PART,
  //                  true,
  //                  spatial_interp,
  //                  false);
  //CHECK_EPS(res.data(), res_expected.data(), res.size());

}

MPI_TEST_CASE("[spatialInter] - 1p - User", 1) {
  std::vector<double> res_expected;
  std::vector<double> res;

  res_expected = {6206.880769, 0., 0., 0., 0., 0., 0., 0., 0., 7768.53956, 0., 0., 0., 0., 0., 0., 0., 0.};
  res = _get_values(test_comm,
                    CWP_COMM_PAR_WITH_PART,
                    CWP_COMM_PAR_WITH_PART,
                    false,
                    CWP_SPATIAL_INTERP_FROM_INTERSECTION,
                    true);
  CHECK_EPS(res.data(), res_expected.data(), res.size());

  res_expected = {4987.122141, 10648.084076, 10648.084076, 8422.621557, 15130.568446, 15130.568446, 8422.621557, 15130.568446, 15130.568446, 6187.854088, 12528.617381, 12528.617381, 17826.422146, 36437.872372, 36437.872372, 17826.422146, 36437.872372, 36437.872372};
  res = _get_values(test_comm,
                    CWP_COMM_PAR_WITH_PART,
                    CWP_COMM_PAR_WITH_PART,
                    false,
                    CWP_SPATIAL_INTERP_FROM_NEAREST_TARGETS_LEAST_SQUARES,
                    true);
  CHECK_EPS(res.data(), res_expected.data(), res.size());

  res_expected = {7969., 18517., 18517., 9375.770798, 22829.288685, 22829.288685, 0., 0., 0., 8599., 12723.407406, 0., 20965., 35057.307515, 0., 20965., 35057.307515, 0.};
  res = _get_values(test_comm,
                    CWP_COMM_PAR_WITH_PART,
                    CWP_COMM_PAR_WITH_PART,
                    false,
                    CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE,
                    true);
  CHECK_EPS(res.data(), res_expected.data(), res.size());

}