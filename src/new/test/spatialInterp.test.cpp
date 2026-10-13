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
                                    static_cast<CWPT_Mesh_nodal_elt_t>(CWP_BLOCK_FACE_QUAD4),
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
          pvtx_coord[i_code][i][3*i_vtx+1] += 0.3;
        }
      }
    }

    if (strcmp(code_names[i_code], "code2") == 0) {
      for (int i = 0; i < n_part; i++) {
        for (int i_vtx = 0; i_vtx < pn_vtx[i_code][i]; i_vtx++) {
          pvtx_coord[i_code][i][3*i_vtx] += 0.2;
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

  SUBCASE("WP - WP - no-shared") {
    res_expected = {2080., 0., 0., 0., 3690., 0., 0., 0., 3690., 2080., 0., 0., 0., 3690., 0., 0., 0., 3690.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITH_PART,
                      CWP_COMM_PAR_WITH_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WoP - WoP - no-shared") {
    res_expected = {2080., 0., 0., 0., 3690., 0., 0., 0., 3690., 2080., 0., 0., 0., 3690., 0., 0., 0., 3690.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITHOUT_PART,
                      CWP_COMM_PAR_WITHOUT_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WP - WoP - no-shared") {
    res_expected = {2080., 0., 0., 0., 3690., 0., 0., 0., 3690., 2080., 0., 0., 0., 3690., 0., 0., 0., 3690.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITH_PART,
                      CWP_COMM_PAR_WITHOUT_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WoP - WP - no-shared") {
    res_expected = {2080., 0., 0., 0., 3690., 0., 0., 0., 3690., 2080., 0., 0., 0., 3690., 0., 0., 0., 3690.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITHOUT_PART,
                      CWP_COMM_PAR_WITH_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

}

MPI_TEST_CASE("[spatialInterp] - 2p - Identity #2", 2) {
  CWP_Spatial_interp_t spatial_interp = CWP_SPATIAL_INTERP_FROM_IDENTITY;

  std::vector<double> res_expected;
  std::vector<double> res;

  SUBCASE("WP - WP - no-shared") {
    res_expected = {2080., 0., 0., 0., 3690., 0., 0., 0., 3690., 2080., 0., 0., 0., 3690., 0., 0., 0., 3690.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITH_PART,
                      CWP_COMM_PAR_WITH_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WP - WoP - no-shared") {
    res_expected = {2080., 0., 0., 0., 3690., 0., 0., 0., 3690., 2080., 0., 0., 0., 3690., 0., 0., 0., 3690.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITH_PART,
                      CWP_COMM_PAR_WITHOUT_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WoP - WoP - no-shared") {
    res_expected = {2080., 0., 0., 0., 3690., 0., 0., 0., 3690., 2080., 0., 0., 0., 3690., 0., 0., 0., 3690.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITHOUT_PART,
                      CWP_COMM_PAR_WITHOUT_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WoP - WP - no-shared") {
    res_expected = {2080., 0., 0., 0., 3690., 0., 0., 0., 3690., 2080., 0., 0., 0., 3690., 0., 0., 0., 3690.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITHOUT_PART,
                      CWP_COMM_PAR_WITH_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WP - WP - share") {
    res_expected = {2080., 0., 0., 0., 4100., 0., 0., 0., 4100., 2080., 0., 0., 0., 4100., 0., 0., 0., 4100.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITH_PART,
                      CWP_COMM_PAR_WITH_PART,
                      true,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WP - WoP - share") {
    res_expected = {4160., 0., 0., 0., 7380., 0., 0., 0., 7380., 2080., 0., 0., 0., 4100., 0., 0., 0., 4100.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITH_PART,
                      CWP_COMM_PAR_WITHOUT_PART,
                      true,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WoP - WP - share") {
    res_expected = {2080., 0., 0., 0., 4100., 0., 0., 0., 4100., 4160., 0., 0., 0., 7380., 0., 0., 0., 7380.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITHOUT_PART,
                      CWP_COMM_PAR_WITH_PART,
                      true,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

}

MPI_TEST_CASE("[spatialInterp] - 1p - Intersection #1", 1) {
  CWP_Spatial_interp_t spatial_interp = CWP_SPATIAL_INTERP_FROM_INTERSECTION;

  std::vector<double> res_expected;
  std::vector<double> res;

  SUBCASE("WP - WP - no-shared") {
    res_expected = {716.2, 0., 0., 0., 0., 0., 0., 0., 0., 1454.2, 0., 0., 0., 0., 0., 0., 0., 0.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITH_PART,
                      CWP_COMM_PAR_WITH_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WoP - WoP - no-shared") {
    res_expected = {716.2, 0., 0., 0., 0., 0., 0., 0., 0., 1454.2, 0., 0., 0., 0., 0., 0., 0., 0.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITHOUT_PART,
                      CWP_COMM_PAR_WITHOUT_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WP - WoP - no-shared") {
    res_expected = {716.2, 0., 0., 0., 0., 0., 0., 0., 0., 1454.2, 0., 0., 0., 0., 0., 0., 0., 0.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITH_PART,
                      CWP_COMM_PAR_WITHOUT_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WoP - WP - no-shared") {
    res_expected = {716.2, 0., 0., 0., 0., 0., 0., 0., 0., 1454.2, 0., 0., 0., 0., 0., 0., 0., 0.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITHOUT_PART,
                      CWP_COMM_PAR_WITH_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

}

MPI_TEST_CASE("[spatialInterp] - 2p - Intersection #2", 2) {
  CWP_Spatial_interp_t spatial_interp = CWP_SPATIAL_INTERP_FROM_INTERSECTION;

  std::vector<double> res_expected;
  std::vector<double> res;

  SUBCASE("WP - WP - no-shared") {
    res_expected = {716.2, 0., 0., 0., 0., 0., 0., 0., 0., 1454.2, 0., 0., 0., 0., 0., 0., 0., 0.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITH_PART,
                      CWP_COMM_PAR_WITH_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WP - WoP - no-shared") {
    res_expected = {716.2, 0., 0., 0., 0., 0., 0., 0., 0., 1454.2, 0., 0., 0., 0., 0., 0., 0., 0.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITH_PART,
                      CWP_COMM_PAR_WITHOUT_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WoP - WoP - no-shared") {
    res_expected = {716.2, 0., 0., 0., 0., 0., 0., 0., 0., 1454.2, 0., 0., 0., 0., 0., 0., 0., 0.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITHOUT_PART,
                      CWP_COMM_PAR_WITHOUT_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WoP - WP - no-shared") {
    res_expected = {716.2, 0., 0., 0., 0., 0., 0., 0., 0., 1454.2, 0., 0., 0., 0., 0., 0., 0., 0.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITHOUT_PART,
                      CWP_COMM_PAR_WITH_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WP - WP - share") {
    res_expected = {716.2, 0., 0., 0., 0., 0., 0., 0., 0., 1520.2, 0., 0., 0., 0., 0., 0., 0., 0.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITH_PART,
                      CWP_COMM_PAR_WITH_PART,
                      true,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WP - WoP - share") {
    res_expected = {1432.4, 0., 0., 0., 0., 0., 0., 0., 0., 1520.2, 0., 0., 0., 0., 0., 0., 0., 0.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITH_PART,
                      CWP_COMM_PAR_WITHOUT_PART,
                      true,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WoP - WP - share") {
    res_expected = {716.2, 0., 0., 0., 0., 0., 0., 0., 0., 2908.4, 0., 0., 0., 0., 0., 0., 0., 0.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITHOUT_PART,
                      CWP_COMM_PAR_WITH_PART,
                      true,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

}

MPI_TEST_CASE("[spatialInterp] - 1p - LS Source #1", 1) {
  CWP_Spatial_interp_t spatial_interp = CWP_SPATIAL_INTERP_FROM_NEAREST_SOURCES_LEAST_SQUARES;

  std::vector<double> res_expected;
  std::vector<double> res;

  SUBCASE("WP - WP - no-shared") {
    res_expected = {1562.547649, -5321.872299, -5321.872299, 1444.805946, 1663.986312, 1663.986312, 1444.805946, 1663.986312, 1663.986312, 9130.97676, 125280.353827, 125280.353827, -308141.3809, -35584.977189, -35584.977189, -308141.3809, -35584.977189, -35584.977189};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITH_PART,
                      CWP_COMM_PAR_WITH_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WoP - WoP - no-shared") {
    res_expected = {1562.547649, -5321.872299, -5321.872299, 1444.805946, 1663.986312, 1663.986312, 1444.805946, 1663.986312, 1663.986312, 9130.97676, 125280.353827, 125280.353827, -308141.3809, -35584.977189, -35584.977189, -308141.3809, -35584.977189, -35584.977189};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITHOUT_PART,
                      CWP_COMM_PAR_WITHOUT_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WP - WoP - no-shared") {
    res_expected = {1562.547649, -5321.872299, -5321.872299, 1444.805946, 1663.986312, 1663.986312, 1444.805946, 1663.986312, 1663.986312, 9130.97676, 125280.353827, 125280.353827, -308141.3809, -35584.977189, -35584.977189, -308141.3809, -35584.977189, -35584.977189};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITH_PART,
                      CWP_COMM_PAR_WITHOUT_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WoP - WP - no-shared") {
    res_expected = {1562.547649, -5321.872299, -5321.872299, 1444.805946, 1663.986312, 1663.986312, 1444.805946, 1663.986312, 1663.986312, 9130.97676, 125280.353827, 125280.353827, -308141.3809, -35584.977189, -35584.977189, -308141.3809, -35584.977189, -35584.977189};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITHOUT_PART,
                      CWP_COMM_PAR_WITH_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

}

MPI_TEST_CASE("[spatialInterp] - 2p - LS Source #2", 2) {
  CWP_Spatial_interp_t spatial_interp = CWP_SPATIAL_INTERP_FROM_NEAREST_SOURCES_LEAST_SQUARES;

  std::vector<double> res_expected;
  std::vector<double> res;

  SUBCASE("WP - WP - no-shared") {
    res_expected = {1562.547649, -5321.872299, -5321.872299, 1444.805946, 1663.986312, 1663.986312, 1444.805946, 1663.986312, 1663.986312, 9130.97676, 125280.353827, 125280.353827, -308141.3809, -35584.977189, -35584.977189, -308141.3809, -35584.977189, -35584.977189};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITH_PART,
                      CWP_COMM_PAR_WITH_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WP - WoP - no-shared") {
    res_expected = {1562.547649, -5321.872299, -5321.872299, 1444.805946, 1663.986312, 1663.986312, 1444.805946, 1663.986312, 1663.986312, 9130.97676, 125280.353827, 125280.353827, -308141.3809, -35584.977189, -35584.977189, -308141.3809, -35584.977189, -35584.977189};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITH_PART,
                      CWP_COMM_PAR_WITHOUT_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WoP - WoP - no-shared") {
    res_expected = {1562.547649, -5321.872299, -5321.872299, 1444.805946, 1663.986312, 1663.986312, 1444.805946, 1663.986312, 1663.986312, 9130.97676, 125280.353827, 125280.353827, -308141.3809, -35584.977189, -35584.977189, -308141.3809, -35584.977189, -35584.977189};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITHOUT_PART,
                      CWP_COMM_PAR_WITHOUT_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WoP - WP - no-shared") {
    res_expected = {1562.547649, -5321.872299, -5321.872299, 1444.805946, 1663.986312, 1663.986312, 1444.805946, 1663.986312, 1663.986312, 9130.97676, 125280.353827, 125280.353827, -308141.3809, -35584.977189, -35584.977189, -308141.3809, -35584.977189, -35584.977189};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITHOUT_PART,
                      CWP_COMM_PAR_WITH_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WP - WP - share") {
    res_expected = {1562.547649, -5121.089646, -5121.089646, 1444.805946, 2083.962699, 2083.962699, 1444.805946, 2083.962699, 2083.962699, 9130.97676, 125280.353827, 125280.353827, -307645.144413, -35017.643386, -35017.643386, -307645.144413, -35017.643386, -35017.643386};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITH_PART,
                      CWP_COMM_PAR_WITH_PART,
                      true,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WP - WoP - share") {
    res_expected = {3125.095297, -10643.744597, -10643.744597, 2889.611892, 3327.972624, 3327.972624, 2889.611892, 3327.972624, 3327.972624, 9130.97676, 125280.353827, 125280.353827, -307645.144413, -35017.643386, -35017.643386, -307645.144413, -35017.643386, -35017.643386};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITH_PART,
                      CWP_COMM_PAR_WITHOUT_PART,
                      true,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WoP - WP - share") {
    res_expected = {1562.547649, -5121.089646, -5121.089646, 1444.805946, 2083.962699, 2083.962699, 1444.805946, 2083.962699, 2083.962699, 18261.953521, 250560.707653, 250560.707653, -616282.761801, -71169.954378, -71169.954378, -616282.761801, -71169.954378, -71169.954378};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITHOUT_PART,
                      CWP_COMM_PAR_WITH_PART,
                      true,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

}

MPI_TEST_CASE("[spatialInterp] - 1p - LS Target #1", 1) {
  CWP_Spatial_interp_t spatial_interp = CWP_SPATIAL_INTERP_FROM_NEAREST_TARGETS_LEAST_SQUARES;

  std::vector<double> res_expected;
  std::vector<double> res;

  SUBCASE("WP - WP - no-shared") {
    res_expected = {760.796303, 733.219025, 733.219025, 1864.466344, 1397.168421, 1397.168421, 1864.466344, 1397.168421, 1397.168421, 1568.417827, 2346.856014, 2346.856014, 1751.720811, 3003.381699, 3003.381699, 1751.720811, 3003.381699, 3003.381699};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITH_PART,
                      CWP_COMM_PAR_WITH_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WoP - WoP - no-shared") {
    res_expected = {760.796303, 733.219025, 733.219025, 1864.466344, 1397.168421, 1397.168421, 1864.466344, 1397.168421, 1397.168421, 1568.417827, 2346.856014, 2346.856014, 1751.720811, 3003.381699, 3003.381699, 1751.720811, 3003.381699, 3003.381699};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITHOUT_PART,
                      CWP_COMM_PAR_WITHOUT_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WP - WoP - no-shared") {
    res_expected = {760.796303, 733.219025, 733.219025, 1864.466344, 1397.168421, 1397.168421, 1864.466344, 1397.168421, 1397.168421, 1568.417827, 2346.856014, 2346.856014, 1751.720811, 3003.381699, 3003.381699, 1751.720811, 3003.381699, 3003.381699};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITH_PART,
                      CWP_COMM_PAR_WITHOUT_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WoP - WP - no-shared") {
    res_expected = {760.796303, 733.219025, 733.219025, 1864.466344, 1397.168421, 1397.168421, 1864.466344, 1397.168421, 1397.168421, 1568.417827, 2346.856014, 2346.856014, 1751.720811, 3003.381699, 3003.381699, 1751.720811, 3003.381699, 3003.381699};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITHOUT_PART,
                      CWP_COMM_PAR_WITH_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

}

MPI_TEST_CASE("[spatialInterp] - 2p - LS Target #2", 2) {
  CWP_Spatial_interp_t spatial_interp = CWP_SPATIAL_INTERP_FROM_NEAREST_TARGETS_LEAST_SQUARES;

  std::vector<double> res_expected;
  std::vector<double> res;

  SUBCASE("WP - WP - no-shared") {
    res_expected = {760.796303, 733.909723, 733.909723, 1863.698903, 1397.168421, 1397.168421, 1863.698903, 1397.168421, 1397.168421, 1568.417827, 2346.856014, 2346.856014, 1749.566865, 3003.381699, 3003.381699, 1749.566865, 3003.381699, 3003.381699};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITH_PART,
                      CWP_COMM_PAR_WITH_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WP - WoP - no-shared") {
    res_expected = {760.796303, 733.909723, 733.909723, 1863.698903, 1397.168421, 1397.168421, 1863.698903, 1397.168421, 1397.168421, 1568.417827, 2346.856014, 2346.856014, 1749.566865, 3003.381699, 3003.381699, 1749.566865, 3003.381699, 3003.381699};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITH_PART,
                      CWP_COMM_PAR_WITHOUT_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WoP - WoP - no-shared") {
    res_expected = {760.796303, 733.909723, 733.909723, 1863.698903, 1397.168421, 1397.168421, 1863.698903, 1397.168421, 1397.168421, 1568.417827, 2346.856014, 2346.856014, 1749.566865, 3003.381699, 3003.381699, 1749.566865, 3003.381699, 3003.381699};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITHOUT_PART,
                      CWP_COMM_PAR_WITHOUT_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WoP - WP - no-shared") {
    res_expected = {760.796303, 733.909723, 733.909723, 1863.698903, 1397.168421, 1397.168421, 1863.698903, 1397.168421, 1397.168421, 1568.417827, 2346.856014, 2346.856014, 1749.566865, 3003.381699, 3003.381699, 1749.566865, 3003.381699, 3003.381699};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITHOUT_PART,
                      CWP_COMM_PAR_WITH_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WP - WP - share") {
    res_expected = {760.196303, 989.828328, 989.828328, 1867.698903, 1598.168421, 1598.168421, 1867.698903, 1598.168421, 1598.168421, 1529.501083, 2300.925165, 2300.925165, 2304.717065, 3378.336462, 3378.336462, 2304.717065, 3378.336462, 3378.336462};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITH_PART,
                      CWP_COMM_PAR_WITH_PART,
                      true,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WP - WoP - share") {
    res_expected = {1521.592605, 1467.819446, 1469.200841, 3727.397805, 2794.336842, 2794.336842, 3727.397805, 2794.336842, 2794.336842, 1529.501083, 2300.925165, 2300.925165, 2304.717065, 3378.336462, 3378.336462, 2304.717065, 3378.336462, 3378.336462};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITH_PART,
                      CWP_COMM_PAR_WITHOUT_PART,
                      true,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WoP - WP - share") {
    res_expected = {760.196303, 989.828328, 989.828328, 1867.698903, 1598.168421, 1598.168421, 1867.698903, 1598.168421, 1598.168421, 3136.835653, 4693.712029, 4693.712029, 3499.13373, 6006.763398, 6006.763398, 3499.13373, 6006.763398, 6006.763398};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITHOUT_PART,
                      CWP_COMM_PAR_WITH_PART,
                      true,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

}

MPI_TEST_CASE("[spatialInterp] - 1p - Location Octree #1", 1) {
  CWP_Spatial_interp_t spatial_interp = CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE;

  std::vector<double> res_expected;
  std::vector<double> res;

  SUBCASE("WP - WP - no-shared") {
    res_expected = {918., 1141., 1141., 1044., 1547., 1547., 0., 0., 0., 1422., 1908., 0., 2044., 2471., 0., 2044., 2471., 0.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITH_PART,
                      CWP_COMM_PAR_WITH_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WoP - WoP - no-shared") {
    res_expected = {918., 1141., 1141., 1044., 1547., 1547., 0., 0., 0., 1422., 1908., 0., 2044., 2471., 0., 2044., 2471., 0.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITHOUT_PART,
                      CWP_COMM_PAR_WITHOUT_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WP - WoP - no-shared") {
    res_expected = {918., 1141., 1141., 1044., 1547., 1547., 0., 0., 0., 1422., 1908., 0., 2044., 2471., 0., 2044., 2471., 0.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITH_PART,
                      CWP_COMM_PAR_WITHOUT_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WoP - WP - no-shared") {
    res_expected = {918., 1141., 1141., 1044., 1547., 1547., 0., 0., 0., 1422., 1908., 0., 2044., 2471., 0., 2044., 2471., 0.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITHOUT_PART,
                      CWP_COMM_PAR_WITH_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

}

MPI_TEST_CASE("[spatialInterp] - 2p - Location Octree #2", 2) {
  CWP_Spatial_interp_t spatial_interp = CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE;

  std::vector<double> res_expected;
  std::vector<double> res;

  SUBCASE("WP - WP - no-shared") {
    res_expected = {918., 1141., 1141., 1044., 1547., 1547., 0., 0., 0., 1422., 1908., 0., 2044., 2471., 0., 2044., 2471., 0.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITH_PART,
                      CWP_COMM_PAR_WITH_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WP - WoP - no-shared") {
    res_expected = {918., 1141., 1141., 1044., 1547., 1547., 0., 0., 0., 1422., 1908., 0., 2044., 2471., 0., 2044., 2471., 0.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITH_PART,
                      CWP_COMM_PAR_WITHOUT_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WoP - WoP - no-shared") {
    res_expected = {918., 1141., 1141., 1044., 1547., 1547., 0., 0., 0., 1422., 1908., 0., 2044., 2471., 0., 2044., 2471., 0.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITHOUT_PART,
                      CWP_COMM_PAR_WITHOUT_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WoP - WP - no-shared") {
    res_expected = {918., 1141., 1141., 1044., 1547., 1547., 0., 0., 0., 1422., 1908., 0., 2044., 2471., 0., 2044., 2471., 0.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITHOUT_PART,
                      CWP_COMM_PAR_WITH_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WP - WP - share") {
    res_expected = {918., 1311., 1311., 1044., 1775., 1775., 0., 0., 0., 1422., 1908., 0., 2329., 2817., 0., 2329., 2817., 0.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITH_PART,
                      CWP_COMM_PAR_WITH_PART,
                      true,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WP - WoP - share") {
    res_expected = {1836., 2282., 2282., 2088., 3094., 3094., 0., 0., 0., 1422., 1908., 0., 2329., 2817., 0., 2329., 2817., 0.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITH_PART,
                      CWP_COMM_PAR_WITHOUT_PART,
                      true,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WoP - WP - share") {
    res_expected = {918., 1311., 1311., 1044., 1775., 1775., 0., 0., 0., 2844., 3816., 0., 4088., 4942., 0., 4088., 4942., 0.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITHOUT_PART,
                      CWP_COMM_PAR_WITH_PART,
                      true,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

}

MPI_TEST_CASE("[spatialInterp] - 1p - Location Boxtree #1", 1) {
  CWP_Spatial_interp_t spatial_interp = CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_BOXTREE;

  std::vector<double> res_expected;
  std::vector<double> res;

  SUBCASE("WP - WP - no-shared") {
    res_expected = {918., 1141., 1141., 1044., 1547., 1547., 0., 0., 0., 1422., 1908., 0., 2044., 2471., 0., 2044., 2471., 0.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITH_PART,
                      CWP_COMM_PAR_WITH_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WoP - WoP - no-shared") {
    res_expected = {918., 1141., 1141., 1044., 1547., 1547., 0., 0., 0., 1422., 1908., 0., 2044., 2471., 0., 2044., 2471., 0.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITHOUT_PART,
                      CWP_COMM_PAR_WITHOUT_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WP - WoP - no-shared") {
    res_expected = {918., 1141., 1141., 1044., 1547., 1547., 0., 0., 0., 1422., 1908., 0., 2044., 2471., 0., 2044., 2471., 0.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITH_PART,
                      CWP_COMM_PAR_WITHOUT_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WoP - WP - no-shared") {
    res_expected = {918., 1141., 1141., 1044., 1547., 1547., 0., 0., 0., 1422., 1908., 0., 2044., 2471., 0., 2044., 2471., 0.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITHOUT_PART,
                      CWP_COMM_PAR_WITH_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

}

MPI_TEST_CASE("[spatialInterp] - 2p - Location Boxtree #2", 2) {
  CWP_Spatial_interp_t spatial_interp = CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_BOXTREE;

  std::vector<double> res_expected;
  std::vector<double> res;

  SUBCASE("WP - WP - no-shared") {
    res_expected = {918., 1141., 1141., 1044., 1547., 1547., 0., 0., 0., 1422., 1908., 0., 2044., 2471., 0., 2044., 2471., 0.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITH_PART,
                      CWP_COMM_PAR_WITH_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WP - WoP - no-shared") {
    res_expected = {918., 1141., 1141., 1044., 1547., 1547., 0., 0., 0., 1422., 1908., 0., 2044., 2471., 0., 2044., 2471., 0.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITH_PART,
                      CWP_COMM_PAR_WITHOUT_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WoP - WoP - no-shared") {
    res_expected = {918., 1141., 1141., 1044., 1547., 1547., 0., 0., 0., 1422., 1908., 0., 2044., 2471., 0., 2044., 2471., 0.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITHOUT_PART,
                      CWP_COMM_PAR_WITHOUT_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WoP - WP - no-shared") {
    res_expected = {918., 1141., 1141., 1044., 1547., 1547., 0., 0., 0., 1422., 1908., 0., 2044., 2471., 0., 2044., 2471., 0.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITHOUT_PART,
                      CWP_COMM_PAR_WITH_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WP - WP - share") {
    res_expected = {918., 1311., 1311., 1044., 1775., 1775., 0., 0., 0., 1422., 1908., 0., 2329., 2817., 0., 2329., 2817., 0.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITH_PART,
                      CWP_COMM_PAR_WITH_PART,
                      true,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WP - WoP - share") {
    res_expected = {1836., 2282., 2282., 2088., 3094., 3094., 0., 0., 0., 1422., 1908., 0., 2329., 2817., 0., 2329., 2817., 0.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITH_PART,
                      CWP_COMM_PAR_WITHOUT_PART,
                      true,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WoP - WP - share") {
    res_expected = {918., 1311., 1311., 1044., 1775., 1775., 0., 0., 0., 2844., 3816., 0., 4088., 4942., 0., 4088., 4942., 0.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITHOUT_PART,
                      CWP_COMM_PAR_WITH_PART,
                      true,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

}

MPI_TEST_CASE("[spatialInterp] - 1p - Location All Tgt #1", 1) {
  CWP_Spatial_interp_t spatial_interp = CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_LOCATE_ALL_TGT;

  std::vector<double> res_expected;
  std::vector<double> res;

  SUBCASE("WP - WP - no-shared") {
    res_expected = {1260., 1788., 1788., 1538.4, 2073.8, 2073.8, 0., 0., 0., 2892., 3864.2, 0., 4058., 5319., 0., 4058., 5319., 0.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITH_PART,
                      CWP_COMM_PAR_WITH_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WoP - WoP - no-shared") {
    res_expected = {1260., 1788., 1788., 1538.4, 2073.8, 2073.8, 0., 0., 0., 2892., 3864.2, 0., 4058., 5319., 0., 4058., 5319., 0.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITHOUT_PART,
                      CWP_COMM_PAR_WITHOUT_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WP - WoP - no-shared") {
    res_expected = {1260., 1788., 1788., 1538.4, 2073.8, 2073.8, 0., 0., 0., 2892., 3864.2, 0., 4058., 5319., 0., 4058., 5319., 0.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITH_PART,
                      CWP_COMM_PAR_WITHOUT_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WoP - WP - no-shared") {
    res_expected = {1260., 1788., 1788., 1538.4, 2073.8, 2073.8, 0., 0., 0., 2892., 3864.2, 0., 4058., 5319., 0., 4058., 5319., 0.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITHOUT_PART,
                      CWP_COMM_PAR_WITH_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

}

MPI_TEST_CASE("[spatialInterp] - 2p - Location All Tgt #2", 2) {
  CWP_Spatial_interp_t spatial_interp = CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_LOCATE_ALL_TGT;

  std::vector<double> res_expected;
  std::vector<double> res;

  SUBCASE("WP - WP - no-shared") {
    res_expected = {1260., 1788., 1788., 1538.4, 2073.8, 2073.8, 0., 0., 0., 2892., 3864.2, 0., 4059., 5315.8, 0., 4059., 5315.8, 0.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITH_PART,
                      CWP_COMM_PAR_WITH_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WP - WoP - no-shared") {
    res_expected = {1260., 1788., 1788., 1538.4, 2073.8, 2073.8, 0., 0., 0., 2892., 3864.2, 0., 4059., 5315.8, 0., 4059., 5315.8, 0.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITH_PART,
                      CWP_COMM_PAR_WITHOUT_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WoP - WoP - no-shared") {
    res_expected = {1260., 1788., 1788., 1538.4, 2073.8, 2073.8, 0., 0., 0., 2892., 3864.2, 0., 4059., 5315.8, 0., 4059., 5315.8, 0.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITHOUT_PART,
                      CWP_COMM_PAR_WITHOUT_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WoP - WP - no-shared") {
    res_expected = {1260., 1788., 1788., 1538.4, 2073.8, 2073.8, 0., 0., 0., 2892., 3864.2, 0., 4059., 5315.8, 0., 4059., 5315.8, 0.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITHOUT_PART,
                      CWP_COMM_PAR_WITH_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WP - WP - share") {
    res_expected = {1260., 1979., 1979., 1538.4, 2285., 2285., 0., 0., 0., 2892., 3864.2, 0., 4518., 5924.6, 0., 4518., 5924.6, 0.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITH_PART,
                      CWP_COMM_PAR_WITH_PART,
                      true,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WP - WoP - share") {
    res_expected = {2520., 3576., 3576., 3076.8, 4147.6, 4147.6, 0., 0., 0., 2892., 3864.2, 0., 4518., 5924.6, 0., 4518., 5924.6, 0.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITH_PART,
                      CWP_COMM_PAR_WITHOUT_PART,
                      true,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WoP - WP - share") {
    res_expected = {1260., 1979., 1979., 1538.4, 2285., 2285., 0., 0., 0., 5784., 7728.4, 0., 8118., 10631.6, 0., 8118., 10631.6, 0.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITHOUT_PART,
                      CWP_COMM_PAR_WITH_PART,
                      true,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

}

MPI_TEST_CASE("[spatialInterp] - 1p - User", 1) {
  std::vector<double> res_expected;
  std::vector<double> res;

  SUBCASE("user-interp INTERSECTION") {
    res_expected = {716.2, 0., 0., 0., 0., 0., 0., 0., 0., 1454.2, 0., 0., 0., 0., 0., 0., 0., 0.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITH_PART,
                      CWP_COMM_PAR_WITH_PART,
                      false,
                      CWP_SPATIAL_INTERP_FROM_INTERSECTION,
                      true);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("user-interp NEAREST_TARGETS_LEAST_SQUARES") {
    res_expected = {781.454105, 718.205851, 718.205851, 1324.454442, 1426.332082, 1426.332082, 1324.454442, 1426.332082, 1426.332082, 1543.056302, 2452.406032, 2452.406032, 1773.036565, 2969.581645, 2969.581645, 1773.036565, 2969.581645, 2969.581645};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITH_PART,
                      CWP_COMM_PAR_WITH_PART,
                      false,
                      CWP_SPATIAL_INTERP_FROM_NEAREST_TARGETS_LEAST_SQUARES,
                      true);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("user-interp LOCATION_OCTREE") {
    res_expected = {918., 1141., 1141., 1044., 1547., 1547., 0., 0., 0., 1422., 1908., 0., 2044., 2471., 0., 2044., 2471., 0.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITH_PART,
                      CWP_COMM_PAR_WITH_PART,
                      false,
                      CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE,
                      true);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

}