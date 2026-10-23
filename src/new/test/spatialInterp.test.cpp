#include "doctest/doctest.h"
#include "doctest/extensions/doctest_mpi.h"
#include "cwp.h"

#include <string.h>
#include <unistd.h>
#include <vector>

static double _max_error = 1e-6;
static std::vector<const char *> _forward  = {"code1", "code2"};
static std::vector<const char *> _backward = {"code2", "code1"};

#define CHECK_EPS(array, array_expected, size) \
  for(size_t i = 0; i < size; ++i){ CHECK(array[i] == doctest::Approx(array_expected[i]).epsilon(_max_error)); }

static std::string
_concat
(
  std::string a,
  std::string b
)
{
  return a+b;
}

struct Cases {

  int n_cases = 0;

  std::vector<int>                idx;                // Index to store results
  std::vector<std::string>        fields_interlaced;  // Interlaced  field names
  std::vector<std::string>        fields_interleaved; // Interleaved field names
  std::vector<CWP_Dof_location_t> dof_location1;      // DOF location for code 1
  std::vector<CWP_Dof_location_t> dof_location2;      // DOF location for code 2

  // Fill structure with spatial interpolation method
  // compatible configurations
  Cases
  (
    CWP_Spatial_interp_t  method,
    bool                  backward,
    int                   offset,
    std::string           suffix
  )
  {
    bool cell_cell = false;
    bool cell_vtx  = false;
    bool cell_usr  = false;
    bool vtx_cell  = false;
    bool vtx_vtx   = false;
    bool vtx_usr   = false;
    bool usr_cell  = false;
    bool usr_vtx   = false;
    bool usr_usr   = false;

    switch (method) {
      case CWP_SPATIAL_INTERP_FROM_IDENTITY:
        cell_cell = true;
        vtx_vtx   = true;
        usr_usr   = true;
        break;
      case CWP_SPATIAL_INTERP_FROM_INTERSECTION:
        cell_cell = true;
        break;
      case CWP_SPATIAL_INTERP_FROM_NEAREST_SOURCES_LEAST_SQUARES:
      case CWP_SPATIAL_INTERP_FROM_NEAREST_TARGETS_LEAST_SQUARES:
        cell_cell = true;
        cell_vtx  = true;
        cell_usr  = true;
        vtx_cell  = true;
        vtx_vtx   = true;
        vtx_usr   = true;
        usr_cell  = true;
        usr_vtx   = true;
        usr_usr   = true;
        break;
      case CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE:
      case CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_BOXTREE:
      case CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_LOCATE_ALL_TGT:
        cell_cell = true;
        cell_vtx  = true;
        cell_usr  = !backward;
        vtx_cell  = true;
        vtx_vtx   = true;
        vtx_usr   = !backward;
        usr_cell  = backward;
        usr_vtx   = backward;
        break;
      }

      if (cell_cell) {
        idx               .push_back(offset + 0);
        fields_interlaced .push_back(_concat("fieldinterlacedcellcell",  suffix));
        fields_interleaved.push_back(_concat("fieldinterleavedcellcell", suffix));
        dof_location1     .push_back(CWP_DOF_LOCATION_CELL_CENTER);
        dof_location2     .push_back(CWP_DOF_LOCATION_CELL_CENTER);
        n_cases++;
      }
      if (cell_vtx) {
        idx               .push_back(offset + 1);
        fields_interlaced .push_back(_concat("fieldinterlacedcellvtx",  suffix));
        fields_interleaved.push_back(_concat("fieldinterleavedcellvtx", suffix));
        dof_location1     .push_back(CWP_DOF_LOCATION_CELL_CENTER);
        dof_location2     .push_back(CWP_DOF_LOCATION_NODE);
        n_cases++;
      }
      if (cell_usr) {
        idx               .push_back(offset + 2);
        fields_interlaced .push_back(_concat("fieldinterlacedcellusr",  suffix));
        fields_interleaved.push_back(_concat("fieldinterleavedcellusr", suffix));
        dof_location1     .push_back(CWP_DOF_LOCATION_CELL_CENTER);
        dof_location2     .push_back(CWP_DOF_LOCATION_USER);
        n_cases++;
      }
      if (vtx_cell) {
        idx               .push_back(offset + 3);
        fields_interlaced .push_back(_concat("fieldinterlacedvtxcell",  suffix));
        fields_interleaved.push_back(_concat("fieldinterleavedvtxcell", suffix));
        dof_location1     .push_back(CWP_DOF_LOCATION_NODE);
        dof_location2     .push_back(CWP_DOF_LOCATION_CELL_CENTER);
        n_cases++;
      }
      if (vtx_vtx) {
        idx               .push_back(offset + 4);
        fields_interlaced .push_back(_concat("fieldinterlacedvtxvtx",  suffix));
        fields_interleaved.push_back(_concat("fieldinterleavedvtxvtx", suffix));
        dof_location1     .push_back(CWP_DOF_LOCATION_NODE);
        dof_location2     .push_back(CWP_DOF_LOCATION_NODE);
        n_cases++;
      }
      if (vtx_usr) {
        idx               .push_back(offset + 5);
        fields_interlaced .push_back(_concat("fieldinterlacedvtxusr",  suffix));
        fields_interleaved.push_back(_concat("fieldinterleavedvtxusr", suffix));
        dof_location1     .push_back(CWP_DOF_LOCATION_NODE);
        dof_location2     .push_back(CWP_DOF_LOCATION_USER);
        n_cases++;
      }
      if (usr_cell) {
        idx               .push_back(offset + 6);
        fields_interlaced .push_back(_concat("fieldinterlacedusrcell",  suffix));
        fields_interleaved.push_back(_concat("fieldinterleavedusrcell", suffix));
        dof_location1     .push_back(CWP_DOF_LOCATION_USER);
        dof_location2     .push_back(CWP_DOF_LOCATION_CELL_CENTER);
        n_cases++;
      }
      if (usr_vtx) {
        idx               .push_back(offset + 7);
        fields_interlaced .push_back(_concat("fieldinterlacedusrvtx",  suffix));
        fields_interleaved.push_back(_concat("fieldinterleavedusrvtx", suffix));
        dof_location1     .push_back(CWP_DOF_LOCATION_USER);
        dof_location2     .push_back(CWP_DOF_LOCATION_NODE);
        n_cases++;
      }
      if (usr_usr) {
        idx               .push_back(offset + 8);
        fields_interlaced .push_back(_concat("fieldinterlacedusrusr",  suffix));
        fields_interleaved.push_back(_concat("fieldinterleavedusrusr", suffix));
        dof_location1     .push_back(CWP_DOF_LOCATION_USER);
        dof_location2     .push_back(CWP_DOF_LOCATION_USER);
        n_cases++;
      }
  }

};

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

  for (int i_part = 0; i_part < n_part; i_part++) {
    free(pface_edge    [i_part]);
    free(pedge_vtx     [i_part]);
    free(pedge_ln_to_gn[i_part]);
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
        int                  *n_part,         // Number of partitions
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
      for (int i_part = 0; i_part < n_part[i_code]; i_part++) {
        if (loc_send == CWP_DOF_LOCATION_CELL_CENTER) {
          for (int i = 0; i < pn_face[i_code][i_part]; i++) {
            if (interleaved) {
              send_val[i_code][i_part][i                          ] = static_cast<double>(     pface_ln_to_gn[i_code][i_part][i]);
              send_val[i_code][i_part][i + pn_face[i_code][i_part]] = static_cast<double>(-2 * pface_ln_to_gn[i_code][i_part][i]);
            }
            else {
              send_val[i_code][i_part][2 * i    ] = static_cast<double>(     pface_ln_to_gn[i_code][i_part][i]);
              send_val[i_code][i_part][2 * i + 1] = static_cast<double>(-2 * pface_ln_to_gn[i_code][i_part][i]);
            }
          }
        }
        else {
          for (int i = 0; i < pn_vtx[i_code][i_part]; i++) {
            if (interleaved) {
              send_val[i_code][i_part][i                         ] = static_cast<double>(     pvtx_ln_to_gn[i_code][i_part][i]);
              send_val[i_code][i_part][i + pn_vtx[i_code][i_part]] = static_cast<double>(-2 * pvtx_ln_to_gn[i_code][i_part][i]);
            }
            else {
              send_val[i_code][i_part][2 * i    ] = static_cast<double>(     pvtx_ln_to_gn[i_code][i_part][i]);
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
      for (int i_part = 0; i_part < n_part[i_code]; i_part++) {
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
        int                  *n_part,        // Number of partitions
        int                 **pn_vtx,        // Number of vertices
        int                 **pn_face,       // Number of faces
        CWP_g_num_t        ***pvtx_ln_to_gn, // Vertex global numbering
        CWP_g_num_t        ***pface_ln_to_gn // Face global numbering
)
{
  // Allocate buffers and set field data
  double ***send_val = static_cast<double ***>(malloc(sizeof(double **) * n_code));
  double ***recv_val = static_cast<double ***>(malloc(sizeof(double **) * n_code));
  for (int i_code = 0; i_code < n_code; i_code++) {
    if (strcmp(code_names[i_code], order[0]) == 0) {
      send_val[i_code] = static_cast<double **>(malloc(sizeof(double *) * n_part[i_code]));
      for (int i_part = 0; i_part < n_part[i_code]; i_part++) {
        if (loc_send == CWP_DOF_LOCATION_CELL_CENTER) {
          send_val[i_code][i_part] = static_cast<double *>(malloc(sizeof(double) * 2 * pn_face[i_code][i_part]));
        }
        else {
          send_val[i_code][i_part] = static_cast<double *>(malloc(sizeof(double) * 2 * pn_vtx[i_code][i_part]));
        }

        CWP_Field_data_set(code_names[i_code],
                           cpl_name,
                           field_name1,
                           i_part,
                           CWP_FIELD_MAP_SOURCE,
                           send_val  [i_code][i_part]);
        CWP_Field_data_set(code_names[i_code],
                           cpl_name,
                           field_name2,
                           i_part,
                           CWP_FIELD_MAP_SOURCE,
                           send_val  [i_code][i_part]);
      }
    }

    if (strcmp(code_names[i_code], order[1]) == 0) {
      recv_val[i_code] = static_cast<double **>(malloc(sizeof(double *) * n_part[i_code]));
      for (int i_part = 0; i_part < n_part[i_code]; i_part++) {
        if (loc_recv == CWP_DOF_LOCATION_CELL_CENTER) {
          recv_val[i_code][i_part] = static_cast<double *>(malloc(sizeof(double) * 2 * pn_face[i_code][i_part]));
        }
        else {
          recv_val[i_code][i_part] = static_cast<double *>(malloc(sizeof(double) * 2 * pn_vtx[i_code][i_part]));
        }

        CWP_Field_data_set(code_names[i_code],
                           cpl_name,
                           field_name1,
                           i_part,
                           CWP_FIELD_MAP_TARGET,
                           recv_val  [i_code][i_part]);
        CWP_Field_data_set(code_names[i_code],
                           cpl_name,
                           field_name2,
                           i_part,
                           CWP_FIELD_MAP_TARGET,
                           recv_val  [i_code][i_part]);
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
      for (int i_part = 0; i_part < n_part[i_code]; i_part++) {
        free(send_val[i_code][i_part]);
      }
      free(send_val[i_code]);
    }
    if (strcmp(code_names[i_code], order[1]) == 0) {
      for (int i_part = 0; i_part < n_part[i_code]; i_part++) {
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
  // Generate 9 cases for send -> recv from code1 to code2
  // and      9 cases for send -> recv from code2 to code1
  // 9 cases to cover all source/target points location
  Cases forward_cases (spatial_interp, false, 0, "_forward");
  Cases backward_cases(spatial_interp, true,  9, "_backward");

  std::vector<double> res_loc = {0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                 0., 0., 0., 0., 0., 0., 0., 0., 0.};
  std::vector<double> res(res_loc.size());

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
  // 2 ranks = dispatch codes to share ranks if share_ranks == true
  else {
    if (share_ranks) {
      code_names = {"code1", "code2"};
    }
    else {
      if (i_rank == 1) {
        code_names.push_back("code1");
      }
      else { // i_rank == 0
        code_names.push_back("code2");
      }
    }
  }
  int n_code = static_cast<int>(code_names.size());

  MPI_Comm *intra_comm = static_cast<MPI_Comm *>(malloc(sizeof(MPI_Comm) * n_code));
  CWP_Init(comm,
           n_code,
           code_names.data(),
           is_active_rank,
           intra_comm);

  // Create coupling
  int        *n_part    = static_cast<int *>(malloc(sizeof(int) * n_code));
  int         n_vtx_seg = 9;
  const char *cpl_name  = "spatialInterpTest";

  for (int i_code = 0; i_code < n_code; i_code++) {

    const char *coupled_code_name;
    CWP_Comm_t  part;

    if (strcmp(code_names[i_code], "code1") == 0) {
      coupled_code_name = "code2";
      part           = part1;
      n_part[i_code] = 1;
    }

    if (strcmp(code_names[i_code], "code2") == 0) {
      coupled_code_name = "code1";
      part           = part2;
      n_part[i_code] = 2;
    }

    CWP_Cpl_create(code_names[i_code],
                   cpl_name,
                   coupled_code_name,
                   CWP_INTERFACE_SURFACE,
                   part,
                   spatial_interp,
                   n_part    [i_code],
                   CWP_DYNAMIC_MESH_DEFORMABLE,
                   CWP_TIME_EXCH_USER_CONTROLLED);

  }

  // Enable visualization
  for (int i_code = 0; i_code < n_code; i_code++) {
    CWP_Visu_set(code_names[i_code],
                 cpl_name,
                 2,
                 CWP_VISU_FORMAT_ENSIGHT,
                 "binary");
  }

  // Define mesh
  int          **pn_face        = static_cast<int          **>(malloc(sizeof(int          *) * n_code));
  int          **pn_vtx         = static_cast<int          **>(malloc(sizeof(int          *) * n_code));
  int         ***pface_vtx_idx  = static_cast<int         ***>(malloc(sizeof(int         **) * n_code));
  int         ***pface_vtx      = static_cast<int         ***>(malloc(sizeof(int         **) * n_code));
  double      ***pvtx_coord     = static_cast<double      ***>(malloc(sizeof(double      **) * n_code));
  CWP_g_num_t ***pface_ln_to_gn = static_cast<CWP_g_num_t ***>(malloc(sizeof(CWP_g_num_t **) * n_code));
  CWP_g_num_t ***pvtx_ln_to_gn  = static_cast<CWP_g_num_t ***>(malloc(sizeof(CWP_g_num_t **) * n_code));

  for (int i_code = 0; i_code < n_code; i_code++) {

    // Generate mesh on intra_comm communicator or local rank communicator
    // if no partitioning is required
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
              n_part         [i_code],
              CWPT_SPLIT_DUAL_WITH_HILBERT,
              n_vtx_seg,
              &pn_face       [i_code],
              &pn_vtx        [i_code],
              &pface_vtx_idx [i_code],
              &pface_vtx     [i_code],
              &pvtx_coord    [i_code],
              &pface_ln_to_gn[i_code],
              &pvtx_ln_to_gn [i_code]);
    MPI_Comm_free(&mesh_comm);

    // Move meshes to have unlocated points
    if (strcmp(code_names[i_code], "code1") == 0) {
      for (int i_part = 0; i_part < n_part[i_code]; i_part++) {
        for (int i_vtx = 0; i_vtx < pn_vtx[i_code][i_part]; i_vtx++) {
          pvtx_coord[i_code][i_part][3*i_vtx+1] += 0.3;
        }
      }
    }

    if (strcmp(code_names[i_code], "code2") == 0) {
      for (int i_part = 0; i_part < n_part[i_code]; i_part++) {
        for (int i_vtx = 0; i_vtx < pn_vtx[i_code][i_part]; i_vtx++) {
          pvtx_coord[i_code][i_part][3*i_vtx] += 0.2;
        }
      }
    }

    // Set mesh
    for (int i_part = 0; i_part < n_part[i_code]; i_part++) {
      CWP_Mesh_interf_vtx_set(code_names   [i_code],
                              cpl_name,
                              i_part,
                              pn_vtx       [i_code][i_part],
                              pvtx_coord   [i_code][i_part],
                              pvtx_ln_to_gn[i_code][i_part]);

      CWP_Mesh_interf_from_facevtx_set(code_names    [i_code],
                                       cpl_name,
                                       i_part,
                                       pn_face       [i_code][i_part],
                                       pface_vtx_idx [i_code][i_part],
                                       pface_vtx     [i_code][i_part],
                                       pface_ln_to_gn[i_code][i_part]);

      CWP_User_tgt_pts_set(code_names   [i_code],
                           cpl_name,
                           i_part,
                           pn_vtx       [i_code][i_part],
                           pvtx_coord   [i_code][i_part],
                           pvtx_ln_to_gn[i_code][i_part]);
    }

    CWP_Mesh_interf_finalize(code_names[i_code], cpl_name);
  }

  // Define fields for each case :
  //  - interlaced and interleaved exchange and interpolation have
  //    to return the same values
  //  - user-defined interpolation functions are set if required
  //  - broadcast is enabled if no partitioning is required
  for (int i_code = 0; i_code < n_code; i_code++) {
    if (strcmp(code_names[i_code], "code1") == 0) {
      for (int i_field = 0; i_field < forward_cases.n_cases; i_field++) {
        CWP_Field_create(code_names[i_code],
                         cpl_name,
                         forward_cases.fields_interlaced[i_field].c_str(),
                         CWP_DOUBLE,
                         CWP_FIELD_STORAGE_INTERLACED,
                         2,
                         forward_cases.dof_location1[i_field],
                         CWP_FIELD_EXCH_SEND,
                         (i_field == 0 ? CWP_STATUS_ON:CWP_STATUS_OFF));
        CWP_Field_create(code_names[i_code],
                         cpl_name,
                         forward_cases.fields_interleaved[i_field].c_str(),
                         CWP_DOUBLE,
                         CWP_FIELD_STORAGE_INTERLEAVED,
                         2,
                         forward_cases.dof_location1[i_field],
                         CWP_FIELD_EXCH_SEND,
                         (i_field == 0 ? CWP_STATUS_ON:CWP_STATUS_OFF));
        if (part1 == CWP_COMM_PAR_WITHOUT_PART) {
          CWP_Involved_srcs_bcast_enable(code_names[i_code],
                                         cpl_name,
                                         forward_cases.fields_interlaced[i_field].c_str());
          CWP_Involved_srcs_bcast_enable(code_names[i_code],
                                         cpl_name,
                                         forward_cases.fields_interleaved[i_field].c_str());
        }
      }
      for (int i_field = 0; i_field < backward_cases.n_cases; i_field++) {
        CWP_Field_create(code_names[i_code],
                         cpl_name,
                         backward_cases.fields_interlaced[i_field].c_str(),
                         CWP_DOUBLE,
                         CWP_FIELD_STORAGE_INTERLACED,
                         2,
                         backward_cases.dof_location1[i_field],
                         CWP_FIELD_EXCH_RECV,
                         (i_field == 0 ? CWP_STATUS_ON:CWP_STATUS_OFF));
        CWP_Field_create(code_names[i_code],
                         cpl_name,
                         backward_cases.fields_interleaved[i_field].c_str(),
                         CWP_DOUBLE,
                         CWP_FIELD_STORAGE_INTERLEAVED,
                         2,
                         backward_cases.dof_location1[i_field],
                         CWP_FIELD_EXCH_RECV,
                         (i_field == 0 ? CWP_STATUS_ON:CWP_STATUS_OFF));
        if (user_interp) {
          CWP_Field_interp_function_set(code_names[i_code],
                                        cpl_name,
                                        backward_cases.fields_interlaced[i_field].c_str(),
                                        user_interpolation);
          CWP_Field_interp_function_set(code_names[i_code],
                                        cpl_name,
                                        backward_cases.fields_interleaved[i_field].c_str(),
                                        user_interpolation);
        }
        if (part1 == CWP_COMM_PAR_WITHOUT_PART) {
          CWP_Computed_tgts_bcast_enable(code_names[i_code],
                                         cpl_name,
                                         backward_cases.fields_interlaced[i_field].c_str());
          CWP_Computed_tgts_bcast_enable(code_names[i_code],
                                         cpl_name,
                                         backward_cases.fields_interleaved[i_field].c_str());
        }
      }
    }

    if (strcmp(code_names[i_code], "code2") == 0) {
      for (int i_field = 0; i_field < forward_cases.n_cases; i_field++) {
        CWP_Field_create(code_names[i_code],
                         cpl_name,
                         forward_cases.fields_interlaced[i_field].c_str(),
                         CWP_DOUBLE,
                         CWP_FIELD_STORAGE_INTERLACED,
                         2,
                         forward_cases.dof_location2[i_field],
                         CWP_FIELD_EXCH_RECV,
                         (i_field == 0 ? CWP_STATUS_ON:CWP_STATUS_OFF));
        CWP_Field_create(code_names[i_code],
                         cpl_name,
                         forward_cases.fields_interleaved[i_field].c_str(),
                         CWP_DOUBLE,
                         CWP_FIELD_STORAGE_INTERLEAVED,
                         2,
                         forward_cases.dof_location2[i_field],
                         CWP_FIELD_EXCH_RECV,
                         (i_field == 0 ? CWP_STATUS_ON:CWP_STATUS_OFF));
        if (user_interp) {
          CWP_Field_interp_function_set(code_names[i_code],
                                        cpl_name,
                                        forward_cases.fields_interlaced[i_field].c_str(),
                                        user_interpolation);
          CWP_Field_interp_function_set(code_names[i_code],
                                        cpl_name,
                                        forward_cases.fields_interleaved[i_field].c_str(),
                                        user_interpolation);
        }
        if (part2 == CWP_COMM_PAR_WITHOUT_PART) {
          CWP_Computed_tgts_bcast_enable(code_names[i_code],
                                         cpl_name,
                                         forward_cases.fields_interlaced[i_field].c_str());
          CWP_Computed_tgts_bcast_enable(code_names[i_code],
                                         cpl_name,
                                         forward_cases.fields_interleaved[i_field].c_str());
        }
      }
      for (int i_field = 0; i_field < backward_cases.n_cases; i_field++) {
        CWP_Field_create(code_names[i_code],
                         cpl_name,
                         backward_cases.fields_interlaced[i_field].c_str(),
                         CWP_DOUBLE,
                         CWP_FIELD_STORAGE_INTERLACED,
                         2,
                         backward_cases.dof_location2[i_field],
                         CWP_FIELD_EXCH_SEND,
                         (i_field == 0 ? CWP_STATUS_ON:CWP_STATUS_OFF));
        CWP_Field_create(code_names[i_code],
                         cpl_name,
                         backward_cases.fields_interleaved[i_field].c_str(),
                         CWP_DOUBLE,
                         CWP_FIELD_STORAGE_INTERLEAVED,
                         2,
                         backward_cases.dof_location2[i_field],
                         CWP_FIELD_EXCH_SEND,
                         (i_field == 0 ? CWP_STATUS_ON:CWP_STATUS_OFF));
        if (part2 == CWP_COMM_PAR_WITHOUT_PART) {
          CWP_Involved_srcs_bcast_enable(code_names[i_code],
                                         cpl_name,
                                         backward_cases.fields_interlaced[i_field].c_str());
          CWP_Involved_srcs_bcast_enable(code_names[i_code],
                                         cpl_name,
                                         backward_cases.fields_interleaved[i_field].c_str());
        }
      }
    }

  }

  // Loop on time to check the update of interpolation weights
  for (int i_step = 0; i_step < 2; i_step++) {

    // Begin time step
    for (int i_code = 0; i_code < n_code; i_code++) {
      CWP_Time_step_beg(code_names[i_code], static_cast<double>(i_step));
    }

    // Compute interpolation weights
    for (int i_code = 0; i_code < n_code; i_code++) {
      CWP_Spatial_interp_weights_compute(code_names[i_code], cpl_name);
    }

    // Get value from forward cases
    for (int i_case = 0; i_case < forward_cases.n_cases; i_case++) {
      res_loc[forward_cases.idx[i_case]] = _get_value_cpl(_forward.data(),
                                                          n_code,
                                                          code_names.data(),
                                                          cpl_name,
                                                          forward_cases.fields_interlaced [i_case].c_str(),
                                                          forward_cases.fields_interleaved[i_case].c_str(),
                                                          forward_cases.dof_location1     [i_case],
                                                          forward_cases.dof_location2     [i_case],
                                                          n_part,
                                                          pn_vtx,
                                                          pn_face,
                                                          pvtx_ln_to_gn,
                                                          pface_ln_to_gn);
    }

    // Get value from backward cases
    for (int i_case = 0; i_case < backward_cases.n_cases; i_case++) {
      res_loc[backward_cases.idx[i_case]] = _get_value_cpl(_backward.data(),
                                                           n_code,
                                                           code_names.data(),
                                                           cpl_name,
                                                           backward_cases.fields_interlaced [i_case].c_str(),
                                                           backward_cases.fields_interleaved[i_case].c_str(),
                                                           backward_cases.dof_location2     [i_case],
                                                           backward_cases.dof_location1     [i_case],
                                                           n_part,
                                                           pn_vtx,
                                                           pn_face,
                                                           pvtx_ln_to_gn,
                                                           pface_ln_to_gn);
    }

    // Sum values across ranks
    MPI_Allreduce(res_loc.data(), res.data(), res.size(), MPI_DOUBLE, MPI_SUM, comm);

    // End time step
    for (int i_code = 0; i_code < n_code; i_code++) {
      CWP_Time_step_end(code_names[i_code]);
    }

  }

  // Check files exist
  MPI_Barrier(comm);
  if (i_rank == 0) {
    const std::vector<std::string> names {
      "code1_code2/CHR.case",
      "code1_code2/chr.geo.00001",
      "code1_code2/chr.partitioning.00001",
      "code1_code2/chr.ranking.00001",
      "code1_code2/chr.r_fieldinterlacedcellcell_backward1.00001",
      "code1_code2/chr.r_fieldinterlacedcellcell_backward2.00001",
      "code1_code2/chr.r_fieldinterlacedcellcell_backward_status.00001",
      "code1_code2/chr.r_fieldinterleavedcellcell_backward1.00001",
      "code1_code2/chr.r_fieldinterleavedcellcell_backward2.00001",
      "code1_code2/chr.r_fieldinterleavedcellcell_backward_status.00001",
      "code1_code2/chr.s_fieldinterlacedcellcell_forward1.00001",
      "code1_code2/chr.s_fieldinterlacedcellcell_forward2.00001",
      "code1_code2/chr.s_fieldinterlacedcellcell_forward_status.00001",
      "code1_code2/chr.s_fieldinterleavedcellcell_forward1.00001",
      "code1_code2/chr.s_fieldinterleavedcellcell_forward2.00001",
      "code1_code2/chr.s_fieldinterleavedcellcell_forward_status.00001",
      "code1_code2/chr.user_tgt_partitioning.00001",
      "code1_code2/chr.user_tgt_ranking.00001",
      "code2_code1/CHR.case",
      "code2_code1/chr.geo.00001",
      "code2_code1/chr.partitioning.00001",
      "code2_code1/chr.ranking.00001",
      "code2_code1/chr.r_fieldinterlacedcellcell_forward1.00001",
      "code2_code1/chr.r_fieldinterlacedcellcell_forward2.00001",
      "code2_code1/chr.r_fieldinterlacedcellcell_forward_status.00001",
      "code2_code1/chr.r_fieldinterleavedcellcell_forward1.00001",
      "code2_code1/chr.r_fieldinterleavedcellcell_forward2.00001",
      "code2_code1/chr.r_fieldinterleavedcellcell_forward_status.00001",
      "code2_code1/chr.s_fieldinterlacedcellcell_backward1.00001",
      "code2_code1/chr.s_fieldinterlacedcellcell_backward2.00001",
      "code2_code1/chr.s_fieldinterlacedcellcell_backward_status.00001",
      "code2_code1/chr.s_fieldinterleavedcellcell_backward1.00001",
      "code2_code1/chr.s_fieldinterleavedcellcell_backward2.00001",
      "code2_code1/chr.s_fieldinterleavedcellcell_backward_status.00001",
      "code2_code1/chr.user_tgt_partitioning.00001",
      "code2_code1/chr.user_tgt_ranking.00001"
    };
    const std::string prefix = "cwipi_writer/spatialInterpTest_";
    for (const auto &n : names) {
      const std::string concat = prefix + n;
      CHECK(access(concat.c_str(), F_OK) == 0);
      remove(concat.c_str());
    }
  }

  // Free all
  for (int i_code = 0; i_code < n_code; i_code++) {
    if (strcmp(code_names[i_code], "code1") == 0) {
      for (int i_field = 0; i_field < forward_cases.n_cases; i_field++) {
        CWP_Field_del(code_names[i_code],
                      cpl_name,
                      forward_cases.fields_interlaced[i_field].c_str());
        CWP_Field_del(code_names[i_code],
                      cpl_name,
                      forward_cases.fields_interleaved[i_field].c_str());
      }
      for (int i_field = 0; i_field < backward_cases.n_cases; i_field++) {
        CWP_Field_del(code_names[i_code],
                      cpl_name,
                      backward_cases.fields_interlaced[i_field].c_str());
        CWP_Field_del(code_names[i_code],
                      cpl_name,
                      backward_cases.fields_interleaved[i_field].c_str());
      }
    }

    if (strcmp(code_names[i_code], "code2") == 0) {
      for (int i_field = 0; i_field < forward_cases.n_cases; i_field++) {
        CWP_Field_del(code_names[i_code],
                      cpl_name,
                      forward_cases.fields_interlaced[i_field].c_str());
        CWP_Field_del(code_names[i_code],
                      cpl_name,
                      forward_cases.fields_interleaved[i_field].c_str());
      }
      for (int i_field = 0; i_field < backward_cases.n_cases; i_field++) {
        CWP_Field_del(code_names[i_code],
                      cpl_name,
                      backward_cases.fields_interlaced[i_field].c_str());
        CWP_Field_del(code_names[i_code],
                      cpl_name,
                      backward_cases.fields_interleaved[i_field].c_str());
      }
    }

    CWP_Mesh_interf_del(code_names[i_code], cpl_name);
    CWP_Cpl_del(code_names[i_code], cpl_name);

  }

  for (int i_code = 0; i_code < n_code; i_code++) {
    for (int i_part = 0; i_part < n_part[i_code]; i_part++) {
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
  free(n_part);

  free(intra_comm);

  CWP_Finalize();

  return res;
}

MPI_TEST_CASE("[spatialInterp] - 1p - Identity #1", 1) {
  CWP_Spatial_interp_t spatial_interp = CWP_SPATIAL_INTERP_FROM_IDENTITY;

  std::vector<double> res_expected;
  std::vector<double> res;

  SUBCASE("WP - WP - no-shared") {
    res_expected = {2080., 0., 0., 0., 3690., 0., 0., 0., 3690., 2080., 0., 0., 0., 3321., 0., 0., 0., 3321.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITH_PART,
                      CWP_COMM_PAR_WITH_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WoP - WoP - no-shared") {
    res_expected = {2080., 0., 0., 0., 3690., 0., 0., 0., 3690., 2080., 0., 0., 0., 3321., 0., 0., 0., 3321.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITHOUT_PART,
                      CWP_COMM_PAR_WITHOUT_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WP - WoP - no-shared") {
    res_expected = {2080., 0., 0., 0., 3690., 0., 0., 0., 3690., 2080., 0., 0., 0., 3321., 0., 0., 0., 3321.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITH_PART,
                      CWP_COMM_PAR_WITHOUT_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WoP - WP - no-shared") {
    res_expected = {2080., 0., 0., 0., 3690., 0., 0., 0., 3690., 2080., 0., 0., 0., 3321., 0., 0., 0., 3321.};
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
    res_expected = {2080., 0., 0., 0., 3690., 0., 0., 0., 3690., 2080., 0., 0., 0., 3321., 0., 0., 0., 3321.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITH_PART,
                      CWP_COMM_PAR_WITH_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WP - WoP - no-shared") {
    res_expected = {2080., 0., 0., 0., 3690., 0., 0., 0., 3690., 2080., 0., 0., 0., 3321., 0., 0., 0., 3321.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITH_PART,
                      CWP_COMM_PAR_WITHOUT_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WoP - WoP - no-shared") {
    res_expected = {2080., 0., 0., 0., 3690., 0., 0., 0., 3690., 2080., 0., 0., 0., 3321., 0., 0., 0., 3321.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITHOUT_PART,
                      CWP_COMM_PAR_WITHOUT_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WoP - WP - no-shared") {
    res_expected = {2080., 0., 0., 0., 3690., 0., 0., 0., 3690., 2080., 0., 0., 0., 3321., 0., 0., 0., 3321.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITHOUT_PART,
                      CWP_COMM_PAR_WITH_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WP - WP - shared") {
    res_expected = {2080., 0., 0., 0., 4100., 0., 0., 0., 4100., 2080., 0., 0., 0., 3690., 0., 0., 0., 3690.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITH_PART,
                      CWP_COMM_PAR_WITH_PART,
                      true,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WP - WoP - shared") {
    res_expected = {4160., 0., 0., 0., 7380., 0., 0., 0., 7380., 2080., 0., 0., 0., 3690., 0., 0., 0., 3690.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITH_PART,
                      CWP_COMM_PAR_WITHOUT_PART,
                      true,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WoP - WP - shared") {
    res_expected = {2080., 0., 0., 0., 4100., 0., 0., 0., 4100., 4160., 0., 0., 0., 6642., 0., 0., 0., 6642.};
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
    res_expected = {716.2, 0., 0., 0., 0., 0., 0., 0., 0., 1411., 0., 0., 0., 0., 0., 0., 0., 0.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITH_PART,
                      CWP_COMM_PAR_WITH_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WoP - WoP - no-shared") {
    res_expected = {716.2, 0., 0., 0., 0., 0., 0., 0., 0., 1411., 0., 0., 0., 0., 0., 0., 0., 0.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITHOUT_PART,
                      CWP_COMM_PAR_WITHOUT_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WP - WoP - no-shared") {
    res_expected = {716.2, 0., 0., 0., 0., 0., 0., 0., 0., 1411., 0., 0., 0., 0., 0., 0., 0., 0.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITH_PART,
                      CWP_COMM_PAR_WITHOUT_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WoP - WP - no-shared") {
    res_expected = {716.2, 0., 0., 0., 0., 0., 0., 0., 0., 1411., 0., 0., 0., 0., 0., 0., 0., 0.};
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
    res_expected = {716.2, 0., 0., 0., 0., 0., 0., 0., 0., 1411., 0., 0., 0., 0., 0., 0., 0., 0.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITH_PART,
                      CWP_COMM_PAR_WITH_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WP - WoP - no-shared") {
    res_expected = {716.2, 0., 0., 0., 0., 0., 0., 0., 0., 1411., 0., 0., 0., 0., 0., 0., 0., 0.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITH_PART,
                      CWP_COMM_PAR_WITHOUT_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WoP - WoP - no-shared") {
    res_expected = {716.2, 0., 0., 0., 0., 0., 0., 0., 0., 1411., 0., 0., 0., 0., 0., 0., 0., 0.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITHOUT_PART,
                      CWP_COMM_PAR_WITHOUT_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WoP - WP - no-shared") {
    res_expected = {716.2, 0., 0., 0., 0., 0., 0., 0., 0., 1411., 0., 0., 0., 0., 0., 0., 0., 0.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITHOUT_PART,
                      CWP_COMM_PAR_WITH_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WP - WP - shared") {
    res_expected = {716.2, 0., 0., 0., 0., 0., 0., 0., 0., 1454.2, 0., 0., 0., 0., 0., 0., 0., 0.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITH_PART,
                      CWP_COMM_PAR_WITH_PART,
                      true,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WP - WoP - shared") {
    res_expected = {1432.4, 0., 0., 0., 0., 0., 0., 0., 0., 1454.2, 0., 0., 0., 0., 0., 0., 0., 0.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITH_PART,
                      CWP_COMM_PAR_WITHOUT_PART,
                      true,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WoP - WP - shared") {
    res_expected = {716.2, 0., 0., 0., 0., 0., 0., 0., 0., 2822., 0., 0., 0., 0., 0., 0., 0., 0.};
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
    res_expected = {1091.571319, 1542.465595, 1542.465595, 1444.805946, 2073.302578, 2073.302578, 1444.805946, 2073.302578, 2073.302578, 3021.93645, 3776.923976, 3776.923976, 3766.380448, 4684.184569, 4684.184569, 3766.380448, 4684.184569, 4684.184569};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITH_PART,
                      CWP_COMM_PAR_WITH_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WoP - WoP - no-shared") {
    res_expected = {1091.571319, 1542.465595, 1542.465595, 1444.805946, 2073.302578, 2073.302578, 1444.805946, 2073.302578, 2073.302578, 3021.93645, 3776.923976, 3776.923976, 3766.380448, 4684.184569, 4684.184569, 3766.380448, 4684.184569, 4684.184569};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITHOUT_PART,
                      CWP_COMM_PAR_WITHOUT_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WP - WoP - no-shared") {
    res_expected = {1091.571319, 1542.465595, 1542.465595, 1444.805946, 2073.302578, 2073.302578, 1444.805946, 2073.302578, 2073.302578, 3021.93645, 3776.923976, 3776.923976, 3766.380448, 4684.184569, 4684.184569, 3766.380448, 4684.184569, 4684.184569};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITH_PART,
                      CWP_COMM_PAR_WITHOUT_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WoP - WP - no-shared") {
    res_expected = {1091.571319, 1542.465595, 1542.465595, 1444.805946, 2073.302578, 2073.302578, 1444.805946, 2073.302578, 2073.302578, 3021.93645, 3776.923976, 3776.923976, 3766.380448, 4684.184569, 4684.184569, 3766.380448, 4684.184569, 4684.184569};
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
    res_expected = {1091.571319, 1542.465595, 1542.465595, 1444.805946, 2073.302578, 2073.302578, 1444.805946, 2073.302578, 2073.302578, 3021.93645, 3776.923976, 3776.923976, 3766.380448, 4684.184569, 4684.184569, 3766.380448, 4684.184569, 4684.184569};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITH_PART,
                      CWP_COMM_PAR_WITH_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WP - WoP - no-shared") {
    res_expected = {1091.571319, 1542.465595, 1542.465595, 1444.805946, 2073.302578, 2073.302578, 1444.805946, 2073.302578, 2073.302578, 3021.93645, 3776.923976, 3776.923976, 3766.380448, 4684.184569, 4684.184569, 3766.380448, 4684.184569, 4684.184569};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITH_PART,
                      CWP_COMM_PAR_WITHOUT_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WoP - WoP - no-shared") {
    res_expected = {1091.571319, 1542.465595, 1542.465595, 1444.805946, 2073.302578, 2073.302578, 1444.805946, 2073.302578, 2073.302578, 3021.93645, 3776.923976, 3776.923976, 3766.380448, 4684.184569, 4684.184569, 3766.380448, 4684.184569, 4684.184569};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITHOUT_PART,
                      CWP_COMM_PAR_WITHOUT_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WoP - WP - no-shared") {
    res_expected = {1091.571319, 1542.465595, 1542.465595, 1444.805946, 2073.302578, 2073.302578, 1444.805946, 2073.302578, 2073.302578, 3021.93645, 3776.923976, 3776.923976, 3766.380448, 4684.184569, 4684.184569, 3766.380448, 4684.184569, 4684.184569};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITHOUT_PART,
                      CWP_COMM_PAR_WITH_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WP - WP - shared") {
    res_expected = {1091.571319, 1728.831409, 1728.831409, 1444.805946, 2316.761505, 2316.761505, 1444.805946, 2316.761505, 2316.761505, 3021.93645, 3776.923976, 3776.923976, 4216.273271, 5232.530723, 5232.530723, 4216.273271, 5232.530723, 5232.530723};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITH_PART,
                      CWP_COMM_PAR_WITH_PART,
                      true,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WP - WoP - shared") {
    res_expected = {2183.142638, 3084.93119, 3084.93119, 2889.611892, 4146.605155, 4146.605155, 2889.611892, 4146.605155, 4146.605155, 3021.93645, 3776.923976, 3776.923976, 4216.273271, 5232.530723, 5232.530723, 4216.273271, 5232.530723, 5232.530723};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITH_PART,
                      CWP_COMM_PAR_WITHOUT_PART,
                      true,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WoP - WP - shared") {
    res_expected = {1091.571319, 1728.831409, 1728.831409, 1444.805946, 2316.761505, 2316.761505, 1444.805946, 2316.761505, 2316.761505, 6043.872901, 7553.847951, 7553.847951, 7532.760897, 9368.369138, 9368.369138, 7532.760897, 9368.369138, 9368.369138};
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
    res_expected = {760.796303, 733.219025, 733.219025, 1326.725839, 1397.168421, 1397.168421, 1326.725839, 1397.168421, 1397.168421, 1532.96506, 2386.106979, 2386.106979, 1510.730113, 2535.464648, 2535.464648, 1510.730113, 2535.464648, 2535.464648};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITH_PART,
                      CWP_COMM_PAR_WITH_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WoP - WoP - no-shared") {
    res_expected = {760.796303, 733.219025, 733.219025, 1326.725839, 1397.168421, 1397.168421, 1326.725839, 1397.168421, 1397.168421, 1532.96506, 2386.106979, 2386.106979, 1510.730113, 2535.464648, 2535.464648, 1510.730113, 2535.464648, 2535.464648};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITHOUT_PART,
                      CWP_COMM_PAR_WITHOUT_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WP - WoP - no-shared") {
    res_expected = {760.796303, 733.219025, 733.219025, 1326.725839, 1397.168421, 1397.168421, 1326.725839, 1397.168421, 1397.168421, 1532.96506, 2386.106979, 2386.106979, 1510.730113, 2535.464648, 2535.464648, 1510.730113, 2535.464648, 2535.464648};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITH_PART,
                      CWP_COMM_PAR_WITHOUT_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WoP - WP - no-shared") {
    res_expected = {760.796303, 733.219025, 733.219025, 1326.725839, 1397.168421, 1397.168421, 1326.725839, 1397.168421, 1397.168421, 1532.96506, 2386.106979, 2386.106979, 1510.730113, 2535.464648, 2535.464648, 1510.730113, 2535.464648, 2535.464648};
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
    res_expected = {760.796303, 733.909723, 733.909723, 1325.958397, 1397.168421, 1397.168421, 1325.958397, 1397.168421, 1397.168421, 1532.96506, 2386.106979, 2386.106979, 1508.576167, 2535.464648, 2535.464648, 1508.576167, 2535.464648, 2535.464648};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITH_PART,
                      CWP_COMM_PAR_WITH_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WP - WoP - no-shared") {
    res_expected = {760.796303, 733.909723, 733.909723, 1325.958397, 1397.168421, 1397.168421, 1325.958397, 1397.168421, 1397.168421, 1532.96506, 2386.106979, 2386.106979, 1508.576167, 2535.464648, 2535.464648, 1508.576167, 2535.464648, 2535.464648};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITH_PART,
                      CWP_COMM_PAR_WITHOUT_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WoP - WoP - no-shared") {
    res_expected = {760.796303, 733.909723, 733.909723, 1325.958397, 1397.168421, 1397.168421, 1325.958397, 1397.168421, 1397.168421, 1532.96506, 2386.106979, 2386.106979, 1508.576167, 2535.464648, 2535.464648, 1508.576167, 2535.464648, 2535.464648};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITHOUT_PART,
                      CWP_COMM_PAR_WITHOUT_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WoP - WP - no-shared") {
    res_expected = {760.796303, 733.909723, 733.909723, 1325.958397, 1397.168421, 1397.168421, 1325.958397, 1397.168421, 1397.168421, 1532.96506, 2386.106979, 2386.106979, 1508.576167, 2535.464648, 2535.464648, 1508.576167, 2535.464648, 2535.464648};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITHOUT_PART,
                      CWP_COMM_PAR_WITH_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WP - WP - shared") {
    res_expected = {760.196303, 989.828328, 989.828328, 1329.958397, 1598.168421, 1598.168421, 1329.958397, 1598.168421, 1598.168421, 1568.417827, 2437.106979, 2437.106979, 1749.566865, 3003.381699, 3003.381699, 1749.566865, 3003.381699, 3003.381699};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITH_PART,
                      CWP_COMM_PAR_WITH_PART,
                      true,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WP - WoP - shared") {
    res_expected = {1521.592605, 1467.819446, 1469.200841, 2651.916794, 2794.336842, 2794.336842, 2651.916794, 2794.336842, 2794.336842, 1568.417827, 2437.106979, 2437.106979, 1749.566865, 3003.381699, 3003.381699, 1749.566865, 3003.381699, 3003.381699};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITH_PART,
                      CWP_COMM_PAR_WITHOUT_PART,
                      true,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WoP - WP - shared") {
    res_expected = {760.196303, 989.828328, 989.828328, 1329.958397, 1598.168421, 1598.168421, 1329.958397, 1598.168421, 1598.168421, 3065.93012, 4772.213957, 4772.213957, 3017.152335, 5070.929296, 5070.929296, 3017.152335, 5070.929296, 5070.929296};
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
    res_expected = {918., 1141., 1141., 1044., 1547., 1547., 0., 0., 0., 1422., 1908., 0., 1680., 2037., 0., 1680., 2037., 0.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITH_PART,
                      CWP_COMM_PAR_WITH_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WoP - WoP - no-shared") {
    res_expected = {918., 1141., 1141., 1044., 1547., 1547., 0., 0., 0., 1422., 1908., 0., 1680., 2037., 0., 1680., 2037., 0.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITHOUT_PART,
                      CWP_COMM_PAR_WITHOUT_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WP - WoP - no-shared") {
    res_expected = {918., 1141., 1141., 1044., 1547., 1547., 0., 0., 0., 1422., 1908., 0., 1680., 2037., 0., 1680., 2037., 0.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITH_PART,
                      CWP_COMM_PAR_WITHOUT_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WoP - WP - no-shared") {
    res_expected = {918., 1141., 1141., 1044., 1547., 1547., 0., 0., 0., 1422., 1908., 0., 1680., 2037., 0., 1680., 2037., 0.};
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
    res_expected = {918., 1141., 1141., 1044., 1547., 1547., 0., 0., 0., 1422., 1908., 0., 1680., 2037., 0., 1680., 2037., 0.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITH_PART,
                      CWP_COMM_PAR_WITH_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WP - WoP - no-shared") {
    res_expected = {918., 1141., 1141., 1044., 1547., 1547., 0., 0., 0., 1422., 1908., 0., 1680., 2037., 0., 1680., 2037., 0.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITH_PART,
                      CWP_COMM_PAR_WITHOUT_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WoP - WoP - no-shared") {
    res_expected = {918., 1141., 1141., 1044., 1547., 1547., 0., 0., 0., 1422., 1908., 0., 1680., 2037., 0., 1680., 2037., 0.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITHOUT_PART,
                      CWP_COMM_PAR_WITHOUT_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WoP - WP - no-shared") {
    res_expected = {918., 1141., 1141., 1044., 1547., 1547., 0., 0., 0., 1422., 1908., 0., 1680., 2037., 0., 1680., 2037., 0.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITHOUT_PART,
                      CWP_COMM_PAR_WITH_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WP - WP - shared") {
    res_expected = {918., 1311., 1311., 1044., 1775., 1775., 0., 0., 0., 1422., 1908., 0., 2044., 2471., 0., 2044., 2471., 0.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITH_PART,
                      CWP_COMM_PAR_WITH_PART,
                      true,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WP - WoP - shared") {
    res_expected = {1836., 2282., 2282., 2088., 3094., 3094., 0., 0., 0., 1422., 1908., 0., 2044., 2471., 0., 2044., 2471., 0.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITH_PART,
                      CWP_COMM_PAR_WITHOUT_PART,
                      true,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WoP - WP - shared") {
    res_expected = {918., 1311., 1311., 1044., 1775., 1775., 0., 0., 0., 2844., 3816., 0., 3360., 4074., 0., 3360., 4074., 0.};
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
    res_expected = {918., 1141., 1141., 1044., 1547., 1547., 0., 0., 0., 1422., 1908., 0., 1680., 2037., 0., 1680., 2037., 0.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITH_PART,
                      CWP_COMM_PAR_WITH_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WoP - WoP - no-shared") {
    res_expected = {918., 1141., 1141., 1044., 1547., 1547., 0., 0., 0., 1422., 1908., 0., 1680., 2037., 0., 1680., 2037., 0.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITHOUT_PART,
                      CWP_COMM_PAR_WITHOUT_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WP - WoP - no-shared") {
    res_expected = {918., 1141., 1141., 1044., 1547., 1547., 0., 0., 0., 1422., 1908., 0., 1680., 2037., 0., 1680., 2037., 0.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITH_PART,
                      CWP_COMM_PAR_WITHOUT_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WoP - WP - no-shared") {
    res_expected = {918., 1141., 1141., 1044., 1547., 1547., 0., 0., 0., 1422., 1908., 0., 1680., 2037., 0., 1680., 2037., 0.};
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
    res_expected = {918., 1141., 1141., 1044., 1547., 1547., 0., 0., 0., 1422., 1908., 0., 1680., 2037., 0., 1680., 2037., 0.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITH_PART,
                      CWP_COMM_PAR_WITH_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WP - WoP - no-shared") {
    res_expected = {918., 1141., 1141., 1044., 1547., 1547., 0., 0., 0., 1422., 1908., 0., 1680., 2037., 0., 1680., 2037., 0.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITH_PART,
                      CWP_COMM_PAR_WITHOUT_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WoP - WoP - no-shared") {
    res_expected = {918., 1141., 1141., 1044., 1547., 1547., 0., 0., 0., 1422., 1908., 0., 1680., 2037., 0., 1680., 2037., 0.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITHOUT_PART,
                      CWP_COMM_PAR_WITHOUT_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WoP - WP - no-shared") {
    res_expected = {918., 1141., 1141., 1044., 1547., 1547., 0., 0., 0., 1422., 1908., 0., 1680., 2037., 0., 1680., 2037., 0.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITHOUT_PART,
                      CWP_COMM_PAR_WITH_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WP - WP - shared") {
    res_expected = {918., 1311., 1311., 1044., 1775., 1775., 0., 0., 0., 1422., 1908., 0., 2044., 2471., 0., 2044., 2471., 0.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITH_PART,
                      CWP_COMM_PAR_WITH_PART,
                      true,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WP - WoP - shared") {
    res_expected = {1836., 2282., 2282., 2088., 3094., 3094., 0., 0., 0., 1422., 1908., 0., 2044., 2471., 0., 2044., 2471., 0.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITH_PART,
                      CWP_COMM_PAR_WITHOUT_PART,
                      true,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WoP - WP - shared") {
    res_expected = {918., 1311., 1311., 1044., 1775., 1775., 0., 0., 0., 2844., 3816., 0., 3360., 4074., 0., 3360., 4074., 0.};
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
    res_expected = {1260., 1788., 1788., 1538.4, 2073.8, 2073.8, 0., 0., 0., 2892., 3864.2, 0., 3612., 4777.2, 0., 3612., 4777.2, 0.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITH_PART,
                      CWP_COMM_PAR_WITH_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WoP - WoP - no-shared") {
    res_expected = {1260., 1788., 1788., 1538.4, 2073.8, 2073.8, 0., 0., 0., 2892., 3864.2, 0., 3612., 4777.2, 0., 3612., 4777.2, 0.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITHOUT_PART,
                      CWP_COMM_PAR_WITHOUT_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WP - WoP - no-shared") {
    res_expected = {1260., 1788., 1788., 1538.4, 2073.8, 2073.8, 0., 0., 0., 2892., 3864.2, 0., 3612., 4777.2, 0., 3612., 4777.2, 0.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITH_PART,
                      CWP_COMM_PAR_WITHOUT_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WoP - WP - no-shared") {
    res_expected = {1260., 1788., 1788., 1538.4, 2073.8, 2073.8, 0., 0., 0., 2892., 3864.2, 0., 3612., 4777.2, 0., 3612., 4777.2, 0.};
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
    res_expected = {1260., 1788., 1788., 1538.4, 2073.8, 2073.8, 0., 0., 0., 2892., 3864.2, 0., 3613., 4774., 0., 3613., 4774., 0.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITH_PART,
                      CWP_COMM_PAR_WITH_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WP - WoP - no-shared") {
    res_expected = {1260., 1788., 1788., 1538.4, 2073.8, 2073.8, 0., 0., 0., 2892., 3864.2, 0., 3613., 4774., 0., 3613., 4774., 0.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITH_PART,
                      CWP_COMM_PAR_WITHOUT_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WoP - WoP - no-shared") {
    res_expected = {1260., 1788., 1788., 1538.4, 2073.8, 2073.8, 0., 0., 0., 2892., 3864.2, 0., 3613., 4774., 0., 3613., 4774., 0.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITHOUT_PART,
                      CWP_COMM_PAR_WITHOUT_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WoP - WP - no-shared") {
    res_expected = {1260., 1788., 1788., 1538.4, 2073.8, 2073.8, 0., 0., 0., 2892., 3864.2, 0., 3613., 4774., 0., 3613., 4774., 0.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITHOUT_PART,
                      CWP_COMM_PAR_WITH_PART,
                      false,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WP - WP - shared") {
    res_expected = {1260., 1979., 1979., 1538.4, 2285., 2285., 0., 0., 0., 2892., 3864.2, 0., 4059., 5315.8, 0., 4059., 5315.8, 0.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITH_PART,
                      CWP_COMM_PAR_WITH_PART,
                      true,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WP - WoP - shared") {
    res_expected = {2520., 3576., 3576., 3076.8, 4147.6, 4147.6, 0., 0., 0., 2892., 3864.2, 0., 4059., 5315.8, 0., 4059., 5315.8, 0.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITH_PART,
                      CWP_COMM_PAR_WITHOUT_PART,
                      true,
                      spatial_interp,
                      false);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("WoP - WP - shared") {
    res_expected = {1260., 1979., 1979., 1538.4, 2285., 2285., 0., 0., 0., 5784., 7728.4, 0., 7226., 9548., 0., 7226., 9548., 0.};
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
    res_expected = {716.2, 0., 0., 0., 0., 0., 0., 0., 0., 1411., 0., 0., 0., 0., 0., 0., 0., 0.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITH_PART,
                      CWP_COMM_PAR_WITH_PART,
                      false,
                      CWP_SPATIAL_INTERP_FROM_INTERSECTION,
                      true);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("user-interp NEAREST_TARGETS_LEAST_SQUARES") {
    res_expected = {781.454105, 718.205851, 718.205851, 1324.454442, 1426.332082, 1426.332082, 1324.454442, 1426.332082, 1426.332082, 1507.209926, 2403.902901, 2403.902901, 1530.680693, 2501.77961, 2501.77961, 1530.680693, 2501.77961, 2501.77961};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITH_PART,
                      CWP_COMM_PAR_WITH_PART,
                      false,
                      CWP_SPATIAL_INTERP_FROM_NEAREST_TARGETS_LEAST_SQUARES,
                      true);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

  SUBCASE("user-interp LOCATION_OCTREE") {
    res_expected = {918., 1141., 1141., 1044., 1547., 1547., 0., 0., 0., 1422., 1908., 0., 1680., 2037., 0., 1680., 2037., 0.};
    res = _get_values(test_comm,
                      CWP_COMM_PAR_WITH_PART,
                      CWP_COMM_PAR_WITH_PART,
                      false,
                      CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE,
                      true);
    CHECK_EPS(res.data(), res_expected.data(), res.size());
  }

}