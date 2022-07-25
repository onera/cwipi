
cdef extern from "pdm_part_extension.h":
  ctypedef struct PDM_part_extension_t:
    pass

  ctypedef enum PDM_extend_type_t:
    PDM_EXTEND_FROM_FACE = 0,
    PDM_EXTEND_FROM_EDGE = 1,
    PDM_EXTEND_FROM_VTX  = 2

  PDM_part_extension_t* PDM_part_extension_create(int                n_domain,
                                                  int               *n_part,
                                                  PDM_extend_type_t  extend_type,
                                                  int                depth,
                                                  PDM_MPI_Comm       comm,
                                                  PDM_ownership_t    owner)

  void PDM_part_extension_compute(PDM_part_extension_t *part_ext)

  void PDM_part_extension_set_part(PDM_part_extension_t *part_ext,
                                   int                   i_domain,
                                   int                   i_part,
                                   int                   n_cell,
                                   int                   n_face,
                                   int                   n_face_part_bound,
                                   int                   n_face_group,
                                   int                   n_edge,
                                   int                   n_vtx,
                                   int                  *cell_face_idx,
                                   int                  *cell_face,
                                   int                  *face_cell,
                                   int                  *face_edge_idx,
                                   int                  *face_edge,
                                   int                  *face_vtx_idx,
                                   int                  *face_vtx,
                                   int                  *edge_vtx,
                                   int                  *face_bound_idx,
                                   int                  *face_bound,
                                   int                  *face_join_idx,
                                   int                  *face_join,
                                   int                  *face_part_bound_proc_idx,
                                   int                  *face_part_bound_part_idx,
                                   int                  *face_part_bound,
                                   int                  *vtx_part_bound_proc_idx,
                                   int                  *vtx_part_bound_part_idx,
                                   int                  *vtx_part_bound,
                                   PDM_g_num_t          *cell_ln_to_gn,
                                   PDM_g_num_t          *face_ln_to_gn,
                                   PDM_g_num_t          *edge_ln_to_gn,
                                   PDM_g_num_t          *vtx_ln_to_gn,
                                   PDM_g_num_t          *face_group_ln_to_gn,
                                   double               *vtx_coord)

  void PDM_part_extension_free(PDM_part_extension_t *part_ext)

  int PDM_part_extension_connectivity_get(PDM_part_extension_t     *part_ext,
                                          int                       i_domain,
                                          int                       i_part,
                                          PDM_connectivity_type_t   connectivity_type,
                                          int                     **connect,
                                          int                     **connect_idx)

  int PDM_part_extension_ln_to_gn_get(PDM_part_extension_t     *part_ext,
                                      int                       i_domain,
                                      int                       i_part,
                                      PDM_mesh_entities_t       connectivity_type,
                                      PDM_g_num_t             **ln_to_gn)

  int PDM_part_extension_interface_get(PDM_part_extension_t     *part_ext,
                                       int                       i_domain,
                                       int                       i_part,
                                       PDM_mesh_entities_t       mesh_entity,
                                       int                     **interface_no);

  int PDM_part_extension_composed_interface_get(PDM_part_extension_t     *part_ext,
                                                int                     **composed_interface_idx,
                                                int                     **composed_interface,
                                                PDM_g_num_t             **composed_ln_to_gn_sorted);


  int PDM_part_extension_group_get(PDM_part_extension_t     *part_ext,
                                   int                       i_domain,
                                   int                       i_part,
                                   PDM_mesh_entities_t       mesh_entity,
                                   int                     **connect,
                                   int                     **connect_idx,
                                   PDM_g_num_t             **ln_to_gn)

  int PDM_part_extension_coord_get(PDM_part_extension_t     *part_ext,
                                   int                       i_domain,
                                   int                       i_part,
                                   double                  **vtx_coord)

  void PDM_part_extension_part_domain_interface_shared_set(PDM_part_extension_t        *part_ext,
                                                           PDM_part_domain_interface_t *pdi);


cdef class PartExtension:
  """
  """
  cdef PDM_part_extension_t* _part_ext
  cdef object pdi
  # ------------------------------------------------------------------
  def __cinit__(self,
                n_domain,
                NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] n_part,
                PDM_extend_type_t                             extend_type,
                depth,
                MPI.Comm                                      comm):
    """
    """
    cdef MPI.MPI_Comm c_comm = comm.ob_mpi

    self._part_ext = PDM_part_extension_create(n_domain,
                                      <int*>   n_part.data,
                                               extend_type,
                                               depth,
                                               PDM_MPI_mpi_2_pdm_mpi_comm (<void *> &c_comm),
                                               PDM_OWNERSHIP_USER)

  # ------------------------------------------------------------------
  def set_part(self,
               int i_domain,
               int i_part,
               int n_cell,
               int n_face,
               int n_face_part_bound,
               int n_face_group,
               int n_edge,
               int n_vtx,
               NPY.ndarray[NPY.int32_t, mode='c', ndim=1]    cell_face_idx,
               NPY.ndarray[NPY.int32_t, mode='c', ndim=1]    cell_face,
               NPY.ndarray[NPY.int32_t, mode='c', ndim=1]    face_cell,
               NPY.ndarray[NPY.int32_t, mode='c', ndim=1]    face_edge_idx,
               NPY.ndarray[NPY.int32_t, mode='c', ndim=1]    face_edge,
               NPY.ndarray[NPY.int32_t, mode='c', ndim=1]    face_vtx_idx,
               NPY.ndarray[NPY.int32_t, mode='c', ndim=1]    face_vtx,
               NPY.ndarray[NPY.int32_t, mode='c', ndim=1]    edge_vtx,
               NPY.ndarray[NPY.int32_t, mode='c', ndim=1]    face_bound_idx,
               NPY.ndarray[NPY.int32_t, mode='c', ndim=1]    face_bound,
               NPY.ndarray[NPY.int32_t, mode='c', ndim=1]    face_join_idx,
               NPY.ndarray[NPY.int32_t, mode='c', ndim=1]    face_join,
               NPY.ndarray[NPY.int32_t, mode='c', ndim=1]    face_part_bound_proc_idx,
               NPY.ndarray[NPY.int32_t, mode='c', ndim=1]    face_part_bound_part_idx,
               NPY.ndarray[NPY.int32_t, mode='c', ndim=1]    face_part_bound,
               NPY.ndarray[NPY.int32_t, mode='c', ndim=1]    vtx_part_bound_proc_idx,
               NPY.ndarray[NPY.int32_t, mode='c', ndim=1]    vtx_part_bound_part_idx,
               NPY.ndarray[NPY.int32_t, mode='c', ndim=1]    vtx_part_bound,
               NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] cell_ln_to_gn,
               NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] face_ln_to_gn,
               NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] edge_ln_to_gn,
               NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] vtx_ln_to_gn,
               NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] face_group_ln_to_gn,
               NPY.ndarray[NPY.double_t  , mode='c', ndim=1] vtx_coord):
    """
    """

    cdef int * cell_face_idx_data
    if (cell_face_idx is None):
      cell_face_idx_data = NULL
    else:
      cell_face_idx_data = <int *> cell_face_idx.data

    cdef int * cell_face_data
    if (cell_face is None):
      cell_face_data = NULL
    else:
      cell_face_data = <int *> cell_face.data

    cdef int * face_cell_data
    if (face_cell is None):
      face_cell_data = NULL
    else:
      face_cell_data = <int *> face_cell.data

    cdef int * face_edge_idx_data
    if (face_edge_idx is None):
      face_edge_idx_data = NULL
    else:
      face_edge_idx_data = <int *> face_edge_idx.data

    cdef int * face_edge_data
    if (face_edge is None):
      face_edge_data = NULL
    else:
      face_edge_data = <int *> face_edge.data

    cdef int * face_vtx_idx_data
    if (face_vtx_idx is None):
      face_vtx_idx_data = NULL
    else:
      face_vtx_idx_data = <int *> face_vtx_idx.data

    cdef int * face_vtx_data
    if (face_vtx is None):
      face_vtx_data = NULL
    else:
      face_vtx_data = <int *> face_vtx.data

    cdef int * edge_vtx_data
    if (edge_vtx is None):
      edge_vtx_data = NULL
    else:
      edge_vtx_data = <int *> edge_vtx.data

    cdef int * face_bound_data
    if (face_bound is None):
      face_bound_data = NULL
    else:
      face_bound_data = <int *> face_bound.data

    cdef int * face_bound_idx_data
    if (face_bound_idx is None):
      face_bound_idx_data = NULL
    else:
      face_bound_idx_data = <int *> face_bound_idx.data

    cdef int * face_join_idx_data
    if (face_join_idx is None):
      face_join_idx_data = NULL
    else:
      face_join_idx_data = <int *> face_join_idx.data

    cdef int * face_join_data
    if (face_join is None):
      face_join_data = NULL
    else:
      face_join_data = <int *> face_join.data

    cdef int * face_part_bound_proc_idx_data
    if (face_part_bound_proc_idx is None):
      face_part_bound_proc_idx_data = NULL
    else:
      face_part_bound_proc_idx_data = <int *> face_part_bound_proc_idx.data

    cdef int * face_part_bound_part_idx_data
    if (face_part_bound_part_idx is None):
      face_part_bound_part_idx_data = NULL
    else:
      face_part_bound_part_idx_data = <int *> face_part_bound_part_idx.data

    cdef int * face_part_bound_data
    if (face_part_bound is None):
      face_part_bound_data = NULL
    else:
      face_part_bound_data = <int *> face_part_bound.data

    cdef int * vtx_part_bound_proc_idx_data
    if (vtx_part_bound_proc_idx is None):
      vtx_part_bound_proc_idx_data = NULL
    else:
      vtx_part_bound_proc_idx_data = <int *> vtx_part_bound_proc_idx.data

    cdef int * vtx_part_bound_part_idx_data
    if (vtx_part_bound_part_idx is None):
      vtx_part_bound_part_idx_data = NULL
    else:
      vtx_part_bound_part_idx_data = <int *> vtx_part_bound_part_idx.data

    cdef int * vtx_part_bound_data
    if (vtx_part_bound is None):
      vtx_part_bound_data = NULL
    else:
      vtx_part_bound_data = <int *> vtx_part_bound.data

    cdef PDM_g_num_t * cell_ln_to_gn_data
    if (cell_ln_to_gn is None):
      cell_ln_to_gn_data = NULL
    else:
      cell_ln_to_gn_data = <PDM_g_num_t *> cell_ln_to_gn.data

    cdef PDM_g_num_t * face_ln_to_gn_data
    if (face_ln_to_gn is None):
      face_ln_to_gn_data = NULL
    else:
      face_ln_to_gn_data = <PDM_g_num_t *> face_ln_to_gn.data

    cdef PDM_g_num_t * edge_ln_to_gn_data
    if (edge_ln_to_gn is None):
      edge_ln_to_gn_data = NULL
    else:
      edge_ln_to_gn_data = <PDM_g_num_t *> edge_ln_to_gn.data

    cdef PDM_g_num_t * vtx_ln_to_gn_data
    if (vtx_ln_to_gn is None):
      vtx_ln_to_gn_data = NULL
    else:
      vtx_ln_to_gn_data = <PDM_g_num_t *> vtx_ln_to_gn.data

    cdef PDM_g_num_t * face_group_ln_to_gn_data
    if (face_group_ln_to_gn is None):
      face_group_ln_to_gn_data = NULL
    else:
      face_group_ln_to_gn_data = <PDM_g_num_t *> face_group_ln_to_gn.data

    PDM_part_extension_set_part(self._part_ext,
                                i_domain,
                                i_part,
                                n_cell,
                                n_face,
                                n_face_part_bound,
                                n_face_group,
                                n_edge,
                                n_vtx,
                                cell_face_idx_data,
                                cell_face_data,
                                face_cell_data,
                                face_edge_idx_data,
                                face_edge_data,
                                face_vtx_idx_data,
                                face_vtx_data,
                                edge_vtx_data,
                                face_bound_idx_data,
                                face_bound_data,
                                face_join_idx_data,
                                face_join_data,
                                face_part_bound_proc_idx_data,
                                face_part_bound_part_idx_data,
                                face_part_bound_data,
                                vtx_part_bound_proc_idx_data,
                                vtx_part_bound_part_idx_data,
                                vtx_part_bound_data,
                                cell_ln_to_gn_data,
                                face_ln_to_gn_data,
                                edge_ln_to_gn_data,
                                vtx_ln_to_gn_data,
                                face_group_ln_to_gn_data,
                     <double*>  vtx_coord.data)

  # ------------------------------------------------------------------
  def compute(self):
    PDM_part_extension_compute(self._part_ext)


  # ------------------------------------------------------------------
  def part_domain_interface_shared_set(self, PartDomainInterface pdi):
    """
    """
    self.pdi = pdi # Keep alive
    PDM_part_extension_part_domain_interface_shared_set(self._part_ext, pdi.pdi)


  # ------------------------------------------------------------------
  def get_connectivity(self,
                       int i_domain,
                       int i_part,
                       PDM_connectivity_type_t   connectivity_type):
    """
    """
    cdef int *connect,
    cdef int *connect_idx,
    cdef int size

    size = PDM_part_extension_connectivity_get(self._part_ext, i_domain, i_part, connectivity_type, &connect, &connect_idx)

    if (connect_idx == NULL) :
      np_connect_idx = None
    else :
      dim = <NPY.npy_intp> (size+1)
      np_connect_idx = NPY.PyArray_SimpleNewFromData(1,
                                                  &dim,
                                                  NPY.NPY_INT32,
                                                  <void *> connect_idx)
      PyArray_ENABLEFLAGS(np_connect_idx, NPY.NPY_OWNDATA)

    if (connect == NULL) :
      np_connect = None
    else :
      dim = <NPY.npy_intp> connect_idx[size]
      np_connect = NPY.PyArray_SimpleNewFromData(1,
                                                  &dim,
                                                  NPY.NPY_INT32,
                                                  <void *> connect)
      PyArray_ENABLEFLAGS(np_connect, NPY.NPY_OWNDATA)

    return (np_connect, np_connect_idx)

  # ------------------------------------------------------------------
  def get_ln_to_gn(self,
                       int i_domain,
                       int i_part,
                       PDM_mesh_entities_t mesh_ety_type):
    """
    """
    cdef PDM_g_num_t *ln_to_gn,
    cdef int size

    size = PDM_part_extension_ln_to_gn_get(self._part_ext, i_domain, i_part, mesh_ety_type, &ln_to_gn)

    if (ln_to_gn == NULL) :
      # np_ln_to_gn = None
      dim = <NPY.npy_intp> 0
      np_ln_to_gn = NPY.PyArray_SimpleNew(1, &dim, PDM_G_NUM_NPY_INT)
    else :
      dim = <NPY.npy_intp> size
      np_ln_to_gn = NPY.PyArray_SimpleNewFromData(1,
                                                  &dim,
                                                  PDM_G_NUM_NPY_INT,
                                                  <void *> ln_to_gn)
      PyArray_ENABLEFLAGS(np_ln_to_gn, NPY.NPY_OWNDATA)

    return np_ln_to_gn

  # ------------------------------------------------------------------
  def get_interface(self,
                    int i_domain,
                    int i_part,
                    PDM_mesh_entities_t mesh_ety_type):
    """
    """
    cdef int *interface_no,
    cdef int size

    size = PDM_part_extension_interface_get(self._part_ext, i_domain, i_part, mesh_ety_type, &interface_no)

    if (interface_no == NULL) :
      np_interface_no = None
    else :
      dim = <NPY.npy_intp> size
      np_interface_no = NPY.PyArray_SimpleNewFromData(1,
                                                  &dim,
                                                  NPY.NPY_INT32,
                                                  <void *> interface_no)
      PyArray_ENABLEFLAGS(np_interface_no, NPY.NPY_OWNDATA)

    return np_interface_no

  # ------------------------------------------------------------------
  def get_composed_interface(self):
    """
    """
    cdef int         *composed_interface_idx,
    cdef int         *composed_interface,
    cdef PDM_g_num_t *composed_ln_to_gn_sorted,
    cdef int          n_composed_interface
    cdef NPY.npy_intp dim

    n_composed_interface = PDM_part_extension_composed_interface_get(self._part_ext,
                                                                     &composed_interface_idx,
                                                                     &composed_interface,
                                                                     &composed_ln_to_gn_sorted)

    if (composed_interface_idx == NULL) :
      np_composed_interface_idx = None
    else :
      dim = <NPY.npy_intp> n_composed_interface+1
      np_composed_interface_idx = NPY.PyArray_SimpleNewFromData(1,
                                                                &dim,
                                                                NPY.NPY_INT32,
                                                                <void *> composed_interface_idx)
      PyArray_ENABLEFLAGS(np_composed_interface_idx, NPY.NPY_OWNDATA)

    if (composed_interface == NULL) :
      np_composed_interface = None
    else :
      dim = <NPY.npy_intp> np_composed_interface_idx[n_composed_interface]
      np_composed_interface = NPY.PyArray_SimpleNewFromData(1,
                                                  &dim,
                                                  NPY.NPY_INT32,
                                                  <void *> composed_interface)
      PyArray_ENABLEFLAGS(np_composed_interface, NPY.NPY_OWNDATA)

    if (composed_ln_to_gn_sorted == NULL) :
      np_composed_ln_to_gn_sorted = None
    else :
      dim = <NPY.npy_intp> n_composed_interface
      np_composed_ln_to_gn_sorted = NPY.PyArray_SimpleNewFromData(1,
                                                  &dim,
                                                  PDM_G_NUM_NPY_INT,
                                                  <void *> composed_ln_to_gn_sorted)
      PyArray_ENABLEFLAGS(np_composed_ln_to_gn_sorted, NPY.NPY_OWNDATA)

    return (np_composed_interface_idx, np_composed_interface, np_composed_ln_to_gn_sorted)

  # ------------------------------------------------------------------
  def get_group(self,
                int i_domain,
                int i_part,
                PDM_mesh_entities_t mesh_ety_type):
    """
    """
    cdef int *entity_group_idx
    cdef int *entity_group
    cdef PDM_g_num_t *group_ln_to_gn
    cdef int size

    n_group = PDM_part_extension_group_get(self._part_ext, i_domain, i_part,
                                        mesh_ety_type,
                                        &entity_group,
                                        &entity_group_idx,
                                        &group_ln_to_gn)

    if (group_ln_to_gn == NULL) :
      np_group_ln_to_gn = None
    else :
      dim = <NPY.npy_intp> entity_group_idx[n_group]
      np_group_ln_to_gn = NPY.PyArray_SimpleNewFromData(1,
                                                  &dim,
                                                  PDM_G_NUM_NPY_INT,
                                                  <void *> group_ln_to_gn)
      PyArray_ENABLEFLAGS(np_group_ln_to_gn, NPY.NPY_OWNDATA)

    if (entity_group_idx == NULL) :
      np_entity_group_idx = None
    else :
      dim = <NPY.npy_intp> (n_group+1)
      np_entity_group_idx = NPY.PyArray_SimpleNewFromData(1,
                                                          &dim,
                                                          NPY.NPY_INT32,
                                                          <void *> entity_group_idx)
      PyArray_ENABLEFLAGS(np_entity_group_idx, NPY.NPY_OWNDATA)

    if (entity_group == NULL) :
      np_entity_group = None
    else :
      dim = <NPY.npy_intp> entity_group_idx[n_group]
      np_entity_group = NPY.PyArray_SimpleNewFromData(1,
                                                      &dim,
                                                      NPY.NPY_INT32,
                                                      <void *> entity_group)
      PyArray_ENABLEFLAGS(np_entity_group, NPY.NPY_OWNDATA)

    return (np_entity_group_idx, np_entity_group, np_group_ln_to_gn)

  # ------------------------------------------------------------------
  def get_coord(self, int i_domain, int i_part):
    """
    """
    cdef double *coord
    cdef int size

    size = PDM_part_extension_coord_get(self._part_ext, i_domain, i_part, &coord)

    if (coord == NULL) :
      np_coord = None
    else :
      dim = <NPY.npy_intp> (3 * size)
      np_coord = NPY.PyArray_SimpleNewFromData(1,
                                                  &dim,
                                                  NPY.NPY_DOUBLE,
                                                  <void *> coord)
      PyArray_ENABLEFLAGS(np_coord, NPY.NPY_OWNDATA)

    return np_coord

  # ------------------------------------------------------------------
  def __dealloc__(self):
    """
    """
    PDM_part_extension_free(self._part_ext)
