# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# > Wrapping of functions
cdef extern from "pdm_extract_part.h":
  ctypedef struct PDM_extract_part_t:
    pass
  PDM_extract_part_t* PDM_extract_part_create(int                    dim,
                                              int                    n_part_in,
                                              int                    n_part_out,
                                              PDM_bool_t             equilibrate,
                                              PDM_split_dual_t       split_dual_method,
                                              PDM_ownership_t        ownership,
                                              PDM_MPI_Comm           comm);
  void PDM_extract_part_compute(PDM_extract_part_t        *extrp);
  void PDM_extract_part_selected_lnum_set(PDM_extract_part_t       *extrp,
                                          int                       i_part,
                                          int                       n_extract,
                                          int                      *extract_lnum);

  void PDM_extract_part_part_set(PDM_extract_part_t        *extrp,
                                 int                       i_part,
                                 int                       n_cell,
                                 int                       n_face,
                                 int                       n_edge,
                                 int                       n_vtx,
                                 int                      *cell_face_idx,
                                 int                      *cell_face,
                                 int                      *face_edge_idx,
                                 int                      *face_edge,
                                 int                      *edge_vtx,
                                 int                      *face_vtx_idx,
                                 int                      *face_vtx,
                                 PDM_g_num_t              *cell_ln_to_gn,
                                 PDM_g_num_t              *face_ln_to_gn,
                                 PDM_g_num_t              *edge_ln_to_gn,
                                 PDM_g_num_t              *vtx_ln_to_gn,
                                 double                   *vtx_coord);

  int PDM_extract_part_n_entity_get(PDM_extract_part_t       *extrp,
                                    int                       i_part_out,
                                    PDM_mesh_entities_t       entity_type);

  int PDM_extract_part_connectivity_get(PDM_extract_part_t        *extrp,
                                        int                        i_part_out,
                                        PDM_connectivity_type_t    connectivity_type,
                                        int                      **connect,
                                        int                      **connect_idx,
                                        PDM_ownership_t           ownership);
  int PDM_extract_part_ln_to_gn_get(PDM_extract_part_t        *extrp,
                                    int                        i_part_out,
                                    PDM_mesh_entities_t        entity_type,
                                    PDM_g_num_t              **pentity_ln_to_gn,
                                    PDM_ownership_t            ownership);

  int PDM_extract_part_parent_ln_to_gn_get(PDM_extract_part_t        *extrp,
                                           int                        i_part_out,
                                           PDM_mesh_entities_t        entity_type,
                                           PDM_g_num_t              **parent_entity_ln_to_gn,
                                           PDM_ownership_t            ownership);

  int PDM_extract_part_vtx_coord_get(PDM_extract_part_t         *extrp,
                                     int                        i_part_out,
                                     double                   **pvtx_coord,
                                     PDM_ownership_t            ownership);

  void PDM_extract_part_free(PDM_extract_part_t  *extrp);

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# > Class definition
cdef class ExtractPart:
  """

  """
  # --------------------------------------------------------------------------
  # > Class attributes
  cdef PDM_extract_part_t* _extrp
  keep_alive = []
  # --------------------------------------------------------------------------

  # ------------------------------------------------------------------
  def __cinit__(self,
                int               dim,
                int               n_part_in,
                int               n_part_out,
                PDM_bool_t        equilibrate,
                PDM_split_dual_t  split_dual_method,
                MPI.Comm          comm):
    """
    Compute the distance from point clouds to a surface
    """
    cdef MPI.MPI_Comm c_comm = comm.ob_mpi
    self._extrp =  PDM_extract_part_create(dim,
                                           n_part_in,
                                           n_part_out,
                                           equilibrate,
                                           split_dual_method,
                                           PDM_OWNERSHIP_USER,
                                           PDM_MPI_mpi_2_pdm_mpi_comm (<void *> &c_comm));

  # ------------------------------------------------------------------
  def selected_lnum_set(self,
                        int i_part,
                        NPY.ndarray[NPY.int32_t, mode='c', ndim=1] extract_lnum):
    """
    """
    cdef int n_extract = extract_lnum.shape[0]
    PDM_extract_part_selected_lnum_set(self._extrp,
                                       i_part,
                                       n_extract,
                              <int* >  extract_lnum.data);


  # ------------------------------------------------------------------
  def part_set(self,
               int                                           i_part,
               int                                           n_cell,
               int                                           n_face,
               int                                           n_edge,
               int                                           n_vtx,
               NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] cell_face_idx,
               NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] cell_face    ,
               NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] face_edge_idx,
               NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] face_edge    ,
               NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] edge_vtx     ,
               NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] face_vtx_idx ,
               NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] face_vtx     ,
               NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] cell_ln_to_gn,
               NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] face_ln_to_gn,
               NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] edge_ln_to_gn,
               NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] vtx_ln_to_gn ,
               NPY.ndarray[NPY.double_t  , mode='c', ndim=1] coords):
    """
    """
    self.keep_alive.append(cell_face_idx)
    self.keep_alive.append(cell_face)
    self.keep_alive.append(face_edge_idx)
    self.keep_alive.append(face_edge)
    self.keep_alive.append(edge_vtx)
    self.keep_alive.append(face_vtx_idx)
    self.keep_alive.append(face_vtx)
    self.keep_alive.append(cell_ln_to_gn)
    self.keep_alive.append(face_ln_to_gn)
    self.keep_alive.append(edge_ln_to_gn)
    self.keep_alive.append(vtx_ln_to_gn)
    self.keep_alive.append(coords)

    PDM_extract_part_part_set(self._extrp,
                              i_part,
                              n_cell,
                              n_face,
                              n_edge,
                              n_vtx,
             <int         *>  cell_face_idx.data,
             <int         *>  cell_face    .data,
             <int         *>  face_edge_idx.data,
             <int         *>  face_edge    .data,
             <int         *>  edge_vtx     .data,
             <int         *>  face_vtx_idx .data,
             <int         *>  face_vtx     .data,
             <PDM_g_num_t *> cell_ln_to_gn.data,
             <PDM_g_num_t *> face_ln_to_gn.data,
             <PDM_g_num_t *> edge_ln_to_gn.data,
             <PDM_g_num_t *> vtx_ln_to_gn .data,
             <double      *> coords       .data)

  # ------------------------------------------------------------------
  def compute(self):
    """
    """
    PDM_extract_part_compute(self._extrp)

  # ------------------------------------------------------------------
  def n_entity_get(self,
                   int ipart,
                   PDM_mesh_entities_t entity_type):
    """
    """
    return PDM_extract_part_n_entity_get(self._extrp, ipart, entity_type)

  # ------------------------------------------------------------------
  def connectivity_get(self,
                       int ipart,
                       PDM_connectivity_type_t connectivity_type):
    """
    """
    cdef int  n_entity
    cdef int *connect
    cdef int *connect_idx

    n_entity = PDM_extract_part_connectivity_get(self._extrp,
                                                 ipart,
                                                 connectivity_type,
                                                 &connect,
                                                 &connect_idx,
                                                 PDM_OWNERSHIP_USER)

    np_connect_idx = None
    np_connect     = None
    if(connect_idx != NULL):
      np_connect_idx = create_numpy_i(connect_idx, n_entity+1           )
      np_connect     = create_numpy_i(connect    , connect_idx[n_entity])
    else:
      np_connect  = create_numpy_i(connect    , 2 * n_entity)

    # return PDM_extract_part_n_entity_get(self._extrp, ipart, entity_type)
    return np_connect_idx, np_connect

  # ------------------------------------------------------------------
  def ln_to_gn_get(self,
                   int ipart,
                   PDM_mesh_entities_t entity_type):
    """
    """
    cdef int  n_entity
    cdef PDM_g_num_t *entity_ln_to_gn

    n_entity = PDM_extract_part_ln_to_gn_get(self._extrp,
                                             ipart,
                                             entity_type,
                                             &entity_ln_to_gn,
                                             PDM_OWNERSHIP_USER)
    return create_numpy_pdm_gnum(entity_ln_to_gn, n_entity)

  # ------------------------------------------------------------------
  def parent_ln_to_gn_get(self,
                          int ipart,
                          PDM_mesh_entities_t entity_type):
    """
    """
    cdef int  n_entity
    cdef PDM_g_num_t *parent_ln_to_gn

    n_entity = PDM_extract_part_parent_ln_to_gn_get(self._extrp,
                                                    ipart,
                                                    entity_type,
                                                    &parent_ln_to_gn,
                                                    PDM_OWNERSHIP_USER)
    return create_numpy_pdm_gnum(parent_ln_to_gn, n_entity)

  # ------------------------------------------------------------------
  def vtx_coord_get(self,
                    int ipart):
    """
    """
    cdef int     n_vtx
    cdef double *pvtx_coord

    n_vtx = PDM_extract_part_vtx_coord_get(self._extrp, ipart, &pvtx_coord, PDM_OWNERSHIP_USER)

    return create_numpy_d(pvtx_coord, 3 * n_vtx)

  # ------------------------------------------------------------------
  def __dealloc__(self):
      PDM_extract_part_free(self._extrp)
