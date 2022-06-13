
cdef extern from "pdm_dmesh.h":
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Wrapping of structure
    ctypedef struct PDM_dmesh_t:
      pass
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Wrapping of function
    PDM_dmesh_t* PDM_dmesh_create(PDM_ownership_t owner,
                                  int             dn_cell,
                                  int             dn_face,
                                  int             dn_edge,
                                  int             dn_vtx,
                                  int             n_bnd,
                                  int             n_join,
                                  PDM_MPI_Comm    comm)

    void PDM_dmesh_set(PDM_dmesh_t  *dm,
                       double       *dvtx_coord,
                       int          *dface_vtx_idx,
                       PDM_g_num_t  *dface_vtx,
                       PDM_g_num_t  *dface_cell,
                       int          *dface_bound_idx,
                       PDM_g_num_t  *dface_bound,
                       int          *join_g_dms,
                       int          *dface_join_idx,
                       PDM_g_num_t  *dface_join)

    void PDM_dmesh_dims_get(PDM_dmesh_t *dm,
                            int         *dn_cell,
                            int         *dn_face,
                            int         *dn_edge,
                            int         *dn_vtx,
                            int         *n_bnd,
                            int         *n_joins)

    void PDM_dmesh_data_get(PDM_dmesh_t   *dm,
                            double       **dvtx_coord,
                            int          **dface_vtx_idx,
                            PDM_g_num_t  **dface_vtx,
                            PDM_g_num_t  **dface_cell,
                            int          **dface_bound_idx,
                            PDM_g_num_t  **dface_bound,
                            int          **join_g_dms,
                            int          **dface_join_idx,
                            PDM_g_num_t  **dface_join)
    int PDM_dmesh_connectivity_get(PDM_dmesh_t              *dmesh,
                                   PDM_connectivity_type_t   connectivity_type,
                                   PDM_g_num_t             **connect,
                                   int                     **connect_idx,
                                   PDM_ownership_t           ownership)
    int PDM_dmesh_distrib_get(PDM_dmesh_t              *dmesh,
                               PDM_mesh_entities_t       entity,
                               PDM_g_num_t             **distrib)
    int  PDM_dmesh_bound_get(PDM_dmesh_t       *dmesh,
                             PDM_bound_type_t   bound_type,
                             PDM_g_num_t      **connect,
                             int              **connect_idx,
                             PDM_ownership_t    ownership)
    void PDM_dmesh_free(PDM_dmesh_t   *dm)
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# ------------------------------------------------------------------
# cdef void PDM_pydmesh_free(object caps):
#   # print("PDM_pydmesh_free", PyCapsule_GetName(caps))
#   cdef PDM_dmesh_t* dm = <PDM_dmesh_t *> PyCapsule_GetPointer(caps, <const char*> PyCapsule_GetName(caps))
#   PDM_dmesh_free(dm);

# ------------------------------------------------------------------
cdef class DistributedMeshCaspule:
  """
  """
  # ************************************************************************
  # > Class attributes
  cdef PDM_dmesh_t* _dm
  # ************************************************************************
  # ------------------------------------------------------------------------
  def __cinit__(self, object caps):
    """
    """
    # print("DistributedMeshCaspule", PyCapsule_GetName(caps))
    cdef PDM_dmesh_t* dm = <PDM_dmesh_t *> PyCapsule_GetPointer(caps, NULL)
    self._dm = dm;

  # ------------------------------------------------------------------------
  def dmesh_connectivity_get(self, PDM_connectivity_type_t connectivity_type):
    """
    """
    return dmesh_connectivity_get(self, connectivity_type)

  # ------------------------------------------------------------------------
  def dmesh_distrib_get(self, PDM_mesh_entities_t entity_type):
    """
    """
    return dmesh_distrib_get(self, entity_type)

  # ------------------------------------------------------------------------
  def dmesh_bound_get(self, PDM_bound_type_t bound_type):
    """
    """
    return dmesh_bound_get(self, bound_type)

  # ------------------------------------------------------------------------
  def __dealloc__(self):
    """
       Use the free method of PDM Lib
    """
    PDM_dmesh_free(self._dm)

# ------------------------------------------------------------------
cdef class DistributedMesh:
  """
     DistributedMesh: Distributed mesh structure
  """
  # ************************************************************************
  # > Class attributes
  cdef PDM_dmesh_t* _dm
  cdef int          n_rank
  # ************************************************************************
  # ------------------------------------------------------------------------
  def __cinit__(self, MPI.Comm comm,
                      dn_cell,
                      dn_face,
                      dn_edge,
                      dn_vtx,
                      n_bnd,
                      n_join):
    """
    TODOUX
    """
    # ************************************************************************
    # > Declaration
    # ************************************************************************
    cdef MPI.MPI_Comm c_comm = comm.ob_mpi
    self.n_rank = comm.Get_size()
    # ::::::::::::::::::::::::::::::::::::::::::::::::::
    self._dm = PDM_dmesh_create(PDM_OWNERSHIP_UNGET_RESULT_IS_FREE,
                                dn_cell,
                                dn_face,
                                dn_edge,
                                dn_vtx,
                                n_bnd,
                                n_join,
                                PDM_MPI_mpi_2_pdm_mpi_comm (<void *> &c_comm))
    # ::::::::::::::::::::::::::::::::::::::::::::::::::

  # ------------------------------------------------------------------------
  def dmesh_set(self, NPY.ndarray[NPY.double_t  , mode='c', ndim=1] dvtx_coord   not None,
                      NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] dface_vtx_idx,
                      NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] dface_vtx    not None,
                      NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] dface_cell,
                      NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] dface_bound_idx,
                      NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] dface_bound,
                      NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] join_g_dms,
                      NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] dface_join_idx,
                      NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] dface_join):
    """
    """
    # print("dvtx_coord", dvtx_coord)
    # print("dface_vtx_idx", dface_vtx_idx)
    # print("dface_vtx", dface_vtx)
    # print("dface_cell", dface_cell)
    # print("dface_bound_idx", dface_bound_idx)
    # print("dface_bound", dface_bound)

    PDM_dmesh_set(self._dm,
                  <double*>      dvtx_coord.data,
                  <int*>         dface_vtx_idx.data,
                  <PDM_g_num_t*> dface_vtx.data,
                  <PDM_g_num_t*> dface_cell.data,
                  <int*>         dface_bound_idx.data,
                  <PDM_g_num_t*> dface_bound.data,
                  <int*>         join_g_dms.data,
                  <int*>         dface_join_idx.data,
                  <PDM_g_num_t*> dface_join.data)
  # ------------------------------------------------------------------------
  def dmesh_connectivity_get(self, PDM_connectivity_type_t connectivity_type):
    """
    """
    return dmesh_connectivity_get(self, connectivity_type)

  # ------------------------------------------------------------------------
  def dmesh_distrib_get(self, PDM_mesh_entities_t entity_type):
    """
    """
    return dmesh_connectivity_get(self, entity_type)
  # ------------------------------------------------------------------------
  def dmesh_bound_get(self, PDM_bound_type_t bound_type):
    """
    """
    return dmesh_bound_get(self, bound_type)

  # ------------------------------------------------------------------------
  def __dealloc__(self):
    """
       Use the free method of PDM Lib
    """
    # ************************************************************************
    # > Declaration
    # ************************************************************************
    # print("DistributedMesh::__dealloc__")
    PDM_dmesh_free(self._dm)
    # print("DistributedMesh::__dealloc__")

ctypedef fused DMesh:
  DistributedMesh
  DistributedMeshCaspule

# ------------------------------------------------------------------------
def dmesh_connectivity_get(DMesh pydm, PDM_connectivity_type_t connectivity_type):
  """
  """
  # ************************************************************************
  # > Declaration
  cdef int          *connect_idx
  cdef PDM_g_num_t  *connect
  cdef int           dn_entity
  cdef NPY.npy_intp  dim
  # ************************************************************************
  dn_entity = PDM_dmesh_connectivity_get(pydm._dm,
                                         connectivity_type,
                                         &connect,
                                         &connect_idx,
                                         PDM_OWNERSHIP_USER)
  if (connect_idx == NULL) :
      np_connect_idx = None
  else :
      dim = <NPY.npy_intp> dn_entity + 1
      np_connect_idx = NPY.PyArray_SimpleNewFromData(1,
                                                     &dim,
                                                     NPY.NPY_INT32,
                                                     <void *> connect_idx)
  PyArray_ENABLEFLAGS(np_connect_idx, NPY.NPY_OWNDATA);

  if (connect == NULL) :
      np_connect = None
  else :
      if(np_connect_idx is not None):
        dim = <NPY.npy_intp> connect_idx[dn_entity]
      else:
        dim = <NPY.npy_intp> 2 * dn_entity # Face cell
      np_connect = NPY.PyArray_SimpleNewFromData(1,
                                                 &dim,
                                                 PDM_G_NUM_NPY_INT,
                                                 <void *> connect)
  PyArray_ENABLEFLAGS(np_connect, NPY.NPY_OWNDATA);
  return (np_connect_idx, np_connect)

# ------------------------------------------------------------------------
def dmesh_distrib_get(DMesh pydm, PDM_mesh_entities_t entity_type):
  """
  """
  # ************************************************************************
  # > Declaration
  cdef PDM_g_num_t  *distrib
  cdef NPY.npy_intp  dim
  cdef int size
  # ************************************************************************
  size = PDM_dmesh_distrib_get(pydm._dm, entity_type, &distrib)

  if (distrib == NULL) :
      np_distrib = None
  else :
      dim = <NPY.npy_intp> size + 1
      np_distrib = NPY.PyArray_SimpleNewFromData(1,
                                                 &dim,
                                                 PDM_G_NUM_NPY_INT,
                                                 <void *> distrib)
  # PyArray_ENABLEFLAGS(np_distrib, NPY.NPY_OWNDATA);
  return np_distrib

# ------------------------------------------------------------------------
def dmesh_bound_get(DMesh pydm, PDM_bound_type_t bound_type):
  """
  """
  # ************************************************************************
  # > Declaration
  cdef int          *connect_idx
  cdef PDM_g_num_t  *connect
  cdef int           dn_entity
  cdef NPY.npy_intp  dim
  cdef int           n_bnd
  # ************************************************************************
  n_bnd = PDM_dmesh_bound_get(pydm._dm,
                              bound_type,
                              &connect,
                              &connect_idx,
                              PDM_OWNERSHIP_USER)

  if (connect_idx == NULL) :
      np_connect_idx = None
  else :
      dim = <NPY.npy_intp> n_bnd + 1
      np_connect_idx = NPY.PyArray_SimpleNewFromData(1,
                                                     &dim,
                                                     NPY.NPY_INT32,
                                                     <void *> connect_idx)
  PyArray_ENABLEFLAGS(np_connect_idx, NPY.NPY_OWNDATA);

  if (connect == NULL) :
      np_connect = None
  else :
      dim = <NPY.npy_intp> connect_idx[n_bnd]
      np_connect = NPY.PyArray_SimpleNewFromData(1,
                                                 &dim,
                                                 PDM_G_NUM_NPY_INT,
                                                 <void *> connect)
  PyArray_ENABLEFLAGS(np_connect, NPY.NPY_OWNDATA);
  return (np_connect_idx, np_connect)
