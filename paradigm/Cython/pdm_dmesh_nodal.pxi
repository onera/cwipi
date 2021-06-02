
cdef extern from "pdm_mesh_nodal.h":
  int PDM_Mesh_nodal_n_vertices_element(PDM_Mesh_nodal_elt_t type,
                                        int            order);

cdef extern from "pdm_dmesh_nodal.h":
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Wrapping of Ppart Structure
    ctypedef struct PDM_dmesh_nodal_t:
      pass
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Wrapping of function
    PDM_dmesh_nodal_t* PDM_DMesh_nodal_create(PDM_MPI_Comm comm,
                                              int          mesh_dimension,
                                              PDM_g_num_t  n_vtx,
                                              PDM_g_num_t  n_cell,
                                              PDM_g_num_t  n_face,
                                              PDM_g_num_t  n_edge)
    void PDM_DMesh_nodal_free(PDM_dmesh_nodal_t* dmn, int partial)

    void PDM_DMesh_nodal_coord_set(PDM_dmesh_nodal_t* dmn, int n_vtx, double* coords, PDM_ownership_t    owner)
    void PDM_DMesh_nodal_section_g_dims_get(PDM_dmesh_nodal_t* dmn,
                                            PDM_g_num_t        *n_cell_abs,
                                            PDM_g_num_t        *n_face_abs,
                                            PDM_g_num_t        *n_edge_abs,
                                            PDM_g_num_t        *n_vtx_abs)

    # PDM_g_num_t *PDM_DMesh_nodal_distrib_vtx_get(PDM_dmesh_nodal_t* dmn)
    # PDM_g_num_t *PDM_DMesh_nodal_distrib_section_get(PDM_dmesh_nodal_t* dmn, int id_section)

    int                  PDM_DMesh_nodal_n_vtx_get(PDM_dmesh_nodal_t* dmn)
    int*                 PDM_DMesh_nodal_sections_id_get(PDM_dmesh_nodal_t* dmn)
    int                  PDM_DMesh_nodal_n_section_get(PDM_dmesh_nodal_t* dmn)
    PDM_Mesh_nodal_elt_t PDM_DMesh_nodal_section_type_get(PDM_dmesh_nodal_t* dmn, int id_section)
    double*              PDM_DMesh_nodal_vtx_get(PDM_dmesh_nodal_t* dmn)
    PDM_g_num_t*         PDM_DMesh_nodal_section_distri_std_copy_get(PDM_dmesh_nodal_t *dmesh_nodal,
                                                                     int                id_section)
    PDM_g_num_t*         PDM_dmesh_nodal_vtx_distrib_copy_get(PDM_dmesh_nodal_t *dmesh_nodal)
    int                  PDM_DMesh_nodal_section_add(PDM_dmesh_nodal_t* dmn, PDM_Mesh_nodal_elt_t t_elt)
    void                 PDM_DMesh_nodal_section_std_set(PDM_dmesh_nodal_t* dmn,
                                                        int                 id_section,
                                                        int                 n_elmts,
                                                        PDM_g_num_t*        connec,
                                                        PDM_ownership_t     owner)

    PDM_g_num_t* PDM_DMesh_nodal_section_std_get(PDM_dmesh_nodal_t* dmn, int id_section)
    int PDM_DMesh_nodal_section_n_elt_get(PDM_dmesh_nodal_t* dmn, int id_section)

    void PDM_DMesh_nodal_section_poly2d_set(PDM_dmesh_nodal_t* dmn, int id_section, PDM_l_num_t n_elt,
                                            PDM_l_num_t        *connec_idx,
                                            PDM_g_num_t        *connec,
                                            PDM_ownership_t     owner)
    void PDM_DMesh_nodal_section_group_elmt_set(PDM_dmesh_nodal_t  *dmesh_nodal,
                                                int                 n_group_elmt,
                                                int                *dgroup_elmt_idx,
                                                PDM_g_num_t        *dgroup_elmt,
                                                PDM_ownership_t     owner)
    void PDM_DMesh_nodal_section_group_elmt_get(PDM_dmesh_nodal_t  *dmesh_nodal,
                                                int                *n_group_elmt,
                                                int               **dgroup_elmt_idx,
                                                PDM_g_num_t       **dgroup_elmt)

    PDM_g_num_t PDM_dmesh_nodal_total_n_cell_get(PDM_dmesh_nodal_t* dmn)
    PDM_g_num_t PDM_dmesh_nodal_total_n_face_get(PDM_dmesh_nodal_t* dmn)
    PDM_g_num_t PDM_dmesh_nodal_total_n_vtx_get(PDM_dmesh_nodal_t* dmn)
    void PDM_dmesh_nodal_generate_distribution(PDM_dmesh_nodal_t* dmn)
    void PDM_DMesh_nodal_cell_face_compute(PDM_dmesh_nodal_t* dmn)
    int PDM_DMesh_nodal_cell_face_get(PDM_dmesh_nodal_t* dmn, PDM_l_num_t** cell_face_idx, PDM_g_num_t **cell_face)
    int PDM_DMesh_nodal_face_cell_get(PDM_dmesh_nodal_t* dmn, PDM_g_num_t** face_cell)
    int PDM_DMesh_nodal_face_vtx_get(PDM_dmesh_nodal_t* dmn, int** dface_vtx_idx, PDM_g_num_t **dface_vtx)

    PDM_g_num_t* PDM_DMesh_nodal_distrib_cell_get(PDM_dmesh_nodal_t* dmn)
    PDM_g_num_t* PDM_DMesh_nodal_distrib_face_get(PDM_dmesh_nodal_t* dmn)
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

cdef extern from "pdm_elt_parent_find.h":
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Wrapping of Ppart Structure
    void PDM_elt_parent_find_from_distrib(PDM_g_num_t  *elt_distrib,
                                          int          *elt_def_idx,
                                          PDM_g_num_t  *elt_def,
                                          PDM_g_num_t  *elt_to_find_distrib,
                                          int          *elt_to_find_def_idx,
                                          PDM_g_num_t  *elt_to_find_def,
                                          PDM_MPI_Comm  comm,
                                          PDM_g_num_t  *parent)

    void PDM_elt_parent_find(int           dnelt,
                             int          *elt_def_idx,
                             PDM_g_num_t  *elt_def,
                             int           dnelt_to_find,
                             int          *elt_to_find_def_idx,
                             PDM_g_num_t  *elt_to_find_def,
                             PDM_MPI_Comm  comm,
                             PDM_g_num_t  *parent)

cdef extern from "pdm_distrib.h":

    void PDM_distrib_compute(int           dnelt,
                             PDM_g_num_t  *elt_distrib,
                             int           offset,
                             PDM_MPI_Comm  comm)

# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# ------------------------------------------------------------------
cdef class DistributedMeshNodal:
    """
       DistributedMeshNodal: Interface to build face from Element->Vtx connectivity
    """
    # ************************************************************************
    # > Class attributes
    cdef PDM_dmesh_nodal_t *dmn
    # cdef int idmesh
    cdef int n_rank
    # ************************************************************************
    # ------------------------------------------------------------------------
    def __cinit__(self, MPI.Comm    comm,
                        PDM_g_num_t n_vtx,
                        PDM_g_num_t n_cell,
                        PDM_g_num_t n_face = -1,
                        PDM_g_num_t n_edge = -1,
                        int         mesh_dimension = 3):
        """
        TODOUX
        """
        # ************************************************************************
        # > Declaration
        # cdef int      nElts
        # cdef int      idx
        # # > Numpy array
        # cdef NPY.ndarray[npy_pdm_gnum_t, ndim=1, mode='fortran'] partLNToGN
        # ************************************************************************

        # ::::::::::::::::::::::::::::::::::::::::::::::::::
        self.n_rank = comm.Get_size()
        # ::::::::::::::::::::::::::::::::::::::::::::::::::

        # ::::::::::::::::::::::::::::::::::::::::::::::::::
        # > Convert mpi4py -> PDM_MPI
        cdef MPI.MPI_Comm c_comm = comm.ob_mpi
        cdef PDM_MPI_Comm PDMC   = PDM_MPI_mpi_2_pdm_mpi_comm(<void *> &c_comm)
        # ::::::::::::::::::::::::::::::::::::::::::::::::::

        # ::::::::::::::::::::::::::::::::::::::::::::::::::
        self.dmn = PDM_DMesh_nodal_create(PDMC, mesh_dimension, n_vtx, n_cell, n_face, n_edge)
        # ::::::::::::::::::::::::::::::::::::::::::::::::::

    # ------------------------------------------------------------------------
    def set_coordinnates(self, NPY.ndarray[NPY.double_t  , mode='c', ndim=1] dvtx_coord):
        """
        """
        # ************************************************************************
        # > Declaration
        cdef int n_vtx
        # ************************************************************************

        # ::::::::::::::::::::::::::::::::::::::::::::::::::
        n_vtx = dvtx_coord.shape[0]//3
        PDM_DMesh_nodal_coord_set(self.dmn, n_vtx, <double *> dvtx_coord.data, PDM_OWNERSHIP_USER)
        # ::::::::::::::::::::::::::::::::::::::::::::::::::

    # ------------------------------------------------------------------------
    def set_sections(self, list elmt_list,
                     NPY.ndarray[NPY.int32_t, mode='c', ndim=1] elmts_type,
                     NPY.ndarray[NPY.int32_t, mode='c', ndim=1] n_elemts):
        """
           TODO : Split function as PDM
        """
        # ************************************************************************
        # > Declaration
        cdef int n_vtx
        cdef NPY.ndarray[npy_pdm_gnum_t, ndim=1, mode='fortran'] connect
        # ************************************************************************

        # ::::::::::::::::::::::::::::::::::::::::::::::::::
        # > Panic assert
        assert(len(elmt_list) == elmts_type.shape[0])
        # ::::::::::::::::::::::::::::::::::::::::::::::::::

        # ::::::::::::::::::::::::::::::::::::::::::::::::::
        for i_elmt, connect in enumerate(elmt_list):
          id_section = PDM_DMesh_nodal_section_add(self.dmn, <PDM_Mesh_nodal_elt_t> elmts_type[i_elmt])
          PDM_DMesh_nodal_section_std_set(self.dmn,
                                          id_section,
                                          n_elemts[i_elmt],
                          <PDM_g_num_t *> connect.data,
                                          PDM_OWNERSHIP_USER)
        # ::::::::::::::::::::::::::::::::::::::::::::::::::

    # ------------------------------------------------------------------------
    def set_group_elmt(self, n_group_elmt,
                       NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] dgroup_elmt_idx,
                       NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] dgroup_elmt):
        """
           TODO : Split function as PDM
        """
        if(dgroup_elmt_idx is None):
          PDM_DMesh_nodal_section_group_elmt_set(self.dmn,
                                                 n_group_elmt,
                                                 NULL, NULL, PDM_OWNERSHIP_USER)
        else:
          PDM_DMesh_nodal_section_group_elmt_set(self.dmn,
                                                 n_group_elmt,
                                          <int*> dgroup_elmt_idx.data,
                                  <PDM_g_num_t*> dgroup_elmt.data,
                                                 PDM_OWNERSHIP_USER)

    # ------------------------------------------------------------------------
    def generate_distribution(self):
        """
        """
        # ************************************************************************
        # > Declaration
        # ************************************************************************
        PDM_dmesh_nodal_generate_distribution(self.dmn)

    # ------------------------------------------------------------------------
    def compute(self):
        """
        """
        # ************************************************************************
        # > Declaration
        # ************************************************************************
        PDM_DMesh_nodal_cell_face_compute(self.dmn)

    # ------------------------------------------------------------------------
    def get_face_cell(self):
        """
        """
        # ************************************************************************
        # > Declaration
        cdef PDM_g_num_t *face_cell
        cdef PDM_g_num_t n_face
        cdef PDM_g_num_t dn_face
        cdef NPY.npy_intp dim
        # ************************************************************************

        # > Get Size
        n_face  = PDM_dmesh_nodal_total_n_face_get(self.dmn)

        # > Get array
        dn_face = PDM_DMesh_nodal_face_cell_get(self.dmn, &face_cell)

        # > Build numpy capsule
        dim = <NPY.npy_intp> 2 * dn_face
        np_face_cell = NPY.PyArray_SimpleNewFromData(1,
                                                   &dim,
                                                   PDM_G_NUM_NPY_INT,
                                                   <void *> face_cell)

        return {'n_face'        : n_face,
                'dn_face'       : dn_face,
                'np_dface_cell' : np_face_cell}

    # ------------------------------------------------------------------------
    def get_cell_face(self):
        """
        """
        # ************************************************************************
        # > Declaration
        cdef PDM_g_num_t *cell_face
        cdef PDM_l_num_t *cell_face_idx
        cdef PDM_g_num_t n_face
        cdef PDM_g_num_t dn_face
        cdef NPY.npy_intp dim
        # ************************************************************************

        # > Get Size
        n_cell  = PDM_dmesh_nodal_total_n_cell_get(self.dmn)
        # n_cell  = 0

        # > Get array
        dn_cell = PDM_DMesh_nodal_cell_face_get(self.dmn,  &cell_face_idx, &cell_face)

        # > Build numpy capsule
        dim = <NPY.npy_intp> dn_cell + 1
        np_cell_face_idx = NPY.PyArray_SimpleNewFromData(1,
                                                         &dim,
                                                         NPY.NPY_INT32,
                                                         <void *> cell_face_idx)
        PyArray_ENABLEFLAGS(np_cell_face_idx, NPY.NPY_OWNDATA);

        dim = <NPY.npy_intp> np_cell_face_idx[np_cell_face_idx.shape[0]-1]
        np_cell_face = NPY.PyArray_SimpleNewFromData(1,
                                                     &dim,
                                                     PDM_G_NUM_NPY_INT,
                                                     <void *> cell_face)
        PyArray_ENABLEFLAGS(np_cell_face, NPY.NPY_OWNDATA);

        return {'n_cell'            : n_cell,
                'dn_cell'           : dn_cell,
                'npdcell_face_idx'  : np_cell_face_idx,
                'npdcell_face'      : np_cell_face}

    # ------------------------------------------------------------------------
    def get_face_vtx(self):
        """
        """
        # ************************************************************************
        # > Declaration
        cdef PDM_g_num_t *dface_vtx
        cdef PDM_l_num_t *dface_vtx_idx
        cdef PDM_g_num_t n_face
        cdef PDM_g_num_t dn_face
        cdef NPY.npy_intp dim
        # ************************************************************************

        # > Get Size
        n_face  = PDM_dmesh_nodal_total_n_face_get(self.dmn)

        # > Get array
        dn_face = PDM_DMesh_nodal_face_vtx_get(self.dmn, &dface_vtx_idx, &dface_vtx)

        # > Build numpy capsule
        dim = <NPY.npy_intp> dn_face + 1
        np_dface_vtx_idx = NPY.PyArray_SimpleNewFromData(1,
                                                         &dim,
                                                         NPY.NPY_INT32,
                                                         <void *> dface_vtx_idx)
        PyArray_ENABLEFLAGS(np_dface_vtx_idx, NPY.NPY_OWNDATA);

        dim = <NPY.npy_intp> np_dface_vtx_idx[np_dface_vtx_idx.shape[0]-1]
        np_face_vtx = NPY.PyArray_SimpleNewFromData(1,
                                                    &dim,
                                                    PDM_G_NUM_NPY_INT,
                                                    <void *> dface_vtx)
        PyArray_ENABLEFLAGS(np_face_vtx, NPY.NPY_OWNDATA);

        return {'n_face'           : n_face,
                'np_dface_vtx_idx' : np_dface_vtx_idx,
                'np_dface_vtx'     : np_face_vtx}

    # ------------------------------------------------------------------------
    def get_distrib_face(self):
        """
        """
        # ************************************************************************
        # > Declaration
        cdef PDM_g_num_t *face_distrib
        cdef NPY.npy_intp dim
        # ************************************************************************

        # > Get array
        face_distrib = PDM_DMesh_nodal_distrib_face_get(self.dmn)

        # > Build numpy capsule
        dim = <NPY.npy_intp> self.n_rank + 1
        np_distrib_face = NPY.PyArray_SimpleNewFromData(1,
                                                        &dim,
                                                        PDM_G_NUM_NPY_INT,
                                                        <void *> face_distrib)
        PyArray_ENABLEFLAGS(np_distrib_face, NPY.NPY_OWNDATA);

        return np_distrib_face

    # ------------------------------------------------------------------------
    def __dealloc__(self):
      """
         Use the free method of PDM Lib
      """
      # ************************************************************************
      # > Declaration
      # ************************************************************************
      # > Free Ppart Structure
      # print('PDM_DMesh_nodal_free')
      PDM_DMesh_nodal_free(self.dmn, 1)

# ------------------------------------------------------------------
cdef class DistributedMeshNodalCaspule:
  """
  """
  # ************************************************************************
  # > Class attributes
  cdef PDM_dmesh_nodal_t* dmn
  # ************************************************************************
  # ------------------------------------------------------------------------
  def __cinit__(self, object caps):
    """
    """
    # print("DistributedMeshNodalCaspule", PyCapsule_GetName(caps))
    cdef PDM_dmesh_nodal_t* casp_dmn = <PDM_dmesh_nodal_t *> PyCapsule_GetPointer(caps, NULL)
    self.dmn = casp_dmn;

  # ------------------------------------------------------------------------
  def dmesh_nodal_get_sections(self, MPI.Comm    comm):
    """
    """
    return dmesh_nodal_get_sections(self, comm)

  # ------------------------------------------------------------------------
  def dmesh_nodal_get_group(self):
    """
    """
    return dmesh_nodal_get_group(self)

  # ------------------------------------------------------------------------
  def dmesh_nodal_get_g_dims(self):
    """
    """
    # print("Wrap dmesh_nodal_get_g_dims")
    return dmesh_nodal_get_g_dims(self)

  # ------------------------------------------------------------------------
  def __dealloc__(self):
    """
       Use the free method of PDM Lib
    """
    # print("DistributedMeshNodalCaspule::__dealloc__")
    PDM_DMesh_nodal_free(self.dmn, 1)
    # print("DistributedMeshNodalCaspule::__dealloc__ end z")

ctypedef fused DMeshNodal:
  DistributedMeshNodal
  DistributedMeshNodalCaspule


def dmesh_nodal_get_g_dims(DMeshNodal pydmn):
  """
  """
  # ************************************************************************
  # > Declaration
  cdef PDM_g_num_t n_cell_abs, n_face_abs, n_edge_abs, n_vtx_abs
  # ************************************************************************

  PDM_DMesh_nodal_section_g_dims_get(pydmn.dmn, &n_cell_abs, &n_face_abs, &n_edge_abs, &n_vtx_abs)
  return {"n_cell_abs" : n_cell_abs,
          "n_face_abs" : n_face_abs,
          "n_edge_abs" : n_edge_abs,
          "n_vtx_abs"  : n_vtx_abs}


def dmesh_nodal_get_sections(DMeshNodal pydmn, MPI.Comm    comm):
  """
  """
  # ************************************************************************
  # > Declaration
  cdef int                   n_section
  cdef int                   n_vtx_per_elmt
  cdef PDM_g_num_t           dn_elmt
  cdef int                  *section_id
  cdef int                  *connect_idx
  cdef PDM_g_num_t          *connect
  cdef PDM_g_num_t          *section_distrib
  cdef double               *vtx_coord
  cdef PDM_Mesh_nodal_elt_t  t_elmt
  cdef NPY.npy_intp          dim
  # ************************************************************************

  n_section  = PDM_DMesh_nodal_n_section_get(pydmn.dmn)
  section_id = PDM_DMesh_nodal_sections_id_get(pydmn.dmn)

  # print("n_section : ", n_section)

  sections = []
  for i_section in range(n_section):
    id_section = section_id[i_section]
    t_elmt = PDM_DMesh_nodal_section_type_get(pydmn.dmn, id_section)
    # > For now only use this interfaces for standard element ... (to be refactor with HO and polyhedra elements)
    assert(t_elmt != PDM_MESH_NODAL_POLY_2D)
    assert(t_elmt != PDM_MESH_NODAL_POLY_3D)

    section_distrib = PDM_DMesh_nodal_section_distri_std_copy_get(pydmn.dmn, id_section)
    connect         = PDM_DMesh_nodal_section_std_get(pydmn.dmn, id_section)

    # > Build numpy capsule
    dim = <NPY.npy_intp> comm.Get_size() + 1
    np_distrib = NPY.PyArray_SimpleNewFromData(1,
                                                    &dim,
                                                    PDM_G_NUM_NPY_INT,
                                                    <void *> section_distrib)
    PyArray_ENABLEFLAGS(np_distrib, NPY.NPY_OWNDATA);

    # > Build numpy capsule
    dn_elmt = np_distrib[comm.Get_rank()+1] - np_distrib[comm.Get_rank()]
    n_vtx_per_elmt = PDM_Mesh_nodal_n_vertices_element(t_elmt, 1)
    dim = <NPY.npy_intp> n_vtx_per_elmt * dn_elmt
    np_connec = NPY.PyArray_SimpleNewFromData(1,
                                               &dim,
                                               PDM_G_NUM_NPY_INT,
                                               <void *> connect)
    PyArray_ENABLEFLAGS(np_connec, NPY.NPY_OWNDATA);

    sections.append({"pdm_type"   : t_elmt,
                     "np_distrib" : np_distrib,
                     "np_connec"  : np_connec})
    # print("t_elmt : ", t_elmt )
    # print("np_distrib : ", np_distrib )
    # print("np_connec : ", np_connec )

  vtx_distrib = PDM_dmesh_nodal_vtx_distrib_copy_get(pydmn.dmn)
  n_vtx = PDM_DMesh_nodal_n_vtx_get(pydmn.dmn);
  vtx_coord = PDM_DMesh_nodal_vtx_get(pydmn.dmn)

  dim = <NPY.npy_intp> 3 * n_vtx
  np_vtx = NPY.PyArray_SimpleNewFromData(1,
                                         &dim,
                                         NPY.NPY_DOUBLE,
                                         <void *> vtx_coord)
  PyArray_ENABLEFLAGS(np_vtx, NPY.NPY_OWNDATA);

  dim = <NPY.npy_intp> comm.Get_size() + 1
  np_vtx_distrib = NPY.PyArray_SimpleNewFromData(1,
                                             &dim,
                                             PDM_G_NUM_NPY_INT,
                                             <void *> vtx_distrib)
  PyArray_ENABLEFLAGS(np_vtx_distrib, NPY.NPY_OWNDATA);

  # print("Return dmesh_nodal_get_sections")
  return {"vtx"      : {"np_vtx" : np_vtx, "np_vtx_distrib" : np_vtx_distrib},
          "sections" : sections}


def dmesh_nodal_get_group(DMeshNodal pydmn):
  """
  """
  # ************************************************************************
  # > Declaration
  cdef int                   n_group
  cdef int                  *dgroup_elmt_idx
  cdef PDM_g_num_t          *dgroup_elmt
  cdef NPY.npy_intp          dim
  # ************************************************************************

  PDM_DMesh_nodal_section_group_elmt_get(pydmn.dmn, &n_group, &dgroup_elmt_idx, &dgroup_elmt);

  dim = <NPY.npy_intp> n_group + 1
  np_dgroup_elmt_idx = NPY.PyArray_SimpleNewFromData(1,
                                                     &dim,
                                                     NPY.NPY_INT32,
                                                     <void *> dgroup_elmt_idx)
  PyArray_ENABLEFLAGS(np_dgroup_elmt_idx, NPY.NPY_OWNDATA);

  dim = <NPY.npy_intp> np_dgroup_elmt_idx[n_group]
  np_dgroup_elmt = NPY.PyArray_SimpleNewFromData(1,
                                                 &dim,
                                                 PDM_G_NUM_NPY_INT,
                                                 <void *> dgroup_elmt)
  PyArray_ENABLEFLAGS(np_dgroup_elmt, NPY.NPY_OWNDATA);

  return {"dgroup_elmt_idx" : np_dgroup_elmt_idx,
          "dgroup_elmt"     : np_dgroup_elmt}


# ------------------------------------------------------------------------
def ElementParentFind(int                                           dnelt,
                      NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] elt_def_idx,
                      NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] elt_def,
                      int                                           dnelt_to_find,
                      NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] elt_to_find_def_idx,
                      NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] elt_to_find_def,
                      MPI.Comm    comm,
                      NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] parent):
    """
    """
    # ************************************************************************
    # > Declaration
    # ************************************************************************

    # ::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Convert mpi4py -> PDM_MPI
    cdef MPI.MPI_Comm c_comm = comm.ob_mpi
    cdef PDM_MPI_Comm PDMC   = PDM_MPI_mpi_2_pdm_mpi_comm(<void *> &c_comm)
    # ::::::::::::::::::::::::::::::::::::::::::::::::::

    # ::::::::::::::::::::::::::::::::::::::::::::::::::
    PDM_elt_parent_find(dnelt,
                        <int *>         elt_def_idx.data,
                        <PDM_g_num_t *> elt_def.data,
                        dnelt_to_find,
                        <int *>         elt_to_find_def_idx.data,
                        <PDM_g_num_t *> elt_to_find_def.data,
                        PDMC,
                        <PDM_g_num_t *> parent.data)
    # ::::::::::::::::::::::::::::::::::::::::::::::::::


# ------------------------------------------------------------------------
def ComputeDistributionFromDelmt(int         dnelt,
                                 MPI.Comm    comm,
                                 int         offset=0):
    """
    """
    # ************************************************************************
    # > Declaration
    cdef NPY.ndarray[npy_pdm_gnum_t, ndim=1, mode='c'] elt_distrib = NPY.empty( comm.Get_size() + 1, dtype=npy_pdm_gnum_dtype, order='C')
    # ************************************************************************

    # ::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Convert mpi4py -> PDM_MPI
    cdef MPI.MPI_Comm c_comm = comm.ob_mpi
    cdef PDM_MPI_Comm PDMC   = PDM_MPI_mpi_2_pdm_mpi_comm(<void *> &c_comm)
    # ::::::::::::::::::::::::::::::::::::::::::::::::::

    # ::::::::::::::::::::::::::::::::::::::::::::::::::
    PDM_distrib_compute(dnelt,
                        <PDM_g_num_t *> elt_distrib.data,
                        offset,
                        PDMC)
    # ::::::::::::::::::::::::::::::::::::::::::::::::::

    return elt_distrib

# ------------------------------------------------------------------------


