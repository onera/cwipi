
cdef extern from "pdm_part_mesh_nodal.h":
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Wrapping of Ppart Structure
    ctypedef struct PDM_part_mesh_nodal_t:
      pass
    ctypedef struct PDM_part_mesh_nodal_elmts_t:
      pass
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    PDM_part_mesh_nodal_t* PDM_part_mesh_nodal_create(int          mesh_dimension,
                                                      int          n_part,
                                                      PDM_MPI_Comm comm)

    void PDM_part_mesh_nodal_coord_set(       PDM_part_mesh_nodal_t *pmn,
                                       int                    id_part,
                                       int                    n_vtx,
                                       double                *coords,
                                       PDM_g_num_t           *numabs,
                                       PDM_ownership_t        owner)

    int PDM_part_mesh_nodal_n_part_get(PDM_part_mesh_nodal_t *pmn)

    int PDM_part_mesh_nodal_n_vtx_get(PDM_part_mesh_nodal_t *pmn,
                                      int                    id_part)

    double* PDM_part_mesh_nodal_vtx_coord_get(PDM_part_mesh_nodal_t *pmn,
                                              int                    id_part)


    PDM_g_num_t* PDM_part_mesh_nodal_vtx_g_num_get(PDM_part_mesh_nodal_t *pmn,
                                                   int                    id_part)

    int PDM_part_mesh_nodal_n_section_get(PDM_part_mesh_nodal_t *pmn,
                                          PDM_geometry_kind_t    geom_kind)

    int * PDM_part_mesh_nodal_sections_id_get(PDM_part_mesh_nodal_t *pmn,
                                              PDM_geometry_kind_t    geom_kind)


    PDM_Mesh_nodal_elt_t PDM_part_mesh_nodal_section_elt_type_get(PDM_part_mesh_nodal_t *pmn,
                                                                  PDM_geometry_kind_t    geom_kind,
                                                                  int                    id_section)

    PDM_Mesh_nodal_elt_t PDM_part_mesh_nodal_section_type_get(PDM_part_mesh_nodal_t *pmn,
                                                              PDM_geometry_kind_t    geom_kind,
                                                              int                    id_section)

    int PDM_part_mesh_nodal_section_add(PDM_part_mesh_nodal_t *pmn,
                                        PDM_geometry_kind_t    geom_kind,
                                        PDM_Mesh_nodal_elt_t   t_elt)

    void PDM_part_mesh_nodal_section_std_set(PDM_part_mesh_nodal_t *pmn,
                                             PDM_geometry_kind_t    geom_kind,
                                             int                    id_block,
                                             int                    id_part,
                                             int                    n_elt,
                                             int                   *connec,
                                             PDM_g_num_t           *numabs,
                                             int                   *parent_num,
                                             PDM_g_num_t           *parent_entity_g_num,
                                             PDM_ownership_t        owner)

    int PDM_part_mesh_nodal_block_n_elt_get(PDM_part_mesh_nodal_t  *pmn,
                                            PDM_geometry_kind_t     geom_kind,
                                            int                     id_block,
                                            int                     id_part)


    void PDM_part_mesh_nodal_block_std_get(PDM_part_mesh_nodal_t  *pmn,
                                           PDM_geometry_kind_t     geom_kind,
                                           int                     id_block,
                                           int                     id_part,
                                           int                   **connec,
                                           PDM_g_num_t           **numabs,
                                           int                   **parent_num,
                                           PDM_g_num_t           **parent_entity_g_num)

    PDM_Mesh_nodal_elt_t PDM_part_mesh_nodal_block_type_get(PDM_part_mesh_nodal_t  *pmn,
                                                            PDM_geometry_kind_t     geom_kind,
                                                            int                     id_block)

    void PDM_part_mesh_nodal_add_part_mesh_nodal_elmts( PDM_part_mesh_nodal_t       *pmn,
                                                       PDM_part_mesh_nodal_elmts_t *pmne,
                                                       PDM_ownership_t              owner)

    void PDM_part_mesh_nodal_free( PDM_part_mesh_nodal_t* pmn)

# ------------------------------------------------------------------
cdef class PartMeshNodal:
    """
       PartMeshNodal: Interface to build face from Element->Vtx connectivity
    """
    # ************************************************************************
    # > Class attributes
    cdef PDM_part_mesh_nodal_t *pmn
    # cdef int idmesh
    cdef int n_rank
    # ************************************************************************
    # ------------------------------------------------------------------------
    def __cinit__(self, MPI.Comm    comm,
                        int         n_part,
                        int         mesh_dimension = 3):
        """
        """
        # ::::::::::::::::::::::::::::::::::::::::::::::::::
        self.n_rank = comm.Get_size()
        # ::::::::::::::::::::::::::::::::::::::::::::::::::

        # ::::::::::::::::::::::::::::::::::::::::::::::::::
        # > Convert mpi4py -> PDM_MPI
        cdef MPI.MPI_Comm c_comm = comm.ob_mpi
        cdef PDM_MPI_Comm PDMC   = PDM_MPI_mpi_2_pdm_mpi_comm(<void *> &c_comm)
        # ::::::::::::::::::::::::::::::::::::::::::::::::::

        # ::::::::::::::::::::::::::::::::::::::::::::::::::
        self.pmn = PDM_part_mesh_nodal_create(mesh_dimension, n_part, PDMC)
        # ::::::::::::::::::::::::::::::::::::::::::::::::::

    # ------------------------------------------------------------------------
    def __dealloc__(self):
      """
      """
      PDM_part_mesh_nodal_free(self.pmn)


# ------------------------------------------------------------------
cdef class PartMeshNodalCaspule:
  """
  """
  # ************************************************************************
  # > Class attributes
  cdef PDM_part_mesh_nodal_t* pmn
  # ************************************************************************
  # ------------------------------------------------------------------------
  def __cinit__(self, object caps):
    """
    """
    # print("DistributedMeshNodalCaspule", PyCapsule_GetName(caps))
    cdef PDM_part_mesh_nodal_t* casp_pmn = <PDM_part_mesh_nodal_t *> PyCapsule_GetPointer(caps, NULL)
    self.pmn = casp_pmn;

  # ------------------------------------------------------------------------
  def part_mesh_nodal_get_sections(self, PDM_geometry_kind_t geom_kind, int i_part):
    """
    """
    return part_mesh_nodal_get_sections(self, geom_kind, i_part)



  # ------------------------------------------------------------------------
  def __dealloc__(self):
    """
    """
    PDM_part_mesh_nodal_free(self.pmn)



ctypedef fused PMeshNodal:
  PartMeshNodal
  PartMeshNodalCaspule

def part_mesh_nodal_get_sections(PMeshNodal pypmn, PDM_geometry_kind_t geom_kind, int i_part):
  """
  """
  # ************************************************************************
  # > Declaration
  cdef int                   n_section
  cdef int                   n_vtx_per_elmt
  cdef int                   n_elmt_in_section
  cdef int                  *section_id
  cdef int                  *parent_num
  cdef int                  *connec
  cdef PDM_g_num_t          *numabs
  cdef PDM_g_num_t          *parent_entity_g_num
  cdef double               *vtx_coord
  cdef PDM_Mesh_nodal_elt_t  t_elmt
  cdef NPY.npy_intp          dim
  # ************************************************************************

  n_section  = PDM_part_mesh_nodal_n_section_get  (pypmn.pmn, geom_kind)
  section_id = PDM_part_mesh_nodal_sections_id_get(pypmn.pmn, geom_kind)

  sections = []
  for i_section in range(n_section):
    id_section = section_id[i_section]

    t_elmt = PDM_part_mesh_nodal_block_type_get(pypmn.pmn, geom_kind, id_section)
    assert(t_elmt != PDM_MESH_NODAL_POLY_2D)
    assert(t_elmt != PDM_MESH_NODAL_POLY_3D)

    n_elmt_in_section = PDM_part_mesh_nodal_block_n_elt_get(pypmn.pmn, geom_kind, id_section, i_part)

    PDM_part_mesh_nodal_block_std_get(pypmn.pmn, geom_kind, id_section, i_part, &connec, &numabs, &parent_num, &parent_entity_g_num)

    n_vtx_per_elmt = PDM_Mesh_nodal_n_vertices_element(t_elmt, 1)

    dim = <NPY.npy_intp> n_elmt_in_section * n_vtx_per_elmt
    np_connec = NPY.PyArray_SimpleNewFromData(1,
                                              &dim,
                                              NPY.NPY_INT32,
                                     <void *> connec)
    PyArray_ENABLEFLAGS(np_connec, NPY.NPY_OWNDATA);

    dim = <NPY.npy_intp> n_elmt_in_section
    np_parent_num = NPY.PyArray_SimpleNewFromData(1,
                                              &dim,
                                              NPY.NPY_INT32,
                                     <void *> parent_num)
    PyArray_ENABLEFLAGS(np_parent_num, NPY.NPY_OWNDATA);

    dim = <NPY.npy_intp> n_elmt_in_section
    np_numabs = NPY.PyArray_SimpleNewFromData(1,
                                              &dim,
                                              PDM_G_NUM_NPY_INT,
                                     <void *> numabs)
    PyArray_ENABLEFLAGS(np_numabs, NPY.NPY_OWNDATA);

    dim = <NPY.npy_intp> n_elmt_in_section
    np_parent_entity_g_num = None
    if(parent_entity_g_num != NULL):
      np_parent_entity_g_num = NPY.PyArray_SimpleNewFromData(1,
                                                             &dim,
                                                             PDM_G_NUM_NPY_INT,
                                                             <void *> parent_entity_g_num)
      PyArray_ENABLEFLAGS(np_parent_entity_g_num, NPY.NPY_OWNDATA);

    sections.append({"pdm_type"               : t_elmt,
                     "np_connec"              : np_connec,
                     "np_parent_num"          : np_parent_num,
                     "np_numabs"              : np_numabs,
                     "np_parent_entity_g_num" : np_parent_entity_g_num})

  return sections
