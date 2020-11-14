
cdef extern from "pdm_dmesh_nodal.h":
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Wrapping of Ppart Structure
    ctypedef struct PDM_DMesh_nodal_t:
      pass

    ctypedef enum PDM_Mesh_nodal_elt_t:
      PDM_MESH_NODAL_POINT    = 1
      PDM_MESH_NODAL_BAR2     = 2
      PDM_MESH_NODAL_TRIA3    = 3
      PDM_MESH_NODAL_QUAD4    = 4
      PDM_MESH_NODAL_POLY_2D  = 5
      PDM_MESH_NODAL_TETRA4   = 6
      PDM_MESH_NODAL_PYRAMID5 = 7
      PDM_MESH_NODAL_PRISM6   = 8
      PDM_MESH_NODAL_HEXA8    = 9
      PDM_MESH_NODAL_POLY_3D  = 10
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Wrapping of function
    int PDM_DMesh_nodal_create(PDM_MPI_Comm comm,
                               int          mesh_dimension,
                               PDM_g_num_t  n_vtx,
                               PDM_g_num_t  n_cell,
                               PDM_g_num_t  n_face,
                               PDM_g_num_t  n_edge)
    void PDM_DMesh_nodal_free(int handle, int partial)

    void PDM_DMesh_nodal_coord_set(int handle, int n_vtx, double* coords)

    PDM_g_num_t *PDM_DMesh_nodal_distrib_vtx_get(int handle)
    PDM_g_num_t *PDM_DMesh_nodal_distrib_section_get(int handle, int id_section)

    int                  PDM_DMesh_nodal_n_vtx_get(int handle)
    int                  PDM_DMesh_nodal_n_sections_get(int handle)
    int*                 PDM_DMesh_nodal_sections_id_get(int handle)
    PDM_Mesh_nodal_elt_t PDM_DMesh_nodal_section_type_get(int handle, int id_section)
    double*              PDM_DMesh_nodal_vtx_get(int handle)

    int                  PDM_DMesh_nodal_section_add(int handle, PDM_Mesh_nodal_elt_t t_elt)
    void                 PDM_DMesh_nodal_section_std_set(int          handle,
                                                        int          id_section,
                                                        int          n_elmts,
                                                        PDM_g_num_t* connec)

    PDM_g_num_t* PDM_DMesh_nodal_section_std_get(int          handle, int id_section)
    int PDM_DMesh_nodal_section_n_elt_get(int          handle, int id_section)

    void PDM_DMesh_nodal_section_poly2d_set(int handle, int id_section, PDM_l_num_t n_elt,
                                            PDM_l_num_t* connec_idx,
                                            PDM_g_num_t   *connec)

    PDM_g_num_t PDM_DMesh_nodal_total_n_cell_get(int handle)
    PDM_g_num_t PDM_DMesh_nodal_total_n_face_get(int handle)
    PDM_g_num_t PDM_DMesh_nodal_total_n_vtx_get(int handle)
    void PDM_DMesh_nodal_cell_face_compute(int handle)
    int PDM_DMesh_nodal_cell_face_get(int handle, PDM_l_num_t** cell_face_idx, PDM_g_num_t **cell_face)
    int PDM_DMesh_nodal_face_cell_get(int handle, PDM_g_num_t** face_cell)
    int PDM_DMesh_nodal_face_vtx_get(int handle, int** dface_vtx_idx, PDM_g_num_t **dface_vtx)

    PDM_g_num_t* PDM_DMesh_nodal_distrib_cell_get(int handle)
    PDM_g_num_t* PDM_DMesh_nodal_distrib_face_get(int handle)
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
    cdef int idmesh
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
        self.idmesh = PDM_DMesh_nodal_create(PDMC, mesh_dimension, n_vtx, n_cell, n_face, n_edge)
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
        n_vtx = dvtx_coord.shape[0]
        PDM_DMesh_nodal_coord_set(self.idmesh, n_vtx, <double *> dvtx_coord.data)
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
          id_section = PDM_DMesh_nodal_section_add(self.idmesh, <PDM_Mesh_nodal_elt_t> elmts_type[i_elmt])
          PDM_DMesh_nodal_section_std_set(self.idmesh, id_section, n_elemts[i_elmt], <PDM_g_num_t *> connect.data)
        # ::::::::::::::::::::::::::::::::::::::::::::::::::

    # ------------------------------------------------------------------------
    def compute(self):
        """
        """
        # ************************************************************************
        # > Declaration
        # ************************************************************************
        PDM_DMesh_nodal_cell_face_compute(self.idmesh)

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
        n_face  = PDM_DMesh_nodal_total_n_face_get(self.idmesh)

        # > Get array
        dn_face = PDM_DMesh_nodal_face_cell_get(self.idmesh, &face_cell)

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
        n_cell  = PDM_DMesh_nodal_total_n_cell_get(self.idmesh)
        # n_cell  = 0

        # > Get array
        dn_cell = PDM_DMesh_nodal_cell_face_get(self.idmesh,  &cell_face_idx, &cell_face)

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
        n_face  = PDM_DMesh_nodal_total_n_face_get(self.idmesh)

        # > Get array
        dn_face = PDM_DMesh_nodal_face_vtx_get(self.idmesh, &dface_vtx_idx, &dface_vtx)

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
        face_distrib = PDM_DMesh_nodal_distrib_face_get(self.idmesh)

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
      print('PDM_DMesh_nodal_free')
      PDM_DMesh_nodal_free(self.idmesh, 1)



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


