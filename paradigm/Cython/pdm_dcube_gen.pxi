cdef extern from "pdm_dcube_gen.h":
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Wrapping of Ppart Structure
    ctypedef struct PDM_dcube_t:
      pass
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    # ------------------------------------------------------------------
    PDM_dcube_t* PDM_dcube_gen_init(PDM_MPI_Comm        comm,
                                    const PDM_g_num_t   n_vtx_seg,
                                    const double        length,
                                    const double        zero_x,
                                    const double        zero_y,
                                    const double        zero_z,
                                    PDM_ownership_t     owner)

    # ------------------------------------------------------------------
    void PDM_dcube_gen_dim_get(PDM_dcube_t        *pdm_dcube,
                               int                *n_face_group,
                               int                *dn_cell,
                               int                *dn_face,
                               int                *dn_vtx,
                               int                *sface_vtx,
                               int                *sface_group)

    # ------------------------------------------------------------------
    void PDM_dcube_gen_data_get(PDM_dcube_t        *pdm_dcube,
                                PDM_g_num_t       **dface_cell,
                                int               **dface_vtx_idx,
                                PDM_g_num_t       **dface_vtx,
                                double            **dvtx_coord,
                                int               **dface_group_idx,
                                PDM_g_num_t       **dface_group)

    # ------------------------------------------------------------------
    void PDM_dcube_gen_free(PDM_dcube_t        *pdm_dcube)
    # ------------------------------------------------------------------

# ------------------------------------------------------------------
cdef class DCubeGenerator:
    """
       DCubeGenerator
    """
    # > For Ppart
    cdef PDM_dcube_t* _dcube
    # ------------------------------------------------------------------
    def __cinit__(self,
                  npy_pdm_gnum_t                         n_vtx_seg,
                  NPY.double_t                           length,
                  NPY.double_t                           zero_x,
                  NPY.double_t                           zero_y,
                  NPY.double_t                           zero_z,
                  MPI.Comm                               comm):

        """
        """
        # ~> Communicator Mpi
        cdef MPI.MPI_Comm c_comm = comm.ob_mpi

        # -> Create dcube
        self._dcube = PDM_dcube_gen_init(PDM_MPI_mpi_2_pdm_mpi_comm (<void *> &c_comm),
                                         n_vtx_seg,
                                         length,
                                         zero_x,
                                         zero_y,
                                         zero_z,
                                         PDM_OWNERSHIP_USER) # Python take owership

    # ------------------------------------------------------------------
    def __dealloc__(self):
        PDM_dcube_gen_free(self._dcube)

    # ------------------------------------------------------------------
    def dcube_dim_get(self):
        """
           Get dcube dimensions
        """
        # ************************************************************************
        # > Declaration
        cdef int n_face_group
        cdef int dn_cell
        cdef int dn_face
        cdef int dn_vtx
        cdef int sface_vtx
        cdef int sface_group
        # ************************************************************************

        PDM_dcube_gen_dim_get(self._dcube,
                              &n_face_group,
                              &dn_cell,
                              &dn_face,
                              &dn_vtx,
                              &sface_vtx,
                              &sface_group)

        return {'n_face_group' : n_face_group,
                'dn_cell'      : dn_cell,
                'dn_face'      : dn_face,
                'dn_vtx'       : dn_vtx,
                'sface_vtx'    : sface_vtx,
                'sface_group'  : sface_group}

    # ------------------------------------------------------------------
    def dcube_val_get(self):
        """
           Get dcube data
        """
        # ************************************************************************
        # > Declaration
        cdef NPY.npy_intp dim
        cdef PDM_g_num_t  *dface_cell
        cdef int          *dface_vtx_idx
        cdef PDM_g_num_t  *dface_vtx
        cdef double       *dvtx_coord
        cdef int          *dface_group_idx
        cdef PDM_g_num_t  *dface_group
        # ************************************************************************

        dims = self.dcube_dim_get()

        PDM_dcube_gen_data_get(self._dcube,
                               &dface_cell,
                               &dface_vtx_idx,
                               &dface_vtx,
                               &dvtx_coord,
                               &dface_group_idx,
                               &dface_group)

        # \param [out]  dface_cell            Face to cell connectivity (size = 2*nFace)
        dim = <NPY.npy_intp> 2*dims['dn_face']
        np_dface_cell = NPY.PyArray_SimpleNewFromData(1,
                                                   &dim,
                                                   PDM_G_NUM_NPY_INT,
                                                   <void *> dface_cell)
        PyArray_ENABLEFLAGS(np_dface_cell, NPY.NPY_OWNDATA);

        # \param [out]  dface_vtx_idx        Face to vtx connectivity indexes (size = nFace+1)
        dim = <NPY.npy_intp> (dims['dn_face'] + 1)
        np_dface_vtx_idx = NPY.PyArray_SimpleNewFromData(1,
                                                     &dim,
                                                     NPY.NPY_INT32,
                                                     <void *> dface_vtx_idx)
        PyArray_ENABLEFLAGS(np_dface_vtx_idx, NPY.NPY_OWNDATA);

        # \param [out]  dface_vtx            Face to vtx connectivity (size = dface_vtx_idx[nFace])
        dim = <NPY.npy_intp> dims['sface_vtx']
        np_dface_vtx = NPY.PyArray_SimpleNewFromData(1,
                                                  &dim,
                                                  PDM_G_NUM_NPY_INT,
                                                  <void *> dface_vtx)
        PyArray_ENABLEFLAGS(np_dface_vtx, NPY.NPY_OWNDATA);

        # \param [out]  dvtx_coords          Vertices coordinates (size = 3*nVtx)
        dim = <NPY.npy_intp> 3*dims['dn_vtx']
        np_dtx_coord = NPY.PyArray_SimpleNewFromData(1,
                                                   &dim,
                                                   NPY.NPY_DOUBLE,
                                                   <void *> dvtx_coord)
        PyArray_ENABLEFLAGS(np_dtx_coord, NPY.NPY_OWNDATA);

        # \param [out]  dface_group_idx       Face group indexes (size = n_face_group + 1)
        dim = <NPY.npy_intp> (dims['n_face_group'] + 1)
        np_dface_group_idx = NPY.PyArray_SimpleNewFromData(1,
                                                       &dim,
                                                       NPY.NPY_INT32,
                                                       <void *> dface_group_idx)
        PyArray_ENABLEFLAGS(np_dface_group_idx, NPY.NPY_OWNDATA);

        # \param [out]  dface_group          Face group (size = dface_group_idx[n_face_group])
        dim = <NPY.npy_intp> dims['sface_group']
        np_dface_group = NPY.PyArray_SimpleNewFromData(1,
                                                    &dim,
                                                    PDM_G_NUM_NPY_INT,
                                                    <void *> dface_group)
        PyArray_ENABLEFLAGS(np_dface_group, NPY.NPY_OWNDATA);

        return {'dface_cell'      : np_dface_cell,
                'dface_vtx_idx'   : np_dface_vtx_idx,
                'dface_vtx'       : np_dface_vtx,
                'dvtx_coord'      : np_dtx_coord,
                'dface_group_idx' : np_dface_group_idx,
                'dface_group'     : np_dface_group}
