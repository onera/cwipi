cdef extern from "pdm_dcube_gen.h":
    # ------------------------------------------------------------------
    void PDM_dcube_gen_init(int                *id,
                            PDM_MPI_Comm        comm,
                            const PDM_g_num_t   n_vtx_seg,
                            const double        length,
                            const double        zero_x,
                            const double        zero_y,
                            const double        zero_z)

    # ------------------------------------------------------------------
    void PDM_dcube_gen_dim_get(int                id,
                               int                *n_face_group,
                               int                *dn_cell,
                               int                *dn_face,
                               int                *dn_vtx,
                               int                *dface_vtxL,
                               int                *dFacegroupL)

    # ------------------------------------------------------------------
    void PDM_dcube_gen_data_get(int                 id,
                                PDM_g_num_t       **dface_cell,
                                int               **dface_vtx_idx,
                                PDM_g_num_t       **dface_vtx,
                                double            **dvtx_coord,
                                int               **dface_group_idx,
                                PDM_g_num_t       **dface_group)

    # ------------------------------------------------------------------
    void PDM_dcube_gen_free(int  id)
    # ------------------------------------------------------------------

# ------------------------------------------------------------------
cdef class DCubeGenerator:
    """
       DCubeGenerator
    """
    # > For Ppart
    cdef int _dcube_id
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
        PDM_dcube_gen_init(&self._dcube_id,
                           PDM_MPI_mpi_2_pdm_mpi_comm (<void *> &c_comm),
                           n_vtx_seg,
                           length,
                           zero_x,
                           zero_y,
                           zero_z)


    # ------------------------------------------------------------------
    def __dealloc__(self):
        PDM_dcube_gen_free(self._dcube_id)

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
        cdef int dface_vtxL
        cdef int dFacegroupL
        # ************************************************************************

        PDM_dcube_gen_dim_get(self._dcube_id,
                              &n_face_group,
                              &dn_cell,
                              &dn_face,
                              &dn_vtx,
                              &dface_vtxL,
                              &dFacegroupL)

        return {'nFaceGroup'    : n_face_group,
                'dnCell'        : dn_cell,
                'dnFace'        : dn_face,
                'dnVtx'         : dn_vtx,
                'dFaceVtxL'     : dface_vtxL,
                'dFaceGroupL'   : dFacegroupL}

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

        PDM_dcube_gen_data_get(self._dcube_id,
                               &dface_cell,
                               &dface_vtx_idx,
                               &dface_vtx,
                               &dvtx_coord,
                               &dface_group_idx,
                               &dface_group)

        # \param [out]  dface_cell            Face to cell connectivity (size = 2*nFace)
        dim = <NPY.npy_intp> 2*dims['dnFace']
        npFaceCell = NPY.PyArray_SimpleNewFromData(1,
                                                   &dim,
                                                   PDM_G_NUM_NPY_INT,
                                                   <void *> dface_cell)

        # \param [out]  dface_vtx_idx        Face to vtx connectivity indexes (size = nFace+1)
        dim = <NPY.npy_intp> (dims['dnFace'] + 1)
        npFaceVtxIdx = NPY.PyArray_SimpleNewFromData(1,
                                                     &dim,
                                                     NPY.NPY_INT32,
                                                     <void *> dface_vtx_idx)

        # \param [out]  dface_vtx            Face to vtx connectivity (size = dface_vtx_idx[nFace])
        dim = <NPY.npy_intp> dims['dFaceVtxL']
        npFaceVtx = NPY.PyArray_SimpleNewFromData(1,
                                                  &dim,
                                                  PDM_G_NUM_NPY_INT,
                                                  <void *> dface_vtx)

        # \param [out]  dvtx_coords          Vertices coordinates (size = 3*nVtx)
        dim = <NPY.npy_intp> 3*dims['dnVtx']
        npVtxCoord = NPY.PyArray_SimpleNewFromData(1,
                                                   &dim,
                                                   NPY.NPY_DOUBLE,
                                                   <void *> dvtx_coord)

        # \param [out]  dface_group_idx       Face group indexes (size = nFaceGroup + 1)
        dim = <NPY.npy_intp> (dims['nFaceGroup'] + 1)
        npFaceGroupIdx = NPY.PyArray_SimpleNewFromData(1,
                                                       &dim,
                                                       NPY.NPY_INT32,
                                                       <void *> dface_group_idx)

        # \param [out]  dface_group          Face group (size = dface_group_idx[nFaceGroup])
        dim = <NPY.npy_intp> dims['dFaceGroupL']
        npFaceGroup = NPY.PyArray_SimpleNewFromData(1,
                                                    &dim,
                                                    PDM_G_NUM_NPY_INT,
                                                    <void *> dface_group)

        return {'dFaceCell'              : npFaceCell,
                'dFaceVtxIdx'            : npFaceVtxIdx,
                'dFaceVtx'               : npFaceVtx,
                'dVtxCoord'              : npVtxCoord,
                'dFaceGroupIdx'          : npFaceGroupIdx,
                'dFaceGroup'             : npFaceGroup}
