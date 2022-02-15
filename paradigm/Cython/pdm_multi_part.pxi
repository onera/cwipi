cdef extern from "pdm_para_graph_dual.h":
    ctypedef enum PDM_split_dual_t:
        PDM_SPLIT_DUAL_WITH_PARMETIS = 1
        PDM_SPLIT_DUAL_WITH_PTSCOTCH = 2
cdef extern from "pdm_multipart.h":
    ctypedef enum PDM_part_size_t:
        PDM_PART_SIZE_HOMONEGEOUS   = 1
        PDM_PART_SIZE_HETEROGENEOUS = 2
    ctypedef struct PDM_multipart_t:
        pass

    # -> PPART bases functions
    # ------------------------------------------------------------------
    # MPI_Comm      comm,
    PDM_multipart_t* PDM_multipart_create(int              n_zone,
                                          int*             n_part,
                                          PDM_bool_t       merge_blocks,
                                          PDM_split_dual_t split_method,
                                          PDM_part_size_t  part_size_method,
                                          double*          part_fraction,
                                          PDM_MPI_Comm     comm,
                                          PDM_ownership_t  owner)

    # ------------------------------------------------------------------
    void PDM_multipart_register_block(PDM_multipart_t *mtp,
                                      int              zone_gid,
                                      PDM_dmesh_t     *dmesh)

    # ------------------------------------------------------------------
    void PDM_multipart_register_dmesh_nodal(PDM_multipart_t   *mtp,
                                            int                zone_gid,
                                            PDM_dmesh_nodal_t *dmesh)

    # ------------------------------------------------------------------
    void PDM_multipart_register_joins(PDM_multipart_t *mtp,
                                      int              n_total_joins,
                                      int*             matching_join_array)

    # ------------------------------------------------------------------
    void PDM_multipart_set_reordering_options(      PDM_multipart_t *mtp,
                                              const int              i_zone,
                                              const char            *renum_cell_method,
                                              const int             *renum_cell_properties,
                                              const char            *renum_face_method)

    # ------------------------------------------------------------------
    void PDM_multipart_run_ppart(PDM_multipart_t *mtp)

    # ------------------------------------------------------------------
    void PDM_multipart_part_dim_get(PDM_multipart_t  *mtp,
                                    int               zone_gid,
                                    int               ipart,
                                    int              *n_section,
                                    int             **n_elt,
                                    int              *n_cell,
                                    int              *n_face,
                                    int              *n_face_part_bound,
                                    int              *nVtx,
                                    int              *n_proc,
                                    int              *nt_part,
                                    int              *scell_face,
                                    int              *sface_vtx,
                                    int              *s_face_bound,
                                    int              *n_face_bound,
                                    int              *s_face_join,
                                    int              *n_face_join)

    # ------------------------------------------------------------------
    void PDM_multipart_part_val_get(PDM_multipart_t   *mtp,
                                    int                zone_gid,
                                    int                ipart,
                                    int             ***elt_vtx_idx,
                                    int             ***elt_vtx,
                                    PDM_g_num_t     ***elt_section_ln_to_gn,
                                    int              **cell_tag,
                                    int              **cell_face_idx,
                                    int              **cell_face,
                                    PDM_g_num_t      **cell_ln_to_gn,
                                    int              **face_tag,
                                    int              **face_cell,
                                    int              **face_vtx_idx,
                                    int              **face_vtx,
                                    PDM_g_num_t      **face_ln_to_gn,
                                    int              **face_part_bound_proc_idx,
                                    int              **face_part_bound_part_idx,
                                    int              **face_part_bound,
                                    int              **vtx_tag,
                                    double           **vtx,
                                    PDM_g_num_t      **vtx_ln_to_gn,
                                    int              **face_bound_idx,
                                    int              **face_bound,
                                    PDM_g_num_t      **face_bound_ln_to_gn,
                                    int              **face_join_idx,
                                    int              **face_join,
                                    PDM_g_num_t      **face_join_ln_to_gn)

    # ------------------------------------------------------------------
    void PDM_multipart_get_part_mesh_nodal(PDM_multipart_t   *mtp,
                                           const int   i_zone,
                                           PDM_part_mesh_nodal_t **pmesh_nodal)

    # ------------------------------------------------------------------
    void PDM_multipart_part_graph_comm_vtx_dim_get(PDM_multipart_t   *mtp,
                                                   int                i_zone,
                                                   int                i_part,
                                                   int               *n_vtx_part_bound)

    # ------------------------------------------------------------------
    void PDM_multipart_part_graph_comm_vtx_data_get(PDM_multipart_t   *mtp,
                                                    int                i_zone,
                                                    int                i_part,
                                                    int              **vtx_part_bound_proc_idx,
                                                    int              **vtx_part_bound_part_idx,
                                                    int              **vtx_part_bound)

    # ------------------------------------------------------------------
    void PDM_multipart_part_color_get(PDM_multipart_t   *mtp,
                                      int                zone_gid,
                                      int                ipart,
                                      int              **cell_color,
                                      int              **face_color,
                                      int              **face_hp_color,
                                      int              **thread_color,
                                      int              **hyper_plane_color)

    # ------------------------------------------------------------------
    void PDM_multipart_part_ghost_infomation_get(PDM_multipart_t *mtp,
                                                 int              zone_gid,
                                                 int              ipart,
                                                 int             **vtx_ghost_information)

    # ------------------------------------------------------------------
    int PDM_multipart_part_connectivity_get(PDM_multipart_t          *mtp,
                                            int                       i_zone,
                                            int                       i_part,
                                            PDM_connectivity_type_t   connectivity_type,
                                            int                     **connect,
                                            int                     **connect_idx,
                                            PDM_ownership_t           ownership);
    # ------------------------------------------------------------------
    int PDM_multipart_part_ln_to_gn_get(PDM_multipart_t      *mtp,
                                        int                   i_zone,
                                        int                   i_part,
                                        PDM_mesh_entities_t   entity_type,
                                        PDM_g_num_t         **entity_ln_to_gn,
                                        PDM_ownership_t       ownership);
    # ------------------------------------------------------------------
    int PDM_multipart_partition_color_get(PDM_multipart_t      *mtp,
                                          int                   i_zone,
                                          int                   i_part,
                                          PDM_mesh_entities_t   entity_type,
                                          int                 **entity_ln_to_gn,
                                          PDM_ownership_t       ownership);

    # ------------------------------------------------------------------
    void PDM_multipart_time_get(PDM_multipart_t      *mtp,
                                int                   zone_gid,
                                double              **elapsed,
                                double              **cpu,
                                double              **cpu_user,
                                double              **cpu_sys)
    # ------------------------------------------------------------------
    void PDM_multipart_free(PDM_multipart_t      *mtp)

# ------------------------------------------------------------------
cdef class MultiPart:
    """
       MultiPpart
    """
    # > For Ppart
    cdef PDM_multipart_t* _mtp
    # ------------------------------------------------------------------
    def __cinit__(self,
                  int                                           n_zone,
                  NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] n_part,
                  int                                           merge_blocks,
                  PDM_split_dual_t                              split_method,
                  PDM_part_size_t                               part_size_method,
                  NPY.ndarray[NPY.double_t  , mode='c', ndim=1] part_fraction,
                  MPI.Comm                                      comm):

        """
        """
        # ~> Communicator Mpi
        cdef MPI.MPI_Comm c_comm = comm.ob_mpi
        cdef double* part_fraction_data

        # print("MultiPart::n_zone -->", n_zone)
        # print("MultiPart::n_part -->", n_part)
        # print("MultiPart::merge_blocks -->", merge_blocks)
        # print("MultiPart::split_method -->", split_method)

        if part_fraction is None:
          part_fraction_data = NULL
        else:
          part_fraction_data = <double *> part_fraction.data

        # -> Create PPART
        self._mtp = PDM_multipart_create(n_zone,
                                              <int*> n_part.data,
                                              <PDM_bool_t> merge_blocks,
                                              split_method,
                                              part_size_method,
                                              part_fraction_data,
                                              PDM_MPI_mpi_2_pdm_mpi_comm (<void *> &c_comm),
                                              PDM_OWNERSHIP_USER) # Python take ownership

    # ------------------------------------------------------------------
    def __dealloc__(self):
        PDM_multipart_free(self._mtp)

    # ------------------------------------------------------------------
    def multipart_register_block(self, int zone_gid,
                                       DMesh dm): # DMesh = DistributedMeshCaspule or DistributedMesh
        """
        """
        PDM_multipart_register_block(self._mtp,
                                     zone_gid,
                                     dm._dm)

    # ------------------------------------------------------------------
    def multipart_register_dmesh_nodal(self, int zone_gid,
                                       DistributedMeshNodal dmn): # DMesh = DistributedMeshCaspule or DistributedMesh
        """
        """
        PDM_multipart_register_dmesh_nodal(self._mtp,
                                           zone_gid,
                                           dmn.dmn)

    # ------------------------------------------------------------------
    def multipart_register_joins(self, int n_total_joins,
                                 NPY.ndarray[NPY.int32_t, mode='c', ndim=1] matching_join):
        """
        """
        PDM_multipart_register_joins(       self._mtp,
                                            n_total_joins,
                                     <int*> matching_join.data)
    # ------------------------------------------------------------------
    def multipart_set_reordering(self, int i_zone,
                                 char *renum_cell_method,
                                 char *renum_face_method,
                                 NPY.ndarray[NPY.int32_t, mode='c', ndim=1] renum_properties_cell):
      """
      """
      cdef int *renum_properties_cell_data
      if (renum_properties_cell is None):
        renum_properties_cell_data = NULL
      else:
        renum_properties_cell_data = <int *> renum_properties_cell.data
      PDM_multipart_set_reordering_options(self._mtp,
                                           i_zone,
                                           renum_cell_method,
                                           renum_properties_cell_data,
                                           renum_face_method)
    # ------------------------------------------------------------------
    def multipart_run_ppart(self):
        """
        """
        PDM_multipart_run_ppart(self._mtp)

    # ------------------------------------------------------------------
    def multipart_dim_get(self, int ipart, int zone_gid):
        """
           Get partition dimensions
        """
        # ************************************************************************
        # > Declaration
        cdef int n_cell
        cdef int n_proc
        cdef int nt_part
        cdef int n_face
        cdef int n_face_part_bound
        cdef int n_vtx
        cdef int scell_face
        cdef int s_face_vtx
        cdef int s_face_bound
        cdef int n_face_bound
        cdef int s_face_join
        cdef int n_face_join

        cdef int *n_elt
        cdef int n_section
        # ************************************************************************

        PDM_multipart_part_dim_get(self._mtp,
                                   zone_gid,
                                   ipart,
                                   &n_section,
                                   &n_elt,
                                   &n_cell,
                                   &n_face,
                                   &n_face_part_bound,
                                   &n_vtx,
                                   &n_proc,
                                   &nt_part,
                                   &scell_face,
                                   &s_face_vtx,
                                   &s_face_bound,
                                   &n_face_bound,
                                   &s_face_join,
                                   &n_face_join)

        cdef NPY.npy_intp dim
        if (n_elt == NULL) :
            np_n_elt = None
        else :
            dim = n_section
            np_n_elt = NPY.PyArray_SimpleNewFromData(1,
                                                     &dim,
                                                     NPY.NPY_INT32,
                                                     <void *> n_elt)
            #PyArray_ENABLEFLAGS(np_n_elt, NPY.NPY_OWNDATA); # well it should be there if PDM_multipart_part_dim_get were to be called once

        return {'n_cell'            : n_cell,
                'n_section'         : n_section,
                'n_elt'             : np_n_elt,
                'ipart'             : ipart,
                'n_face'            : n_face,
                'nt_part'           : nt_part,
                'n_proc'            : n_proc,
                'n_face_part_bound' : n_face_part_bound,
                'n_vtx'             : n_vtx,
                'scell_face'        : scell_face,
                's_face_vtx'        : s_face_vtx,
                's_face_bound'      : s_face_bound,
                'n_face_bound'      : n_face_bound,
                's_face_join'       : s_face_join,
                'n_face_join'       : n_face_join}

    # ------------------------------------------------------------------
    def multipart_val_get(self, int ipart, int zone_gid):
        """
           Get partition dimensions
        """
        # ************************************************************************
        # > Declaration
        cdef int          *cell_tag,
        cdef int          *cell_face_idx,
        cdef int          *cell_face,
        cdef PDM_g_num_t  *cell_ln_to_gn,
        cdef int          *face_tag,
        cdef int          *face_cell,
        cdef int          *face_vtx_idx,
        cdef int          *face_vtx,
        cdef PDM_g_num_t  *face_ln_to_gn,
        cdef int          *face_part_bound,
        cdef int          *face_part_bound_proc_idx,
        cdef int          *face_part_bound_part_idx,
        cdef int          *vtx_tag,
        cdef double       *vtx_coord,
        cdef PDM_g_num_t  *vtx_ln_to_gn,
        cdef int          *face_bound_idx,
        cdef int          *face_bound,
        cdef PDM_g_num_t  *face_bound_ln_to_gn,
        cdef int          *face_join_idx,
        cdef int          *face_join,
        cdef PDM_g_num_t  *face_join_ln_to_gn
        cdef int         **elt_vtx_idx
        cdef int         **elt_vtx
        cdef PDM_g_num_t **elt_section_ln_to_gn
        # ************************************************************************

        # dims = self.part_dim_get(self._mtp, ipart)
        dims = self.multipart_dim_get(ipart, zone_gid)

        # -> Call PPART to get info
        PDM_multipart_part_val_get(self._mtp,
                                   zone_gid,
                                   ipart,
                                   &elt_vtx_idx,
                                   &elt_vtx,
                                   &elt_section_ln_to_gn,
                                   &cell_tag,
                                   &cell_face_idx,
                                   &cell_face,
                                   &cell_ln_to_gn,
                                   &face_tag,
                                   &face_cell,
                                   &face_vtx_idx,
                                   &face_vtx,
                                   &face_ln_to_gn,
                                   &face_part_bound_proc_idx,
                                   &face_part_bound_part_idx,
                                   &face_part_bound,
                                   &vtx_tag,
                                   &vtx_coord,
                                   &vtx_ln_to_gn,
                                   &face_bound_idx,
                                   &face_bound,
                                   &face_bound_ln_to_gn,
                                   &face_join_idx,
                                   &face_join,
                                   &face_join_ln_to_gn)
        # elt sections
        cdef NPY.npy_intp n_section = <NPY.npy_intp> dims['n_section']
        cdef list np_elt_vtx_idx = []
        cdef list np_elt_vtx = []
        cdef list np_elt_section_ln_to_gn = []
        cdef NPY.npy_intp dim_idx
        for i_section in range(n_section):
            dim_idx = dims['n_elt'][i_section] + 1
            np_elt_vtx_idx_i_section = NPY.PyArray_SimpleNewFromData(1,
                                                     &dim_idx,
                                                     NPY.NPY_INT32,
                                                     <void *> elt_vtx_idx[i_section])
            n_elt_i_section = <NPY.npy_intp> dims['n_elt'][i_section]
            np_elt_section_ln_to_gn_i_section = NPY.PyArray_SimpleNewFromData(1,
                                                     &n_elt_i_section,
                                                     PDM_G_NUM_NPY_INT,
                                                     <void *> elt_section_ln_to_gn[i_section])
            dim_elt = <NPY.npy_intp> np_elt_vtx_idx_i_section[-1]
            np_elt_vtx_i_section = NPY.PyArray_SimpleNewFromData(1,
                                                     &dim_elt,
                                                     NPY.NPY_INT32,
                                                     <void *> elt_vtx[i_section])

            PyArray_ENABLEFLAGS(np_elt_vtx_idx_i_section, NPY.NPY_OWNDATA);
            PyArray_ENABLEFLAGS(np_elt_vtx_i_section, NPY.NPY_OWNDATA);
            PyArray_ENABLEFLAGS(np_elt_section_ln_to_gn_i_section, NPY.NPY_OWNDATA);

            np_elt_vtx_idx.append(np_elt_vtx_idx_i_section)
            np_elt_vtx.append(np_elt_vtx_i_section)
            np_elt_section_ln_to_gn.append(np_elt_section_ln_to_gn_i_section)

        # -> Begin
        cdef NPY.npy_intp dim

        # \param [out]  cell_tag            Cell tag (size = n_cell)
        if (cell_tag == NULL) :
            np_cell_tag = None
        else :
            dim = <NPY.npy_intp> dims['n_cell']
            np_cell_tag = NPY.PyArray_SimpleNewFromData(1,
                                                     &dim,
                                                     NPY.NPY_INT32,
                                                     <void *> cell_tag)
            PyArray_ENABLEFLAGS(np_cell_tag, NPY.NPY_OWNDATA);

        # \param [out]  cell_face_idx        Cell to face connectivity index (size = n_cell + 1)
        if (cell_face_idx == NULL) :
            np_cell_face_idx = None
        else :
            dim = <NPY.npy_intp> (dims['n_cell'] + 1)
            np_cell_face_idx = NPY.PyArray_SimpleNewFromData(1,
                                                         &dim,
                                                         NPY.NPY_INT32,
                                                         <void *> cell_face_idx)
            PyArray_ENABLEFLAGS(np_cell_face_idx, NPY.NPY_OWNDATA);

        # \param [out]  cell_face           Cell to face connectivity (size = cell_face_idx[n_cell] = lcell_face)
        if (cell_face == NULL) :
            np_cell_face = None
        else :
            dim = <NPY.npy_intp> dims['scell_face']
            np_cell_face = NPY.PyArray_SimpleNewFromData(1,
                                                      &dim,
                                                      NPY.NPY_INT32,
                                                      <void *> cell_face)
            PyArray_ENABLEFLAGS(np_cell_face, NPY.NPY_OWNDATA);

        # \param [out]  cell_ln_to_gn         Cell local numbering to global numbering (size = n_cell)
        # dim = <NPY.npy_intp> dims['n_cell']
        if (cell_ln_to_gn == NULL) :
            np_cell_ln_to_gn = None
        else :
            dim = <NPY.npy_intp> dims['n_cell']
            np_cell_ln_to_gn = NPY.PyArray_SimpleNewFromData(1,
                                                        &dim,
                                                        PDM_G_NUM_NPY_INT,
                                                        <void *> cell_ln_to_gn)
            PyArray_ENABLEFLAGS(np_cell_ln_to_gn, NPY.NPY_OWNDATA);

        # \param [out]  face_tag            Face tag (size = n_face)
        if (face_tag == NULL) :
            np_face_tag = None
        else :
            dim = <NPY.npy_intp> dims['n_face']
            np_face_tag = NPY.PyArray_SimpleNewFromData(1,
                                                     &dim,
                                                     NPY.NPY_INT32,
                                                     <void *> face_tag)
            PyArray_ENABLEFLAGS(np_face_tag, NPY.NPY_OWNDATA);

        # \param [out]  face_cell           Face to cell connectivity  (size = 2 * n_face)
        if (face_cell == NULL) :
            np_face_cell = None
        else :
            dim = <NPY.npy_intp> (2 * dims['n_face'])
            np_face_cell = NPY.PyArray_SimpleNewFromData(1,
                                                      &dim,
                                                      NPY.NPY_INT32,
                                                      <void *> face_cell)
            PyArray_ENABLEFLAGS(np_face_cell, NPY.NPY_OWNDATA);

        # \param [out]  face_vtx_idx         Face to vtx_coord connectivity index (size = n_face + 1)
        if (face_vtx_idx == NULL) :
            np_face_vtx_idx = None
        else :
            dim = <NPY.npy_intp> (dims['n_face'] + 1)
            np_face_vtx_idx = NPY.PyArray_SimpleNewFromData(1,
                                                           &dim,
                                                           NPY.NPY_INT32,
                                                           <void *> face_vtx_idx)
            PyArray_ENABLEFLAGS(np_face_vtx_idx, NPY.NPY_OWNDATA);

        # \param [out]  face_vtx            Face to vtx_coord connectivity (size = face_vtx_idx[n_face])
        cdef NPY.ndarray[NPY.int32_t, ndim=1] np_face_vtx
        if (face_vtx == NULL) :
            np_face_vtx = None
        else :
            dim = <NPY.npy_intp> dims['s_face_vtx']
            np_face_vtx  = NPY.PyArray_SimpleNewFromData(1,
                                                         &dim,
                                                         NPY.NPY_INT32,
                                                         <void *> face_vtx)
            PyArray_ENABLEFLAGS(np_face_vtx, NPY.NPY_OWNDATA);

        # \param [out]  face_ln_to_gn         Face local numbering to global numbering (size = n_face)
        if (face_ln_to_gn == NULL) :
            np_face_ln_to_gn = None
        else :
            dim = <NPY.npy_intp> dims['n_face']
            np_face_ln_to_gn   = NPY.PyArray_SimpleNewFromData(1,
                                                          &dim,
                                                          PDM_G_NUM_NPY_INT,
                                                          <void *> face_ln_to_gn)
            PyArray_ENABLEFLAGS(np_face_ln_to_gn, NPY.NPY_OWNDATA);

        # \param [out]  face_part_bound      Partitioning boundary faces
        if (face_part_bound == NULL) :
            np_face_part_bound = None
        else :
            dim = <NPY.npy_intp> (4 * dims['n_face_part_bound'])
            np_face_part_bound   = NPY.PyArray_SimpleNewFromData(1,
                                                             &dim,
                                                             NPY.NPY_INT32,
                                                             <void *> face_part_bound)
            PyArray_ENABLEFLAGS(np_face_part_bound, NPY.NPY_OWNDATA);

        # \param [out]  face_part_bound_proc_idx  Partitioning boundary faces block distribution from processus (size = n_proc + 1)
        if (face_part_bound_proc_idx == NULL) :
            np_face_part_bound_proc_idx = None
        else :
            dim = <NPY.npy_intp> ( dims['n_proc'] + 1)
            np_face_part_bound_proc_idx   = NPY.PyArray_SimpleNewFromData(1,
                                                             &dim,
                                                             NPY.NPY_INT32,
                                                             <void *> face_part_bound_proc_idx)
            PyArray_ENABLEFLAGS(np_face_part_bound_proc_idx, NPY.NPY_OWNDATA);

        # \param [out]  face_part_bound_part_idx  Partitioning boundary faces block distribution from partition (size = nt_part + 1)
        if (face_part_bound_part_idx == NULL) :
            np_face_part_bound_part_idx = None
        else :
            dim = <NPY.npy_intp> ( dims['nt_part'] + 1)
            np_face_part_bound_part_idx   = NPY.PyArray_SimpleNewFromData(1,
                                                             &dim,
                                                             NPY.NPY_INT32,
                                                             <void *> face_part_bound_part_idx)
            PyArray_ENABLEFLAGS(np_face_part_bound_part_idx, NPY.NPY_OWNDATA);

        # \param [out]  vtx_tag             vtx_coord tag (size = nVtx)
        if (vtx_tag == NULL) :
            np_vtx_tag = None
        else :
            dim = <NPY.npy_intp> dims['n_vtx']
            np_vtx_tag   = NPY.PyArray_SimpleNewFromData(1,
                                                         &dim,
                                                         NPY.NPY_INT32,
                                                         <void *> vtx_tag)
            PyArray_ENABLEFLAGS(np_vtx_tag, NPY.NPY_OWNDATA);

        # \param [out]  vtx                vtx_coord coordinates (size = 3 * nVtx)
        if (vtx_coord == NULL) :
            np_vtx_coord = None
        else :
            dim = <NPY.npy_intp> (3 * dims['n_vtx'])
            np_vtx_coord  = NPY.PyArray_SimpleNewFromData(1,
                                                     &dim,
                                                     NPY.NPY_DOUBLE,
                                                     <void *> vtx_coord)
            PyArray_ENABLEFLAGS(np_vtx_coord, NPY.NPY_OWNDATA);

        # \param [out]  vtx_ln_to_gn          vtx_coord local numbering to global numbering (size = nVtx)
        if (vtx_ln_to_gn == NULL) :
            np_vtx_ln_to_gn = None
        else :
            dim = <NPY.npy_intp> dims['n_vtx']
            np_vtx_ln_to_gn  = NPY.PyArray_SimpleNewFromData(1,
                                                           &dim,
                                                           PDM_G_NUM_NPY_INT,
                                                           <void *> vtx_ln_to_gn)
            PyArray_ENABLEFLAGS(np_vtx_ln_to_gn, NPY.NPY_OWNDATA);

        # \param [out]  face_bound_idx       face group index (size = n_face_bound + 1)
        if (face_bound_idx == NULL) :
            np_face_bound_idx = None
        else :
            dim = <NPY.npy_intp> dims['n_face_bound'] + 1
            np_face_bound_idx  = NPY.PyArray_SimpleNewFromData(1,
                                                           &dim,
                                                           NPY.NPY_INT32,
                                                           <void *> face_bound_idx)
            PyArray_ENABLEFLAGS(np_face_bound_idx, NPY.NPY_OWNDATA);

        # \param [out]  face_bound          faces for each group (size = face_bound_idx[n_face_bound] = lFace_bound)
        if (face_bound == NULL) :
            np_face_bound = None
        else :
            dim = <NPY.npy_intp> dims['s_face_bound']
            np_face_bound = NPY.PyArray_SimpleNewFromData(1,
                                                       &dim,
                                                       NPY.NPY_INT32,
                                                       <void *> face_bound)
            PyArray_ENABLEFLAGS(np_face_bound, NPY.NPY_OWNDATA);

        # \param [out]  face_bound_ln_to_gn    faces global numbering for each group (size = face_bound_idx[n_face_bound] = lFace_bound)
        if (face_bound_ln_to_gn == NULL) :
            np_face_bound_ln_to_gn = None
        else :
            dim = <NPY.npy_intp> dims['s_face_bound']
            np_face_bound_ln_to_gn = NPY.PyArray_SimpleNewFromData(1,
                                                             &dim,
                                                             PDM_G_NUM_NPY_INT,
                                                             <void *> face_bound_ln_to_gn)
            PyArray_ENABLEFLAGS(np_face_bound_ln_to_gn, NPY.NPY_OWNDATA);

        # \param [out]  face_join_idx       face group index (size = n_face_join + 1)
        if (face_join_idx == NULL) :
            np_face_join_idx = None
        else :
            dim = <NPY.npy_intp> dims['n_face_join'] + 1
            np_face_join_idx  = NPY.PyArray_SimpleNewFromData(1,
                                                          &dim,
                                                          NPY.NPY_INT32,
                                                          <void *> face_join_idx)
            PyArray_ENABLEFLAGS(np_face_join_idx, NPY.NPY_OWNDATA);

        # \param [out]  face_join          faces for each group (size = face_join_idx[n_face_join] = lFace_join)
        if (face_join == NULL) :
            np_face_join = None
        else :
            dim = <NPY.npy_intp> (4 * dims['s_face_join'])
            np_face_join = NPY.PyArray_SimpleNewFromData(1,
                                                      &dim,
                                                      NPY.NPY_INT32,
                                                      <void *> face_join)
            PyArray_ENABLEFLAGS(np_face_join, NPY.NPY_OWNDATA);

        # \param [out]  face_join_ln_to_gn    faces global numbering for each group (size = face_join_idx[n_face_join] = lFace_join)
        if (face_join_ln_to_gn == NULL) :
            np_face_join_ln_to_gn = None
        else :
            dim = <NPY.npy_intp> dims['s_face_join']
            np_face_join_ln_to_gn = NPY.PyArray_SimpleNewFromData(1,
                                                             &dim,
                                                             PDM_G_NUM_NPY_INT,
                                                             <void *> face_join_ln_to_gn)
            PyArray_ENABLEFLAGS(np_face_join_ln_to_gn, NPY.NPY_OWNDATA);

        return {'np_cell_tag'                  : np_cell_tag,
                'np_cell_face_idx'             : np_cell_face_idx,
                'np_cell_face'                 : np_cell_face,
                'np_cell_ln_to_gn'             : np_cell_ln_to_gn,
                'np_face_tag'                  : np_face_tag,
                'np_face_cell'                 : np_face_cell,
                'np_elt_vtx_idx'               : np_elt_vtx_idx,
                'np_elt_vtx'                   : np_elt_vtx,
                'np_elt_section_ln_to_gn'      : np_elt_section_ln_to_gn,
                'np_face_vtx_idx'              : np_face_vtx_idx,
                'np_face_vtx'                  : np_face_vtx,
                'np_face_ln_to_gn'             : np_face_ln_to_gn,
                'np_face_part_bound_proc_idx'  : np_face_part_bound_proc_idx,
                'np_face_part_bound_part_idx'  : np_face_part_bound_part_idx,
                'np_face_part_bound'           : np_face_part_bound,
                'np_vtx_tag'                   : np_vtx_tag,
                'np_vtx_coord'                 : np_vtx_coord,
                'np_vtx_ln_to_gn'              : np_vtx_ln_to_gn,
                'np_face_bound_idx'            : np_face_bound_idx,
                'np_face_bound'                : np_face_bound,
                'np_face_bound_ln_to_gn'       : np_face_bound_ln_to_gn,
                'np_face_join_idx'             : np_face_join_idx,
                'np_face_join'                 : np_face_join,
                'np_face_join_ln_to_gn'        : np_face_join_ln_to_gn}


    # ------------------------------------------------------------------
    def multipart_part_mesh_nodal_get(self, int zone_gid):
        cdef PDM_part_mesh_nodal_t *pmesh_nodal
        PDM_multipart_get_part_mesh_nodal(self._mtp, zone_gid, &pmesh_nodal)
        if pmesh_nodal == NULL:
          return None
        else:
          #See pdm_part_mesh_nodal.pxi
          py_caps = PyCapsule_New(pmesh_nodal, NULL, NULL);
          return PartMeshNodalCaspule(py_caps)

    def multipart_graph_comm_vtx_dim_get(self, int ipart, int zone_gid):
        """
           Get partition dimensions
        """
        # ************************************************************************
        # > Declaration
        cdef int n_vtx_part_bound
        # ************************************************************************

        PDM_multipart_part_graph_comm_vtx_dim_get(self._mtp,
                                                  zone_gid,
                                                  ipart,
                                                  &n_vtx_part_bound)

        return {'n_vtx_part_bound' : n_vtx_part_bound}

    # ------------------------------------------------------------------
    def multipart_graph_comm_vtx_val_get(self, int ipart, int zone_gid):
        """
           Get partition dimensions
        """
        # ************************************************************************
        # > Declaration
        cdef int          *vtx_part_bound
        cdef int          *vtx_part_bound_proc_idx
        cdef int          *vtx_part_bound_part_idx
        # ************************************************************************

        # dims = self.part_dim_get(self._mtp, ipart)
        dims    = self.multipart_dim_get(ipart, zone_gid)
        dims_gc = self.multipart_graph_comm_vtx_dim_get(ipart, zone_gid)

        # -> Call PPART to get info
        PDM_multipart_part_graph_comm_vtx_data_get(self._mtp,
                                                   zone_gid,
                                                   ipart,
                                                   &vtx_part_bound_proc_idx,
                                                   &vtx_part_bound_part_idx,
                                                   &vtx_part_bound)
        # -> Begin
        cdef NPY.npy_intp dim

        # \param [out]  vtx_part_bound      Partitioning boundary vtxs
        if (vtx_part_bound == NULL) :
            np_vtx_part_bound = None
        else :
            dim = <NPY.npy_intp> (4 * dims_gc['n_vtx_part_bound'])
            np_vtx_part_bound   = NPY.PyArray_SimpleNewFromData(1,
                                                                &dim,
                                                                NPY.NPY_INT32,
                                                                <void *> vtx_part_bound)
            PyArray_ENABLEFLAGS(np_vtx_part_bound, NPY.NPY_OWNDATA);

        # \param [out]  vtx_part_bound_proc_idx  Partitioning boundary vtxs block distribution from processus (size = n_proc + 1)
        if (vtx_part_bound_proc_idx == NULL) :
            np_vtx_part_bound_proc_idx = None
        else :
            dim = <NPY.npy_intp> ( dims['n_proc'] + 1)
            np_vtx_part_bound_proc_idx = NPY.PyArray_SimpleNewFromData(1,
                                                                       &dim,
                                                                       NPY.NPY_INT32,
                                                                       <void *> vtx_part_bound_proc_idx)
            PyArray_ENABLEFLAGS(np_vtx_part_bound_proc_idx, NPY.NPY_OWNDATA);

        # \param [out]  vtx_part_bound_part_idx  Partitioning boundary vtxs block distribution from partition (size = nt_part + 1)
        if (vtx_part_bound_part_idx == NULL) :
            np_vtx_part_bound_part_idx = None
        else :
            dim = <NPY.npy_intp> ( dims['nt_part'] + 1)
            np_vtx_part_bound_part_idx = NPY.PyArray_SimpleNewFromData(1,
                                                                       &dim,
                                                                       NPY.NPY_INT32,
                                                                       <void *> vtx_part_bound_part_idx)
            PyArray_ENABLEFLAGS(np_vtx_part_bound_part_idx, NPY.NPY_OWNDATA);

        return {'np_vtx_part_bound_proc_idx'  : np_vtx_part_bound_proc_idx,
                'np_vtx_part_bound_part_idx'  : np_vtx_part_bound_part_idx,
                'np_vtx_part_bound'           : np_vtx_part_bound}

    # ------------------------------------------------------------------
    def multipart_color_get(self, int ipart, int zone_gid):
        """
           Get partition dimensions
        """
        # ************************************************************************
        # > Declaration
        cdef int          *cell_color,
        cdef int          *face_color
        cdef int          *face_hp_color
        cdef int          *thread_color
        cdef int          *hyper_plane_color
        # ************************************************************************

        # dims = self.part_dim_get(self._mtp, ipart)
        dims = self.multipart_dim_get(ipart, zone_gid)

        # -> Call PPART to get info
        PDM_multipart_part_color_get(self._mtp,
                                     zone_gid,
                                     ipart,
                                     &cell_color,
                                     &face_color,
                                     &face_hp_color,
                                     &thread_color,
                                     &hyper_plane_color)
        # -> Begin
        cdef NPY.npy_intp dim

        # \param [out]  cell_color            Cell tag (size = n_cell)
        if (cell_color == NULL):
            np_cell_color = None
        else :
            dim = <NPY.npy_intp> dims['n_cell']
            np_cell_color = NPY.PyArray_SimpleNewFromData(1,
                                                        &dim,
                                                        NPY.NPY_INT32,
                                                        <void *> cell_color)
            PyArray_ENABLEFLAGS(np_cell_color, NPY.NPY_OWNDATA);

        # \param [out]  face_color            Cell tag (size = n_face)
        if (face_color == NULL):
            np_face_color = None
        else :
            dim = <NPY.npy_intp> dims['n_face']
            np_face_color = NPY.PyArray_SimpleNewFromData(1,
                                                        &dim,
                                                        NPY.NPY_INT32,
                                                        <void *> face_color)
            PyArray_ENABLEFLAGS(np_face_color, NPY.NPY_OWNDATA);

        # \param [out]  face_hp_color            Cell tag (size = n_face)
        if (face_hp_color == NULL):
            np_face_hp_color = None
        else :
            dim = <NPY.npy_intp> dims['n_face']
            np_face_hp_color = NPY.PyArray_SimpleNewFromData(1,
                                                        &dim,
                                                        NPY.NPY_INT32,
                                                        <void *> face_hp_color)
            PyArray_ENABLEFLAGS(np_face_hp_color, NPY.NPY_OWNDATA);

        # \param [out]  thread_color            Cell tag (size = n_cell)
        if (thread_color == NULL):
            np_thread_color = None
        else :
            dim = <NPY.npy_intp> dims['n_cell']
            np_thread_color = NPY.PyArray_SimpleNewFromData(1,
                                                          &dim,
                                                          NPY.NPY_INT32,
                                                          <void *> thread_color)
            PyArray_ENABLEFLAGS(np_thread_color, NPY.NPY_OWNDATA);

        # \param [out]  hyper_plane_color            Cell tag (size = n_cell)
        if (hyper_plane_color == NULL):
            np_hyper_plane_color = None
        else :
            dim = <NPY.npy_intp> dims['n_cell']
            np_hyper_plane_color = NPY.PyArray_SimpleNewFromData(1,
                                                              &dim,
                                                              NPY.NPY_INT32,
                                                              <void *> hyper_plane_color)
            PyArray_ENABLEFLAGS(np_hyper_plane_color, NPY.NPY_OWNDATA);

        return {'np_cell_color'        : np_cell_color,
                'np_face_color'        : np_face_color,
                'np_face_hp_color'     : np_face_hp_color,
                'np_thread_color'      : np_thread_color,
                'np_hyper_plane_color' : np_hyper_plane_color}


    # ------------------------------------------------------------------
    def multipart_ghost_information_get(self, int ipart, int zone_gid):
        """
           Get partition ghost information
        """
        # ************************************************************************
        # > Declaration
        cdef int          *vtx_ghost_information
        # ************************************************************************

        # dims = self.part_dim_get(self._mtp, ipart)
        dims = self.multipart_dim_get(ipart, zone_gid)

        # -> Call PPART to get info
        PDM_multipart_part_ghost_infomation_get(self._mtp,
                                                zone_gid,
                                                ipart,
                                                &vtx_ghost_information)
        # -> Begin
        cdef NPY.npy_intp dim

        # \param [out]  cell_color            Cell tag (size = n_cell)
        if (vtx_ghost_information == NULL):
            np_vtx_ghost_information = None
        else :
            dim = <NPY.npy_intp> dims['n_vtx']
            np_vtx_ghost_information = NPY.PyArray_SimpleNewFromData(1,
                                                                     &dim,
                                                                     NPY.NPY_INT32,
                                                                     <void *> vtx_ghost_information)
            PyArray_ENABLEFLAGS(np_vtx_ghost_information, NPY.NPY_OWNDATA);
        return {'np_vtx_ghost_information' : np_vtx_ghost_information}

    # ------------------------------------------------------------------
    def multipart_connectivity_get(self, int ipart, int zone_gid, PDM_connectivity_type_t connectivity_type):
        """
           Get partition ghost information
        """
        # ************************************************************************
        # > Declaration
        cdef int          *entity1_entity2
        cdef int          *entity1_entity2_idx
        # ************************************************************************

        # -> Call PPART to get info
        n_entity1 = PDM_multipart_part_connectivity_get(self._mtp,
                                                        zone_gid,
                                                        ipart,
                                                        connectivity_type,
                                                        &entity1_entity2,
                                                        &entity1_entity2_idx,
                                                        PDM_OWNERSHIP_USER)
        cdef NPY.npy_intp dim

        if (entity1_entity2_idx == NULL):
            np_entity1_entity2_idx = None
            dim = <NPY.npy_intp> 2 * n_entity1
            np_entity1_entity2 = NPY.PyArray_SimpleNewFromData(1,
                                                               &dim,
                                                               NPY.NPY_INT32,
                                                      <void *> entity1_entity2)
            PyArray_ENABLEFLAGS(np_entity1_entity2, NPY.NPY_OWNDATA);
        else :
            dim = <NPY.npy_intp> n_entity1+1
            np_entity1_entity2_idx = NPY.PyArray_SimpleNewFromData(1,
                                                                   &dim,
                                                                   NPY.NPY_INT32,
                                                                   <void *> entity1_entity2_idx)
            PyArray_ENABLEFLAGS(np_entity1_entity2_idx, NPY.NPY_OWNDATA);

            dim = <NPY.npy_intp> np_entity1_entity2_idx[n_entity1]
            np_entity1_entity2 = NPY.PyArray_SimpleNewFromData(1,
                                                               &dim,
                                                               NPY.NPY_INT32,
                                                      <void *> entity1_entity2)
            PyArray_ENABLEFLAGS(np_entity1_entity2, NPY.NPY_OWNDATA);


        return {'np_entity1_entity2'     : np_entity1_entity2,
                'np_entity1_entity2_idx' : np_entity1_entity2_idx}


    # ------------------------------------------------------------------
    def multipart_ln_to_gn_get(self, int ipart, int zone_gid, PDM_mesh_entities_t entity_type):
        """
           Get partition ghost information
        """
        # ************************************************************************
        # > Declaration
        cdef PDM_g_num_t  *entity_ln_to_gn
        # ************************************************************************

        # -> Call PPART to get info
        n_entity1 = PDM_multipart_part_ln_to_gn_get(self._mtp,
                                                    zone_gid,
                                                    ipart,
                                                    entity_type,
                                                    &entity_ln_to_gn,
                                                    PDM_OWNERSHIP_USER)
        cdef NPY.npy_intp dim
        if (entity_ln_to_gn == NULL) :
            np_entity_ln_to_gn = None
        else :
            dim = <NPY.npy_intp> n_entity1
            np_entity_ln_to_gn   = NPY.PyArray_SimpleNewFromData(1,
                                                                 &dim,
                                                                 PDM_G_NUM_NPY_INT,
                                                        <void *> entity_ln_to_gn)
            PyArray_ENABLEFLAGS(np_entity_ln_to_gn, NPY.NPY_OWNDATA);

        return {'np_entity_ln_to_gn'     : np_entity_ln_to_gn}


    # ------------------------------------------------------------------
    def multipart_part_color_get(self, int ipart, int zone_gid, PDM_mesh_entities_t entity_type):
        """
           Get partition ghost information
        """
        # ************************************************************************
        # > Declaration
        cdef int  *entity_color
        # ************************************************************************

        # -> Call PPART to get info
        n_entity1 = PDM_multipart_partition_color_get(self._mtp,
                                                      zone_gid,
                                                      ipart,
                                                      entity_type,
                                                      &entity_color,
                                                      PDM_OWNERSHIP_USER)
        cdef NPY.npy_intp dim
        if (entity_color == NULL) :
            np_entity_color = None
        else :
            dim = <NPY.npy_intp> n_entity1
            np_entity_color   = NPY.PyArray_SimpleNewFromData(1,
                                                                 &dim,
                                                                 NPY.NPY_INT32,
                                                        <void *> entity_color)
            PyArray_ENABLEFLAGS(np_entity_color, NPY.NPY_OWNDATA);

        return {'np_entity_color'     : np_entity_color}

    # ------------------------------------------------------------------
    def multipart_time_get(self, int zone_gid):
        """
        Get times
        """
        # ************************************************************************
        # > Declaration
        cdef double *elapsed
        cdef double *cpu
        cdef double *cpu_user
        cdef double *cpu_sys
        # ************************************************************************

        PDM_multipart_time_get(self._mtp, zone_gid, &elapsed, &cpu, &cpu_user, &cpu_sys)

        d_elapsed = {'total'              : elapsed[0],
                     'building graph'     : elapsed[1],
                     'splitting graph'    : elapsed[2],
                     'building partitions': elapsed[3]}

        d_cpu     = {'total'              : cpu[0],
                     'building graph'     : cpu[1],
                     'splitting graph'    : cpu[2],
                     'building partitions': cpu[3]}

        d_cpu_user = {'total'              : cpu_user[0],
                      'building graph'     : cpu_user[1],
                      'splitting graph'    : cpu_user[2],
                      'building partitions': cpu_user[3]}

        d_cpu_sys = {'total'              : cpu_sys[0],
                     'building graph'     : cpu_sys[1],
                     'splitting graph'    : cpu_sys[2],
                     'building partitions': cpu_sys[3]}

        return {'elapsed' : d_elapsed, 'cpu' : d_cpu, 'cpu_user' : d_cpu_user,  'cpu_sys' : d_cpu_sys}


    # ------------------------------------------------------------------
    # def multipart_stat_get(self, int zone_gid):
    #     """
    #     Get statistics
    #     """
    #     # ************************************************************************
    #     # > Declaration
    #     cdef int      cells_average,
    #     cdef int      cells_median,
    #     cdef double   cells_std_deviation,
    #     cdef int      cells_min,
    #     cdef int      cells_max,
    #     cdef int      bound_part_faces_average,
    #     cdef int      bound_part_faces_median,
    #     cdef double   bound_part_faces_std_deviation,
    #     cdef int      bound_part_faces_min,
    #     cdef int      bound_part_faces_max,
    #     cdef int      bound_part_faces_sum
    #     # ************************************************************************

    #     PDM_multipart_stat_get(self._mtp,
    #                            zone_gid,
    #                            &cells_average,
    #                            &cells_median,
    #                            &cells_std_deviation,
    #                            &cells_min,
    #                            &cells_max,
    #                            &bound_part_faces_average,
    #                            &bound_part_faces_median,
    #                            &bound_part_faces_std_deviation,
    #                            &bound_part_faces_min,
    #                            &bound_part_faces_max,
    #                            &bound_part_faces_sum)

    #     return {'cells_average'                  : cells_average,
    #             'cells_median'                   : cells_median,
    #             'cells_std_deviation'            : cells_std_deviation,
    #             'cells_min'                      : cells_min,
    #             'cells_max'                      : cells_max,
    #             'bound_part_faces_average'       : bound_part_faces_average,
    #             'bound_part_faces_median'        : bound_part_faces_median,
    #             'bound_part_faces_std_deviation' : bound_part_faces_std_deviation,
    #             'bound_part_faces_min'           : bound_part_faces_min,
    #             'bound_part_faces_max'           : bound_part_faces_max,
    #             'bound_part_faces_sum'           : bound_part_faces_sum}
