cdef extern from "pdm_para_graph_dual.h":
    ctypedef enum PDM_split_dual_t:
        PDM_SPLIT_DUAL_WITH_PARMETIS = 1
        PDM_SPLIT_DUAL_WITH_PTSCOTCH = 2
cdef extern from "pdm_multipart.h":
    ctypedef enum PDM_part_size_t:
        PDM_PART_SIZE_HOMONEGEOUS   = 1
        PDM_PART_SIZE_HETEROGENEOUS = 2

    # -> PPART bases functions
    # ------------------------------------------------------------------
    # MPI_Comm      comm,
    int PDM_multipart_create(int              n_zone,
                             int*             n_part,
                             PDM_bool_t       merge_blocks,
                             PDM_split_dual_t split_method,
                             PDM_part_size_t  part_size_method,
                             double*          part_fraction,
                             PDM_MPI_Comm     comm)

    # ------------------------------------------------------------------
    void PDM_multipart_register_block(int        mpart_id,
                                      int        zoneGId,
                                      int        block_data_id)

    # ------------------------------------------------------------------
    void PDM_multipart_register_joins(int        mpart_id,
                                      int        n_total_joins,
                                      int*       matching_join_array)

    # ------------------------------------------------------------------
    void PDM_multipart_set_reordering_options(const int   mpart_id,
                                              const int   i_zone,
                                              const char *renum_cell_method,
                                              const int  *renum_cell_properties,
                                              const char *renum_face_method)

    # ------------------------------------------------------------------
    void PDM_multipart_run_ppart(int id);

    # ------------------------------------------------------------------
    void PDM_multipart_part_dim_get(int  mpartId,
                                    int  zoneGId,
                                    int  ipart,
                                    int *nCell,
                                    int *nFace,
                                    int *nFacePartBound,
                                    int *nVtx,
                                    int *nProc,
                                    int *nTPart,
                                    int *sCellFace,
                                    int *sFaceVtx,
                                    int *sFaceBound,
                                    int *nFaceBound,
                                    int *sFaceJoin,
                                    int *nFaceJoin);

    # ------------------------------------------------------------------
    void PDM_multipart_part_val_get(int            mpartId,
                                    int            zoneGId,
                                    int            ipart,
                                    int          **cellTag,
                                    int          **cellFaceIdx,
                                    int          **cellFace,
                                    PDM_g_num_t  **cellLNToGN,
                                    int          **faceTag,
                                    int          **faceCell,
                                    int          **faceVtxIdx,
                                    int          **faceVtx,
                                    PDM_g_num_t  **faceLNToGN,
                                    int          **facePartBoundProcIdx,
                                    int          **facePartBoundPartIdx,
                                    int          **facePartBound,
                                    int          **vtxTag,
                                    double       **vtx,
                                    PDM_g_num_t  **vtxLNToGN,
                                    int          **faceBoundIdx,
                                    int          **faceBound,
                                    PDM_g_num_t  **faceBoundLNToGN,
                                    int          **faceJoinIdx,
                                    int          **faceJoin,
                                    PDM_g_num_t  **faceJoinLNToGN)
    # ------------------------------------------------------------------
    void PDM_multipart_part_color_get(int            mpartId,
                                      int            zoneGId,
                                      int            ipart,
                                      int          **cellColor,
                                      int          **faceColor,
                                      int          **threadColor,
                                      int          **hyperPlaneColor)

    # ------------------------------------------------------------------
    void PDM_multipart_time_get(int       mpartId,
                                int       zoneGId,
                                double  **elapsed,
                                double  **cpu,
                                double  **cpu_user,
                                double  **cpu_sys);
    # ------------------------------------------------------------------
    void PDM_multipart_free(int id);

# ------------------------------------------------------------------
cdef class MultiPart:
    """
       MultiPpart
    """
    # > For Ppart
    cdef int _mpart_id
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
        self._mpart_id = PDM_multipart_create(n_zone,
                                              <int*> n_part.data,
                                              <PDM_bool_t> merge_blocks,
                                              split_method,
                                              part_size_method,
                                              part_fraction_data,
                                              PDM_MPI_mpi_2_pdm_mpi_comm (<void *> &c_comm))
        print("MultiPart::end ")

    # ------------------------------------------------------------------
    # def multi_part_partial_free(self):
    #     print 'multi_part_partial_free MultiPart'
    #     PDM_multipart_partial_free(self._mpart_id)

    # ------------------------------------------------------------------
    def __dealloc__(self):
        print('__dealloc__ MultiPart ')
        PDM_multipart_free(self._mpart_id)

    # ------------------------------------------------------------------
    def multipart_register_block(self, int zoneGId,
                                       int block_data_id):
        """
        """
        PDM_multipart_register_block(self._mpart_id, zoneGId, block_data_id)

    # ------------------------------------------------------------------
    def multipart_register_joins(self, int n_total_joins,
                                 NPY.ndarray[NPY.int32_t, mode='c', ndim=1] matching_join):
        """
        """
        PDM_multipart_register_joins(       self._mpart_id,
                                            n_total_joins,
                                     <int*> matching_join.data)
    # ------------------------------------------------------------------
    def multipart_set_reordering(self, int i_zone,
                                 char *renum_cell_method, char *renum_face_method,
                                 NPY.ndarray[NPY.int32_t, mode='c', ndim=1] renum_properties_cell):
      """
      """
      cdef int *renum_properties_cell_data
      if (renum_properties_cell is None):
        renum_properties_cell_data = NULL
      else:
        renum_properties_cell_data = <int *> renum_properties_cell.data
      PDM_multipart_set_reordering_options(self._mpart_id,
                                           i_zone,
                                           renum_cell_method,
                                           renum_properties_cell_data,
                                           renum_face_method)
    # ------------------------------------------------------------------
    def multipart_run_ppart(self):
        """
        """
        PDM_multipart_run_ppart(self._mpart_id)

    # ------------------------------------------------------------------
    def multipart_dim_get(self, int ipart, int zoneGId):
        """
           Get partition dimensions
        """
        # ************************************************************************
        # > Declaration
        cdef int nCell
        cdef int nProc
        cdef int nTPart
        cdef int nFace
        cdef int nFacePartBound
        cdef int nVertex
        cdef int sCellFace
        cdef int sFaceVertex
        cdef int sFaceBound
        cdef int nFaceBound
        cdef int sFaceJoin
        cdef int nFaceJoin
        # ************************************************************************

        PDM_multipart_part_dim_get(self._mpart_id,
                                   zoneGId,
                                   ipart,
                                   &nCell,
                                   &nFace,
                                   &nFacePartBound,
                                   &nVertex,
                                   &nProc,
                                   &nTPart,
                                   &sCellFace,
                                   &sFaceVertex,
                                   &sFaceBound,
                                   &nFaceBound,
                                   &sFaceJoin,
                                   &nFaceJoin)

        return {'nCell'          :nCell,
                'ipart'          :ipart,
                'nFace'          :nFace,
                'nTPart'         :nTPart,
                'nProc'          :nProc,
                'nFacePartBound' :nFacePartBound,
                'nVertex'        :nVertex,
                'sCellFace'      :sCellFace,
                'sFaceVertex'    :sFaceVertex,
                'sFaceBound'     :sFaceBound,
                'nFaceBound'     :nFaceBound,
                'sFaceJoin'      :sFaceJoin,
                'nFaceJoin'      :nFaceJoin}

    # ------------------------------------------------------------------
    def multipart_val_get(self, int ipart, int zoneGId):
        """
           Get partition dimensions
        """
        # ************************************************************************
        # > Declaration
        cdef int          *cellTag,
        cdef int          *cellFaceIdx,
        cdef int          *cellFace,
        cdef PDM_g_num_t  *cellLNToGN,
        cdef int          *faceTag,
        cdef int          *faceCell,
        cdef int          *faceVertexIdx,
        cdef int          *faceVertex,
        cdef PDM_g_num_t  *faceLNToGN,
        cdef int          *facePartBound,
        cdef int          *facePartBoundProcIdx,
        cdef int          *facePartBoundPartIdx,
        cdef int          *vertexTag,
        cdef double       *vertex,
        cdef PDM_g_num_t  *vertexLNToGN,
        cdef int          *faceBoundIdx,
        cdef int          *faceBound,
        cdef PDM_g_num_t  *faceBoundLNToGN,
        cdef int          *faceJoinIdx,
        cdef int          *faceJoin,
        cdef PDM_g_num_t  *faceJoinLNToGN
        # ************************************************************************

        # dims = self.part_dim_get(self._mpart_id, ipart)
        dims = self.multipart_dim_get(ipart, zoneGId)

        # -> Call PPART to get info
        PDM_multipart_part_val_get(self._mpart_id,
                                   zoneGId,
                                   ipart,
                                   &cellTag,
                                   &cellFaceIdx,
                                   &cellFace,
                                   &cellLNToGN,
                                   &faceTag,
                                   &faceCell,
                                   &faceVertexIdx,
                                   &faceVertex,
                                   &faceLNToGN,
                                   &facePartBoundProcIdx,
                                   &facePartBoundPartIdx,
                                   &facePartBound,
                                   &vertexTag,
                                   &vertex,
                                   &vertexLNToGN,
                                   &faceBoundIdx,
                                   &faceBound,
                                   &faceBoundLNToGN,
                                   &faceJoinIdx,
                                   &faceJoin,
                                   &faceJoinLNToGN)
        # -> Begin
        cdef NPY.npy_intp dim

        # \param [out]  cellTag            Cell tag (size = nCell)
        if (cellTag == NULL) :
            npCellTag = None
        else :
            dim = <NPY.npy_intp> dims['nCell']
            npCellTag = NPY.PyArray_SimpleNewFromData(1,
                                                     &dim,
                                                     NPY.NPY_INT32,
                                                     <void *> cellTag)

        # \param [out]  cellFaceIdx        Cell to face connectivity index (size = nCell + 1)
        if (cellFaceIdx == NULL) :
            npCellFaceIdx = None
        else :
            dim = <NPY.npy_intp> (dims['nCell'] + 1)
            npCellFaceIdx = NPY.PyArray_SimpleNewFromData(1,
                                                         &dim,
                                                         NPY.NPY_INT32,
                                                         <void *> cellFaceIdx)

        # \param [out]  cellFace           Cell to face connectivity (size = cellFaceIdx[nCell] = lCellFace)
        if (cellFace == NULL) :
            npCellFace = None
        else :
            dim = <NPY.npy_intp> dims['sCellFace']
            npCellFace = NPY.PyArray_SimpleNewFromData(1,
                                                      &dim,
                                                      NPY.NPY_INT32,
                                                      <void *> cellFace)

        # \param [out]  cellLNToGN         Cell local numbering to global numbering (size = nCell)
        # dim = <NPY.npy_intp> dims['nCell']
        if (cellLNToGN == NULL) :
            npCellLNToGN = None
        else :
            dim = <NPY.npy_intp> dims['nCell']
            npCellLNToGN = NPY.PyArray_SimpleNewFromData(1,
                                                        &dim,
                                                        PDM_G_NUM_NPY_INT,
                                                        <void *> cellLNToGN)

        # \param [out]  faceTag            Face tag (size = nFace)
        if (faceTag == NULL) :
            npFaceTag = None
        else :
            dim = <NPY.npy_intp> dims['nFace']
            npFaceTag = NPY.PyArray_SimpleNewFromData(1,
                                                     &dim,
                                                     NPY.NPY_INT32,
                                                     <void *> faceTag)

        # \param [out]  faceCell           Face to cell connectivity  (size = 2 * nFace)
        if (faceCell == NULL) :
            npFaceCell = None
        else :
            dim = <NPY.npy_intp> (2 * dims['nFace'])
            npFaceCell = NPY.PyArray_SimpleNewFromData(1,
                                                      &dim,
                                                      NPY.NPY_INT32,
                                                      <void *> faceCell)

        # \param [out]  faceVtxIdx         Face to Vertex connectivity index (size = nFace + 1)
        if (faceVertexIdx == NULL) :
            npFaceVertexIdx = None
        else :
            dim = <NPY.npy_intp> (dims['nFace'] + 1)
            npFaceVertexIdx = NPY.PyArray_SimpleNewFromData(1,
                                                           &dim,
                                                           NPY.NPY_INT32,
                                                           <void *> faceVertexIdx)

        # \param [out]  faceVtx            Face to Vertex connectivity (size = faceVtxIdx[nFace])
        cdef NPY.ndarray[NPY.int32_t, ndim=1] npFaceVertex
        if (faceVertex == NULL) :
            npFaceVertex = None
        else :
            dim = <NPY.npy_intp> dims['sFaceVertex']
            npFaceVertex  = NPY.PyArray_SimpleNewFromData(1,
                                                         &dim,
                                                         NPY.NPY_INT32,
                                                         <void *> faceVertex)
            # PyArray_ENABLEFLAGS(npFaceVertex, NPY.NPY_OWNDATA)
            # print '*'*1000
            # print 'Take ownership'
            # print '*'*1000

        # \param [out]  faceLNToGN         Face local numbering to global numbering (size = nFace)
        if (faceLNToGN == NULL) :
            npFaceLNToGN = None
        else :
            dim = <NPY.npy_intp> dims['nFace']
            npFaceLNToGN   = NPY.PyArray_SimpleNewFromData(1,
                                                          &dim,
                                                          PDM_G_NUM_NPY_INT,
                                                          <void *> faceLNToGN)

        # \param [out]  facePartBound      Partitioning boundary faces
        if (facePartBound == NULL) :
            npFacePartBound = None
        else :
            dim = <NPY.npy_intp> (4 * dims['nFacePartBound'])
            npFacePartBound   = NPY.PyArray_SimpleNewFromData(1,
                                                             &dim,
                                                             NPY.NPY_INT32,
                                                             <void *> facePartBound)

        # \param [out]  facePartBoundProcIdx  Partitioning boundary faces block distribution from processus (size = nProc + 1)
        if (facePartBoundProcIdx == NULL) :
            npfacePartBoundProcIdx = None
        else :
            dim = <NPY.npy_intp> ( dims['nProc'] + 1)
            npfacePartBoundProcIdx   = NPY.PyArray_SimpleNewFromData(1,
                                                             &dim,
                                                             NPY.NPY_INT32,
                                                             <void *> facePartBoundProcIdx)

        # \param [out]  facePartBoundPartIdx  Partitioning boundary faces block distribution from partition (size = nTPart + 1)
        if (facePartBoundPartIdx == NULL) :
            npfacePartBoundPartIdx = None
        else :
            dim = <NPY.npy_intp> ( dims['nTPart'] + 1)
            npfacePartBoundPartIdx   = NPY.PyArray_SimpleNewFromData(1,
                                                             &dim,
                                                             NPY.NPY_INT32,
                                                             <void *> facePartBoundPartIdx)

        # \param [out]  vtxTag             Vertex tag (size = nVtx)
        if (vertexTag == NULL) :
            npVertexTag = None
        else :
            dim = <NPY.npy_intp> dims['nVertex']
            npVertexTag   = NPY.PyArray_SimpleNewFromData(1,
                                                         &dim,
                                                         NPY.NPY_INT32,
                                                         <void *> vertexTag)

        # \param [out]  vtx                Vertex coordinates (size = 3 * nVtx)
        if (vertex == NULL) :
            npVertex = None
        else :
            dim = <NPY.npy_intp> (3 * dims['nVertex'])
            npVertex  = NPY.PyArray_SimpleNewFromData(1,
                                                     &dim,
                                                     NPY.NPY_DOUBLE,
                                                     <void *> vertex)

        # \param [out]  vtxLNToGN          Vertex local numbering to global numbering (size = nVtx)
        if (vertexLNToGN == NULL) :
            npVertexLNToGN = None
        else :
            dim = <NPY.npy_intp> dims['nVertex']
            npVertexLNToGN  = NPY.PyArray_SimpleNewFromData(1,
                                                           &dim,
                                                           PDM_G_NUM_NPY_INT,
                                                           <void *> vertexLNToGN)

        # \param [out]  faceBoundIdx       face group index (size = nFaceBound + 1)
        if (faceBoundIdx == NULL) :
            npFaceBoundIdx = None
        else :
            dim = <NPY.npy_intp> dims['nFaceBound'] + 1
            npFaceBoundIdx  = NPY.PyArray_SimpleNewFromData(1,
                                                           &dim,
                                                           NPY.NPY_INT32,
                                                           <void *> faceBoundIdx)

        # \param [out]  faceBound          faces for each group (size = faceBoundIdx[nFaceBound] = lFaceBound)
        if (faceBound == NULL) :
            npFaceBound = None
        else :
            dim = <NPY.npy_intp> dims['sFaceBound']
            npFaceBound = NPY.PyArray_SimpleNewFromData(1,
                                                       &dim,
                                                       NPY.NPY_INT32,
                                                       <void *> faceBound)

        # \param [out]  faceBoundLNToGN    faces global numbering for each group (size = faceBoundIdx[nFaceBound] = lFaceBound)
        if (faceBoundLNToGN == NULL) :
            npFaceBoundLNToGN = None
        else :
            dim = <NPY.npy_intp> dims['sFaceBound']
            npFaceBoundLNToGN = NPY.PyArray_SimpleNewFromData(1,
                                                             &dim,
                                                             PDM_G_NUM_NPY_INT,
                                                             <void *> faceBoundLNToGN)

        # \param [out]  faceJoinIdx       face group index (size = nFaceJoin + 1)
        if (faceJoinIdx == NULL) :
            npFaceJoinIdx = None
        else :
            dim = <NPY.npy_intp> dims['nFaceJoin'] + 1
            npFaceJoinIdx  = NPY.PyArray_SimpleNewFromData(1,
                                                          &dim,
                                                          NPY.NPY_INT32,
                                                          <void *> faceJoinIdx)

        # \param [out]  faceJoin          faces for each group (size = faceJoinIdx[nFaceJoin] = lFaceJoin)
        if (faceJoin == NULL) :
            npFaceJoin = None
        else :
            dim = <NPY.npy_intp> (4 * dims['sFaceJoin'])
            npFaceJoin = NPY.PyArray_SimpleNewFromData(1,
                                                      &dim,
                                                      NPY.NPY_INT32,
                                                      <void *> faceJoin)

        # \param [out]  faceJoinLNToGN    faces global numbering for each group (size = faceJoinIdx[nFaceJoin] = lFaceJoin)
        if (faceJoinLNToGN == NULL) :
            npFaceJoinLNToGN = None
        else :
            dim = <NPY.npy_intp> dims['sFaceJoin']
            npFaceJoinLNToGN = NPY.PyArray_SimpleNewFromData(1,
                                                             &dim,
                                                             PDM_G_NUM_NPY_INT,
                                                             <void *> faceJoinLNToGN)

        return {'npCellTag'                  : npCellTag,
                'npCellFaceIdx'              : npCellFaceIdx,
                'npCellFace'                 : npCellFace,
                'npCellLNToGN'               : npCellLNToGN,
                'npFaceTag'                  : npFaceTag,
                'npFaceCell'                 : npFaceCell,
                'npFaceVertexIdx'            : npFaceVertexIdx,
                'npFaceVertex'               : npFaceVertex,
                'npFaceLNToGN'               : npFaceLNToGN,
                'npfacePartBoundProcIdx'     : npfacePartBoundProcIdx,
                'npfacePartBoundPartIdx'     : npfacePartBoundPartIdx,
                'npFacePartBound'            : npFacePartBound,
                'npVertexTag'                : npVertexTag,
                'npVertex'                   : npVertex,
                'npVertexLNToGN'             : npVertexLNToGN,
                'npFaceBoundIdx'             : npFaceBoundIdx,
                'npFaceBound'                : npFaceBound,
                'npFaceBoundLNToGN'          : npFaceBoundLNToGN,
                'npFaceJoinIdx'              : npFaceJoinIdx,
                'npFaceJoin'                 : npFaceJoin,
                'npFaceJoinLNToGN'           : npFaceJoinLNToGN}

    # ------------------------------------------------------------------
    def multipart_color_get(self, int ipart, int zoneGId):
        """
           Get partition dimensions
        """
        # ************************************************************************
        # > Declaration
        cdef int          *cellColor,
        cdef int          *faceColor
        cdef int          *threadColor
        cdef int          *hyperPlaneColor
        # ************************************************************************

        # dims = self.part_dim_get(self._mpart_id, ipart)
        dims = self.multipart_dim_get(ipart, zoneGId)

        # -> Call PPART to get info
        PDM_multipart_part_color_get(self._mpart_id,
                                     zoneGId,
                                     ipart,
                                     &cellColor,
                                     &faceColor,
                                     &threadColor,
                                     &hyperPlaneColor)
        # -> Begin
        cdef NPY.npy_intp dim

        # \param [out]  cellColor            Cell tag (size = nCell)
        if (cellColor == NULL):
            npCellColor = None
        else :
            dim = <NPY.npy_intp> dims['nCell']
            npCellColor = NPY.PyArray_SimpleNewFromData(1,
                                                        &dim,
                                                        NPY.NPY_INT32,
                                                        <void *> cellColor)
        # \param [out]  faceColor            Cell tag (size = nFace)
        if (faceColor == NULL):
            npFaceColor = None
        else :
            dim = <NPY.npy_intp> dims['nFace']
            npFaceColor = NPY.PyArray_SimpleNewFromData(1,
                                                        &dim,
                                                        NPY.NPY_INT32,
                                                        <void *> faceColor)

        # \param [out]  threadColor            Cell tag (size = nCell)
        if (threadColor == NULL):
            npThreadColor = None
        else :
            dim = <NPY.npy_intp> dims['nCell']
            npThreadColor = NPY.PyArray_SimpleNewFromData(1,
                                                          &dim,
                                                          NPY.NPY_INT32,
                                                          <void *> threadColor)

        # \param [out]  hyperPlaneColor            Cell tag (size = nCell)
        if (hyperPlaneColor == NULL):
            npHyperPlaneColor = None
        else :
            dim = <NPY.npy_intp> dims['nCell']
            npHyperPlaneColor = NPY.PyArray_SimpleNewFromData(1,
                                                              &dim,
                                                              NPY.NPY_INT32,
                                                              <void *> hyperPlaneColor)
        return {'npCellColor'       : npCellColor,
                'npFaceColor'       : npFaceColor,
                'npThreadColor'     : npThreadColor,
                'npHyperPlaneColor' : npHyperPlaneColor}

    # ------------------------------------------------------------------
    def multipart_time_get(self, int zoneGId):
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

        PDM_multipart_time_get(self._mpart_id, zoneGId, &elapsed, &cpu, &cpu_user, &cpu_sys)

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
    # def multipart_stat_get(self, int zoneGId):
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

    #     PDM_multipart_stat_get(self._mpart_id,
    #                            zoneGId,
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
