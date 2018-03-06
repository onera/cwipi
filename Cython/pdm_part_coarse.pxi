
# cdef extern from "pdm_part_priv.h":
#     ctypedef struct _part_t:
#        int           nVertex              
#        int           nCell             
#        int           nFace             
#        int           nFacePartBound    
#        int          *cellFaceIdx
#        PDM_g_num_t  *gCellFace        
#        int          *cellFace        
#        PDM_g_num_t  *cellLNToGN       
#        int          *cellTag         
#        int          *faceCell        
#        int          *faceVtxIdx      
#        PDM_g_num_t  *gFaceVtx         
#        int          *faceVtx         
#        PDM_g_num_t  *faceLNToGN       
#        int          *faceTag         
#        int          *facePartBoundProcIdx  
#        int          *facePartBoundPartIdx 
#        int          *facePartBound     
#        int          *faceGroupIdx      
#        int          *faceGroup         
#        PDM_g_num_t  *faceGroupLNToGN   
#        double       *vtx              
#        PDM_g_num_t  *vtxLNToGN        
#        int          *vtxTag           


# cdef extern from "PDM_part_coarse_mesh_priv.h":
#     ctypedef struct _coarse_part_t:
#        _part_t    *part
#        int        *coarseCellCellIdx 
#        int        *coarseCellCell
#        int        *coarseFaceGroupToFineFaceGroup
#        int        *coarseFaceToFineFace
#        int        *coarseVtxToFineVtx

#     ctypedef struct _coarse_mesh_t:
#        pass



cdef extern from "pdm_part_coarse_mesh.h":
    
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Wrapping of Ppart Structure 

    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Wrapping of function 
    void PDM_part_coarse_mesh_create(int           *cmId,
                                     PDM_MPI_Comm  comm,
                                     int           method,
                                     int           nPart, 
                                     int           nTPart,
                                     int           nFaceGroup, 
                                     int           have_cellTag,
                                     int           have_faceTag,
                                     int           have_vtxTag,
                                     int           have_cellWeight,
                                     int           have_faceWeight,
                                     int           have_faceGroup)

    void PDM_part_coarse_mesh_input(int           cmId,
                                    int           iPart,
                                    int           nCoarseCell,
                                    int           nCell,
                                    int           nFace,
                                    int           nVertex,
                                    int           nFaceGroup,
                                    int           nFacePartBound,
                                    int          *cellFaceIdx,
                                    int          *cellFace,
                                    int          *cellTag,
                                    int          *cellWeight,
                                    int          *faceWeight,
                                    PDM_g_num_t  *cellLNToGN,       
                                    int          *faceCell,
                                    int          *faceVtxIdx,
                                    int          *faceVtx,
                                    int          *faceTag,       
                                    PDM_g_num_t  *faceLNToGN,       
                                    double       *vtxCoord,
                                    int          *vtxTag,
                                    PDM_g_num_t  *vtxLNToGN,       
                                    int          *faceGroupIdx,
                                    int          *faceGroup,
                                    PDM_g_num_t  *faceGroupLNToGN,
                                    int          *facePartBoundProcIdx,       
                                    int          *facePartBoundPartIdx,
                                    int          *facePartBound)

    void PDM_part_coarse_mesh_compute(int cmId)

    void PDM_part_coarse_mesh_part_dim_get(int  cmId,
                                           int  iPart,
                                           int *nCell,
                                           int *nFace,
                                           int *nFacePartBound,
                                           int *nVertex,
                                           int *nProc,
                                           int *nTPart,
                                           int *nFaceGroup,
                                           int *sCellFace,
                                           int *sFaceVtx,
                                           int *sFaceGroup,
                                           int *sCoarseCellToFineCell)

    void PDM_part_coarse_mesh_part_get(int            cmId,
                                       int            iPart,       
                                       int          **cellFaceIdx,
                                       int          **cellFace,
                                       int          **cellTag,
                                       PDM_g_num_t  **cellLNToGN,
                                       int          **cellInitCellIdx,                  
                                       int          **cellInitCell,          
                                       int          **faceCell,
                                       int          **faceVtxIdx,
                                       int          **faceVtx,
                                       int          **faceTag,
                                       PDM_g_num_t  **faceLNToGN,  
                                       int          **faceGroupInitFaceGroup,
                                       int          **faceInitFace,          
                                       double       **vtxCoord,
                                       int          **vtxTag,
                                       PDM_g_num_t **vtxLNToGN,        
                                       int          **vtxInitVtx,          
                                       int          **faceGroupIdx,
                                       int          **faceGroup,
                                       PDM_g_num_t  **faceGroupLNToGN,
                                       int          **facePartBoundProcIdx,
                                       int          **facePartBoundPartIdx,
                                       int          **facePartBound)
    
    void PDM_part_coarse_mesh_part_set_anisotropic_info(int   cmId,
                                                        int   iPart,       
                                                        int  *agglomerationLines,
                                                        int  *agglomerationLinesIdx,
                                                        int   agglomerationLinesIdx_size,
                                                        int  *isOnFineBnd)
   
    void PDM_part_coarse_mesh_part_get_anisotropic_info(int   cmId,
                                                        int   iPart,       
                                                        int **agglomerationLines,
                                                        int **agglomerationLinesIdx,
                                                        int   *agglomerationLinesIdx_size,
                                                        int **isOnFineBnd)
    
    void PDM_part_coarse_mesh_add_option_anisotropic(int  cmId, 
                                                     int *anisotropicOption)
    
    void PDM_part_coarse_mesh_free(int cmId)

    void PDM_part_coarse_mesh_time_get(int       cmId,
                                       double  **elapsed,
                                       double  **cpu,
                                       double  **cpu_user,
                                       double  **cpu_sys)

    void PDM_part_coarse_mesh_display(int    cmId)



# ------------------------------------------------------------------
cdef class CoarseMesh:
    """
    Define a coarse mesh partitioning with PDM Library
    """
    cdef int _cmId
    cdef int _nFaceGroup
    # ------------------------------------------------------------------
    def __init__(self,
                  MPI.Comm comm,
                  int      method, 
                  int      nPart, 
                  int      nTPart,
                  int      nFaceGroup,
                  int      have_cellTag    = 0, 
                  int      have_faceTag    = 0, 
                  int      have_vtxTag     = 0, 
                  int      have_cellWeight = 0, 
                  int      have_faceWeight = 0, 
                  int      have_faceGroup  = 0): 
      """
        Create a Coarse Ppart partition from a existing Ppart partition 
        with PDM Library ( Developed at ONERA by Eric Quemerais )

      """

      # ~> Communicator Mpi
      cdef MPI.MPI_Comm c_comm = comm.ob_mpi

      # ~> Set _nFaceGroup
      self._nFaceGroup =  nFaceGroup

      # > Create 
      PDM_part_coarse_mesh_create(&self._cmId, 
                                  PDM_MPI_mpi_2_pdm_mpi_comm (<void *> &c_comm), 
                                  method, 
                                  nPart, 
                                  nTPart, 
                                  self._nFaceGroup, 
                                  have_cellTag,    
                                  have_faceTag,   
                                  have_vtxTag,     
                                  have_cellWeight, 
                                  have_faceWeight, 
                                  have_faceGroup)
      
    # ------------------------------------------------------------------
    def add_option_anisotropic(self, 
                               NPY.ndarray[NPY.int32_t, mode='c', ndim=1] anisotropicOption,
                               ):
        """

        """
        # ************************************************************************
        # > Declaration
        cdef int             *anisotropicOption_data
        # ************************************************************************
        
        # assert(anisotropicOption.shape[0] == 12)
        
        anisotropicOption_data = <int *> anisotropicOption.data
        
        PDM_part_coarse_mesh_add_option_anisotropic(self._cmId, 
                                                    anisotropicOption_data)
        
    # ------------------------------------------------------------------
    def set_mesh_input(self, 
                       int      iPart, 
                       int      nCoarseCell,
                       int      nCell,
                       int      nFace,
                       int      nVertex,
                       int      nFaceGroup,
                       int      nFacePartBound,
                       NPY.ndarray[NPY.int32_t    , mode='c', ndim=1] CellFaceIdx not None,
                       NPY.ndarray[NPY.int32_t    , mode='c', ndim=1] CellFace not None,
                       NPY.ndarray[NPY.int32_t    , mode='c', ndim=1] CellTag,
                       NPY.ndarray[NPY.int32_t    , mode='c', ndim=1] CellWeight,
                       NPY.ndarray[NPY.int32_t    , mode='c', ndim=1] FaceWeight,
                       NPY.ndarray[npy_pdm_gnum_t , mode='c', ndim=1] CellLNToGN not None,       
                       NPY.ndarray[NPY.int32_t    , mode='c', ndim=1] FaceCell not None,
                       NPY.ndarray[NPY.int32_t    , mode='c', ndim=1] FaceVtxIdx not None,
                       NPY.ndarray[NPY.int32_t    , mode='c', ndim=1] FaceVtx not None,
                       NPY.ndarray[NPY.int32_t    , mode='c', ndim=1] FaceTag ,       
                       NPY.ndarray[npy_pdm_gnum_t , mode='c', ndim=1] FaceLNToGN not None,       
                       NPY.ndarray[NPY.double_t   , mode='c', ndim=1] VtxCoord not None,
                       NPY.ndarray[NPY.int32_t    , mode='c', ndim=1] VtxTag,
                       NPY.ndarray[npy_pdm_gnum_t , mode='c', ndim=1] VtxLNToGN not None,       
                       NPY.ndarray[NPY.int32_t    , mode='c', ndim=1] FaceGroupIdx,
                       NPY.ndarray[NPY.int32_t    , mode='c', ndim=1] FaceGroup,
                       NPY.ndarray[npy_pdm_gnum_t , mode='c', ndim=1] FaceGroupLNToGN,
                       NPY.ndarray[NPY.int32_t    , mode='c', ndim=1] FacePartBoundProcIdx,
                       NPY.ndarray[NPY.int32_t    , mode='c', ndim=1] FacePartBoundPartIdx,
                       NPY.ndarray[NPY.int32_t    , mode='c', ndim=1] FacePartBound
                       ):
        """

        """
        # ************************************************************************
        # > Declaration
        # > Cell entity 
        cdef int             *CellFaceIdx_data
        cdef int             *CellFace_data     
        cdef int             *CellTag_data
        cdef int             *CellWeight_data
        cdef int             *FaceWeight_data
        cdef PDM_g_num_t     *CellLNToGN_data
        # > Face entity
        cdef int             *FaceCell_data
        cdef int             *FaceTag_data
        cdef PDM_g_num_t     *FaceLNToGN_data
        cdef int             *FaceVtxIdx_data
        cdef int             *FaceVtx_data
        # > Vertices entity
        cdef double          *VtxCoord_data
        cdef PDM_g_num_t     *VtxLNToGN_data
        cdef int             *VtxTag_data
        # > Boundary face conditions
        cdef int             *FaceGroup_data
        cdef int             *FaceGroupIdx_data
        cdef PDM_g_num_t     *FaceGroupLNToGN_data
        cdef int             *FacePartBoundProcIdx_data
        cdef int             *FacePartBoundPartIdx_data
        cdef int             *FacePartBound_data
        # ************************************************************************

        # :::::::::::::::::::::::::::::::::::::::::::::::::::::
        CellFaceIdx_data = <int *>         CellFaceIdx.data
        CellFace_data    = <int *>         CellFace.data
        CellLNToGN_data  = <PDM_g_num_t *> CellLNToGN.data

        # \param [in]   CellTag       Cell tag (size : nCell) or NULL
        if (CellTag is None):
            CellTag_data = NULL
        else:
            CellTag_data = <int *> CellTag.data

        # \param [in]   CellWeight    Cell weight (size : nCell) or NULL
        if (CellWeight is None):
            CellWeight_data = NULL
        else:
            CellWeight_data = <int *> CellWeight.data

        # :::::::::::::::::::::::::::::::::::::::::::::::::::::
        FaceCell_data   = <int *>         FaceCell.data
        FaceLNToGN_data = <PDM_g_num_t *> FaceLNToGN.data
        FaceVtx_data    = <int *>         FaceVtx.data
        FaceVtxIdx_data = <int *>         FaceVtxIdx.data

        # \param [in]   FaceTag       Distributed face tag
        if (FaceTag is None):
          FaceTag_data = NULL
        else:
          FaceTag_data = <int *> FaceTag.data

        # \param [in]   faceWeight    Face weight (size : nFace) or NULL
        if (FaceWeight is None):
            FaceWeight_data = NULL
        else:
            FaceWeight_data = <int *> FaceWeight.data

        # :::::::::::::::::::::::::::::::::::::::::::::::::::::
        VtxCoord_data  = <double *>      VtxCoord.data
        VtxLNToGN_data = <PDM_g_num_t *> VtxLNToGN.data

        if (VtxTag is None):
          VtxTag_data = NULL
        else:
          VtxTag_data = <int *> VtxTag.data

        # :::::::::::::::::::::::::::::::::::::::::::::::::::::
        if (FaceGroup is None):
          FaceGroup_data = NULL
        else:
          FaceGroup_data            = <int *>         FaceGroup.data
          
        if (FaceGroupIdx is None):
          FaceGroupIdx_data = NULL
        else:
          FaceGroupIdx_data            = <int *>         FaceGroupIdx.data
          
        if (FaceGroupLNToGN is None):
          FaceGroupLNToGN_data = NULL
        else:
          FaceGroupLNToGN_data            = <PDM_g_num_t *> FaceGroupLNToGN.data
                 
        FacePartBoundProcIdx_data = <int *>         FacePartBoundProcIdx.data
        FacePartBoundPartIdx_data = <int *>         FacePartBoundPartIdx.data
        FacePartBound_data        = <int *>         FacePartBound.data

        # :::::::::::::::::::::::::::::::::::::::::::::::::::::
        # > Fill input mesh 
        PDM_part_coarse_mesh_input(self._cmId, 
                                   iPart, 
                                   nCoarseCell, 
                                   nCell, 
                                   nFace, 
                                   nVertex, 
                                   nFaceGroup, 
                                   nFacePartBound,
                                   CellFaceIdx_data,
                                   CellFace_data,
                                   CellTag_data,
                                   CellWeight_data,
                                   FaceWeight_data,
                                   CellLNToGN_data,       
                                   FaceCell_data,
                                   FaceVtxIdx_data,
                                   FaceVtx_data,
                                   FaceTag_data,       
                                   FaceLNToGN_data,       
                                   VtxCoord_data,
                                   VtxTag_data,
                                   VtxLNToGN_data,       
                                   FaceGroupIdx_data,
                                   FaceGroup_data,
                                   FaceGroupLNToGN_data,
                                   FacePartBoundProcIdx_data,       
                                   FacePartBoundPartIdx_data,
                                   FacePartBound_data)

    # ------------------------------------------------------------------
    def set_mesh_input_anisotropic(self, 
                                   int      iPart,
                                   NPY.ndarray[NPY.int32_t, mode='c', ndim=1] agglomerationLinesInit,
                                   NPY.ndarray[NPY.int32_t, mode='c', ndim=1] agglomerationLinesInitIdx,
                                   NPY.ndarray[NPY.int32_t, mode='c', ndim=1] isOnFineBndInit
                                   ):
        """

        """
        # ************************************************************************
        # > Declaration
        # > Cell entity 
        cdef int             *agglomerationLinesInit_data
        cdef int             *agglomerationLinesInitIdx_data
        cdef int             agglomerationLinesInitIdx_size 
        cdef int             *isOnFineBnd_data
        # ************************************************************************
        
        # \param [in]  
        if (agglomerationLinesInit is None):
            agglomerationLinesInit_data = NULL
        else:
            agglomerationLinesInit_data = <int *> agglomerationLinesInit.data
            
        # \param [in]  
        if (agglomerationLinesInitIdx is None):
            agglomerationLinesInitIdx_data  = NULL
            agglomerationLinesInitIdx_size  = 0
        else:
            agglomerationLinesInitIdx_data = <int *> agglomerationLinesInitIdx.data
            agglomerationLinesInitIdx_size  = agglomerationLinesInitIdx.shape[0]
            
        # \param [in]  
        if (isOnFineBndInit is None):
            isOnFineBndInit_data = NULL
        else:
            isOnFineBndInit_data = <int *> isOnFineBndInit.data
            
        print "\t\t\tset_mesh_input_anisotropic agglomerationLinesInitIdx_size=", agglomerationLinesInitIdx_size
        
        # :::::::::::::::::::::::::::::::::::::::::::::::::::::
        # > Fill input mesh for anisotropic
        PDM_part_coarse_mesh_part_set_anisotropic_info(self._cmId,
                                                       iPart,       
                                                       agglomerationLinesInit_data,
                                                       agglomerationLinesInitIdx_data,
                                                       agglomerationLinesInitIdx_size,
                                                       isOnFineBndInit_data)
        # :::::::::::::::::::::::::::::::::::::::::::::::::::::
        
    # ------------------------------------------------------------------
    def computeCoarseMesh(self):
        """
        Effective compute 
        """
        PDM_part_coarse_mesh_compute(self._cmId)

    # ------------------------------------------------------------------
    def __dealloc__(self):
        """
           Free memory in PTJ Lib and structure
        """
        PDM_part_coarse_mesh_free(self._cmId)

    # ------------------------------------------------------------------
    def part_coarse_dim_get(self, int ipart):
        """
           Get partition dimensions
        """
        # ************************************************************************
        # > Declaration
        cdef int nCell,
        cdef int nFace,
        cdef int nFacePartBound,
        cdef int nVertex,
        cdef int nProc,
        cdef int nTPart,
        cdef int nFaceGroup,
        cdef int sCellFace,
        cdef int sFaceVtx,
        cdef int sFaceGroup,
        cdef int sCoarseCellToFineCell
        # ************************************************************************

        PDM_part_coarse_mesh_part_dim_get(self._cmId,
                                          ipart,
                                          &nCell,
                                          &nFace,
                                          &nFacePartBound,
                                          &nVertex,
                                          &nProc,
                                          &nTPart,
                                          &nFaceGroup,
                                          &sCellFace,
                                          &sFaceVtx,
                                          &sFaceGroup,
                                          &sCoarseCellToFineCell)

        return {'nCell'                : nCell,
                'ipart'                : ipart,
                'nCell'                : nCell,
                'nFace'                : nFace,
                'nFacePartBound'       : nFacePartBound,
                'nVertex'              : nVertex,
                'nProc'                : nProc,
                'nTPart'               : nTPart,
                'nFaceGroup'           : nFaceGroup,
                'sCellFace'            : sCellFace,
                'sFaceVtx'             : sFaceVtx,
                'sFaceGroup'           : sFaceGroup,
                'sCoarseCellToFineCel' : sCoarseCellToFineCell
                }

    # ------------------------------------------------------------------
    def part_coarse_val_get(self, int ipart):
        """
           Get partition dimensions
        """
        # ************************************************************************
        # > Declaration
        # > Cell entity 
        cdef int          *cellTag
        cdef int          *cellFaceIdx
        cdef int          *cellFace
        cdef PDM_g_num_t  *cellLNToGN
        cdef int          *cellInitCellIdx                  
        cdef int          *cellInitCell  
        # > Face entity 
        cdef int          *faceTag
        cdef int          *faceCell
        cdef int          *faceVtxIdx
        cdef int          *faceVtx
        cdef PDM_g_num_t  *faceLNToGN
        cdef int          *faceGroupInitFaceGroup
        cdef int          *faceInitFace
        cdef int          *faceGroupIdx
        cdef int          *faceGroup
        cdef PDM_g_num_t  *faceGroupLNToGN
        cdef int          *facePartBound
        cdef int          *facePartBoundProcIdx
        cdef int          *facePartBoundPartIdx
        # > Vertices entity       
        cdef double       *vtxCoord
        cdef int          *vtxTag
        cdef PDM_g_num_t  *vtxLNToGN    
        cdef int          *vtxInitVtx
        # > For numpy capsule
        cdef NPY.npy_intp dim
        # ************************************************************************
        # > Get dim 
        dims = self.part_coarse_dim_get(ipart)

        # > Get array
        PDM_part_coarse_mesh_part_get(self._cmId,
                                      ipart,       
                                      &cellFaceIdx,
                                      &cellFace,
                                      &cellTag,
                                      &cellLNToGN,
                                      &cellInitCellIdx,                  
                                      &cellInitCell,          
                                      &faceCell,
                                      &faceVtxIdx,
                                      &faceVtx,
                                      &faceTag,
                                      &faceLNToGN,  
                                      &faceGroupInitFaceGroup,
                                      &faceInitFace,          
                                      &vtxCoord,
                                      &vtxTag,
                                      &vtxLNToGN,        
                                      &vtxInitVtx,          
                                      &faceGroupIdx,
                                      &faceGroup,
                                      &faceGroupLNToGN,
                                      &facePartBoundProcIdx,
                                      &facePartBoundPartIdx,
                                      &facePartBound)

        # :::::::::::::::::::::::::::::::::::::::::::::::::::::
        # > Translate to numpy capsule (Tout est pas fait encore )
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
        if (cellLNToGN == NULL) :
          npCellLNToGN = None
        else :
          dim = <NPY.npy_intp> dims['nCell']
          npCellLNToGN = NPY.PyArray_SimpleNewFromData(1,
                                                       &dim,
                                                       PDM_G_NUM_NPY_INT,
                                                       <void *> cellLNToGN)

        # \param [out]  cellInitCellIdx        Cell to face connectivity index (size = nCell + 1)
        if (cellInitCellIdx == NULL) :
          npCellInitCellIdx = None
        else :
          dim = <NPY.npy_intp> (dims['nCell'] + 1)
          npCellInitCellIdx = NPY.PyArray_SimpleNewFromData(1,
                                                            &dim,
                                                            NPY.NPY_INT32,
                                                            <void *> cellInitCellIdx)    

        # \param [out]  cellFaceIdx        Cell to face connectivity index (size = nCell + 1)
        if (cellInitCell == NULL) :
          npCellInitCell = None
        else :
          dim = <NPY.npy_intp> dims['sCoarseCellToFineCel'] 
          npCellInitCell = NPY.PyArray_SimpleNewFromData(1,
                                                         &dim,
                                                         NPY.NPY_INT32,
                                                         <void *> cellInitCell)                                                                    
        # :::::::::::::::::::::::::::::::::::::::::::::::::::::
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
        if (faceVtxIdx == NULL) :
          npFaceVertexIdx = None
        else :
          dim = <NPY.npy_intp> (dims['nFace'] + 1)
          npFaceVertexIdx = NPY.PyArray_SimpleNewFromData(1,
                                                         &dim,
                                                         NPY.NPY_INT32,
                                                         <void *> faceVtxIdx)

        # \param [out]  faceVtx            Face to Vertex connectivity (size = faceVtxIdx[nFace])
        if (faceVtx == NULL) :
          npFaceVertex = None
        else :
          dim = <NPY.npy_intp> dims['sFaceVtx']
          npFaceVertex  = NPY.PyArray_SimpleNewFromData(1,
                                                       &dim,
                                                       NPY.NPY_INT32,
                                                       <void *> faceVtx)

        # \param [out]  faceLNToGN         Face local numbering to global numbering (size = nFace)
        if (faceLNToGN == NULL) :
          npFaceLNToGN = None
        else :
          dim = <NPY.npy_intp> dims['nFace']
          npFaceLNToGN   = NPY.PyArray_SimpleNewFromData(1,
                                                         &dim,
                                                         PDM_G_NUM_NPY_INT,
                                                         <void *> faceLNToGN)

        if (faceInitFace == NULL) :
          npFaceInitFace = None
        else:
          # dim = <NPY.npy_intp> dims['nCoarseFace']
          dim = <NPY.npy_intp> dims['nFace']
          npFaceInitFace   = NPY.PyArray_SimpleNewFromData(1,
                                                           &dim,
                                                           NPY.NPY_INT32,
                                                           <void *> faceInitFace)
        if (faceGroupInitFaceGroup == NULL) :
          npfaceGroupInitFaceGroup = None
        else:
          # dim = <NPY.npy_intp> dims['nCoarseFace']
          dim = <NPY.npy_intp> dims['sFaceGroup']
          npfaceGroupInitFaceGroup   = NPY.PyArray_SimpleNewFromData(1,
                                                                     &dim,
                                                                     NPY.NPY_INT32,
                                                                     <void *> faceGroupInitFaceGroup)

        # :::::::::::::::::::::::::::::::::::::::::::::::::::::
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

        # \param [out]  faceGroupIdx       face group index (size = nFaceGroup + 1)
        if (faceGroupIdx == NULL) :
          npFaceGroupIdx = None
        else :
          dim = <NPY.npy_intp> (self._nFaceGroup + 1)
          npFaceGroupIdx  = NPY.PyArray_SimpleNewFromData(1,
                                                         &dim,
                                                         NPY.NPY_INT32,
                                                         <void *> faceGroupIdx)

        # \param [out]  faceGroup          faces for each group (size = faceGroupIdx[nFaceGroup] = lFaceGroup)
        if (faceGroup == NULL) :
          npFaceGroup = None
        else :
          dim = <NPY.npy_intp> dims['sFaceGroup']
          npFaceGroup = NPY.PyArray_SimpleNewFromData(1,
                                                     &dim,
                                                     NPY.NPY_INT32,
                                                     <void *> faceGroup)

        # \param [out]  faceGroupLNToGN    faces global numbering for each group (size = faceGroupIdx[nFaceGroup] = lFaceGroup)
        if (faceGroupLNToGN == NULL) :
          npFaceGroupLNToGN = None
        else :
          dim = <NPY.npy_intp> dims['sFaceGroup']
          npFaceGroupLNToGN = NPY.PyArray_SimpleNewFromData(1,
                                                            &dim,
                                                            PDM_G_NUM_NPY_INT,
                                                            <void *> faceGroupLNToGN)
        # :::::::::::::::::::::::::::::::::::::::::::::::::::::
        # \param [out]  vtxTag             Vertex tag (size = nVertex)
        if (vtxTag == NULL) :
          npVertexTag = None
        else :
          dim = <NPY.npy_intp> dims['nVertex']
          npVertexTag   = NPY.PyArray_SimpleNewFromData(1,
                                                       &dim,
                                                       NPY.NPY_INT32,
                                                       <void *> vtxTag)
        # \param [out]  vtx                Vertex coordinates (size = 3 * nVertex)
        if (vtxCoord == NULL) :
          npVertex = None
        else :
          dim = <NPY.npy_intp> (3 * dims['nVertex'])
          npVertex  = NPY.PyArray_SimpleNewFromData(1,
                                                   &dim,
                                                   NPY.NPY_DOUBLE,
                                                   <void *> vtxCoord)

        # \param [out]  vtxLNToGN          Vertex local numbering to global numbering (size = nVertex)
        if (vtxLNToGN == NULL) :
          npVertexLNToGN = None
        else :
          dim = <NPY.npy_intp> dims['nVertex']
          npVertexLNToGN  = NPY.PyArray_SimpleNewFromData(1,
                                                          &dim,
                                                          PDM_G_NUM_NPY_INT,
                                                          <void *> vtxLNToGN)
          
        # \param [out]  vtxLNToGN          Vertex local numbering to global numbering (size = nVertex)
        if (vtxInitVtx == NULL) :
          npVtxInitVtx = None
        else :
          dim = <NPY.npy_intp> dims['nVertex']
          npVtxInitVtx  = NPY.PyArray_SimpleNewFromData(1,
                                                        &dim,
                                                        NPY.NPY_INT32,
                                                        <void *> vtxInitVtx)
          
        # :::::::::::::::::::::::::::::::::::::::::::::::::::::
        return {'npCellTag'                  : npCellTag,
                'npCellFaceIdx'              : npCellFaceIdx,
                'npCellFace'                 : npCellFace,
                'npCellLNToGN'               : npCellLNToGN,
                'npCellInitCellIdx'          : npCellInitCellIdx,
                'npCellInitCell'             : npCellInitCell,
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
                'npVtxInitVtx'               : npVtxInitVtx,
                'npFaceGroupIdx'             : npFaceGroupIdx,
                'npFaceGroup'                : npFaceGroup,
                'npfaceGroupInitFaceGroup'   : npfaceGroupInitFaceGroup,
                'npFaceGroupLNToGN'          : npFaceGroupLNToGN, 
                'npFaceInitFace'             : npFaceInitFace}
        # :::::::::::::::::::::::::::::::::::::::::::::::::::::

    # ------------------------------------------------------------------
    def part_coarse_val_get_anisotropic(self, int ipart):
        """
           Get partition dimensions
        """
        # ************************************************************************
        # > Declaration
        cdef int          *agglomerationLines
        cdef int          *agglomerationLinesIdx
        cdef int           agglomerationLinesIdx_size
        cdef int          *isOnFineBnd
        # > For numpy capsule
        cdef NPY.npy_intp dim
        # ************************************************************************
        # > Get dim 
        dims = self.part_coarse_dim_get(ipart)

        # > Get array
        PDM_part_coarse_mesh_part_get_anisotropic_info(self._cmId,
                                                       ipart,       
                                                       &agglomerationLines,
                                                       &agglomerationLinesIdx,
                                                       &agglomerationLinesIdx_size,
                                                       &isOnFineBnd)
        print "\t\t\tPXI: \tpart_coarse_val_get_anisotropic: \tagglomerationLinesIdx_size", agglomerationLinesIdx_size
        # :::::::::::::::::::::::::::::::::::::::::::::::::::::
        # > Translate to numpy capsule (Tout est pas fait encore )
        if (agglomerationLines == NULL) :
          npAggloLines = None
        else :
          dim = <NPY.npy_intp> agglomerationLinesIdx[agglomerationLinesIdx_size-1]
          npAggloLines = NPY.PyArray_SimpleNewFromData(1,
                                                       &dim,
                                                       NPY.NPY_INT32,
                                                       <void *> agglomerationLines)

        if (agglomerationLinesIdx == NULL) :
          npAggloLinesIdx = None
        else :
          dim = <NPY.npy_intp> agglomerationLinesIdx_size
          npAggloLinesIdx = NPY.PyArray_SimpleNewFromData(1,
                                                          &dim,
                                                          NPY.NPY_INT32,
                                                          <void *> agglomerationLinesIdx)        

        if (isOnFineBnd == NULL) :
          npIsOnFineBnd = None
        else :
          dim = <NPY.npy_intp> dims['nCell']
          npIsOnFineBnd = NPY.PyArray_SimpleNewFromData(1,
                                                        &dim,
                                                        NPY.NPY_INT32,
                                                        <void *> isOnFineBnd)    

        # :::::::::::::::::::::::::::::::::::::::::::::::::::::
        return {'npAggloLines'    : npAggloLines,
                'npAggloLinesIdx' : npAggloLinesIdx,
                'npIsOnFineBnd'   : npIsOnFineBnd,
                } 
