
cdef extern from "pdm_overlay.h":

  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  # > Wrapping of enum
  ctypedef enum PDM_ol_mesh_t:
    PDM_OL_MESH_A = 0
    PDM_OL_MESH_B = 1

  ctypedef enum PDM_ol_parameter_t:
    PDM_OL_CAR_LENGTH_TOL = 0
    PDM_OL_EXTENTS_TOL    = 1
    PDM_OL_SAME_PLANE_TOL = 2

  ctypedef enum PDM_ol_mv_t:
    PDM_OL_MV_TRANSFORMATION = 0
    PDM_OL_MV_UNKNOWN        = 1
  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  # > Wrapping of function
  int PDM_ol_create(int          n_partMeshA,
                    int          n_partMeshB,
                    double       projectCoeff,
                    PDM_MPI_Comm comm)

  void PDM_ol_parameter_set(int                id,
                            PDM_ol_parameter_t parameter,
                            double             value)

  void PDM_ol_input_mesh_set(int           id,
                             PDM_ol_mesh_t mesh,
                             int           i_part,
                             int           n_face,
                             int          *face_vtx_idx,
                             int          *face_vtx,
                             PDM_g_num_t  *face_ln_to_gn,
                             int           n_vtx,
                             double       *coords,
                             PDM_g_num_t  *vtx_ln_to_gn)

  void PDM_ol_moving_type_set(int           id,
                              PDM_ol_mesh_t mesh,
                              PDM_ol_mv_t   mv)

  void PDM_ol_translation_set(int           id,
                              double       *vect,
                              double       *center)

  void PDM_ol_rotation_set(int      id,
                           double  *direction,
                           double  *center,
                           double   angle)

  void PDM_ol_compute(int          id)

  void PDM_ol_mesh_dim_get(int             id,
                           PDM_ol_mesh_t   mesh,
                           PDM_g_num_t    *nGOlFace,
                           PDM_g_num_t    *nGOlVtx)

  void PDM_ol_part_mesh_dim_get(int            id,
                                PDM_ol_mesh_t  mesh,
                                int            i_part,
                                int           *nOlFace,
                                int           *nOlLinkedFace,
                                int           *nOlVtx,
                                int           *sOlFaceIniVtx,
                                int           *sOlface_vtx,
                                int           *sInit_to_ol_face)

  void PDM_ol_mesh_entities_get(int              id,
                                PDM_ol_mesh_t    mesh,
                                int              i_part,
                                int            **olFaceIniVtxIdx,
                                int            **olFaceIniVtx,
                                int            **olface_vtx_idx,
                                int            **olface_vtx,
                                int            **olLinkedface_procIdx,
                                int            **olLinkedFace,
                                PDM_g_num_t    **olface_ln_to_gn,
                                double         **olCoords,
                                PDM_g_num_t    **olvtx_ln_to_gn,
                                int            **init_to_ol_face_idx,
                                int            **init_to_ol_face)

  void PDM_ol_del(int     id)

  void PDM_ol_dump_times(int     id)
  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# ------------------------------------------------------------------
cdef class Overlay:
  """
  """
  # ************************************************************************
  # > Class attributes
  cdef int _id
  cdef int _size
  cdef int _rank
  cdef int _init_face[2] # Store the init_face for mesh_a / mesh_b
  # ************************************************************************

  # ------------------------------------------------------------------------
  def __init__(self, int          n_partMeshA,
                     int          n_partMeshB,
                     double       projectCoeff,
                     MPI.Comm     comm):
    """
    """

    # ::::::::::::::::::::::::::::::::::::::::::::::::::
    self._rank = comm.Get_rank()
    self._size = comm.Get_size()
    # ::::::::::::::::::::::::::::::::::::::::::::::::::

    # ::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Convert mpi4py -> PDM_MPI
    cdef MPI.MPI_Comm c_comm = comm.ob_mpi
    cdef PDM_MPI_Comm PDMC   = PDM_MPI_mpi_2_pdm_mpi_comm(<void *> &c_comm)
    # ::::::::::::::::::::::::::::::::::::::::::::::::::

    self._id = PDM_ol_create(n_partMeshA,
                             n_partMeshB,
                             projectCoeff,
                             PDMC)

  # ------------------------------------------------------------------------
  def parameter_set(self, PDM_ol_parameter_t parameter,
                          double             value):
    """
    """
    PDM_ol_parameter_set(self._id, parameter, value)

  # ------------------------------------------------------------------------
  def input_mesh_set(self, PDM_ol_mesh_t                                 mesh,
                           int                                           i_part,
                           int                                           n_face,
                           NPY.ndarray[NPY.int32_t, mode='c', ndim=1]    face_vtx_idx,
                           NPY.ndarray[NPY.int32_t, mode='c', ndim=1]    face_vtx,
                           NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] face_ln_to_gn,
                           int                                           n_vtx,
                           NPY.ndarray[NPY.double_t  , mode='c', ndim=1] coords,
                           NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] vtx_ln_to_gn):
    """
    """
    self._init_face[<int>mesh] = n_face
    PDM_ol_input_mesh_set(self._id, mesh,
                          i_part,
                          n_face,
           <int*>         face_vtx_idx.data,
           <int*>         face_vtx.data,
           <PDM_g_num_t*> face_ln_to_gn.data,
                          n_vtx,
           <double*>      coords.data,
           <PDM_g_num_t*> vtx_ln_to_gn.data)

  # ------------------------------------------------------------------------
  def moving_type_set(self, PDM_ol_mesh_t mesh,
                            PDM_ol_mv_t   mv):
    """
    """
    PDM_ol_moving_type_set(self._id, mesh, mv)

  # ------------------------------------------------------------------------
  def translation_set(self, NPY.ndarray[NPY.double_t  , mode='c', ndim=1] vect,
                            NPY.ndarray[NPY.double_t  , mode='c', ndim=1] center):
    """
    """
    PDM_ol_translation_set(self._id,
                 <double*> vect.data,
                 <double*> center.data)

  # ------------------------------------------------------------------------
  def rotation_set(self, NPY.ndarray[NPY.double_t  , mode='c', ndim=1] direction,
                         NPY.ndarray[NPY.double_t  , mode='c', ndim=1] center,
                         double angle):
    """
    """
    PDM_ol_rotation_set(self._id,
              <double*> direction.data,
              <double*> center.data,
                        angle)

  # ------------------------------------------------------------------------
  def compute(self):
    """
    """
    PDM_ol_compute(self._id)

  # ------------------------------------------------------------------------
  def mesh_dim_get(self, PDM_ol_mesh_t mesh):
    """
    """
    # ************************************************************************
    # > Declaration
    cdef PDM_g_num_t nGOlFace
    cdef PDM_g_num_t nGOlVtx
    # ************************************************************************

    PDM_ol_mesh_dim_get(self._id, mesh, &nGOlFace, &nGOlVtx)

    return {'nGOlFace' : nGOlFace,
            'nGOlVtx'  : nGOlVtx
            }

  # ------------------------------------------------------------------------
  def part_mesh_dim_get(self, PDM_ol_mesh_t mesh,
                              int           i_part):
    """
    """
    # ************************************************************************
    # > Declaration
    cdef int nOlFace
    cdef int nOlLinkedFace
    cdef int nOlVtx
    cdef int sOlFaceIniVtx
    cdef int sOlface_vtx
    cdef int sInit_to_ol_face
    # ************************************************************************

    PDM_ol_part_mesh_dim_get(self._id, mesh, i_part,
                             &nOlFace,
                             &nOlLinkedFace,
                             &nOlVtx,
                             &sOlFaceIniVtx,
                             &sOlface_vtx,
                             &sInit_to_ol_face)

    return {'nOlFace'          : nOlFace,
            'nOlLinkedFace'    : nOlLinkedFace,
            'nOlVtx'           : nOlVtx,
            'sOlFaceIniVtx'    : sOlFaceIniVtx,
            'sOlface_vtx'      : sOlface_vtx,
            'sInit_to_ol_face' : sInit_to_ol_face
            }

  # ------------------------------------------------------------------------
  def part_mesh_entities_get(self, PDM_ol_mesh_t mesh,
                                   int           i_part):
    """
    """
    # ************************************************************************
    # > Declaration
    cdef int          nOlFace
    cdef int          nOlLinkedFace
    cdef int          nOlVtx
    cdef int          sOlFaceIniVtx
    cdef int          sOlface_vtx
    cdef int          sInit_to_ol_face
    cdef int         *olFaceIniVtxIdx
    cdef int         *olFaceIniVtx
    cdef int         *olface_vtx_idx
    cdef int         *olface_vtx
    cdef int         *olLinkedface_procIdx
    cdef int         *olLinkedFace
    cdef PDM_g_num_t *olface_ln_to_gn,
    cdef double      *olCoords
    cdef PDM_g_num_t *olvtx_ln_to_gn,
    cdef int         *init_to_ol_face_idx
    cdef int         *init_to_ol_face
    # ************************************************************************

    # > Get dims first
    PDM_ol_part_mesh_dim_get(self._id, mesh, i_part,
                             &nOlFace,
                             &nOlLinkedFace,
                             &nOlVtx,
                             &sOlFaceIniVtx,
                             &sOlface_vtx,
                             &sInit_to_ol_face)

    PDM_ol_mesh_entities_get(self._id, mesh, i_part,
                             &olFaceIniVtxIdx,
                             &olFaceIniVtx,
                             &olface_vtx_idx,
                             &olface_vtx,
                             &olLinkedface_procIdx,
                             &olLinkedFace,
                             &olface_ln_to_gn,
                             &olCoords,
                             &olvtx_ln_to_gn,
                             &init_to_ol_face_idx,
                             &init_to_ol_face)

    cdef NPY.npy_intp dim
    # > Build numpy capsule
    dim = <NPY.npy_intp> sInit_to_ol_face + 1
    np_olFaceIniVtxIdx = NPY.PyArray_SimpleNewFromData(1,
                                                       &dim,
                                                       NPY.NPY_INT32,
                                                       <void *> olFaceIniVtxIdx)

    # dim = <NPY.npy_intp> olFaceIniVtxIdx[sInit_to_ol_face]
    dim = <NPY.npy_intp> sInit_to_ol_face
    np_olFaceIniVtx = NPY.PyArray_SimpleNewFromData(1,
                                                    &dim,
                                                    NPY.NPY_INT32,
                                                    <void *> olFaceIniVtx)

    dim = <NPY.npy_intp> nOlFace + 1
    np_olface_vtx_idx = NPY.PyArray_SimpleNewFromData(1,
                                                      &dim,
                                                      NPY.NPY_INT32,
                                                      <void *> olface_vtx_idx)

    dim = <NPY.npy_intp> olface_vtx_idx[nOlFace]
    np_olface_vtx = NPY.PyArray_SimpleNewFromData(1,
                                                  &dim,
                                                  NPY.NPY_INT32,
                                                  <void *> olface_vtx)

    dim = <NPY.npy_intp> self._size + 1
    np_olLinkedface_procIdx = NPY.PyArray_SimpleNewFromData(1,
                                                            &dim,
                                                            NPY.NPY_INT32,
                                                            <void *> olLinkedface_procIdx)

    dim = <NPY.npy_intp> (4*nOlLinkedFace)
    np_olLinkedFace = NPY.PyArray_SimpleNewFromData(1,
                                                    &dim,
                                                    NPY.NPY_INT32,
                                                    <void *> olLinkedFace)

    dim = <NPY.npy_intp> nOlFace
    np_olface_ln_to_gn = NPY.PyArray_SimpleNewFromData(1,
                                                       &dim,
                                                       PDM_G_NUM_NPY_INT,
                                                       <void *> olface_ln_to_gn)
    dim = <NPY.npy_intp> 3 * nOlVtx
    np_olCoords = NPY.PyArray_SimpleNewFromData(1,
                                                &dim,
                                                NPY.NPY_DOUBLE,
                                                <void *> olCoords)

    dim = <NPY.npy_intp> nOlVtx
    np_olvtx_ln_to_gn = NPY.PyArray_SimpleNewFromData(1,
                                                      &dim,
                                                      PDM_G_NUM_NPY_INT,
                                                      <void *> olvtx_ln_to_gn)

    dim = <NPY.npy_intp> self._init_face[<int>mesh] + 1
    np_init_to_ol_face_idx = NPY.PyArray_SimpleNewFromData(1,
                                                       &dim,
                                                       NPY.NPY_INT32,
                                                       <void *> init_to_ol_face_idx)

    dim = <NPY.npy_intp> sInit_to_ol_face
    np_init_to_ol_face = NPY.PyArray_SimpleNewFromData(1,
                                                       &dim,
                                                       NPY.NPY_INT32,
                                                       <void *> init_to_ol_face)

    return {"olFaceIniVtxIdx"      : np_olFaceIniVtxIdx,
            "olFaceIniVtx"         : np_olFaceIniVtx,
            "olface_vtx_idx"       : np_olface_vtx_idx,
            "olface_vtx"           : np_olface_vtx,
            "olLinkedface_procIdx" : np_olLinkedface_procIdx,
            "olLinkedFace"         : np_olLinkedFace,
            "olface_ln_to_gn,"     : np_olface_ln_to_gn,
            "olCoords"             : np_olCoords,
            "olvtx_ln_to_gn,"      : np_olvtx_ln_to_gn,
            "init_to_ol_face_idx"  : np_init_to_ol_face_idx,
            "init_to_ol_face"      : np_init_to_ol_face,
            }

  # ------------------------------------------------------------------------
  def __dealloc__(self):
    """
       Use the free method of PDM Lib
    """
    # ************************************************************************
    # > Declaration
    # ************************************************************************
    PDM_ol_del(self._id)
