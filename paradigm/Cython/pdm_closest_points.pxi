
cdef extern from "pdm_closest_points.h":
  # > Wrapping of Ppart Structure
  ctypedef struct PDM_closest_point_t:
    pass
  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  PDM_closest_point_t* PDM_closest_points_create(PDM_MPI_Comm    comm,
                                int             n_closest,
                                PDM_ownership_t owner);

  void PDM_closest_points_n_part_cloud_set(PDM_closest_point_t* cls,
                                           int                  n_part_cloud_src,
                                           int                  n_part_cloud_tgt);

  void PDM_closest_points_tgt_cloud_set(PDM_closest_point_t *cls,
                                        int                  i_part,
                                        int                  n_points,
                                        double              *coords,
                                        PDM_g_num_t         *gnum);

  void PDM_closest_points_src_cloud_set(PDM_closest_point_t *cls,
                                        int                  i_part,
                                        int                  n_points,
                                        double              *coords,
                                        PDM_g_num_t         *gnum);

  void PDM_closest_points_compute(PDM_closest_point_t  *cls);

  void PDM_closest_points_get(PDM_closest_point_t  *cls,
                              int                   i_part_tgt,
                              PDM_g_num_t         **closest_src_gnum,
                              double              **closest_src_distance);

  void PDM_closest_points_tgt_in_src_get(PDM_closest_point_t  *cls,
                                         int                   i_part_src,
                                         int                 **tgt_in_src_idx,
                                         PDM_g_num_t         **tgt_in_src);

  void PDM_closest_points_free(PDM_closest_point_t *cls)

  void PDM_closest_points_dump_times(PDM_closest_point_t *cls);


# ------------------------------------------------------------------
cdef class ClosestPoints:
  """
  """
  # ************************************************************************
  # > Class attributes
  cdef PDM_closest_point_t *_cls
  cdef int _size
  cdef int _rank
  cdef int* src_n_points
  cdef int* tgt_n_points
  cdef int n_closest
  # ************************************************************************

  # ------------------------------------------------------------------------
  def __init__(self, MPI.Comm comm,
                     int      n_closest):
    """
    """
    # ::::::::::::::::::::::::::::::::::::::::::::::::::
    self.n_closest = n_closest
    # ::::::::::::::::::::::::::::::::::::::::::::::::::

    # ::::::::::::::::::::::::::::::::::::::::::::::::::
    self._rank = comm.Get_rank()
    self._size = comm.Get_size()
    # ::::::::::::::::::::::::::::::::::::::::::::::::::

    # ::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Convert mpi4py -> PDM_MPI
    cdef MPI.MPI_Comm c_comm = comm.ob_mpi
    cdef PDM_MPI_Comm PDMC   = PDM_MPI_mpi_2_pdm_mpi_comm(<void *> &c_comm)
    # ::::::::::::::::::::::::::::::::::::::::::::::::::

    # ::::::::::::::::::::::::::::::::::::::::::::::::::
    self.src_n_points = NULL
    self.tgt_n_points = NULL
    # ::::::::::::::::::::::::::::::::::::::::::::::::::

    # ::::::::::::::::::::::::::::::::::::::::::::::::::
    self._cls = PDM_closest_points_create(PDMC, n_closest, PDM_OWNERSHIP_USER) # Python take ownership
    # ::::::::::::::::::::::::::::::::::::::::::::::::::

  # ------------------------------------------------------------------------
  def n_part_cloud_set(self, int n_part_cloud_src,
                             int n_part_cloud_tgt):
    """
    """
    assert(self.src_n_points  == NULL)
    assert(self.tgt_n_points  == NULL)
    self.src_n_points = <int *> malloc(sizeof(int *) * n_part_cloud_src )
    self.tgt_n_points = <int *> malloc(sizeof(int *) * n_part_cloud_tgt )
    PDM_closest_points_n_part_cloud_set(self._cls,
                                        n_part_cloud_src,
                                        n_part_cloud_tgt)

  # ------------------------------------------------------------------------
  def tgt_cloud_set(self, int                                           i_part,
                          int                                           n_points,
                          NPY.ndarray[NPY.double_t  , mode='c', ndim=1] coords,
                          NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] gnum):
    """
    """
    self.tgt_n_points[i_part] = n_points
    PDM_closest_points_tgt_cloud_set(self._cls,
                                     i_part,
                                     n_points,
                           <double*> coords.data,
                      <PDM_g_num_t*> gnum.data)


  # ------------------------------------------------------------------------
  def src_cloud_set(self, int                                           i_part,
                          int                                           n_points,
                          NPY.ndarray[NPY.double_t  , mode='c', ndim=1] coords,
                          NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] gnum):
    """
    """
    self.src_n_points[i_part] = n_points
    PDM_closest_points_src_cloud_set(self._cls,
                                     i_part,
                                     n_points,
                           <double*> coords.data,
                      <PDM_g_num_t*> gnum.data)

  # ------------------------------------------------------------------------
  def compute(self):
    """
    """
    PDM_closest_points_compute(self._cls)

  # ------------------------------------------------------------------------
  def points_get(self, int i_part_tgt):
    """
    """
    # ************************************************************************
    # > Declaration
    cdef PDM_g_num_t *closest_src_gnum
    cdef double      *closest_src_distance
    # ************************************************************************

    # > Get
    PDM_closest_points_get(self._cls,
                           i_part_tgt,
                           &closest_src_gnum,
                           &closest_src_distance)

    cdef NPY.npy_intp dim

    # > Build numpy capsule
    dim = <NPY.npy_intp> self.tgt_n_points[i_part_tgt] * self.n_closest
    np_closest_src_gnum = NPY.PyArray_SimpleNewFromData(1,
                                                        &dim,
                                                        PDM_G_NUM_NPY_INT,
                                                        <void *> closest_src_gnum)
    PyArray_ENABLEFLAGS(np_closest_src_gnum, NPY.NPY_OWNDATA);

    np_closest_src_distance = NPY.PyArray_SimpleNewFromData(1,
                                                            &dim,
                                                            NPY.NPY_DOUBLE,
                                                            <void *> closest_src_distance)
    PyArray_ENABLEFLAGS(np_closest_src_distance, NPY.NPY_OWNDATA);

    return {'closest_src_gnum'  : np_closest_src_gnum,
            'closest_src_distance' : np_closest_src_distance
            }

  # ------------------------------------------------------------------------
  def tgt_in_src_get(self, int i_part_src):
    """
    """
    # ************************************************************************
    # > Declaration
    cdef int *tgt_in_src_idx
    cdef PDM_g_num_t *tgt_in_src
    # ************************************************************************

    # > Get
    PDM_closest_points_tgt_in_src_get(self._cls,
                                      i_part_src,
                                      &tgt_in_src_idx,
                                      &tgt_in_src)

    cdef NPY.npy_intp dim

    # > Build numpy capsule
    dim = <NPY.npy_intp> self.src_n_points[i_part_src] + 1
    np_tgt_in_src_idx = NPY.PyArray_SimpleNewFromData(1,
                                                      &dim,
                                                      NPY.NPY_INT32,
                                                      <void *> tgt_in_src_idx)
    PyArray_ENABLEFLAGS(np_tgt_in_src_idx, NPY.NPY_OWNDATA);

    dim = <NPY.npy_intp> np_tgt_in_src_idx[self.src_n_points[i_part_src]]
    np_tgt_in_src = NPY.PyArray_SimpleNewFromData(1,
                                                  &dim,
                                                  PDM_G_NUM_NPY_INT,
                                                  <void *> tgt_in_src)
    PyArray_ENABLEFLAGS(np_tgt_in_src, NPY.NPY_OWNDATA);

    return {'tgt_in_src_idx' : np_tgt_in_src_idx,
            'tgt_in_src'     : np_tgt_in_src
            }

  # ------------------------------------------------------------------------
  def dump_times(self):
    """
    """
    PDM_closest_points_dump_times(self._cls)

  # ------------------------------------------------------------------------
  def __dealloc__(self):
    """
    """
    free(self.tgt_n_points)
    free(self.src_n_points)
    PDM_closest_points_free(self._cls)
