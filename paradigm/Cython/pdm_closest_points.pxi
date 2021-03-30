
cdef extern from "pdm_closest_points.h":

  int PDM_closest_points_create(PDM_MPI_Comm    comm,
                                int             n_closest,
                                PDM_ownership_t owner);

  void PDM_closest_points_n_part_cloud_set(int  id,
                                           int  n_part_cloud_src,
                                           int  n_part_cloud_tgt);

  void PDM_closest_points_tgt_cloud_set(int          id,
                                        int          i_part,
                                        int          n_points,
                                        double      *coords,
                                        PDM_g_num_t *gnum);

  void PDM_closest_points_src_cloud_set(int          id,
                                        int          i_part,
                                        int          n_points,
                                        double      *coords,
                                        PDM_g_num_t *gnum);

  void PDM_closest_points_compute(int id);

  void PDM_closest_points_get(int           id,
                              int           i_part_tgt,
                              PDM_g_num_t **closest_src_gnum,
                              double      **closest_src_distance);

  void PDM_closest_points_free(int id)

  void PDM_closest_points_dump_times(int id);


# ------------------------------------------------------------------
cdef class ClosestPoints:
  """
  """
  # ************************************************************************
  # > Class attributes
  cdef int _id
  cdef int _size
  cdef int _rank
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
    self._id = PDM_closest_points_create(PDMC, n_closest, PDM_OWNERSHIP_USER) # Python take ownership
    # ::::::::::::::::::::::::::::::::::::::::::::::::::

  # ------------------------------------------------------------------------
  def n_part_cloud_set(self, int n_part_cloud_src,
                             int n_part_cloud_tgt):
    """
    """
    self.tgt_n_points = <int *> malloc(sizeof(int *) * n_part_cloud_tgt )
    PDM_closest_points_n_part_cloud_set(self._id,
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
    PDM_closest_points_tgt_cloud_set(self._id,
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
    PDM_closest_points_src_cloud_set(self._id,
                                     i_part,
                                     n_points,
                           <double*> coords.data,
                      <PDM_g_num_t*> gnum.data)

  # ------------------------------------------------------------------------
  def compute(self):
    """
    """
    PDM_closest_points_compute(self._id)

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
    PDM_closest_points_get(self._id,
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
  def dump_times(self):
    """
    """
    PDM_closest_points_dump_times(self._id)

  # ------------------------------------------------------------------------
  def __dealloc__(self):
    """
    """
    free(self.tgt_n_points)
    PDM_closest_points_free(self._id)
