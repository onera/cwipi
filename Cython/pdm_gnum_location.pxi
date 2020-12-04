# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# > Wrapping of functions
cdef extern from "pdm_gnum_location.h":
    int                  PDM_gnum_location_create(const int          n_part_in,
                                                  const int          n_part_out,
                                                  const PDM_MPI_Comm comm);

    void           PDM_gnum_location_elements_set(const int         id,
                                                  const int         i_part_in,
                                                  const int         n_elts_in,
                                                  const PDM_g_num_t *gnum_in);

    void PDM_gnum_location_requested_elements_set(const int         id,
                                                  const int         i_part_out,
                                                  const int         n_elts_out,
                                                  const PDM_g_num_t *gnum_out);

    void                PDM_gnum_location_compute(const int id);

    void                    PDM_gnum_location_get(const int id,
                                                  const int i_part_out,
                                                  int **location_idx,
                                                  int **location);

    void                   PDM_gnum_location_free(const int id,
                                                  const int partial);
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# > Class definition
cdef class GlobalNumberingLocation:
  """
     GlobalNumberingLocation : Interface for pdm_gnum_location.c
  """
  # --------------------------------------------------------------------------
  # > Class attributesint _id
  cdef NPY.int32_t[:] _n_elmts_in_part
  cdef NPY.int32_t[:] _n_elmts_out_part
  # --------------------------------------------------------------------------

  # --------------------------------------------------------------------------
  def __init__(self, int n_part_in, int n_part_out, MPI.Comm comm):
    """
        Init a gnum location structure
    """
    # ************************************************************************
    # > Declaration
    cdef MPI.MPI_Comm c_comm = comm.ob_mpi
    # ************************************************************************
    # > Init private array storing partition sizes
    self._n_elmts_in_part  = NPY.zeros(n_part_in,  dtype=NPY.int32)
    self._n_elmts_out_part = NPY.zeros(n_part_out, dtype=NPY.int32)
    # ************************************************************************

    # ************************************************************************
    # > PDM call
    self._id = PDM_gnum_location_create(n_part_in,
                                        n_part_out,
                                        PDM_MPI_mpi_2_pdm_mpi_comm (<void *> &c_comm))
    # ************************************************************************

  # --------------------------------------------------------------------------
  def gnum_location_elements_set(self,
                                 int i_part_in,
                                 int n_elmts_in,
                                 NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] gnum_in):
    """
       Calls set method for elements location from PDM_gnum_location
    """
    # ************************************************************************
    # > Declaration
    # ************************************************************************

    # ************************************************************************
    # > PDM call
    PDM_gnum_location_elements_set(self._id,
                                   i_part_in,
                                   n_elmts_in,
                                   <PDM_g_num_t*> gnum_in.data)
    # ************************************************************************

    # ************************************************************************
    self._n_elmts_in_part[i_part_in] = n_elmts_in
    # ************************************************************************

  # --------------------------------------------------------------------------
  def gnum_location_requested_elements_set(self,
                                           int i_part_out,
                                           int n_elmts_out,
                                           NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] gnum_out):
    """
       Calls set method for requested elements location from PDM_gnum_location
    """
    # ************************************************************************
    # > Declaration
    # ************************************************************************

    # ************************************************************************
    # > PDM call
    PDM_gnum_location_requested_elements_set(self._id,
                                             i_part_out,
                                             n_elmts_out,
                              <PDM_g_num_t*> gnum_out.data)
    # ************************************************************************

    # ************************************************************************
    self._n_elmts_out_part[i_part_out] = n_elmts_out
    # ************************************************************************

  # --------------------------------------------------------------------------
  def gnum_location_compute(self):
    """
       Calls compute method from PDM_gnum_location
    """
    # ************************************************************************
    # > Declaration
    # ************************************************************************

    # ************************************************************************
    # > PDM call
    PDM_gnum_location_compute(self._id)
    # ************************************************************************

  # --------------------------------------------------------------------------
  def gnum_location_get(self,
                        int i_part_out):
    """
       Calls get method from PDM_gnum_location
    """
    # ************************************************************************
    # > Declaration
    # ************************************************************************
    cdef int          *location_idx
    cdef int          *location
    cdef NPY.npy_intp  dim
    # ************************************************************************
    # > PDM call
    PDM_gnum_location_get(self._id,
                          i_part_out,
                          &location_idx,
                          &location)
    # ************************************************************************

    # ************************************************************************
    if(location_idx == NULL):
      locationIdx = None
    else:
      dim = <NPY.npy_intp> (self._n_elmts_out_part[i_part_out] + 1)
      locationIdx = NPY.PyArray_SimpleNewFromData(1,
                                                  &dim,
                                                  NPY.NPY_INT32,
                                         <void *> location_idx)
    if(location == NULL):
      locationArr = None
    else:
      dim = <NPY.npy_intp> (location_idx[self._n_elmts_out_part[i_part_out]])
      locationArr = NPY.PyArray_SimpleNewFromData(1,
                                                  &dim,
                                                  NPY.NPY_INT32,
                                         <void *> location)
    # ************************************************************************

    # ************************************************************************
    return (locationIdx, locationArr)
    # ************************************************************************

  # --------------------------------------------------------------------------
  def __dealloc__(self):
    """
       Calls the free method of PDM_gnum_location
    """
    # ************************************************************************
    # > Declaration
    # ************************************************************************

    # ************************************************************************
    # > PDM call
    # Todo : tenir compte du partial ?
    PDM_gnum_location_free(self._id, 1)
    # ************************************************************************

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
