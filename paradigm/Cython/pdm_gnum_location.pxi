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
  # > Class attributes
  cdef public int _id
  
  cdef NPY.npy_intp[:] _nEltsInPart
  cdef NPY.npy_intp[:] _nEltsOutPart
  # --------------------------------------------------------------------------

  # --------------------------------------------------------------------------
  def __init__(self, int nPartIn, int nPartOut, MPI.Comm comm):
    """
        Init a gnum location structure
    """
    # ************************************************************************
    # > Declaration
    cdef MPI.MPI_Comm c_comm = comm.ob_mpi
    # ************************************************************************

    # ************************************************************************
    # > Init private array storing partition sizes
    self._nEltsInPart  = NPY.zeros(nPartIn,  dtype=NPY.intp)
    self._nEltsOutPart = NPY.zeros(nPartOut, dtype=NPY.intp)
    # ************************************************************************

    # ************************************************************************
    # > PDM call
    self._id = PDM_gnum_location_create(nPartIn,
                                        nPartOut,
                                        PDM_MPI_mpi_2_pdm_mpi_comm (<void *> &c_comm))
    # ************************************************************************

  # --------------------------------------------------------------------------
  def gnum_location_elements_set(self,
                                 int iPartIn,
                                 int nEltsIn,
                                 NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] gnumIn):
    """
       Calls set method for elements location from PDM_gnum_location
    """
    # ************************************************************************
    # > Declaration
    # ************************************************************************

    # ************************************************************************
    # > PDM call
    PDM_gnum_location_elements_set(self._id,
                                   iPartIn,
                                   nEltsIn,
                                   <PDM_g_num_t*> gnumIn.data)
    # ************************************************************************
    
    # ************************************************************************
    self._nEltsInPart[iPartIn] = nEltsIn
    # ************************************************************************

  # --------------------------------------------------------------------------
  def gnum_location_requested_elements_set(self,
                                           int iPartOut,
                                           int nEltsOut,
                                           NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] gnumOut):
    """
       Calls set method for requested elements location from PDM_gnum_location
    """
    # ************************************************************************
    # > Declaration
    # ************************************************************************

    # ************************************************************************
    # > PDM call
    PDM_gnum_location_requested_elements_set(self._id,
                                             iPartOut,
                                             nEltsOut,
                                             <PDM_g_num_t*> gnumOut.data)
    # ************************************************************************
    
    # ************************************************************************
    self._nEltsOutPart[iPartOut] = nEltsOut
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
                        int iPartOut):
    """
       Calls get method from PDM_gnum_location
    """
    # ************************************************************************
    # > Declaration
    # ************************************************************************
    cdef int *location_idx
    cdef int *location
    cdef int dim
    
    # ************************************************************************
    # > PDM call
    PDM_gnum_location_get(self._id,
                          iPartOut,
                          &location_idx,
                          &location)
    # ************************************************************************

    # ************************************************************************
    
    
    if (location_idx == NULL):
      locationIdx = None
    else:
      dim = <NPY.npy_intp> (self._nEltsOutPart[iPartOut] + 1)
      locationIdx = NPY.PyArray_SimpleNewFromData(1,
                                                  &dim,
                                                  NPY.NPY_INT32,
                                         <void *> location_idx)
      
    if (location == NULL):
      locationArr = None
    else:
      dim = <NPY.npy_intp> (location_idx[self._nEltsOutPart[iPartOut]])
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
