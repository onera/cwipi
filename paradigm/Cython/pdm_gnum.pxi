# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# > Wrapping of functions
cdef extern from "pdm_gnum.h":
  ctypedef struct PDM_gen_gnum_t:
    pass
  PDM_gen_gnum_t* PDM_gnum_create(const int             dim,
                                const int             n_part,
                                const PDM_bool_t      merge,
                                const double          tolerance,
                                const PDM_MPI_Comm    comm,
                                const PDM_ownership_t owner)

  void PDM_gnum_set_from_coords(PDM_gen_gnum_t* gen_gnum,
                                const int     i_part,
                                const int     n_elts,
                                const double *coords,
                                const double *char_length)

  void PDM_gnum_set_from_parents(PDM_gen_gnum_t *gen_gnum,
                                 const int             i_part,
                                 const int             n_elts,
                                 const PDM_g_num_t    *parent_gnum)

  void         PDM_gnum_compute(PDM_gen_gnum_t* gen_gnum)

  PDM_g_num_t*     PDM_gnum_get(PDM_gen_gnum_t* gen_gnum,
                                const int i_part)

  void            PDM_gnum_free(PDM_gen_gnum_t* gen_gnum)
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# > Class definition
cdef class GlobalNumbering:
  """

  """
  # --------------------------------------------------------------------------
  # > Class attributes
  cdef PDM_gen_gnum_t* _gen_gnum
  cdef NPY.npy_intp[:] _n_elem_per_part
  # --------------------------------------------------------------------------

  # --------------------------------------------------------------------------
  def __init__(self, int dim, int n_part, PDM_bool_t merge, double tolerance, MPI.Comm comm):
    """ Init a gnum structure """
    # ************************************************************************
    # > Declaration
    cdef MPI.MPI_Comm c_comm = comm.ob_mpi
    # ************************************************************************

    # ************************************************************************
    # > Init private array storing partition sizes
    self._n_elem_per_part = NPY.zeros(n_part, dtype=NPY.intp)
    # ************************************************************************

    # ************************************************************************
    # > PDM call
    self._gen_gnum = PDM_gnum_create(dim,
                                     n_part,
                                     merge,
                                     tolerance,
                                     PDM_MPI_mpi_2_pdm_mpi_comm (<void *> &c_comm),
                                     PDM_OWNERSHIP_USER) # Python take ownership);
    # ************************************************************************

  # --------------------------------------------------------------------------
  def gnum_set_from_coords(self,
                           int i_part,
                           int n_elts,
                           NPY.ndarray[NPY.double_t  , mode='c', ndim=1] coords not None,
                           NPY.ndarray[NPY.double_t  , mode='c', ndim=1] caracteristic_length):
    """
    """
    # ************************************************************************
    # > Declaration
    cdef double *caracteristic_length_data
    cdef double *coords_data
    # ************************************************************************

    # ************************************************************************
    coords_data = <double *> coords.data
    if (caracteristic_length is None):
      caracteristic_length_data = NULL
    else:
      caracteristic_length_data = <double *> caracteristic_length.data
    # ************************************************************************

    # ************************************************************************
    # > Store size to use it in the get
    self._n_elem_per_part[i_part] = n_elts
    # ************************************************************************

    # ************************************************************************
    # > PDM call
    PDM_gnum_set_from_coords(self._gen_gnum,
                             i_part,
                             n_elts,
                             coords_data,
                             caracteristic_length_data)
    # ************************************************************************

  # --------------------------------------------------------------------------
  def gnum_set_from_parent(self,
                           int i_part,
                           int n_elts,
                           NPY.ndarray[npy_pdm_gnum_t  , mode='c', ndim=1] parent_gnum not None):
    """
    """
    # ************************************************************************
    # > Declaration
    # ************************************************************************

    # ************************************************************************
    # > Store size to use it in the get
    self._n_elem_per_part[i_part] = n_elts
    # ************************************************************************

    # ************************************************************************
    # > PDM call
    PDM_gnum_set_from_parents(self._gen_gnum,
                              i_part,
                              n_elts,
               <PDM_g_num_t*> parent_gnum.data)
    # ************************************************************************

  # --------------------------------------------------------------------------
  def gnum_compute(self):
    """ Calls compute method from PDM_gnum """
    PDM_gnum_compute(self._gen_gnum)

  # --------------------------------------------------------------------------
  def gnum_get(self, int i_part):
    """ Calls the get method from PDM_gnum """
    # ************************************************************************
    # > Declaration
    cdef PDM_g_num_t *gnum_array
    # ************************************************************************

    # ************************************************************************
    # > PDM call
    gnum_array = PDM_gnum_get(self._gen_gnum, i_part)
    # ************************************************************************

    # ************************************************************************
    if (gnum_array == NULL):
      return None
    else:
      dim = <NPY.npy_intp> self._n_elem_per_part[i_part]
      np_gnum_array = NPY.PyArray_SimpleNewFromData(1,
                                                    &dim,
                                                    PDM_G_NUM_NPY_INT,
                                                    <void *> gnum_array)
    PyArray_ENABLEFLAGS(np_gnum_array, NPY.NPY_OWNDATA)
    return {'gnum' : np_gnum_array}
    # ************************************************************************

  # --------------------------------------------------------------------------
  def __dealloc__(self):
    """Calls the free method of PDM_gnum """
    PDM_gnum_free(self._gen_gnum);

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
