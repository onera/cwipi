
cdef extern from "pdm_points_merge.h":
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Wrapping of structure
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Wrapping of function
    int PDM_points_merge_create(int          n_point_cloud,
                                double       tolerance,
                                PDM_MPI_Comm comm);
    void PDM_points_merge_free(int          id);
    void PDM_points_merge_cloud_set(int          id,
                                    int          i_point_cloud,
                                    int          n_points,
                                    double      *coords,
                                    double      *char_length);
    void PDM_points_merge_process(int        id);
    void PDM_points_merge_candidates_get(int    id,
                                         int    i_point_cloud,
                                         int    **candidates_idx,
                                         int    **candidates_desc);


    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# ------------------------------------------------------------------
cdef class PointsMerge:
    """
       PointsMerge: Interface to build connection between multiple cloud in parallel
       Useful for get connection between partiton from faces coordinates
    """
    # ************************************************************************
    # > Class attributes
    cdef int _id
    cdef int _size
    cdef int _rank
    # ************************************************************************
    # ------------------------------------------------------------------------
    def __init__(self, MPI.Comm    comm,
                       int         n_point_cloud,
                       double      relative_tolerance):
        """
        TODOUX
        """
        # ************************************************************************
        # > Declaration
        # cdef int      nElts
        # cdef int      idx
        # # > Numpy array
        # cdef NPY.ndarray[npy_pdm_gnum_t, ndim=1, mode='fortran'] partLNToGN
        # ************************************************************************

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
        self._id = PDM_points_merge_create(n_point_cloud, relative_tolerance, PDMC)
        # ::::::::::::::::::::::::::::::::::::::::::::::::::

    # ------------------------------------------------------------------------
    def cloud_set(self, int i_point_cloud,
                        int n_points,
                        NPY.ndarray[NPY.double_t, mode='c', ndim=1] coords,
                        NPY.ndarray[NPY.double_t, mode='c', ndim=1] char_length):
        """
        """
        # ************************************************************************
        # > Declaration
        # ************************************************************************

        # ::::::::::::::::::::::::::::::::::::::::::::::::::
        PDM_points_merge_cloud_set(self._id,
                                   i_point_cloud,
                                   n_points,
                                   <double *> coords.data,
                                   <double *> char_length.data)
        # ::::::::::::::::::::::::::::::::::::::::::::::::::

    # ------------------------------------------------------------------------
    def compute(self):
        """
        """
        # ************************************************************************
        # > Declaration
        # ************************************************************************
        PDM_points_merge_process(self._id)

    # ------------------------------------------------------------------------
    def get_merge_candidates(self, int i_point_cloud):
        """
        """
        # ************************************************************************
        # > Declaration
        cdef int *candidates_idx
        cdef int *candidates_desc
        # ************************************************************************

        # > Get Size
        PDM_points_merge_candidates_get(self._id,
                                        i_point_cloud,
                                        &candidates_idx,
                                        &candidates_desc)

        # > Build numpy capsule
        # dim = <NPY.npy_intp> 2 * dNface
        # npFaceCell = NPY.PyArray_SimpleNewFromData(1,
        #                                            &dim,
        #                                            PDM_G_NUM_NPY_INT,
        #                                            <void *> faceCell)
        np_candidates_idx  = None
        np_candidates_desc = None

        return {'candidates_idx'  : np_candidates_idx,
                'candidates_desc' : np_candidates_desc
                }

    # ------------------------------------------------------------------------
    def __dealloc__(self):
      """
         Use the free method of PDM Lib
      """
      # ************************************************************************
      # > Declaration
      # ************************************************************************
      print 'PDM_points_merge_free'
      PDM_points_merge_free(self._id)

