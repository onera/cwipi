import warnings

cdef extern from "pdm_part_to_part.h":
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Wrapping of Ppart Structure
    ctypedef struct PDM_part_to_part_t:
        pass
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    PDM_part_to_part_t *PDM_part_to_part_create(PDM_g_num_t   **gnum_elt1,
                                                int            *n_elt1,
                                                int             n_part1,
                                                PDM_g_num_t   **gnum_elt2,
                                                int            *n_elt2,
                                                int             n_part2,
                                                int           **part1_to_part2_idx,
                                                PDM_g_num_t   **part1_to_part2,
                                                PDM_MPI_Comm    comm);
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    void PDM_part_to_part_iexch(PDM_part_to_part_t                *ptp,
                                const PDM_mpi_comm_kind_t          k_comm,
                                const PDM_stride_t                 t_stride,
                                const PDM_part_to_part_data_def_t  t_part1_data_def,
                                const int                          cst_stride,
                                const size_t                       s_data,
                                const int                        **part1_stride,
                                const void                       **part1_data,
                                int                             ***part2_stride,
                                void                            ***part2_data,
                                int                               *request);
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    void PDM_part_to_part_iexch_wait(PDM_part_to_part_t                *ptp,
                                     const int                          request);
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    void PDM_part_to_part_gnum1_come_from_get(PDM_part_to_part_t *ptp,
                                              int              ***gnum1_come_from_idx,
                                              PDM_g_num_t      ***gnum1_come_from);
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    void PDM_part_to_part_ref_lnum2_get(PDM_part_to_part_t   *ptp,
                                        int                 **n_ref_lnum2,
                                        int                ***ref_lnum2);
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    void PDM_part_to_part_reverse_iexch(PDM_part_to_part_t             *ptp,
                                        PDM_mpi_comm_kind_t             k_comm,
                                        PDM_stride_t                    t_stride,
                                        PDM_part_to_part_data_def_t     t_part2_data_def,
                                        int                             cst_stride,
                                        size_t                          s_data,
                                        int                           **part2_stride,
                                        void                          **part2_data,
                                        int                          ***part1_stride,
                                        void                         ***part1_data,
                                        int                            *request);
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    void PDM_part_to_part_reverse_iexch_wait(PDM_part_to_part_t                *ptp,
                                             const int                          request);
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    PDM_part_to_part_t *PDM_part_to_part_free(PDM_part_to_part_t *ptp);
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# ------------------------------------------------------------------
cdef class PartToPart:
    """
       PartToPart: Interface for block_to_part.c
    """
    # ************************************************************************
    # > Class attributes
    cdef public:
        object lpart1_ln_to_gn, lpart2_ln_to_gn, lpart1_to_part2_idx, lpart1_to_part2
    cdef PDM_part_to_part_t         *ptp
    cdef PDM_g_num_t               **gnum_elt1
    cdef int*                        n_elt1
    cdef int                         n_part1
    cdef PDM_g_num_t               **gnum_elt2
    cdef int*                        n_elt2
    cdef int                         n_part2
    cdef MPI.Comm                    py_comm
    cdef PDM_g_num_t               **_part1_ln_to_gn
    cdef PDM_g_num_t               **_part2_ln_to_gn
    cdef int                       **_part1_to_part2_idx
    cdef PDM_g_num_t               **_part1_to_part2
    cdef cpp_map[int, void**]        dict_part_data
    cdef cpp_map[int, int **]        dict_part_stride
    cdef cpp_map[int, size_t]        dict_npy_type
    cdef cpp_map[int, PDM_stride_t]  dict_stride
    cdef cpp_map[int, int         ]  dict_cst_strid
    # ************************************************************************
    # ------------------------------------------------------------------------
    def __cinit__(self, MPI.Comm comm,
                        list part1_ln_to_gn,
                        list part2_ln_to_gn,
                        list part1_to_part2_idx,
                        list part1_to_part2):
      """
      """
      self.n_part1 = len(part1_ln_to_gn)
      self.n_part2 = len(part2_ln_to_gn)

      assert(len(part1_to_part2_idx) == self.n_part1)

      self.lpart1_ln_to_gn     = part1_ln_to_gn
      self.lpart2_ln_to_gn     = part2_ln_to_gn
      self.lpart1_to_part2_idx = part1_to_part2_idx
      self.lpart1_to_part2     = part1_to_part2

      # > Convert input data
      cdef MPI.MPI_Comm c_comm   = comm.ob_mpi
      cdef PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&c_comm)

      self.n_elt1  = list_to_int_pointer([array.size for array in part1_ln_to_gn])
      self.n_elt2  = list_to_int_pointer([array.size for array in part2_ln_to_gn])

      self._part1_ln_to_gn     = np_list_to_gnum_pointers(part1_ln_to_gn)
      self._part2_ln_to_gn     = np_list_to_gnum_pointers(part2_ln_to_gn)
      self._part1_to_part2_idx = np_list_to_int_pointers (part1_to_part2_idx)
      self._part1_to_part2     = np_list_to_gnum_pointers(part1_to_part2)

      self.ptp = PDM_part_to_part_create(<const PDM_g_num_t **> self._part1_ln_to_gn,
                                         <const int *>          self.n_elt1,
                                                                self.n_part1,
                                         <const PDM_g_num_t **> self._part2_ln_to_gn,
                                         <const int *>          self.n_elt2,
                                                                self.n_part2,
                                         <const int **>         self._part1_to_part2_idx,
                                         <const PDM_g_num_t **> self._part1_to_part2,
                                         pdm_comm);


    # ------------------------------------------------------------------------
    def get_referenced_lnum2(self):
      """ Return a copy of the local referenced ids for each part2 partition """
      cdef int  *n_ref_lnum2 = NULL;
      cdef int **ref_lnum2   = NULL;
      PDM_part_to_part_ref_lnum2_get(self.ptp,
                                     &n_ref_lnum2,
                                     &ref_lnum2);
      lnp_ref_lnum2 = list()
      for i_part in range(self.n_part2):
        np_ref_lnum2_id = create_numpy_i(ref_lnum2[i_part], n_ref_lnum2[i_part], flag_owndata=False)
        lnp_ref_lnum2.append(NPY.copy(np_ref_lnum2_id))
      return lnp_ref_lnum2

    # ------------------------------------------------------------------------
    def get_gnum1_come_from(self):
      """ Return a copy of the sending part1 gnum for each part2 partition """
      # Needed to get the sizes n_ref_lnum2
      cdef int  *n_ref_lnum2 = NULL;
      cdef int **ref_lnum2   = NULL;
      PDM_part_to_part_ref_lnum2_get(self.ptp,
                                     &n_ref_lnum2,
                                     &ref_lnum2);

      cdef int         **gnum1_come_from_idx = NULL
      cdef PDM_g_num_t **gnum1_come_from     = NULL
      PDM_part_to_part_gnum1_come_from_get(self.ptp,
                                           &gnum1_come_from_idx,
                                           &gnum1_come_from)
      come_from_l = list()
      for i_part in range(self.n_part2):
        n_ref = n_ref_lnum2[i_part]
        # Idx array is shaped as n_ref_lnum2
        sending_idx = create_numpy_i(gnum1_come_from_idx[i_part], n_ref+1, flag_owndata=False)
        sending_array = create_numpy_pdm_gnum(gnum1_come_from[i_part], gnum1_come_from_idx[i_part][n_ref], flag_owndata=False)
        come_from_l.append({'come_from_idx' : NPY.copy(sending_idx), 'come_from' : NPY.copy(sending_array)})
      return come_from_l

    # ------------------------------------------------------------------------
    def iexch(self,
              PDM_mpi_comm_kind_t         k_comm,
              PDM_part_to_part_data_def_t t_part1_data_def,
              list                        part1_data,
              part1_stride=1,
              bint interlaced_str=True):
      """
      """
      cdef int request_exch
      cdef PDM_stride_t _stride_t

      cdef int   _part1_stride_cst = 0
      cdef int** _part1_stride     = NULL
      if isinstance(part1_stride, int):
        _stride_t = PDM_STRIDE_CST_INTERLACED if interlaced_str else PDM_STRIDE_CST_INTERLEAVED
        _part1_stride_cst = part1_stride
      elif isinstance(part1_stride, list):
        _stride_t = PDM_STRIDE_VAR_INTERLACED
        assert len(part1_stride) == self.n_part1
        for i_part in range(self.n_part1):
          assert_single_dim_np(part1_stride[i_part], NPY.int32, self.n_elt1[i_part])
        _part1_stride = np_list_to_int_pointers(part1_stride)
      else:
        raise ValueError("Invalid stride in PtB exchange")

      cdef void** _part1_data   = np_list_to_void_pointers(part1_data)

      cdef size_t s_data   = part1_data[0].dtype.itemsize
      cdef size_t npy_type = part1_data[0].dtype.num

      cdef int**  _part2_stride = NULL;
      cdef void** _part2_data   = NULL;
      PDM_part_to_part_iexch(self.ptp,
                             k_comm,
                             _stride_t,
                             t_part1_data_def,
                             _part1_stride_cst,
                             s_data,
            <const int ** >  _part1_stride,
            <const void **>  _part1_data,
                             &_part2_stride,
                  <void ***> &_part2_data,
                             &request_exch)

      self.dict_part_data  [request_exch] = _part2_data
      self.dict_part_stride[request_exch] = _part2_stride
      self.dict_npy_type   [request_exch] = npy_type
      self.dict_stride     [request_exch] = _stride_t
      self.dict_cst_strid  [request_exch] = _part1_stride_cst

      if _stride_t == PDM_STRIDE_VAR_INTERLACED:
        free(_part1_stride)
      free(_part1_data)

      return request_exch

    # ------------------------------------------------------------------------
    def wait(self, int request_id):
      """
      Wait for a preceding exchange
      """
      PDM_part_to_part_iexch_wait(self.ptp, request_id)

      cdef int         **gnum1_come_from_idx = NULL
      cdef PDM_g_num_t **gnum1_come_from     = NULL
      PDM_part_to_part_gnum1_come_from_get(self.ptp,
                                           &gnum1_come_from_idx,
                                           &gnum1_come_from)

      cdef int  *n_ref_lnum2 = NULL;
      cdef int **ref_lnum2   = NULL;
      PDM_part_to_part_ref_lnum2_get(self.ptp,
                                     &n_ref_lnum2,
                                     &ref_lnum2);

      cdef void **_part2_data   = self.dict_part_data  [request_id]
      cdef int  **_part2_stride = self.dict_part_stride[request_id]
      cdef NPY.npy_intp dim_np

      # cdef int *__part2_stride = NULL

      lnp_part_strid = list()
      lnp_part_data  = list()
      for i_part in range(self.n_part2):

        if(_part2_stride != NULL):

          strid_size = gnum1_come_from_idx[i_part][n_ref_lnum2[i_part]]

          np_part2_stride = create_numpy_i(_part2_stride[i_part], strid_size)
          dim_np = np_part2_stride.sum()

          np_part2_data = NPY.PyArray_SimpleNewFromData(1, &dim_np, self.dict_npy_type[request_id], <void *> _part2_data[i_part])
          PyArray_ENABLEFLAGS(np_part2_data, NPY.NPY_OWNDATA)

          lnp_part_strid.append(np_part2_stride)
          lnp_part_data .append(np_part2_data)

        else:
          dim_np  = gnum1_come_from_idx[i_part][n_ref_lnum2[i_part]]
          np_part2_data = NPY.PyArray_SimpleNewFromData(1, &dim_np, self.dict_npy_type[request_id], <void *> _part2_data[i_part])
          PyArray_ENABLEFLAGS(np_part2_data, NPY.NPY_OWNDATA)

          lnp_part_data .append(np_part2_data)

      free(self.dict_part_data  [request_id])
      free(self.dict_part_stride[request_id])

      self.dict_part_data  .erase(request_id)
      self.dict_part_stride.erase(request_id)
      self.dict_npy_type   .erase(request_id)
      self.dict_stride     .erase(request_id)
      self.dict_cst_strid  .erase(request_id)

      return lnp_part_strid, lnp_part_data

    # ------------------------------------------------------------------------
    def reverse_iexch(self,
                      PDM_mpi_comm_kind_t         k_comm,
                      PDM_part_to_part_data_def_t t_part2_data_def,
                      list                        part2_data,
                      part2_stride=1,
                      bint interlaced_str=True):
      """
      """
      cdef NPY.ndarray  np_part2_data
      cdef int          request_exch
      cdef PDM_stride_t _stride_t

      cdef int   _part2_stride_cst = 0
      cdef int** _part2_stride = NULL
      if isinstance(part2_stride, int):
        _stride_t = PDM_STRIDE_CST_INTERLACED if interlaced_str else PDM_STRIDE_CST_INTERLEAVED
        _part2_stride_cst = part2_stride
      elif isinstance(part2_stride, list):
        _stride_t = PDM_STRIDE_VAR_INTERLACED
        assert len(part2_stride) == self.n_part2
        for i_part in range(self.n_part2):
          assert_single_dim_np(part2_stride[i_part], NPY.int32, self.n_elt2[i_part])
        _part2_stride = np_list_to_int_pointers(part2_stride)
      else:
        raise ValueError("Invalid stride in PtB exchange")

      cdef void** _part2_data   = np_list_to_void_pointers(part2_data)

      np_part2_data        = part2_data[0]
      cdef size_t s_data   = np_part2_data.dtype.itemsize
      cdef size_t npy_type = part2_data[0].dtype.num

      cdef int**  _part1_stride = NULL;
      cdef void** _part1_data   = NULL;
      PDM_part_to_part_reverse_iexch(self.ptp,
                                     k_comm,
                                     _stride_t,
                                     t_part2_data_def,
                                     _part2_stride_cst,
                                     s_data,
                    <const int ** >  _part2_stride,
                    <const void **>  _part2_data,
                                     &_part1_stride,
                          <void ***> &_part1_data,
                                     &request_exch)

      self.dict_part_data  [request_exch] = _part1_data
      self.dict_part_stride[request_exch] = _part1_stride
      self.dict_npy_type   [request_exch] = npy_type
      self.dict_stride     [request_exch] = _stride_t
      self.dict_cst_strid  [request_exch] = _part2_stride_cst

      if _stride_t == PDM_STRIDE_VAR_INTERLACED:
        free(_part2_stride)
      free(_part2_data)

      return request_exch

    # ------------------------------------------------------------------------
    def reverse_wait(self, int request_id):
      """
      Wait for a preceding exchange
      """
      PDM_part_to_part_reverse_iexch_wait(self.ptp, request_id)

      cdef void **_part1_data   = self.dict_part_data  [request_id]
      cdef int  **_part1_stride = self.dict_part_stride[request_id]
      cdef NPY.npy_intp dim_np

      lnp_part_strid = list()
      lnp_part_data  = list()
      for i_part in range(self.n_part1):

        if(_part1_stride != NULL):

          strid_size = self.n_elt1[i_part]

          np_part1_stride = create_numpy_i(_part1_stride[i_part], strid_size)
          dim_np = np_part1_stride.sum()
          print("dim_np : ", dim_np)
          print("np_part1_stride : ", np_part1_stride)

          np_part1_data = NPY.PyArray_SimpleNewFromData(1, &dim_np, self.dict_npy_type[request_id], <void *> _part1_data[i_part])
          PyArray_ENABLEFLAGS(np_part1_data, NPY.NPY_OWNDATA)

          lnp_part_strid.append(np_part1_stride)
          lnp_part_data .append(np_part1_data)

        else:
          dim_np  = self.n_elt1[i_part] * self.dict_cst_strid[request_id]

          print("dim_np : ", dim_np)

          np_part1_data = NPY.PyArray_SimpleNewFromData(1, &dim_np, self.dict_npy_type[request_id], <void *> _part1_data[i_part])
          PyArray_ENABLEFLAGS(np_part1_data, NPY.NPY_OWNDATA)

          lnp_part_data .append(np_part1_data)

      free(self.dict_part_data  [request_id])

      if(_part1_stride != NULL):
        free(self.dict_part_stride[request_id])

      self.dict_part_data  .erase(request_id)
      self.dict_part_stride.erase(request_id)
      self.dict_npy_type   .erase(request_id)
      self.dict_stride     .erase(request_id)
      self.dict_cst_strid  .erase(request_id)

      return lnp_part_strid, lnp_part_data

    # ------------------------------------------------------------------------
    def __dealloc__(self):
      """ Free structure """
      PDM_part_to_part_free(self.ptp)

      # > Free allocated array
      free(self.n_elt1 )
      free(self.n_elt2 )
      free(self._part1_ln_to_gn    )
      free(self._part2_ln_to_gn    )
      free(self._part1_to_part2_idx)
      free(self._part1_to_part2    )
