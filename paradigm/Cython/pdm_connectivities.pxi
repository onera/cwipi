cdef extern from "pdm_dconnectivity_transform.h":
    # ------------------------------------------------------------------
    void PDM_dconnectivity_transpose(const PDM_MPI_Comm     comm,
                                     const PDM_g_num_t     *entity1_distrib,
                                           PDM_g_num_t     *entity2_distrib,
                                     const int             *dentity1_entity2_idx,
                                     const PDM_g_num_t     *dentity1_entity2,
                                           int              is_signed,
                                           int            **dentity2_entity1_idx,
                                           PDM_g_num_t    **dentity2_entity1)

cdef extern from "pdm_part_connectivity_transform.h":
    # ------------------------------------------------------------------
    void PDM_combine_connectivity(int   n_entity1,
                                  int  *entity1_entity2_idx,
                                  int  *entity1_entity2,
                                  int  *entity2_entity3_idx,
                                  int  *entity2_entity3,
                                  int **entity1_entity3_idx,
                                  int **entity1_entity3)
    # ------------------------------------------------------------------
    void PDM_connectivity_transpose(const int   n_entity1,
                                    const int   n_entity2,
                                          int  *entity1_entity2_idx,
                                          int  *entity1_entity2,
                                          int **entity2_entity1_idx,
                                          int **entity2_entity1)
    # ------------------------------------------------------------------
    void PDM_part_connectivity_transpose(const int    n_part,
                                         const int   *n_entity1,
                                         const int   *n_entity2,
                                               int  **entity1_entity2_idx,
                                               int  **entity1_entity2,
                                               int ***entity2_entity1_idx,
                                               int ***entity2_entity1)
    # ------------------------------------------------------------------
    void PDM_part_connectivity_to_connectity_idx(const int    n_part,
                                                const int   *n_entity1,
                                                      int  **entity1_entity2_in,
                                                      int ***entity1_entity2_idx,
                                                      int ***entity1_entity2)
    # ------------------------------------------------------------------


# ------------------------------------------------------------------------
def dconnectivity_transpose(MPI.Comm comm,
                            NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1]    entity1_distrib,
                            NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1]    entity2_distrib,
                            NPY.ndarray[int           , mode='c', ndim=1]    dentity1_entity2_idx,
                            NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1]    dentity1_entity2,
                            bint                                             is_signed):
                                  
    ## entity2_distrib can be recomputed if allocated but entity2_distrib[0] = -1
    # TODO : manage case entity1_distrib recomputed
    # > Convert mpi4py -> PDM_MPI
    cdef MPI.MPI_Comm c_comm = comm.ob_mpi
    cdef PDM_MPI_Comm PDMC   = PDM_MPI_mpi_2_pdm_mpi_comm(&c_comm)

    cdef PDM_g_num_t * _entity1_distrib
    cdef PDM_g_num_t * _entity2_distrib
    _entity1_distrib = <PDM_g_num_t *> entity1_distrib.data if entity1_distrib is not None else NULL
    _entity2_distrib = <PDM_g_num_t *> entity2_distrib.data if entity2_distrib is not None else NULL

    cdef int *_dentity1_entity2_idx
    _dentity1_entity2_idx = <int *> dentity1_entity2_idx.data if dentity1_entity2_idx is not None else NULL
    cdef PDM_g_num_t *_dentity1_entity2
    _dentity1_entity2 = <PDM_g_num_t *> dentity1_entity2.data if dentity1_entity2 is not None else NULL


    cdef        int * _dentity2_entity1_idx
    cdef PDM_g_num_t* _dentity2_entity1

    PDM_dconnectivity_transpose(PDMC,
                                _entity1_distrib,
                                _entity2_distrib,
                                _dentity1_entity2_idx,
                                _dentity1_entity2,
                         <int>  is_signed,
                              &_dentity2_entity1_idx,
                              &_dentity2_entity1)
    
    dn_entity2 = entity2_distrib[comm.Get_rank()+1] - entity2_distrib[comm.Get_rank()]

    np_dentity2_entity1_idx = create_numpy_i(_dentity2_entity1_idx, dn_entity2 + 1)

    np_dentity2_entity1 = create_numpy_pdm_gnum(_dentity2_entity1, np_dentity2_entity1_idx[dn_entity2])

    return np_dentity2_entity1_idx, np_dentity2_entity1

# ------------------------------------------------------------------------
def combine_connectivity(NPY.ndarray[int, mode='c', ndim=1]    entity1_entity2_idx,
                         NPY.ndarray[int, mode='c', ndim=1]    entity1_entity2,
                         NPY.ndarray[int, mode='c', ndim=1]    entity2_entity3_idx,
                         NPY.ndarray[int, mode='c', ndim=1]    entity2_entity3):

    assert_single_dim_np(entity1_entity2, NPY.int32, entity1_entity2_idx[-1])
    assert_single_dim_np(entity2_entity3, NPY.int32, entity2_entity3_idx[-1])
    
    cdef int n_entity1 = entity1_entity2_idx.size -1
    
    cdef int *_entity1_entity2_idx
    assert entity1_entity2_idx is not None
    _entity1_entity2_idx = <int *> entity1_entity2_idx.data
    
    cdef int *_entity1_entity2
    _entity1_entity2 = <int *> entity1_entity2.data
    
    cdef int *_entity2_entity3_idx
    assert entity1_entity2_idx is not None
    _entity2_entity3_idx = <int *> entity2_entity3_idx.data
    
    cdef int *_entity2_entity3
    _entity2_entity3 = <int *> entity2_entity3.data
    
    cdef int *_entity1_entity3_idx = NULL
    cdef int *_entity1_entity3     = NULL
    
    PDM_combine_connectivity( n_entity1,
                              _entity1_entity2_idx,
                              _entity1_entity2,
                              _entity2_entity3_idx,
                              _entity2_entity3,
                             &_entity1_entity3_idx,
                             &_entity1_entity3)
    
    assert _entity1_entity3_idx != NULL
    
    np_entity1_entity3_idx = create_numpy_i(_entity1_entity3_idx, n_entity1 + 1)
    np_entity1_entity3     = create_numpy_i(_entity1_entity3, np_entity1_entity3_idx[n_entity1])
    
    return np_entity1_entity3_idx, np_entity1_entity3

# ------------------------------------------------------------------------
def connectivity_transpose(NPY.int                               n_entity2, # We have to pass n_entity2 to manage empty tabs and gap numerbering
                           NPY.ndarray[int, mode='c', ndim=1]    entity1_entity2_idx,
                           NPY.ndarray[int, mode='c', ndim=1]    entity1_entity2):
    
    assert_single_dim_np(entity1_entity2, NPY.int32, entity1_entity2_idx[-1])
    
    cdef int n_entity1 = entity1_entity2_idx.size -1
    
    cdef int *_entity1_entity2_idx
    assert entity1_entity2_idx is not None
    _entity1_entity2_idx = <int *> entity1_entity2_idx.data
    
    cdef int *_entity1_entity2
    _entity1_entity2 = <int *> entity1_entity2.data
    
    cdef int *_entity2_entity1_idx = NULL
    cdef int *_entity2_entity1     = NULL
    
    PDM_connectivity_transpose(      n_entity1,
                               <int> n_entity2,
                                     _entity1_entity2_idx,
                                     _entity1_entity2,
                                    &_entity2_entity1_idx,
                                    &_entity2_entity1)
    
    assert _entity2_entity1_idx != NULL
    
    np_entity2_entity1_idx = create_numpy_i(_entity2_entity1_idx, n_entity2 + 1)
    np_entity2_entity1     = create_numpy_i(_entity2_entity1, np_entity2_entity1_idx[n_entity2])
    
    return np_entity2_entity1_idx, np_entity2_entity1

# ------------------------------------------------------------------------
def part_connectivity_transpose(list   n_entity2, # We have to pass n_entity2 to manage empty tabs and gap numerbering for each partition
                                list   l_entity1_entity2_idx,
                                list   l_entity1_entity2):

    assert(len(n_entity2) == len(l_entity1_entity2_idx) == len(l_entity1_entity2))
    
    cdef int n_part = len(n_entity2)
    
    n_entity1 = [entity1_entity2_idx.size -1 for entity1_entity2_idx in l_entity1_entity2_idx]
    cdef int* _n_entity1 = list_to_int_pointer(n_entity1)
    
    cdef int* _n_entity2 = list_to_int_pointer(n_entity2)
    
    cdef int** _entity1_entity2_idx = np_list_to_int_pointers(l_entity1_entity2_idx)
    cdef int** _entity1_entity2     = np_list_to_int_pointers(l_entity1_entity2)

    _entity2_entity1_idx = <int **> malloc(n_part * sizeof(int *))
    _entity2_entity1     = <int **> malloc(n_part * sizeof(int *))
    
    PDM_part_connectivity_transpose( n_part,
                                     _n_entity1,
                                     _n_entity2,
                                     _entity1_entity2_idx,
                                     _entity1_entity2,
                                    &_entity2_entity1_idx,
                                    &_entity2_entity1)
    
    l_np_entity2_entity1_idx = list()
    l_np_entity2_entity1     = list()
    
    for i_part in range(n_part):
        
        assert _entity2_entity1_idx[i_part] != NULL
        
        np_entity2_entity1_idx = create_numpy_i(_entity2_entity1_idx[i_part], _n_entity2[i_part] + 1)
        np_entity2_entity1     = create_numpy_i(_entity2_entity1[i_part], np_entity2_entity1_idx[_n_entity2[i_part]])
        
        l_np_entity2_entity1_idx.append(np_entity2_entity1_idx)
        l_np_entity2_entity1.append(np_entity2_entity1)
    
    free(_n_entity1)
    free(_n_entity2)
    free(_entity1_entity2_idx)
    free(_entity1_entity2)
    free(_entity2_entity1_idx)
    free(_entity2_entity1)
    
    return l_np_entity2_entity1_idx, l_np_entity2_entity1

# ------------------------------------------------------------------------
def part_connectivity_to_connectity_idx(list   n_entity1,
                                        list   l_entity1_entity2_in):

    assert(len(n_entity1) == len(l_entity1_entity2_in))
    
    cdef int n_part = len(n_entity1)

    cdef int* _n_entity1 = list_to_int_pointer(n_entity1)
    
    cdef int** _entity1_entity2_in = np_list_to_int_pointers(l_entity1_entity2_in)

    _entity1_entity2_idx = <int **> malloc(n_part * sizeof(int *))
    _entity1_entity2     = <int **> malloc(n_part * sizeof(int *))
    
    
    PDM_part_connectivity_to_connectity_idx( n_part,
                                             _n_entity1,
                                             _entity1_entity2_in,
                                            &_entity1_entity2_idx,
                                            &_entity1_entity2)
    
    l_np_entity1_entity2_idx = list()
    l_np_entity1_entity2     = list()
    
    for i_part in range(n_part):
        
        assert _entity1_entity2_idx[i_part] != NULL
        
        np_entity1_entity2_idx = create_numpy_i(_entity1_entity2_idx[i_part], _n_entity1[i_part] + 1)
        np_entity1_entity2     = create_numpy_i(_entity1_entity2[i_part], np_entity1_entity2_idx[_n_entity1[i_part]])
        
        l_np_entity1_entity2_idx.append(np_entity1_entity2_idx)
        l_np_entity1_entity2.append(np_entity1_entity2)

    free(_n_entity1)
    free(_entity1_entity2_in)
    free(_entity1_entity2_idx)
    free(_entity1_entity2)
    
    return l_np_entity1_entity2_idx, l_np_entity1_entity2
