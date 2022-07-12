
cdef extern from "pdm_partitioning_algorithm.h":
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    void PDM_part_distgroup_to_partgroup(PDM_MPI_Comm      comm,
                                         PDM_g_num_t      *entity_distribution,
                                         int               n_group,
                                         int              *dgroup_idx,
                                         PDM_g_num_t      *dgroup,
                                         int               n_part,
                                         int              *pn_entity,
                                         PDM_g_num_t     **pentity_ln_to_gn,
                                         int            ***pgroup_idx,
                                         int            ***pgroup,
                                         PDM_g_num_t    ***pgroup_ln_to_gn)
    void PDM_part_dcoordinates_to_pcoordinates(PDM_MPI_Comm    comm,
                                               int             n_part,
                                               PDM_g_num_t    *vertex_distribution,
                                               double         *dvtx_coord,
                                               int            *pn_vtx,
                                               PDM_g_num_t   **pvtx_ln_to_gn,
                                               double       ***pvtx_coord)
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


# ===================================================================================
def part_distgroup_to_partgroup(MPI.Comm                                      comm,
                                NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] entity_distribution,
                                int                                           n_group,
                                NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] dgroup_idx,
                                NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] dgroup,
                                int                                           n_part,
                                list                                          pn_entity,
                                list                                          pentity_ln_to_gn):

    # ::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Convert mpi4py -> PDM_MPI
    cdef MPI.MPI_Comm c_comm = comm.ob_mpi
    cdef PDM_MPI_Comm PDMC   = PDM_MPI_mpi_2_pdm_mpi_comm(&c_comm)

    # ~> \param [in]   entity_distribution   Distributed entity
    cdef PDM_g_num_t * entity_distribution_data
    if (entity_distribution is None):
        entity_distribution_data = NULL
    else:
        entity_distribution_data = <PDM_g_num_t *> entity_distribution.data

    # ~> \param [in]   n_goupe
    cdef int _n_group = n_group

    # ~> \param [in]   dgroup_idx
    cdef int *dgroup_idx_data
    if (dgroup_idx is None):
        dgroup_idx_data = NULL
    else:
        dgroup_idx_data = <int *> dgroup_idx.data

    # ~> \param [in]   distributed group
    cdef PDM_g_num_t * dgroup_data
    if (dgroup is None):
        dgroup_data = NULL
    else:
        dgroup_data = <PDM_g_num_t *> dgroup.data

    # ~> \param [in]   Npart
    cdef int _n_part = n_part

    # ~> \param [in]   pn_entity
    cdef int * _pn_entity
    _pn_entity        = <int *> malloc(sizeof(int) * _n_part )
    for idx, part_pn_entity  in enumerate(pn_entity):
        _pn_entity[idx] = part_pn_entity

    # ~> \param [in]   pentity_ln_to_gn local to global
    cdef NPY.ndarray[npy_pdm_gnum_t, ndim=1, mode='c'] partLNToGN
    cdef PDM_g_num_t ** LNToGN

    LNToGN        = <PDM_g_num_t **> malloc(sizeof(PDM_g_num_t *) * _n_part )
    for idx, partLNToGN in enumerate(pentity_ln_to_gn):
        LNToGN[idx] = <PDM_g_num_t *> partLNToGN.data

    # ~> \param [out]   pgroup_idx
    cdef int **_pgroup_idx

    # ~> \param [out]   pgroup
    cdef int **_pgroup

    # ~> \param [out]   pgroup_ln_to_gn
    cdef PDM_g_num_t **_pgroup_ln_to_gn

    PDM_part_distgroup_to_partgroup(
                                    PDMC,
                                    entity_distribution_data,
                                    _n_group,
                                    dgroup_idx_data,
                                    dgroup_data,
                                    _n_part,
                                    _pn_entity,
             <const PDM_g_num_t **> LNToGN,
                                    &_pgroup_idx,
                                    &_pgroup,
                                    &_pgroup_ln_to_gn
                                    )
    list_group_part = list()

    for i_part in range(_n_part):

        dict_group_part    = dict()
        # ~>  \param [out]  numpy array of pgroup_idx
        if (_pgroup_idx == NULL) :
            npgroup_idx = None
        else :
            dim = <NPY.npy_intp> (_n_group + 1)
            npgroup_idx = NPY.PyArray_SimpleNewFromData(1,
                                                        &dim,
                                                        NPY.NPY_INT32,
                                                        <void *> _pgroup_idx[i_part])
            PyArray_ENABLEFLAGS(npgroup_idx, NPY.NPY_OWNDATA);
        # ~>  \param [out]  numpy array of pgroup
        if (_pgroup == NULL) :
            npgroup = None
        else :
            dim = _pgroup_idx[i_part][_n_group]
            npgroup = NPY.PyArray_SimpleNewFromData(1,
                                                    &dim,
                                                    NPY.NPY_INT32,
                                                    <void *> _pgroup[i_part])
            PyArray_ENABLEFLAGS(npgroup, NPY.NPY_OWNDATA)

        # ~>  \param [out]  numpy array of pgroup
        if (_pgroup_ln_to_gn == NULL) :
            npgroup_ln_to_gn = None
        else :
            dim = _pgroup_idx[i_part][_n_group]
            npgroup_ln_to_gn = NPY.PyArray_SimpleNewFromData(1,
                                                    &dim,
                                                    PDM_G_NUM_NPY_INT,
                                                    <void *> _pgroup_ln_to_gn[i_part])
            PyArray_ENABLEFLAGS(npgroup_ln_to_gn, NPY.NPY_OWNDATA)

        dict_group_part['npZSRGroupIdx']    = npgroup_idx
        dict_group_part['npZSRGroup']       = npgroup
        dict_group_part['npZSRGroupLNToGN'] = npgroup_ln_to_gn

        list_group_part.append(dict_group_part)

    free(_pn_entity)
    free(_pgroup_idx)
    free(_pgroup)
    free(_pgroup_ln_to_gn)
    free(LNToGN)

    return list_group_part

cdef extern from "pdm_closest_points.h":
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    void PDM_transform_to_parent_gnum(int           n_part_initial,
                                      int          *n_elmt_initial,
                                      PDM_g_num_t **child_ln_to_gn,
                                      PDM_g_num_t **parent_ln_to_gn,
                                      int           n_part_to_transform,
                                      int          *n_elmt_to_transform,
                                      PDM_g_num_t **gnum_to_transform,
                                      PDM_MPI_Comm  comm)
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# ===================================================================================
#Why here ? Could be moved in closest_point
def transform_to_parent_gnum(list     gnum_to_transform,
                             list     child_ln_to_gn,
                             list     parent_ln_to_gn,
                             MPI.Comm comm):
  """
  """
  assert len(parent_ln_to_gn) == len(child_ln_to_gn)
  for child_gn, parent_gn in zip(child_ln_to_gn, parent_ln_to_gn):
    assert child_gn.shape == parent_gn.shape

  cdef NPY.ndarray[npy_pdm_gnum_t, ndim=1, mode='c'] np_array
  cdef int n_part_out = len(gnum_to_transform)
  cdef int n_part_ini   = len(child_ln_to_gn)

  cdef PDM_g_num_t **_gnum_to_transform = <PDM_g_num_t **> malloc(sizeof(PDM_g_num_t *) * n_part_out)
  cdef int          *_pn_elt_out        = <int *>          malloc(sizeof(int)           * n_part_out)
  for idx, np_array in enumerate(gnum_to_transform):
      _gnum_to_transform[idx] = <PDM_g_num_t *> np_array.data
      _pn_elt_out[idx] = <int> np_array.shape[0]

  cdef PDM_g_num_t **_child_ln_to_gn  = <PDM_g_num_t **> malloc(sizeof(PDM_g_num_t *) * n_part_ini)
  cdef PDM_g_num_t **_parent_ln_to_gn = <PDM_g_num_t **> malloc(sizeof(PDM_g_num_t *) * n_part_ini)
  cdef int          *_pn_elt_in       = <int *>          malloc(sizeof(int)           * n_part_ini)
  for idx in range(n_part_ini):
      np_array = child_ln_to_gn[idx]
      _child_ln_to_gn[idx] = <PDM_g_num_t *> np_array.data
      np_array = parent_ln_to_gn[idx]
      _parent_ln_to_gn[idx] = <PDM_g_num_t *> np_array.data
      _pn_elt_in[idx] = <int> np_array.shape[0]
    
  cdef MPI.MPI_Comm c_comm = comm.ob_mpi
  cdef PDM_MPI_Comm PDMC   = PDM_MPI_mpi_2_pdm_mpi_comm(&c_comm)
  PDM_transform_to_parent_gnum(                 n_part_ini,
                               <const int          *> _pn_elt_in,
                               <const PDM_g_num_t **> _child_ln_to_gn,
                               <const PDM_g_num_t **> _parent_ln_to_gn,
                                                n_part_out,
                               <const int          *> _pn_elt_out,
                               <PDM_g_num_t **> _gnum_to_transform,
                                                PDMC)

  free(_pn_elt_out)
  free(_pn_elt_in)
  free(_child_ln_to_gn)
  free(_parent_ln_to_gn)
  free(_gnum_to_transform)
# ===================================================================================


# ===================================================================================
def part_dcoordinates_to_pcoordinates(MPI.Comm                                      comm,
                                      NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] vtx_distribution,
                                      NPY.ndarray[NPY.double_t  , mode='c', ndim=1] dvtx_coord   not None,
                                      list                                          l_pvtx_ln_to_gn):
    """
    """
    cdef MPI.MPI_Comm c_comm = comm.ob_mpi
    cdef PDM_MPI_Comm PDMC   = PDM_MPI_mpi_2_pdm_mpi_comm(&c_comm)
    cdef NPY.ndarray[npy_pdm_gnum_t, ndim=1, mode='fortran'] part_ln_to_gn

    cdef int n_part = len(l_pvtx_ln_to_gn)
    cdef int          *pn_vtx        = <int          *> malloc(n_part * sizeof(int          ))
    cdef PDM_g_num_t **pvtx_ln_to_gn = <PDM_g_num_t **> malloc(n_part * sizeof(PDM_g_num_t *))
    cdef double      **pvtx_coord

    for i, part_ln_to_gn in enumerate(l_pvtx_ln_to_gn):
        pn_vtx[i]        = l_pvtx_ln_to_gn[i].shape[0]
        pvtx_ln_to_gn[i] = <PDM_g_num_t *> part_ln_to_gn.data

    PDM_part_dcoordinates_to_pcoordinates(PDMC,
                                          n_part,
                     <const PDM_g_num_t*> vtx_distribution.data,
                     <const double*>      dvtx_coord.data,
                     <const int *>        pn_vtx,
                     <const PDM_g_num_t**>pvtx_ln_to_gn,
                                          &pvtx_coord);

    l_pvtx_coord = list()
    for i in range(n_part):
        dim = <NPY.npy_intp> 3 * pn_vtx[i]
        np_vtx_coord = NPY.PyArray_SimpleNewFromData(1,
                                                   &dim,
                                                   NPY.NPY_DOUBLE,
                                          <void *> pvtx_coord[i])
        PyArray_ENABLEFLAGS(np_vtx_coord, NPY.NPY_OWNDATA);
        l_pvtx_coord.append(np_vtx_coord)

    free(pvtx_coord)
    free(pn_vtx)
    free(pvtx_ln_to_gn)

    return l_pvtx_coord
