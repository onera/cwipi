
cdef extern from "pdm_domain_interface.h":
  ctypedef struct PDM_domain_interface_t:
      pass
  ctypedef enum PDM_domain_interface_mult_t:
      PDM_DOMAIN_INTERFACE_MULT_NO  = 0
      PDM_DOMAIN_INTERFACE_MULT_YES = 1

  PDM_domain_interface_t* PDM_domain_interface_create(const int                   n_interface,
                                                      const int                   n_zone,
                                                      PDM_domain_interface_mult_t multizone_interface,
                                                      PDM_ownership_t             ownership,
                                                      PDM_MPI_Comm                comm)

  void PDM_domain_interface_set(PDM_domain_interface_t *dom_intrf,
                                PDM_bound_type_t        interface_kind,
                                int                    *interface_dn,
                                PDM_g_num_t           **interface_ids,
                                int                   **interface_dom)

  void PDM_domain_interface_get(PDM_domain_interface_t *dom_intrf,
                                PDM_bound_type_t        interface_kind,
                                int                   **interface_dn,
                                PDM_g_num_t          ***interface_ids,
                                int                  ***interface_dom)

  void PDM_domain_interface_translate_face2vtx(PDM_domain_interface_t  *dom_intrf,
                                               int                     *dn_vtx,
                                               int                     *dn_face,
                                               int                    **dface_vtx_idx,
                                               PDM_g_num_t            **dface_vtx)

  void PDM_domain_interface_free(PDM_domain_interface_t *dom_intrf)

# ===================================================================================
def assert_single_dim_np(tab, dtype, size=None):
  assert isinstance(tab, NPY.ndarray) and tab.ndim == 1 and tab.dtype == dtype
  if size is not None:
    assert tab.size == size

def interface_face_to_vertex(int       n_interface,
                             int       n_zone,
                             bint      multizone_interface,
                             list      interface_dn_face,
                             list      interface_ids_face,
                             list      interface_dom_face,
                             list      dn_vtx,
                             list      dn_face,
                             list      dface_vtx_idx,
                             list      dface_vtx,
                             MPI.Comm  comm):
    # Generic declarations
    cdef NPY.ndarray[NPY.int32_t,    ndim=1, mode='c'] numpy_int
    cdef NPY.ndarray[npy_pdm_gnum_t, ndim=1, mode='c'] numpy_gnum

    #Some checks
    assert (len(interface_dn_face) == len(interface_ids_face) == n_interface)
    assert (len(dn_vtx) == len(dn_face) == len(dface_vtx_idx) == len(dface_vtx) == n_zone)
    for i in range(n_zone):
      assert_single_dim_np(dface_vtx_idx[i], NPY.int32)
      assert_single_dim_np(dface_vtx[i], npy_pdm_gnum_dtype)
    for i in range(n_interface):
      assert_single_dim_np(interface_ids_face[i], npy_pdm_gnum_dtype, 2*interface_dn_face[i])
      if multizone_interface:
        assert_single_dim_np(interface_dom_face[i], NPY.int32,          2*interface_dn_face[i])

    #Convert input data
    cdef MPI.MPI_Comm c_comm = comm.ob_mpi
    cdef PDM_MPI_Comm PDMC   = PDM_MPI_mpi_2_pdm_mpi_comm(&c_comm)

    cdef int           *_interface_dn_face = <int          *> malloc(n_interface*sizeof(int));
    cdef PDM_g_num_t **_interface_ids_face = <PDM_g_num_t **> malloc(n_interface*sizeof(PDM_g_num_t*));
    cdef int         **_interface_dom_face = <int         **> malloc(n_interface*sizeof(int*));

    for i in range(n_interface):
      _interface_dn_face[i]  = <int> interface_dn_face[i]
      numpy_gnum = interface_ids_face[i]
      _interface_ids_face[i] = <PDM_g_num_t *> numpy_gnum.data
      if multizone_interface:
        numpy_int = interface_dom_face[i]
        _interface_dom_face[i] = <int *> numpy_int.data
      else:
        _interface_dom_face[i] = <int *> malloc(2*sizeof(int))
        _interface_dom_face[i][0] = interface_dom_face[i][0]
        _interface_dom_face[i][1] = interface_dom_face[i][1]

    cdef int          *_dn_vtx        = <int         *>  malloc(n_zone * sizeof(int));
    cdef int          *_dn_face       = <int         *>  malloc(n_zone * sizeof(int));
    cdef int         **_dface_vtx_idx = <int         **> malloc(n_zone * sizeof(int*));
    cdef PDM_g_num_t **_dface_vtx     = <PDM_g_num_t **> malloc(n_zone * sizeof(PDM_g_num_t*));
    for i in range(n_zone):
      _dn_vtx[i] = <int> dn_vtx[i]
      _dn_face[i] = <int> dn_face[i]
      numpy_int = dface_vtx_idx[i]
      _dface_vtx_idx[i] =  <int *> numpy_int.data
      numpy_gnum = dface_vtx[i]
      _dface_vtx[i] =  <PDM_g_num_t *> numpy_gnum.data

    # Run function
    cdef PDM_domain_interface_t *dom_intrf;
    cdef _multizone_interface = PDM_DOMAIN_INTERFACE_MULT_YES if multizone_interface else PDM_DOMAIN_INTERFACE_MULT_NO
    dom_intrf = PDM_domain_interface_create(n_interface,
                                            n_zone,
                                            _multizone_interface,
                                            PDM_OWNERSHIP_USER,
                                            PDMC)
    PDM_domain_interface_set(dom_intrf,
                             PDM_BOUND_TYPE_FACE,
                             _interface_dn_face,
                             _interface_ids_face,
                             _interface_dom_face)

    PDM_domain_interface_translate_face2vtx(dom_intrf,
                                            _dn_vtx,
                                            _dn_face,
                                            _dface_vtx_idx,
                                            _dface_vtx)

    
    # Convert output data
    cdef int          *_interface_dn_vtx  = NULL;
    cdef PDM_g_num_t **_interface_ids_vtx = NULL;
    cdef int         **_interface_dom_vtx = NULL;
    PDM_domain_interface_get(dom_intrf,
                             PDM_BOUND_TYPE_VTX,
                            &_interface_dn_vtx,
                            &_interface_ids_vtx,
                            &_interface_dom_vtx)

    vtx_interface = list()
    for i in range(n_interface):
      interface_ids_vtx = create_numpy_pdm_gnum(_interface_ids_vtx[i], 2*_interface_dn_vtx[i])
      if interface_ids_vtx is not None:
        PyArray_ENABLEFLAGS(interface_ids_vtx, NPY.NPY_OWNDATA)
      interface_results = {'interface_dn_vtx' : _interface_dn_vtx[i], 'np_interface_ids_vtx' : interface_ids_vtx}
      if multizone_interface: #Return domains only if we had complex interfaces
        interface_dom_vtx = create_numpy_i(_interface_dom_vtx[i], 2*_interface_dn_vtx[i])
        if interface_dom_vtx is not None:
          PyArray_ENABLEFLAGS(interface_dom_vtx, NPY.NPY_OWNDATA)
        interface_results['np_interface_dom_vtx'] = interface_dom_vtx
      vtx_interface.append(interface_results)

    # Free temporary objects and return
    PDM_domain_interface_free(dom_intrf)
    free(_interface_dn_face )
    free(_interface_ids_face)
    if not multizone_interface:
      for i in range(n_interface):
        free(_interface_dom_face[i])

    free(_interface_dom_face)
    free(_dn_vtx       )
    free(_dn_face      )
    free(_dface_vtx_idx)
    free(_dface_vtx    )
    free(_interface_dn_vtx )
    free(_interface_ids_vtx)
    if multizone_interface: #Same than face dom if not multizone_interface
      free(_interface_dom_vtx)

    return vtx_interface

