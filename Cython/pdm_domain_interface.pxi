
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

  int PDM_domain_interface_get_as_graph(PDM_domain_interface_t *dom_intrf,
                                        PDM_bound_type_t        interface_kind,
                                        int                   **interface_graph_idx,
                                        PDM_g_num_t           **interface_graph_ids,
                                        int                   **interface_graph_dom)

  void PDM_domain_interface_translate_face2vtx(PDM_domain_interface_t  *dom_intrf,
                                               int                     *dn_vtx,
                                               int                     *dn_face,
                                               int                    **dface_vtx_idx,
                                               PDM_g_num_t            **dface_vtx)

  void PDM_domain_interface_free(PDM_domain_interface_t *dom_intrf)

# ===================================================================================

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

    # Interfaces data
    cdef int           *_interface_dn_face = list_to_int_pointer(interface_dn_face)
    cdef PDM_g_num_t **_interface_ids_face = np_list_to_gnum_pointers(interface_ids_face)
    cdef int         **_interface_dom_face
    if multizone_interface:
      _interface_dom_face = np_list_to_int_pointers(interface_dom_face)
    else:
      _interface_dom_face = <int **> malloc(n_interface*sizeof(int*));
      for i in range(n_interface):
        _interface_dom_face[i] = <int *> malloc(2*sizeof(int))
        _interface_dom_face[i][0] = interface_dom_face[i][0]
        _interface_dom_face[i][1] = interface_dom_face[i][1]

    # Zone data
    cdef int          *_dn_vtx        = list_to_int_pointer(dn_vtx)
    cdef int          *_dn_face       = list_to_int_pointer(dn_face)
    cdef int         **_dface_vtx_idx = np_list_to_int_pointers(dface_vtx_idx)
    cdef PDM_g_num_t **_dface_vtx     = np_list_to_gnum_pointers(dface_vtx)

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

def interface_to_graph(int       n_interface,
                       bint      multizone_interface,
                       list      interface_dn,
                       list      interface_ids,
                       list      interface_dom,
                       MPI.Comm  comm):


    #Some checks
    assert (len(interface_dn) == len(interface_ids) == len(interface_dom) == n_interface)
    for i in range(n_interface):
      assert_single_dim_np(interface_ids[i], npy_pdm_gnum_dtype, 2*interface_dn[i])
      if multizone_interface:
        assert_single_dim_np(interface_dom[i], NPY.int32, 2*interface_dn[i])

    #Convert input data
    cdef MPI.MPI_Comm c_comm = comm.ob_mpi
    cdef PDM_MPI_Comm PDMC   = PDM_MPI_mpi_2_pdm_mpi_comm(&c_comm)

    cdef int           *_interface_dn = list_to_int_pointer(interface_dn)
    cdef PDM_g_num_t **_interface_ids = np_list_to_gnum_pointers(interface_ids)
    cdef int         **_interface_dom
    if multizone_interface:
      _interface_dom = np_list_to_int_pointers(interface_dom)
    else:
      _interface_dom = <int **> malloc(n_interface*sizeof(int*));
      for i in range(n_interface):
        _interface_dom[i] = <int *> malloc(2*sizeof(int))
        _interface_dom[i][0] = interface_dom[i][0]
        _interface_dom[i][1] = interface_dom[i][1]

    # Run function
    cdef PDM_domain_interface_t *dom_intrf;
    cdef _multizone_interface = PDM_DOMAIN_INTERFACE_MULT_YES if multizone_interface else PDM_DOMAIN_INTERFACE_MULT_NO
    dom_intrf = PDM_domain_interface_create(n_interface,
                                            -1, #Not used
                                            _multizone_interface,
                                            PDM_OWNERSHIP_USER,
                                            PDMC)
    PDM_domain_interface_set(dom_intrf,
                             PDM_BOUND_TYPE_VTX, #Type does not matter
                             _interface_dn,
                             _interface_ids,
                             _interface_dom)

    cdef int*         _interface_graph_idx
    cdef PDM_g_num_t* _interface_graph_ids
    cdef int*         _interface_graph_dom
    
    cdef int graph_dn = PDM_domain_interface_get_as_graph(dom_intrf,
                                                          PDM_BOUND_TYPE_VTX,
                                                         &_interface_graph_idx,
                                                         &_interface_graph_ids,
                                                         &_interface_graph_dom)


    interface_graph_idx = create_numpy_i(_interface_graph_idx, graph_dn+1)
    PyArray_ENABLEFLAGS(interface_graph_idx, NPY.NPY_OWNDATA)

    interface_graph_ids = create_numpy_pdm_gnum(_interface_graph_ids, interface_graph_idx[graph_dn])
    if interface_graph_ids is not None:
      PyArray_ENABLEFLAGS(interface_graph_ids, NPY.NPY_OWNDATA)
    interface_graph_dom = create_numpy_i(_interface_graph_dom, interface_graph_idx[graph_dn])
    if interface_graph_dom is not None:
      PyArray_ENABLEFLAGS(interface_graph_dom, NPY.NPY_OWNDATA)

    # Free temporary objects and return
    PDM_domain_interface_free(dom_intrf)
    free(_interface_dn )
    free(_interface_ids)
    if not multizone_interface:
      for i in range(n_interface):
        free(_interface_dom[i])
    free(_interface_dom)

    return interface_graph_idx, interface_graph_ids, interface_graph_dom
