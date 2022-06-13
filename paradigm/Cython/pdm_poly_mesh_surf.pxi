
cdef extern from "pdm_poly_surf_gen.h":
    # ------------------------------------------------------------------
    void PDM_poly_surf_gen(PDM_MPI_Comm  localComm,
                           double        xmin,
                           double        xmax,
                           double        ymin,
                           double        ymax,
                           int           have_random,
                           int           init_random,
                           PDM_g_num_t   nx,
                           PDM_g_num_t   ny,
                           PDM_g_num_t  *n_face,
                           PDM_g_num_t  *n_vtx,
                           PDM_g_num_t  *n_edge,
                           int          *dn_vtx,
                           double      **dvtx_coord,
                           int          *dn_face,
                           int         **dface_vtx_idx,
                           PDM_g_num_t **dface_vtx,
                           PDM_g_num_t **dface_edge,
                           int          *dn_edge,
                           PDM_g_num_t **dedge_vtx,
                           PDM_g_num_t **dedge_face,
                           int          *n_edge_group,
                           int         **dedge_group_idx,
                           PDM_g_num_t **dedge_group)
    # ------------------------------------------------------------------

# ------------------------------------------------------------------------
def PolyMeshSurf(double      xmin,
                 double      xmax,
                 double      ymin,
                 double      ymax,
                 int         have_random,
                 int         init_random,
                 PDM_g_num_t nx,
                 PDM_g_num_t ny,
                 MPI.Comm    comm):
  """
  """
  # ************************************************************************
  # > Declaration
  cdef PDM_g_num_t  n_face
  cdef PDM_g_num_t  n_vtx
  cdef PDM_g_num_t  n_edge
  cdef int          dn_vtx
  cdef double      *dvtx_coord
  cdef int          dn_face
  cdef int         *dface_vtx_idx
  cdef PDM_g_num_t *dface_vtx
  cdef PDM_g_num_t *dface_edge
  cdef int          dn_edge
  cdef PDM_g_num_t *dedge_vtx
  cdef PDM_g_num_t *dedge_face
  cdef int          n_edge_group
  cdef int         *edge_group
  cdef int         *dedge_group_idx
  cdef PDM_g_num_t *dedge_group
  # ************************************************************************

  cdef MPI.MPI_Comm c_comm = comm.ob_mpi
  cdef PDM_MPI_Comm PDMC   = PDM_MPI_mpi_2_pdm_mpi_comm(<void *> &c_comm)

  PDM_poly_surf_gen(PDMC,
                    xmin,
                    xmax,
                    ymin,
                    ymax,
                    have_random,
                    init_random,
                    nx,
                    ny,
                    &n_face,
                    &n_vtx,
                    &n_edge,
                    &dn_vtx,
                    &dvtx_coord,
                    &dn_face,
                    &dface_vtx_idx,
                    &dface_vtx,
                    &dface_edge,
                    &dn_edge,
                    &dedge_vtx,
                    &dedge_face,
                    &n_edge_group,
                    &dedge_group_idx,
                    &dedge_group)

  # -> Begin
  cdef NPY.npy_intp dim

  # ::::::::::::::::::::::::::::::::::::::::::::::::::
  dim = <NPY.npy_intp> (3 * dn_vtx)
  np_dvtx_coord  = NPY.PyArray_SimpleNewFromData(1,
                                                 &dim,
                                                 NPY.NPY_DOUBLE,
                                                 <void *> dvtx_coord)
  PyArray_ENABLEFLAGS(np_dvtx_coord, NPY.NPY_OWNDATA);
  # ::::::::::::::::::::::::::::::::::::::::::::::::::

  # ::::::::::::::::::::::::::::::::::::::::::::::::::
  dim = <NPY.npy_intp> dn_face + 1
  np_dface_vtx_idx = NPY.PyArray_SimpleNewFromData(1,
                                                   &dim,
                                                   NPY.NPY_INT32,
                                                   <void *> dface_vtx_idx)
  PyArray_ENABLEFLAGS(np_dface_vtx_idx, NPY.NPY_OWNDATA);
  # ::::::::::::::::::::::::::::::::::::::::::::::::::

  # ::::::::::::::::::::::::::::::::::::::::::::::::::
  dim = <NPY.npy_intp> np_dface_vtx_idx[np_dface_vtx_idx.shape[0]-1]
  np_dface_vtx = NPY.PyArray_SimpleNewFromData(1,
                                              &dim,
                                              PDM_G_NUM_NPY_INT,
                                              <void *> dface_vtx)
  PyArray_ENABLEFLAGS(np_dface_vtx, NPY.NPY_OWNDATA);
  # ::::::::::::::::::::::::::::::::::::::::::::::::::

  # ::::::::::::::::::::::::::::::::::::::::::::::::::
  # > In 2D cas number of vertex is the same as number of edge
  dim = <NPY.npy_intp> np_dface_vtx_idx[np_dface_vtx_idx.shape[0]-1]
  np_dface_edge = NPY.PyArray_SimpleNewFromData(1,
                                                &dim,
                                                PDM_G_NUM_NPY_INT,
                                                <void *> dface_edge)
  PyArray_ENABLEFLAGS(np_dface_edge, NPY.NPY_OWNDATA);
  # ::::::::::::::::::::::::::::::::::::::::::::::::::

  # ::::::::::::::::::::::::::::::::::::::::::::::::::
  dim = <NPY.npy_intp> 2 * dn_edge
  np_dedge_vtx= NPY.PyArray_SimpleNewFromData(1,
                                              &dim,
                                              PDM_G_NUM_NPY_INT,
                                              <void *> dedge_vtx)
  PyArray_ENABLEFLAGS(np_dedge_vtx, NPY.NPY_OWNDATA);
  # ::::::::::::::::::::::::::::::::::::::::::::::::::

  # ::::::::::::::::::::::::::::::::::::::::::::::::::
  dim = <NPY.npy_intp> 2 * dn_edge
  np_dedge_face= NPY.PyArray_SimpleNewFromData(1,
                                               &dim,
                                               PDM_G_NUM_NPY_INT,
                                               <void *> dedge_face)
  PyArray_ENABLEFLAGS(np_dedge_face, NPY.NPY_OWNDATA);
  # ::::::::::::::::::::::::::::::::::::::::::::::::::

  # ::::::::::::::::::::::::::::::::::::::::::::::::::
  dim = <NPY.npy_intp> n_edge_group + 1
  np_dedge_group_idx= NPY.PyArray_SimpleNewFromData(1,
                                                    &dim,
                                                    NPY.NPY_INT32,
                                                    <void *> dedge_group_idx)
  PyArray_ENABLEFLAGS(np_dedge_group_idx, NPY.NPY_OWNDATA);
  # ::::::::::::::::::::::::::::::::::::::::::::::::::

  # ::::::::::::::::::::::::::::::::::::::::::::::::::
  dim = <NPY.npy_intp> np_dedge_group_idx[np_dedge_group_idx.shape[0]-1]
  np_dedge_group= NPY.PyArray_SimpleNewFromData(1,
                                                &dim,
                                                PDM_G_NUM_NPY_INT,
                                                <void *> dedge_group)
  PyArray_ENABLEFLAGS(np_dedge_group, NPY.NPY_OWNDATA);
  # ::::::::::::::::::::::::::::::::::::::::::::::::::

  return {'n_face'          : n_face,
          'n_vtx'           : n_vtx,
          'n_edge'          : n_edge,
          'dn_face'         : dn_face,
          'dn_vtx'          : dn_vtx,
          'dn_edge'         : dn_edge,
          'n_edge_group'    : n_edge_group,
          'dvtx_coord'      : np_dvtx_coord,
          'dface_vtx_idx'   : np_dface_vtx_idx,
          'dface_vtx'       : np_dface_vtx,
          'dface_edge'      : np_dface_edge,
          'dedge_vtx'       : np_dedge_vtx,
          'dedge_face'      : np_dedge_face,
          'dedge_group_idx' : np_dedge_group_idx,
          'dedge_group'     : np_dedge_group,
          }
