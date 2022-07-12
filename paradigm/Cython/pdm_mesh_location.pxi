
cdef extern from "pdm_mesh_location.h":
  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  # > Wrapping of Ppart Structure
  ctypedef struct PDM_mesh_location_t:
    pass
  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  # > Wrapping for C
  ctypedef enum PDM_interpolate_kind_t:
    PDM_INTERPOLATE_KIND_FROM_CENTER  = 0
    PDM_INTERPOLATE_KIND_FROM_NODE = 1
  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  PDM_mesh_location_t* PDM_mesh_location_create(PDM_mesh_nature_t mesh_nature,
                               int               _n_point_cloud,
                               PDM_MPI_Comm      comm,
                               PDM_ownership_t owner)
  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  void PDM_mesh_location_n_part_cloud_set(PDM_mesh_location_t* ml,
                                          int                  i_point_cloud,
                                          int                  n_part)
  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  void PDM_mesh_location_cloud_set(PDM_mesh_location_t *ml,
                                   int                  i_point_cloud,
                                   int                  i_part,
                                   int                  n_points,
                                   double              *coords,
                                   PDM_g_num_t         *gnum)

  void PDM_mesh_location_cloud_get (PDM_mesh_location_t  *ml,
                                    int                   i_point_cloud,
                                    int                   i_part,
                                    int                  *n_points,
                                    double              **coords,
                                    PDM_g_num_t         **gnum)

  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  # void PDM_mesh_location_shared_nodal_mesh_set(int  id, PDM_Mesh_nodal_t *mesh_nodal)
  void PDM_mesh_location_mesh_global_data_set (PDM_mesh_location_t  *ml, int  n_part)
  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  void PDM_mesh_location_part_set(PDM_mesh_location_t  *ml,
                                  int                  i_part,
                                  int                  n_cell,
                                  int                 *cell_face_idx,
                                  int                 *cell_face,
                                  PDM_g_num_t         *cell_ln_to_gn,
                                  int                  n_face,
                                  int                 *face_vtx_idx,
                                  int                 *face_vtx,
                                  PDM_g_num_t         *face_ln_to_gn,
                                  int                  n_vtx,
                                  double              *coords,
                                  PDM_g_num_t         *vtx_ln_to_gn)
  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  void PDM_mesh_location_part_set_2d(PDM_mesh_location_t *ml,
                                     int                  i_part,
                                     int                  n_cell,
                                     int                 *cell_edge_idx,
                                     int                 *cell_edge,
                                     PDM_g_num_t         *cell_ln_to_gn,
                                     int                  n_edge,
                                     int                 *edge_vtx_idx,
                                     int                 *edge_vtx,
                                     PDM_g_num_t         *edge_ln_to_gn,
                                     int                  n_vtx,
                                     double              *coords,
                                     PDM_g_num_t         *vtx_ln_to_gn)
  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  void PDM_mesh_location_tolerance_set(PDM_mesh_location_t *ml, double tol)
  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  void PDM_mesh_location_method_set(PDM_mesh_location_t *ml, PDM_mesh_location_method_t method)
  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  void PDM_mesh_location_compute(PDM_mesh_location_t *ml)
  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  int PDM_mesh_location_n_located_get (PDM_mesh_location_t *ml,
                                       int                  i_point_cloud,
                                       int                  i_part)

  int PDM_mesh_location_n_unlocated_get (PDM_mesh_location_t *ml,
                                         int                  i_point_cloud,
                                         int                  i_part)

  int *PDM_mesh_location_unlocated_get (PDM_mesh_location_t *ml,
                                        int                  i_point_cloud,
                                        int                  i_part)


  int *PDM_mesh_location_located_get (PDM_mesh_location_t *ml,
                                      int                  i_point_cloud,
                                      int                  i_part)

  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  void PDM_mesh_location_point_location_get(PDM_mesh_location_t  *ml,
                                            int                   i_point_cloud,
                                            int                   i_part,
                                            PDM_g_num_t         **location,
                                            double              **dist2,
                                            double              **projected_coord)
  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  void PDM_mesh_location_points_in_elt_get (PDM_mesh_location_t  *ml,
                                            int                   i_part,
                                            int                   i_point_cloud,
                                            int                 **elt_pts_inside_idx,
                                            PDM_g_num_t         **points_gnum,
                                            double              **points_coords,
                                            double              **points_uvw,
                                            int                 **points_weights_idx,
                                            double              **points_weights,
                                            double              **points_dist2,
                                            double              **points_projected_coords)

  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  void PDM_mesh_location_free(PDM_mesh_location_t  *ml)
  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  void PDM_mesh_location_reverse_results_enable(PDM_mesh_location_t  *ml)
  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  void PDM_mesh_location_dump_times(PDM_mesh_location_t  *ml)
  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  void PDM_mesh_location_mesh_nodal_id_get(PDM_mesh_location_t  *ml)
  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  int PDM_mesh_location_n_cell_get (PDM_mesh_location_t  *ml, int i_part)

  void PDM_mesh_location_cell_vertex_get (PDM_mesh_location_t  *ml, int i_part, int **cell_vtx_idx, int **cell_vtx)


# ------------------------------------------------------------------
cdef class MeshLocation:
  """
     PointsMerge: Interface to build connection between multiple cloud in parallel
     Useful for get connection between partiton from faces coordinates
  """
  # ************************************************************************
  # > Class attributes
  cdef PDM_mesh_location_t* _ml
  cdef int _n_point_cloud
  cdef int _n_src_part
  cdef int _reverse_enabled
  cdef list _n_tgt_part_per_cloud

  cdef list _np_located, _np_unlocated, _dic_location, _dic_points_in_elt, _dic_cell_vertex
  # ************************************************************************

  # ------------------------------------------------------------------------
  def __init__(self, PDM_mesh_nature_t mesh_nature,
                     int               n_point_cloud,
                     MPI.Comm          comm,
                     bint              enable_reverse=True):
    """
    """

    # ::::::::::::::::::::::::::::::::::::::::::::::::::
    self._n_point_cloud = n_point_cloud
    self._n_tgt_part_per_cloud = [0 for i in range(n_point_cloud)]
    self._reverse_enabled = enable_reverse
    # ::::::::::::::::::::::::::::::::::::::::::::::::::

    # ::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Convert mpi4py -> PDM_MPI
    cdef MPI.MPI_Comm c_comm = comm.ob_mpi
    cdef PDM_MPI_Comm PDMC   = PDM_MPI_mpi_2_pdm_mpi_comm(<void *> &c_comm)
    # ::::::::::::::::::::::::::::::::::::::::::::::::::

    # ::::::::::::::::::::::::::::::::::::::::::::::::::
    self._ml = PDM_mesh_location_create(mesh_nature, n_point_cloud, PDMC, PDM_OWNERSHIP_UNGET_RESULT_IS_FREE)

    if self._reverse_enabled:
      PDM_mesh_location_reverse_results_enable (self._ml)

    # ::::::::::::::::::::::::::::::::::::::::::::::::::

  # ------------------------------------------------------------------------
  def n_part_cloud_set(self, int i_point_cloud,
                             int n_part):
    """
    """
    self._n_tgt_part_per_cloud[i_point_cloud] = n_part
    PDM_mesh_location_n_part_cloud_set(self._ml,
                                       i_point_cloud,
                                       n_part)

  # ------------------------------------------------------------------------
  def cloud_set(self, int i_point_cloud,
                      int i_part,
                      int n_points,
                      NPY.ndarray[NPY.double_t  , mode='c', ndim=1] coords,
                      NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] gnum,
                      ):
    """
    """
    PDM_mesh_location_cloud_set(self._ml,
                                i_point_cloud,
                                i_part,
                                n_points,
                 <double*>      coords.data,
                 <PDM_g_num_t*> gnum.data)

  # ------------------------------------------------------------------------
  # def nodal_mesh_set(self, MeshNodal mesh_nodal):
  #   """
  #   """
  #   PDM_mesh_location_shared_nodal_mesh_set(self._ml, mesh_nodal.mn)

  # ------------------------------------------------------------------------
  def mesh_global_data_set(self, int n_part):
    """
    """
    self._n_src_part = n_part
    PDM_mesh_location_mesh_global_data_set(self._ml, n_part)

  # ------------------------------------------------------------------------
  def part_set(self, int i_part,
                     int n_cell,
                     NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] cell_face_idx,
                     NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] cell_face,
                     NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] cell_ln_to_gn,
                     int n_face,
                     NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] face_vtx_idx,
                     NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] face_vtx,
                     NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] face_ln_to_gn,
                     int n_vtx,
                     NPY.ndarray[NPY.double_t  , mode='c', ndim=1] coords,
                     NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] vtx_ln_to_gn):
    """
    """
    PDM_mesh_location_part_set(self._ml,
                               i_part,
                               n_cell,
                <int*>         cell_face_idx.data,
                <int*>         cell_face.data,
                <PDM_g_num_t*> cell_ln_to_gn.data,
                               n_face,
                <int*>         face_vtx_idx.data,
                <int*>         face_vtx.data,
                <PDM_g_num_t*> face_ln_to_gn.data,
                               n_vtx,
                <double*>      coords.data,
                <PDM_g_num_t*> vtx_ln_to_gn.data)

  # ------------------------------------------------------------------------
  def part_set_2d(self, int i_part,
                     int n_cell,
                     NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] cell_edge_idx,
                     NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] cell_edge,
                     NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] cell_ln_to_gn,
                     int n_face,
                     NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] edge_vtx_idx,
                     NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] edge_vtx,
                     NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] edge_ln_to_gn,
                     int n_vtx,
                     NPY.ndarray[NPY.double_t  , mode='c', ndim=1] coords,
                     NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] vtx_ln_to_gn):
    """
    """
    PDM_mesh_location_part_set_2d(self._ml,
                                  i_part,
                                  n_cell,
                   <int*>         cell_edge_idx.data,
                   <int*>         cell_edge.data,
                   <PDM_g_num_t*> cell_ln_to_gn.data,
                                  n_face,
                   <int*>         edge_vtx_idx.data,
                   <int*>         edge_vtx.data,
                   <PDM_g_num_t*> edge_ln_to_gn.data,
                                  n_vtx,
                   <double*>      coords.data,
                   <PDM_g_num_t*> vtx_ln_to_gn.data)

  # ------------------------------------------------------------------------
  def tolerance_set(self, double tol):
    """
    """
    PDM_mesh_location_tolerance_set(self._ml, tol)

  # ------------------------------------------------------------------------
  def method_set(self, PDM_mesh_location_method_t method):
    """
    """
    PDM_mesh_location_method_set(self._ml, method)

  # ------------------------------------------------------------------------
  def __location_get(self, int i_point_cloud, int i_part):
    """
    """
    # ************************************************************************
    # > Declaration
    cdef int           n_points
    cdef int           n_located
    cdef double       *coords
    cdef PDM_g_num_t  *gnum
    cdef PDM_g_num_t  *location
    cdef double       *dist2
    cdef double       *p_proj_coord
    # ************************************************************************

    # Attention : Mettre les fonction n_located_ et located !!!!
    PDM_mesh_location_cloud_get (self._ml,
                                 i_point_cloud,
                                 i_part,
                                 &n_points,
                                 &coords,
                                 &gnum)

    n_located = PDM_mesh_location_n_located_get(self._ml, i_point_cloud, i_part)

    PDM_mesh_location_point_location_get(self._ml,
                                         i_point_cloud,
                                         i_part,
                                         &location,
                                         &dist2,
                                         &p_proj_coord)

    return {
            'g_num'        : create_numpy_pdm_gnum(gnum,         n_points, flag_owndata=False),
            'location'     : create_numpy_pdm_gnum(location,     n_located),
            'dist2'        : create_numpy_d       (dist2,        n_located),
            'p_proj_coord' : create_numpy_d       (p_proj_coord, 3*n_located)
           }

  def location_get(self, int i_point_cloud, int i_part):
    return self._dic_location[i_point_cloud][i_part]


  def __cell_vertex_get (self, int i_part):

    cdef int *cell_vtx_idx
    cdef int *cell_vtx

    cdef int n_elts = PDM_mesh_location_n_cell_get(self._ml, i_part)
    PDM_mesh_location_cell_vertex_get(self._ml, i_part, &cell_vtx_idx, &cell_vtx)

    return {'cell_vtx_idx'      : create_numpy_i(cell_vtx_idx, n_elts+1),
            'cell_vtx'          : create_numpy_i(cell_vtx, cell_vtx_idx[n_elts])
           }
  def cell_vertex_get (self, int i_part):
    return self._dic_cell_vertex[i_part]



  def __points_in_elt_get(self, int i_part, int i_point_cloud):
    """
    """
    # ************************************************************************
    # > Declaration
    cdef int          *elt_pts_inside_idx
    cdef PDM_g_num_t  *points_gnum
    cdef double       *points_coords
    cdef double       *points_uvw
    cdef int          *points_weights_idx
    cdef double       *points_weights
    cdef double       *points_dist2
    cdef double       *points_projected_coords
    # ************************************************************************

    PDM_mesh_location_points_in_elt_get(self._ml,
                                        i_part,
                                        i_point_cloud,
                                        &elt_pts_inside_idx,
                                        &points_gnum,
                                        &points_coords,
                                        &points_uvw,
                                        &points_weights_idx,
                                        &points_weights,
                                        &points_dist2,
                                        &points_projected_coords)

    cdef int n_elts =  PDM_mesh_location_n_cell_get(self._ml, i_part)
    cdef int s_loc  = elt_pts_inside_idx[n_elts]
    cdef int s_wei  = points_weights_idx[s_loc]

    return {'elt_pts_inside_idx'      : create_numpy_i       (elt_pts_inside_idx,      n_elts+1),
            'points_gnum'             : create_numpy_pdm_gnum(points_gnum,             s_loc   ),
            'points_coords'           : create_numpy_d       (points_coords,           3*s_loc ),
            'points_uvw'              : create_numpy_d       (points_uvw,              3*s_loc ),
            'points_weights_idx'      : create_numpy_i       (points_weights_idx,      s_loc+1 ),
            'points_weights'          : create_numpy_d       (points_weights,          s_wei   ),
            'points_dist2'            : create_numpy_d       (points_dist2,            s_loc   ),
            'points_projected_coords' : create_numpy_d       (points_projected_coords, 3*s_loc )
            }

  def points_in_elt_get(self, int i_part, int i_point_cloud):
    if self._reverse_enabled:
      return self._dic_points_in_elt[i_point_cloud][i_part]
    else:
      raise RuntimeError("Reverse mode was not enabled at class creation")


  def __located_get(self, int i_point_cloud, int i_part):
    """
    """
    cdef int n_located = PDM_mesh_location_n_located_get(self._ml, i_point_cloud, i_part)
    cdef int* located = PDM_mesh_location_located_get(self._ml, i_point_cloud, i_part)

    return create_numpy_i(located, n_located)

  def located_get(self, int i_point_cloud, int i_part):
    return self._np_located[i_point_cloud][i_part]

  def __unlocated_get(self, int i_point_cloud, int i_part):
    """
    """
    cdef int n_unlocated = PDM_mesh_location_n_unlocated_get(self._ml, i_point_cloud, i_part)
    cdef int* unlocated = PDM_mesh_location_unlocated_get(self._ml, i_point_cloud, i_part)

    return create_numpy_i(unlocated, n_unlocated)

  def unlocated_get(self, int i_point_cloud, int i_part):
    return self._np_unlocated[i_point_cloud][i_part]

  # ------------------------------------------------------------------------
  def compute(self):
    """
    """
    # ************************************************************************
    # > Declaration
    # ************************************************************************
    PDM_mesh_location_compute(self._ml)

    #Take ownership of all computed data to avoid memory leaks
    self._np_unlocated      = []
    self._np_located        = []
    self._dic_location      = []
    self._dic_points_in_elt = []
    for i_pt_cloud in range(self._n_point_cloud):
      #Target related data
      n_tgt_part = self._n_tgt_part_per_cloud[i_pt_cloud]
      self._np_unlocated.append([self.__unlocated_get(i_pt_cloud, i_part) for i_part in range(n_tgt_part)])
      self._np_located  .append([self.__located_get  (i_pt_cloud, i_part) for i_part in range(n_tgt_part)])
      self._dic_location.append([self.__location_get (i_pt_cloud, i_part) for i_part in range(n_tgt_part)])
      #Source related data
      if self._reverse_enabled:
        self._dic_points_in_elt.append([self.__points_in_elt_get(i_part, i_pt_cloud) for i_part in range(self._n_src_part)])
    #Source related data
    if self._reverse_enabled:
      self._dic_cell_vertex = [self.__cell_vertex_get(i_part) for i_part in range(self._n_src_part)]

  # ------------------------------------------------------------------------
  def dump_times(self):
    """
    """
    # ************************************************************************
    # > Declaration
    # ************************************************************************
    PDM_mesh_location_dump_times(self._ml)

  # ------------------------------------------------------------------------
  def __dealloc__(self):
    """
       Use the free method of PDM Lib
    """
    # ************************************************************************
    # > Declaration
    # ************************************************************************
    PDM_mesh_location_free(self._ml)
