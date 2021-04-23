
cdef extern from "pdm_mesh_location.h":

  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  # > Wrapping for C
  ctypedef enum PDM_mesh_location_method_t:
    PDM_MESH_LOCATION_OCTREE  = 0
    PDM_MESH_LOCATION_DBBTREE = 1
  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  int PDM_mesh_location_create(PDM_mesh_nature_t mesh_nature,
                               int               n_point_cloud,
                               PDM_MPI_Comm      comm);
  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  void PDM_mesh_location_n_part_cloud_set(int          id,
                                          int          i_point_cloud,
                                          int          n_part);
  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  void PDM_mesh_location_cloud_set(int          id,
                                   int          i_point_cloud,
                                   int          i_part,
                                   int          n_points,
                                   double      *coords,
                                   PDM_g_num_t *gnum);

  void PDM_mesh_location_cloud_get (int           id,
                                    int           i_point_cloud,
                                    int           i_part,
                                    int          *n_points,
                                    double      **coords,
                                    PDM_g_num_t **gnum);

  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  # void PDM_mesh_location_shared_nodal_mesh_set(int  id, PDM_Mesh_nodal_t *mesh_nodal);
  void PDM_mesh_location_mesh_global_data_set (int  id, int  n_part);
  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  void PDM_mesh_location_part_set(int          id,
                                  int          i_part,
                                  int          n_cell,
                                  int         *cell_face_idx,
                                  int         *cell_face,
                                  PDM_g_num_t *cell_ln_to_gn,
                                  int          n_face,
                                  int         *face_vtx_idx,
                                  int         *face_vtx,
                                  PDM_g_num_t *face_ln_to_gn,
                                  int          n_vtx,
                                  double      *coords,
                                  PDM_g_num_t *vtx_ln_to_gn);
  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  void PDM_mesh_location_part_set_2d(int          id,
                                     int          i_part,
                                     int          n_cell,
                                     int         *cell_edge_idx,
                                     int         *cell_edge,
                                     PDM_g_num_t *cell_ln_to_gn,
                                     int          n_edge,
                                     int         *edge_vtx_idx,
                                     int         *edge_vtx,
                                     PDM_g_num_t *edge_ln_to_gn,
                                     int          n_vtx,
                                     double      *coords,
                                     PDM_g_num_t *vtx_ln_to_gn);
  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  void PDM_mesh_location_tolerance_set(int    id, double tol);
  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  void PDM_mesh_location_method_set(int id, PDM_mesh_location_method_t method);
  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  void PDM_mesh_location_compute(int id);
  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  int PDM_mesh_location_n_located_get (int id,
                                       int i_point_cloud,
                                       int i_part);

  int PDM_mesh_location_n_unlocated_get (int  id,
                                         int  i_point_cloud,
                                         int  i_part);

  int *PDM_mesh_location_unlocated_get (int id,
                                        int i_point_cloud,
                                        int i_part);


  int *PDM_mesh_location_located_get (int id,
                                      int i_point_cloud,
                                      int i_part);

  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  void PDM_mesh_location_point_location_get(int           id,
                                            int           i_point_cloud,
                                            int           i_part,
                                            PDM_g_num_t **location,
                                            double      **dist2,
                                            double      **projected_coord);
  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  void PDM_mesh_location_points_in_elt_get (int     id,
                                            int     i_part,
                                            int     i_point_cloud,
                                            int         **elt_pts_inside_idx,
                                            PDM_g_num_t **points_gnum,
                                            double      **points_coords,
                                            double      **points_uvw,
                                            int         **points_weights_idx,
                                            double      **points_weights,
                                            double      **points_dist2,
                                            double      **points_projected_coords);

  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  void PDM_mesh_location_free(int id, int partial);
  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  void PDM_mesh_location_dump_times(int id);
  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  void PDM_mesh_location_mesh_nodal_id_get(int id);
  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  int PDM_mesh_location_n_cell_get (int id, int i_part);

# ------------------------------------------------------------------
cdef class MeshLocation:
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
  def __init__(self, PDM_mesh_nature_t mesh_nature,
                     int               n_point_cloud,
                     MPI.Comm          comm):
    """
    """

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
    self._id = PDM_mesh_location_create(mesh_nature, n_point_cloud, PDMC)
    # ::::::::::::::::::::::::::::::::::::::::::::::::::

  # ------------------------------------------------------------------------
  def n_part_cloud_set(self, int i_point_cloud,
                             int n_part):
    """
    """
    PDM_mesh_location_n_part_cloud_set(self._id,
                                       i_point_cloud,
                                       n_part);

  # ------------------------------------------------------------------------
  def cloud_set(self, int i_point_cloud,
                      int i_part,
                      int n_points,
                      NPY.ndarray[NPY.double_t  , mode='c', ndim=1] coords,
                      NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] gnum,
                      ):
    """
    """
    PDM_mesh_location_cloud_set(self._id,
                                i_point_cloud,
                                i_part,
                                n_points,
                 <double*>      coords.data,
                 <PDM_g_num_t*> gnum.data);

  # ------------------------------------------------------------------------
  # def nodal_mesh_set(self, MeshNodal mesh_nodal):
  #   """
  #   """
  #   PDM_mesh_location_shared_nodal_mesh_set(self._id, mesh_nodal.mn);

  # ------------------------------------------------------------------------
  def mesh_global_data_set(self, int n_part):
    """
    """
    PDM_mesh_location_mesh_global_data_set(self._id, n_part);

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
    PDM_mesh_location_part_set(self._id,
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
                <PDM_g_num_t*> vtx_ln_to_gn.data);

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
    PDM_mesh_location_part_set_2d(self._id,
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
                   <PDM_g_num_t*> vtx_ln_to_gn.data);

  # ------------------------------------------------------------------------
  def tolerance_set(self, double tol):
    """
    """
    PDM_mesh_location_tolerance_set(self._id, tol);

  # ------------------------------------------------------------------------
  def method_set(self, PDM_mesh_location_method_t method):
    """
    """
    PDM_mesh_location_method_set(self._id, method);

  # ------------------------------------------------------------------------
  def location_get(self, int i_point_cloud,
                         int i_part):
    """
    """
    # ************************************************************************
    # > Declaration
    cdef int           n_points
    cdef double       *coords
    cdef PDM_g_num_t  *gnum
    cdef PDM_g_num_t  *location
    cdef double       *dist2
    cdef double       *p_proj_coord
    # ************************************************************************


    # Attention : Mettre les fonction n_located_ et located !!!!  



    PDM_mesh_location_cloud_get (self._id,
                                 i_point_cloud,
                                 i_part,
                                 &n_points,
                                 &coords,
                                 &gnum);

    PDM_mesh_location_point_location_get(self._id,
                                         i_point_cloud,
                                         i_part,
                                         &location,
                                         &dist2,
                                         &p_proj_coord);

    cdef NPY.npy_intp dim
    # > Build numpy capsule
    dim = <NPY.npy_intp> n_points
    np_g_num = NPY.PyArray_SimpleNewFromData(1,
                                             &dim,
                                             PDM_G_NUM_NPY_INT,
                                             <void *> gnum)
    # PyArray_ENABLEFLAGS(np_g_num, NPY.NPY_OWNDATA)

    np_location = NPY.PyArray_SimpleNewFromData(1,
                                             &dim,
                                             PDM_G_NUM_NPY_INT,
                                             <void *> location)
    # PyArray_ENABLEFLAGS(np_location, NPY.NPY_OWNDATA)

    np_dist2 = NPY.PyArray_SimpleNewFromData(1,
                                             &dim,
                                             NPY.NPY_DOUBLE,
                                             <void *> dist2)
    # > Build numpy capsule
    dim = <NPY.npy_intp> 3*n_points
    np_p_proj_coord = NPY.PyArray_SimpleNewFromData(1,
                                              &dim,
                                              NPY.NPY_DOUBLE,
                                              <void *> p_proj_coord)
    # PyArray_ENABLEFLAGS(np_location, NPY.NPY_OWNDATA)

    return {'g_num'        : np_g_num,
            'location'     : np_location,
            'dist2'        : np_dist2,
            'p_proj_coord' : np_p_proj_coord
            }


  def points_in_elt_get(self, int i_part, int i_point_cloud):
    """
    """
    # ************************************************************************
    # > Declaration
    cdef int           n_points
    cdef int           n_elts
    cdef double       *coords
    cdef PDM_g_num_t  *gnum
    cdef PDM_g_num_t  *location
    cdef double       *dist2
    cdef double       *p_proj_coord
    # ************************************************************************

    n_elts =  PDM_mesh_location_n_cell_get(self._id, i_part);


    PDM_mesh_location_points_in_elt_get(self._id, 
                                        i_part,
                                        i_point_cloud,
 int         **elt_pts_inside_idx,
 PDM_g_num_t **points_gnum,
 double      **points_coords,
 double      **points_uvw,
 int         **points_weights_idx,
 double      **points_weights,
 double      **points_dist2,
 double      **points_projected_coords
);



    cdef NPY.npy_intp dim
    # > Build numpy capsule
    dim = <NPY.npy_intp> n_points
    np_g_num = NPY.PyArray_SimpleNewFromData(1,
                                             &dim,
                                             PDM_G_NUM_NPY_INT,
                                             <void *> gnum)
    # PyArray_ENABLEFLAGS(np_g_num, NPY.NPY_OWNDATA)

    np_location = NPY.PyArray_SimpleNewFromData(1,
                                             &dim,
                                             PDM_G_NUM_NPY_INT,
                                             <void *> location)
    # PyArray_ENABLEFLAGS(np_location, NPY.NPY_OWNDATA)

    np_dist2 = NPY.PyArray_SimpleNewFromData(1,
                                             &dim,
                                             NPY.NPY_DOUBLE,
                                             <void *> dist2)
    # > Build numpy capsule
    dim = <NPY.npy_intp> 3*n_points
    np_p_proj_coord = NPY.PyArray_SimpleNewFromData(1,
                                              &dim,
                                              NPY.NPY_DOUBLE,
                                              <void *> p_proj_coord)
    # PyArray_ENABLEFLAGS(np_location, NPY.NPY_OWNDATA)

    return {'g_num'        : np_g_num,
            'location'     : np_location,
            'dist2'        : np_dist2,
            'p_proj_coord' : np_p_proj_coord
            }





  def located_get(self, int i_point_cloud, int i_part):
    """
    """
    cdef int           n_located
    cdef double       *located

    n_located = PDM_mesh_location_n_located_get(self._id,
                                                i_point_cloud,
                                                i_part);

    located = PDM_mesh_location_located_get(self._id,
                                              i_point_cloud,
                                              i_part);
    dim = <NPY.npy_intp> 3*n_located
    np_located = NPY.PyArray_SimpleNewFromData(1,
                                               &dim,
                                               NPY.NPY_INT32,
                                               <void *> located)
    return np_located

  def unlocated_get(self, int i_point_cloud, int i_part):
    """
    """
    cdef int           n_unlocated
    cdef double       *unlocated

    n_unlocated = PDM_mesh_location_n_unlocated_get(self._id,
                                                i_point_cloud,
                                                i_part);

    unlocated = PDM_mesh_location_unlocated_get(self._id,
                                              i_point_cloud,
                                              i_part);
    dim = <NPY.npy_intp> 3*n_unlocated
    np_unlocated = NPY.PyArray_SimpleNewFromData(1,
                                               &dim,
                                               NPY.NPY_INT32,
                                               <void *> unlocated)
    return np_unlocated

  # ------------------------------------------------------------------------
  def compute(self):
    """
    """
    # ************************************************************************
    # > Declaration
    # ************************************************************************
    PDM_mesh_location_compute(self._id)

  # ------------------------------------------------------------------------
  def dump_times(self):
    """
    """
    # ************************************************************************
    # > Declaration
    # ************************************************************************
    PDM_mesh_location_dump_times(self._id)

  # ------------------------------------------------------------------------
  def __dealloc__(self):
    """
       Use the free method of PDM Lib
    """
    # ************************************************************************
    # > Declaration
    # ************************************************************************
    print('PDM_mesh_location_free')
    PDM_mesh_location_free(self._id, 0)
