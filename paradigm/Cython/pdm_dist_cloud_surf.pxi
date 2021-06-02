
cdef extern from "pdm_dist_cloud_surf.h":
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Wrapping of Ppart Structure
    ctypedef struct PDM_dist_cloud_surf_t:
        pass
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Wrapping of function
    PDM_dist_cloud_surf_t* PDM_dist_cloud_surf_create(PDM_mesh_nature_t mesh_nature,
                                                      int               n_point_cloud,
                                                      PDM_MPI_Comm      comm,
                                                      PDM_ownership_t   owner)

    void PDM_dist_cloud_surf_n_part_cloud_set(PDM_dist_cloud_surf_t *dist,
                                              int                    i_point_cloud,
                                              int                    n_part)

    void PDM_dist_cloud_surf_cloud_set(PDM_dist_cloud_surf_t *dist,
                                       int                    i_point_cloud,
                                       int                    i_part,
                                       int                    n_points,
                                       double                *coords,
                                       PDM_g_num_t           *gnum)

#    void PDM_dist_cloud_surf_nodal_mesh_set(PDM_dist_cloud_surf_t *dist,
#                                            int                    mesh_nodal_id)

#    void PDM_dist_cloud_surf_surf_mesh_map(PDM_dist_cloud_surf_t *dist,
#                                           PDM_surf_mesh_t       *surf_mesh)

    void PDM_dist_cloud_surf_surf_mesh_global_data_set(PDM_dist_cloud_surf_t *dist,
                                                       PDM_g_num_t            n_g_face,
                                                       PDM_g_num_t            n_g_vtx,
                                                       int                    n_part)

    void PDM_dist_cloud_surf_surf_mesh_part_set(PDM_dist_cloud_surf_t *dist,
                                                int                    i_part,
                                                int                     n_face,
                                                int                   *face_vtx_idx,
                                                int                   *face_vtx,
                                                PDM_g_num_t           *face_ln_to_gn,
                                                int                    n_vtx,
                                                double                *coords,
                                                PDM_g_num_t           *vtx_ln_to_gn)

    void PDM_dist_cloud_surf_compute(PDM_dist_cloud_surf_t *dist)

    void PDM_dist_cloud_surf_get(PDM_dist_cloud_surf_t   *dist,
                                 int                      i_point_cloud,
                                 int                      i_part,
                                 double                **closest_elt_distance,
                                 double                **closest_elt_projected,
                                 PDM_g_num_t           **closest_elt_gnum)

    void PDM_dist_cloud_surf_free(PDM_dist_cloud_surf_t   *dist)

    void PDM_dist_cloud_surf_dump_times(PDM_dist_cloud_surf_t   *dist)


    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# ------------------------------------------------------------------
cdef class DistCloudSurf:
    """
    Define a method to compute the distance from a cloud to a surface
    """
    # ************************************************************************
    # > Class attributes
    cdef PDM_dist_cloud_surf_t*  _dist
    cdef int  _n_point_cloud
    cdef int  **_nb_pts
    cdef int  _part_n
    # ************************************************************************
    # ------------------------------------------------------------------
    def __cinit__(self,
                  PDM_mesh_nature_t mesh_nature,
                  int n_point_cloud,
                  MPI.Comm  comm):
        """
        Compute the distance from point clouds to a surface
        """
        # ::::::::::::::::::::::::::::::::::::::::::::::::::
        # > Convert mpi4py -> PDM_MPI
        cdef MPI.MPI_Comm c_comm   = comm.ob_mpi
        cdef PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&c_comm)

        self._dist = PDM_dist_cloud_surf_create(mesh_nature,
                                                n_point_cloud,
                                                pdm_comm,
                                                PDM_OWNERSHIP_USER) # Python take ownership)
        self._n_point_cloud = n_point_cloud
        self._nb_pts = <int **>  malloc(sizeof(int *) * n_point_cloud)

    # ------------------------------------------------------------------
    def n_part_cloud_set(self,
                         int i_point_cloud,
                         int n_part):
        """
        Give the number of partitions of a point cloud
        """
        PDM_dist_cloud_surf_n_part_cloud_set(self._dist, i_point_cloud, n_part)
        self._nb_pts[i_point_cloud] =  <int *>  malloc(sizeof(int) * n_part)
        self._part_n = n_part

    # ------------------------------------------------------------------
    def cloud_set(self,
                  int i_point_cloud,
                  int i_part,
                  int n_points,
                  NPY.ndarray[NPY.double_t   , mode='c', ndim=1] coords not None,
                  NPY.ndarray[npy_pdm_gnum_t , mode='c', ndim=1] gnum not None):
        """
        Give the properties of a partition of a point cloud
        """
        PDM_dist_cloud_surf_cloud_set(self._dist,
                                      i_point_cloud,
                                      i_part,
                                      n_points,
                                      <double *> coords.data,
                                      <PDM_g_num_t *> gnum.data)
        self._nb_pts[i_point_cloud][i_part] = n_points

    # ------------------------------------------------------------------
    def surf_mesh_global_data_set(self,
                                  PDM_g_num_t n_g_face,
                                  PDM_g_num_t n_g_vtx,
                                  int n_part):
        """
        Give the global properties of the surface mesh
        """
        PDM_dist_cloud_surf_surf_mesh_global_data_set(self._dist,
                                                      n_g_face,
                                                      n_g_vtx,
                                                      n_part)

    # ------------------------------------------------------------------
    def surf_mesh_part_set(self,
                           i_part,
                           n_face,
                           NPY.ndarray[NPY.int32_t    , mode='c', ndim=1] face_vtx_idx not None,
                           NPY.ndarray[NPY.int32_t    , mode='c', ndim=1] face_vtx not None,
                           NPY.ndarray[npy_pdm_gnum_t , mode='c', ndim=1] face_ln_to_gn not None,
                           n_vtx,
                           NPY.ndarray[NPY.double_t   , mode='c', ndim=1] coords not None,
                           NPY.ndarray[npy_pdm_gnum_t , mode='c', ndim=1] vtx_ln_to_gn not None):
        """
        Give the properties of a partition of the surface mesh
        """
        PDM_dist_cloud_surf_surf_mesh_part_set (self._dist,
                                                i_part,
                                                n_face,
                                                <int *> face_vtx_idx.data,
                                                <int *> face_vtx.data,
                                                <PDM_g_num_t *> face_ln_to_gn.data,
                                                n_vtx,
                                                <double *> coords.data,
                                                <PDM_g_num_t *> vtx_ln_to_gn.data)

    # ------------------------------------------------------------------
    def compute(self):
        """
        Compute distance
        """
        PDM_dist_cloud_surf_compute (self._dist)

    # ------------------------------------------------------------------
    def get(self,
           int i_point_cloud,
           int i_part):
        """
        Return the properties of the closest surface element
        (distance, projected point coordinates and global numbering)
        """

        cdef double *closest_elt_distance
        cdef double *closest_elt_projected
        cdef PDM_g_num_t *closest_elt_gnum

        PDM_dist_cloud_surf_get (self._dist,
                                 i_point_cloud,
                                 i_part,
                                 &closest_elt_distance,
                                 &closest_elt_projected,
                                 &closest_elt_gnum)

        if (closest_elt_distance == NULL) :
            npClosestEltDistance = None
        else :
            dim = <NPY.npy_intp> self._nb_pts[i_point_cloud][i_part]
            npClosestEltDistance = NPY.PyArray_SimpleNewFromData(1,
                                                                 &dim,
                                                                 NPY.NPY_DOUBLE,
                                                                 <void *> closest_elt_distance)
            PyArray_ENABLEFLAGS(npClosestEltDistance, NPY.NPY_OWNDATA);

        if (closest_elt_projected == NULL) :
            npClosestEltProjected = None
        else :
            dim = <NPY.npy_intp> (3 * self._nb_pts[i_point_cloud][i_part])
            npClosestEltProjected = NPY.PyArray_SimpleNewFromData(1,
                                                                  &dim,
                                                                  NPY.NPY_DOUBLE,
                                                                  <void *> closest_elt_projected)
            PyArray_ENABLEFLAGS(npClosestEltProjected, NPY.NPY_OWNDATA);

        if (closest_elt_gnum == NULL) :
            npClosestEltGnum = None
        else :
            dim = <NPY.npy_intp> self._nb_pts[i_point_cloud][i_part]
            npClosestEltGnum = NPY.PyArray_SimpleNewFromData(1,
                                                             &dim,
                                                             PDM_G_NUM_NPY_INT,
                                                             <void *> closest_elt_gnum)
            PyArray_ENABLEFLAGS(npClosestEltGnum, NPY.NPY_OWNDATA);

        return {'ClosestEltDistance'    : npClosestEltDistance,
                'ClosestEltProjected'   : npClosestEltProjected,
                'ClosestEltGnum'        : npClosestEltGnum}

    # ------------------------------------------------------------------
    def dump_times(self):
        """
        Print elapsed time
        """
        PDM_dist_cloud_surf_dump_times(self._dist)

    # ------------------------------------------------------------------
    def __dealloc__(self):
        """
        """
        for idx in xrange(self._part_n):
            free (self._nb_pts[idx])
        free (self._nb_pts)
        PDM_dist_cloud_surf_free(self._dist)
