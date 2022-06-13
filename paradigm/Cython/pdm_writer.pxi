cdef extern from "pdm_io.h":

  ctypedef enum PDM_io_kind_t:
  
    PDM_IO_KIND_MPIIO_EO   = 0
    PDM_IO_KIND_MPIIO_IP   = 1
    PDM_IO_KIND_MPI_SIMPLE = 2
    PDM_IO_KIND_SEQ        = 3

cdef extern from "pdm_mesh_nodal.h":

  ctypedef struct PDM_Mesh_nodal_t:
    pass

cdef extern from "pdm_writer.h":
  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  # > Wrapping of Ppart Structure
  ctypedef struct PDM_writer_t:
    pass
  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  # > Wrapping for C


  ctypedef enum PDM_writer_status_t:
    PDM_WRITER_OFF = 0
    PDM_WRITER_ON  = 1


  ctypedef enum PDM_writer_elt_geom_t:
    PDM_WRITER_POINT    = PDM_MESH_NODAL_POINT
    PDM_WRITER_BAR2     = PDM_MESH_NODAL_BAR2
    PDM_WRITER_TRIA3    = PDM_MESH_NODAL_TRIA3
    PDM_WRITER_QUAD4    = PDM_MESH_NODAL_QUAD4
    PDM_WRITER_POLY_2D  = PDM_MESH_NODAL_POLY_2D
    PDM_WRITER_TETRA4   = PDM_MESH_NODAL_TETRA4
    PDM_WRITER_PYRAMID5 = PDM_MESH_NODAL_PYRAMID5
    PDM_WRITER_PRISM6   = PDM_MESH_NODAL_PRISM6
    PDM_WRITER_HEXA8    = PDM_MESH_NODAL_HEXA8
    PDM_WRITER_POLY_3D  = PDM_MESH_NODAL_POLY_3D


  ctypedef enum PDM_writer_fmt_fic_t: 
    PDM_WRITER_FMT_BIN   = 0
    PDM_WRITER_FMT_ASCII = 1



  ctypedef enum PDM_writer_var_dim_t:
    PDM_WRITER_VAR_CST           = 0
    PDM_WRITER_VAR_SCALAR        = 1
    PDM_WRITER_VAR_VECTOR        = 3
    PDM_WRITER_VAR_TENSOR_SYM    = 6
    PDM_WRITER_VAR_TENSOR_ASYM   = 9

  ctypedef enum PDM_writer_var_loc_t:
    PDM_WRITER_VAR_VERTICES     = 0
    PDM_WRITER_VAR_ELEMENTS     = 1
    PDM_WRITER_VAR_PARTICLES    = 2

  ctypedef enum PDM_writer_topology_t:
    PDM_WRITER_TOPO_CST         = 0
    PDM_WRITER_TOPO_DEFORMABLE  = 1
    PDM_WRITER_TOPO_VARIABLE    = 2


  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  
  PDM_writer_t * PDM_writer_create(char                   *fmt,
                                   PDM_writer_fmt_fic_t    fmt_fic,
                                   PDM_writer_topology_t  topologie,
                                   PDM_writer_status_t     st_reprise,
                                   char                   *rep_sortie,
                                   char                   *nom_sortie,
                                   PDM_MPI_Comm            pdm_mpi_comm,
                                   PDM_io_kind_t          acces,
                                   double                  prop_noeuds_actifs,
                                   char                   *options)
  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  void PDM_writer_free(PDM_writer_t *cs)

  void PDM_writer_step_beg(PDM_writer_t  *cs, double   physical_time)

  void PDM_writer_step_end(PDM_writer_t  *cs)

  int PDM_writer_geom_create(PDM_writer_t               *cs,
                             char                 *nom_geom,
                             int                   n_part)

  int PDM_writer_geom_create_from_mesh_nodal(PDM_writer_t        *cs,
                                            char                *nom_geom,
                                            PDM_Mesh_nodal_t    *mesh)


  void PDM_writer_geom_coord_set(PDM_writer_t   *cs,
                                 int             id_geom,
                                 int             id_part,
                                 int             n_som,
                                 PDM_real_t      *coords,
                                 PDM_g_num_t     *numabs,
                                 PDM_ownership_t owner)


  void PDM_writer_geom_coord_from_parent_set(PDM_writer_t *cs,
                                             int           id_geom,
                                             int           id_part,
                                             int           n_som,
                                             int           n_som_parent,
                                             PDM_g_num_t  *numabs,
                                             int          *num_parent,
                                             PDM_real_t   *coords_parent,
                                             PDM_g_num_t  *numabs_parent)


  int PDM_writer_geom_bloc_add(PDM_writer_t          *cs,
                               int                    id_geom,
                               PDM_writer_status_t    st_free_data,
                               PDM_writer_elt_geom_t  t_elt)
  

  void PDM_writer_geom_bloc_std_set(PDM_writer_t  *cs,
                                    int      id_geom,
                                    int      id_bloc,
                                    int      id_part,
                                    int      n_elt,
                                    PDM_l_num_t   *connec,
                                    PDM_g_num_t   *numabs)
  

  void PDM_writer_geom_bloc_poly2d_set(PDM_writer_t  *cs,
                                       int            id_geom,
                                       int            id_bloc,
                                       int            id_part,
                                       PDM_l_num_t    n_elt,
                                       PDM_l_num_t   *connec_idx,
                                       PDM_l_num_t   *connec,
                                       PDM_g_num_t   *numabs)
  
  void PDM_writer_geom_bloc_poly3d_set(PDM_writer_t  *cs,
                                       int            id_geom,
                                       int            id_bloc,
                                       int            id_part,
                                       PDM_l_num_t    n_elt,
                                       PDM_l_num_t    n_face,
                                       PDM_l_num_t   *facsom_idx,
                                       PDM_l_num_t   *facsom,
                                       PDM_l_num_t   *cellfac_idx,
                                       PDM_l_num_t   *cellfac,
                                       PDM_g_num_t   *numabs)
  
  void PDM_writer_geom_cell3d_cellface_add(PDM_writer_t *cs,
                                           int           id_geom,
                                           int           id_part,
                                           int           n_cell,
                                           int           n_face,
                                           PDM_l_num_t  *face_som_idx,
                                           PDM_l_num_t  *face_som_nb,
                                           PDM_l_num_t  *face_som,
                                           PDM_l_num_t  *cell_face_idx,
                                           PDM_l_num_t  *cell_face_nb,
                                           PDM_l_num_t  *cell_face,
                                           PDM_g_num_t  *numabs)
  

  void PDM_writer_geom_cell2d_cellface_add(PDM_writer_t *cs,
                                           int           id_geom,
                                           int           id_part,
                                           int           n_cell,
                                           int           n_face,
                                           PDM_l_num_t  *face_som_idx,
                                           PDM_l_num_t  *face_som_nb,
                                           PDM_l_num_t  *face_som,
                                           PDM_l_num_t  *cell_face_idx,
                                           PDM_l_num_t  *cell_face_nb,
                                           PDM_l_num_t  *cell_face,
                                           PDM_g_num_t  *numabs)

  
  void PDM_writer_geom_faces_facesom_add(PDM_writer_t *cs,
                                         int     id_geom,
                                         int     id_part,
                                         int     n_face,
                                         PDM_l_num_t  *face_som_idx,
                                         PDM_l_num_t  *face_som_nb,
                                         PDM_l_num_t  *face_som,
                                         PDM_g_num_t  *numabs)

  
  void PDM_writer_geom_write(PDM_writer_t *cs, int id_geom)

  
  void PDM_writer_geom_data_free(PDM_writer_t *cs, int id_geom)
  

  void PDM_writer_geom_free(PDM_writer_t *cs, int id_geom)

  
  int PDM_writer_var_create(PDM_writer_t         *cs,
                            PDM_writer_status_t   st_dep_tps,
                            PDM_writer_var_dim_t  dim,
                            PDM_writer_var_loc_t  loc,
                            char                 *nom_var)
  
  void PDM_writer_name_map_add(PDM_writer_t *cs,
                               char   *public_name,
                               char   *private_name)
  
  void PDM_writer_var_write(PDM_writer_t *cs, int id_var)
  
  void PDM_writer_var_set(PDM_writer_t  *cs,
                          int           id_var,
                          int           id_geom,
                          int           id_part,
                          PDM_real_t   *val)
  
  void PDM_writer_var_data_free(PDM_writer_t *cs,
                                int           id_var)
  
  void PDM_writer_var_free(PDM_writer_t *cs, int id_var)
  
  # void PDM_writer_fmt_add(char                  *name,            
  #                         PDM_writer_fct_t       create_fct,      
  #                         PDM_writer_fct_t       free_fct,        
  #                         PDM_writer_fct_t       beg_step_fct,    
  #                         PDM_writer_fct_t       end_step_fct,    
  #                         PDM_writer_geom_fct_t  geom_create_fct, 
  #                         PDM_writer_geom_fct_t  geom_write_fct,  
  #                         PDM_writer_geom_fct_t  geom_free_fct,   
  #                         PDM_writer_var_fct_t   var_create_fct,  
  #                         PDM_writer_var_fct_t   var_write_fct,   
  #                         PDM_writer_var_fct_t   var_free_fct)
  
  # void PDM_writer_fmt_free ()
  
  void PDM_writer_geom_data_reset(PDM_writer_t *cs, int     id_geom)
 
# ------------------------------------------------------------------

cdef class Writer:
  """
  """

  cdef PDM_writer_t* _wt

  def __cinit__(self,                  char  *fmt,
                PDM_writer_fmt_fic_t   fmt_fic,
                PDM_writer_topology_t  topologie,
                PDM_writer_status_t    st_reprise,
                char                  *rep_sortie,
                char                  *nom_sortie,
                MPI.Comm               comm,
                PDM_io_kind_t         acces,
                NPY.double_t           prop_noeuds_actifs,
                char                  *options):
      """
      """

      cdef MPI.MPI_Comm c_comm = comm.ob_mpi
      cdef PDM_MPI_Comm PDMC   = PDM_MPI_mpi_2_pdm_mpi_comm(<void *> &c_comm)

      self._wt = PDM_writer_create(fmt,
                                   fmt_fic,
                                   topologie,
                                   st_reprise,
                                   rep_sortie,
                                   nom_sortie,
                                   PDMC,
                                   acces,
                                   prop_noeuds_actifs,
                                   options)

  def geom_create(self,
                  char  *nom_geom,
                  int    n_part):
      """
      """

      return PDM_writer_geom_create(self._wt,
                                    nom_geom,
                                    n_part)


  # def geom_create_from_mesh_nodal(self,   
  #                                 char                *nom_geom,
  #                                 PDM_Mesh_nodal_t    *mesh):
  #     """
  #     """

  #     return PDM_writer_geom_create_from_mesh_nodal(self._wt,
  #                                                   nom_geom,
  #                                                   mesh)

  def geom_cell2d_cellface_add(self,
                               int id_geom,
                               int id_part,
                               int n_cell,
                               int n_face,
                               NPY.ndarray[NPY.int32_t  , mode='c', ndim=1] face_som_idx,
                               NPY.ndarray[NPY.int32_t  , mode='c', ndim=1] face_som_nb,
                               NPY.ndarray[NPY.int32_t  , mode='c', ndim=1] face_som,
                               NPY.ndarray[NPY.int32_t  , mode='c', ndim=1] cell_face_idx,
                               NPY.ndarray[NPY.int32_t  , mode='c', ndim=1] cell_face_nb,
                               NPY.ndarray[NPY.int32_t  , mode='c', ndim=1] cell_face,
                               NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] numabs):
      """
      """

      PDM_writer_geom_cell2d_cellface_add(self._wt,
                                          id_geom,
                                          id_part,
                                          n_cell,
                                          n_face,
                                          <int *> face_som_idx.data,
                                          <int *> face_som_nb.data,
                                          <int *> face_som.data,
                                          <int *> cell_face_idx.data,
                                          <int *> cell_face_nb.data,
                                          <int *> cell_face.data,
                                          <PDM_g_num_t*> numabs.data)

  def geom_cell3d_cellface_add(self,
                               int id_geom,
                               int id_part,
                               int n_cell,
                               int n_face,
                               NPY.ndarray[NPY.int32_t  , mode='c', ndim=1] face_som_idx,
                               NPY.ndarray[NPY.int32_t  , mode='c', ndim=1] face_som_nb,
                               NPY.ndarray[NPY.int32_t  , mode='c', ndim=1] face_som,
                               NPY.ndarray[NPY.int32_t  , mode='c', ndim=1] cell_face_idx,
                               NPY.ndarray[NPY.int32_t  , mode='c', ndim=1] cell_face_nb,
                               NPY.ndarray[NPY.int32_t  , mode='c', ndim=1] cell_face,
                               NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] numabs):
      """
      """

      PDM_writer_geom_cell3d_cellface_add(self._wt,
                                          id_geom,
                                          id_part,
                                          n_cell,
                                          n_face,
                                          <int *> face_som_idx.data,
                                          <int *> face_som_nb.data,
                                          <int *> face_som.data,
                                          <int *> cell_face_idx.data,
                                          <int *> cell_face_nb.data,
                                          <int *> cell_face.data,
                                          <PDM_g_num_t*> numabs.data)

  def geom_coord_set(self,
                     int            id_geom,
                     int            id_part,
                     int            n_som,
                     NPY.ndarray[NPY.double_t, mode='c', ndim=1] coords,
                     NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] numabs):
      """
      """

      PDM_writer_geom_coord_set(self._wt,
                               id_geom,
                               id_part,
                               n_som,
                               <PDM_real_t *> coords.data,
                               <PDM_g_num_t *> numabs.data,
                               PDM_OWNERSHIP_USER)
  

  def geom_faces_facesom_add(self,
                             int     id_geom,
                             int     id_part,
                             int     n_face,
                             NPY.ndarray[NPY.int_t  , mode='c', ndim=1] face_som_idx,
                             NPY.ndarray[NPY.int_t  , mode='c', ndim=1] face_som_nb,
                             NPY.ndarray[NPY.int_t  , mode='c', ndim=1] face_som,
                             NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] numabs):
      """
      """

      PDM_writer_geom_faces_facesom_add(self._wt,
                                        id_geom,
                                        id_part,
                                        n_face,
                                        <int *> face_som_idx.data,
                                        <int *> face_som_nb.data,
                                        <int *> face_som.data,
                                        <PDM_g_num_t*> numabs.data)

  def geom_write(self, int id_geom):
      """
      """

      PDM_writer_geom_write(self._wt, id_geom)

  
  def geom_data_free(self, int id_geom):
      """
      """
  
      PDM_writer_geom_data_free(self._wt, id_geom)


  def geom_free(self, int id_geom):
      """
      """
      PDM_writer_geom_free(self._wt, id_geom)


  def var_create(self,
                 PDM_writer_status_t   st_dep_tps,
                 PDM_writer_var_dim_t  dim,
                 PDM_writer_var_loc_t  loc,
                 char                 *nom_var):
      """
      """
  
      return PDM_writer_var_create(self._wt,
                                   st_dep_tps,
                                   dim,
                                   loc,
                                   nom_var)

  def name_map_add(self,
                   char   *public_name,
                   char   *private_name):
      """
      """
  
      PDM_writer_name_map_add(self._wt,
                              public_name,
                              private_name)

  def var_write(self, int id_var):
      """
      """
  
      PDM_writer_var_write(self._wt, id_var)


  def var_set(self,
              int           id_var,
              int           id_geom,
              int           id_part,
              NPY.ndarray[NPY.double_t  , mode='c', ndim=1] val):
      """
      """
  
      PDM_writer_var_set(self._wt,
                         id_var,
                         id_geom,
                         id_part,
                         <double *> val.data)

  def var_data_free(self,
                    int           id_var):
  
      """
      """
      PDM_writer_var_data_free(self._wt,
                              id_var)


  def var_free(self, int id_var):
      """
      """

      PDM_writer_var_free(self._wt, id_var)

  def step_beg(self, double t):
    PDM_writer_step_beg(self._wt, t)

  def step_end(self):
    PDM_writer_step_end(self._wt);

  # ------------------------------------------------------------------------

  def __dealloc__(self):
    """
       Use the free method of PDM Lib
    """
    # ************************************************************************
    # > Declaration
    # ************************************************************************

    PDM_writer_free(self._wt)
