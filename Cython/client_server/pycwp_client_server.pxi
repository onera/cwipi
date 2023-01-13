cdef extern from "client.h":

  # --> general functions
  void CWP_client_Init(MPI.MPI_Comm       global_comm,
                       int                n_code,
                       char             **code_names,
                       CWP_Status_t      *is_active_rank,
                       double            *time_init)
  void CWP_client_Finalize()

  # --> functions about current code properties
  void CWP_client_State_update(char* local_code_name,
                               CWP_State_t state)
  void CWP_client_Time_update(char* local_code_name,
                              double current_time)
  void CWP_client_Output_file_set(FILE *output_file)

  # --> functions about other code properties
  CWP_State_t CWP_client_State_get(char    *code_name)
  int CWP_client_Codes_nb_get()
  char **CWP_client_Codes_list_get()
  int CWP_client_Loc_codes_nb_get()
  char **CWP_client_Loc_codes_list_get()

  # --> functions about properties
  void CWP_client_Properties_dump()

  # --> general functions about coupling
  void CWP_client_Cpl_create(char                *local_code_name,
                             char                *cpl_id,
                             char                *coupled_code_name,
                             CWP_Interface_t      entities_dim,
                             CWP_Comm_t           comm_type,
                             CWP_Spatial_interp_t spatial_interp,
                             int                  n_part,
                             CWP_Dynamic_mesh_t   displacement,
                             CWP_Time_exch_t      recv_freq_type)
  void CWP_client_Cpl_del(char *local_code_name,
                          char *cpl_id)
  int CWP_client_N_uncomputed_tgts_get(char *local_code_name,
                                       char *cpl_id,
                                       char *field_id,
                                       int   i_part)
  int *CWP_client_Uncomputed_tgts_get(char *local_code_name,
                                      char *cpl_id,
                                      char *field_id,
                                      int   i_part)
  int CWP_client_N_computed_tgts_get(char *local_code_name,
                                     char *cpl_id,
                                     char *field_id,
                                     int   i_part)
  int *CWP_client_Computed_tgts_get(char *local_code_name,
                                    char *cpl_id,
                                    char *field_id,
                                    int   i_part)
  int CWP_client_N_involved_srcs_get(char *local_code_name,
                                     char *cpl_id,
                                     char *field_id,
                                     int   i_part)
  int *CWP_client_Involved_srcs_get(char *local_code_name,
                                    char *cpl_id,
                                    char *field_id,
                                    int   i_part)

  # --> functions about spatial interpolation
  void CWP_client_Spatial_interp_weights_compute(char     *local_code_name,
                                                 char     *cpl_id)
  void CWP_client_Spatial_interp_property_set(char     *local_code_name,
                                              char     *cpl_id,
                                              char     *property_name,
                                              char     *property_type,
                                              char     *property_value)

  # --> functions about visualization
  void CWP_client_Visu_set(char                 *local_code_name,
                           char                 *cpl_id,
                           int                   freq,
                           CWP_Visu_format_t     visu_format,
                           char                 *format_option)

  # --> functions about User target points
  void CWP_client_User_tgt_pts_set(char    *local_code_name,
                                   char    *cpl_id,
                                   int      i_part,
                                   int      n_pts,
                                   double   *coord,
                                   long     *global_num)

  # --> functions about Mesh
  void CWP_client_Mesh_interf_finalize(char         *local_code_name,
                                       char         *cpl_id)
  void CWP_client_Mesh_interf_vtx_set(char           *local_code_name,
                                      char           *cpl_id,
                                      int             i_part,
                                      int             n_pts,
                                      double          *coord,
                                      long            *global_num)
  int CWP_client_Mesh_interf_block_add(char           *local_code_name,
                                       char           *cpl_id,
                                       CWP_Block_t     block_type)
  void CWP_client_Mesh_interf_block_std_set(char        *local_code_name,
                                            char        *cpl_id,
                                            int          i_part,
                                            int          block_id,
                                            int          n_elts,
                                            int          *connec,
                                            long         *global_num)
  void CWP_client_Mesh_interf_block_std_get(char         *local_code_name,
                                            char         *cpl_id,
                                            int           i_part,
                                            int           block_id,
                                            int          *n_elts,
                                            int         **connec,
                                            long        **global_num)
  CWP_Block_t CWP_client_std_block_type_get(char             *local_code_name,
                                            char             *cpl_id,
                                            int               block_id)
  void CWP_client_Mesh_interf_f_poly_block_set(char             *local_code_name,
                                               char             *cpl_id,
                                               int               i_part,
                                               int               block_id,
                                               int               n_elts,
                                               int              *connec_idx,
                                               int              *connec,
                                               long             *global_num)
  void CWP_client_Mesh_interf_f_poly_block_get(char             *local_code_name,
                                               char             *cpl_id,
                                               int               i_part,
                                               int               block_id,
                                               int              *n_elts,
                                               int             **connec_idx,
                                               int             **connec,
                                               long            **global_num)
  void CWP_client_Mesh_interf_c_poly_block_set(char           *local_code_name,
                                               char           *cpl_id,
                                               int             i_part,
                                               int             block_id,
                                               int             n_elts,
                                               int             n_faces,
                                               int            *connec_faces_idx,
                                               int            *connec_faces,
                                               int            *connec_cells_idx,
                                               int            *connec_cells,
                                               long           *global_num)
  void CWP_client_Mesh_interf_c_poly_block_get(char           *local_code_name,
                                               char           *cpl_id,
                                               int             i_part,
                                               int             block_id,
                                               int              *n_elts,
                                               int              *n_faces,
                                               int             **connec_faces_idx,
                                               int             **connec_faces,
                                               int             **connec_cells_idx,
                                               int             **connec_cells,
                                               long            **global_num)
  void CWP_client_Mesh_interf_del(char *local_code_name,
                                  char *cpl_id)
  void CWP_client_Mesh_interf_from_cellface_set(char           *local_code_name,
                                                char           *cpl_id,
                                                int             i_part,
                                                int             n_cells,
                                                int            *cell_face_idx,
                                                int            *cell_face,
                                                int             n_faces,
                                                int            *face_vtx_idx,
                                                int            *face_vtx,
                                                long           *global_num)
  void CWP_client_Mesh_interf_from_faceedge_set(char           *local_code_name,
                                                char           *cpl_id,
                                                int             i_part,
                                                int             n_faces,
                                                int            *face_edge_idx,
                                                int            *face_edge,
                                                int             n_edges,
                                                int            *edge_vtx_idx,
                                                int            *edge_vtx,
                                                long           *global_num)

  # --> functions about field
  void CWP_client_Field_create(char                  *local_code_name,
                               char                  *cpl_id,
                               char                  *field_id,
                               CWP_Type_t             data_type,
                               CWP_Field_storage_t    storage,
                               int                    n_component,
                               CWP_Dof_location_t     target_location,
                               CWP_Field_exch_t       exch_type,
                               CWP_Status_t           visu_status)
  void CWP_client_Field_data_set(char              *local_code_name,
                                 char              *cpl_id,
                                 char              *field_id,
                                 int                i_part,
                                 CWP_Field_map_t    map_type,
                                 int                n_entities,
                                 double            *data)
  int CWP_client_Field_n_component_get(char      *local_code_name,
                                       char      *cpl_id,
                                       char      *field_id)
  CWP_Dof_location_t CWP_client_Field_target_dof_location_get(char      *local_code_name,
                                                              char      *cpl_id,
                                                              char      *field_id)
  CWP_Field_storage_t CWP_client_Field_storage_get(char      *local_code_name,
                                                   char      *cpl_id         ,
                                                   char      *field_id)
  void CWP_client_Field_del(char      *local_code_name,
                            char      *cpl_id         ,
                            char      *field_id)

  # --> functions about exchange
  void CWP_client_Field_issend(char     *local_code_name,
                               char     *cpl_id,
                               char     *src_field_id)
  void CWP_client_Field_irecv(char        *local_code_name,
                              char        *cpl_id,
                              char        *tgt_field_id)
  void CWP_client_Field_wait_issend(char  *local_code_name,
                                    char  *cpl_id,
                                    char  *src_field_id)
  void CWP_client_Field_wait_irecv(char  *local_code_name,
                                   char  *cpl_id,
                                   char  *tgt_field_id)

  # --> functions about control parameters
  void CWP_client_Param_add(char        *local_code_name,
                            char        *param_name,
                            CWP_Type_t   data_type,
                            void        *initial_value)
  void CWP_client_Param_set(char             *local_code_name,
                            char             *param_name,
                            CWP_Type_t        data_type,
                            void             *value)
  void CWP_client_Param_del(char       *local_code_name,
                            char       *param_name,
                            CWP_Type_t  data_type)

  # --> functions about all code parameters
  int CWP_client_Param_n_get(char             *code_name,
                             CWP_Type_t        data_type)
  void CWP_client_Param_list_get(char             *code_name,
                                 CWP_Type_t        data_type,
                                 int              *nParam,
                                 char           ***paramNames)
  int CWP_client_Param_is(char             *code_name,
                          char             *param_name,
                          CWP_Type_t        data_type)
  void CWP_client_Param_get(char       *code_name,
                            char       *param_name,
                            CWP_Type_t  data_type,
                            void       *value)
  void CWP_client_Param_reduce(CWP_Op_t    op,
                               char       *param_name,
                               CWP_Type_t  data_type,
                               void       *res,
                               int         nCode,
                               char      **code_names)
  void CWP_client_Param_lock(char *code_name)
  void CWP_client_Param_unlock(char *code_name)
