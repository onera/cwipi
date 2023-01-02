# cython: c_string_type=str, c_string_encoding=ascii
#-----------------------------------------------------------------------------
# This file is part of the CWIPI library.
#
# Copyright (C) 2011  ONERA
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 3 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library. If not, see <http://www.gnu.org/licenses/>.
#-----------------------------------------------------------------------------

"""
cwipi - Coupling With Interpolation Parallel Interface library.
"""

#-----------------------------------------------------------------------------
# IMPORTS

# --> Python
import numpy     as NPY
from mpi4py import MPI

# --> Cython
from libc.stdlib       cimport malloc, free
from cpython.ref       cimport PyObject, Py_INCREF, Py_DECREF, Py_XDECREF
from cpython.pycapsule cimport PyCapsule_New, PyCapsule_GetPointer, PyCapsule_IsValid, PyCapsule_GetName
from libc.stdio cimport FILE, fdopen
from cpython.object cimport PyObject_AsFileDescriptor
cimport numpy as NPY
cimport mpi4py.MPI as MPI

# initialize the numpy C API
NPY.import_array()

#-----------------------------------------------------------------------------
# EXTERN

cdef extern from "cwp.h":

  # TYPES

  ctypedef enum CWP_Type_t:
      CWP_DOUBLE
      CWP_INT
      CWP_CHAR

  ctypedef enum CWP_Visu_format_t:
      CWP_VISU_FORMAT_ENSIGHT

  ctypedef enum CWP_Comm_t:
      CWP_COMM_PAR_WITH_PART
      CWP_COMM_PAR_WITHOUT_PART
      CWP_COMM_SEQ

  ctypedef enum CWP_Time_exch_t:
      CWP_TIME_EXCH_USER_CONTROLLED

  ctypedef enum CWP_Dof_location_t:
      CWP_DOF_LOCATION_UNDEF
      CWP_DOF_LOCATION_CELL_CENTER
      CWP_DOF_LOCATION_NODE
      CWP_DOF_LOCATION_USER

  ctypedef enum CWP_Field_exch_t:
      CWP_FIELD_EXCH_SEND
      CWP_FIELD_EXCH_RECV
      CWP_FIELD_EXCH_SENDRECV

  ctypedef enum CWP_Field_map_t:
      CWP_FIELD_MAP_SOURCE
      CWP_FIELD_MAP_TARGET

  ctypedef enum CWP_Field_storage_t:
      CWP_FIELD_STORAGE_INTERLACED
      CWP_FIELD_STORAGE_INTERLEAVED

  ctypedef enum CWP_Status_t:
      CWP_STATUS_OFF
      CWP_STATUS_ON

  ctypedef enum CWP_Err_t:
      CWP_ERR_NO_ERROR
      CWP_ERR_DEFAULT

  ctypedef enum CWP_Block_t:
      CWP_BLOCK_NODE
      CWP_BLOCK_EDGE2
      CWP_BLOCK_FACE_TRIA3
      CWP_BLOCK_FACE_QUAD4
      CWP_BLOCK_FACE_POLY
      CWP_BLOCK_CELL_TETRA4
      CWP_BLOCK_CELL_HEXA8
      CWP_BLOCK_CELL_PRISM6
      CWP_BLOCK_CELL_PYRAM5
      CWP_BLOCK_CELL_POLY

  ctypedef enum CWP_Dynamic_mesh_t:
      CWP_DYNAMIC_MESH_STATIC
      CWP_DYNAMIC_MESH_DEFORMABLE
      CWP_DYNAMIC_MESH_VARIABLE

  ctypedef enum CWP_Spatial_interp_t:
      CWP_SPATIAL_INTERP_FROM_CLOSEST_POINT_LEAST_SQUARES
      CWP_SPATIAL_INTERP_FROM_INTERSECTION
      CWP_SPATIAL_INTERP_FROM_LOCATION_DIST_CLOUD_SURF
      CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE
      CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_DBBTREE

  ctypedef enum CWP_Interface_t:
      CWP_INTERFACE_POINT
      CWP_INTERFACE_LINEAR
      CWP_INTERFACE_SURFACE
      CWP_INTERFACE_VOLUME

  ctypedef enum CWP_State_t:
      CWP_STATE_IN_PROGRESS
      CWP_STATE_END
      CWP_STATE_OUTPUT_ERROR

  ctypedef enum CWP_Op_t:
      CWP_OP_MIN
      CWP_OP_MAX
      CWP_OP_SUM

  ctypedef void (*CWP_interp_from_location_t) (int                  interface_type,
                                               char                *code_name,
                                               int                  src_n_block,
                                               CWP_Block_t         *src_blocks_type,
                                               int                  src_i_part,
                                               int                  src_n_vtx,
                                               double              *src_vtx_coords,
                                               long                *src_vtx_global_num,
                                               int                  src_n_elts,
                                               int                 *src_id_block,
                                               int                 *src_elt_in_block,
                                               int                 *src_elt_vtx_idx,
                                               int                 *src_elt_vtx,
                                               long                *src_elts_global_num,
                                               int                  tgt_n_pts,
                                               int                 *tgt_pts_elt_idx,
                                               double              *tgt_pts_coords,
                                               double              *tgt_pts_dist,
                                               double              *tgt_pts_uvw,
                                               int                 *tgt_pts_weights_idx,
                                               double              *tgt_pts_weights,
                                               int                  stride,
                                               CWP_Dof_location_t   src_field_dof_location,
                                               void                *src_field,
                                               void                *tgt_field)

  ctypedef void (*CWP_interp_from_intersect_t) (int interface_type)

  ctypedef void (*CWP_interp_from_closest_pts_t) (int interface_type)

  #-----------------------------------------------------------------------------
  # CWIPI
  # --> general functions
  void CWP_Init(MPI.MPI_Comm       global_comm,
                int                n_code,
                char             **code_names,
                CWP_Status_t      *is_active_rank,
                double            *time_init,
                MPI.MPI_Comm      *intra_comms)
  void CWP_Finalize()

  # --> functions about current code properties
  void CWP_State_update(char* local_code_name,
                        CWP_State_t state)
  void CWP_Time_update(char* local_code_name,
                       double current_time)
  void CWP_Output_file_set(FILE *output_file)
  void CWP_User_structure_set(char* local_code_name,
                              void* user_structure)
  void *CWP_User_structure_get(char* local_code_name)

  # --> functions about other code properties
  CWP_State_t CWP_State_get(char    *code_name)
  int CWP_Codes_nb_get()
  char **CWP_Codes_list_get()
  int CWP_Loc_codes_nb_get()
  char **CWP_Loc_codes_list_get()

  # --> functions about properties
  void CWP_Properties_dump()

  # --> general functions about coupling
  void CWP_Cpl_create(char                *local_code_name,
                      char                *cpl_id,
                      char                *coupled_code_name,
                      CWP_Interface_t      entities_dim,
                      CWP_Comm_t           comm_type,
                      CWP_Spatial_interp_t spatial_interp,
                      int                  n_part,
                      CWP_Dynamic_mesh_t   displacement,
                      CWP_Time_exch_t      recv_freq_type)
  void CWP_Cpl_del(char *local_code_name,
                   char *cpl_id)
  int CWP_N_uncomputed_tgts_get(char *local_code_name,
                                char *cpl_id,
                                char *field_id,
                                int   i_part)
  int *CWP_Uncomputed_tgts_get(char *local_code_name,
                               char *cpl_id,
                               char *field_id,
                               int   i_part)
  int CWP_N_computed_tgts_get(char *local_code_name,
                              char *cpl_id,
                              char *field_id,
                              int   i_part)
  int *CWP_Computed_tgts_get(char *local_code_name,
                             char *cpl_id,
                             char *field_id,
                             int   i_part)
  int CWP_N_involved_srcs_get(char *local_code_name,
                              char *cpl_id,
                              char *field_id,
                              int   i_part)
  int *CWP_Involved_srcs_get(char *local_code_name,
                             char *cpl_id,
                             char *field_id,
                             int   i_part)

  # --> functions about spatial interpolation
  void CWP_Spatial_interp_weights_compute(char     *local_code_name,
                                          char     *cpl_id)
  void CWP_Spatial_interp_property_set(char     *local_code_name,
                                       char     *cpl_id,
                                       char     *property_name,
                                       char     *property_type,
                                       char     *property_value)

  # --> functions about visualization
  void CWP_Visu_set(char                 *local_code_name,
                    char                 *cpl_id,
                    int                   freq,
                    CWP_Visu_format_t     visu_format,
                    char                 *format_option)

  # --> functions about User target points
  void CWP_User_tgt_pts_set(char    *local_code_name,
                            char    *cpl_id,
                            int      i_part,
                            int      n_pts,
                            double   *coord,
                            long     *global_num)

  # --> functions about Mesh
  void CWP_Mesh_interf_finalize(char         *local_code_name,
                                char         *cpl_id)
  void CWP_Mesh_interf_vtx_set(char           *local_code_name,
                               char           *cpl_id,
                               int             i_part,
                               int             n_pts,
                               double          *coord,
                               long            *global_num)
  int CWP_Mesh_interf_block_add(char           *local_code_name,
                                char           *cpl_id,
                                CWP_Block_t     block_type)
  void CWP_Mesh_interf_block_std_set(char        *local_code_name,
                                     char        *cpl_id,
                                     int          i_part,
                                     int          block_id,
                                     int          n_elts,
                                     int          *connec,
                                     long         *global_num)
  void CWP_Mesh_interf_block_std_get(char         *local_code_name,
                                     char         *cpl_id,
                                     int           i_part,
                                     int           block_id,
                                     int          *n_elts,
                                     int         **connec,
                                     long        **global_num)
  void CWP_Mesh_interf_f_poly_block_set(char             *local_code_name,
                                        char             *cpl_id,
                                        int               i_part,
                                        int               block_id,
                                        int               n_elts,
                                        int              *connec_idx,
                                        int              *connec,
                                        long             *global_num)
  void CWP_Mesh_interf_f_poly_block_get(char             *local_code_name,
                                        char             *cpl_id,
                                        int               i_part,
                                        int               block_id,
                                        int              *n_elts,
                                        int             **connec_idx,
                                        int             **connec,
                                        long            **global_num)
  void CWP_Mesh_interf_c_poly_block_set(char           *local_code_name,
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
  void CWP_Mesh_interf_c_poly_block_get(char           *local_code_name,
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
  void CWP_Mesh_interf_del(char *local_code_name,
                           char *cpl_id)
  void CWP_Mesh_interf_from_cellface_set(char           *local_code_name,
                                         char           *cpl_id,
                                         int             i_part,
                                         int             n_cells,
                                         int            *cell_face_idx,
                                         int            *cell_face,
                                         int             n_faces,
                                         int            *face_vtx_idx,
                                         int            *face_vtx,
                                         long           *global_num)
  void CWP_Mesh_interf_from_faceedge_set(char           *local_code_name,
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
  void CWP_Field_create(char                  *local_code_name,
                        char                  *cpl_id,
                        char                  *field_id,
                        CWP_Type_t             data_type,
                        CWP_Field_storage_t    storage,
                        int                    n_component,
                        CWP_Dof_location_t     target_location,
                        CWP_Field_exch_t       exch_type,
                        CWP_Status_t           visu_status)
  void CWP_Field_data_set(char              *local_code_name,
                          char              *cpl_id,
                          char              *field_id,
                          int                i_part,
                          CWP_Field_map_t    map_type,
                          double            *data)
  int CWP_Field_n_component_get(char      *local_code_name,
                                char      *cpl_id,
                                char      *field_id)
  CWP_Dof_location_t CWP_Field_target_dof_location_get(char      *local_code_name,
                                                       char      *cpl_id,
                                                       char      *field_id)
  CWP_Field_storage_t CWP_Field_storage_get(char      *local_code_name,
                                            char      *cpl_id         ,
                                            char      *field_id)
  void CWP_Field_del(char      *local_code_name,
                     char      *cpl_id         ,
                     char      *field_id)

  # --> functions about exchange
  void CWP_Field_issend(char     *local_code_name,
                        char     *cpl_id,
                        char     *src_field_id)
  void CWP_Field_irecv(char        *local_code_name,
                       char        *cpl_id,
                       char        *tgt_field_id)
  void CWP_Field_wait_issend(char  *local_code_name,
                             char  *cpl_id,
                             char  *src_field_id)
  void CWP_Field_wait_irecv(char  *local_code_name,
                            char  *cpl_id,
                            char  *tgt_field_id)

  # --> functions about user interpolation
  void CWP_Interp_from_location_unset(char                 *local_code_name,
                                      char                 *cpl_id,
                                      char                 *src_field_id)
  void CWP_Interp_from_location_set(char                       *local_code_name,
                                    char                       *cpl_id,
                                    char                       *src_field_id,
                                    CWP_interp_from_location_t  fct)

  # --> functions about control parameters
  void CWP_Param_add(char        *local_code_name,
                     char        *param_name,
                     CWP_Type_t   data_type,
                     void        *initial_value)
  void CWP_Param_set(char             *local_code_name,
                     char             *param_name,
                     CWP_Type_t        data_type,
                     void             *value)
  void CWP_Param_del(char       *local_code_name,
                     char       *param_name,
                     CWP_Type_t  data_type)

  # --> functions about all code parameters
  int CWP_Param_n_get(char             *code_name,
                      CWP_Type_t        data_type)
  void CWP_Param_list_get(char             *code_name,
                          CWP_Type_t        data_type,
                          int              *nParam,
                          char           ***paramNames)
  int CWP_Param_is(char             *code_name,
                   char             *param_name,
                   CWP_Type_t        data_type)
  void CWP_Param_get(char       *code_name,
                     char       *param_name,
                     CWP_Type_t  data_type,
                     void       *value)
  void CWP_Param_reduce(CWP_Op_t    op,
                        char       *param_name,
                        CWP_Type_t  data_type,
                        void       *res,
                        int         nCode,
                        char      **code_names)
  void CWP_Param_lock(char *code_name)
  void CWP_Param_unlock(char *code_name)

  # --> client-server TO DO

#-----------------------------------------------------------------------------
# CALLBACK
g_interp_fct = {}
current_cpl_id = -1

cdef void interp_callback(const int                  interface_type,
                          const char                *code_name,
                          const int                  src_n_block,
                          const CWP_Block_t          src_blocks_type[],
                          const int                  src_i_part,
                          const int                  src_n_vtx,
                          const double               src_vtx_coords[],
                          const long                 src_vtx_global_num[], # TO DO replace by CWP_g_num_t ?
                          const int                  src_n_elts,
                          const int                  src_id_block[],
                          const int                  src_elt_in_block[],
                          const int                  src_elt_vtx_idx[],
                          const int                  src_elt_vtx[],
                          const long                 src_elts_global_num[], # TO DO replace by CWP_g_num_t ?
                          const int                  tgt_n_pts,
                          const int                  tgt_pts_elt_idx[],
                          const double               tgt_pts_coords[],
                          const double               tgt_pts_dist[],
                          const double               tgt_pts_uvw[],
                          const int                  tgt_pts_weights_idx[],
                          const double               tgt_pts_weights[],
                          const int                  stride,
                          const CWP_Dof_location_t   src_field_dof_location,
                          const void                *src_field,
                          void                      *tgt_field):
    """
    User interpolation function interface from location into a mesh.

    Parameters:
      global_comm (MPI.Comm) : MPI global communicator
      interface_type              Interface type
      code_name                   Name of code
      src_n_block                 Number of blocks
      src_block_type              Block types (size = n_block)
      src_i_part                  Part id
      src_n_vtx                   Number of vertices
      src_vtx_coords              Coordinates of vertices (size = 3 * src_n_vtx)
      src_vtx_global_num          Global number of vertices (size = src_n_vtx)
      src_n_elts                  Number of elements
      src_id_block                 block id of the element
                                  (size = src_n_elts)
      src_elt_in_block            Element number of the elements into it block
                                  (size = src_n_elts)
      src_elt_vtx_idx              Element to vertex index
                                  (src_elt_vtx_idx[0] = 0 and
                                  size = src_n_elts + 1)
      src_elt_vtx                  Element to vertex connectivity.
                                  (size = src_elt_vtx_idx[src_n_elts])
      src_elts_global_num         Global number of elements (size = src_n_elts)
      tgt_n_pts                   Number of target points
      tgt_pts_elt_idx             The list of target points located in each element
                                  (size = src_n_elts + 1)
      tgt_pts_coords              Target points coordinates
                                  (size = 3 * tgt_n_pts)
      tgt_pts_dist                target points distance to location element
                                  (size = tgt_n_pts)
      tgt_pts_uvw                 Parametric coordinates of target points in the elements
                                  ( 0 <= u <= 1, 0 <= v <= 1, -1 : for polydra and polygons)
                                  (size = dim_interface * tgt_n_pts)
      tgt_pts_weights_idx         Index of Barycentric coordinates target points
                                  in location element
                                  (tgt_pts_bary_coords_idx[0] = 0 and
                                  size = n_tgt_pts + 1)
      tgt_pts_weights             Barycentric coordinates target points
                                  in location element
                                  (size = tgt_pts_weights_idx[n_tgt_pts])
      stride                      Number of field components
      src_field_dof_location      source field location
      src_field                   source field
                                  (size depends on field type and stride)

    Returns:
      tgt_field                   target field
                                  (size = stride * n_tgt_pts)
    """

    global g_interp_fct
    global current_cpl_id
    cdef NPY.npy_intp dims = 0

    if (src_n_block == 0 or src_blocks_type == NULL):
      src_blocks_type_a = None
    else:
      dims = <NPY.npy_intp>(src_n_block)
      src_blocks_type_a = NPY.PyArray_SimpleNewFromData(1,
                                                        &dims,
                                                        NPY.NPY_LONG, # TO DO: ok for CWP_Block_t?
                                                        <void *> src_blocks_type)

    if (src_n_vtx == 0 or src_vtx_coords == NULL):
      src_vtx_coords_a = None
    else:
      dims = <NPY.npy_intp>(3 * src_n_vtx)
      src_vtx_coords_a = NPY.PyArray_SimpleNewFromData(1,
                                                       &dims,
                                                       NPY.NPY_DOUBLE,
                                                       <void *> src_vtx_coords)

    if (src_n_vtx == 0 or src_vtx_global_num == NULL):
      src_vtx_global_num_a = None
    else:
      dims = <NPY.npy_intp>(src_n_vtx)
      src_vtx_global_num_a = NPY.PyArray_SimpleNewFromData(1,
                                                           &dims,
                                                           NPY.NPY_LONG,
                                                           <void *> src_vtx_global_num)

    if (src_n_elts == 0 or src_id_block == NULL):
      src_id_block_a = None
    else:
      dims = <NPY.npy_intp>(src_n_elts)
      src_id_block_a = NPY.PyArray_SimpleNewFromData(1,
                                                     &dims,
                                                     NPY.NPY_INT,
                                                     <void *> src_id_block)

    if (src_n_elts == 0 or src_elt_in_block == NULL):
      src_elt_in_block_a = None
    else:
      dims = <NPY.npy_intp>(src_n_elts)
      src_elt_in_block_a = NPY.PyArray_SimpleNewFromData(1,
                                                     &dims,
                                                     NPY.NPY_INT,
                                                     <void *> src_elt_in_block)

    if (src_n_elts == 0 or src_elt_vtx_idx == NULL):
      src_elt_vtx_idx_a = None
      src_elt_vtx_a     = None
    else:
      dims = <NPY.npy_intp>(src_n_elts+1)
      src_elt_vtx_idx_a = NPY.PyArray_SimpleNewFromData(1,
                                                        &dims,
                                                        NPY.NPY_INT,
                                                        <void *> src_elt_vtx_idx)
      if (src_elt_vtx_idx[src_n_elts] == 0 or src_elt_vtx == NULL):
        src_elt_vtx_a = None
      else:
        dims = <NPY.npy_intp>(src_elt_vtx[src_n_elts])
        src_elt_vtx_a = NPY.PyArray_SimpleNewFromData(1,
                                                      &dims,
                                                      NPY.NPY_INT,
                                                      <void *> src_elt_vtx)

    if (src_n_elts == 0 or src_elts_global_num == NULL):
      src_elts_global_num_a = None
    else:
      dims = <NPY.npy_intp>(src_n_elts)
      src_elts_global_num_a = NPY.PyArray_SimpleNewFromData(1,
                                                            &dims,
                                                            NPY.NPY_LONG,
                                                            <void *> src_elts_global_num)

    if (tgt_n_pts == 0 or tgt_pts_elt_idx == NULL):
      tgt_pts_elt_idx_a = None
    else:
      dims = <NPY.npy_intp>(tgt_n_pts+1)
      tgt_pts_elt_idx_a = NPY.PyArray_SimpleNewFromData(1,
                                                        &dims,
                                                        NPY.NPY_INT,
                                                        <void *> tgt_pts_elt_idx)

    if (tgt_n_pts == 0 or tgt_pts_coords == NULL):
      tgt_pts_coords_a = None
    else:
      dims = <NPY.npy_intp>(3*tgt_n_pts)
      tgt_pts_coords_a = NPY.PyArray_SimpleNewFromData(1,
                                                       &dims,
                                                       NPY.NPY_DOUBLE,
                                                       <void *> tgt_pts_coords)

    if (tgt_n_pts == 0 or tgt_pts_dist == NULL):
      tgt_pts_dist_a = None
    else:
      dims = <NPY.npy_intp>(tgt_n_pts)
      tgt_pts_dist_a = NPY.PyArray_SimpleNewFromData(1,
                                                     &dims,
                                                     NPY.NPY_DOUBLE,
                                                     <void *> tgt_pts_dist)

    if (tgt_n_pts == 0 or tgt_pts_uvw == NULL):
      tgt_pts_uvw_a = None
    else:
      dim_interface = 3 # tmp
      dims = <NPY.npy_intp>(dim_interface * tgt_n_pts) # TO DO ?
      tgt_pts_uvw_a = NPY.PyArray_SimpleNewFromData(1,
                                                    &dims,
                                                    NPY.NPY_DOUBLE,
                                                    <void *> tgt_pts_uvw)

    if (tgt_n_pts == 0 or tgt_pts_weights_idx == NULL):
      tgt_pts_weights_idx_a = None
      tgt_pts_weights_a     = None
    else:
      dims = <NPY.npy_intp>(tgt_n_pts+1)
      tgt_pts_weights_idx_a = NPY.PyArray_SimpleNewFromData(1,
                                                    &dims,
                                                    NPY.NPY_INT,
                                                    <void *> tgt_pts_weights_idx)
      if (tgt_pts_weights_idx[tgt_n_pts] == 0 or tgt_pts_weights == NULL):
        tgt_pts_weights_a = None
      else:
        dims = <NPY.npy_intp>(tgt_pts_weights[tgt_n_pts])
        tgt_pts_weights_a = NPY.PyArray_SimpleNewFromData(1,
                                                         &dims,
                                                         NPY.NPY_DOUBLE,
                                                         <void *> tgt_pts_weights)

    if (src_n_elts == 0 or src_field == NULL):
      src_field_a = None
    else:
      size = stride*tgt_n_pts # tmp
      dims = <NPY.npy_intp>(size) # TO DO ?
      src_field_a = NPY.PyArray_SimpleNewFromData(1,
                                                  &dims,
                                                  NPY.NPY_DOUBLE,
                                                  <void *> src_field)

    if (tgt_n_pts == 0 or tgt_field == NULL):
      tgt_field_a = None
    else:
      dims = <NPY.npy_intp>(stride*tgt_n_pts)
      tgt_field_a = NPY.PyArray_SimpleNewFromData(1,
                                                  &dims,
                                                  NPY.NPY_DOUBLE,
                                                  <void *> tgt_field)

    (<object> g_interp_fct[current_cpl_id]) (interface_type,
                                                            code_name,
                                                            src_n_block,
                                                            src_blocks_type_a,
                                                            src_i_part,
                                                            src_n_vtx,
                                                            src_vtx_coords_a,
                                                            src_vtx_global_num_a,
                                                            src_n_elts,
                                                            src_id_block_a,
                                                            src_elt_in_block_a,
                                                            src_elt_vtx_idx_a,
                                                            src_elt_vtx_a,
                                                            src_elts_global_num_a,
                                                            tgt_n_pts,
                                                            tgt_pts_elt_idx_a,
                                                            tgt_pts_coords_a,
                                                            tgt_pts_dist_a,
                                                            tgt_pts_uvw_a,
                                                            tgt_pts_weights_idx_a,
                                                            tgt_pts_weights_a,
                                                            stride,
                                                            src_field_dof_location,
                                                            src_field_a,
                                                            tgt_field_a)

#-----------------------------------------------------------------------------
# UTILS

# > Enable transfer of ownership between C and Python
cdef extern from "numpy/arrayobject.h":
    void PyArray_ENABLEFLAGS(NPY.ndarray arr, int flags)

cdef create_numpy_i(int* array, int size, bint flag_owndata=True):
    dim = <NPY.npy_intp> size
    nparray = NPY.PyArray_SimpleNewFromData(1, &dim, NPY.NPY_INT32, <void *> array)
    if flag_owndata:
        PyArray_ENABLEFLAGS(nparray, NPY.NPY_OWNDATA)
    return nparray

cdef create_numpy_l(long int* array, int size, bint flag_owndata=True):
    dim = <NPY.npy_intp> size
    nparray = NPY.PyArray_SimpleNewFromData(1, &dim, NPY.NPY_INT64, <void *> array)
    if flag_owndata:
        PyArray_ENABLEFLAGS(nparray, NPY.NPY_OWNDATA)
    return nparray

#-----------------------------------------------------------------------------
# DEFINITION
# --> general functions

def init(MPI.Comm global_comm,
         int      n_code,
         NPY.ndarray[NPY.str] code_names):
  """
  Initialize CWIPI library and create intra-communicator.

  Parameters:
    global_comm (MPI.Comm) : MPI global communicator
    n_code      (int)      : Number of codes on the current rank
    code_names  (char**)   : Names of codes on the current rank (size = n_code)

  Returns:
    is_active_rank (int)      : Is current rank have to be used by CWIPI (size = n_code)
    time_init      (double)   : Initial time (size = n_code)
    intra_comms    (MPI.Comm) : MPI intra communicators of each code (size = n_code)
  """

  cdef char **c_code_names = <char**> malloc(sizeof(char*) * n_code)
  for i in range(n_code):
    c_code_names[i] = <char *> malloc(sizeof(char) * len(code_names[i]))

  for i in range(n_code):
    py_byte_code_name = code_names[i].encode('UTF-8')
    c_code_names[i] = py_byte_code_name

  cdef CWP_Status_t is_active_rank
  cdef double       time_init

  cdef MPI.MPI_Comm c_global_comm = global_comm.ob_mpi

  cdef MPI.MPI_Comm* intra_comms = <MPI.MPI_Comm *> malloc(n_code * sizeof(MPI.MPI_Comm))
  cdef py_intra_comms = [MPI.Comm()]*n_code

  cdef MPI.MPI_Comm c_comm
  cdef MPI.Comm comm

  for i_code in range(n_code):
    comm = py_intra_comms[i_code]
    c_comm = <MPI.MPI_Comm> comm.ob_mpi
    intra_comms[i_code] = c_comm

  CWP_Init(c_global_comm,
           n_code,
           c_code_names,
           &(is_active_rank),
           &(time_init),
           intra_comms)

  free(intra_comms)

  return {
          'is_active_rank' : is_active_rank,
          'time_init'      : time_init,
          'intra_comms'    : py_intra_comms
         }

def finalize():
  """
  Finalize this module.

  After this call, no other function from this module can be called.
  It is a synchronization point between all applications.
  """

  CWP_Finalize()

# --> functions about current code properties
def state_update(char*       local_code_name,
                 CWP_State_t state):
    """
    Update code state.

    Parameters:
      local_code_name (char*)       : Local code name
      state           (CWP_State_t) : State
    """

    CWP_State_update(local_code_name,
                     state)

def time_update(char*  local_code_name,
                double current_time):
  """
  Update code time.

  Parameters:
    local_code_name (char*)  : Local code name
    current_time    (double) : Current time
  """

  CWP_Time_update(local_code_name,
                  current_time)

def output_file_set(output_file):
  """
  Define output file.

  Parameters:
    output_file(file):  Output file

  Note:
    If this file is written by both cwipi and python, I/O buffering may mix up the output.
    In this case, calling output_listing.flush() before any cwipi call may help.
  """

  cdef int fd
  cdef FILE* c_file
  cdef str mode

  # next line will most likely fail if output_file is not
  # an io.TextIOWrapper(in py3) or a file(in py2)
  mode = output_file.mode

  fd = PyObject_AsFileDescriptor(output_file)
  c_file = fdopen(fd, mode)

  CWP_Output_file_set(c_file)

def user_structure_set(char*     local_code_name,
                       object user_structure): # TO DO: object ?
  """
  Define a user structure associated to a code.

  Parameters:
    local_code_name  (char*)      : Local code name
    user_structure   (PyObject *) : User structure
  """

  cdef void *c_user_structure = <void *> PyCapsule_GetPointer(user_structure, NULL)

  CWP_User_structure_set(local_code_name,
                         c_user_structure)

def user_structure_get(local_code_name):
  """
  Return the user structure associated.

  Parameters:
    local_code_name  (char*) : Local code name

  Returns:
    user_structure   (PyObject *) : User structure
  """

  cdef void *c_user_structure = CWP_User_structure_get(local_code_name)

  if c_user_structure: # not Null
    user_structure = PyCapsule_New(c_user_structure, NULL, NULL)

  return user_structure

# --> functions about other code properties

def state_get(char* code_name):
  """
  Get code state.

  Parameters:
    local_code_name (char*)       : Local code name
  """

  CWP_State_get(code_name)

def codes_nb_get():
  """
  Return the number of codes known by CWIPI.

  Returns:
    n_code (int) : Number of codes known by CWIPI
  """

  cdef int n_code = CWP_Codes_nb_get()

  return n_code

def codes_list_get():
  """
  Return the list of code names known by CWIPI.

  Returns:
    code_list (char**) : Number of codes known by CWIPI
  """

  cdef char** c_code_list = CWP_Codes_list_get()

  cdef int n_code = CWP_Codes_nb_get()
  code_list = NPY.array(size=n_code, dtype=str)
  for i in range(n_code):
    code_list[i] = c_code_list[i]

  return code_list

def loc_codes_nb_get():
  """
  Return the number of local codes known by CWIPI.

  Returns:
    n_loc_code (int) : Number of codes known by CWIPI
  """

  cdef int n_loc_code = CWP_Loc_codes_nb_get()

  return n_loc_code

def loc_codes_list_get():
  """
  Return the list of code names known by CWIPI.

  Returns:
    loc_code_list (char**) : Number of local codes known by CWIPI
  """

  cdef char** c_loc_code_list = CWP_Loc_codes_list_get()

  cdef int n_loc_code = CWP_Loc_codes_nb_get()
  loc_code_list = NPY.array(size=n_loc_code, dtype=str)
  for i in range(n_loc_code):
    loc_code_list[i] = c_loc_code_list[i]

  return loc_code_list

# --> functions about properties
def properties_dump():
  """
  Dump application properties.
  """

  CWP_Properties_dump()

# --> functions about control parameters
def param_add_int(char        *local_code_name,
                  char        *param_name,
                  int          initial_value):
  """
  Add a new interger parameter and intialize it.

  Parameters:
    local_code_name  (char*) : Local code name
    param_name       (char*) : Parameter name
    initial_value    (int)   : Initial value
  """

  CWP_Param_add(local_code_name,
                param_name,
                CWP_INT,
                &(initial_value))

def param_add_dbl(char        *local_code_name,
                  char        *param_name,
                  int          initial_value):
  """
  Add a new double parameter and intialize it.

  Parameters:
    local_code_name  (char*)  : Local code name
    param_name       (char*)  : Parameter name
    initial_value    (double) : Initial value
  """

  CWP_Param_add(local_code_name,
                param_name,
                CWP_DOUBLE,
                &(initial_value))

def param_add_str(char        *local_code_name,
                  char        *param_name,
                  char        *initial_value):
  """
  Add a new string parameter and intialize it.

  Parameters:
    local_code_name  (char*) : Local code name
    param_name       (char*) : Parameter name
    initial_value    (char*) : Initial value
  """

  CWP_Param_add(local_code_name,
                param_name,
                CWP_CHAR,
                &(initial_value))

def param_set_int(char        *local_code_name,
                  char        *param_name,
                  int          initial_value):
  """
  Set an integer parameter.

  Parameters:
    local_code_name  (char*) : Local code name
    param_name       (char*) : Parameter name
    initial_value    (int)   : Initial value
  """

  CWP_Param_set(local_code_name,
                param_name,
                CWP_INT,
                &(initial_value))

def param_set_dbl(char        *local_code_name,
                  char        *param_name,
                  double       initial_value):
  """
  Set a double parameter.

  Parameters:
    local_code_name  (char*)  : Local code name
    param_name       (char*)  : Parameter name
    initial_value    (double) : Initial value
  """

  CWP_Param_set(local_code_name,
                param_name,
                CWP_DOUBLE,
                &(initial_value))

def param_set_str(char        *local_code_name,
                  char        *param_name,
                  char        *initial_value):
  """
  Set a string parameter.

  Parameters:
    local_code_name  (char*) : Local code name
    param_name       (char*) : Parameter name
    initial_value    (char*) : Initial value
  """

  CWP_Param_set(local_code_name,
                param_name,
                CWP_CHAR,
                &(initial_value))

def param_del(char        *local_code_name,
              char        *param_name,
              CWP_Type_t   data_type):
  """
  Delete a parameter.

  Parameters:
    local_code_name  (char*)      : Local code name
    param_name       (char*)      : Parameter name
    data_type        (CWP_Type_t) : Parameter type
  """

  CWP_Param_del(local_code_name,
                param_name,
                data_type)

# -> functions about all code parameters
def param_n_get(char        *code_name,
                CWP_Type_t   data_type):
  """
  Return the number of parameters for the code code_name.

  Parameters:
    code_name        (char*)      : Code name
    data_type        (CWP_Type_t) : Parameter type

  Returns:
    n_param  (int) : Number of parameters for the code code_name
  """

  cdef int n_param = CWP_Param_n_get(code_name,
                                     data_type)

  return n_param

def param_list_get(char        *code_name,
                   CWP_Type_t   data_type):
  """
  Return the number of parameters for the code code_name.

  Parameters:
    code_name        (char*)      : Code name
    data_type        (CWP_Type_t) : Parameter type

  Returns:
    n_param     (int)    : Number of parameters for the code code_name
    param_names (char**) : Parameter names
  """

  cdef int n_param = -1;
  cdef char **c_param_names = NULL

  CWP_Param_list_get(code_name,
                     data_type,
                     &(n_param),
                     &(c_param_names))

  param_names = NPY.array(size=n_param, dtype=str)
  for i in range(n_param):
    param_names[i] = c_param_names[i]

  return {
          'n_param'     : n_param,
          'param_names' : param_names
         }

def param_is(char        *code_name,
             char        *param_name,
             CWP_Type_t   data_type):
  """
  Is this code_name a parameter ?

  Parameters:
    code_name        (char*)      : Code name
    param_name       (char*)      : Parameter name
    data_type        (CWP_Type_t) : Parameter type

  Returns:
   is_param (int) : 1 : true / 0 : false
  """

  cdef int is_param = CWP_Param_is(code_name,
                              param_name,
                              data_type)

  return is_param

def param_get(char        *code_name,
              char        *param_name,
              CWP_Type_t   data_type):
  """
  Return the parameter value of param_name on code_name.

  Parameters:
    code_name        (char*)      : Code name
    param_name       (char*)      : Parameter name
    data_type        (CWP_Type_t) : Parameter type

  Returns:
   value (data_type) : Parameter value
  """

  cdef int int_value = -1
  cdef double dbl_value = -1
  cdef char *str_value = NULL

  if (data_type == CWP_INT):
    CWP_Param_get(code_name,
                  param_name,
                  data_type,
                  &(int_value))
    return int_value

  elif (data_type == CWP_DOUBLE):
    CWP_Param_get(code_name,
                  param_name,
                  data_type,
                  &(dbl_value))
    return dbl_value

  elif (data_type == CWP_CHAR):
    CWP_Param_get(code_name,
                  param_name,
                  data_type,
                  &(str_value))
    return str_value

def param_reduce(CWP_Op_t    op,
                 char       *param_name,
                 CWP_Type_t  data_type,
                 int         nCode,
                 NPY.ndarray[NPY.str] code_names):
  """
  Return the result of a reduce operation about a parameter.
  The parameter name has to be the same for all codes.

  Parameters:
    op         (CWP_Op_t)   :Operation
    param_name (char*)      :Parameter name
    data_type  (CWP_Type_t) :Parameter type
    nCode      (int)        :Number of codes
    code_names (char**)     :Codes name

  Returns:
    res (data_type) : Result
  """

  cdef char **c_code_names = <char**> malloc(sizeof(char*) * nCode)
  for i in range(nCode):
    c_code_names[i] = <char *> malloc(sizeof(char) * len(code_names[i]))

  for i in range(nCode):
    py_byte_code_name = code_names[i].encode('UTF-8')
    c_code_names[i] = py_byte_code_name

  cdef int int_res = -1
  cdef double dbl_res = -1
  cdef char *str_res = NULL

  if (data_type == CWP_INT):
    CWP_Param_reduce(op,
                     param_name,
                     data_type,
                     &(int_res),
                     nCode,
                     c_code_names)
    return int_res

  elif (data_type == CWP_DOUBLE):
    CWP_Param_reduce(op,
                     param_name,
                     data_type,
                     &(dbl_res),
                     nCode,
                     c_code_names)
    return dbl_res

  elif (data_type == CWP_CHAR):
    CWP_Param_reduce(op,
                     param_name,
                     data_type,
                     &(str_res),
                     nCode,
                     c_code_names)
    return str_res

def param_lock(char        *code_name):
  """
  Lock access to local parameters from a distant code.

  Parameters:
    code_name (char*) : Code name
  """

  CWP_Param_lock(code_name)

def param_unlock(char        *code_name):
  """
  Unlock access to local parameters from a distant code.

  Parameters:
    code_name (char*) : Code name
  """

  CWP_Param_unlock(code_name)

# --> functions about exchange
def field_issend(char     *local_code_name,
                 char     *cpl_id,
                 char     *src_field_id):
  """
  Send a spatially interpolated field to the coupled code with
  nonblocking communications.

  This function is independant of CWP_Time_exch_t mode. The user has to
  manually check the consistency of the exchanges.

  Parameters:
    local_code_name (char*) : Local code name
    cpl_id          (char*) : Coupling identifier
    src_field_id    (char*) : Source field identifier
  """

  CWP_Field_issend(local_code_name,
                   cpl_id,
                   src_field_id)

def field_irecv(char     *local_code_name,
                char     *cpl_id,
                char     *tgt_field_id):
  """
  Wait the end of an exchange related to request from CWP_Field_issend.

  Parameters:
    local_code_name (char*) : Local code name
    cpl_id          (char*) : Coupling identifier
    tgt_field_id    (char*) : Target field identifier
  """

  CWP_Field_irecv(local_code_name,
                  cpl_id,
                  tgt_field_id)

def field_wait_issend(char     *local_code_name,
                      char     *cpl_id,
                      char     *src_field_id):
  """
  Wait the end of an exchange related to request from CWP_Field_issend.

  Parameters:
    local_code_name (char*) : Local code name
    cpl_id          (char*) : Coupling identifier
    src_field_id    (char*) : Source field identifier
  """

  CWP_Field_wait_issend(local_code_name,
                        cpl_id,
                        src_field_id)

def field_wait_irecv(char     *local_code_name,
                     char     *cpl_id,
                     char     *tgt_field_id):
  """
  Wait the end of an exchange related to request from CWP_Field_irecv.

  This function waits the end of exchange related to request
  from CWP_Field_irecv

  Parameters:
    local_code_name (char*) : Local code name
    cpl_id          (char*) : Coupling identifier
    tgt_field_id    (char*) : Target field identifier
  """

  CWP_Field_wait_irecv(local_code_name,
                       cpl_id,
                       tgt_field_id)

#-----------------------------------------------------------------------------
# FIELD CLASS
cdef class Field (object):

    """
    Create a CWIPI field object.
    """

    def __init__(self,
                 char                  *local_code_name,
                 char                  *cpl_id,
                 char                  *field_id,
                 CWP_Type_t             data_type,
                 CWP_Field_storage_t    storage,
                 int                    n_component,
                 CWP_Dof_location_t     target_location,
                 CWP_Field_exch_t       exch_type,
                 CWP_Status_t           visu_status):
      """
      Initialize field object.

      Parameters:
        local_code_name (char*)               : Local code name
        cpl_id          (char*)               : Field id
        field_id        (char*)               : Field id
        data_type       (CWP_Type_t)          : Data type
        storage         (CWP_Field_storage_t) : Storage type
        n_component     (int)                 : Number of component
        target_location (CWP_Dof_location_t)  : Target location
        exch_type       (CWP_Field_exch_t)    : Exchange type
        visu_status     (CWP_Status_t)        : Visualization status
      """

      self.local_code_name = local_code_name
      self.cpl_id = cpl_id
      self.id = field_id

      CWP_Field_create(self.local_code_name,
                       self.cpl_id,
                       self.id,
                       data_type,
                       storage,
                       n_component,
                       target_location,
                       exch_type,
                       visu_status)

    def __del__(self):
      """
      Delete a field object.
      """

      CWP_Field_del(self.local_code_name,
                    self.cpl_id,
                    self.id)

    # --> DATA
    def data_set(self,
                 int                i_part,
                 CWP_Field_map_t    map_type,
                 NPY.ndarray[NPY.double_t, mode='c', ndim=1] data):
      """
      Set field data.

      Parameters:
        i_part    (int)                       : Current partition
        map_type  (CWP_Field_map_t)           : Choice if data is setted for the source or the target
        data      (NPY.ndarray[NPY.double_t]) : Storage array (Mapping)
      """

      CWP_Field_data_set(self.local_code_name,
                         self.cpl_id,
                         self.id,
                         i_part,
                         map_type,
              <double *> data.data)

    def n_component_get(self):
      """
      Get number of field components.

      Returns:
        n_comp (int) : Number of field components
      """

      cdef int n_comp = CWP_Field_n_component_get(self.local_code_name,
                                                  self.cpl_id,
                                                  self.id)

      return n_comp

    def target_dof_location_get(self):
      """
      Get target degrees of freedom location.

      Returns:
        dof_loc (CWP_Dof_location_t) : Location of degrees of freedom
      """

      cdef CWP_Dof_location_t dof_loc = CWP_Field_target_dof_location_get(self.local_code_name,
                                                                          self.cpl_id,
                                                                          self.id)

      return dof_loc

    def storage_get(self):
      """
      Get field storage type.

      Returns:
        storage (CWP_Field_storage_t) : Field storage type
      """

      cdef CWP_Field_storage_t storage = CWP_Field_storage_get(self.local_code_name,
                                                               self.cpl_id,
                                                               self.id)

      return storage

#-----------------------------------------------------------------------------
# COUPLING CLASS
cdef class Cpl (object):

    """
    Create a CWIPI coupling object.
    """

    def __init__(self,
                 char                *local_code_name,
                 char                *cpl_id,
                 char                *coupled_code_name,
                 CWP_Interface_t      entities_dim,
                 CWP_Comm_t           comm_type,
                 CWP_Spatial_interp_t spatial_interp,
                 int                  n_part,
                 CWP_Dynamic_mesh_t   displacement,
                 CWP_Time_exch_t      recv_freq_type):
      """
      Initialize coupling object.

      Parameters:
        local_code_name   (char*)                : Local code name
        cpl_id            (char*)                : Coupling identifier
        coupled_code_name (char*)                : Distant or local coupled code name
        entities_dim      (CWP_Interface_t)      : ?
        comm_type         (CWP_Comm_t)           : Communication type
        spatial_interp    (CWP_Spatial_interp_t) : Spatial interpolation method
        n_part            (int)                  : Number of interface partition
        displacement      (CWP_Dynamic_mesh_t)   : Mesh moving status
        recv_freq_type    (CWP_Time_exch_t)      : Type of receiving frequency
      """

      self.local_code_name = local_code_name
      self.id = cpl_id
      self.fields = {}

      CWP_Cpl_create(local_code_name,
                     cpl_id,
                     coupled_code_name,
                     entities_dim,
                     comm_type,
                     spatial_interp,
                     n_part,
                     displacement,
                     recv_freq_type)

    def __del__(self):
      """
      Delete a coupling object.
      """

      CWP_Cpl_del(self.local_code_name, self.id)

    # --> FIELD
    def field_create(self,
                     char                  *field_id,
                     CWP_Type_t             data_type,
                     CWP_Field_storage_t    storage,
                     int                    n_component,
                     CWP_Dof_location_t     target_location,
                     CWP_Field_exch_t       exch_type,
                     CWP_Status_t           visu_status):
      """
      Create field object and add in self.fields dictionnary.

      Parameters:
        field_id        (char*)               : Field id
        data_type       (CWP_Type_t)          : Data type
        storage         (CWP_Field_storage_t) : Storage type
        n_component     (int)                 : Number of component
        target_location (CWP_Dof_location_t)  : Target location
        exch_type       (CWP_Field_exch_t)    : Exchange type
        visu_status     (CWP_Status_t)        : Visualization status
      """
      self.fields[field_id] = Field(self.local_code_name,
                                    self.id,
                                    field_id,
                                    data_type,
                                    storage,
                                    n_component,
                                    target_location,
                                    exch_type,
                                    visu_status)

    def field_set(self,
                  char                  *field_id,
                  int                i_part,
                  CWP_Field_map_t    map_type,
                  NPY.ndarray[NPY.double_t, mode='c', ndim=1] data):
      """
      Set field data for a given Field object in the dictionnary.

      Parameters:
        field_id        (char*)               : Field id
        i_part    (int)                       : Current partition
        data_type (CWP_Field_map_t)           : Choice if data is setted for the source or the target
        data      (NPY.ndarray[NPY.double_t]) : Storage array (Mapping)
      """
      self.fields[field_id].data_set(i_part,
                                     map_type,
                                     data)

    def field_get(self,
                  char                  *field_id):
      """
      Get field properties.

      Parameters:
        field_id        (char*)               : Field id

      Returns:
        n_comp (int) : Number of field components
        dof_loc (CWP_Dof_location_t) : Location of degrees of freedom
        storage (CWP_Field_storage_t) : Field storage type
      """

      cdef int n_comp = self.fields[field_id].n_component_get()
      cdef CWP_Dof_location_t dof_loc = self.fields[field_id].target_dof_location_get()
      cdef CWP_Field_storage_t storage = self.fields[field_id].storage_get()

      return {
              'n_comp' : n_comp,
              'dof_loc': dof_loc,
              'storage': storage
             }

    def field_del(self,
                  char                  *field_id):
      """
      Delete Field object and field form self.fiels dictionnary

      Parameters:
        field_id        (char*)               : Field id
      """
      del self.fields[field_id]
      # same as:
      # field = self.fields.pop(field_id)
      # del(field)

    # --> MESH
    def mesh_interf_finalize(self):
      """
      Finalize interface mesh.
      This function computes the global numbers of mesh entities if they are
      not provided.
      """

      CWP_Mesh_interf_finalize(self.local_code_name, self.id)

    def mesh_interf_vtx_set(self,
                            int i_part,
                            int n_pts,
                            NPY.ndarray[NPY.double_t, mode='c', ndim=1] coord not None,
                            NPY.ndarray[NPY.long, mode='c', ndim=1] global_num):
      """
      Set vertices.

      Parameters:
        i_part     (int)                       : Current partition
        n_pts      (int)                       : Number of points
        coord      (NPY.ndarray[NPY.double_t]) : Coordinates (size = 3 * n_pts)
        global_num (NPY.ndarray[NPY.long])     : Pointer to parent element number (or None)
      """

      assert 3 * n_pts <= coord.size
      assert n_pts <= global_num.size

      CWP_Mesh_interf_vtx_set(self.local_code_name,
                              self.id,
                              i_part,
                              n_pts,
                   <double *> coord.data,
                     <long *> global_num.data)

    def mesh_interf_block_add(self,
                              CWP_Block_t block_type):
      """
      Add a connectivity block to the interface mesh.

      Parameters:
        block_type (CWP_Block_t) : Block type

      Returns:
        block_id (int) : Block identifier
      """

      cdef int block_id = CWP_Mesh_interf_block_add(self.local_code_name,
                                                    self.id,
                                                    block_type)

      return block_id

    def mesh_interf_block_std_set(self,
                                  i_part,
                                  block_id,
                                  n_elts,
                                  NPY.ndarray[NPY.int32_t, mode='c', ndim=1] connec not None,
                                  NPY.ndarray[NPY.long, mode='c', ndim=1]    global_num):
      """
      Set a standard block to the interface mesh.
      This function adds a connectivity block to the interface mesh.

      Parameters:
        i_part     (int)                      : Partition identifier
        block_id   (int)                      : Block identifier
        n_elts     (int)                      : Number of elements
        connec     (NPY.ndarray[NPY.int32_t]) : Connectivity (size = n_vertex_elt * n_elts)
        global_num (NPY.ndarray[NPY.long])    : Pointer to parent element number (or None)
      """

      CWP_Mesh_interf_block_std_set(self.local_code_name,
                                    self.id,
                                    i_part,
                                    block_id,
                                    n_elts,
                            <int *> connec.data,
                           <long *> global_num.data)

    def mesh_interf_block_std_get(self,
                                  i_part,
                                  block_id):
      """
      Get the properties of a standard block of the interface mesh.

      Parameters:
        i_part     (int) : Partition identifier
        block_id   (int) : Block identifier

      Returns:
        n_elts     (int)                      : Number of elements
        connec     (NPY.ndarray[NPY.int32_t]) : Connectivity (size = n_vertex_elt * n_elts)
        global_num (NPY.ndarray[NPY.long])    : Pointer to parent element number (or None)
      """

      cdef int n_elts = -1
      cdef int  *connec
      cdef long *global_num

      CWP_Mesh_interf_block_std_get(self.local_code_name,
                                    self.id,
                                    i_part,
                                    block_id,
                                    &(n_elts),
                                    &(connec),
                                    &(global_num))

      cdef int size = n_elts * 3 # TO DO: tmp value
      return {
              'n_elts'     : n_elts,
              'connec'     : create_numpy_i(connec, size), # TO DO: get n_vtx_elt
              'global_num' : create_numpy_l(global_num, n_elts)
             }

    def mesh_interf_f_poly_block_set(self,
                                     i_part,
                                     block_id,
                                     n_elts,
                                     NPY.ndarray[NPY.int32_t, mode='c', ndim=1] connec_idx not None,
                                     NPY.ndarray[NPY.int32_t, mode='c', ndim=1] connec not None,
                                     NPY.ndarray[NPY.long, mode='c', ndim=1]    global_num):
      """
      Set the connectivity of a polygon block in a interface mesh partition.

      Parameters:
        i_part     (int)                      : Partition identifier
        block_id   (int)                      : Block identifier
        n_elts     (int)                      : Number of elements
        connec_idx (NPY.ndarray[NPY.int32_t]) : Connectivity index ( connec_id[0] = 0 and size = n_elts + 1)
        connec     (NPY.ndarray[NPY.int32_t]) : Connectivity (size = connec_idx[n_elts])
        global_num (NPY.ndarray[NPY.long])    : Pointer to global element number (or None)
      """

      CWP_Mesh_interf_f_poly_block_set(self.local_code_name,
                                       self.id,
                                       i_part,
                                       block_id,
                                       n_elts,
                               <int *> connec_idx.data,
                               <int *> connec.data,
                              <long *> global_num.data)

    def mesh_interf_f_poly_block_set(self,
                                     i_part,
                                     block_id):
      """
      Get the connectivity of a polygon block in a interface mesh partition.

      Parameters:
        i_part     (int)                      : Partition identifier
        block_id   (int)                      : Block identifier

      Returns:
        n_elts     (int)                      : Number of elements
        connec_idx (NPY.ndarray[NPY.int32_t]) : Connectivity index ( connec_id[0] = 0 and size = n_elts + 1)
        connec     (NPY.ndarray[NPY.int32_t]) : Connectivity (size = connec_idx[n_elts])
        global_num (NPY.ndarray[NPY.long])    : Pointer to global element number (or None)
      """

      cdef int n_elts = -1
      cdef int  *connec
      cdef int  *connec_idx
      cdef long *global_num

      CWP_Mesh_interf_f_poly_block_get(self.local_code_name,
                                       self.id,
                                       i_part,
                                       block_id,
                                       &(n_elts),
                                       &(connec_idx),
                                       &(connec),
                                       &(global_num))

      return {
              'n_elts'     : n_elts,
              'connec_idx' : create_numpy_i(connec_idx, n_elts+1),
              'connec'     : create_numpy_i(connec, connec_idx[n_elts]), # TO DO: correct to do this connec_idx[n_elts] ?
              'global_num' : create_numpy_l(global_num, n_elts)
             }

    def mesh_interf_c_poly_block_set(self,
                                     i_part,
                                     block_id,
                                     n_elts,
                                     n_faces,
                                     NPY.ndarray[NPY.int32_t, mode='c', ndim=1] connec_faces_idx,
                                     NPY.ndarray[NPY.int32_t, mode='c', ndim=1] connec_faces,
                                     NPY.ndarray[NPY.int32_t, mode='c', ndim=1] connec_cells_idx,
                                     NPY.ndarray[NPY.int32_t, mode='c', ndim=1] connec_cells,
                                     NPY.ndarray[NPY.long, mode='c', ndim=1]    global_num):
      """
      Set the properties of a polyhedron block in a interface mesh partition.

      Parameters:
        i_part           (int)                      : Partition identifier
        block_id         (int)                      : Block identifier
        n_elts           (int)                      : Number of elements
        n_faces          (int)                      : Number of faces
        connec_faces_idx (NPY.ndarray[NPY.int32_t]) : Polyhedron face to vertex index (face_vertex_idx[0] = 0 and size = max(cell_face_connec) + 1)
        connec_faces     (NPY.ndarray[NPY.int32_t]) : Polyhedron face to vertex connectivity (size = face_vertex_idx[n_elts])
        connec_cells_idx (NPY.ndarray[NPY.int32_t]) : Polyhedron to face index (cell_face_idx[0] = 0 and size = n_elts + 1)
        connec_cells     (NPY.ndarray[NPY.int32_t]) : Polyhedron to face connectivity (size = cell_face_idx[n_elts])
        global_num       (NPY.ndarray[NPY.long])    : Pointer to global element number (or None)
      """

      CWP_Mesh_interf_c_poly_block_set(self.local_code_name,
                                       self.id,
                                       i_part,
                                       block_id,
                                       n_elts,
                                       n_faces,
                               <int *> connec_faces_idx.data,
                               <int *> connec_faces.data,
                               <int *> connec_cells_idx.data,
                               <int *> connec_cells.data,
                              <long *> global_num.data)

    def mesh_interf_c_poly_block_get(self,
                                     i_part,
                                     block_id):
      """
      Get the properties of a polyhedron block of the interface mesh partition.

      Parameters:
        i_part           (int)                      : Partition identifier
        block_id         (int)                      : Block identifier

      Returns:
        n_elts           (int)                      : Number of elements
        n_faces          (int)                      : Number of faces
        connec_faces_idx (NPY.ndarray[NPY.int32_t]) : Polyhedron face to vertex index (face_vertex_idx[0] = 0 and size = max(cell_face_connec) + 1)
        connec_faces     (NPY.ndarray[NPY.int32_t]) : Polyhedron face to vertex connectivity (size = face_vertex_idx[n_elts])
        connec_cells_idx (NPY.ndarray[NPY.int32_t]) : Polyhedron to face index (cell_face_idx[0] = 0 and size = n_elts + 1)
        connec_cells     (NPY.ndarray[NPY.int32_t]) : Polyhedron to face connectivity (size = cell_face_idx[n_elts])
        global_num       (NPY.ndarray[NPY.long])    : Pointer to global element number (or None)
      """

      cdef int n_elts  = -1
      cdef int n_faces = -1
      cdef int  *connec_faces_idx
      cdef int  *connec_faces
      cdef int  *connec_cells_idx
      cdef int  *connec_cells
      cdef long *global_num

      CWP_Mesh_interf_c_poly_block_get(self.local_code_name,
                                       self.id,
                                       i_part,
                                       block_id,
                                       &(n_elts),
                                       &(n_faces),
                                       &(connec_faces_idx),
                                       &(connec_faces),
                                       &(connec_cells_idx),
                                       &(connec_cells),
                                       &(global_num))

      return {
              'n_elts'           : n_elts,
              'n_faces'          : n_faces,
              'connec_faces_idx' : create_numpy_i(connec_faces_idx, n_faces+1),
              'connec_faces'     : create_numpy_i(connec_faces, connec_faces_idx[n_faces]), # TO DO: correct to do this connec_faces_idx[n_faces] ?
              'connec_cells_idx' : create_numpy_i(connec_cells_idx, n_elts+1),
              'connec_cells'     : create_numpy_i(connec_cells, connec_cells_idx[n_elts]), # TO DO: correct to do this connec_faces_idx[n_faces] ?
              'global_num'       : create_numpy_l(global_num, n_elts)
             }

    def mesh_interf_del(self):
      """
      Delete interface mesh.
      """
      CWP_Mesh_interf_del(self.local_code_name,
                          self.id)

    def mesh_interf_from_cellface_set(self,
                                      i_part,
                                      n_cells,
                                      NPY.ndarray[NPY.int32_t, mode='c', ndim=1] cell_face_idx,
                                      NPY.ndarray[NPY.int32_t, mode='c', ndim=1] cell_face,
                                      n_faces,
                                      NPY.ndarray[NPY.int32_t, mode='c', ndim=1] face_vtx_idx,
                                      NPY.ndarray[NPY.int32_t, mode='c', ndim=1] face_vtx,
                                      NPY.ndarray[NPY.long, mode='c', ndim=1]    global_num):
      """
      Define the interface mesh from a cell to face connectivity.

      Parameters:
        i_part        (int)                      : Current partition
        n_cells       (int)                      : Number of cells
        cell_face_idx (NPY.ndarray[NPY.int32_t]) : Polyhedron to face index (cell_face_idx[0] = 0 and size = n_elts + 1)
        cell_face     (NPY.ndarray[NPY.int32_t]) : Cell to face connectivity (size = cell_face_idx[n_elts])
        n_faces       (int)                      : Number of faces
        face_vtx_idx  (NPY.ndarray[NPY.int32_t]) : Polyhedron face to vertex index (face_vtx_idx[0] = 0 and size = n_faces + 1)
        face_vtx      (NPY.ndarray[NPY.int32_t]) : Face to vertex connectivity (size = face_vtx_idx[n_elts])
        global_num    (NPY.ndarray[NPY.long])    : Pointer to parent element number (or None)
      """

      CWP_Mesh_interf_from_cellface_set(self.local_code_name,
                                        self.id,
                                        i_part,
                                        n_cells,
                                <int *> cell_face_idx.data,
                                <int *> cell_face.data,
                                        n_faces,
                                <int *> face_vtx_idx.data,
                                <int *> face_vtx.data,
                               <long *> global_num.data)

    def mesh_interf_from_faceedge_set(self,
                                      i_part,
                                      n_faces,
                                      NPY.ndarray[NPY.int32_t, mode='c', ndim=1] face_edge_idx,
                                      NPY.ndarray[NPY.int32_t, mode='c', ndim=1] face_edge,
                                      n_edges,
                                      NPY.ndarray[NPY.int32_t, mode='c', ndim=1] edge_vtx_idx,
                                      NPY.ndarray[NPY.int32_t, mode='c', ndim=1] edge_vtx,
                                      NPY.ndarray[NPY.long, mode='c', ndim=1]    global_num):
      """
      Define the surface interface mesh from a face to edge connectivity.

      Parameters:
        i_part        (int)                      : Current partition
        n_faces       (int)                      : Number of faces
        face_edge_idx (NPY.ndarray[NPY.int32_t]) : Polygon to edge index (face_edge_idx[0] = 0 and size =  n_faces + 1)
        face_edge     (NPY.ndarray[NPY.int32_t]) : Face to edge connectivity (size = face_edge_idx[n_faces])
        n_edges       (int)                      : Number of edges
        edge_vtx_idx  (NPY.ndarray[NPY.int32_t]) : Polyhedron edge to vertex index  (edge_vtx_idx[0] = 0 and size = n_edges + 1)
        edge_vtx      (NPY.ndarray[NPY.int32_t]) : Face to vertex connectivity (size = edge_vtx_idx[n_edges])
        global_num    (NPY.ndarray[NPY.long])    : Pointer to parent element number (or None)
      """

      CWP_Mesh_interf_from_faceedge_set(self.local_code_name,
                                        self.id,
                                        i_part,
                                        n_faces,
                                <int *> face_edge_idx.data,
                                <int *> face_edge.data,
                                        n_edges,
                                <int *> edge_vtx_idx.data,
                                <int *> edge_vtx.data,
                               <long *> global_num.data)

    # --> SPATIAL INTERPOLATION
    def spatial_interp_weights_compute(self):
      """
      Compute spatial interpolation weights.
      """

      CWP_Spatial_interp_weights_compute(self.local_code_name,
                                         self.id)

    def spatial_interp_property_set(self,
                                    char *property_name,
                                    char *property_type,
                                    char *property_value):
      """
      Set a property of the spatial interpolation algorithm.

      Parameters:
        property_name  (char*) : Name of the property
        property_type  (char*) : Type of the property ("double" or "int")
        property_value (char*) : Value of the property
      """

      CWP_Spatial_interp_property_set(self.local_code_name,
                                      self.id,
                                      property_name,
                                      property_type,
                                      property_value)

    # --> VISU
    def visu_set(self,
                 int                   freq,
                 CWP_Visu_format_t     visu_format,
                 char                 *format_option):
      """
      Enable visualization output.

      Parameters:
        freq            (int)               : Output frequency
        visu_format     (CWP_Visu_format_t) : Output format to visualize exchanged fieldsDouble
        format_option   (char*)             : Output options "opt1, opt2, ..."
      """

      CWP_Visu_set(self.local_code_name,
                   self.id,
                   freq,
                   visu_format,
                   format_option)

    # --> USER TGT
    def user_tgt_pts_set(self,
                         int i_part,
                         int n_pts,
                         NPY.ndarray[NPY.double_t, mode='c', ndim=1] coord,
                         NPY.ndarray[NPY.long, mode='c', ndim=1] global_num):
      """
      This function must be called if the degrees of freedom locations are CWP_DOF_LOCATION_USER

      Parameters:
         i_part     (int)                      : Current partition
         n_pts      (int)                      : Number of points
         coord      (NPY.ndarray[NPY.double_t) : Coordinates (size = 3 * n_pts)
         global_num (NPY.ndarray[NPY.long])    : Global number or NUL (size = n_pts)
      """

      CWP_User_tgt_pts_set(self.local_code_name,
                           self.id,
                           i_part,
                           n_pts,
                <double *> coord.data,
                  <long *> global_num.data)

    # --> USER INTERPOLATION
    def interp_from_location_unset(self,
                                   char  *src_field_id):
      """
      Unsetting of an user interpolation.

      Parameters:
        src_field_id (char*) : Source field identifier
      """

      CWP_Interp_from_location_unset(self.local_code_name,
                                     self.id,
                                     src_field_id)

    def interp_from_location_set(self,
                                 const char                 *src_field_id,
                                 fct):
      """
      Setting of an user interpolation from location.
      This function takes into account an user interpolation function written with
      void (*\ref CWP_Interp_from_location_t) interface.

      Parameters:
        src_field_id (char*)                      : Source field id
        fct                                       : Function
      """
      global current_cpl_id
      global g_interp_fct
      current_cpl_id = self.id
      g_interp_fct[current_cpl_id] = fct

      CWP_Interp_from_location_set(self.local_code_name,
                                   self.id,
                                   src_field_id,
                                   interp_callback)
      current_cpl_id = -1
