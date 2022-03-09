#include "pdm_configf.h"
#ifdef PDM_LONG_G_NUM
  integer, parameter :: pdm_g_num_s = 8
#else
  integer, parameter :: pdm_g_num_s = 4
#endif
  integer, parameter :: pdm_l_num_s = 4

  integer, parameter :: pdm_max_char_length = 100

  integer, parameter :: PDM_STRIDE_CST_INTERLACED  = 0
  integer, parameter :: PDM_STRIDE_CST_INTERLEAVED = 1
  integer, parameter :: PDM_STRIDE_VAR_INTERLACED  = 2

  integer, parameter :: PDM_OWNERSHIP_KEEP                 = 0
  integer, parameter :: PDM_OWNERSHIP_USER                 = 1
  integer, parameter :: PDM_OWNERSHIP_UNGET_RESULT_IS_FREE = 2

  integer, parameter :: PDM_MESH_ENTITY_CELL    = 0  !< Cell entity
  integer, parameter :: PDM_MESH_ENTITY_FACE    = 1  !< Face entity
  integer, parameter :: PDM_MESH_ENTITY_EDGE    = 2  !< Edge entity
  integer, parameter :: PDM_MESH_ENTITY_VERTEX  = 3   !< Vertex entity

  integer, parameter :: PDM_MESH_NATURE_NODAL_SHARED   = 0  !< PDM_mesh_nodal
  integer, parameter :: PDM_MESH_NATURE_MESH_SETTED    = 1  !< PDm_surface_mesh

  integer, parameter :: PDM_MPI_COMM_KIND_P2P                                = 0
  integer, parameter :: PDM_MPI_COMM_KIND_COLLECTIVE                         = 1
  integer, parameter :: PDM_MPI_COMM_KIND_NEIGHBOR_COLLECTIVE                = 2
  integer, parameter :: PDM_MPI_COMM_KIND_WIN_SHARED_AND_P2P                 = 3
  integer, parameter :: PDM_MPI_COMM_KIND_WIN_SHARED_AND_COLLECTIVE          = 4
  integer, parameter :: PDM_MPI_COMM_KIND_WIN_SHARED_AND_NEIGHBOR_COLLECTIVE = 5
  integer, parameter :: PDM_MPI_COMM_KIND_WIN_RMA                            = 6
