#ifndef __PDM_H__
#define __PDM_H__

#include <stdio.h>
#include "pdm_config.h"
#include "pdm_mpi.h"

/*=============================================================================
 * Macro definitions
 *============================================================================*/

#if defined (__uxpv__)
#define ARGF_SUPP_CHAINE
#else
#define ARGF_SUPP_CHAINE , ...
#endif

#ifdef PDM_LONG_G_NUM
#define PDM_FMT_G_NUM "%ld"
#else
#define PDM_FMT_G_NUM "%d"
#endif

#define PDM_FMT_L_NUM "%d"

#define PDM_MAX_CHAR_LENGTH 100

#define PDM_UNUSED(x) (void)(x)

/**
 * \brief Interface to for hexa manipulation
 * \param [in]       Hex config flags
 * \param [in]       Flag to test
 * \param [out]      An integer, 1 if the flag is set, 0 if unset
 */

#define PDM_HASFLAG( config, flag )   ((config & flag) == flag)
/**
 * \brief Interface to for hexa manipulation
 * \param Arg:       Hex config flags
 * \param Arg:       Flag to test
 * \param [out]      The flags as integer
 */
#define PDM_SETFLAG( config, flag )   ( config |= flag)

/*
 * \brief Interface to for hexa manipulation
 * \param [in]       Hex config flags
 * \param [in]       flag to test
 * \param [out]      The flags as integer
 */
#define PDM_UNSETFLAG( config, flag ) ( config &= ~flag)

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Type
 *============================================================================*/

#if defined(_WIN64)
# define __int64 long long
#endif


/**
 * \enum PDM_g_num_t
 * \brief Long int in pdm
 *
 */

#ifdef PDM_LONG_G_NUM
typedef  long PDM_g_num_t;
#define  PDM__MPI_G_NUM MPI_LONG
#define  PDM__PDM_MPI_G_NUM PDM_MPI_LONG
#else
typedef int PDM_g_num_t;
#define  PDM__MPI_G_NUM MPI_INT
#define  PDM__PDM_MPI_G_NUM PDM_MPI_INT
#endif

typedef  double PDM_real_t;
#define  PDM__MPI_REAL MPI_DOUBLE
#define  PDM__PDM_MPI_REAL PDM_MPI_DOUBLE

/**
 * \enum PDM_l_num_t
 * \brief Long int in pdm
 *
 */

typedef int PDM_l_num_t;
#define  PDM__MPI_L_NUM MPI_INT
#define  PDM__PDM_MPI_L_NUM PDM_MPI_INT


/**
 * \enum PDM_data_t
 * \brief Type of data
 *
 */

typedef enum {

  PDM_INT    = 0,  /*!< Integer */
  PDM_DOUBLE = 1   /*!< Double */

} PDM_data_t;

/**
 * \enum PDM_part_to_block_distrib_t
 * \brief Type of distribution
 *
 */

typedef enum {

  PDM_STRIDE_CST_INTERLACED = 0, /*!< Constant stride interlaced */
  PDM_STRIDE_CST_INTERLEAVED= 1, /*!< Constant stride interleaved */
  PDM_STRIDE_VAR_INTERLACED = 2  /*!< Variable stride interlaced */

} PDM_stride_t;


/**
 * \enum PDM_mpi_comm_t
 * \brief Framework used for MPI communications
 *
 */

typedef enum {

  PDM_MPI_COMM_KIND_P2P                                = 0, /*!< Peer ro peer (MPI_issend/MPI_irecv) */
  PDM_MPI_COMM_KIND_COLLECTIVE                         = 1, /*!< Collective communications (MPI_Ialltoall, ...) */
  PDM_MPI_COMM_KIND_NEIGHBOR_COLLECTIVE                = 2, /*!< Neighborhood communcations (MPI_I_neighbor_alltoall, ...) */
  PDM_MPI_COMM_KIND_WIN_SHARED_AND_P2P                 = 3, /*!< Shared windows (MPI_Put, MPI_GET, ...) */
  PDM_MPI_COMM_KIND_WIN_SHARED_AND_COLLECTIVE          = 4, /*!< Shared windows (MPI_Put, MPI_GET, ...) */
  PDM_MPI_COMM_KIND_WIN_SHARED_AND_NEIGHBOR_COLLECTIVE = 5, /*!< Shared windows (MPI_Put, MPI_GET, ...) */
  PDM_MPI_COMM_KIND_WIN_RMA                            = 6  /*!< RMA windows (MPI_Put, MPI_GET, ...) */
  
} PDM_mpi_comm_kind_t;


/**
 * \enum PDM_bool_t
 * \brief Bool type
 *
 */

typedef enum {

  PDM_FALSE = 0, /*!< False */
  PDM_TRUE  = 1  /*!< True  */

} PDM_bool_t;


/**
 * \enum PDM_mesh_entities_t
 * \brief Mesh entities
 *
 */

typedef enum {

  PDM_MESH_ENTITY_CELL    = 0,  /*!< Cell entity  */
  PDM_MESH_ENTITY_FACE    = 1,  /*!< Face entity  */
  PDM_MESH_ENTITY_EDGE    = 2,  /*!< Edge entity  */
  PDM_MESH_ENTITY_VERTEX  = 3,  /*!< Vertex entity  */
  PDM_MESH_ENTITY_MAX     = 4   /*!<  */

} PDM_mesh_entities_t;


/**
 * \enum PDM_mesh_nature_t
 * \brief Mesh natures
 *
 */

typedef enum {

  PDM_MESH_NATURE_NODAL_SHARED   = 0,  /*!< Shared PDM_mesh_nodal  */
  PDM_MESH_NATURE_MESH_SETTED    = 1   /*!< PDm_surface_mesh  */

} PDM_mesh_nature_t;


/**
 * \enum PDM_ownership_t
 * \brief Bool type
 *
 */

typedef enum {

  PDM_OWNERSHIP_KEEP                 = 0, /*!< paradigm will free results */
  PDM_OWNERSHIP_USER                 = 1, /*!< Ownership is gives         */
  PDM_OWNERSHIP_UNGET_RESULT_IS_FREE = 2, /*!< Free all memory that not be getted by user */
  PDM_OWNERSHIP_BAD_VALUE            = 3  /*!< Wrong  */

} PDM_ownership_t;

typedef enum {
  PDM_CONNECTIVITY_TYPE_CELL_ELMT   = 0,
  PDM_CONNECTIVITY_TYPE_CELL_CELL   = 1,
  PDM_CONNECTIVITY_TYPE_CELL_FACE   = 2,
  PDM_CONNECTIVITY_TYPE_CELL_EDGE   = 3,
  PDM_CONNECTIVITY_TYPE_CELL_VTX    = 4,
  PDM_CONNECTIVITY_TYPE_FACE_ELMT   = 5,
  PDM_CONNECTIVITY_TYPE_FACE_CELL   = 6,
  PDM_CONNECTIVITY_TYPE_FACE_FACE   = 7,
  PDM_CONNECTIVITY_TYPE_FACE_EDGE   = 8,
  PDM_CONNECTIVITY_TYPE_FACE_VTX    = 9,
  PDM_CONNECTIVITY_TYPE_EDGE_ELMT   = 10,
  PDM_CONNECTIVITY_TYPE_EDGE_CELL   = 11,
  PDM_CONNECTIVITY_TYPE_EDGE_FACE   = 12,
  PDM_CONNECTIVITY_TYPE_EDGE_EDGE   = 13,
  PDM_CONNECTIVITY_TYPE_EDGE_VTX    = 14,
  PDM_CONNECTIVITY_TYPE_VTX_ELMT    = 15,
  PDM_CONNECTIVITY_TYPE_VTX_CELL    = 16,
  PDM_CONNECTIVITY_TYPE_VTX_FACE    = 17,
  PDM_CONNECTIVITY_TYPE_VTX_EDGE    = 18,
  PDM_CONNECTIVITY_TYPE_VTX_VTX     = 19,
  PDM_CONNECTIVITY_TYPE_ELMT_CELL   = 20,
  PDM_CONNECTIVITY_TYPE_ELMT_FACE   = 21,
  PDM_CONNECTIVITY_TYPE_ELMT_EDGE   = 22,
  PDM_CONNECTIVITY_TYPE_ELMT_VTX    = 23,
  PDM_CONNECTIVITY_TYPE_MAX         = 24
} PDM_connectivity_type_t;

typedef enum {
  PDM_BOUND_TYPE_ELMT   = 0,
  PDM_BOUND_TYPE_CELL   = 1,
  PDM_BOUND_TYPE_FACE   = 2,
  PDM_BOUND_TYPE_EDGE   = 3,
  PDM_BOUND_TYPE_VTX    = 4,
  PDM_BOUND_TYPE_MAX    = 5,
} PDM_bound_type_t;

typedef enum {
  PDM_GEOMETRY_KIND_VOLUMIC  = 0,
  PDM_GEOMETRY_KIND_SURFACIC = 1,
  PDM_GEOMETRY_KIND_RIDGE    = 2,
  PDM_GEOMETRY_KIND_CORNER   = 3,
  PDM_GEOMETRY_KIND_MAX      = 4,
} PDM_geometry_kind_t;

typedef enum {

  PDM_MESH_LOCATION_OCTREE,
  PDM_MESH_LOCATION_DBBTREE,

} PDM_mesh_location_method_t;


typedef enum {
  PDM_VTX_KIND_NONE               = 0x0000000, /*=* empty flag, all values set to false */
  PDM_VTX_KIND_NULL               = 0xFFFFFFF, /*=* unsignificant flags value           */
  PDM_VTX_KIND_ON_VOLUME          = 0x0000001, /*=*  0                                  */
  PDM_VTX_KIND_ON_SURFACE         = 0x0000002, /*=*  1                                  */
  PDM_VTX_KIND_ON_RIDGE           = 0x0000004, /*=*  2                                  */
  PDM_VTX_KIND_ON_CORNER          = 0x0000008, /*=*  3                                  */
} PDM_vtx_kind;

typedef enum {
  PDM_REDUCE_OP_MIN,
  PDM_REDUCE_OP_MAX,
  PDM_REDUCE_OP_SUM//,
  // PDM_REDUCE_OP_MEAN
} PDM_reduce_op_t;


typedef enum {

  PDM_DOMAIN_INTERFACE_MULT_NO  = 0,  /*!< Each interface involves only 2 domains */
  PDM_DOMAIN_INTERFACE_MULT_YES = 1,  /*!< Each interface involves several domains */

} PDM_domain_interface_mult_t;


/**
 * \enum PDM_part_split_t
 * \brief Split method
 *
 */

typedef enum {
  PDM_SPLIT_DUAL_WITH_PARMETIS = 1,
  PDM_SPLIT_DUAL_WITH_PTSCOTCH = 2,
  PDM_SPLIT_DUAL_WITH_HILBERT  = 3
} PDM_split_dual_t;


/**
 * \enum PDM_iso_surface_kind_t
 * \brief Type of iso surface
 *
 */
typedef enum {

  PDM_ISO_SURFACE_KIND_PLANE  = 0,
  PDM_ISO_SURFACE_KIND_SPHERE = 1,
  PDM_ISO_SURFACE_KIND_FIELD  = 2,

} PDM_iso_surface_kind_t;




/*=============================================================================
 * Public function prototypes
 *============================================================================*/


/**
 * \brief Finalize PDM
 *
 * This function frees all allocated global variables
 *
 */

void
PDM_Finalize
(
void
);

#ifdef __cplusplus
}
#endif /* __cplusplus */


#endif  /* __PDM_H__ */
