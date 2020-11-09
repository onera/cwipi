#ifndef __PDM_H__
#define __PDM_H__

#include <stdio.h>
#include "pdm_config.h"
#include "pdm_mpi.h"

/*=============================================================================
 * Macro definitions
 *============================================================================*/

#if !defined (__hpux) && !defined (_AIX)
#define PROCF(x, y) x##_
#else
#define PROCF(x, y) x
#endif

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

  PDM_STRIDE_CST = 0,  /*!< Constant stride element */
  PDM_STRIDE_VAR = 1   /*!< Variable stride element */

} PDM_stride_t;

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
  PDM_MESH_ENTITY_VERTEX  = 3   /*!< Vertex entity  */

} PDM_mesh_entities_t;


/**
 * \enum PDM_mesh_nature_t
 * \brief Mesh natures
 *
 */

typedef enum {

  PDM_MESH_NATURE_NODAL_SHARED   = 0,  /*!< Shared PDM_mesh_nodal  */
  PDM_MESH_NATURE_MESH_SETTED    = 1,  /*!< PDm_surface_mesh  */

} PDM_mesh_nature_t;


/**
 * \enum PDM_ownership_t
 * \brief Bool type
 *
 */

typedef enum {

  PDM_OWNERSHIP_KEEP                 = 0, /*!< paradigm will free results */
  PDM_OWNERSHIP_USER                 = 1, /*!< Ownership is give          */
  PDM_OWNERSHIP_UNGET_RESULT_IS_FREE = 2  /*!< Free all memory that not be getted by user */

} PDM_ownership_t;


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
