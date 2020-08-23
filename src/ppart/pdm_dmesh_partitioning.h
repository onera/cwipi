#ifndef __PDM_DMESH_PARTITIONING_H__
#define __PDM_DMESH_PARTITIONING_H__

#include <stdio.h>
#include "pdm.h"
#include "pdm_mpi.h"

/*=============================================================================
 * Macro definitions
 *============================================================================*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Type
 *============================================================================*/

typedef enum {
  PDM_PART_NONE               = 0x0000000, /*=* empty flag, all values set to false */
  PDM_PART_NULL               = 0xFFFFFFF, /*=* unsignificant flags value           */
  PDM_PART_CELL_VTX           = 0x0000001, /*=*  0 */
  PDM_PART_FACE_CELL          = 0x0000002, /*=*  1 */
  PDM_PART_CELL_FACE          = 0x0000004, /*=*  2 */
  PDM_PART_EDGE_VTX           = 0x0000008, /*=*  3 */
  PDM_PART_FACE_VTX           = 0x0000010, /*=*  4 */
  PDM_PART_CELL_WEIGHT        = 0x0000020, /*=*  5 */
  PDM_PART_EDGE_WEIGHT        = 0x0000040, /*=*  6 */
  PDM_PART_CELL_TAG           = 0x0000080, /*=*  7 */
  PDM_PART_FACE_TAG           = 0x0000100, /*=*  8 */
  PDM_PART_VTX_TAG            = 0x0000200, /*=*   */
  PDM_PART_FACE_GROUP         = 0x0000400, /*=*   */
  PDM_PART_CELL_GROUP         = 0x0000800, /*=*   */
  PDM_PART_VTX_GROUP          = 0x0001000, /*=*   */
  PDM_PART_EDGE_GROUP         = 0x0002000, /*=*   */
  PDM_PART_CELL_PART          = 0x0004000, /*=*   */
  PDM_PART_VTX_PART           = 0x0008000, /*=*   */
  PDM_PART_GRAPH_COMM_FACE    = 0x0010000, /*=*   */
  PDM_PART_GRAPH_COMM_VTX     = 0x0020000, /*=*   */
  PDM_PART_GRAPH_COMM_EDGE    = 0x0040000, /*=*   */
  PDM_PART_VTX_COORD          = 0x0080000, /*=*   */
  PDM_PART_OWNDATA            = 0x0100000, /*=*   */
  PDM_PART_UNUSED1            = 0x0200000, /*=*   */
  PDM_PART_UNUSED2            = 0x0400000, /*=*   */
  PDM_PART_OWNDAT3            = 0x0800000, /*=*   */
  PDM_PART_UNUSED4            = 0x1000000, /*=*   */
  PDM_PART_UNUSED5            = 0x2000000, /*=*   */
  PDM_PART_UNUSED6            = 0x4000000, /*=*   */
  PDM_PART_UNUSED7            = 0x8000000, /*=*   */
} PDM_partitioning_option_t;

/**
 * \enum PDM_part_split_t
 * \brief Split method
 *
 */

typedef enum {
  PDM_PARTITIONING_WITH_PARMETIS = 1,
  PDM_PARTITIONING_WITH_PTSCOTCH = 2,
  PDM_PARTITIONING_WITH_HILBERT  = 3
} PDM_partitioning_method_t;

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/**
 *
 * \brief Return _dmesh_partitioning_t object from it identifier
 *
 * \param [in]   Comm         MPI Communicator
 * \return       dmpartitioning_id
 */
int
PDM_dmesh_partitioning_create
(
 const PDM_MPI_Comm              comm,
 const PDM_partitioning_method_t split_method
);

/**
 *
 * \brief Return _dmesh_partitioning_t object from it identifier
 *
 * \param [in]   Comm         MPI Communicator
 * \return       dmpartitioning_id
 */
void
PDM_dmesh_partitioning_compute
(
 const int  dmesh_partitioning_id,
 const int  input_flags,
 const int  queries_flags
);

/**
 *
 * \brief Setup structure with a dmesh
 *
 * \param [in]   dmpartitioning_id        ppart identifier
 *
 */
void
PDM_dmesh_partitioning_set_from_dmesh
(
 const int dmesh_partitioning_id,
 const int dmesh_id
);

/**
 *
 * \brief Setup structure with a dmesh
 *
 * \param [in]   dmpartitioning_id        ppart identifier
 *
 */
void
PDM_dmesh_partitioning_part_get
(
 const int  dmesh_partitioning_id,
 const int  part_id,
 const int  input_field_key,
       int *field
);

/**
 *
 * \brief Setup structure with a dmesh
 *
 * \param [in]   dmpartitioning_id        ppart identifier
 *
 */
void
PDM_dmesh_partitioning_get
(
 const int   dmesh_partitioning_id,
 const int   input_field_key,
       int **field
);
/**
 *
 * \brief Free
 *
 * \param [in]   id           Identifier
 *
 */
void
PDM_dmesh_partitioning_free
(
 const int id
);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif  /* __PDM_DMESH_PARTITIONING_H__ */
