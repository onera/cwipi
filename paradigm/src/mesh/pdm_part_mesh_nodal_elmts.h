#ifndef __PDM_PART_MESH_NODAL_ELMTS_H__
#define __PDM_PART_MESH_NODAL_ELMTS_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"
#include "pdm_io.h"
#include "pdm_mesh_nodal.h"

/*=============================================================================
 * Macro definition
 *============================================================================*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Types definition
 *============================================================================*/

typedef struct _pdm_part_mesh_nodal_elmts_t PDM_part_mesh_nodal_elmts_t;

/*=============================================================================
 * Global variables
 *============================================================================*/

/*=============================================================================
 * Public function interfaces
 *============================================================================*/

/**
 * \brief Create a Mesh nodal structure
 *
 * \param [in]   n_part   Number of partition on the current process
 * \param [in]   comm     MPI communicator
 *
 * \return       New mesh nodal handle
 *
 */
PDM_part_mesh_nodal_elmts_t*
PDM_part_mesh_nodal_elmts_create
(
const int          mesh_dimension,
const int          n_part,
const PDM_MPI_Comm comm
);

int
PDM_part_mesh_nodal_elmts_add
(
      PDM_part_mesh_nodal_elmts_t *pmne,
const PDM_Mesh_nodal_elt_t         t_elt
);

int
PDM_part_mesh_nodal_elmts_ho_add
(
      PDM_part_mesh_nodal_elmts_t *pmne,
const PDM_Mesh_nodal_elt_t         t_elt,
const int                          order,
const char                        *ho_ordering
);

void
PDM_part_mesh_nodal_elmts_free
(
 PDM_part_mesh_nodal_elmts_t *pmne
);

void
PDM_part_mesh_nodal_elmts_std_set
(
      PDM_part_mesh_nodal_elmts_t *pmne,
const int                          id_block,
const int                          id_part,
const int                          n_elt,
const int                         *connec,
const PDM_g_num_t                 *numabs,
const int                         *parent_num,
const PDM_g_num_t                 *parent_entity_g_num,
      PDM_ownership_t              owner
);

/**
 * \brief Define a polygon block
 *
 * \param [in]  mesh           Pointer to \ref PDM_Mesh_nodal object
 * \param [in]  id_block       Block identifier
 * \param [in]  id_part        Partition identifier
 * \param [in]  n_elt          Number of elements
 * \param [in]  connect_idx    Connectivity index (size = \ref n_elt + 1)
 * \param [in]  connect        Connectivity (size = \ref connect_idx[\ref n_elt])
 * \param [in]  numabs         Global numbering
 * \param [in]  parent_num     Parent numbering or NULL
 *
 */

void
PDM_part_mesh_nodal_elmts_block_poly2d_set
(
      PDM_part_mesh_nodal_elmts_t *pmne,
const int                          id_block,
const int                          id_part,
const int                          n_elt,
const int                         *connec_idx,
const int                         *connec,
const PDM_g_num_t                 *numabs,
const int                         *parent_num,
      PDM_ownership_t              owner
);

/**
 * \brief Define a polyhedra block
 *
 * \param [in]  mesh           Pointer to \ref PDM_part_mesh_nodal_elmts_t object
 * \param [in]  id_block       Block identifier
 * \param [in]  id_part        Partition identifier
 * \param [in]  n_elt          Number of polyhedra
 * \param [in]  n_face         Number of faces used to describe polyhedra
 * \param [in]  facvtx_idx     Index of face vertex connectivity
 * \param [in]  facvtx         Face vertex connectivity
 * \param [in]  cellfac_idx    Index of cell face connectivity
 * \param [in]  cellfac        Cell face connectivity
 * \param [in]  numabs         Global numbering
 * \param [in]  parent_num     Parent numbering or NULL
 *
 */

void
PDM_part_mesh_nodal_elmts_block_poly3d_set
(
      PDM_part_mesh_nodal_elmts_t *pmne,
const int                          id_block,
const int                          id_part,
const int                          n_elt,
const int                          n_face,
const int                         *facvtx_idx,
const int                         *facvtx,
const PDM_g_num_t                 *face_ln_to_gn,
const int                         *cellfac_idx,
const int                         *cellfac,
const PDM_g_num_t                 *numabs,
const int                         *parent_num,
      PDM_ownership_t              owner
);



void
PDM_part_mesh_nodal_elmts_block_std_get
(
      PDM_part_mesh_nodal_elmts_t  *pmne,
const int                           id_block,
const int                           id_part,
      int                         **connec,
      PDM_g_num_t                 **numabs,
      int                         **parent_num,
      PDM_g_num_t                 **parent_entity_g_num
);

void
PDM_part_mesh_nodal_elmts_block_std_ho_get
(
      PDM_part_mesh_nodal_elmts_t  *pmne,
const int                           id_block,
const int                           id_part,
      int                         **connec,
      PDM_g_num_t                 **numabs,
      int                         **parent_num,
      PDM_g_num_t                 **parent_entity_g_num,
      int                          *order,
const char                        **ho_ordering
);

/**
 * \brief Return a polygon block description
 *
 * \param [in]  pmne           Pointer to \ref PDM_part_mesh_nodal_elmts object
 * \param [in]  id_block       Block identifier
 * \param [in]  id_part        Partition identifier
 * \param [out] connect_idx    Connectivity index (size = \ref n_elt + 1)
 * \param [out] connect        Connectivity (size = \ref connect_idx[\ref n_elt])
 *
 */

void
PDM_part_mesh_nodal_elmts_block_poly2d_get
(
       PDM_part_mesh_nodal_elmts_t  *pmne,
 const int                           id_block,
 const int                           id_part,
       int                         **connec_idx,
       int                         **connec
);

/**
 * \brief Get the cell-vertex connectivity of a polyhedra block
 *
 * \param [in]  pmne           Pointer to \ref PDM_part_mesh_nodal_elmts object
 * \param [in]  id_block       Block identifier
 * \param [in]  id_part        Partition identifier
 * \param [out] cell_vtx_idx   Index of cell vertex connectivity
 * \param [out] cell_vtx       Cell vertex connectivity
 *
 */

void
PDM_part_mesh_nodal_elmts_block_poly3d_cell_vtx_connect_get
(
       PDM_part_mesh_nodal_elmts_t  *pmne,
 const int                           id_block,
 const int                           id_part,
       int                         **cell_vtx_idx,
       int                         **cell_vtx
);


void
PDM_part_mesh_nodal_elmts_block_poly3d_get
(
      PDM_part_mesh_nodal_elmts_t  *pmne,
const int                           id_block,
const int                           id_part,
      int                          *n_face,
      PDM_g_num_t                 **face_ln_to_gn,
      int                         **face_vtx_idx,
      int                         **face_vtx,
      PDM_g_num_t                 **numabs,
      int                         **cell_face_idx,
      int                         **cell_face,
      int                         **parent_num,
      PDM_g_num_t                 **parent_entity_g_num
);

int
PDM_part_mesh_nodal_elmts_block_n_elt_get
(
      PDM_part_mesh_nodal_elmts_t *pmne,
const int                          id_block,
const int                          id_part
);

int
PDM_part_mesh_nodal_elmts_n_section_get
(
 PDM_part_mesh_nodal_elmts_t *pmne
);

int *
PDM_part_mesh_nodal_elmts_sections_id_get
(
  PDM_part_mesh_nodal_elmts_t *pmne
);

PDM_Mesh_nodal_elt_t
PDM_part_mesh_nodal_elmts_block_type_get
(
      PDM_part_mesh_nodal_elmts_t *mesh,
const int                          id_block
);

PDM_part_mesh_nodal_elmts_t*
PDM_part_mesh_nodal_create_from_part3d
(
  const int                n_part,
  const int               *n_cell,
  const int               *n_face,
  const int              **face_vtx_idx,
  const int              **face_vtx,
  const PDM_g_num_t      **face_ln_to_gn,
  const int              **cell_face_idx,
  const int              **cell_face,
  const double           **vtx_coord,
  const PDM_g_num_t      **numabs,
        PDM_MPI_Comm       comm
);

int *
PDM_part_mesh_nodal_elmts_parent_num_get
(
      PDM_part_mesh_nodal_elmts_t *pmne,
const int                          id_block,
const int                          id_part
);


PDM_part_mesh_nodal_elmts_t*
PDM_part_mesh_nodal_create_from_part2d
(
  const int                n_part,
  const int               *n_face,
  const int               *n_edge,
  const int               *n_vtx,
  const int              **edge_vtx_idx,
  const int              **edge_vtx,
  const int              **face_edge_idx,
  const int              **face_edge,
  const PDM_g_num_t      **numabs,
        PDM_MPI_Comm       comm
);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_PART_MESH_NODAL_ELMTS_H__ */
