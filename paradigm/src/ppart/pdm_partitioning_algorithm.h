#ifndef __PDM_PARTITIONING_ALGORITHM_H__
#define __PDM_PARTITIONING_ALGORITHM_H__

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

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/**
 *  \brief Gather the entities splitted by the partitioner
 *   (usually cells) to their attributed partition, using the array mapping
 *   entities id to their assigned partition number.
*/
int
PDM_part_assemble_partitions
(
 const PDM_MPI_Comm    comm,
       PDM_g_num_t    *part_distribution,
 const PDM_g_num_t    *entity_distribution,
 const int            *dentity_to_part,
       int           **pn_entity,
       PDM_g_num_t  ***pentity_ln_to_gn
);

/**
 *  \brief Construct the face->cell connectivity from the cell->face connectivity
 */
void
PDM_part_reverse_pcellface
(
  const int         n_part,
  const int        *n_cell,
  const int        *n_face,
  const int       **pcell_face_idx,
  const int       **pcell_face,
        int      ***pface_cell
);

/**
 *  \brief Reorient the boundary faces such that they have a outward normal for the boundary cell.
 */
void
PDM_part_reorient_bound_faces
(
  const int         n_part,
  const int        *np_face,
        int       **pface_cell,
  const int       **pcell_face_idx,
        int       **pcell_face,
  const int       **pface_vtx_idx,
        int       **pface_vtx,
        int       **pface_edge_idx,
        int       **pface_edge
);

/**
 *  \brief Recover partitioned entity groups (cell, face, vertex) from distributed
 *   entity groups.
 */
void
PDM_part_distgroup_to_partgroup
(
 const PDM_MPI_Comm      comm,
 const PDM_g_num_t      *entity_distribution,
 const int               n_group,
 const int              *dgroup_idx,
 const PDM_g_num_t      *dgroup,
 const int               n_part,
 const int              *pn_entity,
 const PDM_g_num_t     **pentity_ln_to_gn,
       int            ***pgroup_idx,
       int            ***pgroup,
       PDM_g_num_t    ***pgroup_ln_to_gn
);

/**
 *  \brief Generated the partitioned connectivity (entity->child_elements) associated
 *   to the given distributed connectivity, using element distribution and element local
 *   to global numbering. In addition, return the partitioned number of unique child_element
 *   and the corresponding local to global numbering for the child elements.
 */
void
PDM_part_dconnectivity_to_pconnectivity_sort
(
 const PDM_MPI_Comm    comm,
 const PDM_g_num_t    *entity_distribution,
 const int            *dconnectivity_idx,
 const PDM_g_num_t    *dconnectivity,
 const int             n_part,
 const int            *pn_entity,
 const PDM_g_num_t   **pentity_ln_to_gn,
       int           **pn_child_entity,
       PDM_g_num_t  ***pchild_ln_to_gn,
       int          ***pconnectivity_idx,
       int          ***pconnectivity
);

void
PDM_part_dconnectivity_to_pconnectivity_sort_single_part
(
 const PDM_MPI_Comm    comm,
 const PDM_g_num_t    *entity_distribution,
 const int            *dconnectivity_idx,
 const PDM_g_num_t    *dconnectivity,
 const int             pn_entity,
 const PDM_g_num_t    *pentity_ln_to_gn,
       int            *pn_child_entity,
       PDM_g_num_t   **pchild_ln_to_gn,
       int           **pconnectivity_idx,
       int           **pconnectivity
);

/**
 *  \brief Generated the partitioned connectivity (entity->child_elements) associated
 *   to the given distributed connectivity, using element distribution and element local
 *   to global numbering. In addition, return the partitioned number of unique child_element
 *   and the corresponding local to global numbering for the child elements.
 *
 *   -- multi-section version
 */
void
PDM_part_multi_dconnectivity_to_pconnectivity_sort
(
 const PDM_MPI_Comm    comm,
 const int             n_part,
 const int             n_section,
 const int            *section_idx,
       PDM_g_num_t   **entity_distribution,
       int            *dconnectivity_idx,
       PDM_g_num_t    *dconnectivity,
       int           **pn_entity,
       PDM_g_num_t  ***pentity_ln_to_gn,
       int           **pn_child_entity,
       PDM_g_num_t  ***pchild_ln_to_gn,
       int         ****pconnectivity_idx,
       int         ****pconnectivity
);

/**
 *  \brief Generated the partitioned connectivity (entity->child_elements) associated
 *   to the given distributed connectivity, using element distribution and element local
 *   to global numbering. In addition, return the partitioned number of unique child_element
 *   and the corresponding local to global numbering for the child elements.
 */
void
PDM_part_dconnectivity_to_pconnectivity_hash
(
 const PDM_MPI_Comm    comm,
 const PDM_g_num_t    *entity_distribution,
 const int            *dconnectivity_idx,
 const PDM_g_num_t    *dconnectivity,
 const int             n_part,
 const int            *pn_entity,
 const PDM_g_num_t   **pentity_ln_to_gn,
       int           **pn_child_entity,
       PDM_g_num_t  ***pchild_ln_to_gn,
       int          ***pconnectivity_idx,
       int          ***pconnectivity
);


/**
 *  \brief Generated the communication information at the partition interfaces for the
 *   given entity. The communication data associates to
 *   each partitioned entity belonging to an (internal) interface the 4-tuple
 *   (local id, opposite proc number, opposite part number on this proc, local id in the
 *   opposite partition).
 */
void
PDM_part_generate_entity_graph_comm
(
 const PDM_MPI_Comm   comm,
 const PDM_g_num_t   *part_distribution,
 const PDM_g_num_t   *entity_distribution,
 const int            n_part,
 const int           *pn_entity,
 const PDM_g_num_t  **pentity_ln_to_gn,
 const int          **pentity_hint,
       int         ***pproc_bound_idx,
       int         ***ppart_bound_idx,
       int         ***pentity_bound,
       int         ***pentity_priority
);

/**
 *  \brief Recover partitioned coordinates from distributed coordinates and
 *   vertex ln_to_gn indirection.
 *   This function basically calls PDM_block_to_part on to exchange vertex coordinates.
 *
 */
void
PDM_part_dcoordinates_to_pcoordinates
(
  const PDM_MPI_Comm    comm,
  const int             n_part,
  const PDM_g_num_t    *vertex_distribution,
  const double         *dvtx_coord,
  const int            *pn_vtx,
  const PDM_g_num_t   **pvtx_ln_to_gn,
        double       ***pvtx_coord
);

void
PDM_part_dfield_to_pfield
(
  const PDM_MPI_Comm    comm,
  const int             n_part,
  size_t                s_data,
  const PDM_g_num_t    *field_distribution,
  const unsigned char  *dfield,
  const int            *pn_field,
  const PDM_g_num_t   **pfield_ln_to_gn,
        unsigned char ***pfield
);


void
PDM_part_dfield_to_pfield2
(
  const PDM_MPI_Comm     comm,
  const int              n_part,
  size_t                 s_data,
  PDM_stride_t           t_stride,
  const PDM_g_num_t     *field_distribution,
  const int             *dfield_stri,
  const unsigned char   *dfield,
  const int             *pn_field,
  const PDM_g_num_t    **pfield_ln_to_gn,
  int                 ***pfield_stride,
        unsigned char ***pfield
);

void
PDM_extend_mesh
(
 const PDM_MPI_Comm    comm,
 const PDM_g_num_t    *part_distribution,
 const PDM_g_num_t    *entity_distribution,
 const int            *dentity_to_part,
 const int             n_part,
 const PDM_g_num_t    *dual_graph_idx,
 const PDM_g_num_t    *dual_graph,
 const int            *pn_entity,
       PDM_g_num_t   **pentity_ln_to_gn,
       int           **pn_entity_extented,
       PDM_g_num_t  ***pentity_ln_to_gn_extended
);


void
PDM_part_dentity_group_to_pentity_group
(
  const PDM_MPI_Comm     comm,
  const int              n_part,
  const PDM_g_num_t     *entity_distribution,
  const int             *dentity_group_idx,
  const int             *dentity_group,
  const int             *pn_entity,
  const PDM_g_num_t    **pentity_ln_to_gn,
  int                 ***pentity_group_idx,
  int                 ***pentity_group
);

void
PDM_setup_connectivity_idx
(
  int           dn_entity1,
  int           stride,
  PDM_g_num_t  *dentity1_dentity2,
  int         **dentity1_dentity2_idx,
  PDM_g_num_t **dentity1_dentity2_new
);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif  /* __PDM_PARTITIONING_ALGORITHM_H__ */
