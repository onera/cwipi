#ifndef __PDM_PART_EXTENSION_H__
#define __PDM_PART_EXTENSION_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_mpi.h"
#include "pdm.h"
#include "pdm_part_domain_interface.h"

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

typedef struct _pdm_part_extension_t PDM_part_extension_t;


typedef enum {

  PDM_EXTEND_FROM_FACE = 0,
  PDM_EXTEND_FROM_EDGE = 1,
  PDM_EXTEND_FROM_VTX  = 2

} PDM_extend_type_t;


/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/


PDM_part_extension_t*
PDM_part_extension_create
(
 const int                n_domain,
 const int               *n_part,
       PDM_extend_type_t  extend_type,
       int                depth,
 const PDM_MPI_Comm       comm,
 const PDM_ownership_t    owner
);

/**
 *
 * \brief Compute a part extension structure
 *
 * \param [in]   part_ext          PDM_part_extension_t
 *
 */
void
PDM_part_extension_compute
(
  PDM_part_extension_t *part_ext
);



// void
// PDM_part_extension_set_part_part_bound
// (
//   PDM_part_extension_t *part_ext,
//   int                   i_domain,
//   int                   i_part,
//   PDM_bound_type_t      bound_type,
//   int                  *bound_proc_idx,
//   int                  *bound_part_idx,
//   int                  *bound
// )
// {
//   // Caratère optionelle sympa, si l'utilisateur

// }

// void
// PDM_part_extension_set_part
// (
//   PDM_part_extension_t *part_ext,
//   int                   i_domain,
//   int                   i_part,
//   int                   n_cell,
//   int                   n_face,
//   int                   n_face_part_bound,
//   int                   n_face_group,
//   int                   n_edge,
//   int                   n_vtx,
//   int                  *cell_face_idx,
//   int                  *cell_face,
//   int                  *face_edge_idx,
//   int                  *face_edge,
//   int                  *face_vtx_idx,
//   int                  *face_vtx,
//   int                  *edge_vtx,
//   int                  *face_bound_idx,
//   int                  *face_bound,
//   int                  *edge_bound_idx,
//   int                  *edge_bound,
//   int                  *vtx_bound_idx,
//   int                  *vtx_bound,
//   PDM_g_num_t          *cell_ln_to_gn,
//   PDM_g_num_t          *face_ln_to_gn,
//   PDM_g_num_t          *edge_ln_to_gn,
//   PDM_g_num_t          *vtx_ln_to_gn,
//   PDM_g_num_t          *face_group_ln_to_gn,
//   double               *vtx_coord
// );


void
PDM_part_extension_set_part
(
  PDM_part_extension_t *part_ext,
  int                   i_domain,
  int                   i_part,
  int                   n_cell,
  int                   n_face,
  int                   n_face_part_bound,
  int                   n_face_group,
  int                   n_edge,
  int                   n_vtx,
  int                  *cell_face_idx,
  int                  *cell_face,
  int                  *face_cell,
  int                  *face_edge_idx,
  int                  *face_edge,
  int                  *face_vtx_idx,
  int                  *face_vtx,
  int                  *edge_vtx,
  int                  *face_bound_idx,
  int                  *face_bound,
  int                  *face_join_idx,
  int                  *face_join,
  int                  *face_part_bound_proc_idx,
  int                  *face_part_bound_part_idx,
  int                  *face_part_bound,
  int                  *vtx_part_bound_proc_idx,
  int                  *vtx_part_bound_part_idx,
  int                  *vtx_part_bound,
  PDM_g_num_t          *cell_ln_to_gn,
  PDM_g_num_t          *face_ln_to_gn,
  PDM_g_num_t          *edge_ln_to_gn,
  PDM_g_num_t          *vtx_ln_to_gn,
  PDM_g_num_t          *face_group_ln_to_gn,
  double               *vtx_coord
);


void
PDM_part_extension_part_domain_interface_shared_set
(
  PDM_part_extension_t        *part_ext,
  PDM_part_domain_interface_t *pdi
);

// void
// PDM_part_extension_n_domain_interface_set
// (
//   PDM_part_extension_t         *part_ext,
//   const int                     n_interface,
//   PDM_domain_interface_mult_t   multidomain_interface
// );

// void
// PDM_part_extension_domain_interface_set
// (
//   PDM_part_extension_t         *part_ext,
//   int                           i_interface,
//   PDM_bound_type_t              interface_kind,
//   int                           n_pair,
//   PDM_g_num_t                  *interface_ids,
//   int                          *interface_dom
// );

// void
// PDM_part_extension_domain_interface_translation_set
// (
//         PDM_part_extension_t *part_ext,
//         int                   i_interface,
//   const double               *vect
// );

// void
// PDM_part_extension_domain_interface_rotation_set
// (
//         PDM_part_extension_t    *part_ext,
//   const int                      i_interface,
//   const double                  *direction,
//   const double                  *center,
//   const double                   angle
// );


void
PDM_part_extension_free
(
 PDM_part_extension_t *part_ext
);


/**
 *
 * \brief Get connectivity
 *
 * \param [in]  part_ext     Pointer to \ref PDM_part_extension_t object
 * \param [in]  i_domain     Id of current domain
 * \param [in]  i_part       Id of current partition
 * \param [in]  mesh_entity  Type of mesh entity
 * \param [out] connect      Entity->group graph (size = \ref connect_idx[\ref n_elt])
 * \param [out] connect_idx  Index for entity->group graph (size = \ref n_elt + 1)
 *
 * \return  n_elt  Number of elements
 *
 */

int
PDM_part_extension_connectivity_get
(
 PDM_part_extension_t     *part_ext,
 int                       i_domain,
 int                       i_part,
 PDM_connectivity_type_t   connectivity_type,
 int                     **connect,
 int                     **connect_idx
);


/**
 *
 * \brief Get global ids
 *
 * \param [in]  part_ext     Pointer to \ref PDM_part_extension_t object
 * \param [in]  i_domain     Id of current domain
 * \param [in]  i_part       Id of current partition
 * \param [in]  mesh_entity  Type of mesh entity
 * \param [out] ln_to_gn     Global ids (size = \ref n_elt)
 *
 * \return  n_elt  Number of elements
 *
 */

int
PDM_part_extension_ln_to_gn_get
(
 PDM_part_extension_t     *part_ext,
 int                       i_domain,
 int                       i_part,
 PDM_mesh_entities_t       mesh_entity,
 PDM_g_num_t             **ln_to_gn
);

int
PDM_part_extension_interface_get
(
 PDM_part_extension_t     *part_ext,
 int                       i_domain,
 int                       i_part,
 PDM_mesh_entities_t       mesh_entity,
 int                     **interface_no
);

/**
 *
 * \brief Get groups
 *
 * \param [in]  part_ext     Pointer to \ref PDM_part_extension_t object
 * \param [in]  i_domain     Id of current domain
 * \param [in]  i_part       Id of current partition
 * \param [in]  mesh_entity  Type of mesh entity
 * \param [out] connect      Entity->group graph (size = \ref connect_idx[\ref n_elt])
 * \param [out] connect_idx  Index for entity->group graph (size = \ref n_elt + 1)
 * \param [out] ln_to_gn     Global ids (size = \ref connect_idx[\ref n_elt])
 *
 * \return  n_elt  Number of elements
 *
 */

int
PDM_part_extension_group_get
(
 PDM_part_extension_t     *part_ext,
 int                       i_domain,
 int                       i_part,
 PDM_mesh_entities_t       mesh_entity,
 int                     **connect,
 int                     **connect_idx,
 PDM_g_num_t             **ln_to_gn
);


/**
 *
 * \brief Get vertex coordinates
 *
 * \param [in]  part_ext     Pointer to \ref PDM_part_extension_t object
 * \param [in]  i_domain     Id of current domain
 * \param [in]  i_part       Id of current partition
 * \param [out] vtx_coord    Vertex coordinates (size = \ref n_vtx * 3)
 *
 * \return  n_vtx  Number of vertices
 *
 */

int
PDM_part_extension_coord_get
(
 PDM_part_extension_t     *part_ext,
 int                       i_domain,
 int                       i_part,
 double                  **vtx_coord
);


int
PDM_part_extension_composed_interface_get
(
 PDM_part_extension_t     *part_ext,
 int                     **composed_interface_idx,
 int                     **composed_interface,
 PDM_g_num_t             **composed_ln_to_gn_sorted
);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_PART_EXTENSION_H__ */
