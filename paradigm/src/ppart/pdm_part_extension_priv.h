#ifndef __PDM_PART_EXTENSION_PRIV_H__
#define __PDM_PART_EXTENSION_PRIV_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

// #include "pdm_multipart.h"
#include "pdm_part_priv.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation csback to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/**
 * \struct _pdm_part_extension_t
 * \brief  Distributed cube
 *
 * _dcube_t define a distributed mesh of a cube
 *
 */

struct _pdm_part_extension_t {
  PDM_MPI_Comm      comm;            /*!< MPI communicator                          */
  PDM_ownership_t   owner;           /*!< Which have the responsabilities of results*/

  PDM_extend_type_t extend_type;
  int               depth;

  int             n_domain;
  const int      *n_part;
  int            *n_part_idx;

  _part_t  **parts;

  /* Store for each depth / each domain / each part */
  int **neighbor_idx;
  int **neighbor_desc;
  int  *n_entity_bound;
  int  *n_cell;
  int  *n_cell_border;

  /* Graph of cell */
  int **dist_neighbor_cell_n;
  int **dist_neighbor_cell_idx;
  int **dist_neighbor_cell_desc;

  int ***cell_cell_idx;
  int ***cell_cell;

  /* This one is only on the border and contains only border cells */
  int ***cell_cell_extended_idx;
  int ***cell_cell_extended_n;
  int ***cell_cell_extended;
  int  **border_cell_list;

  int  **cell_cell_extended_pruned_idx;
  int  **cell_cell_extended_pruned;

  int  *n_tot_part_by_domain;

  int **entity_cell_idx;
  int **entity_cell_n;
  int **entity_cell;

  int **entity_cell_opp_idx;
  int **entity_cell_opp_n;
  int **entity_cell_opp;

  /* Internal results */
  int **face_face_extended_idx;
  int **face_face_extended;

  int **edge_edge_extended_idx;
  int **edge_edge_extended;

  int **vtx_vtx_extended_idx;
  int **vtx_vtx_extended;

  /* Results */
  int **border_cell_face_idx;
  int **border_cell_face;

  int **border_face_edge_idx;
  int **border_face_edge;

  int **border_edge_vtx_idx;
  int **border_edge_vtx;

  int **border_face_vtx_idx;
  int **border_face_vtx;

  // int  *n_face_group;
  int **border_face_group_idx;
  int **border_face_group;

  PDM_g_num_t **border_cell_ln_to_gn;
  PDM_g_num_t **border_face_ln_to_gn;
  PDM_g_num_t **border_edge_ln_to_gn;
  PDM_g_num_t **border_vtx_ln_to_gn;

  PDM_g_num_t **border_face_group_ln_to_gn;

  double **border_vtx;


};

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_PART_EXTENSION_PRIV_H__ */
