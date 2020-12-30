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

  int             n_domain;
  const int      *n_part;

  _part_t  **parts;

  /* Store for each depth / each domain / each part */
  int **neighbor_idx;
  int **neighbor_desc;
  int  *n_entity_bound;

  /* Graph of cell */
  int **dist_neighbor_cell_n;
  int **dist_neighbor_cell_idx;
  int **dist_neighbor_cell_desc;

  int ****graph_comm_cell_opp;
  int ****graph_comm_cell_opp_idx;
  int ****graph_comm_cell;
  int ****cell_list;
  int  ***n_cell_per_bound;

  int ***cell_cell_idx;
  int ***cell_cell;

  /* This one is only on the border and contains only border cells */
  int ****cell_cell_extended_idx;
  int ****cell_cell_extended;

  int  *n_tot_part_by_domain;

  int **entity_cell_idx;
  int **entity_cell_n;
  int **entity_cell;

  int **entity_cell_opp_idx;
  int **entity_cell_opp_n;
  int **entity_cell_opp;

};

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_PART_EXTENSION_PRIV_H__ */
