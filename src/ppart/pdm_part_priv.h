#ifndef __PDM_PART_PRIV_H__
#define __PDM_PART_PRIV_H__

#include <stdlib.h>

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_timer.h"
#include "pdm_part.h"
#include "pdm_mpi.h"

/*============================================================================
 * Type definitions
 *============================================================================*/


/**
 * \struct _subpartlayout_t
 * \brief  Partition object
 *
 * _subpartlayout_t define a mesh partition layouts base on sudomaine
 *
 */

typedef struct  _subpartlayout_t {

  int           n_sdom;                  /*!< Number of subDomain                     */
  int           n_face_int;               /*!< Number of Interior face                 */
  int           n_face_ext;               /*!< Number of Exterior face                 */

  /* Idx array of displacement */
  int*          cell_tile_idx;           /*!< Cell Tile Index     (Size = n_sdom + 1)   */
  int*          face_tile_idx;           /*!< Face Tile Index     (Size = n_sdom + 1)   */
  int*          face_bnd_tile_idx;        /*!< Face Bnd Tile Index (Size = n_sdom + 1)   */

  /* Idx array of displacement */
  int*          mask_tile_idx;           /*!< Mask Tile Index   (Size = n_sdom + 1)     */
  int*          cell_vect_tile_idx;       /*!< Cell Tile Index   (Size = n_sdom + 1)     */
  int*          mask_tile_n;             /*!< Mask Tile number  (Size = n_sdom + 1)     */
  int*          cell_vect_tile_n;         /*!< Cell Tile number  (Size = n_sdom + 1)     */
  int*          mask_tile;              /*!< Mask Tile number                         */


} _subpartlayout_t;


/**
 * \struct _part_t
 * \brief  Partition object
 *
 * _part_t define a mesh partition
 *
 */

typedef struct  _part_t {
  int           n_vtx;               /*!< Number of vertices                   */
  int           n_cell;              /*!< Number of cells                      */
  int           n_face;              /*!< Number of faces                      */
  int           n_face_part_bound;     /*!< Number of partitioning boundary faces*/
  int           n_face_group;         /*!< Number of boundary faces             */

  int          *cell_face_idx;        /*!< Cell face connectivity index
                                      (size = n_cell + 1)                      */
  PDM_g_num_t *gcell_face;          /*!< Global numbering cell face connectivity
                                      (size = cell_face_idx[n_cell])             */
  int          *cell_face;           /*!< Cell face connectivity
                                      (size = cell_face_idx[n_cell])             */
  PDM_g_num_t *cell_ln_to_gn;         /*!< Local to global cell numbering
                                      (size = n_cell)                          */
  int          *cell_tag;            /*!< Cell tag
                                      (size = n_cell)                          */

  int          *face_cell;           /*!< Face cell connectivity
                                      (size = 2 * n_face)                      */
  int          *face_vtx_idx;         /*!< Face vertex connectivity index
                                      (size = n_face + 1)                      */
  PDM_g_num_t  *gface_vtx;           /*!< Global numbering face vtx connectivity
                                      (size = face_vtx_idx[n_face])              */
  int          *face_vtx;            /*!< Face vertex connectivity
                                      (size = face_vtx_idx[n_face])              */
  PDM_g_num_t *face_ln_to_gn ;         /*!< Local to global cell numbering
                                      (size = n_face)                          */
  int          *face_tag;            /*!< Tag face connectivity index
                                      (size = n_face)                          */

  int          *face_part_bound_proc_idx;   /*!< Partitioning boundary bloc distribution
                                      (size = n_ranks + 1)                     */
  int          *face_part_bound_part_idx;   /*!< Partitioning boundary bloc distribution
                                      (size = n_partTotal + 1)                 */
  int          *face_part_bound;      /*!< Partitioning boundary faces sorted by
                                         proc, sorted by part in the proc, and
                                         sorted by absolute face number in the part
                                         For each face :
                                          - Face local number
                                          - Connected process
                                          - Connected Partition
                                            on the connected process
                                          - Connected face local number
                                            in the connected partition
                                      (size = 4* n_face_part_bound)              */

  int          *face_group_idx;       /*!< Face group index
                                      (size = n_face_group + 1)                     */
  int          *face_group;          /*!< Face group index
                                      (size = face_group_idx[n_face_group])       */

  PDM_g_num_t *face_group_ln_to_gn;    /*!< Local to global boundary face numbering
                                      (size = face_group_idx[n_face_group])           */

  double       *vtx;                /*!< Vertex coordinates
                                      (size = 3 * n_vtx)                           */
  PDM_g_num_t *vtx_ln_to_gn;          /*!< Local to global vertex numbering
                                      (size = n_vtx)                           */
  int          *vtx_tag;             /*!< Tag vertex
                                      (size = n_vtx)                           */


  const int          *cell_weight;             /*!< Cell weight - For coarse mesh
                                            (size = n_cell)                           */
  const int          *face_weight;             /*!< Face weight - For coarse mesh
                                            (size = nFac)                           */

  int          *cell_color;             /*!< Cell color - For cache blocking
                                            (size = n_cell)                           */
  int          *face_color;             /*!< Face color - For cache blocking
                                            (size = nFac)                           */
  int          *thread_color;             /*!< Thread color - For cache blocking
                                            (size = nThread)                        */
  int          *hyperplane_color;         /*!< Thread color - For cache blocking
                                            (size = nThread)                         */

  int          *new_to_old_order_cell;   /*!< Cell reordering
                                         (size = n_cell)                           */
  int          *new_to_old_order_face;   /*!< Face reordering
                                            (size = nFac)                        */

  _subpartlayout_t *subpartlayout;    /*!< Layouts of subdomain                     */

} _part_t;

/**
 * \struct _PDM_part_t
 * \brief   PPART object
 *
 * _part_t define a parallel mesh partitioning
 *
 */

typedef struct _PDM_part_t {

  /* Local dimensions */

  int                 dn_vtx;         /*!< Number of distributed vertices      */
  int                 dn_cell;        /*!< Number of distributed cells         */
  int                 dn_face;        /*!< Number of distributed faces         */
  int                 n_face_group;    /*!< Number of boundaries                */

  /* Cell definitions */

  const int          *_dcell_face_idx;  /*!< Cell-face connectivity of distributed
                                       cells (size = dcell_face_idx[dn_cell],
                                       shared array) (computed)               */
  const PDM_g_num_t *_dcell_face;    /*!< Tag of distributed cells
                                       (size = dn_cell, shared array)          */
  const int          *_dcell_tag;     /*!< Tag of distributed cells
                                       (size = dn_cell, shared array)          */
  const int          *_dcell_weight;  /*!< Weight of distributed cells
                                       (size = dn_cell, shared array)          */
  const int          *_dcell_part;    /*!< Partitioning of distributed cells
                                       (size = dn_cell, shared array)          */
  int                *dcell_face_idx;  /*!< Cell-face connectivity index of
                                       distributed cells
                                       (size = dn_cell + 1)
                                        computed                             */
  PDM_g_num_t       *dcell_face;     /*!< Cell-face connectivity of distributed
                                       cells (size = dcell_face_idx[dn_cell],
                                       computed                               */
  PDM_g_num_t       *dcell_proc;     /*!< Initial cell distribution on processes
                                       (size = n_rank + 1) (computed)          */

  /* Face definitions */

  const int          *_dface_tag;     /*!< Tag of distributed face
                                       (size = dn_face, shared array)          */
  const PDM_g_num_t *_dface_cell;    /*!< Face-cell connectivity of distributed
                                       faces (size = 2 * dn_face, shared array)
                                       if iface is a boundary face,
                                       _dface_cell[2*iface + 1] = 0           */
  const int          *_dface_vtx_idx;  /*!< Face-vertex connectivity index of
                                       distributed faces (size = dn_face + 1,
                                       shared array)                          */
  const PDM_g_num_t *_dface_vtx;     /*!< Face-vertex connectivity of
                                       distributed faces
                                       (size = dface_vtx_idx[dn_face],shared array)
                                     */
  PDM_g_num_t       *dface_proc;     /*!< Initial face distribution on processes
                                       (size = n_rank + 1)                     */

  PDM_g_num_t       *dface_cell;    /*!< Face-cell connectivity of distributed
                                       faces (size = 2 * dn_face,) computed
                                       if iface is a boundary face,
                                       _dface_cell[2*iface + 1] = 0           */
  /* Vertex definitions */

  const double       *_dvtx_coord;    /*!< Coordinates of ditributed vertices
                                       (size = 3 * dn_vtx, shared array)       */
  const int          *_dvtx_tag;      /*!< Tag of distributed vertices
                                       (size = dn_vtx, shared array)           */
  PDM_g_num_t       *dvtx_proc;      /*!< Initial vertex distribution
                                       (size = n_rank + 1                      */

  /* Face group */

  const int          *_dface_group_idx; /*!< Index of distributed faces list of
                                        each boundary (size = nBound + 1)
                                        or NULL                               */
  const PDM_g_num_t *_dface_group    ;/*!< Distributed faces list of each
                                       boundary (size = dfaceBoundIdx[nBound])
                                       or NULL                                */

  /* Partitioning boundary faces */

  int          *dpart_bound;          /*!< Partitioning boundary faces definition
                                       For each face :
                                       (rank1, part1, lFace1, rank2, part2, lFace2)
                                       or -1 size = 6 * dn_face  */

  /* Dual graph */

  PDM_g_num_t *ddual_graph_idx;      /*!< Dual graph index
                                      (size = dn_cell + 1)                     */
  PDM_g_num_t *ddual_graph;         /*!< Dual graph
                                      (size = dualGraphIdx[dn_cell])           */

  PDM_timer_t *timer;             /*!< Timer */

  double times_elapsed [4];          /*!< Elapsed times :
                                      - Total,
                                      - build dualgraph,
                                      - split graph
                                      - build meshes partition */

  double times_cpu[4];             /*!< CPU times :
                                      - Total,
                                      - build dualgraph,
                                      - split graph
                                      - build meshes partition */

  double times_cpu_u[4];           /*!< User CPU times :
                                      - Total,
                                      - build dualgraph,
                                      - split graph
                                      - build meshes partition */

  double times_cpu_s[4];          /*!< Systeme CPU times :
                                      - Total,
                                      - build dualgraph,
                                      - split graph
                                      - build meshes partition */

  /* Communicator */

  PDM_MPI_Comm      comm;               /*!< Communicator */

  /* Partitions */

  PDM_part_split_t split_method;             /*!< Partitioning method */

  int renum_face_method;                     /*!< Renumbering face method */

  int renum_cell_method;                     /*!< Renumbering cell method */

  int  n_property_cell;                         /*!< Size of cells properties      */
  int  n_property_face;                         /*!< Size of faces properties      */
  const int* renum_properties_cell;           /*!< Renumbering cells properties  */
  const int* renum_properties_face;           /*!< Renumbering faces properties  */


  int          n_part;               /*!< Number of partitions to define
                                      on this process */
  int          tn_part;              /*!< Total number of partitions           */
  int          mn_part;              /*!< Maximum number of partitions to define
                                      on one process */
  int         *dpart_proc;           /*!< Initial cell distribution on processes
                                      (size = n_rank + 1)                      */
  int         *gpart_to_lproc_part;    /*!< For each lobal part number :
                                          - process storing this partition
                                          - local number partition on this
                                            partition
                                         (size = 2*tn_part)                    */
 _part_t     **mesh_parts;           /*!< Partitions built ont this process
                                      (size = n_part)                          */

} _PDM_part_t;

/**
 *
 * \brief Return an initialized part object
 *
 */

static inline _part_t *
_part_create
(
void
 )
{
  _part_t *part = (_part_t *) malloc(sizeof(_part_t));

  part->n_vtx = 0;
  part->n_cell = 0;
  part->n_face = 0;
  part->n_face_part_bound = 0;
  part->cell_face_idx = NULL;
  part->gcell_face = NULL;
  part->cell_face = NULL;
  part->cell_ln_to_gn = NULL;
  part->cell_tag = NULL;
  part->face_cell = NULL;
  part->face_vtx_idx = NULL;
  part->gface_vtx = NULL;
  part->face_vtx = NULL;
  part->face_ln_to_gn = NULL;
  part->face_tag = NULL;
  part->face_part_bound_proc_idx = NULL;
  part->face_part_bound_part_idx = NULL;
  part->face_part_bound = NULL;
  part->face_group_idx = NULL;
  part->face_group = NULL;
  part->face_group_ln_to_gn = NULL;
  part->vtx = NULL;
  part->vtx_ln_to_gn = NULL;
  part->vtx_tag = NULL;
  part->cell_color = NULL;
  part->face_color = NULL;
  part->thread_color = NULL;
  part->hyperplane_color = NULL;
  part->new_to_old_order_cell = NULL;
  part->new_to_old_order_face = NULL;
  part->subpartlayout = NULL;
  return part;
}

#endif
