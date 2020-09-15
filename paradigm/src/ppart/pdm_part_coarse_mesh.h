/*
 * File:   pdm_part_coarse_mesh.h
 * Author: jmagnene
 *
 * Created on July 8, 2016, 9:29 AM
 */

#ifndef PDM_PART_COARSE_MESH_H
#define	PDM_PART_COARSE_MESH_H


/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_part.h"


#if !defined (__hpux) && !defined (_AIX)
#define PROCF(x, y) x##_
#else
#define PROCF(x, y) x
#endif


#ifdef	__cplusplus
extern "C" {
#endif


/**
 *
 * \brief Return an initialized coarse mesh object
 *
 * \param [out]  cmId              Coarse mesh identifier
 *
 * \param [in]   pt_comm           Communicator
 * \param [in]   method            Choice between (1 for ParMETIS or 2 for PT-Scotch)
 * \param [in]   n_part             Number of partitions
 * \param [in]   n_total_part            Total number of partitions
 * \param [in]   have_cell_tag      Presence d'un tableau de tags pour les cellules
 * \param [in]   have_face_tag      Presence d'un tableau de tags pour les faces
 * \param [in]   have_vtx_tag       Presence d'un tableau de tags pour les sommets
 * \param [in]   have_cell_weight   Presence d'un tableau de poids pour les cellules
 * \param [in]   have_face_weight   Presence d'un tableau de poids pour les faces
 * \param [in]   have_face_group    Presence des tableaux de groupes de faces
 */

void
PDM_part_coarse_mesh_create
(
 int                *cmId,
 PDM_MPI_Comm        comm,
 const char*         method,
 const char         *renum_cell_method,
 const char         *renum_face_method,
 const int           n_property_cell,
 const int          *renum_properties_cell,
 const int           n_property_face,
 const int          *renum_properties_face,
 const int           n_part,
 const int           n_total_part,
 const int           n_face_group,
 const int           have_cell_tag,
 const int           have_face_tag,
 const int           have_vtx_tag,
 const int           have_cell_weight,
 const int           have_face_weight,
 const int           have_face_group
);

void
PROCF (pdm_part_coarse_mesh_create_cf, PDM_PART_COARSE_MESH_CREATE_CF)
(
 int                *cmId,
 PDM_MPI_Fint       *fcomm,
 const char         *method,
 const int          *l_method,
 const char         *renum_cell_method,
 const int          *l_renum_cell_method,
 const char         *renum_face_method,
 const int          *l_renum_face_method,
 const int          *n_property_cell,
 const int          *renum_properties_cell,
 const int          *n_property_face,
 const int          *renum_properties_face,
 const int          *n_part,
 const int          *n_total_part,
 const int          *n_face_group,
 const int          *have_cell_tag,
 const int          *have_face_tag,
 const int          *have_vtx_tag,
 const int          *have_cell_weight,
 const int          *have_face_weight,
 const int          *have_face_group
 );


/**
 *
 * \brief Build a coarse mesh
 *
 * \param [in]  cmId               Coarse mesh identifier
 * \param [in]  i_part              Partition identifier
 * \param [in]  n_coarse_cell        Number of cells in the coarse grid
 * \param [in]  n_cell              Number of cells
 * \param [in]  n_face              Number of faces
 * \param [in]  n_face_part_bound     Number of partitioning boundary faces
 * \param [in]  n_vtx               Number of vertices
 * \param [in]  n_face_group         Number of face groups
 * \param [in]  cell_face_idx        Cell to face connectivity index (size = n_cell + 1, numbering : 0 to n-1)
 * \param [in]  cell_face           Cell to face connectivity (size = cell_face_idx[n_cell] = lcell_face
 *                                                             numbering : 1 to n)
 * \param [in]  cell_tag            Cell tag (size = n_cell)
 * \param [in]  cell_ln_to_gn         Cell local numbering to global numbering (size = n_cell, numbering : 1 to n)
 * \param [in]  cell_weight         Cell weight (size = n_cell)
 * \param [in]  face_weight         Face weight (size = n_face)
 * \param [in]  face_cell           Face to cell connectivity  (size = 2 * n_face, numbering : 1 to n)
 * \param [in]  face_vtx_idx         Face to Vertex connectivity index (size = n_face + 1, numbering : 0 to n-1)
 * \param [in]  face_vtx            Face to Vertex connectivity (size = faceVertexIdx[n_face], numbering : 1 to n)
 * \param [in]  face_tag            Face tag (size = n_face)
 * \param [in]  face_ln_to_gn         Face local numbering to global numbering (size = n_face, numbering : 1 to n)
 * \param [in]  vtxCoord           Vertex coordinates (size = 3 * nVertex)
 * \param [in]  vtx_tag             Vertex tag (size = nVertex)
 * \param [in]  vtx_ln_to_gn          Vertex local numbering to global numbering (size = n_vtx, numbering : 1 to n)
 * \param [in]  face_group_idx       Face group index (size = n_face_group + 1, numbering : 1 to n-1)
 * \param [in]  face_group          faces for each group (size = face_group_idx[n_face_group] = lFaceGroup, numbering : 1 to n)
 * \param [in]  face_group_ln_to_gn    Faces global numbering for each group
 *                                  (size = face_group_idx[n_face_group] = lFaceGroup, numbering : 1 to n)
 * \param [in]  face_part_bound_proc_idx  Partitioning boundary faces block distribution from processus (size = n_proc + 1)
 * \param [in]  face_part_bound_part_idx  Partitioning boundary faces block distribution from partition (size = n_total_part + 1)
 * \param [in]  face_part_bound      Partitioning boundary faces (size = 4 * n_face_part_bound)
 *                                       sorted by processus, sorted by partition in each processus, and
 *                                       sorted by absolute face number in each partition
 *                                   For each face :
 *                                        - Face local number (numbering : 1 to n)
 *                                        - Connected process (numbering : 0 to n-1)
 *                                        - Connected Partition
 *                                          on the connected process (numbering :1 to n)
 *                                        - Connected face local number
 *                                          in the connected partition (numbering :1 to n)
 */

void
PDM_part_coarse_mesh_input
(
 int                 cmId,
 int                 i_part,
 const int           n_coarse_cell,
 const int           n_cell,
 const int           n_face,
 const int           n_vtx,
 const int           n_face_group,
 const int           n_face_part_bound,
 const int          *cell_face_idx,
 const int          *cell_face,
 const int          *cell_tag,
 const int          *cell_weight,
 const int          *face_weight,
 const PDM_g_num_t  *cell_ln_to_gn,
 const int          *face_cell,
 const int          *face_vtx_idx,
 const int          *face_vtx,
 const int          *face_tag,
 const PDM_g_num_t  *face_ln_to_gn,
 const double       *vtxCoord,
 const int          *vtx_tag,
 const PDM_g_num_t  *vtx_ln_to_gn,
 const int          *face_group_idx,
 const int          *face_group,
 const PDM_g_num_t  *face_group_ln_to_gn,
 const int          *face_part_bound_proc_idx,
 const int          *face_part_bound_part_idx,
 const int          *face_part_bound
);

void
PROCF (pdm_part_coarse_mesh_input, PDM_PART_COARSE_MESH_INPUT)
(
 int                *cmId,
 int                *i_part,
 const int          *n_coarse_cell,
 const int          *n_cell,
 const int          *n_face,
 const int          *n_vtx,
 const int          *n_face_group,
 const int          *n_face_part_bound,
 const int          *cell_face_idx,
 const int          *cell_face,
 const int          *have_cell_tag,
 const int          *cell_tag,
 const int          *have_cell_weight,
 const int          *cell_weight,
 const int          *have_face_weight,
 const int          *face_weight,
 const PDM_g_num_t *cell_ln_to_gn,
 const int          *face_cell,
 const int          *face_vtx_idx,
 const int          *face_vtx,
 const int          *have_face_tag,
 const int          *face_tag,
 const PDM_g_num_t *face_ln_to_gn,
 const double       *vtxCoord,
 const int          *have_vtx_tag,
 const int          *vtx_tag,
 const PDM_g_num_t *vtx_ln_to_gn,
 const int          *have_face_group,
 const int          *face_group_idx,
 const int          *face_group,
 const PDM_g_num_t *face_group_ln_to_gn,
 const int          *face_part_bound_proc_idx,
 const int          *face_part_bound_part_idx,
 const int          *face_part_bound
);

/**
 *
 * \brief Updates all the arrays dealing with MPI exchanges
 *
 * \param [in] cmId               Coarse mesh identifier
 */

void
PDM_part_coarse_mesh_compute
(
  int cmId
);

void
PROCF (pdm_part_coarse_mesh_compute, PDM_PART_COARSE_MESH_COMPUTE)
(
  int *cmId
);

/**
 *
 * \brief Return a coarse mesh partition dimensions
 *
 * \param [in]   cmId               Coarse mesh identifier
 * \param [in]   i_part              Current partition
 *
 * \param [out]  n_cell              Number of cells
 * \param [out]  n_face              Number of faces
 * \param [out]  n_face_part_bound     Number of partitioning boundary faces
 * \param [out]  n_vtx               Number of vertices
 * \param [out]  n_proc              Number of processus
 * \param [out]  n_total_part             Number of partitions
 * \param [out]  scell_face          Size of cell-face connectivity
 * \param [out]  sface_vtx           Size of face-vertex connectivity
 * \param [out]  sFacePartBound     Size of face_part_bound array
 * \param [out]  sface_group         Size of face_group array
 * \param [out]  sCoarseCellToFineCell  Size of coarseCellToFineCell array
 *
 */

void
PDM_part_coarse_mesh_part_dim_get
(
 const int      cmId,
 int            i_part,
 int           *n_cell,
 int           *n_face,
 int           *n_face_part_bound,
 int           *n_vtx,
 int           *n_proc,
 int           *n_total_part,
 int           *n_face_group,
 int           *scell_face,
 int           *sface_vtx,
 int           *sface_group,
 int           *sCoarseCellToFineCell
);

void
PROCF (pdm_part_coarse_mesh_part_dim_get, PDM_PART_COARSE_MESH_PART_DIM_GET)
(
 const int     *cmId,
 int           *i_part,
 int           *n_cell,
 int           *n_face,
 int           *n_face_part_bound,
 int           *n_vtx,
 int           *n_proc,
 int           *n_total_part,
 int           *n_face_group,
 int           *scell_face,
 int           *sface_vtx,
 int           *sface_group,
 int           *sCoarseCellToFineCell
);

/**
 *
 * \brief Return a mesh partition
 *
 * \param [in]   cmId               Coarse mesh identifier
 * \param [in]   i_part              Current partition
 *
 * \param [out]  cell_face_idx        Cell to face connectivity index (size = n_cell + 1, numbering : 0 to n-1)
 * \param [out]  cell_face           Cell to face connectivity (size = cell_face_idx[n_cell] = lcell_face
 *                                                             numbering : 1 to n)
 * \param [out]  cell_tag            Cell tag (size = n_cell)
 * \param [out]  cell_ln_to_gn         Cell local numbering to global numbering (size = n_cell, numbering : 1 to n)
 * \param [out]  cellInitCellIdx    Array of indexes of the connected partitions (size : n_coarse_cell + 1)
 * \param [out]  cellInitCell       Partitioning array (size : cellInitCellIdx[n_coarse_cell])
 *
 * \param [out]  face_cell           Face to cell connectivity  (size = 2 * n_face, numbering : 1 to n)
 * \param [out]  face_vtx_idx         Face to Vertex connectivity index (size = n_face + 1, numbering : 0 to n-1)
 * \param [out]  face_vtx            Face to Vertex connectivity (size = faceVertexIdx[n_face], numbering : 1 to n)
 * \param [out]  face_tag            Face tag (size = n_face)
 * \param [out]  face_ln_to_gn         Face local numbering to global numbering (size = n_face, numbering : 1 to n)
 * \param [out]  faceInitFace       Coarse face - fine face connectivity (size = nCoarseFace)
 *
 * \param [out]  vtxCoord           Vertex coordinates (size = 3 * n_vtx)
 * \param [out]  vtx_tag             Vertex tag (size = n_vtx)
 * \param [out]  vtx_ln_to_gn          Vertex local numbering to global numbering (size = n_vtx, numbering : 1 to n)
 * \param [out]  vtxInitVtx         Coarse vertex - fine vertex connectivity (size = nCoarseVtx)
 *
 * \param [out]  face_group_idx       Face group index (size = n_face_group + 1, numbering : 1 to n-1)
 * \param [out]  face_group          faces for each group (size = face_group_idx[n_face_group] = lFaceGroup, numbering : 1 to n)
 * \param [out]  face_group_ln_to_gn    Faces global numbering for each group
 *                                  (size = face_group_idx[n_face_group] = lFaceGroup, numbering : 1 to n)
 *
 * \param [out]  face_part_bound_proc_idx  Partitioning boundary faces block distribution from processus (size = n_proc + 1)
 * \param [out]  face_part_bound_part_idx  Partitioning boundary faces block distribution from partition (size = n_total_part + 1)
 * \param [out]  face_part_bound      Partitioning boundary faces (size = 4 * n_face_part_bound)
 *                                       sorted by processus, sorted by partition in each processus, and
 *                                       sorted by absolute face number in each partition
 *                                   For each face :
 *                                        - Face local number (numbering : 1 to n)
 *                                        - Connected process (numbering : 0 to n-1)
 *                                        - Connected Partition
 *                                          on the connected process (numbering :1 to n)
 *                                        - Connected face local number
 *                                          in the connected partition (numbering :1 to n)
 *
 */

void
PDM_part_coarse_mesh_part_get
(
 const int    cmId,
 const int    i_part,
 int          **cell_face_idx,
 int          **cell_face,
 int          **cell_tag,
 PDM_g_num_t **cell_ln_to_gn,
 int          **cellInitCellIdx,
 int          **cellInitCell,
 int          **face_cell,
 int          **face_vtx_idx,
 int          **face_vtx,
 int          **face_tag,
 PDM_g_num_t **face_ln_to_gn,
 int          **faceGroupInitFaceGroup,
 int          **faceInitFace,
 double       **vtxCoord,
 int          **vtx_tag,
 PDM_g_num_t **vtx_ln_to_gn,
 int          **vtxInitVtx,
 int          **face_group_idx,
 int          **face_group,
 PDM_g_num_t **face_group_ln_to_gn,
 int          **face_part_bound_proc_idx,
 int          **face_part_bound_part_idx,
 int          **face_part_bound
);

void
PROCF (pdm_part_coarse_mesh_part_get, PDM_PART_COARSE_MESH_PART_GET)
(
 int          *cmId,
 int          *i_part,
 int          *cell_face_idx,
 int          *cell_face,
 int          *cell_tag,
 PDM_g_num_t  *cell_ln_to_gn,
 int          *cellInitCellIdx,
 int          *cellInitCell,
 int          *face_cell,
 int          *face_vtx_idx,
 int          *face_vtx,
 int          *face_tag,
 PDM_g_num_t *face_ln_to_gn,
 int          *faceGroupInitFaceGroup,
 int          *faceInitFace,
 double       *vtxCoord,
 int          *vtx_tag,
 PDM_g_num_t *vtx_ln_to_gn,
 int          *vtxInitVtx,
 int          *face_group_idx,
 int          *face_group,
 PDM_g_num_t *face_group_ln_to_gn,
 int          *face_part_bound_proc_idx,
 int          *face_part_bound_part_idx,
 int          *face_part_bound
);

/**
 *
 * \brief Return a mesh partition
 *
 * \param [in]   ppart_id               ppart identifier
 * \param [in]   i_part                 Current partition
 * \param [out]  cell_color             Cell Color (size = n_cell)
 * \param [out]  face_color             Face Color (size = n_face)
 */

void PDM_part_coarse_color_get
(
 const int   cmId,
 const int   i_part,
       int **cell_color,
       int **face_color,
       int **thread_color,
       int **hyperplane_color
);

/**
 *
 * \brief Free coarse mesh
 *
 * \param [in]   cmId        Coarse mesh identifier
 *
 */

void
PDM_part_coarse_mesh_free
(
 const int    cmId
);

void
PROCF (pdm_part_coarse_mesh_free, PDM_PART_COARSE_MESH_FREE)
(
  int *cmId
);

/**
 *
 * \brief Return times
 *
 * \param [in]   cmId        coarse mesh identifier
 * \param [out]  elapsed     elapsed times (size = 18)
 * \param [out]  cpu         cpu times (size = 18)
 * \param [out]  cpu_user    user cpu times (size = 18)
 * \param [out]  cpu_sys     system cpu times (size = 18)
 *
 */

void PDM_part_coarse_mesh_time_get
(
 int       cmId,
 double  **elapsed,
 double  **cpu,
 double  **cpu_user,
 double  **cpu_sys
);

void
PROCF (pdm_part_coarse_mesh_time_get, PDM_PART_COARSE_MESH_TIME_GET)
(
 int      *cmId,
 double   *elapsed,
 double   *cpu,
 double   *cpu_user,
 double   *cpu_sys
 );

/**
 *
 * \brief Displays all the arrays of a coarse mesh
 *
 * \param [in]   cmId        Coarse mesh identifier
 *
 */

void
PDM_part_coarse_mesh_display
(
 const int    cmId
);

void
PROCF (pdm_part_coarse_mesh_display, PDM_PART_COARSE_MESH_DISPLAY)
(
  int *cmId
);


#ifdef	__cplusplus
}
#endif

#endif	/* PDM_PART_COARSE_MESH_H */

