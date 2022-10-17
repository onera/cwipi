#ifndef PDM_MESH_LOCATION_H
#define PDM_MESH_LOCATION_H

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"
#include "pdm_mesh_nodal.h"
#include "pdm_part_to_part.h"

#ifdef	__cplusplus
extern "C" {
#endif

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

typedef struct _pdm_mesh_location_t PDM_mesh_location_t;

/*============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief Create a structure to compute the location of point clouds inta a mesh
 *
 * \param [in]   mesh_nature    Nature of the mesh
 * \param [in]   n_point_cloud  Number of point cloud
 * \param [in]   comm           MPI communicator
 * \param [in]   owner          Ownership
 *
 * \return     Pointer to \ref PDM_mesh_location object
 *
 */

PDM_mesh_location_t*
PDM_mesh_location_create
(
 const PDM_mesh_nature_t mesh_nature,
 const int               n_point_cloud,
 const PDM_MPI_Comm      comm,
 const PDM_ownership_t   owner
);


/**
 *
 * \brief Set the number of partitions of a point cloud
 *
 * \param [in]   id              Pointer to \ref PDM_mesh_location object
 * \param [in]   i_point_cloud   Index of point cloud
 * \param [in]   n_part          Number of partitions
 *
 */

void
PDM_mesh_location_n_part_cloud_set
(
       PDM_mesh_location_t *ml,
 const int                  i_point_cloud,
 const int                  n_part
);



/**
 *
 * \brief Set a point cloud
 *
 * \param [in]   id              Pointer to \ref PDM_mesh_location object
 * \param [in]   i_point_cloud   Index of point cloud
 * \param [in]   i_part          Index of partition
 * \param [in]   n_points        Number of points
 * \param [in]   coords          Point coordinates
 * \param [in]   gnum            Point global number
 *
 */

void
PDM_mesh_location_cloud_set
(
       PDM_mesh_location_t *ml,
 const int                  i_point_cloud,
 const int                  i_part,
 const int                  n_points,
       double              *coords,
       PDM_g_num_t         *gnum
);


/**
 *
 * \brief Get a point cloud
 *
 * \param [in]   id              Pointer to \ref PDM_mesh_location object
 * \param [in]   i_point_cloud   Index of point cloud
 * \param [in]   i_part          Index of partition
 * \param [out]  n_points        Number of points
 * \param [out]  coords          Point coordinates
 * \param [out]  gnum            Point global number
 *
 */

void
PDM_mesh_location_cloud_get
(
       PDM_mesh_location_t  *ml,
 const int                   i_point_cloud,
 const int                   i_part,
       int                  *n_points,
       double              **coords,
       PDM_g_num_t         **gnum
);


/**
 *
 * \brief Get the number of located points
 *
 * \param [in]   id              Pointer to \ref PDM_mesh_location object
 * \param [in]   i_point_cloud   Index of point cloud
 * \param [in]   i_part          Index of partition
 *
 * \return     The number of located points
 *
 */

int
PDM_mesh_location_n_located_get
(
       PDM_mesh_location_t *ml,
 const int                  i_point_cloud,
 const int                  i_part
);


/**
 *
 * \brief Get the number of unlocated points
 *
 * \param [in]   id              Pointer to \ref PDM_mesh_location object
 * \param [in]   i_point_cloud   Index of point cloud
 * \param [in]   i_part          Index of partition
 *
 * \return     The number of unlocated points
 *
 */

int
PDM_mesh_location_n_unlocated_get
(
       PDM_mesh_location_t *ml,
 const int                  i_point_cloud,
 const int                  i_part
);


/**
 *
 * \brief Get the list of unlocated points
 *
 * \param [in]   id              Pointer to \ref PDM_mesh_location object
 * \param [in]   i_point_cloud   Index of point cloud
 * \param [in]   i_part          Index of partition
 *
 * \return     The list of unlocated points
 *
 */

int *
PDM_mesh_location_unlocated_get
(
       PDM_mesh_location_t *ml,
 const int                  i_point_cloud,
 const int                  i_part
);


/**
 *
 * \brief Get the list of located points
 *
 * \param [in]   id              Pointer to \ref PDM_mesh_location object
 * \param [in]   i_point_cloud   Index of point cloud
 * \param [in]   i_part          Index of partition
 *
 * \return     The list of located points
 *
 */

int *
PDM_mesh_location_located_get
(
       PDM_mesh_location_t *ml,
 const int                  i_point_cloud,
 const int                  i_part
);



/**
 *
 * \brief Set the mesh nodal
 *
 * \param [in]   id             Pointer to \ref PDM_mesh_location object
 * \param [in]   mesh_nodal_id  Mesh nodal Pointer to \ref PDM_mesh_location object
 *
 */

void
PDM_mesh_location_shared_nodal_mesh_set
(
 PDM_mesh_location_t *ml,
 PDM_Mesh_nodal_t    *mesh_nodal
);


/**
 *
 * \brief Set global data of a mesh
 *
 * \param [in]   id             Pointer to \ref PDM_mesh_location object
 * \param [in]   n_g_cell       Global number of cells
 * \param [in]   n_g_face       Global number of faces
 * \param [in]   n_g_vtx        Global number of vertices
 * \param [in]   n_part         Number of partition
 *
 */

void
PDM_mesh_location_mesh_global_data_set
(
       PDM_mesh_location_t *ml,
 const int                  n_part
);


/**
 *
 * \brief Set a part of a mesh
 *
 * \param [in]   id            Pointer to \ref PDM_mesh_location object
 * \param [in]   i_part        Partition to define
 * \param [in]   n_cell        Number of cells
 * \param [in]   cell_face_idx Index in the cell -> face connectivity
 * \param [in]   cell_face     cell -> face connectivity
 * \param [in]   cell_ln_to_gn Local cell numbering to global cel numbering
 * \param [in]   n_face        Number of faces
 * \param [in]   face_vtx_idx  Index in the face -> vertex connectivity
 * \param [in]   face_vtx      face -> vertex connectivity
 * \param [in]   face_ln_to_gn Local face numbering to global face numbering
 * \param [in]   n_vtx         Number of vertices
 * \param [in]   coords        Coordinates
 * \param [in]   vtx_ln_to_gn  Local vertex numbering to global vertex numbering
 *
 */

void
PDM_mesh_location_part_set
(
       PDM_mesh_location_t *ml,
 const int                  i_part,
 const int                  n_cell,
 const int                 *cell_face_idx,
 const int                 *cell_face,
 const PDM_g_num_t         *cell_ln_to_gn,
 const int                  n_face,
 const int                 *face_vtx_idx,
 const int                 *face_vtx,
 const PDM_g_num_t         *face_ln_to_gn,
 const int                  n_vtx,
 const double              *coords,
 const PDM_g_num_t         *vtx_ln_to_gn
);

/**
 *
 * \brief Set a part of a mesh
 *
 * \param [in]   id                     Pointer to \ref PDM_mesh_location object
 * \param [in]   i_part                 Partition to define
 * \param [in]   n_cell                 Number of cells
 * \param [in]   is_elmt_select_by_user Flag to determine if user want or no to extract current cell
 *
 */
void
PDM_mesh_location_user_extract_set
(
       PDM_mesh_location_t *ml,
 const int                  i_part,
 const int                 *is_elmt_select_by_user
);



/**
 *
 * \brief Set a part of a mesh (2d version)
 *
 * \param [in]   id            Pointer to \ref PDM_mesh_location object
 * \param [in]   i_part        Partition to define
 * \param [in]   n_cell        Number of cells
 * \param [in]   cell_edge_idx Index in the cell -> edge connectivity
 * \param [in]   cell_edge     cell -> edge connectivity
 * \param [in]   cell_ln_to_gn Local cell numbering to global cel numbering
 * \param [in]   n_edge        Number of edges
 * \param [in]   edge_vtx_idx  Index in the edge -> vertex connectivity
 * \param [in]   edge_vtx      edge -> vertex connectivity
 * \param [in]   edge_ln_to_gn Local edge numbering to global edge numbering
 * \param [in]   n_vtx         Number of vertices
 * \param [in]   coords        Coordinates
 * \param [in]   vtx_ln_to_gn  Local vertex numbering to global vertex numbering
 *
 */

void
PDM_mesh_location_part_set_2d
(
       PDM_mesh_location_t *ml,
 const int                  i_part,
 const int                  n_cell,
 const int                 *cell_edge_idx,
 const int                 *cell_edge,
 const PDM_g_num_t         *cell_ln_to_gn,
 const int                  n_edge,
 const int                 *edge_vtx_idx,
 const int                 *edge_vtx,
 const PDM_g_num_t         *edge_ln_to_gn,
 const int                  n_vtx,
 const double              *coords,
 const PDM_g_num_t         *vtx_ln_to_gn
);

/**
 *
 * \brief Set the tolerance for bounding boxes
 *
 * \param [in]   id              Pointer to \ref PDM_mesh_location object
 * \param [in]   tol             Tolerance
 *
 */

void
PDM_mesh_location_tolerance_set
(
       PDM_mesh_location_t *ml,
 const double               tol
);


/**
 *
 * \brief Set the method for computing location
 *
 * \param [in]   id              Pointer to \ref PDM_mesh_location object
 * \param [in]   method          Method
 *
 */

void
PDM_mesh_location_method_set
(
       PDM_mesh_location_t        *ml,
 const PDM_mesh_location_method_t  method
);


/**
 *
 * \brief Compute point location
 *
 * \param [in]   id  Pointer to \ref PDM_mesh_location object
 *
 */

void
PDM_mesh_location_compute
(
PDM_mesh_location_t        *ml
);


/**
 *
 * \brief Compute point location optim
 *
 * \param [in]   id  Pointer to \ref PDM_mesh_location object
 *
 */

void
PDM_mesh_location_compute_optim
(
 PDM_mesh_location_t        *ml
);


/**
 *
 * \brief Get point location
 *
 * \param [in]   id                    Pointer to \ref PDM_mesh_location object
 * \param [in]   i_point_cloud         Current cloud
 * \param [in]   i_part                Index of partition of the cloud
 * \param [out]  n_points              Number of points in point cloud
 * \param [out]  coord                 Coordinates of points in point cloud
 * \param [out]  location              The global number of the closest element if the point is located,
 *                                     -1 otherwise
 *
 */

void
PDM_mesh_location_point_location_get
(
       PDM_mesh_location_t  *ml,
 const int                   i_point_cloud,
 const int                   i_part,
       PDM_g_num_t         **location,
       double              **dist2,
       double              **projected_coord
);


/**
 *
 * \brief get cell vertex connectivity
 *
 * \param [in]   id                    Pointer to \ref PDM_mesh_location object
 * \param [in]   i_part                Index of partition of the cloud
 * \param [out]  cell_vtx_idx          Index in (size = n_elt + 1)
 * \param [out]  cell_vtx              Cell vertex connectivity
 *
 */

void
PDM_mesh_location_cell_vertex_get
(
       PDM_mesh_location_t  *ml,
 const int                   i_part,
       int                 **cell_vtx_idx,
       int                 **cell_vtx
);


/**
 *
 * \brief Get point list located in elements
 *
 * \param [in]   id                      Pointer to \ref PDM_mesh_location object
 * \param [in]   i_part                  Index of partition of the mesh
 * \param [in]   i_point_cloud           Index of cloud
 * \param [out]  elt_pts_inside_idx      Points index (size = n_elt + 1)
 * \param [out]  points_gnum             Points global number
 * \param [out]  points_coords           Points coordinates
 * \param [out]  points_uvw              Points parametric coordinates in elements
 * \param [out]  points_weights_idx      Interpolation weights index (size = elt_pts_inside_idx[n_elt] + 1)
 * \param [out]  points_weights          Interpolation weights
 * \param [out]  points_dist2            Distance element-points (dist < 0 if the point is inside)
 * \param [out]  points_projected_coords Point projection on element if the point is outside
 *
 */

void
PDM_mesh_location_points_in_elt_get
(
       PDM_mesh_location_t  *ml,
 const int                   i_part,
 const int                   i_point_cloud,
       int                 **elt_pts_inside_idx,
       PDM_g_num_t         **points_gnum,
       double              **points_coords,
       double              **points_uvw,
       int                 **points_weights_idx,
       double              **points_weights,
       double              **points_dist2,
       double              **points_projected_coords
);


/**
 *
 * \brief Free a mesh location structure
 *
 * \param [in]  ml       Pointer to \ref PDM_mesh_location object
 *
 */

void
PDM_mesh_location_free
(
 PDM_mesh_location_t  *ml
);


/**
 *
 * \brief Get the number of cells
 *
 * \param [in]  id       Pointer to \ref PDM_mesh_location object
 * \param [in]  i_part   Index of partition of the mesh
 *
 * \return Number of cells
 */

int
PDM_mesh_location_n_cell_get
(
       PDM_mesh_location_t *ml,
 const int                  i_part
);


/**
 *
 * \brief  Dump elapsed an CPU time
 *
 * \param [in]  id       Pointer to \ref PDM_mesh_location object
 *
 */

void
PDM_mesh_location_dump_times
(
PDM_mesh_location_t *ml
);

PDM_Mesh_nodal_t*
PDM_mesh_location_mesh_nodal_get
(
PDM_mesh_location_t *ml
);


/**
 * Enable reverse results computation (To call PDM_mesh_location_points_in_elt_get)
 */

void
PDM_mesh_location_reverse_results_enable
(
PDM_mesh_location_t *ml
);


/**
 * \brief Get part_to_part object to exchange data between
 * the source mesh and a target point cloud (both in user frame)
 *
 * \param [in ] ml         Pointer to \ref PDM_mesh_location_t object
 * \param [in ] icloud     Point cloud ID
 * \param [out] ptp        Pointer to \ref PDM_part_to_part_t object
 * \param [in ] ownership  Ownership for ptp
 *
 */

void
PDM_mesh_location_part_to_part_get
(
       PDM_mesh_location_t  *ml,
 const int                   icloud,
       PDM_part_to_part_t  **ptp,
       PDM_ownership_t       ownership
 );

#ifdef	__cplusplus
}
#endif

#endif // PDM_MESH_LOCATION_H
