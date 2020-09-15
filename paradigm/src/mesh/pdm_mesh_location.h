#ifndef PDM_MESH_LOCATION_H
#define PDM_MESH_LOCATION_H

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"

#ifdef	__cplusplus
extern "C" {
#endif

  /*============================================================================
   * Macro definitions
   *============================================================================*/

  /*============================================================================
   * Type definitions
   *============================================================================*/
  typedef enum {

    PDM_MESH_LOCATION_OCTREE,
    PDM_MESH_LOCATION_DBBTREE,

  } PDM_mesh_location_method_t;

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
   *
   * \return     Identifier
   *
   */

  int
  PDM_mesh_location_create
  (
   const PDM_mesh_nature_t mesh_nature,
   const int n_point_cloud,
   const PDM_MPI_Comm comm
   );

  void
  PDM_mesh_location_create_cf
  (
   const PDM_mesh_nature_t mesh_nature,
   const int n_point_cloud,
   const PDM_MPI_Fint comm,
   int *id
   );


  /**
   *
   * \brief Set the number of partitions of a point cloud
   *
   * \param [in]   id              Identifier
   * \param [in]   i_point_cloud   Index of point cloud
   * \param [in]   n_part          Number of partitions
   *
   */

  void
  PDM_mesh_location_n_part_cloud_set
  (
   const int          id,
   const int          i_point_cloud,
   const int          n_part
   );



  /**
   *
   * \brief Set a point cloud
   *
   * \param [in]   id              Identifier
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
   const int          id,
   const int          i_point_cloud,
   const int          i_part,
   const int          n_points,
   double      *coords,
   PDM_g_num_t *gnum
   );


  /**
   *
   * \brief Set the mesh nodal
   *
   * \param [in]   id             Identifier
   * \param [in]   mesh_nodal_id  Mesh nodal identifier
   *
   */

  void
  PDM_mesh_location_shared_nodal_mesh_set
  (
   const int  id,
   const int  mesh_nodal_id
   );


  /**
   *
   * \brief Set global data of a mesh
   *
   * \param [in]   id             Identifier
   * \param [in]   n_g_cell       Global number of cells
   * \param [in]   n_g_face       Global number of faces
   * \param [in]   n_g_vtx        Global number of vertices
   * \param [in]   n_part         Number of partition
   *
   */

  void
  PDM_mesh_location_mesh_global_data_set
  (
   const int         id,
   const int         n_part
   );


  /**
   *
   * \brief Set a part of a mesh
   *
   * \param [in]   id            Identifier
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
   const int          id,
   const int          i_part,
   const int          n_cell,
   const int         *cell_face_idx,
   const int         *cell_face,
   const PDM_g_num_t *cell_ln_to_gn,
   const int          n_face,
   const int         *face_vtx_idx,
   const int         *face_vtx,
   const PDM_g_num_t *face_ln_to_gn,
   const int          n_vtx,
   const double      *coords,
   const PDM_g_num_t *vtx_ln_to_gn
   );

  /**
 *
 * \brief Set a part of a mesh (2d version)
 *
 * \param [in]   id            Identifier
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
 const int          id,
 const int          i_part,
 const int          n_cell,
 const int         *cell_edge_idx,
 const int         *cell_edge,
 const PDM_g_num_t *cell_ln_to_gn,
 const int          n_edge,
 const int         *edge_vtx_idx,
 const int         *edge_vtx,
 const PDM_g_num_t *edge_ln_to_gn,
 const int          n_vtx,
 const double      *coords,
 const PDM_g_num_t *vtx_ln_to_gn
 );

  /**
   *
   * \brief Set the tolerance for bounding boxes
   *
   * \param [in]   id              Identifier
   * \param [in]   tol             Tolerance
   *
   */

  void
  PDM_mesh_location_tolerance_set
  (
   const int    id,
   const double tol
   );


  /**
   *
   * \brief Set the method for computing location
   *
   * \param [in]   id              Identifier
   * \param [in]   method          Method
   *
   */

  void
  PDM_mesh_location_method_set
  (
   const int                        id,
   const PDM_mesh_location_method_t method
   );


  /**
   *
   * \brief Compute point location
   *
   * \param [in]   id  Identifier
   *
   */

  void
  PDM_mesh_location_compute
  (
   const int id
   );

  /**
   *
   * \brief Get mesh location
   *
   * \param [in]   id                    Identifier
   * \param [in]   i_point_cloud         Current cloud
   * \param [in]   i_part                Index of partition of the cloud
   * \param [out]  n_points              Number of points in point cloud
   * \param [out]  coord                 Coordinates of points in point cloud
   * \param [out]  g_num                 Global numbers of points in point cloud
   * \param [out]  location              The global number of the closest element if the point is located,
   *                                     -1 otherwise
   *
   */

  void
  PDM_mesh_location_get
  (
   const int     id,
   const int     i_point_cloud,
   const int     i_part,
   int          *n_points,
   double      **coord,
   PDM_g_num_t **g_num,
   PDM_g_num_t **location,
   int         **weights_idx,
   double      **weights
   );

  /**
   *
   * \brief Free a locationd mesh structure
   *
   * \param [in]  id       Identifier
   * \param [in]  partial  if partial is equal to 0, all data are removed.
   *                       Otherwise, results are kept.
   *
   */

  void
  PDM_mesh_location_free
  (
   const int id,
   const int partial
   );


  /**
   *
   * \brief  Dump elapsed an CPU time
   *
   * \param [in]  id       Identifier
   *
   */

  void
  PDM_mesh_location_dump_times
  (
   const int id
   );

  int
  PDM_mesh_location_mesh_nodal_id_get
  (
   const int id
   );

#ifdef	__cplusplus
}
#endif

#endif // PDM_MESH_LOCATION_H
