/*----------------------------------------------------------------------------
 * System headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"
#include "pdm_mesh_dist.h"

/*----------------------------------------------------------------------------*/

#ifdef	__cplusplus
extern "C" {
#endif

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief Create a structure to compute distance to a mesh nodal
 *
 * \param [in]   mesh_nodal_id  Mesh nodal identifier
 * \param [in]   n_point_cloud  Number of point cloud
 * \param [in]   comm           MPI communicator
 *
 * \return     Identifier
 */

int
PDM_mesh_dist_create
(
 const int mesh_nodal_id,
 const int n_point_cloud,
 const PDM_MPI_Comm comm
)
{
  return 0;
}


/**
 *
 * \brief Set a point cloud
 *
 * \param [in]   id              Identifier
 * \param [in]   i_point_cloud   Index of point cloud
 * \param [in]   i_part          Index of partition
 * \param [in]   n_points        Number of points
 * \param [in]   coords          Point coordinates
 *
 */

void
PDM_mesh_dist_cloud_set
(
 const int          id,
 const int          i_point_cloud,
 const int          i_part,
 const int          n_points,
 const double      *coords
)
{
}


/**
 *
 * \brief Set a point cloud with initial distance
 *
 * \param [in]   id              Identifier
 * \param [in]   i_point_cloud   Index of point cloud
 * \param [in]   i_part          Index of partition
 * \param [in]   n_points        Number of points
 * \param [in]   initial_dist    Initial distance  
 * \param [in]   coords          Point coordinates
 *
 */

void
PDM_mesh_dist_cloud_with_initial_set
(
 const int          id,
 const int          i_point_cloud,
 const int          i_part,
 const int          n_points,
 const double      *initial_dist,
 const double      *coords
)
{
}


/**
 *
 * \brief Set normal surface mesh
 *
 * \param [in]   id              Identifier
 * \param [in]   i_part          Index of partition
 * \param [in]   normal          Normal
 *
 */

void
PDM_mesh_dist_normal_set
(
 const int          id,
 const int          i_part,
 const double      *normal
)
{
}

  

/**
 *
 * \brief Set normal surface mesh
 *
 * \param [in]   id              Identifier
 * \param [in]   i_part          Index of partition
 * \param [in]   normal          Normal
 *
 */

void
PDM_mesh_dist_center_set
(
 const int          id,
 const int          i_part,
 const double      *center
)
{
}


/**
 *
 * \brief Process merge points
 *
 * \param [in]   id  Identifier
 *
 */

void
PDM_mesh_dist_process
(
 const int id
)
{

  /* 
   * Construction octree distribue avec les sommets de la surface (ou centre face a voir) 
   */

  /*
   *  Pour chaque point recherche du sommet le plus proche (initialisation du bbtree) 
   */

  /*
   *  Construction du dbbtree 
   */
  
  /* 
   * Pour chaque point determination des boites situee a une plus courte distance que le maxima produit par l'octree 
   */

  /* 
   * Repartition des sommets suivant la numerotation absolue en fonction (poids sur le nombre de candidats)
   *  Necessite de faire une block_to_part avec poids
   * 
   * Il faut envoyer les coordonnees des sommets de chaque triangle ou quadrangle ou polygone
   *
   */

  /* 
   * Calcul des distances pour chaque candidat 
   */

  /* 
   * Envoi du resultat selon la repartition initiale des points  
   */
  
  
  


  
}


/**
 *
 * \brief Get mesh distance
 *
 * \param [in]   id              Identifier
 * \param [in]   i_point_cloud   Current cloud
 * \param [in]   i_part          Index of partition
 * \param [out]  dist            Distance
 * \param [out]  proj            Projected point coordinates
 * \param [out]  closest_part    Closest partition
 * \param [out]  closest_elt     Closest element
 *
 */

void
PDM_mesh_dist_get
(
 const int       id,
 const int       i_point_cloud,
       double  **dist,
       double  **proj,
       int     **closest_part,
       int     **closest_elt
)
{
}

/**
 *
 * \brief Free a distance mesh structure
 *
 * \param [in]  id  Identifier
 *
 * \return     Identifier
 */

int
PDM_mesh_dist_free
(
 const int id
)
{
  return 0;
}

  
#ifdef	__cplusplus
}
#endif
