/* 
 * File:   pdm_part_coarse_mesh_priv.h
 * Author: jmagnene
 *
 * Created on July 12, 2016, 10:54 AM
 */

#include "pdm_part_priv.h"
#include "pdm_timer.h"
#include "pdm_mpi.h"

#ifndef PDM_PART_COARSE_MESH_PRIV_H
#define	PDM_PART_COARSE_MESH_PRIV_H

#ifdef	__cplusplus
extern "C" {
#endif

/*============================================================================
 * Type definitions
 *============================================================================*/

    

/**
 * \struct _coarse_part_t
 * \brief  Coarse partition object
 * 
 * _part_t defines a coarse partition 
 *
 */

typedef struct  {
  
  _part_t      *part;        //Coarse mesh
  
  int *coarseCellCellIdx;    //Array of indexes of the connected partitions (size : nCoarseCell + 1)
  
  int *coarseCellCell;       //Partitioning array (size : coarseCellCellIdx[nCoarseCell]) 
  
  int *coarseFaceGroupToFineFaceGroup; //Coarse face group - fine face group connectivity (size = faceGroupIdx[nFaceGroup])
  
  int *coarseFaceToFineFace; //Coarse face - fine face connectivity (size = nCoarseFace)

  int *coarseVtxToFineVtx;   //Coarse vertex - fine vertex connectivity (size = nCoarseVtx)

} _coarse_part_t;


/**
 * \struct _coarse_part_t
 * \brief  Coarse partition object
 * 
 * _coarse_mesh_t defines a coarse mesh 
 *
 */

typedef struct  {

  /* Partitions */          
    
  int nPart;        /*!< Number of partitions to define
                                      on this process */ 
  
  int method;       /*!< Partitioning method */
  int nTPart;       /*!< Total number of partitions */
  int nFaceGroup;   /*!< Number of boundaries */
  
  int have_cellTag;
  int have_faceTag;
  int have_vtxTag;
  int have_cellWeight;
  int have_faceWeight;
  int have_faceGroup;

  //TIMER
  
  PDM_timer_t *timer;             /*!< Timer */ 

  double times_elapsed [18];          /*!< Elapsed times :
                                      - Total,
                                      - build dualgraph,
                                      - split graph
                                      - build meshes partition */

  double times_cpu[18];             /*!< CPU times :
                                      - Total,
                                      - build dualgraph,
                                      - split graph
                                      - build meshes partition */

  double times_cpu_u[18];           /*!< User CPU times :
                                      - Total,
                                      - build dualgraph,
                                      - split graph
                                      - build meshes partition */

  double times_cpu_s[18];          /*!< Systeme CPU times :
                                      - Total,
                                      - build dualgraph,
                                      - split graph
                                      - build meshes partition */

  
  /* Communicator */
  
  PDM_MPI_Comm  comm;   /*!< Communicator */
    
  _part_t        **part_ini;               /*!< Partition: fine mesh                            */
  
  _coarse_part_t **part_res;               //Coarse mesh
  
} _coarse_mesh_t;
    

/*============================================================================
 * Private function definitions
 *============================================================================*/


/**
 *
 * \brief Return an initialized coarse part object
 *
 */

static inline _coarse_part_t * 
_coarse_part_create
(
void
 )
{
  _coarse_part_t *cp = (_coarse_part_t *) malloc(sizeof(_coarse_part_t));
  cp->part = _part_create();
  
  cp->coarseCellCell = NULL;
  
  cp->coarseCellCellIdx = NULL;
  
  cp->coarseFaceGroupToFineFaceGroup = NULL;
  
  cp->coarseFaceToFineFace = NULL; 

  cp->coarseVtxToFineVtx = NULL;
  
  return cp;
  
}

    
/**
 *
 * \brief Return an initialized coarse part object
 * 
 * \param [in]   pt_comm           Communicator
 * \param [in]   method            Choice between (1 for ParMETIS or 2 for PT-Scotch)
 * \param [in]   nPart             Number of partitions
 * \param [in]   nTPart            Total number of partitions
 * \param [in]   nFaceGroup        Number of boundaries
 * \param [in]   have_cellTag      Presence d'un tableau de tags pour les cellules
 * \param [in]   have_faceTag      Presence d'un tableau de tags pour les faces
 * \param [in]   have_vtxTag       Presence d'un tableau de tags pour les sommets
 * \param [in]   have_cellWeight   Presence d'un tableau de poids pour les cellules
 * \param [in]   have_faceWeight   Presence d'un tableau de poids pour les faces
 * \param [in]   have_faceGroup    Presence des tableaux de groupes de faces
 */

static inline _coarse_mesh_t * 
_coarse_mesh_create
(
 const PDM_MPI_Comm  comm,        
 const int           method,
 const int           nPart,
 const int           nTPart,
 const int           nFaceGroup,
 const int           have_cellTag,
 const int           have_faceTag,
 const int           have_vtxTag,
 const int           have_cellWeight,
 const int           have_faceWeight,
 const int           have_faceGroup

 )
{     
   _coarse_mesh_t *cm = (_coarse_mesh_t *) malloc(sizeof(_coarse_mesh_t));

   cm->nPart = nPart;
   cm->comm = comm; 
   
   cm->method = method;
   cm->nTPart = nTPart;
   cm->nFaceGroup = nFaceGroup;

   cm->have_cellTag = have_cellTag;
   cm->have_faceTag = have_faceTag;
   cm->have_vtxTag = have_vtxTag;
   cm->have_cellWeight = have_cellWeight;
   cm->have_faceWeight = have_faceWeight;
   cm->have_faceGroup = have_faceGroup;
   
   cm->part_ini = malloc(sizeof(_part_t *) * nPart); //On déclare un tableau de partitions
   
   cm->part_res = malloc(sizeof(_coarse_part_t *) * nPart);
   
   for (int i = 0; i < nPart; i++) {
     cm->part_ini[i] = _part_create(); 
     
     cm->part_res[i] = _coarse_part_create();     
     
   }   
    
   return cm;
}

#ifdef	__cplusplus
}
#endif

#endif	/* PDM_PART_COARSE_MESH_PRIV_H */
