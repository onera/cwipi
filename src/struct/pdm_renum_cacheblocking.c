/*============================================================================
 * Hilbert encoding for 2D or 3D coordinates.
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_part.h"
#include "pdm_renum_cacheblocking.h"
#include "pdm_part_graph.h"
#include "pdm_part_renum.h"
#include "pdm_printf.h"
#include "pdm_error.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

#define _MIN(a,b)   ((a) < (b) ?  (a) : (b))  /* Minimum of a et b */

#define _MAX(a,b)   ((a) > (b) ?  (a) : (b))  /* Maximum of a et b */

/*=============================================================================
 * Static global variables
 *============================================================================*/


/*============================================================================
 * Private function definitions
 *============================================================================*/

/**
 *
 * \brief Compute Bandwidth of a graph 
 *
 * \param [in,out]  node_num            The number of nodes.
 * \param [in,out]  adj_row[ADJ_NUM]    Information about row I is stored in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ
 * \param [in,out]  adj                 The adjacency structure. For each row, it contains the column indices of the nonzero entries.
 * \param [out]     ADJ_BANDWIDTH,      The bandwidth of the adjacency matrix.
 */


/*============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief Perform a cells renumbering with cache blocking (Synchrone)
 *
 * \param [in]  part ppart structure
 */

void 
PDM_renum_cacheblocking
(
 _part_t     *part,
int           split_method,  
int           nCellPerCacheWanted, 
int           isAsynchrone,  
int           isVectorisation 
)
{
  // @ Eric :  Deplacer depuis pdm_part_renum ???
  // Il faut deplcer et rendre public : _order_faceCell, renum_cell, renum_face ...
      
  /* Get nFac and nCel */
  const int nCell = part->nCell;
  const int nFace = part->nFace;
  
  /** Allocate reoerdering/permutation array **/
  int *CellOrder = (int *) malloc (sizeof(int) * nCell);
  int *FaceOrder = (int *) malloc (sizeof(int) * nFace);
  
  /*
   * I/ Coloring the cells with classic graph library
   */
  int *CellCellIdx = NULL;
  int *CellCell    = NULL;   

  /* Compute graph associate to mesh */
  PDM_compute_graph_from_face_cell(part,
                                   (int **) &CellCellIdx,
                                   (int **) &CellCell); 
  
  /* 
   * Determine the optimal size for cache blocking 
   *   -> The user specify the desired number of cells on each subdomain
   *   -> An other idea is to specify directly the number of block he want
   */
  int nBlkCacheWanted;
  if(nCellPerCacheWanted == 0){nBlkCacheWanted = 1;}
  else                        {nBlkCacheWanted = PDM_MAX(nCellPerCacheWanted,1);}
  
  /* Split the graph */
  PDM_split_graph(split_method,
                  nBlkCacheWanted,
                  part,
                  CellCellIdx, 
                  CellCell, 
                  (int *) NULL,
                  (int *) NULL,
                  (int **) &part->cellColor);
  
  /*
   * II/ Create a proper cells order :
   *       -> [ [SubDom1] , [SubDom2], ..., [SubDomN] ]
   */
  
  /* Allocate */
  int *nCellPerCache    = (int *) malloc( nBlkCacheWanted      * sizeof(int));
  int *nCellPerCacheBeg = (int *) malloc( nBlkCacheWanted      * sizeof(int));
  int* partCellIdx      = (int *) malloc((nBlkCacheWanted + 1) * sizeof(int));
  
  /* II-Bis - Asynchrone : We change only number of each sub-domain 
   * in order to have interior sub-domain first and exterior after 
   */
  if(isAsynchrone == 1)
  {
    /* Allocate */
    int *flagHaveBnd     = (int *) malloc(nBlkCacheWanted * sizeof(int));
    int *syncToAsynColor = (int *) malloc(nBlkCacheWanted * sizeof(int));
    
    /* All sub-domain is initialize */
    for(int i = 0; i < nBlkCacheWanted; i++) {
      flagHaveBnd[i] = 0;
    }
    
    /* 
     * Loop on faces : if at least one is border, flagHaveBnd is set to 1 
     */
    for (int i = 0; i < part->nFace; i++) {
      int iCell1 = part->faceCell[2*i    ];
      int iCell2 = part->faceCell[2*i + 1];
      int color1 = part->cellColor[iCell1-1];

      /* Solution 1 */        
      // if(iCell2 > 0 ){
      //   int color2 = part->cellColor[iCell2-1];
      //   int color  = _PDM_part_MIN(color1, color2);
      // }
      // else
      // {
      //   flagHaveBnd[color1] = 1;
      // }
      
      /* Solution 2 */     
      if(iCell2 == 0){flagHaveBnd[color1] = 1;}
      
    }
    
    /*
     * Built the syncToAsynColor array
     */
    int nSdomB = 0;                     /* Number of   Blocking domain */
    int nSdomU = 0;                     /* Number of Unblocking domain */
    int nextColorB = nBlkCacheWanted-1; /* Next Blocking   color       */
    int nextColorU = 0;                 /* Next UnBlocking color       */ 
    
    for(int i = 0; i < nBlkCacheWanted; i++) {
      if(flagHaveBnd[i] == 0)
      {
        syncToAsynColor[i] = nextColorU;
        nextColorU++;
        nSdomU++;
      }
      else
      {
        syncToAsynColor[i] = nextColorB;
        nextColorB--;
        nSdomB++;
      }
    }
    
    /* Dramatic verbose */
    if(0 == 1)
    {
      printf(" ----------- : %i \n", nBlkCacheWanted);
      for (int i = 0; i < nBlkCacheWanted; i++){
        printf("SyncToAsynColor[%i] = %i\n", i, syncToAsynColor[i]);
      }
    }
    
    /*
     * Apply the computed array on cellPart 
     */
    for(int i = 0; i < part->nCell; i++) {
      part->cellColor[i] = syncToAsynColor[part->cellColor[i]];
    }
    
    
    /* Free memory */
    free(flagHaveBnd);
    free(syncToAsynColor);
    
  } /* End asynchrone */

  
  /* Init to Zero */
  for (int i = 0; i < nBlkCacheWanted; i++){
    nCellPerCache[i] = 0;
  }
  
  /* First loop to identify the size of each color */
  for (int i = 0; i < nCell; i++){
    int color = part->cellColor[i]; // A color is a number of partition (output of Metis or Scotch)
    nCellPerCache[color]++;
  }

  /* Compute the index of each Sub-Domain **/
  partCellIdx[0] = 0;
  for (int i = 0; i < nBlkCacheWanted; i++){
    partCellIdx[i + 1] = partCellIdx[i] + nCellPerCache[i];
  }
  

  if(isVectorisation == 0){
    /** First method - Standard **/
    
    /* Reset array */
    for (int i = 0; i < nBlkCacheWanted; i++){
      nCellPerCache[i] = 0;
    }
    
    for (int i = 0; i < nCell; i++){
      int color  = part->cellColor[i]; 
      int idx    = partCellIdx[color] + nCellPerCache[color];
      
      CellOrder[idx] = i;
      nCellPerCache[color]++;
    }
  }
  else
  {
    /** Second method - Vectorisation on cell on current domain **/
    
    /* Store the current Index to add a vectorisable cell */
    for (int i = 0; i < nBlkCacheWanted; i++){
      nCellPerCacheBeg[i] = 0;
    }
    
    /* Loop on cells : 
     *     -> If current cell is adjacent to another color : add to end 
     *     -> Else add to begin
     * In fact this function sort the interior cell and exterior cell leads to the following numbering : 
     *         [ [ Interior Cells] , [ Exterior Cells]]
     */
    for (int i = 0; i < nCell; i++){
      
      /* Get color and prepare flag */
      int color = part->cellColor[i];
      int flag  = -1;
      for(int j = CellCellIdx[i]; j < CellCellIdx[i+1]; j++){
        int iCell = CellCell[j];
        if(part->cellColor[iCell] != color){flag = 1;} // Alors on est au bord d'un sous domaine !
      }
      
      if(flag == -1){  // Cell is interior : add to begin
        
        int idx = partCellIdx[color] + nCellPerCacheBeg[color];
        
        CellOrder[idx] = i;
        nCellPerCacheBeg[color]++;
        
      }
      else{ // Cell is exterior : add to end
        
        int idx = partCellIdx[color] + nCellPerCache[color]-1;
        
        CellOrder[idx] = i;   
        nCellPerCache[color]--;
        
      }
      
      /* Panic verbose */
      // if(0 == 1){
      //   PDM_printf("Begin : %i %i %i \n", color, idx, nCellPerCacheBeg[color]);
      // }
      
    }
  }
  
  /* 
   * Verbose 
   */
  if(0 == 1)
  {
    printf(" ----------- : %i \n", nBlkCacheWanted);
    for (int i = 0; i < nCell; i++){
      printf("~> %i\n", i);
      printf(" =====> CellOrder[%i]  = %i\n", i, CellOrder[i]);
    }
  }
  
  /*
   * We apply renumbering here because we need it for faces renumbering 
   */
  PDM_part_reorder_cell(part, CellOrder);
  
  /*
   * III/ Create a proper faces order :
   *       -> [ [SubDom1Int/SubDom1Ext] , [SubDom2Int/SubDom2Ext], ..., [SubDomNInt/SubDomNExt],  ]
   *      Renumbering the faces according to the sub block new ordering 
   */

  /* Allocate */    
  int *nFacePerCache    = (int *) malloc(nBlkCacheWanted * sizeof(int));
  int *nFaceBndPerCache = (int *) malloc(nBlkCacheWanted * sizeof(int));
  
  /* Init array */
  for (int i = 0; i < nBlkCacheWanted; i++){
    nFacePerCache[i]    = 0;
    nFaceBndPerCache[i] = 0;
  }
  
  /* 
   * First pass : 
   *      -> the face is associate with the subdomain that have the lowest color number 
   */
  for (int i = 0; i < part->nFace; i++) {
    int iCell1 = part->faceCell[2*i    ];
    int iCell2 = part->faceCell[2*i + 1];
    
    int color1 = part->cellColor[iCell1-1];
    if(iCell2 > 0 ){
      int color2 = part->cellColor[iCell2-1];
      int color  = PDM_MIN(color1, color2);
      nFacePerCache[color]++;
    }
    else
    {
      nFaceBndPerCache[color1]++;
    }
    
    /* Dramatic test */
    if(iCell1 < 0 ){ 
      printf("PPART internal error \n");
      exit(1); 
    }
  }
   
  /* Second pass : 
   *      -> the face is associate with the subdomain that have the lowest color number 
   */
  
  /* Allocate */
  int* partFaceIdx    = (int *) malloc((nBlkCacheWanted + 1) * sizeof(int));
  int* partFaceBndIdx = (int *) malloc((nBlkCacheWanted + 1) * sizeof(int));
  
  /* Initialise */
  partFaceIdx[0]    = 0;
  partFaceBndIdx[0] = 0;
  for (int i = 0; i < nBlkCacheWanted; i++){
    partFaceIdx[i + 1]    = partFaceIdx[i]    + nFacePerCache[i];
    partFaceBndIdx[i + 1] = partFaceBndIdx[i] + nFaceBndPerCache[i];
  }
  
  for (int i = 0; i < nBlkCacheWanted+1; i++){
    partFaceBndIdx[i] = partFaceBndIdx[i] + partFaceIdx[nBlkCacheWanted];
  }
  
  /* 
   * Verbose 
   */
  if(0 == 1)
  {
    printf(" ----------- : %i \n", nBlkCacheWanted);
    for (int i = 0; i < nBlkCacheWanted; i++){
      printf("~> %i\n", i);
      printf(" =====> partFaceIdx    %i\n", partFaceIdx[i + 1]);
      printf(" =====> partFaceBndIdx %i\n", partFaceBndIdx[i + 1]);
      
      printf(" =====> nFacePerCache    %i\n", nFacePerCache[i]);
      printf(" =====> nFaceBndPerCache %i\n", nFaceBndPerCache[i]);
    }
  }
      
  /* Reset Idx */
  for (int i = 0; i < nBlkCacheWanted; i++){
    nFacePerCache[i]    = 0;
    nFaceBndPerCache[i] = 0;
  }

  /* Determine for each subDomain the associate faces */    
  for (int i = 0; i < part->nFace; i++) {
    int iCell1 = part->faceCell[2*i    ];
    int iCell2 = part->faceCell[2*i + 1];
    
    int color1 = part->cellColor[iCell1-1];
    
    if(iCell2 > 0 ){
      int color2 = part->cellColor[iCell2-1];
      int color  = PDM_MIN(color1, color2);
      int idx    = partFaceIdx[color]    + nFacePerCache[color];
      
      FaceOrder[idx] = i;
      nFacePerCache[color]++;
    }
    else
    {
      int color  = color1;
      int idxBnd = partFaceBndIdx[color] + nFaceBndPerCache[color];
      
      FaceOrder[idxBnd] = i;
      nFaceBndPerCache[color]++;
    }
    
    if(iCell1 < 0 ){ 
      printf("PPART internal error \n");
      exit(1); 
    }
    
  } /* End file FaceOrder */
  
  /* Dramatic verbose */
  if(0 == 1)
  {
    for (int i = 0; i < part->nFace; i++) {
      printf("FaceOrder[%i] = %i \n", i, FaceOrder[i]);
    }
  }
  
  /*
   * Apply renumbering
   */
  PDM_part_reorder_face(part, FaceOrder);
  
  /*
   * Save in array the color of each faces 
   */
  part->faceColor = (int *) malloc (sizeof(int) * nFace);
  
  for (int i = 0; i < nBlkCacheWanted; i++){
    for (int iface = partFaceIdx[i]; iface < partFaceIdx[i+1]; iface++){
      part->faceColor[iface] = i;
    }
    for (int iface = partFaceBndIdx[i]; iface < partFaceBndIdx[i+1]; iface++){
      part->faceColor[iface] = i;
    }
  }
  
      
  /*
   * IV/ Create a proper faces order for vectorisation :
   *         -> The idea is to found a pool of faces that not pointed same cells : "independant faces"
   *         Rmk : Boundary faces are not vecotrised
   */
  for (int i = 0; i < part->nFace; i++) {
    FaceOrder[i] = i;
  }
  
  /* Allocate */
  int *flagCell = (int *) malloc (sizeof(int) * nCell);
  int *flagFace = (int *) malloc (sizeof(int) * nFace);
  
  /* Core loop of reordenencing */
  for (int i = 0; i < nBlkCacheWanted; i++){
    
    /* Get local size of SubDomain */
    int nFacSub = nFacePerCache[i];
    // int nCelSub = nCellPerCache[i];
    
    /* Reset for one pacquet all flag */
    for (int iface = partFaceIdx[i]; iface < partFaceIdx[i+1]; iface++){
      flagFace[iface] = -1;
    }
    
    int nFacTreated = 0;
    // printf("nFacSub : %i \n", nFacSub);
    while(nFacTreated != nFacSub){
      
      /* Reset for one pacquet all flag */
      // for (int iCel = 0; iCel < nCell; iCel++){
      //   flagCell[iCel] = -1;
      // }
      // TEST : A remplacer par : (Beacoup moins couteux)
      for (int iface = partFaceIdx[i]; iface < partFaceIdx[i+1]; iface++){
        int iCell1 = part->faceCell[2*iface  ];
        int iCell2 = part->faceCell[2*iface+1];
        flagCell[iCell1-1] = -1;
        flagCell[iCell2-1] = -1;
      }
      
      /* Dramatic verbose */
      if(0 == 1){
        printf("nFacTreated/nFacSub : %i -> %i \n", nFacTreated, nFacSub);
      }
      
      /* Loop on face */
      for (int iface = partFaceIdx[i]; iface < partFaceIdx[i+1]; iface++){
        
        if(flagFace[iface] == -1){
          int iCell1 = part->faceCell[2*iface  ];
          int iCell2 = part->faceCell[2*iface+1];
          
          int t1 = flagCell[iCell1-1];
          int t2 = flagCell[iCell2-1];
          
          // printf("t1/t2 : %i/%i \n", t1, t2);
          if( (t1 == -1) && (t2 == -1)){
            flagCell[iCell1-1] = 1;
            flagCell[iCell2-1] = 1;
            
            int bFac = partFaceIdx[i];
            
            FaceOrder[bFac+nFacTreated] = iface;
            
            flagFace[iface] = 1;
            
            nFacTreated++;
          }
        } /* End Face loop */
      } /* End While */
    } /* End SubDom loop */
     
  }
  
  /*
   * Apply renumbering -> Attention au tableau en plus à reordencer !
   */
  PDM_part_reorder_face(part, FaceOrder);
  
  /* Dramatic verbose */
  if(0 == 1)
  {
    for (int i = 0; i < part->nFace; i++) {
      printf("FaceOrderVect[%i] = %i \n", i, FaceOrder[i]);
    }
  }
  
  /* Free memory */
  free(CellOrder);
  free(FaceOrder);
  free(CellCellIdx);
  free(CellCell);
  free(nCellPerCacheBeg);
  free(nCellPerCache);
  free(nFacePerCache);
  free(nFaceBndPerCache);
  free(partFaceIdx);
  free(partFaceBndIdx);
  free(flagCell);
  free(flagFace);
  
  

}

#ifdef  __cplusplus
}
#endif
