
/*============================================================================
 * Parallel partitioning
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  System headers
 *----------------------------------------------------------------------------*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include <string.h>

#include "pdm.h"
#include "pdm_mpi.h"
#include "pdm_priv.h"
#include "pdm_config.h"

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm_part.h"
#include "pdm_part_priv.h"
#include "pdm_timer.h"

#include "pdm_part_geom.h"
#include "pdm_part_renum.h"
#include "pdm_hilbert.h"
#include "pdm_handles.h"
#include "pdm_geom_elem.h"
#include "pdm_part_graph.h"
#include "pdm_sort.h"
#include "pdm_order.h"
#include "pdm_cuthill.h"
#include "pdm_printf.h"
#include "pdm_error.h"

#ifdef PDM_IN_PDMA
#include "pdm_renum_cacheblocking.h"
#endif

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Fortran function header
 *============================================================================*/


/*============================================================================
 * Local macro definitions
 *============================================================================*/


/*============================================================================
 * Type definitions
 *============================================================================*/

/**
 * \struct PDM_writer_fmt_t
 * \brief  Writer format
 *
 */

typedef struct _renum_method_t {

  char                  *name; /*!< Name of method */
  PDM_part_renum_fct_t   fct;  /*!< Renumbering function */

} _renum_method_t;

/*============================================================================
 * Global variable
 *============================================================================*/

/**
 * Storage of face renumbering methods   
 */

static PDM_Handles_t *face_methods = NULL; 

/**
 * Storage of cell renumbering methods   
 */

static PDM_Handles_t *cell_methods = NULL; 

/*============================================================================
 * Private function definitions
 *============================================================================*/


/**
 *
 * \brief Random order 
 *
 * \param [in]   nElt    Number of elements
 * \param [out]  order   Random order
 *
 */

static void 
_random_order
(
const int nElt,
      int *order
)
{
  int *tmpArray = (int *) malloc (sizeof(int) * nElt);
  
  time_t _seed = time(NULL);
  srand(( unsigned int) _seed);
  for (int i = 0; i < nElt; i++) {
    tmpArray[i] = rand()%nElt;
    order[i] = i;
  }
  
  PDM_sort_int (tmpArray, order, nElt);
}

/**
 * \brief Renumber face to cell connectivity 
 * 
 * \param [in]      nCell        Number of cells
 * \param [in]      nFace        Number of faces
 * \param [in]      cellFaceIdx  Cell face connectivity Index
 * \param [in]      cellFace     Cell face connectivity
 * \param [in, out] faceCell     Face cell connectivity
 *
 */

static void 
_renum_faceCell 
(
const int  nCell,
const int  nFace,
const int *cellFaceIdx, 
const int *cellFace, 
      int *faceCell 
)
{

  for (int i = 0; i < 2 * nFace; i++) {
    faceCell[i] = 0;
  }

  for (int i = 0; i < nCell; ++i) {
    for (int j = cellFaceIdx[i]; j < cellFaceIdx[i+1]; j++) {
      int idx = 2 * (PDM_ABS(cellFace[j])-1);
      if (faceCell[idx] == 0) { 
        faceCell[idx] = i + 1;
        if (cellFace[j] < 0) {
          faceCell[idx] = -faceCell[idx];
        }
      }
      else { 
        faceCell[idx + 1] = i + 1;
        if (cellFace[j] < 0) {
          faceCell[idx+1] = -faceCell[idx+1];
        }
      }
    }
  }
 
}

/**
 * \brief Order faceCell array
 * 
 * \param [in]      sizeArray       Number of elements
 * \param [in]      newToOldOrder   New order (size = \ref nElt
 * \param [in, out] faceCell        Array to order
 *
 */

static void 
_order_faceCell 
(
int          nFace,
int         *newToOldOrder,
int         *faceCell 
)
{
  int *oldFaceCell = 
          (int *) malloc (nFace * 2 * sizeof(int));
  for(int i = 0; i < nFace * 2; ++i) {
    oldFaceCell [i] = faceCell [i];
  }
  
  for(int i = 0; i < nFace; ++i) {
    faceCell[i*2+0] = oldFaceCell[newToOldOrder[i]*2+0];
    faceCell[i*2+1] = oldFaceCell[newToOldOrder[i]*2+1];
  }
  
  free(oldFaceCell);
}


/**
 * \brief Order an array
 * 
 * \param [in]      sizeArray       Number of elements
 * \param [in]      newToOldOrder        New order (size = \ref nElt
 * \param [in, out] Array         	Array to renumber
 *
 */

static void 
_renum_array 
(
const int  sizeArray,
const int *olToNewOrder,
int       *array
)
{
  int *oldArray = (int *) malloc (sizeof(int) * sizeArray);
  
  for (int i = 0; i < sizeArray; ++i) {
    oldArray[i] = array[i];
  }
  
  for (int i = 0; i < sizeArray; ++i) {
    array[i] = olToNewOrder[oldArray[i]-1] + 1;
  }
  
  free(oldArray);
}


/**
 * \brief Renumber connectivities 
 * 
 * \param [in]      nElt            Number of elements
 * \param [in]      newToOldOrder        New order (size = \ref nElt
 * \param [in, out] connectivityIdx	Connectivity index
 * \param [in, out] connectivities	Element connectivities
 *
 */
 
static void 
_renum_connectivities 
(
const int nElt,
const int *newToOldOrder,
int       *connectivityIdx,
int       *connectivities
)
{
  
  int *oldConnectivities = (int *) malloc (connectivityIdx[nElt] * sizeof(int));
  
  for (int i = 0; i < connectivityIdx[nElt]; ++i) { 
    oldConnectivities[i] = connectivities[i];
  }
    
  int *oldConnectivityIdx = (int *) malloc ((nElt + 1) * sizeof(int));
  
  for (int i = 0; i < nElt+1; ++i) { 
    oldConnectivityIdx[i] = connectivityIdx[i];
  }

  for (int elem = 0; elem < nElt; ++elem) {
    int nbSSElem = oldConnectivityIdx[newToOldOrder[elem] + 1] - oldConnectivityIdx[newToOldOrder[elem]];
    connectivityIdx[elem+1] = nbSSElem;
  }
  
  connectivityIdx[0] = 0;
  for (int elem = 1; elem < nElt+1; ++elem) {
    connectivityIdx[elem] += connectivityIdx[elem-1];
  }
    
  for (int elem = 0; elem < nElt; ++elem) {
    int nbSSElem = oldConnectivityIdx[newToOldOrder[elem] + 1] - oldConnectivityIdx[newToOldOrder[elem]];
    
    for (int ssElem = 0; ssElem < nbSSElem; ++ssElem) {
      connectivities[connectivityIdx[elem] + ssElem] = oldConnectivities[oldConnectivityIdx[newToOldOrder[elem]]+ssElem];
    }
  }
  
  free(oldConnectivities);
  free(oldConnectivityIdx);
  
}


/**
 *
 * \brief Perform a cells renumbering from a Hilbert curve
 *
 * \param [in]  part      Current PPART structure
 * \param [out] cellCenter Cell center
 *
 */

static void 
_compute_cellCenter 
(
_part_t *part,
double  *cellCenter
)
{
  if (part->nCell <= 0) {
    return;
  }
  
  int isPoly3D = (part->faceVtxIdx[1] > 2);

  if (isPoly3D) {
    const int isOriented = 0;
    double *volume = (double *) malloc (part->nCell * sizeof(double));
    int isDegenerated;  

    if(1 == 0){
      PDM_geom_elem_polyhedra_properties (isOriented,
                                          part->nCell,
                                          part->nFace,
                                          part->faceVtxIdx,
                                          part->faceVtx,
                                          part->cellFaceIdx,
                                          part->cellFace,
                                          part->nVtx,
                                          part->vtx,
                                          volume,
                                          cellCenter,
                                          NULL,
                                          &isDegenerated);
    }
    else /*Trash patch */
    {
      /* Allocate */
      double *cellWeight = (double *) malloc (part->nCell * sizeof(double));
      
      /* Nulliffy cellCenterArray */
      for(int iCell = 0; iCell < part->nCell; iCell++) {
        cellCenter[3*iCell  ] = 0.;
        cellCenter[3*iCell+1] = 0.;
        cellCenter[3*iCell+2] = 0.;
        cellWeight[iCell]     = 0.;
      }
      
      /* Compute */
      for(int iCell = 0; iCell < part->nCell; iCell++) {
        
        /* Cellule composé de nFace */
        int aFac = part->cellFaceIdx[iCell];
        int nFac = part->cellFaceIdx[iCell+1] - aFac;
  
        for(int iFac = 0; iFac < nFac; iFac++) {
          
          /* Face composé de nVtx */
          int lFac = PDM_ABS(part->cellFace[aFac + iFac]) - 1;
  
          int aVtx = part->faceVtxIdx[lFac];
          int nVtx = part->faceVtxIdx[lFac+1] - aVtx;
          
          for(int iVtx = 0; iVtx < nVtx; iVtx++) {
  
            /* Face composé de nVtx */
            int lVtx = part->faceVtx[aVtx + iVtx] - 1;

            /* Add to current cell and stack weight */
            cellCenter[3*iCell  ] += part->vtx[3*lVtx  ];
            cellCenter[3*iCell+1] += part->vtx[3*lVtx+1];
            cellCenter[3*iCell+2] += part->vtx[3*lVtx+2];
  
            cellWeight[iCell] += 1.;
          }
        }
      }   
  
      /* Nulliffy cellCenterArray */
      for(int iCell = 0; iCell < part->nCell; iCell++) {
        cellCenter[3*iCell  ] = cellCenter[3*iCell  ]/cellWeight[iCell];
        cellCenter[3*iCell+1] = cellCenter[3*iCell+1]/cellWeight[iCell];
        cellCenter[3*iCell+2] = cellCenter[3*iCell+2]/cellWeight[iCell];
      }
  
      /* Verbose */
      if(0 == 1){
        for(int iCell = 0; iCell < part->nCell; iCell++) {
          PDM_printf("cellCenter (X,Y,Z) : %f - %f - %f \n", cellCenter[3*iCell  ], cellCenter[3*iCell+1], cellCenter[3*iCell+2]);
          PDM_printf("cellWeight         : %f  \n", cellWeight[iCell  ]);
        }
      }

      /* Free */
      free(cellWeight);

    }

    /* Free */
    free (volume);
  }  
  else {   /* isPoly3D */
    double *surfaceVector = (double * ) malloc( sizeof(double) * 3 * part->nCell); 
    int isDegenerated;  
    
    int *connectivity = (int *) malloc (part->cellFaceIdx[part->nCell] 
                        * sizeof(int));

    
    int idx = 0;    
    for (int icell = 0; icell < part->nCell; icell++) {
      int faceCurr = PDM_ABS(part->cellFace[part->cellFaceIdx[icell]]) - 1;
      int deb  = part->faceVtx[part->faceVtxIdx[faceCurr]];
      int next = part->faceVtx[part->faceVtxIdx[faceCurr]+1];
      connectivity[idx++] = deb;
      connectivity[idx++] = next;
      while (next != deb) {
        for (int j = part->cellFaceIdx[icell]; j < part->cellFaceIdx[icell+1]; ++j) {
          int face = PDM_ABS(part->cellFace[j]) - 1;
          if (faceCurr != face) {
            int s1 = part->faceVtx[part->faceVtxIdx[face]  ];
            int s2 = part->faceVtx[part->faceVtxIdx[face]+1];
						
            if ((s1 == next) || (s2 == next)) { 
              if (s1 == next) {
                next = s2;
              }
              else if (s2 == next) {
                next = s1;
              }
              if (next != deb) {
                connectivity[idx++] = next;
              }
              faceCurr = face;
              break;
            }              
          }
        }
      }
    }

		assert (idx == part->cellFaceIdx[part->nCell]);
		
    PDM_geom_elem_polygon_properties (part->nCell,   
                                      part->cellFaceIdx, 
                                      connectivity,
                                      part->vtx,
                                      surfaceVector,
                                      cellCenter,
                                      NULL,
                                      &isDegenerated);

    free (surfaceVector);
    free (connectivity);

  }

}

/**
 *
 * \brief Perform a cells renumbering from a Hilbert curve
 *
 * \param [in,out]  ppart    Current PPART structure
 *
 */

static void 
_renum_cells_hilbert 
(
_PDM_part_t* ppart
)
{
  
  for(int ipart = 0; ipart < ppart->nPart; ++ipart) {
    _part_t *part = ppart->meshParts[ipart];       
    double *cellCenter = 
        (double *) malloc (part->nCell * 3 * sizeof(double ));
    PDM_hilbert_code_t *hilbertCodes = 
        (PDM_hilbert_code_t *) malloc (part->nCell * sizeof(PDM_hilbert_code_t));
      
    /** Barycentre computation **/
    
    _compute_cellCenter (part, cellCenter);
        
    double extents[3 * 2]; 
    
    /** Get EXTENTS LOCAL **/
    
    PDM_hilbert_get_coord_extents_seq(3, part->nCell, cellCenter, extents);
    
    /** Hilbert Coordinates Computation **/
    
    PDM_hilbert_encode_coords(3, PDM_HILBERT_CS, extents, part->nCell, cellCenter, hilbertCodes);
    
    /** CHECK H_CODES **/
    
    free(cellCenter);
    
    int *newToOldOrder = (int *) malloc (part->nCell * sizeof(int));
    for(int i = 0; i < part->nCell; ++i) {
      newToOldOrder [i] = i;
    }
      
    PDM_sort_double (hilbertCodes, newToOldOrder, part->nCell);
	  
    PDM_part_reorder_cell(part, newToOldOrder);
          
    free (hilbertCodes);
    free (newToOldOrder);
    
  }
}

/**
 *
 * \brief Perform a cells renumbering reverse CutHill Mac-Kee
 *
 * \param [in,out]  ppart    Current PPART structure
 *
 */

static void
_renum_cells_cuthill
(
_PDM_part_t* ppart
)
{

  /** Loop over all part of the current process **/
  for(int ipart = 0; ipart < ppart->nPart; ++ipart) {
    /** Get current part id **/
    _part_t *part = ppart->meshParts[ipart];
    const int nCell = part->nCell;
    
    /** Allocate reoerdering/permutation array **/
    int *order = (int *) malloc (sizeof(int) * nCell);

    /** Verbose bandwidth **/
    // dualBandWidth = PDM_checkbandwidth(part);
    // PDM_printf("Bandwidth of graph before reordering : %d \n", dualBandWidth);

    /** Compute reordering **/
    PDM_cuthill_generate(part, order);
  
    /** Apply renumbering **/
    PDM_part_reorder_cell(part, order);

    /** Verbose bandwidth **/
    // dualBandWidth = PDM_checkbandwidth(part);
    // PDM_printf("Bandwidth of graph after reordering : %d \n", dualBandWidth);

    /** Free memory **/
    free(order);
  }
}


/**
 *
 * \brief Perform a cells renumbering from a Hilbert curve
 *
 * \param [in,out]  ppart    Current PPART structure
 *
 */

static void 
_renum_cells_random 
(
_PDM_part_t* ppart
)
{
  for(int ipart = 0; ipart < ppart->nPart; ++ipart) {
    _part_t *part = ppart->meshParts[ipart];       
    const int nCell = part->nCell;
    
    int *order = (int *) malloc (sizeof(int) * nCell);
    
    _random_order (nCell, order);
    
    PDM_part_reorder_cell(part, order);
      
    free (order);
  }
}


/**
 *
 * \brief Perform a face random renumbering 
 *
 * \param [in,out]  ppart    Current PPART structure
 *
 */

static void 
_renum_faces_random 
(
_PDM_part_t* ppart
)
{
  for(int ipart = 0; ipart < ppart->nPart; ++ipart) {
    _part_t *part = ppart->meshParts[ipart];       
    const int nFace = part->nFace;
    
    int *order = (int *) malloc (sizeof(int) * nFace);
    
    _random_order (nFace, order);
    
    PDM_part_reorder_face(part, order);
      
    free (order);
  }
}

/**
 *
 * \brief Perform a face random renumbering
 *
 * \param [in,out]  ppart    Current PPART structure
 *
 */

static void
_renum_faces_lexicographic
(
_PDM_part_t* ppart
)
{
  for(int ipart = 0; ipart < ppart->nPart; ++ipart) {
    _part_t *part = ppart->meshParts[ipart];
    const int nFace = part->nFace;

    int *order = (int *) malloc (sizeof(int) * nFace);

    /** Build a pre-array face cell ordered */
    int *faceCellTmp = (int *) malloc(2*nFace * sizeof(int)); 

    for(int i = 0; i < nFace; i++) {
       int iL = PDM_ABS (part->faceCell[2*i  ]);
       int iR = PDM_ABS (part->faceCell[2*i+1]);
       if(iL < iR )
       {
          faceCellTmp[2*i  ] = iR;
          faceCellTmp[2*i+1] = iL;
       }
       else
       {
          faceCellTmp[2*i  ] = iL;
          faceCellTmp[2*i+1] = iR;
       }
    }
  
    /** Reorder lexicographicly the array */
    PDM_order_lnum_s (faceCellTmp, 2, order, nFace);

    /** Update face array with the new array **/
    PDM_part_reorder_face(part, order);

    /** Free memory **/
    free (order);
    free (faceCellTmp);
  }
}


/*=============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief Purge renumbering methods 
 *
 */

void 
PDM_part_renum_purge
(
 void
)        
{
  if (face_methods != NULL) {

    const int *index =  PDM_Handles_idx_get (face_methods);
    int n_methods = PDM_Handles_n_get (face_methods);
    
    while (n_methods > 0) {
      int idx = index[0];
      _renum_method_t *method_ptr = 
              (_renum_method_t *) PDM_Handles_get (face_methods, idx);
      free (method_ptr->name);
      PDM_Handles_handle_free (face_methods, idx, PDM_TRUE);
      n_methods = PDM_Handles_n_get (face_methods);
    }

    face_methods = PDM_Handles_free (face_methods);
    
  }

  if (cell_methods != NULL) {

    const int *index =  PDM_Handles_idx_get (cell_methods);
    int n_methods = PDM_Handles_n_get (cell_methods);
    
    while (n_methods > 0) {
      int idx = index[0];
      _renum_method_t *method_ptr = 
              (_renum_method_t *) PDM_Handles_get (cell_methods, idx);
      free (method_ptr->name);
      PDM_Handles_handle_free (cell_methods, idx, PDM_TRUE);
      n_methods = PDM_Handles_n_get (cell_methods);
    }

    cell_methods = PDM_Handles_free (cell_methods);
    
  }
}


/**
 *
 * \brief Get index of a renumbering cell method
 * 
 * \param [in]  name   Name of the method
 * 
 * \return Index (-1 if not found)
 */

int 
PDM_part_renum_cell_method_idx_get
(
const char *name
)        
{
  if (cell_methods == NULL) {
    PDM_part_renum_load_local();
  }
  int idx = -1;
  
  if (cell_methods != NULL) {
    int n_methods = PDM_Handles_n_get (cell_methods);
    const int *index =  PDM_Handles_idx_get (cell_methods);

    for (int i = 0; i < n_methods; i++) {
      _renum_method_t *method_ptr = 
              (_renum_method_t *) PDM_Handles_get (cell_methods, index[i]);
      if (!strcmp(method_ptr->name, name)) {
        idx = index[i];
        break;
      }      
    }
  }
  return idx;

}

/**
 *
 * \brief Get index of a renumbering face method
 * 
 * \param [in]  name   Name of the method
 * 
 * \return Index (-1 if not found)
 */

int 
PDM_part_renum_face_method_idx_get
(
const char *name
)        
{
  if (face_methods == NULL) {
    PDM_part_renum_load_local();
  }
  int idx = -1;
  int n_methods = PDM_Handles_n_get (face_methods);
  const int *index =  PDM_Handles_idx_get (face_methods);

  for (int i = 0; i < n_methods; i++) {
    _renum_method_t *method_ptr = 
            (_renum_method_t *) PDM_Handles_get (face_methods, index[i]);
    if (!strcmp(method_ptr->name, name)) {
      idx = index[i];
      break;
    }      
  }
  return idx;
}



/**
 *
 * \brief Get name of the cell renumbering method 
 * 
 * \param [in]  idx     Index of the method
 * 
 * \return Name of the method
 *
 */

const char * 
PDM_part_renum_cell_method_name_get
(
const int idx
)
{
  if (cell_methods == NULL) {
    PDM_part_renum_load_local();
  }

  int n_methods = PDM_Handles_n_get (cell_methods);

  if (idx >= n_methods) {
    return NULL;
  }

  const int *index =  PDM_Handles_idx_get (cell_methods);

  _renum_method_t *method_ptr = 
            (_renum_method_t *) PDM_Handles_get (cell_methods, index[idx]);

  return method_ptr->name;
}


/**
 *
 * \brief Get name of the face renumbering method 
 * 
 * \param [in]  idx     Index of the method
 * 
 * \return Name of the method (NULL otherwise)
 *
 */

const char * 
PDM_part_renum_face_method_name_get
(
const int idx
)
{
  if (face_methods == NULL) {
    PDM_part_renum_load_local();
  }

  int n_methods = PDM_Handles_n_get (face_methods);

  if (idx >= n_methods) {
    return NULL;
  }

  const int *index =  PDM_Handles_idx_get (face_methods);

  _renum_method_t *method_ptr = 
            (_renum_method_t *) PDM_Handles_get (face_methods, index[idx]);

  return method_ptr->name;
}



/**
 *
 * \brief Get the number of renumbering face methods 
 * 
 * \return Name of the method
 *
 */

int  
PDM_part_n_renum_cell_method_get
(
void 
)
{
  if (cell_methods == NULL) {
    PDM_part_renum_load_local();
  }
  
  return PDM_Handles_n_get (cell_methods);
  
}



/**
 *
 * \brief Get the number of renumbering face methods 
 * 
 * \return Name of the method
 *
 */

int  
PDM_part_n_renum_face_method_get
(
void 
)
{
  if (face_methods == NULL) {
    PDM_part_renum_load_local();
  }

  return PDM_Handles_n_get (face_methods);

}

/**
 *
 * \brief Add a new method for cell renumbering 
 *
 * \param [in]      name           Mesh entity to renumber
 * \param [in]      renum_fct      Renumbering function
 *
 */

int 
PDM_part_renum_cell_add
(
 const char                 *name,     /*!< Name          */ 
 const PDM_part_renum_fct_t  renum_fct /*!< Customize \ref PDM_part_renum_cell function for the format */             
)
{
  _renum_method_t *method_ptr = malloc (sizeof(_renum_method_t));
  
  int idx = PDM_Handles_store  (cell_methods, method_ptr);

  method_ptr->name = malloc (sizeof(char) * (strlen(name) + 1));
  strcpy (method_ptr->name, name);
  
  method_ptr->fct = renum_fct;

  return idx;
}

/**
 *
 * \brief Add a new method for face renumbering 
 *
 * \param [in]      name           Mesh entity to renumber
 * \param [in]      renum_fct      Renumbering function
 *
 */

int 
PDM_part_renum_face_add
(
 const char                 *name,     /*!< Name          */ 
 const PDM_part_renum_fct_t  renum_fct /*!< Customize \ref PDM_part_renum_face function for the format */             
)
{
  _renum_method_t *method_ptr = malloc (sizeof(_renum_method_t));
  
  int idx = PDM_Handles_store  (face_methods, method_ptr);

  method_ptr->name = malloc (sizeof(char) * (strlen(name) + 1));
  strcpy (method_ptr->name, name);
  
  method_ptr->fct = renum_fct;

  return idx;
}


/**
 *
 * \brief Load local renumbering methods 
 *
 */

void 
PDM_part_renum_load_local
(
void
)
{
  if (cell_methods == NULL)  {
    
    const int n_default_methods = 4;
    cell_methods = PDM_Handles_create (n_default_methods);
    
    PDM_part_renum_cell_add ("PDM_PART_RENUM_CELL_NONE", 
                             NULL);
    PDM_part_renum_cell_add ("PDM_PART_RENUM_CELL_RANDOM", 
                             _renum_cells_random);
    PDM_part_renum_cell_add ("PDM_PART_RENUM_CELL_HILBERT",
                             _renum_cells_hilbert);
    PDM_part_renum_cell_add ("PDM_PART_RENUM_CELL_CUTHILL",
                             _renum_cells_cuthill);
  }

  if (face_methods == NULL)  {
    const int n_default_methods = 3;
    face_methods = PDM_Handles_create (n_default_methods);

    PDM_part_renum_face_add ("PDM_PART_RENUM_FACE_NONE", 
                             NULL);
    PDM_part_renum_face_add ("PDM_PART_RENUM_CELL_RANDOM", 
                             _renum_faces_random);
    PDM_part_renum_face_add ("PDM_PART_RENUM_CELL_LEXICOGRAPHIC",
                             _renum_faces_lexicographic);
  }
  
}

  
/**
 *
 * \brief Perform cell renumbering
 *
 * \param [in,out]  part       part structure
 *
 */

void 
PDM_part_renum_cell
(
  _PDM_part_t   *ppart                
)        
{
  
  if (cell_methods == NULL)  {
    PDM_part_renum_load_local ();
  }
  
  PDM_part_renum_fct_t fct =
          (PDM_part_renum_fct_t) PDM_Handles_get (cell_methods, 
                                                  ppart->renum_cell_method);
  
  if (fct != NULL) {
    (fct) (ppart);
  }
  
}


//  case PDM_PART_RENUM_CELL_CACHEBLOCKING_SYNC :
//    _renum_cells_cacheblocking(ppart); 
//    break;
//  case PDM_PART_RENUM_CELL_CACHEBLOCKING_ASYNC :
//    _renum_cells_cacheblocking(ppart); 
//    break;
//}



/**
 *
 * \brief Perform mesh entities renumbering
 *
 * \param [in,out]  part       part structure
 *
 */

void 
PDM_part_renum_face
(
 _PDM_part_t           *ppart                
)
{
  
  if (cell_methods == NULL)  {
    PDM_part_renum_load_local ();
  }
  
  PDM_part_renum_fct_t fct =
          (PDM_part_renum_fct_t) PDM_Handles_get (cell_methods, 
                                                  ppart->renum_face_method);
  
  if (fct != NULL) {
    (fct) (ppart);
  }
}

/**
 *
 * \brief Perform cells renumbering from a new order
 *
 * \param [in,out]  part        Current partition
 * \param [in]      newToOldOrder    NewOrder
 *
 */
void 
PDM_part_reorder_cell
(
 _part_t *part, 
 int     *newToOldOrder
)
{
  /*
   * Cell Renumbering
   */
   
  _renum_connectivities (part->nCell,
                         newToOldOrder,
                         part->cellFaceIdx, 
                         part->cellFace); 
  
  if (part->cellTag != NULL) {
    PDM_order_array (part->nCell,
                     sizeof(int),
                     newToOldOrder,
                     part->cellTag);
  }
  
  if (part->cellColor != NULL) {
    PDM_order_array (part->nCell,
                     sizeof(int),
                     newToOldOrder,
                     part->cellColor);
  }
   
  PDM_order_array (part->nCell,
                   sizeof(PDM_g_num_t),
                   newToOldOrder,
                   part->cellLNToGN); 
   
  _renum_faceCell (part->nCell,
                   part->nFace,
                   part->cellFaceIdx, 
                   part->cellFace, 
                   part->faceCell); 
   
}


/**
 *
 * \brief Perform faces renumbering from a new order
 *
 * \param [in,out]  part        Current partition
 * \param [in]      newToOldOrder    NewOrder
 *
 */
void 
PDM_part_reorder_face
(
_part_t *part, 
int     *newToOldOrder
)
{
  
  /** Renum FaceVtx / FaceVtxIdx **/

  _renum_connectivities (part->nFace, 
                         newToOldOrder, 
                         part->faceVtxIdx, 
                         part->faceVtx);

  /** CellFace **/
  
   int *oldToNewOrder = (int *) malloc (part->nFace * sizeof(int));
  
   for(int i = 0; i < part->nFace; i++) {
    oldToNewOrder[newToOldOrder[i]] = i;
   }
 
  _renum_array (part->cellFaceIdx[part->nCell], 
                oldToNewOrder,
                part->cellFace);
  
  free (oldToNewOrder);
    
  /** FaceTag **/
  if (part->faceTag != NULL) {
    PDM_order_array (part->nFace,
                     sizeof(int),
                     newToOldOrder,
                     part->faceTag); 
  }
  
  /** FaceColor **/
  if (part->faceColor != NULL) {
    PDM_order_array (part->nFace,
                     sizeof(int),
                     newToOldOrder,
                     part->faceColor); 
  }
    
   /** FaceLNToGN **/

  PDM_order_array (part->nFace,
                   sizeof(PDM_g_num_t),
                   newToOldOrder,
                   part->faceLNToGN); // OK
    
  /** FaceCell Face **/
  
  _order_faceCell (part->nFace, 
                   newToOldOrder,
                   part->faceCell);  

}

#ifdef __cplusplus
}
#endif /* __cplusplus */

