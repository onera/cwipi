
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


/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"
#include "pdm_fortran_to_c_string.h"
#include "pdm_priv.h"
#include "pdm_config.h"
#include "pdm_part.h"
#include "pdm_part_priv.h"
#include "pdm_timer.h"

#include "pdm_part_geom.h"
#include "pdm_part_renum.h"
#include "pdm_hilbert.h"
#include "pdm_handles.h"
#include "pdm_geom_elem.h"
#include "pdm_sort.h"
#include "pdm_cuthill.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_order.h"

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
 * \struct _renum_method_t
 * \brief Renumbering method
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

/**
 * Storage of vtx renumbering methods
 */

static PDM_Handles_t *edge_methods = NULL;

/**
 * Storage of vtx renumbering methods
 */

static PDM_Handles_t *vtx_methods = NULL;

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
 * \param [in]      n_cell        Number of cells
 * \param [in]      n_face        Number of faces
 * \param [in]      cell_face_idx  Cell face connectivity Index
 * \param [in]      cell_face     Cell face connectivity
 * \param [in, out] face_cell     Face cell connectivity
 *
 */

static void
_renum_faceCell
(
const int  n_cell,
const int  n_face,
      int *face_cell,
      int *new_to_old_order
)
{

  int *old_to_new_order = (int *) malloc (n_cell * sizeof(int));

  for(int i = 0; i < n_cell; i++) {
   old_to_new_order[new_to_old_order[i]] = i;
  }

  PDM_part_renum_array_face_cell (n_face,
                                  old_to_new_order,
                                  face_cell);
  free (old_to_new_order);

}

/**
 * \brief Order face_cell array
 *
 * \param [in]      sizeArray       Number of elements
 * \param [in]      new_to_old_order   New order (size = \ref nElt
 * \param [in, out] face_cell        Array to order
 *
 */

static void
_order_face_cell
(
int          n_face,
int         *new_to_old_order,
int         *face_cell
)
{
  int *oldface_cell = (int *) malloc (n_face * 2 * sizeof(int));
  for(int i = 0; i < n_face * 2; ++i) {
    oldface_cell [i] = face_cell [i];
  }

  for(int i = 0; i < n_face; ++i) {
    face_cell[i*2+0] = oldface_cell[new_to_old_order[i]*2+0];
    face_cell[i*2+1] = oldface_cell[new_to_old_order[i]*2+1];
  }

  free(oldface_cell);
}


/**
 * \brief Order an array
 *
 * \param [in]      sizeArray       Number of elements
 * \param [in]      new_to_old_order        New order (size = \ref nElt
 * \param [in, out] Array         	Array to renumber
 *
 */

void
PDM_part_renum_array
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
    int old_idx = PDM_ABS(oldArray[i]);
    int sign    = PDM_SIGN(oldArray[i]);
    array[i] = sign * (olToNewOrder[old_idx-1] + 1);
  }

  free(oldArray);
}


/**
 * \brief Order an array
 *
 * \param [in]      sizeArray       Number of elements
 * \param [in]      new_to_old_order   New order (size = \ref nElt
 * \param [in, out] Array           Array to renumber
 *
 */

void
PDM_part_renum_array_face_cell
(
const int  n_face,
const int *olToNewOrder,
int       *array
)
{
  int *oldArray = (int *) malloc (sizeof(int) * 2 * n_face);

  for (int i = 0; i < 2*n_face; ++i) {
    oldArray[i] = array[i];
  }

  for (int i = 0; i < 2*n_face; ++i) {
    if(PDM_ABS(oldArray[i]) == 0 ){
      array[i] = 0;
      // printf("[%i] BND \n", i);
    }
    else{
      int signe = (oldArray[i] < 0 ) ? -1 : 1 ;
      int iCel  = PDM_ABS(oldArray[i])-1;
      // array[i] = olToNewOrder[oldArray[i]-1] + 1;
      // printf("[%i] = %i/%i - %i  \n", i,signe, iCel, oldArray[i]);
      array[i] = signe * (olToNewOrder[iCel] + 1);
    }
  }

  free(oldArray);
}

/**
 * \brief Renumber connectivities
 *
 * \param [in]      nElt            Number of elements
 * \param [in]      new_to_old_order        New order (size = \ref nElt
 * \param [in, out] connectivityIdx	Connectivity index
 * \param [in, out] connectivities	Element connectivities
 *
 */

void
PDM_part_renum_connectivities
(
const int nElt,
const int *new_to_old_order,
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
    int nbSSElem = oldConnectivityIdx[new_to_old_order[elem] + 1] - oldConnectivityIdx[new_to_old_order[elem]];
    connectivityIdx[elem+1] = nbSSElem;
  }

  connectivityIdx[0] = 0;
  for (int elem = 1; elem < nElt+1; ++elem) {
    connectivityIdx[elem] += connectivityIdx[elem-1];
  }

  for (int elem = 0; elem < nElt; ++elem) {
    int nbSSElem = oldConnectivityIdx[new_to_old_order[elem] + 1] - oldConnectivityIdx[new_to_old_order[elem]];

    for (int ssElem = 0; ssElem < nbSSElem; ++ssElem) {
      connectivities[connectivityIdx[elem] + ssElem] = oldConnectivities[oldConnectivityIdx[new_to_old_order[elem]]+ssElem];
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
  if (part->n_cell <= 0) {
    return;
  }

  int isPoly3D = (part->face_vtx_idx[1] > 2);

  if (isPoly3D) {
    const int isOriented = 0;
    double *volume = (double *) malloc (part->n_cell * sizeof(double));
    int isDegenerated;

    if(1 == 0){
      PDM_geom_elem_polyhedra_properties (isOriented,
                                          part->n_cell,
                                          part->n_face,
                                          part->face_vtx_idx,
                                          part->face_vtx,
                                          part->cell_face_idx,
                                          part->cell_face,
                                          part->n_vtx,
                                          part->vtx,
                                          volume,
                                          cellCenter,
                                          NULL,
                                          &isDegenerated);
    }
    else /*Trash patch */
    {
      /* Allocate */
      double *cell_weight = (double *) malloc (part->n_cell * sizeof(double));

      /* Nulliffy cellCenterArray */
      for(int iCell = 0; iCell < part->n_cell; iCell++) {
        cellCenter[3*iCell  ] = 0.;
        cellCenter[3*iCell+1] = 0.;
        cellCenter[3*iCell+2] = 0.;
        cell_weight[iCell]     = 0.;
      }

      /* Compute */
      for(int iCell = 0; iCell < part->n_cell; iCell++) {

        /* Cellule composé de n_face */
        int aFac = part->cell_face_idx[iCell];
        int nFac = part->cell_face_idx[iCell+1] - aFac;

        for(int iFac = 0; iFac < nFac; iFac++) {

          /* Face composé de n_vtx */
          int lFac = PDM_ABS(part->cell_face[aFac + iFac]) - 1;

          int aVtx = part->face_vtx_idx[lFac];
          int n_vtx = part->face_vtx_idx[lFac+1] - aVtx;

          for(int iVtx = 0; iVtx < n_vtx; iVtx++) {

            /* Face composé de n_vtx */
            int lVtx = part->face_vtx[aVtx + iVtx] - 1;

            /* Add to current cell and stack weight */
            cellCenter[3*iCell  ] += part->vtx[3*lVtx  ];
            cellCenter[3*iCell+1] += part->vtx[3*lVtx+1];
            cellCenter[3*iCell+2] += part->vtx[3*lVtx+2];

            cell_weight[iCell] += 1.;
          }
        }
      }

      /* Nulliffy cellCenterArray */
      for(int iCell = 0; iCell < part->n_cell; iCell++) {
        cellCenter[3*iCell  ] = cellCenter[3*iCell  ]/cell_weight[iCell];
        cellCenter[3*iCell+1] = cellCenter[3*iCell+1]/cell_weight[iCell];
        cellCenter[3*iCell+2] = cellCenter[3*iCell+2]/cell_weight[iCell];
      }

      /* Verbose */
      if(0 == 1){
        for(int iCell = 0; iCell < part->n_cell; iCell++) {
          PDM_printf("cellCenter (X,Y,Z) : %f - %f - %f \n", cellCenter[3*iCell  ], cellCenter[3*iCell+1], cellCenter[3*iCell+2]);
          PDM_printf("cell_weight         : %f  \n", cell_weight[iCell  ]);
        }
      }

      /* Free */
      free(cell_weight);

    }

    /* Free */
    free (volume);
  }
  else {   /* isPoly3D */
    double *surface_vector = (double * ) malloc( sizeof(double) * 3 * part->n_cell);
    int isDegenerated;

    int *connectivity = (int *) malloc (part->cell_face_idx[part->n_cell]
                        * sizeof(int));


    int idx = 0;
    for (int icell = 0; icell < part->n_cell; icell++) {
      int faceCurr = PDM_ABS(part->cell_face[part->cell_face_idx[icell]]) - 1;
      int deb  = part->face_vtx[part->face_vtx_idx[faceCurr]];
      int next = part->face_vtx[part->face_vtx_idx[faceCurr]+1];
      connectivity[idx++] = deb;
      connectivity[idx++] = next;
      while (next != deb) {
        for (int j = part->cell_face_idx[icell]; j < part->cell_face_idx[icell+1]; ++j) {
          int face = PDM_ABS(part->cell_face[j]) - 1;
          if (faceCurr != face) {
            int s1 = part->face_vtx[part->face_vtx_idx[face]  ];
            int s2 = part->face_vtx[part->face_vtx_idx[face]+1];

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

		assert (idx == part->cell_face_idx[part->n_cell]);

    PDM_geom_elem_polygon_properties (part->n_cell,
                                      part->cell_face_idx,
                                      connectivity,
                                      part->vtx,
                                      surface_vector,
                                      cellCenter,
                                      NULL,
                                      &isDegenerated);

    free (surface_vector);
    free (connectivity);

  }

}

/**
 *
 * \brief Quick sort
 *
 * \param [inout]   a     Array to sort
 * \param [in]      l     First element
 * \param [in]      r     Last  element
 *
 */

static void
_quickSort_int
(
 int a[],
 int l,
 int r
)
{
  if (l < r) {
    int j = r+1;
    int t;
    int pivot = a[l];
    int i = l;

    while(1) {
      do ++i; while (a[i] <= pivot && i < r);
      do --j; while (a[j] > pivot);
      if (i >= j) break;

      t    = a[i];
      a[i] = a[j];
      a[j] = t;

    }
    t    = a[l];
    a[l] = a[j];
    a[j] = t;

    _quickSort_int(a, l  , j-1);
    _quickSort_int(a, j+1,   r);
  }
}

/**
 *
 * \brief Builds dual graph face cell connectivity
 *
 * \param [inout] part_ini                 Part object - fine mesh partition
 *
 * \param [inout] cell_cell_comp_idx    Array of indexes of the dual graph
 * \param [inout] cell_cell_comp       Dual graph
 *
 */

static void
_dual_graph_firstrank
(
  _part_t        *part_ini,
  int           **cell_cell_comp_idx,
  int           **cell_cell_comp
)
{
  //cell_cell_n: array of counters of the numbers of connectivities
  //cell_cell: dual graph to be built
  //cell_cell_idx: array of indexes of the dual graph (same as cell_face_idx)

  int *cell_cell_n = (int *) malloc(part_ini->n_cell * sizeof(int));
  for (int i = 0; i < part_ini->n_cell; i++) {
    cell_cell_n[i] = 0;
  }

  int *cell_cell = (int *) malloc(part_ini->cell_face_idx[part_ini->n_cell] * sizeof(int));
  for (int i = 0; i < part_ini->cell_face_idx[part_ini->n_cell]; i++) {
    cell_cell[i] = -1;
  }

  int *cell_cell_idx = (int *) malloc((part_ini->n_cell + 1) * sizeof(int));
  for(int i = 0; i < part_ini->n_cell + 1; i++) {
    cell_cell_idx[i] = part_ini->cell_face_idx[i];
  }

  for (int i = 0; i < part_ini->n_face; i++) {
    int i_cell1 = PDM_ABS (part_ini->face_cell[2*i    ]) - 1;
    int i_cell2 = PDM_ABS (part_ini->face_cell[2*i + 1]) - 1;
    //Only the non-boundary faces are stored
    if (i_cell2 > 0) {
      int idx1 = cell_cell_idx[i_cell1] + cell_cell_n[i_cell1];
      cell_cell[idx1] = i_cell2 + 1;
      cell_cell_n[i_cell1] += 1;

      int idx2 = cell_cell_idx[i_cell2] + cell_cell_n[i_cell2];
      cell_cell[idx2] = i_cell1 + 1;
      cell_cell_n[i_cell2] += 1;
    }
  }

  if (0 == 1) {
    PDM_printf("Content of cell_cell_n after looping over cell_face: ");
    for(int i = 0; i < part_ini->n_cell; i++) {
      PDM_printf(" %d ", cell_cell_n[i]);
    }
    PDM_printf("\n");

    PDM_printf("Content of cell_cell after looping over cell_face: ");
    for(int i = 0; i < part_ini->cell_face_idx[part_ini->n_cell]; i++) {
      PDM_printf(" %d ", cell_cell[i]);
    }
    PDM_printf("\n");
  }

  // cell_cell_idx is rebuilt
  *cell_cell_comp_idx = malloc((part_ini->n_cell + 1) * sizeof(int));

  (*cell_cell_comp_idx)[0] = 0;
  for(int i = 0; i < part_ini->n_cell; i++) {
    (*cell_cell_comp_idx)[i + 1] = (*cell_cell_comp_idx)[i] + cell_cell_n[i];
  }

  //We compress the dual graph since cell_cell_idx was built from cell_face_idx
  //We have then n_face elements in cell_cell whereas it needs to be composed of n_cell elements

  //    PDM_printf("(*cell_cell_comp_idx)[part_ini->n_cell] : %d \n", (*cell_cell_comp_idx)[part_ini->n_cell]);
  *cell_cell_comp = malloc((*cell_cell_comp_idx)[part_ini->n_cell] * sizeof(int));

  int cpt_cell_cell_comp = 0;
  for(int i = 0; i < part_ini->cell_face_idx[part_ini->n_cell]; i++) {
    //        PDM_printf("I am testing a value for the %d time! \n", i);

    //We have an information to store when a neighboring cell exists
    if(cell_cell[i] > -1){
      //            PDM_printf("I am storing a value for the %d time! \n", i);
      //We add a -1 to have the graph vertices numbered from 0 to n (C numbering)
      (*cell_cell_comp)[cpt_cell_cell_comp++] = cell_cell[i] - 1;
      //            PDM_printf("Valeur stockee : %d \n ", (*cell_cell_comp)[cpt_cell_cell_comp - 1]);
    }
  }

  if( 0 == 1) {
    PDM_printf("Content of cell_cell_comp after compression and renumbering: ");
    for(int i = 0; i < (*cell_cell_comp_idx)[part_ini->n_cell]; i++) {
      PDM_printf(" %d ", (*cell_cell_comp)[i]);
    }
    PDM_printf("\n");
  }

  /* Free temporary arrays*/

  free(cell_cell_n);
  free(cell_cell);
  free(cell_cell_idx);

  //Remove duplicate cells of the dual graph
  //We use the following scheme:
  //We loop over the indexes for the whole array to subdivide it into subarrays
  //We sort locally each subarray (determined thanks to cell_cell_comp_idx)
  //We loop over each subarray
  //We store the first value of each subarray anyway
  //We store each non-duplicated value and increment the writing index
  //We update the index array at each iteration

  int idx_write = 0;
  int tab_idx_temp = 0;

  for (int i = 0; i < part_ini->n_cell; i++) {
    _quickSort_int((*cell_cell_comp), tab_idx_temp, (*cell_cell_comp_idx)[i + 1] - 1);

    int last_value = -1;

    for (int j = tab_idx_temp; j < (*cell_cell_comp_idx)[i + 1]; j++) {
      //We need to have a local index (between 0 and n_face)
      //If the value is different from the previous one (higher than is the same as different since the array is sorted)

      if(last_value != (*cell_cell_comp)[j]) {
        (*cell_cell_comp)[idx_write++] = (*cell_cell_comp)[j];
        last_value = (*cell_cell_comp)[j];
      }
    }

    if (0 == 1) {
      PDM_printf("\n Contenu de cell_cell_comp apres reecriture: \n");
      for(int i1 = 0; i1 < (*cell_cell_comp_idx)[part_ini->n_cell]; i1++) {
        PDM_printf(" %d ", (*cell_cell_comp)[i1]);
      }
      PDM_printf("\n");
    }

    tab_idx_temp = (*cell_cell_comp_idx)[i + 1];
    (*cell_cell_comp_idx)[i + 1] = idx_write;

    if (0 == 1) {
      PDM_printf("\n Contenu de cell_cell_comp_idx apres reecriture: \n");
      for(int i1 = 0; i1 < part_ini->n_cell + 1; i1++) {
        PDM_printf(" %d ", (*cell_cell_comp_idx)[i1]);
      }
      PDM_printf("\n");
    }
  }

  if (0 == 1) {
    PDM_printf("Content of cell_cell_comp_idx after compression: ");
    for(int i1 = 0; i1 < part_ini->n_cell + 1; i1++) {
      PDM_printf(" %d ", (*cell_cell_comp_idx)[i1]);
    }
    PDM_printf("\n");

    PDM_printf("Content of cell_cell_comp after compression: ");
    for(int i1 = 0; i1 < (*cell_cell_comp_idx)[part_ini->n_cell]; i1++) {
      PDM_printf(" %d ", (*cell_cell_comp)[i1]);
    }
    PDM_printf("\n");
  }

  //We reallocate the memory in case of duplicated values removed
  //The new array size is idx_write (stored in (*cell_cell_comp_idx)[part_ini->n_cell])
  *cell_cell_comp = realloc(*cell_cell_comp,
                                (*cell_cell_comp_idx)[part_ini->n_cell] * sizeof(int));

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
 _part_t **mesh_parts,
 int       n_part,
 void     *specific_data
)
{
  PDM_UNUSED(specific_data);

  for(int i_part = 0; i_part < n_part; ++i_part) {
    _part_t *part = mesh_parts[i_part];
    double *cellCenter =
        (double *) malloc (part->n_cell * 3 * sizeof(double ));
    PDM_hilbert_code_t *hilbert_codes =
        (PDM_hilbert_code_t *) malloc (part->n_cell * sizeof(PDM_hilbert_code_t));

    /** Barycentre computation **/

    _compute_cellCenter (part, cellCenter);

    double extents[3 * 2];

    /** Get EXTENTS LOCAL **/

    PDM_hilbert_get_coord_extents_seq(3, part->n_cell, cellCenter, extents);

    /** Hilbert Coordinates Computation **/

    PDM_hilbert_encode_coords(3, PDM_HILBERT_CS, extents, part->n_cell, cellCenter, hilbert_codes);

    /** CHECK H_CODES **/

    free(cellCenter);

    int *new_to_old_order = (int *) malloc (part->n_cell * sizeof(int));
    for(int i = 0; i < part->n_cell; ++i) {
      new_to_old_order [i] = i;
    }

    PDM_sort_double (hilbert_codes, new_to_old_order, part->n_cell);

    PDM_part_reorder_cell (part, new_to_old_order);

    free (hilbert_codes);
    free (new_to_old_order);

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
 _part_t **mesh_parts,
 int       n_part,
 void     *specific_data
)
{
  PDM_UNUSED(specific_data);

  /** Loop over all part of the current process **/
  for(int i_part = 0; i_part < n_part; ++i_part) {
    /** Get current part id **/
    _part_t *part = mesh_parts[i_part];
    const int n_cell = part->n_cell;

    /** Allocate reoerdering/permutation array **/
    int *order = (int *) malloc (sizeof(int) * n_cell);

    /** Graph computation (in the new partition ) **/
    int *dual_graph_idx = NULL;
    int *dual_graph     = NULL;
    _dual_graph_firstrank(part,
                         (int **) &dual_graph_idx,
                         (int **) &dual_graph);

    /** Verbose bandwidth **/
    // int dualBandWidth;
    // dualBandWidth = PDM_cuthill_checkbandwidth(part->n_cell, dual_graph_idx, dual_graph);
    // PDM_printf("Bandwidth of graph before reordering : %d \n", dualBandWidth);
    // PDM_printf("Bandwidth of graph before reordering \n");

    /** Compute reordering **/
    PDM_cuthill_generate(part->n_cell, dual_graph_idx, dual_graph, order);

    /** Apply renumbering **/
    PDM_part_reorder_cell(part, order);

    /** Verbose bandwidth **/
    // dualBandWidth = PDM_cuthill_checkbandwidth(part->n_cell, dual_graph_idx, dual_graph);
    // PDM_printf("Bandwidth of graph after reordering : %d \n", dualBandWidth);

    /* Copy in partition */
    if(part->new_to_old_order_cell == NULL){
      part->new_to_old_order_cell = (int *) malloc (sizeof(int) * n_cell);
      for (int i = 0; i < n_cell; i++){
        part->new_to_old_order_cell[i] = order[i];
      }
    }

    /** Free memory **/
    free(order);
    free(dual_graph_idx);
    free(dual_graph);
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
 _part_t **mesh_parts,
 int       n_part,
 void     *specific_data
)
{
  PDM_UNUSED(specific_data);

  for(int i_part = 0; i_part < n_part; ++i_part) {
    _part_t *part = mesh_parts[i_part];
    const int n_cell = part->n_cell;

    int *order = (int *) malloc (sizeof(int) * n_cell);

    _random_order (n_cell, order);

    PDM_part_reorder_cell (part, order);

    /* Copy in partition */
    if(part->new_to_old_order_cell == NULL){
      part->new_to_old_order_cell = (int *) malloc (sizeof(int) * n_cell);
      for (int i = 0; i < n_cell; i++){
        part->new_to_old_order_cell[i] = order[i];
      }
    }

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
 _part_t **mesh_parts,
 int       n_part,
 void     *specific_data
)
{
  PDM_UNUSED(specific_data);

  for(int i_part = 0; i_part < n_part; ++i_part) {
    _part_t *part = mesh_parts[i_part];
    const int n_face = part->n_face;

    int *order = (int *) malloc (sizeof(int) * n_face);

    _random_order (n_face, order);

    PDM_part_reorder_face (part, order);

    /* Copy in partition */
    if(part->new_to_old_order_face == NULL){
      part->new_to_old_order_face = (int *) malloc (sizeof(int) * n_face);
      for (int i = 0; i < n_face; i++){
        part->new_to_old_order_face[i] = order[i];
      }
    }

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
 _part_t** mesh_parts,
 int       n_part,
 void     *specific_data
)
{
  PDM_UNUSED(specific_data);

  for(int i_part = 0; i_part < n_part; ++i_part) {
    _part_t *part = mesh_parts[i_part];
    const int n_face = part->n_face;

    int *order = (int *) malloc (sizeof(int) * n_face);

    /** Build a pre-array face cell ordered */
    int *face_cell_tmp = (int *) malloc(2*n_face * sizeof(int));

    for(int i = 0; i < n_face; i++) {
       int icell1 = PDM_ABS (part->face_cell[2*i  ]);
       int icell2 = PDM_ABS (part->face_cell[2*i+1]);
       if(icell1 < icell2 )
       {
          face_cell_tmp[2*i  ] = icell2;
          face_cell_tmp[2*i+1] = icell1;
       }
       else
       {
          face_cell_tmp[2*i  ] = icell1;
          face_cell_tmp[2*i+1] = icell2;
       }
    }

    /** Reorder lexicographicly the array */
    PDM_order_lnum_s (face_cell_tmp, 2, order, n_face);

    /** Update face array with the new array **/
    PDM_part_reorder_face (part, order);

    if(part->new_to_old_order_face == NULL){
      part->new_to_old_order_face = (int *) malloc (sizeof(int) * n_face);
      for (int i = 0; i < n_face; i++){
        part->new_to_old_order_face[i] = order[i];
      }
    }

    /** Free memory **/
    free (order);
    free (face_cell_tmp);
  }
}


static void
_renum_vtx_sort_int_ext
(
 _part_t **mesh_parts,
 int       n_part,
 void     *specific_data
)
{
  assert(specific_data != NULL);

  /* Specific data in this context is a array for each part of priority give by graph_comm */
  int** pvtx_priority = (int ** ) specific_data;

  for(int i_part = 0; i_part < n_part; ++i_part) {
    _part_t *part = mesh_parts[i_part];
    const int n_vtx = part->n_vtx;

    int *order    = (int *) malloc (sizeof(int) * n_vtx);

    // _random_order (n_vtx, order);

    /* By default the partition own all vertex */
    int n_kind = 3;
    int* type_idx = (int * ) malloc( (n_kind+1) * sizeof(int));
    int* type_n   = (int * ) malloc( (n_kind+1) * sizeof(int));
    for (int i = 0; i < n_kind + 1; i++){
      type_idx[i] = 0;
      type_n  [i] = 0;
    }

    for (int i = 0; i < n_vtx; i++){
      type_idx[pvtx_priority[i_part][i]+1]++;
    }

    for (int i = 0; i < n_kind; i++){
      type_idx[i+1] += type_idx[i];
    }

    /* Setup reordering */
    for (int i = 0; i < n_vtx; i++){
      int priority = pvtx_priority[i_part][i];
      int idx = type_idx[priority] + type_n[priority];
      order[idx] = i;
      type_n[priority]++;
    }
    free(type_n);
    free(type_idx);

    /*
     * Verbose
     */
    if(0 == 1) {
      printf(" vtx_order : ");
      for (int i = 0; i < n_vtx; i++){
         printf(" %i", order[i]);
      }
      printf("\n");
    }

    PDM_part_reorder_vtx (part, order);

    /* Copy in partition */
    if(part->new_to_old_order_vtx == NULL){
      part->new_to_old_order_vtx = (int *) malloc (sizeof(int) * n_vtx);
      for (int i = 0; i < n_vtx; i++){
        part->new_to_old_order_vtx[i] = order[i];
      }
    }

    free (order   );
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
PDM_part_renum_method_purge
(
 void
)
{
  if (face_methods != NULL) {
    // printf("PDM_part_renum_method_purge: face_methods\n");

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

    // printf("PDM_part_renum_method_purge: cell_methods\n");
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

void
PROCF (pdm_part_renum_method_cell_idx_get_cf, PDM_PART_RENUM_METHOD_CELL_IDX_GET_CF)
(
 char *name,
 int  *l_name,
 int  *idx
 )
{
  char *_name = PDM_fortran_to_c_string (name, *l_name);

  *idx = PDM_part_renum_method_cell_idx_get (_name);

  free (_name);

}


int
PDM_part_renum_method_cell_idx_get
(
const char *name
)
{
  if (cell_methods == NULL) {
    PDM_part_renum_method_load_local();
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

void
PROCF (pdm_part_renum_method_face_idx_get_cf, PDM_PART_RENUM_METHOD_FACE_IDX_GET_CF)
(
 char *name,
 int  *l_name,
 int  *idx
 )
{
  char *_name = PDM_fortran_to_c_string (name, *l_name);

  *idx = PDM_part_renum_method_face_idx_get (_name);

  free (_name);

}

int
PDM_part_renum_method_face_idx_get
(
const char *name
)
{
  if (face_methods == NULL) {
    PDM_part_renum_method_load_local();
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

void
PROCF (pdm_part_renum_method_cell_name_get_cf, PDM_PART_RENUM_METHOD_CELL_NAME_GET_CF)
(
 char *name,
 int  *l_name,
 int  *idx
 )
{
  const char *_name = PDM_part_renum_method_cell_name_get (*idx);

  const int _l_name = strlen(_name);

  *l_name = PDM_MAX (_l_name, PDM_MAX_CHAR_LENGTH);

  strncpy (name, _name, *l_name);
}

const char *
PDM_part_renum_method_cell_name_get
(
const int idx
)
{
  if (cell_methods == NULL) {
    PDM_part_renum_method_load_local();
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
 * \brief Get the number of renumbering cell methods
 *
 * \return Number of methods
 *
 */

void
PROCF (pdm_part_n_renum_method_cell_get, PDM_PART_N_RENUM_METHOD_CELL_GET)
(
 int  *n_method
 )
{
  *n_method = PDM_part_n_renum_method_cell_get ();
}

int
PDM_part_n_renum_method_cell_get
(
void
)
{
  if (cell_methods == NULL) {
    PDM_part_renum_method_load_local();
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

void
PROCF (pdm_part_n_renum_method_face_get, PDM_PART_N_RENUM_METHOD_FACE_GET)
(
 int  *n_method
 )
{
  *n_method = PDM_part_n_renum_method_face_get ();
}

int
PDM_part_n_renum_method_face_get
(
void
)
{
  if (face_methods == NULL) {
    PDM_part_renum_method_load_local();
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
PDM_part_renum_method_cell_add
(
 const char                 *name,     /*!< Name          */
 const PDM_part_renum_fct_t  renum_fct /*!< Customize \ref PDM_part_renum_cell function for the format */
)
{
  if (cell_methods == NULL) {
    PDM_part_renum_method_load_local ();
  }

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
PDM_part_renum_method_face_add
(
 const char                 *name,     /*!< Name          */
 const PDM_part_renum_fct_t  renum_fct /*!< Customize \ref PDM_part_renum_face function for the format */
)
{
  if (face_methods == NULL) {
    PDM_part_renum_method_load_local ();
  }

  _renum_method_t *method_ptr = malloc (sizeof(_renum_method_t));

  int idx = PDM_Handles_store  (face_methods, method_ptr);

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
PDM_part_renum_method_vtx_add
(
 const char                 *name,     /*!< Name          */
 const PDM_part_renum_fct_t  renum_fct /*!< Customize \ref PDM_part_renum_face function for the format */
)
{
  if (vtx_methods == NULL) {
    PDM_part_renum_method_load_local ();
  }

  _renum_method_t *method_ptr = malloc (sizeof(_renum_method_t));

  int idx = PDM_Handles_store  (vtx_methods, method_ptr);

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
PDM_part_renum_method_load_local
(
void
)
{
  if (cell_methods == NULL)  {

    const int n_default_methods = 4;
    cell_methods = PDM_Handles_create (n_default_methods);

    PDM_part_renum_method_cell_add ("PDM_PART_RENUM_CELL_NONE",
                             NULL);
    PDM_part_renum_method_cell_add ("PDM_PART_RENUM_CELL_RANDOM",
                             _renum_cells_random);
    PDM_part_renum_method_cell_add ("PDM_PART_RENUM_CELL_HILBERT",
                             _renum_cells_hilbert);
    PDM_part_renum_method_cell_add ("PDM_PART_RENUM_CELL_CUTHILL",
                             _renum_cells_cuthill);
  }

  if (face_methods == NULL)  {
    const int n_default_methods = 3;
    face_methods = PDM_Handles_create (n_default_methods);

    PDM_part_renum_method_face_add ("PDM_PART_RENUM_FACE_NONE",
                                    NULL);
    PDM_part_renum_method_face_add ("PDM_PART_RENUM_FACE_RANDOM",
                                    _renum_faces_random);
    PDM_part_renum_method_face_add ("PDM_PART_RENUM_FACE_LEXICOGRAPHIC",
                                    _renum_faces_lexicographic);
  }

  if (vtx_methods == NULL)  {
    const int n_default_methods = 2;
    vtx_methods = PDM_Handles_create (n_default_methods);

    PDM_part_renum_method_vtx_add ("PDM_PART_RENUM_VTX_NONE",
                                    NULL);
    PDM_part_renum_method_vtx_add ("PDM_PART_RENUM_VTX_SORT_INT_EXT",
                                    _renum_vtx_sort_int_ext);
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
 _part_t **mesh_parts,
 int       n_part,
 int       renum_cell_method,
 void     *specific_data
)
{

  if (cell_methods == NULL)  {
    PDM_part_renum_method_load_local ();
  }

  const _renum_method_t *method_ptr = (const _renum_method_t *)
                                    PDM_Handles_get (cell_methods, renum_cell_method);
                                    // PDM_Handles_get (cell_methods, ppart->renum_cell_method);

  PDM_part_renum_fct_t fct = method_ptr->fct;

  if (fct != NULL) {
    (fct) (mesh_parts, n_part, specific_data);
  }

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

void
PROCF (pdm_part_renum_method_face_name_get_cf, PDM_PART_RENUM_METHOD_FACE_NAME_GET_CF)
(
 char *name,
 int  *l_name,
 int  *idx
 )
{
  const char *_name = PDM_part_renum_method_face_name_get (*idx);

  const int _l_name = strlen(_name);

  *l_name = PDM_MAX (_l_name, PDM_MAX_CHAR_LENGTH);

  strncpy (name, _name, *l_name);
}


const char *
PDM_part_renum_method_face_name_get
(
const int idx
)
{
  if (face_methods == NULL) {
    PDM_part_renum_method_load_local();
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
 * \brief Perform mesh entities renumbering
 *
 * \param [in,out]  part       part structure
 *
 */

void
PDM_part_renum_face
(
 _part_t **mesh_parts,
 int       n_part,
 int       renum_face_method,
 void     *specific_data
)
{
  if (face_methods == NULL)  {
    PDM_part_renum_method_load_local ();
  }

  const _renum_method_t *method_ptr = (const _renum_method_t *)
                                      PDM_Handles_get (face_methods, renum_face_method);
                                      // PDM_Handles_get (face_methods, ppart->renum_face_method);

  PDM_part_renum_fct_t fct = method_ptr->fct;

  if (fct != NULL) {
    (fct) (mesh_parts, n_part, specific_data);
  }
}

/**
 *
 * \brief Perform mesh entities renumbering
 *
 * \param [in,out]  part       part structure
 *
 */

void
PDM_part_renum_edge
(
 _part_t **mesh_parts,
 int       n_part,
 int       renum_edge_method,
 void     *specific_data
)
{
  if (edge_methods == NULL)  {
    PDM_part_renum_method_load_local ();
  }

  const _renum_method_t *method_ptr = (const _renum_method_t *)
                                      PDM_Handles_get (edge_methods, renum_edge_method);

  PDM_part_renum_fct_t fct = method_ptr->fct;

  if (fct != NULL) {
    (fct) (mesh_parts, n_part, specific_data);
  }
}

/**
 *
 * \brief Perform mesh entities renumbering
 *
 * \param [in,out]  part       part structure
 *
 */

void
PDM_part_renum_vtx
(
 _part_t **mesh_parts,
 int       n_part,
 int       renum_vtx_method,
 void     *specific_data
)
{
  if (vtx_methods == NULL)  {
    PDM_part_renum_method_load_local ();
  }

  const _renum_method_t *method_ptr = (const _renum_method_t *)
                                      PDM_Handles_get (vtx_methods, renum_vtx_method);

  PDM_part_renum_fct_t fct = method_ptr->fct;

  if (fct != NULL) {
    (fct) (mesh_parts, n_part, specific_data);
  }
}


/**
 *
 * \brief Perform cells renumbering from a new order
 *
 * \param [in,out]  part        Current partition
 * \param [in]      new_to_old_order    NewOrder
 *
 */
void
PDM_part_reorder_cell
(
 _part_t *part,
 int     *new_to_old_order
)
{
  /*
   * Cell Renumbering
   */

  PDM_part_renum_connectivities (part->n_cell,
                                 new_to_old_order,
                                 part->cell_face_idx,
                                 part->cell_face);

  if (part->cell_tag != NULL) {
    PDM_order_array (part->n_cell,
                     sizeof(int),
                     new_to_old_order,
                     part->cell_tag);
  }

  if (part->cell_color != NULL) {
    PDM_order_array (part->n_cell,
                     sizeof(int),
                     new_to_old_order,
                     part->cell_color);
  }

  if (part->thread_color != NULL) {
    PDM_order_array (part->n_cell,
                     sizeof(int),
                     new_to_old_order,
                     part->thread_color);
  }

  if (part->hyperplane_color != NULL) {
    PDM_order_array (part->n_cell,
                     sizeof(int),
                     new_to_old_order,
                     part->hyperplane_color);
  }

  if (part->new_to_old_order_cell != NULL) {
    // printf("PDM_order_array :new_to_old_order_cell \n");
    PDM_order_array (part->n_cell,
                     sizeof(int),
                     new_to_old_order,
                     part->new_to_old_order_cell);
  }

  PDM_order_array (part->n_cell,
                   sizeof(PDM_g_num_t),
                   new_to_old_order,
                   part->cell_ln_to_gn);

  _renum_faceCell (part->n_cell,
                   part->n_face,
                   part->face_cell,
                   new_to_old_order);

}


/**
 *
 * \brief Perform faces renumbering from a new order
 *
 * \param [in,out]  part        Current partition
 * \param [in]      new_to_old_order    NewOrder
 *
 */
void
PDM_part_reorder_face
(
_part_t *part,
int     *new_to_old_order
)
{

  /** Renum face_vtx / face_vtx_idx **/
  PDM_part_renum_connectivities (part->n_face,
                                 new_to_old_order,
                                 part->face_vtx_idx,
                                 part->face_vtx);

  /** Renum face_edge / face_edge_idx **/
  if(part->face_edge != NULL) {
    PDM_part_renum_connectivities (part->n_face,
                                   new_to_old_order,
                                   part->face_edge_idx,
                                   part->face_edge);
  }

  /** cell_face **/
  int *old_to_new_order = (int *) malloc (part->n_face * sizeof(int));
  for(int i = 0; i < part->n_face; i++) {
   old_to_new_order[new_to_old_order[i]] = i;
  }

  PDM_part_renum_array (part->cell_face_idx[part->n_cell],
                        old_to_new_order,
                        part->cell_face);

  /** face_tag **/
  if (part->face_tag != NULL) {
    PDM_order_array (part->n_face,
                     sizeof(int),
                     new_to_old_order,
                     part->face_tag);
  }

  /** face_color **/
  if (part->face_color != NULL) {
    PDM_order_array (part->n_face,
                     sizeof(int),
                     new_to_old_order,
                     part->face_color);
  }

  /** face_color **/
  if (part->new_to_old_order_face != NULL) {
    // printf("PDM_order_array :new_to_old_order_face \n");
    PDM_order_array (part->n_face,
                     sizeof(int),
                     new_to_old_order,
                     part->new_to_old_order_face);
  }

   /** face_ln_to_gn **/
  PDM_order_array (part->n_face,
                   sizeof(PDM_g_num_t),
                   new_to_old_order,
                   part->face_ln_to_gn); // OK

  /** face_group **/
  if (part->face_group != NULL) {
    PDM_part_renum_array (part->face_group_idx[part->n_face_group],
                          old_to_new_order,
                          part->face_group); // OK
  }

  /** face_cell Face **/
  _order_face_cell (part->n_face,
                    new_to_old_order,
                    part->face_cell);

  /* Free */
  free (old_to_new_order);

}

/**
 *
 * \brief Perform vtx renumbering from a new order
 *
 * \param [in,out]  part        Current partition
 * \param [in]      new_to_old_order    NewOrder
 *
 */
void
PDM_part_reorder_edge
(
_part_t *part,
int     *new_to_old_order
)
{
  /** face_vtx **/
  int *old_to_new_order = (int *) malloc (part->n_edge * sizeof(int));
  for(int i = 0; i < part->n_edge; i++) {
   old_to_new_order[new_to_old_order[i]] = i;
  }

  if (part->face_edge != NULL) {
    PDM_part_renum_array (part->face_edge_idx[part->n_face],
                          old_to_new_order,
                          part->face_edge);
  }

  if(part->edge_face != NULL) {
    PDM_part_renum_connectivities (part->n_face,
                                   new_to_old_order,
                                   part->edge_face_idx,
                                   part->edge_face);
  }

  if( part->edge_vtx != NULL) {
    _order_face_cell(part->n_edge,
                     new_to_old_order,
                     part->edge_vtx);
  }

  /** edge_tag **/
  if (part->edge_tag != NULL) {
    PDM_order_array (part->n_edge,
                     sizeof(int),
                     new_to_old_order,
                     part->edge_tag);
  }

  /** edge_color **/
  if (part->edge_color != NULL) {
    PDM_order_array (part->n_edge,
                     sizeof(int),
                     new_to_old_order,
                     part->edge_color);
  }

  if (part->new_to_old_order_edge != NULL) {
    // printf("PDM_order_array :new_to_old_order_edge \n");
    PDM_order_array (part->n_edge,
                     sizeof(int),
                     new_to_old_order,
                     part->new_to_old_order_edge);
  }

  if( part->edge_ln_to_gn != NULL) {
    PDM_order_array (part->n_edge,
                     sizeof(PDM_g_num_t),
                     new_to_old_order,
                     part->edge_ln_to_gn); // OK
  }

  free (old_to_new_order);
}

/**
 *
 * \brief Perform vtx renumbering from a new order
 *
 * \param [in,out]  part        Current partition
 * \param [in]      new_to_old_order    NewOrder
 *
 */
void
PDM_part_reorder_vtx
(
_part_t *part,
int     *new_to_old_order
)
{
  /** face_vtx **/
  int *old_to_new_order = (int *) malloc (part->n_vtx * sizeof(int));
  for(int i = 0; i < part->n_vtx; i++) {
   old_to_new_order[new_to_old_order[i]] = i;
  }

  if (part->n_face!=0) {
    PDM_part_renum_array (part->face_vtx_idx[part->n_face],
                          old_to_new_order,
                          part->face_vtx);
  } else { // the mesh is supposed to be described by elements
    int n_section = part->n_section;
    assert(n_section!=0);
    for (int i_section=0; i_section<n_section; ++i_section) {
      int n_elt = part->n_elt[i_section];
      int n_elt_vtx = part->elt_vtx_idx[i_section][n_elt];
      PDM_part_renum_array (n_elt_vtx,
                            old_to_new_order,
                            part->elt_vtx[i_section]);
    }
    PDM_order_array (part->n_vtx,
                     sizeof(int),
                     new_to_old_order,
                     part->vtx_ghost_information);
  }

  if(part->edge_vtx != NULL) {
    printf("PDM_part_reorder_vtx TO DEBUG \n");
    abort();
    PDM_part_renum_array (2 * part->n_edge,
                          old_to_new_order,
                          part->edge_vtx);
  }

  /** vtx **/
  PDM_order_array (part->n_vtx,
                   3*sizeof(double),
                   new_to_old_order,
                   part->vtx);

  /** vtx_tag **/
  if (part->vtx_tag != NULL) {
    PDM_order_array (part->n_vtx,
                     sizeof(int),
                     new_to_old_order,
                     part->vtx_tag);
  }

  /** vtx_color **/
  // if (part->vtx_color != NULL) {
  //   PDM_order_array (part->n_vtx,
  //                    sizeof(int),
  //                    new_to_old_order,
  //                    part->vtx_color);
  // }

  /** vtx_color **/
  if (part->new_to_old_order_vtx != NULL) {
    // printf("PDM_order_array :new_to_old_order_vtx \n");
    PDM_order_array (part->n_vtx,
                     sizeof(int),
                     new_to_old_order,
                     part->new_to_old_order_vtx);
  }

   /** vtx_ln_to_gn **/
  PDM_order_array (part->n_vtx,
                   sizeof(PDM_g_num_t),
                   new_to_old_order,
                   part->vtx_ln_to_gn); // OK

  /* Free */
  free (old_to_new_order);

}

#ifdef __cplusplus
}
#endif /* __cplusplus */
