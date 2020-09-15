
/*----------------------------------------------------------------------------
 *  Standar headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_part_geom.h"
#include "pdm_hilbert.h"
#include "pdm_sort.h"
#include "pdm_geom_elem.h"
#include "pdm_binary_search.h"
#include "pdm_block_to_part.h"
#include "pdm_printf.h"


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
 * Private function definitions
 *============================================================================*/


static void
_prepare_connectivity
(
const PDM_MPI_Comm    comm,
const PDM_g_num_t     dn_cell,
const int            *dcell_face_idx,
const PDM_g_num_t    *dcell_face,
const int            *dface_vtx_idx,
const PDM_g_num_t    *dface_vtx,
const PDM_g_num_t    *dface_proc,
int                **face_vtx,
int                **face_vtx_idx,
int                  *sizeface_vtx_idx,
const PDM_g_num_t    *dvtx_proc,
const double         *dvtx_coord,
      double       **lvtx_coord
)
{
  int i_rank;
  int n_rank;

  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  /* Offset Distribution to work with block_to_part */
  PDM_g_num_t* dface_proc_loc = (PDM_g_num_t *) malloc( sizeof(PDM_g_num_t) * (n_rank + 1));
  PDM_g_num_t* dvtx_proc_loc  = (PDM_g_num_t *) malloc( sizeof(PDM_g_num_t) * (n_rank + 1));
  for (int i = 0; i < n_rank+1; i++) {
    dface_proc_loc[i] = dface_proc[i]-1;
    dvtx_proc_loc[i]  = dvtx_proc[i]-1;
  }


  /* Verbose */
  if(0 == 1){
    PDM_printf("dface_proc : \n");
    for (int i = 0; i < n_rank+1; i++) {
      PDM_printf("[%i] -> %i \n ",i,  dface_proc[i]);
      PDM_printf("[%i] -> %i \n ",i,  dface_proc_loc[i]);
    }
    for (int i = 0; i < dn_cell; i++) {
      int nFac = dcell_face_idx[i+1] - dcell_face_idx[i];
      printf("nFac[%i] :  %i \n ", i, nFac);
      for (int iVtx = dcell_face_idx[i]; iVtx < dcell_face_idx[i+1]; iVtx++) {
        printf("dcell_face[%i] :  "PDM_FMT_G_NUM" \n ", iVtx, dcell_face[iVtx]);
      }
    }
  }


  /* -------------------------------------------------------------- */
  PDM_block_to_part_t *ptb = PDM_block_to_part_create (dface_proc_loc,
                                                       &dcell_face,
                                                       &dcell_face_idx[dn_cell],
                                                        1,
                                                        comm);

  /* Echange sur les nombre de sommet de chaque face */
  int dn_face = dface_proc[i_rank+1] - dface_proc[i_rank];

  int* n_vtx_face_loc = (int *) malloc (sizeof(int) * dn_face              );
  int* n_vtx_face    = (int *) malloc (sizeof(int) * dcell_face_idx[dn_cell]);

  for (int i = 0; i < dn_face; i++) {
    n_vtx_face_loc[i] = dface_vtx_idx[i+1] - dface_vtx_idx[i];
    for (int iVtx = dface_vtx_idx[i]; iVtx < dface_vtx_idx[i+1]; iVtx++) {
    }
  }

  int stride_one = 1;

  PDM_block_to_part_exch (ptb,
                          sizeof(PDM_l_num_t),
                          PDM_STRIDE_CST,
                          &stride_one,
                          (void *) n_vtx_face_loc,
                          NULL,
                          (void **) &n_vtx_face);

  /* Verbose */
  if(0 == 1){
    for (int i = 0; i < dcell_face_idx[dn_cell]; i++) {
      printf("n_vtx_face :  %i \n ", n_vtx_face[i]);
    }
  }

  int ndcell_face_tot = 0;
  for (int i = 0; i < dcell_face_idx[dn_cell]; i++) {
    ndcell_face_tot += n_vtx_face[i];
  }
  // printf("ndcell_face_tot : %i \n", ndcell_face_tot);

  /* Alloc field */
  PDM_g_num_t* lface_vtx    = (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) * ndcell_face_tot);

  PDM_block_to_part_exch (           ptb,
                                     sizeof(PDM_g_num_t),
                                     PDM_STRIDE_VAR,
                                     n_vtx_face_loc,
                          (void *)   dface_vtx,
                                    &n_vtx_face,
                          (void **) &lface_vtx);

  /* Verbose */
  if(0 == 1){
    for (int i = 0; i < dcell_face_idx[dn_cell]; i++) {
      printf("n_vtx_face :  %i \n ", n_vtx_face[i]);
    }
    for (int i = 0; i < ndcell_face_tot; i++) {
      printf("lface_vtx[%i] :  "PDM_FMT_G_NUM" \n ",i, lface_vtx[i]);
    }
  }

  PDM_block_to_part_free(ptb);
  free (dface_proc_loc);
  /* -------------------------------------------------------------- */


  /* -------------------------------------------------------------- */
  ptb = PDM_block_to_part_create (dvtx_proc_loc,
                                  (const PDM_g_num_t **) &lface_vtx,
                                  &ndcell_face_tot,
                                  1,
                                  comm);

  int stride = 3;


  // int dn_vtx = dvtx_proc[i_rank+1] - dvtx_proc[i_rank];
  // for (int i = 0; i < dn_vtx; i++) {
  //   printf("dvtx_coord[%i] :  %12.5e \n ",i, dvtx_coord[i]);
  // }

  *lvtx_coord = (double *) malloc (sizeof(double) * 3 * ndcell_face_tot);
  PDM_block_to_part_exch (ptb,
                          sizeof(double),
                          PDM_STRIDE_CST,
                          &stride,
                          (void *) dvtx_coord,
                          NULL,
                          (void **) lvtx_coord);

  /* Verbose */
  if(0 == 1){
    for (int i = 0; i < ndcell_face_tot; i++) {
      printf("lvtx_coord[%i] :  %12.5e \n ",i, (*lvtx_coord)[i]);
    }
  }

  PDM_block_to_part_free(ptb);
  free (dvtx_proc_loc);
  /* -------------------------------------------------------------- */

  *sizeface_vtx_idx = dcell_face_idx[dn_cell] + 1 ;

  *face_vtx    = (int *) malloc( sizeof(int *) *   ndcell_face_tot              );
  *face_vtx_idx = (int *) malloc( sizeof(int *) * ( dcell_face_idx[dn_cell] + 1 ) );

  for (int i = 0; i < ndcell_face_tot; i++) {
    (*face_vtx)[i] = i+1;
  }

  (*face_vtx_idx)[0] = 0;
  for (int i = 0; i < dcell_face_idx[dn_cell]; i++) {
    (*face_vtx_idx)[i+1] = (*face_vtx_idx)[i] + n_vtx_face[i];
  }

  /* Free memory */
  free(n_vtx_face);
  free(n_vtx_face_loc);
  free(lface_vtx);

}
/**
 *
 * \brief Compute cell center of elements
 *
 * \param [in]  comm         MPI Communicator
 * \param [in]  dn_cell       Number of cells in the current process
 * \param [in]  dn_face       Number of faces in the current process
 * \param [in]  dn_vtx        Number of vertices in the current process
 * \param [in]  dcell_face_idx Index of cell_face
 * \param [in]  dcell_face    cell face connectivity in the current process
 * \param [in]  dface_vtx_idx  Index of face_vtx
 * \param [in]  dface_vtx     face vertex connectivity in the current process
 * \param [in]  dface_proc    face distribution
 * \param [in]  dvtx_coord    coordinates of vertices
 * \param [in]  dvtx_proc     Vertex distribution
 *
 * \param [out] cellCenter   Cell centers
 *
 */

static void
_compute_cell_center
(
  const PDM_MPI_Comm  comm,
  const int           dn_cell,
  const int          *dcell_face_idx,
  const PDM_g_num_t  *dcell_face,
  const int          *dface_vtx_idx,
  const PDM_g_num_t  *dface_vtx,
  const PDM_g_num_t  *dface_proc,
  const double       *dvtx_coord,
  const PDM_g_num_t  *dvtx_proc,
  double             *cellCenter
)
{

  int i_rank;
  int n_rank;

  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  /* Local connectivity (Allocate in _prepareconnectivity */
  int sizeface_vtx_idx = 0;

  int *face_vtx;
  int *face_vtx_idx;

  /* Local coordinates */
  double *lvtx_coord;

  /*
   * Get face connectivities from other process
   */
  _prepare_connectivity( comm,
                         dn_cell,
                         dcell_face_idx,
                         dcell_face,
                         dface_vtx_idx,
                         dface_vtx,
                         dface_proc,
                        &face_vtx,
                        &face_vtx_idx,
                        &sizeface_vtx_idx,
                         dvtx_proc,
                         dvtx_coord,
                        &lvtx_coord);


  /* Verbose */
  if(0 == 1 ){
    PDM_printf("face_vtx : \n");
    for (int i1 = 0; i1 < sizeface_vtx_idx; i1++) {
      for (int j = face_vtx_idx[i1]; j < face_vtx_idx[i1+1]; j++)
        PDM_printf(" %i", face_vtx[j]);
      PDM_printf("\n");
    }
  }
  free (face_vtx);

  /*
   * Compute cell centers
   */
  int iFace = 0;
  for (int iCell = 0; iCell < dn_cell; iCell++) {

    int iBeg = dcell_face_idx[iCell  ];
    int iEnd = dcell_face_idx[iCell+1];

    cellCenter[3*iCell    ] = 0.;
    cellCenter[3*iCell + 1] = 0.;
    cellCenter[3*iCell + 2] = 0.;

    int n_vtxOn_cell = 0;

    /* Loop on all faces */
    for (int iFaceG = iBeg; iFaceG < iEnd; iFaceG++) {

      int iBegVtx = face_vtx_idx[iFace  ];
      int iEndVtx = face_vtx_idx[iFace+1];
      n_vtxOn_cell += iEndVtx - iBegVtx;

      /* Loop on all Vtx of current faces */
      for (int iVtx = iBegVtx; iVtx < iEndVtx; iVtx++){

        cellCenter[3*iCell    ] += lvtx_coord[3*iVtx    ];
        cellCenter[3*iCell + 1] += lvtx_coord[3*iVtx + 1];
        cellCenter[3*iCell + 2] += lvtx_coord[3*iVtx + 2];

      }

      /* Go to next face in local numerotation */
      iFace++;

    }

    /* Ponderate */
    cellCenter[3*iCell    ] = cellCenter[3*iCell    ]/n_vtxOn_cell;
    cellCenter[3*iCell + 1] = cellCenter[3*iCell + 1]/n_vtxOn_cell;
    cellCenter[3*iCell + 2] = cellCenter[3*iCell + 2]/n_vtxOn_cell;
  }
  free (lvtx_coord);
  free (face_vtx_idx);

  /* Verbose */
  if(0 == 1){
    for (int iCell = 0; iCell < dn_cell; iCell++) {
      PDM_printf(" --------- %i / %i \n ", iCell);
      PDM_printf(" %12.5e %12.5e %12.5e \n", cellCenter[3*iCell    ], cellCenter[3*iCell + 1], cellCenter[3*iCell + 2]);
    }
  }

}

/*=============================================================================
 * Public function definitions
 *============================================================================*/


/**
 *
 * \brief Perform geomtric patitionning
 *
 * \param [in]   method         Geometric method
 * \param [in]   n_part          Number of partition to build on this process
 * \param [in]   comm           Communicator
 * \param [in]   dn_cell         Number of distributed cells
 * \param [in]   dn_face         Number of distributed faces
 * \param [in]   dn_vtx          Number of distributed vertices
 * \param [in]   dcell_face_idx   Distributed cell face connectivity index or NULL
 *                              (size : dn_cell + 1, numbering : 0 to n-1)
 * \param [in]   dcell_face      Distributed cell face connectivity or NULL
 *                              (size : dface_vtx_idx[dn_cell], numbering : 1 to n)
 * \param [in]   dcell_weight    Cell weight (size : n_cell) or NULL
 * \param [in]   dface_vtx_idx    Distributed face to vertex connectivity index
 *                              (size : dn_face + 1, numbering : 0 to n-1)
 * \param [in]   dface_vtx       Distributed face to vertex connectivity
 *                              (size : dface_vtx_idx[dn_face], numbering : 1 to n)
 * \param [in]   dvtx_coord      Distributed vertex coordinates
 *                              (size : 3*dn_vtx)
 * \param [inout]   dcell_part      Distributed cell partitioning
 *                              (size = dn_cell)
 *
 */

void
PDM_part_geom
(
 PDM_part_geom_t     method,
 const int           n_part,
 const PDM_MPI_Comm  comm,
 const int           dn_cell,
 const int          *dcell_face_idx,
 const PDM_g_num_t  *dcell_face,
 const int          *dcell_weight,
 const int          *dface_vtx_idx,
 const PDM_g_num_t  *dface_vtx,
 const PDM_g_num_t  *dface_proc,
 const double        *dvtx_coord,
 const PDM_g_num_t  *dvtx_proc,
 int                *dcell_part
)
{

  assert (method == PDM_PART_GEOM_HILBERT);

  const int dim = 3;

  double *barycenterCoords = (double *) malloc (dn_cell * 3 * sizeof(double ));

  PDM_hilbert_code_t *hilbert_codes     = (PDM_hilbert_code_t *) malloc (dn_cell * sizeof(PDM_hilbert_code_t));
  PDM_hilbert_code_t *tmp_hilbert_codes = (PDM_hilbert_code_t *) malloc (dn_cell * sizeof(PDM_hilbert_code_t));

  /*
   * cell center computation
   */
  _compute_cell_center (comm,
                        dn_cell,
                        dcell_face_idx,
                        dcell_face,
                        dface_vtx_idx,
                        dface_vtx,
                        dface_proc,
                        dvtx_coord,
                        dvtx_proc,
                        barycenterCoords);


  /** TRAITEMENT HILBERT FVM **/

	/** Initialisation **/

  double extents[2*dim]; /** DIM x 2**/

	/** Get EXTENTS **/
  PDM_hilbert_get_coord_extents_par(dim, dn_cell, barycenterCoords, extents, comm);

	/** Hilbert Coordinates Computation **/
  PDM_hilbert_encode_coords(dim, PDM_HILBERT_CS, extents, dn_cell, barycenterCoords, hilbert_codes);

  for (int i = 0; i < dn_cell; ++i) {
    tmp_hilbert_codes [i] = hilbert_codes [i];
	}

  ///** Calcul des index des codes Hilbert **/

  int * hilbert_order = (int * ) malloc (dn_cell * sizeof(int));

  for (int i = 0; i < dn_cell; ++i) {
    hilbert_order [i] = i;
  }

  assert (sizeof(double) == sizeof(PDM_hilbert_code_t));
  PDM_sort_double (tmp_hilbert_codes, hilbert_order, dn_cell);

  free(tmp_hilbert_codes);

  int n_rank;
  PDM_MPI_Comm_size (comm, &n_rank);

  PDM_hilbert_code_t *hilbert_codesIdx = (PDM_hilbert_code_t *) malloc ((n_rank+1)*n_part * sizeof(PDM_hilbert_code_t));

  int * weight = (int *) malloc (dn_cell * sizeof(int));
  if (dcell_weight != NULL) {
    for(int i = 0; i < dn_cell; ++i) {
		  weight [i] = dcell_weight [i];
    }
  }
  else {
    for(int i = 0; i < dn_cell; ++i) {
		  weight [i] = 1;
    }
  }

  int n_total_part;
  PDM_MPI_Allreduce ((void *) &n_part, &n_total_part, 1, PDM_MPI_INT, PDM_MPI_SUM, comm);

  PDM_hilbert_build_rank_index (dim,
                                n_total_part,
                                dn_cell,
                                hilbert_codes,
                                weight,
                                hilbert_order,
                                hilbert_codesIdx,
                                comm);

  free(weight);

  /** Remplissage de cell_parts -> en fct des codes Hilbert **/

  for(int i = 0; i < dn_cell; ++i) {
    size_t quantile = PDM_hilbert_quantile_search(n_total_part,
                                                hilbert_codes[i],
                                                hilbert_codesIdx);
    dcell_part [i] = (int) quantile;

  }

  free(barycenterCoords);
  free(hilbert_codesIdx);
  free(hilbert_order);
  free(hilbert_codes);

}


#ifdef __cplusplus
}
#endif /* __cplusplus */

