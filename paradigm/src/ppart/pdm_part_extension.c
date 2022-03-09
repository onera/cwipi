/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <float.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_mpi.h"
#include "pdm.h"
#include "pdm_distant_neighbor.h"
#include "pdm_logging.h"
#include "pdm_unique.h"
#include "pdm_binary_search.h"
#include "pdm_order.h"
#include "pdm_error.h"
#include "pdm_part_extension.h"
#include "pdm_part_extension_priv.h"
#include "pdm_part_connectivity_transform.h"
#include "pdm_array.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type
 *============================================================================*/

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Static function definitions
 *============================================================================*/

static inline
int
_is_same_triplet
(
int iproc1, int ipart1, int ielt1,
int iproc2, int ipart2, int ielt2
)
{
  if(iproc1 == iproc2){
    if(ipart1 == ipart2){
      if(ielt1 == ielt2){
        return 1;
      }
    }
  }
  return 0;
}

static
int
_setup_unique_order_triplet
(
const int    n_entity,
const int   *unique_neighbor_entity_idx,
const int   *unique_neighbor_entity,
      int  **unique_order_neighbor_entity
)
{
  int* order        = malloc( unique_neighbor_entity_idx[n_entity] * sizeof(int));
  int* unique_order = malloc( unique_neighbor_entity_idx[n_entity] * sizeof(int));

  PDM_order_lnum_s(unique_neighbor_entity,
                   3,
                   order,
                   unique_neighbor_entity_idx[n_entity]);

  int idx_unique = -1;
  int last_proc  = -1;
  int last_part  = -1;
  int last_elmt  = -1;

  for(int i = 0; i < unique_neighbor_entity_idx[n_entity]; i++){

    int old_order = order[i];
    int curr_proc = unique_neighbor_entity[3*old_order  ];
    int curr_part = unique_neighbor_entity[3*old_order+1];
    int curr_elmt = unique_neighbor_entity[3*old_order+2];
    int is_same = _is_same_triplet(last_proc, last_part, last_elmt,
                                   curr_proc, curr_part, curr_elmt);
    // log_trace(" curr:: ( %d / %d / %d ) | last:: ( %d / %d / %d ) \n",
    //             curr_proc, curr_part, curr_elmt,
    //             last_proc, last_part, last_elmt);

    if(is_same == 0){ // N'est pas le meme
      idx_unique++;
      last_proc = curr_proc;
      last_part = curr_part;
      last_elmt = curr_elmt;
    }
    // order_unique[i] = idx_unique; // NO NO NO
    unique_order[old_order] = idx_unique; // NO NO NO
    // unique_order[old_order] = idx_unique;
    // log_trace("[%d] = %d --> %d \n", i, is_same, idx_unique);
    // new_order[order[i]] = i;
  }

  *unique_order_neighbor_entity = unique_order;
  free(order);
  return idx_unique+1;
}


static
void
_unique_triplet
(
  int   n_entity,
  int  *neighbor_entity_idx,
  int  *neighbor_entity,
  int **unique_neighbor_entity_idx,
  int **unique_neighbor_entity_n,
  int **unique_neighbor_entity
)
{

  int* _unique_neighbor_entity_idx = malloc( (n_entity + 1) * sizeof(int));
  int* _unique_neighbor_entity_n   = malloc( (n_entity    ) * sizeof(int));
  int* _unique_neighbor_entity     = malloc( 3 * neighbor_entity_idx[n_entity] * sizeof(int));
  int* order                       = malloc(     neighbor_entity_idx[n_entity] * sizeof(int)); // Suralloc

  _unique_neighbor_entity_idx[0] = 0;
  for(int i_entity = 0; i_entity < n_entity; ++i_entity) {

    int beg       = neighbor_entity_idx[i_entity];
    int n_connect = neighbor_entity_idx[i_entity+1] - beg;

    PDM_order_lnum_s(&neighbor_entity[3*beg], 3, order, n_connect);

    _unique_neighbor_entity_n  [i_entity  ] = 0;
    _unique_neighbor_entity_idx[i_entity+1] = _unique_neighbor_entity_idx[i_entity];

    int last_proc  = -1;
    int last_part  = -1;
    int last_elmt  = -1;
    for(int i = 0; i < n_connect; ++i) {
      int old_order   = order[i];
      int curr_proc   = neighbor_entity[3*(beg+old_order)  ];
      int curr_part   = neighbor_entity[3*(beg+old_order)+1];
      int curr_entity = neighbor_entity[3*(beg+old_order)+2];
      int is_same  = _is_same_triplet(last_proc, last_part, last_elmt,
                                      curr_proc, curr_part, curr_entity);

      if(is_same == 0){ // N'est pas le meme
        // idx_unique++;
        last_proc = curr_proc;
        last_part = curr_part;
        last_elmt = curr_entity;

        int beg_write = 3 * _unique_neighbor_entity_idx[i_entity+1];
        // printf("beg_write = %i | curr_proc = %i | curr_part = %i | curr_entity = %i \n", beg_write, curr_proc, curr_part, curr_entity);
        _unique_neighbor_entity[beg_write  ] = curr_proc;
        _unique_neighbor_entity[beg_write+1] = curr_part;
        _unique_neighbor_entity[beg_write+2] = curr_entity;

        /* Increment the new counter */
        _unique_neighbor_entity_idx[i_entity+1]++;
        _unique_neighbor_entity_n  [i_entity  ]++;
      }
    }
  }

  _unique_neighbor_entity = realloc(_unique_neighbor_entity, 3 * neighbor_entity_idx[n_entity] * sizeof(int));

  *unique_neighbor_entity_idx = _unique_neighbor_entity_idx;
  *unique_neighbor_entity_n   = _unique_neighbor_entity_n;
  *unique_neighbor_entity     = _unique_neighbor_entity;
  free(order);
}

static
void
_create_cell_cell_graph
(
  PDM_part_extension_t *part_ext,
  PDM_extend_type_t     extend_type
)
{
  // assert(extend_type == PDM_EXTEND_FROM_FACE);

  part_ext->cell_cell_idx = malloc( part_ext->n_domain * sizeof(int ***));
  part_ext->cell_cell     = malloc( part_ext->n_domain * sizeof(int ***));

  if(extend_type == PDM_EXTEND_FROM_FACE ){

    /* In order to call generic fonction we tranform face_cell with idx */
    int ***face_cell_idx = malloc( part_ext->n_domain * sizeof(int ***));
    int ***face_cell     = malloc( part_ext->n_domain * sizeof(int ***));

    int shift_part = 0;
    for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {

      int  *pn_face        = (int * ) malloc( part_ext->n_part[i_domain] * sizeof(int  ));
      int  *pn_cell        = (int * ) malloc( part_ext->n_part[i_domain] * sizeof(int  ));
      int **pface_cell     = (int **) malloc( part_ext->n_part[i_domain] * sizeof(int *));
      int **pcell_face     = (int **) malloc( part_ext->n_part[i_domain] * sizeof(int *));
      int **pcell_face_idx = (int **) malloc( part_ext->n_part[i_domain] * sizeof(int *));

      for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {
        pn_face       [i_part] = part_ext->parts[i_domain][i_part].n_face;
        pn_cell       [i_part] = part_ext->parts[i_domain][i_part].n_cell;
        pface_cell    [i_part] = part_ext->parts[i_domain][i_part].face_cell;
        pcell_face_idx[i_part] = part_ext->parts[i_domain][i_part].cell_face_idx;
        pcell_face    [i_part] = part_ext->parts[i_domain][i_part].cell_face;
      }

      PDM_part_connectivity_to_connectity_idx(part_ext->n_part[i_domain],
                                              pn_face,
                                              pface_cell,
                                              &face_cell_idx[i_domain],
                                              &face_cell[i_domain]);

      part_ext->cell_cell_idx[i_domain] = (int **) malloc( part_ext->n_part[i_domain] * sizeof(int *));
      part_ext->cell_cell    [i_domain] = (int **) malloc( part_ext->n_part[i_domain] * sizeof(int *));

      for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {

        // PDM_log_trace_array_int(face_cell_idx[i_domain][i_part], pn_face[i_part], "face_cell_idx::");
        // PDM_log_trace_array_int(face_cell    [i_domain][i_part], face_cell_idx[i_domain][i_part][pn_face[i_part]], "face_cell::");
        // printf("[%i] face_cell_idx -> \n ", i_part);
        // for(int i_face = 0; i_face < pn_face[i_part]; ++i_face) {
        //   printf(" (%i) = ", i_face);
        //   for(int idx_cell = face_cell_idx[i_domain][i_part][i_face]; idx_cell < face_cell_idx[i_domain][i_part][i_face+1]; ++idx_cell) {
        //     printf(" %i", face_cell[i_domain][i_part][idx_cell]-1);
        //   }
        //   printf("\n");
        // }

        // On fait égalemnt le cell_cell
        // printf("PDM_combine_connectivity[%i] -> %i %i \n", i_part, pn_cell[i_part], pn_face[i_part]);
        PDM_combine_connectivity(pn_cell[i_part],
                                 pcell_face_idx[i_part],
                                 pcell_face[i_part],
                                 face_cell_idx[i_domain][i_part],
                                 face_cell[i_domain][i_part],
                                 &part_ext->cell_cell_idx[i_domain][i_part],
                                 &part_ext->cell_cell[i_domain][i_part]);

        // Remove sign for cell_cell
        for(int i = 0; i < part_ext->cell_cell_idx[i_domain][i_part][pn_cell[i_part]]; ++i) {
          part_ext->cell_cell[i_domain][i_part][i] = PDM_ABS(part_ext->cell_cell[i_domain][i_part][i]);
        }

        /*
         * Setup shortcut and free useless memory
         */
        part_ext->entity_cell_idx[i_part+shift_part] = face_cell_idx[i_domain][i_part];
        part_ext->entity_cell    [i_part+shift_part] = face_cell    [i_domain][i_part];
        part_ext->entity_cell_n  [i_part+shift_part] = (int * ) malloc( pn_face[i_part] * sizeof(int));
        for(int i_face = 0; i_face < pn_face[i_part]; ++i_face) {
          int n_connected = part_ext->entity_cell_idx[i_part+shift_part][i_face+1] - part_ext->entity_cell_idx[i_part+shift_part][i_face];
          part_ext->entity_cell_n[i_part+shift_part][i_face] = n_connected;
        }
      }

      free(pn_face);
      free(pn_cell);
      free(pface_cell);
      free(pcell_face);
      free(pcell_face_idx);

      free(face_cell_idx[i_domain]);
      free(face_cell[i_domain]);

      shift_part += part_ext->n_part[i_domain];
    }

    free(face_cell_idx);
    free(face_cell);
  } else if(extend_type == PDM_EXTEND_FROM_VTX) {

    // Check

    // Compute cell_cell = cell_face + face_edge + edge_vtx then transpose

    int shift_part = 0;
    for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {

      part_ext->cell_cell_idx[i_domain] = (int **) malloc( part_ext->n_part[i_domain] * sizeof(int *));
      part_ext->cell_cell    [i_domain] = (int **) malloc( part_ext->n_part[i_domain] * sizeof(int *));

      for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {
        int pn_cell = part_ext->parts[i_domain][i_part].n_cell;
        int pn_edge = part_ext->parts[i_domain][i_part].n_edge;
        int pn_vtx  = part_ext->parts[i_domain][i_part].n_vtx;

        int* edge_vtx_idx = (int *) malloc( (pn_edge+1) * sizeof(int));
        edge_vtx_idx[0] = 0;
        for(int i_edge = 0; i_edge < pn_edge; ++i_edge) {
          edge_vtx_idx[i_edge+1] = edge_vtx_idx[i_edge] + 2;
        }

        /* Compute cell_edge */
        int* cell_edge_idx = NULL;
        int* cell_edge     = NULL;
        PDM_combine_connectivity(part_ext->parts[i_domain][i_part].n_cell,
                                 part_ext->parts[i_domain][i_part].cell_face_idx,
                                 part_ext->parts[i_domain][i_part].cell_face,
                                 part_ext->parts[i_domain][i_part].face_edge_idx,
                                 part_ext->parts[i_domain][i_part].face_edge,
                                 &cell_edge_idx,
                                 &cell_edge);

        /* Compute cell_vtx */
        int* cell_vtx_idx = NULL;
        int* cell_vtx     = NULL;
        PDM_combine_connectivity(part_ext->parts[i_domain][i_part].n_cell,
                                 cell_edge_idx,
                                 cell_edge,
                                 edge_vtx_idx,
                                 part_ext->parts[i_domain][i_part].edge_vtx,
                                 &cell_vtx_idx,
                                 &cell_vtx);
        free(cell_edge_idx);
        free(cell_edge);
        free(edge_vtx_idx);

        int *vtx_cell     = NULL;
        int *vtx_cell_idx = NULL;
        PDM_connectivity_transpose(part_ext->parts[i_domain][i_part].n_cell,
                                   part_ext->parts[i_domain][i_part].n_vtx,
                                   cell_vtx_idx,
                                   cell_vtx,
                                   &vtx_cell_idx,
                                   &vtx_cell);

        PDM_combine_connectivity(part_ext->parts[i_domain][i_part].n_cell,
                                 cell_vtx_idx,
                                 cell_vtx,
                                 vtx_cell_idx,
                                 vtx_cell,
                                 &part_ext->cell_cell_idx[i_domain][i_part],
                                 &part_ext->cell_cell[i_domain][i_part]);

        // Remove sign for cell_cell
        for(int i = 0; i < part_ext->cell_cell_idx[i_domain][i_part][pn_cell]; ++i) {
          part_ext->cell_cell[i_domain][i_part][i] = PDM_ABS(part_ext->cell_cell[i_domain][i_part][i]);
        }
        free(cell_vtx_idx);
        free(cell_vtx);

        /*
         * Setup shortcut and free useless memory
         */
        part_ext->entity_cell_idx[i_part+shift_part] = vtx_cell_idx;
        part_ext->entity_cell    [i_part+shift_part] = vtx_cell    ;
        part_ext->entity_cell_n  [i_part+shift_part] = (int * ) malloc( pn_vtx * sizeof(int));
        for(int i_vtx = 0; i_vtx < pn_vtx; ++i_vtx) {
          int n_connected = part_ext->entity_cell_idx[i_part+shift_part][i_vtx+1] - part_ext->entity_cell_idx[i_part+shift_part][i_vtx];
          part_ext->entity_cell_n[i_part+shift_part][i_vtx] = n_connected;
        }

      } /* End i_part */
      shift_part += part_ext->n_part[i_domain];
    }


  } else {
    abort();
  }
}

static
void
_create_cell_graph_comm
(
  PDM_part_extension_t *part_ext
)
{
  /*
   * The idea is to use the first communicator graph to exchange directly the entity_cell
   *     -> We deduced first from this graph the data necessary for distant_neigbor
   *     -> We apply exchange with the connectivity entitiy_cell
   *        The information in graph comm is usefull also for directly adressing the results
   */
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(part_ext->comm, &i_rank);
  PDM_MPI_Comm_size(part_ext->comm, &n_rank);

  assert(part_ext->n_domain == 1);
  /* Si multidomain on fait un shift et tt roule */

  int n_tot_all_domain = 0;
  int n_part_loc_all_domain = 0;
  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
    n_tot_all_domain      += part_ext->n_tot_part_by_domain[i_domain];
    n_part_loc_all_domain += part_ext->n_part[i_domain];
  }
  // printf(" n_tot_all_domain      = %i \n", n_tot_all_domain);
  // printf(" n_part_loc_all_domain = %i \n", n_part_loc_all_domain);

  /* We flat all partition */
  assert(part_ext->neighbor_idx   == NULL);
  assert(part_ext->neighbor_idx   == NULL);
  assert(part_ext->n_entity_bound == NULL);

  part_ext->neighbor_idx   = (int **) malloc( n_part_loc_all_domain * sizeof(int *) );
  part_ext->neighbor_desc  = (int **) malloc( n_part_loc_all_domain * sizeof(int *) );
  part_ext->n_entity_bound = (int  *) malloc( n_part_loc_all_domain * sizeof(int  ) );

  part_ext->entity_cell_opp_idx = (int  **) malloc( n_part_loc_all_domain * sizeof(int  *) );

  part_ext->dist_neighbor_cell_n              = (int **) malloc( n_part_loc_all_domain * sizeof(int *) );
  part_ext->dist_neighbor_cell_idx            = (int **) malloc( n_part_loc_all_domain * sizeof(int *) );
  part_ext->dist_neighbor_cell_desc           = (int **) malloc( n_part_loc_all_domain * sizeof(int *) );
  part_ext->unique_order_dist_neighbor_cell   = (int **) malloc( n_part_loc_all_domain * sizeof(int *) );
  part_ext->n_unique_order_dist_neighbor_cell = (int  *) malloc( n_part_loc_all_domain * sizeof(int  ) );

  // Begin with exchange by the connectivity the cell opposite number
  int shift_part = 0;
  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {

    int n_part_total = part_ext->n_tot_part_by_domain[i_domain];

    /* First loop to count */
    for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {

      int n_part_entity_bound_tot = -1;

      // printf(" i_part+shift_part = %i \n", i_part+shift_part);
      int* _entity_part_bound = NULL;
      if(part_ext->extend_type == PDM_EXTEND_FROM_FACE) {
        // part_ext->n_entity_bound[i_part+shift_part] = n_part_entity_bound_tot;
        n_part_entity_bound_tot = part_ext->parts[i_domain][i_part].face_part_bound_part_idx[n_part_total];
        part_ext->n_entity_bound[i_part+shift_part] = part_ext->parts[i_domain][i_part].n_face;
        _entity_part_bound = part_ext->parts[i_domain][i_part].face_part_bound;
      } else if (part_ext->extend_type == PDM_EXTEND_FROM_VTX){
        n_part_entity_bound_tot = part_ext->parts[i_domain][i_part].vtx_part_bound_part_idx[n_part_total];
        part_ext->n_entity_bound[i_part+shift_part] = part_ext->parts[i_domain][i_part].n_vtx;
        _entity_part_bound = part_ext->parts[i_domain][i_part].vtx_part_bound;
      } else {
        PDM_error(__FILE__, __LINE__, 0, "PDM_part_extension_compute wrong extend_type \n");
      }

      part_ext->neighbor_idx       [i_part+shift_part] = (int *) malloc(    (part_ext->n_entity_bound[i_part+shift_part]+1) * sizeof(int) );

      int n_cell = part_ext->parts[i_domain][i_part].n_cell;
      part_ext->n_cell[i_part+shift_part] = n_cell;

      int* _neighbor_idx  = part_ext->neighbor_idx [i_part+shift_part];

      int* _neighbor_n = PDM_array_zeros_int(part_ext->n_entity_bound[i_part+shift_part]);

      /* Just copy to remove the 4 indices to 3 indices */
      for(int idx_entity = 0; idx_entity < n_part_entity_bound_tot; ++idx_entity) {
        int i_entity = _entity_part_bound[4*idx_entity]-1;
        _neighbor_n[i_entity] += 1;
      }

      /* Compute index */
      _neighbor_idx[0] = 0;
      for(int i_entity = 0; i_entity < part_ext->n_entity_bound[i_part+shift_part]; ++i_entity) {
        _neighbor_idx[i_entity+1] = _neighbor_idx[i_entity] + _neighbor_n[i_entity];
        _neighbor_n[i_entity] = 0;
      }

      /* Ici il faut faire les raccords entre domaine ---> Count */
      part_ext->neighbor_desc[i_part+shift_part] = (int *) malloc( 3 * _neighbor_idx[part_ext->n_entity_bound[i_part+shift_part]] * sizeof(int) );
      int* _neighbor_desc = part_ext->neighbor_desc[i_part+shift_part];

      /* Second loop for fill */
      /* Just copy to remove the 4 indices to 3 indices */
      for(int idx_entity = 0; idx_entity < n_part_entity_bound_tot; ++idx_entity) {
        int i_entity = _entity_part_bound[4*idx_entity]-1;
        int idx_write = _neighbor_idx[i_entity] + _neighbor_n[i_entity]++;
        _neighbor_desc[3*idx_write  ] = _entity_part_bound[4*idx_entity+1];   // i_proc_opp;
        _neighbor_desc[3*idx_write+1] = _entity_part_bound[4*idx_entity+2]-1; // i_part_opp
        _neighbor_desc[3*idx_write+2] = _entity_part_bound[4*idx_entity+3]-1; // i_entity_opp
      }

      free(_neighbor_n);

      // PDM_log_trace_array_int(_neighbor_idx , part_ext->n_entity_bound[i_part+shift_part]+1, "_neighbor_idx::");
      // PDM_log_trace_array_int(_neighbor_desc , 3 * _neighbor_idx[part_ext->n_entity_bound[i_part+shift_part]], "_neighbor_desc::");

      /* Ici il faut faire les raccords entre domaine ---> Fill - Count _distant_neighbor_cell also */

    }

    shift_part += part_ext->n_part[i_domain];
  }


  /*
   * All partition of all domain is treated in the same time
   *    -> All partition are flatten along i_domain
   */
  PDM_distant_neighbor_t* dn = PDM_distant_neighbor_create(part_ext->comm,
                                                           n_part_loc_all_domain,
                                                           part_ext->n_entity_bound,
                                                           part_ext->neighbor_idx,
                                                           part_ext->neighbor_desc);


  PDM_distant_neighbor_exch(dn,
                            sizeof(int),
                            PDM_STRIDE_VAR_INTERLACED,
                            -1,
                            part_ext->entity_cell_n,
                  (void **) part_ext->entity_cell,
                           &part_ext->entity_cell_opp_n,
                 (void ***)&part_ext->entity_cell_opp);
  // PDM_distant_neighbor_exch_int(dn,
  //                           sizeof(int),
  //                           PDM_STRIDE_VAR,
  //                           -1,
  //                           part_ext->entity_cell_n,
  //                 (void **) part_ext->entity_cell,
  //                          &part_ext->entity_cell_opp_n,
  //                (void ***)&part_ext->entity_cell_opp);

  /* Compute idx */
  shift_part = 0;
  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
    for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {

      int n_cell = part_ext->parts[i_domain][i_part].n_cell;
      int* _neighbor_idx        = part_ext->neighbor_idx       [i_part+shift_part];
      int* _neighbor_desc       = part_ext->neighbor_desc      [i_part+shift_part];

      int* _entity_cell_opp_n   = part_ext->entity_cell_opp_n  [i_part+shift_part];
      int* _entity_cell_opp     = part_ext->entity_cell_opp    [i_part+shift_part];

      int* _entity_cell_idx = part_ext->entity_cell_idx[i_part+shift_part];
      int* _entity_cell     = part_ext->entity_cell    [i_part+shift_part];

      int n_entity = part_ext->n_entity_bound[i_part+shift_part]; // n_face / n_edge / n_vtx
      int n_elmt   = _neighbor_idx[n_entity];

      // int n_cell = part_ext->parts[i_domain][i_part].n_cell;
      part_ext->entity_cell_opp_idx   [i_part+shift_part] = PDM_array_new_idx_from_sizes_int(_entity_cell_opp_n, n_elmt);
      part_ext->dist_neighbor_cell_n  [i_part+shift_part] = PDM_array_zeros_int(n_cell);
      part_ext->dist_neighbor_cell_idx[i_part+shift_part] = (int *) malloc( (n_cell+1) * sizeof(int) );

      int* _entity_cell_opp_idx    = part_ext->entity_cell_opp_idx   [i_part+shift_part];
      int* _dist_neighbor_cell_n   = part_ext->dist_neighbor_cell_n  [i_part+shift_part];
      int* _dist_neighbor_cell_idx = part_ext->dist_neighbor_cell_idx[i_part+shift_part];


      // part_ext->border_cell_list[i_part+shift_part] = malloc( n_elmt * sizeof(int));
      part_ext->border_cell_list[i_part+shift_part] = malloc( n_cell * sizeof(int));
      // for(int i = 0; i < n_elmt; ++i) {
      //   part_ext->border_cell_list[i_part+shift_part][i] = -1000;
      // }

      /* For each border count number of opposite */
      for(int i_entity = 0; i_entity < n_entity; ++i_entity) {
        for(int idx_entity = _neighbor_idx[i_entity]; idx_entity < _neighbor_idx[i_entity+1]; ++idx_entity) {
          for(int idx_cell = _entity_cell_idx[i_entity]; idx_cell < _entity_cell_idx[i_entity+1]; ++idx_cell) {
            int i_cell = PDM_ABS(_entity_cell[idx_cell])-1;
            _dist_neighbor_cell_n[i_cell] += _entity_cell_opp_n[idx_entity];
          }
        }
      }


      // PDM_log_trace_array_int(_dist_neighbor_cell_n , n_cell, "_dist_neighbor_cell_n::");
      // PDM_log_trace_array_int(_entity_cell_opp_n , n_part_entity_bound_tot, "_entity_cell_opp_n::");
      // PDM_log_trace_array_int(_entity_cell_opp , n_part_entity_bound_tot, "_entity_cell_opp::");

      /* Ici il faut faire les raccords entre domaine ---> Count */

      /* Compute idx */
      int idx_indic = 0;
      _dist_neighbor_cell_idx[0] = 0;
      for(int i_cell = 0; i_cell < n_cell; ++i_cell) {
        _dist_neighbor_cell_idx[i_cell+1] = _dist_neighbor_cell_idx[i_cell] + _dist_neighbor_cell_n[i_cell];
        if(_dist_neighbor_cell_n[i_cell] > 0){
          part_ext->border_cell_list[i_part+shift_part][idx_indic++] = i_cell;
        }
      }
      part_ext->border_cell_list[i_part+shift_part] = realloc(part_ext->border_cell_list[i_part+shift_part], idx_indic * sizeof(int));
      /* Because cell can be associated twice */
      part_ext->n_cell_border[i_part+shift_part] = idx_indic;

      // printf(" idx_indic = %i\n", idx_indic);
      // PDM_log_trace_array_int(part_ext->border_cell_list[i_part+shift_part]  , idx_indic  , "border_cell_list::");

      /* Reset */
      PDM_array_reset_int(_dist_neighbor_cell_n, n_cell, 0);

      /* Allocate */
      part_ext->dist_neighbor_cell_desc[i_part+shift_part] = (int * ) malloc( 3 * _dist_neighbor_cell_idx[n_cell] * sizeof(int));
      int* _dist_neighbor_cell_desc = part_ext->dist_neighbor_cell_desc[i_part+shift_part];

      for(int i_entity = 0; i_entity < n_entity; ++i_entity) {
        for(int idx_entity = _neighbor_idx[i_entity]; idx_entity < _neighbor_idx[i_entity+1]; ++idx_entity) {
          int i_proc_opp   = _neighbor_desc[3*idx_entity  ];
          int i_part_opp   = _neighbor_desc[3*idx_entity+1];
          // int i_entity_opp = _neighbor_desc[3*idx_entity+2];
          for(int idx_cell = _entity_cell_idx[i_entity]; idx_cell < _entity_cell_idx[i_entity+1]; ++idx_cell) {
            int i_cell = PDM_ABS(_entity_cell[idx_cell])-1;

            for(int idx_cell_opp = _entity_cell_opp_idx[idx_entity]; idx_cell_opp < _entity_cell_opp_idx[idx_entity+1]; ++idx_cell_opp) {
              int idx_write = _dist_neighbor_cell_idx[i_cell] + _dist_neighbor_cell_n[i_cell]++;
              int i_opp_cell = PDM_ABS(_entity_cell_opp[idx_cell_opp])-1; // Commence à zero
              // printf("[%i] - _entity_cell_opp[%i] = %i | idx_write = %i  \n", i_part, idx_cell_opp,  _entity_cell_opp[idx_cell_opp], idx_write);
              _dist_neighbor_cell_desc[3*idx_write  ] = i_proc_opp;
              _dist_neighbor_cell_desc[3*idx_write+1] = i_part_opp;
              _dist_neighbor_cell_desc[3*idx_write+2] = i_opp_cell;
              // printf("[%i][%i] _dist_neighbor_cell[%i] = %i %i %i - i_entity = %i from ii_part_opp = %i and  i_entity_opp = %i \n",
              //        i_part, i_cell, idx_write, i_proc_opp, i_part_opp, i_opp_cell, i_entity, i_part_opp, i_entity_opp);
            }
          }
        }
      }

      int* _unique_dist_neighbor_cell_idx = NULL;
      int* _unique_dist_neighbor_cell_n   = NULL;
      int* _unique_dist_neighbor_cell     = NULL;
      _unique_triplet(n_cell,
                      _dist_neighbor_cell_idx,
                      _dist_neighbor_cell_desc,
                      &_unique_dist_neighbor_cell_idx,
                      &_unique_dist_neighbor_cell_n,
                      &_unique_dist_neighbor_cell);
      free(_dist_neighbor_cell_idx);
      free(_dist_neighbor_cell_desc);
      free(_dist_neighbor_cell_n);
      part_ext->dist_neighbor_cell_idx         [i_part+shift_part] = _unique_dist_neighbor_cell_idx;
      part_ext->dist_neighbor_cell_n           [i_part+shift_part] = _unique_dist_neighbor_cell_n;
      part_ext->dist_neighbor_cell_desc        [i_part+shift_part] = _unique_dist_neighbor_cell;

      int *_unique_order_dist_neighbor_cell = NULL;
      int n_unique = _setup_unique_order_triplet(n_cell,
                                                 _unique_dist_neighbor_cell_idx,
                                                 _unique_dist_neighbor_cell,
                                                 &_unique_order_dist_neighbor_cell);
      part_ext->unique_order_dist_neighbor_cell  [i_part+shift_part] = _unique_order_dist_neighbor_cell;
      part_ext->n_unique_order_dist_neighbor_cell[i_part+shift_part] = n_unique;

    }
    shift_part += part_ext->n_part[i_domain];
  }

  PDM_distant_neighbor_free(dn);

  /*
   * Panic verbose
   */
  if(0 == 1) {
    shift_part = 0;
    for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {

      for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {
        int* _neighbor_idx = part_ext->neighbor_idx                [i_part+shift_part];
        int n_elmt         = _neighbor_idx[part_ext->n_entity_bound[i_part+shift_part]];
        int n_data         = part_ext->entity_cell_opp_idx[i_part+shift_part][n_elmt];

        PDM_log_trace_array_int(part_ext->entity_cell_opp_n  [i_part+shift_part], n_elmt  , "entity_cell_opp_n::");
        PDM_log_trace_array_int(part_ext->entity_cell_opp_idx[i_part+shift_part], n_elmt+1, "entity_cell_opp_idx::");
        PDM_log_trace_array_int(part_ext->entity_cell_opp    [i_part+shift_part], n_data  , "entity_cell_opp::");

        int n_cell = part_ext->parts[i_domain][i_part].n_cell;
        int* _dist_neighbor_cell_idx  = part_ext->dist_neighbor_cell_idx [i_part+shift_part];
        int* _dist_neighbor_cell_desc = part_ext->dist_neighbor_cell_desc[i_part+shift_part];

        printf(" _dist_neighbor_cell_idx ---------------- \n");
        for(int i_cell = 0; i_cell < n_cell; ++i_cell) {
          printf("[%i] i_cell -> %i -->  ", i_part, i_cell);
          for(int idx = _dist_neighbor_cell_idx[i_cell]; idx < _dist_neighbor_cell_idx[i_cell+1]; ++idx) {
            printf("(%i, %i) ", _dist_neighbor_cell_desc[3*idx+1], _dist_neighbor_cell_desc[3*idx+2]);
          }
          printf("\n");
        }
        printf(" _dist_neighbor_cell_idx ---------------- END \n");

        PDM_log_trace_array_int(_dist_neighbor_cell_idx , n_cell+1                       , "_dist_neighbor_cell_idx::");
        PDM_log_trace_array_int(_dist_neighbor_cell_desc, 3 * _dist_neighbor_cell_idx[n_cell], "_dist_neighbor_cell_desc::");
      }
      shift_part += part_ext->n_part[i_domain];
    }
  }
  // exit(1);

  /*
   *   Pour la reconstruction des ghost cell
   *     -> On peut à mon avis faire une astuce pour faire 1 distant_neighbor avec
   *         "1 voisins = 1 rangées "
   *     -> Assez pratique car par partitione on va récupérer une liste de rangés également -_-
   *     -> Pb possible quand une cellule apparait dans une partition de rang 2 et dans une autre en rang 1
   */

}


static
void
_compute_dual_graph
(
  PDM_part_extension_t   *part_ext,
  PDM_distant_neighbor_t *dn,
  int                     i_depth
)
{
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(part_ext->comm, &i_rank);
  PDM_MPI_Comm_size(part_ext->comm, &n_rank);

  // printf("_compute_dual_graph : %i  \n", i_depth);

  assert(i_depth > 0);
  int** prev_cell_cell_extended_idx = part_ext->cell_cell_extended_idx[i_depth-1];
  int** prev_cell_cell_extended_n   = part_ext->cell_cell_extended_n  [i_depth-1];
  int** prev_cell_cell_extended     = part_ext->cell_cell_extended    [i_depth-1];

  // int** next_cell_cell_extended_idx = NULL;
  int** next_cell_cell_extended_n   = NULL;
  int** next_cell_cell_extended     = NULL;

  /*
   * Exchange of the previous rank
   */
  PDM_distant_neighbor_exch(dn,
                            3 * sizeof(int),
                            PDM_STRIDE_VAR_INTERLACED,
                            -1,
                            prev_cell_cell_extended_n,
                  (void **) prev_cell_cell_extended,
                           &next_cell_cell_extended_n,
                 (void ***)&next_cell_cell_extended);

  int shift_part = 0;
  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
    for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {

      int n_cell        = part_ext->n_cell       [i_part+shift_part];
      int n_cell_border = part_ext->n_cell_border[i_part+shift_part];

      int* _border_cell_cell_extended_n   = next_cell_cell_extended_n[i_part+shift_part];
      int* _border_cell_cell_extended     = next_cell_cell_extended  [i_part+shift_part];
      int* _border_cell_cell_extended_idx = (int * ) malloc( (n_cell_border+1) * sizeof(int) );

      if( 0 == 1) {
        printf("prev_cell_cell_extended :: --------------------- \n");
        for(int i_cell = 0; i_cell < n_cell; ++i_cell) {
          if( prev_cell_cell_extended_idx[i_part+shift_part][i_cell+1] > prev_cell_cell_extended_idx[i_part+shift_part][i_cell]){
            printf("[%i] i_cell -> %i -->  ", i_part, i_cell);
          }
          for(int idx = prev_cell_cell_extended_idx[i_part+shift_part][i_cell]; idx < prev_cell_cell_extended_idx[i_part+shift_part][i_cell+1]; ++idx) {
            printf("(%i, %i) ", prev_cell_cell_extended[i_part+shift_part][3*idx+1], prev_cell_cell_extended[i_part+shift_part][3*idx+2]);
          }
          if( prev_cell_cell_extended_idx[i_part+shift_part][i_cell+1] > prev_cell_cell_extended_idx[i_part+shift_part][i_cell]){
            printf("\n");
          }
        }
        printf("prev_cell_cell_extended :: --------------------- END \n");
      }

      // _border_cell_cell_extended_idx[0] = 0;
      // for(int i = 0; i < n_cell_border; ++i) {
      //   _border_cell_cell_extended_idx[i+1] = _border_cell_cell_extended_idx[i] + _border_cell_cell_extended_n[i];
      // }

      int lidx_write = 0;
      _border_cell_cell_extended_idx[0] = 0;
      for(int i_cell = 0; i_cell < n_cell; ++i_cell) {
        int n_neight = 0;
        for(int idx_neight = part_ext->dist_neighbor_cell_idx [i_part+shift_part][i_cell]; idx_neight < part_ext->dist_neighbor_cell_idx[i_part+shift_part][i_cell+1]; ++idx_neight){
          n_neight += _border_cell_cell_extended_n[idx_neight];
        }
        // printf(" beg = %i | end = %i | n_neight = %i \n", part_ext->dist_neighbor_cell_idx [i_part+shift_part][i_cell], part_ext->dist_neighbor_cell_idx [i_part+shift_part][i_cell+1], n_neight );
        if(n_neight > 0) {
          _border_cell_cell_extended_idx[lidx_write+1] = _border_cell_cell_extended_idx[lidx_write] + n_neight;
          lidx_write++;
        }
      }

      for(int i = 0; i < n_cell_border; ++i) {
        _border_cell_cell_extended_n[i] = _border_cell_cell_extended_idx[i+1] - _border_cell_cell_extended_idx[i];
      }

      int* _unique_border_cell_cell_extended_n   = NULL;
      int* _unique_border_cell_cell_extended     = NULL;
      int* _unique_border_cell_cell_extended_idx = NULL;
      _unique_triplet(n_cell_border,
                      _border_cell_cell_extended_idx,
                      _border_cell_cell_extended,
                      &_unique_border_cell_cell_extended_idx,
                      &_unique_border_cell_cell_extended_n,
                      &_unique_border_cell_cell_extended);

      free(_border_cell_cell_extended_n   );
      free(_border_cell_cell_extended     );
      free(_border_cell_cell_extended_idx );

      _border_cell_cell_extended_n   = _unique_border_cell_cell_extended_n  ;
      _border_cell_cell_extended     = _unique_border_cell_cell_extended    ;
      _border_cell_cell_extended_idx = _unique_border_cell_cell_extended_idx;
      next_cell_cell_extended_n[i_part+shift_part] = _border_cell_cell_extended_n;
      next_cell_cell_extended  [i_part+shift_part] = _border_cell_cell_extended;

      int *_unique_order_border_cell_cell_extended = NULL;
      int n_unique_border = _setup_unique_order_triplet(n_cell_border,
                                                        _unique_border_cell_cell_extended_idx,
                                                        _unique_border_cell_cell_extended,
                                                        &_unique_order_border_cell_cell_extended);

      if(0 == 1) {
        PDM_log_trace_array_int(_border_cell_cell_extended_idx, n_cell_border+1, "next_border_cell_cell_extended_idx::");
        PDM_log_trace_array_int(_border_cell_cell_extended_n  , n_cell_border,   "next_border_cell_cell_extended_n::");
        PDM_log_trace_array_int(_border_cell_cell_extended    , 3 * _border_cell_cell_extended_idx[n_cell_border], "next_border_cell_cell_extended::");

        printf("_border_cell_cell_extended :: --------------------- \n");
        for(int i = 0; i < n_cell_border; ++i) {
          int i_cell = part_ext->border_cell_list[i_part+shift_part][i];
          if( _border_cell_cell_extended_idx[i+1] > _border_cell_cell_extended_idx[i]){
            printf("i_cell -> %i -->  ", i_cell);
          }
          for(int idx = _border_cell_cell_extended_idx[i]; idx < _border_cell_cell_extended_idx[i+1]; ++idx) {
            printf("(%i, %i) ", _border_cell_cell_extended[3*idx+1], _border_cell_cell_extended[3*idx+2]);
          }
          if( _border_cell_cell_extended_idx[i+1] > _border_cell_cell_extended_idx[i]){
            printf("\n");
          }
        }
        printf("_border_cell_cell_extended :: --------------------- END \n");

      }

      /* Now we have the opposite border_cell_cell_extended
       *   -> We need to reput all together
       *   -> The current cell has to be also extendted but caution !!!
       *      the extension depend also of border
       */

      /*
       * Generate tag to know is a local cell is a border or not
       */
      int* idx_border_cell = PDM_array_const_int(n_cell, -1);

      for(int idx_cell = 0; idx_cell < n_cell_border; ++idx_cell) {
        int i_cell = part_ext->border_cell_list[i_part+shift_part][idx_cell];
        idx_border_cell[i_cell] = idx_cell; /* Keep idx_cell for imediate loop */
      }

      /* Allocate current depth and shorcut */
      part_ext->cell_cell_extended_idx[i_depth][i_part+shift_part] = (int *) malloc( (n_cell + 1 ) * sizeof(int));
      part_ext->cell_cell_extended_n  [i_depth][i_part+shift_part] = (int *) malloc( (n_cell     ) * sizeof(int));

      int* _prev_cell_cell_extended_n   = prev_cell_cell_extended_n  [i_part+shift_part];
      int* _prev_cell_cell_extended_idx = prev_cell_cell_extended_idx[i_part+shift_part];
      int* _prev_cell_cell_extended     = prev_cell_cell_extended    [i_part+shift_part];

      int* _cell_cell_extended_idx = part_ext->cell_cell_extended_idx[i_depth][i_part+shift_part];
      int* _cell_cell_extended_n   = part_ext->cell_cell_extended_n  [i_depth][i_part+shift_part];

      int i_depth_cur = 0;
      int* _cell_cell_idx = part_ext->cell_cell_idx[i_depth_cur][i_part+shift_part];
      int* _cell_cell     = part_ext->cell_cell    [i_depth_cur][i_part+shift_part];

      for(int i_cell = 0; i_cell < n_cell; ++i_cell) {
        _cell_cell_extended_idx[i_cell] = 0;
        _cell_cell_extended_n  [i_cell] = 0;
      }
      _cell_cell_extended_idx[n_cell] = 0;

      /* Weak hash table */
      int n_unique_dist_neighbor_cell = part_ext->n_unique_order_dist_neighbor_cell[i_part+shift_part];
      int n_unique_cell_cell_extended = part_ext->n_unique_order_cell_cell_extended[i_depth-1][i_part+shift_part];

      int *_unique_prev_cell_cell_extended = part_ext->unique_order_cell_cell_extended[i_depth-1][i_part+shift_part];

      // printf("n_unique_dist_neighbor_cell = %i \n", n_unique_dist_neighbor_cell);
      // printf("n_unique_cell_cell_extended = %i \n", n_unique_cell_cell_extended);

      int *tag_border_cell_cell_extended = malloc( n_unique_border             * sizeof(int));
      int *tag_dist_neighbor_cell        = malloc( n_unique_dist_neighbor_cell * sizeof(int));
      int *tag_prev_cell_cell_extended   = malloc( n_unique_cell_cell_extended * sizeof(int));
      int *tag_interior_cell             = malloc( n_cell                      * sizeof(int));

      int n_work = n_unique_border + n_unique_dist_neighbor_cell + n_unique_cell_cell_extended;
      int *icell_to_reset = malloc( n_work * sizeof(int)); // Suralloc but don't really have choice

      int *icell_border_to_reset = malloc( n_unique_border * sizeof(int)); // Suralloc but don't really have choice

      int n_work2 = _cell_cell_idx[n_cell];
      int *icell_to_reset_interior = malloc( n_work2 * sizeof(int)); // Suralloc but don't really have choice

      for(int i = 0; i < n_unique_border; ++i) {
        tag_border_cell_cell_extended[i] = 0;
      }
      for(int i = 0; i < n_unique_dist_neighbor_cell; ++i) {
        tag_dist_neighbor_cell[i] = 0;
      }
      for(int i = 0; i < n_unique_cell_cell_extended; ++i) {
        tag_prev_cell_cell_extended[i] = 0;
      }

      for(int i = 0; i < n_cell; ++i) {
        tag_interior_cell[i] = 0;
      }

      /* Let's go - First pass to count */
      for(int idx_cell = 0; idx_cell < n_cell_border; ++idx_cell) {
        int i_cell = part_ext->border_cell_list[i_part+shift_part][idx_cell];
        int n_unique_loc = 0;
        int n_interior_unique_loc = 0;
        int n_border_unique_loc = 0;

        assert(_cell_cell_extended_n[i_cell] == 0); // All are sorted before

        /* From interior - we add the previous rank */
        _cell_cell_extended_n[i_cell] = _prev_cell_cell_extended_n[i_cell];
        for(int idx_neight = _prev_cell_cell_extended_idx[i_cell]; idx_neight < _prev_cell_cell_extended_idx[i_cell+1]; ++idx_neight) {
          int i_unique2 = _unique_prev_cell_cell_extended[idx_neight]; // Pas la cellule le tableau est unique !!!!!
          tag_prev_cell_cell_extended[i_unique2 ] = 1;
          icell_to_reset[n_unique_loc++ ] = i_unique2;
        }

        /* From border */
        // printf("_prev_cell_cell_extended_n[%i][%i] = %i\n", i_part, i_cell, _prev_cell_cell_extended_n[i_cell]);
        // printf("_border_cell_cell_extended_n[%i][%i] = %i\n", i_part, idx_cell, _border_cell_cell_extended_n[idx_cell]);
        assert(_border_cell_cell_extended_n[idx_cell] > 0);
        // _cell_cell_extended_n[i_cell] += _border_cell_cell_extended_n[idx_cell];

        for(int idx_neight = _border_cell_cell_extended_idx[idx_cell]; idx_neight < _border_cell_cell_extended_idx[idx_cell+1]; ++idx_neight) {
          int i_unique2 = _unique_order_border_cell_cell_extended[idx_neight];
          int is_treat2 = tag_border_cell_cell_extended          [i_unique2 ];
          if(is_treat2 == 0) {
            _cell_cell_extended_n[i_cell] += 1;
            tag_border_cell_cell_extended[i_unique2]  = 1;
            icell_border_to_reset[n_border_unique_loc++] = i_unique2;
          }
        }

        /* Now we have to extend the interior */
        for(int idx_neight = _prev_cell_cell_extended_idx[i_cell]; idx_neight < _prev_cell_cell_extended_idx[i_cell+1]; ++idx_neight) {
          int i_rank_neight = _prev_cell_cell_extended[3*idx_neight  ];
          int i_part_neight = _prev_cell_cell_extended[3*idx_neight+1];
          /* We add stencil only if it's local */
          if(i_part+shift_part == i_part_neight && i_rank == i_rank_neight) {
            int i_cell_neight = _prev_cell_cell_extended[3*idx_neight+2];

            int i_unique = _unique_prev_cell_cell_extended[idx_neight]; // Pas la cellule le tableau est unique !!!!!
            // int is_treat = tag_prev_cell_cell_extended    [i_unique  ];

            // printf("i_cell_neight = %i | tag_prev_cell_cell_extended[%i] = %i \n", i_cell_neight, i_unique, is_treat);

            /* From interior */
            // _cell_cell_extended_n[i_cell] += _prev_cell_cell_extended_n[i_cell_neight];

            for(int idx_neight2 = _prev_cell_cell_extended_idx[i_cell_neight]; idx_neight2 < _prev_cell_cell_extended_idx[i_cell_neight+1]; ++idx_neight2) {

              int i_unique2 = _unique_prev_cell_cell_extended[idx_neight2]; // Pas la cellule le tableau est unique !!!!!
              // printf("[%i] [neighbor = %i] _unique_prev_cell_cell_extended[%i] = %i \n", i_cell, i_cell_neight, idx_neight2, i_unique2);
              int is_treat2 = tag_prev_cell_cell_extended    [i_unique2  ];
              if(is_treat2 == 0) {
                _cell_cell_extended_n      [i_cell   ] += 1;
                tag_prev_cell_cell_extended[i_unique2]  = 1;
                icell_to_reset[n_unique_loc++] = i_unique2;
              }
            }

            int idx_border_neight = idx_border_cell[i_cell_neight];
            if(idx_border_neight != -1) {
              // Il faut rajouter les voisins aussi
              // _cell_cell_extended_n[i_cell] += _border_cell_cell_extended_n[idx_border_neight];
              for(int idx_neight2 = _border_cell_cell_extended_idx[idx_border_neight]; idx_neight2 < _border_cell_cell_extended_idx[idx_border_neight+1]; ++idx_neight2) {
                int i_unique2 = _unique_order_border_cell_cell_extended[idx_neight2];
                int is_treat2 = tag_border_cell_cell_extended          [i_unique2 ];
                if(is_treat2 == 0) {
                  _cell_cell_extended_n[i_cell] += 1;
                  tag_border_cell_cell_extended[i_unique2]  = 1;
                  icell_border_to_reset[n_border_unique_loc++] = i_unique2;
                }
              }
            }

            /* Rajout du vrai intérieur */
            // _cell_cell_extended_n[i_cell] += _cell_cell_idx[i_cell_neight+1] - _cell_cell_idx[i_cell_neight];

            for(int idx_neight2 = _cell_cell_idx[i_cell_neight]; idx_neight2 < _cell_cell_idx[i_cell_neight+1]; ++idx_neight2 ) {
              int i_cell_neight2 = _cell_cell[idx_neight2];
              if(tag_interior_cell[i_cell_neight2-1] == 0) {

                _cell_cell_extended_n[i_cell] += 1;

                icell_to_reset_interior[n_interior_unique_loc++] = i_cell_neight2-1;
                tag_interior_cell[i_cell_neight2-1] = 1;
              }
            }

            tag_prev_cell_cell_extended[i_unique] = 1;

          } /* End if same part and same proc */
        } /* End loop neighbor */

        // Reset tag
        for(int j = 0; j < n_unique_loc; ++j) {
          tag_prev_cell_cell_extended[icell_to_reset[j]] = 0;
        }
        n_unique_loc = 0;
        for(int j = 0; j < n_interior_unique_loc; ++j) {
          tag_interior_cell[icell_to_reset_interior[j]] = 0;
        }
        n_interior_unique_loc = 0;
        for(int j = 0; j < n_border_unique_loc; ++j) {
          tag_border_cell_cell_extended[icell_border_to_reset[j]] = 0;
        }
        n_border_unique_loc = 0;


        /* Now we have to extend the interior */
        // for(int idx_neight = _prev_cell_cell_extended_idx[i_cell]; idx_neight < _prev_cell_cell_extended_idx[i_cell+1]; ++idx_neight) {
        //   int i_rank_neight = _prev_cell_cell_extended[3*idx_neight  ];
        //   int i_part_neight = _prev_cell_cell_extended[3*idx_neight+1];
        //   /* We add stencil only if it's local */
        //   if(i_part+shift_part == i_part_neight && i_rank == i_rank_neight) {
        //     int i_cell_neight = _prev_cell_cell_extended[3*idx_neight+2];

        //     int i_unique = _unique_prev_cell_cell_extended[i_cell_neight];
        //     int is_treat = tag_prev_cell_cell_extended    [i_unique];

        //     printf("i_cell_neight = %i | tag_prev_cell_cell_extended[%i] = %i \n", i_cell_neight, i_unique, is_treat);

        //     /* From interior */
        //     _cell_cell_extended_n[i_cell] += _prev_cell_cell_extended_n[i_cell_neight];

        //     int idx_border_neight = idx_border_cell[i_cell_neight];
        //     if(idx_border_neight != -1) {
        //       // Il faut rajouter les voisins aussi
        //       _cell_cell_extended_n[i_cell] += _border_cell_cell_extended_n[idx_border_neight];
        //     }

        //     /* Rajout du vrai intérieur */
        //     _cell_cell_extended_n[i_cell] += _cell_cell_idx[i_cell_neight+1] - _cell_cell_idx[i_cell_neight];

        //     tag_prev_cell_cell_extended[i_unique] = 1;

        //   } /* End if same part and same proc */
        // } /* End loop neighbor */

      } /* End loop border */

      /* Setup idx and reset */
      int max_neight = 0;
      _cell_cell_extended_idx[0] = 0;
      for(int i_cell = 0; i_cell < n_cell; ++i_cell) {
        _cell_cell_extended_idx[i_cell+1] = _cell_cell_extended_idx[i_cell] + _cell_cell_extended_n[i_cell];
        max_neight = PDM_MAX(max_neight, _cell_cell_extended_n[i_cell]);
        _cell_cell_extended_n[i_cell] = 0;
      }

      /* Allocate */
      part_ext->cell_cell_extended[i_depth][i_part+shift_part] = (int *) malloc( 3 * _cell_cell_extended_idx[n_cell] * sizeof(int));
      int* _cell_cell_extended = part_ext->cell_cell_extended[i_depth][i_part+shift_part];


      /* Let's go - Second pass to fill */
      for(int idx_cell = 0; idx_cell < n_cell_border; ++idx_cell) {
        int i_cell = part_ext->border_cell_list[i_part+shift_part][idx_cell];
        int n_unique_loc = 0;
        int n_interior_unique_loc = 0;
        int n_border_unique_loc = 0;
        assert(_cell_cell_extended_n[i_cell] == 0); // All are sorted before

        /* From interior - we add the previous rank */
        for(int idx_neight = _prev_cell_cell_extended_idx[i_cell]; idx_neight < _prev_cell_cell_extended_idx[i_cell+1]; ++idx_neight) {
          int idx_write = _cell_cell_extended_idx[i_cell] + _cell_cell_extended_n[i_cell]++;
          _cell_cell_extended[3*idx_write  ] = _prev_cell_cell_extended[3*idx_neight  ];
          _cell_cell_extended[3*idx_write+1] = _prev_cell_cell_extended[3*idx_neight+1];
          _cell_cell_extended[3*idx_write+2] = _prev_cell_cell_extended[3*idx_neight+2];

          int i_unique2 = _unique_prev_cell_cell_extended[idx_neight]; // Pas la cellule le tableau est unique !!!!!
          tag_prev_cell_cell_extended[i_unique2 ] = 1;
          icell_to_reset[n_unique_loc++] = i_unique2;
        }

        /* From border */
        assert(_border_cell_cell_extended_n[idx_cell] > 0);
        for(int idx_neight = _border_cell_cell_extended_idx[idx_cell]; idx_neight < _border_cell_cell_extended_idx[idx_cell+1]; ++idx_neight) {
          int i_unique2 = _unique_order_border_cell_cell_extended[idx_neight];
          int is_treat2 = tag_border_cell_cell_extended          [i_unique2 ];
          if(is_treat2 == 0) {

            int idx_write = _cell_cell_extended_idx[i_cell] + _cell_cell_extended_n[i_cell]++;
            _cell_cell_extended[3*idx_write  ] = _border_cell_cell_extended[3*idx_neight  ];
            _cell_cell_extended[3*idx_write+1] = _border_cell_cell_extended[3*idx_neight+1];
            _cell_cell_extended[3*idx_write+2] = _border_cell_cell_extended[3*idx_neight+2];

            tag_border_cell_cell_extended[i_unique2]  = 1;
            icell_border_to_reset[n_border_unique_loc++] = i_unique2;
          }
        }

        /* Now we have to extend the interior */
        for(int idx_neight = _prev_cell_cell_extended_idx[i_cell]; idx_neight < _prev_cell_cell_extended_idx[i_cell+1]; ++idx_neight) {
          int i_rank_neight = _prev_cell_cell_extended[3*idx_neight  ];
          int i_part_neight = _prev_cell_cell_extended[3*idx_neight+1];
          /* We add stencil only if it's local */
          if(i_part+shift_part == i_part_neight && i_rank == i_rank_neight) {

            int i_cell_neight = _prev_cell_cell_extended[3*idx_neight+2];

            /* From interior */
            for(int idx_neight2 = _prev_cell_cell_extended_idx[i_cell_neight]; idx_neight2 < _prev_cell_cell_extended_idx[i_cell_neight+1]; ++idx_neight2) {

              int i_unique2 = _unique_prev_cell_cell_extended[idx_neight2]; // Pas la cellule le tableau est unique !!!!!
              int is_treat2 = tag_prev_cell_cell_extended    [i_unique2  ];
              if(is_treat2 == 0) {
                int idx_write = _cell_cell_extended_idx[i_cell] + _cell_cell_extended_n[i_cell]++;
                _cell_cell_extended[3*idx_write  ] = _prev_cell_cell_extended[3*idx_neight2  ];
                _cell_cell_extended[3*idx_write+1] = _prev_cell_cell_extended[3*idx_neight2+1];
                _cell_cell_extended[3*idx_write+2] = _prev_cell_cell_extended[3*idx_neight2+2];

                // Update tag
                tag_prev_cell_cell_extended[i_unique2]  = 1;
                icell_to_reset[n_unique_loc++] = i_unique2;
              }
            }

            int idx_border_neight = idx_border_cell[i_cell_neight];
            if(idx_border_neight != -1) {
              // Il faut rajouter les voisins aussi
              for(int idx_neight2 = _border_cell_cell_extended_idx[idx_border_neight]; idx_neight2 < _border_cell_cell_extended_idx[idx_border_neight+1]; ++idx_neight2) {

                int i_unique2 = _unique_order_border_cell_cell_extended[idx_neight2];
                int is_treat2 = tag_border_cell_cell_extended          [i_unique2 ];
                if(is_treat2 == 0) {
                  int idx_write = _cell_cell_extended_idx[i_cell] + _cell_cell_extended_n[i_cell]++;
                  _cell_cell_extended[3*idx_write  ] = _border_cell_cell_extended[3*idx_neight2  ];
                  _cell_cell_extended[3*idx_write+1] = _border_cell_cell_extended[3*idx_neight2+1];
                  _cell_cell_extended[3*idx_write+2] = _border_cell_cell_extended[3*idx_neight2+2];

                  tag_border_cell_cell_extended[i_unique2]  = 1;
                  icell_border_to_reset[n_border_unique_loc++] = i_unique2;

                }
              }
            }

            /* Rajout du vrai intérieur */
            for(int idx_neight2 = _cell_cell_idx[i_cell_neight]; idx_neight2 < _cell_cell_idx[i_cell_neight+1]; ++idx_neight2 ) {
              int i_cell_neight2 = _cell_cell[idx_neight2];

              if(tag_interior_cell[i_cell_neight2-1] == 0) {

                int idx_write =  _cell_cell_extended_idx[i_cell] + _cell_cell_extended_n[i_cell]++;
                _cell_cell_extended[3*idx_write  ] = i_rank;
                _cell_cell_extended[3*idx_write+1] = i_part;
                _cell_cell_extended[3*idx_write+2] = i_cell_neight2-1;

                icell_to_reset_interior[n_interior_unique_loc++] = i_cell_neight2-1;
                tag_interior_cell[i_cell_neight2-1] = 1;
              }
            }

          } /* End if same part and same proc */

        } /* End loop neighbor */

        // Reset tag for next cell
        for(int j = 0; j < n_unique_loc; ++j) {
          tag_prev_cell_cell_extended[icell_to_reset[j]] = 0;
        }
        for(int j = 0; j < n_interior_unique_loc; ++j) {
          tag_interior_cell[icell_to_reset_interior[j]] = 0;
        }
        for(int j = 0; j < n_border_unique_loc; ++j) {
          tag_border_cell_cell_extended[icell_border_to_reset[j]] = 0;
        }

      } /* End loop border */
      // free(_tag_cell_is_treated);
      free(tag_border_cell_cell_extended);
      free(tag_dist_neighbor_cell);
      free(tag_prev_cell_cell_extended);
      free(tag_interior_cell);
      free(icell_to_reset);
      free(icell_border_to_reset);
      free(icell_to_reset_interior);

      /* The _cell_cell_extended need to be sorted because many entry is duplicated */
      int* order = malloc( max_neight * sizeof(int));
      int* _ncell_cell_extended_n   = malloc( (n_cell    ) * sizeof(int));
      int* _ncell_cell_extended_idx = malloc( (n_cell + 1) * sizeof(int));
      int* _ncell_cell_extended     = malloc( 3 * _cell_cell_extended_idx[n_cell] * sizeof(int));
      _ncell_cell_extended_idx[0] = 0;
      for(int i_cell = 0; i_cell < n_cell; ++i_cell) {

        int beg = _cell_cell_extended_idx[i_cell];
        int n_connect = _cell_cell_extended_idx[i_cell+1] - beg;
        assert(n_connect == _cell_cell_extended_n[i_cell]);

        PDM_order_lnum_s(&_cell_cell_extended[3*beg], 3, order, n_connect);

        // log_trace("i_cell = %i | n_cell = %i | n_connect = %i \n", i_cell, n_cell, n_connect);
        // PDM_log_trace_array_int(order, )
        // printf(" order[%i] = ", i_cell);
        // for(int i = 0; i < n_connect; ++i) {
        //   printf(" %i", order[i]);
        // }
        // printf("\n");

        _ncell_cell_extended_n  [i_cell  ] = 0;
        _ncell_cell_extended_idx[i_cell+1] = _ncell_cell_extended_idx[i_cell];

        // int idx_unique = -1;
        int last_proc  = -1;
        int last_part  = -1;
        int last_elmt  = -1;
        for(int i = 0; i < n_connect; ++i) {
          int old_order = order[i];
          int curr_proc = _cell_cell_extended[3*(beg+old_order)  ];
          int curr_part = _cell_cell_extended[3*(beg+old_order)+1];
          int curr_cell = _cell_cell_extended[3*(beg+old_order)+2];
          int is_same  = _is_same_triplet(last_proc, last_part, last_elmt,
                                          curr_proc, curr_part, curr_cell);

          if(is_same == 0){ // N'est pas le meme
            // idx_unique++;
            last_proc = curr_proc;
            last_part = curr_part;
            last_elmt = curr_cell;

            int beg_write = 3 * _ncell_cell_extended_idx[i_cell+1];
            // printf(" write in = %i  | beg_write = %i | idx_unique = %i\n", beg_write + idx_unique, beg_write, idx_unique);
            _ncell_cell_extended[beg_write  ] = curr_proc;
            _ncell_cell_extended[beg_write+1] = curr_part;
            _ncell_cell_extended[beg_write+2] = curr_cell;

            /* Increment the new counter */
            _ncell_cell_extended_idx[i_cell+1]++;
            _ncell_cell_extended_n  [i_cell  ]++;
          }
        }
      }

      /* Free old ptr and assign the sort one */
      free(part_ext->cell_cell_extended_idx[i_depth][i_part+shift_part]);
      free(part_ext->cell_cell_extended_n  [i_depth][i_part+shift_part]);
      free(part_ext->cell_cell_extended    [i_depth][i_part+shift_part]);

      _ncell_cell_extended = realloc(_ncell_cell_extended, 3 * _ncell_cell_extended_idx[n_cell] * sizeof(int));

      part_ext->cell_cell_extended_idx[i_depth][i_part+shift_part] = _ncell_cell_extended_idx;
      part_ext->cell_cell_extended_n  [i_depth][i_part+shift_part] = _ncell_cell_extended_n;
      part_ext->cell_cell_extended    [i_depth][i_part+shift_part] = _ncell_cell_extended;

      int *_unique_order_cell_cell_extended = NULL;
      int n_unique = _setup_unique_order_triplet(n_cell,
                                                 _ncell_cell_extended_idx,
                                                 _ncell_cell_extended,
                                                 &_unique_order_cell_cell_extended);
      part_ext->unique_order_cell_cell_extended  [i_depth][i_part+shift_part] = _unique_order_cell_cell_extended;
      part_ext->n_unique_order_cell_cell_extended[i_depth][i_part+shift_part] = n_unique;

      if(0 == 1) {
        PDM_log_trace_array_int(_ncell_cell_extended_idx, n_cell+1, "_ncell_cell_extended_idx:: ");
        PDM_log_trace_array_int(_ncell_cell_extended_n  , n_cell  , "_ncell_cell_extended_n:: ");
        PDM_log_trace_array_int(_ncell_cell_extended    , 3 * _ncell_cell_extended_idx[n_cell], "_ncell_cell_extended:: ");

        printf(" --------------------------------------------------------------- \n");
        printf("_ncell_cell_extended :: --------------------- \n");
        for(int i_cell = 0; i_cell < n_cell; ++i_cell) {
          if( _ncell_cell_extended_idx[i_cell+1] > _ncell_cell_extended_idx[i_cell]){
            printf("[%i] i_cell -> %i -->  ", i_part, i_cell);
          }
          for(int idx = _ncell_cell_extended_idx[i_cell]; idx < _ncell_cell_extended_idx[i_cell+1]; ++idx) {
            printf("(%i, %i) ", _ncell_cell_extended[3*idx+1], _ncell_cell_extended[3*idx+2]);
          }
          if( _ncell_cell_extended_idx[i_cell+1] > _ncell_cell_extended_idx[i_cell]){
            printf("\n");
          }
        }
        printf("_ncell_cell_extended :: --------------------- END \n");
        printf(" --------------------------------------------------------------- \n");
      }

      /* Free */
      free(idx_border_cell);
      free(_border_cell_cell_extended_idx);
      free(order);

      /* Free allocated memory in distant neigbor exhange */
      free(next_cell_cell_extended_n[i_part+shift_part]);
      free(next_cell_cell_extended  [i_part+shift_part]);
      free(_unique_order_border_cell_cell_extended);

    }
    shift_part += part_ext->n_part[i_domain];
  }

  free(next_cell_cell_extended_n);
  free(next_cell_cell_extended);

}



static
void
_compute_first_extended_cell_graph
(
 PDM_part_extension_t *part_ext
)
{
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(part_ext->comm, &i_rank);
  PDM_MPI_Comm_size(part_ext->comm, &n_rank);

  int i_depth_cur = 0;

  int shift_part = 0;
  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
    for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {

      int n_cell = part_ext->n_cell[i_part+shift_part];
      part_ext->cell_cell_extended_idx[i_depth_cur][i_part+shift_part] = (int *) malloc( (n_cell + 1 ) * sizeof(int));
      part_ext->cell_cell_extended_n  [i_depth_cur][i_part+shift_part] = (int *) malloc( (n_cell     ) * sizeof(int));

      int* _dist_neighbor_cell_n    = part_ext->dist_neighbor_cell_n   [i_part+shift_part];
      int* _dist_neighbor_cell_idx  = part_ext->dist_neighbor_cell_idx [i_part+shift_part];
      int* _dist_neighbor_cell_desc = part_ext->dist_neighbor_cell_desc[i_part+shift_part];

      /* Uniquement besoin du graph sur les cellules de bords */
      // int n_cell_border = _dist_neighbor_cell_idx[n_cell];
      int n_cell_border = part_ext->n_cell_border[i_part+shift_part];

      int* _cell_cell_extended_idx = part_ext->cell_cell_extended_idx[i_depth_cur][i_part+shift_part];
      int* _cell_cell_extended_n   = part_ext->cell_cell_extended_n  [i_depth_cur][i_part+shift_part];

      int* _cell_cell_idx = part_ext->cell_cell_idx[i_depth_cur][i_part+shift_part];
      int* _cell_cell     = part_ext->cell_cell    [i_depth_cur][i_part+shift_part];

      for(int i_cell = 0; i_cell < n_cell; ++i_cell) {
        _cell_cell_extended_idx[i_cell] = 0;
        _cell_cell_extended_n  [i_cell] = 0;
      }
      _cell_cell_extended_idx[n_cell] = 0;

      // printf(" n_cell_border = %i \n", n_cell_border);

      /* First pass to count */
      for(int idx_cell = 0; idx_cell < n_cell_border; ++idx_cell) {
        int i_cell = part_ext->border_cell_list[i_part+shift_part][idx_cell];
        assert(_cell_cell_extended_n[i_cell] == 0); // All are sorted before

        /* From interior */
        _cell_cell_extended_n[i_cell] = _cell_cell_idx[i_cell+1] - _cell_cell_idx[i_cell];

        /* From border */
        assert(_dist_neighbor_cell_n[i_cell] > 0);
        _cell_cell_extended_n[i_cell] += _dist_neighbor_cell_n[i_cell];
      }

      _cell_cell_extended_idx[0] = 0;
      for(int i_cell = 0; i_cell < n_cell; ++i_cell) {
        _cell_cell_extended_idx[i_cell+1] = _cell_cell_extended_idx[i_cell] + _cell_cell_extended_n[i_cell];
        _cell_cell_extended_n[i_cell] = 0;
      }

      part_ext->cell_cell_extended[i_depth_cur][i_part+shift_part] = (int *) malloc( 3 * _cell_cell_extended_idx[n_cell] * sizeof(int));
      int* _cell_cell_extended = part_ext->cell_cell_extended[i_depth_cur][i_part+shift_part];

      /* Second pass to fill */
      for(int idx_cell = 0; idx_cell < n_cell_border; ++idx_cell) {
        int i_cell = part_ext->border_cell_list[i_part+shift_part][idx_cell];

        /* From interior */
        for(int idx_neight = _cell_cell_idx[i_cell]; idx_neight < _cell_cell_idx[i_cell+1]; ++idx_neight ) {
          int i_cell_neight = _cell_cell[idx_neight];
          int idx_write =  _cell_cell_extended_idx[i_cell] + _cell_cell_extended_n[i_cell]++;
          _cell_cell_extended[3*idx_write  ] = i_rank;
          _cell_cell_extended[3*idx_write+1] = i_part;
          _cell_cell_extended[3*idx_write+2] = i_cell_neight-1;
        }

        /* From border */
        for(int idx_neight = _dist_neighbor_cell_idx[i_cell]; idx_neight < _dist_neighbor_cell_idx[i_cell+1]; ++idx_neight) {
          int i_rank_neight = _dist_neighbor_cell_desc[3*idx_neight  ];
          int i_part_neight = _dist_neighbor_cell_desc[3*idx_neight+1];
          int i_cell_neight = _dist_neighbor_cell_desc[3*idx_neight+2];
          int idx_write =  _cell_cell_extended_idx[i_cell] + _cell_cell_extended_n[i_cell]++;
          _cell_cell_extended[3*idx_write  ] = i_rank_neight;
          _cell_cell_extended[3*idx_write+1] = i_part_neight;
          _cell_cell_extended[3*idx_write+2] = i_cell_neight;
        }

        // printf("[%i] _cell_cell_extended_n[%i] = %i\n", i_part, i_cell, _cell_cell_extended_n[i_cell]);
      }

      int* _unique_cell_cell_extended_idx = NULL;
      int* _unique_cell_cell_extended_n   = NULL;
      int* _unique_cell_cell_extended     = NULL;

      /* Sort and unique inside connectivity */
      _unique_triplet(n_cell,
                      _cell_cell_extended_idx,
                      _cell_cell_extended,
                      &_unique_cell_cell_extended_idx,
                      &_unique_cell_cell_extended_n,
                      &_unique_cell_cell_extended);

      free(part_ext->cell_cell_extended_idx[i_depth_cur][i_part+shift_part]);
      free(part_ext->cell_cell_extended_n  [i_depth_cur][i_part+shift_part]);
      free(part_ext->cell_cell_extended    [i_depth_cur][i_part+shift_part]);

      part_ext->cell_cell_extended_idx[i_depth_cur][i_part+shift_part] = _unique_cell_cell_extended_idx;
      part_ext->cell_cell_extended_n  [i_depth_cur][i_part+shift_part] = _unique_cell_cell_extended_n;
      part_ext->cell_cell_extended    [i_depth_cur][i_part+shift_part] = _unique_cell_cell_extended;

      /*
       * Setup a unique among all connecetvity to avoid adding same entity multiple times
       */
      int *_unique_order_cell_cell_extended = NULL;
      int n_unique = _setup_unique_order_triplet(n_cell,
                                                 _unique_cell_cell_extended_idx,
                                                 _unique_cell_cell_extended,
                                                 &_unique_order_cell_cell_extended);

      part_ext->unique_order_cell_cell_extended  [i_depth_cur][i_part+shift_part] = _unique_order_cell_cell_extended;
      part_ext->n_unique_order_cell_cell_extended[i_depth_cur][i_part+shift_part] = n_unique;

      // PDM_log_trace_array_int(_unique_order_cell_cell_extended, _unique_cell_cell_extended_idx[n_cell]  , "_unique_order_cell_cell_extended::");
      if(0 == 1) {
        // PDM_log_trace_array_int(_cell_cell_extended_n  , n_cell  , "t_cell_cell_extended_n::");
        // PDM_log_trace_array_int(_cell_cell_extended_idx, n_cell+1, "t_cell_cell_extended_idx::");
        // PDM_log_trace_array_int(_cell_cell_extended, 3 * _cell_cell_extended_idx[n_cell]  , "t_cell_cell_extended::");
        PDM_log_trace_array_int(part_ext->cell_cell_extended_n  [i_depth_cur][i_part+shift_part], n_cell  , "t_cell_cell_extended_n::");
        PDM_log_trace_array_int(part_ext->cell_cell_extended_idx[i_depth_cur][i_part+shift_part], n_cell+1, "t_cell_cell_extended_idx::");
        PDM_log_trace_array_int(part_ext->cell_cell_extended    [i_depth_cur][i_part+shift_part], 3 * _cell_cell_extended_idx[n_cell]  , "t_cell_cell_extended::");
      }

    }
    shift_part += part_ext->n_part[i_domain];
  }
}


static
void
_prune_cell_cell_extented
(
  PDM_part_extension_t *part_ext,
  int i_depth
)
{
  // printf("_prune_cell_cell_extented \n");
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(part_ext->comm, &i_rank);
  PDM_MPI_Comm_size(part_ext->comm, &n_rank);

  int shift_part = 0;
  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
    for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {

      int* _cell_cell_extended_idx = part_ext->cell_cell_extended_idx[i_depth][i_part+shift_part];
      int* _cell_cell_extended     = part_ext->cell_cell_extended    [i_depth][i_part+shift_part];

      int n_cell      = part_ext->parts[i_domain][i_part].n_cell;
      int s_tot       = _cell_cell_extended_idx[n_cell];

      if( 1 == 0) {
        PDM_log_trace_array_int(_cell_cell_extended_idx, n_cell+1, "_cell_cell_extended_idx::");
        PDM_log_trace_array_int(_cell_cell_extended    , 3 * _cell_cell_extended_idx[n_cell], "_cell_cell_extended::");
      }

      int* order = (int * ) malloc( s_tot * sizeof(int));
      PDM_order_lnum_s(_cell_cell_extended, 3, order, s_tot);

      part_ext->cell_cell_extended_pruned    [i_part+shift_part] = (int * ) malloc(  3 * s_tot * sizeof(int));
      part_ext->cell_cell_extended_pruned_idx[i_part+shift_part] = (int * ) malloc( (n_cell+1) * sizeof(int));

      int* _cell_cell_extended_pruned     = part_ext->cell_cell_extended_pruned    [i_part+shift_part];
      int* _cell_cell_extended_pruned_idx = part_ext->cell_cell_extended_pruned_idx[i_part+shift_part];

      int idx_unique = 0;
      int last_proc  = -1;
      int last_part  = -1;
      int last_elmt  = -1;
      _cell_cell_extended_pruned_idx[0] = 0;
      _cell_cell_extended_pruned_idx[1] = 0;
      for(int i = 0; i < s_tot; ++i) {
        int old_order = order[i];
        int curr_proc = _cell_cell_extended[3*old_order  ];
        int curr_part = _cell_cell_extended[3*old_order+1];
        int curr_cell = _cell_cell_extended[3*old_order+2];
        int is_same  = _is_same_triplet(last_proc, last_part, last_elmt,
                                        curr_proc, curr_part, curr_cell);

        // On peut également trie les locaux qui ne serve à rien
        int is_local = (curr_proc == i_rank) && (curr_part == i_part+shift_part);
        if(is_same == 0 && !is_local){ // N'est pas le meme
          // idx_unique++;
          last_proc = curr_proc;
          last_part = curr_part;
          last_elmt = curr_cell;

          _cell_cell_extended_pruned[idx_unique++] = curr_proc;
          _cell_cell_extended_pruned[idx_unique++] = curr_part;
          _cell_cell_extended_pruned[idx_unique++] = curr_cell;

          /* Increment the new counter */
          _cell_cell_extended_pruned_idx[1]++;
        }
      }

      /* On considère que une cellule est connecté a toutes les autres */
      for(int i = 1; i < n_cell; ++i) {
        _cell_cell_extended_pruned_idx[i+1] = _cell_cell_extended_pruned_idx[i];
      }

      part_ext->cell_cell_extended_pruned[i_part+shift_part] = realloc(part_ext->cell_cell_extended_pruned[i_part+shift_part], idx_unique * sizeof(int));

      if(0 == 1) {
        _cell_cell_extended_pruned     = part_ext->cell_cell_extended_pruned    [i_part+shift_part];
        PDM_log_trace_array_int(_cell_cell_extended_pruned_idx, n_cell+1, "_cell_cell_extended_pruned_idx:: ");
        PDM_log_trace_array_int(_cell_cell_extended_pruned, 3 * _cell_cell_extended_pruned_idx[1], "_cell_cell_extended_pruned:: ");
      }
      free(order);
    }
    shift_part += part_ext->n_part[i_domain];
  }

  // printf("_prune_cell_cell_extented end \n");
}


static
void
_generate_extended_partition_connectivity
(
 int           n_entity1,
 int           n_entity2,
 int          *entity1_entity1_extended_idx,
 int          *entity1_entity1_extended,
 PDM_g_num_t  *entity2_ln_to_gn,
 int          *entity1_entity2_idx,
 int          *entity1_entity2,
 int          *border_gentity1_entity2_n,
 PDM_g_num_t  *border_gentity1_entity2,
 int          *border_part_and_proc_id,
 int          *border_lentity1_entity2,
 int         **entity2_entity2_extended_idx,
 int         **entity2_entity2_extended,
 PDM_g_num_t **border_entity2_ln_to_gn
)
{
  PDM_UNUSED(entity1_entity1_extended);

  PDM_g_num_t* gentity1_entity2 = (PDM_g_num_t *) malloc( entity1_entity2_idx[n_entity1] * sizeof(PDM_g_num_t));
  for(int i = 0; i < entity1_entity2_idx[n_entity1]; ++i) {
    gentity1_entity2[i] = entity2_ln_to_gn[PDM_ABS(entity1_entity2[i])-1];
  }

  int n_neight_tot = entity1_entity1_extended_idx[n_entity1];
  int *_border_gentity1_entity2_idx = (int * ) malloc( (entity1_entity1_extended_idx[n_entity1]+1) * sizeof(int) );

  _border_gentity1_entity2_idx[0] = 0;
  int s_tot = 0;
  for(int i = 0; i < n_neight_tot; ++i) {
    s_tot += border_gentity1_entity2_n[i];
    _border_gentity1_entity2_idx[i+1] = _border_gentity1_entity2_idx[i] + border_gentity1_entity2_n[i];
  }

  if(0 == 1) {
    PDM_log_trace_array_int (_border_gentity1_entity2_idx, n_neight_tot+1, "_border_gcell_face_idx::");
    PDM_log_trace_array_int (border_gentity1_entity2_n  , n_neight_tot  , "border_gentity1_entity2_n::");
    PDM_log_trace_array_long(border_gentity1_entity2, s_tot, "border_gentity1_entity2::");
  }

  PDM_g_num_t* _border_entity2_ln_to_gn = (PDM_g_num_t * ) malloc( s_tot * sizeof(PDM_g_num_t));

  /*
   * Prepare and order the current entity ln_to_gn
   *   Cause each entity1 can have connectivity of interior or a new entitity2 (from neightborood)
   */
  PDM_g_num_t* _sorted_entity2_ln_to_gn = (PDM_g_num_t * ) malloc( n_entity2 * sizeof(PDM_g_num_t));
  for(int i_entity2 = 0; i_entity2 < n_entity2; ++i_entity2 ) {
    _sorted_entity2_ln_to_gn[i_entity2] = entity2_ln_to_gn[i_entity2];
  }

  int* order                 = (int *) malloc( n_entity2                      * sizeof(int));
  int* order_entity1_entity2 = (int *) malloc( entity1_entity2_idx[n_entity1] * sizeof(int));
  for(int i = 0; i < n_entity2; ++i) {
    order[i] = i;
  }
  for(int i = 0; i < entity1_entity2_idx[n_entity1]; ++i) {
    order_entity1_entity2[i] = i;
  }

  PDM_sort_long(_sorted_entity2_ln_to_gn, order                , n_entity2                     );
  PDM_sort_long(gentity1_entity2        , order_entity1_entity2, entity1_entity2_idx[n_entity1]);
  // abort(); // Il faut trier le cell_face !!!!! --> Permet de prendre le bon signe aprés !

  if(0 == 1) {
    PDM_log_trace_array_long(gentity1_entity2, entity1_entity2_idx[n_entity1], "gentity1_entity2::");
    PDM_log_trace_array_int (order_entity1_entity2, entity1_entity2_idx[n_entity1], "order_entity1_entity2::");
  }

  /*
   * Do the same but for the boundary limit
   */
  int* border_order = (int * ) malloc( s_tot * sizeof(int));
  for(int i = 0; i < s_tot; ++i) {
    _border_entity2_ln_to_gn[i] = PDM_ABS(border_gentity1_entity2[i]);
  }

  // Ce qui compte ici c'est le order de la cellule !!!!
  for(int i = 0; i < n_neight_tot; ++i) {
    for(int j = _border_gentity1_entity2_idx[i]; j < _border_gentity1_entity2_idx[i+1]; ++j) {
      border_order[j] = j;
    }
  }

  if(0 == 1) {
    PDM_log_trace_array_int (border_order, s_tot, "border_order::");
    PDM_log_trace_array_long(_border_entity2_ln_to_gn, s_tot, "_border_entity2_ln_to_gn (avant tri) ::");
  }

  int n_entity2_unique = PDM_inplace_unique_long(_border_entity2_ln_to_gn, border_order, 0, s_tot-1);
  _border_entity2_ln_to_gn = realloc(_border_entity2_ln_to_gn, n_entity2_unique * sizeof(PDM_g_num_t));
  *border_entity2_ln_to_gn = _border_entity2_ln_to_gn;

  int *old_to_new_order_border = (int *) malloc (s_tot * sizeof(int));
  for(int i = 0; i < s_tot; i++) {
   old_to_new_order_border[border_order[i]] = i;
  }

  int* border_entity1_order  = (int *) malloc( s_tot * sizeof(int));
  int* border_entity2_unique = (int *) malloc( s_tot * sizeof(int));

  // Indices of the first unique in original array in order to find out the new
  int* border_entity2_first_unique = (int *) malloc( n_entity2_unique * sizeof(int));

  PDM_g_num_t last  = 0; // Gn is never 0 because start as 1
  int idx_unique = -1;
  int idx_first  = 0;
  for(int i = 0; i < s_tot; ++i) {
    int old_order = border_order[i];
    int pos = PDM_binary_search_gap_int(old_order, _border_gentity1_entity2_idx, n_neight_tot+1);
    border_entity1_order[i] = pos;

    if(last != PDM_ABS(border_gentity1_entity2[old_order])) {
      idx_unique = i;
      last = PDM_ABS(border_gentity1_entity2[old_order]);
      border_entity2_first_unique[idx_first++] = i;
    }
    border_entity2_unique[old_order] = idx_unique;

    // printf(" Search idx -> border_order[%i] = %i  --> idx_unique = %i | last = %i \n", i, old_order, idx_unique, (int)last);
    // border_entity1_order[i] = -1;
    // printf(" Associated cell = %i \n", pos);
  }

  if(0 == 1) {
    // printf(" --------------------------------------------  \n");
    // for(int i = 0; i < s_tot; ++i) {
    //   printf("border_entity2_unique[%i] = %i -> %i \n", i, border_entity2_unique[i], border_entity2_unique[border_order[i]]);
    // }
    // printf(" -------------------------------------------- \n");
    PDM_log_trace_array_long(_border_entity2_ln_to_gn, n_entity2_unique, "_border_entity2_ln_to_gn::");
    PDM_log_trace_array_int(border_order, s_tot, "border_order::");
    PDM_log_trace_array_int(border_entity1_order, s_tot, "border_entity1_order::");
    PDM_log_trace_array_int(border_entity2_unique, s_tot, "border_entity2_unique::");
    PDM_log_trace_array_int(border_entity2_first_unique, n_entity2_unique, "border_entity2_first_unique::");
  }
  // exit(1);

  /* Pour chaque elements on chercher si il est dans les entity2_ln_to_gn
   *   Si ce n'est pas le cas, c'est un nouvelle element, on parcours la liste unique des bords
   *   donc les eléments des bordes sont également unique
   *   En même temps on réalise la construction du nouveau graph d'échanges pour entity2
   */
  *entity2_entity2_extended_idx = (int * ) malloc( (    n_entity2 + 1   ) * sizeof(int));
  *entity2_entity2_extended     = (int * ) malloc( (3 * n_entity2_unique) * sizeof(int));
  int* _entity2_entity2_extended_idx = *entity2_entity2_extended_idx;
  int* _entity2_entity2_extended     = *entity2_entity2_extended;

  PDM_g_num_t *entity2_extended_gnum = (PDM_g_num_t * ) malloc( n_entity2_unique * sizeof(PDM_g_num_t));
  int n_entity2_extended = 0;
  _entity2_entity2_extended_idx[0] = 0;
  _entity2_entity2_extended_idx[1] = 0;
  int idx_write = 0;

  for(int i_entity2 = 0; i_entity2 < n_entity2_unique; ++i_entity2) {
    PDM_g_num_t g_entity2 = _border_entity2_ln_to_gn[i_entity2];
    int pos = PDM_binary_search_long(g_entity2, _sorted_entity2_ln_to_gn, n_entity2);
    if(pos == -1) {
      entity2_extended_gnum[n_entity2_extended++] = g_entity2;
      // printf(" [%i] found [%i] = %i\n", i_part+shift_part, i_entity2, pos);

      _entity2_entity2_extended_idx[1]++;

      // int old_entity1_order = border_entity1_order[i_entity2];
      // int old_entity2_order = border_order        [i_entity2];
      // int opp_proc    = entity1_entity1_extended[3*old_entity1_order  ];
      // int opp_part    = entity1_entity1_extended[3*old_entity1_order+1];
      // int opp_entity1 = entity1_entity1_extended[3*old_entity1_order+2];

      int pos_first_unique = border_entity2_first_unique[i_entity2];
      int old_order        = border_order[pos_first_unique];
      int gopp_entity2     = border_gentity1_entity2[old_order];
      int opp_entity2      = border_lentity1_entity2[old_order];
      int opp_proc    = border_part_and_proc_id[2*old_order  ];
      int opp_part    = border_part_and_proc_id[2*old_order+1];

      // printf(" border_gentity1_entity2[%i] = %i for g_entity2 = %i \n", old_order, border_gentity1_entity2[old_order], g_entity2);

      // printf(" old_entity1_order = %i | old_entity2_order = %i | new_entity2 = %i \n", old_entity1_order, old_entity2_order, new_entity2);
      // printf(" [pos = %i] [opp_proc = %i | opp_part = %i | opp_entity1 = %i | opp_entity2 = %i] | g_entity2 = %i | gopp_entity2 = %i\n", pos, opp_proc, opp_part, opp_entity1, PDM_ABS(opp_entity2), (int)g_entity2, gopp_entity2);
      // printf(" pos_first_unique = %i | old_order = %i | g_num = %i | g_num_check = %i \n",
      //        pos_first_unique, old_order, g_entity2, gopp_entity2);
      assert(g_entity2 == PDM_ABS(gopp_entity2));

      // printf("_entity2_entity2_extended[%i] = [%i/%i/%i] \n", idx_write, opp_proc, opp_part, PDM_ABS(opp_entity2)-1);

      _entity2_entity2_extended[3*idx_write  ] = opp_proc;
      _entity2_entity2_extended[3*idx_write+1] = opp_part;
      _entity2_entity2_extended[3*idx_write+2] = PDM_ABS(opp_entity2)-1;
      idx_write++;
    }

  }
  // printf(" ----- \n");

  for(int i = 1; i < n_entity2; ++i) {
    _entity2_entity2_extended_idx[i+1] = _entity2_entity2_extended_idx[i];
  }
  *entity2_entity2_extended = realloc(*entity2_entity2_extended,  (3 * n_entity2_extended) * sizeof(int));

  if(0 == 1) {
    _entity2_entity2_extended = *entity2_entity2_extended;
    PDM_log_trace_array_long(entity2_extended_gnum       , n_entity2_extended    , "entity2_extended_gnum::");
    PDM_log_trace_array_int(_entity2_entity2_extended_idx, n_entity2+1           , "_entity2_entity2_extended_idx::");
    PDM_log_trace_array_int(_entity2_entity2_extended    , 3 * n_entity2_extended, "_entity2_entity2_extended::");
  }

  free(old_to_new_order_border);
  // exit(1);

  /*
   * Reconstruction de la connectivité de bord
   *   On écrase border_lentity1_entity2
   */
  int idx = 0;
  int i_entity2_extented = 0;
  for(int i = 0; i < s_tot; ++i) {
    PDM_g_num_t g_entity2 = PDM_ABS(border_gentity1_entity2[i]);

    /* On cherche d'abord dans le bord - face_extended_gnum is sort by construction */
    int pos = PDM_binary_search_long(g_entity2, entity2_extended_gnum, n_entity2_extended);

    if(pos != -1) {
      int sgn    = PDM_SIGN(border_lentity1_entity2[i]); // A aller cherche dans le cell_face de depart

      // On doit chercher le sign de l'entity2 qu'on garde, car elle impose le signe
      // int i_unique  = border_entity2_unique[i];
      // int old_order = border_order[i_unique];
      // int g_sgn     = PDM_SIGN(border_gentity1_entity2[old_order]);
      // PDM_g_num_t g_num_check = border_gentity1_entity2[old_order];

      // printf(" Border face comming for other proc : pos = %i | idx = %i | old_order = %i | g_sgn = %i | g_num_check = %i | g_entity2 = %i\n", pos, idx, old_order, g_sgn, (int)g_num_check, (int)g_entity2);
      // printf(" Rebuild from exterior [%i] with gnum = "PDM_FMT_G_NUM" and pos : %i - new numbering %i \n ", i, g_entity2, pos, ( pos + n_entity2 + 1 ));

      border_lentity1_entity2[idx++] = sgn * ( pos + n_entity2 + 1 ); // Car on shift
      i_entity2_extented++;
    } else {

      int pos_interior2 = PDM_binary_search_long(g_entity2, _sorted_entity2_ln_to_gn, n_entity2);
      int pos_interior = PDM_binary_search_long(g_entity2, gentity1_entity2, entity1_entity2_idx[n_entity1]);
      // printf(" Border face comming from interior %i - %i \n", pos_interior, idx);
      assert(pos_interior  != -1);
      assert(pos_interior2 != -1);

      // int old_order_entity1_entity2 = order_entity1_entity2[pos_interior];

      // Keep it to check
      // int sgn       = PDM_SIGN(gentity1_entity2[old_order_entity1_entity2]);
      // int i_entity2 = PDM_ABS ( entity1_entity2[old_order_entity1_entity2]);
      // int i_unique  = border_entity2_unique[i];
      // int old_order = border_order[i_unique];
      // int g_sgn     = PDM_SIGN(border_gentity1_entity2[old_order]);
      // PDM_g_num_t g_num_check = border_gentity1_entity2[old_order];
      // printf(" Border face comming from interior %i - %i | i_unique = %i | old_order = %i | g_sgn = %i | sgn = %i | g_num_check = %i | g_entity2 = %i\n", pos, idx, i_unique, old_order, g_sgn, sgn, (int)g_num_check, (int)g_entity2);
      // printf("OKOK border_gentity1_entity2[%i] = %i for g_entity2 = %i \n", old_order, border_gentity1_entity2[old_order], g_entity2);

      // border_lentity1_entity2[idx++] = sgn * ( order[pos_interior] + 1 ); // Car le tableau est trié pas comme la partition
      // ON MET FORCEMENT L'OPPOSE
      // printf(" Rebuild from interior [%i] with gnum = "PDM_FMT_G_NUM" and pos : %i \n ", i, g_entity2, pos_interior);
      // OLD -> FOnctionne pas pour le edge_vtx
      // border_lentity1_entity2[idx++] = - entity1_entity2[old_order_entity1_entity2]; // Car le tableau est trié pas comme la partition

      int sgn    = PDM_SIGN(border_lentity1_entity2[i]); // A aller cherche dans le cell_face de depart
      int old_pos = order[pos_interior2];
      border_lentity1_entity2[idx++] = sgn * ( old_pos + 1 );

      // PDM_g_num_t old_g_num = _sorted_entity2_ln_to_gn[pos_interior2];
      // printf("Cas 1 : border_lentity1_entity2 [%i] | Cas 2  : %i \n", - entity1_entity2[old_order_entity1_entity2], sgn * ( old_pos + 1 ));
      // printf("Cas 1 : border_lentity1_entity2 [%i] \n", sgn * old_g_num);

    }
  }
  // exit(1);

  // Mise à jour
  _border_entity2_ln_to_gn = realloc(_border_entity2_ln_to_gn, n_entity2_extended * sizeof(PDM_g_num_t));
  *border_entity2_ln_to_gn = _border_entity2_ln_to_gn;
  for(int i = 0; i < n_entity2_extended; ++i){
    _border_entity2_ln_to_gn[i] = entity2_extended_gnum[i];
  }

  /*
   * Free
   */
  free(gentity1_entity2);
  free(order_entity1_entity2);
  free(entity2_extended_gnum);
  free(border_entity1_order);
  free(border_entity2_unique);
  free(border_entity2_first_unique);
  free(border_order);
  free(order);
  free(_sorted_entity2_ln_to_gn);
  // free(_border_entity2_ln_to_gn); //
  free(_border_gentity1_entity2_idx);
}




// static
// void
// _rebuild_connectivity_cell_face_debug
// (
//   PDM_part_extension_t *part_ext
// )
// {
//   // printf("_rebuild_connectivity_cell_face \n");

//   int n_tot_all_domain = 0;
//   int n_part_loc_all_domain = 0;
//   for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
//     n_tot_all_domain      += part_ext->n_tot_part_by_domain[i_domain];
//     n_part_loc_all_domain += part_ext->n_part[i_domain];
//   }

//   PDM_distant_neighbor_t* dn = PDM_distant_neighbor_create(part_ext->comm,
//                                                            n_part_loc_all_domain,
//                                                            part_ext->n_cell,
//                                                            part_ext->cell_cell_extended_pruned_idx,
//                                                            part_ext->cell_cell_extended_pruned    );

//   /* On doit échanger toutes les connectivités en un seul coup et de manière descendante */
//   /* Donc par exemple cell_face + face_vtx
//    * Ou cell_face + face_edge + edge_vtx
//    * La deduction des autres se fait par transitivité local
//    * Il faut également gerer les conditions limites
//    */

//   /* Prepare */
//   int         **cell_face_n   = (int         **) malloc( n_part_loc_all_domain * sizeof(int         *));
//   PDM_g_num_t **gcell_face    = (PDM_g_num_t **) malloc( n_part_loc_all_domain * sizeof(PDM_g_num_t *));
//   int         **lcell_face    = (int         **) malloc( n_part_loc_all_domain * sizeof(int         *));
//   PDM_g_num_t **cell_ln_to_gn = (PDM_g_num_t **) malloc( n_part_loc_all_domain * sizeof(PDM_g_num_t *));
//   // PDM_g_num_t **cell_flags    = (PDM_g_num_t **) malloc( n_part_loc_all_domain * sizeof(PDM_g_num_t *));

//   int shift_part = 0;
//   for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
//     for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {

//       int* cell_face_idx =  part_ext->parts[i_domain][i_part].cell_face_idx;
//       int* cell_face     =  part_ext->parts[i_domain][i_part].cell_face;

//       lcell_face[i_part+shift_part] = cell_face;

//       int n_cell      = part_ext->parts[i_domain][i_part].n_cell;
//       // int n_face      = part_ext->parts[i_domain][i_part].n_face;
//       int s_cell_face = cell_face_idx[n_cell];

//       cell_face_n[i_part+shift_part] = (int         *) malloc( n_cell      * sizeof(int        ));
//       gcell_face [i_part+shift_part] = (PDM_g_num_t *) malloc( s_cell_face * sizeof(PDM_g_num_t));

//       PDM_g_num_t* face_ln_to_gn = part_ext->parts[i_domain][i_part].face_ln_to_gn;

//       cell_ln_to_gn[i_part+shift_part] = part_ext->parts[i_domain][i_part].cell_ln_to_gn;

//       for(int i_cell = 0; i_cell < n_cell; ++i_cell) {
//         cell_face_n[i_part+shift_part][i_cell] = cell_face_idx[i_cell+1] - cell_face_idx[i_cell];
//         for(int idx_face = cell_face_idx[i_cell]; idx_face < cell_face_idx[i_cell+1]; ++idx_face) {
//           int sgn    = PDM_SIGN(cell_face[idx_face]);
//           int i_face = PDM_ABS (cell_face[idx_face])-1;
//           gcell_face[i_part+shift_part][idx_face] = sgn * face_ln_to_gn[i_face];
//           // printf("gcell_face[%i][%i] = %i \n", i_part+shift_part, idx_face, i_part);
//         }
//       }

//     }
//     shift_part += part_ext->n_part[i_domain];
//   }

//   /* Exchange */
//   int         **border_gcell_face_n;
//   PDM_g_num_t **border_gcell_face;
//   PDM_distant_neighbor_exch(dn,
//                             sizeof(PDM_g_num_t),
//                             PDM_STRIDE_VAR,
//                             -1,
//                             cell_face_n,
//                  (void **)  gcell_face,
//                            &border_gcell_face_n,
//                 (void ***) &border_gcell_face);

//   int         **border_lcell_face_n;
//   int         **border_lcell_face;
//   PDM_distant_neighbor_exch(dn,
//                             sizeof(int),
//                             PDM_STRIDE_VAR,
//                             -1,
//                             cell_face_n,
//                  (void **)  lcell_face,
//                            &border_lcell_face_n,
//                 (void ***) &border_lcell_face);

//   /* On fait le cell_ln_to_gn par la même occasion */
//   PDM_distant_neighbor_exch(dn,
//                             sizeof(PDM_g_num_t),
//                             PDM_STRIDE_CST,
//                             1,
//                             NULL,
//                  (void **)  cell_ln_to_gn,
//                             NULL,
//                 (void ***) &part_ext->border_cell_ln_to_gn);

//   free(lcell_face);
//   free(cell_ln_to_gn);

//   /* Post treatment */
//   shift_part = 0;
//   for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
//     for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {

//       int n_cell        = part_ext->n_cell[i_part+shift_part];
//       int n_face        = part_ext->parts[i_domain][i_part].n_face;

//       int* face_face_extended_pruned_idx;
//       int* face_face_extended_pruned;
//       PDM_g_num_t* face_ln_to_gn;
//       _generate_extended_partition_connectivity(n_cell,
//                                                 n_face,
//                                                 part_ext->cell_cell_extended_pruned_idx[i_part+shift_part],
//                                                 part_ext->cell_cell_extended_pruned    [i_part+shift_part],
//                                                 part_ext->parts[i_domain][i_part].face_ln_to_gn,
//                                                 part_ext->parts[i_domain][i_part].cell_face_idx,
//                                                 part_ext->parts[i_domain][i_part].cell_face,
//                                                 border_gcell_face_n    [i_part+shift_part],
//                                                 border_gcell_face      [i_part+shift_part],
//                                                 border_lcell_face      [i_part+shift_part],
//                                                 border_part_and_proc_id[i_part+shift_part],
//                                                 NULL,
//                                                &face_face_extended_pruned_idx,
//                                                &face_face_extended_pruned,
//                                                &face_ln_to_gn);
//       free(face_face_extended_pruned_idx);
//       free(face_face_extended_pruned);
//       free(face_ln_to_gn);
//     }
//     shift_part += part_ext->n_part[i_domain];
//   }

//   /* Pour les faces group on peut faire aussi le gnum location --> Marche pas en multidomain (ou il faut shifter )*/
//   shift_part = 0;
//   for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
//     for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {
//       free(cell_face_n[i_part+shift_part]);
//       free(gcell_face[i_part+shift_part]);
//       free(border_gcell_face_n[i_part+shift_part]);
//       free(border_lcell_face_n[i_part+shift_part]);
//       free(border_lcell_face[i_part+shift_part]);
//       free(border_gcell_face[i_part+shift_part]);
//     }
//     shift_part += part_ext->n_part[i_domain];
//   }

//   PDM_distant_neighbor_free(dn);
//   free(cell_face_n);
//   free(gcell_face);
//   free(border_gcell_face_n);
//   free(border_gcell_face);
//   free(border_lcell_face_n);
//   free(border_lcell_face);
//   // printf("_rebuild_connectivity end \n");
// }


static
void
_rebuild_connectivity
(
  PDM_part_extension_t  *part_ext,
  int                   *n_entity1,
  int                   *n_entity2,
  int                  **entity1_entity1_extended_idx,
  int                  **entity1_entity1_extended,
  int                  **entity1_entity2_idx,
  int                  **entity1_entity2,
  PDM_g_num_t          **entity2_ln_to_gn,
  int                 ***border_lentity1_entity2_idx,
  int                 ***border_lentity1_entity2,
  int                 ***entity2_entity2_extended_idx,
  int                 ***entity2_entity2_extended,
  PDM_g_num_t         ***border_entity2_ln_to_gn
)
{
  // printf("_rebuild_connectivity \n");

  int i_rank = -1;
  PDM_MPI_Comm_rank(part_ext->comm, &i_rank);

  int n_tot_all_domain = 0;
  int n_part_loc_all_domain = 0;
  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
    n_tot_all_domain      += part_ext->n_tot_part_by_domain[i_domain];
    n_part_loc_all_domain += part_ext->n_part[i_domain];
  }

  *entity2_entity2_extended_idx = (int         ** ) malloc( n_part_loc_all_domain * sizeof(int         *));
  *entity2_entity2_extended     = (int         ** ) malloc( n_part_loc_all_domain * sizeof(int         *));
  *border_entity2_ln_to_gn      = (PDM_g_num_t ** ) malloc( n_part_loc_all_domain * sizeof(PDM_g_num_t *));
  *border_lentity1_entity2_idx  = (int         ** ) malloc( n_part_loc_all_domain * sizeof(int         *));

  PDM_distant_neighbor_t* dn = PDM_distant_neighbor_create(part_ext->comm,
                                                           n_part_loc_all_domain,
                                                           n_entity1,
                                                           entity1_entity1_extended_idx,
                                                           entity1_entity1_extended);

  /* Prepare */
  int         **entity1_entity2_n = (int         **) malloc( n_part_loc_all_domain * sizeof(int         *));
  PDM_g_num_t **gentity1_entity2  = (PDM_g_num_t **) malloc( n_part_loc_all_domain * sizeof(PDM_g_num_t *));

  int         **part_and_proc_id  = (int         **) malloc( n_part_loc_all_domain * sizeof(int         *));

  int shift_part = 0;
  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
    for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {

      int* _pentity1_entity2_idx = entity1_entity2_idx[i_part+shift_part];
      int* _pentity1_entity2     = entity1_entity2    [i_part+shift_part];

      int pn_entity1        = n_entity1[i_part+shift_part];
      int s_entity1_entity2 = _pentity1_entity2_idx[pn_entity1];

      entity1_entity2_n[i_part+shift_part] = (int         *) malloc( pn_entity1        * sizeof(int        ));
      gentity1_entity2 [i_part+shift_part] = (PDM_g_num_t *) malloc( s_entity1_entity2 * sizeof(PDM_g_num_t));

      part_and_proc_id [i_part+shift_part] = (int         *) malloc( 2 * s_entity1_entity2 * sizeof(int));

      PDM_g_num_t* _pentity2_ln_to_gn = entity2_ln_to_gn[i_part+shift_part];

      for(int i1 = 0; i1 < pn_entity1; ++i1) {
        entity1_entity2_n[i_part+shift_part][i1] = _pentity1_entity2_idx[i1+1] - _pentity1_entity2_idx[i1];
        for(int idx_entity2 = _pentity1_entity2_idx[i1]; idx_entity2 < _pentity1_entity2_idx[i1+1]; ++idx_entity2) {
          int sgn    = PDM_SIGN(_pentity1_entity2[idx_entity2]);
          int i_face = PDM_ABS (_pentity1_entity2[idx_entity2])-1;
          gentity1_entity2[i_part+shift_part][idx_entity2] = sgn * _pentity2_ln_to_gn[i_face];
          part_and_proc_id[i_part+shift_part][2*idx_entity2  ] = i_rank;
          part_and_proc_id[i_part+shift_part][2*idx_entity2+1] = i_part+shift_part;
        }
      }

    }
    shift_part += part_ext->n_part[i_domain];
  }

  /* Exchange */
  int         **border_gentity1_entity2_n;
  PDM_g_num_t **border_gentity1_entity2;
  PDM_distant_neighbor_exch(dn,
                            sizeof(PDM_g_num_t),
                            PDM_STRIDE_VAR_INTERLACED,
                            -1,
                            entity1_entity2_n,
                 (void **)  gentity1_entity2,
                           &border_gentity1_entity2_n,
                (void ***) &border_gentity1_entity2);

  int         **border_lentity1_entity2_n;
  PDM_distant_neighbor_exch(dn,
                            sizeof(int),
                            PDM_STRIDE_VAR_INTERLACED,
                            -1,
                            entity1_entity2_n,
                 (void **)  entity1_entity2,
                           &border_lentity1_entity2_n,
                (void ***)  border_lentity1_entity2);

  int         **border_part_and_proc_id;
  int         **border_lentity1_entity2_n_tmp;
  PDM_distant_neighbor_exch(dn,
                            2 * sizeof(int),
                            PDM_STRIDE_VAR_INTERLACED,
                            -1,
                            entity1_entity2_n,
                 (void **)  part_and_proc_id,
                           &border_lentity1_entity2_n_tmp,
                (void ***) &border_part_and_proc_id);

  shift_part = 0;
  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
    for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {
      free(border_lentity1_entity2_n_tmp[i_part+shift_part]);
      free(part_and_proc_id[i_part+shift_part]);
    }
    shift_part += part_ext->n_part[i_domain];
  }
  free(border_lentity1_entity2_n_tmp);
  free(part_and_proc_id);

  /* Post treatment */
  shift_part = 0;
  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
    for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {

      int pn_entity1 = n_entity1[i_part+shift_part];
      int pn_entity2 = n_entity2[i_part+shift_part];

      int         *pborder_lentity1_entity2      = (*border_lentity1_entity2)[i_part+shift_part];
      int         *pentity2_entity2_extended_idx = NULL;
      int         *pentity2_entity2_extended     = NULL;
      PDM_g_num_t *pentity2_ln_to_gn             = NULL;
      _generate_extended_partition_connectivity(pn_entity1,
                                                pn_entity2,
                                                entity1_entity1_extended_idx[i_part+shift_part],
                                                entity1_entity1_extended    [i_part+shift_part],
                                                entity2_ln_to_gn            [i_part+shift_part],
                                                entity1_entity2_idx         [i_part+shift_part],
                                                entity1_entity2             [i_part+shift_part],
                                                border_gentity1_entity2_n   [i_part+shift_part],
                                                border_gentity1_entity2     [i_part+shift_part],
                                                border_part_and_proc_id     [i_part+shift_part],
                                                pborder_lentity1_entity2,
                                               &pentity2_entity2_extended_idx,
                                               &pentity2_entity2_extended,
                                               &pentity2_ln_to_gn);

      (*entity2_entity2_extended_idx)[i_part+shift_part] = pentity2_entity2_extended_idx;
      (*entity2_entity2_extended    )[i_part+shift_part] = pentity2_entity2_extended;
      (*border_entity2_ln_to_gn     )[i_part+shift_part] = pentity2_ln_to_gn;

      /* On refait l'index */
      int n_neight_tot = entity1_entity1_extended_idx[i_part+shift_part][pn_entity1];
      (*border_lentity1_entity2_idx)[i_part+shift_part] = PDM_array_new_idx_from_sizes_int(border_lentity1_entity2_n[i_part+shift_part], n_neight_tot);
    }
    shift_part += part_ext->n_part[i_domain];
  }

  /* Pour les faces group on peut faire aussi le gnum location --> Marche pas en multidomain (ou il faut shifter )*/
  shift_part = 0;
  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
    for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {
      free(entity1_entity2_n[i_part+shift_part]);
      free(gentity1_entity2[i_part+shift_part]);
      free(border_gentity1_entity2_n[i_part+shift_part]);
      free(border_lentity1_entity2_n[i_part+shift_part]);
      // free(border_lentity1_entity2[i_part+shift_part]);
      free(border_gentity1_entity2[i_part+shift_part]);
      free(border_part_and_proc_id[i_part+shift_part]);
    }
    shift_part += part_ext->n_part[i_domain];
  }

  PDM_distant_neighbor_free(dn);
  free(entity1_entity2_n);
  free(gentity1_entity2);
  free(border_gentity1_entity2_n);
  free(border_gentity1_entity2);
  free(border_lentity1_entity2_n);
  free(border_part_and_proc_id);
  // printf("_rebuild_connectivity end \n");

}



static
void
_rebuild_connectivity_cell_face
(
  PDM_part_extension_t *part_ext
)
{

  int n_tot_all_domain = 0;
  int n_part_loc_all_domain = 0;
  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
    n_tot_all_domain      += part_ext->n_tot_part_by_domain[i_domain];
    n_part_loc_all_domain += part_ext->n_part[i_domain];
  }

  /* Cell face */
  int          *n_cell        = (int         * ) malloc( n_part_loc_all_domain * sizeof(int          ));
  int          *n_face        = (int         * ) malloc( n_part_loc_all_domain * sizeof(int          ));
  int         **cell_face_idx = (int         **) malloc( n_part_loc_all_domain * sizeof(int         *));
  int         **cell_face     = (int         **) malloc( n_part_loc_all_domain * sizeof(int         *));
  PDM_g_num_t **face_ln_to_gn = (PDM_g_num_t **) malloc( n_part_loc_all_domain * sizeof(PDM_g_num_t *));
  PDM_g_num_t **cell_ln_to_gn = (PDM_g_num_t **) malloc( n_part_loc_all_domain * sizeof(PDM_g_num_t *));

  int shift_part = 0;
  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
    for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {
      n_cell       [i_part+shift_part] = part_ext->parts[i_domain][i_part].n_cell;
      n_face       [i_part+shift_part] = part_ext->parts[i_domain][i_part].n_face;
      cell_face_idx[i_part+shift_part] = part_ext->parts[i_domain][i_part].cell_face_idx;
      cell_face    [i_part+shift_part] = part_ext->parts[i_domain][i_part].cell_face;
      face_ln_to_gn[i_part+shift_part] = part_ext->parts[i_domain][i_part].face_ln_to_gn;
      cell_ln_to_gn[i_part+shift_part] = part_ext->parts[i_domain][i_part].cell_ln_to_gn;
    }
    shift_part += part_ext->n_part[i_domain];
  }


  PDM_distant_neighbor_t* dn = PDM_distant_neighbor_create(part_ext->comm,
                                                           n_part_loc_all_domain,
                                                           part_ext->n_cell,
                                                           part_ext->cell_cell_extended_pruned_idx,
                                                           part_ext->cell_cell_extended_pruned    );
  assert(part_ext->border_cell_ln_to_gn == NULL);
  /* On fait le cell_ln_to_gn par la même occasion */
  PDM_distant_neighbor_exch(dn,
                            sizeof(PDM_g_num_t),
                            PDM_STRIDE_CST_INTERLACED,
                            1,
                            NULL,
                 (void **)  cell_ln_to_gn,
                            NULL,
                (void ***) &part_ext->border_cell_ln_to_gn);

  PDM_distant_neighbor_free(dn);

  // int **border_lcell_face_idx;
  // int **border_lcell_face;
  // int **face_face_extended_idx;
  // int **face_face_extended;
  assert(part_ext->face_face_extended_idx == NULL);
  assert(part_ext->face_face_extended     == NULL);
  _rebuild_connectivity(part_ext,
                        n_cell,
                        n_face,
                        part_ext->cell_cell_extended_pruned_idx,
                        part_ext->cell_cell_extended_pruned,
                        cell_face_idx,
                        cell_face,
                        face_ln_to_gn,
                       &part_ext->border_cell_face_idx,
                       &part_ext->border_cell_face,
                       &part_ext->face_face_extended_idx,
                       &part_ext->face_face_extended,
                       &part_ext->border_face_ln_to_gn);

  free(n_cell       );
  free(n_face       );
  free(cell_face_idx);
  free(cell_face    );
  free(face_ln_to_gn);
  free(cell_ln_to_gn);
}



static
void
_rebuild_connectivity_face_vtx
(
  PDM_part_extension_t *part_ext
)
{

  int n_tot_all_domain = 0;
  int n_part_loc_all_domain = 0;
  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
    n_tot_all_domain      += part_ext->n_tot_part_by_domain[i_domain];
    n_part_loc_all_domain += part_ext->n_part[i_domain];
  }

  /* Cell face */
  int          *n_face        = (int         * ) malloc( n_part_loc_all_domain * sizeof(int          ));
  int          *n_vtx         = (int         * ) malloc( n_part_loc_all_domain * sizeof(int          ));
  int         **face_vtx_idx  = (int         **) malloc( n_part_loc_all_domain * sizeof(int         *));
  int         **face_vtx      = (int         **) malloc( n_part_loc_all_domain * sizeof(int         *));
  PDM_g_num_t **vtx_ln_to_gn  = (PDM_g_num_t **) malloc( n_part_loc_all_domain * sizeof(PDM_g_num_t *));

  int shift_part = 0;
  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
    for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {
      n_face       [i_part+shift_part] = part_ext->parts[i_domain][i_part].n_face;
      n_vtx        [i_part+shift_part] = part_ext->parts[i_domain][i_part].n_vtx;
      face_vtx_idx [i_part+shift_part] = part_ext->parts[i_domain][i_part].face_vtx_idx;
      face_vtx     [i_part+shift_part] = part_ext->parts[i_domain][i_part].face_vtx;
      vtx_ln_to_gn [i_part+shift_part] = part_ext->parts[i_domain][i_part].vtx_ln_to_gn;
    }
    shift_part += part_ext->n_part[i_domain];
  }

  assert(part_ext->vtx_vtx_extended_idx == NULL);
  assert(part_ext->vtx_vtx_extended     == NULL);
  _rebuild_connectivity(part_ext,
                        n_face,
                        n_vtx,
                        part_ext->face_face_extended_idx,
                        part_ext->face_face_extended,
                        face_vtx_idx,
                        face_vtx,
                        vtx_ln_to_gn,
                       &part_ext->border_face_vtx_idx,
                       &part_ext->border_face_vtx,
                       &part_ext->vtx_vtx_extended_idx,
                       &part_ext->vtx_vtx_extended,
                       &part_ext->border_vtx_ln_to_gn);

  free(n_face      );
  free(n_vtx       );
  free(face_vtx_idx);
  free(face_vtx    );
  free(vtx_ln_to_gn);
}


static
void
_rebuild_connectivity_face_edge
(
  PDM_part_extension_t *part_ext
)
{

  int n_tot_all_domain = 0;
  int n_part_loc_all_domain = 0;
  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
    n_tot_all_domain      += part_ext->n_tot_part_by_domain[i_domain];
    n_part_loc_all_domain += part_ext->n_part[i_domain];
  }

  /* Cell face */
  int          *n_face         = (int         * ) malloc( n_part_loc_all_domain * sizeof(int          ));
  int          *n_edge         = (int         * ) malloc( n_part_loc_all_domain * sizeof(int          ));
  int         **face_edge_idx  = (int         **) malloc( n_part_loc_all_domain * sizeof(int         *));
  int         **face_edge      = (int         **) malloc( n_part_loc_all_domain * sizeof(int         *));
  PDM_g_num_t **edge_ln_to_gn  = (PDM_g_num_t **) malloc( n_part_loc_all_domain * sizeof(PDM_g_num_t *));

  int shift_part = 0;
  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
    for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {
      n_face        [i_part+shift_part] = part_ext->parts[i_domain][i_part].n_face;
      n_edge        [i_part+shift_part] = part_ext->parts[i_domain][i_part].n_edge;
      face_edge_idx [i_part+shift_part] = part_ext->parts[i_domain][i_part].face_edge_idx;
      face_edge     [i_part+shift_part] = part_ext->parts[i_domain][i_part].face_edge;
      edge_ln_to_gn [i_part+shift_part] = part_ext->parts[i_domain][i_part].edge_ln_to_gn;
    }
    shift_part += part_ext->n_part[i_domain];
  }

  assert(part_ext->edge_edge_extended_idx == NULL);
  assert(part_ext->edge_edge_extended     == NULL);
  _rebuild_connectivity(part_ext,
                        n_face,
                        n_edge,
                        part_ext->face_face_extended_idx,
                        part_ext->face_face_extended,
                        face_edge_idx,
                        face_edge,
                        edge_ln_to_gn,
                       &part_ext->border_face_edge_idx,
                       &part_ext->border_face_edge,
                       &part_ext->edge_edge_extended_idx,
                       &part_ext->edge_edge_extended,
                       &part_ext->border_edge_ln_to_gn);

  free(n_face       );
  free(n_edge        );
  free(face_edge_idx);
  free(face_edge    );
  free(edge_ln_to_gn );
}


static
void
_rebuild_connectivity_edge_vtx
(
  PDM_part_extension_t *part_ext
)
{

  int n_tot_all_domain = 0;
  int n_part_loc_all_domain = 0;
  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
    n_tot_all_domain      += part_ext->n_tot_part_by_domain[i_domain];
    n_part_loc_all_domain += part_ext->n_part[i_domain];
  }

  /* Cell face */
  int          *n_edge        = (int         * ) malloc( n_part_loc_all_domain * sizeof(int          ));
  int          *n_vtx         = (int         * ) malloc( n_part_loc_all_domain * sizeof(int          ));
  int         **edge_vtx_idx  = (int         **) malloc( n_part_loc_all_domain * sizeof(int         *));
  int         **edge_vtx      = (int         **) malloc( n_part_loc_all_domain * sizeof(int         *));
  PDM_g_num_t **vtx_ln_to_gn  = (PDM_g_num_t **) malloc( n_part_loc_all_domain * sizeof(PDM_g_num_t *));

  int shift_part = 0;
  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
    for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {
      n_edge      [i_part+shift_part] = part_ext->parts[i_domain][i_part].n_edge;
      n_vtx       [i_part+shift_part] = part_ext->parts[i_domain][i_part].n_vtx;
      edge_vtx_idx[i_part+shift_part] = (int *) malloc((n_edge[i_part+shift_part]+1) * sizeof(int));
      edge_vtx    [i_part+shift_part] = part_ext->parts[i_domain][i_part].edge_vtx;
      vtx_ln_to_gn[i_part+shift_part] = part_ext->parts[i_domain][i_part].vtx_ln_to_gn;

      edge_vtx_idx[i_part+shift_part][0] = 0;
      for(int i_edge = 0; i_edge < n_edge[i_part+shift_part]; ++i_edge ){
        edge_vtx_idx[i_part+shift_part][i_edge+1] = edge_vtx_idx[i_part+shift_part][i_edge] + 2;
      }

    }
    shift_part += part_ext->n_part[i_domain];
  }

  assert(part_ext->vtx_vtx_extended_idx == NULL);
  assert(part_ext->vtx_vtx_extended     == NULL);
  _rebuild_connectivity(part_ext,
                        n_edge,
                        n_vtx,
                        part_ext->edge_edge_extended_idx,
                        part_ext->edge_edge_extended,
                        edge_vtx_idx,
                        edge_vtx,
                        vtx_ln_to_gn,
                       &part_ext->border_edge_vtx_idx,
                       &part_ext->border_edge_vtx,
                       &part_ext->vtx_vtx_extended_idx,
                       &part_ext->vtx_vtx_extended,
                       &part_ext->border_vtx_ln_to_gn);

  shift_part = 0;
  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
    for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {
      free(edge_vtx_idx  [i_part+shift_part]);

      /* Generic algorithm setup a sign on edge_vtx - We remove it by swap */
      // int pn_edge = part_ext->parts[i_domain][i_part].n_edge;
      // int pn_edge_extented = part_ext->edge_edge_extended_idx[shift_part+i_part][pn_edge];
      // for(int i_edge = 0; i_edge < pn_edge_extented; ++i_edge) {
      //   int i_vtx1 = part_ext->border_edge_vtx[i_part+shift_part][2*i_edge  ];
      //   int i_vtx2 = part_ext->border_edge_vtx[i_part+shift_part][2*i_edge+1];
      //   if(i_vtx1 < 0) {
      //     part_ext->border_edge_vtx[i_part+shift_part][2*i_edge  ] = PDM_ABS(i_vtx2);
      //     part_ext->border_edge_vtx[i_part+shift_part][2*i_edge+1] = PDM_ABS(i_vtx1);
      //   } else {
      //     part_ext->border_edge_vtx[i_part+shift_part][2*i_edge  ] = PDM_ABS(i_vtx1);
      //     part_ext->border_edge_vtx[i_part+shift_part][2*i_edge+1] = PDM_ABS(i_vtx2);
      //   }
      // }

    }
    shift_part += part_ext->n_part[i_domain];
  }
  free(n_edge      );
  free(n_vtx       );
  free(edge_vtx_idx);
  free(edge_vtx    );
  free(vtx_ln_to_gn);
}

static
void
_rebuild_face_group
(
  PDM_part_extension_t *part_ext
)
{
  int n_tot_all_domain = 0;
  int n_part_loc_all_domain = 0;
  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
    n_tot_all_domain      += part_ext->n_tot_part_by_domain[i_domain];
    n_part_loc_all_domain += part_ext->n_part[i_domain];
  }

  /* Cell face */
  int          *n_face              = (int         * ) malloc( n_part_loc_all_domain * sizeof(int          ));
  int         **face_group_idg      = (int         **) malloc( n_part_loc_all_domain * sizeof(int         *));
  int         **face_group_n        = (int         **) malloc( n_part_loc_all_domain * sizeof(int         *));
  PDM_g_num_t **face_group_ln_to_gn = (PDM_g_num_t **) malloc( n_part_loc_all_domain * sizeof(PDM_g_num_t *));
  PDM_g_num_t **face_ln_to_gn_check = (PDM_g_num_t **) malloc( n_part_loc_all_domain * sizeof(PDM_g_num_t *));

  int shift_part = 0;
  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
    for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {
      n_face       [i_part+shift_part] = part_ext->parts[i_domain][i_part].n_face;

      int* _pface_group_idx = part_ext->parts[i_domain][i_part].face_bound_idx;
      int* _pface_group     = part_ext->parts[i_domain][i_part].face_bound;

      PDM_g_num_t* _pface_group_ln_to_gn = part_ext->parts[i_domain][i_part].face_group_ln_to_gn;
      PDM_g_num_t* _pface_ln_to_gn       = part_ext->parts[i_domain][i_part].face_ln_to_gn;

      /* Creation d'un champs de face contenant les id de group */
      int n_face_group = part_ext->parts[i_domain][i_part].n_face_group;
      int pn_face      = part_ext->parts[i_domain][i_part].n_face;

      int* face_group_idx = malloc( ( pn_face + 1 ) * sizeof(int));
      face_group_n[i_part+shift_part] = malloc( ( pn_face ) * sizeof(int));

      for(int i_face = 0; i_face < pn_face; ++i_face) {
        face_group_n  [i_part+shift_part][i_face] = 0;
      }

      // printf("n_face_group = %i \n", n_face_group);
      for(int i_group = 0; i_group < n_face_group; ++i_group) {
        // printf("_pface_group_idx[%i] = %i --> %i \n", i_group, _pface_group_idx[i_group], _pface_group_idx[i_group+1]);
        for(int idx_face = _pface_group_idx[i_group]; idx_face < _pface_group_idx[i_group+1]; ++idx_face) {
          int i_face = _pface_group[idx_face];
          // printf("[%i] - iface = %i \n ", i_group, i_face);
          face_group_n[i_part+shift_part][i_face-1]++;
        }
      }

      face_group_idx[0] = 0;
      for(int i_face = 0; i_face < pn_face; ++i_face) {
        // printf(" face_group_n[%i] = %i \n", i_face, face_group_n[i_part+shift_part][i_face]);
        face_group_idx[i_face+1] = face_group_idx[i_face] + face_group_n[i_part+shift_part][i_face];
        face_group_n[i_part+shift_part][i_face] = 0;
      }

      face_group_idg     [i_part+shift_part] = malloc( face_group_idx[pn_face] * sizeof(int        ));
      face_group_ln_to_gn[i_part+shift_part] = malloc( face_group_idx[pn_face] * sizeof(PDM_g_num_t));
      face_ln_to_gn_check[i_part+shift_part] = malloc( face_group_idx[pn_face] * sizeof(PDM_g_num_t));

      int         *_face_group_idg      = face_group_idg     [i_part+shift_part];
      PDM_g_num_t *_face_group_ln_to_gn = face_group_ln_to_gn[i_part+shift_part];
      PDM_g_num_t *_face_ln_to_gn_check = face_ln_to_gn_check[i_part+shift_part];

      // PDM_log_trace_array_long(_pface_ln_to_gn, pn_face, "_pface_ln_to_gn::");

      for(int i_group = 0; i_group < n_face_group; ++i_group) {
        for(int idx_face = _pface_group_idx[i_group]; idx_face < _pface_group_idx[i_group+1]; ++idx_face) {
          int i_face    = _pface_group[idx_face];
          int idx_write = face_group_idx[i_face-1] + face_group_n[i_part+shift_part][i_face-1]++;
          _face_group_idg     [idx_write] = i_group;
          _face_group_ln_to_gn[idx_write] = _pface_group_ln_to_gn[idx_face];
          // printf("[%i] _face_ln_to_gn_check[%i] = %i \n", i_part+shift_part, idx_write, (int)_pface_ln_to_gn[i_face-1]);
          _face_ln_to_gn_check[idx_write] = _pface_ln_to_gn[i_face-1];
        }
      }

      free(face_group_idx);
    }
    shift_part += part_ext->n_part[i_domain];
  }


  PDM_distant_neighbor_t* dn = PDM_distant_neighbor_create(part_ext->comm,
                                                           n_part_loc_all_domain,
                                                           n_face,
                                                           part_ext->face_face_extended_idx,
                                                           part_ext->face_face_extended);

  int** border_face_group_idg_n;
  int** border_face_group_idg;
  PDM_distant_neighbor_exch(dn,
                            sizeof(int),
                            PDM_STRIDE_VAR_INTERLACED,
                            -1,
                            face_group_n,
                  (void **) face_group_idg,
                           &border_face_group_idg_n,
                 (void ***)&border_face_group_idg);

  int         **border_face_group_ln_to_gn_n;
  PDM_g_num_t **border_face_group_ln_to_gn;
  PDM_distant_neighbor_exch(dn,
                            sizeof(PDM_g_num_t),
                            PDM_STRIDE_VAR_INTERLACED,
                            -1,
                            face_group_n,
                  (void **) face_group_ln_to_gn,
                           &border_face_group_ln_to_gn_n,
                 (void ***)&border_face_group_ln_to_gn);

  int         **border_face_ln_to_gn_check_n;
  PDM_g_num_t **border_face_ln_to_gn_check;
  PDM_distant_neighbor_exch(dn,
                            sizeof(PDM_g_num_t),
                            PDM_STRIDE_VAR_INTERLACED,
                            -1,
                            face_group_n,
                  (void **) face_ln_to_gn_check,
                           &border_face_ln_to_gn_check_n,
                 (void ***)&border_face_ln_to_gn_check);

  PDM_distant_neighbor_free(dn);

  /* Post treatment */
  // TODO MANAGEMENT of multiple domain
  assert(part_ext->n_domain == 1);

  part_ext->border_face_group_idx      = malloc( n_part_loc_all_domain * sizeof(int         *));
  part_ext->border_face_group          = malloc( n_part_loc_all_domain * sizeof(int         *));
  part_ext->border_face_group_ln_to_gn = malloc( n_part_loc_all_domain * sizeof(PDM_g_num_t *));

  shift_part = 0;
  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
    for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {
      int pn_face = part_ext->parts[i_domain][i_part].n_face;
      int n_face_border = part_ext->face_face_extended_idx[shift_part+i_part][pn_face];

      // En tout rigeur le nombre de group peut être different entre 2 domaines
      int n_face_group = part_ext->parts[i_domain][i_part].n_face_group;
      part_ext->border_face_group_idx[shift_part+i_part] = malloc( (n_face_group+1) * sizeof(int));
      int *_pborder_face_group_idx = part_ext->border_face_group_idx[shift_part+i_part];

      PDM_array_reset_int(_pborder_face_group_idx, n_face_group+1, 0);

      int idx = 0;
      for(int i = 0; i < n_face_border; ++i) {
        for(int j = 0; j < border_face_group_idg_n[shift_part+i_part][i]; ++j) {
          int i_group = border_face_group_idg[shift_part+i_part][idx++];
          _pborder_face_group_idx[i_group+1]++;
        }
      }

      for(int i_group = 1; i_group < n_face_group; ++i_group) {
        _pborder_face_group_idx[i_group+1] += _pborder_face_group_idx[i_group];
      }

      // printf(" _pborder_face_group_idx[%i] = %i\n", n_face_group, _pborder_face_group_idx[n_face_group]);

      part_ext->border_face_group         [shift_part+i_part] = malloc( _pborder_face_group_idx[n_face_group] * sizeof(int        ));
      part_ext->border_face_group_ln_to_gn[shift_part+i_part] = malloc( _pborder_face_group_idx[n_face_group] * sizeof(PDM_g_num_t));

      int         *_pborder_face_group          = part_ext->border_face_group         [shift_part+i_part];
      PDM_g_num_t *_pborder_face_group_ln_to_gn = part_ext->border_face_group_ln_to_gn[shift_part+i_part];

      int* pborder_face_group_n = PDM_array_zeros_int(n_face_group);

      // PDM_log_trace_array_long(part_ext->border_face_ln_to_gn[shift_part+i_part], n_face_border, "border_face_ln_to_gn::");

      idx = 0;
      for(int i = 0; i < n_face_border; ++i) {
        for(int j = 0; j < border_face_group_idg_n[shift_part+i_part][i]; ++j) {
          int i_group = border_face_group_idg[shift_part+i_part][idx];
          int idx_write = _pborder_face_group_idx[i_group] + pborder_face_group_n[i_group]++;
          // _pborder_face_group         [idx_write] = pn_face+idx; // NON car si on a plusieurs group
          _pborder_face_group         [idx_write] = pn_face+i+1;
          _pborder_face_group_ln_to_gn[idx_write] = border_face_group_ln_to_gn[shift_part+i_part][idx];

          PDM_g_num_t g_num_face = border_face_ln_to_gn_check[shift_part+i_part][idx];
          int pos = PDM_binary_search_long(g_num_face, part_ext->border_face_ln_to_gn[shift_part+i_part], n_face_border);

          // printf("Find "PDM_FMT_G_NUM" pos = %i\n", g_num_face, pos);
          // printf(" check_ln_to_gn[%i] = "PDM_FMT_G_NUM" | "PDM_FMT_G_NUM" \n", idx, g_num_face, part_ext->border_face_ln_to_gn[shift_part+i_part][pos] );
          assert(g_num_face == part_ext->border_face_ln_to_gn[shift_part+i_part][pos]);
          idx++;
        }
      }

      free(pborder_face_group_n);
    }
    shift_part += part_ext->n_part[i_domain];
  }

  for(int i = 0; i < n_part_loc_all_domain; ++i) {
    free(border_face_group_idg_n[i]);
    free(border_face_group_idg[i]);
    free(border_face_group_ln_to_gn_n[i]);
    free(border_face_group_ln_to_gn[i]);
    free(border_face_ln_to_gn_check_n[i]);
    free(border_face_ln_to_gn_check[i]);
    free(face_group_n[i]);
    free(face_group_idg[i]);
    free(face_group_ln_to_gn[i]);
    free(face_ln_to_gn_check[i]);
  }
  free(border_face_group_idg_n);
  free(border_face_group_idg);
  free(border_face_group_ln_to_gn_n);
  free(border_face_group_ln_to_gn);
  free(border_face_ln_to_gn_check_n);
  free(border_face_ln_to_gn_check);


  free(n_face      );
  free(face_group_idg);
  free(face_group_n);
  free(face_group_ln_to_gn);
  free(face_ln_to_gn_check);

}


/*=============================================================================
 * Public function definitions
 *============================================================================*/


/**
 *
 * \brief Create a part extension structure
 *
 * \param [in]   comm           Communicator
 *
 */

PDM_part_extension_t*
PDM_part_extension_create
(
 const int                n_domain,
 const int               *n_part,
       PDM_extend_type_t  extend_type,
       int                depth,
 const PDM_MPI_Comm       comm,
 const PDM_ownership_t    owner
)
{
  PDM_part_extension_t *part_ext = (PDM_part_extension_t *) malloc(sizeof(PDM_part_extension_t));

  part_ext->n_domain    = n_domain;
  part_ext->n_part      = malloc( n_domain * sizeof(int)); // Make a copy to avoid pb in cython
  for(int i = 0; i < part_ext->n_domain; ++i) {
    part_ext->n_part[i] = n_part[i];
  }
  part_ext->comm        = comm;
  part_ext->owner       = owner;
  part_ext->extend_type = extend_type;
  part_ext->depth       = depth;

  part_ext->n_part_idx  = (int * ) malloc( (n_domain + 1) * sizeof(int));
  part_ext->parts = malloc(n_domain * sizeof(_part_t *));
  part_ext->n_part_idx[0] = 0;
  for(int i_domain = 0; i_domain < n_domain; ++i_domain) {
    part_ext->parts[i_domain] = malloc( n_part[i_domain] * sizeof(_part_t));
    part_ext->n_part_idx[i_domain+1] = part_ext->n_part_idx[i_domain] + part_ext->n_part[i_domain];
  }

  part_ext->neighbor_idx     = NULL;
  part_ext->neighbor_desc    = NULL;
  part_ext->n_entity_bound   = NULL;
  part_ext->border_cell_list = NULL;

  part_ext->dist_neighbor_cell_n              = NULL;
  part_ext->dist_neighbor_cell_idx            = NULL;
  part_ext->dist_neighbor_cell_desc           = NULL;
  part_ext->unique_order_dist_neighbor_cell   = NULL;
  part_ext->n_unique_order_dist_neighbor_cell = NULL;

  part_ext->n_tot_part_by_domain = NULL;

  part_ext->entity_cell_idx     = NULL;
  part_ext->entity_cell         = NULL;
  part_ext->entity_cell_opp_idx = NULL;
  part_ext->entity_cell_opp     = NULL;

  part_ext->cell_cell_extended_idx            = NULL;
  part_ext->cell_cell_extended_n              = NULL;
  part_ext->cell_cell_extended                = NULL;
  part_ext->unique_order_cell_cell_extended   = NULL;
  part_ext->n_unique_order_cell_cell_extended = NULL;
  part_ext->cell_cell_extended_pruned_idx     = NULL;
  part_ext->cell_cell_extended_pruned         = NULL;

  part_ext->face_face_extended_idx        = NULL;
  part_ext->face_face_extended            = NULL;

  part_ext->edge_edge_extended_idx        = NULL;
  part_ext->edge_edge_extended            = NULL;

  part_ext->vtx_vtx_extended_idx          = NULL;
  part_ext->vtx_vtx_extended              = NULL;

  part_ext->border_cell_face_idx          = NULL;
  part_ext->border_cell_face              = NULL;

  part_ext->border_face_edge_idx          = NULL;
  part_ext->border_face_edge              = NULL;

  part_ext->border_edge_vtx_idx           = NULL;
  part_ext->border_edge_vtx               = NULL;

  part_ext->border_face_vtx_idx           = NULL;
  part_ext->border_face_vtx               = NULL;

  part_ext->border_face_group_idx         = NULL;
  part_ext->border_face_group             = NULL;

  part_ext->border_vtx                    = NULL;

  part_ext->border_cell_ln_to_gn          = NULL;
  part_ext->border_face_ln_to_gn          = NULL;
  part_ext->border_edge_ln_to_gn          = NULL;
  part_ext->border_vtx_ln_to_gn           = NULL;
  part_ext->border_face_group_ln_to_gn    = NULL;

  return part_ext;
}

void
PDM_part_extension_set_part
(
  PDM_part_extension_t *part_ext,
  int                   i_domain,
  int                   i_part,
  int                   n_cell,
  int                   n_face,
  int                   n_face_part_bound,
  int                   n_face_group,
  int                   n_edge,
  int                   n_vtx,
  int                  *cell_face_idx,
  int                  *cell_face,
  int                  *face_cell,
  int                  *face_edge_idx,
  int                  *face_edge,
  int                  *face_vtx_idx,
  int                  *face_vtx,
  int                  *edge_vtx,
  int                  *face_bound_idx,
  int                  *face_bound,
  int                  *face_join_idx,
  int                  *face_join,
  int                  *face_part_bound_proc_idx,
  int                  *face_part_bound_part_idx,
  int                  *face_part_bound,
  int                  *vtx_part_bound_proc_idx,
  int                  *vtx_part_bound_part_idx,
  int                  *vtx_part_bound,
  PDM_g_num_t          *cell_ln_to_gn,
  PDM_g_num_t          *face_ln_to_gn,
  PDM_g_num_t          *edge_ln_to_gn,
  PDM_g_num_t          *vtx_ln_to_gn,
  PDM_g_num_t          *face_group_ln_to_gn,
  double               *vtx_coord
)
{

  part_ext->parts[i_domain][i_part].n_cell            = n_cell;
  part_ext->parts[i_domain][i_part].n_face            = n_face;
  part_ext->parts[i_domain][i_part].n_face_part_bound = n_face_part_bound;
  part_ext->parts[i_domain][i_part].n_face_group      = n_face_group;
  part_ext->parts[i_domain][i_part].n_edge            = n_edge;
  part_ext->parts[i_domain][i_part].n_vtx             = n_vtx;

  part_ext->parts[i_domain][i_part].cell_face_idx = cell_face_idx;
  part_ext->parts[i_domain][i_part].cell_face     = cell_face;
  part_ext->parts[i_domain][i_part].face_cell     = face_cell;

  part_ext->parts[i_domain][i_part].face_edge_idx = face_edge_idx;
  part_ext->parts[i_domain][i_part].face_edge     = face_edge;

  part_ext->parts[i_domain][i_part].face_vtx_idx  = face_vtx_idx;
  part_ext->parts[i_domain][i_part].face_vtx      = face_vtx;

  part_ext->parts[i_domain][i_part].edge_vtx      = edge_vtx;

  part_ext->parts[i_domain][i_part].cell_ln_to_gn = cell_ln_to_gn;
  part_ext->parts[i_domain][i_part].face_ln_to_gn = face_ln_to_gn;
  part_ext->parts[i_domain][i_part].edge_ln_to_gn = edge_ln_to_gn;
  part_ext->parts[i_domain][i_part].vtx_ln_to_gn  = vtx_ln_to_gn;

  part_ext->parts[i_domain][i_part].face_part_bound_proc_idx = face_part_bound_proc_idx;
  part_ext->parts[i_domain][i_part].face_part_bound_part_idx = face_part_bound_part_idx;
  part_ext->parts[i_domain][i_part].face_part_bound          = face_part_bound;

  part_ext->parts[i_domain][i_part].vtx_part_bound_proc_idx = vtx_part_bound_proc_idx;
  part_ext->parts[i_domain][i_part].vtx_part_bound_part_idx = vtx_part_bound_part_idx;
  part_ext->parts[i_domain][i_part].vtx_part_bound          = vtx_part_bound;

  part_ext->parts[i_domain][i_part].face_bound_idx      = face_bound_idx;
  part_ext->parts[i_domain][i_part].face_bound          = face_bound;
  part_ext->parts[i_domain][i_part].face_bound_ln_to_gn = face_group_ln_to_gn;

  part_ext->parts[i_domain][i_part].face_group_ln_to_gn = face_group_ln_to_gn;

  part_ext->parts[i_domain][i_part].face_join_idx       = face_join_idx;
  part_ext->parts[i_domain][i_part].face_join           = face_join;

  part_ext->parts[i_domain][i_part].vtx = vtx_coord;
}


/**
 *
 * \brief Compute a part extension structure
 *
 * \param [in]   part_ext          PDM_part_extension_t
 *
 */
void
PDM_part_extension_compute
(
  PDM_part_extension_t *part_ext
)
{
  /*
   *  A prevoir : reconstruire les "ghost" sans toute la topologie
   *              par exemple pour l'algébre linéaire
   */
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(part_ext->comm, &i_rank);
  PDM_MPI_Comm_size(part_ext->comm, &n_rank);

  int depth = part_ext->depth;

  part_ext->n_tot_part_by_domain = (int *) malloc( part_ext->n_domain * sizeof(int));
  int n_part_loc_all_domain = 0;
  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {

    part_ext->n_tot_part_by_domain[i_domain] = -1;
    int n_part_loc = part_ext->n_part[i_domain];
    PDM_MPI_Allreduce(&n_part_loc, &part_ext->n_tot_part_by_domain[i_domain], 1, PDM_MPI_INT, PDM_MPI_SUM, part_ext->comm);

    n_part_loc_all_domain += part_ext->n_part[i_domain];
  }

  part_ext->entity_cell_n    = (int **) malloc( n_part_loc_all_domain * sizeof(int *));
  part_ext->entity_cell_idx  = (int **) malloc( n_part_loc_all_domain * sizeof(int *));
  part_ext->entity_cell      = (int **) malloc( n_part_loc_all_domain * sizeof(int *));
  part_ext->n_cell           = (int  *) malloc( n_part_loc_all_domain * sizeof(int  ));
  part_ext->n_cell_border    = (int  *) malloc( n_part_loc_all_domain * sizeof(int  ));
  part_ext->border_cell_list = (int **) malloc( n_part_loc_all_domain * sizeof(int *));

  assert(part_ext != NULL);
  // printf(" PDM_part_extension_compute : depth = %i | extend_type = %i \n", depth, part_ext->extend_type);

  // assert(part_ext->extend_type == PDM_EXTEND_FROM_FACE);

  _create_cell_cell_graph(part_ext, part_ext->extend_type);

  part_ext->cell_cell_extended_idx            = (int *** ) malloc( (depth + 1) * sizeof(int  **));
  part_ext->cell_cell_extended_n              = (int *** ) malloc( (depth + 1) * sizeof(int  **));
  part_ext->cell_cell_extended                = (int *** ) malloc( (depth + 1) * sizeof(int  **));
  part_ext->unique_order_cell_cell_extended   = (int *** ) malloc( (depth + 1) * sizeof(int  **));
  part_ext->n_unique_order_cell_cell_extended = (int  ** ) malloc( (depth + 1) * sizeof(int   *));

  for(int i_depth = 0; i_depth < depth+1; ++i_depth) {
    part_ext->cell_cell_extended_idx           [i_depth] = (int **) malloc( n_part_loc_all_domain * sizeof(int *));
    part_ext->cell_cell_extended_n             [i_depth] = (int **) malloc( n_part_loc_all_domain * sizeof(int *));
    part_ext->cell_cell_extended               [i_depth] = (int **) malloc( n_part_loc_all_domain * sizeof(int *));
    part_ext->unique_order_cell_cell_extended  [i_depth] = (int **) malloc( n_part_loc_all_domain * sizeof(int *));
    part_ext->n_unique_order_cell_cell_extended[i_depth] = (int  *) malloc( n_part_loc_all_domain * sizeof(int  ));
  }

  part_ext->cell_cell_extended_pruned_idx = (int **) malloc( n_part_loc_all_domain * sizeof(int *));
  part_ext->cell_cell_extended_pruned     = (int **) malloc( n_part_loc_all_domain * sizeof(int *));

  // TODO : vtx_cell
  _create_cell_graph_comm(part_ext);

  /*
   * Create for first level the proper graph
   */
  _compute_first_extended_cell_graph(part_ext);


  /*
   *   Step 3 : Compute the graph cell with triplet
   *      -> Now we have all the things ok to exchange direclty on cells
   *      -> In order to build the next depth of ghost cell
   *         we need to prepare an extended graph with triplet
   *      -> With the same distant neigbor we exchange for each depth this graph (containing the information of surrounding cells)
   *   --> Init the distant_neightbor from cell
   */
  PDM_distant_neighbor_t* dn = PDM_distant_neighbor_create(part_ext->comm,
                                                           n_part_loc_all_domain,
                                                           part_ext->n_cell,
                                                           part_ext->dist_neighbor_cell_idx,
                                                           part_ext->dist_neighbor_cell_desc);

  for(int i_depth = 1; i_depth < depth+1; ++i_depth) {
    /* Graph compute = local + distant */
    _compute_dual_graph(part_ext, dn, i_depth);
  }

  PDM_distant_neighbor_free(dn);

  /*
   * At this step we have for each level the opposite connected cell
   *   In order to setup all ghost cell we need to deduced all descending connectivity
   *       cell -> face -> edge -> vtx
   * Another step is to compute the ln_to_gn in mulitpart context -> (i_domain, ln_to_gn)
   */

  // _prune_cell_cell_extented(part_ext, depth);
  _prune_cell_cell_extented(part_ext, depth-1);

  // if( 0 == 1) {
  //   _rebuild_connectivity_cell_face_debug(part_ext);
  // }
  if(part_ext->extend_type == PDM_EXTEND_FROM_FACE) {
    _rebuild_connectivity_cell_face(part_ext);
    _rebuild_connectivity_face_vtx(part_ext);
    // _rebuild_connectivity_face_edge(part_ext);
    // _rebuild_connectivity_edge_vtx (part_ext);
  } else if (part_ext->extend_type == PDM_EXTEND_FROM_VTX) {
    _rebuild_connectivity_cell_face(part_ext);
    _rebuild_connectivity_face_edge(part_ext);
    _rebuild_connectivity_edge_vtx (part_ext);
  } else {
    PDM_error(__FILE__, __LINE__, 0, "PDM_part_extension_compute wrong extend_type \n");
  }

  int     *n_vtx     = (int     * ) malloc( n_part_loc_all_domain * sizeof(int     ));
  double **vtx_coord = (double ** ) malloc( n_part_loc_all_domain * sizeof(double *));

  int shift_part = 0;
  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
    for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {
      n_vtx        [i_part+shift_part] = part_ext->parts[i_domain][i_part].n_vtx;
      vtx_coord    [i_part+shift_part] = part_ext->parts[i_domain][i_part].vtx;
    }
    shift_part += part_ext->n_part[i_domain];
  }

  /* Finalize with vertex */
  PDM_distant_neighbor_t* dn_vtx = PDM_distant_neighbor_create(part_ext->comm,
                                                               n_part_loc_all_domain,
                                                               n_vtx,
                                                               part_ext->vtx_vtx_extended_idx,
                                                               part_ext->vtx_vtx_extended);


  assert(part_ext->border_vtx == NULL);


  PDM_distant_neighbor_exch(dn_vtx,
                            sizeof(double),
                            PDM_STRIDE_CST_INTERLACED,
                            3,
                            NULL,
                 (void **)  vtx_coord,
                            NULL,
                (void ***) &part_ext->border_vtx);


  // // Debug
  // PDM_g_num_t **gnum_check = (PDM_g_num_t ** ) malloc( n_part_loc_all_domain * sizeof(PDM_g_num_t *));
  // PDM_g_num_t **stride  = (int ** ) malloc( n_part_loc_all_domain * sizeof(int *));

  // shift_part = 0;
  // for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
  //   for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {
  //     n_vtx        [i_part+shift_part] = part_ext->parts[i_domain][i_part].n_vtx;
  //     gnum_check   [i_part+shift_part] = part_ext->parts[i_domain][i_part].vtx_ln_to_gn;

  //     stride    [i_part+shift_part] = (int *) malloc( part_ext->parts[i_domain][i_part].n_vtx * sizeof(int ));

  //     for( int i_vtx = 0; i_vtx < part_ext->parts[i_domain][i_part].n_vtx; i_vtx++) {
  //        stride    [i_part+shift_part][i_vtx] = 1;
  //     }


  //     PDM_log_trace_array_long(gnum_check   [i_part+shift_part], part_ext->parts[i_domain][i_part].n_vtx        , "vtx_ln_to_gn::");
  //     PDM_log_trace_array_int(part_ext->vtx_vtx_extended_idx[shift_part+i_part], part_ext->parts[i_domain][i_part].n_vtx        , "vtx_vtx_extended_idx::");
  //     PDM_log_trace_array_int(part_ext->vtx_vtx_extended    [shift_part+i_part], 3 * part_ext->vtx_vtx_extended_idx[shift_part+i_part][part_ext->parts[i_domain][i_part].n_vtx], "vtx_vtx_extended_idx::");
  //   }
  //   shift_part += part_ext->n_part[i_domain];
  // }

  // PDM_g_num_t **pgnum_check;
  // int  **recv_strid;
  // // PDM_distant_neighbor_exch(dn_vtx,
  // //                           sizeof(PDM_g_num_t),
  // //                           PDM_STRIDE_CST,
  // //                           1,
  // //                           NULL,
  // //                (void **)  gnum_check,
  // //                           NULL,
  // //               (void ***) &pgnum_check);
  // PDM_distant_neighbor_exch(dn_vtx,
  //                           sizeof(PDM_g_num_t),
  //                           PDM_STRIDE_VAR,
  //                           -1,
  //                           stride,
  //                (void **)  gnum_check,
  //                           &recv_strid,
  //               (void ***) &pgnum_check);

  // shift_part = 0;
  // for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
  //   for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {
  //     for(int i_vtx = 0; i_vtx < part_ext->vtx_vtx_extended_idx[shift_part+i_part][part_ext->parts[i_domain][i_part].n_vtx]; ++i_vtx){
  //       printf(" [%i] pgnum_check[%i] = %i \n", i_part, i_vtx, pgnum_check[i_part][i_vtx]);
  //     }
  //   }
  // }


  free(n_vtx);
  free(vtx_coord);

  PDM_distant_neighbor_free(dn_vtx);
  // printf(" PDM_part_extension_compute end \n");

  /* Condition limite - Face uniquement pour l'instant */
  _rebuild_face_group(part_ext);

}


/**
 *
 * \brief Free a part extension structure
 *
 * \param [in]   part_ext          PDM_part_extension_t
 *
 */
void
PDM_part_extension_free
(
 PDM_part_extension_t *part_ext
)
{

  if(part_ext->n_tot_part_by_domain != NULL) {
    free(part_ext->n_tot_part_by_domain);
    part_ext->n_tot_part_by_domain = NULL;
  }

  if(part_ext->neighbor_idx != NULL) {
    int shift_part = 0;
    for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
      for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {
        free(part_ext->neighbor_idx       [i_part+shift_part]);
        free(part_ext->neighbor_desc      [i_part+shift_part]);

        free(part_ext->dist_neighbor_cell_n           [i_part+shift_part]);
        free(part_ext->dist_neighbor_cell_idx         [i_part+shift_part]);
        free(part_ext->dist_neighbor_cell_desc        [i_part+shift_part]);
        free(part_ext->unique_order_dist_neighbor_cell[i_part+shift_part]);

        free(part_ext->entity_cell_opp_idx[i_part+shift_part]);
        free(part_ext->entity_cell_opp_n  [i_part+shift_part]);
        free(part_ext->entity_cell_opp    [i_part+shift_part]);

        free(part_ext->border_cell_list    [i_part+shift_part]);

        // if(part_ext->extend_type == PDM_EXTEND_FROM_FACE) {
        free(part_ext->entity_cell_idx[i_part+shift_part]);
        free(part_ext->entity_cell_n  [i_part+shift_part]);
        free(part_ext->entity_cell    [i_part+shift_part]);
        // }

        if(part_ext->face_face_extended_idx != NULL) {
          free(part_ext->face_face_extended_idx[i_part+shift_part]);
        }

        if(part_ext->face_face_extended != NULL) {
          free(part_ext->face_face_extended[i_part+shift_part]);
        }

        if(part_ext->edge_edge_extended_idx != NULL) {
          free(part_ext->edge_edge_extended_idx[i_part+shift_part]);
        }

        if(part_ext->edge_edge_extended != NULL) {
          free(part_ext->edge_edge_extended[i_part+shift_part]);
        }

        if(part_ext->vtx_vtx_extended_idx != NULL) {
          free(part_ext->vtx_vtx_extended_idx[i_part+shift_part]);
        }

        if(part_ext->vtx_vtx_extended != NULL) {
          free(part_ext->vtx_vtx_extended[i_part+shift_part]);
        }

        if(part_ext->owner == PDM_OWNERSHIP_KEEP) {

          if(part_ext->border_cell_face_idx != NULL) {
            free(part_ext->border_cell_face_idx[i_part+shift_part]);
            free(part_ext->border_cell_face    [i_part+shift_part]);
          }

          if(part_ext->border_face_edge_idx != NULL) {
            free(part_ext->border_face_edge_idx[i_part+shift_part]);
            free(part_ext->border_face_edge    [i_part+shift_part]);
          }

          if(part_ext->border_edge_vtx_idx != NULL) {
            free(part_ext->border_edge_vtx_idx[i_part+shift_part]);
            free(part_ext->border_edge_vtx    [i_part+shift_part]);
          }

          if(part_ext->border_face_vtx_idx != NULL) {
            free(part_ext->border_face_vtx_idx[i_part+shift_part]);
            free(part_ext->border_face_vtx    [i_part+shift_part]);
          }

          if(part_ext->border_face_group_idx != NULL) {
            free(part_ext->border_face_group_idx[i_part+shift_part]);
            free(part_ext->border_face_group    [i_part+shift_part]);
          }

          if(part_ext->border_vtx != NULL) {
            free(part_ext->border_vtx[i_part+shift_part]);
          }

          if(part_ext->border_cell_ln_to_gn != NULL) {
            free(part_ext->border_cell_ln_to_gn[i_part+shift_part]);
          }

          if(part_ext->border_face_ln_to_gn != NULL) {
            free(part_ext->border_face_ln_to_gn[i_part+shift_part]);
          }

          if(part_ext->border_edge_ln_to_gn != NULL) {
            free(part_ext->border_edge_ln_to_gn[i_part+shift_part]);
          }

          if(part_ext->border_vtx_ln_to_gn != NULL) {
            free(part_ext->border_vtx_ln_to_gn[i_part+shift_part]);
          }

          if(part_ext->border_face_group_ln_to_gn != NULL) {
            free(part_ext->border_face_group_ln_to_gn[i_part+shift_part]);
          }

        }

        free(part_ext->cell_cell_idx[i_domain][i_part]);
        free(part_ext->cell_cell    [i_domain][i_part]);

      }
      free(part_ext->cell_cell_idx[i_domain]);
      free(part_ext->cell_cell    [i_domain]);

      shift_part += part_ext->n_part[i_domain];
    }
    free(part_ext->n_entity_bound);
  }
  free(part_ext->neighbor_idx);
  free(part_ext->neighbor_desc);
  free(part_ext->n_cell);
  free(part_ext->n_cell_border);
  free(part_ext->border_cell_list);

  free(part_ext->cell_cell_idx);
  free(part_ext->cell_cell);

  if(part_ext->face_face_extended_idx != NULL) {
    free(part_ext->face_face_extended_idx);
    free(part_ext->face_face_extended);
  }

  if(part_ext->edge_edge_extended_idx != NULL) {
    free(part_ext->edge_edge_extended_idx);
    free(part_ext->edge_edge_extended);
  }

  if(part_ext->vtx_vtx_extended_idx != NULL) {
    free(part_ext->vtx_vtx_extended_idx);
    free(part_ext->vtx_vtx_extended);
  }

  /* Peu import l'ownership on free car on rend à l'utilisateur l'interface i_domain / i_part */
  if(part_ext->border_cell_face_idx != NULL) {
    free(part_ext->border_cell_face_idx);
    free(part_ext->border_cell_face);
  }

  if(part_ext->border_face_edge_idx != NULL) {
    free(part_ext->border_face_edge_idx);
    free(part_ext->border_face_edge);
  }

  if(part_ext->border_edge_vtx_idx != NULL) {
    free(part_ext->border_edge_vtx_idx);
    free(part_ext->border_edge_vtx);
  }

  if(part_ext->border_face_vtx_idx != NULL) {
    free(part_ext->border_face_vtx_idx);
    free(part_ext->border_face_vtx);
  }

  if(part_ext->border_face_group_idx != NULL) {
    free(part_ext->border_face_group_idx);
    free(part_ext->border_face_group);
  }

  if(part_ext->border_vtx != NULL) {
    free(part_ext->border_vtx);
  }

  if(part_ext->border_cell_ln_to_gn != NULL) {
    free(part_ext->border_cell_ln_to_gn);
  }

  if(part_ext->border_face_ln_to_gn != NULL) {
    free(part_ext->border_face_ln_to_gn);
  }

  if(part_ext->border_edge_ln_to_gn != NULL) {
    free(part_ext->border_edge_ln_to_gn);
  }

  if(part_ext->border_vtx_ln_to_gn != NULL) {
    free(part_ext->border_vtx_ln_to_gn);
  }

  if(part_ext->border_face_group_ln_to_gn != NULL) {
    free(part_ext->border_face_group_ln_to_gn);
  }

  free(part_ext->n_part_idx);

  for(int i_depth = 0; i_depth < part_ext->depth+1; ++i_depth) {

    int shift_part = 0;
    for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
      for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {
        free(part_ext->cell_cell_extended_idx         [i_depth][i_part+shift_part]);
        free(part_ext->cell_cell_extended_n           [i_depth][i_part+shift_part]);
        free(part_ext->cell_cell_extended             [i_depth][i_part+shift_part]);
        free(part_ext->unique_order_cell_cell_extended[i_depth][i_part+shift_part]);
      }
      shift_part += part_ext->n_part[i_domain];
    }
    free(part_ext->cell_cell_extended_idx           [i_depth]);
    free(part_ext->cell_cell_extended_n             [i_depth]);
    free(part_ext->cell_cell_extended               [i_depth]);
    free(part_ext->unique_order_cell_cell_extended  [i_depth]);
    free(part_ext->n_unique_order_cell_cell_extended[i_depth]);
  }

  free(part_ext->cell_cell_extended_idx);
  free(part_ext->cell_cell_extended_n);
  free(part_ext->cell_cell_extended);
  free(part_ext->unique_order_cell_cell_extended);
  free(part_ext->n_unique_order_cell_cell_extended);

  /* Only shortcut of user data */
  free(part_ext->entity_cell_idx    );
  free(part_ext->entity_cell        );
  free(part_ext->entity_cell_n);

  free(part_ext->dist_neighbor_cell_n   );
  free(part_ext->dist_neighbor_cell_idx );
  free(part_ext->dist_neighbor_cell_desc);
  free(part_ext->unique_order_dist_neighbor_cell);
  free(part_ext->n_unique_order_dist_neighbor_cell);

  /* Allocated by distant neightbor */
  free(part_ext->entity_cell_opp_idx);
  free(part_ext->entity_cell_opp_n  );
  free(part_ext->entity_cell_opp    );

  int shift_part = 0;
  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
    for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {
      free(part_ext->cell_cell_extended_pruned_idx[i_part+shift_part]);
      free(part_ext->cell_cell_extended_pruned    [i_part+shift_part]);
    }
    shift_part += part_ext->n_part[i_domain];
  }
  free(part_ext->cell_cell_extended_pruned_idx);
  free(part_ext->cell_cell_extended_pruned    );

  free(part_ext->n_part);
  part_ext->n_part = NULL;

  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
    free(part_ext->parts[i_domain]);
  }
  free(part_ext->parts);


  free(part_ext);
}


int
PDM_part_extension_connectivity_get
(
 PDM_part_extension_t     *part_ext,
 int                       i_domain,
 int                       i_part,
 PDM_connectivity_type_t   connectivity_type,
 int                     **connect,
 int                     **connect_idx
)
{
  int n_entity = -1;
  int shift_part = part_ext->n_part_idx[i_domain];
  switch(connectivity_type)
  {
    case PDM_CONNECTIVITY_TYPE_CELL_FACE:
    {
      int n_cell   = part_ext->parts[i_domain][i_part].n_cell;
      n_entity     = part_ext->cell_cell_extended_pruned_idx[shift_part+i_part][n_cell];
      *connect_idx = part_ext->border_cell_face_idx         [shift_part+i_part];
      *connect     = part_ext->border_cell_face             [shift_part+i_part];
    }
    break;

    case PDM_CONNECTIVITY_TYPE_FACE_EDGE:
    {
      int n_face   = part_ext->parts[i_domain][i_part].n_face;
      n_entity     = part_ext->face_face_extended_idx[shift_part+i_part][n_face];
      *connect_idx = part_ext->border_face_edge_idx  [shift_part+i_part];
      *connect     = part_ext->border_face_edge      [shift_part+i_part];
    }
    break;

    case PDM_CONNECTIVITY_TYPE_FACE_VTX:
    {
      int n_face   = part_ext->parts[i_domain][i_part].n_face;
      n_entity     = part_ext->face_face_extended_idx[shift_part+i_part][n_face];
      *connect_idx = part_ext->border_face_vtx_idx   [shift_part+i_part];
      *connect     = part_ext->border_face_vtx       [shift_part+i_part];
    }
    break;

    case PDM_CONNECTIVITY_TYPE_EDGE_VTX:
    {
      int n_edge   = part_ext->parts[i_domain][i_part].n_edge;
      n_entity     = part_ext->edge_edge_extended_idx[shift_part+i_part][n_edge];
      *connect_idx = part_ext->border_edge_vtx_idx   [shift_part+i_part];
      *connect     = part_ext->border_edge_vtx       [shift_part+i_part];
    }
    break;

  default :
    PDM_error(__FILE__, __LINE__, 0, "Unknown connectivity_type \n");
    break;

  }

  return n_entity;
}


int
PDM_part_extension_ln_to_gn_get
(
 PDM_part_extension_t     *part_ext,
 int                       i_domain,
 int                       i_part,
 PDM_mesh_entities_t       mesh_entity,
 PDM_g_num_t             **ln_to_gn
)
{
  int n_entity = -1;
  int shift_part = part_ext->n_part_idx[i_domain];
  switch(mesh_entity)
  {
    case PDM_MESH_ENTITY_CELL:
    {
      int n_cell   = part_ext->parts[i_domain][i_part].n_cell;
      n_entity     = part_ext->cell_cell_extended_pruned_idx[shift_part+i_part][n_cell];
      *ln_to_gn    = part_ext->border_cell_ln_to_gn         [shift_part+i_part];
    }
    break;

    case PDM_MESH_ENTITY_FACE:
    {
      int n_face   = part_ext->parts[i_domain][i_part].n_face;
      n_entity     = part_ext->face_face_extended_idx[shift_part+i_part][n_face];
      *ln_to_gn    = part_ext->border_face_ln_to_gn  [shift_part+i_part];
    }
    break;

    case PDM_MESH_ENTITY_EDGE:
    {
      int n_edge   = part_ext->parts[i_domain][i_part].n_edge;
      n_entity     = part_ext->edge_edge_extended_idx[shift_part+i_part][n_edge];
      *ln_to_gn    = part_ext->border_edge_ln_to_gn  [shift_part+i_part];
    }
    break;

    case PDM_MESH_ENTITY_VERTEX:
    {
      int n_vtx   = part_ext->parts[i_domain][i_part].n_vtx;
      n_entity     = part_ext->vtx_vtx_extended_idx[shift_part+i_part][n_vtx];
      *ln_to_gn    = part_ext->border_vtx_ln_to_gn  [shift_part+i_part];
    }
    break;

  default :
    PDM_error(__FILE__, __LINE__, 0, "Unknown connectivity_type \n");
    break;

  }

  return n_entity;
}

int
PDM_part_extension_group_get
(
 PDM_part_extension_t     *part_ext,
 int                       i_domain,
 int                       i_part,
 PDM_mesh_entities_t       mesh_entity,
 int                     **connect,
 int                     **connect_idx,
 PDM_g_num_t             **ln_to_gn
)
{
  int n_entity = -1;
  int shift_part = part_ext->n_part_idx[i_domain];
  switch(mesh_entity)
  {
    case PDM_MESH_ENTITY_CELL:
    {
      abort();
    }
    break;

    case PDM_MESH_ENTITY_FACE:
    {
      int n_face_group = part_ext->parts[i_domain][i_part].n_face_group;
      n_entity     = n_face_group;
      *connect_idx = part_ext->border_face_group_idx[shift_part+i_part];
      *connect     = part_ext->border_face_group    [shift_part+i_part];
      *ln_to_gn    = part_ext->border_face_group_ln_to_gn [shift_part+i_part];
    }
    break;

    case PDM_MESH_ENTITY_EDGE:
    {
      abort();
    }
    break;

    case PDM_MESH_ENTITY_VERTEX:
    {
      abort();
    }
    break;

  default :
    PDM_error(__FILE__, __LINE__, 0, "Unknown connectivity_type \n");
    break;

  }

  return n_entity;
}

int
PDM_part_extension_coord_get
(
 PDM_part_extension_t     *part_ext,
 int                       i_domain,
 int                       i_part,
 double                  **vtx_coord
)
{
  int shift_part     = part_ext->n_part_idx[i_domain];
  int n_vtx          = part_ext->parts[i_domain][i_part].n_vtx;
  int n_vtx_extended = part_ext->vtx_vtx_extended_idx[shift_part+i_part][n_vtx];
  *vtx_coord = part_ext->border_vtx[shift_part+i_part];

  return n_vtx_extended;
}

#ifdef __cplusplus
}
#endif /* __cplusplus */

