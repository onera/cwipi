
/*----------------------------------------------------------------------------
 *  System headers
 *----------------------------------------------------------------------------*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_error.h"
#include "pdm_binary_search.h"
#include "pdm_dmesh_nodal_to_dmesh_priv.h"
#include "pdm_dmesh_nodal_priv.h"
#include "pdm_dconnectivity_transform.h"
#include "pdm_dmesh_nodal_to_dmesh.h"
#include "pdm_dmesh_nodal_elements_utils.h"
#include "pdm_part_to_block.h"
#include "pdm_block_to_part.h"
#include "pdm_logging.h"
#include "pdm_dmesh.h"
#include "pdm_distrib.h"
#include "pdm_quick_sort.h"
#include "pdm_para_graph_dual.h"

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

/*============================================================================
 * Global variable
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

static
int
_is_a_3d_elment
(
 PDM_Mesh_nodal_elt_t t_elt
)
{
  switch (t_elt) {
  case PDM_MESH_NODAL_TETRA4   :
  case PDM_MESH_NODAL_PYRAMID5 :
  case PDM_MESH_NODAL_PRISM6   :
  case PDM_MESH_NODAL_HEXA8    :
    {

      return 1;
    }
    break;

  case PDM_MESH_NODAL_POLY_2D  :
    {
      return 0;
    }
    break;

  case PDM_MESH_NODAL_POLY_3D  :
    {
      return 1;
    }
    break;


  default :
    return 0;
    break;
  }

}

static
int
_is_a_2d_elment
(
 PDM_Mesh_nodal_elt_t t_elt
)
{
  switch (t_elt) {
  case PDM_MESH_NODAL_TRIA3    :
  case PDM_MESH_NODAL_QUAD4    :
    {
      return 1;
    }
    break;

  case PDM_MESH_NODAL_POLY_2D  :
    {
      return 1;
    }
    break;

  case PDM_MESH_NODAL_POLY_3D  :
    {
      return 0;
    }
    break;

  default :
    return 0;
    break;
  }

}

/**
 * \def _compute_keys
 */
static
void
_compute_keys
(
const int          n_face_elt_tot,
const int         *delmt_face_vtx_idx,
const PDM_g_num_t *delmt_face_vtx,
      PDM_g_num_t *ln_to_gn,
      PDM_g_num_t  key_mod
)
{
  for(int i_face = 0; i_face < n_face_elt_tot; ++i_face ) {
    PDM_g_num_t key = 0;
    for(int idx = delmt_face_vtx_idx[i_face]; idx < delmt_face_vtx_idx[i_face+1]; ++idx) {
      key += delmt_face_vtx[idx];
    }
    // min_vtx =
    ln_to_gn[i_face] = key % key_mod + 1;
  }
}

static
PDM_g_num_t*
_make_absolute_entity_numbering
(
       int          dn_entity,
 const PDM_MPI_Comm comm
)
{
  int n_rank;

  PDM_MPI_Comm_size(comm, &n_rank);

  PDM_g_num_t n_entity_proc = dn_entity;
  PDM_g_num_t beg_num_abs;

  PDM_MPI_Scan(&n_entity_proc, &beg_num_abs, 1, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, comm);
  beg_num_abs -= n_entity_proc;

  /** Compute the distribution of elements amont proc **/
  PDM_g_num_t *entity_distrib = (PDM_g_num_t *) malloc((n_rank+1) * sizeof(PDM_g_num_t));
  PDM_g_num_t _dn_face = (PDM_g_num_t) dn_entity;
  PDM_MPI_Allgather((void *) &_dn_face,
                    1,
                    PDM__PDM_MPI_G_NUM,
                    (void *) (&entity_distrib[1]),
                    1,
                    PDM__PDM_MPI_G_NUM,
                    comm);

  // entity_distrib[0] = 1;
  entity_distrib[0] = 0;

  for (int i = 1; i < n_rank+1; i++) {
    entity_distrib[i] +=  entity_distrib[i-1];
  }

  if (0 == 1) {
    printf("beg_num_abs::Face : "PDM_FMT_G_NUM" \n", beg_num_abs);
    printf("entity_distrib : "PDM_FMT_G_NUM,  entity_distrib[0]);
    for (int i = 1; i < n_rank+1; i++) {
      printf(" "PDM_FMT_G_NUM, entity_distrib[i]);
    }
    printf("\n");
  }
  return entity_distrib;
}

void
PDM_generate_entitiy_connectivity
(
PDM_MPI_Comm   comm,
PDM_g_num_t    n_vtx_abs,
int            n_entity_elt_tot,
PDM_g_num_t   *delmt_entity,
int           *delmt_entity_vtx_idx,
PDM_g_num_t   *delmt_entity_vtx,
int           *dn_entity,
PDM_g_num_t  **entity_distrib,
int          **dentity_vtx_idx,
PDM_g_num_t  **dentity_vtx,
int          **dentity_elmt_idx,
PDM_g_num_t  **dentity_elmt
)
{
  /*
   * We are now all information flatten - we only need to compute hash_keys for each entitys
   */
  PDM_g_num_t* ln_to_gn = (PDM_g_num_t*) malloc( n_entity_elt_tot * sizeof(PDM_g_num_t));
  PDM_g_num_t key_mod = 4 * n_vtx_abs;

  // Option des clés
  _compute_keys(n_entity_elt_tot,
                delmt_entity_vtx_idx,
                delmt_entity_vtx,
                ln_to_gn,
                key_mod);

  if(0 == 1) {
    log_trace("n_entity_elt_tot = %i \n", n_entity_elt_tot);
    PDM_log_trace_array_long(ln_to_gn, n_entity_elt_tot , "ln_to_gn:: ");
    PDM_log_trace_array_int(delmt_entity_vtx_idx  , n_entity_elt_tot+1 , "delmt_entity_vtx_idx:: ");
    PDM_log_trace_array_long(delmt_entity_vtx, delmt_entity_vtx_idx[n_entity_elt_tot] , "delmt_entity_vtx:: ");
  }

  /*
   * Prepare exchange by computing stride
   */
  int* delmt_entity_vtx_n = (int        *) malloc( n_entity_elt_tot * sizeof(int        ));
  int* stride_one         = (int        *) malloc( n_entity_elt_tot * sizeof(int        ));
  for(int i_entity = 0; i_entity < n_entity_elt_tot; ++i_entity) {
    delmt_entity_vtx_n[i_entity] = delmt_entity_vtx_idx[i_entity+1] - delmt_entity_vtx_idx[i_entity];
    stride_one[i_entity]       = 1;
  }
  free(delmt_entity_vtx_idx);

  // PDM_log_trace_array_int(delmt_entity_vtx_n  , n_entity_elt_tot , "delmt_entity_vtx_n:: ");

  // for(int i = 0; i < n_entity_elt_tot; ++i) {
  //   log_trace("delmt_entity_vtx_n[%i] = %i \n", i, delmt_entity_vtx_n[i]);
  // }
  /*
   * Setup part_to_block to filter all keys
   */
  double* weight = (double *) malloc( n_entity_elt_tot * sizeof(double));
  for(int i = 0; i < n_entity_elt_tot; ++i) {
    weight[i] = 1.;
  }

  PDM_part_to_block_t *ptb = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                      PDM_PART_TO_BLOCK_POST_MERGE,
                                                      1.,
                                                      &ln_to_gn,
                                                      &weight,
                                                      &n_entity_elt_tot,
                                                      1,
                                                      comm);
  free(weight);
  /*
   * Exchange data
   */
  int*         blk_tot_entity_vtx_n = NULL;
  PDM_g_num_t* blk_tot_entity_vtx   = NULL;

  int blk_tot_entity_vtx_size = PDM_part_to_block_exch(         ptb,
                                                              sizeof(PDM_g_num_t),
                                                              PDM_STRIDE_VAR,
                                                              -1,
                                                              &delmt_entity_vtx_n,
                                                    (void **) &delmt_entity_vtx,
                                                              &blk_tot_entity_vtx_n,
                                                    (void **) &blk_tot_entity_vtx);

  int* blk_n_entity_per_key = NULL;
  int* blk_entity_vtx_n     = NULL;
  int blk_entity_vtx_n_size = PDM_part_to_block_exch(         ptb,
                                                            sizeof(int),
                                                            PDM_STRIDE_VAR,
                                                            -1,
                                                            &stride_one,
                                                  (void **) &delmt_entity_vtx_n,
                                                            &blk_n_entity_per_key,
                                                  (void **) &blk_entity_vtx_n);
  free(delmt_entity_vtx_n);

  int*         blk_elmt_entity_elmt_stri = NULL;
  PDM_g_num_t* blk_elmt_entity_elmt      = NULL;
  int blk_entity_elmt_size = PDM_part_to_block_exch(         ptb,
                                                           sizeof(PDM_g_num_t),
                                                           PDM_STRIDE_VAR,
                                                           -1,
                                                           &stride_one,
                                                 (void **) &delmt_entity,
                                                           &blk_elmt_entity_elmt_stri,
                                                 (void **) &blk_elmt_entity_elmt);

  /*
   *  Get the size of the current process bloc
   */
  int blk_size = PDM_part_to_block_n_elt_block_get(ptb);

  if( 0 == 1 ) {
    PDM_log_trace_array_int(blk_tot_entity_vtx_n, blk_size             , "blk_tot_entity_vtx_n:: ");
    PDM_log_trace_array_long(blk_tot_entity_vtx , blk_tot_entity_vtx_size, "blk_tot_entity_vtx:: "  );

    PDM_log_trace_array_int(blk_n_entity_per_key, blk_size         , "blk_n_entity_per_key:: ");
    PDM_log_trace_array_int(blk_entity_vtx_n    , blk_entity_vtx_n_size, "blk_entity_vtx_n:: ");

    PDM_log_trace_array_long(blk_elmt_entity_elmt, blk_entity_elmt_size, "blk_elmt_entity_elmt:: ");
  }

  PDM_part_to_block_free(ptb);
  free(stride_one);
  free(blk_elmt_entity_elmt_stri); // Same as blk_n_entity_per_key
  free(blk_tot_entity_vtx_n);

  /*
   * Get the max number of vertex of entitys
   */
  int* blk_entity_vtx_idx  = (int        *) malloc( (blk_entity_vtx_n_size+1) * sizeof(int        ));
  int n_max_entity_per_key = -1;
  int n_tot_entity_per_key = 0;
  for(int i_entity = 0; i_entity < blk_size; ++i_entity) {
    n_max_entity_per_key = PDM_MAX(n_max_entity_per_key, blk_n_entity_per_key[i_entity]);
    n_tot_entity_per_key += blk_n_entity_per_key[i_entity];
  }

  int n_max_vtx       = -1;
  blk_entity_vtx_idx[0] = 0;
  for(int i_entity = 0; i_entity < blk_entity_vtx_n_size; ++i_entity) {
    n_max_vtx          = PDM_MAX(n_max_vtx         , blk_entity_vtx_n    [i_entity]);
    blk_entity_vtx_idx[i_entity+1] = blk_entity_vtx_idx[i_entity] + blk_entity_vtx_n[i_entity];
  }

  // PDM_log_trace_array_int(blk_entity_vtx_idx, blk_entity_vtx_n_size, "blk_entity_vtx_idx:: ");
  /*
   * We need to identify each uniques entitys :
   *      - We have multiple packet to treat
   *      - The connectivity can be sorted in place
   *      - Multiple case can occur :
   *           - Alone entity normaly boundary
   *           - Multiple entitys
   *           - Same entity, we remove and replace by the first
   */
  PDM_g_num_t* loc_entity_vtx_1 = (PDM_g_num_t *) malloc(  n_max_vtx               * sizeof(PDM_g_num_t) );
  PDM_g_num_t* loc_entity_vtx_2 = (PDM_g_num_t *) malloc(  n_max_vtx               * sizeof(PDM_g_num_t) );
  int*         already_treat    = (int         *) malloc(  n_max_entity_per_key    * sizeof(int        ) );
  int*         same_entity_idx  = (int         *) malloc( (n_max_entity_per_key+1) * sizeof(int        ) );
  int*         sens_entity      = (int         *) malloc(  n_max_entity_per_key    * sizeof(int        ) );

  /*
   * Allocate Memory - entity_vtx - entity_elmt
   */
  *dentity_vtx      = (PDM_g_num_t *) malloc( sizeof(PDM_g_num_t) *  blk_tot_entity_vtx_size  );
  *dentity_vtx_idx  = (int         *) malloc( sizeof(int        ) * (blk_entity_vtx_n_size+1 ));
  *dentity_elmt     = (PDM_g_num_t *) malloc( sizeof(PDM_g_num_t) *  n_tot_entity_per_key );
  *dentity_elmt_idx = (int         *) malloc( sizeof(int)         * (blk_entity_elmt_size+1)  );

  PDM_g_num_t *_dentity_vtx      = *dentity_vtx;
  int         *_dentity_vtx_idx  = *dentity_vtx_idx;
  PDM_g_num_t *_dentity_elmt     = *dentity_elmt;
  int         *_dentity_elmt_idx = *dentity_elmt_idx;

  // printf("blk_entity_elmt_size::%i\n", blk_entity_elmt_size);
  // printf("n_tot_entity_per_key::%i\n", n_tot_entity_per_key);

  /*
   * Init global numbering
   */
  int i_abs_entity = 0;
  _dentity_vtx_idx[0] = 0;

  int idx = 0;
  int idx_entity_vtx = 0;
  _dentity_elmt_idx[0] = 0;
  for(int i_key = 0; i_key < blk_size; ++i_key) {
    // printf(" --- Number of conflicting keys :: %i \n", blk_n_entity_per_key[i_key]);

    int n_conflict_entitys = blk_n_entity_per_key[i_key];

    /* Reset */
    for(int j = 0; j < n_conflict_entitys; ++j) {
      already_treat[j] = -1;
    }

    /* Loop over all entitys in conflict */
    for(int i_entity = 0; i_entity < n_conflict_entitys; ++i_entity) {
      // printf("Number of vtx on entitys %i :: %i with index [%i] \n", i_entity, blk_entity_vtx_n[idx+i_entity], idx+i_entity);

      int n_vtx_entity_1 = blk_entity_vtx_n  [idx+i_entity];
      int beg_1          = blk_entity_vtx_idx[idx+i_entity];
      int idx_next_same_entity = 0;
      sens_entity[idx_next_same_entity] = 1;
      same_entity_idx[idx_next_same_entity++] = i_entity;

      if(already_treat[i_entity] != 1) {

        PDM_g_num_t key_1 = 0;
        int idx_min_1 = -1;
        PDM_g_num_t min_1 = n_vtx_abs+1;
        for(int j = 0; j < n_vtx_entity_1; ++j) {
          loc_entity_vtx_1[j] = blk_tot_entity_vtx[beg_1+j];
          key_1 += loc_entity_vtx_1[j];
          if(loc_entity_vtx_1[j] < min_1) {
            min_1 = loc_entity_vtx_1[j];
            idx_min_1 = j;
          };
        }
        PDM_quick_sort_long(loc_entity_vtx_1, 0, n_vtx_entity_1-1);

        for(int i_entity_next = i_entity+1; i_entity_next < n_conflict_entitys; ++i_entity_next) {

          int n_vtx_entity_2 = blk_entity_vtx_n[idx+i_entity_next];

          if(n_vtx_entity_1 == n_vtx_entity_2 ) {

            int beg_2 = blk_entity_vtx_idx[idx+i_entity_next];
            PDM_g_num_t key_2 = 0;
            int idx_min_2 = -1;
            PDM_g_num_t min_2 = n_vtx_abs+1;
            for(int j = 0; j < n_vtx_entity_1; ++j) {
              loc_entity_vtx_2[j] = blk_tot_entity_vtx[beg_2+j];
              key_2 += loc_entity_vtx_2[j];
              if(loc_entity_vtx_2[j] < min_2) {
                min_2 = loc_entity_vtx_2[j];
                idx_min_2 = j;
              };
            }
            PDM_quick_sort_long(loc_entity_vtx_2, 0, n_vtx_entity_2-1);

            assert(key_1 == key_2);

            int is_same_entity = 1;
            for(int i_vtx = 0; i_vtx < n_vtx_entity_1; ++i_vtx) {
              if(loc_entity_vtx_1[i_vtx] != loc_entity_vtx_2[i_vtx]) {
                is_same_entity = -1;
                break;
              }
            }

            if(is_same_entity == 1 ){
              // printf("idx_min_1 = %i | idx_min_2 = %i \n", idx_min_1, idx_min_2);
              // printf(" test1 :: ");
              // for(int i_vtx = 0; i_vtx < n_vtx_entity_1; ++i_vtx) {
              //   printf(" %i", (int)blk_tot_entity_vtx[beg_1+i_vtx]);
              // }
              // printf(" \n");
              // printf(" test2 :: ");
              // for(int i_vtx = 0; i_vtx < n_vtx_entity_1; ++i_vtx) {
              //   printf(" %i", (int)blk_tot_entity_vtx[beg_2+i_vtx]);
              // }
              // printf(" \n");

              if(n_vtx_entity_1 == 2) {
                PDM_g_num_t i1 = blk_tot_entity_vtx[beg_1  ];
                PDM_g_num_t i2 = blk_tot_entity_vtx[beg_1+1];

                PDM_g_num_t j1 = blk_tot_entity_vtx[beg_2];
                PDM_g_num_t j2 = blk_tot_entity_vtx[beg_2+1];
                if(i1 == j1) {
                  sens_entity[idx_next_same_entity] = 1;
                } else {
                  assert(i1 == j2);
                  assert(i2 == j1);
                  sens_entity[idx_next_same_entity] = -1;
                }

              } else {
                // Determine the sens
                PDM_g_num_t i1 = blk_tot_entity_vtx[beg_1 +  idx_min_1                    ];
                PDM_g_num_t i2 = blk_tot_entity_vtx[beg_1 + (idx_min_1+1) % n_vtx_entity_1];

                PDM_g_num_t j1 = blk_tot_entity_vtx[beg_2 +  idx_min_2                    ];
                PDM_g_num_t j2 = blk_tot_entity_vtx[beg_2 + (idx_min_2+1) % n_vtx_entity_1];
                // printf(" i1 = %i | i2 = %i | j1 = %i | j2 = %i\n", i1, i2, j1, j2);
                assert(i1 == j1); // Panic
                if(i2 != j2) {
                  sens_entity[idx_next_same_entity] = -1;
                  // printf(" idx3 = %i \n", (idx_min_1+n_vtx_entity_1) % n_vtx_entity_1);
                  PDM_g_num_t i3 = blk_tot_entity_vtx[beg_1 + (idx_min_1+n_vtx_entity_1-1) % n_vtx_entity_1];
                  // printf(" i1 = %i | i2 = %i | i3 = %i | j1 = %i | j2 = %i\n", i1, i2, i3, j1, j2);
                  assert(i3 == j2);
                } else {
                  sens_entity[idx_next_same_entity] = 1;
                  assert(i2 == j2);
                }
              }

              // Check if same sens
              same_entity_idx[idx_next_same_entity++] = i_entity_next;
            }
          } /* End if same number of vertex */
        } /* End loop next entity */

         // printf("[%i] - same_entity_idx::", idx_next_same_entity);
         // for(int ii = 0; ii < idx_next_same_entity; ++ii) {
         //   printf(" %i", same_entity_idx[ii]);
         // }
         // printf("\n");

        _dentity_vtx_idx[i_abs_entity+1] = _dentity_vtx_idx[i_abs_entity] + n_vtx_entity_1;
        for(int i_vtx = 0; i_vtx < n_vtx_entity_1; ++i_vtx) {
          _dentity_vtx[idx_entity_vtx++] = blk_tot_entity_vtx[beg_1+i_vtx];
          // Ecriture à partir du minumun
          // _dentity_vtx[idx_entity_vtx++] = blk_tot_entity_vtx[beg_1+i_vtx];
        }

        _dentity_elmt_idx[i_abs_entity+1] = _dentity_elmt_idx[i_abs_entity];
        for(int i = 0; i < idx_next_same_entity; ++i) {
          int i_same_entity = same_entity_idx[i];
          int sign = sens_entity[i];
          // Signe à faire
          _dentity_elmt[_dentity_elmt_idx[i_abs_entity+1]++] = sign*blk_elmt_entity_elmt[idx+i_same_entity];
          already_treat[i_same_entity] = 1;
        }
        i_abs_entity++;

      } /* End loop already treated */

    } /* End loop entity in conflict */

    idx += n_conflict_entitys;

  }

  // printf(" realloc dentity_elmt : %i --> %i \n", n_tot_entity_per_key, _dentity_elmt_idx[i_abs_entity]);
  *dentity_elmt = realloc(*dentity_elmt, sizeof(PDM_g_num_t) *  _dentity_elmt_idx[i_abs_entity] );
  _dentity_elmt = *dentity_elmt;

  /*
   * Free all unused structure
   */
  free(loc_entity_vtx_1);
  free(loc_entity_vtx_2);
  free(blk_entity_vtx_idx);
  free(already_treat);
  free(same_entity_idx);
  free(sens_entity);
  free(delmt_entity);
  free(delmt_entity_vtx);
  free(ln_to_gn);
  free(blk_tot_entity_vtx);
  free(blk_n_entity_per_key);
  free(blk_entity_vtx_n);
  free(blk_elmt_entity_elmt);

  if( 0 == 1 ){
    printf("i_abs_entity::%i \n", i_abs_entity+1);
    PDM_log_trace_array_int(_dentity_vtx_idx, i_abs_entity+1                   , "_dentity_vtx_idx:: " );
    PDM_log_trace_array_long(_dentity_vtx   , _dentity_vtx_idx[i_abs_entity]   , "_dentity_vtx:: "     );
    PDM_log_trace_array_int(_dentity_elmt_idx, i_abs_entity+1                  , "_dentity_elmt_idx:: ");
    PDM_log_trace_array_long(_dentity_elmt    , _dentity_elmt_idx[i_abs_entity], "_dentity_elmt:: "    );
  }

  /*
   * Fill up structure
   */
  *dn_entity = i_abs_entity;
  int _dn_entity = *dn_entity;

  /*
   * Realloc
   */
  *dentity_vtx_idx = (int *        ) realloc(*dentity_vtx_idx, (_dn_entity + 1) * sizeof(int * ) );
  _dentity_vtx_idx = *dentity_vtx_idx;

  *dentity_vtx     = (PDM_g_num_t *) realloc(*dentity_vtx    , _dentity_vtx_idx[_dn_entity] * sizeof(PDM_g_num_t * ));
  _dentity_vtx     = *dentity_vtx;

  /*
   * Generate absolute numerotation of entitys
   */
  *entity_distrib = _make_absolute_entity_numbering(_dn_entity, comm);
}


static
void
_generate_faces_from_dmesh_nodal
(
  _pdm_link_dmesh_nodal_to_dmesh_t *link
)
{

  PDM_dmesh_nodal_t* dmesh_nodal = link->dmesh_nodal;

  int n_face_elt_tot     = 0;
  int n_sum_vtx_face_tot = 0;

  PDM_dmesh_nodal_decompose_faces_get_size(dmesh_nodal, &n_face_elt_tot, &n_sum_vtx_face_tot);

  PDM_g_num_t* delmt_face         = (PDM_g_num_t*) malloc(  n_face_elt_tot     * sizeof(PDM_g_num_t));
  int*         delmt_face_vtx_idx = (int        *) malloc( (n_face_elt_tot +1) * sizeof(int        ));
  PDM_g_num_t* delmt_face_vtx     = (PDM_g_num_t*) malloc(  n_sum_vtx_face_tot * sizeof(PDM_g_num_t));

  delmt_face_vtx_idx[0] = 0;
  PDM_sections_decompose_faces(dmesh_nodal,
                               delmt_face_vtx_idx,
                               delmt_face_vtx,
                               delmt_face,
                               NULL, NULL);

  /*
   *  Create empty dmesh
   */
  PDM_dmesh_t* dm = PDM_dmesh_create(link->owner, -1, -1, -1, -1, -1, -1, dmesh_nodal->pdm_mpi_comm);
  assert(link->dmesh == NULL);
  link->dmesh = dm;

  PDM_generate_entitiy_connectivity(dmesh_nodal->pdm_mpi_comm,
                                    dmesh_nodal->n_vtx_abs,
                                    n_face_elt_tot,
                                    delmt_face,
                                    delmt_face_vtx_idx,
                                    delmt_face_vtx,
                                    &dm->dn_face,
                                    &dm->face_distrib,
                                    &dm->_dface_vtx_idx,
                                    &dm->_dface_vtx,
                                    &link->_dface_elmt_idx,
                                    &link->_dface_elmt);

  dm->is_owner_connectivity[PDM_CONNECTIVITY_TYPE_FACE_VTX] = PDM_TRUE;
  dm->dconnectivity_idx    [PDM_CONNECTIVITY_TYPE_FACE_VTX] = dm->_dface_vtx_idx;
  dm->dconnectivity        [PDM_CONNECTIVITY_TYPE_FACE_VTX] = dm->_dface_vtx;

  int is_signed = 1;
  assert(link->elmt_distrib == NULL);
  link->elmt_distrib = (PDM_g_num_t * ) malloc( (dmesh_nodal->n_rank + 1 ) * sizeof(PDM_g_num_t));
  link->elmt_distrib[0] = -1;
  PDM_dconnectivity_transpose(dmesh_nodal->pdm_mpi_comm,
                              dm->face_distrib,
                              link->elmt_distrib,
                              link->_dface_elmt_idx,
                              link->_dface_elmt,
                              is_signed,
                              &link->_delmt_face_idx,
                              &link->_delmt_face);
  // PDM_log_trace_array_long(link->elmt_distrib  , dmesh_nodal->n_rank+1, "elmt_distrib:: ");
  // free(elmt_distrib);

  int dn_elmt = link->elmt_distrib[dmesh_nodal->i_rank+1] - link->elmt_distrib[dmesh_nodal->i_rank];
  link->dn_elmt = dn_elmt;

  if( 0 == 1 ){
    printf("dn_elmt ::%i\n", dn_elmt );
    PDM_log_trace_array_int (link->_delmt_face_idx, dn_elmt+1                     , "link->_delmt_face_idx:: ");
    PDM_log_trace_array_long(link->_delmt_face    , link->_delmt_face_idx[dn_elmt], "link->_delmt_face:: ");
  }

  // Count the number of edges
  int n_edge_elt_tot = dm->_dface_vtx_idx[dm->dn_face];

  int*         dface_edge_vtx_idx = (int         *) malloc( ( n_edge_elt_tot + 1) * sizeof(int        ) );
  PDM_g_num_t* dface_edge         = (PDM_g_num_t *) malloc(     n_edge_elt_tot    * sizeof(PDM_g_num_t) );
  PDM_g_num_t* dface_edge_vtx     = (PDM_g_num_t *) malloc( 2 * n_edge_elt_tot    * sizeof(PDM_g_num_t) );
  int idx = 0;
  int i_edge = 0;
  dface_edge_vtx_idx[0] = 0;
  for(int i_face = 0; i_face < dm->dn_face; ++i_face ) {
    int beg        = dm->_dface_vtx_idx[i_face  ];
    int n_vtx_elmt = dm->_dface_vtx_idx[i_face+1] - beg;

    for(int ivtx = 0; ivtx < n_vtx_elmt; ++ivtx ) {
      int inext = (ivtx + 1) % n_vtx_elmt;

      dface_edge_vtx[idx++] = dm->_dface_vtx[beg+ivtx ];
      dface_edge_vtx[idx++] = dm->_dface_vtx[beg+inext];

      dface_edge_vtx_idx[i_edge+1] = dface_edge_vtx_idx[i_edge] + 2;
      dface_edge[i_edge] = (PDM_g_num_t) i_face + dm->face_distrib[dmesh_nodal->i_rank] + 1;
      i_edge++;
    }
  }

  if( 0 == 1 ){
    printf("n_edge_elt_tot ::%i\n", n_edge_elt_tot );
    PDM_log_trace_array_int (dface_edge_vtx_idx, n_edge_elt_tot+1              , "dface_edge_vtx_idx:: ");
    PDM_log_trace_array_long(dface_edge_vtx    , dface_edge_vtx_idx[n_edge_elt_tot], "dface_edge_vtx:: ");
    PDM_log_trace_array_long(dface_edge        , n_edge_elt_tot                , "dface_edge:: ");
  }

  PDM_generate_entitiy_connectivity(dmesh_nodal->pdm_mpi_comm,
                                    dmesh_nodal->n_vtx_abs,
                                    n_edge_elt_tot,
                                    dface_edge,
                                    dface_edge_vtx_idx,
                                    dface_edge_vtx,
                                    &dm->dn_edge,
                                    &dm->edge_distrib,
                                    &dm->_dedge_vtx_idx,
                                    &dm->_dedge_vtx,
                                    &dm->_dedge_face_idx,
                                    &dm->_dedge_face);

  dm->is_owner_connectivity[PDM_CONNECTIVITY_TYPE_EDGE_VTX ] = PDM_TRUE;
  dm->is_owner_connectivity[PDM_CONNECTIVITY_TYPE_EDGE_FACE] = PDM_TRUE;
  dm->dconnectivity_idx[PDM_CONNECTIVITY_TYPE_EDGE_VTX]  = dm->_dedge_vtx_idx;
  dm->dconnectivity    [PDM_CONNECTIVITY_TYPE_EDGE_VTX]  = dm->_dedge_vtx;
  dm->dconnectivity_idx[PDM_CONNECTIVITY_TYPE_EDGE_FACE] = dm->_dedge_face_idx;
  dm->dconnectivity    [PDM_CONNECTIVITY_TYPE_EDGE_FACE] = dm->_dedge_face;

  // &dmesh_nodal->dface_edge_idx,
  // &dmesh_nodal->dface_edge,
  if( 0 == 1 ){
    printf("dmesh_nodal->dn_edge ::%i\n", dm->dn_edge );
    PDM_log_trace_array_int (dm->_dedge_vtx_idx, dm->dn_edge+1                  , "dm->_dedge_vtx_idx:: ");
    PDM_log_trace_array_long(dm->_dedge_vtx    , dm->_dedge_vtx_idx[dm->dn_edge], "dm->_dedge_vtx:: ");

    // PDM_log_trace_array_int (dm->dface_edge_idx, dm->dn_face+1                           , "dm->dface_edge_idx:: ");
    // PDM_log_trace_array_long(dm->dface_edge    , dm->dface_edge_idx[dm->dn_face], "dm->dface_edge:: ");

    PDM_log_trace_array_int (dm->_dedge_face_idx, dm->dn_edge+1                   , "dm->_dedge_face_idx:: ");
    PDM_log_trace_array_long(dm->_dedge_face    , dm->_dedge_face_idx[dm->dn_edge], "dm->_dedge_face:: ");

  }

}



static
void
_generate_edges_from_dmesh_nodal
(
  _pdm_link_dmesh_nodal_to_dmesh_t *link
)
{
  PDM_dmesh_nodal_t* dmesh_nodal = link->dmesh_nodal;

  int n_edge_elt_tot     = 0;
  int n_sum_vtx_edge_tot = 0;

  PDM_dmesh_nodal_decompose_edges_get_size(dmesh_nodal, &n_edge_elt_tot, &n_sum_vtx_edge_tot);

  PDM_g_num_t* delmt_edge         = (PDM_g_num_t*) malloc(  n_edge_elt_tot     * sizeof(PDM_g_num_t));
  int*         delmt_edge_vtx_idx = (int        *) malloc( (n_edge_elt_tot +1) * sizeof(int        ));
  PDM_g_num_t* delmt_edge_vtx     = (PDM_g_num_t*) malloc(  n_sum_vtx_edge_tot * sizeof(PDM_g_num_t));

  delmt_edge_vtx_idx[0] = 0;
  PDM_sections_decompose_edges(dmesh_nodal,
                               delmt_edge_vtx_idx,
                               delmt_edge_vtx,
                               delmt_edge,
                               NULL, NULL);

  /*
   *  Create empty dmesh
   */
  PDM_dmesh_t* dm = PDM_dmesh_create(link->owner, -1, -1, -1, -1, -1, -1, dmesh_nodal->pdm_mpi_comm);
  assert(link->dmesh == NULL);
  link->dmesh = dm;

  PDM_generate_entitiy_connectivity(dmesh_nodal->pdm_mpi_comm,
                                    dmesh_nodal->n_vtx_abs,
                                    n_edge_elt_tot,
                                    delmt_edge,
                                    delmt_edge_vtx_idx,
                                    delmt_edge_vtx,
                                    &dm->dn_edge,
                                    &dm->edge_distrib,
                                    &dm->_dedge_vtx_idx,
                                    &dm->_dedge_vtx,
                                    &link->_dedge_elmt_idx,
                                    &link->_dedge_elmt);

  dm->is_owner_connectivity[PDM_CONNECTIVITY_TYPE_EDGE_VTX] = PDM_TRUE;
  dm->dconnectivity_idx    [PDM_CONNECTIVITY_TYPE_EDGE_VTX] = dm->_dedge_vtx_idx;
  dm->dconnectivity        [PDM_CONNECTIVITY_TYPE_EDGE_VTX] = dm->_dedge_vtx;

  int is_signed = 1;
  assert(link->elmt_distrib == NULL);
  link->elmt_distrib = (PDM_g_num_t * ) malloc( (dmesh_nodal->n_rank + 1 ) * sizeof(PDM_g_num_t));
  link->elmt_distrib[0] = -1;
  PDM_dconnectivity_transpose(dmesh_nodal->pdm_mpi_comm,
                              dm->edge_distrib,
                              link->elmt_distrib,
                              link->_dedge_elmt_idx,
                              link->_dedge_elmt,
                              is_signed,
                              &link->_delmt_edge_idx,
                              &link->_delmt_edge);

  int dn_elmt = link->elmt_distrib[dmesh_nodal->i_rank+1] - link->elmt_distrib[dmesh_nodal->i_rank];
  link->dn_elmt = dn_elmt;

  if( 0 == 1 ){
    printf("dn_elmt ::%i\n", dn_elmt );
    PDM_log_trace_array_int (link->_delmt_edge_idx, dn_elmt+1                     , "link->_delmt_edge_idx:: ");
    PDM_log_trace_array_long(link->_delmt_edge    , link->_delmt_edge_idx[dn_elmt], "link->_delmt_edge:: ");
  }

}


static
void
_translate_element_group_to_entity
(
 PDM_MPI_Comm  comm,
 PDM_g_num_t  *entity_distrib,
 PDM_g_num_t  *dgroup_elmt,
 int          *dgroup_elmt_idx,
 int           n_group_elmt,
 PDM_g_num_t  *delmt_entity,
 int          *delmt_entity_idx,
 PDM_g_num_t **dentity_bound,
 int         **dentity_bound_idx
)
{
  int i_rank;
  PDM_MPI_Comm_rank (comm, &i_rank);
  int n_rank;
  PDM_MPI_Comm_size (comm, &n_rank);

  // if(1 == 1) {
  //   PDM_log_trace_array_int (dgroup_elmt_idx, n_group_elmt+1         , "dgroup_elmt_idx:: ");
  //   PDM_log_trace_array_long(dgroup_elmt    , dgroup_elmt_idx[n_group_elmt], "dgroup_elmt:: ");
  // }
  /*
   * Prepare exchange protocol
   */
  PDM_block_to_part_t* btp = PDM_block_to_part_create(entity_distrib,
                               (const PDM_g_num_t **) &dgroup_elmt,
                                                      &dgroup_elmt_idx[n_group_elmt],
                                                      1,
                                                      comm);

  /*
   * Exchange
   */
  int dn_entity = entity_distrib[i_rank+1] - entity_distrib[i_rank];
  int* delmt_entity_n = (int *) malloc( dn_entity * sizeof(int));
  for(int i = 0; i < dn_entity; ++i) {
    delmt_entity_n[i] = delmt_entity_idx[i+1] - delmt_entity_idx[i];
  }

  if(0 == 1) {
    PDM_log_trace_array_long (entity_distrib, n_rank+1, "entity_distrib:: ");
    PDM_log_trace_array_int (delmt_entity_n, dn_entity, "delmt_entity_n:: ");
    PDM_log_trace_array_long(delmt_entity  , delmt_entity_idx[dn_entity], "delmt_entity:: ");
  }

  int**         part_group_stri;
  PDM_g_num_t** part_group_data;
  PDM_block_to_part_exch2(btp,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_VAR,
                          delmt_entity_n,
             (void *  )   delmt_entity,
             (int  ***)  &part_group_stri,
             (void ***)  &part_group_data);
  free(delmt_entity_n);

  int*         _part_group_stri = part_group_stri[0];
  PDM_g_num_t* _part_group_data = part_group_data[0];

  *dentity_bound_idx = (int * ) malloc( (n_group_elmt+1) * sizeof(int) );
  int* _dentity_bound_idx = *dentity_bound_idx;

  int idx_data = 0;
  int idx_stri = 0;
  _dentity_bound_idx[0] = 0;
  for(int i_group = 0; i_group < n_group_elmt; ++i_group) {
    _dentity_bound_idx[i_group+1] = _dentity_bound_idx[i_group];
    for(int ielmt = dgroup_elmt_idx[i_group]; ielmt < dgroup_elmt_idx[i_group+1]; ++ielmt) {
      _dentity_bound_idx[i_group+1] += _part_group_stri[idx_stri++];
      _part_group_data[idx_data] = PDM_ABS(_part_group_data[idx_data]);
      idx_data++;
    }
  }
  *dentity_bound = _part_group_data;

  free(part_group_stri);
  free(_part_group_stri);
  free(part_group_data);

  if(0 == 1) {
    PDM_log_trace_array_int (_dentity_bound_idx, n_group_elmt+1                  , "_dentity_bound_idx:: ");
    PDM_log_trace_array_long(_part_group_data  , _dentity_bound_idx[n_group_elmt], "_dentity_bound:: ");
  }
  PDM_block_to_part_free(btp);
}

static
void
_translate_element_group_to_faces
(
 _pdm_link_dmesh_nodal_to_dmesh_t* link
)
{
  PDM_dmesh_nodal_t *dmesh_nodal = link->dmesh_nodal;
  PDM_dmesh_t       *dm          = link->dmesh;

  // printf("_translate_element_group_to_faces \n");
  assert(dmesh_nodal != NULL);
  assert(dm          != NULL);

  PDM_g_num_t *dface_bound;
  int         *dface_bound_idx;

  _translate_element_group_to_entity(dmesh_nodal->pdm_mpi_comm,
                                     link->elmt_distrib,
                                     dmesh_nodal->dgroup_elmt,
                                     dmesh_nodal->dgroup_elmt_idx,
                                     dmesh_nodal->n_group_elmt,
                                     link->_delmt_face,
                                     link->_delmt_face_idx,
                                     &dface_bound,
                                     &dface_bound_idx);

  dm->is_owner_bound[PDM_BOUND_TYPE_FACE] = PDM_TRUE;
  dm->dbound_idx    [PDM_BOUND_TYPE_FACE] = dface_bound_idx;
  dm->dbound        [PDM_BOUND_TYPE_FACE] = dface_bound;
  dm->n_bnd                               = dmesh_nodal->n_group_elmt;

  // Par recursion on peut avoir les group de vertex ou de edge

}

static
void
_translate_element_group_to_edges
(
 _pdm_link_dmesh_nodal_to_dmesh_t* link
)
{
  PDM_dmesh_nodal_t *dmesh_nodal = link->dmesh_nodal;
  PDM_dmesh_t       *dm          = link->dmesh;

  PDM_g_num_t *dedge_bound;
  int         *dedge_bound_idx;

  _translate_element_group_to_entity(dmesh_nodal->pdm_mpi_comm,
                                     link->elmt_distrib,
                                     dmesh_nodal->dgroup_elmt,
                                     dmesh_nodal->dgroup_elmt_idx,
                                     dmesh_nodal->n_group_elmt,
                                     link->_delmt_edge,
                                     link->_delmt_edge_idx,
                                     &dedge_bound,
                                     &dedge_bound_idx);

  dm->is_owner_bound[PDM_BOUND_TYPE_EDGE] = PDM_TRUE;
  dm->dbound_idx    [PDM_BOUND_TYPE_EDGE] = dedge_bound_idx;
  dm->dbound        [PDM_BOUND_TYPE_EDGE] = dedge_bound;
  dm->n_bnd                               = dmesh_nodal->n_group_elmt;

  // Par recursion on peut avoir les group de vertex

}

static
void
_to_coherent_2d
(
 PDM_MPI_Comm                      comm,
 _pdm_link_dmesh_nodal_to_dmesh_t *link,
 PDM_dmesh_nodal_t                *dmesh_nodal,
 PDM_dmesh_t                      *dmesh
)
{
  int n_rank;
  int i_rank;

  PDM_MPI_Comm_size (comm, &n_rank);
  PDM_MPI_Comm_rank (comm, &i_rank);

  PDM_g_num_t *section_distribution = dmesh_nodal->section_distribution;

  // 1. delmt_edge --> dface_edge (because 2D)
  PDM_g_num_t *dface_edge     = (PDM_g_num_t *) malloc( link->_delmt_edge_idx[link->dn_elmt] * sizeof(PDM_g_num_t));
  int         *dface_edge_idx = (int         *) malloc( (link->dn_elmt + 1 )                 * sizeof(int        ));

  if( 0 == 1 ){
    PDM_log_trace_array_long(link->elmt_distrib , n_rank+1, "link->elmt_distrib");
    PDM_log_trace_array_long(dmesh->edge_distrib, n_rank+1, "link->edge_distrib");
    PDM_log_trace_array_long(section_distribution, dmesh_nodal->n_section_tot+1, "section_distribution");
  }

  // Find out in dmesh_nodal where the last element of current proc
  PDM_g_num_t first_elmt = link->elmt_distrib[i_rank  ];
  PDM_g_num_t last_elmt  = link->elmt_distrib[i_rank+1]-1;
  int first_section = PDM_binary_search_gap_long(first_elmt, section_distribution, dmesh_nodal->n_section_tot+1);
  int last_section  = PDM_binary_search_gap_long(last_elmt , section_distribution, dmesh_nodal->n_section_tot+1);

  printf("[%i] first_section = %i \n", i_rank, first_section);
  printf("[%i] last_section  = %i \n", i_rank, last_section );

  int idx   = 0;
  int dn_face = 0;
  dface_edge_idx[0] = 0;
  PDM_g_num_t cur_elmt = link->elmt_distrib[dmesh_nodal->i_rank];
  for(int i_section = first_section; i_section <= last_section; ++i_section) {

    int id_section = dmesh_nodal->sections_id[i_section];
    PDM_Mesh_nodal_elt_t t_elt   = PDM_DMesh_nodal_section_type_get   (dmesh_nodal, id_section);
    // const PDM_g_num_t*   distrib = PDM_DMesh_nodal_distrib_section_get(dmesh_nodal, id_section);

    if( _is_a_2d_elment(t_elt) ){

      printf(" filter face \n");

      // log_trace("section_distribution[%i] = %i \n", i_section, section_distribution[i_section  ]);
      // log_trace("section_distribution[%i] = %i \n", i_section+1, section_distribution[i_section+1]);
      // log_trace("elmt_distrib[%i] = %i \n", dmesh_nodal->i_rank, link->elmt_distrib[dmesh_nodal->i_rank]);
      // log_trace("elmt_distrib[%i] = %i \n", dmesh_nodal->i_rank+1,  link->elmt_distrib[dmesh_nodal->i_rank+1]);

      int beg_elmt = -1;
      int end_elmt = -1;

      if( section_distribution[i_section] < link->elmt_distrib[dmesh_nodal->i_rank]) {
        // log_trace("[begin] cas 1 : section_distribution[i_section] < link->elmt_distrib[dmesh_nodal->i_rank] \n");
        beg_elmt = 0; // On est dans la section mais on commence
      } else if( section_distribution[i_section] == cur_elmt ) {
        // log_trace("[begin] cas 2 : section_distribution[i_section] == cur_elmt \n");
        beg_elmt = dn_face; // On continue car c'est la suite
      }

      if( section_distribution[i_section+1] > link->elmt_distrib[dmesh_nodal->i_rank+1]) {
        // log_trace("[end] cas 1 : section_distribution[i_section+1] > link->elmt_distrib[dmesh_nodal->i_rank+1] \n");
        end_elmt = link->elmt_distrib[dmesh_nodal->i_rank+1] - link->elmt_distrib[dmesh_nodal->i_rank]; // Si section trop grande on s'arrete à cell du proca
      } else {
        // log_trace("[end] cas 2 : %i ->  %i \n",section_distribution[i_section+1], link->elmt_distrib[dmesh_nodal->i_rank] );
        end_elmt = section_distribution[i_section+1] - link->elmt_distrib[dmesh_nodal->i_rank];
      }

      // printf("[%i] - beg_elmt = %i | end_elmt = %i \n", i_rank, beg_elmt, end_elmt);
      // log_trace(" ---------------------  \n");
      // log_trace("beg_elmt::%i | end_elmt::%i \n", beg_elmt, end_elmt);

      for( int ielmt = beg_elmt; ielmt < end_elmt; ++ielmt ) {
        // dcell_elmt[dn_face] = ielmt;
        // log_trace("ielmt::%i --> ", ielmt);
        // log_trace("idxedge = %i -> %i \n ", link->_delmt_edge_idx[ielmt], link->_delmt_edge_idx[ielmt+1]);
        dface_edge_idx[dn_face+1] = dface_edge_idx[dn_face];
        for(int iedge = link->_delmt_edge_idx[ielmt]; iedge < link->_delmt_edge_idx[ielmt+1]; ++iedge ){
          dface_edge[idx++] = link->_delmt_edge[iedge];
          dface_edge_idx[dn_face+1]++;
        }
        dn_face  += 1;
        cur_elmt += 1;
      }
    }
  }

  // printf(" dn_face = %i\n", dn_face);
  dface_edge_idx = (int         *) realloc(dface_edge_idx, (dn_face+1)             * sizeof(int        ));
  dface_edge     = (PDM_g_num_t *) realloc(dface_edge    , dface_edge_idx[dn_face] * sizeof(PDM_g_num_t));

  if(0 == 1) {
    PDM_log_trace_array_int (dface_edge_idx, dn_face+1              , "dface_edge_idx");
    PDM_log_trace_array_long(dface_edge    , dface_edge_idx[dn_face], "dface_edge");
    log_trace("dn_face::%i --> %p", dn_face, dface_edge);
  }

  assert(dmesh->dn_face == -1);
  dmesh->dn_face = dn_face;

  dmesh->is_owner_connectivity[PDM_CONNECTIVITY_TYPE_FACE_EDGE] = PDM_TRUE;
  dmesh->dconnectivity_idx    [PDM_CONNECTIVITY_TYPE_FACE_EDGE] = dface_edge_idx;
  dmesh->dconnectivity        [PDM_CONNECTIVITY_TYPE_FACE_EDGE] = dface_edge;

  // 2. dedge_elmt --> dedge_face
  // Caution : if edge_face is in the bad sens we need to adapt edge_vtx (see partitioning with flip)
  int         *edge_face_idx;
  PDM_g_num_t *edge_face_tmp;
  dmesh->face_distrib = PDM_compute_entity_distribution(comm, dn_face); // Begin a 1

  if( 0 == 1 ){
    PDM_log_trace_array_long(dmesh->face_distrib, n_rank+1, "dmesh->face_distrib::");
    PDM_log_trace_array_long(dmesh->edge_distrib, n_rank+1, "dmesh->edge_distrib::");
  }

  PDM_dconnectivity_transpose(dmesh_nodal->pdm_mpi_comm,
                              dmesh->face_distrib,
                              dmesh->edge_distrib,
                              dface_edge_idx,
                              dface_edge,
                              1, // is_signed
                              &edge_face_idx,
                              &edge_face_tmp);
  int dn_edge = dmesh->edge_distrib[i_rank+1] - dmesh->edge_distrib[i_rank];
  assert(dn_edge == dmesh->dn_edge);

  if(0 == 1) {
    PDM_log_trace_array_int (edge_face_idx, dn_edge+1             , "edge_face_idx::");
    PDM_log_trace_array_long(edge_face_tmp, edge_face_idx[dn_edge], "edge_face_tmp::");
  }

  // Post_treat
  PDM_g_num_t *edge_face = (PDM_g_num_t *) malloc( 2 * dn_edge * sizeof(PDM_g_num_t));;
  for(int i_edge = 0; i_edge < dn_edge; ++i_edge) {

    int beg = edge_face_idx[i_edge];
    int n_connect_face = edge_face_idx[i_edge+1] - beg;
    if(n_connect_face == 1) {
      edge_face[2*i_edge  ] = PDM_ABS(edge_face_tmp[beg]); // Attention on peut être retourner !!!!
      edge_face[2*i_edge+1] = 0;
    } else {
      assert(n_connect_face == 2);
      edge_face[2*i_edge  ] = PDM_ABS(edge_face_tmp[beg  ]);
      edge_face[2*i_edge+1] = PDM_ABS(edge_face_tmp[beg+1]);
    }

    // Flip if the edge is in the other sens
    if( PDM_SIGN(edge_face_tmp[beg]) == -1 ) {
      int beg_edge_vtx = dmesh->_dedge_vtx_idx[i_edge  ];
      int end_edge_vtx = dmesh->_dedge_vtx_idx[i_edge+1];

      PDM_g_num_t tmp_swap;
      tmp_swap = dmesh->_dedge_vtx[beg_edge_vtx];
      dmesh->_dedge_vtx[beg_edge_vtx  ] = dmesh->_dedge_vtx[end_edge_vtx-1];
      dmesh->_dedge_vtx[end_edge_vtx-1] = tmp_swap;

    }
  }

  free(edge_face_idx);
  free(edge_face_tmp);

  dmesh->is_owner_connectivity[PDM_CONNECTIVITY_TYPE_EDGE_FACE] = PDM_TRUE;
  dmesh->dconnectivity_idx    [PDM_CONNECTIVITY_TYPE_EDGE_FACE] = NULL;
  dmesh->dconnectivity        [PDM_CONNECTIVITY_TYPE_EDGE_FACE] = edge_face;

  if(0 == 1) {
    PDM_log_trace_array_long(edge_face, 2 * dn_edge, "edge_face::");
  }

}

static
void
_to_coherent_3d
(
 PDM_MPI_Comm                      comm,
 _pdm_link_dmesh_nodal_to_dmesh_t *link,
 PDM_dmesh_nodal_t                *dmesh_nodal,
 PDM_dmesh_t                      *dmesh
)
{
  int n_rank;
  int i_rank;

  PDM_MPI_Comm_size (comm, &n_rank);
  PDM_MPI_Comm_rank (comm, &i_rank);

  PDM_g_num_t *section_distribution = dmesh_nodal->section_distribution;

  // 1. delmt_face --> dcell_face
  PDM_g_num_t *dcell_face     = (PDM_g_num_t *) malloc( link->_delmt_face_idx[link->dn_elmt] * sizeof(PDM_g_num_t));
  int         *dcell_face_idx = (int         *) malloc( (link->dn_elmt + 1 )                 * sizeof(int        ));

  if( 0 == 1 ){
    PDM_log_trace_array_long(link->elmt_distrib , n_rank+1, "link->elmt_distrib");
    PDM_log_trace_array_long(dmesh->face_distrib, n_rank+1, "link->face_distrib");
    PDM_log_trace_array_long(section_distribution, dmesh_nodal->n_section_tot+1, "section_distribution");
  }

  // Find out in dmesh_nodal where the last element of current proc
  PDM_g_num_t first_elmt = link->elmt_distrib[i_rank  ];
  PDM_g_num_t last_elmt  = link->elmt_distrib[i_rank+1]-1;
  int first_section = PDM_binary_search_gap_long(first_elmt, section_distribution, dmesh_nodal->n_section_tot+1);
  int last_section  = PDM_binary_search_gap_long(last_elmt , section_distribution, dmesh_nodal->n_section_tot+1);

  // printf("first_section = %i \n", first_section);
  // printf("last_section  = %i \n", last_section);

  int idx     = 0;
  int dn_cell = 0;
  dcell_face_idx[0] = 0;
  PDM_g_num_t cur_elmt = link->elmt_distrib[dmesh_nodal->i_rank];
  for(int i_section = first_section; i_section <= last_section; ++i_section) {

    int id_section = dmesh_nodal->sections_id[i_section];
    PDM_Mesh_nodal_elt_t t_elt   = PDM_DMesh_nodal_section_type_get   (dmesh_nodal, id_section);
    // const PDM_g_num_t*   distrib = PDM_DMesh_nodal_distrib_section_get(dmesh_nodal, id_section);

    if( _is_a_3d_elment(t_elt) ){

      printf(" filter cell \n");

      // log_trace("section_distribution[%i] = %i \n", i_section, section_distribution[i_section  ]);
      // log_trace("section_distribution[%i] = %i \n", i_section+1, section_distribution[i_section+1]);
      // log_trace("elmt_distrib[%i] = %i \n", dmesh_nodal->i_rank, link->elmt_distrib[dmesh_nodal->i_rank]);
      // log_trace("elmt_distrib[%i] = %i \n", dmesh_nodal->i_rank+1,  link->elmt_distrib[dmesh_nodal->i_rank+1]);

      int beg_elmt = -1;
      int end_elmt = -1;

      if( section_distribution[i_section] < link->elmt_distrib[dmesh_nodal->i_rank]) {
        // log_trace("[begin] cas 1 : section_distribution[i_section] < link->elmt_distrib[dmesh_nodal->i_rank] \n");
        beg_elmt = 0; // On est dans la section mais on commence
      } else if( section_distribution[i_section] == cur_elmt ) {
        // log_trace("[begin] cas 2 : section_distribution[i_section] == cur_elmt \n");
        beg_elmt = dn_cell; // On continue car c'est la suite
      }

      if( section_distribution[i_section+1] > link->elmt_distrib[dmesh_nodal->i_rank+1]) {
        // log_trace("[end] cas 1 : section_distribution[i_section+1] > link->elmt_distrib[dmesh_nodal->i_rank+1] \n");
        end_elmt = link->elmt_distrib[dmesh_nodal->i_rank+1] - link->elmt_distrib[dmesh_nodal->i_rank]; // Si section trop grande on s'arrete à cell du proc
      } else {
        // log_trace("[end] cas 2 : %i ->  %i \n",section_distribution[i_section+1], link->elmt_distrib[dmesh_nodal->i_rank] );
        end_elmt = section_distribution[i_section+1] - link->elmt_distrib[dmesh_nodal->i_rank];
      }

      // int beg_elmt = PDM_MAX(section_distribution[i_section  ] - link->elmt_distrib[dmesh_nodal->i_rank], 0);
      // int beg_elmt = section_distribution[i_section]+link->elmt_distrib[dmesh_nodal->i_rank]-cur_elmt;
      // int beg_elmt = section_distribution[i_section]-cur_elmt;
      // int end_elmt = PDM_MIN(section_distribution[i_section+1], link->elmt_distrib[dmesh_nodal->i_rank+1]) - cur_elmt;

      // log_trace("PDM_MIN(section_distribution[i_section  ], link->elmt_distrib[dmesh_nodal->i_rank  ]) :: %i \n", PDM_MIN(section_distribution[i_section ], link->elmt_distrib[dmesh_nodal->i_rank]));
      // log_trace("PDM_MIN(section_distribution[i_section+1], link->elmt_distrib[dmesh_nodal->i_rank+1]) :: %i \n", PDM_MIN(section_distribution[i_section+1], link->elmt_distrib[dmesh_nodal->i_rank+1]));
      // log_trace("cur_elmt :: %i \n", cur_elmt);
      // log_trace("[%i] - beg_elmt = %i | end_elmt = %i \n", i_rank, beg_elmt, end_elmt);
      // printf("[%i] - beg_elmt = %i | end_elmt = %i \n", i_rank, beg_elmt, end_elmt);

      for( int ielmt = beg_elmt; ielmt < end_elmt; ++ielmt ) {
        // log_trace("ielmt::%i \n", ielmt);
        // dcell_elmt[dn_cell] = (PDM_g_num_t) section_distribution[i_section] + ielmt;
        dcell_face_idx[dn_cell+1] = dcell_face_idx[dn_cell];
        for(int iface = link->_delmt_face_idx[ielmt]; iface < link->_delmt_face_idx[ielmt+1]; ++iface ){
          dcell_face[idx++] = link->_delmt_face[iface];
          dcell_face_idx[dn_cell+1]++;
        }
        dn_cell += 1;
        cur_elmt += 1;
      }
    }
  }

  // printf(" dn_cell = %i\n", dn_cell);
  dcell_face_idx = (int         *) realloc(dcell_face_idx, (dn_cell+1)             * sizeof(int        ));
  dcell_face     = (PDM_g_num_t *) realloc(dcell_face    , dcell_face_idx[dn_cell] * sizeof(PDM_g_num_t));

  if(0 == 1) {
    PDM_log_trace_array_int (dcell_face_idx, dn_cell+1              , "dcell_face_idx");
    PDM_log_trace_array_long(dcell_face    , dcell_face_idx[dn_cell], "dcell_face");
  }

  assert(dmesh->dn_cell == -1);
  dmesh->dn_cell = dn_cell;

  dmesh->is_owner_connectivity[PDM_CONNECTIVITY_TYPE_CELL_FACE] = PDM_TRUE;
  dmesh->dconnectivity_idx    [PDM_CONNECTIVITY_TYPE_CELL_FACE] = dcell_face_idx;
  dmesh->dconnectivity        [PDM_CONNECTIVITY_TYPE_CELL_FACE] = dcell_face;

  // 2. dface_elmt --> dface_cell
  // Caution : if face_cell is in the bad sens we need to adapt face_vtx (see partitioning with flip)
  int         *face_cell_idx;
  PDM_g_num_t *face_cell_tmp;
  dmesh->cell_distrib = PDM_compute_entity_distribution(comm, dn_cell); // Begin a 1
  // PDM_log_trace_array_long(dmesh->cell_distrib, n_rank+1, "dmesh->cell_distrib::");
  // PDM_log_trace_array_long(dmesh->face_distrib, n_rank+1, "dmesh->face_distrib::");
  PDM_dconnectivity_transpose(dmesh_nodal->pdm_mpi_comm,
                              dmesh->cell_distrib,
                              dmesh->face_distrib,
                              dcell_face_idx,
                              dcell_face,
                              1, // is_signed
                              &face_cell_idx,
                              &face_cell_tmp);
  int dn_face = dmesh->face_distrib[i_rank+1] - dmesh->face_distrib[i_rank];
  assert(dn_face == dmesh->dn_face);
  if(0 == 1) {
    PDM_log_trace_array_int (face_cell_idx, dn_face+1             , "face_cell_idx::");
    PDM_log_trace_array_long(face_cell_tmp, face_cell_idx[dn_face], "face_cell_tmp::");
  }

  // log_trace("face_distrib[i_rank  ] = %i \n", dmesh->face_distrib[i_rank]);
  // log_trace("face_distrib[i_rank+1] = %i \n", dmesh->face_distrib[i_rank+1]);
  // printf("dcell_face_idx[dn_cell] = %i \n", dcell_face_idx[dmesh->dn_cell]);
  // printf("face_cell_idx[dn_face] = %i \n", face_cell_idx[dmesh->dn_face]);
  // log_trace("dcell_face_idx[dn_cell] = %i \n", dcell_face_idx[dmesh->dn_cell]);
  // log_trace("face_cell_idx[dn_face] = %i \n", face_cell_idx[dmesh->dn_face]);
  // log_trace("face_distrib[i_rank+1] = %i \n", dmesh->face_distrib[i_rank+1]);
  // log_trace("dn_face = %i \n", dn_face);

  // Post_treat
  PDM_g_num_t *face_cell = (PDM_g_num_t *) malloc( 2 * dn_face * sizeof(PDM_g_num_t));;
  for(int i_face = 0; i_face < dn_face; ++i_face) {

    int beg = face_cell_idx[i_face];
    int n_connect_cell = face_cell_idx[i_face+1] - beg;
    if(n_connect_cell == 1) {
      face_cell[2*i_face  ] = PDM_ABS(face_cell_tmp[beg]); // Attention on peut être retourner !!!!
      face_cell[2*i_face+1] = 0;
    } else {
      assert(n_connect_cell == 2);
      face_cell[2*i_face  ] = PDM_ABS(face_cell_tmp[beg  ]);
      face_cell[2*i_face+1] = PDM_ABS(face_cell_tmp[beg+1]);
    }

    // Flip if the face is in the other sens
    if( PDM_SIGN(face_cell_tmp[beg]) == -1 ) {
      int beg_face_vtx = dmesh->_dface_vtx_idx[i_face  ];
      int end_face_vtx = dmesh->_dface_vtx_idx[i_face+1];

      int offset = beg_face_vtx + 1;
      int n_vtx_on_face = end_face_vtx - beg_face_vtx - 1;
      int end_index = n_vtx_on_face - 1;
      PDM_g_num_t tmp_swap;
      for(int i_vtx = 0; i_vtx < n_vtx_on_face/2; ++i_vtx) {
        tmp_swap = dmesh->_dface_vtx[offset+end_index];
        dmesh->_dface_vtx[offset+end_index] = dmesh->_dface_vtx[offset+i_vtx];
        dmesh->_dface_vtx[offset+i_vtx] = tmp_swap;
        end_index--;

      }

    }
  }

  free(face_cell_tmp);
  free(face_cell_idx);

  dmesh->is_owner_connectivity[PDM_CONNECTIVITY_TYPE_FACE_CELL] = PDM_TRUE;
  dmesh->dconnectivity_idx    [PDM_CONNECTIVITY_TYPE_FACE_CELL] = NULL;
  dmesh->dconnectivity        [PDM_CONNECTIVITY_TYPE_FACE_CELL] = face_cell;

  if(0 == 1) {
    PDM_log_trace_array_long(face_cell, 2 * dn_face, "face_cell::");
  }


}


static
_pdm_link_dmesh_nodal_to_dmesh_t*
_link_dmesh_nodal_to_dmesh_init
(
 PDM_ownership_t owner
)
{
  _pdm_link_dmesh_nodal_to_dmesh_t* link = malloc( sizeof(_pdm_link_dmesh_nodal_to_dmesh_t));

  link->dmesh_nodal       = NULL;
  link->dmesh             = NULL;

  link->owner             = owner;
  link->results_is_getted = PDM_FALSE;

  link->dn_elmt           = 0;
  link->elmt_distrib      = NULL;

  link->_delmt_face       = NULL;
  link->_delmt_face_idx   = NULL;

  link->_dface_elmt       = NULL;
  link->_dface_elmt_idx   = NULL;

  link->_delmt_edge       = NULL;
  link->_delmt_edge_idx   = NULL;

  link->_dedge_elmt       = NULL;
  link->_dedge_elmt_idx   = NULL;

  return link;
}


static
void
_link_dmesh_nodal_to_dmesh_free
(
 _pdm_link_dmesh_nodal_to_dmesh_t* link
)
{
  link->dmesh_nodal      = NULL; /* On a pas l'onwership de cette structure */

  if(( link->owner == PDM_OWNERSHIP_KEEP ) ||
     ( link->owner == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE && !link->results_is_getted)){
    PDM_dmesh_free(link->dmesh);
  }

  if(link->elmt_distrib != NULL) {
    free(link->elmt_distrib);
    link->elmt_distrib = NULL;
  }

  if(link->_delmt_face != NULL) {
    free(link->_delmt_face);
    link->_delmt_face = NULL;
  }

  if(link->_delmt_face_idx != NULL) {
    free(link->_delmt_face_idx);
    link->_delmt_face_idx = NULL;
  }

  if(link->_dface_elmt != NULL) {
    free(link->_dface_elmt);
    link->_dface_elmt = NULL;
  }

  if(link->_dface_elmt_idx != NULL) {
    free(link->_dface_elmt_idx);
    link->_dface_elmt_idx = NULL;
  }

  if(link->_dedge_elmt != NULL) {
    free(link->_dedge_elmt);
    link->_dedge_elmt = NULL;
  }

  if(link->_dedge_elmt_idx != NULL) {
    free(link->_dedge_elmt_idx);
    link->_dedge_elmt_idx = NULL;
  }

  free(link);

}

/*=============================================================================
 * Public function definitions
 *============================================================================*/

PDM_dmesh_nodal_to_dmesh_t*
PDM_dmesh_nodal_to_dmesh_create
(
const int             n_mesh,
const PDM_MPI_Comm    comm,
const PDM_ownership_t owner
)
{
  PDM_dmesh_nodal_to_dmesh_t *dmesh_nodal_to_dm = (PDM_dmesh_nodal_to_dmesh_t *) malloc(sizeof(PDM_dmesh_nodal_to_dmesh_t));

  dmesh_nodal_to_dm->comm              = comm;
  dmesh_nodal_to_dm->owner             = owner;
  dmesh_nodal_to_dm->results_is_getted = PDM_FALSE;
  dmesh_nodal_to_dm->n_mesh            = n_mesh;
  dmesh_nodal_to_dm->link              = malloc( n_mesh * sizeof(_pdm_link_dmesh_nodal_to_dmesh_t*) );
  for(int i_mesh = 0; i_mesh < n_mesh; ++i_mesh) {
    dmesh_nodal_to_dm->link[i_mesh] = _link_dmesh_nodal_to_dmesh_init(owner);
  }

  return (PDM_dmesh_nodal_to_dmesh_t *) dmesh_nodal_to_dm;
}

/**
 * \brief  Add dmesh_nodal
 */
void
PDM_dmesh_nodal_to_dmesh_add_dmesh_nodal
(
        PDM_dmesh_nodal_to_dmesh_t *dmesh_nodal_to_dm,
  const int                         i_mesh,
        PDM_dmesh_nodal_t          *dmesh_nodal
)
{
  dmesh_nodal_to_dm->link[i_mesh]->dmesh_nodal = dmesh_nodal;
}


void
PDM_dmesh_nodal_to_dmesh_get_dmesh
(
        PDM_dmesh_nodal_to_dmesh_t  *dmesh_nodal_to_dm,
  const int                          i_mesh,
        PDM_dmesh_t                **dm
)
{
  *dm = dmesh_nodal_to_dm->link[i_mesh]->dmesh;
  dmesh_nodal_to_dm->results_is_getted = PDM_TRUE;
}


/**
 * \brief  Free
 */
void
PDM_dmesh_nodal_to_dmesh_free
(
  PDM_dmesh_nodal_to_dmesh_t* dmesh_nodal_to_dm
)
{
  // printf("PDM_dmesh_nodal_to_dmesh_free\n");

  if(( dmesh_nodal_to_dm->owner == PDM_OWNERSHIP_KEEP ) ||
     ( dmesh_nodal_to_dm->owner == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE && !dmesh_nodal_to_dm->results_is_getted)){

    for(int i_mesh = 0; i_mesh < dmesh_nodal_to_dm->n_mesh; ++i_mesh) {
      _link_dmesh_nodal_to_dmesh_free(dmesh_nodal_to_dm->link[i_mesh]);
    }
  } else {
    for(int i_mesh = 0; i_mesh < dmesh_nodal_to_dm->n_mesh; ++i_mesh) {
      _link_dmesh_nodal_to_dmesh_free(dmesh_nodal_to_dm->link[i_mesh]);
    }
  }

  free(dmesh_nodal_to_dm->link);

  free(dmesh_nodal_to_dm);

  // printf("PDM_dmesh_nodal_to_dmesh_free end \n");
}



void
PDM_dmesh_nodal_to_dmesh_compute
(
        PDM_dmesh_nodal_to_dmesh_t                 *dmesh_nodal_to_dm,
  const PDM_dmesh_nodal_to_dmesh_transform_t        transform_kind,
  const PDM_dmesh_nodal_to_dmesh_translate_group_t  transform_group_kind
)
{

  for(int i_mesh = 0; i_mesh < dmesh_nodal_to_dm->n_mesh; ++i_mesh) {

    if(dmesh_nodal_to_dm->link[i_mesh]->dmesh_nodal->mesh_dimension == 2) {
      assert(transform_kind != PDM_DMESH_NODAL_TO_DMESH_TRANSFORM_TO_FACE);
    }

    switch (transform_kind) {

      case PDM_DMESH_NODAL_TO_DMESH_TRANSFORM_TO_FACE:
        {
          _generate_faces_from_dmesh_nodal(dmesh_nodal_to_dm->link[i_mesh]);
        }
        break;

      case PDM_DMESH_NODAL_TO_DMESH_TRANSFORM_TO_EDGE:
        {
          _generate_edges_from_dmesh_nodal(dmesh_nodal_to_dm->link[i_mesh]);
        }
        break;
    }
  }

  // Boundary management
  for(int i_mesh = 0; i_mesh < dmesh_nodal_to_dm->n_mesh; ++i_mesh) {
    if(dmesh_nodal_to_dm->link[i_mesh]->dmesh_nodal->dgroup_elmt != NULL) {
      switch (transform_group_kind) {
        case PDM_DMESH_NODAL_TO_DMESH_TRANSLATE_GROUP_NONE:
          break;
        case PDM_DMESH_NODAL_TO_DMESH_TRANSLATE_GROUP_TO_FACE:
          {
            assert(transform_kind == PDM_DMESH_NODAL_TO_DMESH_TRANSFORM_TO_FACE);
            _translate_element_group_to_faces(dmesh_nodal_to_dm->link[i_mesh]);
          }
          break;
        case PDM_DMESH_NODAL_TO_DMESH_TRANSLATE_GROUP_TO_EDGE:
          {
            _translate_element_group_to_edges(dmesh_nodal_to_dm->link[i_mesh]);
            // PDM_error (__FILE__, __LINE__, 0, "PDM_DMESH_NODAL_TO_DMESH_TRANSLATE_GROUP_TO_EDGE not implemented \n");
          }
          break;
        case PDM_DMESH_NODAL_TO_DMESH_TRANSLATE_GROUP_TO_VTX:
          {
            PDM_error (__FILE__, __LINE__, 0, "PDM_DMESH_NODAL_TO_DMESH_TRANSLATE_GROUP_TO_EDGE not implemented \n");
          }
          break;
      }
    }
  }

  // Join management


  // Keep all link between elements and generated entity


}


// API change : PDM_dmesh_nodal_to_dmesh_transform_extract_dmesh
void
PDM_dmesh_nodal_to_dmesh_transform_to_coherent_dmesh
(
        PDM_dmesh_nodal_to_dmesh_t *dmesh_nodal_to_dm,
  const int                         extract_dim
)
{
  assert(extract_dim >= 2);
  // assert(extract_dim == 3); // Not implemented

  for(int i_mesh = 0; i_mesh < dmesh_nodal_to_dm->n_mesh; ++i_mesh) {

    _pdm_link_dmesh_nodal_to_dmesh_t *link        = dmesh_nodal_to_dm->link[i_mesh];
    PDM_dmesh_nodal_t                *dmesh_nodal = link->dmesh_nodal;
    PDM_dmesh_t                      *dmesh       = link->dmesh;

    if(extract_dim == 3) {
      _to_coherent_3d(dmesh_nodal_to_dm->comm, link, dmesh_nodal, dmesh);
    } else if (extract_dim == 2) {
      _to_coherent_2d(dmesh_nodal_to_dm->comm, link, dmesh_nodal, dmesh);
    }
  }
}
