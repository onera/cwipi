#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "pdm.h"
#include "pdm_config.h"
#include "pdm_mpi.h"
#include "pdm_part.h"
#include "pdm_dcube_gen.h"
#include "pdm_logging.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_priv.h"
#include "pdm_part_to_block.h"
#include "pdm_block_to_part.h"
#include "pdm_multi_block_to_part.h"
#include "pdm_distrib.h"
#include "pdm_vtk.h"
#include "pdm_sort.h"
#include "pdm_array.h"
#include "pdm_partitioning_algorithm.h"
#include "pdm_dconnectivity_transform.h"
#include "pdm_dmesh_nodal_elements_utils.h"
#include "pdm_dmesh_nodal_to_dmesh.h"

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

static
void
_exchange_point_list
(
 int            n_group_join,
 int           *group_join_to_join_opp,
 int           *dface_join_idx,
 PDM_g_num_t   *dface_join,
 PDM_g_num_t  **dface_join_opp,
 PDM_MPI_Comm   comm
)
{
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  /*
   * We have for all extraction zone, we need to exchange id
   */
  PDM_g_num_t **distrib_join   = malloc(n_group_join * sizeof(PDM_g_num_t *));
  // PDM_g_num_t **dface_join_opp = malloc(n_group_join * sizeof(PDM_g_num_t *));
  for(int i_group_join = 0; i_group_join < n_group_join; ++i_group_join) {
    int dn_face_join = dface_join_idx[i_group_join+1] - dface_join_idx[i_group_join];
    distrib_join[i_group_join] = PDM_compute_entity_distribution(comm, dn_face_join);
  }

  *dface_join_opp = malloc(dface_join_idx[n_group_join] * sizeof(PDM_g_num_t));
  PDM_g_num_t* _dface_join_opp = *dface_join_opp;

  for(int i_group_join = 0; i_group_join < n_group_join; ++i_group_join) {

    int i_group_join_opp = group_join_to_join_opp[i_group_join];
    int dn_face_join = dface_join_idx[i_group_join+1] - dface_join_idx[i_group_join];
    PDM_g_num_t* distrib_join_cur = distrib_join[i_group_join    ];
    PDM_g_num_t* distrib_join_opp = distrib_join[i_group_join_opp];

    PDM_g_num_t *join_ln_to_gn = malloc(dn_face_join * sizeof(PDM_g_num_t));
    for(int i = 0; i < dn_face_join; ++i) {
      join_ln_to_gn[i] = distrib_join_cur[i_rank] + i + 1;
    }

    /*
     * Exchange
     */
    PDM_g_num_t *blk_dface_join_cur = &dface_join[dface_join_idx[i_group_join    ]];
    PDM_g_num_t *blk_dface_join_opp = &dface_join[dface_join_idx[i_group_join_opp]];

    PDM_block_to_part_t *btp = PDM_block_to_part_create(distrib_join_opp,
                                 (const PDM_g_num_t **) &join_ln_to_gn,
                                                        &dn_face_join,
                                                        1,
                                                        comm);

    int cst_stride = 1;
    PDM_g_num_t* sub_dface_join_opp = &_dface_join_opp[dface_join_idx[i_group_join]];
    PDM_block_to_part_exch(btp,
                           sizeof(PDM_g_num_t),
                           PDM_STRIDE_CST,
                           &cst_stride,
                  (void *) blk_dface_join_opp,
                           NULL,
               (void ** ) &sub_dface_join_opp);

    if(0 == 1) {
      PDM_log_trace_array_long(blk_dface_join_cur, dn_face_join, "dface_join_cur :: ");
      PDM_log_trace_array_long(sub_dface_join_opp, dn_face_join, "dface_join_opp :: ");
    }

    PDM_block_to_part_free(btp);
    free(join_ln_to_gn);
  }

  for(int i_group_join = 0; i_group_join < n_group_join; ++i_group_join) {
    free(distrib_join[i_group_join]);
  }
  free(distrib_join);
}






static
void
_deduce_descending_join
(
 int            n_zone,
 int            n_group_join,
 int           *group_join_to_zone_cur,
 int           *group_join_to_zone_opp,
 int           *group_join_to_join_opp,
 int           *dface_join_idx,
 PDM_g_num_t   *dface_join,
 PDM_g_num_t   *dface_join_opp,
 int          **dface_vtx_idx,
 PDM_g_num_t  **dface_vtx,
 PDM_g_num_t  **extract_face_distribution,
 PDM_g_num_t  **extract_vtx_distribution,
 int          **dextract_face_vtx_idx,
 PDM_g_num_t  **dextract_face_vtx,
 int          **dextract_face_group_id,
 int          **dextract_face_group_sens,
 PDM_g_num_t  **dextract_face_join,
 PDM_g_num_t  **dextract_face_join_opp,
 PDM_g_num_t  **dparent_face_g_num,
 PDM_g_num_t  **dparent_vtx_g_num,
 PDM_g_num_t  **pextract_old_to_new,
 double       **dextract_vtx_coord,
 PDM_MPI_Comm   comm
)
{
  PDM_UNUSED(dface_join_idx);
  PDM_UNUSED(dface_join);
  PDM_UNUSED(dface_vtx_idx);
  PDM_UNUSED(dface_vtx);
  PDM_UNUSED(extract_face_distribution);
  PDM_UNUSED(extract_vtx_distribution);
  PDM_UNUSED(dextract_face_vtx_idx);
  PDM_UNUSED(dextract_face_vtx);
  PDM_UNUSED(dparent_face_g_num);
  PDM_UNUSED(dparent_vtx_g_num);
  PDM_UNUSED(pextract_old_to_new);
  PDM_UNUSED(dextract_vtx_coord);

  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);


  PDM_g_num_t **dedge_distrib  = malloc(n_zone * sizeof(PDM_g_num_t *));
  int         **dedge_vtx_idx  = malloc(n_zone * sizeof(int         *));
  PDM_g_num_t **dedge_vtx      = malloc(n_zone * sizeof(PDM_g_num_t *));
  int         **dedge_face_idx = malloc(n_zone * sizeof(int         *));
  PDM_g_num_t **dedge_face     = malloc(n_zone * sizeof(PDM_g_num_t *));

  PDM_g_num_t **key_ln_to_gn    = (PDM_g_num_t **) malloc( n_zone * sizeof(PDM_g_num_t*));

  /* Connectivity */
  PDM_g_num_t **data_send_connect    = (PDM_g_num_t **) malloc( n_zone * sizeof(PDM_g_num_t*));
  PDM_g_num_t **data_send_edge_g_num = (PDM_g_num_t **) malloc( n_zone * sizeof(PDM_g_num_t*));
  PDM_g_num_t **data_send_group      = (PDM_g_num_t **) malloc( n_zone * sizeof(PDM_g_num_t*));
  PDM_g_num_t **data_send_sens       = (PDM_g_num_t **) malloc( n_zone * sizeof(PDM_g_num_t*));
  PDM_g_num_t **data_send_face_g_num = (PDM_g_num_t **) malloc( n_zone * sizeof(PDM_g_num_t*));

  int         **stride_one      = (int         **) malloc( n_zone * sizeof(int        *));
  int         **stride_two      = (int         **) malloc( n_zone * sizeof(int        *));
  int         **stride_four     = (int         **) malloc( n_zone * sizeof(int        *));
  int         **zone_id         = (int         **) malloc( n_zone * sizeof(int        *));

  int          *dn_edge       = (int          *) malloc( n_zone * sizeof(int         ));
  int          *dn_internal_edge       = (int          *) malloc( n_zone * sizeof(int         ));
  int          *dn_external_edge       = (int          *) malloc( n_zone * sizeof(int         ));

  // PDM_g_num_t key_mod = ; // Comment definir un bon key_mod ?
  for(int i_zone = 0; i_zone < n_zone; ++i_zone) {
    /*
     * Generate edge numbering
     */
    int dn_face = extract_face_distribution[i_zone][i_rank+1] - extract_face_distribution[i_zone][i_rank];
    int n_edge_elt_tot = dextract_face_vtx_idx[i_zone][dn_face];

    PDM_g_num_t* tmp_dface_edge         = (PDM_g_num_t *) malloc(     n_edge_elt_tot    * sizeof(PDM_g_num_t) );
    int*         tmp_parent_elmt_pos    = (int         *) malloc(     n_edge_elt_tot    * sizeof(int        ) );
    int*         tmp_dface_edge_vtx_idx = (int         *) malloc( ( n_edge_elt_tot + 1) * sizeof(int        ) );
    PDM_g_num_t* tmp_dface_edge_vtx     = (PDM_g_num_t *) malloc( 2 * n_edge_elt_tot    * sizeof(PDM_g_num_t) );

    int n_elmt_current = 0;
    int n_edge_current = 0;
    tmp_dface_edge_vtx_idx[0] = 0;
    PDM_poly2d_decomposes_edges(dn_face,
                                &n_elmt_current,
                                &n_edge_current,
                                extract_face_distribution[i_zone][i_rank],
                                -1,
                                dextract_face_vtx[i_zone],
                                dextract_face_vtx_idx[i_zone],
                                tmp_dface_edge_vtx_idx,
                                tmp_dface_edge_vtx,
                                tmp_dface_edge,
                                NULL,
                                NULL,
                                tmp_parent_elmt_pos);
    assert(n_edge_current == n_edge_elt_tot);

    /*
     *  Compute edges connectivity
     */
    PDM_generate_entitiy_connectivity_raw(comm,
                                          extract_vtx_distribution[i_zone][n_rank],
                                          n_edge_elt_tot,
                                          tmp_dface_edge,
                                          tmp_dface_edge_vtx_idx,
                                          tmp_dface_edge_vtx,
                                          &dn_edge[i_zone],
                                          &dedge_distrib [i_zone],
                                          &dedge_vtx_idx [i_zone],
                                          &dedge_vtx     [i_zone],
                                          &dedge_face_idx[i_zone],
                                          &dedge_face    [i_zone]);
    free(tmp_parent_elmt_pos    );

    // Count the number of internal & external edges
    dn_internal_edge[i_zone] = 0;
    for(int i_edge = 0; i_edge < dn_edge[i_zone]; ++i_edge) {
      int n_face_this_edge = dedge_face_idx[i_zone][i_edge+1] - dedge_face_idx[i_zone][i_edge];
      dn_internal_edge[i_zone] += (int) (n_face_this_edge > 1);
    }
    dn_external_edge[i_zone] = dn_edge[i_zone] - dn_internal_edge[i_zone];
    log_trace("dn internal edges for zone %i is %i \n", i_zone, dn_internal_edge[i_zone]);
    log_trace("dn external edges for zone %i is %i \n", i_zone, dn_external_edge[i_zone]);

    if(1 == 1) {
      PDM_log_trace_array_long(dedge_vtx_idx[i_zone], dn_edge[i_zone]+1, "dedge_vtx_idx :: ");
      PDM_log_trace_array_long(dedge_vtx   [i_zone], dedge_vtx_idx [i_zone][dn_edge[i_zone]], "dedge_vtx :: ");
      PDM_log_trace_array_long(dedge_face_idx[i_zone], dn_edge[i_zone]+1, "dedge_face_idx :: ");
      PDM_log_trace_array_long(dedge_face  [i_zone], dedge_face_idx[i_zone][dn_edge[i_zone]], "dedge_face :: ");
      // PDM_log_trace_array_long(key_ln_to_gn[i_zone], dn_face, "key_ln_to_gn :: ");
    }

    /*
     * Echange des données pour construire la clé
     *   - exchange throw dedge_face
     */
    PDM_g_num_t *dedge_face_abs = (PDM_g_num_t *) malloc(dedge_face_idx[i_zone][dn_edge[i_zone]] * sizeof(PDM_g_num_t));
    int         *dedge_face_sgn = (int         *) malloc(dedge_face_idx[i_zone][dn_edge[i_zone]] * sizeof(int        ));
    for(int i = 0; i < dedge_face_idx[i_zone][dn_edge[i_zone]]; ++i) {
      dedge_face_abs[i] = PDM_ABS (dedge_face[i_zone][i]);
      dedge_face_sgn[i] = PDM_SIGN(dedge_face[i_zone][i]);
    }

    PDM_block_to_part_t *btp = PDM_block_to_part_create(extract_face_distribution[i_zone],
                                 (const PDM_g_num_t **) &dedge_face_abs,
                                                        &dedge_face_idx[i_zone][dn_edge[i_zone]],
                                                        1,
                                                        comm);

    /*
     * Exchange
     */
    int* dedge_face_group_id   = malloc(dedge_face_idx[i_zone][dn_edge[i_zone]] * sizeof(int        ));
    int* dedge_face_group_sens = malloc(dedge_face_idx[i_zone][dn_edge[i_zone]] * sizeof(int        ));
    int* dedge_face_join       = malloc(dedge_face_idx[i_zone][dn_edge[i_zone]] * sizeof(PDM_g_num_t));
    int* dedge_face_join_opp   = malloc(dedge_face_idx[i_zone][dn_edge[i_zone]] * sizeof(PDM_g_num_t));

    int cst_stride = 1;
    PDM_block_to_part_exch(btp,
                           sizeof(int),
                           PDM_STRIDE_CST,
                           &cst_stride,
                  (void *) dextract_face_group_id[i_zone],
                           NULL,
               (void ** ) &dedge_face_group_id);

    PDM_block_to_part_exch(btp,
                           sizeof(int),
                           PDM_STRIDE_CST,
                           &cst_stride,
                  (void *) dextract_face_group_sens[i_zone],
                           NULL,
               (void ** ) &dedge_face_group_sens);

    PDM_block_to_part_exch(btp,
                           sizeof(PDM_g_num_t),
                           PDM_STRIDE_CST,
                           &cst_stride,
                  (void *) dextract_face_join[i_zone],
                           NULL,
               (void ** ) &dedge_face_join);

    PDM_block_to_part_exch(btp,
                           sizeof(PDM_g_num_t),
                           PDM_STRIDE_CST,
                           &cst_stride,
                  (void *) dextract_face_join_opp[i_zone],
                           NULL,
               (void ** ) &dedge_face_join_opp);

    PDM_block_to_part_free(btp);

    free(dedge_face_abs);
    free(dedge_face_sgn);

    /*
     * Remplissage des clés + préparation buffer d'envoi
     *    - Une clé par edge pour unifier !!
     *    - Maybe besoin de signé le dgroup_face pour identifié les vtx aprés !!!
     *  Attention les l'ordre des faces dans le dfaces_group n'est pas le même que dans l'extraction !!!!!
     *  Il faut appliqué le pextract_old_to_new
     *
     *   -> Le mieux est peut-être de prétraité pour ne plus avoir le pb aprés ...
     */
    key_ln_to_gn[i_zone] = malloc(dn_internal_edge[i_zone] * sizeof(PDM_g_num_t)); // Toutes les edges ont une clé car tout vient de l'extraction

    stride_one          [i_zone] = malloc( (    dn_internal_edge[i_zone]                       ) * sizeof(int        ));
    stride_two          [i_zone] = malloc( (    dn_internal_edge[i_zone]                       ) * sizeof(int        ));
    stride_four         [i_zone] = malloc( (    dn_internal_edge[i_zone]                       ) * sizeof(int        ));
    zone_id             [i_zone] = malloc( (    dn_internal_edge[i_zone]                       ) * sizeof(int        ));
    data_send_connect   [i_zone] = malloc( (4 * dn_internal_edge[i_zone]) * sizeof(PDM_g_num_t));
    data_send_edge_g_num[i_zone] = malloc( (    dn_internal_edge[i_zone]                        ) * sizeof(PDM_g_num_t));
    data_send_group     [i_zone] = malloc( (4 * dn_internal_edge[i_zone]) * sizeof(PDM_g_num_t));
    data_send_sens      [i_zone] = malloc( (    2*dn_internal_edge[i_zone]) * sizeof(PDM_g_num_t));
    data_send_face_g_num[i_zone] = malloc( (    2*dn_internal_edge[i_zone]) * sizeof(PDM_g_num_t));

    /*
     *  Let's build key
     */
    // int idx_write     = 0;
    int idx_write_connect    = 0;
    int idx_write_group      = 0;
    int idx_write_sens       = 0;
    int idx_write_face_g_num = 0;
    int idx_write_key        = 0;

    // Je pense qu'il faut speraer le traitement bord + interne
    //  --> On peut le faire directiement ici !!


    int i_int_edge = 0;
    for(int i_edge = 0; i_edge < dn_edge[i_zone]; ++i_edge) {

      int n_face_this_edge = dedge_face_idx[i_zone][i_edge+1] - dedge_face_idx[i_zone][i_edge];
      if(dedge_face_idx[i_zone][i_edge+1] - dedge_face_idx[i_zone][i_edge] == 1) {
        continue;
      }
      assert (n_face_this_edge == 2);

      stride_one          [i_zone][i_int_edge] = 1;
      stride_two          [i_zone][i_int_edge] = 2;
      stride_four         [i_zone][i_int_edge] = 4;
      zone_id             [i_zone][i_int_edge] = i_zone;

      data_send_edge_g_num[i_zone][i_int_edge] = dedge_distrib [i_zone][i_rank] + i_edge + 1;

      PDM_g_num_t key = 0;
      for(int j = dedge_face_idx[i_zone][i_edge]; j < dedge_face_idx[i_zone][i_edge+1]; ++j) {
        key += (dedge_face_join[j] + dedge_face_join_opp[j]);
        data_send_connect[i_zone][idx_write_connect++] = dedge_face_join    [j];
        data_send_connect[i_zone][idx_write_connect++] = dedge_face_join_opp[j];
      }

      /* Send group id */
      for(int j = dedge_face_idx[i_zone][i_edge]; j < dedge_face_idx[i_zone][i_edge+1]; ++j) {
        data_send_group[i_zone][idx_write_group++] = dedge_face_group_id[j];
        int group_join_opp = group_join_to_join_opp[dedge_face_group_id[j]];
        data_send_group[i_zone][idx_write_group++] = group_join_opp;
      }

      /* Send face_g_num */
      for(int j = dedge_face_idx[i_zone][i_edge]; j < dedge_face_idx[i_zone][i_edge+1]; ++j) {
        data_send_face_g_num[i_zone][idx_write_face_g_num++] = dedge_face[i_zone][j];
      }

      /* Send group sens - Caution / 2 */
      for(int j = dedge_face_idx[i_zone][i_edge]; j < dedge_face_idx[i_zone][i_edge+1]; ++j) {
        data_send_sens[i_zone][idx_write_sens++] = dedge_face_group_sens[j];
      }

      // key_ln_to_gn[i_zone][i_edge] = key % key_mod + 1; // TODO
      key_ln_to_gn[i_zone][i_int_edge] = key;
      ++i_int_edge;

    }

    free(dedge_face_group_id);
    free(dedge_face_join);
    free(dedge_face_join_opp);
    free(dedge_face_group_sens);

    if(1 == 1) {
      PDM_log_trace_array_long(key_ln_to_gn[i_zone], dn_internal_edge[i_zone], "key_ln_to_gn :: ");
    }
  }


  /*
   * Distributed hash table
   */
  PDM_part_to_block_t* ptb = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                      PDM_PART_TO_BLOCK_POST_MERGE,
                                                      1.,
                                                      key_ln_to_gn,
                                                      NULL,
                                                      dn_internal_edge,
                                                      n_zone,
                                                      comm);

  // Get protocol data
  int blk_size = PDM_part_to_block_n_elt_block_get(ptb);
  // This one will contain 2 for non conflicting hash key, more otherwise
  int *gnum_n_occurences = PDM_part_to_block_block_gnum_count_get(ptb);

  int gnum_n_occurences_tot = 0;
  for (int k = 0; k < blk_size; k++) {
    gnum_n_occurences_tot += gnum_n_occurences[k];
  }

  int *blk_entity_idx    = (int *) malloc((gnum_n_occurences_tot + 1)*sizeof(int));
  int *blk_data_face_idx = (int *) malloc((gnum_n_occurences_tot + 1)*sizeof(int));
  for (int k = 0; k < gnum_n_occurences_tot + 1; k++) {
    blk_entity_idx[k]    = 4*k;
    blk_data_face_idx[k] = 2*k;
  }


  // Send data : zone_id, edge gnum, face & facedonor id, face & facedonor group id, face_sens and faces_gnum

  int *unused_recv_stride = NULL;
  int *blk_zone_id   = NULL;
  int exch_size = PDM_part_to_block_exch(ptb,
                                         sizeof(int),
                                         PDM_STRIDE_VAR,
                                         -1,
                                         stride_one,
                               (void **) zone_id,
                                         &unused_recv_stride,
                               (void **) &blk_zone_id);
  free(unused_recv_stride); // Same as gnum_n_occurences 
  assert (exch_size == gnum_n_occurences_tot);

  PDM_g_num_t* blk_edge_g_num           = NULL;
  exch_size = PDM_part_to_block_exch(ptb,
                                     sizeof(PDM_g_num_t),
                                     PDM_STRIDE_VAR,
                                     -1,
                                     stride_one,
                           (void **) data_send_edge_g_num,
                                     &unused_recv_stride,
                           (void **) &blk_edge_g_num);
  free(unused_recv_stride); // Same as gnum_n_occurences 
  assert (exch_size == gnum_n_occurences_tot);

  PDM_g_num_t* blk_data_sens   = NULL;
  exch_size = PDM_part_to_block_exch(ptb,
                                     sizeof(int),
                                     PDM_STRIDE_VAR,
                                     -1,
                                     stride_two,
                           (void **) data_send_sens,
                                     &unused_recv_stride,
                           (void **) &blk_data_sens);
  free(unused_recv_stride); // Same as 2*gnum_n_occurences 
  assert (exch_size == 2*gnum_n_occurences_tot);

  PDM_g_num_t* blk_data_face_g_num   = NULL;
  exch_size = PDM_part_to_block_exch(ptb,
                                     sizeof(PDM_g_num_t),
                                     PDM_STRIDE_VAR,
                                     -1,
                                     stride_two,
                           (void **) data_send_face_g_num,
                                     &unused_recv_stride,
                           (void **) &blk_data_face_g_num);
  free(unused_recv_stride); // Same as 2*gnum_n_occurences
  assert (exch_size == 2*gnum_n_occurences_tot);

  PDM_g_num_t* blk_data_connect   = NULL;
  exch_size = PDM_part_to_block_exch(ptb,
                                     sizeof(PDM_g_num_t),
                                     PDM_STRIDE_VAR,
                                     -1,
                                    stride_four,
                          (void **) data_send_connect,
                                    &unused_recv_stride,
                          (void **) &blk_data_connect);
  free(unused_recv_stride); // Same as 4*gnum_n_occurences 
  assert (exch_size == 4*gnum_n_occurences_tot);

  PDM_g_num_t* blk_data_group   = NULL;
  exch_size = PDM_part_to_block_exch(ptb,
                                     sizeof(int),
                                     PDM_STRIDE_VAR,
                                     -1,
                                     stride_four,
                           (void **) data_send_group,
                                     &unused_recv_stride,
                           (void **) &blk_data_group);
  free(unused_recv_stride); // Same as 4*gnum_n_occurences 
  assert (exch_size == 4*gnum_n_occurences_tot);

  /*
   * Free
   */
  for(int i_zone = 0; i_zone < n_zone; ++i_zone) {
    free(key_ln_to_gn        [i_zone]);
    free(data_send_connect   [i_zone]);
    free(data_send_group     [i_zone]);
    free(data_send_edge_g_num[i_zone]);
    free(data_send_sens      [i_zone]);
    free(stride_one          [i_zone]);
    free(stride_two          [i_zone]);
    free(stride_four         [i_zone]);
    free(data_send_face_g_num[i_zone]);
    free(zone_id             [i_zone]);
  }
  free(key_ln_to_gn        );
  free(data_send_connect   );
  free(data_send_group     );
  free(data_send_edge_g_num);
  free(data_send_sens      );
  free(data_send_face_g_num);
  free(stride_one          );
  free(stride_two          );
  free(stride_four         );
  free(zone_id             );


  if(1 == 1) {
    PDM_log_trace_array_int(gnum_n_occurences   , blk_size               , "gnum_n_occurences   :: ");

    PDM_log_trace_array_int(blk_zone_id         , gnum_n_occurences_tot  , "blk_zone_id         :: ");
    PDM_log_trace_array_long(blk_edge_g_num     , gnum_n_occurences_tot  , "blk_edge_g_num      :: ");
    PDM_log_trace_array_long(blk_data_face_g_num, 2*gnum_n_occurences_tot, "blk_data_face_g_num :: ");
    PDM_log_trace_array_long(blk_data_sens      , 2*gnum_n_occurences_tot, "blk_data_sens       :: ");
    PDM_log_trace_array_long(blk_data_connect   , 4*gnum_n_occurences_tot, "blk_data_connect    :: ");
    PDM_log_trace_array_long(blk_data_group     , 4*gnum_n_occurences_tot, "blk_data_group      :: ");
  }

  /*
   * Post-treatment
   */
  int n_max_entity_per_key = 0;
  for(int i = 0; i < blk_size; ++i) {
    n_max_entity_per_key = PDM_MAX(gnum_n_occurences   [i], n_max_entity_per_key);
  }
  int n_max_connec = 4*n_max_entity_per_key;


  //Number of edge to treat on each zone
  int* zone_id_n = PDM_array_zeros_int(n_zone);
  for (int i_edge = 0; i_edge < gnum_n_occurences_tot; i_edge++) {
    zone_id_n[blk_zone_id[i_edge]]++;
  }

  if(1 == 1) {
    PDM_log_trace_array_int(zone_id_n        , n_zone                , "zone_id_n         :: ");
    PDM_log_trace_array_int(blk_entity_idx   , gnum_n_occurences_tot+1   , "blk_entity_idx    :: ");
    PDM_log_trace_array_int(blk_data_face_idx, gnum_n_occurences_tot+1, "blk_data_face_idx :: ");
  }

  printf("n_max_entity_per_key = %i \n", n_max_entity_per_key);
  printf("n_max_connec         = %i \n", n_max_connec);

  // PDM_g_num_t** loc_connec       = (PDM_g_num_t *) malloc(  n_max_entity_per_key * sizeof(PDM_g_num_t *) );
  int*          already_treat    = (int         *) malloc(  n_max_entity_per_key * sizeof(int          ) );
  int*          same_entity_idx  = (int         *) malloc(  n_max_entity_per_key * sizeof(int          ) );
  int*          sens_entity      = (int         *) malloc(  n_max_entity_per_key * sizeof(int          ) );

  PDM_g_num_t **ln_to_gn_res     = malloc(n_zone * sizeof(PDM_g_num_t *));
  PDM_g_num_t **results_edge     = malloc(n_zone * sizeof(PDM_g_num_t *));
  PDM_g_num_t **results_edge_opp = malloc(n_zone * sizeof(PDM_g_num_t *));
  int         **results_zone_opp = malloc(n_zone * sizeof(int         *));
  for(int i = 0; i < n_zone; ++i) {
    ln_to_gn_res    [i] = malloc(40 * zone_id_n[i] * sizeof(PDM_g_num_t));
    results_edge    [i] = malloc(zone_id_n[i] * sizeof(PDM_g_num_t));
    results_edge_opp[i] = malloc(zone_id_n[i] * sizeof(PDM_g_num_t));
    results_zone_opp[i] = malloc(zone_id_n[i] * sizeof(int        ));
  }

  /* Reset to fill */
  PDM_array_reset_int(zone_id_n, n_zone, 0);

  int idx  = 0;
  for(int i_key = 0; i_key < blk_size; ++i_key) {

    int n_matching_edge = gnum_n_occurences[i_key];

    log_trace(" i_key = %i | n_matching_edge = %i \n", i_key, n_matching_edge);

    /* Reset */
    PDM_array_reset_int(already_treat, n_matching_edge, -1);

    /* Loop over all entitys in conflict and sort all */
    for(int i_entity = 0; i_entity < n_matching_edge; ++i_entity) {

      //Each internal edge comes with 2 faces -> 4 data
      int beg    = 4*(idx+i_entity);

      // Caution inplace sort !!!!
      PDM_sort_long(&blk_data_connect[beg], NULL, 4);
      PDM_sort_int (&blk_data_group  [beg], NULL, 4);

      if(1 == 1) {
        PDM_log_trace_array_long(&blk_data_connect[beg], 4, "blk_data_connect (sort) :: ");
        PDM_log_trace_array_int (&blk_data_group  [beg], 4, "blk_data_group   (sort) :: ");
      }
    }

    /*
     *  Identify pair or invalid other
     */
    for(int i_entity = 0; i_entity < n_matching_edge; ++i_entity) {
      int beg1    = 4*(idx+i_entity);

      int n_same         = 0;
      int i_entity2_same = -1;

      if(already_treat[i_entity] == 1) {
        continue;
      }

      for(int i_entity2 = i_entity+1; i_entity2 < n_matching_edge; ++i_entity2) {

        if(already_treat[i_entity2] == -1) {
          int beg2    = 4*(idx+i_entity2);

          if(!PDM_array_are_equal_int(&blk_data_group[beg1], &blk_data_group[beg2], 4)) {
            continue;
          }

          if(!PDM_array_are_equal_int(&blk_data_connect[beg1], &blk_data_connect[beg2], 4)) {
            continue;
          }

          already_treat[i_entity2] = 1;
          same_entity_idx[n_same++] = i_entity2;
          i_entity2_same = i_entity2;

        }
      } /* End for i_entity2 */
      assert(n_same == 1);

      log_trace("i_entity = %i | i_entity2_same = %i | n_same = %i \n", i_entity, i_entity2_same, n_same);

      /*
       * Renvoie des resultats :
       *    - Par edge --> Il faut donc également le numero de zones
       */
      int edge_idx     = (idx+i_entity);
      int edge_idx_opp = (idx+i_entity2_same);

      int i_zone_cur = blk_zone_id[edge_idx];
      int i_zone_opp = blk_zone_id[edge_idx_opp];

      // Set data for edge
      results_edge    [i_zone_cur][zone_id_n[i_zone_cur]  ] = blk_edge_g_num[edge_idx];
      results_edge_opp[i_zone_cur][zone_id_n[i_zone_cur]  ] = blk_edge_g_num[edge_idx_opp];
      results_zone_opp[i_zone_cur][zone_id_n[i_zone_cur]++] = blk_zone_id[edge_idx_opp];

      // Set data for opposite edge
      results_edge    [i_zone_opp][zone_id_n[i_zone_opp]  ] = blk_edge_g_num[edge_idx_opp];
      results_edge_opp[i_zone_opp][zone_id_n[i_zone_opp]  ] = blk_edge_g_num[edge_idx];
      results_zone_opp[i_zone_opp][zone_id_n[i_zone_opp]++] = blk_zone_id[edge_idx];

      // Renvoi de tout les edges candidats à travers les faces ???
      already_treat[i_entity] = 1;
    }

    idx  += n_matching_edge;
  }

  /*
   *  Part to block avec face_extraction_dsitrib
   */
  for(int i = 0; i < n_zone; ++i) {
    log_trace("Results for zone %i\n", i);
    // PDM_log_trace_array_long(ln_to_gn_res    [i], zone_id_n[i], "ln_to_gn_res     :: ");
    PDM_log_trace_array_long(results_edge    [i], zone_id_n[i], "results_edge     :: ");
    PDM_log_trace_array_long(results_edge_opp[i], zone_id_n[i], "results_edge_opp :: ");
    PDM_log_trace_array_int (results_zone_opp[i], zone_id_n[i], "results_zone_opp :: ");
  }
  PDM_part_to_block_free(ptb);

  // Construction du dface_edge et du dface_opp_edge
  //First build dface to edge connectivity -- only face having undermined edge are needed but do all for now
  int **dface_edge_idx = (int **) malloc(n_zone*sizeof(int *));
  PDM_g_num_t **dface_edge = (PDM_g_num_t **) malloc(n_zone*sizeof(PDM_g_num_t *));
  for (int i_zone = 0; i_zone < n_zone; i_zone++) {
    PDM_dconnectivity_transpose(comm,
                                dedge_distrib[i_zone],
                                extract_face_distribution[i_zone], // face distribution
                                dedge_face_idx[i_zone],
                                dedge_face[i_zone],
                                1,
                               &dface_edge_idx[i_zone],
                               &dface_edge[i_zone]);

    int dn_face = extract_face_distribution[i_zone][i_rank+1] - extract_face_distribution[i_zone][i_rank];
    /*PDM_g_num_t *dface_edge_abs = (PDM_g_num_t *) malloc(dface_edge_idx[dn_face] * sizeof(PDM_g_num_t));*/
    /*for (int i = 0; i < dface_edge_idx[dn_face]; i++) {*/
      /*dface_edge_abs[i] = PDM_ABS(dface_edge[i]);*/
    /*}*/

    PDM_log_trace_array_long(extract_face_distribution[i_zone], n_rank +1, "face distribution:: ");
    PDM_log_trace_connectivity_long(dface_edge_idx[i_zone], dface_edge[i_zone], dn_face, "dface_edge :: ");
  }

  //Reconstruire un dface_edge total

  PDM_g_num_t *multi_distri_idx = (PDM_g_num_t *) malloc((n_zone+1) * sizeof(PDM_g_num_t));
  multi_distri_idx[0] = 0;
  int n_face_tot = 0;
  for (int i_zone = 0; i_zone < n_zone; i_zone++) {
    n_face_tot += extract_face_distribution[i_zone][i_rank+1] - extract_face_distribution[i_zone][i_rank]; //For this proc
    multi_distri_idx[i_zone+1] = multi_distri_idx[i_zone] + extract_face_distribution[i_zone][n_rank];  //Global
  }
  PDM_g_num_t *multi_gnum        = malloc(n_face_tot * sizeof(PDM_g_num_t));

  n_face_tot = 0;
  for (int i_zone = 0; i_zone < n_zone; i_zone++) {
    int dn_face = extract_face_distribution[i_zone][i_rank+1] - extract_face_distribution[i_zone][i_rank];

    for (int i_face = 0; i_face < dn_face; i_face++) {
      PDM_g_num_t opp_face_loc = dextract_face_join_opp[i_zone][pextract_old_to_new[i_zone][i_face]-1];
      int         opp_zone_id  = group_join_to_zone_opp[dextract_face_group_id[i_zone][i_face]];
      log_trace("Zone %i face %i : loc opp is %d, opp id is %d \n", i_zone, i_face, opp_face_loc, opp_zone_id);
      multi_gnum[n_face_tot++] = i_face + multi_distri_idx[opp_zone_id] + 1; //Numéro de la face opposée shifté p/ au nombre de block
    }

  }

  PDM_log_trace_array_long(multi_distri_idx, n_zone+1, "multi_distri_idx :: ");
  PDM_log_trace_array_long(multi_gnum, n_face_tot, "multi_gnum :: ");

  PDM_multi_block_to_part_t *mptb = PDM_multi_block_to_part_create(multi_distri_idx, //todo
                                                                   n_zone,
                                            (const PDM_g_num_t **) extract_face_distribution,
                                            (const PDM_g_num_t **)&multi_gnum,
                                                                  &n_face_tot,
                                                                   1,
                                                                   comm);
  PDM_multi_block_to_part_free(mptb);
  free(multi_distri_idx);
  free(multi_gnum);


  // On considère les edges unifiés comme une partition, que l'on renvoie vers la distribution
  // de edge
  for (int i_zone = 0; i_zone < n_zone; i_zone++) {
    // Part to block avec la distribution des edges dedge_distrib [i_zone],
    PDM_part_to_block_t* ptb_z = PDM_part_to_block_create2(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                           PDM_PART_TO_BLOCK_POST_NOTHING,
                                                           1.,
                                                          &(results_edge[i_zone]),
                                                           dedge_distrib[i_zone],
                                                          &(zone_id_n[i_zone]),
                                                           1,
                                                           comm);

    PDM_g_num_t *dedge_gnum     = malloc(PDM_part_to_block_n_elt_block_get(ptb_z) * sizeof(PDM_g_num_t));
    PDM_g_num_t *dedge_gnum_opp = NULL;

    PDM_g_num_t *dedge_gnum_tmp = PDM_part_to_block_block_gnum_get(ptb_z);
    memcpy(dedge_gnum, dedge_gnum_tmp, PDM_part_to_block_n_elt_block_get(ptb_z) * sizeof(PDM_g_num_t));

    PDM_part_to_block_exch(ptb_z,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_CST,
                          1,
                          NULL,
               (void **) &(results_edge_opp[i_zone]),
                          NULL,
                (void **) &dedge_gnum_opp);

    PDM_part_to_block_free(ptb_z);

    // Rebuild idx array (no data on external edges)
    int *dedge_gnum_idx = (int *) malloc ((dn_edge[i_zone] + 1) * sizeof(int)); //This one will be useless probably
    int count = 0;
    dedge_gnum_idx[0] = 0;
    for (int i = 0; i < dn_edge[i_zone]; i++) {
      if (i + dedge_distrib[i_zone][i_rank] + 1 == dedge_gnum[count]) {
        count++;
      }
      dedge_gnum_idx[i+1] = count;
    }
    int *dedge_gnum_n = (int *) malloc ((dn_edge[i_zone]) * sizeof(int));
    for (int i = 0; i < dn_edge[i_zone]; i ++)
      dedge_gnum_n[i] = dedge_gnum_idx[i+1] - dedge_gnum_idx[i];

    PDM_log_trace_array_long(dedge_gnum_idx, dn_edge[i_zone]+1, "dedge gnum  idx:: ");
    PDM_log_trace_array_long(dedge_gnum, 4, "dedge gnum :: ");
    PDM_log_trace_array_long(dedge_gnum_opp, 4, "dedge gnum_opp :: ");



    //Il faut obtenir le dface_opp_edge

    // Request opposite edge ids
    /*
    PDM_block_to_part_t *btp_z = PDM_block_to_part_create(dedge_distrib[i_zone],
                                  (const PDM_g_num_t **) &dface_edge_abs,
                                                         &dface_edge_idx[dn_face],
                                                         1,
                                                         comm);

    int         **pedge_gnum_n   = NULL;
    PDM_g_num_t **pedge_gnum_opp = NULL;
    PDM_block_to_part_exch2(btp_z,
                            sizeof(PDM_g_num_t),
                            PDM_STRIDE_VAR,
                            dedge_gnum_n,
                            dedge_gnum_opp,
                           &pedge_gnum_n,
                (void ***) &pedge_gnum_opp);
    PDM_block_to_part_free(btp_z);

    int n_data_recv = 0;
    for (int i = 0; i < dface_edge_idx[dn_face]; i++)
      n_data_recv += pedge_gnum_n[0][i];
    
    PDM_log_trace_array_int(pedge_gnum_n[0], dface_edge_idx[dn_face], "pedge_gnum_n::");
    PDM_log_trace_array_int(pedge_gnum_opp[0], n_data_recv, "pedge_gnum_opp::");
    */
  }



  /*
   * A verifier que les edges qui sont indeterminé pointe sur une seule face !!!
   *   Quand on decompose il faut également envoyer le numero de la face (de l'extraction)
   */
  free(zone_id_n);

  for(int i = 0; i < n_zone; ++i) {
    free(ln_to_gn_res    [i]);
    free(results_edge    [i]);
    free(results_edge_opp[i]);
    free(results_zone_opp[i]);
  }
  free(ln_to_gn_res    );
  free(results_edge    );
  free(results_edge_opp);
  free(results_zone_opp);
  free(dn_internal_edge);
  free(dn_external_edge);

  // free(loc_entity_1    );
  // free(loc_entity_2    );
  free(already_treat   );
  free(same_entity_idx );
  free(sens_entity     );


  free(blk_data_connect  );
  free(blk_entity_idx);
  free(blk_data_face_idx);

  free(blk_data_group  );
  free(blk_data_sens);
  free(blk_edge_g_num);
  free(blk_zone_id);
  free(blk_data_face_g_num);


  PDM_part_to_block_free(ptb);

  /*
   * Renvoie des resultats
   *   btp avec les key_ln_to_gn
   *   Attention aux block partiels !!!
   */

  /*
   * Pour les edges recalitrant
   *   --> On extract les faces dont les edges_opp ne sont pas referencé
   *   --> Normalement il y au moins un edge par face qui est connecté
   *   --> On échange le dface_edge + l'information du sens devrait permmettre de tout reorienté
   *       ==> On trouve le edge commun puis on tourne dans le sens de l'orientation respective
   */


  /*
   * Pour les coins coins il faut probablement un echange en plus entre les graph 2 à 2 puis unifier
   * Maybe utiliser le PDM_part_dentity_group_to_pentity_group
   * A la fin de l'algo un vtx est potentiellement lié à plusieurs groupes.
   * On peut merger avec un part_to_block tout les tags d'un vertex en appliquant sur les PL et PLD
   * --> Ca gère le coin
   */
  for(int i_zone = 0; i_zone < n_zone; ++i_zone) {
    free(dedge_distrib  [i_zone]);
    free(dedge_vtx_idx  [i_zone]);
    free(dedge_vtx      [i_zone]);
    free(dedge_face_idx [i_zone]);
    free(dedge_face     [i_zone]);
  }

  free(dedge_distrib);
  free(dedge_vtx_idx);
  free(dedge_vtx);
  free(dedge_face_idx);
  free(dedge_face);
  free(dn_edge);
}

static void
merge_zone
(
 int            n_zone,
 int            n_group_join,
 int           *group_join_to_zone_cur,
 int           *group_join_to_zone_opp,
 int           *group_join_to_join_opp,
 int           *dface_join_idx,
 PDM_g_num_t   *dface_join,
 PDM_g_num_t   *dface_join_opp,
 int          **dface_vtx_idx,
 PDM_g_num_t  **dface_vtx,
 PDM_g_num_t  **dface_cell,
 PDM_MPI_Comm   comm
)
{

  /*
   * On a plusieurs connectivités distribué  :
   *   -> On veut unifier ce qui est commun
   *   -> Et actualiser les connectivités
   *
   *  PDM_multi_block_to_part
   *    -> Avec face_ln_to_gn = implicite face_ln_to_gn - face à supprimer
   *          --> PDM_redistribute
   *
   */


  /*
   * Par connectivité ascendante ou descendane ?
   *    --> face_cell
   *    --> Par ascendance --> Trop dure je pense
   *   Dans pdm_mesh_adpation c'est en ascendant mais on part des vertex .
   *   Subtile on reconstruit dans ce sens --> vtx -> edge -> face -> cell
   *   Mais on construit les connectivités descendantes --> edge_vtx -> face_edge...
   */


  /*
   * Par descendance :
   *   --> Calcul du dcell_face --> Attention au signe !!!!!!!
   *   --> Quand on transpose la connectvity il faut d'aborder transformé le dface_cell en orienté
   *   --> This one : PDM_setup_connectivity_idx
   *
   *   PMD_multiblock_to_part avec dcell_face + cell_ln_to_gn concateante
   *     --> Donc un dcell_face re-repartie mais faux en face
   *     --> Si on est des ouf --> Reorder_block_hilbert
   *   --> Pour update le dcell_face
   *         --> mbtp avec ln_to_gn = face_to_keep
   *   --> On utilise le dface_jon + dface_join_opp
   *   Attention au desequilibre de charge possibleeeeeeeee
   *   Pour l'update des numero de face --> block_to_part avec un block variable de numero de faces old_to_new
   *     --> Si la face est supprimé strid = 0
   *
   */






}



/**
 *
 * \brief  Usage
 *
 */

static void
_usage(int exit_code)
{
  PDM_printf
    ("\n"
     "  Usage: \n\n"
     "  -n      <level>  Number of vertices on the cube side.\n\n"
     "  -l      <level>  Cube length.\n\n"
     "  -h               This message.\n\n");

  exit(exit_code);
}

/**
 *
 * \brief  Read arguments from the command line
 *
 * \param [in]      argc     Number of arguments
 * \param [in]      argv     Arguments
 * \param [inout]   n_vtx_seg  Number of vertices on the cube side
 * \param [inout]   length   Cube length
 * \param [inout]   n_part   Number of partitions par process
 * \param [inout]   post     Ensight outputs status
 * \param [inout]   method   Partitioner (1 ParMETIS, 2 Pt-Scotch)
 *
 */

static void
_read_args(int            argc,
           char         **argv,
           PDM_g_num_t  *n_vtx_seg,
           double        *length)
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {

    if (strcmp(argv[i], "-h") == 0)
      _usage(EXIT_SUCCESS);

    else if (strcmp(argv[i], "-n") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        long _n_vtx_seg = atol(argv[i]);
        *n_vtx_seg = (PDM_g_num_t) _n_vtx_seg;
      }
    }
    else if (strcmp(argv[i], "-l") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *length = atof(argv[i]);
    }
    else
      _usage(EXIT_FAILURE);
    i++;
  }
}


/**
 *
 * \brief  Main
 *
 */

int main(int argc, char *argv[])
{
  PDM_g_num_t        n_vtx_seg = 10;
  double             length    = 1.;
  int                n_zone    = 2;


  _read_args(argc,
             argv,
             &n_vtx_seg,
             &length);

  int i_rank;
  int n_rank;

  PDM_MPI_Init(&argc, &argv);
  PDM_MPI_Comm     comm = PDM_MPI_COMM_WORLD;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  /* Alloc for distributed mesh */
  int *dn_cell       = (int *) malloc(n_zone * sizeof(int));
  int *dn_face       = (int *) malloc(n_zone * sizeof(int));
  int *dn_vtx        = (int *) malloc(n_zone * sizeof(int));
  int *n_face_group  = (int *) malloc(n_zone * sizeof(int));
  int *dface_vtx_s   = (int *) malloc(n_zone * sizeof(int));
  int *dface_group_s = (int *) malloc(n_zone * sizeof(int));

  PDM_g_num_t  **dface_cell      = (PDM_g_num_t **) malloc(n_zone * sizeof(PDM_g_num_t *));
  int          **dface_vtx_idx   = (int         **) malloc(n_zone * sizeof(int         *));
  PDM_g_num_t  **dface_vtx       = (PDM_g_num_t **) malloc(n_zone * sizeof(PDM_g_num_t *));
  double       **dvtx_coord      = (double      **) malloc(n_zone * sizeof(double      *));
  int          **dface_group_idx = (int         **) malloc(n_zone * sizeof(int         *));
  PDM_g_num_t  **dface_group     = (PDM_g_num_t **) malloc(n_zone * sizeof(PDM_g_num_t *));
  int          **dface_bnd_idx   = (int         **) malloc(n_zone * sizeof(int         *));
  PDM_g_num_t  **dface_bnd       = (PDM_g_num_t **) malloc(n_zone * sizeof(PDM_g_num_t *));

  PDM_g_num_t  **extract_face_distribution = (PDM_g_num_t **) malloc( n_zone * sizeof(PDM_g_num_t *));
  PDM_g_num_t  **extract_vtx_distribution  = (PDM_g_num_t **) malloc( n_zone * sizeof(PDM_g_num_t *));
  int          **dextract_face_vtx_idx     = (int         **) malloc( n_zone * sizeof(int         *));
  PDM_g_num_t  **dextract_face_vtx         = (PDM_g_num_t **) malloc( n_zone * sizeof(PDM_g_num_t *));
  PDM_g_num_t  **dparent_face_g_num        = (PDM_g_num_t **) malloc( n_zone * sizeof(PDM_g_num_t *));
  PDM_g_num_t  **dparent_vtx_g_num         = (PDM_g_num_t **) malloc( n_zone * sizeof(PDM_g_num_t *));
  PDM_g_num_t  **pextract_old_to_new       = (PDM_g_num_t **) malloc( n_zone * sizeof(PDM_g_num_t *));
  double       **dextract_vtx_coord        = (double      **) malloc( n_zone * sizeof(double      *));

  int n_group_join = 2*(n_zone-1);
  int *group_join_to_zone_cur = (int *) malloc( n_group_join * sizeof(int));
  int *group_join_to_zone_opp = (int *) malloc( n_group_join * sizeof(int));
  int *group_join_to_join_opp = (int *) malloc( n_group_join * sizeof(int));

  int          *dface_join_idx  = (int         *) malloc((n_group_join + 1) * sizeof(int        ));
  PDM_g_num_t  *dface_join      = NULL; // A allouer propremet

  int tmp_i_zone = 0;
  for (int i_join = 0; i_join < n_group_join; i_join++) {
    if (i_join % 2 == 0) {
      group_join_to_join_opp[i_join] = i_join + 1;
      group_join_to_zone_opp[i_join] = tmp_i_zone + 1;
    } else {
      group_join_to_join_opp[i_join] = i_join - 1;
      group_join_to_zone_opp[i_join] = tmp_i_zone++;
    }
  }


  for(int i_group_join = 0; i_group_join < n_group_join+1; ++i_group_join) {
    dface_join_idx[i_group_join] = 0;
  }

  PDM_dcube_t **dcube = (PDM_dcube_t **) malloc(n_zone * sizeof(PDM_dcube_t *));
  int tmp_i_group_join = 0;
  for (int i_zone = 0; i_zone < n_zone; i_zone++) {

    dcube[i_zone] = PDM_dcube_gen_init(comm, n_vtx_seg, length, i_zone, 0., 0., PDM_OWNERSHIP_KEEP);
    PDM_dcube_gen_dim_get(dcube         [i_zone],
                          &n_face_group [i_zone],
                          &dn_cell      [i_zone],
                          &dn_face      [i_zone],
                          &dn_vtx       [i_zone],
                          &dface_vtx_s  [i_zone],
                          &dface_group_s[i_zone]);

    PDM_dcube_gen_data_get(dcube          [i_zone],
                          &dface_cell     [i_zone],
                          &dface_vtx_idx  [i_zone],
                          &dface_vtx      [i_zone],
                          &dvtx_coord     [i_zone],
                          &dface_group_idx[i_zone],
                          &dface_group    [i_zone]);

    /*
     * Les faces groups du dcube sont : zmin, zmax, xmin, xmax, ymin, ymax
     * Il faut les séparer en faces de bords et faces raccord, sachant que
     * les zones sont alignées selon X
     */
    int n_bnd = 4;
    int n_jn  = 2;
    if (i_zone == 0){
      n_bnd++;
      n_jn-- ;
    }
    if (i_zone == n_zone-1){
      n_bnd++;
      n_jn-- ;
    }

    // Join numbering (left to right, increasing i_zone)

    dface_bnd_idx [i_zone] = (int *) malloc((n_bnd        + 1) * sizeof(int));

    // First pass to count and allocate
    int i_bnd = 1;
    int i_jn  = tmp_i_group_join+1;
    dface_bnd_idx[i_zone][0]  = 0;
    // dface_join_idx[0] = 0;
    for (int igroup = 0; igroup < n_face_group[i_zone]; igroup++) {
      int copy_to_bnd = (igroup != 2 || (igroup == 2 && i_zone == 0)) && (igroup != 3 || (igroup == 3 && i_zone == n_zone-1));
      int group_size = dface_group_idx[i_zone][igroup+1] - dface_group_idx[i_zone][igroup];
      if (copy_to_bnd) { //Its a boundary
        dface_bnd_idx[i_zone][i_bnd++] = group_size;
      } else { //Its a join
        group_join_to_zone_cur[i_jn-1] = i_zone;
        dface_join_idx[i_jn++] = group_size;
      }
    }
    for (int i = 0; i < n_bnd; i++) {
      dface_bnd_idx[i_zone][i+1] = dface_bnd_idx[i_zone][i+1] + dface_bnd_idx[i_zone][i];
    }


    for (int i = tmp_i_group_join; i < i_jn-1; i++) {
      dface_join_idx[i+1] = dface_join_idx[i+1] + dface_join_idx[i];
    }

    /* A bit stupid but complicated to made it in ohter way for a clear test */
    dface_join = realloc(dface_join, dface_join_idx[i_jn-1] * sizeof(PDM_g_num_t));

    // Second pass to copy
    dface_bnd [i_zone] = (PDM_g_num_t *) malloc(dface_bnd_idx [i_zone][n_bnd        ] * sizeof(PDM_g_num_t));
    i_bnd = 0;
    i_jn  = tmp_i_group_join;
    for (int igroup = 0; igroup < n_face_group[i_zone]; igroup++) {
      int copy_to_bnd = (igroup != 2 || (igroup == 2 && i_zone == 0)) && (igroup != 3 || (igroup == 3 && i_zone == n_zone-1));
      if (copy_to_bnd){ //Its a boundary
        for (int i = dface_group_idx[i_zone][igroup]; i < dface_group_idx[i_zone][igroup+1]; i++) {
          dface_bnd[i_zone][i_bnd++] = dface_group[i_zone][i];
        }
      } else { //Its a join
        int k = 0;
        for (int i = dface_group_idx[i_zone][igroup]; i < dface_group_idx[i_zone][igroup+1]; i++) {
          dface_join[dface_join_idx[i_jn]+k++] = dface_group[i_zone][i];
        }
        i_jn++;
      }
    }

    /*
     *  Go to nexts join
     */
    tmp_i_group_join += n_jn;
  }

  /*
   * Setup dface_join_opp + group_id in current layout
   */
  PDM_g_num_t *dface_join_opp = NULL;
  _exchange_point_list(n_group_join,
                       group_join_to_join_opp,
                       dface_join_idx,
                       dface_join,
                       &dface_join_opp,
                       comm);

  log_trace("Global join data (%d)\n", n_group_join);
  PDM_log_trace_array_int(group_join_to_zone_cur, n_group_join, "group_join_to_zone_cur :: ");
  PDM_log_trace_array_int(group_join_to_join_opp, n_group_join, "group_join_to_join_opp :: ");
  PDM_log_trace_array_int(group_join_to_zone_opp, n_group_join, "group_join_to_zone_opp :: ");

  PDM_log_trace_array_int(dface_join_idx, n_group_join+1, "dface_join_idx :: ");
  PDM_log_trace_array_long(dface_join    , dface_join_idx[n_group_join], "dface_join     :: ");
  PDM_log_trace_array_long(dface_join_opp, dface_join_idx[n_group_join], "dface_join_opp :: ");

  // New version begins

  // Extract all the jn faces
  
  PDM_g_num_t *face_per_block_offset = (PDM_g_num_t *) malloc((n_zone+1) * sizeof(PDM_g_num_t));
  PDM_g_num_t *vtx_per_block_offset  = (PDM_g_num_t *) malloc((n_zone+1) * sizeof(PDM_g_num_t));
  face_per_block_offset[0] = 0;
  vtx_per_block_offset[0] = 0;
  PDM_MPI_Allreduce(dn_face, &face_per_block_offset[1], n_zone, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, comm);
  PDM_MPI_Allreduce(dn_vtx , &vtx_per_block_offset[1] , n_zone, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, comm);
  PDM_array_accumulate_gnum(face_per_block_offset, n_zone+1);
  PDM_array_accumulate_gnum(vtx_per_block_offset , n_zone+1);

  //int n_face_join = dface_join_idx[n_group_join]; //
  int n_face_join = 0; // Put faces only once
  for (int i_join = 0; i_join < n_group_join; i_join++) {
    if (i_join <= group_join_to_join_opp[i_join])
      n_face_join += 2*(dface_join_idx[i_join+1] - dface_join_idx[i_join]);
  }

  PDM_g_num_t *multi_gnum        = malloc(n_face_join * sizeof(PDM_g_num_t));

  int idx = 0;
  for (int i_join = 0; i_join < n_group_join; i_join++) {
    int i_join_opp = group_join_to_join_opp[i_join];
    int i_zone_cur = group_join_to_zone_cur[i_join];
    int i_zone_opp = group_join_to_zone_opp[i_join];

    log_trace("Treat jn %i\n", i_join);
    if (i_join <= i_join_opp) {
      for (int i_face_jn = dface_join_idx[i_join]; i_face_jn < dface_join_idx[i_join+1]; i_face_jn++) {
        multi_gnum[idx++] = dface_join[i_face_jn] + face_per_block_offset[i_zone_cur];
        multi_gnum[idx++] = dface_join_opp[i_face_jn] + face_per_block_offset[i_zone_opp];
      }
      /*for (int i_face_jn = dface_join_idx[i_join]; i_face_jn < dface_join_idx[i_join+1]; i_face_jn++) {*/
        /*multi_gnum[idx++] = dface_join_opp[i_face_jn] + face_per_block_offset[i_zone_opp];*/
      /*}*/
    }
  }
  assert (idx == n_face_join);

  PDM_log_trace_array_long(face_per_block_offset, n_zone+1, "face_per_block_offset :: ");
  PDM_log_trace_array_long(vtx_per_block_offset,  n_zone+1, "vtx_per_block_offset :: ");
  PDM_log_trace_array_long(multi_gnum, n_face_join, "multi_gnum :: ");

  PDM_g_num_t **all_face_distribution = malloc(n_zone * sizeof(PDM_g_num_t*));
  for (int i_zone = 0; i_zone < n_zone; i_zone++) {
    all_face_distribution[i_zone] = PDM_compute_entity_distribution(comm, dn_face[i_zone]);
  }
  PDM_multi_block_to_part_t *mptb = PDM_multi_block_to_part_create(face_per_block_offset,
                                                                   n_zone,
                                            (const PDM_g_num_t **) all_face_distribution,
                                            (const PDM_g_num_t **)&multi_gnum,
                                                                  &n_face_join,
                                                                   1,
                                                                   comm);
  //Prepare data to send : face -> vtx connectivity 
  int **face_vtx_n       = malloc(n_zone * sizeof(int*));
  PDM_g_num_t **face_vtx_shifted = malloc(n_zone * sizeof(int*));
  for (int i_zone = 0; i_zone < n_zone; i_zone++) {
    face_vtx_n[i_zone]       = malloc(dn_face[i_zone] * sizeof(int));
    face_vtx_shifted[i_zone] = malloc(dface_vtx_idx[i_zone][dn_face[i_zone]] * sizeof(PDM_g_num_t));

    for (int i_face = 0; i_face < dn_face[i_zone]; i_face++) {
      face_vtx_n[i_zone][i_face] = dface_vtx_idx[i_zone][i_face+1] - dface_vtx_idx[i_zone][i_face];
      for (int j = dface_vtx_idx[i_zone][i_face]; j < dface_vtx_idx[i_zone][i_face+1]; j++) {
        face_vtx_shifted[i_zone][j] = dface_vtx[i_zone][j] + vtx_per_block_offset[i_zone];
      }
    }
    //PDM_log_trace_array_long(face_vtx_shifted[i_zone], dface_vtx_idx[i_zone][dn_face[i_zone]], "face_vtx_shifted :: ");
  }
  //int **face_vtx = malloc(n_zone * sizeof(int*));



  int         **part_stride = NULL;
  PDM_g_num_t **part_data   = NULL;
  PDM_multi_block_to_part_exch2(mptb,
                                sizeof(PDM_g_num_t),
                                PDM_STRIDE_VAR,
                                face_vtx_n,
                 (void **)      face_vtx_shifted,
                               &part_stride,
                 (void ***)    &part_data);

  int         *face_vtx_both_n   = part_stride[0];
  int         *face_vtx_both_idx = PDM_array_new_idx_from_sizes_int(face_vtx_both_n, n_face_join);
  PDM_g_num_t *face_vtx_both     = part_data[0];

  free(part_data);
  free(part_stride);
  PDM_multi_block_to_part_free(mptb);

  int n_recv = 0;
  for (int i = 0; i < n_face_join; i++)
    n_recv += face_vtx_both_n[i];

  log_trace("Face vtx received after MBTP\n");
  PDM_log_trace_array_int (face_vtx_both_idx, n_face_join+1, "face_vtx_idx :: ");
  PDM_log_trace_array_long(face_vtx_both, n_recv, "face_vtx :: ");


  int *face_vtx_idx = malloc(((n_face_join / 2 )+1) * sizeof(int));
  PDM_g_num_t *face_vtx = malloc((n_recv / 2) * sizeof(PDM_g_num_t));
  PDM_g_num_t *face_vtx_opp = malloc((n_recv / 2 ) * sizeof(PDM_g_num_t));
  face_vtx_idx[0] = 0;

  idx = 0;
  int recv_idx = 0;
  int data_idx = 0;
  int recv_data_idx = 0;
  for (int i_join = 0; i_join < n_group_join; i_join++) {
    int i_join_opp = group_join_to_join_opp[i_join];
    int n_face_this_jn = dface_join_idx[i_join+1] - dface_join_idx[i_join];
    if (i_join <= i_join_opp) {
      for (int k = dface_join_idx[i_join]; k < dface_join_idx[i_join+1]; k++) {
        int n_edge = face_vtx_both_n[2*recv_idx]; //Jump opposite face stride
        // Update idx array
        face_vtx_idx[idx+1] = face_vtx_idx[idx] + n_edge;
        recv_idx++;
        idx++;
        // Updata data
        memcpy(&face_vtx    [data_idx], &face_vtx_both[recv_data_idx], n_edge * sizeof(PDM_g_num_t));
        memcpy(&face_vtx_opp[data_idx], &face_vtx_both[recv_data_idx+n_edge], n_edge * sizeof(PDM_g_num_t));
        recv_data_idx += 2*n_edge;
        data_idx += n_edge;
      }
    }
  }

  //Old version, when jn were not interlaced
  
  /*
  //We have jn_idx, prepare jn_data_idx to allow easier access
  int *jn_data_idx = malloc(n_group_join * sizeof(int));
  jn_data_idx[0] = 0;
  for (int i_join = 0; i_join < n_group_join; i_join++) {
    jn_data_idx[i_join+1] = 0;
    for (int k = dface_join_idx[i_join]; k < dface_join_idx[i_join+1]; k++) {
      jn_data_idx[i_join+1] += part_stride[0][k];
    }
    jn_data_idx[i_join+1] += jn_data_idx[i_join];
  }
  PDM_log_trace_array_int(dface_join_idx, n_group_join+1, "jn_idx :: ");
  PDM_log_trace_array_int(jn_data_idx, n_group_join+1, "jn_data_idx :: ");

  //Reorder all received face->edge connectivity to have pairs of faces
  int *face_edge_idx = malloc(((n_face_join / 2 )+1) * sizeof(int));
  PDM_g_num_t *face_edge = malloc((n_recv / 2) * sizeof(PDM_g_num_t));
  PDM_g_num_t *face_edge_opp = malloc((n_recv / 2 ) * sizeof(PDM_g_num_t));
  face_edge_idx[0] = 0;
  idx = 0;
  int data_idx = 0;
  for (int i_join = 0; i_join < n_group_join; i_join++) {
    int i_join_opp = group_join_to_join_opp[i_join];
    int n_face_this_jn = dface_join_idx[i_join+1] - dface_join_idx[i_join];

    log_trace("Post treat jn %i\n", i_join);
    if (i_join <= i_join_opp) {
      //Copy idx
      for (int k = dface_join_idx[i_join]; k < dface_join_idx[i_join+1]; k++) {
        face_edge_idx[idx+1] = face_edge_idx[idx] + part_stride[0][k];
        idx++;
      }
      //Copy data
      int data_size = jn_data_idx[i_join+1] - jn_data_idx[i_join];
      log_trace("data size is %d\n", data_size);
      assert (jn_data_idx[i_join_opp+1] - jn_data_idx[i_join_opp] == data_size);
      memcpy(&face_edge[data_idx], &part_data[0][jn_data_idx[i_join]], data_size * sizeof(PDM_g_num_t));
      memcpy(&face_edge_opp[data_idx], &part_data[0][jn_data_idx[i_join_opp]], data_size * sizeof(PDM_g_num_t));
      data_idx += data_size;
    }
  }
  */
  assert (idx == n_face_join / 2 );
  assert (data_idx == n_recv / 2);

  log_trace("Split face_vtx & face_vtx donor \n");
  PDM_log_trace_array_int(face_vtx_idx, idx+1, "face_vtx_idx :: ");
  PDM_log_trace_array_long(face_vtx, n_recv/2, "face_vtx :: ");
  PDM_log_trace_array_long(face_vtx_opp, n_recv/2, "face_vtx_opp :: ");

  //Now we have some pairs of faces (each pair appears only one) + face_vtx for this pairs
  
  //To get the number of unique vertex, we can do a part to block using face_vtx_both as gnum
  PDM_part_to_block_t *ptb = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                      PDM_PART_TO_BLOCK_POST_MERGE,
                                                      1.,
                                                     &face_vtx_both,
                                                      NULL,
                                                     &face_vtx_both_idx[n_face_join],
                                                      1,
                                                      comm);
  int n_vtx;
  int n_vtx_loc = PDM_part_to_block_n_elt_block_get(ptb);
  PDM_MPI_Allreduce(&n_vtx_loc, &n_vtx, 1, PDM_MPI_INT, PDM_MPI_SUM, comm);
  PDM_part_to_block_free(ptb);

  //Decompose into edges
  int ext_dn_face = n_face_join;
  int n_edge_elt = face_vtx_both_idx[ext_dn_face];
  int n_elmt_current = 0;
  int n_edge_current = 0;
  PDM_g_num_t *face_distri = PDM_compute_entity_distribution(comm, ext_dn_face);

  PDM_g_num_t* tmp_dface_edge         = (PDM_g_num_t *) malloc(n_edge_elt       * sizeof(PDM_g_num_t) );
  int*         tmp_parent_elmt_pos    = (int         *) malloc(n_edge_elt       * sizeof(int        ) );
  int*         tmp_dface_edge_vtx_idx = (int         *) malloc((n_edge_elt + 1) * sizeof(int        ) );
  PDM_g_num_t* tmp_dface_edge_vtx     = (PDM_g_num_t *) malloc(2*n_edge_elt     * sizeof(PDM_g_num_t) );

  tmp_dface_edge_vtx_idx[0] = 0;
  PDM_poly2d_decomposes_edges(ext_dn_face,
                              &n_elmt_current,
                              &n_edge_current,
                              face_distri[i_rank],
                              -1,
                              face_vtx_both,
                              face_vtx_both_idx,
                              tmp_dface_edge_vtx_idx, //Numéro de sommet des edges
                              tmp_dface_edge_vtx,
                              tmp_dface_edge,  // Numéro des edges pour une face
                              NULL,
                              NULL,
                              tmp_parent_elmt_pos);
  assert(n_edge_current == n_edge_elt);

  /*PDM_log_trace_array_long(tmp_dface_edge, n_edge_elt, "dface_edge :: ");*/
  /*PDM_log_trace_connectivity_long(tmp_dface_edge_vtx_idx, tmp_dface_edge_vtx, n_edge_elt, "dface_edge :: ");*/

  PDM_g_num_t dn_edge;
  PDM_g_num_t *dedge_distrib  = NULL;
  int         *dedge_vtx_idx  = NULL;
  PDM_g_num_t *dedge_vtx      = NULL;
  int         *dedge_face_idx = NULL;
  PDM_g_num_t *dedge_face     = NULL;
  PDM_generate_entitiy_connectivity_raw(comm,
                                        n_vtx, // n vtx tot
                                        n_edge_elt,
                                        tmp_dface_edge,
                                        tmp_dface_edge_vtx_idx,
                                        tmp_dface_edge_vtx,
                                        &dn_edge,
                                        &dedge_distrib ,
                                        &dedge_vtx_idx ,
                                        &dedge_vtx     ,
                                        &dedge_face_idx,
                                        &dedge_face    );
  free(tmp_parent_elmt_pos    );

  log_trace("Edges rebuild\n");
  //PDM_log_trace_array_long(dedge_vtx, 2*dn_edge, "dedge_vtx ::");
  //PDM_log_trace_connectivity_long(dedge_face_idx, dedge_face, dn_edge, "dedge_face ::");

  // Count the number of internal & external edges
  int dn_internal_edge = 0;
  for(int i_edge = 0; i_edge < dn_edge; ++i_edge) {
    int n_face_this_edge = dedge_face_idx[i_edge+1] - dedge_face_idx[i_edge];
    dn_internal_edge += (int) (n_face_this_edge > 1);
  }
  int dn_external_edge = dn_edge - dn_internal_edge;
  log_trace("dn internal edges is %i \n", dn_internal_edge);
  log_trace("dn external edges is %i \n", dn_external_edge);



  // Transport data to edges from face to build the key


  //Prepare numbering
  PDM_g_num_t *dedge_face_abs = (PDM_g_num_t *) malloc(dedge_face_idx[dn_edge] * sizeof(PDM_g_num_t));
  int         *dedge_face_sgn = (int         *) malloc(dedge_face_idx[dn_edge] * sizeof(int        ));
  for(int i = 0; i < dedge_face_idx[dn_edge]; ++i) {
    dedge_face_abs[i] = PDM_ABS (dedge_face[i]);
    dedge_face_sgn[i] = PDM_SIGN(dedge_face[i]);
  }

  //Prepare data on multi_gnum (use same ordering than before) to transfert to edges
  int *dextract_face_group_idNEW = malloc(ext_dn_face*sizeof(int));
  PDM_g_num_t *dextract_face_joinNEW     = malloc(ext_dn_face*sizeof(PDM_g_num_t));
  PDM_g_num_t *dextract_face_join_oppNEW = malloc(ext_dn_face*sizeof(PDM_g_num_t));
  idx = 0;
  for (int i_join = 0; i_join < n_group_join; i_join++) {
    int i_join_opp = group_join_to_join_opp[i_join];
    int i_zone_cur = group_join_to_zone_cur[i_join];
    int i_zone_opp = group_join_to_zone_opp[i_join];

    if (i_join <= i_join_opp) {
      for (int i_face_jn = dface_join_idx[i_join]; i_face_jn < dface_join_idx[i_join+1]; i_face_jn++) {
        //Face data
        dextract_face_joinNEW    [idx] = dface_join[i_face_jn]; //Here unshifted
        dextract_face_join_oppNEW[idx] = dface_join_opp[i_face_jn];
        dextract_face_group_idNEW[idx++] = i_join;
        //Opp face data
        dextract_face_joinNEW    [idx] = dface_join_opp[i_face_jn];
        dextract_face_join_oppNEW[idx] = dface_join[i_face_jn];
        dextract_face_group_idNEW[idx++] = i_join_opp;
      }
    }
  }
  assert (idx == n_face_join);


  PDM_g_num_t* dedge_face_join       = malloc(dedge_face_idx[dn_edge] * sizeof(PDM_g_num_t));
  PDM_g_num_t* dedge_face_join_opp   = malloc(dedge_face_idx[dn_edge] * sizeof(PDM_g_num_t));
  int        * dedge_face_group_id   = malloc(dedge_face_idx[dn_edge] * sizeof(int        ));
  PDM_block_to_part_t *btp = PDM_block_to_part_create(face_distri,
                               (const PDM_g_num_t **) &dedge_face_abs,
                                                      &dedge_face_idx[dn_edge],
                                                      1,
                                                      comm);
  int cst_stride = 1;
  PDM_block_to_part_exch(btp,
                         sizeof(int),
                         PDM_STRIDE_CST,
                         &cst_stride,
                (void *) dextract_face_group_idNEW,
                         NULL,
             (void ** ) &dedge_face_group_id);
  PDM_block_to_part_exch(btp,
                         sizeof(PDM_g_num_t),
                         PDM_STRIDE_CST,
                         &cst_stride,
                (void *) dextract_face_joinNEW,
                         NULL,
             (void ** ) &dedge_face_join);

  PDM_block_to_part_exch(btp,
                         sizeof(PDM_g_num_t),
                         PDM_STRIDE_CST,
                         &cst_stride,
                (void *) dextract_face_join_oppNEW,
                         NULL,
             (void ** ) &dedge_face_join_opp);

  PDM_block_to_part_free(btp);

  log_trace("Transport data on edges\n");
  PDM_log_trace_array_int (dedge_face_group_id, dedge_face_idx[dn_edge], "dedge_face_group_id ::");
  PDM_log_trace_array_long(dedge_face_join    , dedge_face_idx[dn_edge], "dedge_face_join     ::");
  PDM_log_trace_array_long(dedge_face_join_opp, dedge_face_idx[dn_edge], "dedge_face_join_opp ::");

  // Begin hash table
  PDM_g_num_t *key_ln_to_gn = malloc(dn_internal_edge * sizeof(PDM_g_num_t)); 
  int         *stride_one   = malloc(dn_internal_edge * sizeof(int        ));
  int         *stride_two   = malloc(dn_internal_edge * sizeof(int        ));
  int         *stride_four  = malloc(dn_internal_edge * sizeof(int        ));

  int         *zone_id              = malloc(  dn_internal_edge * sizeof(int        ));
  PDM_g_num_t *data_send_connect    = malloc(4*dn_internal_edge * sizeof(PDM_g_num_t));
  PDM_g_num_t *data_send_edge_g_num = malloc(  dn_internal_edge * sizeof(PDM_g_num_t));
  PDM_g_num_t *data_send_group      = malloc(4*dn_internal_edge * sizeof(PDM_g_num_t));
  PDM_g_num_t *data_send_sens       = malloc(2*dn_internal_edge * sizeof(PDM_g_num_t));
  PDM_g_num_t *data_send_face_g_num = malloc(2*dn_internal_edge * sizeof(PDM_g_num_t));

  int i_int_edge = 0;
  int idx_write2 = 0;
  int idx_write4 = 0;
  for(int i_edge = 0; i_edge < dn_edge; ++i_edge) {

    int n_face_this_edge = dedge_face_idx[i_edge+1] - dedge_face_idx[i_edge];
    if (n_face_this_edge == 1) {
      continue;
    }
    assert (n_face_this_edge == 2);

    stride_one [i_int_edge] = 1;
    stride_two [i_int_edge] = 2;
    stride_four[i_int_edge] = 4;
    //Retrive zone id using group id of any of two faces
    data_send_edge_g_num[i_int_edge] = dedge_distrib[i_rank] + i_edge + 1;
    int group_id = dedge_face_group_id[dedge_face_idx[i_edge]];
    zone_id[i_int_edge] = group_join_to_zone_cur[group_id];

    int key = 0;
    for(int j = dedge_face_idx[i_edge]; j < dedge_face_idx[i_edge+1]; ++j) { //Do it for the two faces data
      key += (dedge_face_join[j] + dedge_face_join_opp[j]);

      data_send_face_g_num[idx_write2]   = dedge_face[j];
      idx_write2++;
      //data_send_sens      [idx_write2++] = dedge_face_group_sens[j];

      data_send_connect[idx_write4] = dedge_face_join    [j];
      data_send_group[idx_write4++] = dedge_face_group_id[j];

      data_send_connect[idx_write4] = dedge_face_join_opp[j];
      data_send_group[idx_write4++] = group_join_to_join_opp[dedge_face_group_id[j]]; // This is group join opp
    }
    key_ln_to_gn[i_int_edge] = key;

    i_int_edge++;
  }
  assert(idx_write2 == 2*dn_internal_edge);
  assert(idx_write4 == 4*dn_internal_edge);

  PDM_log_trace_array_long(key_ln_to_gn, dn_internal_edge, "key_ln_to_gn :: ");

  //Attention, pb d'équilibrage car les clés sont réparties vers la fin ... Un proc risque
  // de se retrouver avec tt les clés
                       ptb = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                      PDM_PART_TO_BLOCK_POST_MERGE,
                                                      1.,
                                    (PDM_g_num_t **) &key_ln_to_gn,
                                                      NULL,
                                                     &dn_internal_edge,
                                                      1,
                                                      comm);
  // Get protocol data
  int blk_size = PDM_part_to_block_n_elt_block_get(ptb);
  // This one will contain 2 for non conflicting hash key, more otherwise
  int *gnum_n_occurences = PDM_part_to_block_block_gnum_count_get(ptb);

  int gnum_n_occurences_tot = 0;
  for (int k = 0; k < blk_size; k++) {
    gnum_n_occurences_tot += gnum_n_occurences[k];
  }
  int *blk_entity_idx    = (int *) malloc((gnum_n_occurences_tot + 1)*sizeof(int));
  int *blk_data_face_idx = (int *) malloc((gnum_n_occurences_tot + 1)*sizeof(int));
  for (int k = 0; k < gnum_n_occurences_tot + 1; k++) {
    blk_entity_idx[k]    = 4*k;
    blk_data_face_idx[k] = 2*k;
  }
  
  // Exch data : zone_id, edge gnum, face & facedonor id, face & facedonor group id, face_sens and faces_gnum

  int *unused_recv_stride = NULL;
  int *blk_zone_id   = NULL;
  int exch_size = PDM_part_to_block_exch(ptb,
                                         sizeof(int),
                                         PDM_STRIDE_VAR,
                                         -1,
                               (int **)  &stride_one,
                               (void **) &zone_id,
                                         &unused_recv_stride,
                               (void **) &blk_zone_id);
  free(unused_recv_stride); // Same as gnum_n_occurences 
  assert (exch_size == gnum_n_occurences_tot);

  PDM_g_num_t* blk_edge_g_num           = NULL;
  exch_size = PDM_part_to_block_exch(ptb,
                                     sizeof(PDM_g_num_t),
                                     PDM_STRIDE_VAR,
                                     -1,
                           (int **)  &stride_one,
                           (void **) &data_send_edge_g_num,
                                     &unused_recv_stride,
                           (void **) &blk_edge_g_num);
  free(unused_recv_stride); // Same as gnum_n_occurences 
  assert (exch_size == gnum_n_occurences_tot);

  /*
  PDM_g_num_t* blk_data_sens   = NULL;
  exch_size = PDM_part_to_block_exch(ptb,
                                     sizeof(int),
                                     PDM_STRIDE_VAR,
                                     -1,
                           (int **)  &stride_two,
                           (void **) &data_send_sens,
                                     &unused_recv_stride,
                           (void **) &blk_data_sens);
  free(unused_recv_stride); // Same as 2*gnum_n_occurences 
  assert (exch_size == 2*gnum_n_occurences_tot);
  */

  PDM_g_num_t* blk_data_face_g_num   = NULL;
  exch_size = PDM_part_to_block_exch(ptb,
                                     sizeof(PDM_g_num_t),
                                     PDM_STRIDE_VAR,
                                     -1,
                           (int **)  &stride_two,
                           (void **) &data_send_face_g_num,
                                     &unused_recv_stride,
                           (void **) &blk_data_face_g_num);
  free(unused_recv_stride); // Same as 2*gnum_n_occurences
  assert (exch_size == 2*gnum_n_occurences_tot);

  PDM_g_num_t* blk_data_connect   = NULL;
  exch_size = PDM_part_to_block_exch(ptb,
                                     sizeof(PDM_g_num_t),
                                     PDM_STRIDE_VAR,
                                     -1,
                          (int **)  &stride_four,
                          (void **) &data_send_connect,
                                    &unused_recv_stride,
                          (void **) &blk_data_connect);
  free(unused_recv_stride); // Same as 4*gnum_n_occurences 
  assert (exch_size == 4*gnum_n_occurences_tot);

  PDM_g_num_t* blk_data_group   = NULL;
  exch_size = PDM_part_to_block_exch(ptb,
                                     sizeof(int),
                                     PDM_STRIDE_VAR,
                                     -1,
                           (int **)  &stride_four,
                           (void **) &data_send_group,
                                     &unused_recv_stride,
                           (void **) &blk_data_group);
  free(unused_recv_stride); // Same as 4*gnum_n_occurences 
  assert (exch_size == 4*gnum_n_occurences_tot);


  PDM_log_trace_array_int(gnum_n_occurences   , blk_size               , "gnum_n_occurences   :: ");
  PDM_log_trace_array_int(blk_zone_id         , gnum_n_occurences_tot  , "blk_zone_id         :: ");
  PDM_log_trace_array_long(blk_edge_g_num     , gnum_n_occurences_tot  , "blk_edge_g_num      :: ");
  PDM_log_trace_array_long(blk_data_face_g_num, 2*gnum_n_occurences_tot, "blk_data_face_g_num :: ");
  //PDM_log_trace_array_long(blk_data_sens      , 2*gnum_n_occurences_tot, "blk_data_sens       :: ");
  PDM_log_trace_array_long(blk_data_connect   , 4*gnum_n_occurences_tot, "blk_data_connect    :: ");
  PDM_log_trace_array_long(blk_data_group     , 4*gnum_n_occurences_tot, "blk_data_group      :: ");


  free(key_ln_to_gn        );
  free(data_send_connect   );
  free(data_send_group     );
  free(data_send_edge_g_num);
  free(data_send_sens      );
  free(stride_one          );
  free(stride_two          );
  free(stride_four         );
  free(data_send_face_g_num);
  free(zone_id             );

  // Post treatemement : resolve conflicting keys
  int n_max_entity_per_key = 0;
  for(int i = 0; i < blk_size; ++i) {
    n_max_entity_per_key = PDM_MAX(gnum_n_occurences   [i], n_max_entity_per_key);
  }
  int n_max_connec = 4*n_max_entity_per_key;


  //Number of edge to treat on each zone (Usefull ?)
  int* zone_id_n = PDM_array_zeros_int(n_zone);
  for (int i_edge = 0; i_edge < gnum_n_occurences_tot; i_edge++) {
    zone_id_n[blk_zone_id[i_edge]]++;
  }

  PDM_log_trace_array_int(zone_id_n        , n_zone                , "zone_id_n         :: ");
  log_trace("n_max_entity_per_key = %i \n", n_max_entity_per_key);
  log_trace("n_max_connec         = %i \n", n_max_connec);
  int*          already_treat    = (int         *) malloc(  n_max_entity_per_key * sizeof(int          ) );
  int*          same_entity_idx  = (int         *) malloc(  n_max_entity_per_key * sizeof(int          ) );
  //int*          sens_entity      = (int         *) malloc(  n_max_entity_per_key * sizeof(int          ) );
  PDM_g_num_t *results_edge     = malloc(gnum_n_occurences_tot * sizeof(PDM_g_num_t));
  PDM_g_num_t *results_edge_opp = malloc(gnum_n_occurences_tot * sizeof(PDM_g_num_t));

  /* Reset to fill */
  PDM_array_reset_int(zone_id_n, n_zone, 0);

  idx  = 0;
  int idx_w = 0;
  for(int i_key = 0; i_key < blk_size; ++i_key) {

    int n_matching_edge = gnum_n_occurences[i_key];

    log_trace(" i_key = %i | n_matching_edge = %i \n", i_key, n_matching_edge);

    /* Reset */
    PDM_array_reset_int(already_treat, n_matching_edge, -1);

    /* Loop over all entitys in conflict and sort all */
    for(int i_entity = 0; i_entity < n_matching_edge; ++i_entity) {

      //Each internal edge comes with 2 faces -> 4 data
      int beg    = 4*(idx+i_entity);

      // Caution inplace sort !!!!
      PDM_sort_long(&blk_data_connect[beg], NULL, 4);
      PDM_sort_int (&blk_data_group  [beg], NULL, 4);

      if(0 == 1) {
        PDM_log_trace_array_long(&blk_data_connect[beg], 4, "blk_data_connect (sort) :: ");
        PDM_log_trace_array_int (&blk_data_group  [beg], 4, "blk_data_group   (sort) :: ");
      }
    }

    /*
     *  Identify pair or invalid other
     */
    for(int i_entity = 0; i_entity < n_matching_edge; ++i_entity) {
      int beg1    = 4*(idx+i_entity);

      int n_same         = 0;
      int i_entity2_same = -1;

      if(already_treat[i_entity] == 1) {
        continue;
      }

      for(int i_entity2 = i_entity+1; i_entity2 < n_matching_edge; ++i_entity2) {

        if(already_treat[i_entity2] == -1) {
          int beg2    = 4*(idx+i_entity2);

          if(!PDM_array_are_equal_int(&blk_data_group[beg1], &blk_data_group[beg2], 4)) {
            continue;
          }

          if(!PDM_array_are_equal_int(&blk_data_connect[beg1], &blk_data_connect[beg2], 4)) {
            continue;
          }

          already_treat[i_entity2] = 1;
          same_entity_idx[n_same++] = i_entity2;
          i_entity2_same = i_entity2;

        }
      } /* End for i_entity2 */
      assert(n_same == 1);

      //log_trace("i_entity = %i | i_entity2_same = %i | n_same = %i \n", i_entity, i_entity2_same, n_same);

      /*
       * Renvoie des resultats :
       *    - Par edge --> Il faut donc également le numero de zones
       */
      int edge_idx     = (idx+i_entity);
      int edge_idx_opp = (idx+i_entity2_same);

      int i_zone_cur = blk_zone_id[edge_idx];
      int i_zone_opp = blk_zone_id[edge_idx_opp];

      // Set data for edge
      results_edge    [idx_w] = blk_edge_g_num[edge_idx];
      results_edge_opp[idx_w++] = blk_edge_g_num[edge_idx_opp];
      //results_edge    [i_zone_cur][zone_id_n[i_zone_cur]  ] = blk_edge_g_num[edge_idx];
      //results_edge_opp[i_zone_cur][zone_id_n[i_zone_cur]  ] = blk_edge_g_num[edge_idx_opp];
      //results_zone_opp[i_zone_cur][zone_id_n[i_zone_cur]++] = blk_zone_id[edge_idx_opp];

      // Set data for opposite edge
      results_edge    [idx_w] = blk_edge_g_num[edge_idx_opp];
      results_edge_opp[idx_w++] = blk_edge_g_num[edge_idx];
      //results_edge    [i_zone_opp][zone_id_n[i_zone_opp]  ] = blk_edge_g_num[edge_idx_opp];
      //results_edge_opp[i_zone_opp][zone_id_n[i_zone_opp]  ] = blk_edge_g_num[edge_idx];
      //results_zone_opp[i_zone_opp][zone_id_n[i_zone_opp]++] = blk_zone_id[edge_idx];

      // Renvoi de tout les edges candidats à travers les faces ???
      already_treat[i_entity] = 1;
    }

    idx  += n_matching_edge;
  }
  assert (idx_w == gnum_n_occurences_tot);

  PDM_part_to_block_free(ptb); // Needed for gnum_n_occurences

  log_trace("Conflict resolved, gnum are\n");
  PDM_log_trace_array_long(results_edge, gnum_n_occurences_tot, "edge gnum ::");
  PDM_log_trace_array_long(results_edge_opp, gnum_n_occurences_tot, "edge gnum opp::");

  free(zone_id_n);


  // Send back result on edge distribution (we were on key distribution)
                       ptb = PDM_part_to_block_create2(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                       PDM_PART_TO_BLOCK_POST_NOTHING,
                                                       1.,
                                                      &results_edge,
                                                       dedge_distrib,
                                                      &gnum_n_occurences_tot,
                                                       1,
                                                       comm);

  assert (PDM_part_to_block_n_elt_block_get(ptb) == dn_internal_edge);
  PDM_g_num_t *dedge_gnum     = malloc(PDM_part_to_block_n_elt_block_get(ptb) * sizeof(PDM_g_num_t));
  PDM_g_num_t *dedge_gnum_opp = NULL;

  PDM_g_num_t *dedge_gnum_tmp = PDM_part_to_block_block_gnum_get(ptb);
  memcpy(dedge_gnum, dedge_gnum_tmp, PDM_part_to_block_n_elt_block_get(ptb) * sizeof(PDM_g_num_t));

  PDM_part_to_block_exch(ptb,
                        sizeof(PDM_g_num_t),
                        PDM_STRIDE_CST,
                        1,
                        NULL,
             (void **) &(results_edge_opp),
                        NULL,
              (void **) &dedge_gnum_opp);
  PDM_part_to_block_free(ptb);
  
  // Rebuild idx array (no data on external edges)
  /*int *dedge_gnum_idx = (int *) malloc ((dn_edge + 1) * sizeof(int)); //This one will be useless probably*/
  /*int count = 0;*/
  /*dedge_gnum_idx[0] = 0;*/
  /*for (int i = 0; i < dn_edge[i_zone]; i++) {*/
    /*if (i + dedge_distrib[i_zone][i_rank] + 1 == dedge_gnum[count]) {*/
      /*count++;*/
    /*}*/
    /*dedge_gnum_idx[i+1] = count;*/
  /*}*/
  log_trace("Internal edge matches after conflict resolution \n");
  PDM_log_trace_array_long(dedge_gnum, dn_internal_edge, "dedge gnum :: ");
  PDM_log_trace_array_long(dedge_gnum_opp, dn_internal_edge, "dedge gnum_opp :: ");


  free(dextract_face_joinNEW);
  free(dextract_face_join_oppNEW);
  free(dextract_face_group_idNEW);

  free(dedge_face_abs);
  free(dedge_face_sgn);

  free(dedge_face_join);
  free(dedge_face_join_opp);
  free(dedge_face_group_id);


  free(face_distri);

  free(face_vtx_both_n);
  free(face_vtx_both_idx);
  free(face_vtx_both);

  free(face_per_block_offset);
  free(vtx_per_block_offset);
  free(multi_gnum);

  for (int i_zone = 0; i_zone < n_zone; i_zone++) {
    free(all_face_distribution[i_zone]);
    free(face_vtx_n[i_zone]);
    free(face_vtx_shifted[i_zone]);
  }
  free(all_face_distribution);
  free(face_vtx_n);
  free(face_vtx_shifted);

  free(face_vtx_idx);
  free(face_vtx);
  free(face_vtx_opp);

  free(dedge_distrib);
  free(dedge_vtx_idx);
  free(dedge_vtx);
  free(dedge_face_idx);
  free(dedge_face);


  exit(0);
  /*
   * Algorithm begin - Extract faces
   */
  int         **dextract_face_group_id   = malloc(n_zone * sizeof(int         *));
  int         **dextract_face_group_sens = malloc(n_zone * sizeof(int         *));
  PDM_g_num_t **dextract_face_join       = malloc(n_zone * sizeof(PDM_g_num_t *));
  PDM_g_num_t **dextract_face_join_opp   = malloc(n_zone * sizeof(PDM_g_num_t *));

  tmp_i_group_join = 0;
  for (int i_zone = 0; i_zone < n_zone; i_zone++) {
    int n_jn  = 2;
    if (i_zone == 0       ){n_jn--;}
    if (i_zone == n_zone-1){n_jn--;}

    /*
     *  Now we have all joins create we need to extract them
     */
    PDM_g_num_t* face_distribution = PDM_compute_entity_distribution(comm, dn_face[i_zone]);
    PDM_g_num_t* vtx_distribution  = PDM_compute_entity_distribution(comm, dn_vtx [i_zone]);

    /*
     *  By construction of actual join, each join of same zone are concatenate
     */
    int dn_l_face_join = dface_join_idx[tmp_i_group_join+n_jn] - dface_join_idx[tmp_i_group_join];
    PDM_g_num_t* l_dface_join = &dface_join[dface_join_idx[tmp_i_group_join]];

    PDM_dconnectivity_to_extract_dconnectivity(comm,
                                               dn_l_face_join,
                                               l_dface_join,
                                               face_distribution,
                                               dface_vtx_idx[i_zone],
                                               dface_vtx[i_zone],
                                               &extract_face_distribution[i_zone],
                                               &extract_vtx_distribution[i_zone],
                                               &dextract_face_vtx_idx[i_zone],
                                               &dextract_face_vtx[i_zone],
                                               &dparent_face_g_num[i_zone],
                                               &dparent_vtx_g_num[i_zone],
                                               &pextract_old_to_new[i_zone]);

    /*
     * Equilibrate the pextract_old_to_new
     */
    PDM_g_num_t* dface_group_init_distrib = PDM_compute_entity_distribution(comm, dn_l_face_join);
    int dn_extract_face = extract_face_distribution[i_zone][i_rank+1] - extract_face_distribution[i_zone][i_rank];
    PDM_g_num_t* extract_face_ln_to_gn = malloc(dn_extract_face * sizeof(PDM_g_num_t));
    for(int i = 0; i < dn_extract_face; ++i) {
      extract_face_ln_to_gn[i] = extract_face_distribution[i_zone][i_rank] + i + 1;
    }

    PDM_g_num_t** tmp_dextract_old_to_new = NULL;
    PDM_part_dfield_to_pfield(comm,
                              1,
                              sizeof(PDM_g_num_t),
                              dface_group_init_distrib,
      (unsigned char    *)    pextract_old_to_new[i_zone],
                              &dn_extract_face,
      (const PDM_g_num_t **)  &extract_face_ln_to_gn,
      (unsigned char ***)     &tmp_dextract_old_to_new);
    PDM_g_num_t *dextract_old_to_new = tmp_dextract_old_to_new[0];
    free(tmp_dextract_old_to_new);

    int dn_extract_vtx  = extract_vtx_distribution[i_zone][i_rank+1] - extract_vtx_distribution[i_zone][i_rank];

    log_trace("Extracted mesh for zone %i\n", i_zone);
    if (1 == 1) {
      PDM_log_trace_array_long(pextract_old_to_new[i_zone], dface_group_init_distrib[i_rank+1] - dface_group_init_distrib[i_rank], "pextract_old2new :: ");
      PDM_log_trace_array_long(dparent_face_g_num[i_zone], dn_extract_face, "dparent_face_g_num :: ");
      PDM_log_trace_array_long(dparent_vtx_g_num[i_zone], dn_extract_vtx, "dparent_vtx_g_num :: ");
      PDM_log_trace_array_long(dextract_old_to_new, dn_extract_face, "dextract_old_to_new :: ");
    }

    double** tmp_dextract_vtx_coord = NULL;
    PDM_part_dcoordinates_to_pcoordinates(comm,
                                          1,
                                          vtx_distribution,
                                          dvtx_coord[i_zone],
                                          &dn_extract_vtx,
                   (const PDM_g_num_t **) &dparent_vtx_g_num[i_zone],
                                          &tmp_dextract_vtx_coord);

    /*
     *  Pour chaque face extraite :
     *     - Echange du numero de face/face_opposé (Attention on est obligé de le faire car PDM_dconnectivity_to_extract_dconnectivity reordonne tout)
     *     - Exchange i_group also
     */
    int* dface_group_tag  = malloc(dn_l_face_join * sizeof(int));
    int* dface_group_sens = malloc(dn_l_face_join * sizeof(int)); // 1 = Exterior / -1 = Interior
    int k = 0;
    for(int i_group = tmp_i_group_join; i_group < tmp_i_group_join+n_jn; ++i_group) {
      for(int i_face = dface_join_idx[i_group]; i_face < dface_join_idx[i_group+1]; ++i_face) {
        dface_group_tag [k  ] = i_group;
        dface_group_sens[k++] = 1;
      }
    }

    if(0 == 1) {
      log_trace("dn_extract_face = %i \n", dn_extract_face);
      log_trace("dn_l_face_join  = %i \n", dn_l_face_join);
    }

    /*
     * Exchange tag
     */
    int** tmp_dextract_face_group_id = NULL;
    PDM_part_dfield_to_pfield(comm,
                              1,
                              sizeof(int),
                              dface_group_init_distrib,
      (unsigned char    *)    dface_group_tag,
                              &dn_extract_face,
      (const PDM_g_num_t **)  &dextract_old_to_new,
      (unsigned char ***)     &tmp_dextract_face_group_id);
    dextract_face_group_id[i_zone] = tmp_dextract_face_group_id[0];
    free(tmp_dextract_face_group_id);

    /*
     * Exchange sens
     */
    int** tmp_dextract_face_group_sens = NULL;
    PDM_part_dfield_to_pfield(comm,
                              1,
                              sizeof(int),
                              dface_group_init_distrib,
      (unsigned char    *)    dface_group_sens,
                              &dn_extract_face,
      (const PDM_g_num_t **)  &dextract_old_to_new,
      (unsigned char ***)     &tmp_dextract_face_group_sens);
    dextract_face_group_sens[i_zone] = tmp_dextract_face_group_sens[0];
    free(tmp_dextract_face_group_sens);

    /*
     *  Exchange dface_join
     */
    PDM_g_num_t* sub_dface_join = &dface_join[dface_join_idx[tmp_i_group_join]];
    int** tmp_dextract_face_join = NULL;
    PDM_part_dfield_to_pfield(comm,
                              1,
                              sizeof(PDM_g_num_t),
                              dface_group_init_distrib,
      (unsigned char    *)    sub_dface_join,
                              &dn_extract_face,
      (const PDM_g_num_t **)  &dextract_old_to_new,
      (unsigned char ***)     &tmp_dextract_face_join);
    dextract_face_join[i_zone] = tmp_dextract_face_join[0];
    free(tmp_dextract_face_join);

    /*
     *  Exchange dface_join_opp
     */
    PDM_g_num_t* sub_dface_join_opp = &dface_join_opp[dface_join_idx[tmp_i_group_join]];
    int** tmp_dextract_face_join_opp = NULL;
    PDM_part_dfield_to_pfield(comm,
                              1,
                              sizeof(PDM_g_num_t),
                              dface_group_init_distrib,
      (unsigned char    *)    sub_dface_join_opp,
                              &dn_extract_face,
      (const PDM_g_num_t **)  &dextract_old_to_new,
      (unsigned char ***)     &tmp_dextract_face_join_opp);
    dextract_face_join_opp[i_zone] = tmp_dextract_face_join_opp[0];
    free(tmp_dextract_face_join_opp);


    free(extract_face_ln_to_gn);

    dextract_vtx_coord[i_zone] = tmp_dextract_vtx_coord[0];
    free(tmp_dextract_vtx_coord);

    free(face_distribution);
    free(vtx_distribution);
    free(dface_group_init_distrib);
    free(dface_group_tag);
    free(dface_group_sens);
    free(dextract_old_to_new); // Maybe to keep ??

    /*
     * To keep
     */
    if(1 == 1) {
      PDM_log_trace_array_long(dextract_face_group_id  [i_zone], dn_extract_face, "dextract_face_group_id   :: ");
      PDM_log_trace_array_long(dextract_face_group_sens[i_zone], dn_extract_face, "dextract_face_group_sens :: ");
      PDM_log_trace_array_long(dextract_face_join      [i_zone], dn_extract_face, "dextract_face_join       :: ");
      PDM_log_trace_array_long(dextract_face_join_opp  [i_zone], dn_extract_face, "dextract_face_join_opp   :: ");
    }

    /*
     *  Go to nexts join
     */
    tmp_i_group_join += n_jn;

  }

  /*
   * Real algorithm
   */
  _deduce_descending_join(n_zone,
                          n_group_join,
                          group_join_to_zone_cur,
                          group_join_to_zone_opp,
                          group_join_to_join_opp,
                          dface_join_idx,
                          dface_join,
                          dface_join_opp,
                          dface_vtx_idx,
                          dface_vtx,
                          extract_face_distribution,
                          extract_vtx_distribution,
                          dextract_face_vtx_idx,
                          dextract_face_vtx,
                          dextract_face_group_id,
                          dextract_face_group_sens,
                          dextract_face_join,
                          dextract_face_join_opp,
                          dparent_face_g_num,
                          dparent_vtx_g_num,
                          pextract_old_to_new,
                          dextract_vtx_coord,
                          comm);

  if(0 == 1) {
    merge_zone(n_zone,
               n_group_join,
               group_join_to_zone_cur,
               group_join_to_zone_opp,
               group_join_to_join_opp,
               dface_join_idx,
               dface_join,
               dface_join_opp,
               dface_vtx_idx,
               dface_vtx,
               dface_cell,
               comm);
  }


  /* Free memory */
  for (int i_zone = 0; i_zone < n_zone; i_zone++) {
    free(dface_bnd_idx [i_zone]);
    free(dface_bnd     [i_zone]);
    free(extract_face_distribution[i_zone]);
    free(extract_vtx_distribution [i_zone]);
    free(dextract_face_vtx_idx    [i_zone]);
    free(dextract_face_vtx        [i_zone]);
    free(dparent_face_g_num       [i_zone]);
    free(dparent_vtx_g_num        [i_zone]);
    free(pextract_old_to_new      [i_zone]);
    free(dextract_vtx_coord       [i_zone]);
    free(dextract_face_group_id   [i_zone]);
    free(dextract_face_group_sens [i_zone]);
    free(dextract_face_join       [i_zone]);
    free(dextract_face_join_opp   [i_zone]);
    PDM_dcube_gen_free(dcube[i_zone]);
  }
  free(dcube);
  free(dn_cell);
  free(dn_face);
  free(dn_vtx);
  free(n_face_group);
  free(dface_group_s);
  free(dface_vtx_s);
  free(dface_cell);
  free(dface_vtx_idx);
  free(dface_vtx);
  free(dvtx_coord);
  free(dface_group_idx);
  free(dface_group);
  free(dface_bnd_idx);
  free(dface_bnd);
  free(dface_join_idx);
  free(dface_join);
  free(group_join_to_zone_cur);
  free(group_join_to_zone_opp);
  free(group_join_to_join_opp);
  free(dextract_face_group_id);
  free(dextract_face_group_sens);
  free(dextract_face_join    );
  free(dextract_face_join_opp);
  free(dface_join_opp);

  free(extract_face_distribution);
  free(extract_vtx_distribution );
  free(dextract_face_vtx_idx    );
  free(dextract_face_vtx        );
  free(dparent_face_g_num       );
  free(dparent_vtx_g_num        );
  free(pextract_old_to_new      );
  free(dextract_vtx_coord       );

  PDM_MPI_Finalize();

  return 0;
}
