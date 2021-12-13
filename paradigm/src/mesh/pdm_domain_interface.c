/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_config.h"
#include "pdm_mpi.h"
#include "pdm_error.h"
#include "pdm_logging.h"
#include "pdm_array.h"
#include "pdm_sort.h"
#include "pdm_unique.h"
#include "pdm_distrib.h"
#include "pdm_part_to_block.h"
#include "pdm_block_to_block.h"
#include "pdm_block_to_part.h"
#include "pdm_multi_block_to_part.h"
#include "pdm_dconnectivity_transform.h"
#include "pdm_dmesh_nodal_to_dmesh.h"
#include "pdm_dmesh_nodal_elements_utils.h"

#include "pdm_domain_interface.h"
#include "pdm_domain_interface_priv.h"

/*=============================================================================
 * Macro definitions
 *============================================================================*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Private function definition
 *============================================================================*/
static int _unique_pairs(int n_pairs, PDM_g_num_t *ids, int *dom_ids) {
  //Rebuild n_occurences
  int n_read = 0;
  PDM_g_num_t last = -1;
  for (int i = 0; i < n_pairs; i++) {
    if (ids[2*i] != last) {
      n_read++;
      last = ids[2*i];
    }
  }
  int *n_occurences = PDM_array_const_int(n_read, 1);
  n_read = -1;
  last   = -1;
  for (int i = 0; i < n_pairs; i++) {
    if (ids[2*i] != last) {
      n_read++;
      last = ids[2*i];
    }
    else {
      n_occurences[n_read]++;
    }
  }
  n_read++; //Add Last
  int max_occur = 0;
  for (int j = 0; j < n_read; j++)
    max_occur = PDM_MAX(max_occur, n_occurences[j]);

  PDM_g_num_t *working_array     = (PDM_g_num_t *) malloc(max_occur*sizeof(PDM_g_num_t));
  PDM_g_num_t *backup_array_gnum = (PDM_g_num_t *) malloc(2*max_occur*sizeof(PDM_g_num_t));
  int         *backup_array_int  = (int *)         malloc(2*max_occur*sizeof(int));
  int         *ordering_array    = (int *)         malloc(max_occur*sizeof(int));

  int start = 0;
  int rewrite_start = 0;
  for (int j = 0; j < n_read; j++) {
    //Fill working array
    for (int k = 0; k < n_occurences[j]; k++)
      working_array[k] = ids[2*(start+k)+1];
    int n_unique = PDM_inplace_unique_long2(working_array, ordering_array, 0, n_occurences[j]-1);
    //Copy into array
    memcpy(backup_array_gnum, &ids[2*start],     2*n_occurences[j]*sizeof(PDM_g_num_t));
    memcpy(backup_array_int,  &dom_ids[2*start], 2*n_occurences[j]*sizeof(int));
    for (int k = 0; k < n_occurences[j]; k++) {
      int pos = ordering_array[k];
      ids[2*(rewrite_start+pos)]   = backup_array_gnum[2*k];
      ids[2*(rewrite_start+pos)+1] = backup_array_gnum[2*k+1];
      dom_ids[2*(rewrite_start+pos)]   = backup_array_int[2*k];
      dom_ids[2*(rewrite_start+pos)+1] = backup_array_int[2*k+1];
    }
    start += n_occurences[j];
    rewrite_start += n_unique;
  }

  free(n_occurences);
  free(working_array);
  free(ordering_array);
  free(backup_array_gnum);
  free(backup_array_int);
  return rewrite_start;
}

static PDM_g_num_t* _per_block_offset(int n_block, int *sizes, PDM_MPI_Comm comm) {
  PDM_g_num_t *sizes_as_gn = (PDM_g_num_t *) malloc(n_block*sizeof(PDM_g_num_t));
  for (int i = 0; i < n_block; i ++)
    sizes_as_gn[i] = sizes[i];
  PDM_g_num_t *per_block_offset = (PDM_g_num_t *) malloc((n_block+1) * sizeof(PDM_g_num_t));
  per_block_offset[0] = 0;
  PDM_MPI_Allreduce(sizes_as_gn, &per_block_offset[1], n_block, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, comm);
  PDM_array_accumulate_gnum(per_block_offset, n_block+1);
  free(sizes_as_gn);
  return per_block_offset;
}

static int _extract_and_shift_jn_faces
(
 int           n_zone,
 int          *dn_face,
 int          *dn_vtx,
 int           n_interface,
 int          *interfaces_size,
 PDM_g_num_t **interface_face_ids,
 int         **interface_domains_ids,
 int         **dface_vtx_idx,
 PDM_g_num_t **dface_vtx,
 int         **face_vtx_both_idx,
 PDM_g_num_t **face_vtx_both,
 int         **dextract_face_group_id,
 int         **dextract_face_dom_id,
 PDM_g_num_t **dextract_face_id,
 PDM_MPI_Comm  comm
)
{
  int n_rank;
  int i_rank;
  PDM_MPI_Comm_size(comm, &n_rank);
  PDM_MPI_Comm_rank(comm, &i_rank);

  PDM_g_num_t *face_per_block_offset = _per_block_offset(n_zone, dn_face, comm);
  PDM_g_num_t *vtx_per_block_offset  = _per_block_offset(n_zone, dn_vtx,  comm);
  
  int n_face_join = 0; // Each interface comes with a pair of faces
  for (int i_interface = 0; i_interface < n_interface; i_interface++) {
    n_face_join += 2*interfaces_size[i_interface];
  }

  PDM_g_num_t *_dextract_face_id_tmp       = malloc(n_face_join * sizeof(PDM_g_num_t));
  // Also transport some data to the extracted faces
  int         *_dextract_face_group_id_tmp = malloc(n_face_join*sizeof(int));
  int         *_dextract_face_dom_id_tmp   = malloc(n_face_join*sizeof(int));

  int idx = 0;
  for (int i_interface = 0; i_interface < n_interface; i_interface++) {

    PDM_g_num_t *_interface_ids = interface_face_ids[i_interface]; //Shortcuts
    int         *_interface_dom = interface_domains_ids[i_interface];
    for (int i_pair = 0; i_pair < interfaces_size[i_interface]; i_pair++) {
      int i_zone_cur = _interface_dom[2*i_pair];
      int i_zone_opp = _interface_dom[2*i_pair+1];

      _dextract_face_id_tmp[idx]         = _interface_ids[2*i_pair] + face_per_block_offset[i_zone_cur];
      _dextract_face_dom_id_tmp  [idx]   = i_zone_cur;
      _dextract_face_group_id_tmp[idx++] = i_interface;

      _dextract_face_id_tmp[idx]         = _interface_ids[2*i_pair+1] + face_per_block_offset[i_zone_opp];
      _dextract_face_dom_id_tmp  [idx]   = i_zone_opp;
      _dextract_face_group_id_tmp[idx++] = i_interface;
    }
  }
  assert (idx == n_face_join);
  
  // Multi gnum is not equilibrated, we have to redistribute it but we want to keep the face/face_opp groups
  PDM_g_num_t *cur_distri   = PDM_compute_entity_distribution(comm, n_face_join/2);
  PDM_g_num_t *ideal_distri = PDM_compute_uniform_entity_distribution(comm, cur_distri[n_rank]);
  PDM_block_to_block_t *btb = PDM_block_to_block_create(cur_distri, ideal_distri, comm);

  //dextract_face_id will be used to get the face->vtx of extracted faces, and also
  //to build the hash table later so return it
  PDM_block_to_block_exch(btb,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_CST,
                          2,
                          NULL,
                          _dextract_face_id_tmp,
                          NULL,
              (void **)   dextract_face_id);
  PDM_block_to_block_exch(btb,
                          sizeof(int),
                          PDM_STRIDE_CST,
                          2,
                          NULL,
                          _dextract_face_group_id_tmp,
                          NULL,
              (void **)  dextract_face_group_id);
  PDM_block_to_block_exch(btb,
                          sizeof(int),
                          PDM_STRIDE_CST,
                          2,
                          NULL,
                          _dextract_face_dom_id_tmp,
                          NULL,
              (void **)  dextract_face_dom_id);

  
  PDM_block_to_block_free(btb);
  // Update n_face_join before freeing distribution
  n_face_join = 2*(ideal_distri[i_rank+1]-ideal_distri[i_rank]);
  free(ideal_distri);
  free(cur_distri);
  free(_dextract_face_id_tmp);
  free(_dextract_face_group_id_tmp);
  free(_dextract_face_dom_id_tmp);

  if (0 == 1) {
    PDM_log_trace_array_long(face_per_block_offset, n_zone+1, "face_per_block_offset :: ");
    PDM_log_trace_array_long(vtx_per_block_offset,  n_zone+1, "vtx_per_block_offset :: ");
    PDM_log_trace_array_long(*dextract_face_id, n_face_join, "dextract_face_id :: ");
  }

  PDM_g_num_t **all_face_distribution = malloc(n_zone * sizeof(PDM_g_num_t*));
  for (int i_zone = 0; i_zone < n_zone; i_zone++) {
    all_face_distribution[i_zone] = PDM_compute_entity_distribution(comm, dn_face[i_zone]);
  }

  PDM_multi_block_to_part_t *mptb = PDM_multi_block_to_part_create(face_per_block_offset,
                                                                   n_zone,
                                            (const PDM_g_num_t **) all_face_distribution,
                                            (const PDM_g_num_t **) dextract_face_id,
                                                                  &n_face_join,
                                                                   1,
                                                                   comm);
  //Prepare data to send : face -> vtx connectivity 
  int         **face_vtx_n       = malloc(n_zone * sizeof(int*));
  PDM_g_num_t **face_vtx_shifted = malloc(n_zone * sizeof(PDM_g_num_t*));
  for (int i_zone = 0; i_zone < n_zone; i_zone++) {
    face_vtx_n[i_zone]       = malloc(dn_face[i_zone] * sizeof(int));
    face_vtx_shifted[i_zone] = malloc(dface_vtx_idx[i_zone][dn_face[i_zone]] * sizeof(PDM_g_num_t));

    for (int i_face = 0; i_face < dn_face[i_zone]; i_face++) {
      face_vtx_n[i_zone][i_face] = dface_vtx_idx[i_zone][i_face+1] - dface_vtx_idx[i_zone][i_face];
      for (int j = dface_vtx_idx[i_zone][i_face]; j < dface_vtx_idx[i_zone][i_face+1]; j++) {
        face_vtx_shifted[i_zone][j] = dface_vtx[i_zone][j] + vtx_per_block_offset[i_zone];
      }
    }
  }

  int         **part_stride = NULL;
  PDM_g_num_t **part_data   = NULL;
  PDM_multi_block_to_part_exch2(mptb,
                                sizeof(PDM_g_num_t),
                                PDM_STRIDE_VAR,
                                face_vtx_n,
                 (void **)      face_vtx_shifted,
                               &part_stride,
                 (void ***)    &part_data);

  *face_vtx_both_idx = PDM_array_new_idx_from_sizes_int(part_stride[0], n_face_join);
  *face_vtx_both     = part_data[0];

  free(part_data);
  free(part_stride[0]);
  free(part_stride);

  PDM_multi_block_to_part_free(mptb);

  for (int i_zone = 0; i_zone < n_zone; i_zone++) {
    free(all_face_distribution[i_zone]);
    free(face_vtx_n[i_zone]);
    free(face_vtx_shifted[i_zone]);
  }
  free(all_face_distribution);
  free(face_vtx_n);
  free(face_vtx_shifted);

  free(face_per_block_offset);
  free(vtx_per_block_offset);

  return n_face_join;
}
static int _generate_edge_face
(
 int            n_face,
 int           *face_vtx_idx,
 PDM_g_num_t   *face_vtx,
 PDM_g_num_t  **dedge_distrib,
 int          **dedge_vtx_idx,
 PDM_g_num_t  **dedge_vtx,
 int          **dedge_face_idx,
 PDM_g_num_t  **dedge_face,
 PDM_MPI_Comm  comm
)
{

  int i_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  
  // 1. Get the number of unique vertex
  PDM_part_to_block_t *ptb = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                      PDM_PART_TO_BLOCK_POST_MERGE,
                                                      1.,
                                                     &face_vtx,
                                                      NULL,
                                                     &face_vtx_idx[n_face],
                                                      1,
                                                      comm);
  int n_vtx;
  int n_vtx_loc = PDM_part_to_block_n_elt_block_get(ptb);
  PDM_MPI_Allreduce(&n_vtx_loc, &n_vtx, 1, PDM_MPI_INT, PDM_MPI_SUM, comm);
  PDM_part_to_block_free(ptb);

  // 2. Get (unmerged) edges
  int n_edge_unmerged = face_vtx_idx[n_face];
  int n_elmt_current  = 0;
  int n_edge_current  = 0;

  PDM_g_num_t *face_distri         = PDM_compute_entity_distribution(comm, n_face);
  PDM_g_num_t *tmp_edge_face       = (PDM_g_num_t *) malloc(n_edge_unmerged       * sizeof(PDM_g_num_t));
  int         *tmp_parent_elmt_pos = (int         *) malloc(n_edge_unmerged       * sizeof(int        ));
  int         *tmp_edge_vtx_idx    = (int         *) malloc((n_edge_unmerged + 1) * sizeof(int        ));
  PDM_g_num_t *tmp_edge_vtx        = (PDM_g_num_t *) malloc(2*n_edge_unmerged     * sizeof(PDM_g_num_t));

  tmp_edge_vtx_idx[0] = 0;
  PDM_poly2d_decomposes_edges(n_face,
                             &n_elmt_current,
                             &n_edge_current,
                              face_distri[i_rank],
                              -1,
                              face_vtx,
                              face_vtx_idx,
                              tmp_edge_vtx_idx, //Numéro de sommet des edges
                              tmp_edge_vtx,
                              tmp_edge_face,         // Numéro des edges pour une face
                              NULL,
                              NULL,
                              tmp_parent_elmt_pos);
  assert(n_edge_current == n_edge_unmerged);

  /*PDM_log_trace_array_long(tmp_dface_edge, n_edge_unmerged, "dface_edge :: ");*/
  /*PDM_log_trace_connectivity_long(tmp_dface_edge_vtx_idx, tmp_dface_edge_vtx, n_edge_unmerged, "dface_edge :: ");*/

  // 3. Merge shared edges
  int dn_edge;
  PDM_generate_entitiy_connectivity_raw(comm,
                                        n_vtx,
                                        n_edge_unmerged,
                                        tmp_edge_face,
                                        tmp_edge_vtx_idx,
                                        tmp_edge_vtx,
                                       &dn_edge,
                                        dedge_distrib,
                                        dedge_vtx_idx,
                                        dedge_vtx,
                                        dedge_face_idx,
                                        dedge_face);
  free(tmp_parent_elmt_pos);
  free(face_distri);
  return dn_edge;
}

static int _match_internal_edges
(
 int            dn_edge,
 PDM_g_num_t   *dedge_distrib,
 int           *dedge_face_idx,
 PDM_g_num_t   *dedge_face,
 PDM_g_num_t   *dedge_face_join,
 PDM_g_num_t   *dedge_face_join_opp,
 PDM_g_num_t  **dedge_gnum,
 PDM_g_num_t  **dedge_gnum_opp,
 PDM_MPI_Comm  comm
)
{
  PDM_UNUSED(dedge_face);
  int i_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  
  // 0. Count the number of internal edges
  int dn_internal_edge = 0;
  for(int i_edge = 0; i_edge < dn_edge; ++i_edge) {
    int n_face_this_edge = dedge_face_idx[i_edge+1] - dedge_face_idx[i_edge];
    dn_internal_edge += (int) (n_face_this_edge > 1);
  }
  if (0 == 1) {
    log_trace("dn internal edges is %i \n", dn_internal_edge);
    log_trace("dn external edges is %i \n", dn_edge - dn_internal_edge);
  }

  // 1. Build hash keys
  PDM_g_num_t *key_ln_to_gn = malloc(dn_internal_edge * sizeof(PDM_g_num_t)); 
  int         *stride_one   = malloc(dn_internal_edge * sizeof(int        ));
  int         *stride_two   = malloc(dn_internal_edge * sizeof(int        ));
  int         *stride_four  = malloc(dn_internal_edge * sizeof(int        ));
  double      *weight       = malloc(dn_internal_edge * sizeof(double     ));

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
    weight[i_int_edge] = 1.;
    //Retrive zone id using group id of any of two faces
    data_send_edge_g_num[i_int_edge] = dedge_distrib[i_rank] + i_edge + 1;

    int key = 0;
    for(int j = dedge_face_idx[i_edge]; j < dedge_face_idx[i_edge+1]; ++j) { //Do it for the two faces data
      key += (dedge_face_join[j] + dedge_face_join_opp[j]);

      //data_send_face_g_num[idx_write2]   = dedge_face[j];
      //data_send_sens      [idx_write2++] = dedge_face_group_sens[j];
      idx_write2++;

      data_send_connect[idx_write4++] = dedge_face_join    [j];

      data_send_connect[idx_write4++] = dedge_face_join_opp[j];
    }
    key_ln_to_gn[i_int_edge] = key;

    i_int_edge++;
  }
  assert(idx_write2 == 2*dn_internal_edge);
  assert(idx_write4 == 4*dn_internal_edge);

  // 2. Exchange data over hash key
  //Attention, pb d'équilibrage car les clés sont réparties vers la fin ... Un proc risque
  // de se retrouver avec tt les clés
  PDM_part_to_block_t *ptb = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                      PDM_PART_TO_BLOCK_POST_MERGE,
                                                      1.,
                                    (PDM_g_num_t **) &key_ln_to_gn,
                                                     &weight,
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
  
  int         *unused_recv_stride = NULL;
  PDM_g_num_t *blk_edge_g_num     = NULL;
  int exch_size = PDM_part_to_block_exch(ptb,
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

  /*
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
  */

  PDM_g_num_t *blk_data_connect = NULL;
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

  if (0 == 1) {
    PDM_log_trace_array_long(key_ln_to_gn, dn_internal_edge, "key_ln_to_gn :: ");
    PDM_log_trace_array_int(gnum_n_occurences   , blk_size               , "gnum_n_occurences   :: ");
    PDM_log_trace_array_long(blk_edge_g_num     , gnum_n_occurences_tot  , "blk_edge_g_num      :: ");
    //PDM_log_trace_array_long(blk_data_face_g_num, 2*gnum_n_occurences_tot, "blk_data_face_g_num :: ");
    //PDM_log_trace_array_long(blk_data_sens      , 2*gnum_n_occurences_tot, "blk_data_sens       :: ");
    PDM_log_trace_array_long(blk_data_connect   , 4*gnum_n_occurences_tot, "blk_data_connect    :: ");
  }


  free(key_ln_to_gn        );
  free(data_send_connect   );
  free(data_send_group     );
  free(data_send_sens      );
  free(data_send_face_g_num);
  free(data_send_edge_g_num);
  free(stride_one          );
  free(stride_two          );
  free(stride_four         );
  free(weight              );



  // 3. Post treatemement : resolve conflicting keys
  PDM_g_num_t *results_edge     = malloc(gnum_n_occurences_tot * sizeof(PDM_g_num_t));
  PDM_g_num_t *results_edge_opp = malloc(gnum_n_occurences_tot * sizeof(PDM_g_num_t));

  int n_max_entity_per_key = 0;
  for(int i = 0; i < blk_size; ++i) {
    n_max_entity_per_key = PDM_MAX(gnum_n_occurences   [i], n_max_entity_per_key);
  }

  int *already_treat   = (int *) malloc(n_max_entity_per_key * sizeof(int));
  int *same_entity_idx = (int *) malloc(n_max_entity_per_key * sizeof(int));
  // int *sens_entity     = (int *) malloc(n_max_entity_per_key * sizeof(int));


  int idx  = 0;
  int idx_w = 0;
  for(int i_key = 0; i_key < blk_size; ++i_key) {

    int n_matching_edge = gnum_n_occurences[i_key];

    /* Reset */
    PDM_array_reset_int(already_treat, n_matching_edge, -1);

    /* Loop over all entitys in conflict and sort all */
    for(int i_entity = 0; i_entity < n_matching_edge; ++i_entity) {

      //Each internal edge comes with 2 faces -> 4 data
      int beg    = 4*(idx+i_entity);

      // Caution inplace sort !!!!
      PDM_sort_long(&blk_data_connect[beg], NULL, 4);

      if(0 == 1) {
        PDM_log_trace_array_long(&blk_data_connect[beg], 4, "blk_data_connect (sort) :: ");
      }
    }

    /*
     *  Identify pair or invalid other //todo : shortcut only 2 occurences or let for debug
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

          if(!PDM_array_are_equal_gnum(&blk_data_connect[beg1], &blk_data_connect[beg2], 4)) {
            continue;
          }

          already_treat[i_entity2] = 1;
          same_entity_idx[n_same++] = i_entity2;
          i_entity2_same = i_entity2;

        }
      } /* End for i_entity2 */
      assert(n_same <= 1);


      if (n_same == 1) {
        int edge_idx     = (idx+i_entity);
        int edge_idx_opp = (idx+i_entity2_same);

        // Set data for edge
        results_edge    [idx_w] = blk_edge_g_num[edge_idx];
        results_edge_opp[idx_w++] = blk_edge_g_num[edge_idx_opp];

        // Set data for opposite edge
        results_edge    [idx_w] = blk_edge_g_num[edge_idx_opp];
        results_edge_opp[idx_w++] = blk_edge_g_num[edge_idx];
      }

      already_treat[i_entity] = 1;
    }
    idx  += n_matching_edge;
  }
  free(already_treat);
  free(same_entity_idx);
  // Some pairs can be still unresolved, eg if a edge is internal from one interface point of view but
  // external for the other
  int rsvd_gnum_n_occurences_tot = idx_w;
  results_edge     = realloc(results_edge,     rsvd_gnum_n_occurences_tot*sizeof(PDM_g_num_t));
  results_edge_opp = realloc(results_edge_opp, rsvd_gnum_n_occurences_tot*sizeof(PDM_g_num_t));


  free(blk_edge_g_num);
  free(blk_data_connect);
  PDM_part_to_block_free(ptb); // Needed until here for gnum_n_occurences

  if (0 == 1) {
    log_trace("Conflict resolved, gnum are\n");
    PDM_log_trace_array_long(results_edge, rsvd_gnum_n_occurences_tot, "edge gnum ::");
    PDM_log_trace_array_long(results_edge_opp, rsvd_gnum_n_occurences_tot, "edge gnum opp::");
  }



  // 4. Send back result on edge distribution
                       ptb = PDM_part_to_block_create2(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                       PDM_PART_TO_BLOCK_POST_NOTHING,
                                                       1.,
                                                      &results_edge,
                                                       dedge_distrib,
                                                      &rsvd_gnum_n_occurences_tot,
                                                       1,
                                                       comm);

  int resolved_dn_internal_edge = PDM_part_to_block_n_elt_block_get(ptb);
  *dedge_gnum = malloc(resolved_dn_internal_edge * sizeof(PDM_g_num_t));

  PDM_g_num_t *dedge_gnum_tmp = PDM_part_to_block_block_gnum_get(ptb);
  memcpy(*dedge_gnum, dedge_gnum_tmp, resolved_dn_internal_edge*sizeof(PDM_g_num_t));

  PDM_part_to_block_exch(ptb,
                        sizeof(PDM_g_num_t),
                        PDM_STRIDE_CST,
                        1,
                        NULL,
             (void **) &(results_edge_opp),
                        NULL,
              (void **) dedge_gnum_opp);
  PDM_part_to_block_free(ptb);

  free(results_edge);
  free(results_edge_opp);

  return resolved_dn_internal_edge;
}

static void _match_all_edges_from_faces
(
  int          dn_face,
  int         *face_edge_idx,
  PDM_g_num_t *face_edge,
  PDM_g_num_t *face_edge_wopp,
  PDM_g_num_t *pedge_vtx,
  PDM_g_num_t *p_all_vtx,
  PDM_g_num_t *p_all_vtx_opp
)
{
  // Avoid multiple malloc using max size
  int max_face_len = 0;
  for (int i_face = 0; i_face < dn_face/2; i_face++) {
    max_face_len = PDM_MAX(max_face_len, face_edge_idx[2*i_face+1] - face_edge_idx[2*i_face]);
  }
  PDM_g_num_t *ordered_edge     = malloc(max_face_len * sizeof(PDM_g_num_t));
  PDM_g_num_t *ordered_edge_opp = malloc(max_face_len * sizeof(PDM_g_num_t));
  PDM_g_num_t *ordered_vtx      = malloc(max_face_len * sizeof(PDM_g_num_t));
  PDM_g_num_t *ordered_vtx_opp  = malloc(max_face_len * sizeof(PDM_g_num_t));

  // Avec la construction des faces de bord, on a des paires faces / face opp
  assert (dn_face%2 == 0);
  assert (face_edge_idx[dn_face]%2 == 0);
  int glob_idx     = 0;
  for (int i_face = 0; i_face < dn_face/2; i_face++) {
    int start_idx     = face_edge_idx[2*i_face];
    int start_idx_opp = face_edge_idx[2*i_face+1];
    int face_len = face_edge_idx[2*i_face+1] - face_edge_idx[2*i_face];
    if (0 == 1) {
      log_trace("iface %d :\n", i_face);
      log_trace("face data\n");
      // PDM_log_trace_array_long(&face_vtx[face_edge_idx[i_face]], face_len, "face_vtx");
      PDM_log_trace_array_long(&face_edge[start_idx], face_len, "face_edge");
      PDM_log_trace_array_long(&pedge_vtx[2*start_idx], 2*face_len, "edge vertices");
      PDM_log_trace_array_long(&face_edge_wopp[start_idx], face_len, "face_edge_wopp");
      log_trace("face opp data\n");
      // PDM_log_trace_array_long(&face_vtx_opp[face_edge_idx[i_face]], face_len, "face_vtx");
      PDM_log_trace_array_long(&face_edge[start_idx_opp], face_len, "face_opp_edge_opp");
      PDM_log_trace_array_long(&pedge_vtx[2*start_idx_opp], 2*face_len, "edge vertices opp");
    }

    // Search any received edge (we should have at least one)
    PDM_g_num_t opp_edge_key = 0;
    int idx = -1;
    while (opp_edge_key == 0)
      opp_edge_key = face_edge_wopp[start_idx + 1 + idx++];
    assert (idx < face_len);

    //Search idx of opposite edge in opposite 
    PDM_g_num_t candidate = 0;
    int opp_idx = -1;
    while (candidate != opp_edge_key) {
      candidate = face_edge[start_idx_opp + 1 + opp_idx++];
      candidate = PDM_ABS(candidate);
    }
    assert (opp_idx < face_len);

    // Find starting point using the two edge indices and signs
    int sign     = PDM_SIGN(face_edge[start_idx + idx]);
    int opp_sign = PDM_SIGN(face_edge[start_idx_opp + opp_idx]);

    int next_vtx, next_vtx_opp, cur_vtx, cur_vtx_opp;
    if (sign == 1) {
      cur_vtx  = pedge_vtx[2*(start_idx + idx)];
      next_vtx = pedge_vtx[2*(start_idx + idx)+1];
    }
    else {
      cur_vtx  = pedge_vtx[2*(start_idx + idx)+1];
      next_vtx = pedge_vtx[2*(start_idx + idx)];
    }
    if (opp_sign == 1) {//Invert looping order for opposite face
      cur_vtx_opp  = pedge_vtx[2*(start_idx_opp + opp_idx)+1];
      next_vtx_opp = pedge_vtx[2*(start_idx_opp + opp_idx)];
    }
    else {
      cur_vtx_opp  = pedge_vtx[2*(start_idx_opp + opp_idx)];
      next_vtx_opp = pedge_vtx[2*(start_idx_opp + opp_idx)+1];
    }

    for (int i = 0; i < face_len; i++) {
      //Fill
      // log_trace("Cur vtx %d and opp %d \n", cur_vtx, cur_vtx_opp);
      ordered_edge[i]     = PDM_ABS(face_edge[start_idx+idx]);
      ordered_edge_opp[i] = PDM_ABS(face_edge[start_idx_opp+opp_idx]);
      ordered_vtx[i]      = cur_vtx;
      ordered_vtx_opp[i]  = cur_vtx_opp;
      //This is for face
      for (int j = 0; j < face_len; j++) {
        if (j != idx) {
          int vtx1 = pedge_vtx[2*(start_idx + j)];
          int vtx2 = pedge_vtx[2*(start_idx + j)+1];
          if (vtx1 == next_vtx) {
            idx = j;
            cur_vtx  = vtx1;
            next_vtx = vtx2;
            break;
          }
          else if (vtx2 == next_vtx) {
            idx = j;
            cur_vtx  = vtx2;
            next_vtx = vtx1;
            break;
          }
        }
      }
      //This is for opposite face
      for (int j = 0; j < face_len; j++) {
        if (j != opp_idx) {
          int vtx1 = pedge_vtx[2*(start_idx_opp+j)];
          int vtx2 = pedge_vtx[2*(start_idx_opp+j)+1];
          if (vtx1 == next_vtx_opp) {
            opp_idx = j;
            cur_vtx_opp  = vtx1;
            next_vtx_opp = vtx2;
            break;
          }
          else if (vtx2 == next_vtx_opp) {
            opp_idx = j;
            cur_vtx_opp  = vtx2;
            next_vtx_opp = vtx1;
            break;
          }
        }
      }
    }
    if (0 == 1) {
      PDM_log_trace_array_long(ordered_edge,     face_len, "ordered edges");
      PDM_log_trace_array_long(ordered_edge_opp, face_len, "ordered edges_opp");
      PDM_log_trace_array_long(ordered_vtx,     face_len, "ordered vtx");
      PDM_log_trace_array_long(ordered_vtx_opp, face_len, "ordered vtx_opp");
    }
    //Copy results for this face
    memcpy(&p_all_vtx[glob_idx],      ordered_vtx,      face_len*sizeof(PDM_g_num_t));
    memcpy(&p_all_vtx_opp[glob_idx],  ordered_vtx_opp,  face_len*sizeof(PDM_g_num_t));
    //memcpy(&p_all_edge[glob_idx],     ordered_edge      face_len*sizeof(PDM_g_num_t));
    //memcpy(&p_all_edge_opp[glob_idx], ordered_edge_opp, face_len*sizeof(PDM_g_num_t));
    glob_idx     += face_len;
  }
  free(ordered_edge);
  free(ordered_edge_opp);
  free(ordered_vtx);
  free(ordered_vtx_opp);
}

static void
_create_vtx_join
(
int            n_interface,
int            p_all_vtx_n,
PDM_g_num_t   *p_all_vtx,
PDM_g_num_t   *p_all_vtx_opp,
int           *p_all_vtx_group,
int           *p_all_vtx_dom_id,
int           *p_all_vtx_domopp_id,
int           *vtx_interface_size,
PDM_g_num_t  **interface_vtx_ids,
int          **interface_vtx_dom_ids,
PDM_MPI_Comm   comm
)
{
  //Todo : we could shift back to position of vtx in extraction to have a better
  //balance of edge distribution
  int *stride_one = PDM_array_const_int(p_all_vtx_n, 1);
  //Merge vertices using gnum
  PDM_part_to_block_t *ptb = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                      PDM_PART_TO_BLOCK_POST_MERGE,
                                                      1.,
                                                     &p_all_vtx,
                                                      NULL,
                                                     &p_all_vtx_n,
                                                      1,
                                                      comm);

  int blk_size = PDM_part_to_block_n_elt_block_get(ptb);
  PDM_g_num_t *dall_vtx   = PDM_part_to_block_block_gnum_get(ptb);

  int         *recv_stride      = NULL;
  int         *dall_vtx_dom     = NULL;
  int         *dall_vtx_dom_opp = NULL;
  int         *dall_vtx_group   = NULL;
  PDM_g_num_t *dall_vtx_opp     = NULL;
  int exch_size = PDM_part_to_block_exch(ptb,
                                     sizeof(PDM_g_num_t),
                                     PDM_STRIDE_VAR,
                                     -1,
                           (int  **) &stride_one,
                           (void **) &p_all_vtx_opp,
                                     &recv_stride,
                           (void **) &dall_vtx_opp);

  int *unused_recv_stride = NULL;
  exch_size = PDM_part_to_block_exch(ptb,
                                     sizeof(int),
                                     PDM_STRIDE_VAR,
                                     -1,
                           (int  **) &stride_one,
                           (void **) &p_all_vtx_group,
                                     &unused_recv_stride, //Same  than recv stride
                           (void **) &dall_vtx_group);
  free(unused_recv_stride);

  exch_size = PDM_part_to_block_exch(ptb,
                                     sizeof(int),
                                     PDM_STRIDE_VAR,
                                     -1,
                           (int  **) &stride_one,
                           (void **) &p_all_vtx_dom_id,
                                     &unused_recv_stride, //Same  than recv stride
                           (void **) &dall_vtx_dom);
  free(unused_recv_stride);

  exch_size = PDM_part_to_block_exch(ptb,
                                     sizeof(int),
                                     PDM_STRIDE_VAR,
                                     -1,
                           (int  **) &stride_one,
                           (void **) &p_all_vtx_domopp_id,
                                     &unused_recv_stride, //Same  than recv stride
                           (void **) &dall_vtx_dom_opp);
  free(unused_recv_stride);
  free(stride_one);

  if (0 == 1) {
    PDM_log_trace_array_long(dall_vtx,       blk_size,  "dall_vtx            :");
    PDM_log_trace_array_int (recv_stride,    blk_size,  "recv stride         :");
    PDM_log_trace_array_long(dall_vtx_opp,   exch_size, "recv dall_vtx_opp   :");
    PDM_log_trace_array_int (dall_vtx_group, exch_size, "recv dall_vtx_group :");
  }

  //First, dispatch vertices depending of the original interface
  PDM_array_reset_int(vtx_interface_size, n_interface, 0);
  int start_vtx = 0;
  for (int i_vtx =  0; i_vtx < blk_size; i_vtx++) {
    int n_recv = recv_stride[i_vtx];
    for (int i = 0; i < n_recv; i++) {
      vtx_interface_size[dall_vtx_group[start_vtx + i]]++;
    }
    start_vtx += n_recv;
  }
  for (int i_interface = 0; i_interface < n_interface; i_interface++) {
    interface_vtx_ids[i_interface]     = malloc(2*vtx_interface_size[i_interface]*sizeof(PDM_g_num_t));
    interface_vtx_dom_ids[i_interface] = malloc(2*vtx_interface_size[i_interface]*sizeof(int));
  }
  PDM_array_reset_int(vtx_interface_size, n_interface, 0);
  start_vtx = 0;
  for (int i_vtx =  0; i_vtx < blk_size; i_vtx++) {
    int n_recv = recv_stride[i_vtx];
    for (int i = 0; i < n_recv; i++) {
      int i_interface = dall_vtx_group[start_vtx + i];
      interface_vtx_ids[i_interface][2*vtx_interface_size[i_interface]]   = dall_vtx[i_vtx];
      interface_vtx_ids[i_interface][2*vtx_interface_size[i_interface]+1] = dall_vtx_opp[start_vtx+i];
      interface_vtx_dom_ids[i_interface][2*vtx_interface_size[i_interface]]   = dall_vtx_dom    [start_vtx+i];
      interface_vtx_dom_ids[i_interface][2*vtx_interface_size[i_interface]+1] = dall_vtx_dom_opp[start_vtx+i];
      vtx_interface_size[i_interface]++;
    }
    start_vtx += n_recv;
  }

  //Then, for each interface, eliminate pairs of vertices occuring more than once
  //If dom_id & dom_opp_id differs, vtx_id & opp should also differ because of the shift
  for (int i_interface = 0; i_interface < n_interface; i_interface++) {
    
    int n_pairs_u = _unique_pairs(vtx_interface_size[i_interface],
                                  interface_vtx_ids[i_interface],
                                  interface_vtx_dom_ids[i_interface]);

    //Update
    vtx_interface_size[i_interface] = n_pairs_u;
    interface_vtx_ids[i_interface] = realloc(interface_vtx_ids[i_interface], 2*n_pairs_u*sizeof(PDM_g_num_t));
    interface_vtx_dom_ids[i_interface] = realloc(interface_vtx_dom_ids[i_interface], 2*n_pairs_u*sizeof(int));
  }

  PDM_part_to_block_free(ptb);
  free(recv_stride);
  free(dall_vtx_dom);
  free(dall_vtx_dom_opp);
  free(dall_vtx_opp);
  free(dall_vtx_group);
}

static void _domain_interface_face_to_vertex
(
 int            n_interface,             /* Total number of interfaces */
 int           *interfaces_size,         /* Number of face pairs in each interface */
 PDM_g_num_t  **interface_face_ids,      /* For each interface, list of pairs face,face_opp */
 int          **interface_domains_ids,   /* For each interface, list of domains dom,dom_opp */
 int            n_zone,                  /* Number of zones */
 int           *dn_vtx,                  /* Number of vertex in each zone (distributed) */
 int           *dn_face,                 /* Number of face in each zone (distributed) */
 int          **dface_vtx_idx,           /* Face->vertex connectivity for each domain */
 PDM_g_num_t  **dface_vtx,
 int           *vtx_interface_size,      /* [OUT] Number of vtx pairs in each interface */
 PDM_g_num_t  **interface_vtx_ids,       /* [OUT] For each interface, list of pairs vtx,vtx_opp */
 int          **interface_vtx_dom_ids,   /* [OUT] For each interface, list of domains dom,dom_opp */
 PDM_MPI_Comm   comm                     /* Mpi comunicator */
)
{
  int n_rank;
  int i_rank;
  PDM_MPI_Comm_size(comm, &n_rank);
  PDM_MPI_Comm_rank(comm, &i_rank);

  PDM_g_num_t *face_per_block_offset = _per_block_offset(n_zone, dn_face, comm);
  PDM_g_num_t *vtx_per_block_offset  = _per_block_offset(n_zone, dn_vtx,  comm);

  int         *face_vtx_both_idx = NULL;
  PDM_g_num_t *face_vtx_both     = NULL;
  int         *dextract_face_dom_id   = NULL;
  int         *dextract_face_group_id = NULL;
  PDM_g_num_t *dextract_face_join     = NULL;
  int n_extr_face = _extract_and_shift_jn_faces(n_zone,
                                                dn_face,
                                                dn_vtx,
                                                n_interface,
                                                interfaces_size,
                                                interface_face_ids,
                                                interface_domains_ids,
                                                dface_vtx_idx,
                                                dface_vtx,
                                               &face_vtx_both_idx,
                                               &face_vtx_both,
                                               &dextract_face_group_id,
                                               &dextract_face_dom_id,
                                               &dextract_face_join,
                                                comm);

  PDM_g_num_t *extracted_face_distri = PDM_compute_entity_distribution(comm, n_extr_face);

  //Duplicate this data for easier send to edges
  PDM_g_num_t *dextract_face_join_opp   = malloc(n_extr_face*sizeof(PDM_g_num_t));
  for (int i = 0; i < n_extr_face/2; i++) {
    dextract_face_join_opp[2*i]     = dextract_face_join[2*i+1];
    dextract_face_join_opp[2*i+1]   = dextract_face_join[2*i];
  }

  if (0 == 1) {
    log_trace("Face vtx received after MBTP\n");
    PDM_log_trace_array_int (face_vtx_both_idx, n_extr_face+1, "face_vtx_idx :: ");
    PDM_log_trace_array_long(face_vtx_both, face_vtx_both_idx[n_extr_face], "face_vtx :: ");
  }


  //Now we have some pairs of faces (each pair appears only one) + face_vtx for this pairs

  //Generate edge numbering

  PDM_g_num_t *dedge_distrib  = NULL;
  int         *dedge_vtx_idx  = NULL;
  PDM_g_num_t *dedge_vtx      = NULL;
  int         *dedge_face_idx = NULL;
  PDM_g_num_t *dedge_face     = NULL;
  int dn_edge = _generate_edge_face(n_extr_face,
                                    face_vtx_both_idx,
                                    face_vtx_both,
                                   &dedge_distrib,
                                   &dedge_vtx_idx,
                                   &dedge_vtx,
                                   &dedge_face_idx,
                                   &dedge_face,
                                    comm);
  
  if (0 == 1) {
    log_trace("Edges rebuild\n");
    PDM_log_trace_array_long(dedge_distrib, n_rank+1, "dedge_distri ::");
    //PDM_log_trace_array_long(dedge_vtx, 2*dn_edge, "dedge_vtx ::");
    //PDM_log_trace_connectivity_long(dedge_face_idx, dedge_face, dn_edge, "dedge_face ::");
  }

  // Transport face data to edges

  //Prepare numbering
  PDM_g_num_t *dedge_face_abs = (PDM_g_num_t *) malloc(dedge_face_idx[dn_edge] * sizeof(PDM_g_num_t));
  int         *dedge_face_sgn = (int         *) malloc(dedge_face_idx[dn_edge] * sizeof(int        ));
  for(int i = 0; i < dedge_face_idx[dn_edge]; ++i) {
    dedge_face_abs[i] = PDM_ABS (dedge_face[i]);
    dedge_face_sgn[i] = PDM_SIGN(dedge_face[i]);
  }


  PDM_g_num_t *dedge_face_join         = malloc(dedge_face_idx[dn_edge] * sizeof(PDM_g_num_t));
  PDM_g_num_t *dedge_face_join_opp     = malloc(dedge_face_idx[dn_edge] * sizeof(PDM_g_num_t));
  PDM_block_to_part_t *btp = PDM_block_to_part_create(extracted_face_distri,
                               (const PDM_g_num_t **) &dedge_face_abs,
                                                      &dedge_face_idx[dn_edge],
                                                      1,
                                                      comm);
  int cst_stride = 1;
  PDM_block_to_part_exch(btp,
                         sizeof(PDM_g_num_t),
                         PDM_STRIDE_CST,
                         &cst_stride,
                (void *) dextract_face_join,
                         NULL,
             (void ** ) &dedge_face_join);

  PDM_block_to_part_exch(btp,
                         sizeof(PDM_g_num_t),
                         PDM_STRIDE_CST,
                         &cst_stride,
                (void *) dextract_face_join_opp,
                         NULL,
             (void ** ) &dedge_face_join_opp);

  PDM_block_to_part_free(btp);
  free(dedge_face_abs);
  free(dedge_face_sgn);
  free(dextract_face_join_opp);


  if (0 == 1) {
    log_trace("Transport data on edges\n");
    PDM_log_trace_array_int (dedge_face_idx,        dn_edge+1,               "dedge_face_idx        ::");
    PDM_log_trace_array_long(dedge_face_join    ,   dedge_face_idx[dn_edge], "dedge_face_join       ::");
    PDM_log_trace_array_long(dedge_face_join_opp,   dedge_face_idx[dn_edge], "dedge_face_join_opp   ::");
  }

  // Match internal edges
  PDM_g_num_t *dedge_gnum     = NULL;
  PDM_g_num_t *dedge_gnum_opp = NULL;
  int dn_internal_edge = _match_internal_edges(dn_edge,
                                               dedge_distrib,
                                               dedge_face_idx,
                                               dedge_face,
                                               dedge_face_join,
                                               dedge_face_join_opp,
                                              &dedge_gnum,
                                              &dedge_gnum_opp,
                                               comm);

  if (0 == 1) {
    log_trace("Internal edge matches after conflict resolution \n");
    PDM_log_trace_array_long(dedge_gnum, dn_internal_edge, "dedge gnum :: ");
    PDM_log_trace_array_long(dedge_gnum_opp, dn_internal_edge, "dedge gnum_opp :: ");
  }

  // To match external edges, we will work on face distribution. We need face->edge connectivity

  int         *dface_edge_idx = NULL;
  PDM_g_num_t *dface_edge     = NULL;
  PDM_dconnectivity_transpose(comm,
                              dedge_distrib,
                              extracted_face_distri,
                              dedge_face_idx,
                              dedge_face,
                              1,
                             &dface_edge_idx,
                             &dface_edge);

  if (0 == 1) {
    log_trace("Generate dface->edge from dedge->face\n");
    PDM_log_trace_connectivity_long(dface_edge_idx, dface_edge, n_extr_face, "dface_edge :: ");
  }

  //Transfert some data from the edges to the edges know by the faces

  //Prepare gnum and stride
  PDM_g_num_t *dface_edge_abs = (PDM_g_num_t *) malloc(dface_edge_idx[n_extr_face] * sizeof(PDM_g_num_t));
  for (int i = 0; i < dface_edge_idx[n_extr_face]; i++) {
    dface_edge_abs[i] = PDM_ABS(dface_edge[i]);
  }
  int idx = 0;
  int *dedge_gnum_n = PDM_array_zeros_int(dn_edge);
  for (int i = 0; i < dn_edge; i++) {
    if (i + dedge_distrib[i_rank] + 1 == dedge_gnum[idx]) {
      dedge_gnum_n[i] = 1;
      idx++;
    }
    if (idx == dn_internal_edge) { //End of internal edge reached, no more comparaison is needed
      break;
    }
  }

                       btp = PDM_block_to_part_create(dedge_distrib,
                               (const PDM_g_num_t **) &dface_edge_abs,
                                                      &dface_edge_idx[n_extr_face],
                                                      1,
                                                      comm);


  int         **recv_stride_tmp = NULL;
  PDM_g_num_t **recv_data_tmp   = NULL;
  PDM_block_to_part_exch2(btp,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_VAR,
                          dedge_gnum_n,
                          dedge_gnum_opp,
                         &recv_stride_tmp,
              (void ***) &recv_data_tmp);
  int         *pedge_gnum_n   = recv_stride_tmp[0];
  PDM_g_num_t *pedge_gnum_opp = recv_data_tmp[0];
  free(recv_stride_tmp);
  free(recv_data_tmp);

  int stride2 = 2;
  PDM_block_to_part_exch2(btp,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_CST,
                         &stride2,
                          dedge_vtx,
                          NULL,
              (void ***) &recv_data_tmp);
  PDM_g_num_t *pedge_vtx = recv_data_tmp[0];
  free(recv_data_tmp);
  //PDM_log_trace_array_long(pedge_vtx, 2*dface_edge_idx[n_extr_face], "pedge_vtx");



  //Attention, on devrait pouvoir travailler sur face externes uniquement (filtre le dface_edge_abs)
  PDM_g_num_t *face_edge_wopp = malloc(dface_edge_idx[n_extr_face]*sizeof(PDM_g_num_t));
  idx = 0;
  for (int i = 0; i < dface_edge_idx[n_extr_face]; i++) {
    if (pedge_gnum_n[i] == 1)
      face_edge_wopp[i] = pedge_gnum_opp[idx++];
    else
      face_edge_wopp[i] = 0;
  }

  //PDM_log_trace_connectivity_long(dface_edge_idx, face_edge_wopp, n_extr_face, "dface_edge :: ");

  free(dface_edge_abs);
  free(dedge_gnum_n);
  free(pedge_gnum_n);
  free(pedge_gnum_opp);
  PDM_block_to_part_free(btp);

  //Match external edges
  assert (dface_edge_idx[n_extr_face] % 2 == 0);
  int n_vtx_interface_tot = dface_edge_idx[n_extr_face] / 2;
  PDM_g_num_t *p_all_vtx      = malloc(n_vtx_interface_tot * sizeof(PDM_g_num_t));
  PDM_g_num_t *p_all_vtx_opp  = malloc(n_vtx_interface_tot * sizeof(PDM_g_num_t));
  /*PDM_g_num_t *p_all_edge_gnum     = malloc(dface_edge_idx[n_extr_face] * sizeof(PDM_g_num_t));*/
  /*PDM_g_num_t *p_all_edge_gnum_opp = malloc(dface_edge_idx[n_extr_face] * sizeof(PDM_g_num_t));*/
  _match_all_edges_from_faces(n_extr_face,
                              dface_edge_idx,
                              dface_edge,
                              face_edge_wopp,
                              pedge_vtx,
                              p_all_vtx,
                              p_all_vtx_opp);

  //Copy group from face to vertices
  int *p_all_vtx_group     = malloc(n_vtx_interface_tot * sizeof(int));
  int *p_all_vtx_dom_id    = malloc(n_vtx_interface_tot * sizeof(int));
  int *p_all_vtx_domopp_id = malloc(n_vtx_interface_tot * sizeof(int));
  int glob_idx = 0;
  for (int i_face = 0; i_face < n_extr_face/2; i_face++) {
    int face_len = dface_edge_idx[2*i_face+1] - dface_edge_idx[2*i_face];
    for (int i = 0; i < face_len; i++) {
      p_all_vtx_group[glob_idx]     = dextract_face_group_id[2*i_face];
      p_all_vtx_dom_id[glob_idx]    = dextract_face_dom_id  [2*i_face];
      p_all_vtx_domopp_id[glob_idx] = dextract_face_dom_id  [2*i_face+1];
      glob_idx ++;
    }
  }


  free(face_edge_wopp);
  if (0 == 1) {
    log_trace("Vtx matching on face distribution\n");
    /*PDM_log_trace_array_long(p_all_edge_gnum,     dface_edge_idx[n_extr_face], "p_all_edge_gnum     ::");*/
    /*PDM_log_trace_array_long(p_all_edge_gnum_opp, dface_edge_idx[n_extr_face], "p_all_edge_gnum_opp ::");*/
    PDM_log_trace_array_long(p_all_vtx,     n_vtx_interface_tot, "p_all_vtx     ::");
    PDM_log_trace_array_long(p_all_vtx_opp, n_vtx_interface_tot, "p_all_vtx_opp ::");
    PDM_log_trace_array_int (p_all_vtx_group, n_vtx_interface_tot, "p_all_vtx_group ::");
    PDM_log_trace_array_int (p_all_vtx_dom_id, n_vtx_interface_tot, "p_all_vtx_dom_id ::");
    PDM_log_trace_array_int (p_all_vtx_domopp_id, n_vtx_interface_tot, "p_all_vtx_domopp_id ::");
  }

  _create_vtx_join(n_interface,
                   n_vtx_interface_tot,
                   p_all_vtx,
                   p_all_vtx_opp,
                   p_all_vtx_group,
                   p_all_vtx_dom_id,
                   p_all_vtx_domopp_id,
                   vtx_interface_size,
                   interface_vtx_ids,
                   interface_vtx_dom_ids,
                   comm);


  //Ultimate step : go back to original vtx numbering. All we have to do is retrieve zone
  // and substract zone offset
  for (int i_interface = 0; i_interface < n_interface; i_interface++) {
    for (int i_vtx = 0; i_vtx < 2*vtx_interface_size[i_interface]; i_vtx++) {
      interface_vtx_ids[i_interface][i_vtx] -= vtx_per_block_offset[interface_vtx_dom_ids[i_interface][i_vtx]];
    }
  }
  if (0 == 1) {
    PDM_log_trace_array_int(vtx_interface_size, n_interface, "Vtx interfaces sizes");
    for (int i = 0; i < n_interface; i++) {
      log_trace("Vtx & vtx opp in vtx global numbering for interface %d \n", i);
      PDM_log_trace_array_long(interface_vtx_ids[i], 2*vtx_interface_size[i], "vertex ids     ::");
      PDM_log_trace_array_int(interface_vtx_dom_ids[i], 2*vtx_interface_size[i], "vertex doms     ::");
    }
  }



  free(dface_edge_idx);
  free(dface_edge);

  free(dextract_face_join);
  free(dextract_face_group_id);
  free(dextract_face_dom_id);

  free(dedge_face_join);
  free(dedge_face_join_opp);


  free(extracted_face_distri);

  free(face_vtx_both_idx);
  free(face_vtx_both);

  free(dedge_distrib);
  free(dedge_vtx_idx);
  free(dedge_vtx);
  free(dedge_face_idx);
  free(dedge_face);

  free(dedge_gnum);
  free(dedge_gnum_opp);

  free(p_all_vtx);
  free(p_all_vtx_opp);
  free(p_all_vtx_group);
  free(p_all_vtx_dom_id);
  free(p_all_vtx_domopp_id);

  free(pedge_vtx);
  free(face_per_block_offset);
  free(vtx_per_block_offset);
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

PDM_domain_interface_t* PDM_domain_interface_create
(
const int                   n_interface,
const int                   n_zone,
PDM_domain_interface_mult_t multizone_interface,
PDM_ownership_t             ownership,
PDM_MPI_Comm                comm
)
{
  PDM_domain_interface_t *dom_intrf = malloc (sizeof(PDM_domain_interface_t));
  dom_intrf->n_interface       = n_interface;
  dom_intrf->n_zone            = n_zone;
  dom_intrf->multizone_intrf   = multizone_interface;
  dom_intrf->ownership         = ownership;
  dom_intrf->comm              = comm;

  dom_intrf->interface_dn_face   = NULL;
  dom_intrf->interface_ids_face  = NULL;
  dom_intrf->interface_dom_face  = NULL;
  dom_intrf->interface_dn_vtx    = NULL;
  dom_intrf->interface_ids_vtx   = NULL;
  dom_intrf->interface_dom_vtx   = NULL;

  for (int i = 0; i < PDM_BOUND_TYPE_MAX; i++)
    dom_intrf->is_result[i] = 0;

  return dom_intrf;
}

void PDM_domain_interface_set
(
 PDM_domain_interface_t *dom_intrf,
 PDM_bound_type_t        interface_kind,
 int                    *interface_dn,
 PDM_g_num_t           **interface_ids,
 int                   **interface_dom
)
{
  assert (dom_intrf != NULL);
  if (interface_kind != PDM_BOUND_TYPE_FACE) {
    PDM_error(__FILE__, __LINE__, 0, "Only face domain connectivity is currently supported\n");
  }
  else {
    dom_intrf->interface_dn_face  = interface_dn;
    dom_intrf->interface_ids_face = interface_ids;
    dom_intrf->interface_dom_face = interface_dom;
  }

}

void PDM_domain_interface_translate_face2vtx
(
 PDM_domain_interface_t  *dom_intrf,
 int                     *dn_vtx,
 int                     *dn_face,
 int                    **dface_vtx_idx,
 PDM_g_num_t            **dface_vtx
)
{
  assert (dom_intrf != NULL);
  assert (dom_intrf->interface_dn_face != NULL);
  assert (dom_intrf->interface_dn_vtx  == NULL);
  dom_intrf->interface_dn_vtx  = (int *)          malloc(dom_intrf->n_interface * sizeof(int));
  dom_intrf->interface_ids_vtx = (PDM_g_num_t **) malloc(dom_intrf->n_interface * sizeof(PDM_g_num_t*));
  dom_intrf->interface_dom_vtx = (int         **) malloc(dom_intrf->n_interface * sizeof(int*));

  // Simple case is not yet managed, copy to go back to full case
  int **_interface_dom_face = NULL;
  if (dom_intrf->multizone_intrf == PDM_DOMAIN_INTERFACE_MULT_NO) {
    _interface_dom_face = (int **) malloc(dom_intrf->n_interface*sizeof(int*));
    for (int i_intrf = 0; i_intrf < dom_intrf->n_interface; i_intrf++) {
      _interface_dom_face[i_intrf] = (int *) malloc(2*dom_intrf->interface_dn_face[i_intrf]*sizeof(int));
      for (int j = 0; j < dom_intrf->interface_dn_face[i_intrf]; j++) {
        _interface_dom_face[i_intrf][2*j]   = dom_intrf->interface_dom_face[i_intrf][0];
        _interface_dom_face[i_intrf][2*j+1] = dom_intrf->interface_dom_face[i_intrf][1];
      }
    }
  }
  else {
    _interface_dom_face = dom_intrf->interface_dom_face;
  }

  _domain_interface_face_to_vertex(dom_intrf->n_interface,
                                   dom_intrf->interface_dn_face,
                                   dom_intrf->interface_ids_face,
                                   _interface_dom_face,
                                   dom_intrf->n_zone,
                                   dn_vtx,
                                   dn_face,
                                   dface_vtx_idx,
                                   dface_vtx,
                                   dom_intrf->interface_dn_vtx,
                                   dom_intrf->interface_ids_vtx,
                                   dom_intrf->interface_dom_vtx,
                                   dom_intrf->comm);

  dom_intrf->is_result[PDM_BOUND_TYPE_VTX] = 1;

  // Simple case is not yet managed, free working arrays
  if (dom_intrf->multizone_intrf == PDM_DOMAIN_INTERFACE_MULT_NO) {
    for (int i_intrf = 0; i_intrf < dom_intrf->n_interface; i_intrf++) {
      for (int j = 0; j < dom_intrf->interface_dn_vtx[i_intrf]; j++) {
        assert(dom_intrf->interface_dom_vtx[i_intrf][2*j]   == dom_intrf->interface_dom_face[i_intrf][0]);
        assert(dom_intrf->interface_dom_vtx[i_intrf][2*j+1] == dom_intrf->interface_dom_face[i_intrf][1]);
      }
      free(dom_intrf->interface_dom_vtx[i_intrf]);
      free(_interface_dom_face[i_intrf]);
    }
    free(_interface_dom_face);
    free(dom_intrf->interface_dom_vtx);
    dom_intrf->interface_dom_vtx = dom_intrf->interface_dom_face;
  }
}

void PDM_domain_interface_get
(
 PDM_domain_interface_t *dom_intrf,
 PDM_bound_type_t        interface_kind,
 int                   **interface_dn,
 PDM_g_num_t          ***interface_ids,
 int                  ***interface_dom
)
{
  assert (dom_intrf != NULL);
  if (interface_kind == PDM_BOUND_TYPE_FACE) {
    assert (dom_intrf->interface_dn_face != NULL);
    *interface_dn  = dom_intrf->interface_dn_face;
    *interface_ids = dom_intrf->interface_ids_face;
    *interface_dom = dom_intrf->interface_dom_face;
  }
  else if (interface_kind == PDM_BOUND_TYPE_VTX) {
    assert (dom_intrf->interface_dn_vtx != NULL);
    *interface_dn  = dom_intrf->interface_dn_vtx;
    *interface_ids = dom_intrf->interface_ids_vtx;
    *interface_dom = dom_intrf->interface_dom_vtx;
  }
  else  {
    PDM_error(__FILE__, __LINE__, 0, "This kind of entity is not yet supported\n");
  }
}

void PDM_domain_interface_free
(
 PDM_domain_interface_t *dom_intrf
)
{
  assert (dom_intrf != NULL);
  if (dom_intrf->ownership == PDM_OWNERSHIP_KEEP) {
    if (dom_intrf->is_result[PDM_BOUND_TYPE_VTX]) {
      for (int i_interface = 0; i_interface < dom_intrf->n_interface; i_interface++) {
        free(dom_intrf->interface_ids_vtx[i_interface]);
        if (dom_intrf->multizone_intrf == PDM_DOMAIN_INTERFACE_MULT_YES)
          free(dom_intrf->interface_dom_vtx[i_interface]);
      }
      free(dom_intrf->interface_dn_vtx);
      free(dom_intrf->interface_ids_vtx);
      if (dom_intrf->multizone_intrf == PDM_DOMAIN_INTERFACE_MULT_YES)
        free(dom_intrf->interface_dom_vtx);
    }
    if (dom_intrf->is_result[PDM_BOUND_TYPE_FACE]) {
      for (int i_interface = 0; i_interface < dom_intrf->n_interface; i_interface++) {
        free(dom_intrf->interface_ids_face[i_interface]);
        if (dom_intrf->multizone_intrf == PDM_DOMAIN_INTERFACE_MULT_YES)
          free(dom_intrf->interface_dom_face[i_interface]);
      }
      free(dom_intrf->interface_dn_face);
      free(dom_intrf->interface_ids_face);
      if (dom_intrf->multizone_intrf == PDM_DOMAIN_INTERFACE_MULT_YES)
        free(dom_intrf->interface_dom_face);
    }
  }

  free(dom_intrf);
}


#ifdef __cplusplus
}
#endif /* __cplusplus */
