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
#include "pdm_binary_search.h"
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

// Translate interfaces list to a graph representation
static int _interface_to_graph
(
  const int                    n_interface,
  PDM_domain_interface_mult_t  multidomain_intrf,
  int                         *interface_dn,
  PDM_g_num_t                **interface_ids,
  int                        **interface_dom,
  int                        **graph_idx,
  PDM_g_num_t                **graph_ids,
  int                        **graph_dom,
  PDM_MPI_Comm                 comm
)
{
  int n_rank;
  PDM_MPI_Comm_size (comm, &n_rank);

  // Step 0 : retrieve some data. We need (to offset gnums)
  //   - the number of involved blocks
  //   - the max id occuring in each block
  int n_domain = -1;
  if (multidomain_intrf == PDM_DOMAIN_INTERFACE_MULT_YES) {
    int max_domain_loc = 0;
    for (int itrf = 0; itrf < n_interface; itrf++) {
      for (int k = 0; k < 2*interface_dn[itrf]; k++) {
        max_domain_loc = PDM_MAX(max_domain_loc, interface_dom[itrf][k]);
      }
    }
    PDM_MPI_Allreduce(&max_domain_loc, &n_domain, 1, PDM_MPI_INT, PDM_MPI_MAX, comm);
  }
  else {
    for (int itrf = 0; itrf < n_interface; itrf++) {
      n_domain = PDM_MAX(n_domain, interface_dom[itrf][0]);
      n_domain = PDM_MAX(n_domain, interface_dom[itrf][1]);
    }
  }
  n_domain++; //Because domain numbering start at 0

  PDM_g_num_t *max_per_domain_loc = PDM_array_const_gnum(n_domain, 0);
  PDM_g_num_t *max_per_domain     = (PDM_g_num_t *) malloc((n_domain+1) * sizeof(PDM_g_num_t));
  for (int itrf = 0; itrf < n_interface; itrf++) {
    int dom, domopp;
    if (multidomain_intrf == PDM_DOMAIN_INTERFACE_MULT_NO) {
      dom    = interface_dom[itrf][0];
      domopp = interface_dom[itrf][1];
    }
    for (int k = 0; k < interface_dn[itrf]; k++) {
      if (multidomain_intrf == PDM_DOMAIN_INTERFACE_MULT_YES) {
        dom    = interface_dom[itrf][2*k];
        domopp = interface_dom[itrf][2*k+1];
      }
      max_per_domain_loc[dom] = PDM_MAX(max_per_domain_loc[dom], interface_ids[itrf][2*k]);
      max_per_domain_loc[domopp] = PDM_MAX(max_per_domain_loc[domopp], interface_ids[itrf][2*k+1]);
    }
  }
  max_per_domain[0] = 0;
  PDM_MPI_Allreduce(max_per_domain_loc, &max_per_domain[1], n_domain, PDM__PDM_MPI_G_NUM, PDM_MPI_MAX, comm);
  PDM_array_accumulate_gnum(max_per_domain, n_domain+1);
  if (0 == 1)
    PDM_log_trace_array_long(max_per_domain, n_domain+1, "max per domain");
  free(max_per_domain_loc);
  
  // Prepare first PtB with multiple partitions.
  // Use (shifted) ids as gnum and send tuple (shited) id, opp_id
  PDM_g_num_t **interface_ids_shifted = (PDM_g_num_t **) malloc(n_interface * sizeof(PDM_g_num_t*));
  PDM_g_num_t **send_data             = (PDM_g_num_t **) malloc(n_interface * sizeof(PDM_g_num_t*));
  double      **weight                = (double      **) malloc(n_interface * sizeof(double*     ));
  int         **stride_one            = (int         **) malloc(n_interface * sizeof(int*        ));
  int          *interface_dn_twice    = (int          *) malloc(n_interface * sizeof(int         ));
  for (int itrf = 0; itrf < n_interface; itrf++) {
    stride_one[itrf]            = (int         *) malloc(2*interface_dn[itrf]*sizeof(int        ));
    interface_ids_shifted[itrf] = (PDM_g_num_t *) malloc(2*interface_dn[itrf]*sizeof(PDM_g_num_t));
    send_data[itrf]             = (PDM_g_num_t *) malloc(2*interface_dn[itrf]*sizeof(PDM_g_num_t));
    weight[itrf]                = (double      *) malloc(2*interface_dn[itrf]*sizeof(double     ));
    interface_dn_twice[itrf]    = 2*interface_dn[itrf];
    int dom, domopp;
    if (multidomain_intrf == PDM_DOMAIN_INTERFACE_MULT_NO) {
      dom    = interface_dom[itrf][0];
      domopp = interface_dom[itrf][1];
    }
    for (int k = 0; k < interface_dn[itrf]; k++) {
      if (multidomain_intrf == PDM_DOMAIN_INTERFACE_MULT_YES) {
        dom    = interface_dom[itrf][2*k];
        domopp = interface_dom[itrf][2*k+1];
      }
      interface_ids_shifted[itrf][2*k]   = interface_ids[itrf][2*k] + max_per_domain[dom];
      interface_ids_shifted[itrf][2*k+1] = interface_ids[itrf][2*k+1] + max_per_domain[domopp];
      send_data[itrf][2*k]   = interface_ids[itrf][2*k+1] + max_per_domain[domopp];
      send_data[itrf][2*k+1] = interface_ids[itrf][2*k] + max_per_domain[dom];
      weight[itrf][2*k]   = 1.;
      weight[itrf][2*k+1] = 1.;
      stride_one[itrf][2*k]   = 1;
      stride_one[itrf][2*k+1] = 1;
    }
    if (0 == 1) {
      log_trace("Interface %d\n", itrf);
      PDM_log_trace_array_long(interface_ids_shifted[itrf], 2*interface_dn[itrf], "  shifted gnum");
      PDM_log_trace_array_long(send_data[itrf], 2*interface_dn[itrf], "  send");
    }
  }
  
  PDM_part_to_block_t *ptb = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                      PDM_PART_TO_BLOCK_POST_MERGE,
                                                      1.,
                                                      interface_ids_shifted,
                                                      weight,
                                                      interface_dn_twice,
                                                      n_interface,
                                                      comm);
  // Save distribution & gnum from first PtB. We will use it for following PtBs
  int n_gnum = PDM_part_to_block_n_elt_block_get(ptb);
  PDM_g_num_t *distri = (PDM_g_num_t *) malloc((n_rank+1) * sizeof(PDM_g_num_t));
  PDM_g_num_t *gnum   = (PDM_g_num_t *) malloc(n_gnum    * sizeof(PDM_g_num_t));
  memcpy(gnum,   PDM_part_to_block_block_gnum_get(ptb),        n_gnum*sizeof(PDM_g_num_t));
  memcpy(distri, PDM_part_to_block_distrib_index_get(ptb), (n_rank+1)*sizeof(PDM_g_num_t));

  int         *recv_stride = NULL;
  PDM_g_num_t *recv_data   = NULL;
  int n_connected_l = PDM_part_to_block_exch(ptb,
                                             sizeof(PDM_g_num_t),
                                             PDM_STRIDE_VAR_INTERLACED,
                                             -1,
                                             stride_one,
                                   (void **) send_data,
                                             &recv_stride,
                                   (void **) &recv_data);
  if (0 == 1) {
    PDM_log_trace_array_long(PDM_part_to_block_block_gnum_get(ptb), n_gnum, "gnum");
    PDM_log_trace_array_int(recv_stride, n_gnum, "recv stride");
    PDM_log_trace_array_long(recv_data, n_connected_l, "recv data");
  }

  PDM_part_to_block_free(ptb);
  for (int itrf = 0; itrf < n_interface; itrf++) {
    free(stride_one           [itrf]);
    free(send_data            [itrf]);
    free(weight               [itrf]);
    free(interface_ids_shifted[itrf]);
  }
  free(stride_one);
  free(weight);
  free(send_data);
  free(interface_dn_twice);
  free(interface_ids_shifted);

  int n_connected;
  PDM_MPI_Allreduce(&n_connected_l, &n_connected, 1, PDM_MPI_INT, PDM_MPI_SUM, comm);
  // log_trace("Initial size of graph : %d \n", n_connected);


  /* After first exchange, we received for each gnum a list (usually of size one) of connected
   * gnum. The idea is to consider this received list as a partition, and to send to each
   * entity of this partition : the neighbors in the received list (if stride was > 1) + the original gnum.
   *
   * After some iterations this will group the related ids together
  */

  int n_connected_prev = 0;
  while(n_connected_prev != n_connected) {
    // log_trace("\nSize of graph : %d (prev %d), start new it\n", n_connected, n_connected_prev);
    n_connected_prev = n_connected;

    ptb = PDM_part_to_block_create2(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                    PDM_PART_TO_BLOCK_POST_MERGE,
                                    1.,
                                   &recv_data,
                                    distri,
                                   &n_connected_l,
                                    1,
                                    comm);

    
    int *send_stride = (int *) malloc(n_connected_l*sizeof(int));
    int w_idx = 0;
    int n_data = 0;
    for (int k = 0; k < n_gnum; k++) {
      for (int j = 0; j < recv_stride[k]; j++) {
        send_stride[w_idx++] = recv_stride[k];
        n_data += recv_stride[k];
      }
    }
    assert (w_idx == n_connected_l);
    PDM_g_num_t *send_data2 = (PDM_g_num_t *) malloc(n_data*sizeof(PDM_g_num_t));
    w_idx = 0;
    int r_idx = 0;
    for (int k = 0; k < n_gnum; k++) {
      for (int j = 0; j < recv_stride[k]; j++) {
        //Add gnum
        send_data2[w_idx++] = gnum[k];
        //Add others
        for (int i = 0; i < j; i++)
          send_data2[w_idx++] = recv_data[r_idx + i];
        for (int i = j+1; i < recv_stride[k]; i++)
          send_data2[w_idx++] = recv_data[r_idx + i];
      }
      r_idx += recv_stride[k];
    }
    assert (r_idx == n_connected_l);
    assert (w_idx == n_data);
    if (0 == 1) {
      PDM_log_trace_array_int(send_stride, n_connected_l, "  send stride");
      PDM_log_trace_array_long(send_data2, n_data, "  send_data2");
    }

    int         *recv_stride_next = NULL;
    PDM_g_num_t *recv_data_next   = NULL;
    PDM_part_to_block_exch(ptb,
                           sizeof(PDM_g_num_t),
                           PDM_STRIDE_VAR_INTERLACED,
                           -1,
                           &send_stride,
                 (void **) &send_data2,
                           &recv_stride_next,
                 (void **) &recv_data_next);

    free(send_stride);
    free(send_data2);
    PDM_part_to_block_free(ptb);
    
    // Post treat recv data to remove duplicated per gnum and count size of graph
    int start = 0;
    n_connected_l = 0;
    for (int i = 0; i < n_gnum; i ++) {
      int n_unique = PDM_inplace_unique_long(recv_data_next, NULL, start, start+recv_stride_next[i]-1);
      //Compress array at the same time (let meanless data at the end of array)
      memcpy(&recv_data_next[n_connected_l], &recv_data_next[start], n_unique*sizeof(PDM_g_num_t));
      start += recv_stride_next[i];
      recv_stride_next[i] = n_unique;
      n_connected_l += n_unique;
    }

    PDM_MPI_Allreduce(&n_connected_l, &n_connected, 1, PDM_MPI_INT, PDM_MPI_SUM, comm);
    if (0 == 1) {
      PDM_log_trace_array_long(gnum, n_gnum, "  gnum");
      PDM_log_trace_array_int(recv_stride_next, n_gnum, "  recv stride");
      PDM_log_trace_array_long(recv_data_next, n_connected_l, "  recv data");
      log_trace("  Total size of graph : %d \n", n_connected);
    }

    free(recv_stride);
    free(recv_data); // To free after PTB because it was used as lngn
    recv_data = recv_data_next;
    recv_stride = recv_stride_next;
  }
  
  // When iteration are completed, all the connections are known by every id.
  // Last step is to compress the graph and to redistribute it
  // To do that we take for each group of related id the min of it as lngn
  int n_keys = 0;
  n_connected_l = 0;
  int r_idx = 0;
  /* is_key_gr : -1 if gnum is not a key (only min of group is the key),
   *              1 if gnum is a key and gnum is not included in received ids
   *              0 otherwise */
  int *is_key_gr = PDM_array_const_int(n_gnum, 1);
  for (int k = 0; k < n_gnum; k++) {
    for (int j = 0; j < recv_stride[k]; j++) {
      if (recv_data[r_idx+j] == gnum[k]) {
        is_key_gr[k] = 0;
      }
      else if (recv_data[r_idx+j] < gnum[k]) {
        is_key_gr[k] = -1;
        break;
      }
    }
    if (is_key_gr[k] != -1) {
      n_keys++;
      n_connected_l += recv_stride[k] + is_key_gr[k];
    }
    r_idx  += recv_stride[k];
  }

  PDM_g_num_t *lngn_gr        = (PDM_g_num_t *) malloc(n_keys*sizeof(PDM_g_num_t));
  int         *send_stride_gr = (int         *) malloc(n_keys*sizeof(int        ));
  double      *weight_gr      = (double      *) malloc(n_keys*sizeof(double     ));
  PDM_g_num_t *send_data_gr   = (PDM_g_num_t *) malloc(n_connected_l *sizeof(PDM_g_num_t));
  int w_idx = 0;
  int w_idx2 = 0;
  r_idx = 0;
  for (int k = 0; k < n_gnum; k++) {
    if (is_key_gr[k] != -1) {
      lngn_gr[w_idx]        = gnum[k];
      send_stride_gr[w_idx] = recv_stride[k] + is_key_gr[k]; //Include gnum in send data (if needed) so we have directly graph
      weight_gr[w_idx]      = (double) (recv_stride[k] + is_key_gr[k]);
      w_idx++;
      if (is_key_gr[k] == 1) {
        send_data_gr[w_idx2++] = gnum[k];
      }
      memcpy(&send_data_gr[w_idx2], &recv_data[r_idx], recv_stride[k]*sizeof(PDM_g_num_t));
      w_idx2 += recv_stride[k];
    }
    r_idx += recv_stride[k];
  }
  if (0 == 1) {
    log_trace("Build graph\n");
    PDM_log_trace_array_long(lngn_gr, n_keys, "  keys graph");
    PDM_log_trace_array_int(send_stride_gr, n_keys, "  send stride");
    PDM_log_trace_array_long(send_data_gr, n_connected_l, "  send data");
  }
  
  // Data of previous iteration is not usefull anymore
  free(gnum);
  free(distri);
  free(recv_stride);
  free(recv_data);
  
  //TODO In fact we just want to do a block to block, but PTB + weights compute distribution for us
  //and we are lazy
  ptb = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                 PDM_PART_TO_BLOCK_POST_MERGE,
                                 1.,
                                &lngn_gr,
                                &weight_gr,
                                &n_keys,
                                 1,
                                 comm);
  int graph_dn = PDM_part_to_block_n_elt_block_get(ptb);

  int         *graph_size = NULL;
  PDM_g_num_t *graph_gnum = NULL;
  PDM_part_to_block_exch(ptb,
                         sizeof(PDM_g_num_t),
                         PDM_STRIDE_VAR_INTERLACED,
                         -1,
                         &send_stride_gr,
               (void **) &send_data_gr,
                         &graph_size,
               (void **) &graph_gnum);

  PDM_part_to_block_free(ptb);
  free(send_stride_gr);
  free(send_data_gr);
  free(weight_gr);
  free(lngn_gr);
  free(is_key_gr);

  int* _graph_idx = PDM_array_new_idx_from_sizes_int(graph_size, graph_dn);
  int *_graph_dom = (int *) malloc(_graph_idx[graph_dn]*sizeof(int));

  if (0 == 1) {
    PDM_log_trace_array_int(graph_size, graph_dn, "  recv stride");
    PDM_log_trace_array_long(graph_gnum, _graph_idx[graph_dn], "  recv data");
  }

  // Retrieve domain and local gnum in domain
  for (int i = 0; i < _graph_idx[graph_dn]; i++) {
    _graph_dom[i] = PDM_binary_search_gap_long(graph_gnum[i]-1, max_per_domain, n_domain+1);
    graph_gnum[i] -= max_per_domain[_graph_dom[i]];
  }
  free(graph_size);
  free(max_per_domain);

  *graph_idx = _graph_idx;
  *graph_ids =  graph_gnum;
  *graph_dom = _graph_dom;
  return graph_dn;
}

static int _extract_and_shift_jn_faces
(
 int           n_domain,
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

  PDM_g_num_t *face_per_block_offset = _per_block_offset(n_domain, dn_face, comm);
  PDM_g_num_t *vtx_per_block_offset  = _per_block_offset(n_domain, dn_vtx,  comm);
  
  int n_face_join = 0; // Each interface comes with a pair of faces
  for (int i_interface = 0; i_interface < n_interface; i_interface++) {
    n_face_join += 2*interfaces_size[i_interface];
  }

  PDM_g_num_t *_dextract_face_id_tmp       = (PDM_g_num_t *) malloc(n_face_join * sizeof(PDM_g_num_t));
  // Also transport some data to the extracted faces
  int         *_dextract_face_group_id_tmp = (int *) malloc(n_face_join*sizeof(int));
  int         *_dextract_face_dom_id_tmp   = (int *) malloc(n_face_join*sizeof(int));

  int idx = 0;
  for (int i_interface = 0; i_interface < n_interface; i_interface++) {

    PDM_g_num_t *_interface_ids = interface_face_ids[i_interface]; //Shortcuts
    int         *_interface_dom = interface_domains_ids[i_interface];
    for (int i_pair = 0; i_pair < interfaces_size[i_interface]; i_pair++) {
      int i_domain_cur = _interface_dom[2*i_pair];
      int i_domain_opp = _interface_dom[2*i_pair+1];

      _dextract_face_id_tmp[idx]         = _interface_ids[2*i_pair] + face_per_block_offset[i_domain_cur];
      _dextract_face_dom_id_tmp  [idx]   = i_domain_cur;
      _dextract_face_group_id_tmp[idx++] = i_interface;

      _dextract_face_id_tmp[idx]         = _interface_ids[2*i_pair+1] + face_per_block_offset[i_domain_opp];
      _dextract_face_dom_id_tmp  [idx]   = i_domain_opp;
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
                          PDM_STRIDE_CST_INTERLACED,
                          2,
                          NULL,
                          _dextract_face_id_tmp,
                          NULL,
              (void **)   dextract_face_id);
  PDM_block_to_block_exch(btb,
                          sizeof(int),
                          PDM_STRIDE_CST_INTERLACED,
                          2,
                          NULL,
                          _dextract_face_group_id_tmp,
                          NULL,
              (void **)  dextract_face_group_id);
  PDM_block_to_block_exch(btb,
                          sizeof(int),
                          PDM_STRIDE_CST_INTERLACED,
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
    PDM_log_trace_array_long(face_per_block_offset, n_domain+1, "face_per_block_offset :: ");
    PDM_log_trace_array_long(vtx_per_block_offset,  n_domain+1, "vtx_per_block_offset :: ");
    PDM_log_trace_array_long(*dextract_face_id, n_face_join, "dextract_face_id :: ");
  }

  PDM_g_num_t **all_face_distribution = (PDM_g_num_t **) malloc(n_domain * sizeof(PDM_g_num_t*));
  for (int i_domain = 0; i_domain < n_domain; i_domain++) {
    all_face_distribution[i_domain] = PDM_compute_entity_distribution(comm, dn_face[i_domain]);
  }

  PDM_multi_block_to_part_t *mptb = PDM_multi_block_to_part_create(face_per_block_offset,
                                                                   n_domain,
                                            (const PDM_g_num_t **) all_face_distribution,
                                            (const PDM_g_num_t **) dextract_face_id,
                                                                  &n_face_join,
                                                                   1,
                                                                   comm);
  //Prepare data to send : face -> vtx connectivity 
  int         **face_vtx_n       = (int         **) malloc(n_domain * sizeof(int*));
  PDM_g_num_t **face_vtx_shifted = (PDM_g_num_t **) malloc(n_domain * sizeof(PDM_g_num_t*));
  for (int i_domain = 0; i_domain < n_domain; i_domain++) {
    face_vtx_n[i_domain]       = (int         *) malloc(dn_face[i_domain] * sizeof(int));
    face_vtx_shifted[i_domain] = (PDM_g_num_t *) malloc(dface_vtx_idx[i_domain][dn_face[i_domain]] * sizeof(PDM_g_num_t));

    for (int i_face = 0; i_face < dn_face[i_domain]; i_face++) {
      face_vtx_n[i_domain][i_face] = dface_vtx_idx[i_domain][i_face+1] - dface_vtx_idx[i_domain][i_face];
      for (int j = dface_vtx_idx[i_domain][i_face]; j < dface_vtx_idx[i_domain][i_face+1]; j++) {
        face_vtx_shifted[i_domain][j] = dface_vtx[i_domain][j] + vtx_per_block_offset[i_domain];
      }
    }
  }

  int         **part_stride = NULL;
  PDM_g_num_t **part_data   = NULL;
  PDM_multi_block_to_part_exch2(mptb,
                                sizeof(PDM_g_num_t),
                                PDM_STRIDE_VAR_INTERLACED,
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

  for (int i_domain = 0; i_domain < n_domain; i_domain++) {
    free(all_face_distribution[i_domain]);
    free(face_vtx_n[i_domain]);
    free(face_vtx_shifted[i_domain]);
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
  PDM_g_num_t *key_ln_to_gn = (PDM_g_num_t *) malloc(dn_internal_edge * sizeof(PDM_g_num_t)); 
  int         *stride_one   = (int         *) malloc(dn_internal_edge * sizeof(int        ));
  int         *stride_two   = (int         *) malloc(dn_internal_edge * sizeof(int        ));
  int         *stride_four  = (int         *) malloc(dn_internal_edge * sizeof(int        ));
  double      *weight       = (double      *) malloc(dn_internal_edge * sizeof(double     ));

  PDM_g_num_t *data_send_connect    = (PDM_g_num_t *) malloc(4*dn_internal_edge * sizeof(PDM_g_num_t));
  PDM_g_num_t *data_send_edge_g_num = (PDM_g_num_t *) malloc(  dn_internal_edge * sizeof(PDM_g_num_t));
  PDM_g_num_t *data_send_group      = (PDM_g_num_t *) malloc(4*dn_internal_edge * sizeof(PDM_g_num_t));
  PDM_g_num_t *data_send_sens       = (PDM_g_num_t *) malloc(2*dn_internal_edge * sizeof(PDM_g_num_t));
  PDM_g_num_t *data_send_face_g_num = (PDM_g_num_t *) malloc(2*dn_internal_edge * sizeof(PDM_g_num_t));

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
    //Retrive domain id using group id of any of two faces
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
                                     PDM_STRIDE_VAR_INTERLACED,
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
                                     PDM_STRIDE_VAR_INTERLACED,
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
                                     PDM_STRIDE_VAR_INTERLACED,
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
                                     PDM_STRIDE_VAR_INTERLACED,
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
  PDM_g_num_t *results_edge     = (PDM_g_num_t *) malloc(gnum_n_occurences_tot * sizeof(PDM_g_num_t));
  PDM_g_num_t *results_edge_opp = (PDM_g_num_t *) malloc(gnum_n_occurences_tot * sizeof(PDM_g_num_t));

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
  results_edge     = (PDM_g_num_t *) realloc(results_edge,     rsvd_gnum_n_occurences_tot*sizeof(PDM_g_num_t));
  results_edge_opp = (PDM_g_num_t *) realloc(results_edge_opp, rsvd_gnum_n_occurences_tot*sizeof(PDM_g_num_t));


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
  *dedge_gnum = (PDM_g_num_t *) malloc(resolved_dn_internal_edge * sizeof(PDM_g_num_t));

  PDM_g_num_t *dedge_gnum_tmp = PDM_part_to_block_block_gnum_get(ptb);
  memcpy(*dedge_gnum, dedge_gnum_tmp, resolved_dn_internal_edge*sizeof(PDM_g_num_t));

  PDM_part_to_block_exch(ptb,
                        sizeof(PDM_g_num_t),
                        PDM_STRIDE_CST_INTERLACED,
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
  PDM_g_num_t *ordered_edge     = (PDM_g_num_t *) malloc(max_face_len * sizeof(PDM_g_num_t));
  PDM_g_num_t *ordered_edge_opp = (PDM_g_num_t *) malloc(max_face_len * sizeof(PDM_g_num_t));
  PDM_g_num_t *ordered_vtx      = (PDM_g_num_t *) malloc(max_face_len * sizeof(PDM_g_num_t));
  PDM_g_num_t *ordered_vtx_opp  = (PDM_g_num_t *) malloc(max_face_len * sizeof(PDM_g_num_t));

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
                                     PDM_STRIDE_VAR_INTERLACED,
                                     -1,
                           (int  **) &stride_one,
                           (void **) &p_all_vtx_opp,
                                     &recv_stride,
                           (void **) &dall_vtx_opp);

  int *unused_recv_stride = NULL;
  exch_size = PDM_part_to_block_exch(ptb,
                                     sizeof(int),
                                     PDM_STRIDE_VAR_INTERLACED,
                                     -1,
                           (int  **) &stride_one,
                           (void **) &p_all_vtx_group,
                                     &unused_recv_stride, //Same  than recv stride
                           (void **) &dall_vtx_group);
  free(unused_recv_stride);

  exch_size = PDM_part_to_block_exch(ptb,
                                     sizeof(int),
                                     PDM_STRIDE_VAR_INTERLACED,
                                     -1,
                           (int  **) &stride_one,
                           (void **) &p_all_vtx_dom_id,
                                     &unused_recv_stride, //Same  than recv stride
                           (void **) &dall_vtx_dom);
  free(unused_recv_stride);

  exch_size = PDM_part_to_block_exch(ptb,
                                     sizeof(int),
                                     PDM_STRIDE_VAR_INTERLACED,
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
    interface_vtx_ids[i_interface]     = (PDM_g_num_t *) malloc(2*vtx_interface_size[i_interface]*sizeof(PDM_g_num_t));
    interface_vtx_dom_ids[i_interface] = (int         *) malloc(2*vtx_interface_size[i_interface]*sizeof(int));
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
    interface_vtx_ids[i_interface] = (PDM_g_num_t *) realloc(interface_vtx_ids[i_interface], 2*n_pairs_u*sizeof(PDM_g_num_t));
    interface_vtx_dom_ids[i_interface] = (int *) realloc(interface_vtx_dom_ids[i_interface], 2*n_pairs_u*sizeof(int));
  }

  PDM_part_to_block_free(ptb);
  free(recv_stride);
  free(dall_vtx_dom);
  free(dall_vtx_dom_opp);
  free(dall_vtx_opp);
  free(dall_vtx_group);
}


static void _connect_additional_edges
(
 int           n_extr_face,
 int          *face_vtx_both_idx,
 PDM_g_num_t  *face_vtx_both,
 PDM_g_num_t  *dface_edge,
 PDM_g_num_t  *dextract_face_join,
 PDM_g_num_t  *pedge_vtx,
 int          *face_status,
 PDM_g_num_t  *face_edge_wopp,
 PDM_MPI_Comm  comm
)
{
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  //TODO : equilibrate with weights
  //1. Use face_vtx connectivity to send to the vertex the id of face to which they belong
  PDM_part_to_block_t *ptb_vtx = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                          PDM_PART_TO_BLOCK_POST_MERGE,
                                                          1.,
                                                         &face_vtx_both,
                                                          NULL,
                                                         &face_vtx_both_idx[n_extr_face],
                                                          1,
                                                          comm);

  PDM_g_num_t *send_data = (PDM_g_num_t *) malloc(2*face_vtx_both_idx[n_extr_face]*sizeof(PDM_g_num_t));
  for (int i_face = 0; i_face < n_extr_face / 2; i_face++) { //Split loop because of face,face_opp order
    for (int i_vtx = face_vtx_both_idx[2*i_face]; i_vtx < face_vtx_both_idx[2*i_face+1]; i_vtx++) {
      send_data[2*i_vtx] = dextract_face_join[2*i_face];
      send_data[2*i_vtx+1] = dextract_face_join[2*i_face+1];
    }
    for (int i_vtx = face_vtx_both_idx[2*i_face+1]; i_vtx < face_vtx_both_idx[2*(i_face+1)]; i_vtx++) {
      send_data[2*i_vtx] = dextract_face_join[2*i_face+1];
      send_data[2*i_vtx+1] = dextract_face_join[2*i_face];
    }
  }
  //PDM_log_trace_array_long(send_data, 2*face_vtx_both_idx[n_extr_face], "send data");
  int *stride_two = PDM_array_const_int(face_vtx_both_idx[n_extr_face], 2);

  int *vtx_stride_two = NULL;
  PDM_g_num_t *vtx_face_ids = NULL;
  int n_recv = PDM_part_to_block_exch(ptb_vtx,
                                      sizeof(PDM_g_num_t),
                                      PDM_STRIDE_VAR_INTERLACED,
                                      -1,
                            (int **)  &stride_two,
                           (void **)  &send_data,
                                      &vtx_stride_two,
                            (void **) &vtx_face_ids);

  int n_vtx_blk = PDM_part_to_block_n_elt_block_get(ptb_vtx);
  PDM_g_num_t *vtx_distri = (PDM_g_num_t *) malloc((n_rank+1)*sizeof(PDM_g_num_t));
  memcpy(vtx_distri, PDM_part_to_block_distrib_index_get(ptb_vtx), (n_rank+1)*sizeof(PDM_g_num_t));

  PDM_g_num_t *vtx_gnum = PDM_part_to_block_block_gnum_get(ptb_vtx);
  
  if (0 == 1) {
    PDM_log_trace_array_long(vtx_gnum, n_vtx_blk, "block gnum");
    PDM_log_trace_array_int(vtx_stride_two, n_vtx_blk, "recv stride 2");
    PDM_log_trace_array_long(vtx_face_ids, n_recv, "recv data");
  }

  int *vtx_face_n = (int *) malloc(n_vtx_blk * sizeof(int)); //We received 2 data per face connected to each vtx
  for (int i=0; i < n_vtx_blk; i++) {
    vtx_face_n[i] = vtx_stride_two[i] / 2;
  }

  free(stride_two);
  free(send_data);

  // 2. Build key and send face ids & gnum using key numbering
  PDM_g_num_t *vtx_key = PDM_array_const_gnum(n_vtx_blk, 0);
  int read_idx = 0;
  for (int j=0; j < n_vtx_blk; j++) {
    for (int k = 0; k < vtx_stride_two[j]; k++)
      vtx_key[j] += vtx_face_ids[read_idx++];
  }

  // PDM_log_trace_array_long(vtx_key, n_vtx_blk, "vtx_key");

  PDM_part_to_block_t *ptb_key = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                          PDM_PART_TO_BLOCK_POST_MERGE,
                                                          1.,
                                        (PDM_g_num_t **) &vtx_key,
                                                          NULL,
                                                         &n_vtx_blk,
                                                          1,
                                                          comm);

  int *stride_one = PDM_array_const_int(n_vtx_blk, 1);

  int         *unused_recv_stride = NULL;
  PDM_g_num_t *key_vtx_gnum       = NULL;
  // n_key_vtx is the number of vertices involved in keys (counted multiple times) managed by this proc
  int n_key_vtx =  PDM_part_to_block_exch(ptb_key,
                                     sizeof(PDM_g_num_t),
                                     PDM_STRIDE_VAR_INTERLACED,
                                     -1,
                           (int **)  &stride_one,
                           (void **) &vtx_gnum,
                                     &unused_recv_stride, //Same as key count
                           (void **) &key_vtx_gnum);
  free(unused_recv_stride);
  int *key_recv_face_n = NULL; // For each key (unmerged), number of face related to key
  PDM_part_to_block_exch(ptb_key,
                         sizeof(int),
                         PDM_STRIDE_VAR_INTERLACED,
                         -1,
               (int **)  &stride_one,
               (void **) &vtx_face_n,
                         &unused_recv_stride,  //Same as key count
               (void **) &key_recv_face_n);
  free(unused_recv_stride);

  int         *key_recv_stride = NULL;
  PDM_g_num_t *key_recv_data = NULL; //For each key (unmerged), tuples (face/face_opp)  * nb of face related to key
  PDM_part_to_block_exch(ptb_key,
                         sizeof(PDM_g_num_t),
                         PDM_STRIDE_VAR_INTERLACED,
                         -1,
               (int **)  &vtx_stride_two,
               (void **) &vtx_face_ids,
                         &key_recv_stride,
               (void **) &key_recv_data);


  free(stride_one);
  PDM_part_to_block_free(ptb_vtx); // Needed until here for vtx gnum

  int n_keys = PDM_part_to_block_n_elt_block_get(ptb_key);
  PDM_g_num_t *keys_ids = PDM_part_to_block_block_gnum_get(ptb_key);
  int         *keys_cnt = PDM_part_to_block_block_gnum_count_get(ptb_key);

  if (0 == 1) {
    PDM_log_trace_array_long(keys_ids, n_keys, "key to treat");
    PDM_log_trace_array_int(keys_cnt, n_keys, "n recept");
    PDM_log_trace_array_int(key_recv_stride, n_keys, "key recv stride");
    PDM_log_trace_array_int(key_recv_face_n, n_key_vtx, "key recv facen"); 
    //PDM_log_trace_array_long(key_recv_data, n_recv, "key recv data"); 
    PDM_log_trace_array_long(key_vtx_gnum, n_key_vtx, "key recv gnum");
  }

  free(vtx_face_ids);
  free(vtx_face_n);
  free(vtx_stride_two);
  free(vtx_key);
  
  //3. Match data on key distribution
  PDM_g_num_t *key_vtx_gnum_opp = PDM_array_const_gnum(n_key_vtx, 0);
  int count_idx = 0; //Start of data in key_recv_face_n
  int data_idx = 0; //Start of data in key_recv_data
  for (int i_key = 0; i_key < n_keys; i_key++) {
    int n_vtx_this_key = keys_cnt[i_key];
    // Each key occurs n_vtx_this_key times, each time with 2 face ids * nb of face connected to vertex
    // First we sort these sections
    for (int k=0; k < n_vtx_this_key; k++) {
      PDM_sort_long(&key_recv_data[data_idx], NULL, 2*key_recv_face_n[count_idx]);
      data_idx += 2*key_recv_face_n[count_idx];
      count_idx++;
    }
    // Now search matches
    int idx1 = data_idx - key_recv_stride[i_key]; //Reset idx1 for comparaison
    count_idx -= n_vtx_this_key;
    for (int k = 0; k < n_vtx_this_key; k++) {
      int n_match = 0;
      int i_match;
      int idx2 = data_idx - key_recv_stride[i_key]; //Reset idx2 for comparaison

      for (int k2 = 0; k2 < n_vtx_this_key; k2++) {
        if (k2 != k) { //Skip myself
          if (key_recv_face_n[count_idx+k] == key_recv_face_n[count_idx+k2]) {
            if (PDM_array_are_equal_gnum(&key_recv_data[idx1], &key_recv_data[idx2], 2*key_recv_face_n[count_idx+k])) {
              n_match++;
              i_match = k2;
            }
          }
        }
        idx2 += 2*key_recv_face_n[count_idx+k2];
      }
      if (n_match == 1) { //Register match
        key_vtx_gnum_opp[count_idx + k] = key_vtx_gnum[count_idx+i_match];
      }
      idx1 += 2*key_recv_face_n[count_idx+k];
    }
    count_idx += n_vtx_this_key;
  }

  PDM_part_to_block_free(ptb_key);
  free(key_recv_stride);
  free(key_recv_data);
  free(key_recv_face_n);

  // 4. We send back the matches to vertex distribution to have block property
  PDM_part_to_block_t *ptb_vtx2 = PDM_part_to_block_create2(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                            PDM_PART_TO_BLOCK_POST_NOTHING,
                                                            1.,
                                          (PDM_g_num_t **) &key_vtx_gnum,
                                                            vtx_distri,
                                                           &n_key_vtx,
                                                            1,
                                                            comm);
  assert (PDM_part_to_block_n_elt_block_get(ptb_vtx2) == n_vtx_blk);
  PDM_g_num_t *matched_gnum = (PDM_g_num_t *) malloc(n_vtx_blk*sizeof(PDM_g_num_t));
  memcpy(matched_gnum, PDM_part_to_block_block_gnum_get(ptb_vtx2), n_vtx_blk*sizeof(PDM_g_num_t));

  PDM_g_num_t *matched_gnum_opp = NULL;
  PDM_part_to_block_exch(ptb_vtx2,
                         sizeof(PDM_g_num_t),
                         PDM_STRIDE_CST_INTERLACED,
                         1,
                         NULL,
               (void **) &key_vtx_gnum_opp,
                         NULL,
               (void **) &matched_gnum_opp);

  PDM_part_to_block_free(ptb_vtx2);
  free(key_vtx_gnum);
  free(key_vtx_gnum_opp);


  // 5a. Prepare edge matching : get the vtx_gnum_opp only for the untreated faces
  int unsolvable_edge = 0;
  for (int i_face=0; i_face < n_extr_face; i_face++) {
    if (face_status[i_face] == 0)
      unsolvable_edge += face_vtx_both_idx[i_face+1] - face_vtx_both_idx[i_face];
  }

  PDM_g_num_t *requested_gnum = (PDM_g_num_t *) malloc(unsolvable_edge*sizeof(PDM_g_num_t));
  int idx = 0;
  for (int i_face = 0; i_face < n_extr_face; i_face++) {
    if (face_status[i_face] == 0) {
      int n_vtx_face = face_vtx_both_idx[i_face+1] - face_vtx_both_idx[i_face];
      memcpy(&requested_gnum[idx], &face_vtx_both[face_vtx_both_idx[i_face]], n_vtx_face*sizeof(PDM_g_num_t));
      idx += n_vtx_face;
    }
  }
  PDM_block_to_part_t *btp = PDM_block_to_part_create(vtx_distri,
                               (const PDM_g_num_t **) &requested_gnum,
                                                      &unsolvable_edge,
                                                      1,
                                                      comm);
  int *blk_stride = PDM_array_zeros_int(vtx_distri[i_rank+1] - vtx_distri[i_rank]);
  for (int i = 0; i < n_vtx_blk; i++) {
    blk_stride[matched_gnum[i] - vtx_distri[i_rank] - 1] = 1;
  }
  int         **recv_stride;
  PDM_g_num_t **recv_data;
  PDM_block_to_part_exch2(btp,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_VAR_INTERLACED,
                          blk_stride,
                          matched_gnum_opp, 
                         &recv_stride,
              (void ***) &recv_data);
  free(recv_stride[0]);
  free(recv_stride);
  PDM_g_num_t *requested_gnum_opp = recv_data[0];
  free(recv_data);
  free(blk_stride);
  PDM_block_to_part_free(btp);

  free(matched_gnum_opp);
  free(matched_gnum);
  free(vtx_distri);

  //5b. Now we can match edges and update face_edge_wopp !
  int extracted_face_id = 0;
  for (int i_face = 0; i_face < n_extr_face/2; i_face++) {
    if (face_status[2*i_face] == 0) {
      int edge_start     = face_vtx_both_idx[2*i_face];
      int edge_opp_start = face_vtx_both_idx[2*i_face+1];
      int n_vtx_face = face_vtx_both_idx[2*i_face+1] - face_vtx_both_idx[2*i_face];
      
      //Find any opp vtx gnum to have a gnum / gnum opp couple
      PDM_g_num_t opp_vtx_gnum = 0;
      int         opp_vtx_pos  = 0;
      for (int j=0; j < n_vtx_face; j++) {
        opp_vtx_gnum = requested_gnum_opp[extracted_face_id + j];
        if (opp_vtx_gnum != 0) {
          opp_vtx_pos = j;
          break;
        }
      }
      assert (opp_vtx_gnum != 0);
      PDM_g_num_t my_vtx_gnum = requested_gnum[extracted_face_id+opp_vtx_pos];

      // Search my vtx gnum in edges
      int edge_pos = -1;
      int edge_sens;
      for (int j = 0; j < n_vtx_face; j++) {
        if (pedge_vtx[2*edge_start + 2*j] == my_vtx_gnum) {
          edge_pos = j;
          edge_sens = 0;
          break;
        }
        else if (pedge_vtx[2*edge_start + 2*j + 1] == my_vtx_gnum) {
          edge_pos = j;
          edge_sens = 1;
          break;
        }
      }
      assert (edge_pos >= 0);
      
      //Search opp vtx in edge opp, reverse order
      int opp_pos = -1;
      int opp_sens = 1 - edge_sens;
      for (int j = 0; j < n_vtx_face; j++) {
        if (pedge_vtx[2*edge_opp_start + 2*j + opp_sens] == opp_vtx_gnum) {
          opp_pos = j;
          break;
        }
      }
      assert (opp_pos >= 0);

      //Complete face_edge_wopp for face & opp at the same time
      face_edge_wopp[edge_start + edge_pos] = dface_edge[edge_opp_start + opp_pos];
      face_edge_wopp[edge_opp_start + opp_pos] = dface_edge[edge_start + edge_pos];

      face_status[2*i_face]   = 1;
      face_status[2*i_face+1] = 1;
      extracted_face_id += 2*n_vtx_face;
    }
  }

  free(requested_gnum);
  free(requested_gnum_opp);
}



static void _domain_interface_face_to_vertex
(
 int            n_interface,             /* Total number of interfaces */
 int           *interfaces_size,         /* Number of face pairs in each interface */
 PDM_g_num_t  **interface_face_ids,      /* For each interface, list of pairs face,face_opp */
 int          **interface_domains_ids,   /* For each interface, list of domains dom,dom_opp */
 int            n_domain,                  /* Number of domains */
 int           *dn_vtx,                  /* Number of vertex in each domain (distributed) */
 int           *dn_face,                 /* Number of face in each domain (distributed) */
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

  PDM_g_num_t *face_per_block_offset = _per_block_offset(n_domain, dn_face, comm);
  PDM_g_num_t *vtx_per_block_offset  = _per_block_offset(n_domain, dn_vtx,  comm);

  int         *face_vtx_both_idx = NULL;
  PDM_g_num_t *face_vtx_both     = NULL;
  int         *dextract_face_dom_id   = NULL;
  int         *dextract_face_group_id = NULL;
  PDM_g_num_t *dextract_face_join     = NULL;
  int n_extr_face = _extract_and_shift_jn_faces(n_domain,
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
  PDM_g_num_t *dextract_face_join_opp   = (PDM_g_num_t *) malloc(n_extr_face*sizeof(PDM_g_num_t));
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


  PDM_g_num_t *dedge_face_join         = (PDM_g_num_t *) malloc(dedge_face_idx[dn_edge] * sizeof(PDM_g_num_t));
  PDM_g_num_t *dedge_face_join_opp     = (PDM_g_num_t *) malloc(dedge_face_idx[dn_edge] * sizeof(PDM_g_num_t));
  PDM_block_to_part_t *btp = PDM_block_to_part_create(extracted_face_distri,
                               (const PDM_g_num_t **) &dedge_face_abs,
                                                      &dedge_face_idx[dn_edge],
                                                      1,
                                                      comm);
  int cst_stride = 1;
  PDM_block_to_part_exch(btp,
                         sizeof(PDM_g_num_t),
                         PDM_STRIDE_CST_INTERLACED,
                         &cst_stride,
                (void *) dextract_face_join,
                         NULL,
             (void ** ) &dedge_face_join);

  PDM_block_to_part_exch(btp,
                         sizeof(PDM_g_num_t),
                         PDM_STRIDE_CST_INTERLACED,
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
                          PDM_STRIDE_VAR_INTERLACED,
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
                          PDM_STRIDE_CST_INTERLACED,
                         &stride2,
                          dedge_vtx,
                          NULL,
              (void ***) &recv_data_tmp);
  PDM_g_num_t *pedge_vtx = recv_data_tmp[0];
  free(recv_data_tmp);
  //PDM_log_trace_array_long(pedge_vtx, 2*dface_edge_idx[n_extr_face], "pedge_vtx");



  //Attention, on devrait pouvoir travailler sur face externes uniquement (filtre le dface_edge_abs)
  PDM_g_num_t *face_edge_wopp = (PDM_g_num_t *) malloc(dface_edge_idx[n_extr_face]*sizeof(PDM_g_num_t));
  idx = 0;
  for (int i = 0; i < dface_edge_idx[n_extr_face]; i++) {
    if (pedge_gnum_n[i] == 1)
      face_edge_wopp[i] = pedge_gnum_opp[idx++];
    else
      face_edge_wopp[i] = 0;
  }

  //PDM_log_trace_connectivity_long(dface_edge_idx, face_edge_wopp, n_extr_face, "dface_edge :: ");
  int *face_status = (int *) malloc(n_extr_face*sizeof(int));
  for (int i_face=0; i_face < n_extr_face; i_face++) {
    int n_treated_edge = 0;
    for (int i_edge = dface_edge_idx[i_face]; i_edge < dface_edge_idx[i_face+1]; i_edge++)
      n_treated_edge += pedge_gnum_n[i_edge];
    face_status[i_face] = (int) (n_treated_edge > 0);
  }


  free(dface_edge_abs);
  free(dedge_gnum_n);
  free(pedge_gnum_n);
  free(pedge_gnum_opp);
  PDM_block_to_part_free(btp);

  int need_more_edge_l = 0;
  int need_more_edge;
  for (int i_face=0; i_face < n_extr_face; i_face++) {
    if (face_status[i_face] == 0) {
      need_more_edge_l = 1; 
      break;
    }
  }
  PDM_MPI_Allreduce(&need_more_edge_l, &need_more_edge, 1, PDM_MPI_INT, PDM_MPI_SUM, comm);
  if (need_more_edge > 0) {
    log_trace("Warning -- Some face have not shared edges. Try to retrieve it using vertex connectivity\n");
    _connect_additional_edges(n_extr_face,
                              face_vtx_both_idx,
                              face_vtx_both,
                              dface_edge,
                              dextract_face_join,
                              pedge_vtx,
                              face_status,
                              face_edge_wopp,
                              comm);
  }
  free(face_status);

  //Match external edges
  assert (dface_edge_idx[n_extr_face] % 2 == 0);
  int n_vtx_interface_tot = dface_edge_idx[n_extr_face] / 2;
  PDM_g_num_t *p_all_vtx      = (PDM_g_num_t *) malloc(n_vtx_interface_tot * sizeof(PDM_g_num_t));
  PDM_g_num_t *p_all_vtx_opp  = (PDM_g_num_t *) malloc(n_vtx_interface_tot * sizeof(PDM_g_num_t));
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
  int *p_all_vtx_group     = (int *) malloc(n_vtx_interface_tot * sizeof(int));
  int *p_all_vtx_dom_id    = (int *) malloc(n_vtx_interface_tot * sizeof(int));
  int *p_all_vtx_domopp_id = (int *) malloc(n_vtx_interface_tot * sizeof(int));
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


  //Ultimate step : go back to original vtx numbering. All we have to do is retrieve domain
  // and substract domain offset
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
const int                   n_domain,
PDM_domain_interface_mult_t multidomain_interface,
PDM_ownership_t             ownership,
PDM_MPI_Comm                comm
)
{
  PDM_domain_interface_t *dom_intrf = (PDM_domain_interface_t *) malloc (sizeof(PDM_domain_interface_t));
  dom_intrf->n_interface       = n_interface;
  dom_intrf->n_domain            = n_domain;
  dom_intrf->multidomain_intrf   = multidomain_interface;
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
  if (interface_kind == PDM_BOUND_TYPE_VTX) {
    dom_intrf->interface_dn_vtx  = interface_dn;
    dom_intrf->interface_ids_vtx = interface_ids;
    dom_intrf->interface_dom_vtx = interface_dom;
  }
  else if (interface_kind == PDM_BOUND_TYPE_FACE) {
    dom_intrf->interface_dn_face  = interface_dn;
    dom_intrf->interface_ids_face = interface_ids;
    dom_intrf->interface_dom_face = interface_dom;
  }
  else {
    PDM_error(__FILE__, __LINE__, 0, "Kind of interface not supported\n");
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
  if (dom_intrf->multidomain_intrf == PDM_DOMAIN_INTERFACE_MULT_NO) {
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
                                   dom_intrf->n_domain,
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
  if (dom_intrf->multidomain_intrf == PDM_DOMAIN_INTERFACE_MULT_NO) {
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

int PDM_domain_interface_get_as_graph
(
 PDM_domain_interface_t *dom_intrf,
 PDM_bound_type_t        interface_kind,
 int                   **interface_graph_idx,
 PDM_g_num_t           **interface_graph_ids,
 int                   **interface_graph_dom
)
{
  assert (dom_intrf != NULL);
  int          *interface_dn = NULL;
  PDM_g_num_t **interface_ids = NULL;
  int         **interface_dom = NULL;
  if (interface_kind == PDM_BOUND_TYPE_FACE) {
    assert (dom_intrf->interface_dn_face != NULL);
    interface_dn  = dom_intrf->interface_dn_face;
    interface_ids = dom_intrf->interface_ids_face;
    interface_dom = dom_intrf->interface_dom_face;
  }
  else if (interface_kind == PDM_BOUND_TYPE_VTX) {
    assert (dom_intrf->interface_dn_vtx != NULL);
    interface_dn  = dom_intrf->interface_dn_vtx;
    interface_ids = dom_intrf->interface_ids_vtx;
    interface_dom = dom_intrf->interface_dom_vtx;
  }
  else  {
    PDM_error(__FILE__, __LINE__, 0, "This kind of entity is not yet supported\n");
  }

  int graph_dn = _interface_to_graph(dom_intrf->n_interface,
                                     dom_intrf->multidomain_intrf,
                                     interface_dn,
                                     interface_ids,
                                     interface_dom,
                                     interface_graph_idx,
                                     interface_graph_ids,
                                     interface_graph_dom,
                                     dom_intrf->comm);
  return graph_dn;
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
        if (dom_intrf->multidomain_intrf == PDM_DOMAIN_INTERFACE_MULT_YES)
          free(dom_intrf->interface_dom_vtx[i_interface]);
      }
      free(dom_intrf->interface_dn_vtx);
      free(dom_intrf->interface_ids_vtx);
      if (dom_intrf->multidomain_intrf == PDM_DOMAIN_INTERFACE_MULT_YES)
        free(dom_intrf->interface_dom_vtx);
    }
    if (dom_intrf->is_result[PDM_BOUND_TYPE_FACE]) {
      for (int i_interface = 0; i_interface < dom_intrf->n_interface; i_interface++) {
        free(dom_intrf->interface_ids_face[i_interface]);
        if (dom_intrf->multidomain_intrf == PDM_DOMAIN_INTERFACE_MULT_YES)
          free(dom_intrf->interface_dom_face[i_interface]);
      }
      free(dom_intrf->interface_dn_face);
      free(dom_intrf->interface_ids_face);
      if (dom_intrf->multidomain_intrf == PDM_DOMAIN_INTERFACE_MULT_YES)
        free(dom_intrf->interface_dom_face);
    }
  }

  free(dom_intrf);
}


#ifdef __cplusplus
}
#endif /* __cplusplus */
