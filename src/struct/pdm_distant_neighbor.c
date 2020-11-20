/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_priv.h"
#include "pdm_mpi.h"
#include "pdm_logging.h"
#include "pdm_printf.h"
#include "pdm_distant_neighbor_priv.h"
#include "pdm_distant_neighbor.h"
#include "pdm_error.h"
#include "pdm_order.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Type
 *============================================================================*/

/*============================================================================
 * Local structure definitions
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
void
_compute_unique_idx
(
 int order[],
 int order_unique[],
 int connect_triplet[],
 const int nb_ent
)
{
  int idx_unique = -1;
  int last_proc  = -1;
  int last_part  = -1;
  int last_elmt  = -1;

  for(int i = 0; i < nb_ent; i++){

    int old_order = order[i];
    int curr_proc = connect_triplet[3*old_order  ];
    int curr_part = connect_triplet[3*old_order+1];
    int curr_elmt = connect_triplet[3*old_order+2];
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
    order_unique[i] = idx_unique;
    // log_trace("[%d] = %d --> %d \n", i, is_same, idx_unique);
  }


  if(0 == 1){
    log_trace("order_unique:: \n");
    for(int i = 0; i < nb_ent; i++){
      log_trace(" -------------------------- \n");
      // int pos_unique = order_unique_j1[i];
      // int curr_idx   = order_j1[pos_unique];
      int pos_unique = order_unique[i];
      int curr_idx   = order[i];

      int curr_proc  = connect_triplet[3*curr_idx  ];
      int curr_part  = connect_triplet[3*curr_idx+1];
      int curr_elmt  = connect_triplet[3*curr_idx+2];

      log_trace("\t pos_unique :: %d \n", pos_unique);
      log_trace("\t curr_idx   :: %d \n", curr_idx  );
      log_trace("\t triplet    :: ( %d / %d / %d ) \n", curr_proc, curr_part, curr_elmt);

    }
    log_trace("\n");
  }
}

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/**
 * \brief Return an initialized \ref PDM_distant_neighbor_t structure
 *
 * This function returns an initialized \ref PDM_distant_neighbor_t structure
 *
 * \param [in]   comm          MPI communicator
 * \param [in]   n_part        Number of partitions
 * \param [in]   n_entity      Number of entities for each partition
 * \param [out]  neighbor_idx  Indexes of candidate for each current part point
 *                              (size = number of entities in the current part + 1)
 * \param [out]  neighbor_desc Candidates description (process,
 *                                                     part in the process,
 *                                                     entitiy in the part)
 *
 * \return      A new initialized \ref PDM_distant_neighbor structure
 *
 */
PDM_distant_neighbor_t*
PDM_distant_neighbor_create
(
const PDM_MPI_Comm   comm,
const int            n_part,
const int           *n_entity,
      int          **neighbor_idx,
      int          **neighbor_desc
)
{
  PDM_distant_neighbor_t *dn = (PDM_distant_neighbor_t *) malloc(sizeof(PDM_distant_neighbor_t));

  dn->comm          = comm;

  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(dn->comm, &i_rank);
  PDM_MPI_Comm_size(dn->comm, &n_rank);

  dn->n_part               = n_part;
  dn->n_entity             = n_entity;
  dn->neighbor_idx         = neighbor_idx;
  dn->neighbor_desc        = neighbor_desc;
  dn->order                = (int **) malloc(   dn->n_part       * sizeof(int **));
  dn->order_unique         = (int **) malloc(   dn->n_part       * sizeof(int **));
  dn->requested_data_n     = (int * ) malloc( ( n_rank          ) * sizeof(int * ));
  dn->requested_data_idx   = (int * ) malloc( ( n_rank + 1      ) * sizeof(int * ));
  dn->distributed_part_n   = (int * ) malloc( ( dn->n_part     ) * sizeof(int * ));
  dn->distributed_part_idx = (int * ) malloc( ( dn->n_part + 1 ) * sizeof(int * ));

  /*
   * Init the requested counter
   */
  for (int i = 0; i < n_rank; i++) {
    dn->requested_data_idx[i] = 0;
    dn->requested_data_n  [i] = 0;
  }

  /*
   * Sort/unique the triplet (iproc, i_part, ientity)
   * Attention on doit faire le tri pour chaque partttion
   *    --> Dans le cas des vertex on aurai 2 joins qui demande 2 fois le même vertex et le tri
   *       --> Il faudrait copier dans un tableau 1D les pattions les une à la suite des autres !
   *           On change l'interface ?
   *       --> En concatemant avant le MPI est plus simple car on a pas à creer un
   *           tableau d'index sur chaque proc de deplacement entre les différentes partitions
   *       Les 2 sont pratiques car on pourrait imaginer avoir des traitement locaux differents non ?
   */

  for(int i_part = 0; i_part < dn->n_part; i_part++){

    int *_part_neighbor_idx  = dn->neighbor_idx [i_part];
    int *_part_neighbor_desc = dn->neighbor_desc[i_part];

    // printf("[%i] - n_entity:: %d\n", i_part, n_entity[i_part]);

    dn->order       [i_part] = (int *) malloc( _part_neighbor_idx[n_entity[i_part]] * sizeof(int *));
    dn->order_unique[i_part] = (int *) malloc( _part_neighbor_idx[n_entity[i_part]] * sizeof(int *));

    // Sort
    PDM_order_lnum_s(dn->neighbor_desc[i_part],
                     3,
                     dn->order[i_part],
                     _part_neighbor_idx[n_entity[i_part]]);

    // Compute the unique idx from sort
    _compute_unique_idx(dn->order[i_part],
                        dn->order_unique[i_part],
                        dn->neighbor_desc[i_part],
                        _part_neighbor_idx[n_entity[i_part]]);

    // Il faut connaitre le nombre d'occurence une fois trié --> Taille du buffer d'envoie
    int lastidx = -1;
    for(int i_entity = 0; i_entity < _part_neighbor_idx[n_entity[i_part]]; i_entity++){
      int u_entity = dn->order_unique[i_part][i_entity];
      int s_entity = dn->order[i_part][i_entity];
      // printf("[%d] - order:: %d | unique:: %d | lastidx:: %d \n", i_entity, s_entity, u_entity, lastidx);
      if(lastidx != u_entity){
        int opp_proc = _part_neighbor_desc[3*s_entity];
        dn->requested_data_n[opp_proc]++;
        lastidx = u_entity;
      }
    }
  }

  for (int i = 0; i < n_rank; i++) {
    dn->requested_data_idx[i+1] = dn->requested_data_idx[i] +
                                   dn->requested_data_n[i];
  }

  /*
   * Compute size and reset counter
   */
  int s_requested_data = dn->requested_data_idx[n_rank-1] + dn->requested_data_n[n_rank-1];

  // printf("s_requested_data:: %d \n", s_requested_data);
  for (int i = 0; i < n_rank; i++) {
    dn->requested_data_n[i] = 0;
  }

  int *requested_data = malloc (sizeof(int) *  2 * s_requested_data); // Store i_part/ientity

  for(int i_part = 0; i_part < dn->n_part; i_part++){

    int *_part_neighbor_idx  = dn->neighbor_idx[i_part];
    int *_part_neighbor_desc = dn->neighbor_desc[i_part];
    dn->distributed_part_n  [i_part] = 0;
    dn->distributed_part_idx[i_part] = 0;

    int lastidx = -1;
    for(int i_entity = 0; i_entity < _part_neighbor_idx[n_entity[i_part]]; i_entity++){
      int u_entity = dn->order_unique[i_part][i_entity];
      int s_entity = dn->order[i_part][i_entity];
      // printf("[%d] - order:: %d | unique:: %d | lastidx:: %d \n", i_entity, s_entity, u_entity, lastidx);
      if(lastidx != u_entity){

        int opp_proc = _part_neighbor_desc[3*s_entity  ];
        int opp_part = _part_neighbor_desc[3*s_entity+1];
        int opp_etty = _part_neighbor_desc[3*s_entity+2];

        int idx = dn->requested_data_idx[opp_proc] + dn->requested_data_n[opp_proc]++;

        requested_data[2*idx  ] = opp_part;
        requested_data[2*idx+1] = opp_etty;

        dn->distributed_part_n[i_part]++;

        lastidx = u_entity;
      }
    }
  }

  // Each pacquet have 2 value so we multiply by 2 temporary
  for (int i = 0; i < n_rank; i++) {
    dn->requested_data_n[i] = 2*dn->requested_data_n[i];
  }

  for (int i = 0; i < dn->n_part; i++) {
    dn->distributed_part_idx[i+1] = dn->distributed_part_n[i] +
                                     dn->distributed_part_idx[i];
  }

  if(0 == 1){
    log_trace("PDM_distant_neighbor_create::requested_data :: --> ");
    for(int i = 0; i < s_requested_data; ++i){
      log_trace("[%d/%d] ", requested_data[2*i], requested_data[2*i+1]);
    }
    log_trace("\n");
  }

  if(0 == 1){
    log_trace("PDM_distant_neighbor_create::distributed_part_n :: --> ");
    for(int i = 0; i < dn->n_part; ++i){
      log_trace("%d ", dn->distributed_part_n[i]);
    }
    log_trace("\n");
    log_trace("PDM_distant_neighbor_create::distributed_part_idx :: --> ");
    for(int i = 0; i < dn->n_part+1; ++i){
      log_trace("%d ", dn->distributed_part_idx[i]);
    }
    log_trace("\n");
  }

  /*
   * Exchange the requested data
   */
  dn->distributed_data_n = malloc (sizeof(int) * n_rank);

  PDM_MPI_Alltoall (dn->requested_data_n,   1, PDM_MPI_INT,
                    dn->distributed_data_n, 1, PDM_MPI_INT,
                    dn->comm);

  if(0 == 1){
    log_trace("PDM_distant_neighbor_create::distributed_data_n :: --> ");
    for(int i = 0; i < n_rank; ++i){
      log_trace("%d ",  dn->distributed_data_n[i]);
    }
    log_trace("\n");
  }

  dn->distributed_data_idx = malloc (sizeof(int) * (n_rank + 1));
  dn->distributed_data_idx[0] = 0;

  for (int i = 0; i < n_rank; i++) {
    dn->distributed_data_idx[i+1] = dn->distributed_data_n[i] +
                                     dn->distributed_data_idx[i];
  }

  if(0 == 1){
    log_trace("PDM_distant_neighbor_create::distributed_data_idx :: --> ");
    for(int i = 0; i < n_rank+1; ++i){
      log_trace("%d ",  dn->distributed_data_idx[i]);
    }
    log_trace("\n");
  }

  dn->distributed_data = malloc (sizeof(int) * 2 * dn->distributed_data_idx[n_rank]);

  PDM_MPI_Alltoallv (requested_data,
                     dn->requested_data_n,
                     dn->requested_data_idx,
                     PDM_MPI_INT,
                     dn->distributed_data,
                     dn->distributed_data_n,
                     dn->distributed_data_idx,
                     PDM_MPI_INT,
                     dn->comm);

  if(0 == 1){
    log_trace("PDM_distant_neighbor_create::distributed_data :: --> ");
    for(int i = 0; i < dn->distributed_data_idx[n_rank]/2; ++i){
    // for(int i = 0; i < dn->distributed_data_idx[n_rank]; ++i){
      log_trace("[%d/%d] ", dn->distributed_data[2*i], dn->distributed_data[2*i+1]);
    }
    log_trace("\n");
  }

  /*
   * Store in structure the ordering in the recv buffer
   */

  /*
   * Re setup size of each member of struct
   */
  dn->requested_data_idx[0] = 0;
  dn->distributed_data_idx[0] = 0;
  for (int i = 0; i < n_rank; i++) {
    dn->requested_data_n[i]   = dn->requested_data_n  [i]/2;
    dn->distributed_data_n[i] = dn->distributed_data_n[i]/2;

    dn->requested_data_idx  [i+1] = dn->requested_data_idx[i] + dn->requested_data_n    [i];
    dn->distributed_data_idx[i+1] = dn->distributed_data_n[i] + dn->distributed_data_idx[i];
  }

  if(0 == 1){
    log_trace("Re-Setup --- ");
    log_trace("PDM_distant_neighbor_create::distributed_data :: --> ");
    // for(int i = 0; i < dn->distributed_data_idx[n_rank]/2; ++i){
    for(int i = 0; i < dn->distributed_data_idx[n_rank]; ++i){
      log_trace("[%d/%d] ", dn->distributed_data[2*i], dn->distributed_data[2*i+1]);
    }
    log_trace("\n");

    log_trace("PDM_distant_neighbor_create::distributed_data_idx :: --> ");
    for(int i = 0; i < n_rank+1; ++i){
      log_trace("%d ",  dn->distributed_data_idx[i]);
    }
    log_trace("\n");

    log_trace("PDM_distant_neighbor_create::distributed_data_n :: --> ");
    for(int i = 0; i < n_rank; ++i){
      log_trace("%d ",  dn->distributed_data_n[i]);
    }
    log_trace("\n");

    log_trace("PDM_distant_neighbor_create::requested_data_idx :: --> ");
    for(int i = 0; i < n_rank+1; ++i){
      log_trace("%d ",  dn->requested_data_idx[i]);
    }
    log_trace("\n");

    log_trace("PDM_distant_neighbor_create::requested_data_n :: --> ");
    for(int i = 0; i < n_rank; ++i){
      log_trace("%d ",  dn->requested_data_n[i]);
    }
    log_trace("\n");

    log_trace("PDM_distant_neighbor_create::order/order_unique :: --> ");
    for(int i_part = 0; i_part < dn->n_part; i_part++){
      int *_part_neighbor_idx  = dn->neighbor_idx[i_part];
      for(int i_entity = 0; i_entity < _part_neighbor_idx[n_entity[i_part]]; i_entity++){
        int u_entity = dn->order_unique[i_part][i_entity];
        int s_entity = dn->order[i_part][i_entity];
        log_trace("[%d] - order:: %d | unique:: %d | \n", i_entity, s_entity, u_entity);
      }
    }

  }

  /*
   * Free
   */
  free(requested_data);

  return (PDM_distant_neighbor_t* ) dn;
}


/**
 * \brief Exchange data between \ref PDM_distant_neighbor_t structure
 * \param [in]   id          identifier of internal structre
 *  NB : On va commencer par des entiers en stride constantes
 *
 */
void
PDM_distant_neighbor_exch
(
 PDM_distant_neighbor_t   *dn,
 size_t                    s_data,
 PDM_stride_t              t_stride,
 int                       cst_stride,
 int                     **send_entity_stride,
 void                    **send_entity_data,
 int                    ***recv_entity_stride,
 void                   ***recv_entity_data
)
{
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(dn->comm, &i_rank);
  PDM_MPI_Comm_size(dn->comm, &n_rank);

  int s_distributed_data = dn->distributed_data_idx[n_rank];
  int s_requested_data   = dn->requested_data_idx[n_rank];

  size_t *i_send_buffer = (size_t *) malloc (sizeof(size_t) * n_rank);
  size_t *i_recv_buffer = (size_t *) malloc (sizeof(size_t) * n_rank);
  int *n_send_buffer = (int *) malloc (sizeof(int) * n_rank);
  int *n_recv_buffer = (int *) malloc (sizeof(int) * n_rank);

  size_t s_send_buffer = 0;
  size_t s_recv_buffer = 0;

  for (int i = 0; i < n_rank; i++) {
    n_send_buffer[i] = 0;
    n_recv_buffer[i] = 0;
    i_send_buffer[i] = 0;
    i_recv_buffer[i] = 0;
  }

  unsigned char *send_buffer = NULL;
  unsigned char *recv_buffer = NULL;

  unsigned char **_send_entity_data = (unsigned char **) send_entity_data;
  // unsigned char **_recv_entity_data;

  int *recv_stride = NULL;
  int *recv_stride_idx = NULL;

  int s_block_unit = -1;

  if (t_stride == PDM_STRIDE_VAR) {

    int *send_stride = (int *) malloc (sizeof(int) *   s_distributed_data    );
    recv_stride      = (int *) malloc (sizeof(int) *   s_requested_data      );
    recv_stride_idx  = (int *) malloc (sizeof(int) * ( s_requested_data + 1) );

    /*
     * Prepare send stride
     */
    int idx_send = 0;
    for (int i = 0; i < s_distributed_data; i++) {
      int i_part = dn->distributed_data[2*i  ];
      int ienty = dn->distributed_data[2*i+1];
      // log_trace("send_stride[%d/%d] --> [%d,%d] -> %d \n", idx_send, s_distributed_data, i_part, ienty, send_entity_stride[i_part][ienty]);
      send_stride[idx_send++] = send_entity_stride[i_part][ienty];
    }

    PDM_MPI_Alltoallv (send_stride,
                       dn->distributed_data_n,
                       dn->distributed_data_idx,
                       PDM_MPI_INT,
                       recv_stride,
                       dn->requested_data_n,
                       dn->requested_data_idx,
                       PDM_MPI_INT,
                       dn->comm);

    /*
     * Fill the recv stride for all parts / entity
     */
    *recv_entity_stride = malloc( dn->n_part * sizeof(int *) );
    int **_recv_entity_stride = (*(int ***) recv_entity_stride);

    for(int i_part = 0; i_part < dn->n_part; i_part++){
      int *_part_neighbor_idx  = dn->neighbor_idx[i_part];
      _recv_entity_stride[i_part] = (int *) malloc( _part_neighbor_idx[dn->n_entity[i_part]] * sizeof(int *));

      for(int i_entity = 0; i_entity < _part_neighbor_idx[dn->n_entity[i_part]]; i_entity++){
        int idx = dn->distributed_part_idx[i_part] + dn->order_unique[i_part][i_entity];
        // log_trace("recv strid::[%d/%d] --> [%d,%d] -> %d \n", idx, _part_neighbor_idx[dn->n_entity[i_part]], i_part, i_entity, recv_stride[idx]);
        _recv_entity_stride[i_part][i_entity] = recv_stride[idx];
      }
    }

    /*
     * Build buffer
     */
    for (int i = 0; i < n_rank; i++) {

      /*
       * Setup send data
       */

      int iBeg = dn->distributed_data_idx[i];
      int iEnd = dn->distributed_data_idx[i] + dn->distributed_data_n[i];

      n_send_buffer[i] = 0;
      for (int k = iBeg; k < iEnd; k++)  {
        n_send_buffer[i] += send_stride[k];
      }

      n_send_buffer[i] *= (int) s_data;

      if (i > 0) {
        i_send_buffer[i] = i_send_buffer[i-1] + n_send_buffer[i-1];
      }
      else {
        i_send_buffer[i] = 0;
      }

      /*
       * Setup recv data
       */
      iBeg = dn->requested_data_idx[i];
      iEnd = dn->requested_data_idx[i] + dn->requested_data_n[i];

      n_recv_buffer[i] = 0;
      for (int k = iBeg; k < iEnd; k++) {
        n_recv_buffer[i] += recv_stride[k];
      }

      n_recv_buffer[i] *= (int) s_data;

      if (i > 0) {
        i_recv_buffer[i] = i_recv_buffer[i-1] + n_recv_buffer[i-1];
      }
      else {
        i_recv_buffer[i] = 0;
      }
    } /* End i_rank loop */

    s_send_buffer = i_send_buffer[n_rank-1] + n_send_buffer[n_rank-1];
    s_recv_buffer = i_recv_buffer[n_rank-1] + n_recv_buffer[n_rank-1];

    send_buffer = (unsigned char *) malloc(sizeof(unsigned char) * s_send_buffer);
    recv_buffer = (unsigned char *) malloc(sizeof(unsigned char) * s_recv_buffer);

    // log_trace("PDM_distant_neighbor_exch::s_send_buffer :: %d --> \n ", s_send_buffer);
    // log_trace("PDM_distant_neighbor_exch::s_recv_buffer :: %d --> \n ", s_recv_buffer);

    recv_stride_idx[0] = 0;
    for (int i = 0; i < s_requested_data; i++) {
      recv_stride_idx[i+1] = recv_stride_idx[i] + recv_stride[i];
    }

    /*
     * Compute stride for each part (in the order of send buffer )
     */
    int** stride_idx = (int**) malloc( dn->n_part * sizeof(int **));
    for(int i_part = 0; i_part < dn->n_part; i_part++){
      stride_idx[i_part] = (int*) malloc( (dn->n_entity[i_part] + 1)  * sizeof(int*));
      stride_idx[i_part][0] = 0;
      for(int i_entity = 0; i_entity < dn->n_entity[i_part]; i_entity++){
        stride_idx[i_part][i_entity+1] = stride_idx[i_part][i_entity] + send_entity_stride[i_part][i_entity];
      }
    }

    /*
     * Fill the send_buffer
     */
    // log_trace("PDM_distant_neighbor_exch::send_buffer :: --> \n ");
    int idx1 = 0;
    for (int i = 0; i < s_distributed_data; i++) {
      int i_part = dn->distributed_data[2*i  ];
      int ienty = dn->distributed_data[2*i+1];
      for(int idata = 0; idata < send_entity_stride[i_part][ienty] * (int) s_data ; idata++){
        int idxdata = stride_idx[i_part][ienty] * s_data + idata;
        // log_trace("send[%d/%d] --> [%d,%d] -> %d \n", idx1, s_distributed_data, i_part, ienty, _send_entity_data[i_part][idxdata]);
        send_buffer[idx1++] = _send_entity_data[i_part][idxdata];
      }
    }
    // log_trace("PDM_distant_neighbor_exch::send_buffer END \n ");

    for(int i_part = 0; i_part < dn->n_part; i_part++){
      free (stride_idx[i_part]);
    }
    free (stride_idx);
    free (send_stride);

  } else if (t_stride == PDM_STRIDE_CST) {

    s_block_unit = cst_stride * (int) s_data;

    // log_trace("PDM_distant_neighbor_exch::requested_data :: --> \n ");

    for (int i = 0; i < n_rank; i++) {

      i_send_buffer[i] = dn->distributed_data_idx[i] * s_block_unit;
      i_recv_buffer[i] = dn->requested_data_idx[i] * s_block_unit;

      n_send_buffer[i] = dn->distributed_data_n[i] * s_block_unit;
      n_recv_buffer[i] = dn->requested_data_n[i] * s_block_unit;

      // log_trace("[%d] send::[%d/%d] | recv::[%d/%d] \n ", i, i_send_buffer[i], n_send_buffer[i],
      //                                                        i_recv_buffer[i], n_recv_buffer[i]);

    }

    s_send_buffer = i_send_buffer[n_rank-1] + n_send_buffer[n_rank-1];
    s_recv_buffer = i_recv_buffer[n_rank-1] + n_recv_buffer[n_rank-1];

    // log_trace("PDM_distant_neighbor_exch::s_send_buffer :: %d --> \n ", s_send_buffer);
    // log_trace("PDM_distant_neighbor_exch::s_recv_buffer :: %d --> \n ", s_recv_buffer);

    send_buffer = (unsigned char *) malloc(sizeof(unsigned char) * s_send_buffer);
    recv_buffer = (unsigned char *) malloc(sizeof(unsigned char) * s_recv_buffer);

    // log_trace("PDM_distant_neighbor_exch::send_buffer :: --> \n ");
    int idx1 = 0;
    for (int i = 0; i < s_distributed_data; i++) {
      int i_part = dn->distributed_data[2*i  ];
      int ienty = dn->distributed_data[2*i+1];
      for(int idata = 0; idata < s_block_unit; idata++){
        send_buffer[idx1++] = _send_entity_data[i_part][s_block_unit*ienty+idata];
      }
    }
  }

  /*
   * Exchange
   */
  PDM_MPI_Alltoallv_l(send_buffer,
                      n_send_buffer,
                      i_send_buffer,
                      PDM_MPI_BYTE,
                      recv_buffer,
                      n_recv_buffer,
                      i_recv_buffer,
                      PDM_MPI_BYTE,
                      dn->comm);

  free(send_buffer);
  free(n_send_buffer);
  free(n_recv_buffer);
  free(i_recv_buffer);

  *recv_entity_data = malloc( dn->n_part * sizeof(unsigned char *) );
  unsigned char **_recv_entity_data = (*(unsigned char ***) recv_entity_data);


  if (t_stride == PDM_STRIDE_VAR) {

    /*
     * Compute the recv stride index for each entity
     */
    int** _recv_entity_stride = (*(int ***) recv_entity_stride);
    int** _recv_entity_stride_idx = (int**) malloc( dn->n_part * sizeof(int **));
    for(int i_part = 0; i_part < dn->n_part; i_part++){
      int *_part_neighbor_idx  = dn->neighbor_idx[i_part];
      _recv_entity_stride_idx[i_part] = (int*) malloc( (dn->n_entity[i_part] + 1)  * sizeof(int*));
      _recv_entity_stride_idx[i_part][0] = 0;
      // for(int i_entity = 0; i_entity < dn->n_entity[i_part]; i_entity++){
      for(int i_entity = 0; i_entity < _part_neighbor_idx[dn->n_entity[i_part]]; i_entity++){
        _recv_entity_stride_idx[i_part][i_entity+1] = _recv_entity_stride_idx[i_part][i_entity] + _recv_entity_stride[i_part][i_entity];
      }
    }

    if(0 == 1){
      log_trace("PDM_distant_neighbor_exch::recv_stride_idx :: --> ");
      for (int i = 0; i < s_requested_data+1; i++) {
        log_trace("%d ", recv_stride_idx[i]);
      }
      log_trace("\n");
    }

    /*
     * Copy buffer into the recv_data
     */
    for(int i_part = 0; i_part < dn->n_part; i_part++){
      int *_part_neighbor_idx  = dn->neighbor_idx[i_part];

      int recv_part_size = 0;
      for(int i_entity = 0; i_entity < _part_neighbor_idx[dn->n_entity[i_part]]; i_entity++){
        int u_enty = dn->distributed_part_idx[i_part]+dn->order_unique[i_part][i_entity];
        recv_part_size += recv_stride[u_enty];
      }

      // log_trace("PDM_distant_neighbor_exch::recv_part_size :: --> %d --> %d ( Octet ) \n ", recv_part_size, recv_part_size * s_data);
      _recv_entity_data[i_part] = (unsigned char *) malloc( recv_part_size * s_data * sizeof(unsigned char *) );

      for(int i_entity = 0; i_entity < _part_neighbor_idx[dn->n_entity[i_part]]; i_entity++){
        int u_enty   = dn->distributed_part_idx[i_part]+dn->order_unique[i_part][i_entity];
        int s_entity = dn->order[i_part][i_entity]; // On doit remettre dans l'ordre initiale !
        int idx      = recv_stride_idx[u_enty] * s_data;

        for(int idata = 0; idata < _recv_entity_stride[i_part][i_entity] * (int) s_data; idata++){
          int idxdata = _recv_entity_stride_idx[i_part][s_entity] * s_data + idata;
          // log_trace("_recv_entity_data[%d,%d] = %d \n", i_part, idxdata, recv_buffer[idx+idata]);
          _recv_entity_data[i_part][idxdata] = recv_buffer[idx+idata];
        }
      } /* End ientity */
    } /* End part */


    /*
     * Free
     */
    for(int i_part = 0; i_part < dn->n_part; i_part++){
     free(_recv_entity_stride_idx[i_part]);
    }
    free(_recv_entity_stride_idx);

  } else if (t_stride == PDM_STRIDE_CST) {

    for(int i_part = 0; i_part < dn->n_part; i_part++){
      int *_part_neighbor_idx  = dn->neighbor_idx[i_part];

      _recv_entity_data[i_part] = (unsigned char *) malloc( _part_neighbor_idx[dn->n_entity[i_part]] * s_block_unit * sizeof(unsigned char *));

      // log_trace("PDM_distant_neighbor_exch::size :: --> %d \n ", _part_neighbor_idx[dn->n_entity[i_part]] * s_block_unit);
      // log_trace("PDM_distant_neighbor_exch::recv_buffer :: --> \n ");
      for(int i_entity = 0; i_entity < _part_neighbor_idx[dn->n_entity[i_part]]; i_entity++){
        int idx      = dn->distributed_part_idx[i_part] + dn->order_unique[i_part][i_entity];
        int s_entity = dn->order[i_part][i_entity]; // On doit remettre dans l'ordre initiale !
        for(int idata = 0; idata < s_block_unit; idata++){
          _recv_entity_data[i_part][s_block_unit*s_entity+idata] = recv_buffer[s_block_unit*idx+idata];
        }
      }
    }

  }

  /*
   * Free
   */
  free(i_send_buffer);
  free(recv_buffer);
  if(recv_stride_idx != NULL){
    free(recv_stride_idx);
  }
  if(recv_stride != NULL){
    free(recv_stride);
  }

}


/**
 * \brief Exchange data between \ref PDM_distant_neighbor_t structure
 * \param [in]   id          identifier of internal structre
 *  NB : On va commencer par des entiers en stride constantes
 *
 */
void
PDM_distant_neighbor_exch_int
(
 PDM_distant_neighbor_t    *dn,
 size_t                    s_data,
 PDM_stride_t              t_stride,
 int                       cst_stride,
 int                     **send_entity_stride,
 int                     **send_entity_data,
 int                    ***recv_entity_stride,
 int                    ***recv_entity_data
)
{
  PDM_UNUSED(s_data);

  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(dn->comm, &i_rank);
  PDM_MPI_Comm_size(dn->comm, &n_rank);

  int s_distributed_data = dn->distributed_data_idx[n_rank];
  int s_requested_data   = dn->requested_data_idx[n_rank];

  size_t *i_send_buffer = (size_t *) malloc (sizeof(size_t) * n_rank);
  size_t *i_recv_buffer = (size_t *) malloc (sizeof(size_t) * n_rank);
  int *n_send_buffer = (int *) malloc (sizeof(int) * n_rank);
  int *n_recv_buffer = (int *) malloc (sizeof(int) * n_rank);

  size_t s_send_buffer = 0;
  size_t s_recv_buffer = 0;

  for (int i = 0; i < n_rank; i++) {
    n_send_buffer[i] = 0;
    n_recv_buffer[i] = 0;
    i_send_buffer[i] = 0;
    i_recv_buffer[i] = 0;
  }

  // unsigned char *send_buffer = NULL;
  // unsigned char *recv_buffer = NULL;
  int *send_buffer = NULL;
  int *recv_buffer = NULL;

  int *recv_stride = NULL;
  int *recv_stride_idx = NULL;
  if (t_stride == PDM_STRIDE_VAR) {

    int *send_stride = (int *) malloc (sizeof(int) *   s_distributed_data    );
    recv_stride      = (int *) malloc (sizeof(int) *   s_requested_data      );
    recv_stride_idx  = (int *) malloc (sizeof(int) * ( s_requested_data + 1) );

    /*
     * Prepare send stride
     */
    int idx_send = 0;
    for (int i = 0; i < s_distributed_data; i++) {
      int i_part = dn->distributed_data[2*i  ];
      int ienty = dn->distributed_data[2*i+1];
      // log_trace("send_stride[%d/%d] --> [%d,%d] -> %d \n", idx_send, s_distributed_data, i_part, ienty, send_entity_stride[i_part][ienty]);
      send_stride[idx_send++] = send_entity_stride[i_part][ienty];
    }

    PDM_MPI_Alltoallv (send_stride,
                       dn->distributed_data_n,
                       dn->distributed_data_idx,
                       PDM_MPI_INT,
                       recv_stride,
                       dn->requested_data_n,
                       dn->requested_data_idx,
                       PDM_MPI_INT,
                       dn->comm);

    /*
     * Fill the recv stride for all parts / entity
     */
    *recv_entity_stride = malloc( dn->n_part * sizeof(int *) );
    int **_recv_entity_stride = (*(int ***) recv_entity_stride);

    for(int i_part = 0; i_part < dn->n_part; i_part++){
      int *_part_neighbor_idx  = dn->neighbor_idx[i_part];
      _recv_entity_stride[i_part] = (int *) malloc( _part_neighbor_idx[dn->n_entity[i_part]] * sizeof(int *));

      for(int i_entity = 0; i_entity < _part_neighbor_idx[dn->n_entity[i_part]]; i_entity++){
        int idx = dn->distributed_part_idx[i_part] + dn->order_unique[i_part][i_entity];
        // log_trace("recv strid::[%d/%d] --> [%d,%d] -> %d \n", idx, _part_neighbor_idx[dn->n_entity[i_part]], i_part, i_entity, recv_stride[idx]);
        _recv_entity_stride[i_part][i_entity] = recv_stride[idx];
      }
    }

    /*
     * Build buffer
     */
    for (int i = 0; i < n_rank; i++) {

      /*
       * Setup send data
       */

      int iBeg = dn->distributed_data_idx[i];
      int iEnd = dn->distributed_data_idx[i] + dn->distributed_data_n[i];

      n_send_buffer[i] = 0;
      for (int k = iBeg; k < iEnd; k++)  {
        n_send_buffer[i] += send_stride[k];
      }

      // n_send_buffer[i] *= (int) s_data;

      if (i > 0) {
        i_send_buffer[i] = i_send_buffer[i-1] + n_send_buffer[i-1];
      }
      else {
        i_send_buffer[i] = 0;
      }

      /*
       * Setup recv data
       */
      iBeg = dn->requested_data_idx[i];
      iEnd = dn->requested_data_idx[i] + dn->requested_data_n[i];

      n_recv_buffer[i] = 0;
      for (int k = iBeg; k < iEnd; k++) {
        n_recv_buffer[i] += recv_stride[k];
      }

      // n_recv_buffer[i] *= (int) s_data;

      if (i > 0) {
        i_recv_buffer[i] = i_recv_buffer[i-1] + n_recv_buffer[i-1];
      }
      else {
        i_recv_buffer[i] = 0;
      }
    } /* End i_rank loop */

    s_send_buffer = i_send_buffer[n_rank-1] + n_send_buffer[n_rank-1];
    s_recv_buffer = i_recv_buffer[n_rank-1] + n_recv_buffer[n_rank-1];

    // send_buffer = (unsigned char *) malloc(sizeof(unsigned char) * s_send_buffer);
    // recv_buffer = (unsigned char *) malloc(sizeof(unsigned char) * s_recv_buffer);
    // log_trace("PDM_distant_neighbor_exch::s_send_buffer :: %d --> \n ", s_send_buffer);
    // log_trace("PDM_distant_neighbor_exch::s_recv_buffer :: %d --> \n ", s_recv_buffer);

    send_buffer = (int *) malloc(sizeof(int) * s_send_buffer);
    recv_buffer = (int *) malloc(sizeof(int) * s_recv_buffer);

    recv_stride_idx[0] = 0;
    for (int i = 0; i < s_requested_data; i++) {
      recv_stride_idx[i+1] = recv_stride_idx[i] + recv_stride[i];
    }

    /*
     * Compute stride for each part (in the order of send buffer )
     */
    int** stride_idx = (int**) malloc( dn->n_part * sizeof(int **));
    for(int i_part = 0; i_part < dn->n_part; i_part++){
      stride_idx[i_part] = (int*) malloc( (dn->n_entity[i_part] + 1)  * sizeof(int*));
      stride_idx[i_part][0] = 0;
      for(int i_entity = 0; i_entity < dn->n_entity[i_part]; i_entity++){
        stride_idx[i_part][i_entity+1] = stride_idx[i_part][i_entity] + send_entity_stride[i_part][i_entity];
      }
    }

    /*
     * Fill the send_buffer
     */
    // log_trace("PDM_distant_neighbor_exch::send_buffer :: --> \n ");
    int idx1 = 0;
    for (int i = 0; i < s_distributed_data; i++) {
      int i_part = dn->distributed_data[2*i  ];
      int ienty = dn->distributed_data[2*i+1];
      for(int idata = 0; idata < send_entity_stride[i_part][ienty]; idata++){
        int idxdata = stride_idx[i_part][ienty] + idata;
        // log_trace("send[%d/%d] --> [%d,%d] -> %d \n", idx1, s_distributed_data, i_part, ienty, send_entity_data[i_part][idxdata]);
        send_buffer[idx1++] = send_entity_data[i_part][idxdata];
      }
    }
    // log_trace("PDM_distant_neighbor_exch::send_buffer END \n ");

    for(int i_part = 0; i_part < dn->n_part; i_part++){
      free (stride_idx[i_part]);
    }
    free (stride_idx);
    free (send_stride);

  } else if (t_stride == PDM_STRIDE_CST) {

    // int s_block_unit = cst_stride * (int) s_data;

    // log_trace("PDM_distant_neighbor_exch::requested_data :: --> \n ");

    for (int i = 0; i < n_rank; i++) {

      // i_send_buffer[i] = dn->distributed_data_idx[i] * cst_stride * (int) s_data;
      // i_recv_buffer[i] = dn->requested_data_idx[i] * cst_stride * (int) s_data;

      // n_send_buffer[i] = dn->distributed_data_n[i] * cst_stride * (int) s_data;
      // n_recv_buffer[i] = dn->requested_data_n[i] * cst_stride * (int) s_data;

      i_send_buffer[i] = dn->distributed_data_idx[i] * cst_stride;
      i_recv_buffer[i] = dn->requested_data_idx[i] * cst_stride;

      n_send_buffer[i] = dn->distributed_data_n[i] * cst_stride;
      n_recv_buffer[i] = dn->requested_data_n[i] * cst_stride;

      // log_trace("[%d] send::[%d/%d] | recv::[%d/%d] \n ", i, i_send_buffer[i], n_send_buffer[i],
      //                                                        i_recv_buffer[i], n_recv_buffer[i]);

    }

    s_send_buffer = i_send_buffer[n_rank-1] + n_send_buffer[n_rank-1];
    s_recv_buffer = i_recv_buffer[n_rank-1] + n_recv_buffer[n_rank-1];

    // send_buffer = (unsigned char *) malloc(sizeof(unsigned char) * s_send_buffer);
    // recv_buffer = (unsigned char *) malloc(sizeof(unsigned char) * s_recv_buffer);
    // log_trace("PDM_distant_neighbor_exch::s_send_buffer :: %d --> \n ", s_send_buffer);
    // log_trace("PDM_distant_neighbor_exch::s_recv_buffer :: %d --> \n ", s_recv_buffer);
    send_buffer = (int *) malloc(sizeof(int) * s_send_buffer);
    recv_buffer = (int *) malloc(sizeof(int) * s_recv_buffer);

    // log_trace("PDM_distant_neighbor_exch::send_buffer :: --> \n ");
    int idx1 = 0;
    for (int i = 0; i < s_distributed_data; i++) {
      int i_part = dn->distributed_data[2*i  ];
      int ienty = dn->distributed_data[2*i+1];
      for(int idata = 0; idata < cst_stride; idata++){
        // log_trace("send[%d/%d] --> [%d,%d] -> %d \n", idx1, s_distributed_data, i_part, ienty, send_entity_data[i_part][cst_stride*ienty+idata]);
        send_buffer[idx1++] = send_entity_data[i_part][cst_stride*ienty+idata];
      }
    }
    // log_trace("PDM_distant_neighbor_exch::send_buffer END \n ");
  }


  /*
   * Exchange
   */
  // PDM_MPI_Alltoallv_l(send_buffer,
  //                     n_send_buffer,
  //                     i_send_buffer,
  //                     PDM_MPI_BYTE,
  //                     recv_buffer,
  //                     n_recv_buffer,
  //                     i_recv_buffer,
  //                     PDM_MPI_BYTE,
  //                     dn->comm);
  PDM_MPI_Alltoallv_l(send_buffer,
                      n_send_buffer,
                      i_send_buffer,
                      PDM_MPI_INT,
                      recv_buffer,
                      n_recv_buffer,
                      i_recv_buffer,
                      PDM_MPI_INT,
                      dn->comm);


  free(send_buffer);
  free(n_send_buffer);
  free(n_recv_buffer);
  free(i_recv_buffer);

  /*
   * Une seule valeur est echangé mais plusieurs occurence peuvent exister donc on passe du buffer MPI
   * au donné sans le sort/unique
   */
  *recv_entity_data = malloc( dn->n_part * sizeof(int *) );
  int **_recv_entity_data = (*(int ***) recv_entity_data);

  if (t_stride == PDM_STRIDE_VAR) {

    /*
     * Compute the recv stride index for each entity
     */
    int** _recv_entity_stride = (*(int ***) recv_entity_stride);
    int** _recv_entity_stride_idx = (int**) malloc( dn->n_part * sizeof(int **));
    for(int i_part = 0; i_part < dn->n_part; i_part++){
      int *_part_neighbor_idx  = dn->neighbor_idx[i_part];
      _recv_entity_stride_idx[i_part] = (int*) malloc( (dn->n_entity[i_part] + 1)  * sizeof(int*));
      _recv_entity_stride_idx[i_part][0] = 0;
      // for(int i_entity = 0; i_entity < dn->n_entity[i_part]; i_entity++){
      for(int i_entity = 0; i_entity < _part_neighbor_idx[dn->n_entity[i_part]]; i_entity++){
        _recv_entity_stride_idx[i_part][i_entity+1] = _recv_entity_stride_idx[i_part][i_entity] + _recv_entity_stride[i_part][i_entity];
      }
    }

    if(0 == 1){
      log_trace("PDM_distant_neighbor_exch::recv_stride_idx :: --> ");
      for (int i = 0; i < s_requested_data+1; i++) {
        log_trace("%d ", recv_stride_idx[i]);
      }
      log_trace("\n");
    }

    /*
     * Copy buffer into the recv_data
     */
    for(int i_part = 0; i_part < dn->n_part; i_part++){
      int *_part_neighbor_idx  = dn->neighbor_idx[i_part];

      int recv_part_size = 0;
      for(int i_entity = 0; i_entity < _part_neighbor_idx[dn->n_entity[i_part]]; i_entity++){
        int u_enty = dn->distributed_part_idx[i_part]+dn->order_unique[i_part][i_entity];
        recv_part_size += recv_stride[u_enty];
      }

      // log_trace("PDM_distant_neighbor_exch::recv_part_size :: --> %d \n ", recv_part_size);
      _recv_entity_data[i_part] = (int *) malloc( recv_part_size * sizeof(int *));

      for(int i_entity = 0; i_entity < _part_neighbor_idx[dn->n_entity[i_part]]; i_entity++){

        int u_enty   = dn->distributed_part_idx[i_part]+dn->order_unique[i_part][i_entity];
        int s_entity = dn->order[i_part][i_entity]; // On doit remettre dans l'ordre initiale !
        int idx      = recv_stride_idx[u_enty];

        // log_trace(" debug:: [%d] - [u_enty::%d] | [unique::%d] | [recv_stride_idx::%d] | [shiftpart::%d]\n",
        //           i_entity, u_enty, dn->order_unique[i_part][i_entity], recv_stride_idx[u_enty], dn->distributed_part_idx[i_part]);

        for(int idata = 0; idata < _recv_entity_stride[i_part][i_entity]; idata++){
          int idxdata = _recv_entity_stride_idx[i_part][s_entity] + idata;
          // log_trace("recv ::[%d/%d] --> [%d,%d] -> %d \n", idx+idata, _part_neighbor_idx[dn->n_entity[i_part]], i_part, i_entity, recv_buffer[idx+idata]);
          // log_trace("_recv_entity_data[%d,%d] = %d \n", i_part, idxdata, recv_buffer[idx+idata]);
          _recv_entity_data[i_part][idxdata] = recv_buffer[idx+idata];
        }
      } /* End ientity */
    } /* End part */


    /*
     * Free
     */
    for(int i_part = 0; i_part < dn->n_part; i_part++){
     free(_recv_entity_stride_idx[i_part]);
    }
    free(_recv_entity_stride_idx);


  } else if (t_stride == PDM_STRIDE_CST) {

    // Shift is not good because the buffer contains only one occurence of each elements !!!
    for(int i_part = 0; i_part < dn->n_part; i_part++){
      int *_part_neighbor_idx  = dn->neighbor_idx[i_part];
      _recv_entity_data[i_part] = (int *) malloc( _part_neighbor_idx[dn->n_entity[i_part]] * cst_stride * sizeof(int *));

      // log_trace("PDM_distant_neighbor_exch::size :: --> %d \n ", _part_neighbor_idx[dn->n_entity[i_part]] * cst_stride);
      // log_trace("PDM_distant_neighbor_exch::recv_buffer :: --> \n ");
      for(int i_entity = 0; i_entity < _part_neighbor_idx[dn->n_entity[i_part]]; i_entity++){
        int idx      = dn->distributed_part_idx[i_part] + dn->order_unique[i_part][i_entity];
        int s_entity = dn->order[i_part][i_entity]; // On doit remettre dans l'ordre initiale !
        for(int idata = 0; idata < cst_stride; idata++){
          // log_trace("recv ::[%d/%d] --> [%d,%d] -> %d \n", idx, _part_neighbor_idx[dn->n_entity[i_part]], i_part, i_entity, recv_buffer[idx]);
          // log_trace("recv ::[%d/%d] --> [%d,%d] -> %d \n", cst_stride*idx+idata, _part_neighbor_idx[dn->n_entity[i_part]], i_part, i_entity, recv_buffer[cst_stride*idx+idata]);
          // log_trace("recv ::[%d/%d] --> [%d,%d] -> %d \n", cst_stride*i_entity+idata, _part_neighbor_idx[dn->n_entity[i_part]], i_part, i_entity, recv_buffer[cst_stride*idx+idata]);
          _recv_entity_data[i_part][cst_stride*s_entity+idata] = recv_buffer[cst_stride*idx+idata];
        }
      }
    }
  }


  /*
   * Free
   */
  free(i_send_buffer);
  free(recv_buffer);
  if(recv_stride_idx != NULL){
    free(recv_stride_idx);
  }
  if(recv_stride != NULL){
    free(recv_stride);
  }

}


/**
 *
 * \brief Free an distant negihtbor structure
 *
 * \param [in]   id                 Identifier
 *
 */

void
PDM_distant_neighbor_free
(
 PDM_distant_neighbor_t* dn
)
{

  for(int i_part = 0; i_part < dn->n_part; i_part++){
    free(dn->order[i_part]);
    free(dn->order_unique[i_part]);
  }
  free(dn->order);
  free(dn->order_unique);

  free(dn->requested_data_n);
  free(dn->requested_data_idx);
  free(dn->distributed_part_idx);
  free(dn->distributed_part_n);

  free(dn->distributed_data);
  free(dn->distributed_data_n);
  free(dn->distributed_data_idx);

  free (dn);

}

#ifdef __cplusplus
}
#endif /* __cplusplus */
