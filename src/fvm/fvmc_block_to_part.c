/*============================================================================
 * Convert between block distribution and general domain partition.
 *============================================================================*/

/*
  This file is part of the "Finite Volume Mesh" library, intended to provide
  finite volume mesh and associated fields I/O and manipulation services.

  Copyright (C) 2011       ONERA

  Copyright (C) 2008  EDF

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

#include <bftc_mem.h>
#include <bftc_error.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "fvmc_defs.h"
#include "fvmc_config_defs.h"

#include "fvmc_block_dist.h"
#include "fvmc_order.h"
#include "fvmc_parall.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "fvmc_block_to_part.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*=============================================================================
 * Local type definitions
 *============================================================================*/

/* Structure used to redistribute data */

#if defined(FVMC_HAVE_MPI)

struct _fvmc_block_to_part_t {

  MPI_Comm     comm;         /* Associated MPI communicator */

  int          n_ranks;      /* Number of ranks associated with
                                communicator */

  size_t       n_block_ents; /* Number of entities to send (block) */
  size_t       n_part_ents;  /* Number of entities to receive (partition) */
  size_t       send_size;    /* Size of send buffer for MPI_Alltoall
                                (recv_size not necessary, as recv_size
                                should always be equal to n_recv_ents,
                                though elements may arrive in a different
                                order) */

  int         *send_count;   /* Send counts for MPI_Alltoall */
  int         *recv_count;   /* Receive counts for MPI_Alltoall */
  int         *send_displ;   /* Send displs for MPi_Alltoall */
  int         *recv_displ;   /* Receive displs for MPi_Alltoall */

  fvmc_lnum_t  *send_list;    /* List of entities to send in rank order */
  fvmc_lnum_t  *recv_order;   /* Ordering of received entities by
                                increasing global number (duplicates removed) */

  const fvmc_gnum_t  *recv_global_num;  /* Possibly shared global numbers */
  fvmc_gnum_t        *_recv_global_num; /* Private global entity numbers
                                          (NULL if shared); */

};

#endif /* defined(FVMC_HAVE_MPI) */

/*============================================================================
 * Local function defintions
 *============================================================================*/

#if defined(FVMC_HAVE_MPI)

/*----------------------------------------------------------------------------
 * Compute rank displacement based on count.
 *
 * arguments:
 *   n_ranks <-- number of ranks
 *   count   <-- number of entities per rank (size: n_ranks)
 *   displ   --> entity displacement in cumulative array (size: n_ranks)
 *
 * returns:
 *   cumulative count for all ranks
 *----------------------------------------------------------------------------*/

static fvmc_lnum_t
_compute_displ(int        n_ranks,
               const int  count[],
               int        displ[])
{
  int i;
  fvmc_lnum_t total_count = 0;

  displ[0] = 0;

  for (i = 1; i < n_ranks; i++)
    displ[i] = displ[i-1] + count[i-1];

  total_count = displ[n_ranks-1] + count[n_ranks-1];

  return total_count;
}

/*----------------------------------------------------------------------------
 * Build ordered list of entities based on their global number, ignoring
 * duplicates.
 *
 * The caller is responsible for freing the ordered_ent_list allocated by
 * this function when no longer needed.
 *
 * arguments:
 *   n_ents         <-- number of entities
 *   ent_global_num <-- global entity numbers
 *   n_ordered_ents --> number of distinct ordered entities
 *   ordered_ent    --> list of distinct ordered entities
 *----------------------------------------------------------------------------*/

static void
_ordered_list(size_t              n_ents,
              const fvmc_gnum_t    ent_global_num[],
              size_t             *n_ordered_ents,
              fvmc_lnum_t        **ordered_ent)
{
  size_t i, j;
  size_t _n_ordered_ents = 1;

  fvmc_lnum_t *order = NULL;
  fvmc_lnum_t *_ordered_ent = NULL;

  if (n_ents == 0)
    return;

  /* Sort global numbers */

  order = fvmc_order_local(NULL, ent_global_num, n_ents);

  /* Count number of distinct global entities */

  for (i = 1; i < n_ents; i++) {
    if (ent_global_num[order[i]] > ent_global_num[order[i-1]])
      _n_ordered_ents += 1;
  }

  /* Build retrieval list, counting multiply-encountered entities once */

  if (_n_ordered_ents == n_ents) {

    _ordered_ent = order;
    order = NULL;

  }
  else {

    BFTC_MALLOC(_ordered_ent, _n_ordered_ents, fvmc_lnum_t);

    _ordered_ent[0] = order[0];
    for (i = 1, j = 1; i < n_ents; i++) {
      if (ent_global_num[order[i]] > ent_global_num[order[i-1]])
        _ordered_ent[j++] = order[i];
    }

    BFTC_FREE(order);
  }

  *n_ordered_ents = _n_ordered_ents;
  *ordered_ent = _ordered_ent;
}

/*----------------------------------------------------------------------------
 * Initialize distributor with block data
 *
 * arguments:
 *   d  <-> distribution helper
 *   bi <-- block to partition range and size info
 *----------------------------------------------------------------------------*/

static void
_init_global_num(fvmc_block_to_part_t       *d,
                 fvmc_block_to_part_info_t   bi)
{
  size_t j;

  size_t recv_size = 0;

  fvmc_gnum_t *send_global_num = NULL;
  fvmc_gnum_t *recv_global_num = NULL;

  /* Build temporay global numbers */

  BFTC_MALLOC(send_global_num, d->send_size, fvmc_gnum_t);

  for (j = 0; j < d->send_size; j++)
    send_global_num[j] = d->send_list[j] + bi.gnum_range[0];

  /* Exchange global numbers */

  BFTC_MALLOC(recv_global_num, d->n_part_ents, fvmc_gnum_t);

  MPI_Alltoallv(send_global_num, d->send_count, d->send_displ, FVMC_MPI_GNUM,
                recv_global_num, d->recv_count, d->recv_displ, FVMC_MPI_GNUM,
                d->comm);

  /* Count number of distinct global entities and build retrieval index,
     counting multiply-encountered entities once */

  _ordered_list(d->n_part_ents,
                recv_global_num,
                &recv_size,
                &(d->recv_order));

  if (d->n_part_ents != recv_size)
    bftc_error
      (__FILE__, __LINE__, 0,
       _("inconsistent sizes computed for a block to partition distributor\n"
         "(%lu expected, %lu determined)."),
       (unsigned long)(d->n_part_ents), (unsigned long)recv_size);

  /* Now build global number list */

  BFTC_MALLOC(d->_recv_global_num, d->n_part_ents, fvmc_gnum_t);
  d->recv_global_num = d->_recv_global_num;

  for (j = 0; j < d->n_part_ents; j++)
    d->_recv_global_num[j] = recv_global_num[d->recv_order[j]];

  BFTC_FREE(recv_global_num);
  BFTC_FREE(send_global_num);
}

/*----------------------------------------------------------------------------
 * Create distribution helper structure.
 *
 * Send and receive counts and displacements are allocated, but not
 * fully initialized at this point: only the send count is set to zero.
 *
 * arguments:
 *   comm <-- communicator
 *
 * returns:
 *   empty communicator structure
 *----------------------------------------------------------------------------*/

static fvmc_block_to_part_t *
_block_to_part_create(MPI_Comm comm)
{
  int i;

  fvmc_block_to_part_t *d;

  BFTC_MALLOC(d, 1, fvmc_block_to_part_t);

  d->comm = comm;

  MPI_Comm_size(comm, &(d->n_ranks));

  d->n_block_ents = 0;
  d->n_part_ents = 0;
  d->send_size = 0;

  BFTC_MALLOC(d->send_count, d->n_ranks, int);
  BFTC_MALLOC(d->recv_count, d->n_ranks, int);
  BFTC_MALLOC(d->send_displ, d->n_ranks, int);
  BFTC_MALLOC(d->recv_displ, d->n_ranks, int);

  for (i = 0; i < d->n_ranks; i++)
    d->send_count[i] = 0;

  d->send_list = NULL;
  d->recv_order = NULL;

  d->recv_global_num = NULL;
  d->_recv_global_num = NULL;

  return d;
}

#endif /* defined(FVMC_HAVE_MPI) */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Compute block size and rank info for use with a block distribution.
 *
 * arguments:
 *   rank_id        <-- id of local rank (ignored in serial mode)
 *   n_ranks        <-- number of associated ranks
 *   min_rank_step  <-- minimum rank step between blocks
 *   min_block_size <-- minimum number of entities per block
 *   n_g_ents       <-- total number of associated entities
 *
 * returns:
 *   block size and range info structure
 *----------------------------------------------------------------------------*/

fvmc_block_to_part_info_t
fvmc_block_to_part_compute_sizes(int         rank_id,
                                int         n_ranks,
                                int         min_rank_step,
                                fvmc_lnum_t  min_block_size,
                                fvmc_gnum_t  n_g_ents)
{
  fvmc_block_to_part_info_t bi;

  fvmc_block_dist_info_t _bi= fvmc_block_dist_compute_sizes(rank_id,
                                                          n_ranks,
                                                          min_rank_step,
                                                          min_block_size,
                                                          n_g_ents);

  bi.gnum_range[0] = _bi.gnum_range[0];
  bi.gnum_range[1] = _bi.gnum_range[1];
  bi.n_ranks = _bi.n_ranks;
  bi.rank_step = _bi.rank_step;
  bi.block_size = _bi.block_size;

  return bi;
}

#if defined(FVMC_HAVE_MPI)

/*----------------------------------------------------------------------------
 * Create block to partition distributor using entity destination rank array.
 *
 * arguments:
 *   comm     <-- communicator
 *   block    <-- block to partition range and size info
 *   ent_rank <-- destination rank for each entity
 *                (size: block.gnum_range[1] - block.gnum_range[0])
 *
 * returns:
 *   initialized block to partition distributor
 *----------------------------------------------------------------------------*/

fvmc_block_to_part_t *
fvmc_block_to_part_create_by_rank(MPI_Comm                  comm,
                                 fvmc_block_to_part_info_t  block,
                                 int                       ent_rank[])
{
  int i;
  size_t j;

  fvmc_block_to_part_t *d = _block_to_part_create(comm);

  const int n_ranks = d->n_ranks;

  d->n_block_ents = block.gnum_range[1] - block.gnum_range[0];

  /* Count values to send and receive */

  for (j = 0; j < d->n_block_ents; j++)
    d->send_count[ent_rank[j]] += 1;

  MPI_Alltoall(d->send_count, 1, MPI_INT, d->recv_count, 1, MPI_INT, comm);

  d->send_size = _compute_displ(n_ranks, d->send_count, d->send_displ);
  d->n_part_ents = _compute_displ(n_ranks, d->recv_count, d->recv_displ);

  /* Prepare send list (using send_displ for insertion positions) */

  BFTC_MALLOC(d->send_list, d->send_size, fvmc_lnum_t);

  for (j = 0; j < d->send_size; j++) {
    int send_rank = ent_rank[j];
    d->send_list[d->send_displ[send_rank]] = j;
    d->send_displ[send_rank] += 1;
  }

  /* Reset send_displ */

  for (i = 0; i < n_ranks; i++)
    d->send_displ[i] -= d->send_count[i];

  /* Build global numbering and retrieval index */

  _init_global_num(d, block);

  /* Return initialized structure */

  return d;
}

/*----------------------------------------------------------------------------
 * Initialize block to partition distributor with block data using
 * strided adjacency array.
 *
 * The adjacency array uses 1-n based global numbers. 0 values are
 * allowed and may be used to represent empty adjacencies.
 *
 * For example, in a face -> element adjacency relation, each face
 * is adjacent to 2 elements (thus a stride of 2), except for
 * boundary faces which are adjacent to only 1 element; in this case,
 * the adjacent element number for the exterior side of the face is 0.
 *
 * arguments:
 *   comm              <-- communicator
 *   block             <-- block size and range info
 *   adjacent_block    <-- block info for adjacent entities
 *   stride            <-- stride of adjacency array
 *   adjacency         <-- entity adjacency (1 to n numbering)
 *   adjacent_ent_rank <-- destination rank for adjacent entities, or
 *                         NULL if based on block size and range only.
 *
 * returns:
 *   initialized block to partition distributor
 *----------------------------------------------------------------------------*/

fvmc_block_to_part_t *
fvmc_block_to_part_create_by_adj_s(MPI_Comm                  comm,
                                  fvmc_block_to_part_info_t  block,
                                  fvmc_block_to_part_info_t  adjacent_block,
                                  int                       stride,
                                  fvmc_gnum_t                adjacency[],
                                  int                       adjacent_ent_rank[])
{
  int i, k;
  fvmc_lnum_t j;

  fvmc_lnum_t   adj_send_size = 0, adj_recv_size = 0;
  int         *adj_send_count = NULL, *adj_recv_count = NULL;
  int         *adj_send_displ = NULL, *adj_recv_displ = NULL;
  fvmc_lnum_t  *rank_flag = NULL;
  fvmc_gnum_t  *adj_send_num = NULL, *adj_recv_num = NULL;

  fvmc_block_to_part_t *d = _block_to_part_create(comm);

  int rank = -1;
  const int n_ranks = d->n_ranks;

  const fvmc_lnum_t n_ents = block.gnum_range[1] - block.gnum_range[0];

  MPI_Comm_rank(comm, &rank);

  /* Determine an entity's adjacent entities, and use their
     global number to determine by which ranks they were read;
     using this knowledge, we can query these ranks to know to
     which ranks these adjacent entities were sent
     --------------------------------------------------------- */

  /* Count values to send and receive to each processor */

  BFTC_MALLOC(adj_send_count, n_ranks, int);
  BFTC_MALLOC(adj_recv_count, n_ranks, int);

  BFTC_MALLOC(adj_send_displ, n_ranks, int);
  BFTC_MALLOC(adj_recv_displ, n_ranks, int);

  for (i = 0; i < n_ranks; i++)
    adj_send_count[i] = 0;

  for (j = 0; j < n_ents; j++) {
    for (k = 0; k < stride; k++) {
      fvmc_gnum_t adj_g_num = adjacency[j*stride + k];
      if (adj_g_num > 0) {
        int adj_ent_rank =   ((adj_g_num-1) / adjacent_block.block_size)
                           * adjacent_block.rank_step;
        adj_send_count[adj_ent_rank] += 1;
      }
    }
  }

  MPI_Alltoall(adj_send_count, 1, MPI_INT, adj_recv_count, 1, MPI_INT, comm);

  adj_send_size = _compute_displ(n_ranks, adj_send_count, adj_send_displ);
  adj_recv_size = _compute_displ(n_ranks, adj_recv_count, adj_recv_displ);

  /* Prepare destination rank request (temporarily modifying adj_send_displ) */

  BFTC_MALLOC(adj_send_num, adj_send_size, fvmc_gnum_t);
  BFTC_MALLOC(adj_recv_num, adj_recv_size, fvmc_gnum_t);

  for (j = 0; j < n_ents; j++) {
    for (k = 0; k < stride; k++) {
      fvmc_gnum_t adj_g_num = adjacency[j*stride + k];
      if (adj_g_num > 0) {
        int adj_ent_rank =   ((adj_g_num-1) / adjacent_block.block_size)
                           * adjacent_block.rank_step;
        adj_send_num[adj_send_displ[adj_ent_rank]] = adj_g_num;
        adj_send_displ[adj_ent_rank] += 1;
      }
    }
  }

  for (i = 0; i < n_ranks; i++)
    adj_send_displ[i] -= adj_send_count[i];

  /* Adj_Send destination rank request */

  MPI_Alltoallv(adj_send_num, adj_send_count, adj_send_displ, FVMC_MPI_GNUM,
                adj_recv_num, adj_recv_count, adj_recv_displ, FVMC_MPI_GNUM,
                d->comm);

  /* Reply by indicating to which ranks indicated adjacent entities
     are assigned (reusing adj_recv_num) */

  if (adjacent_ent_rank != NULL) {
    for (j = 0; j < adj_recv_size; j++) {
      fvmc_lnum_t adj_l_id = (adj_recv_num[j] - 1) % adjacent_block.block_size;
      adj_recv_num[j] = adjacent_ent_rank[adj_l_id];
    }
  }

  else {
    for (j = 0; j < adj_recv_size; j++)
      adj_recv_num[j] = rank;
  }

  /* Send and receive arguments are inverted as this is a "reply" */

  MPI_Alltoallv(adj_recv_num, adj_recv_count, adj_recv_displ, FVMC_MPI_GNUM,
                adj_send_num, adj_send_count, adj_send_displ, FVMC_MPI_GNUM,
                d->comm);

  BFTC_FREE(adj_recv_num);
  BFTC_FREE(adj_recv_count);
  BFTC_FREE(adj_recv_displ);

  /* We now need to extract the ranks to which each entity will be sent,
     based on where its adjacent entities were sent (and thus where
     it will be needed).
     -------------------------------------------------------------------- */

  BFTC_MALLOC(rank_flag, n_ranks, int);

  for (i = 0; i < n_ranks; i++)
    rank_flag[i] = -1;

  for (j = 0; j < n_ents; j++) {
    for (k = 0; k < stride; k++) {
      fvmc_gnum_t adj_g_num = adjacency[j*stride + k];
      if (adj_g_num > 0) {
        int adj_ent_rank =   ((adj_g_num-1) / adjacent_block.block_size)
                           * adjacent_block.rank_step;
        int send_rank = adj_send_num[adj_send_displ[adj_ent_rank]];
        if (rank_flag[send_rank] < j) {
          d->send_count[send_rank] += 1;
          rank_flag[send_rank] = j;
        }
        adj_send_displ[adj_ent_rank] += 1;
      }
    }
  }

  for (i = 0; i < n_ranks; i++)
    adj_send_displ[i] -= adj_send_count[i];

  /* Exchange send count */

  MPI_Alltoall(d->send_count, 1, MPI_INT, d->recv_count, 1, MPI_INT, comm);

  d->send_size = _compute_displ(n_ranks, d->send_count, d->send_displ);
  d->n_part_ents = _compute_displ(n_ranks, d->recv_count, d->recv_displ);

  /* Prepare send list (using send_displ for insertion positions) */

  BFTC_MALLOC(d->send_list, d->send_size, fvmc_lnum_t);

  for (i = 0; i < n_ranks; i++)
    rank_flag[i] = -1;

  for (j = 0; j < n_ents; j++) {
    for (k = 0; k < stride; k++) {
      fvmc_gnum_t adj_g_num = adjacency[j*stride + k];
      if (adj_g_num > 0) {
        int adj_ent_rank =   ((adj_g_num-1) / adjacent_block.block_size)
                           * adjacent_block.rank_step;
        int send_rank = adj_send_num[adj_send_displ[adj_ent_rank]];
        if (rank_flag[send_rank] < j) {
          d->send_list[d->send_displ[send_rank]] = j;
          d->send_displ[send_rank] += 1;
          rank_flag[send_rank] = j;
        }
        adj_send_displ[adj_ent_rank] += 1;
      }
    }
  }

  for (i = 0; i < n_ranks; i++)
    d->send_displ[i] -= d->send_count[i];

  BFTC_FREE(rank_flag);

  BFTC_FREE(adj_send_num);
  BFTC_FREE(adj_send_count);
  BFTC_FREE(adj_send_displ);

  /* Build global numbering and retrieval index */

  _init_global_num(d, block);

  /* Return initialized structure */

  return d;
}

/*----------------------------------------------------------------------------
 * Initialize block to partition distributor for entities adjacent to
 * already distributed entities.
 *
 * arguments:
 *   comm           <-- communicator
 *   block          <-- block size and range info
 *   adjacent_block <-- block info for adjacent entities
 *   adjacency      <-- entity adjacency (1 to n numbering)
 *
 * returns:
 *   initialized block to partition distributor
 *----------------------------------------------------------------------------*/

fvmc_block_to_part_t *
fvmc_block_to_part_create_adj(MPI_Comm                  comm,
                             fvmc_block_to_part_info_t  adjacent_block,
                             size_t                    adjacency_size,
                             const fvmc_gnum_t          adjacency[])
{
  int i;
  size_t j;

  size_t recv_size = 0;
  fvmc_lnum_t *adj_list = NULL, *_adj_list = NULL;
  fvmc_gnum_t *send_num = NULL, *recv_num = NULL;

  fvmc_block_to_part_t *d = _block_to_part_create(comm);

  const int n_ranks = d->n_ranks;

  /* Sort adjacency list so as to remove duplicates */

  _ordered_list(adjacency_size,
                adjacency,
                &d->n_part_ents,
                &_adj_list);

  /* Use adjacent global number to determine to which ranks they will
     be sent. Send and receive counts are exchanged here, as ranks
     requesting (i.e receiving entities) are determined first */

  for (i = 0; i < d->n_ranks; i++)
    d->recv_count[i] = 0;

 /* Ignore possible 0 values in global numbering */

  if (d->n_part_ents > 0) {
    if (adjacency[_adj_list[0]] == 0) {
      d->n_part_ents -= 1;
      adj_list = _adj_list + 1;
    }
    else
      adj_list = _adj_list;
  }

  for (j = 0; j < d->n_part_ents; j++) {
    fvmc_gnum_t adj_g_id = adjacency[adj_list[j]] - 1;
    int adj_ent_rank =   (adj_g_id / adjacent_block.block_size)
                       * adjacent_block.rank_step;
    d->recv_count[adj_ent_rank] += 1;
  }

  MPI_Alltoall(d->recv_count, 1, MPI_INT, d->send_count, 1, MPI_INT, comm);

  d->send_size = _compute_displ(n_ranks, d->send_count, d->send_displ);
  recv_size = _compute_displ(n_ranks, d->recv_count, d->recv_displ);

  if (d->n_part_ents != recv_size)
    bftc_error
      (__FILE__, __LINE__, 0,
       _("inconsistent sizes computed for a block to partition distributor\n"
         "(%lu expected, %lu determined)."),
       (unsigned long)(d->n_part_ents), (unsigned long)recv_size);

  /* Allocate distributor arrays */

  BFTC_MALLOC(d->send_list, d->send_size, fvmc_lnum_t);
  BFTC_MALLOC(d->recv_order, d->n_part_ents, fvmc_lnum_t);

  BFTC_MALLOC(d->_recv_global_num, d->n_part_ents, fvmc_gnum_t);
  d->recv_global_num = d->_recv_global_num;

  /* We already have all the necessary info to build the global numbering */

  for (j = 0; j < d->n_part_ents; j++)
    d->_recv_global_num[j] = adjacency[adj_list[j]];

  /* Prepare destination rank request and receive_order at the same time
     (temporarily modifying d->recv_displ) */

  BFTC_MALLOC(send_num, d->send_size, fvmc_gnum_t);
  BFTC_MALLOC(recv_num, d->n_part_ents, fvmc_gnum_t);

  for (j = 0; j < d->n_part_ents; j++) {
    fvmc_gnum_t adj_g_num = adjacency[adj_list[j]];
    int adj_ent_rank =   ((adj_g_num-1) / adjacent_block.block_size)
                       * adjacent_block.rank_step;
    recv_num[d->recv_displ[adj_ent_rank]] = adj_g_num;
    d->recv_order[j] = d->recv_displ[adj_ent_rank];
    d->recv_displ[adj_ent_rank] += 1;
  }

  for (i = 0; i < n_ranks; i++)
    d->recv_displ[i] -= d->recv_count[i];

  BFTC_FREE(_adj_list);
  adj_list = NULL;

  MPI_Alltoallv(recv_num, d->recv_count, d->recv_displ, FVMC_MPI_GNUM,
                send_num, d->send_count, d->send_displ, FVMC_MPI_GNUM,
                d->comm);

  BFTC_FREE(recv_num);

  /* Now prepare send list */

  for (j = 0; j < d->send_size; j++)
    d->send_list[j] = send_num[j] - adjacent_block.gnum_range[0];

  BFTC_FREE(send_num);

  /* Return initialized structure */

  return d;
}

/*----------------------------------------------------------------------------
 * Initialize block to partition distributor based global element numbers
 * for partitioned data.
 *
 * This function is similar to fvmc_block_to_part_create_adj(), but is
 * intended to be used with already partitioned data, which may have
 * been renumbered; as such, it does not sort or remove duplicates from the
 * list of global elements, and global entity numbers must be nonzero.
 *
 * arguments:
 *   comm           <-- communicator
 *   bi             <-- block size and range info
 *   n_part_ents    <-- number of elements in partition
 *   global_ent_num <-- global entity numbers (in partition)
 *
 * returns:
 *   initialized partition to block distributor
 *----------------------------------------------------------------------------*/

fvmc_block_to_part_t *
fvmc_block_to_part_create_by_gnum(MPI_Comm                   comm,
                                 fvmc_block_to_part_info_t   bi,
                                 fvmc_lnum_t                 n_part_ents,
                                 const fvmc_gnum_t           global_ent_num[])
{
  int i;
  size_t j;

  size_t recv_size = 0;
  fvmc_gnum_t *send_num = NULL, *recv_num = NULL;

  fvmc_block_to_part_t *d = _block_to_part_create(comm);

  const int n_ranks = d->n_ranks;

  d->n_part_ents = n_part_ents;

  /* Use adjacent global number to determine to which ranks they will
     be sent. Send and receive counts are exchanged here, as ranks
     requesting (i.e receiving entities) are determined first */

  for (i = 0; i < d->n_ranks; i++)
    d->recv_count[i] = 0;

  for (j = 0; j < d->n_part_ents; j++) {
    fvmc_gnum_t ent_g_id = global_ent_num[j] - 1;
    int ent_ent_rank =   (ent_g_id / bi.block_size)
                       * bi.rank_step;
    d->recv_count[ent_ent_rank] += 1;
  }

  MPI_Alltoall(d->recv_count, 1, MPI_INT, d->send_count, 1, MPI_INT, comm);

  d->send_size = _compute_displ(n_ranks, d->send_count, d->send_displ);
  recv_size = _compute_displ(n_ranks, d->recv_count, d->recv_displ);

  if (d->n_part_ents != recv_size)
    bftc_error
      (__FILE__, __LINE__, 0,
       _("inconsistent sizes computed for a block to partition distributor\n"
         "(%lu expected, %lu determined)."),
       (unsigned long)(d->n_part_ents), (unsigned long)recv_size);

  /* Allocate distributor arrays */

  BFTC_MALLOC(d->send_list, d->send_size, fvmc_lnum_t);
  BFTC_MALLOC(d->recv_order, d->n_part_ents, fvmc_lnum_t);

  BFTC_MALLOC(d->_recv_global_num, d->n_part_ents, fvmc_gnum_t);
  d->recv_global_num = d->_recv_global_num;

  /* We already have all the necessary info to build the global numbering */

  for (j = 0; j < d->n_part_ents; j++)
    d->_recv_global_num[j] = global_ent_num[j];

  /* Prepare destination rank request and receive_order at the same time
     (temporarily modifying d->recv_displ) */

  BFTC_MALLOC(send_num, d->send_size, fvmc_gnum_t);
  BFTC_MALLOC(recv_num, d->n_part_ents, fvmc_gnum_t);

  for (j = 0; j < d->n_part_ents; j++) {
    fvmc_gnum_t ent_g_num = global_ent_num[j];
    int ent_ent_rank =   ((ent_g_num-1) / bi.block_size)
                       * bi.rank_step;
    recv_num[d->recv_displ[ent_ent_rank]] = ent_g_num;
    d->recv_order[j] = d->recv_displ[ent_ent_rank];
    d->recv_displ[ent_ent_rank] += 1;
  }

  for (i = 0; i < n_ranks; i++)
    d->recv_displ[i] -= d->recv_count[i];

  MPI_Alltoallv(recv_num, d->recv_count, d->recv_displ, FVMC_MPI_GNUM,
                send_num, d->send_count, d->send_displ, FVMC_MPI_GNUM,
                d->comm);

  BFTC_FREE(recv_num);

  /* Now prepare send list */

  for (j = 0; j < d->send_size; j++)
    d->send_list[j] = send_num[j] - bi.gnum_range[0];

  BFTC_FREE(send_num);

  /* Return initialized structure */

  return d;
}

/*----------------------------------------------------------------------------
 * Destroy a block to partition distributor structure.
 *
 * arguments:
 *   d <-> pointer to block to partition distributor structure pointer
 *----------------------------------------------------------------------------*/

void
fvmc_block_to_part_destroy(fvmc_block_to_part_t **d)
{
  fvmc_block_to_part_t *_d = *d;

  BFTC_FREE(_d->send_count);
  BFTC_FREE(_d->recv_count);
  BFTC_FREE(_d->send_displ);
  BFTC_FREE(_d->recv_displ);

  BFTC_FREE(_d->send_list);
  BFTC_FREE(_d->recv_order);

  BFTC_FREE(_d->_recv_global_num);

  BFTC_FREE(*d);
}

/*----------------------------------------------------------------------------
 * Return number of entities associated with local partition
 *
 * arguments:
 *   d <-- distribtor helper
 *
 * returns:
 *   number of entities associated with distribution receive
 *----------------------------------------------------------------------------*/

fvmc_lnum_t
fvmc_block_to_part_get_n_part_ents(fvmc_block_to_part_t *d)
{
  fvmc_lnum_t retval = 0;

  if (d != NULL)
    retval = d->n_part_ents;

  return retval;
}

/*----------------------------------------------------------------------------
 * Transfer a block to partition distributor's associated global numbering.
 *
 * The pointer to the global number array is returned, and ownership
 * of this array is given to the caller.
 *
 * arguments:
 *   d <-> pointer to block to partition distributor structure pointer
 *
 * returns:
 *   pointer to receiver global numbering, or NULL if the block to
 *   domain partition distributor was not the owner of this array.
 *----------------------------------------------------------------------------*/

fvmc_gnum_t *
fvmc_block_to_part_transfer_gnum(fvmc_block_to_part_t *d)
{
  fvmc_gnum_t *retval = d->_recv_global_num;

  d->_recv_global_num = NULL;

  return retval;
}

/*----------------------------------------------------------------------------
 * Copy array data from block distribution to general domain partition.
 *
 * arguments:
 *   d            <-- block to partition distributor
 *   datatype     <-- type of data considered
 *   stride       <-- number of values per entity (interlaced)
 *   block_values --> values in block distribution
 *   part_values  --> values in general domain partition
 *----------------------------------------------------------------------------*/

void
fvmc_block_to_part_copy_array(fvmc_block_to_part_t   *d,
                             fvmc_datatype_t         datatype,
                             int                    stride,
                             const void            *block_values,
                             void                  *part_values)
{
  int        i;
  size_t     j, k;

  unsigned char *send_buf = NULL;
  unsigned char *recv_buf = NULL;

  size_t stride_size = fvmc_datatype_size[datatype]*stride;
  MPI_Datatype mpi_type = fvmc_datatype_to_mpi[datatype];

  const unsigned char *_block_values = (const unsigned char *) block_values;
  unsigned char *_part_values = (unsigned char *) part_values;

  const int n_ranks = d->n_ranks;
  const size_t _send_size = d->send_size;
  const size_t _n_recv_ents = d->n_part_ents;

  /* Adjust send and receive dimensions */
  if (stride > 1) {
    for (i = 0; i < n_ranks; i++) {
      d->send_count[i] *= stride;
      d->recv_count[i] *= stride;
      d->send_displ[i] *= stride;
      d->recv_displ[i] *= stride;
    }
  }

  /* Prepare MPI buffers */

  BFTC_MALLOC(send_buf, _send_size*stride_size, unsigned char);

  for (j = 0; j < _send_size; j++) {

    size_t w_displ = j*stride_size;
    size_t r_displ = d->send_list[j]*stride_size;

    for (k = 0; k < stride_size; k++)
      send_buf[w_displ + k] = _block_values[r_displ + k];
  }

  BFTC_MALLOC(recv_buf, d->n_part_ents*stride_size, unsigned char);

  /* Exchange values */

  MPI_Alltoallv(send_buf, d->send_count, d->send_displ, mpi_type,
                recv_buf, d->recv_count, d->recv_displ, mpi_type,
                d->comm);

  /* Distribute received values */

  for (j = 0; j < _n_recv_ents; j++) {

    size_t w_displ = j*stride_size;
    size_t r_displ = d->recv_order[j]*stride_size;

    for (k = 0; k < stride_size; k++)
      _part_values[w_displ + k] = recv_buf[r_displ + k];
  }

  BFTC_FREE(recv_buf);
  BFTC_FREE(send_buf);

  /* Reset send and receive dimensions */

  if (stride > 1) {
    for (i = 0; i < n_ranks; i++) {
      d->send_count[i] /= stride;
      d->recv_count[i] /= stride;
      d->send_displ[i] /= stride;
      d->recv_displ[i] /= stride;
    }
  }

}

/*----------------------------------------------------------------------------
 * Copy a local index from block distribution to general domain partition.
 *
 * This is useful for distribution of entity connectivity information.
 *
 * arguments:
 *   d          <-- block to partition distributor
 *   send_index <-- local index in block distribution
 *   recv_index --> local index in general partition distribution
 *                  (size: n_part_entities + 1)
 *----------------------------------------------------------------------------*/

void
fvmc_block_to_part_copy_index(fvmc_block_to_part_t  *d,
                             const fvmc_lnum_t     *block_index,
                             fvmc_lnum_t           *part_index)

{
  size_t i;

  fvmc_lnum_t *send_recv_size = NULL;
  fvmc_lnum_t *send_ent_size = NULL, *recv_ent_size = NULL;

  /* Convert send index to count, then exchange */

  BFTC_MALLOC(send_recv_size, d->send_size + d->n_part_ents, fvmc_lnum_t);
  send_ent_size = send_recv_size;
  recv_ent_size = send_recv_size + d->send_size;

  for (i = 0; i < d->send_size; i++) {
    size_t ent_id = d->send_list[i];
    send_ent_size[i] = block_index[ent_id+1] - block_index[ent_id];
  }

  /* Exchange entity sizes */

  MPI_Alltoallv(send_ent_size, d->send_count, d->send_displ, FVMC_MPI_LNUM,
                recv_ent_size, d->recv_count, d->recv_displ, FVMC_MPI_LNUM,
                d->comm);

  send_ent_size = NULL;

  /* Build received index */

  if (part_index != NULL) {
    part_index[0] = 0;
    for (i = 0; i < d->n_part_ents; i++)
      part_index[i+1] = part_index[i] + recv_ent_size[d->recv_order[i]];
  }

  /* Free temporary memory */

  recv_ent_size = NULL;
  BFTC_FREE(send_recv_size);
}

/*----------------------------------------------------------------------------
 * Copy indexed data from block distribution to general domain partition.
 *
 * This is useful for distribution of entity connectivity information.
 *
 * arguments:
 *   d           <-- block to partition distributor
 *   datatype    <-- type of data considered
 *   block_index <-- local index in block distribution
 *   block_val   <-- values in block distribution
 *                   (size: send_index[n_block_ents])
 *   part_index  --> local index in general distribution
 *   part_val    --> numbers in general  distribution
 *                   (size: recv_index[n_part_ents])
 *----------------------------------------------------------------------------*/

void
fvmc_block_to_part_copy_indexed(fvmc_block_to_part_t   *d,
                               fvmc_datatype_t         datatype,
                               const fvmc_lnum_t      *send_index,
                               const void            *send_val,
                               const fvmc_lnum_t      *recv_index,
                               void                  *recv_val)
{
  int    i;
  size_t j, k, w_displ, r_displ;

  size_t send_size = 0;
  size_t recv_size = 0;

  int  *send_count = NULL;
  int  *recv_count = NULL;
  int  *send_displ = NULL;
  int  *recv_displ = NULL;

  fvmc_lnum_t  *inv_order = NULL;
  size_t  *recv_val_index = NULL;

  unsigned char *send_buf = NULL;
  unsigned char *recv_buf = NULL;

  const unsigned char *_send_val = (const unsigned char *) send_val;
  unsigned char *_recv_val = (unsigned char *) recv_val;

  size_t type_size = fvmc_datatype_size[datatype];
  MPI_Datatype mpi_type = fvmc_datatype_to_mpi[datatype];

  const int n_ranks = d->n_ranks;

  /* Build send and receive counts */
  /*-------------------------------*/

  BFTC_MALLOC(send_count, n_ranks, int);
  BFTC_MALLOC(recv_count, n_ranks, int);
  BFTC_MALLOC(send_displ, n_ranks, int);
  BFTC_MALLOC(recv_displ, n_ranks, int);

  for (i = 0; i < n_ranks; i++) {
    send_count[i] = 0;
    recv_count[i] = 0;
  }

  for (i = 0; i < n_ranks; i++) {
    const size_t start_id = d->send_displ[i];
    const size_t end_id = start_id + d->send_count[i];
    for (j = start_id; j < end_id; j++) {
      const size_t send_id = d->send_list[j];
      send_count[i] += send_index[send_id+1] - send_index[send_id];
    }
  }

  BFTC_MALLOC(inv_order, d->n_part_ents, fvmc_lnum_t);

  for (j = 0; j < d->n_part_ents; j++)
    inv_order[d->recv_order[j]] = j;

  for (i = 0; i < n_ranks; i++) {
    const size_t start_id = d->recv_displ[i];
    const size_t end_id = start_id + d->recv_count[i];
    for (j = start_id; j < end_id; j++) {
      const size_t recv_id = inv_order[j];
      recv_count[i] += recv_index[recv_id+1] - recv_index[recv_id];
    }
  }

  BFTC_FREE(inv_order);

  send_size = _compute_displ(n_ranks, send_count, send_displ);
  recv_size = _compute_displ(n_ranks, recv_count, recv_displ);

  /* Build send and receive buffers */
  /*--------------------------------*/

  BFTC_MALLOC(send_buf, send_size * type_size, unsigned char);
  BFTC_MALLOC(recv_buf, recv_size * type_size, unsigned char);

  w_displ = 0;
  r_displ = 0;

  for (j = 0; j < d->send_size; j++) {

    size_t ent_id = d->send_list[j];
    size_t ent_size =  (send_index[ent_id + 1] - send_index[ent_id])
                      * type_size;
    r_displ = send_index[ent_id] * type_size;

    for (k = 0; k < ent_size; k++)
      send_buf[w_displ++] = _send_val[r_displ++];
  }

  assert(w_displ == send_size*type_size);

  /* Exchange values */

  MPI_Alltoallv(send_buf, send_count, send_displ, mpi_type,
                recv_buf, recv_count, recv_displ, mpi_type, d->comm);

  BFTC_FREE(send_buf);
  BFTC_FREE(send_count);
  BFTC_FREE(send_displ);

  BFTC_FREE(recv_count);
  BFTC_FREE(recv_displ);

  /* Retrieve values */
  /*-----------------*/

  BFTC_MALLOC(recv_val_index, d->n_part_ents + 1, size_t);

  /* Build receive buffer index */

  recv_val_index[0] = 0;

  for (j = 0; j < d->n_part_ents; j++)
    recv_val_index[d->recv_order[j] + 1] = (  recv_index[j+1]
                                            - recv_index[j])*type_size;

  for (j = 0; j < d->n_part_ents; j++)
    recv_val_index[j+1] += recv_val_index[j];

  w_displ = 0;

  for (j = 0; j < d->n_part_ents; j++) {
    const size_t recv_id = d->recv_order[j];
    const size_t ent_size =   (recv_index[recv_id+1] - recv_index[recv_id])
                            * type_size;
    r_displ = recv_val_index[recv_id];
    for (k = 0; k < ent_size; k++)
      _recv_val[w_displ++] = recv_buf[r_displ++];
  }

  BFTC_FREE(recv_buf);
  BFTC_FREE(recv_val_index);
}

#endif /* defined(FVMC_HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Determine local references from references to global numbers.
 *
 * This is based on finding the local id of a given global number
 * in a sorted global list using a binary search.
 *
 * Global numbers use a 1 to n numbering, while local numbers use a
 * 0+base to n-1+base numbering. If an entity's global number does not
 * appear in the global list, base-1 is assigned for that entity's
 * local list.
 *
 * If the sorted list contains duplicate values, any local id having
 * a multiple global number (i.e not necessarily the smallest one)
 * may be assigned to the corresponding local_number[] entry.
 *
 * arguments:
 *   n_ents           <-- number of entities
 *   base             <-- base numbering (typically 0 or 1)
 *   global_list_size <-- size of global entity list
 *   global_list      <-- global entity list
 *   global_number    <-- entity global numbers
 *                        (size: n_ents)
 *   local_number     --> entity local numbers
 *                        (size: n_ents)
 *----------------------------------------------------------------------------*/

void
fvmc_block_to_part_global_to_local(fvmc_lnum_t        n_ents,
                                  fvmc_lnum_t        base,
                                  fvmc_lnum_t        global_list_size,
                                  const fvmc_gnum_t  global_list[],
                                  const fvmc_gnum_t  global_number[],
                                  fvmc_lnum_t        local_number[])
{
  fvmc_lnum_t i;

  if (n_ents == 0)
    return;

 #if defined(DEBUG) && !defined(NDEBUG)
  for (i = 1; i < global_list_size; i++)
    assert(global_list[i] > global_list[i-1]);
#endif

  for (i = 0; i < n_ents; i++) {

    fvmc_lnum_t start_id = 0;
    fvmc_lnum_t end_id = global_list_size;

    const fvmc_gnum_t num_1 = global_number[i];

    /* Use binary search */

    while (start_id < end_id) {
      fvmc_lnum_t mid_id = start_id + ((end_id - start_id) / 2);
      if (global_list[mid_id] < num_1)
        start_id = mid_id + 1;
      else
        end_id = mid_id;  /* Can't be end_id = mid_id -1;
                             global_list[mid_id] >= num_1, so
                             end_id must not be < mid_id in case
                             global_list[mid_id] == num_1 */
    }

    /* start_id == end_id at this stage; */

    if (start_id < global_list_size && global_list[start_id] == num_1)
      local_number[i] = start_id + base;
    else
      local_number[i] = base - 1;
  }
}

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

