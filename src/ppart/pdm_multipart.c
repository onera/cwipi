/*
  This file is part of the ParaDiGM library.

  Copyright (C) 2017       ONERA

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

/*============================================================================
 * TODO : write module description here
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"
#include "pdm_part.h"
#include "pdm_handles.h"
#include "pdm_dmesh.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_binary_search.h"
#include "pdm_part_to_block.h"
#include "pdm_block_to_part.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_multipart.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Local structure definitions
 *============================================================================*/

typedef struct  {
  int  n_bound;
  int  n_join;
  int *joins_ids;
  int *face_bound_idx;
  int *face_join_idx;
  int *face_bound;
  int *face_join;
  PDM_g_num_t *face_bound_ln_to_gn;
  PDM_g_num_t *face_join_ln_to_gn;

} _bounds_and_joins_t;

/**
 * \struct _pdm_multipart_t
 * \brief  This structure describe a multipart. In addition to splitting
 *         parameters, it stores the multiples blocks and part as well as
 *         the global numbering.
 *
 */

typedef struct  {

  int               n_zone;           /*!< Number of initial zones */
  const int        *n_part;          /*!< Number of partitions per proc in each zone */
  PDM_bool_t        merge_blocks;     /*!< Merge before partitionning or not */
  PDM_part_split_t  split_method;     /*!< Partitioning method */
  PDM_MPI_Comm      comm;             /*!< MPI communicator */
  int               *part_zone_distribution;  /*!< Number of part in each zone (distribution,
                                           size = n_zone + 1)                         */
  int               *gpart_to_proc;     /*!< For each global part id, proc storing this
                                           part and localId of part in this process   */
  int               *dmeshes_ids;      /*!< Ids of distributed blocks (size = n_zone)  */
  int               *part_ids;         /*!< Ids of partitions built on each block of this
                                           process (size = n_zone)                    */
  int               *n_bounds_and_joins; /*!< Number of boundaries and joins in each zone
                                           (size = 2*n_zone, global data)             */
  _bounds_and_joins_t **pbounds_and_joins;/*!< partitionned boundary and join data in each
                                           zone/part                                  */
  int               n_total_joins;      /*Total number of joins between zones
                                          (nb : each counts twice)                    */
  int               *join_to_opposite;  /*For each global joinId, give the globalId of
                                          the opposite join (size = n_total_joins)    */
} _pdm_multipart_t;

/*============================================================================
 * Global variable
 *============================================================================*/

static PDM_Handles_t *_multiparts   = NULL;

/*=============================================================================
 * Private function definitions
 *============================================================================*/

/**
 *
 * \brief Return multipart object from it identifier
 *
 * \param [in]   multipartId    multipart identifier
 *
 */

static _pdm_multipart_t *
_get_from_id
(
 int  id
)
{

  _pdm_multipart_t *multipart = (_pdm_multipart_t *) PDM_Handles_get (_multiparts, id);

  if (multipart == NULL) {
    PDM_error(__FILE__, __LINE__, 0, "PDM_multipart error : Bad identifier\n");
  }

  return multipart;
}

// static int
//  PDM_search_rank
// (
//  PDM_g_num_t   elt,
//  PDM_g_num_t  *array,
//  int            id1,
//  int            id2
// )
// {
//   if (elt >= array[id2]) {
//     PDM_printf("PPART error : Element not in initial distributed array "
//            PDM_FMT_G_NUM" "PDM_FMT_G_NUM" "PDM_FMT_G_NUM"\n",
//            elt, array[id1], array[id2]);
//     abort();
//   }

//   if (elt < array[id1]) {
//     PDM_printf("PPART error : Element not in initial distributed array "
//            PDM_FMT_G_NUM" "PDM_FMT_G_NUM" "PDM_FMT_G_NUM"\n",
//            elt, array[id1], array[id2]);
//     abort();
//   }

//   if (id2 == id1 + 1) {
//     return id1;
//   }

//   else {

//     while(array[id1] == array[id1+1]) id1++;

//     int midId = (id2 + id1) / 2;

//     if (elt == array[id1])
//       return id1;
//     else if (elt < array[midId])
//       return PDM_search_rank(elt, array, id1, midId);
//     else if (elt >= array[midId])
//       return PDM_search_rank(elt, array, midId, id2);
//   }
//   return -1;
// }
static void
_build_join_uface_distribution
(
 _pdm_multipart_t *_multipart,
 int              *join_to_ref_join,   //Size n_total_join
 int              *face_in_join_distri //Size n_unique_joins + 1
)
{

  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(_multipart->comm, &i_rank);
  PDM_MPI_Comm_size(_multipart->comm, &n_rank);

  PDM_printf("pdm::_build_join_uface_distribution\n");
  int n_total_joins  = _multipart->n_total_joins;
  int n_unique_joins = n_total_joins/2;

  //Build join_to_ref_join : we want the join and opposite join to have the same shift index,
  // so we take the smaller join global id as the reference
  int ref_join_gid = 0;
  for (int ijoin = 0; ijoin < n_total_joins; ijoin++)
  {
    int opp_join = _multipart->join_to_opposite[ijoin];
    if (ijoin < opp_join)
    {
      join_to_ref_join[ijoin] = ref_join_gid;
      join_to_ref_join[opp_join] = ref_join_gid;
      ref_join_gid ++;
    }
  }
  /*
  PDM_printf("Join to reference join :");
  for (int ijoin = 0; ijoin < n_total_joins; ijoin++)
   PDM_printf(" %d ", join_to_ref_join[ijoin]);
  PDM_printf("\n");
  */

  //Count faces in joins
  int *nb_face_in_joins = (int *) malloc(n_unique_joins * sizeof(int));
  for (int i = 0; i < n_unique_joins; i++)
    nb_face_in_joins[i] = 0;

  for (int izone = 0; izone < _multipart->n_zone; izone++)
  {
    int block_id = _multipart->dmeshes_ids[izone];
    int dn_cell  = 0;
    int dn_face  = 0;
    int dn_vtx   = 0;
    int n_bnd    = 0;
    int n_join   = 0;
    const double       *dvtx_coord;
    const int          *dface_vtx_idx;
    const PDM_g_num_t  *dface_vtx;
    const PDM_g_num_t  *dface_cell;
    const int          *dface_bound_idx;
    const PDM_g_num_t  *dface_bound;
    const int          *djoin_gids;
    const int          *dface_join_idx;
    const PDM_g_num_t  *dface_join;

    PDM_dmesh_dims_get(block_id, &dn_cell, &dn_face, &dn_vtx, &n_bnd, &n_join);
    PDM_dmesh_data_get(block_id, &dvtx_coord, &dface_vtx_idx, &dface_vtx, &dface_cell,
                       &dface_bound_idx, &dface_bound, &djoin_gids, &dface_join_idx, &dface_join);
    for (int ijoin=0; ijoin < n_join; ijoin ++)
    {
      int join_gid = djoin_gids[2*ijoin];
      int join_opp_gid = djoin_gids[2*ijoin + 1];
      //Paired joins must be counted only once
      if (join_gid < join_opp_gid)
        nb_face_in_joins[join_to_ref_join[join_gid-1]] = dface_join_idx[ijoin+1] - dface_join_idx[ijoin];
    }
  }
  /*
  PDM_printf("[%d] nb_face_joins : ", i_rank);
  for (int i = 0; i < n_unique_joins ; i++)
    PDM_printf(" %d ", nb_face_in_joins[i]);
  PDM_printf("\n");
  */

  //Sum faces and build distribution
  PDM_MPI_Allreduce(nb_face_in_joins, &face_in_join_distri[1], n_unique_joins, PDM_MPI_INT, PDM_MPI_SUM, _multipart->comm);

  face_in_join_distri[0] = 0;
  for (int i=0; i < n_unique_joins; i++)
    face_in_join_distri[i+1] = face_in_join_distri[i+1] + face_in_join_distri[i];

  /*
  PDM_printf("[%d] face_in_join_distri : ", i_rank);
  for (int i = 0; i < n_unique_joins + 1; i++)
    PDM_printf(" %d ", face_in_join_distri[i]);
  PDM_printf("\n");
  */

  free(nb_face_in_joins);
}

static void
_set_dface_tag_from_joins
(
 const int          dn_face,
 const int          n_join,
 const int         *djoin_gids,
 const int         *dface_join_idx,
 const PDM_g_num_t *dface_join,
 int               *dface_tag,
 const PDM_MPI_Comm  comm
)
{
  int i_rank;
  int n_rank;

  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  // 1. Construct face distribution -- will be needed to find owner of faces
  PDM_g_num_t * dface_proc = (PDM_g_num_t *) malloc((n_rank+1) * sizeof(PDM_g_num_t));
  int *dn_face_proc = (int *) malloc((n_rank) * sizeof(int));

  PDM_MPI_Allgather((void *) &dn_face, 1, PDM_MPI_INT, (void *) dn_face_proc, 1, PDM_MPI_INT, comm);

  dface_proc[0] = 1;
  for (int i = 1; i < n_rank+1; i++)
    dface_proc[i] = (PDM_g_num_t) dn_face_proc[i-1] + dface_proc[i-1];
  free(dn_face_proc);

  int  nData = 2; //Face Id, JoinId
  int *face_to_send_n   = (int *) malloc(n_rank * sizeof(int));
  int *face_to_send_idx = (int *) malloc((n_rank+1) * sizeof(int));
  for (int i = 0; i < n_rank; i++)
    face_to_send_n[i] = 0;

  // 2. Prepare and send data
  //Count faces to send
  for (int ijoin = 0; ijoin < n_join; ijoin++) {
    for (int iface = dface_join_idx[ijoin]; iface < dface_join_idx[ijoin+1]; iface++) {
      int rank = PDM_search_rank(dface_join[iface], dface_proc, 0, n_rank);
      face_to_send_n[rank] += nData;
    }
  }
  //Prepare variable stride
  face_to_send_idx[0] = 0;
  for (int i = 1; i < n_rank + 1; i++) {
    face_to_send_idx[i] = face_to_send_idx[i-1] + face_to_send_n[i-1];
    face_to_send_n[i-1] = 0;
  }
  //Prepare data
  PDM_g_num_t *face_to_send = (PDM_g_num_t *) malloc(face_to_send_idx[n_rank] * sizeof(PDM_g_num_t));
  for (int ijoin = 0; ijoin < n_join; ijoin++) {
    for (int iface = dface_join_idx[ijoin]; iface < dface_join_idx[ijoin+1]; iface++) {
      int rank = PDM_search_rank(dface_join[iface], dface_proc, 0, n_rank);
      int idx   = face_to_send_idx[rank] + face_to_send_n[rank];
      face_to_send[idx  ]   = dface_join[iface];
      face_to_send[idx+1]   = djoin_gids[2*ijoin];
      face_to_send_n[rank] += nData;
    }
  }
  //Exchange sizes
  int *face_to_recv_n   = (int *) malloc(n_rank * sizeof(int));
  PDM_MPI_Alltoall(face_to_send_n, 1, PDM_MPI_INT, face_to_recv_n, 1, PDM_MPI_INT, comm);
  int *face_to_recv_idx = (int *) malloc((n_rank+1) * sizeof(int));
  face_to_recv_idx[0] = 0;
  for(int i = 1; i < (n_rank+1); i++) {
    face_to_recv_idx[i] = face_to_recv_idx[i-1] + face_to_recv_n[i-1];
  }
  //Exchange data
  PDM_g_num_t *face_to_recv = (PDM_g_num_t *) malloc(face_to_recv_idx[n_rank]*sizeof(PDM_g_num_t));
  PDM_MPI_Alltoallv(face_to_send,
                    face_to_send_n,
                    face_to_send_idx,
                    PDM__PDM_MPI_G_NUM,
                    face_to_recv,
                    face_to_recv_n,
                    face_to_recv_idx,
                    PDM__PDM_MPI_G_NUM,
                    comm);
  int nRecv = face_to_recv_idx[n_rank]/nData;

  free(face_to_send_n);
  free(face_to_send_idx);
  free(face_to_send);
  free(face_to_recv_n);
  free(face_to_recv_idx);

  // 3. Process received data : go back to local numerotation and flag received faces
  for (int iface = 0; iface < dn_face; iface++){
    dface_tag[iface] = -1;
  }
  for (int i=0; i<nRecv; i++) {
    int lfaceId = face_to_recv[nData*i] - dface_proc[i_rank];
    dface_tag[lfaceId] = face_to_recv[nData*i + 1];
  }
  free(face_to_recv);
}


static void
_split_bounds_and_joins
(
 _pdm_multipart_t *_multipart
)
{

  printf("_split_bounds_and_joins::\n");

  //Set structure : we need a to retrive pboundsAndJoin for a given zone/part
  int *bounds_and_joins_idx = (int *) malloc((_multipart->n_zone + 1) * sizeof(int));
  bounds_and_joins_idx[0] = 0;
  for (int i = 0; i < _multipart->n_zone; i++) {
    bounds_and_joins_idx[i + 1] = _multipart->n_part[i] + bounds_and_joins_idx[i];
  }

  _multipart->pbounds_and_joins = (_bounds_and_joins_t **)
  malloc(bounds_and_joins_idx[_multipart->n_zone] * sizeof(_bounds_and_joins_t *));

  // Loop over zones and part to get data
  for (int zone_gid = 0; zone_gid<_multipart->n_zone; zone_gid++) {
    for (int i_part = 0; i_part < _multipart->n_part[zone_gid]; i_part++) {
      int n_cell, n_face, n_face_part_bound, n_vtx, n_proc, n_total_part, scell_face, sface_vtx, sface_group, n_face_group;
      PDM_part_part_dim_get(_multipart->part_ids[zone_gid],
                        i_part,
                        &n_cell,
                        &n_face,
                        &n_face_part_bound,
                        &n_vtx,
                        &n_proc,
                        &n_total_part,
                        &scell_face,
                        &sface_vtx,
                        &sface_group,
                        &n_face_group);

      int n_bound = _multipart->n_bounds_and_joins[2*zone_gid];
      int n_join  = _multipart->n_bounds_and_joins[2*zone_gid+1];
      assert(n_face_group == n_bound + n_join);

      int          *cell_tag;
      int          *cell_face_idx;
      int          *cell_face;
      PDM_g_num_t *cell_ln_to_gn;
      int          *face_tag;
      int          *face_cell;
      int          *face_vtx_idx;
      int          *face_vtx;
      PDM_g_num_t *face_ln_to_gn;
      int          *face_part_bound_proc_idx;
      int          *face_part_bound_part_idx;
      int          *face_part_bound;
      int          *vtx_tag;
      double       *vtx;
      PDM_g_num_t  *vtx_ln_to_gn;
      int          *face_group_idx;
      int          *face_group;
      PDM_g_num_t *face_group_ln_to_gn;
      PDM_part_part_val_get(_multipart->part_ids[zone_gid],
                        i_part,
                        &cell_tag,
                        &cell_face_idx,
                        &cell_face,
                        &cell_ln_to_gn,
                        &face_tag,
                        &face_cell,
                        &face_vtx_idx,
                        &face_vtx,
                        &face_ln_to_gn,
                        &face_part_bound_proc_idx,
                        &face_part_bound_part_idx,
                        &face_part_bound,
                        &vtx_tag,
                        &vtx,
                        &vtx_ln_to_gn,
                        &face_group_idx,
                        &face_group,
                        &face_group_ln_to_gn
                        );

      //Retrieve boundaries and joins from face_group
      int *pface_bound_idx = (int *) malloc((n_bound+1) * sizeof(int));
      int *pface_join_idx  = (int *) malloc((n_join +1) * sizeof(int));
      for (int i = 0; i < n_bound + 1; i++)
        pface_bound_idx[i] = face_group_idx[i];
      pface_join_idx[0] = 0;
      for (int i = n_bound + 1; i < n_bound + n_join + 1; i++)
        pface_join_idx[i-n_bound] = face_group_idx[i] - face_group_idx[n_bound];

      int *pface_bound = (int *) malloc(pface_bound_idx[n_bound] * sizeof(int));
      int *pface_join  = (int *) malloc(4*pface_join_idx[n_join]   * sizeof(int));
      for (int i = 0; i < pface_bound_idx[n_bound]; i++)
        pface_bound[i] = face_group[i];
      for (int i = pface_bound_idx[n_bound]; i < face_group_idx[n_face_group]; i++)
        pface_join[4*(i - pface_bound_idx[n_bound])] = face_group[i];

      PDM_g_num_t *pface_bound_ln_to_gn = (PDM_g_num_t *) malloc(pface_bound_idx[n_bound] * sizeof(PDM_g_num_t));
      PDM_g_num_t *pface_join_ln_to_gn  = (PDM_g_num_t *) malloc(pface_join_idx[n_join]   * sizeof(PDM_g_num_t));
      for (int i = 0; i < pface_bound_idx[n_bound]; i++)
        pface_bound_ln_to_gn[i] = face_group_ln_to_gn[i];
      for (int i = pface_bound_idx[n_bound]; i < face_group_idx[n_face_group]; i++)
        pface_join_ln_to_gn[i - pface_bound_idx[n_bound]] = face_group_ln_to_gn[i];

      // Retrieve joinId using facetag
      int *pjoins_ids      = (int *) malloc(n_join * sizeof(int)); //Faut il init ?
      for (int ijoin = 0; ijoin < n_join; ijoin++) {
        if (pface_join_idx[ijoin] != pface_join_idx[ijoin + 1]) {
          //TODO : coherence des ordinaux joins (démarrage a 0 pour joinOpp, mais à 1 dans le dmesh ...)
          int first_face_lid = pface_join[4*pface_join_idx[ijoin]];
          int joinId         = face_tag[first_face_lid-1];
          pjoins_ids[ijoin] = joinId - 1;
        }
      }

      // Store data in pbounds_and_joins
      int idx = bounds_and_joins_idx[zone_gid] + i_part;
      _multipart->pbounds_and_joins[idx] = malloc(sizeof(_bounds_and_joins_t));
      _multipart->pbounds_and_joins[idx]->n_bound         = n_bound;
      _multipart->pbounds_and_joins[idx]->n_join          = n_join;
      _multipart->pbounds_and_joins[idx]->face_bound_idx    = pface_bound_idx;
      _multipart->pbounds_and_joins[idx]->face_join_idx     = pface_join_idx;
      _multipart->pbounds_and_joins[idx]->face_bound       = pface_bound;
      _multipart->pbounds_and_joins[idx]->face_join        = pface_join;
      _multipart->pbounds_and_joins[idx]->face_bound_ln_to_gn = pface_bound_ln_to_gn;
      _multipart->pbounds_and_joins[idx]->face_join_ln_to_gn  = pface_join_ln_to_gn;
      _multipart->pbounds_and_joins[idx]->joins_ids  = pjoins_ids;

    }
  }
  free(bounds_and_joins_idx);
}
static void
_search_matching_joins
(
 _pdm_multipart_t *_multipart
)
{
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(_multipart->comm, &i_rank);
  PDM_MPI_Comm_size(_multipart->comm, &n_rank);

  int *bounds_and_joins_idx = (int *) malloc((_multipart->n_zone + 1) * sizeof(int));
  bounds_and_joins_idx[0] = 0;
  for (int i = 0; i < _multipart->n_zone; i++){
    bounds_and_joins_idx[i + 1] = _multipart->n_part[i] + bounds_and_joins_idx[i];
  }

  //Construction of (unique) join distribution
  int *join_to_ref_join = (int *) malloc(_multipart->n_total_joins * sizeof(int));
  int *face_in_join_distri = (int *) malloc(((_multipart->n_total_joins)/2+1) * sizeof(int));
  _build_join_uface_distribution(_multipart, join_to_ref_join, face_in_join_distri);

  //Count total nb of join_faces
  int nb_of_joins      = 0;
  for (int izone = 0 ; izone < _multipart->n_zone; izone ++) {
    for (int i_part = 0; i_part < _multipart->n_part[izone]; i_part++) {
      int idx = bounds_and_joins_idx[izone] + i_part; //TO CHECK
      nb_of_joins += _multipart->pbounds_and_joins[idx]->n_join;
    }
  }

  // Prepare lntogn numbering and partitioned data
  PDM_g_num_t **shifted_lntogn = (PDM_g_num_t **) malloc(nb_of_joins * sizeof(PDM_g_num_t*));
  int              **part_data = (int **) malloc(nb_of_joins * sizeof(int *));
  int       *nb_face_per_join  = (int *) malloc(nb_of_joins * sizeof(int));

  int ijoin_pos  = 0;
  for (int izone = 0 ; izone < _multipart->n_zone; izone ++) {
    for (int i_part = 0; i_part < _multipart->n_part[izone]; i_part++) {
      int idx = bounds_and_joins_idx[izone] + i_part; //TO CHECK
      int  n_join           = _multipart->pbounds_and_joins[idx]->n_join;
      int *face_join_idx    = _multipart->pbounds_and_joins[idx]->face_join_idx;
      int *face_join        = _multipart->pbounds_and_joins[idx]->face_join;
      int *face_join_lntogn = _multipart->pbounds_and_joins[idx]->face_join_ln_to_gn;
      for (int ijoin = 0; ijoin < n_join; ijoin++) {
        int join_size = face_join_idx[ijoin + 1] - face_join_idx[ijoin];
        nb_face_per_join[ijoin_pos] = join_size;
        PDM_g_num_t *shifted_lntogn_loc = (PDM_g_num_t *) malloc(join_size * sizeof(PDM_g_num_t));
        int         *part_data_loc      = (int *) malloc(3 * join_size * sizeof(int));
        if (join_size != 0)
        {
          //Attention à la cohérence des gids, qui démarrent parfois à 1 parfois à 0... ici à 0
          //Get shift value from join unique distribution
          int join_gid    = _multipart->pbounds_and_joins[idx]->joins_ids[ijoin];
          int shift_value = face_in_join_distri[join_to_ref_join[join_gid]];
          int j = 0;
          //Prepare partitioned data : (PL, i_rank, i_part)
          for (int iface = face_join_idx[ijoin]; iface < face_join_idx[ijoin + 1]; iface ++) {
            shifted_lntogn_loc[j] = (PDM_g_num_t) shift_value + face_join_lntogn[iface];
            part_data_loc[3*j]    = face_join[4*iface];
            part_data_loc[3*j+1]  = i_rank;
            part_data_loc[3*j+2]  = i_part;
            j++;
          }
        }
        shifted_lntogn[ijoin_pos] = shifted_lntogn_loc;
        part_data[ijoin_pos] = part_data_loc;
        ijoin_pos += 1;
      }
    }
  }
  /*
  PDM_printf("[%d] nb_face_per_join : ", i_rank);
  for (int i = 0; i < nb_of_joins; i++)
    PDM_printf(" %d ", nb_face_per_join[i]);
  PDM_printf("\n");
  */

  //Now exchange join information using part_to_block / block_to_part
  PDM_part_to_block_t *ptb = PDM_part_to_block_create (PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                       PDM_PART_TO_BLOCK_POST_NOTHING,
                                                       1.,
                                                       shifted_lntogn,
                                                       NULL,
                                                       nb_face_per_join,
                                                       nb_of_joins,
                                                       _multipart->comm);

  PDM_g_num_t *distrib_index = PDM_part_to_block_distrib_index_get(ptb);

  /*
  PDM_printf("[%d] PTB distri : ", i_rank);
  for (int i=0; i < n_rank + 1; i++)
    PDM_printf(" %d ", distrib_index[i]);
  PDM_printf("\n");
  */

  int         *block_stride;
  int         *block_data;
  PDM_part_to_block_exch(ptb,
                         sizeof(int),
                         PDM_STRIDE_CST,
                         3,
                         NULL,
                         (void **) part_data,
                         &block_stride,
                         (void **) &block_data);

  /*
  int n_elt_block = PDM_part_to_block_n_elt_block_get(ptb);
  PDM_printf("[%d] PTB nb_elem : %d\n", i_rank, n_elt_block);
  if (i_rank == 1)
  {
    PDM_g_num_t *glob_num = PDM_part_to_block_block_gnum_get(ptb);
    PDM_printf("[%d] PTB globnum : ", i_rank);
    for (int i = 0; i < n_elt_block; i++)
      printf(" %d ", glob_num[i]);
    PDM_printf("\n");
    PDM_printf("[%d] PTB data : ", i_rank);
    for (int i = 0; i < n_elt_block; i++)
      printf(" (%d %d %d) ", block_data[3*i], block_data[3*i+1], block_data[3*i+2]);
    PDM_printf("\n");
  }
  */

  // Don't free ptb now since we need the distribution and the block_data

  PDM_block_to_part_t *btp = PDM_block_to_part_create(distrib_index,
                                                      (const PDM_g_num_t **) shifted_lntogn,
                                                      nb_face_per_join,
                                                      nb_of_joins,
                                                      _multipart->comm);

  int **new_part_data = (int **) malloc(nb_of_joins * sizeof(int *));
  for (int ijoin = 0; ijoin < nb_of_joins; ijoin ++){
    new_part_data[ijoin] = (int *) malloc(6*nb_face_per_join[ijoin] * sizeof(int));
  }
  int cst_stride = 6;

  PDM_block_to_part_exch(btp,
                         sizeof(int),
                         PDM_STRIDE_CST,
                         &cst_stride,
                         (void *) block_data,
                         NULL,
                         (void **) new_part_data);

  PDM_part_to_block_free(ptb);
  PDM_block_to_part_free(btp);

  /*
  if (i_rank == 0)
  {
    PDM_printf("[%d] BTP data : \n",  i_rank);
    for (int ijoin = 0; ijoin < nb_of_joins; ijoin++)
    {
      PDM_printf("  ijoin %d(%d) :", ijoin, nb_face_per_join[ijoin]);
      for (int iface = 0; iface < nb_face_per_join[ijoin]; iface++)
        PDM_printf(" (%d %d %d %d %d %d) ", new_part_data[ijoin][6*iface], new_part_data[ijoin][6*iface+1], new_part_data[ijoin][6*iface+2],
                                            new_part_data[ijoin][6*iface+3], new_part_data[ijoin][6*iface+4], new_part_data[ijoin][6*iface+5]);
      PDM_printf("\n");
    }
  }
  */


  //Process received data
  ijoin_pos = 0;
  for (int izone = 0 ; izone < _multipart->n_zone; izone ++) {
    for (int i_part = 0; i_part < _multipart->n_part[izone]; i_part++) {
      int idx = bounds_and_joins_idx[izone] + i_part; //TO CHECK
      int  n_join        = _multipart->pbounds_and_joins[idx]->n_join;
      int *face_join_idx = _multipart->pbounds_and_joins[idx]->face_join_idx;
      int *face_join     = _multipart->pbounds_and_joins[idx]->face_join;
      for (int ijoin = 0; ijoin < n_join; ijoin++) {
        int join_size = face_join_idx[ijoin + 1] - face_join_idx[ijoin];
        int *part_data_loc = new_part_data[ijoin_pos];
        for (int i = 0; i < join_size; i++) {
          int opp_proc = -1;
          int opp_part = -1;
          int opp_pl   = -1;
          if (part_data_loc[6*i + 1] != i_rank)
          {
            opp_proc = part_data_loc[6*i + 1];
            opp_part = part_data_loc[6*i + 2];
            opp_pl   = part_data_loc[6*i + 0];
          }
          else if (part_data_loc[6*i + 4] != i_rank)
          {
            opp_proc = part_data_loc[6*i + 4];
            opp_part = part_data_loc[6*i + 5];
            opp_pl   = part_data_loc[6*i + 3];
          }
          // The two joins are on the same proc, look at the parts
          else
          {
            opp_proc = i_rank;
            if (part_data_loc[6*i + 2] != i_part)
            {
              opp_part = part_data_loc[6*i + 2];
              opp_pl   = part_data_loc[6*i + 0];
            }
            else if (part_data_loc[6*i + 5] != i_part)
            {
              opp_part = part_data_loc[6*i + 4];
              opp_pl   = part_data_loc[6*i + 3];
            }
            // The two joins have the same proc id / part id, we need to check original pl
            else
            {
              opp_part = i_part;
              int original_pl = face_join[4*(face_join_idx[ijoin] + i)];
              if (part_data_loc[6*i] != original_pl)
                opp_pl = part_data_loc[6*i];
              else
                opp_pl = part_data_loc[6*i+3];
            }
          }
          //Fill values opp_proc, opp_part, opp_plvalue
          face_join[4*(face_join_idx[ijoin] + i) + 1] = opp_proc;
          face_join[4*(face_join_idx[ijoin] + i) + 2] = opp_part;
          face_join[4*(face_join_idx[ijoin] + i) + 3] = opp_pl;
        }
        ijoin_pos += 1;
      }
    }
  }

  //Deallocate
  for (int i = 0; i < nb_of_joins; i++) {
    free(shifted_lntogn[i]);
    free(part_data[i]);
    free(new_part_data[i]);
  }
  free(shifted_lntogn);
  free(part_data);
  free(new_part_data);
  free(face_in_join_distri);
  free(nb_face_per_join);
  free(bounds_and_joins_idx);
}

static void
_search_matching_joins_old
(
 _pdm_multipart_t *_multipart
)
{

  // Implementation avec le alltoall + tri
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(_multipart->comm, &i_rank);
  PDM_MPI_Comm_size(_multipart->comm, &n_rank);

  int *bounds_and_joins_idx = (int *) malloc((_multipart->n_zone + 1) * sizeof(int));
  bounds_and_joins_idx[0] = 0;
  for (int i = 0; i < _multipart->n_zone; i++){
    bounds_and_joins_idx[i + 1] = _multipart->n_part[i] + bounds_and_joins_idx[i];
  }



  // Step 0 : for a given join unique id, we want the list of process having a least one zone/part that
  // includes this join (after partitionning).
  int n_total_joins = _multipart->n_total_joins;

  // Array of size njoin indicating if join is known or not by this proc
  int *proc_knows_join = (int *) malloc(n_total_joins * sizeof(int));
  for (int i = 0; i < n_total_joins; i++)
    proc_knows_join[i] = 0;
  for (int izone = 0 ; izone < _multipart->n_zone; izone ++) {
    for (int i_part = 0; i_part < _multipart->n_part[izone]; i_part++) {
      int idx = bounds_and_joins_idx[izone] + i_part; //TO CHECK
      int *face_join_idx = _multipart->pbounds_and_joins[idx]->face_join_idx;
      int *join_gid      = _multipart->pbounds_and_joins[idx]->joins_ids;
      for (int ijoin = 0; ijoin < _multipart->n_bounds_and_joins[2*izone+1]; ijoin++) {
        if (face_join_idx[ijoin+1] != face_join_idx[ijoin])
          proc_knows_join[join_gid[ijoin]] = 1;
      }
    }
  }
  // To exchange data, we move to a representation proc->list of known joins
  int nb_joins_on_proc = 0;
  for (int i = 0; i < n_total_joins; i++) {
    nb_joins_on_proc += proc_knows_join[i];
  }
  int *joins_on_proc = (int *) malloc(nb_joins_on_proc * sizeof(int));
  int k = 0;
  for (int i = 0; i < n_total_joins; i++) {
    if (proc_knows_join[i] == 1) {
      joins_on_proc[k] = i;
      k++;
    }
  }
  free(proc_knows_join);

  int *nb_joins_on_procs = (int *) malloc((n_rank) * sizeof(int));
  PDM_MPI_Allgather(&nb_joins_on_proc, 1, PDM_MPI_INT, nb_joins_on_procs, 1, PDM_MPI_INT, _multipart->comm);

  int *proc_to_joins_idx = (int *) malloc((n_rank + 1) * sizeof(int));
  proc_to_joins_idx[0] = 0;
  for (int i=0; i < n_rank; i++) {
    proc_to_joins_idx[i+1] = proc_to_joins_idx[i] + nb_joins_on_procs[i];
  }

  PDM_printf("[%d] proc_to_joins_idx are :", i_rank);
  for (int i=0; i < n_rank + 1; i ++){
    PDM_printf(" %d ", proc_to_joins_idx[i]);
  }
  PDM_printf("\n");

  int *proc_to_joins = (int *) malloc(proc_to_joins_idx[n_rank] * sizeof(int));
  PDM_MPI_Allgatherv(joins_on_proc, nb_joins_on_proc, PDM_MPI_INT,
                     proc_to_joins, nb_joins_on_procs, proc_to_joins_idx, PDM_MPI_INT,
                     _multipart->comm);

  free(joins_on_proc);
  PDM_printf("[%d] proc_to_joins is :", i_rank);
  for (int i=0; i < proc_to_joins_idx[n_rank]; i ++){
    PDM_printf(" %d ", proc_to_joins[i]);
  }
  PDM_printf("\n");

  // We have proc -> join, move back to join -> proc
  int *join_to_procs_idx = (int *) malloc((n_total_joins + 1) * sizeof(int));
  for (int i=0; i < n_total_joins + 1; i++){
    join_to_procs_idx[i] = 0;
  }
  for (int iproc = 0; iproc < n_rank; iproc++) {
    for (int ijoin = proc_to_joins_idx[iproc]; ijoin < proc_to_joins_idx[iproc+1]; ijoin++) {
      int join_id = proc_to_joins[ijoin];
      join_to_procs_idx[join_id + 1] ++;
    }
  }
  for (int i=1; i < n_total_joins + 1; i++){
    join_to_procs_idx[i] = join_to_procs_idx[i] + join_to_procs_idx[i-1];
  }

  int *join_to_procs = (int *) malloc(join_to_procs_idx[n_total_joins] * sizeof(int));
  int *ngiven_on_join = (int *) malloc(n_total_joins * sizeof(int));
  for (int i=0; i < n_total_joins; i++){
    ngiven_on_join[i] = 0;
  }
  for (int iproc = 0; iproc < n_rank; iproc++) {
    for (int ijoin = proc_to_joins_idx[iproc]; ijoin < proc_to_joins_idx[iproc+1]; ijoin++) {
      int join_id = proc_to_joins[ijoin];
      join_to_procs[join_to_procs_idx[join_id] + ngiven_on_join[join_id]] = iproc;
      ngiven_on_join[join_id] ++;
    }
  }
  free(ngiven_on_join);
  free(proc_to_joins_idx);
  free(proc_to_joins);

  PDM_printf("[%d] join_to_procs_idx are :", i_rank);
  for (int i=0; i < n_total_joins + 1; i ++) {
    PDM_printf(" %d ", join_to_procs_idx[i]);
  }
  PDM_printf("\n");

  PDM_printf("[%d] join_to_procs is :", i_rank);
  for (int i=0; i < join_to_procs_idx[n_total_joins]; i ++) {
    PDM_printf(" %d ", join_to_procs[i]);
  }
  PDM_printf("\n");

  // Step 1. Count data
  int *data_to_send_n = (int *) malloc(n_rank * sizeof(int));
  for (int i=0; i < n_rank; i++) {
    data_to_send_n[i] = 0;
  }

  for (int izone = 0; izone < _multipart->n_zone; izone++) {
    for (int i_part = 0; i_part < _multipart->n_part[izone]; i_part++) {
      int idx = bounds_and_joins_idx[izone] + i_part; //TO CHECK
      int *face_join_idx = _multipart->pbounds_and_joins[idx]->face_join_idx;
      int *joins_ids     = _multipart->pbounds_and_joins[idx]->joins_ids;
      for (int ijoin = 0; ijoin < _multipart->n_bounds_and_joins[2*izone+1]; ijoin++) {
        int n_face_on_join = face_join_idx[ijoin+1] - face_join_idx[ijoin];
        if (n_face_on_join > 0) {
          // Get destination and deduce procs that could require this data
          int join_gid = joins_ids[ijoin];
          int opp_join_gid = _multipart->join_to_opposite[join_gid];
          PDM_printf("[%i] Zone %i, i_part %i, ijoin %i (gid %i) : joinopp %i --> receiving parts are",
                     i_rank, izone, i_part, ijoin, join_gid, opp_join_gid);
          for (int i = join_to_procs_idx[opp_join_gid]; i < join_to_procs_idx[opp_join_gid+1]; i++) {
            int destProc = join_to_procs[i];

            PDM_printf(" proc %d ", destProc);
            //We have the destination, exchanged data is 5 times the lenght of point list
            // (pl value, LNToGN value, join_gid value, giving proc, giving part)
            data_to_send_n[destProc] += 5*n_face_on_join;
          }
        }
        PDM_printf("\n");

      }
    }
  }
  // Step 2. Prepare data and performs alltoall
  // Prepare stride
  int *data_to_send_idx = (int *) malloc((n_rank+1) * sizeof(int));
  data_to_send_idx[0] = 0;
  for (int i = 1; i < n_rank + 1; i++) {
    data_to_send_idx[i] = data_to_send_idx[i-1] + data_to_send_n[i-1];
    data_to_send_n[i-1] = 0;
  }
  //Prepare data
  PDM_g_num_t *data_to_send = (PDM_g_num_t *) malloc(data_to_send_idx[n_rank] * sizeof(PDM_g_num_t));
  for (int izone = 0; izone < _multipart->n_zone; izone++) {
    for (int i_part = 0; i_part < _multipart->n_part[izone]; i_part++) {
      int idx = bounds_and_joins_idx[izone] + i_part; //TO CHECK
      int *face_join_idx    = _multipart->pbounds_and_joins[idx]->face_join_idx;
      int *joins_ids     = _multipart->pbounds_and_joins[idx]->joins_ids;
      int *face_join       = _multipart->pbounds_and_joins[idx]->face_join;
      PDM_g_num_t *face_join_ln_to_gn = _multipart->pbounds_and_joins[idx]->face_join_ln_to_gn;
      for (int ijoin = 0; ijoin < _multipart->n_bounds_and_joins[2*izone+1]; ijoin++) {
        int n_face_on_join = face_join_idx[ijoin+1] - face_join_idx[ijoin];
        if (n_face_on_join > 0) {
          int join_gid = joins_ids[ijoin];
          int opp_join_gid = _multipart->join_to_opposite[join_gid];
          for (int i = join_to_procs_idx[opp_join_gid]; i < join_to_procs_idx[opp_join_gid+1]; i++) {
            int destProc = join_to_procs[i];
            int idx2 = data_to_send_idx[destProc] + data_to_send_n[destProc];
            k = 0;
            for (int iface = face_join_idx[ijoin]; iface < face_join_idx[ijoin+1]; iface++) {
              data_to_send[idx2 + 5*k    ] = face_join[4*iface];
              data_to_send[idx2 + 5*k + 1] = face_join_ln_to_gn[iface];
              data_to_send[idx2 + 5*k + 2] = opp_join_gid;
              data_to_send[idx2 + 5*k + 3] = i_rank;
              data_to_send[idx2 + 5*k + 4] = i_part;
              k += 1;
            }
            data_to_send_n[destProc] += 5*k;
          }
        }
      }
    }
  }
  //Exchange sizes
  int *data_to_recv_n   = (int *) malloc(n_rank * sizeof(int));
  PDM_MPI_Alltoall(data_to_send_n, 1, PDM_MPI_INT, data_to_recv_n, 1, PDM_MPI_INT, _multipart->comm);
  int *data_to_recv_idx = (int *) malloc((n_rank+1) * sizeof(int));
  data_to_recv_idx[0] = 0;
  for(int i = 1; i < (n_rank+1); i++) {
    data_to_recv_idx[i] = data_to_recv_idx[i-1] + data_to_recv_n[i-1];
  }
  //Exchange data
  PDM_g_num_t *data_to_recv = (PDM_g_num_t *) malloc(data_to_recv_idx[n_rank]*sizeof(PDM_g_num_t));
  PDM_MPI_Alltoallv(data_to_send,
                    data_to_send_n,
                    data_to_send_idx,
                    PDM__PDM_MPI_G_NUM,
                    data_to_recv,
                    data_to_recv_n,
                    data_to_recv_idx,
                    PDM__PDM_MPI_G_NUM,
                    _multipart->comm);

  int nRecv = data_to_recv_idx[n_rank]/5;
  free(data_to_send);
  free(data_to_send_n);
  free(data_to_send_idx);
  free(data_to_recv_n);
  free(data_to_recv_idx);

  // Step 3. Search in received data the matching faces -> To optimize
  for (int izone = 0; izone < _multipart->n_zone; izone++) {
    for (int ipart = 0; ipart < _multipart->n_part[izone]; ipart++) {
      int idx = bounds_and_joins_idx[izone] + ipart; //TO CHECK
      int *face_join_idx    = _multipart->pbounds_and_joins[idx]->face_join_idx;
      int *face_join        = _multipart->pbounds_and_joins[idx]->face_join;
      int *joins_ids     = _multipart->pbounds_and_joins[idx]->joins_ids;
      int *face_join_ln_to_gn = _multipart->pbounds_and_joins[idx]->face_join_ln_to_gn;
      // A optimiser
      for (int ijoin = 0; ijoin < _multipart->n_bounds_and_joins[2*izone + 1]; ijoin++) {
        int join_gid = joins_ids[ijoin];
        for (int i = 0; i < nRecv; i++) {
          if (data_to_recv[5*i + 2] == join_gid)
          {
            for (int j = face_join_idx[ijoin]; j < face_join_idx[ijoin+1]; j++) {
              if (data_to_recv[5*i+1] == face_join_ln_to_gn[j])
              {
                face_join[4*j+1] = data_to_recv[5*i+3]; //Proc
                face_join[4*j+2] = data_to_recv[5*i+4]; //Part
                face_join[4*j+3] = data_to_recv[5*i];
              }
            }
          }
        }
      }
    }
  }

  free(bounds_and_joins_idx);
  free(join_to_procs_idx);
  free(join_to_procs);
}

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief Build a multipart structure
 *
 * \param [in]   n_zone       Number of zones in the original mesh
 * \param [in]   n_part       Number of partition per proc in each zone
 * \param [in]   merge_blocks Merge or not the zones before splitting
 * \param [in]   split_method Choice of library used to split the mesh
 * \param [in]   comm         PDM_MPI communicator
 *
 * \return     Identifier
 */

int
PDM_multipart_create
(
 const int              n_zone,
 const int             *n_part,
 const PDM_bool_t       merge_blocks,
 const PDM_part_split_t split_method,
 const PDM_MPI_Comm     comm
)
{
  printf("PDM_multipart_create::n_zone:: %d \n", n_zone);
  printf("PDM_multipart_create::n_part:: %d \n", n_part[0]);
  printf("PDM_multipart_create::split_method:: %d \n", split_method);
  /*
   * Search a ppart free id
   */

  if (_multiparts == NULL) {
    _multiparts = PDM_Handles_create (4);
  }

  _pdm_multipart_t *_multipart = (_pdm_multipart_t *) malloc(sizeof(_pdm_multipart_t));
  int id = PDM_Handles_store (_multiparts, _multipart);

  _multipart->n_zone      = n_zone;
  _multipart->n_part      = n_part;
  _multipart->merge_blocks= merge_blocks;
  _multipart->split_method= split_method;
  _multipart->comm        = comm;

  _multipart->dmeshes_ids      = (int *) malloc(_multipart->n_zone * sizeof(int));
  _multipart->part_ids         = (int *) malloc(_multipart->n_zone * sizeof(int));
  _multipart->n_bounds_and_joins = (int *) malloc(_multipart->n_zone * 2 * sizeof(int));

  for (int izone = 0; izone < _multipart->n_zone; izone++) {
    _multipart->dmeshes_ids[izone] = -1;
    _multipart->part_ids   [izone] = -1;
    _multipart->n_bounds_and_joins[2*izone]   = -1;
    _multipart->n_bounds_and_joins[2*izone+1] = -1;
  }

  _multipart->n_total_joins = 0;
  _multipart->join_to_opposite = NULL;

  int n_rank;
  PDM_MPI_Comm_size(comm, &n_rank);

  // Number of partitions in each zone (distribution)
  _multipart->part_zone_distribution = (int *) malloc((n_zone + 1) * sizeof(int));
  int *part_zone_distribution = _multipart->part_zone_distribution;
  part_zone_distribution[0] = 0;

  // For each zone (slot of n_rank + 1 in array), number of part per proc (distribution)
  int *dpart_proc = (int *) malloc(n_zone*(n_rank + 1) * sizeof(int));
  for (int izone = 0; izone < _multipart->n_zone; izone++) {
    dpart_proc[izone*(n_rank + 1)] = 0;
    PDM_MPI_Allgather((void *) &n_part[izone],
                      1,
                      PDM_MPI_INT,
                      (void *) (&dpart_proc[izone*(n_rank+1) + 1]),
                      1,
                      PDM_MPI_INT,
                      _multipart->comm);

    for (int i = 1; i < n_rank+1; i++) {
      dpart_proc[izone*(n_rank+1) + i] = dpart_proc[izone*(n_rank+1) + i] + dpart_proc[izone*(n_rank+1) + i-1];
    }

    part_zone_distribution[izone+1] = part_zone_distribution[izone] + dpart_proc[izone*(n_rank+1) + n_rank];
  }
  // For each global part number, owner proc and i_part in proc
  _multipart->gpart_to_proc = (int *) malloc(2*part_zone_distribution[n_zone] * sizeof(int));
  for (int izone = 0; izone < _multipart->n_zone; izone++){
    int zshift = part_zone_distribution[izone];
    for (int i = 0; i < n_rank; i++) {
      for (int j = dpart_proc[izone*(n_rank+1) + i]; j < dpart_proc[izone*(n_rank+1) + i+1]; j++) {
        _multipart->gpart_to_proc[2*(zshift + j)] = i;
        _multipart->gpart_to_proc[2*(zshift + j) + 1] = j - dpart_proc[izone*(n_rank+1) + i];

      }
    }
  }
  free(dpart_proc);

  PDM_printf("Created from PDM_multipart_create. You requested a multipart with %d zones \n", n_zone);
  return id;

}

/* TODO : copy doc of the function */
void PDM_multipart_register_block
(
 const int        mpart_id,
 const int        zone_gid,
 const int        block_data_id
)
{
  PDM_printf("In multipart %d, set zone n°%d using blockdata %d \n",
             mpart_id, zone_gid, block_data_id);

  _pdm_multipart_t *_multipart = _get_from_id (mpart_id);
  assert(zone_gid < _multipart->n_zone);
  _multipart->dmeshes_ids[zone_gid] = block_data_id;
}

void PDM_multipart_register_joins
(
 const int        mpart_id,
 const int        n_total_joins,
 int             *matching_join_array
)
{
  _pdm_multipart_t *_multipart = _get_from_id (mpart_id);

  _multipart->n_total_joins    = n_total_joins;
  _multipart->join_to_opposite = matching_join_array;
}

void
PDM_multipart_run_ppart
(
 const int id
)
{
  _pdm_multipart_t *_multipart = _get_from_id (id);

  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(_multipart->comm, &i_rank);
  PDM_MPI_Comm_size(_multipart->comm, &n_rank);

  if (_multipart->merge_blocks)
  {
    // 1. Generate global numerotation using all blocks
    // 2. Call the partitionner once on the global numbering
  }
  else
  {
    // 2. Loop over the blocks and call the partitionner
    for (int zone_gid = 0; zone_gid < _multipart->n_zone; zone_gid++) {
      PDM_printf("You requested no merge : partitionning zone %d/%d \n", zone_gid+1, _multipart->n_zone);
      int block_id = _multipart->dmeshes_ids[zone_gid];
      PDM_printf("block id for zone %d is %d\n", zone_gid, block_id);
      int dn_cell  = 0;
      int dn_face  = 0;
      int dn_vtx   = 0;
      int n_bnd    = 0;
      int n_join   = 0;
      const double       *dvtx_coord;
      const int          *dface_vtx_idx;
      const PDM_g_num_t  *dface_vtx;
      const PDM_g_num_t  *dface_cell;
      const int          *dface_bound_idx;
      const PDM_g_num_t  *dface_bound;
      const int          *djoin_gids;
      const int          *dface_join_idx;
      const PDM_g_num_t  *dface_join;

      int n_face_group = 0;
      int          *dface_group_idx = NULL;
      PDM_g_num_t  *dface_group    = NULL;
      int          *dface_tag      = NULL;

      PDM_dmesh_dims_get(block_id, &dn_cell, &dn_face, &dn_vtx, &n_bnd, &n_join);
      PDM_dmesh_data_get(block_id, &dvtx_coord, &dface_vtx_idx, &dface_vtx, &dface_cell,
                         &dface_bound_idx, &dface_bound, &djoin_gids, &dface_join_idx, &dface_join);

      //Merge face_bounds and face_joins into face_group
      n_face_group = n_bnd + n_join;
      dface_group_idx = (int *)         malloc((n_face_group + 1) * sizeof(int));
      dface_group     = (PDM_g_num_t *) malloc((dface_bound_idx[n_bnd] + dface_join_idx[n_join]) * sizeof(PDM_g_num_t));

      for (int i=0; i < n_bnd + 1; i++)
        dface_group_idx[i] = dface_bound_idx[i];
      for (int i=0; i < dface_bound_idx[n_bnd]; i++)
        dface_group[i] = dface_bound[i];

      for (int i=1; i < n_join + 1; i++)
        dface_group_idx[n_bnd + i] = dface_bound_idx[n_bnd] + dface_join_idx[i];
      for (int i=0; i < dface_join_idx[n_join]; i++)
        dface_group[dface_bound_idx[n_bnd] + i] = dface_join[i];

      //Store number of bounds and joins in the structure
      _multipart->n_bounds_and_joins[2*zone_gid    ] = n_bnd;
      _multipart->n_bounds_and_joins[2*zone_gid + 1] = n_join;

      dface_tag = (int *) malloc((dn_face) * sizeof(int));
      _set_dface_tag_from_joins(dn_face, n_join, djoin_gids, dface_join_idx, dface_join, dface_tag, _multipart->comm);

      int ppart_id = 0;
      int have_dcell_part = 0;
      int *dcell_part = (int *) malloc(dn_cell*sizeof(int));
      // We probably should create a subcomm to imply only the procs sharing this block
      PDM_part_create(&ppart_id,
                      _multipart->comm,
                      _multipart->split_method,
                      "PDM_PART_RENUM_CELL_NONE",
                      "PDM_PART_RENUM_FACE_NONE",
                      0,                          // n_property_cell
                      NULL,                       // renum_properties_cell
                      0,                          // n_property_face
                      NULL,                       // renum_properties_face
                      _multipart->n_part[zone_gid],
                      dn_cell,
                      dn_face,
                      dn_vtx,
                      n_face_group,
                      NULL,                       // dcell_faceIdx
                      NULL,                       // dcell_face
                      NULL,                       // dcell_tag
                      NULL,                       // dcell_weight
                      have_dcell_part,
                      dcell_part,                  // dcell_part
                      dface_cell,
                      dface_vtx_idx,
                      dface_vtx,
                      dface_tag,
                      dvtx_coord,
                      NULL,                       // dvtx_tag
                      dface_group_idx,
                      dface_group);
      PDM_printf("Partitionning done, ppardId is %d \n", ppart_id);
      //Store the partition id for future access
      _multipart->part_ids[zone_gid] = ppart_id;

      free(dcell_part);
      free(dface_group_idx);
      free(dface_group);
      free(dface_tag);
    }
    // Now separate joins and boundaries and we rebuild joins over the zones
    _split_bounds_and_joins(_multipart);
    _search_matching_joins(_multipart);
  }
}

void
PDM_multipart_part_dim_get
(
const   int  mpart_id,
const   int  zone_gid,
const   int  i_part,
 int        *n_cell,
 int        *n_face,
 int        *n_face_part_bound,
 int        *n_vtx,
 int        *n_proc,
 int        *n_total_part,
 int        *scell_face,
 int        *sface_vtx,
 int        *sface_bound,
 int        *n_face_bound,
 int        *sface_join,
 int        *n_face_join
)
{
  _pdm_multipart_t *_multipart = _get_from_id (mpart_id);

  assert(zone_gid < _multipart->n_zone && i_part < _multipart->n_part[zone_gid]);
  int ppart_id = _multipart->part_ids[zone_gid];

  PDM_part_part_dim_get(ppart_id,
                        i_part,
                        n_cell,
                        n_face,
                        n_face_part_bound,
                        n_vtx,
                        n_proc,
                        n_total_part,
                        scell_face,
                        sface_vtx,
                        sface_bound,
                        n_face_bound);

  // Get boundary and join data from pbounds_and_joins
  int idx = 0;
  for (int i = 0; i < zone_gid; i++)
    idx += _multipart->n_part[i];
  idx += i_part;

  //Attention au cas ou pas de face de bord
  *n_face_bound = _multipart->pbounds_and_joins[idx]->n_bound;
  *sface_bound = _multipart->pbounds_and_joins[idx]->face_bound_idx[*n_face_bound];
  *n_face_join  = _multipart->pbounds_and_joins[idx]->n_join;
  *sface_join  = _multipart->pbounds_and_joins[idx]->face_join_idx[*n_face_join];

}

void
PDM_multipart_part_val_get
(
const int            mpart_id,
const int            zone_gid,
const int            i_part,
      int          **cell_tag,
      int          **cell_face_idx,
      int          **cell_face,
      PDM_g_num_t  **cell_ln_to_gn,
      int          **face_tag,
      int          **face_cell,
      int          **face_vtx_idx,
      int          **face_vtx,
      PDM_g_num_t  **face_ln_to_gn,
      int          **face_part_bound_proc_idx,
      int          **face_part_bound_part_idx,
      int          **face_part_bound,
      int          **vtx_tag,
      double       **vtx,
      PDM_g_num_t  **vtx_ln_to_gn,
      int          **face_bound_idx,
      int          **face_bound,
      PDM_g_num_t  **face_bound_ln_to_gn,
      int          **face_join_idx,
      int          **face_join,
      PDM_g_num_t  **face_join_ln_to_gn
)
{
   _pdm_multipart_t *_multipart = _get_from_id (mpart_id);

  assert(zone_gid < _multipart->n_zone && i_part < _multipart->n_part[zone_gid]);
  int ppart_id = _multipart->part_ids[zone_gid];

  PDM_part_part_val_get(ppart_id,
                        i_part,
                        cell_tag,
                        cell_face_idx,
                        cell_face,
                        cell_ln_to_gn,
                        face_tag,
                        face_cell,
                        face_vtx_idx,
                        face_vtx,
                        face_ln_to_gn,
                        face_part_bound_proc_idx,
                        face_part_bound_part_idx,
                        face_part_bound,
                        vtx_tag,
                        vtx,
                        vtx_ln_to_gn,
                        face_bound_idx,
                        face_bound,
                        face_bound_ln_to_gn
                        );

  // Get boundary and join data from pbounds_and_joins
  int idx = 0;
  for (int i = 0; i < zone_gid; i++)
    idx += _multipart->n_part[i];
  idx += i_part;
  //Attention au cas ou pas de face de bord
  *face_bound_idx       = _multipart->pbounds_and_joins[idx]->face_bound_idx;
  *face_bound           = _multipart->pbounds_and_joins[idx]->face_bound;
  *face_bound_ln_to_gn  = _multipart->pbounds_and_joins[idx]->face_bound_ln_to_gn;
  *face_join_idx        = _multipart->pbounds_and_joins[idx]->face_join_idx;
  *face_join            = _multipart->pbounds_and_joins[idx]->face_join;
  *face_join_ln_to_gn   = _multipart->pbounds_and_joins[idx]->face_join_ln_to_gn;

}

void
PDM_multipart_part_color_get
(
const int            mpart_id,
const int            zone_gid,
const int            i_part,
      int          **cell_color,
      int          **face_color,
      int          **thread_color,
      int          **hyperplane_color
)
{
  _pdm_multipart_t *_multipart = _get_from_id (mpart_id);

  assert(zone_gid < _multipart->n_zone && i_part < _multipart->n_part[zone_gid]);
  int ppart_id = _multipart->part_ids[zone_gid];

  PDM_part_part_color_get(ppart_id,
                          i_part,
                          cell_color,
                          face_color,
                          thread_color,
                          hyperplane_color);

}

void
PDM_multipart_time_get
(
const int       mpart_id,
const int       zone_gid,
      double  **elapsed,
      double  **cpu,
      double  **cpu_user,
      double  **cpu_sys
)
{
  _pdm_multipart_t *_multipart = _get_from_id (mpart_id);
  assert(zone_gid < _multipart->n_zone);
  int ppart_id = _multipart->part_ids[zone_gid];

  PDM_part_time_get(ppart_id,
                    elapsed,
                    cpu,
                    cpu_user,
                    cpu_sys);

}

/**
 *
 * \brief Free
 *
 * \param [in]   id           Identifier
 *
 */

void
PDM_multipart_free
(
 const int id
)
{
  _pdm_multipart_t *_multipart = _get_from_id (id);

  free(_multipart->dmeshes_ids);
  free(_multipart->n_bounds_and_joins);
  free(_multipart->part_zone_distribution);
  free(_multipart->gpart_to_proc);

  for (int izone = 0; izone<_multipart->n_zone; izone++)
    PDM_part_free(_multipart->part_ids[izone]);
  free(_multipart->part_ids);

  free (_multipart);

  PDM_Handles_handle_free (_multiparts, id, PDM_FALSE);

  const int n_multipart = PDM_Handles_n_get (_multiparts);

  if (n_multipart == 0) {
    _multiparts = PDM_Handles_free (_multiparts);
  }
  PDM_printf("Cleaned from PDM_multipart_free\n");
}


/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
