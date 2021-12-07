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
#include "pdm_distrib.h"
#include "pdm_timer.h"
#include "pdm_partitioning_algorithm.h"
#include "pdm_para_graph_dual.h"
#include "pdm_part_renum.h"
#include "pdm_handles.h"
#include "pdm_dmesh.h"
#include "pdm_dmesh_nodal.h"
#include "pdm_mesh_nodal.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_binary_search.h"
#include "pdm_part_to_block.h"
#include "pdm_block_to_part.h"
#include "pdm_multi_block_to_part.h"
#include "pdm_distrib.h"
#include "pdm_array.h"
#include "pdm_logging.h"
#include "pdm_para_graph_dual.h"
#include "pdm_dconnectivity_transform.h"
#include "pdm_unique.h"
#include "pdm_partitioning_nodal_algorithm.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_multipart_priv.h"

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

/**
 *
 * \brief Map each pair of (join, opposite join) to a same global id and count
 *        the total number of faces in this unified join. Return a distribution.
 *        Arrays are allocated in this function.
 *
 * \param [in]   _multipart          multipart object
 * \param [out]  join_to_ref_join    Unique join id associated to each join
 *                                     (size = n_total_join)
 * \param [out]  face_in_join_distri Distribution of join faces over the ref
 *                                   join ids (size = n_unique_joins+1)
 */
static void
_build_join_uface_distribution
(
 _pdm_multipart_t  *_multipart,
 int              **join_to_ref_join,
 int              **face_in_join_distri
)
{

  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(_multipart->comm, &i_rank);
  PDM_MPI_Comm_size(_multipart->comm, &n_rank);

  // PDM_printf("pdm::_build_join_uface_distribution\n");
  int n_total_joins  = _multipart->n_total_joins;
  int n_unique_joins = n_total_joins/2;
  *join_to_ref_join    = (int *) malloc(n_total_joins  * sizeof(int));
  *face_in_join_distri = (int *) malloc((n_unique_joins+1) * sizeof(int));
  int* _face_in_join_distri = *face_in_join_distri;
  int* _join_to_ref_join    = *join_to_ref_join;

  //Build join_to_ref_join : we want the join and opposite join to have the same shift index,
  // so we take the smaller join global id as the reference
  int ref_join_gid = 0;
  for (int ijoin = 0; ijoin < n_total_joins; ijoin++)
  {
    int opp_join = _multipart->join_to_opposite[ijoin];
    if (ijoin < opp_join)
    {
      _join_to_ref_join[ijoin] = ref_join_gid;
      _join_to_ref_join[opp_join] = ref_join_gid;
      ref_join_gid ++;
    }
  }
  /*
  PDM_printf("Join to reference join :");
  for (int ijoin = 0; ijoin < n_total_joins; ijoin++)
   PDM_printf(" %d ", _join_to_ref_join[ijoin]);
  PDM_printf("\n");
  */

  //Count faces in joins
  int *nb_face_in_joins = PDM_array_zeros_int(n_unique_joins);

  for (int izone = 0; izone < _multipart->n_zone; izone++){
    for (int i_part = 0; i_part < _multipart->n_part[izone]; i_part++){
      int *pface_join_idx = _multipart->pmeshes[izone].parts[i_part]->face_join_idx;
      for (int ijoin=0; ijoin < _multipart->pmeshes[izone].n_joins; ijoin ++){
        int join_gid     = _multipart->pmeshes[izone].joins_ids[ijoin];
        int join_opp_gid = _multipart->join_to_opposite[join_gid];
        //Paired joins must be counted only once
        if (join_gid < join_opp_gid)
          nb_face_in_joins[_join_to_ref_join[join_gid]] += pface_join_idx[ijoin+1] - pface_join_idx[ijoin];
      }
    }
  }
  /*
  PDM_printf("[%d] nb_face_joins : ", i_rank);
  for (int i = 0; i < n_unique_joins ; i++)
    PDM_printf(" %d ", nb_face_in_joins[i]);
  PDM_printf("\n");
  */

  //Sum faces and build distribution
  PDM_MPI_Allreduce(nb_face_in_joins, &_face_in_join_distri[1], n_unique_joins,
                    PDM_MPI_INT, PDM_MPI_SUM, _multipart->comm);

  _face_in_join_distri[0] = 0;
  PDM_array_accumulate_int(_face_in_join_distri, n_unique_joins+1);

  /*
  PDM_printf("[%d] _face_in_join_distri : ", i_rank);
  for (int i = 0; i < n_unique_joins + 1; i++)
    PDM_printf(" %d ", _face_in_join_distri[i]);
  PDM_printf("\n");
  */

  free(nb_face_in_joins);
  // PDM_printf("pdm::_build_join_uface_distribution end \n");
}

/**
 *
 * \brief Complete join data, which originally only contains face local id,
 *        with the connecting data opp proc, opp part, opp face local id.
 *
 * \param [inout]   _multipart          multipart object
 */
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

  //Construction of (unique) join distribution
  int *join_to_ref_join;
  int *face_in_join_distri;
  _build_join_uface_distribution(_multipart, &join_to_ref_join, &face_in_join_distri);

  //Count total nb of join_faces
  int nb_of_joins = 0;
  for (int izone = 0 ; izone < _multipart->n_zone; izone ++) {
    nb_of_joins += _multipart->n_part[izone] * _multipart->pmeshes[izone].n_joins;
  }

  // Prepare lntogn numbering and partitioned data
  PDM_g_num_t **shifted_lntogn = (PDM_g_num_t **) malloc(nb_of_joins * sizeof(PDM_g_num_t*));
  int              **part_data = (int **)         malloc(nb_of_joins * sizeof(int *));
  int        *nb_face_per_join = (int *)          malloc(nb_of_joins * sizeof(int));

  int ijoin_pos  = 0;
  for (int izone = 0 ; izone < _multipart->n_zone; izone ++) {
    int n_join            = _multipart->pmeshes[izone].n_joins;
    _part_mesh_t _pmeshes = _multipart->pmeshes[izone];
    for (int i_part = 0; i_part < _multipart->n_part[izone]; i_part++) {
      int         *face_join_idx    = _pmeshes.parts[i_part]->face_join_idx;
      int         *face_join        = _pmeshes.parts[i_part]->face_join;
      PDM_g_num_t *face_join_lntogn = _pmeshes.parts[i_part]->face_join_ln_to_gn;
      for (int ijoin = 0; ijoin < n_join; ijoin++) {
        int join_size = face_join_idx[ijoin + 1] - face_join_idx[ijoin];
        nb_face_per_join[ijoin_pos] = join_size;
        PDM_g_num_t *shifted_lntogn_loc = (PDM_g_num_t *) malloc(join_size * sizeof(PDM_g_num_t));
        int         *part_data_loc      = (int *)         malloc(3 * join_size * sizeof(int));
        //Get shift value from join unique distribution
        int join_gid    = _multipart->pmeshes[izone].joins_ids[ijoin];
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
        shifted_lntogn[ijoin_pos] = shifted_lntogn_loc;
        part_data[ijoin_pos]      = part_data_loc;
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

  int         *block_data;
  PDM_part_to_block_exch(ptb,
                         sizeof(int),
                         PDM_STRIDE_CST,
                         3,
                         NULL,
                         (void **) part_data,
                         NULL,
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
      printf(" (%d %d %d) ", block_data[3*i],
                             block_data[3*i+1],
                             block_data[3*i+2]);
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
    new_part_data[ijoin] = (int *) malloc(6*nb_face_per_join[ijoin]*sizeof(int));
  }
  int cst_stride = 6;

  PDM_block_to_part_exch(btp,
                         sizeof(int),
                         PDM_STRIDE_CST,
                         &cst_stride,
                         (void *) block_data,
                         NULL,
                         (void **) new_part_data);

  free(block_data);
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
        PDM_printf(" (%d %d %d %d %d %d) ", new_part_data[ijoin][6*iface],
                                            new_part_data[ijoin][6*iface+1],
                                            new_part_data[ijoin][6*iface+2],
                                            new_part_data[ijoin][6*iface+3],
                                            new_part_data[ijoin][6*iface+4],
                                            new_part_data[ijoin][6*iface+5]);
      PDM_printf("\n");
    }
  }
  */


  //Process received data
  ijoin_pos = 0;
  for (int izone = 0 ; izone < _multipart->n_zone; izone ++) {
    int n_join = _multipart->pmeshes[izone].n_joins;
    _part_mesh_t _pmeshes = _multipart->pmeshes[izone];
    for (int i_part = 0; i_part < _multipart->n_part[izone]; i_part++) {
      int *face_join_idx = _pmeshes.parts[i_part]->face_join_idx;
      int *face_join     = _pmeshes.parts[i_part]->face_join;
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
              opp_part = part_data_loc[6*i + 5];
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
  free(join_to_ref_join);
  free(face_in_join_distri);
  free(nb_face_per_join);
}

/**
 *
 * \brief Free the memory occuped by a partition structure
 *
 * \param [inout]   part          _part_t object
 */
static void
_part_free
(
 _part_t         *part,
 PDM_ownership_t  owner
)
{
  if(owner == PDM_OWNERSHIP_KEEP){
    if (part->vtx != NULL)
      free(part->vtx);
    part->vtx = NULL;

    if (part->face_vtx_idx != NULL)
      free(part->face_vtx_idx);
    part->face_vtx_idx = NULL;

    if (part->face_vtx != NULL)
      free(part->face_vtx);
    part->face_vtx = NULL;

    if (part->gface_vtx != NULL)
      free(part->gface_vtx);
    part->gface_vtx = NULL;

    if (part->cell_face_idx != NULL)
      free(part->cell_face_idx);
    part->cell_face_idx = NULL;

    if (part->cell_face != NULL)
      free(part->cell_face);
    part->cell_face = NULL;

    if (part->gcell_face != NULL)
      free(part->gcell_face);
    part->gcell_face = NULL;

    if (part->face_cell != NULL)
      free(part->face_cell);
    part->face_cell = NULL;

    if (part->face_group_idx != NULL)
      free(part->face_group_idx);
    part->face_group_idx = NULL;

    if (part->face_group != NULL)
      free(part->face_group);
    part->face_group = NULL;

    if (part->face_part_bound_proc_idx != NULL)
      free(part->face_part_bound_proc_idx);
    part->face_part_bound_proc_idx = NULL;

    if (part->face_part_bound_part_idx != NULL)
      free(part->face_part_bound_part_idx);
    part->face_part_bound_part_idx = NULL;

    if (part->face_part_bound != NULL)
      free(part->face_part_bound);
    part->face_part_bound = NULL;

    if (part->vtx_part_bound_proc_idx != NULL)
      free(part->vtx_part_bound_proc_idx);
    part->vtx_part_bound_proc_idx = NULL;

    if (part->vtx_part_bound_part_idx != NULL)
      free(part->vtx_part_bound_part_idx);
    part->vtx_part_bound_part_idx = NULL;

    if (part->vtx_part_bound != NULL)
      free(part->vtx_part_bound);
    part->vtx_part_bound = NULL;

    if(part->vtx_ghost_information != NULL)
      free(part->vtx_ghost_information);
    part->vtx_ghost_information = NULL;

    if (part->face_bound_idx != NULL)
      free(part->face_bound_idx);
    part->face_bound_idx = NULL;

    if (part->face_bound != NULL)
      free(part->face_bound);
    part->face_bound = NULL;

    if (part->face_join_idx != NULL)
      free(part->face_join_idx);
    part->face_join_idx = NULL;

    if (part->face_join != NULL)
      free(part->face_join);
    part->face_join = NULL;

    if (part->vtx_ln_to_gn != NULL)
      free(part->vtx_ln_to_gn);
    part->vtx_ln_to_gn = NULL;

    if (part->face_ln_to_gn != NULL)
      free(part->face_ln_to_gn);
    part->face_ln_to_gn = NULL;

    if (part->cell_ln_to_gn != NULL)
      free(part->cell_ln_to_gn);
    part->cell_ln_to_gn = NULL;

    if (part->face_group_ln_to_gn != NULL)
      free(part->face_group_ln_to_gn);
    part->face_group_ln_to_gn = NULL;

    if (part->face_bound_ln_to_gn != NULL)
      free(part->face_bound_ln_to_gn);
    part->face_bound_ln_to_gn = NULL;

    if (part->face_join_ln_to_gn != NULL)
      free(part->face_join_ln_to_gn);
    part->face_join_ln_to_gn = NULL;

    if (part->cell_tag != NULL)
      free(part->cell_tag);
    part->cell_tag = NULL;

    if (part->face_tag != NULL)
    free(part->face_tag);
    part->face_tag = NULL;

    if (part->vtx_tag != NULL)
      free(part->vtx_tag);
    part->vtx_tag = NULL;

    if (part->cell_color != NULL)
      free(part->cell_color);
    part->cell_color = NULL;

    if (part->face_color != NULL)
      free(part->face_color);
    part->face_color = NULL;

    if (part->thread_color != NULL)
      free(part->thread_color);
    part->thread_color = NULL;

    if (part->hyperplane_color != NULL)
      free(part->hyperplane_color);
    part->hyperplane_color = NULL;

  } /* End owner */

  /* Following is not results but internal array */
  if (part->new_to_old_order_cell != NULL)
    free(part->new_to_old_order_cell);
  part->new_to_old_order_cell = NULL;

  if (part->new_to_old_order_face != NULL)
    free(part->new_to_old_order_face);
  part->new_to_old_order_face = NULL;

  if (part->new_to_old_order_vtx != NULL)
    free(part->new_to_old_order_vtx);
  part->new_to_old_order_vtx = NULL;

  if(part->subpartlayout != NULL){
    if(part->subpartlayout->cell_tile_idx!= NULL)
      free(part->subpartlayout->cell_tile_idx);
    if(part->subpartlayout->face_tile_idx!= NULL)
      free(part->subpartlayout->face_tile_idx);
    if(part->subpartlayout->face_bnd_tile_idx!= NULL)
      free(part->subpartlayout->face_bnd_tile_idx);
    if(part->subpartlayout->mask_tile_idx!= NULL)
      free(part->subpartlayout->mask_tile_idx);
    if(part->subpartlayout->cell_vect_tile_idx!= NULL)
      free(part->subpartlayout->cell_vect_tile_idx);
    if(part->subpartlayout->mask_tile_n!= NULL)
      free(part->subpartlayout->mask_tile_n);
    if(part->subpartlayout->cell_vect_tile_n!= NULL)
      free(part->subpartlayout->cell_vect_tile_n);
    if(part->subpartlayout->mask_tile!= NULL)
      free(part->subpartlayout->mask_tile);
    free(part->subpartlayout);
  }


  free(part);
}


// static inline
// int
// _vtx_is_in_connectivity
// (
//   PDM_g_num_t  vtx,
//   PDM_g_num_t *first_vtx,
//   int          n_vtx
// )
// {
//   for (int i=0; i<n_vtx; ++i) {
//     if (first_vtx[i]==vtx) return 1;
//   }
//   return 0;
// }

// static
// int
// _is_parent
// (
//   PDM_g_num_t *first_vtx,
//   int          n_vtx,
//   PDM_g_num_t *first_face_vtx,
//   int          n_face_vtx
// )
// {
//   // WARNING: quadratic but n_face_vtx should be 3 or 4
//   for (int i=0; i<n_face_vtx; ++i) {
//     if (!(_vtx_is_in_connectivity(first_face_vtx[i],first_vtx,n_vtx))) {
//       return 0;
//     }
//   }
//   return 1;
// }

static void
_setup_ghost_information
(
const int   i_rank,
const int   n_part,
const int  *pn_vtx,
      int **pinternal_vtx_priority
)
{
  /* 0 : Interior / 1 : owner join (at least one) / 2 : not owner */
  for (int ipart = 0; ipart < n_part; ipart++) {
    for(int ivtx = 0; ivtx < pn_vtx[ipart]; ++ivtx) {
      if(pinternal_vtx_priority[ipart][ivtx] == i_rank){
        pinternal_vtx_priority[ipart][ivtx] = 1;
      } else if(pinternal_vtx_priority[ipart][ivtx] == -1) {
        pinternal_vtx_priority[ipart][ivtx] = 0;
      } else { /* Pas owner / pas interieur donc sur 1 autre proc */
         pinternal_vtx_priority[ipart][ivtx] = 2;
      }
    }
  }
}

static void
_compute_part_mesh_nodal_3d
(
 PDM_dmesh_t       *dmesh,
 PDM_dmesh_nodal_t *dmn,
 _part_mesh_t      *pm,
 int                n_part
)
{
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(dmn->comm, &i_rank);
  PDM_MPI_Comm_size(dmn->comm, &n_rank);

  /*
   * Rebuild the volumic part from cell
   */
  PDM_g_num_t  **pcell_ln_to_gn = (PDM_g_num_t ** ) malloc( n_part * sizeof(PDM_g_num_t *));
  PDM_g_num_t  **pvtx_ln_to_gn  = (PDM_g_num_t ** ) malloc( n_part * sizeof(PDM_g_num_t *));
  int           *pn_cell        = (int *  )         malloc( n_part * sizeof(int          ));
  int           *pn_vtx         = (int *  )         malloc( n_part * sizeof(int          ));
  for(int i_part = 0; i_part < n_part; ++i_part){
    pcell_ln_to_gn[i_part] = pm->parts[i_part]->cell_ln_to_gn;
    pn_cell       [i_part] = pm->parts[i_part]->n_cell;

    pvtx_ln_to_gn[i_part] = pm->parts[i_part]->vtx_ln_to_gn;
    pn_vtx       [i_part] = pm->parts[i_part]->n_vtx;
  }

  PDM_part_mesh_nodal_elmts_t* pmn_vol = PDM_dmesh_nodal_elmts_to_part_mesh_nodal_elmts(dmn->volumic,
                                                                                        n_part,
                                                                                        pn_vtx,
                                                                                        pvtx_ln_to_gn,
                                                                                        pn_cell,
                                                                                        pcell_ln_to_gn);

  PDM_g_num_t  **pface_ln_to_gn = (PDM_g_num_t ** ) malloc( n_part * sizeof(PDM_g_num_t *));
  int           *pn_face        = (int *  )         malloc( n_part * sizeof(int          ));
  for(int i_part = 0; i_part < n_part; ++i_part) {
    pface_ln_to_gn[i_part] = pm->parts[i_part]->face_ln_to_gn;
    pn_face       [i_part] = pm->parts[i_part]->n_face;
  }
  int          *pn_surf;
  PDM_g_num_t **psurf_gnum;
  PDM_g_num_t **psurf_to_face_g_num;
  PDM_reverse_dparent_gnum(dmn->surfacic->dparent_gnum,
                           NULL, // dparent_sign
                           dmesh->face_distrib,
                           dmn->surfacic->delmt_child_distrib,
                           n_part,
                           pn_face,
                           pface_ln_to_gn,
                          &pn_surf,
                          &psurf_gnum,
                          &psurf_to_face_g_num,
                           NULL, // pchild_parent_sign
                           dmn->comm);

  PDM_part_mesh_nodal_elmts_t* pmn_surf = PDM_dmesh_nodal_elmts_to_part_mesh_nodal_elmts(dmn->surfacic,
                                                                                         n_part,
                                                                                         pn_vtx,
                                                                                         pvtx_ln_to_gn,
                                                                                         pn_surf,
                                                                                         psurf_gnum);

  for(int i_part = 0; i_part < n_part; ++i_part){
    free(psurf_gnum[i_part]);
    free(psurf_to_face_g_num[i_part]);
  }
  free(pn_surf);
  free(psurf_to_face_g_num);

  /* Create top structure */
  PDM_part_mesh_nodal_t* pmn = PDM_part_mesh_nodal_create(dmn->mesh_dimension,
                                                          n_part,
                                                          dmn->comm);

  PDM_part_mesh_nodal_add_part_mesh_nodal_elmts(pmn, pmn_vol , PDM_OWNERSHIP_KEEP);
  PDM_part_mesh_nodal_add_part_mesh_nodal_elmts(pmn, pmn_surf, PDM_OWNERSHIP_KEEP);
  // PDM_part_mesh_nodal_add_part_mesh_nodal_elmts(pmn, pmn_ridge, PDM_OWNERSHIP_KEEP);
  PDM_part_mesh_nodal_free(pmn);

  free(pface_ln_to_gn);
  free(pn_face);
  free(pcell_ln_to_gn);
  free(pn_cell);
  free(pvtx_ln_to_gn);
  free(pn_vtx);
}

// LET IT COMMENT BUT REPLACE BY GENERIC PDM_reverse_dparent_gnum
/* Rebuild surface from volumic */
// int dn_surf_elmt = dmn->surfacic->delmt_child_distrib[i_rank+1] - dmn->surfacic->delmt_child_distrib[i_rank];
// PDM_part_to_block_t* ptb = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
//                                                     PDM_PART_TO_BLOCK_POST_MERGE,
//                                                     1.,
//                                                     &dmn->surfacic->dparent_gnum,
//                                                     NULL,
//                                                     &dn_surf_elmt,
//                                                     1,
//                                                     dmn->comm);

// int         *pblk_surf_n    = (int         *) malloc( dn_surf_elmt * sizeof(int        ));
// PDM_g_num_t *pblk_surf_gnum = (PDM_g_num_t *) malloc( dn_surf_elmt * sizeof(PDM_g_num_t));
// for(int i = 0; i < dn_surf_elmt; ++i) {
//   pblk_surf_n   [i] = 1;
//   pblk_surf_gnum[i] = dmn->surfacic->delmt_child_distrib[i_rank] + 1; // Donc le gnum de surf ...
// }

// int         *blk_surf_n = NULL;
// PDM_g_num_t *blk_surf_gnum   = NULL;
// PDM_part_to_block_exch(ptb,
//                        sizeof(PDM_g_num_t),
//                        PDM_STRIDE_VAR,
//                        -1,
//                       &pblk_surf_n,
//            (void **)  &pblk_surf_gnum,
//                       &blk_surf_n,
//            (void **)  &blk_surf_gnum);

// PDM_g_num_t n_g_face = dmesh->face_distrib[n_rank]+1;
// PDM_g_num_t* block_distrib_tmp_idx = PDM_part_to_block_adapt_partial_block_to_block(ptb, &blk_surf_n, n_g_face);

// PDM_part_to_block_free(ptb);
// free(pblk_surf_n   );
// free(pblk_surf_gnum);

// PDM_block_to_part_t* btp = PDM_block_to_part_create(block_distrib_tmp_idx,
//                             (const PDM_g_num_t **)  pface_ln_to_gn,
//                                                     pn_face,
//                                                     n_part,
//                                                     dmn->comm);

// int         **psurf_n    = NULL;
// PDM_g_num_t **psurf_gnum = NULL;
// PDM_block_to_part_exch2(btp,
//                         sizeof(PDM_g_num_t),
//                         PDM_STRIDE_VAR,
//                         blk_surf_n,
//                         blk_surf_gnum,
//                        &psurf_n,
//             (void ***) &psurf_gnum);

// PDM_block_to_part_free(btp);

// /*
//  * At this stage we have for each partition the number AND the gnum of surfs inside
//  *          -> We sort psurf_gnum but normally it's unique
//  *
//  */
// int* pn_surf = (int *) malloc(n_part * sizeof(int));
// for(int i_part = 0; i_part < n_part; ++i_part) {

//   int pn_surf_tmp = 0;
//   for(int i_face = 0; i_face < pn_face[i_part]; ++i_face) {
//     pn_surf_tmp += psurf_n[i_part][i_face];
//     assert(psurf_n[i_part][i_face] <= 1); // DOnc soit 0 soit 1
//   }
//   pn_surf   [i_part] = PDM_inplace_unique_long(psurf_gnum[i_part], NULL, 0, pn_surf_tmp-1);
//   psurf_gnum[i_part] = realloc(psurf_gnum[i_part], pn_surf[i_part] * sizeof(PDM_g_num_t));
// }


static void
_compute_part_mesh_nodal_2d
(
 PDM_dmesh_t       *dmesh,
 PDM_dmesh_nodal_t *dmn,
 _part_mesh_t      *pm,
 int                n_part
)
{
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(dmn->comm, &i_rank);
  PDM_MPI_Comm_size(dmn->comm, &n_rank);

  /*
   * Rebuild the volumic part from cell
   */
  PDM_g_num_t  **pface_ln_to_gn = (PDM_g_num_t ** ) malloc( n_part * sizeof(PDM_g_num_t *));
  int           *pn_face        = (int *  )         malloc( n_part * sizeof(int          ));

  PDM_g_num_t  **pvtx_ln_to_gn  = (PDM_g_num_t ** ) malloc( n_part * sizeof(PDM_g_num_t *));
  int           *pn_vtx         = (int *  )         malloc( n_part * sizeof(int          ));
  for(int i_part = 0; i_part < n_part; ++i_part){
    pface_ln_to_gn[i_part] = pm->parts[i_part]->face_ln_to_gn;
    pn_face       [i_part] = pm->parts[i_part]->n_face;

    pvtx_ln_to_gn[i_part] = pm->parts[i_part]->vtx_ln_to_gn;
    pn_vtx       [i_part] = pm->parts[i_part]->n_vtx;
  }

  PDM_part_mesh_nodal_elmts_t* pmn_surf = PDM_dmesh_nodal_elmts_to_part_mesh_nodal_elmts(dmn->surfacic,
                                                                                         n_part,
                                                                                         pn_vtx,
                                                                                         pvtx_ln_to_gn,
                                                                                         pn_face,
                                                                                         pface_ln_to_gn);

  PDM_g_num_t  **pedge_ln_to_gn = (PDM_g_num_t ** ) malloc( n_part * sizeof(PDM_g_num_t *));
  int           *pn_edge        = (int *  )         malloc( n_part * sizeof(int          ));
  for(int i_part = 0; i_part < n_part; ++i_part) {
    pedge_ln_to_gn[i_part] = pm->parts[i_part]->edge_ln_to_gn;
    pn_edge       [i_part] = pm->parts[i_part]->n_edge;
  }
  int          *pn_surf;
  PDM_g_num_t **psurf_gnum;
  PDM_g_num_t **psurf_to_face_g_num;
  PDM_reverse_dparent_gnum(dmn->surfacic->dparent_gnum,
                           NULL, // dparent_sign
                           dmesh->edge_distrib,
                           dmn->surfacic->delmt_child_distrib,
                           n_part,
                           pn_edge,
                           pedge_ln_to_gn,
                          &pn_surf,
                          &psurf_gnum,
                          &psurf_to_face_g_num,
                           NULL, // pchild_parent_sign
                           dmn->comm);

  PDM_part_mesh_nodal_elmts_t* pmn_ridge = PDM_dmesh_nodal_elmts_to_part_mesh_nodal_elmts(dmn->ridge,
                                                                                          n_part,
                                                                                          pn_vtx,
                                                                                          pvtx_ln_to_gn,
                                                                                          pn_surf,
                                                                                          psurf_gnum);

  for(int i_part = 0; i_part < n_part; ++i_part){
    free(psurf_gnum[i_part]);
    free(psurf_to_face_g_num[i_part]);
  }
  free(pn_surf);
  free(psurf_to_face_g_num);

  /* Create top structure */
  PDM_part_mesh_nodal_t* pmn = PDM_part_mesh_nodal_create(dmn->mesh_dimension,
                                                          n_part,
                                                          dmn->comm);

  PDM_part_mesh_nodal_add_part_mesh_nodal_elmts(pmn, pmn_surf , PDM_OWNERSHIP_KEEP);
  PDM_part_mesh_nodal_add_part_mesh_nodal_elmts(pmn, pmn_ridge, PDM_OWNERSHIP_KEEP);
  PDM_part_mesh_nodal_free(pmn);

  free(pedge_ln_to_gn);
  free(pn_edge);
  free(pface_ln_to_gn);
  free(pn_face);
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
 * \param [in]   part_size_method Choice of homogeneous or heterogeneous partitions
 * \param [in]   part_weight  Weight (in %) of each partition in heterogeneous case
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
 const PDM_split_dual_t split_method,
 const PDM_part_size_t  part_size_method,
 const double          *part_fraction,
 const PDM_MPI_Comm     comm,
 const PDM_ownership_t  owner
)
{
  // printf("PDM_multipart_create::n_zone:: %d \n", n_zone);
  // printf("PDM_multipart_create::n_part:: %d \n", n_part[0]);
  // printf("PDM_multipart_create::split_method:: %d \n", split_method);

  /*
   * Search a ppart free id
   */

  if (_multiparts == NULL) {
    _multiparts = PDM_Handles_create (4);
  }

  _pdm_multipart_t *_multipart = (_pdm_multipart_t *) malloc(sizeof(_pdm_multipart_t));
  int id = PDM_Handles_store (_multiparts, _multipart);

  _multipart->n_zone           = n_zone;
  _multipart->n_part           = n_part;
  _multipart->merge_blocks     = merge_blocks;
  _multipart->split_method     = split_method;
  _multipart->part_size_method = part_size_method;
  _multipart->part_fraction    = part_fraction;
  _multipart->comm             = comm;
  _multipart->owner            = owner;

  _multipart->n_total_joins    = 0;
  _multipart->join_to_opposite = NULL;

  // _multipart->dmeshes_ids = (int *) malloc(_multipart->n_zone * sizeof(int));

  _multipart->dmeshes       = (PDM_dmesh_t       **) malloc(_multipart->n_zone * sizeof(PDM_dmesh_t       *));
  _multipart->dmeshes_nodal = (PDM_dmesh_nodal_t **) malloc(_multipart->n_zone * sizeof(PDM_dmesh_nodal_t *));
  for (int i=0; i<_multipart->n_zone; ++i) {
    _multipart->dmeshes_nodal[i] = NULL;
    _multipart->dmeshes[i]       = NULL;
  }

  _multipart->pmeshes       = (_part_mesh_t *) malloc(_multipart->n_zone * sizeof(_part_mesh_t ));

  int _renum_cell_method = PDM_part_renum_method_cell_idx_get("PDM_PART_RENUM_CELL_NONE");
  int _renum_face_method = PDM_part_renum_method_face_idx_get("PDM_PART_RENUM_FACE_NONE");
  for (int izone = 0; izone < _multipart->n_zone; izone++) {
    _multipart->pmeshes[izone].renum_cell_method = _renum_cell_method;
    _multipart->pmeshes[izone].renum_face_method = _renum_face_method;
    _multipart->pmeshes[izone].renum_cell_properties = NULL;
    _multipart->pmeshes[izone].joins_ids = NULL;
    _multipart->pmeshes[izone].parts     = NULL;
  }

  return id;
}

/**
 *
 * \brief Set distributed mesh data for the input zone
 *
 * \param [in]   mpart_id       Multipart structure id
 * \param [in]   zone_id        Global zone id
 * \param [in]   dmesh_id       Id of the distributed mesh structure to use
 */
void PDM_multipart_register_block
(
 const int          mpart_id,
 const int          zone_id,
       PDM_dmesh_t *dmesh
)
{
  _pdm_multipart_t *_multipart = _get_from_id (mpart_id);
  assert(zone_id < _multipart->n_zone);
  _multipart->dmeshes[zone_id] = dmesh;
}

/**
 *
 * \brief Set distributed mesh data for the input zone
 *
 * \param [in]   mpart_id       Multipart structure id
 * \param [in]   zone_id        Global zone id
 * \param [in]   dmesh_id       Id of the distributed mesh structure to use
 */
void PDM_multipart_register_dmesh_nodal
(
 const int                mpart_id,
 const int                zone_id,
       PDM_dmesh_nodal_t *dmesh_nodal
)
{
  _pdm_multipart_t *_multipart = _get_from_id (mpart_id);
  assert(zone_id < _multipart->n_zone);
  assert(_multipart->dmeshes_nodal[zone_id] == NULL);
  _multipart->dmeshes_nodal[zone_id] = dmesh_nodal;
}

/**
 *
 * \brief Set connecting data between all the zones
 *
 * \param [in]   mpart_id          Multipart structure id
 * \param [in]   n_total_joins     Total number of interfaces
 * \param [in]   join_to_opposite  For each global join id, give the global id
 *                                   of the opposite join (size = n_total_joins)
 *
 * \note Join global id numbering must start at 0 and be continuous.
 */
void PDM_multipart_register_joins
(
 const int        mpart_id,
 const int        n_total_joins,
 const int       *join_to_opposite
)
{
  _pdm_multipart_t *_multipart = _get_from_id (mpart_id);

  _multipart->n_total_joins    = n_total_joins;
  _multipart->join_to_opposite = join_to_opposite;
}

/**
 *
 * \brief Set the reordering methods to be used after partitioning
 *
 * \param [in]   mpart_id           Multipart structure id
 * \param [in]   i_zone             Id of zone which parameters apply (or -1 for all zones)
 * \param [in]   renum_cell_method  Choice of renumbering method for cells
 * \param [in]   renum_cell_properties Parameters used by cacheblocking method :
 *                                     [n_cell_per_cache_wanted, is_asynchrone, is_vectorisation,
                                        n_vect_face, split_method]
 * \param [in]   renum_face_method  Choice of renumbering method for faces
 *
 */
void PDM_multipart_set_reordering_options
(
 const int        mpart_id,
 const int        i_zone,
 const char      *renum_cell_method,
 const int       *renum_cell_properties,
 const char      *renum_face_method
)
{
  _pdm_multipart_t *_multipart = _get_from_id (mpart_id);

  int _renum_cell_method = PDM_part_renum_method_cell_idx_get(renum_cell_method);
  int _renum_face_method = PDM_part_renum_method_face_idx_get(renum_face_method);
  if (_renum_cell_method == -1) {
    PDM_error (__FILE__, __LINE__, 0, "'%s' is an unknown renumbering cell method\n", renum_cell_method);
  }
  if (_renum_face_method == -1) {
    PDM_error (__FILE__, __LINE__, 0, "'%s' is an unknown renumbering face method\n", renum_face_method);
  }

  if(i_zone < 0) {
    for (int izone = 0; izone < _multipart->n_zone; izone++) {
      _multipart->pmeshes[izone].renum_cell_method = _renum_cell_method;
      _multipart->pmeshes[izone].renum_face_method = _renum_face_method;
      _multipart->pmeshes[izone].renum_cell_properties = renum_cell_properties;
    }
  }
  else {
    assert(i_zone < _multipart->n_zone);
    _multipart->pmeshes[i_zone].renum_cell_method = _renum_cell_method;
    _multipart->pmeshes[i_zone].renum_face_method = _renum_face_method;
    _multipart->pmeshes[i_zone].renum_cell_properties = renum_cell_properties;
  }
}


static
void
_run_ppart_zone2
(
PDM_dmesh_t       *dmesh,
PDM_dmesh_nodal_t *dmesh_nodal,
_part_mesh_t      *pmeshes,
int                n_part,
PDM_split_dual_t   split_method,
PDM_part_size_t    part_size_method,
const double*      part_fraction,
PDM_MPI_Comm       comm
)
{
  PDM_UNUSED(dmesh_nodal);
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  // Get distributed mesh data for this zone
  int dn_cell, dn_face, dn_vtx, dn_edge, n_bnd, n_join;
  const double       *dvtx_coord;
  const int          *dface_vtx_idx;
  const PDM_g_num_t  *dface_vtx;
  PDM_g_num_t        *dface_cell;
  int                *dface_cell_idx;
  PDM_g_num_t        *dface_cell_tmp;

  PDM_g_num_t        *dual_graph_idx;
  int                *dcell_face_idx;
  PDM_g_num_t        *dual_graph, *dcell_face;
  const int          *dface_bound_idx;
  const PDM_g_num_t  *dface_bound;
  const int          *joins_ids;
  const int          *dface_join_idx;
  const PDM_g_num_t  *dface_join;
  PDM_dmesh_dims_get(dmesh, &dn_cell, &dn_face, &dn_edge, &dn_vtx, &n_bnd, &n_join);

  printf(" dn_cell = %i \n", dn_cell);
  printf(" dn_face = %i \n", dn_face);
  printf(" dn_edge = %i \n", dn_edge);
  printf(" dn_vtx  = %i \n", dn_vtx );
  printf(" n_bnd   = %i \n", n_bnd  );
  if(dmesh->dn_cell == 0){

    PDM_dmesh_connectivity_get(dmesh, PDM_CONNECTIVITY_TYPE_EDGE_FACE,
                               &dface_cell_tmp,
                               &dface_cell_idx,
                               PDM_OWNERSHIP_KEEP);

    PDM_dmesh_connectivity_get(dmesh, PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                               &dcell_face,
                               &dcell_face_idx,
                               PDM_OWNERSHIP_KEEP);


  } else {
    abort();
    PDM_dmesh_data_get(dmesh, &dvtx_coord, &dface_vtx_idx, &dface_vtx, (const PDM_g_num_t**)&dface_cell,
                       &dface_bound_idx, &dface_bound, &joins_ids, &dface_join_idx, &dface_join);

  }

  // This will store all the partitions created by this proc on this zone
  // Copy number of bounds and joins (global data) in the part structure
  // printf("pmeshes->n_bounds :: %i \n", n_bnd);
  // pmeshes->n_bounds  = n_bnd;
  // pmeshes->n_joins   = n_join;
  // pmeshes->joins_ids = (int *) malloc(n_join * sizeof(int));
  // for (int i_join = 0; i_join < n_join; i_join++){
  //   pmeshes->joins_ids[i_join] = joins_ids[i_join];
  // }

  // Compute total number of partitions for this zone
  int tn_part;
  PDM_MPI_Allreduce(&n_part, &tn_part, 1, PDM_MPI_INT, PDM_MPI_SUM, comm);
  pmeshes->tn_part = tn_part;

  // PDM_log_trace_array_int(dcell_face_idx, dn_face+1, "dcell_face_idx :: ");
  // PDM_log_trace_array_long(dcell_face, dcell_face_idx[dn_face], "dcell_face :: ");

  // PDM_log_trace_array_int(dface_cell_idx, dn_edge+1, "dface_cell_idx :: ");
  // PDM_log_trace_array_long(dface_cell, dface_cell_idx[dn_edge], "dface_cell :: ");

  // PDM_log_trace_array_long(dface_cell, 2*dn_edge, "dface_cell :: ");

  //Construct graph and split -- no interface for now, do it by hand
  PDM_g_num_t *cell_distri = NULL;
  PDM_g_num_t *face_distri = NULL;
  PDM_g_num_t *vtx_distri  = NULL;
  PDM_g_num_t *part_distri = NULL;
  if(dmesh->dn_cell != 0) {
    cell_distri = PDM_compute_entity_distribution(comm, dn_cell);
    face_distri = PDM_compute_entity_distribution(comm, dn_face);
    assert(dface_cell_idx == NULL);
    PDM_setup_connectivity_idx(dn_face,
                               2,
                               dface_cell_tmp,
                               &dface_cell_idx,
                               &dface_cell);
    // if(dface_cell_idx == NULL) {
    //   dface_cell_idx = (int *) malloc( sizeof(int) * (dn_edge+1));
    //   dface_cell_idx[0] = 0;
    //   for(int i = 0; i < dn_edge; ++i){
    //     dface_cell_idx[i+1] = dface_cell_idx[i] + 2;
    //   }
    // }
  } else {
    cell_distri = PDM_compute_entity_distribution(comm, dn_face);
    face_distri = PDM_compute_entity_distribution(comm, dn_edge);
    assert(dface_cell_idx == NULL);
    PDM_setup_connectivity_idx(dn_edge,
                               2,
                               dface_cell_tmp,
                               &dface_cell_idx,
                               &dface_cell);
    // if(dface_cell_idx == NULL) {
    //   dface_cell_idx = (int *) malloc( sizeof(int) * (dn_edge+1));
    //   dface_cell_idx[0] = 0;
    //   for(int i = 0; i < dn_edge; ++i){
    //     dface_cell_idx[i+1] = dface_cell_idx[i] + 2;
    //   }
    // }
  }
  vtx_distri = PDM_compute_entity_distribution(comm, dn_vtx );
  part_distri = PDM_compute_entity_distribution(comm, n_part );

  PDM_UNUSED(vtx_distri);

  PDM_log_trace_array_int (dface_cell_idx, dn_edge+1, "dface_cell_idx :: ");
  PDM_log_trace_array_long(dface_cell    , dface_cell_idx[dn_edge], "dface_cell :: ");

  PDM_deduce_combine_connectivity_dual(comm,
                                       cell_distri,
                                       face_distri,
                                       dcell_face_idx,
                                       dcell_face,
                                       dface_cell_idx,
                                       dface_cell,
                                       1, // is signed
                                       &dual_graph_idx,
                                       &dual_graph);
  free(dface_cell_idx);
  free(dface_cell);

  PDM_log_trace_array_long(dual_graph_idx, dn_face+1, "dual_graph_idx :: ");
  PDM_log_trace_array_long(dual_graph, dual_graph_idx[dn_face], "dual_graph :: ");

  double *part_fractions = NULL;
  if (part_size_method == PDM_PART_SIZE_HETEROGENEOUS){
    int *n_part_per_rank = (int    *) malloc( n_rank * sizeof(int   ));
    int *displ           = (int    *) malloc( n_rank * sizeof(int   ));
    part_fractions       = (double *) malloc(tn_part * sizeof(double));
    for (int i =0; i < n_rank; i++){
      n_part_per_rank[i] = part_distri[i+1] - part_distri[i];
      displ[i] = part_distri[i];
    }

    PDM_MPI_Allgatherv((void*) part_fraction,
                       n_part,
                       PDM_MPI_DOUBLE,
                       part_fractions,
                       n_part_per_rank,
                       displ,
                       PDM_MPI_DOUBLE,
                       comm);
    free(n_part_per_rank);
    free(displ);
  }
  int *cell_part = NULL;
  if(dmesh->dn_cell != 0) {
    cell_part = (int *) malloc(dn_cell * sizeof(int));
  } else {
    cell_part = (int *) malloc(dn_face * sizeof(int));
  }
  PDM_para_graph_split (split_method,
                        cell_distri,
                        dual_graph_idx,
                        dual_graph,
                        NULL, NULL,
                        tn_part,
                        part_fractions,
                        cell_part,
                        comm);

  // PDM_log_trace_array_int (cell_part, dn_face, "cell_part :: ");


  free(dual_graph_idx);
  free(dual_graph);
  if (part_size_method == PDM_PART_SIZE_HETEROGENEOUS) {
    free(part_fractions);
  }

  // Partitioning algorithm except to work on 2d arrays (external size = n_part) whereas
  // _part_t struture stores n_part * 1d arrays so we have to use tmp pointers
  // int          *pn_vtx                        = NULL;
  int          *pn_cell                       = NULL;
  // int          *pn_face                       = NULL;
  PDM_g_num_t **pcell_ln_to_gn                = NULL;
  // PDM_g_num_t **pface_ln_to_gn                = NULL;

  PDM_part_assemble_partitions(comm,
                               part_distri,
                               cell_distri,
                               cell_part,
                              &pn_cell,
                              &pcell_ln_to_gn);
  free(cell_part);

  for(int i_part = 0; i_part < n_part; ++i_part) {
    PDM_log_trace_array_long(pcell_ln_to_gn[i_part], pn_cell[i_part], "pcell_ln_to_gn :: ");
  }


  // Fill _part_t structures with temporary arrays
  pmeshes->parts = (_part_t **) malloc(pmeshes->tn_part*sizeof(_part_t*));
  for (int ipart = 0; ipart < n_part; ipart++) {
    pmeshes->parts[ipart] = _part_create();
    pmeshes->parts[ipart]->n_cell = pn_cell[ipart];
    pmeshes->parts[ipart]->cell_ln_to_gn = pcell_ln_to_gn[ipart];
  }

  free(pn_cell);
  free(pcell_ln_to_gn);

  free(vtx_distri);
  free(part_distri);
  free(cell_distri);
  free(face_distri);
}

static
void
_run_ppart_zone
(
PDM_dmesh_t      *dmesh,
_part_mesh_t     *pmeshes,
int               n_part,
PDM_split_dual_t  split_method,
PDM_part_size_t   part_size_method,
const double*     part_fraction,
PDM_MPI_Comm      comm
)
{
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  // Get distributed mesh data for this zone
  int dn_cell, dn_face, dn_vtx, dn_edge, n_bnd, n_join;
  const double       *dvtx_coord;
  const int          *dface_vtx_idx;
  const PDM_g_num_t  *dface_vtx;
  const PDM_g_num_t  *dface_cell;
  const int          *dface_bound_idx;
  const PDM_g_num_t  *dface_bound;
  const int          *joins_ids;
  const int          *dface_join_idx;
  const PDM_g_num_t  *dface_join;
  PDM_dmesh_dims_get(dmesh, &dn_cell, &dn_face, &dn_edge, &dn_vtx, &n_bnd, &n_join);
  PDM_dmesh_data_get(dmesh, &dvtx_coord, &dface_vtx_idx, &dface_vtx, &dface_cell,
                     &dface_bound_idx, &dface_bound, &joins_ids, &dface_join_idx, &dface_join);

  printf(" dn_cell = %i \n", dn_cell);
  printf(" dn_face = %i \n", dn_face);
  printf(" dn_edge = %i \n", dn_edge);
  printf(" dn_vtx  = %i \n", dn_vtx );
  printf(" n_bnd   = %i \n", n_bnd  );

  // This will store all the partitions created by this proc on this zone
  // Copy number of bounds and joins (global data) in the part structure
  // printf("pmeshes->n_bounds :: %i \n", n_bnd);
  pmeshes->n_bounds  = n_bnd;
  pmeshes->n_joins   = n_join;
  pmeshes->joins_ids = (int *) malloc(n_join * sizeof(int));
  for (int i_join = 0; i_join < n_join; i_join++){
    pmeshes->joins_ids[i_join] = joins_ids[i_join];
  }

  // Compute total number of partitions for this zone
  int tn_part;
  PDM_MPI_Allreduce(&n_part, &tn_part, 1, PDM_MPI_INT, PDM_MPI_SUM, comm);
  pmeshes->tn_part = tn_part;

  //Construct graph and split -- no interface for now, do it by hand
  PDM_g_num_t *cell_distri = PDM_compute_entity_distribution(comm, dn_cell);
  PDM_g_num_t *face_distri = PDM_compute_entity_distribution(comm, dn_face);
  PDM_g_num_t *vtx_distri  = PDM_compute_entity_distribution(comm, dn_vtx );
  PDM_g_num_t *part_distri = PDM_compute_entity_distribution(comm, n_part );

  PDM_g_num_t *dual_graph_idx;
  int         *dcell_face_idx;
  PDM_g_num_t *dual_graph, *dcell_face;

  PDM_para_graph_dual_from_arc2node(comm,
                                    cell_distri,
                                    face_distri,
                                    dface_cell,
                                   &dual_graph_idx,
                                   &dual_graph,
                                    1,
                                   &dcell_face_idx,
                                   &dcell_face);

  // free(dual_graph_idx);
  // free(dual_graph);
  // PDM_para_graph_dual_from_combine_connectivity(comm,
  //                                               cell_distri,
  //                                               face_distri,
  //                                               vtx_distri,
  //                                               dcell_face_idx,
  //                                               dcell_face,
  //                                               dface_vtx_idx,
  //                                               dface_vtx,
  //                              (PDM_g_num_t**) &dual_graph_idx,
  //                              (PDM_g_num_t**) &dual_graph);

  double *part_fractions = NULL;
  if (part_size_method == PDM_PART_SIZE_HETEROGENEOUS){
    int *n_part_per_rank = (int    *) malloc( n_rank * sizeof(int   ));
    int *displ           = (int    *) malloc( n_rank * sizeof(int   ));
    part_fractions       = (double *) malloc(tn_part * sizeof(double));
    for (int i =0; i < n_rank; i++){
      n_part_per_rank[i] = part_distri[i+1] - part_distri[i];
      displ[i] = part_distri[i];
    }

    PDM_MPI_Allgatherv((void*) part_fraction,
                       n_part,
                       PDM_MPI_DOUBLE,
                       part_fractions,
                       n_part_per_rank,
                       displ,
                       PDM_MPI_DOUBLE,
                       comm);
    free(n_part_per_rank);
    free(displ);
  }
  int *cell_part = (int *) malloc(dn_cell * sizeof(int));
  PDM_para_graph_split (split_method,
                        cell_distri,
                        dual_graph_idx,
                        dual_graph,
                        NULL, NULL,
                        tn_part,
                        part_fractions,
                        cell_part,
                        comm);

  // free(dual_graph_idx);
  // free(dual_graph);
  if (part_size_method == PDM_PART_SIZE_HETEROGENEOUS) {
    free(part_fractions);
  }

  // Partitioning algorithm except to work on 2d arrays (external size = n_part) whereas
  // _part_t struture stores n_part * 1d arrays so we have to use tmp pointers
  int          *pn_vtx                        = NULL;
  int          *pn_cell                       = NULL;
  int          *pn_face                       = NULL;
  double      **pvtx_coord                    = NULL;
  int         **pface_vtx_idx                 = NULL;
  int         **pface_vtx                     = NULL;
  int         **pcell_face_idx                = NULL;
  int         **pcell_face                    = NULL;
  int         **pface_cell                    = NULL;
  int         **pface_bound_idx               = NULL;
  int         **pface_bound                   = NULL;
  int         **pface_join_idx                = NULL;
  int         **pinternal_face_bound_proc_idx = NULL;
  int         **pinternal_face_bound_part_idx = NULL;
  int         **pinternal_face_bound          = NULL;
  int         **pinternal_vtx_bound_proc_idx  = NULL;
  int         **pinternal_vtx_bound_part_idx  = NULL;
  int         **pinternal_vtx_bound           = NULL;
  int         **pinternal_vtx_priority        = NULL;
  PDM_g_num_t **pcell_ln_to_gn                = NULL;
  PDM_g_num_t **pface_ln_to_gn                = NULL;
  PDM_g_num_t **pvtx_ln_to_gn                 = NULL;
  PDM_g_num_t **pface_bound_ln_to_gn          = NULL;
  PDM_g_num_t **pface_join_ln_to_gn           = NULL;

  PDM_part_assemble_partitions(comm,
                               part_distri,
                               cell_distri,
                               cell_part,
                              &pn_cell,
                              &pcell_ln_to_gn);

  free(dual_graph_idx);
  free(dual_graph);
  free(cell_part);

  PDM_part_dconnectivity_to_pconnectivity_sort(comm,
                                               cell_distri,
                                               dcell_face_idx,
                                               dcell_face,
                                               n_part,
                                               pn_cell,
                        (const PDM_g_num_t **) pcell_ln_to_gn,
                                              &pn_face,
                                              &pface_ln_to_gn,
                                              &pcell_face_idx,
                                              &pcell_face);

  PDM_part_dconnectivity_to_pconnectivity_sort(comm,
                                               face_distri,
                                               dface_vtx_idx,
                                               dface_vtx,
                                               n_part,
                                               pn_face,
                        (const PDM_g_num_t **) pface_ln_to_gn,
                                              &pn_vtx,
                                              &pvtx_ln_to_gn,
                                              &pface_vtx_idx,
                                              &pface_vtx);
  free(dcell_face_idx);
  free(dcell_face);

  PDM_part_reverse_pcellface(n_part, pn_cell, pn_face,
                             (const int **) pcell_face_idx, (const int **) pcell_face,
                            &pface_cell);

  PDM_part_reorient_bound_faces(n_part,
                                pn_face,
                                pface_cell,
                 (const int **) pcell_face_idx,
                                pcell_face,
                 (const int **) pface_vtx_idx,
                                pface_vtx,
                                NULL,  // pface_edge_idx
                                NULL); // pface_edge

  PDM_part_dcoordinates_to_pcoordinates(comm,
                                        n_part,
                                        vtx_distri,
                                        dvtx_coord,
                                        pn_vtx,
                 (const PDM_g_num_t **) pvtx_ln_to_gn,
                                       &pvtx_coord);

  // Fill _part_t structures with temporary arrays
  pmeshes->parts = (_part_t **) malloc(pmeshes->tn_part*sizeof(_part_t*));
  for (int ipart = 0; ipart < n_part; ipart++) {
    pmeshes->parts[ipart] = _part_create();

    pmeshes->parts[ipart]->n_cell = pn_cell[ipart];
    pmeshes->parts[ipart]->n_face = pn_face[ipart];
    pmeshes->parts[ipart]->n_vtx  = pn_vtx[ipart];

    pmeshes->parts[ipart]->vtx           = pvtx_coord[ipart];
    pmeshes->parts[ipart]->face_vtx_idx  = pface_vtx_idx[ipart];
    pmeshes->parts[ipart]->face_vtx      = pface_vtx[ipart];
    pmeshes->parts[ipart]->cell_face_idx = pcell_face_idx[ipart];
    pmeshes->parts[ipart]->cell_face     = pcell_face[ipart];
    pmeshes->parts[ipart]->face_cell     = pface_cell[ipart];

    pmeshes->parts[ipart]->cell_ln_to_gn = pcell_ln_to_gn[ipart];
    pmeshes->parts[ipart]->face_ln_to_gn = pface_ln_to_gn[ipart];
    pmeshes->parts[ipart]->vtx_ln_to_gn  = pvtx_ln_to_gn[ipart];
  }

  // Use the structure to call renumbering procedures
  PDM_part_renum_cell(pmeshes->parts, n_part, pmeshes->renum_cell_method,
                      (void *) pmeshes->renum_cell_properties);
  PDM_part_renum_face(pmeshes->parts, n_part, pmeshes->renum_face_method, NULL);
  // PDM_part_renum_edge(pmeshes->parts, n_part, pmeshes->renum_face_method, NULL);
  // PDM_part_renum_vtx(pmeshes->parts, n_part, pmeshes->renum_face_method, NULL);

  /* Let's call graph communication to do the reordering - Temporary  */
  PDM_part_generate_entity_graph_comm(comm,
                                      part_distri,
                                      vtx_distri,
                                      n_part,
                                      pn_vtx,
               (const PDM_g_num_t **) pvtx_ln_to_gn,
                                      NULL,
                                     &pinternal_vtx_bound_proc_idx,
                                     &pinternal_vtx_bound_part_idx,
                                     &pinternal_vtx_bound,
                                     &pinternal_vtx_priority);

  _setup_ghost_information(i_rank, n_part, pn_vtx, pinternal_vtx_priority);
  PDM_part_renum_vtx(pmeshes->parts, n_part, 1, (void *) pinternal_vtx_priority);

  /* Free in order to be correclty */
  for (int ipart = 0; ipart < n_part; ipart++) {
    free(pinternal_vtx_bound_proc_idx[ipart]);
    free(pinternal_vtx_bound_part_idx[ipart]);
    free(pinternal_vtx_bound[ipart]);
    free(pinternal_vtx_priority[ipart]);
  }
  free(pinternal_vtx_bound_proc_idx);
  free(pinternal_vtx_bound_part_idx);
  free(pinternal_vtx_bound);
  free(pinternal_vtx_priority);


  // Now genererate bounds and comm data -- we need update pface_ln_to_gn which has been modified
  for (int ipart = 0; ipart < n_part; ipart++) {
    pface_ln_to_gn[ipart] = pmeshes->parts[ipart]->face_ln_to_gn;
  }

  PDM_part_distgroup_to_partgroup(comm,
                                  face_distri,
                                  n_bnd,
                                  dface_bound_idx,
                                  dface_bound,
                                  n_part,
                                  pn_face,
           (const PDM_g_num_t **) pface_ln_to_gn,
                                 &pface_bound_idx,
                                 &pface_bound,
                                 &pface_bound_ln_to_gn);
  int **pface_join_tmp = NULL;
  PDM_part_distgroup_to_partgroup(comm,
                                  face_distri,
                                  n_join,
                                  dface_join_idx,
                                  dface_join,
                                  n_part,
                                  pn_face,
           (const PDM_g_num_t **) pface_ln_to_gn,
                                 &pface_join_idx,
                                 &pface_join_tmp,
                                 &pface_join_ln_to_gn);

  PDM_part_generate_entity_graph_comm(comm,
                                      part_distri,
                                      face_distri,
                                      n_part,
                                      pn_face,
               (const PDM_g_num_t **) pface_ln_to_gn,
                                      NULL,
                                     &pinternal_face_bound_proc_idx,
                                     &pinternal_face_bound_part_idx,
                                     &pinternal_face_bound,
                                      NULL);

  PDM_part_generate_entity_graph_comm(comm,
                                      part_distri,
                                      vtx_distri,
                                      n_part,
                                      pn_vtx,
               (const PDM_g_num_t **) pvtx_ln_to_gn,
                                      NULL,
                                     &pinternal_vtx_bound_proc_idx,
                                     &pinternal_vtx_bound_part_idx,
                                     &pinternal_vtx_bound,
                                     &pinternal_vtx_priority); // Egalemet possible de permeuter dans PMD_part_renum

  // Re setup the array properlly to have the good output
  _setup_ghost_information(i_rank, n_part, pn_vtx, pinternal_vtx_priority);

  // Finally complete parts structure with internal join data and bounds
  for (int ipart = 0; ipart < n_part; ipart++) {
    pmeshes->parts[ipart]->face_bound_idx      = pface_bound_idx[ipart];
    pmeshes->parts[ipart]->face_bound          = pface_bound[ipart];
    pmeshes->parts[ipart]->face_bound_ln_to_gn = pface_bound_ln_to_gn[ipart];

    /* For face_join, the function only returns local id of face in join, we have to
       allocate to set up expected size (4*nb_face_join) */
    pmeshes->parts[ipart]->face_join_idx = pface_join_idx[ipart];
    int s_face_join = pmeshes->parts[ipart]->face_join_idx[n_join];
    pmeshes->parts[ipart]->face_join = (int *) malloc( 4 * s_face_join * sizeof(int));
    for (int i_face = 0; i_face < s_face_join; i_face++){
      pmeshes->parts[ipart]->face_join[4*i_face] = pface_join_tmp[ipart][i_face];
    }
    pmeshes->parts[ipart]->face_join_ln_to_gn = pface_join_ln_to_gn[ipart];
    free(pface_join_tmp[ipart]);

    pmeshes->parts[ipart]->face_part_bound_proc_idx = pinternal_face_bound_proc_idx[ipart];
    pmeshes->parts[ipart]->face_part_bound_part_idx = pinternal_face_bound_part_idx[ipart];
    pmeshes->parts[ipart]->face_part_bound          = pinternal_face_bound[ipart];

    pmeshes->parts[ipart]->vtx_part_bound_proc_idx = pinternal_vtx_bound_proc_idx[ipart];
    pmeshes->parts[ipart]->vtx_part_bound_part_idx = pinternal_vtx_bound_part_idx[ipart];
    pmeshes->parts[ipart]->vtx_part_bound          = pinternal_vtx_bound[ipart];

    pmeshes->parts[ipart]->vtx_ghost_information   = pinternal_vtx_priority[ipart];
  }

  // Free temporary arrays (top-level)
  free(pn_vtx);
  free(pn_cell);
  free(pn_face);
  free(pvtx_coord);
  free(pface_vtx_idx);
  free(pface_vtx);
  free(pcell_face_idx);
  free(pcell_face);
  free(pface_cell);
  free(pface_bound_idx);
  free(pface_bound);
  free(pface_join_idx);
  free(pinternal_face_bound_proc_idx);
  free(pinternal_face_bound_part_idx);
  free(pinternal_face_bound);
  free(pinternal_vtx_bound_proc_idx);
  free(pinternal_vtx_bound_part_idx);
  free(pinternal_vtx_bound);
  free(pcell_ln_to_gn);
  free(pface_ln_to_gn);
  free(pvtx_ln_to_gn);
  free(pface_bound_ln_to_gn);
  free(pface_join_ln_to_gn);
  free(pinternal_vtx_priority);
  free(pface_join_tmp);
  free(cell_distri);
  free(face_distri);
  free(vtx_distri);
  free(part_distri);
}

/*
 * \brief give the elt_part of a 2d element as the one of its parent
 *
 * \param [in]  n_elt    number of elements of the searched 2d element section
 * \param [in]  dn_elt   number of elements in all sections
 * \param [in]  elt_part partition number for each element
 * \param [out] elt2d_part partition number for face of the section
 */
// static
// void
// face_part_from_parent
// (
//   int           n_elt,
//   int           dn_elt,
//   int           section_idx,
//   PDM_g_num_t  *elt_dist,
//   int          *elt_part,
//   int          *delt_vtx_idx,
//   PDM_g_num_t  *delt_vtx,
//   PDM_g_num_t  *elt_elt_idx,
//   PDM_g_num_t  *elt_elt,
//   int          *elt2d_part,
//   PDM_MPI_Comm  comm
// )
// {
//   int i_rank;
//   int n_rank;
//   PDM_MPI_Comm_rank(comm, &i_rank);
//   PDM_MPI_Comm_size(comm, &n_rank);

//   PDM_g_num_t* face_neighbors_idx = (PDM_g_num_t*) malloc ((n_elt+1) * sizeof(PDM_g_num_t));
//   int pos = 0;
//   for (int i=0; i<n_elt; ++i) {
//     int elt_idx = section_idx + i;
//     int n_neighbor = elt_elt_idx[elt_idx+1]-elt_elt_idx[elt_idx];
//     face_neighbors_idx[i] = pos;
//     pos += n_neighbor;
//   }
//   face_neighbors_idx[n_elt] = pos;

//   int n_face_neighbor = face_neighbors_idx[n_elt];
//   PDM_g_num_t* face_neighbors = (PDM_g_num_t*) malloc (n_face_neighbor * sizeof(PDM_g_num_t));
//   for (int i=0; i<n_elt; ++i) {
//     int elt_idx = section_idx + i;
//     int n_neighbor = elt_elt_idx[elt_idx+1]-elt_elt_idx[elt_idx];
//     for (int j=0; j<n_neighbor; ++j) {
//       face_neighbors[face_neighbors_idx[i]+j] = elt_elt[elt_elt_idx[elt_idx]+j]+1; // +1 because block_to_part uses 1-indexed ln_to_gn
//     }
//   }

//   PDM_block_to_part_t* btp = PDM_block_to_part_create(elt_dist,(const PDM_g_num_t**)&face_neighbors,&n_face_neighbor,1,comm);

//   int stride_one = 1;
//   int* block_data1 = elt_part;
//   int** neighbor_part;
//   PDM_block_to_part_exch2(btp,sizeof(int),PDM_STRIDE_CST,&stride_one,block_data1,NULL,(void***)&neighbor_part);

//   int* block_idx2 = malloc(dn_elt * sizeof(int));
//   for (int i=0; i<dn_elt; ++i) {
//     block_idx2[i] = delt_vtx_idx[i+1] - delt_vtx_idx[i];
//   }
//   PDM_g_num_t* block_data2 = delt_vtx;
//   int** neighbor_vtx_stri;
//   PDM_g_num_t** neighbor_vtx;
//   PDM_block_to_part_exch2(btp,sizeof(PDM_g_num_t),PDM_STRIDE_VAR,block_idx2,block_data2,&neighbor_vtx_stri,(void***)&neighbor_vtx);

//   int pos2 = 0;
//   int neighbor_vtx_cur_idx = 0;
//   for (int i=0; i<n_elt; ++i) {
//     int face_vtx_idx = delt_vtx_idx[section_idx+i];
//     int n_face_vtx = delt_vtx_idx[section_idx+i+1] - delt_vtx_idx[section_idx+i];
//     PDM_g_num_t* first_face_vtx = delt_vtx+face_vtx_idx;

//     int n_neighbor = face_neighbors_idx[i+1] - face_neighbors_idx[i];
//     int parent_found = 0;
//     for (int j=0; j<n_neighbor; ++j) {
//       int n_vtx = neighbor_vtx_stri[0][pos2];
//       PDM_g_num_t* first_vtx = neighbor_vtx[0] + neighbor_vtx_cur_idx;
//       neighbor_vtx_cur_idx += n_vtx;
//       if (_is_parent(first_vtx,n_vtx,first_face_vtx,n_face_vtx)) {
//         elt2d_part[i] = neighbor_part[0][pos2];
//         parent_found = 1;
//       }
//       ++pos2;
//     }
//     assert(parent_found);
//   }

//   free(face_neighbors);
// }

// void
// _run_ppart_zone_nodal
// (
//   PDM_dmesh_nodal_t *dmesh_nodal,
//   _part_mesh_t      *pmesh,
//   PDM_split_dual_t   split_method,
//   int                dn_part,
//   PDM_MPI_Comm       comm
// )
// {
//   // printf("run_ppart_zone_nodal\n");
//   // TODO: joins
//   int i_rank;
//   int n_rank;
//   PDM_MPI_Comm_rank(comm, &i_rank);
//   PDM_MPI_Comm_size(comm, &n_rank);

//   // -1. misc
//   int n_bnd = 0; // TODO
//   int n_join = 0; // TODO
//   int* joins_ids = NULL; // TODO
//   pmesh->n_bounds  = n_bnd;
//   pmesh->n_joins   = n_join;
//   pmesh->joins_ids = (int *) malloc(n_join * sizeof(int));
//   for (int i_join = 0; i_join < n_join; i_join++){
//     pmesh->joins_ids[i_join] = joins_ids[i_join];
//   }
//   int tn_part;
//   PDM_MPI_Allreduce(&dn_part, &tn_part, 1, PDM_MPI_INT, PDM_MPI_SUM, comm);
//   pmesh->tn_part = tn_part;

//   // 0. concat sections
//   /// 0.0. all elts
//   int* section_idx;
//   int* delt_vtx_idx;
//   PDM_g_num_t* delt_vtx;
//   int n_section = PDM_concat_elt_sections(dmesh_nodal,&section_idx,&delt_vtx_idx,&delt_vtx);
//   /// 0.1. cells only
//   int* cell_section_idx;
//   int* dcell_vtx_idx;
//   PDM_g_num_t* dcell_vtx;
//   int n_cell_section = PDM_concat_cell_sections(dmesh_nodal,&cell_section_idx,&dcell_vtx_idx,&dcell_vtx);

//   // 1. distributions
//   int dn_vtx = dmesh_nodal->vtx->n_vtx;
//   PDM_g_num_t* vtx_dist = PDM_compute_entity_distribution(comm, dn_vtx);

//   int dn_elt = section_idx[n_section];
//   PDM_g_num_t* elt_dist = PDM_compute_entity_distribution(comm, dn_elt);

//   int dn_cell = cell_section_idx[n_cell_section];
//   PDM_g_num_t* cell_dist = PDM_compute_entity_distribution(comm, dn_cell);

//   // 2. elt_elt and cell_cell graph
//   // 2.0. elt_elt
//   PDM_g_num_t* elt_elt_idx;
//   PDM_g_num_t* elt_elt;
//   PDM_dmesh_nodal_dual_graph(vtx_dist,elt_dist,delt_vtx_idx,delt_vtx,&elt_elt_idx,&elt_elt,comm);
//   // 2.1. cell_cell
//   PDM_g_num_t* cell_cell_idx;
//   PDM_g_num_t* cell_cell;
//   PDM_dmesh_nodal_dual_graph(vtx_dist,cell_dist,dcell_vtx_idx,dcell_vtx,&cell_cell_idx,&cell_cell,comm);

//   // 3. partitioning
//   // 3.0 call graph partitionning on cells
//   double* part_fractions = NULL;
//   int* elt_part = (int*) malloc(dn_elt * sizeof(int));
//   for (int i=0; i<elt_elt_idx[dn_elt]; ++i) {
//     elt_elt[i]--;
//   }
//   int* cell_part = (int*) malloc(dn_cell * sizeof(int));
//   for (int i=0; i<cell_cell_idx[dn_cell]; ++i) {
//     cell_cell[i]--;
//   }
//   PDM_para_graph_split(split_method,
//                        cell_dist,
//                        cell_cell_idx,
//                        cell_cell,
//                        NULL, NULL,
//                        tn_part,
//                        part_fractions,
//                        cell_part,
//                        comm);


//   // 3.1 create the complete elt_part
//   // 3.1.0 Fill for 3D elements
//   for (int i_section=0; i_section<n_section; ++i_section) {
//     int* elt_section_part = elt_part + section_idx[i_section];

//     PDM_Mesh_nodal_elt_t type = PDM_DMesh_nodal_section_elt_type_get(dmesh_nodal,i_section);
//     int counter=0;
//     if (PDM_Mesh_nodal_is_3D_element(type)) {
//       int n_elt = section_idx[i_section+1] - section_idx[i_section];
//       for (int i=0; i<n_elt; ++i) {
//         elt_section_part[i] = cell_part[counter];
//         ++counter;
//       }
//     }
//   }
//   // 3.1.1 2D elements must go to the partition of their 3D parent
//   for (int i_section=0; i_section<n_section; ++i_section) {
//     int* elt_section_part = elt_part + section_idx[i_section];

//     PDM_Mesh_nodal_elt_t type = PDM_DMesh_nodal_section_elt_type_get(dmesh_nodal,i_section);
//     if (PDM_Mesh_nodal_is_2D_element(type)) {
//       int n_elt = section_idx[i_section+1] - section_idx[i_section];
//       face_part_from_parent(n_elt,
//                             dn_elt,
//                             section_idx[i_section],
//                             elt_dist,
//                             elt_part,
//                             delt_vtx_idx,
//                             delt_vtx,
//                             elt_elt_idx,
//                             elt_elt,
//                             elt_section_part,
//                             comm);
//     }
//   }

//   // 4. ln_to_gn for elt sections and cells
//   PDM_g_num_t* part_distri = PDM_compute_entity_distribution(comm, dn_part );
//   // 4.0 de-concetenate
//   int** pn_elt_section = (int**)malloc(n_section * sizeof(int*));
//   PDM_g_num_t*** pelt_section_ln_to_gn = (PDM_g_num_t***)malloc(n_section * sizeof(PDM_g_num_t**));
//   PDM_g_num_t** elt_section_distri = (PDM_g_num_t**)malloc(n_section * sizeof(PDM_g_num_t*));
//   for (int i_section=0; i_section<n_section; ++i_section) {
//     elt_section_distri[i_section] = PDM_DMesh_nodal_section_distri_std_get(dmesh_nodal,i_section);
//   }
//   // 4.1 ln_to_gn for elt sections
//   for (int i_section=0; i_section<n_section; ++i_section) {
//     int* elt_section_part = elt_part + section_idx[i_section];
//     PDM_part_assemble_partitions(comm,
//                                  part_distri,
//                                  elt_section_distri[i_section],
//                                  elt_section_part,
//                                 &pn_elt_section[i_section],
//                                 &pelt_section_ln_to_gn[i_section]);
//   }

//   // 4.2 ln_to_gn for cells
//   int* pn_cell = (int*)malloc(dn_part * sizeof(int));
//   for (int i_part=0; i_part<dn_part; ++i_part) {
//     pn_cell[i_part] = 0;
//     for (int i_section=0; i_section<n_section; ++i_section) {
//       PDM_Mesh_nodal_elt_t type = PDM_DMesh_nodal_section_elt_type_get(dmesh_nodal,i_section);
//       if (PDM_Mesh_nodal_is_3D_element(type)) { // only counting cells
//         pn_cell[i_part] += pn_elt_section[i_section][i_part];
//       }
//     }
//   }

//   PDM_g_num_t** pcell_ln_to_gn = (PDM_g_num_t**)malloc(dn_part * sizeof(PDM_g_num_t*));
//   for (int i_part=0; i_part<dn_part; ++i_part) {
//     pcell_ln_to_gn[i_part] = (PDM_g_num_t*)malloc(pn_cell[i_part]* sizeof(PDM_g_num_t));
//   }

//   for (int i_part=0; i_part<dn_part; ++i_part) {
//     int pos = 0;
//     int offset = 0;
//     for (int i_section=0; i_section<n_section; ++i_section) {
//       PDM_Mesh_nodal_elt_t type = PDM_DMesh_nodal_section_elt_type_get(dmesh_nodal,i_section);
//       if (PDM_Mesh_nodal_is_3D_element(type)) { // only counting cells
//         int pn_elt = pn_elt_section[i_section][i_part];
//         for (int i=0; i<pn_elt; ++i) {
//           pcell_ln_to_gn[i_part][pos] = pelt_section_ln_to_gn[i_section][i_part][i] + offset;
//           ++pos;
//         }
//         int n_elt_section = elt_section_distri[i_section][n_rank];
//         offset += n_elt_section;
//       }
//     }
//   }
//   // 4.3 ln_to_gn for elts
//   int* pn_elt = (int*)malloc(dn_part * sizeof(int));
//   for (int i_part=0; i_part<dn_part; ++i_part) {
//     pn_elt[i_part] = 0;
//     for (int i_section=0; i_section<n_section; ++i_section) {
//       pn_elt[i_part] += pn_elt_section[i_section][i_part];
//     }
//   }
//   PDM_g_num_t** pelt_ln_to_gn = (PDM_g_num_t**)malloc(dn_part * sizeof(PDM_g_num_t*));
//   for (int i_part=0; i_part<dn_part; ++i_part) {
//     pelt_ln_to_gn[i_part] = (PDM_g_num_t*)malloc(pn_elt[i_part]* sizeof(PDM_g_num_t));
//   }
//   for (int i_part=0; i_part<dn_part; ++i_part) {
//     int pos = 0;
//     int offset = 0;
//     for (int i_section=0; i_section<n_section; ++i_section) {
//       int part_n_elt = pn_elt_section[i_section][i_part];
//       for (int i=0; i<part_n_elt; ++i) {
//         pelt_ln_to_gn[i_part][pos] = pelt_section_ln_to_gn[i_section][i_part][i] + offset;
//         ++pos;
//       }
//       int n_elt_section = elt_section_distri[i_section][n_rank];
//       offset += n_elt_section;
//     }
//   }
//   free(elt_part);

//   // 5. reconstruct elts on partitions
//   int* pn_vtx;
//   PDM_g_num_t** pvtx_ln_to_gn;
//   int*** pelt_vtx_idx;
//   int*** pelt_vtx;
//   PDM_part_multi_dconnectivity_to_pconnectivity_sort(comm,
//                                                      dn_part,
//                                                      n_section,
//                                                      section_idx,
//                                                      elt_section_distri,
//                                                      delt_vtx_idx,
//                                                      delt_vtx,
//                                                      pn_elt_section,
//                                                      pelt_section_ln_to_gn,
//                                                     &pn_vtx,
//                                                     &pvtx_ln_to_gn,
//                                                     &pelt_vtx_idx,
//                                                     &pelt_vtx);
//   free(elt_section_distri);

//   // 6. coordinates
//   const double* dvtx_coord = dmesh_nodal->vtx->_coords;
//   double** pvtx_coord = NULL;
//   PDM_part_dcoordinates_to_pcoordinates(comm,
//                                         dn_part,
//                                         vtx_dist,
//                                         dvtx_coord,
//                                         pn_vtx,
//                  (const PDM_g_num_t **) pvtx_ln_to_gn,
//                                        &pvtx_coord);

//   // 7. Fill _part_t structures
//   pmesh->parts = (_part_t **) malloc(dn_part*sizeof(_part_t*));
//   for (int i_part = 0; i_part < dn_part; ++i_part) {
//     pmesh->parts[i_part] = _part_create();

//     pmesh->parts[i_part]->n_cell = pn_cell[i_part];
//     pmesh->parts[i_part]->n_vtx  = pn_vtx[i_part];

//     pmesh->parts[i_part]->vtx = pvtx_coord[i_part];

//     pmesh->parts[i_part]->cell_ln_to_gn = pcell_ln_to_gn[i_part];
//     pmesh->parts[i_part]->vtx_ln_to_gn  = pvtx_ln_to_gn[i_part];

//     pmesh->parts[i_part]->n_section            = n_section;
//     pmesh->parts[i_part]->n_elt                = (int          *) malloc(n_section * sizeof(int          ));
//     pmesh->parts[i_part]->elt_section_ln_to_gn = (PDM_g_num_t **) malloc(n_section * sizeof(PDM_g_num_t *));
//     pmesh->parts[i_part]->elt_vtx_idx          = (int         **) malloc(n_section * sizeof(int         *));
//     pmesh->parts[i_part]->elt_vtx              = (int         **) malloc(n_section * sizeof(int         *));
//     for (int i_section=0; i_section<n_section; ++i_section) {
//       pmesh->parts[i_part]->n_elt               [i_section] = pn_elt_section[i_section][i_part];
//       pmesh->parts[i_part]->elt_section_ln_to_gn[i_section] = pelt_section_ln_to_gn[i_section][i_part];
//       pmesh->parts[i_part]->elt_vtx_idx         [i_section] = pelt_vtx_idx[i_section][i_part];
//       pmesh->parts[i_part]->elt_vtx             [i_section] = pelt_vtx[i_section][i_part];
//     }
//   }
//   free(pcell_ln_to_gn);
//   free(pelt_section_ln_to_gn);


//   /* Let's call graph communication to do the reordering - Temporary  */
//   int** pinternal_vtx_bound_proc_idx  = NULL;
//   int** pinternal_vtx_bound_part_idx  = NULL;
//   int** pinternal_vtx_bound           = NULL;
//   int** pinternal_vtx_priority        = NULL;
//   PDM_part_generate_entity_graph_comm(comm,
//                                       part_distri,
//                                       vtx_dist,
//                                       dn_part,
//                                       pn_vtx,
//                (const PDM_g_num_t **) pvtx_ln_to_gn,
//                                       NULL,
//                                      &pinternal_vtx_bound_proc_idx,
//                                      &pinternal_vtx_bound_part_idx,
//                                      &pinternal_vtx_bound,
//                                      &pinternal_vtx_priority);


//   for (int i_part = 0; i_part<dn_part; ++i_part) {
//     pmesh->parts[i_part]->vtx_ghost_information = (int*) malloc(pn_vtx[i_part]*sizeof(int));
//     for (int i=0; i<pn_vtx[i_part]; ++i) {
//       pmesh->parts[i_part]->vtx_ghost_information[i] = pinternal_vtx_priority[i_part][i];
//     }
//   }

//   _setup_ghost_information(i_rank, dn_part, pn_vtx, pinternal_vtx_priority);
//   PDM_part_renum_vtx(pmesh->parts, dn_part, 1, (void *) pinternal_vtx_priority);
//   // reorder vtx priority itself
//   for (int i_part=0; i_part<dn_part; ++i_part) {
//     int n_kind = 3;
//     int* priority_count = PDM_array_zeros_int(n_kind);

//     for (int i=0; i<pn_vtx[i_part]; ++i) {
//       ++priority_count[pinternal_vtx_priority[i_part][i]];
//     }
//     // reorder
//     int* first = pinternal_vtx_priority[i_part];
//     for (int k=0; k<n_kind; ++k) {
//       for (int i=0; i<priority_count[k]; ++i) {
//         *first++ = k;
//       }
//     }
//     free(pinternal_vtx_priority[i_part]);
//     free(priority_count);
//   }
//   free(pinternal_vtx_priority);


//   /* Free in order to be correct */
//   for (int ipart = 0; ipart < dn_part; ipart++) {
//     free(pinternal_vtx_bound_proc_idx[ipart]);
//     free(pinternal_vtx_bound_part_idx[ipart]);
//     free(pinternal_vtx_bound[ipart]);
//   }
//   free(pinternal_vtx_bound_proc_idx);
//   free(pinternal_vtx_bound_part_idx);
//   free(pinternal_vtx_bound);


//   //// Now genererate bounds and comm data -- we need update pface_ln_to_gn which has been modified
//   //for (int ipart = 0; ipart < n_part; ipart++) {
//   //  pface_ln_to_gn[ipart] = pmeshes->parts[ipart]->face_ln_to_gn;
//   //}

//   int           n_group_elmt    = dmesh_nodal->n_group_elmt;
//   int        *  dgroup_elmt_idx = dmesh_nodal->dgroup_elmt_idx;
//   PDM_g_num_t*  dgroup_elmt     = dmesh_nodal->dgroup_elmt;
//   int        ** pface_bound_idx      = NULL;
//   int        ** pface_bound          = NULL;
//   PDM_g_num_t** pface_bound_ln_to_gn = NULL;
//   PDM_part_distgroup_to_partgroup(comm,
//                                   elt_dist,
//                                   n_group_elmt,
//                                   dgroup_elmt_idx,
//                                   dgroup_elmt,
//                                   dn_part,
//                                   pn_elt,
//            (const PDM_g_num_t **) pelt_ln_to_gn,
//                                  &pface_bound_idx,
//                                  &pface_bound,
//                                  &pface_bound_ln_to_gn);
//   //int **pface_join_tmp = NULL;
//   //PDM_part_distgroup_to_partgroup(comm,
//   //                                face_distri,
//   //                                n_join,
//   //                                dface_join_idx,
//   //                                dface_join,
//   //                                n_part,
//   //                                pn_face,
//   //         (const PDM_g_num_t **) pface_ln_to_gn,
//   //                               &pface_join_idx,
//   //                               &pface_join_tmp,
//   //                               &pface_join_ln_to_gn);

//   //PDM_part_generate_entity_graph_comm(comm,
//   //                                    part_distri,
//   //                                    face_distri,
//   //                                    n_part,
//   //                                    pn_face,
//   //             (const PDM_g_num_t **) pface_ln_to_gn,
//   //                                    NULL,
//   //                                   &pinternal_face_bound_proc_idx,
//   //                                   &pinternal_face_bound_part_idx,
//   //                                   &pinternal_face_bound,
//   //                                    NULL);
//   PDM_part_generate_entity_graph_comm(comm,
//                                       part_distri,
//                                       vtx_dist,
//                                       dn_part,
//                                       pn_vtx,
//                (const PDM_g_num_t **) pvtx_ln_to_gn,
//                                       NULL,
//                                      &pinternal_vtx_bound_proc_idx,
//                                      &pinternal_vtx_bound_part_idx,
//                                      &pinternal_vtx_bound,
//                                       NULL); // do NOT compute a new priority, else a new reordering would be needed

//   // Finally complete parts structure with internal join data and bounds
//   pmesh->n_bounds = n_group_elmt;
//   for (int ipart = 0; ipart < dn_part; ++ipart) {
//     pmesh->parts[ipart]->face_bound_idx      = pface_bound_idx[ipart];
//     pmesh->parts[ipart]->face_bound          = pface_bound[ipart];
//     pmesh->parts[ipart]->face_bound_ln_to_gn = pface_bound_ln_to_gn[ipart];

//     ///* For face_join, the function only returns local id of face in join, we have to
//     //   allocate to set up expected size (4*nb_face_join) */
//     //pmesh->parts[ipart]->face_join_idx = pface_join_idx[ipart];
//     //int s_face_join   pmesh->parts[ipart]->face_join_idx[n_join];
//     //pmesh->parts[ipart]->face_join = (int *) malloc( 4 * s_face_join * sizeof(int));
//     //for (int i_face = 0; i_face < s_face_join; i_face++){
//     //  pmesh->parts[ipart]->face_join[4*i_face] = pface_join_tmp[ipart][i_face];
//     //}
//     //pmesh->parts[ipart]->face_join_ln_to_gn = pface_join_ln_to_gn[ipart];
//     //free(pface_join_tmp[ipart]);

//     //pmesh->parts[ipart]->face_part_bound_proc_idx = pinternal_face_bound_proc_idx[ipart];
//     //pmesh->parts[ipart]->face_part_bound_part_idx = pinternal_face_bound_part_idx[ipart];
//     //pmesh->parts[ipart]->face_part_bound          = pinternal_face_bound[ipart];

//     pmesh->parts[ipart]->vtx_part_bound_proc_idx = pinternal_vtx_bound_proc_idx[ipart];
//     pmesh->parts[ipart]->vtx_part_bound_part_idx = pinternal_vtx_bound_part_idx[ipart];
//     pmesh->parts[ipart]->vtx_part_bound          = pinternal_vtx_bound[ipart];
//   }

//   free(section_idx);
//   free(delt_vtx_idx);
//   free(delt_vtx);

//   free(elt_dist);
//   free(vtx_dist);
//   free(elt_elt_idx);
//   free(elt_elt);

//   free(pn_elt);
//   for(int i_part=0; i_part<dn_part; ++i_part) {
//     free(pelt_ln_to_gn[i_part]);
//   }
//   free(pelt_ln_to_gn);
// }

// static
// void
// _run_ppart_zone_nodal2
// (
//   PDM_dmesh_nodal_t *dmesh_nodal,
//   _part_mesh_t      *pmesh,
//   PDM_split_dual_t   split_method,
//   int                dn_part,
//   PDM_MPI_Comm       comm
// )
// {
//   PDM_UNUSED(dmesh_nodal);
//   PDM_UNUSED(pmesh);
//   PDM_UNUSED(split_method);
//   PDM_UNUSED(dn_part);
//   PDM_UNUSED(comm);
//   // printf("_run_ppart_zone_nodal2\n");
//   // TODO: joins
//   // int i_rank;
//   // int n_rank;
//   // PDM_MPI_Comm_rank(comm, &i_rank);
//   // PDM_MPI_Comm_size(comm, &n_rank);

//   // PDM_dmesh_nodal_to_dmesh_t* dmn_to_dm = PDM_dmesh_nodal_to_dmesh_create(1, comm, PDM_OWNERSHIP_KEEP);
//   // PDM_dmesh_nodal_to_dmesh_add_dmesh_nodal(dmn_to_dm, 0, dmn);

//   // if(dmn->mesh_dimension == 2){
//   //   PDM_dmesh_nodal_to_dmesh_compute2(dmn_to_dm,
//   //                                     PDM_DMESH_NODAL_TO_DMESH_TRANSFORM_TO_FACE,
//   //                                     PDM_DMESH_NODAL_TO_DMESH_TRANSLATE_GROUP_TO_FACE);
//   // } else if(dmn->mesh_dimension == 2 ) {
//   //   PDM_dmesh_nodal_to_dmesh_compute2(dmn_to_dm,
//   //                                     PDM_DMESH_NODAL_TO_DMESH_TRANSFORM_TO_EDGE,
//   //                                     PDM_DMESH_NODAL_TO_DMESH_TRANSLATE_GROUP_TO_EDGE);
//   // } else {
//   //   PDM_error(__FILE__, __LINE__, 0, "PDM_multipart error : Bad dmesh_nodal mesh_dimension \n");
//   // }
// }

/**

 * \brief Construct the partitioned meshes on every zones
 *
 * \param [in]   mpart_id          Multipart structure id
 */
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


  /*
   * Step 1 : Split the graph (If we have a dmesh_nodal and prepare the dmesh before the treatment for faces and elements are the same)
   * Step 2 : Rebuild all connectivity in coherent manner
   * Step 3 : Apply all reordering
   * Step 4 : Deduce all mesh_nodal connectivity
   */

  if (_multipart->merge_blocks)
  {
    // 1. Generate global numerotation using all blocks
    // 2. Call the partitionner once on the global numbering
  } else {
    PDM_timer_t *timer = PDM_timer_create();
    double cum_elapsed_time = 0;
    int *starting_part_idx =  PDM_array_new_idx_from_sizes_int(_multipart->n_part, _multipart->n_zone);

    int is_by_elt = 0;
    for (int i_zone = 0; i_zone < _multipart->n_zone; ++i_zone) {
      PDM_dmesh_nodal_t* dmesh_nodal = _multipart->dmeshes_nodal[i_zone];
      if (dmesh_nodal != NULL) { // element representation
        is_by_elt = 1;
        // PDM_printf("Partitionning elt zone %d/%d \n", i_zone+1, _multipart->n_zone);
        PDM_MPI_Comm comm = _multipart->comm;
        PDM_split_dual_t split_method = _multipart->split_method;
        int n_part = _multipart->n_part[i_zone];
        _part_mesh_t* pmesh = &(_multipart->pmeshes[i_zone]);

        // _run_ppart_zone_nodal(dmesh_nodal,pmesh,split_method,n_part,comm);

        const double* part_fraction      = &_multipart->part_fraction[starting_part_idx[i_zone]];
        PDM_part_size_t part_size_method = _multipart->part_size_method;
        PDM_dmesh_t  *_dmeshes           = _multipart->dmeshes[i_zone];
        _run_ppart_zone2(_dmeshes, dmesh_nodal, pmesh, n_part, split_method, part_size_method, part_fraction, comm);

      } else { // face representation
        // PDM_printf("Partitionning face zone %d/%d \n", i_zone+1, _multipart->n_zone);

        PDM_MPI_Comm comm = _multipart->comm;

        PDM_split_dual_t split_method    = _multipart->split_method;
        PDM_part_size_t part_size_method = _multipart->part_size_method;

        const double* part_fraction = &_multipart->part_fraction[starting_part_idx[i_zone]];

        PDM_dmesh_t  *_dmeshes =   _multipart->dmeshes[i_zone];
        _part_mesh_t *_pmeshes = &(_multipart->pmeshes[i_zone]);

        int n_part = _multipart->n_part[i_zone];


        if (i_rank == 0)
          PDM_printf("Running partitioning for block %i...\n", i_zone+1);
        PDM_timer_resume(timer);
        _run_ppart_zone(_dmeshes, _pmeshes, n_part, split_method, part_size_method, part_fraction, comm);
        PDM_timer_hang_on(timer);
        if (i_rank == 0)
          PDM_printf("...completed (elapsed time : %f)\n", PDM_timer_elapsed(timer) - cum_elapsed_time);
        cum_elapsed_time = PDM_timer_elapsed(timer);
      }
    }
    PDM_timer_free(timer);

    free(starting_part_idx);
    // Now rebuild joins over the zones
    if (!is_by_elt) { // WARNING TODO also implement if element representation
      _search_matching_joins(_multipart);
    }
  }
}

/*
 * \brief Rebuild part_mesh_nodal from dmesh_nodal after the compute
 */
void
PDM_multipart_compute_part_mesh_nodal
(
const int   id
)
{
  _pdm_multipart_t *_multipart = _get_from_id (id);
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(_multipart->comm, &i_rank);
  PDM_MPI_Comm_size(_multipart->comm, &n_rank);

  for (int i_zone = 0; i_zone < _multipart->n_zone; ++i_zone) {
    _part_mesh_t      *pmesh       = &(_multipart->pmeshes    [i_zone]);
    PDM_dmesh_t       *dmesh       = _multipart->dmeshes      [i_zone];
    PDM_dmesh_nodal_t *dmesh_nodal = _multipart->dmeshes_nodal[i_zone];

    int n_part = _multipart->n_part[i_zone];
    if(dmesh_nodal->mesh_dimension == 3){
      _compute_part_mesh_nodal_3d(dmesh, dmesh_nodal, pmesh, n_part);
    } else if(dmesh_nodal->mesh_dimension == 2){
      _compute_part_mesh_nodal_2d(dmesh, dmesh_nodal, pmesh, n_part);
    } else {
      PDM_error(__FILE__, __LINE__, 0, "PDM_multipart_compute_part_mesh_nodal error : Bad dmesh_nodal dimension \n");
    }
  }
}

/**
 *
 * \brief Returns the dimensions of a given partition
 */
void
PDM_multipart_part_dim_get
(
const int   mpart_id,
const int   i_zone,
const int   i_part,
      int  *n_section,
      int **n_elt,
      int  *n_cell,
      int  *n_face,
      int  *n_face_part_bound,
      int  *n_vtx,
      int  *n_proc,
      int  *n_total_part,
      int  *s_cell_face,
      int  *s_face_vtx,
      int  *s_face_bound,
      int  *n_bound_groups,
      int  *s_face_join,
      int  *n_join_groups
)
{
  _pdm_multipart_t *_multipart = _get_from_id (mpart_id);

  assert(i_zone < _multipart->n_zone && i_part < _multipart->n_part[i_zone]);
  _part_mesh_t _pmeshes = _multipart->pmeshes[i_zone];

  *n_cell = _pmeshes.parts[i_part]->n_cell;
  *n_face = _pmeshes.parts[i_part]->n_face;
  *n_vtx  = _pmeshes.parts[i_part]->n_vtx;

  PDM_MPI_Comm_size(_multipart->comm, n_proc);
  *n_total_part = _pmeshes.tn_part;

  *n_section = _pmeshes.parts[i_part]->n_section;
  *n_elt = _pmeshes.parts[i_part]->n_elt;
  if (*n_section>0) {
    *s_cell_face = -1;
    *s_face_vtx  = -1;

    *n_face_part_bound = -1;

    *n_join_groups  = -1;
    *s_face_join    = -1;
  } else {
    *s_cell_face = _pmeshes.parts[i_part]->cell_face_idx[*n_cell];
    *s_face_vtx  = _pmeshes.parts[i_part]->face_vtx_idx[*n_face];

    *n_face_part_bound = _pmeshes.parts[i_part]->face_part_bound_part_idx[*n_total_part];

    *n_join_groups  = _pmeshes.n_joins;
    *s_face_join    = _pmeshes.parts[i_part]->face_join_idx[*n_join_groups];
  }
  *n_bound_groups = _pmeshes.n_bounds;
  *s_face_bound   = _pmeshes.parts[i_part]->face_bound_idx[*n_bound_groups];
}


void
PDM_multipart_part_graph_comm_vtx_dim_get
(
 const int   mpart_id,
 const int   i_zone,
 const int   i_part,
       int  *n_vtx_part_bound
)
{
  _pdm_multipart_t *_multipart = _get_from_id (mpart_id);

  assert(i_zone < _multipart->n_zone && i_part < _multipart->n_part[i_zone]);
  _part_mesh_t _pmeshes = _multipart->pmeshes[i_zone];

  // printf(" n_vtx_part_bound = %i \n", _pmeshes.parts[i_part]->vtx_part_bound_part_idx[_pmeshes.tn_part]);
  // printf(" _pmeshes.tn_part = %i \n", _pmeshes.tn_part);

  *n_vtx_part_bound = _pmeshes.parts[i_part]->vtx_part_bound_part_idx[_pmeshes.tn_part];
}

/**
 *
 * \brief Returns the data arrays of a given partition
 */
void
PDM_multipart_part_val_get
(
const int            mpart_id,
const int            i_zone,
const int            i_part,
      int         ***elt_vtx_idx,
      int         ***elt_vtx,
      PDM_g_num_t ***elt_section_ln_to_gn,
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

  assert(i_zone < _multipart->n_zone && i_part < _multipart->n_part[i_zone]);
  _part_mesh_t _pmeshes = _multipart->pmeshes[i_zone];

  *cell_tag = NULL;
  *face_tag = NULL;
  *vtx_tag  = NULL;

  *cell_ln_to_gn = _pmeshes.parts[i_part]->cell_ln_to_gn;
  *face_ln_to_gn = _pmeshes.parts[i_part]->face_ln_to_gn;
  *vtx_ln_to_gn  = _pmeshes.parts[i_part]->vtx_ln_to_gn;

  *cell_face_idx = _pmeshes.parts[i_part]->cell_face_idx;
  *cell_face     = _pmeshes.parts[i_part]->cell_face;
  *face_cell     = _pmeshes.parts[i_part]->face_cell;
  *face_vtx_idx  = _pmeshes.parts[i_part]->face_vtx_idx;
  *face_vtx      = _pmeshes.parts[i_part]->face_vtx;

  *vtx           = _pmeshes.parts[i_part]->vtx;

  *face_part_bound_proc_idx = _pmeshes.parts[i_part]->face_part_bound_proc_idx;
  *face_part_bound_part_idx = _pmeshes.parts[i_part]->face_part_bound_part_idx;
  *face_part_bound          = _pmeshes.parts[i_part]->face_part_bound;

  *face_bound_idx       = _pmeshes.parts[i_part]->face_bound_idx;
  *face_bound           = _pmeshes.parts[i_part]->face_bound;
  *face_bound_ln_to_gn  = _pmeshes.parts[i_part]->face_bound_ln_to_gn;
  *face_join_idx        = _pmeshes.parts[i_part]->face_join_idx;
  *face_join            = _pmeshes.parts[i_part]->face_join;
  *face_join_ln_to_gn   = _pmeshes.parts[i_part]->face_join_ln_to_gn;

  *elt_vtx_idx = _pmeshes.parts[i_part]->elt_vtx_idx;
  *elt_vtx     = _pmeshes.parts[i_part]->elt_vtx;
  *elt_section_ln_to_gn = _pmeshes.parts[i_part]->elt_section_ln_to_gn;
}


// void PDM_multipart_part_get
// (
//  const int            mpart_id,
//  const int            i_zone,
//  const int            i_part,
//  const result_type    res,
//                    **vtx_part_bound
// )
// {
//   if(res == VTX_COMM_GRAPH)
//   {
//     *vtx_part_bound_proc_idx = _pmeshes.parts[i_part]->vtx_part_bound_proc_idx;
//     *vtx_part_bound_part_idx = _pmeshes.parts[i_part]->vtx_part_bound_part_idx;
//     *vtx_part_bound          = _pmeshes.parts[i_part]->vtx_part_bound;
//   }
// }

void
PDM_multipart_part_graph_comm_vtx_data_get
(
const int            mpart_id,
const int            i_zone,
const int            i_part,
      int          **vtx_part_bound_proc_idx,
      int          **vtx_part_bound_part_idx,
      int          **vtx_part_bound
)
{
   _pdm_multipart_t *_multipart = _get_from_id (mpart_id);

   // table_res[] = is_getted;

  assert(i_zone < _multipart->n_zone && i_part < _multipart->n_part[i_zone]);
  _part_mesh_t _pmeshes = _multipart->pmeshes[i_zone];

  *vtx_part_bound_proc_idx = _pmeshes.parts[i_part]->vtx_part_bound_proc_idx;
  *vtx_part_bound_part_idx = _pmeshes.parts[i_part]->vtx_part_bound_part_idx;
  *vtx_part_bound          = _pmeshes.parts[i_part]->vtx_part_bound;
}


void
PDM_multipart_part_color_get
(
const int            mpart_id,
const int            i_zone,
const int            i_part,
      int          **cell_color,
      int          **face_color,
      int          **face_hp_color,
      int          **thread_color,
      int          **hyperplane_color
)
{
  _pdm_multipart_t *_multipart = _get_from_id (mpart_id);

  assert(i_zone < _multipart->n_zone && i_part < _multipart->n_part[i_zone]);
  _part_mesh_t _pmeshes = _multipart->pmeshes[i_zone];
  *cell_color       = _pmeshes.parts[i_part]->cell_color;
  *face_color       = _pmeshes.parts[i_part]->face_color;
  *face_hp_color    = _pmeshes.parts[i_part]->face_hp_color;
  *thread_color     = _pmeshes.parts[i_part]->thread_color;
  *hyperplane_color = _pmeshes.parts[i_part]->hyperplane_color;

}

void
PDM_multipart_part_ghost_infomation_get
(
const int            mpart_id,
const int            i_zone,
const int            i_part,
      int          **vtx_ghost_information
)
{
  _pdm_multipart_t *_multipart = _get_from_id (mpart_id);

  assert(i_zone < _multipart->n_zone && i_part < _multipart->n_part[i_zone]);
  _part_mesh_t _pmeshes = _multipart->pmeshes[i_zone];

  *vtx_ghost_information = _pmeshes.parts[i_part]->vtx_ghost_information;
}

void
PDM_multipart_time_get
(
const int       mpart_id,
const int       i_zone,
      double  **elapsed,
      double  **cpu,
      double  **cpu_user,
      double  **cpu_sys
)
{
  _pdm_multipart_t *_multipart = _get_from_id (mpart_id);
  assert(i_zone < _multipart->n_zone);

  // PDM_printf("PDM_multipart_time_get: Not implemented\n");
  *elapsed  = NULL;
  *cpu      = NULL;
  *cpu_user = NULL;
  *cpu_sys  = NULL;

}

/**
 *
 * \brief Free the structure
 *
 * \param [in]   mpart_id  Multipart structure id
 */

void
PDM_multipart_free
(
 const int mpart_id
)
{
  _pdm_multipart_t *_multipart = _get_from_id (mpart_id);

  // free(_multipart->dmeshes_ids);
  for (int izone = 0; izone < _multipart->n_zone; izone++) {
    if (_multipart->pmeshes[izone].joins_ids != NULL)
      free(_multipart->pmeshes[izone].joins_ids);
    if (_multipart->pmeshes[izone].parts != NULL){
      for (int ipart = 0; ipart < _multipart->n_part[izone]; ipart++) {
        _part_free(_multipart->pmeshes[izone].parts[ipart], _multipart->owner);
      }
      free(_multipart->pmeshes[izone].parts);
    }
  }
  free(_multipart->pmeshes);
  free(_multipart->dmeshes);
  free(_multipart->dmeshes_nodal);

  //PDM_part_renum_method_purge();
  free (_multipart);

  PDM_Handles_handle_free (_multiparts, mpart_id, PDM_FALSE);

  const int n_multipart = PDM_Handles_n_get (_multiparts);

  if (n_multipart == 0) {
    _multiparts = PDM_Handles_free (_multiparts);
  }
  // PDM_printf("Cleaned from PDM_multipart_free\n");
}


/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
