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
#include "pdm_partitioning_algorithm.h"
#include "pdm_para_graph_dual.h"
#include "pdm_part_priv.h"
#include "pdm_part_renum.h"
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

/**
 * \struct _part_mesh_t
 * \brief  This private structure stores partionned meshes obtained on a given
 *         zone.
 */
typedef struct {
  /* Data shared between all the partitions */
  int  tn_part;          // Total number of partitions created for this mesh
  int  n_bounds;         // Global number of boundary groups
  int  n_joins;          // Global number of interface groups
  int *joins_ids;        // Global id of each interface (size=n_joins)

  /* Partitions -- see pdm_part_priv.h for struct definition */
  _part_t  **parts;
  int        renum_cell_method;     // Choice of renumbering method for cells
  const int *renum_cell_properties; // Parameters used by some renumbering methods
  int        renum_face_method;     // Choice of renumbering method for faces
} _part_mesh_t;

/**
 * \struct _pdm_multipart_t
 * \brief  This structure describe a multipart.
 *         It includes distributed meshes, partioned meshes and
 *         partitioning parameters.
 */

typedef struct  {
  /* Multipart description */
  int           n_zone;              // Number of initial zones
  int          *dmeshes_ids;         // Ids of dmesh structure storing
                                     //   distributed meshes (size = n_zone)
  int           n_total_joins;       // Total number of joins between zones (each counts twice)
  const int    *join_to_opposite;    // For each global joinId, give the globalId of
                                     //   the opposite join (size = n_total_joins)
  PDM_MPI_Comm comm;                 // MPI communicator

  /* Partitioning parameters */
  PDM_bool_t       merge_blocks;     // Merge before partitionning or not
  PDM_split_dual_t split_method;     // Partitioning method (Metis or Scotch)
  PDM_part_size_t  part_size_method; // Procude homogeneous or heterogeneous partitions
  const int       *n_part;           // Number of wanted partitions per proc
                                     // in each zone (size = n_zone)
  const double    *part_fraction;    // Weight (in %) of each partition, in each zone
                                     //   (size = sum n_part[i]), if heterogeneous
  /* Partitioned meshes */
  _part_mesh_t *pmeshes;             // Partitioned meshes structures (size=n_zone)

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

  PDM_printf("pdm::_build_join_uface_distribution\n");
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
  int *nb_face_in_joins = (int *) malloc(n_unique_joins * sizeof(int));
  for (int i = 0; i < n_unique_joins; i++)
    nb_face_in_joins[i] = 0;

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
  for (int i=0; i < n_unique_joins; i++)
    _face_in_join_distri[i+1] = _face_in_join_distri[i+1] + _face_in_join_distri[i];

  /*
  PDM_printf("[%d] _face_in_join_distri : ", i_rank);
  for (int i = 0; i < n_unique_joins + 1; i++)
    PDM_printf(" %d ", _face_in_join_distri[i]);
  PDM_printf("\n");
  */

  free(nb_face_in_joins);
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
  free(block_stride);
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
 _part_t *part
)
{
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

  if (part->new_to_old_order_cell != NULL)
    free(part->new_to_old_order_cell);
  part->new_to_old_order_cell = NULL;

  if (part->new_to_old_order_face != NULL)
    free(part->new_to_old_order_face);
  part->new_to_old_order_face = NULL;

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

  _multipart->n_zone           = n_zone;
  _multipart->n_part           = n_part;
  _multipart->merge_blocks     = merge_blocks;
  _multipart->split_method     = split_method;
  _multipart->part_size_method = part_size_method;
  _multipart->part_fraction    = part_fraction;
  _multipart->comm             = comm;

  _multipart->n_total_joins    = 0;
  _multipart->join_to_opposite = NULL;

  _multipart->dmeshes_ids  = (int *) malloc(_multipart->n_zone * sizeof(int));
  _multipart->pmeshes = (_part_mesh_t *) malloc(_multipart->n_zone * sizeof(_part_mesh_t));

  int _renum_cell_method = PDM_part_renum_method_cell_idx_get("PDM_PART_RENUM_CELL_NONE");
  int _renum_face_method = PDM_part_renum_method_face_idx_get("PDM_PART_RENUM_FACE_NONE");
  for (int izone = 0; izone < _multipart->n_zone; izone++) {
    _multipart->dmeshes_ids[izone] = -1;
    _multipart->pmeshes[izone].renum_cell_method = _renum_cell_method;
    _multipart->pmeshes[izone].renum_face_method = _renum_face_method;
    _multipart->pmeshes[izone].renum_cell_properties = NULL;
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
 const int        mpart_id,
 const int        zone_id,
 const int        dmesh_id
)
{
  _pdm_multipart_t *_multipart = _get_from_id (mpart_id);
  assert(zone_id < _multipart->n_zone);
  _multipart->dmeshes_ids[zone_id] = dmesh_id;
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

/**
 *
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

  if (_multipart->merge_blocks)
  {
    // 1. Generate global numerotation using all blocks
    // 2. Call the partitionner once on the global numbering
  }
  else
  {
    // 2. Loop over the blocks and call the partitionner
    for (int i_zone = 0; i_zone < _multipart->n_zone; i_zone++) {
      PDM_printf("Partitionning zone %d/%d \n", i_zone+1, _multipart->n_zone);

      // Get distributed mesh data for this zone
      int dmesh_id = _multipart->dmeshes_ids[i_zone];
      int dn_cell, dn_face, dn_vtx, n_bnd, n_join;
      const double       *dvtx_coord;
      const int          *dface_vtx_idx;
      const PDM_g_num_t  *dface_vtx;
      const PDM_g_num_t  *dface_cell;
      const int          *dface_bound_idx;
      const PDM_g_num_t  *dface_bound;
      const int          *joins_ids;
      const int          *dface_join_idx;
      const PDM_g_num_t  *dface_join;
      PDM_dmesh_dims_get(dmesh_id, &dn_cell, &dn_face, &dn_vtx, &n_bnd, &n_join);
      PDM_dmesh_data_get(dmesh_id, &dvtx_coord, &dface_vtx_idx, &dface_vtx, &dface_cell,
                         &dface_bound_idx, &dface_bound, &joins_ids, &dface_join_idx, &dface_join);

      // This will store all the partitions created by this proc on this zone
      _part_mesh_t *_pmeshes = &(_multipart->pmeshes[i_zone]);

      //Copy number of bounds and joins (global data) in the part structure
      _pmeshes->n_bounds = n_bnd;
      _pmeshes->n_joins  = n_join;
      _pmeshes->joins_ids = (int *) malloc(n_join*sizeof(int));
      for (int i_join = 0; i_join < n_join; i_join++)
        _pmeshes->joins_ids[i_join] = joins_ids[i_join];

      // Compute total number of partitions for this zone
      int n_part = _multipart->n_part[i_zone];
      int tn_part;
      PDM_MPI_Allreduce(&n_part, &tn_part, 1, PDM_MPI_INT, PDM_MPI_SUM, _multipart->comm);
      _pmeshes->tn_part = tn_part;

      //Construct graph and split -- no interface for now, do it by hand
      PDM_g_num_t *cell_distri = PDM_compute_entity_distribution(_multipart->comm, dn_cell);
      PDM_g_num_t *face_distri = PDM_compute_entity_distribution(_multipart->comm, dn_face);
      PDM_g_num_t *vtx_distri  = PDM_compute_entity_distribution(_multipart->comm, dn_vtx);
      PDM_g_num_t *part_distri = PDM_compute_entity_distribution(_multipart->comm, n_part);

      PDM_g_num_t *dual_graph_idx;
      int         *dcell_face_idx;
      PDM_g_num_t *dual_graph, *dcell_face;

      PDM_para_graph_dual_from_arc2node(_multipart->comm,
                                        cell_distri,
                                        face_distri,
                                        dface_cell,
                                       &dual_graph_idx,
                                       &dual_graph,
                                        1,
                                       &dcell_face_idx,
                                       &dcell_face);

      double *part_fraction = NULL;
      if (_multipart->part_size_method == PDM_PART_SIZE_HETEROGENEOUS){
        int *n_part_per_rank = (int *)    malloc(n_rank * sizeof(int));
        int *displ           = (int *)    malloc(n_rank * sizeof(int));
        part_fraction        = (double *) malloc(tn_part * sizeof(double));
        int starting_part_idx = 0;
        for (int i =0; i < n_rank; i++){
          n_part_per_rank[i] = part_distri[i+1] - part_distri[i];
          displ[i] = part_distri[i]-1;
        }
        for (int jzone = 0; jzone < i_zone; jzone++)
          starting_part_idx += _multipart->n_part[jzone];

        PDM_MPI_Allgatherv((void*) &_multipart->part_fraction[starting_part_idx],
                           n_part,
                           PDM_MPI_DOUBLE,
                           part_fraction,
                           n_part_per_rank,
                           displ,
                           PDM_MPI_DOUBLE,
                           _multipart->comm);
        free(n_part_per_rank);
        free(displ);
      }
      int *cell_part = (int *) malloc(dn_cell * sizeof(int));
      PDM_para_graph_split (_multipart->split_method,
                            cell_distri,
                            dual_graph_idx,
                            dual_graph,
                            NULL, NULL,
                            tn_part,
                            part_fraction,
                            cell_part,
                            _multipart->comm);

      free(dual_graph_idx);
      free(dual_graph);
      if (_multipart->part_size_method == PDM_PART_SIZE_HETEROGENEOUS)
        free(part_fraction);

      //Partitioning algorithm except to work on 2d arrays (external size = n_part) whereas
      //_part_t struture stores n_part * 1d arrays so we have to use tmp pointers
      int          *pn_vtx                       = NULL;
      int          *pn_cell                      = NULL;
      int          *pn_face                      = NULL;
      double      **pvtx_coord                   = NULL;
      int         **pface_vtx_idx                = NULL;
      int         **pface_vtx                    = NULL;
      int         **pcell_face_idx               = NULL;
      int         **pcell_face                   = NULL;
      int         **pface_cell                   = NULL;
      int         **pface_bound_idx              = NULL;
      int         **pface_bound                  = NULL;
      int         **pface_join_idx               = NULL;
      int         **pinternal_face_bound_procidx = NULL;
      int         **pinternal_face_bound_partidx = NULL;
      int         **pinternal_face_bound         = NULL;
      PDM_g_num_t **pcell_ln_to_gn               = NULL;
      PDM_g_num_t **pface_ln_to_gn               = NULL;
      PDM_g_num_t **pvtx_ln_to_gn                = NULL;
      PDM_g_num_t **pface_bound_ln_to_gn         = NULL;
      PDM_g_num_t **pface_join_ln_to_gn          = NULL;

      PDM_part_assemble_partitions(_multipart->comm,
                                   part_distri,
                                   cell_distri,
                                   cell_part,
                                  &pn_cell,
                                  &pcell_ln_to_gn);
      free(cell_part);
      PDM_part_dconnectivity_to_pconnectivity_sort(_multipart->comm,
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
      PDM_part_dconnectivity_to_pconnectivity_sort(_multipart->comm,
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
      PDM_part_reorient_bound_faces(n_part, pn_face, pface_cell,
                                    (const int **) pcell_face_idx, pcell_face,
                                    (const int **) pface_vtx_idx, pface_vtx);

      PDM_part_dcoordinates_to_pcoordinates(_multipart->comm,
                                            n_part,
                                            vtx_distri,
                                            dvtx_coord,
                                            pn_vtx,
                                            (const PDM_g_num_t **) pvtx_ln_to_gn,
                                           &pvtx_coord);
      PDM_part_distgroup_to_partgroup(_multipart->comm,
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
      PDM_part_distgroup_to_partgroup(_multipart->comm,
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


      //Fill _part_t structures with temporary arrays
      _pmeshes->parts = (_part_t **) malloc(_pmeshes->tn_part*sizeof(_part_t*));
      for (int ipart = 0; ipart < n_part; ipart++)
      {
        _pmeshes->parts[ipart] = _part_create();

        _pmeshes->parts[ipart]->n_cell = pn_cell[ipart];
        _pmeshes->parts[ipart]->n_face = pn_face[ipart];
        _pmeshes->parts[ipart]->n_vtx    = pn_vtx[ipart];

        _pmeshes->parts[ipart]->vtx      = pvtx_coord[ipart];
        _pmeshes->parts[ipart]->face_vtx_idx = pface_vtx_idx[ipart];
        _pmeshes->parts[ipart]->face_vtx = pface_vtx[ipart];
        _pmeshes->parts[ipart]->cell_face_idx = pcell_face_idx[ipart];
        _pmeshes->parts[ipart]->cell_face = pcell_face[ipart];
        _pmeshes->parts[ipart]->face_cell = pface_cell[ipart];

        _pmeshes->parts[ipart]->cell_ln_to_gn = pcell_ln_to_gn[ipart];
        _pmeshes->parts[ipart]->face_ln_to_gn = pface_ln_to_gn[ipart];
        _pmeshes->parts[ipart]->vtx_ln_to_gn = pvtx_ln_to_gn[ipart];
        _pmeshes->parts[ipart]->face_bound_ln_to_gn = pface_bound_ln_to_gn[ipart];
        _pmeshes->parts[ipart]->face_join_ln_to_gn = pface_join_ln_to_gn[ipart];

        _pmeshes->parts[ipart]->face_bound_idx = pface_bound_idx[ipart];
        _pmeshes->parts[ipart]->face_bound = pface_bound[ipart];

        /* For face_join, the function only returns local id of face in join, we have to
         allocate to set up expected size (4*nb_face_join) */
        _pmeshes->parts[ipart]->face_join_idx = pface_join_idx[ipart];
        int s_face_join = _pmeshes->parts[ipart]->face_join_idx[n_join];
        _pmeshes->parts[ipart]->face_join = (int *) malloc(4*s_face_join*sizeof(int));
        for (int i_face = 0; i_face < s_face_join; i_face++){
          _pmeshes->parts[ipart]->face_join[4*i_face] = pface_join_tmp[ipart][i_face];
        }
        free(pface_join_tmp[ipart]);
      }

      //Use the structure to call renumbering procedures
      PDM_part_renum_cell(_pmeshes->parts, n_part, _pmeshes->renum_cell_method,
                          (void *) _pmeshes->renum_cell_properties);
      PDM_part_renum_face(_pmeshes->parts, n_part, _pmeshes->renum_face_method, NULL);
      //Now genererate comm data -- we need update pface_ln_to_gn which has been modified
      for (int ipart = 0; ipart < n_part; ipart++) {
        pface_ln_to_gn[ipart] = _pmeshes->parts[ipart]->face_ln_to_gn;
      }
      PDM_part_generate_entity_graph_comm(_multipart->comm,
                                          part_distri,
                                          face_distri,
                                          n_part,
                                          pn_face,
                                          (const PDM_g_num_t **) pface_ln_to_gn,
                                          NULL,
                                         &pinternal_face_bound_procidx,
                                         &pinternal_face_bound_partidx,
                                         &pinternal_face_bound);
      //Finally complete parts structure with internal join data
      for (int ipart = 0; ipart < n_part; ipart++)
      {
        _pmeshes->parts[ipart]->face_part_bound_proc_idx = pinternal_face_bound_procidx[ipart];
        _pmeshes->parts[ipart]->face_part_bound_part_idx = pinternal_face_bound_partidx[ipart];
        _pmeshes->parts[ipart]->face_part_bound = pinternal_face_bound[ipart];
      }

      //Free temporary arrays (top-level)
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
      free(pinternal_face_bound_procidx);
      free(pinternal_face_bound_partidx);
      free(pinternal_face_bound);
      free(pcell_ln_to_gn);
      free(pface_ln_to_gn);
      free(pvtx_ln_to_gn);
      free(pface_bound_ln_to_gn);
      free(pface_join_ln_to_gn);

      free(pface_join_tmp);
      free(cell_distri);
      free(face_distri);
      free(vtx_distri);
      free(part_distri);
    }
    // Now rebuild joins over the zones
    _search_matching_joins(_multipart);
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
      int  *n_cell,
      int  *n_face,
      int  *n_face_part_bound,
      int  *n_vtx,
      int  *n_proc,
      int  *n_total_part,
      int  *scell_face,
      int  *sface_vtx,
      int  *sface_bound,
      int  *n_bound_groups,
      int  *sface_join,
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

  *scell_face = _pmeshes.parts[i_part]->cell_face_idx[*n_cell];
  *sface_vtx  = _pmeshes.parts[i_part]->face_vtx_idx[*n_face];

  *n_face_part_bound = _pmeshes.parts[i_part]->face_part_bound_part_idx[*n_total_part];

  *n_bound_groups = _pmeshes.n_bounds;
  *n_join_groups  = _pmeshes.n_joins;
  *sface_bound    = _pmeshes.parts[i_part]->face_bound_idx[*n_bound_groups];
  *sface_join     = _pmeshes.parts[i_part]->face_join_idx[*n_join_groups];
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

}

void
PDM_multipart_part_color_get
(
const int            mpart_id,
const int            i_zone,
const int            i_part,
      int          **cell_color,
      int          **face_color,
      int          **thread_color,
      int          **hyperplane_color
)
{
  _pdm_multipart_t *_multipart = _get_from_id (mpart_id);

  assert(i_zone < _multipart->n_zone && i_part < _multipart->n_part[i_zone]);

  PDM_printf("PDM_multipart_part_color_get: Not implemented\n");
  *cell_color       = NULL;
  *face_color       = NULL;
  *thread_color     = NULL;
  *hyperplane_color = NULL;

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

  PDM_printf("PDM_multipart_time_get: Not implemented\n");
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

  free(_multipart->dmeshes_ids);
  for (int izone = 0; izone < _multipart->n_zone; izone++) {
    free(_multipart->pmeshes[izone].joins_ids);
    for (int ipart = 0; ipart < _multipart->n_part[izone]; ipart++)
      _part_free(_multipart->pmeshes[izone].parts[ipart]);
    free(_multipart->pmeshes[izone].parts);
  }
  free(_multipart->pmeshes);

  PDM_part_renum_method_purge();
  free (_multipart);

  PDM_Handles_handle_free (_multiparts, mpart_id, PDM_FALSE);

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
