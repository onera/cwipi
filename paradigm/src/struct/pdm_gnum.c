/*============================================================================
 * Main structure for an I/O numbering scheme associated with mesh entities
 * (such as cells, faces, and vertices);
 *
 * In parallel mode, such a scheme is important so as to redistribute
 * locally numbered entities on n processes to files written by p
 * processes, with p <= n.
 *
 * Only the case where p = 1 is presently implemented, so the numbering
 * scheme is simply based on entity's global labels.
 *
 * For p > 1, it would probably be necessary to extend the numbering
 * schemes so as to account for the fact that a given entity may have
 * a main index on its main associated domain, but may be present
 * as a ghost entity with another index on neighboring domains.
 *============================================================================*/

/*
  This file is part of the ParaDiGM library, intended to provide
  finite volume mesh and associated fields I/O and manipulation services.

  Copyright (C) 2017       ONERA
  Copyright (C) 2004-2009  EDF

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
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_config.h"
#include "pdm_morton.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_gnum.h"

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

/*----------------------------------------------------------------------------
 * Structure defining an I/O numbering scheme
 *----------------------------------------------------------------------------*/

/**
 * \struct _pdm_gnum_t
 * \brief  Define a global numberring 
 * 
 */

typedef struct  {

  int          n_part;      /*!< Number of partitions */
  int          dim;         /*!< Spatial dimension */
  PDM_MPI_Comm comm;        /*!< MPI communicator */
  PDM_g_num_t  n_g_elt;     /*!< Global number of elements */
  int          *n_elts;     /*!< Number of elements in partitions */
  PDM_g_num_t **g_nums;     /*!< Global number of elements */ 
  double      **coords;     /*!< Coordinates of elements */

} _pdm_gnum_t;

/*============================================================================
 * Global variable
 *============================================================================*/

static _pdm_gnum_t **_gnums   = NULL;
static int         _l_gnums   = 10;

/*=============================================================================
 * Private function definitions
 *============================================================================*/

/**
 *
 * \brief Return ppart object from it identifier
 *
 * \param [in]   ppartId        ppart identifier
 *
 */

static _pdm_gnum_t *
_get_from_id
(
 int  id
)
{
  if (id >= _l_gnums) {
    printf("PDM_gnum error : Bad identifier\n");
    exit(1);
  }
    
  if (_gnums[id] == NULL) {
    printf("PDM_gnum error : Bad identifier\n");
    exit(1);
  }

  return _gnums[id];
}


/**
 * \brief Use bubble sort on an expectedly short sequence of coordinates
 * to ensure lexicographical ordering.
 *
 *  \param[in]      dim        <-- spatial dimension
 *  \param[in]      start_id   <-- start id in array
 *  \param[in]      end_id     <-- past-the-end id in array
 *  \param[in]      coords     <-- pointer to entity coordinates (interlaced)
 *  \param[in, out] order      <-> ordering array base on Morton encoding, or
 *                                 lexicographical coordinate ordering for ties
 */

/*
  This file is part of the "Finite Volume Mesh" library, intended to provide
  finite volume mesh and associated fields I/O and manipulation services.

  Copyright (C) 2004-2009  EDF

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

inline static void
_reorder_coords_lexicographic
(
int                dim,
size_t             start_id,
size_t             end_id,
const double  coords[],
PDM_l_num_t         order[]
)
{
  size_t  i;
  bool g_swap;

  do {

    g_swap = false;

    for (i = start_id + 1; i < end_id; i++) {

      size_t j_prev = order[i-1], j = order[i];
      bool l_swap = false;

      if (dim == 3) {
        if (coords[j_prev*3] < coords[j*3])
          continue;
        else if (coords[j_prev*3] > coords[j*3])
          l_swap = true;
        else if (coords[j_prev*3 + 1] < coords[j*3 + 1])
          continue;
        else if (   coords[j_prev*3 + 1] > coords[j*3 + 1]
                 || coords[j_prev*3 + 2] > coords[j*3 + 2])
          l_swap = true;
      }
      else if (dim == 2) {
        if (coords[j_prev*2] < coords[j*2 + 1])
          continue;
        else if (   coords[j_prev*2]     > coords[j*2]
                 || coords[j_prev*2 + 1] > coords[j*2 + 1])
          l_swap = true;
      }
      else { /* if (dim == 1) */
        if (coords[j_prev] > coords[j])
          l_swap = true;
      }

      if (l_swap) {
        PDM_l_num_t o_save = order[i-1];
        order[i-1] = order[i];
        order[i] = o_save;
        g_swap = true;
      }
    }

  } while (g_swap);
}


/**
 * 
 * \brief Creation of an I/O numbering structure based on coordinates.
 *
 * The ordering is based on a Morton code, and it is expected that
 * entities are unique (i.e. not duplicated on 2 or more ranks).
 * In the case that 2 entities have a same Morton code, their global
 * number will be different, but their order is undetermined.
 *
 * parameters:
 *  \param [in]      dim        spatial dimension
 *  \param [in]      n_entities number of entities considered
 *  \param [in]      coords     pointer to entity coordinates (interlaced)
 *  \param [in]      m_code     Morton code associated with each entity
 *  \param [in, out] order      ordering array base on Morton encoding, or
 *                              lexicographical coordinate ordering for ties
 *
 */

/*
  This file is part of the "Finite Volume Mesh" library, intended to provide
  finite volume mesh and associated fields I/O and manipulation services.

  Copyright (C) 2004-2009  EDF

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

static void
_check_morton_ordering
(
int                      dim,
size_t                   n_entities,
const double        coords[],
const PDM_morton_code_t  m_code[],
PDM_l_num_t               order[]
)
{
  size_t  i_prev = 0, i = 1;

  if (n_entities == 0)
    return;

  /* Check ordering; if two entities have the same Morton codes,
     use lexicographical coordinates ordering to ensure the
     final order is deterministic. */

  for (i = 1; i < n_entities; i++) {

    size_t j_prev = order[i_prev], j = order[i];

    if (   m_code[j_prev].X[0] != m_code[j].X[0]
        || m_code[j_prev].X[1] != m_code[j].X[1]
        || m_code[j_prev].X[2] != m_code[j].X[2]) {

      /* If successive values have the same Morton code,
         order them lexicographically */
      if (i_prev < i - 1)
        _reorder_coords_lexicographic(dim, i_prev, i-1, coords, order);

    }
    i_prev = i;
  }

  if (i_prev < n_entities - 1)
    _reorder_coords_lexicographic(dim, i_prev, n_entities - 1, coords, order);
}


/*=============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief Build a global numbering 
 *
 * \param [in]   dim          Spatial dimension 
 * \param [in]   n_part       Number of local partitions 
 * \param [in]   comm         PDM_MPI communicator
 *
 * \return     Identifier    
 */

int
PDM_gnum_create
(
 const int dim,
 const int n_part,
 const PDM_MPI_Comm comm
)
{
  int id;
  
  /*
   * Search a ppart free id
   */

  if (_gnums == NULL) {
    _gnums = (_pdm_gnum_t **) malloc(_l_gnums * sizeof(_pdm_gnum_t *));
    for (int i = 0; i < _l_gnums; i++)
      _gnums[i] = NULL;
  }

  _pdm_gnum_t * ppart = NULL;
  for (int i = 0; i < _l_gnums; i++) {
    if (_gnums[i] == NULL) {
      _gnums[i] = (_pdm_gnum_t *) malloc(sizeof(_pdm_gnum_t));
      ppart = _gnums[i];
      id = i;
      break;
    }
  }

  if (ppart == NULL) {
    int l_gnums_old = _l_gnums;
    _l_gnums = 2 * _l_gnums;

    _gnums = (_pdm_gnum_t **) realloc(_gnums, _l_gnums * sizeof(_pdm_gnum_t *));
    for (int i = l_gnums_old; i < _l_gnums; i++)
      _gnums[i] = NULL;

    _gnums[l_gnums_old] = (_pdm_gnum_t *) malloc(sizeof(_pdm_gnum_t));
    ppart = _gnums[l_gnums_old];
    id = l_gnums_old;
  }

  _pdm_gnum_t *_gnum = _gnums[id];

  _gnum->n_part    = n_part;
  _gnum->dim       = dim;  
  _gnum->comm      = comm;
  _gnum->n_g_elt   = -1;
  _gnum->g_nums = (PDM_g_num_t **) malloc (sizeof(PDM_g_num_t * ) * n_part); 
  _gnum->coords = (double **) malloc (sizeof(double * ) * n_part);
  _gnum->n_elts = (int *) malloc (sizeof(int) * n_part);
  
  
  for (int i = 0; i < n_part; i++) {
    _gnum->g_nums[i] = NULL;
    _gnum->coords[i] = NULL;
  }

  return id;

}

void
PROCF (pdm_gnum_create, PDM_GNUM_CREATE)
(
 const int *dim,
 const int *n_part,
 const PDM_MPI_Fint *fcomm,
       int *id
)
{
  const PDM_MPI_Comm c_comm = PDM_MPI_Comm_f2c (*fcomm);

  *id = PDM_gnum_create (*dim, *n_part, c_comm);
}


/**
 *
 * \brief Set from coordinates
 *
 * The ordering is based on a Morton code, and it is expected that
 * entities are unique (i.e. not duplicated on 2 or more ranks).
 * In the case that 2 entities have a same Morton code, their global
 * number will be determined by lexicographical ordering of coordinates.
 *
 * \param [in]   id           Identifier
 * \param [in]   i_part       Current partition
 * \param [in]   n_elts       Number of elements
 * \param [in]   coords       Coordinates (size = 3 * \ref n_elts)
 *
 */

void
PDM_gnum_set_from_coords
(
 const int id,
 const int i_part,
 const int n_elts,
 const double *coords
)
{
  _pdm_gnum_t *_gnum = _get_from_id (id);
  
  _gnum->coords[i_part] = (double *) coords;
  _gnum->n_elts[i_part] = n_elts;
  _gnum->g_nums[i_part] = NULL;
  
}     

void
PROCF (pdm_gnum_set_from_coords, PDM_GNUM_SET_FROM_COORDS)
(
 const int *id,
 const int *i_part,
 const int *n_elts,
 const double *coords
)
{
  PDM_gnum_set_from_coords (*id, *i_part, *n_elts, coords);
}


/**
 *
 * \brief Compute
 *
 * \param [in]   id           Identifier
 *
 */

/*
  This file is part of the "Finite Volume Mesh" library, intended to provide
  finite volume mesh and associated fields I/O and manipulation services.

  Copyright (C) 2004-2009  EDF

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

void
PDM_gnum_compute
(
 const int id
)
{
  _pdm_gnum_t *_gnum = _get_from_id (id);

  size_t  i;
  double extents[6];
  PDM_l_num_t  *order = NULL;
  PDM_morton_code_t *m_code = NULL;

  PDM_MPI_Comm comm = _gnum->comm;

  const int level = sizeof(PDM_morton_int_t)*8 - 1;
  
  int n_ranks;
  PDM_MPI_Comm_size (comm, &n_ranks);
  
  int n_entities = 0;
  
  for (int ipart = 0; ipart < _gnum->n_part; ipart++) {
    n_entities += _gnum->n_elts[ipart];
  }

  double *coords = malloc (sizeof(double) * _gnum->dim * n_entities);
  int k = 0;
  
  for (int ipart = 0; ipart < _gnum->n_part; ipart++) {
    for (int j = 0; j < _gnum->dim * _gnum->n_elts[ipart]; j++) {
      coords[k++] = _gnum->coords[ipart][j];
    }
  }
  
  /* Build Morton encoding and order it */

  PDM_morton_get_coord_extents(_gnum->dim, n_entities, coords, extents, comm);

  m_code = malloc (n_entities * sizeof(PDM_morton_code_t));
  order = malloc (n_entities * sizeof(PDM_l_num_t));

  PDM_morton_encode_coords(_gnum->dim, level, extents, n_entities, coords, m_code);

  PDM_morton_local_order(n_entities, m_code, order);

  if (n_ranks > 1) {

    int rank_id;
    PDM_l_num_t j, shift;

    size_t n_block_ents = 0;
    PDM_g_num_t current_global_num = 0, global_num_shift = 0;

    int *c_rank = NULL;
    int *send_count = NULL, *send_shift = NULL;
    int *recv_count = NULL, *recv_shift = NULL;
    double *send_coords = NULL, *recv_coords = NULL;
    PDM_l_num_t *weight = NULL;
    PDM_g_num_t *block_global_num = NULL, *part_global_num = NULL;
    PDM_morton_code_t *morton_index = NULL;

    weight = malloc (n_entities * sizeof(PDM_l_num_t));
    morton_index = malloc ((n_ranks + 1) * sizeof(PDM_morton_code_t));

    for (i = 0; i < n_entities; i++)
      weight[i] = 1;

    PDM_morton_build_rank_index(_gnum->dim,
                                level,
                                n_entities,
                                m_code,
                                weight,
                                order,
                                morton_index,
                                comm);

    free(order);
    free(weight);
    c_rank = malloc (n_entities * sizeof(int));

    for (i = 0; i < n_entities; i++) {
      size_t _c_rank = PDM_morton_quantile_search((size_t) n_ranks,
                                                   m_code[i],
                                                   morton_index);
      c_rank[i] = (int) _c_rank; 
    }
      

    free(morton_index);
    free(m_code);

    /* Build send_buf, send_count and send_shift
       to build a rank to coords indexed list */

    send_count = malloc (n_ranks * sizeof (int));
    recv_count = malloc (n_ranks * sizeof (int));
    send_shift = malloc ((n_ranks + 1) * sizeof (int));
    recv_shift = malloc ((n_ranks + 1) * sizeof (int));

    for (rank_id = 0; rank_id < n_ranks; rank_id++)
      send_count[rank_id] = 0;

    for (i = 0; i < n_entities; i++)
      send_count[c_rank[i]] += _gnum->dim;

    /* Exchange number of coords to send to each process */

    PDM_MPI_Alltoall(send_count, 1, PDM_MPI_INT, recv_count, 1, PDM_MPI_INT, comm);

    send_shift[0] = 0;
    recv_shift[0] = 0;
    for (rank_id = 0; rank_id < n_ranks; rank_id++) {
      send_shift[rank_id + 1] = send_shift[rank_id] + send_count[rank_id];
      recv_shift[rank_id + 1] = recv_shift[rank_id] + recv_count[rank_id];
    }

    /* Build send and receive buffers */

    send_coords = malloc (send_shift[n_ranks] * sizeof(double));

    for (rank_id = 0; rank_id < n_ranks; rank_id++)
      send_count[rank_id] = 0;

    for (i = 0; i < n_entities; i++) {
      rank_id = c_rank[i];
      shift = send_shift[rank_id] + send_count[rank_id];
      for (j = 0; j < _gnum->dim; j++)
        send_coords[shift + j] = coords[i*_gnum->dim + j];
      send_count[rank_id] += _gnum->dim;
    }

    recv_coords = malloc (recv_shift[n_ranks] * sizeof(double));

    /* Exchange coords between processes */

    PDM_MPI_Alltoallv(send_coords, send_count, send_shift, PDM_MPI_DOUBLE,
                      recv_coords, recv_count, recv_shift, PDM_MPI_DOUBLE,
                      comm);

    free(send_coords);

    /* Now re-build Morton codes on block distribution */

    n_block_ents = recv_shift[n_ranks] / _gnum->dim;

    m_code = malloc (n_block_ents * sizeof(PDM_morton_code_t));
    order = malloc (n_block_ents * sizeof(PDM_l_num_t));

    PDM_morton_encode_coords(_gnum->dim,
                             level,
                             extents,
                             n_block_ents,
                             recv_coords,
                             m_code);

    PDM_morton_local_order((int) n_block_ents, m_code, order);

    /* Check ordering; if two entities have the same Morton codes,
       use lexicographical coordinates ordering to ensure the
       final order is deterministic. */

    _check_morton_ordering(_gnum->dim, n_block_ents, recv_coords, m_code, order);

    /* Determine global order; requires ordering to loop through buffer by
       increasing number (slice blocks associated with each process are
       already sorted, but the whole "gathered" slice is not).
       We build an initial global order based on the initial global numbering,
       such that for each slice, the global number of an entity is equal to
       the cumulative number of sub-entities */

    free(m_code);
    free(recv_coords);
    block_global_num = malloc(n_block_ents * sizeof(PDM_g_num_t));

    for (i = 0; i < n_block_ents; i++)
      block_global_num[order[i]] = (PDM_g_num_t) i + 1;

    free(order);

    current_global_num = (PDM_g_num_t) n_block_ents;

    /* At this stage, block_global_num[] is valid for this process, and
       current_global_num indicates the total number of entities handled
       by this process; we must now shift global numberings on different
       processes by the cumulative total number of entities handled by
       each process */

    PDM_MPI_Scan(&current_global_num, &global_num_shift, 1, PDM__PDM_MPI_G_NUM,
                 PDM_MPI_SUM, comm);
    global_num_shift -= current_global_num;

    for (i = 0; i < n_block_ents; i++)
      block_global_num[i] += global_num_shift;

    /* Return global order to all processors */

    for (rank_id = 0; rank_id < n_ranks; rank_id++) {
      send_count[rank_id] /= _gnum->dim;
      recv_count[rank_id] /= _gnum->dim;
      send_shift[rank_id] /= _gnum->dim;
      recv_shift[rank_id] /= _gnum->dim;
    }

    send_shift[n_ranks] /= _gnum->dim;

    part_global_num = malloc (send_shift[n_ranks] * sizeof(PDM_g_num_t));

    PDM_MPI_Alltoallv(block_global_num, recv_count, recv_shift, PDM__PDM_MPI_G_NUM,
                      part_global_num, send_count, send_shift, PDM__PDM_MPI_G_NUM,
                      comm);

    for (rank_id = 0; rank_id < n_ranks; rank_id++)
      send_count[rank_id] = 0;

    k = 0;
    PDM_g_num_t _max_loc = -1;
    for (int ipart = 0; ipart < _gnum->n_part; ipart++) {
      _gnum->g_nums[ipart] = malloc (sizeof(PDM_g_num_t) * _gnum->n_elts[ipart]);
      for (int j1 = 0; j1 < _gnum->n_elts[ipart]; j1++) {
        rank_id = c_rank[k++];
        shift = send_shift[rank_id] + send_count[rank_id];
        _gnum->g_nums[ipart][j1] = part_global_num[shift];
        _max_loc = PDM_MAX (_max_loc, part_global_num[shift]);
        send_count[rank_id] += 1;
      }  
    }
    
    /* Free memory */

    free(c_rank);

    free(block_global_num);
    free(part_global_num);

    free(send_count);
    free(recv_count);
    free(send_shift);
    free(recv_shift);

    /* Get final maximum global number value */

    PDM_MPI_Allreduce (&_max_loc, 
                       &_gnum->n_g_elt, 
                       1, 
                       PDM__PDM_MPI_G_NUM, 
                       PDM_MPI_MAX, 
                       comm);

  }

  else if (n_ranks == 1) {

    _check_morton_ordering(_gnum->dim, n_entities, coords, m_code, order);

    free(m_code);
    
    PDM_g_num_t *tmp_gnum = malloc (sizeof(PDM_g_num_t) * n_entities);
    
    for (i = 0; i < n_entities; i++)
      tmp_gnum[order[i]] = (PDM_g_num_t) i+1;
    
    k = 0;
    for (int ipart = 0; ipart < _gnum->n_part; ipart++) {
      _gnum->g_nums[ipart] = malloc (sizeof(PDM_g_num_t) * _gnum->n_elts[ipart]);
      for (int j1 = 0; j1 <  _gnum->n_elts[ipart]; j1++) {
        _gnum->g_nums[ipart][j1] = tmp_gnum[k++];
      }  
    }

    free(order);
    free(tmp_gnum);
    
    _gnum->n_g_elt = n_entities;

  }

  free (coords);
}

void
PROCF (pdm_gnum_compute, PDM_GNUM_COMPUTE)
(
 const int *id
)
{
  PDM_gnum_compute (*id);
}


/**
 *
 * \brief Set from coordinates
 *
 * \param [in]   id           Identifier
 * \param [in]   i_part       Current partition
 * \param [in]   n_elts       Number of elements
 * \param [in]   coords       Coordinates (size = 3 * \ref n_elts)
 *
 */

const PDM_g_num_t *
PDM_gnum_get
(
 const int id,
 const int i_part
)
{
  _pdm_gnum_t *_gnum = _get_from_id (id);

  return _gnum->g_nums[i_part];
}

void
PROCF (pdm_gnum_get, PDM_GNUM_GET)
(
 const int *id,
 const int *i_part,
 PDM_g_num_t *gnum
)
{
  _pdm_gnum_t *_gnum = _get_from_id (*id);
  
  const PDM_g_num_t *tmp = PDM_gnum_get (*id, *i_part);
  for (int i = 0; i < _gnum->n_elts[*i_part]; i++) {
    gnum[i] = tmp[i];    
  }
}


/**
 *
 * \brief Free
 *
 * \param [in]   id           Identifier
 *
 */

void
PDM_gnum_free
(
 const int id,
 const int partial       
)
{
  _pdm_gnum_t *_gnum = _get_from_id (id);

  free (_gnum->coords);  
  
  if (partial != 1) {
    for (int i = 0; i < _gnum->n_part; i++) {
      free (_gnum->g_nums[i]);
    }
  }
  
  free (_gnum->g_nums);
  
  free(_gnums[id]);
  _gnums[id] = NULL;

  int toDel = 1;
  for (int i = 0; i < _l_gnums; i++) { 
    if (_gnums[i] != NULL) {
      toDel = 0;
      break;
    }
  }
  if (toDel) {
    free(_gnums);
    _gnums = NULL;
  }

}

void
PROCF (pdm_gnum_free, PDM_GNUM_FREE)
(
 const int *id,
 const int *partial       
)
{
  PDM_gnum_free (*id, *partial);
}


/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
