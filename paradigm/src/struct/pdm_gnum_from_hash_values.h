#ifndef __PDM_GNUM_FROM_HASH_VALUES_H__
#define __PDM_GNUM_FROM_HASH_VALUES_H__

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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"

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

/*============================================================================
 * Type definitions
 *============================================================================*/
typedef int (*gnum_from_hv_compare)(const void* a, const void* b, void* );
typedef int (*gnum_from_hv_equal  )(const void* a, const void* b, void* );

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/**
 *
 * \brief Build a global numbering structure
 *
 * \param [in]   n_part       Number of local partitions
 * \param [in]   equilibrate  Use algorithm to equilibrate the block treatment (hash value is not a priori equi-reparti)
 * \param [in]   tolerance    Geometric tolerance (used if merge double points is activated)
 * \param [in]   comm         PDM_MPI communicator
 *
 * \return     Identifier
 */
int
PDM_gnum_from_hash_values_create
(
 const int            n_part,
 const PDM_bool_t     equilibrate,
 const size_t         s_data,
 gnum_from_hv_compare fcompare,
 gnum_from_hv_equal   fequal,
 const PDM_MPI_Comm   comm
);

// void
// PROCF (pdm_gnum_from_hash_values_create, PDM_GNUM_FROM_HVALUES_CREATE)
// (
//  const int          *n_part,
//  const int          *equilibrate,
//  const size_t       *s_data,
//  const PDM_MPI_Fint *fcomm,
//        int          *id
// );

/**
 *
 * \brief Set hash values for one partition
 *
 * \param [in]   id           Identifier
 * \param [in]   i_part       Current partition
 * \param [in]   n_elts       Number of elements
 * \param [in]   part_hkey    For each elements the hash value associated
 * \param [in]   part_strid   Stride between each data in part_hdata
 * \param [in]   part_hdata   Partition data which compute the hash value, we need it to setup in a block way
 *
 */

void
PDM_gnum_set_hash_values
(
 const int            id,
 const int            i_part,
 const int            n_elts,
 const size_t        *part_hkeys,
 const int           *part_hstri,
 const unsigned char *part_hdata
);

void
PROCF (pdm_gnum_set_hash_values, PDM_GNUM_SET_FROM_HASH_VALUES)
(
 const int           *id,
 const int           *i_part,
 const int           *n_elts,
 const size_t        *part_hkeys,
 const int           *part_hstri,
 const unsigned char *part_hdata
);

/**
 *
 * \brief Compute
 *
 * \param [in]   id           Identifier
 *
 */

void
PDM_gnum_from_hv_compute
(
 const int id
);


void
PROCF (PDM_gnum_from_hv_compute, PDM_GNUM_FROM_HV_COMPUTE)
(
 const int *id
);


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


PDM_g_num_t *
PDM_gnum_from_hv_get
(
 const int id,
 const int i_part
);

void
PROCF (pdm_gnum_from_hv_get, PDM_GNUM_FROM_HV_GET)
(
 const int *id,
 const int *i_part,
 PDM_g_num_t *gnum
);

/**
 *
 * \brief Dump elapsed an CPU time
 *
 * \param [in]   id           Identifier
 *
 */

void
PDM_gnum_from_hv_dump_times
(
 const int id
);

/**
 *
 * \brief Free
 *
 * \param [in]   id           Identifier
 *
 */

void
PDM_gnum_from_hv_free
(
 const int id,
 const int partial
);

void
PROCF (pdm_gnum_from_hv_free, PDM_GNUM_FROM_HV_FREE)
(
 const int *id,
 const int *partial
);


/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_GNUM_FROM_HASH_VALUES_H__ */
