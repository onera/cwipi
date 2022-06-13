#ifndef __PDM_GNUM_LOCATION_H__
#define __PDM_GNUM_LOCATION_H__

/*
  This file is part of the ParaDiGM library.

  Copyright (C) 2019       ONERA

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
 * Look for the location of a global numbering element. The location has three
 * propertie : process, partition, number of element in this partition.
 * A global numbering can be located in multiple partitions
 *============================================================================*/

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

typedef struct _pdm_gnum_location_t PDM_gnum_location_t;

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/**
 *
 * \brief Build a global numbering location structure
 *
 * \param [in]   n_part_in      Number of local partitions for elements
 * \param [in]   n_part_out     Number of local partitions for requested locations
 * \param [in]   comm           PDM_MPI communicator
 * \param [in]   owner          Owner
 *
 * \return     Pointer to \ref PDM_gnum_locaion object
 */

PDM_gnum_location_t*
PDM_gnum_location_create
(
 const int             n_part_in,
 const int             n_part_out,
 const PDM_MPI_Comm    comm,
 const PDM_ownership_t owner
);


/**
 *
 * \brief Set global numbering
 *
 * \param [in]   gnum_loc    Pointer to \ref PDM_gnum_locaion object
 * \param [in]   i_part_in   Current partition
 * \param [in]   n_elts_in   Number of elements
 * \param [in]   gnum_in     Global numbering
 *
 */

void
PDM_gnum_location_elements_set
(
       PDM_gnum_location_t *gnum_loc,
 const int                  i_part_in,
 const int                  n_elts_in,
 const PDM_g_num_t         *gnum_in
);


/**
 *
 * \brief Set requested elements
 *
 * \param [in]   gnum_loc     Pointer to \ref PDM_gnum_locaion object
 * \param [in]   i_part_out   Current partition
 * \param [in]   n_elts_out   Number of elements
 * \param [in]   gnum_out     Global numbering
 *
 */

void
PDM_gnum_location_requested_elements_set
(
       PDM_gnum_location_t *gnum_loc,
 const int                  i_part_out,
 const int                  n_elts_out,
 const PDM_g_num_t         *gnum_out
);


/**
 *
 * \brief Compute the location (processus, partittion, local number in the partition)
 *
 * \param [in]   gnum_loc     Pointer to \ref PDM_gnum_locaion object
 *
 */

void
PDM_gnum_location_compute
(
 PDM_gnum_location_t *gnum_loc
);


/**
 *
 * \brief Get localtion
 *
 * \param [in]    gnum_loc       Pointer to \ref PDM_gnum_locaion object
 * \param [in]    i_part_out     Current partition
 * \param [out]   location_idx   Index in the location arrays (size = 3 * \ref n_elts + 1)
 * \param [out]   location       Locations of each element
 *                                (Three informations : process, partition, element)
 *
 */

void
PDM_gnum_location_get
(
       PDM_gnum_location_t  *gnum_loc,
 const int                   i_part_out,
       int                 **location_idx,
       int                 **location
);


/**
 *
 * \brief Free
 *
 * \param [in]   gnum_loc      Pointer to \ref PDM_gnum_locaion object
 *
 */

void
PDM_gnum_location_free
(
  PDM_gnum_location_t *gnum_loc
);


/**
 *
 * \brief Get the number of requested elements in a given partition
 *
 * \param [in]  gnum_loc      Pointer to \ref PDM_gnum_locaion object
 * \param [in]  i_part_out    Id of current partition
 *
 * \return  Number of requested elements in the current partition
 *
 */

int
PDM_gnum_location_n_requested_elt_get
(
       PDM_gnum_location_t *gnum_loc,
 const int                  i_part_out
 );

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_GNUM_H__ */
