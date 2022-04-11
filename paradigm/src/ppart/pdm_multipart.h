#ifndef __PDM_MULTIPART_H__
#define __PDM_MULTIPART_H__

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
#include "pdm_dmesh.h"
#include "pdm_dmesh_nodal.h"
#include "pdm_part_mesh_nodal.h"
#include "pdm_domain_interface.h"

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

typedef struct _pdm_multipart_t PDM_multipart_t;

/**
 * \enum PDM_part_size_t
 * \brief Use homogeneous or heterogeneous partition sizes
 */
typedef enum {
  PDM_PART_SIZE_HOMOGENEOUS   = 1,
  PDM_PART_SIZE_HETEROGENEOUS = 2,
} PDM_part_size_t;
/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
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
 * \return     Pointer to a new \ref PDM_multipart_t object
 */

PDM_multipart_t *
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
);


/**
 *
 * \brief Set distributed mesh data for the input zone
 *
 * \param [in]   multipart      Pointer to \ref PDM_multipart_t object
 * \param [in]   zone_id        Global zone id
 * \param [in]   dmesh_id       Id of the distributed mesh structure to use
 */

void PDM_multipart_register_block
(
 PDM_multipart_t   *multipart,
 const int          zone_id,
       PDM_dmesh_t *dmesh
);

/**
 *
 * \brief Set distributed mesh data for the input zone
 *
 * \param [in]   multipart      Pointer to \ref PDM_multipart_t object
 * \param [in]   zone_id        Global zone id
 * \param [in]   dmesh_id       Id of the distributed mesh structure to use
 */

void PDM_multipart_register_dmesh_nodal
(
 PDM_multipart_t         *multipart,
 const int                zone_id,
       PDM_dmesh_nodal_t *dmesh_nodal
);


/**
 *
 * \brief Set connecting data between all the zones
 *
 * \param [in]   multipart         Pointer to \ref PDM_multipart_t object
 * \param [in]   n_total_joins     Total number of interfaces
 * \param [in]   join_to_opposite  For each global join id, give the global id
 *                                   of the opposite join (size = n_total_joins)
 *
 * \note Join global id numbering must start at 0 and be continuous.
 */

void PDM_multipart_register_joins
(
 PDM_multipart_t *multipart,
 const int        n_total_joins,
 const int       *join_to_opposite
);

/**
 *
 * \brief Set the reordering methods to be used after partitioning
 *
 * \param [in]   multipart             Pointer to \ref PDM_multipart_t object
 * \param [in]   i_zone                Id of zone which parameters apply (or -1 for all zones)
 * \param [in]   renum_cell_method     Choice of renumbering method for cells
 * \param [in]   renum_cell_properties Parameters used by cacheblocking method :
 *                                     [n_cell_per_cache_wanted, is_asynchrone, is_vectorisation,
                                        n_vect_face, split_method]
 * \param [in]   renum_face_method     Choice of renumbering method for faces
 *
 */
void PDM_multipart_set_reordering_options
(
 PDM_multipart_t *multipart,
 const int        i_zone,
 const char      *renum_cell_method,
 const int       *renum_cell_properties,
 const char      *renum_face_method
);


/**
 *
 * \brief Construct the partitioned meshes on every zones
 *
 * \param [in]   multipart             Pointer to \ref PDM_multipart_t object
 */
void
PDM_multipart_run_ppart
(
 PDM_multipart_t *multipart
);


/**
 * \brief ???
 *
 * \param [in]  multipart             Pointer to \ref PDM_multipart_t object
 * \param [in]  i_zone                Id of zone
 * \param [out] pmesh_nodal           ?
 *
 */
void
PDM_multipart_get_part_mesh_nodal
(
       PDM_multipart_t        *multipart,
 const int                     i_zone,
       PDM_part_mesh_nodal_t **pmesh_nodal
);

/**
 * \brief Set number of element in the block entity
 *
 * \param [in]  multipart             Pointer to \ref PDM_multipart_t object
 * \param [in]  i_zone                Id of zone
 * \param [in]  entity_type           Type of entity (can be cell/face/edge/vtx)
 * \param [in]  dn_entity             Distributed number of entity in current process
 *
 */
void
PDM_multipart_dn_entity_set
(
       PDM_multipart_t     *multipart,
 const int                  i_zone,
       PDM_mesh_entities_t  entity_type,
       int                  dn_entity
);

/**
 * \brief Set number connectivity for current block
 *
 * \param [in]  multipart             Pointer to \ref PDM_multipart_t object
 * \param [in]  i_zone                Id of zone
 * \param [in]  connectivity_type     Type of connectivity
 * \param [in]  connect               connectivity (size = connect_idx[dn_entity] )
 * \param [in]  connect_idx           Index of connectivity or NULL if face_cell for example  (size = dn_entity )
 *
 */
void
PDM_multipart_dconnectivity_set
(
       PDM_multipart_t         *multipart,
 const int                      i_zone,
       PDM_connectivity_type_t  connectivity_type,
       PDM_g_num_t             *dconnect,
       int                     *dconnect_idx
);

/**
 * \brief Set group connectivity by kind
 *
 * \param [in]  multipart             Pointer to \ref PDM_multipart_t object
 * \param [in]  i_zone                Id of zone
 * \param [in]  bound_type            Type of bound
 * \param [in]  connect               connectivity (size = connect_idx[dn_entity] )
 * \param [in]  connect_idx           Index of connectivity or NULL if face_cell for example  (size = dn_entity )
 *
 */
void
PDM_multipart_dgroup_set
(
       PDM_multipart_t          *multipart,
 const int                       i_zone,
       PDM_bound_type_t          bound_type,
       PDM_g_num_t              *dconnect,
       int                      *dconnect_idx
);


/**
 * \brief Set group connectivity by kind
 *
 * \param [in]  multipart             Pointer to \ref PDM_multipart_t object
 * \param [in]  i_zone                Id of zone
 * \param [in]  dvtx_coord            Mesh coordinates (size = 3 * dn_vtx)
 */
void
PDM_multipart_dvtx_coord_set
(
       PDM_multipart_t *multipart,
 const int              i_zone,
 const double          *dvtx_coord
);


void
PDM_multipart_domain_interface_shared_set
(
  PDM_multipart_t        *multipart,
  PDM_domain_interface_t *ditrf
);


/**
 *
 * \brief Construct the partitioned meshes on every zones
 *
 * \param [in]   multipart      Pointer to \ref PDM_multipart_t object
 */
// void
// PDM_multipart_vtx_graph_comm_compute
// (
//  const int id
// );

/**
 *
 * \brief Returns the dimensions of a given partition
 */
void
PDM_multipart_part_dim_get
(
PDM_multipart_t *multipart,
const int        i_zone,
const int        i_part,
      int       *n_section,
      int      **n_elt,
      int       *n_cell,
      int       *n_face,
      int       *n_face_part_bound,
      int       *n_vtx,
      int       *n_proc,
      int       *n_total_part,
      int       *s_cell_face,
      int       *s_face_vtx,
      int       *s_face_bound,
      int       *n_bound_groups,
      int       *s_face_join,
      int       *n_join_groups
);

/**
 *
 * \brief Returns the dimensions of a given partition
 */
void
PDM_multipart_part_graph_comm_vtx_dim_get
(
 PDM_multipart_t *multipart,
 const int        i_zone,
 const int        i_part,
       int       *n_vtx_part_bound
);

/**
 *
 * \brief Returns the data arrays of a given partition
 */
void
PDM_multipart_part_val_get
(
PDM_multipart_t     *multipart,
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
);

/**
 *
 * \brief Returns the data arrays of a given partition
 */
int
PDM_multipart_part_connectivity_get
(
PDM_multipart_t                *multipart,
const int                       i_zone,
const int                       i_part,
      PDM_connectivity_type_t   connectivity_type,
      int                     **connect,
      int                     **connect_idx,
      PDM_ownership_t           ownership
);

int
PDM_multipart_part_ln_to_gn_get
(
PDM_multipart_t            *multipart,
const int                   i_zone,
const int                   i_part,
      PDM_mesh_entities_t   entity_type,
      PDM_g_num_t         **entity_ln_to_gn,
      PDM_ownership_t       ownership
);

int
PDM_multipart_partition_color_get
(
PDM_multipart_t            *multipart,
const int                   i_zone,
const int                   i_part,
      PDM_mesh_entities_t   entity_type,
      int                 **entity_color,
      PDM_ownership_t       ownership
);

void
PDM_multipart_part_graph_comm_vtx_data_get
(
PDM_multipart_t     *multipart,
const int            i_zone,
const int            i_part,
      int          **vtx_part_bound_proc_idx,
      int          **vtx_part_bound_part_idx,
      int          **vtx_part_bound
);

void
PDM_multipart_part_color_get
(
PDM_multipart_t     *multipart,
const int            i_zone,
const int            i_part,
      int          **cell_color,
      int          **face_color,
      int          **face_hp_color,
      int          **thread_color,
      int          **hyperplane_color
);

void
PDM_multipart_part_ghost_infomation_get
(
PDM_multipart_t  *multipart,
const int         i_zone,
const int         i_part,
      int       **vtx_ghost_information
);


/**
 *
 * \brief Return times for a given zone
 * (NOT IMPLEMENTED)
 *
 * \param [in]   multipart      Pointer to \ref PDM_multipart_t object
 * \param [in]   i_zone         Id of current zone
 * \param [out]  elapsed        Elapsed time
 * \param [out]  cpu            CPU time
 * \param [out]  cpu_user       User CPU time
 * \param [out]  cpu_sys        System CPU time
 *
 */

void
PDM_multipart_time_get
(
 PDM_multipart_t *multipart,
 const int        i_zone,
 double         **elapsed,
 double         **cpu,
 double         **cpu_user,
 double         **cpu_sys
);


/**
 *
 * \brief Free the structure
 *
 * \param [in]   multipart      Pointer to \ref PDM_multipart_t object
 */

void
PDM_multipart_free
(
 PDM_multipart_t *multipart
);


int
PDM_multipart_part_vtx_coord_get
(
PDM_multipart_t                *multipart,
const int                       i_zone,
const int                       i_part,
      double                  **vtx_coord,
      PDM_ownership_t           ownership
);


/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_MULTIPART_H__ */
