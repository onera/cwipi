#ifndef __CWIPI_NEWTOOLD_H__
#define __CWIPI_NEWTOOLD_H__
/*
  This file is part of the CWIPI library.

  Copyright (C) 2011-2017  ONERA

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library. If not, see <http://www.gnu.org/licenses/>.
*/


/** \file cwipi_newtoold.h
  * \brief CWIPI API faking CWIPI 0.x API but using CWIPI 1.0
  *
  */

#include <stdio.h>
/**
 * \cond
 */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

#if !defined (__hpux) && !defined (_AIX)
#define PROCF(x, y) x##_
#else
#define PROCF(x, y) x
#endif

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */



typedef enum {

  CWIPI_STATIC_MESH,
  CWIPI_MOBILE_MESH,
  CWIPI_CYCLIC_MESH

} cwipi_mesh_type_t;


//TODO: Look at the precisely the equivalence for MOBILE and CYCLIC MESH
const std::map<cwipi_mesh_type_t,CWP_Dynamic_mesh_t> mesh_type_conv = { {CWIPI_STATIC_MESH, CWP_DYNAMIC_MESH_STATIC},
                                                                  {CWIPI_MOBILE_MESH, CWP_DYNAMIC_MESH_VARIABLE},
                                                                  {CWIPI_CYCLIC_MESH, CWP_DYNAMIC_MESH_VARIABLE}
                                                                };



/*----------------------------------------------------------------------------
 * MPI ranks used for the coupling
 *----------------------------------------------------------------------------*/

typedef enum {

  CWIPI_COUPLING_PARALLEL_WITH_PARTITIONING,
  CWIPI_COUPLING_PARALLEL_WITHOUT_PARTITIONING,
  CWIPI_COUPLING_SEQUENTIAL

} cwipi_coupling_type_t;

const std::map<cwipi_coupling_type_t,CWP_Comm_t> coupling_type_conv = { {CWIPI_COUPLING_PARALLEL_WITH_PARTITIONING   , CWP_COMM_PAR_WITH_PART   },
                                                                        {CWIPI_COUPLING_PARALLEL_WITHOUT_PARTITIONING, CWP_COMM_PAR_WITHOUT_PART},
                                                                        {CWIPI_COUPLING_SEQUENTIAL                   , CWP_COMM_SEQ             }
                                                                      };


/*----------------------------------------------------------------------------
 * Solver type
 *----------------------------------------------------------------------------*/

typedef enum {

  CWIPI_SOLVER_CELL_CENTER,
  CWIPI_SOLVER_CELL_VERTEX

} cwipi_solver_type_t;


/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 *
 * Initialize the cwipi library.
 * Redirect outputs in a file (Standard output with output_listing = NULL or
 * output_logical_unit = -1)
 * Create the current communicator application from 'common_comm'.
 *
 * parameters:
 *   common_comm       <-- Common MPI communicator
 *   application_name  <-- Current application name
 *   application_comm  --> Internal MPI communicator for the current
 *                         application
 *
 * It is a synchronization point between all applications
 *----------------------------------------------------------------------------*/

void cwipi_init_new_to_old
(const MPI_Comm                           global_comm,
 const char                               *code_name,
 MPI_Comm                                 *local_comm
);


/*----------------------------------------------------------------------------
 *
 * Create a coupling object
 *
 * parameters:
 *   coupling_name           <-- Coupling identifier
 *   coupling_type           <-- Coupling type
 *   coupled_application     <-- Coupled application name
 *   entitiesDim             <-- Mesh entities dimension (1, 2 or 3)
 *   tolerance               <-- Geometric tolerance to locate
 *   mesh_type               <-- CWIPI_STATIC_MESH
 *                               CWIPI_MOBILE_MESH (not implemented yet)
 *                               CWIPI_CYCLIC_MESH
 *   solver_type             <-- CWIPI_SOLVER_CELL_CENTER
 *                               CWIPI_SOLVER_CELL_VERTEX
 *   output_frequency        <-- Output frequency
 *   output_format           <-- Output format to visualize exchanged fields
 *                               on the coupled mesh. Choice between :
 *                                 - "EnSight Gold"
 *                                 - "MED_fichier"
 *                                 - "CGNS"
 *   output_format_option    <-- Output options "opt1, opt2,
 *                             text             output text files
 *                             binary              output binary files (default)
 *                             big_endian          force binary files
 *                                                 to big-endian
 *                             discard_polygons    do not output polygons
 *                                                 or related values
 *                             discard_polyhedra   do not output polyhedra
 *                                                 or related values
 *                             divide_polygons     tesselate polygons
 *                                                 with triangles
 *                             divide_polyhedra    tesselate polyhedra
 *                                                 with tetrahedra and pyramids
 *                                                 (adding a vertex near
 *                                                 each polyhedron's center)
 *   nblocations             <-- Number of possible localisations with
 *                               CWIPI_CYCLIC_MESH, optional
 *
 *----------------------------------------------------------------------------*/

void cwipi_create_coupling_new_to_old
( const char  *coupling_name,
  const cwipi_coupling_type_t coupling_type,
  const char  *coupled_application,
  const int    entitiesDim,
  const double tolerance,
  const cwipi_mesh_type_t mesh_type,
  const cwipi_solver_type_t solver_type,
  const int    output_frequency,
  const char  *output_format,
  const char  *output_format_option,
  ...);





#ifdef __cplusplus
}
#endif /* __cplusplus */

/**
 * \endcond
 */

#endif /* __CWIPI_H__ */
