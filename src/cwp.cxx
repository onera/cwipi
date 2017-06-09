
/*
  This file is part of the CWIPI library. 

  Copyright (C) 2011  ONERA

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

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <cstring>
#include <cstdlib>

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

#include <bftc_mem.h>
#include <bftc_printf.h>

/*----------------------------------------------------------------------------
 * FVM library headers
 *----------------------------------------------------------------------------*/
#include <fvmc_parall.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cwp.h"
#include "cwipi_config.h"
#include "factory.hpp"
#include "codeProperties.hxx"
#include "codePropertiesDB.hxx"
#include "codePropertiesDB_i.hxx"
#include "couplingDB.hxx"
#include "couplingDB_i.hxx"
#include "coupling.hxx"
#include "coupling_i.hxx"
#include "commWithPart.hxx"
#include "commWithoutPart.hxx"
#include "commSeq.hxx"
// #include "geometry.hxx"
// #include "location.hxx"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

using namespace std;

/*============================================================================
 * Fortran function header
 *============================================================================*/

/*============================================================================
 * Local macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Global variable
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Output listing File (C printing)
 *----------------------------------------------------------------------------*/

static FILE* _cwipi_output_listing;

/*============================================================================
 * Private function definitions
 *============================================================================*/

/**
 *
 * \brief Intermediate function to write an output in the C file output
 * 
 * \param [in] format     Format (as printf)
 * \param [in] arg_ptr    values to write (as printf)
 *
 * \return                writting status
 */

static int 
_cwipi_print_with_c
( 
 const char     *const format,
 va_list         arg_ptr
)
{
  return vfprintf(_cwipi_output_listing, format, arg_ptr);
}

/**
 *
 * \brief Return Coupling instance from it identifier
 * 
 * \param [in] cpl_id     Coupling identifier
 * \param [in] fct        Function
 *
 * \return                Coupling instance from it identifier
 */

static cwipi::Coupling&
_cpl_get
(
 const char     *cpl_id
 )
{
  cwipi::CouplingDB & couplingDB =
    cwipi::CouplingDB::getInstance();
  
  const string &cpl_name_str = cpl_id;
  return couplingDB.couplingGet(cpl_name_str);
}


/*=============================================================================
 * Public function prototypes - with xml data
 *============================================================================*/

/*============================================================================*
 *                                                                            *
 *                     Public function prototypes                             *
 *                     --------------------------                             *
 *                          without xml                                       *
 *                          -----------                                       *
 *                                                                            *
 *============================================================================*/

/*----------------------------------------------------------------------------*
 * General functions                                                          *
 *----------------------------------------------------------------------------*/

/*============================================================================*
 *                                                                            *
 *                     Public function prototypes                             *
 *                     --------------------------                             *
 *                        without xml data                                    *
 *                        ----------------                                    *
 *                                                                            *
 *============================================================================*/

/*----------------------------------------------------------------------------*
 * General functions                                                          *
 *----------------------------------------------------------------------------*/

/**
 * \brief Initialize CWIPI.
 *
 * This function create the MPI intra communicator for this code from
 * the MPI inter communicator that contains all code process. It is a
 * synchronization point between all codes
 *
 * \param [in]  global_comm       MPI global communicator
 * \param [in]  n_code            Number of codes on the current rank
 * \param [in]  is_coupled_rank   Is current rank used for coupling (size = \ref n_code)
 * \param [in]  code_name         Names of codes on the current rank (size = \ref n_code)
 * \param [in]  time_init         Time init (size = \ref n_code)
 * \param [out] intra_comm        MPI intra communicators of each code
 *
 */

void 
CWP_Init
(
 const MPI_Comm           global_comm,
 const int                n_code,
 const char             **code_names,
 const CWP_Status_t      *is_coupled_rank,
 const double            *time_init,
 MPI_Comm                *intra_comms
)

{
  const int n_param_max_default = 100;
  const int str_size_max_default = 80;
  
  int n_param_max = n_param_max_default;
  int str_size_max = str_size_max_default;
      
  char* pPath;
  pPath = getenv ("CWP_N_PARAM_MAX");
  if (pPath!=NULL) {
    n_param_max = atoi(pPath);
  }

  if (pPath != NULL) {  
   printf ("environment variable 'CWP_N_PARAM_MAX'"
           "is define to: %s\n",pPath);
  }
  pPath = getenv ("CWP_STR_SIZE_MAX");
  if (pPath!=NULL) {
    str_size_max = atoi(pPath);
  }
  if (pPath != NULL) {  
   printf ("environment variable 'CWP_STR_SIZE_MAX'"
           "is define to: %s\n",pPath);
  }
  
  /*
   * Get application properties
   */

  cwipi::CodePropertiesDB & properties =
    cwipi::CodePropertiesDB::getInstance();

  bftc_printf("\ncwipi " CWIPI_VERSION " initializing\n");
  bftc_printf("------------------------\n\n");

  /*
   * Builds application communicator
   */

  properties.init (global_comm,
                   n_code,
                   code_names,
                   is_coupled_rank,
                   n_param_max,
                   str_size_max,
                   intra_comms);

  /*
   * Create default parameters
   */
  
  MPI_Barrier(global_comm);
  
  for (int i = 0; i < n_code; i++) {
    const string &codeNameStr = code_names[i]; 
    properties.ctrlParamAdd <double> (codeNameStr, "time", time_init[i]);
    properties.ctrlParamAdd <int> (codeNameStr, "state", CWP_STATE_IN_PROGRESS);
  }

  MPI_Barrier(global_comm);

  /*
   * Create communication abstract factory 
   */

  cwipi::Factory<cwipi::Communication, CWP_Comm_t> &factoryComm = 
    cwipi::Factory<cwipi::Communication, CWP_Comm_t>::getInstance();

  factoryComm.Register<cwipi::CommWithPart>(CWP_COMM_PAR_WITH_PART);
  factoryComm.Register<cwipi::CommWithoutPart>(CWP_COMM_PAR_WITHOUT_PART);
  factoryComm.Register<cwipi::CommSeq>(CWP_COMM_SEQ);

  /*
   * Create geometry abstract factory 
   */

  // Factory<Geometry, CWP_geom_t> &factoryGeom = 
  //   cwipi::Factory<Geometry, CWP_geom_t>::getInstance();

  // factoryGeom.Register<Location>(CWP_GEOM_LOCATION);
  // factoryGeom.Register<Intersection>(CWP_GEOM_INTERSECTION);
  // factoryGeom.Register<ClosestPoint>(CWP_GEOM_CLOSEST_POINT);

  /*
   * Create support abstract factory 
   */

  // Factory<Support, CWP_support_t> &factorySupport = 
  //   cwipi::Factory<Support, CWP_support_t>::getInstance();

  // factorySupport.Register<Mesh>(CWP_SUPPORT_MESH);
  // factorySupport.Register<PointCloud>(CWP_SUPPORT_POINT_CLOUD);

  /*
   * Create support abstract factory 
   */

  // Factory<Block, CWP_block_t> &factoryBlock = 
  //   cwipi::Factory<Block, CWP_block_t>::getInstance();

  // factoryBlock.Register<BlockNode>(CWP_BLOCK_NODE);
  // factoryBlock.Register<BlockEdge2>(CWP_BLOCK_FACE_EDGE2);
  // factoryBlock.Register<BlockFaceTria3>(CWP_BLOCK_FACE_TRIA3);
  // factoryBlock.Register<BlockFacepoly>(CWP_BLOCK_FACE_POLY);
  // factoryBlock.Register<BlockCellTetra4>(CWP_BLOCK_CELL_TERTRA4);
  // factoryBlock.Register<BlockCellHexa8>(CWP_BLOCK_CELL_HEXA8);
  // factoryBlock.Register<BlockCellPrism6>(CWP_BLOCK_CELL_PRISM6);
  // factoryBlock.Register<BlockCellPyram5>(CWP_BLOCK_CELL_PYRAM5);
  // factoryBlock.Register<BlockCellPoly>(CWP_BLOCK_CELL_POLY);

  MPI_Barrier(global_comm);

}

/*============================================================================*
 *                                                                            *
 *                   Other public function prototypes                         *
 *                   --------------------------------                         *
 *                                                                            *
 *============================================================================*/

/*----------------------------------------------------------------------------*
 * General functions                                                          *
 *----------------------------------------------------------------------------*/

/**
 *
 * \brief CWIPI completion
 *
 * This function finalize CWP
 *
 */

void 
CWP_Finalize
(
)
{
  int flag = 0;
  MPI_Initialized(&flag);

  cwipi::CodePropertiesDB & properties =
    cwipi::CodePropertiesDB::getInstance();

  const MPI_Comm globalComm = properties.globalCommGet();
  printf("CWP_Finalize\n");
  fflush(stdout);
  if (flag != 0) {
    bftc_printf_flush();
    MPI_Barrier(globalComm);
    MPI_Comm oldFVMComm = fvmc_parall_get_mpi_comm();
  }

  properties.kill();

}

/*----------------------------------------------------------------------------*
 * Functions about current application properties                             *
 *----------------------------------------------------------------------------*/

/**
 * \brief Update code state.
 *
 * This function set the code state.
 *
 * \param [in] local_code_name    Local code name
 * \param [in] state              State
 *
 */

void 
CWP_State_update
(
 const char* local_code_name,
 const CWP_State_t state
)
{    
  cwipi::CodePropertiesDB & properties =
    cwipi::CodePropertiesDB::getInstance();
  
  properties.ctrlParamSet<int>(string(local_code_name), "state", state);
}

/**
 * \brief Get current code state.
 *
 * This function return the state of the current code.
 *
 * \param [in]  code_name    Local or distant code name
 *
 * \return state    State
 *
 */

CWP_State_t
CWP_State_get
(
 const char* code_name
)
{
  cwipi::CodePropertiesDB & properties =
    cwipi::CodePropertiesDB::getInstance();

  return static_cast<CWP_State_t> (
      properties.ctrlParamGet<int>(string(code_name), "state"));
}

/**
 * \brief Update application time
 *
 * This function update the application current time.
 *
 * \param [in] local_code_name  Local code name
 * \param [in]  current_time    Current time
 *
 */

void 
CWP_Time_update
(
 const char* local_code_name,
 const double current_time
)
{
  cwipi::CodePropertiesDB & properties =
    cwipi::CodePropertiesDB::getInstance();
  properties.ctrlParamSet<double>(string(local_code_name),"time", current_time);
}

/**
 * \brief Writing output to file.
 *
 * This function set the file for writing output.
 *
 * \param [in] output_file    Output file
 *
 */

void 
CWP_Output_file_set
(
 FILE *output_file
)
{
  _cwipi_output_listing = output_file;
  bftc_printf_proxy_set(_cwipi_print_with_c);
}

/*----------------------------------------------------------------------------*
 * Functions about properties                                                 *
 *----------------------------------------------------------------------------*/

/**
 * \brief Dump application properties.
 *
 * This function dump application properties.
 *
 */

void 
CWP_Properties_dump
(
 void
)
{
  cwipi::CodePropertiesDB & properties =
    cwipi::CodePropertiesDB::getInstance();
  properties.dump();
}

/*----------------------------------------------------------------------------*
 * General functions about coupling                                           *
 *----------------------------------------------------------------------------*/

/**
 * \brief Creating a coupling object.
 *
 * This function creates a coupling object and defines its properties.
 *
  * \param [in]  cpl_id              Coupling identifier
 * \param [in]  local_code_name     Local code name
 * \param [in]  coupled_code_name   Distant or local coupled code name
 * \param [in]  comm_type           Communication type
 * \param [in]  geom_algo           Geometric algorithm
 * \param [in]  support_type        Support type
 * \param [in]  n_part              Number of interface partition 
 * \param [in]  moving_status       Support moving status
 * \param [in]  recv_freq_type      Type of receiving frequency
 *
 */

void 
CWP_Cpl_create
(
 const char               *cpl_id,
 const char               *local_code_name,
 const char               *coupled_code_name,
 const CWP_Comm_t          comm_type, 
 const CWP_Geom_t          geom_algo,
 const CWP_Support_t       support_type,
 const int                 n_part,
 const CWP_Displacement_t  displacement,   
 const CWP_Freq_t          recv_freq_type 
)
{
  cwipi::CouplingDB & couplingDB =
    cwipi::CouplingDB::getInstance();

  cwipi::CodePropertiesDB & properties =
    cwipi::CodePropertiesDB::getInstance();

  const string &coupling_name_str = cpl_id;
  const string &coupled_application_str = coupled_code_name;
  const string &local_application_str = local_code_name;

  couplingDB.couplingCreate(coupling_name_str,
                            comm_type,
                            properties.codePropertiesGet(local_application_str),
                            properties.codePropertiesGet(coupled_application_str),
                            geom_algo,
                            support_type,
                            n_part,
                            displacement,
                            recv_freq_type);
}

/**
 *
 * \brief Removing a coupling object
 *
 * This function delete a coupling abject
 * 
 * \param [in] cpl_id     Coupling identifier
 *
 */

void 
CWP_Cpl_del
(const char *cpl_id
)
{
  cwipi::CouplingDB & couplingDB =
    cwipi::CouplingDB::getInstance();

  const string &cpl_id_str = cpl_id;

  couplingDB.couplingDel(cpl_id_str);
}

/**
 * \brief data exchange  <b>(Not implemented yet)</b> 
 *
 * This function exchanges for each coupling depending on exchange frequency
 * 
 * \param [in] cpl_id     Coupling identifier
 */

// void
// CWP_Exch
// (
//  const char *cpl_id
// )
// {
//   //TODO: Voir comment enchainer les appels, voir comment prendre en compte l'interpolation temporelle 
//   cwipi::Coupling& cpl = _cpl_get(cpl_id);
//   cpl.exchange();
// }

/**
 *
 * \brief Return the number of uncomputed targets
 * 
 * \param [in] cpl_id     Coupling identifier
 *
 * \return                Number of uncomputed targets
 */

// int 
// CWP_n_uncomputed_tgts_get
// (
//  const char *cpl_id
// )
// {
//   cwipi::Coupling& cpl = _cpl_get(cpl_id);

//   return cpl.nUncomputedTargetsGet();
// }

/**
 *
 * \brief Return uncomputed targets
 * 
 * \param [in] cpl_id     Coupling identifier
 *
 * \return                Uncomputed targets
 */

// const int *
// CWP_uncomputed_tgts_get
// (
//  const char *cpl_id
// )
// {
//   cwipi::Coupling& cpl = _cpl_get(cpl_id);
//   return cpl.uncomputedTargetsGet();
// }

/**
 *
 * \brief Return the number of computed targets
 * 
 * \param [in] cpl_id     Coupling identifier
 *
 * \return                Number of computed targets
 */

// int 
// CWP_n_computed_tgts_get
// (
//  const char *cpl_id
// )
// {
//   cwipi::Coupling& cpl = _cpl_get(cpl_id);
//   return cpl.nComputedTargetsGet();
// }

/**
 *
 * \brief Return computed targets
 * 
 * \param [in] cpl_id     Coupling identifier
 *
 * \return                Computed targets
 */

// const int *
// CWP_computed_tgts_get
// (
//  const char *cpl_id
// )
// {
//   cwipi::Coupling& cpl = _cpl_get(cpl_id);
//   return cpl.getLocatedPoint();
// }

/**
 * \brief Return distance from each target to the geometric interface                 
 *
 * \param [in]  cpl_id   Coupling identifier
 *
 * \return               Distance
 *
 */

// const double *
// CWP_computed_tgts_dist_to_geom_get
// (
//  const char *cpl_id
// )
// {
//   TODO
// } 


/*----------------------------------------------------------------------------*
 * Functions about exchange frequency                                         *
 *----------------------------------------------------------------------------*/

/**
 * \brief Setting receiving frequency.
 *
 * This function set receiving frequency. It must be used when
 * the type of receiving frequency is \ref FREQ_RELATED_N_TIME_STEP
 *
 * \param [in]  cpl_id     Coupling identifier
 * \param [in]  n_step     Frequency in steps number
 *
 */

// void 
// CWP_recv_freq_set
// (
//  const char                 *cpl_id,
//  const int                   n_step
// )
// {
//   cwipi::Coupling& cpl = _cpl_get(cpl_id);
//   cpl.recvFreqSet(n_step);
// }

/**
 * \brief Setting the next receiving time.
 *
 * This function set the next receiving time. It must be used when
 * the type of receiving frequency is \ref FREQ_ASYNCHRONOUS
 *
 * \param [in]  cpl_id        Coupling identifier
 * \param [in]  next_time     Next receiving time
 *
 */

// void 
// CWP_next_recv_time_set
// (
//  const char                 *cpl_id,
//  const double                next_time
// )
// {
//   cwipi::Coupling& cpl = _cpl_get(cpl_id);
//   cpl.recvNextTimeSet(next_time);
// }

/*----------------------------------------------------------------------------*
 * Functions about geometry                                                   *
 *----------------------------------------------------------------------------*/

/**
 * \brief Computation geometry                                  
 *
 * This function compute geometry 
 *
 * \param [in]  cpl_id              Coupling identifier
 * \param [in]  storage_activation  Storage activation of geometric results
 * \param [out] n_uncomputed_tgt    Number of uncomputed target
 * \param [out] storage_id          Storage identifier              
 *
 */

//void 
//CWP_geom_compute
//(
// const char     *cpl_id,
// CWP_status_t  storage_activation,
// int            *n_uncomputed_tgt,
// int            *storage_id
//)
// {
//   cwipi::Coupling& cpl = _cpl_get(cpl_id);
//
//   cpl.geomCompute(n_uncomputed_tgt);
// }

/**
 * \brief Stored geometric results activation                 
 *
 * \param [in]  cpl_id      Coupling identifier
 * \param [in] storage_id   Storage identifier              
 *
 */

// void 
// CWP_geom_update
// (
//  const char     *cpl_id,
//  const int       storage_id
// )
// {
// TODO
// }

/**
 * \brief Setting geometry properties
 *
 * This function set the geometric algorithm properties.
 *
 * \param [in]  cpl_id        Coupling identifier
 * \param [in]  fmt           Format with the syntax : "prop1, prop2, ..."
 * \param       ...           Values of each properties
 *
 */

// void 
// CWP_geom_properties_set
// (
//  const char     *cpl_id,
//  const char     *fmt,
//                   ...
// )
// {
//   cwipi::Coupling& cpl = _cpl_get(cpl_id);

//   va_list pa;
//   va_start(pa, fmt);
//   cpl.geomPropertiesSet(fmt, &pa);
//   va_end(pa);
// }

/*----------------------------------------------------------------------------*
 * Functions about visualization                                              *
 *----------------------------------------------------------------------------*/

/**
 * \brief Enable visualization output
 *
 * This function enable visualization output.
 *
 * \param [in]  cpl_id           Coupling identifier
 * \param [in]  freq             Output frequency
 * \param [in]  format           Output format to visualize exchanged fields
 *                               on the coupled mesh. Choice between :
 *                               - "EnSight Gold"
 *                               - "MED_ficher"
 *                               - "CGNS"
 *                               .
 * \param [in]  format_option   Output options "opt1, opt2, ..."
 *                               - text               output text files
 *                               - binary             output binary files (default)
 *                               - big_endian         force binary files
 *                                                    to big-endian
 *                               - discard_polygons   do not output polygons
 *                                                    or related values
 *                               - discard_polyhedra  do not output polyhedra
 *                                                    or related values
 *                               - divide_polygons    tesselate polygons
 *                                                    with triangles
 *                               - divide_polyhedra   tesselate polyhedra
 *                                                    with tetrahedra and pyramids
 *                                                    (adding a vertex near
 *                                                    each polyhedron's center)
 *                               .
 *
 */

// void 
// CWP_visu_set
// (
//  const char                 *cpl_id,
//  const int                   freq,
//  const char                 *format,
//  const char                 *format_option
// )
// {
//   cwipi::Coupling& cpl = _cpl_get(cpl_id);
//   cpl.visuSet(freq, format, format_option);
// }

/*----------------------------------------------------------------------------*
 * Functions about User target points                                         *
 *----------------------------------------------------------------------------*/

/**
 * \brief Setting user target points
 *
 * This function must be called if the nature of receiving fields 
 * is \ref CWP_FIELD_NATURE_USER
 *
 * \param [in]  cpl_id  Coupling identifier
 * \param [in]  n_pts   Number of points
 * \param [in]  coord   Coordinates (size = 3 * \ref n_pts)          
 *
 */

// void 
// CWP_user_tgt_pts_set
// (
//  const char                 *cpl_id,
//  const int                   n_pts,
//  double                      coord[]
// )
// {
//   cwipi::Coupling& cpl = _cpl_get(cpl_id);
//   cpl.userTgtPtsSet(n_pts, coord);
// }

/*----------------------------------------------------------------------------*
 * Functions about Support                                                    *
 *----------------------------------------------------------------------------*/

/**
 * \brief Setting vertices
 *
 * This function set partition vertices
 *
 * \param [in]  cpl_id      Coupling identifier
 * \param [in]  i_part      Current partition
 * \param [in]  n_pts       Number of points
 * \param [in]  coord       Coordinates (size = 3 * \ref n_pts)          
 * \param [in]  parent_num  Pointer to parent element number (or NULL)
 *
 */

// void 
// CWP_support_vtcs_set
// (
//  const char                 *cpl_id,
//  const int                   i_part,
//  const int                   n_pts,
//  const double                coord[],
//  const CWP_long_t          parent_num[]
// )
// {
//   cwipi::Coupling& cpl = _cpl_get(cpl_id);
//   cpl.supportVtcsSet(i_part,
//                      n_pts,
//                      coord[],
//                      parent_num[]);
// }

/**
 * \brief End setting support
 *
 * This function finalizes the support building
 *
 * \param [in]  cpl_id      Coupling identifier
 *
 */

// void 
// CWP_support_end_set
// (
//  const char                 *cpl_id
// )
// {
//   cwipi::Coupling& cpl = _cpl_get(cpl_id);
//   cpl.supportEndSet();
// }

/**
 * \brief Adding a connectivity block to the geometric support
 *
 * This function adds a connectivity block to the geometric support for
 * \ref CWP_SUPPORT_MESH support type. Definition of element connectivity is :
 *
 *  - edge (\ref CWP_BLOCK_EDGE2) :
 *
 *   \code
 *       1 x-------x 2
 *   \endcode
 *
 *  - triangle (\ref CWP_BLOCK_FACE_TRIA3):
 *
 *   \code 
 *       1 x-------x 3
 *          \     /
 *           \   /
 *            \ /
 *             x 2
 *   \endcode
 *
 *  - quadrangle (\ref CWP_BLOCK_FACE_QUAD4) :
 *
 *   \code
 *          4 x-------x 3
 *           /       /
 *          /       /
 *       1 x-------x2
 *   \endcode
 *
 *   - tetrahedron (\ref CWP_BLOCK_CELL_TETRA4) : 
 *
 *   \code
 *             x 4
 *            /|\
 *           / | \
 *          /  |  \
 *       1 x- -|- -x 3
 *          \  |  /
 *           \ | /
 *            \|/
 *             x 2
 *   \endcode
 *
 *   - pyramid (\ref CWP_BLOCK_CELL_PYRAM5) :
 *
 *   \code
 *              5 x
 *               /|\
 *              //| \
 *             // |  \
 *          4 x/--|---x 3
 *           //   |  /
 *          //    | /
 *       1 x-------x 2
 *   \endcode
 *
 *  - prism (\ref CWP_BLOCK_CELL_PRISM6) :
 *
 *   \code
 *       4 x-------x 6
 *         |\     /|
 *         | \   / |
 *       1 x- \-/ -x 3
 *          \ 5x  /
 *           \ | /
 *            \|/
 *             x 2
 *   \endcode
 *
 *  -  hexaedron (\ref CWP_BLOCK_CELL_HEXA8) :
 *
 *   \code
 *          8 x-------x 7
 *           /|      /|
 *          / |     / |
 *       5 x-------x6 |
 *         | 4x----|--x 3
 *         | /     | /
 *         |/      |/
 *       1 x-------x 2
 *   \endcode
 *
 * \param [in]  cpl_id      Coupling identifier
 * \param [in]  i_part      Current partition
 * \param [in]  block_type  Block type
 * \param [in]  n_elts      Number of elements
 * \param [in]  connec      Connectivity (size = n_vertex_elt * \ref n_elts)          
 * \param [in]  parent_num  Pointer to parent element number (or NULL)
 *
 */

// void 
// CWP_support_block_add
// (
//  const char                 *cpl_id,
//  const int                   i_part,
//  const CWP_block_t         block_type,
//  const int                   n_elts,
//  const int                   connec[],
//  const CWP_long_t          parent_num[])
// {
//   cwipi::Coupling& cpl = _cpl_get(cpl_id);
//   cpl.supportBlockAdd(i_part,
//                       block_type,
//                       n_elts,
//                       connec,
//                       parent_num);
// }

/**
 * \brief Adding a polygon connectivity block to the geometric support
 *
 * This function adds a polygon connectivity block to the geometric support for
 * \ref CWP_SUPPORT_MESH support type.
 *
 * \param [in]  cpl_id      Coupling identifier
 * \param [in]  i_part      Current partition
 * \param [in]  n_elts      Number of elements
 * \param [in]  connec_idx  Connectivity index (connec_id[0] = 0 and 
 *                          size = \ref n_elts + 1)          
 * \param [in]  connec      Connectivity (size = connec_id[n_elts] * \ref n_elts)
 * \param [in]  parent_num  Pointer to parent element number (or NULL)
 *
 */

// void 
// CWP_support_f_poly_block_add
// (
//  const char                 *cpl_id,
//  const int                   i_part,
//  const CWP_block_t         block_type,
//  const int                   n_elts,
//  const int                   connec[],
//  const CWP_long_t          parent_num[]
// )
// {
//   cwipi::Coupling& cpl = _cpl_get(cpl_id);
//   cpl.supportFPolyBlockAdd(i_part,
//                            block_type,
//                            n_elts,
//                            connec,
//                            parent_num);
// }

/**
 * \brief Adding a polyhedron connectivity block to the geometric support
 *
 * This function add a connectivity block to the geometric support if support
 * type is only \ref CWP_SUPPORT_MESH. Definition of element connectivity is :
 *
 * \param [in]  cpl_id            Coupling identifier
 * \param [in]  i_part            Current partition
 * \param [in]  n_elts            Number of elements
 * \param [in]  cell_face_idx     Polyhedron to face index 
 *                                (src_poly_cell_face_idx[0] = 0 and
 *                                 size = n_elts + 1)
 * \param [in]  cell_face_connec  Polyhedron to face connectivity 
 *                                (size = cell_face_idx[n_elts])
 * \param [in]  n_faces           Number of faces      
 * \param [in]  face_vtx_idx      Polyhedron face to vertex index 
 *                                (face_vertex_idx[0] = 0 and
 *                                 size_idx = max(cell_face_connec) + 1)
 * \param [in]  face_vtx_connec   Polyhedron face to vertex connectivity
 *                                (size = face_vertex_idx[size_idx - 1])
 * \param [in]  parent_num        Pointer to parent element number (or NULL)
 *
 */

// void 
// CWP_support_c_poly_block_add
// (
//  const char           *cpl_id,
//  const int             i_part,
//  const int             n_elts,
//  const int             cell_face_idx[],
//  const int             cell_face[],
//  const int             n_faces,
//  const int             face_vtx_idx[],
//  const int             face_vtx[],
//  const CWP_long_t    parent_num[]
// )
// {
//   cwipi::Coupling& cpl = _cpl_get(cpl_id);
//   cpl.supportCPolyBlockAdd(i_part,
//                            n_elts,
//                            cell_face_idx,
//                            cell_face,
//                            n_faces,
//                            face_vtx_idx,
//                            face_vtx,
//                            parent_num);
// }

/**
 * \brief Geometric support removal                                  
 *
 * This function delete the geometric support  
 *
 * \param [in]  cpl_id    Coupling identifier
 *
 */

// void 
// CWP_support_del
// (
//  const char *cpl_id
// )
// {
//   cwipi::Coupling& cpl = _cpl_get(cpl_id);
//   cpl.supportDel();
// }

/**
 * \brief Map a fvm nodal as support mesh                                  
 *
 * This function  map a fvm nodal as support mesh
 *
 * \param [in]  cpl_id            Coupling identifier
 * \param [in]  i_part            Current partition
 * \param [in]  fvmc_nodal        fvm nodal mes     
 *
 */

// void 
// CWP_support_shared_fvm_nodal
// (
//  const char    *cpl_id,
//  const int      i_part,
//  void          *fvmc_nodal
// )
// {
//   cwipi::Coupling& cpl = _cpl_get(cpl_id);
//   cpl.fvmcNodalShared(i_part,
//                       fvmc_nodal);
// }

/*----------------------------------------------------------------------------*
 * Functions about field                                                      *
 *----------------------------------------------------------------------------*/

/**
 *
 * \brief Creating a new field
 * 
 * \param [in]  cpl_id         Coupling identifier
 * \param [in]  field_id       Field id
 * \param [in]  data_type      Data type          
 * \param [in]  storage        Storage type          
 * \param [in]  n_component    Number of componenent
 * \param [in]  nature         Nature
 * \param [in]  exch_type      Exchange type
 * \param [in]  visu_status    Visualization status
 * 
 */

// void
// CWP_field_create
// (
//  const char                  *cpl_id,
//  const char                  *field_id,
//  const CWP_type_t     data_type,
//  const CWP_field_storage_t  storage,
//  const int                    n_component,
//  const CWP_Field_location_t   nature,
//  const CWP_field_exch_t     exchange_type,
//  const CWP_status_t         visu_status)
// {
//   cwipi::Coupling& cpl = _cpl_get(cpl_id);
//   cpl.fieldCreate(field_id,
//                   data_type,
//                   storage,
//                   n_component,
//                   nature,
//                   exchange_type, 
//                   visu_status);
// }

/**
 *
 * \brief Set data mapping
 * 
 * \param [in]  cpl_id         Coupling identifier
 * \param [in]  field_id       Field identifier
 * \param [in]  data           Storage array (Mapping)
 * 
 */

// void
// CWP_field_mapping_set
// (
//  const char                  *cpl_id,
//  const char                  *field_id,
//  double                       data[]   
//  )
// {
//   cwipi::Coupling& cpl = _cpl_get(cpl_id);
//   cpl.fieldMappingSet(field_id,
//                       data);   
// }

/**
 *
 * \brief Get nunmber of field components
 * 
 * \param [in]   cpl_id         Coupling identifier
 * \param [in]   field_id       Field identifier
 *
 * \return                      number of field components
 * 
 */

// int
// CWP_field_n_component_get
// (
//  const char                  *cpl_id,
//  const char                  *field_id
// )
// {
//   cwipi::Coupling& cpl = _cpl_get(cpl_id);
//   return cpl.fieldNComponentGet(field_id);
// }

/**
 *
 * \brief Get field nature
 * 
 * \param [in]   cpl_id         Coupling identifier
 * \param [in]   field_id       Field identifier
 *
 * \return                      Field nature
 * 
 */

// CWP_Field_location_t
// CWP_Field_location_get
// (
//  const char                  *cpl_id,
//  const char                  *field_id
// )
// {
//   cwipi::Coupling& cpl = _cpl_get(cpl_id);
//   return cpl.fieldNatureGet(field_id);
// }

/**
 *
 * \brief Get field data type
 * 
 * \param [in]   cpl_id         Coupling identifier
 * \param [in]   field_id       Field identifier
 *
 * \return                      Field data type
 * 
 */

// CWP_type_t
// CWP_field_type_get
// (
//  const char                  *cpl_id,
//  const char                  *field_id
// )
// {
//   cwipi::Coupling& cpl = _cpl_get(cpl_id);
//   return cpl.fieldTypeGet(field_id);
// }

/**
 *
 * \brief Get field storage type
 * 
 * \param [in]   cpl_id         Coupling identifier
 * \param [in]   field_id       Field identifier
 * 
 */

// CWP_field_storage_t
// CWP_field_storage_get
// (
//  const char                  *cpl_id,
//  const char                  *field_id
//  )
// {
//   cwipi::Coupling& cpl = _cpl_get(cpl_id);
//   return cpl.fieldStorageGet(field_id);
// }

/**
 *
 * \brief Removing a field
 * 
 * \param [in]  cpl_id         Coupling identifier
 * \param [in]  field_id       Field identifier
 * 
 */

// void
// CWP_field_del
// (
//  const char *cpl_id,
//  const char *field_id
// )
// {
//   cwipi::Coupling& cpl = _cpl_get(cpl_id);
//   cpl.fieldDelete();
// }

/*----------------------------------------------------------------------------*
 * Functions about exchange                                                   *
 *----------------------------------------------------------------------------*/

/**
 *
 * \brief Exchange data field with the coupled application with blocking 
 *        communications.
 *
 * This function exchanges interpolated fields between coupled codes. 
 * 
 * \warning  The size of \ref tgt_field_id size is n_computed_tgt. 
 *           If \f$ n\_uncomputed\_tgt \ne n\_tgt\_pts \f$,
 *           user himself must set values for uncomputed target points.
 *
 * \param [in]  cpl_id              Coupling identifier
 * \param [in]  src_field_id        Source field id (0 -> no sending)
 * \param [in]  tgt_field_id        Target field id (0 -> no receiving)
 * \param [out] n_uncomputed_tgt    Number of uncomputed target
 *
 * \return                          Exchange status
 *
 */

// CWP_error_t 
// CWP_sendrecv
// (
//  const char   *cpl_id,
//  const char   *src_field_id,
//  const char   *tgt_field_id,
//  int          *n_uncomputed_tgt
// )
// {
//   cwipi::Coupling& cpl = _cpl_get(cpl_id);
//   const string &src_field_id_str = src_field_id;
//   const string &tgt_field_id_str = tgt_field_id;

//   return cpl.sendRecv(src_field_id_str,
//                       tgt_field_id_str,
//                       NULL,
//                       n_uncomputed_tgt);
// }

/**
 *
 * \brief Sending of data field to the coupled application with nonblocking 
 *        communications.
 *
 * This function sends interpolated field to the coupled code. 
 * 
 * \param [in]  cpl_id          Coupling identifier
 * \param [in]  src_field_id    Source field id
 *
 * \param [out] request         Request to call by \ref CWP_wait_issend 
 *                              to wait the end of exchange
 *
 */

// void 
// CWP_issend
// (const char     *cpl_id,
//  const char     *src_field_id,
//  int            *request)
// {
//   cwipi::Coupling& cpl = _cpl_get(cpl_id);
//   const string &src_field_id_str = src_field_id;

//   cpl.issend(src_field_id_str,
//              NULL,
//              request);
// }

/**
 *
 * \brief Receiving of Data field from the coupled application with nonblocking 
 *        communications.
 *
 * This function receives interpolated field from the coupled code 
 * 
 * \param [in]  cpl_id          Coupling identifier
 * \param [in]  tgt_field_id    Target field id
 *
 * \param [out] request         Request to call by \ref CWP_wait_irecv  
 *                              to wait the end of exchange
 *
 */

// void 
// CWP_irecv
// (const char   *cpl_id,
//  const char   *tgt_field_id,
//  int          *request)
// {
//   cwipi::Coupling& cpl = _cpl_get(cpl_id);
//   const string &tgt_field_id_str = tgt_field_id;

//   cpl.irecv(tgt_field_id_str,
//             request);
// }

/**
 *
 * \brief Waiting of the end of exchange related to \ref request.
 *
 * This function waits the end of exchange related to \ref request
 * from \ref CWP_issend
 * 
 * \param [in] cpl_id     Coupling identifier
 * \param [in] request    Request to wait the end of exchange
 *
 */

// void 
// CWP_wait_issend
// (const char  *cpl_id,
//  int          request)
// {
//   cwipi::Coupling& cpl = _cpl_get(cpl_id);
//   cpl.waitIssend(request);
// }

/**
 *
 * \brief Waiting of the end of exchange related to \ref request.
 *
 * This function waits the end of exchange related to \ref request 
 * from \ref CWP_irecv
 * 
 * \param [in] cpl_id     Coupling identifier
 * \param [in] request    Request to wait the end of exchange
 *
 */

// void 
// CWP_wait_irecv
// (const char  *cpl_id,
//  int          request)
// {
//   cwipi::Coupling& cpl = _cpl_get(cpl_id);
//   cpl.waitIrecv(request);
// }


/*----------------------------------------------------------------------------*
 * Functions about user interpolation                                         *
 *----------------------------------------------------------------------------*/

/**
 *
 * \brief Setting of an user interpolation from location.
 *
 * This function takes into account an user interpolation function written with
 * \ref void (*CWP_interp_from_location_t) interface.
 * 
 * \param [in] cpl_id     Coupling identifier
 * \param [in] fct        Function
 *
 */

// void 
// CWP_interp_from_loc_set
// (
//  const char                  *cpl_id,
//  CWP_interp_from_location_t fct
// )
// {
//   cwipi::Coupling& cpl = _cpl_get(cpl_id);
//   cpl.interpFromLocSet(fct);
// }

/**
 *
 * \brief Setting of a FORTRAN user interpolation from location.
 *
 * This function takes into account an user interpolation function written
 * in FORTRAN .
 * 
 * \param [in] cpl_id     Coupling identifier
 * \param [in] fct        Function
 *
 */

// void 
// CWP_interp_from_loc_set_f
// (
//  const char *cpl_id,
//  void       *fct
// )
// {
//   cwipi::Coupling& cpl = _cpl_get(cpl_id);
//   cpl.interpFromLocSetF(fct);
// }

/**
 *
 * \brief Setting of an user interpolation from intersection.
 *
 * This function takes into account an user interpolation function written with
 * \ref void (*CWP_interp_from_intersec_t) interface.
 * 
 * \param [in] cpl_id     Coupling identifier
 * \param [in] fct        Function
 *
 */

// void 
// CWP_interp_from_inter_set
// (
//  const char                  *cpl_id,
//  CWP_interp_from_intersec_t fct
// )
// {
//   cwipi::Coupling& cpl = _cpl_get(cpl_id);
//   cpl.interpFromInterSet(fct);
// }

/**
 *
 * \brief Setting of a FORTRAN user interpolation from intersection.
 *
 * This function takes into account an user interpolation function written
 * in FORTRAN .
 * 
 * \param [in] cpl_id     Coupling identifier
 * \param [in] fct        Function
 *
 */

// void 
// CWP_interp_from_inter_set_f
// (
//  const char *cpl_id,
//  void       *fct
// )
// {
//   cwipi::Coupling& cpl = _cpl_get(cpl_id);
//   cpl.interpFromInterSetF(fct);
// }

/**
 *
 * \brief Setting of an user interpolation from closest points
 *
 * This function takes into account an user interpolation function written with
 * \ref void (*CWP_interp_from_closest_pts_t) interface.
 * 
 * \param [in] cpl_id     Coupling identifier
 * \param [in] fct        Function
 *
 */

// void 
// CWP_interp_from_closest_set
// (
//  const char                     *cpl_id,
//  CWP_interp_from_closest_pts_t fct
// )
// {
//   cwipi::Coupling& cpl = _cpl_get(cpl_id);
//   cpl.interpFromClosestSet(fct);
// }

/**
 *
 * \brief Setting of a FORTRAN user interpolation from closest points
 *
 * This function takes into account an user interpolation function written
 * in FORTRAN .
 * 
 * \param [in] cpl_id     Coupling identifier
 * \param [in] fct        Function
 *
 */

// void 
// CWP_interp_from_closest_set_f
// (
//  const char *cpl_id,
//  void       *fct
// )
// {
//   cwipi::Coupling& cpl = _cpl_get(cpl_id);
//   cpl.interpFromClosestSetF(fct);
// }

/*----------------------------------------------------------------------------*
 * Functions about current application control parameters                     *
 *----------------------------------------------------------------------------*/

/**
 *
 * \brief Add a control parameter
 * 
 * \param [in] local_code_name  Local code name
 * \param [in] param_name       Parameter name
 * \param [in] data_type        Parameter type
 * \param [in] init_value       Initial value
 *
 */

void 
CWP_Param_add
(
 const char        *local_code_name, 
 const char        *param_name, 
 const CWP_Type_t  data_type,
 void              *initial_value
)
{
  cwipi::CodePropertiesDB & properties =
    cwipi::CodePropertiesDB::getInstance();

  const string &codeNameStr = local_code_name;
  const string &nameStr = param_name;

  switch(data_type) {
  case CWP_INT :
    properties.ctrlParamAdd<int>(codeNameStr,nameStr, *(int *)initial_value);
    break;
  case CWP_DOUBLE :
    properties.ctrlParamAdd<double>(codeNameStr,nameStr, *(double *)initial_value);
    break;
  case CWP_CHAR :
    properties.ctrlParamAdd<char *>(codeNameStr,nameStr, *(char **)initial_value);
    break;
  default :
    bftc_error(__FILE__, __LINE__, 0,
               "Not yet implemented data type\n");
  }
}


/**
 *
 * \brief Set a local control parameter
 * 
 * \param [in] local_code_name  Local code name
 * \param [in] param_name       Parameter name
 * \param [in] data_type        Parameter type
 * \param [in] value          Value
 *
 */

void 
CWP_Param_set
(
 const char             *local_code_name,
 const char             *param_name,
 const CWP_Type_t        data_type,
 void                   *value
)
{
  cwipi::CodePropertiesDB & properties =
    cwipi::CodePropertiesDB::getInstance();

  const string &codeNameStr = local_code_name;
  const string &nameStr = param_name;
  switch(data_type) {
  case CWP_INT :
    properties.ctrlParamSet<int>(codeNameStr, nameStr, * (int *) value);
    break;
  case CWP_DOUBLE :
    properties.ctrlParamSet<double>(codeNameStr, nameStr, * (double *) value);
    break;
  case CWP_CHAR :
    properties.ctrlParamSet<char *>(codeNameStr, nameStr, * (char **) value);
    break;
  default :
    bftc_error(__FILE__, __LINE__, 0,
               "Not yet implemented data type\n");
  }
}

/**
 *
 * \brief Removing a local int control parameter
 * 
 * \param [in] local_code_name  Local code name
 * \param [in] param_name       Parameter name
 * \param [in] data_type        Parameter type,
 *
 */

void 
CWP_Param_del
(
 const char       *local_code_name,
 const char       *param_name,
 const CWP_Type_t  data_type
)
{
  cwipi::CodePropertiesDB & properties =
    cwipi::CodePropertiesDB::getInstance();

  const string &codeNameStr = local_code_name;
  const string &nameStr = param_name;
  switch(data_type) {
  case CWP_INT :
    properties.ctrlParamCancel<int>(codeNameStr, nameStr);
    break;
  case CWP_DOUBLE :
    properties.ctrlParamCancel<double>(codeNameStr, nameStr);
    break;
  case CWP_CHAR :
    properties.ctrlParamCancel<string>(codeNameStr, nameStr);
    break;
  default :
    bftc_error(__FILE__, __LINE__, 0,
               "Not yet implemented data type\n");
  }
}

/*----------------------------------------------------------------------------*
 * Functions about other application control parameters                       *
 *----------------------------------------------------------------------------*/

/**
 *
 * \brief Return the number of parameters of a code
 * 
 * 
 * \param [in] code_name       Local or distant code name
 * \param [in] data_type       Parameter type,
 *
 * \return  Number of parameters
 *
 */

int 
CWP_Param_n_get
(
 const char             *code_name,
 const CWP_Type_t        data_type
)
{
  cwipi::CodePropertiesDB & properties =
    cwipi::CodePropertiesDB::getInstance();

  const string &codeNameStr = code_name;
  int nParam;
  switch(data_type) {
  case CWP_INT :
    nParam = properties.ctrlParamNGet<int>(codeNameStr);
    break;
  case CWP_DOUBLE :
    nParam = properties.ctrlParamNGet<double>(codeNameStr);
    break;
  case CWP_CHAR :
    nParam = properties.ctrlParamNGet<string>(codeNameStr);
    break;
  default :
    bftc_error(__FILE__, __LINE__, 0,
                "Not yet implemented data type\n");
  }
  return nParam;
}


/**
 *
 * \brief Return the parameter list of a code
 * 
 * \param [in]  code_name      Local or distant code name
 * \param [in]  data_type      Parameter type,
 * \param [out] nParam         Number of parameters
 * \param [out] paramNames     Parameter names
 *
 *
 */

void
CWP_Param_list_get
(
 const char             *code_name,
 const CWP_Type_t        data_type,
 int                    *nParam,
 char                 ***paramNames   
)
{
  cwipi::CodePropertiesDB & properties =
    cwipi::CodePropertiesDB::getInstance();

  const string &codeNameStr = code_name;

  switch(data_type) {
  case CWP_INT :
    properties.ctrlParamListGet<int>(codeNameStr, nParam, paramNames);
    break;
  case CWP_DOUBLE :
    properties.ctrlParamListGet<double>(codeNameStr, nParam, paramNames);
    break;
  case CWP_CHAR :
    properties.ctrlParamListGet<string>(codeNameStr, nParam, paramNames);
    break;
  default :
    bftc_error(__FILE__, __LINE__, 0,
               "Not yet implemented data type\n");
  }
}


/**
 *
 * \brief Is a parameter ?
 * 
 * \param [in] code_name      Local or distant code name
 * \param [in] param_name     Parameter name
 * \param [in] data_type      Parameter type,
 *
 * \return  1 : true / 0 : false
 *
 */

int 
CWP_Param_is
(
 const char             *code_name,
 const char             *param_name,
 const CWP_Type_t        data_type
)
{
  cwipi::CodePropertiesDB & properties =
    cwipi::CodePropertiesDB::getInstance();

  const string &codeNameStr = code_name;
  const string &nameStr     = param_name;
  int isParam;
  switch(data_type) {
  case CWP_INT :
    isParam = properties.ctrlParamIs<int>(codeNameStr, nameStr);
    break;
  case CWP_DOUBLE :
    isParam = properties.ctrlParamIs<double>(codeNameStr, nameStr);
    break;
  case CWP_CHAR :
    isParam = properties.ctrlParamIs<string>(codeNameStr, nameStr);
    break;
  default :
    bftc_error(__FILE__, __LINE__, 0,
               "Not yet implemented data type\n");
  }
  return isParam;
}


/**
 *
 * \brief Return a control parameter
 * 
 * \param [in]  code_name  Local or distant code name
 * \param [in]  param_name Parameter name
 * \param [in]  data_type  Parameter type
 * \param [out] value      Parameter value       
 *
 */

void
CWP_Param_get
(
 const char       *code_name,
 const char       *param_name,
 const CWP_Type_t  data_type,
 void             *value
)
{
  cwipi::CodePropertiesDB & properties =
    cwipi::CodePropertiesDB::getInstance();

  const string &codeNameStr = code_name;
  const string &nameStr     = param_name;
  switch(data_type) {
  case CWP_INT : {
    int *intValue = (int *) value;
    *intValue = properties.ctrlParamGet<int>(codeNameStr, nameStr);
    break;
  }
  case CWP_DOUBLE : {
    double *dblValue = (double *) value;
    *dblValue = properties.ctrlParamGet<double>(codeNameStr, nameStr);
    break;
  }
  case CWP_CHAR : {
    char **charValue = (char **) value;
    *charValue = properties.ctrlParamGet<char *>(codeNameStr, nameStr);
    break;
  }
  default : {
    bftc_error(__FILE__, __LINE__, 0,
               "Not yet implemented data type\n");
  }
  }
}


/**
 *
 * \brief Return the operation result on a int control parameter 
 * (same name in two codes)
 * 
 * \param [in]  op     Operation
 * \param [in]  name   Parameter name
 * \param [in]  data_type  Parameter type
 * \param [out] res    Result   
 * \param [in]  nCode  Number of codes
 * \param       ...    Codes name
 *
 */

void
CWP_Param_reduce
(
 const CWP_Op_t    op,
 const char       *name,
 const CWP_Type_t  data_type,
 void             *res,
 const int         nCode,     
 ...
 )
{
  cwipi::CodePropertiesDB & properties =
    cwipi::CodePropertiesDB::getInstance();
  
  const string &nameStr = name;
  va_list pa;
  
  va_start(pa, nCode);
  switch(data_type) {
  case CWP_INT : {
    int *intRes = (int *) res;
    properties.ctrlParamReduce<int>(op, nameStr, intRes, nCode, &pa);
    break;
  }
  case CWP_DOUBLE : {
    double *doubleRes = (double *) res;
    properties.ctrlParamReduce<double>(op, nameStr, doubleRes, nCode, &pa);
    break;
  }
  default :
    bftc_error(__FILE__, __LINE__, 0,
               "Not yet implemented data type\n");
    
    va_end(pa);
  }
}


/**
 *
 * \brief Lock access to local parameters from a distant code 
 *
 * \param [in]  code_name  Code to lock
 * 
 */

void          
CWP_Param_lock
(
const char *code_name
)
{
  const string &nameStr = code_name;

  cwipi::CodePropertiesDB & properties =
    cwipi::CodePropertiesDB::getInstance();

  properties.lock(nameStr);
}


/**
 *
 * \brief unlock access to local parameters from a distant code 
 *
 * \param [in]  code_name  Code to lock
 * 
 */

void          
CWP_Param_unlock
(
const char *code_name
)
{
  const string &nameStr = code_name;

  cwipi::CodePropertiesDB & properties =
    cwipi::CodePropertiesDB::getInstance();

  properties.unLock(nameStr);
}

/*-----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
