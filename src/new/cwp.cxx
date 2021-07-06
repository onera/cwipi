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

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <cstring>
#include <cstdlib>

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * FVM library headers
 *----------------------------------------------------------------------------*/
#include <fvmc_parall.h>
#include <fvmc_nodal.h>

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
#include "field.hxx"
#include "pdm.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "surfMeshGenerator.hxx"
#include "surfMeshGeneratorDB.hxx"
#include "spatialInterpIntersection.hxx"
#include "spatialInterpClosestPoint.hxx"
#include "spatialInterpLocationDistSurf.hxx"
#include "spatialInterpLocationMeshLocation.hxx"

 #include "mesh.hxx"
 #include "block.hxx"
 #include "blockStd.hxx"
 #include "blockFP.hxx"
 #include "blockCP.hxx"
 #include <algorithm>
 #include <vector>
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
#include <stdlib.h>


static bool
_is_coupled
(
 const char *local_code_name
 )
{
   cwipi::CodePropertiesDB & propertiesDB = cwipi::CodePropertiesDB::getInstance();
   const cwipi::CodeProperties &properties = propertiesDB.codePropertiesGet(local_code_name);
   return properties.isCoupledRank();
}

static bool
_cpl_exist
(
 const char *local_code_name,
 const char *cpl_id
 )
{
  cwipi::CouplingDB & couplingDB =
    cwipi::CouplingDB::getInstance();

  cwipi::CodePropertiesDB & properties =
    cwipi::CodePropertiesDB::getInstance();

   const string &cpl_name_str = cpl_id;
   return couplingDB.couplingIs(properties.codePropertiesGet(string(local_code_name)),
                             cpl_name_str);
}

static cwipi::Coupling&
_cpl_get
(
 const char *local_code_name,
 const char *cpl_id
 )
{
  cwipi::CouplingDB & couplingDB =
    cwipi::CouplingDB::getInstance();

  cwipi::CodePropertiesDB & properties =
    cwipi::CodePropertiesDB::getInstance();

   const string &cpl_name_str = cpl_id;
   return couplingDB.couplingGet (properties.codePropertiesGet(string(local_code_name)),
                                 cpl_name_str);
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
 * This function creates the MPI intra communicators of the codes from
 * the \p global_comm MPI communicator that contains all code ranks. This
 * function has to be called from all ranks contained in the \p global_comm.
 *
 * \param [in]  global_comm    MPI global communicator
 * \param [in]  n_code         Number of codes on the current rank
 * \param [in]  code_names     Names of codes on the current rank (size = \p n_code)
 * \param [in]  is_active_rank Is current rank have to be used by CWIPI (size = \p n_code)
 * \param [in]  time_init      Initial time (size = \p n_code)
 * \param [out] intra_comms    MPI intra communicators of each code (size = \p n_code)
 *
 */

void
CWP_Init
(
 const MPI_Comm           global_comm,
 const int                n_code,
 const char             **code_names,
 const CWP_Status_t      *is_active_rank,
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

  //PDM_printf("\ncwipi " CWIPI_VERSION " initializing\n");
  //PDM_printf("------------------------\n\n");

  /*
   * Builds application communicator
   */

  properties.init (global_comm,
                   n_code,
                   code_names,
                   is_active_rank,
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
   * Create spatial interpolation abstract factory
   */

  cwipi::Factory<cwipi::SpatialInterp, CWP_Spatial_interp_t> &factorySpatialInterp =
    cwipi::Factory<cwipi::SpatialInterp, CWP_Spatial_interp_t>::getInstance();

  factorySpatialInterp.Register<cwipi::SpatialInterpLocationDistSurf>(CWP_SPATIAL_INTERP_FROM_LOCATION_DIST_CLOUD_SURF);
  factorySpatialInterp.Register<cwipi::SpatialInterpLocationMeshLocationOctree>(CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE);
  factorySpatialInterp.Register<cwipi::SpatialInterpLocationMeshLocationDbbtree>(CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_DBBTREE);
  factorySpatialInterp.Register<cwipi::SpatialInterpIntersection>(CWP_SPATIAL_INTERP_FROM_INTERSECTION);
  factorySpatialInterp.Register<cwipi::SpatialInterpClosestPoint>(CWP_SPATIAL_INTERP_FROM_CLOSEST_POINT_LEAST_SQUARES);

  /*
   * Create block abstract factory
   */

  cwipi::Factory<cwipi::Block, CWP_Block_t> &factoryBlock =
    cwipi::Factory<cwipi::Block, CWP_Block_t>::getInstance();

  factoryBlock.Register<cwipi::BlockStd>(CWP_BLOCK_NODE);
  factoryBlock.Register<cwipi::BlockStd>(CWP_BLOCK_EDGE2);
  factoryBlock.Register<cwipi::BlockStd>(CWP_BLOCK_FACE_TRIA3);
  factoryBlock.Register<cwipi::BlockStd>(CWP_BLOCK_FACE_QUAD4);
  factoryBlock.Register<cwipi::BlockFP >(CWP_BLOCK_FACE_POLY);
  factoryBlock.Register<cwipi::BlockStd>(CWP_BLOCK_CELL_TETRA4);
  factoryBlock.Register<cwipi::BlockStd>(CWP_BLOCK_CELL_HEXA8);
  factoryBlock.Register<cwipi::BlockStd>(CWP_BLOCK_CELL_PRISM6);
  factoryBlock.Register<cwipi::BlockStd>(CWP_BLOCK_CELL_PYRAM5);
  factoryBlock.Register<cwipi::BlockCP > (CWP_BLOCK_CELL_POLY);

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
 * \brief Finalize CWIPI.
 *
 */

void
CWP_Finalize
(
 void
)
{
  int flag = 0;

  MPI_Initialized(&flag);

  cwipi::CodePropertiesDB & properties =
    cwipi::CodePropertiesDB::getInstance();

  const MPI_Comm globalComm = properties.globalCommGet();

 // PDM_printf("CWP_Finalize\n");
  fflush(stdout);
  if (flag != 0) {
    PDM_printf_flush();
    MPI_Barrier(globalComm);
//    MPI_Comm oldFVMComm = fvmc_parall_get_mpi_comm();
  }



  properties.kill();

}

/*----------------------------------------------------------------------------*
 * Functions about current application properties                             *
 *----------------------------------------------------------------------------*/

/**
 * \brief Update code state.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] state            State
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
 * \brief Return code state.
 *
 * \param [in]  code_name    Code name
 *
 * \return      Code state
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
 * \brief Return the number of codes known by CWIPI.
 *
 * \return Number of codes
 *
 */

int
CWP_Codes_nb_get
(void
)
{
  cwipi::CodePropertiesDB & properties =
    cwipi::CodePropertiesDB::getInstance();

  return properties.codesNbGet();
}

/**
 * \brief Return the list of code names known by CWIPI.
 *
 * \return list of codes.
 */

const char **
CWP_Codes_list_get
(void
)
{
  cwipi::CodePropertiesDB & properties =
    cwipi::CodePropertiesDB::getInstance();
  return properties.codesListGet();
}

/**
 * \brief Return the number of local codes known by CWIPI.
 *
 * \return number of local codes.
 */

int
CWP_Loc_codes_nb_get
(void
)
{
  cwipi::CodePropertiesDB & properties =
    cwipi::CodePropertiesDB::getInstance();

  return properties.localCodesNbGet();
}


/**
 * \brief Return the list of local code names known by CWIPI.
 *
 * \return list of local codes.
 */

const char **
CWP_Loc_codes_list_get
(void
)
{
  cwipi::CodePropertiesDB & properties =
    cwipi::CodePropertiesDB::getInstance();
  return properties.localCodesListGet();
}


/**
 * \brief Update code time.
 *
 * \param [in] local_code_name  Local code name
 * \param [in]  current_time Current time
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
 * \brief Define output file.
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
  PDM_printf_proxy_set(_cwipi_print_with_c);
}

//
//**
// * \brief Writing output to fortran file.
// *
// * This function set the file fortran logical unit for writing output.
// *
// * \param [in]  iunit        File fortan logical unit
// *
// */
//
//void
//PROCF (cwp_output_fortran_unit_set, CWP_OUTPUT_FORTRAN_UNIT_SET)
//(
// int *iunit
//);

/*----------------------------------------------------------------------------*
 * Functions about properties                                                 *
 *----------------------------------------------------------------------------*/


/**
 * \brief Dump code properties.
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
 * \brief Create a coupling object and define its properties.
 *
 * \param [in]  local_code_name     Local code name
 * \param [in]  cpl_id              Coupling identifier
 * \param [in]  coupled_code_name   Distant or local coupled code name
 * \param [in]  entities_dim        Dimensions of the entities
 * \param [in]  comm_type           Communication type
 * \param [in]  spatial_interp      Spatial interpolation method
 * \param [in]  n_part              Number of interface partition
 * \param [in]  displacement        Mesh moving status
 * \param [in]  recv_freq_type      Type of receiving frequency
 *
 */

void
CWP_Cpl_create
(
 const char                 *local_code_name,
 const char                 *cpl_id,
 const char                 *coupled_code_name,
 CWP_Interface_t             entities_dim,
 const CWP_Comm_t            comm_type,
 const CWP_Spatial_interp_t  spatial_interp,
 const int                   n_part,
 const CWP_Dynamic_mesh_t    displacement,
 const CWP_Time_exch_t       recv_freq_type
)
{
  cwipi::CouplingDB & couplingDB =
    cwipi::CouplingDB::getInstance();

  cwipi::CodePropertiesDB & properties =
    cwipi::CodePropertiesDB::getInstance();

  const string &coupling_name_str = cpl_id;
  const string &coupled_application_str = coupled_code_name;
  const string &local_application_str = local_code_name;


 // printf("couplingDB.couplingCreate(proper\n");
  couplingDB.couplingCreate(properties.codePropertiesGet(local_application_str),
                            coupling_name_str,
                            properties.codePropertiesGet(coupled_application_str),
                            entities_dim,
                            comm_type,
                            spatial_interp,
                            n_part,
                            displacement,
                            recv_freq_type);
}

/**
 *
 * \brief Delete a coupling object.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 *
 */

void
CWP_Cpl_del
(
const char *local_code_name,
const char *cpl_id
)
{
  cwipi::CouplingDB & couplingDB =
    cwipi::CouplingDB::getInstance();

  cwipi::CodePropertiesDB & properties =
    cwipi::CodePropertiesDB::getInstance();

  const string &cpl_id_str = cpl_id;

  couplingDB.couplingDel(properties.codePropertiesGet(string(local_code_name)),
                         cpl_id_str);
}


/**
 * \brief Exchange spatially interpolated fields. <b>(Not implemented yet)</b>
 *
 * This function exchanges the interpolated fields for each coupling depending
 * on mode of time exchange \ref CWP_Time_exch_t.
 *
 * \param [in] local_code_name      Local code name
 * \param [in] cpl_id               Coupling identifier
 *
 */

void
CWP_Field_exch
(
 const char *local_code_name,
 const char *cpl_id
)
{
  //TODO: Voir comment enchainer les appels, voir comment prendre en compte l'interpolation temporelle
  cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
  // cpl.exchange();
}


/**
 *
 * \brief Return the number of uncomputed targets.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 * \param [in] target_location  Target location
 * \param [in] i_part           Current partition
 *
 * \return                Number of uncomputed targets
 */


int
CWP_N_uncomputed_tgts_get
(
 const char               *local_code_name,
 const char               *cpl_id,
 const CWP_Dof_location_t  target_location,
 const int                 i_part
)
{
   cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);

   return cpl.nUncomputedTargetsGet(target_location,i_part);
}


const int *
CWP_Uncomputed_tgts_get
(
  const char *local_code_name,
  const char *cpl_id
)
{
  // TODO
//  cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
//  return cpl.uncomputedTargetsGet();
}


/**
 *
 * \brief Return the number of computed targets. <b>(Not implemented yet)</b>
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 *
 * \return                Number of computed targets
 */

int
CWP_N_computed_tgts_get
(
  const char *local_code_name,
  const char *cpl_id
)
{
  // TODO
//  cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
//  return cpl.nComputedTargetsGet();
}



const int *
CWP_Computed_tgts_get
(
  const char *local_code_name,
  const char *cpl_id
)
{
  // TODO
//  cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
//  return cpl.getLocatedPoint();
}



const double *
CWP_Computed_tgts_dist_to_spatial_interp_get
(
  const char *local_code_name,
  const char *cpl_id
)
{
  // TODO
}


/*----------------------------------------------------------------------------*
 * Functions about exchange frequency                                         *
 *----------------------------------------------------------------------------*/



// void
// CWP_recv_freq_set
// (
//  const char                 *cpl_id,
//  const int                   n_step
// )
// {
//   cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
//   cpl.recvFreqSet(n_step);
// }


/**
 * \brief Set the next receiving time.
 *
 * It must be used when the type of receiving frequency is
 * \ref CWP_TIME_EXCH_ASYNCHRONOUS
 *
 * \param [in]  local_code_name  Local code name
 * \param [in]  cpl_id           Coupling identifier
 * \param [in]  next_time        Next receiving time
 *
 */

 void
 CWP_next_recv_time_set
 (const char     *local_code_name,
  const char     *cpl_id,
  const double    next_time
 )
 {
   cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
   cpl.recvNextTimeSet(next_time);
 }

/*----------------------------------------------------------------------------*
 * Functions about spatial interpolation                                      *
 *----------------------------------------------------------------------------*/

/**
 * \brief Compute spatial interpolation weights.
 *
 * \param [in]  local_code_name     Local code name
 * \param [in]  cpl_id              Coupling identifier
 *
 */

void
CWP_Spatial_interp_weights_compute
(
 const char     *local_code_name,
 const char     *cpl_id
)
{

  typedef struct field_exch_type {
    CWP_Dof_location_t loc ;
    CWP_Field_exch_t  exch;
  };

  cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);

  /* Only if the coupling exists */
  map <string, cwipi::Field *>* fields = cpl.fieldsGet();

  std::map <std::string, field_exch_type> field_exch_type_map;
  std::string fieldName="";
  vector<int> fieldNameIdx;
  fieldNameIdx.push_back(0);
  std::vector<CWP_Field_exch_t> fieldExch;
  std::vector<CWP_Dof_location_t> fieldLocationV;
  std::map <std::string, cwipi::Field *>::iterator it = fields -> begin();

  while(it != fields -> end()){
    cwipi::Field* field = it -> second;
    CWP_Dof_location_t fieldLocation = field -> typeGet();
    CWP_Field_exch_t  exchangeType  = field -> exchangeTypeGet();
    field_exch_type field_eT;
    field_eT.loc  = fieldLocation;
    field_eT.exch = exchangeType ;
    field_exch_type_map.insert( std::pair<string,field_exch_type>(it->first, field_eT) );
    fieldName += it->first;
    fieldNameIdx.push_back(fieldNameIdx[fieldNameIdx.size() - 1] + it->first.size());
    fieldLocationV.push_back(fieldLocation);
    fieldExch.push_back(exchangeType);

    it++;
  }

  int nb_field = fieldNameIdx.size() -1;
  vector<int              > fieldNameIdx_cpl (nb_field+1,0);
  vector<CWP_Field_exch_t > fieldExch_cpl    (nb_field  );
  vector<CWP_Dof_location_t> fieldLocationV_cpl(nb_field  );
  string fieldName_cpl;

  MPI_Comm unionComm = cpl.communicationGet() -> unionCommGet();
  int unionCommCplCodeRootRank = cpl.communicationGet() -> unionCommCplCodeRootRanksGet();
  int unionCommLocCodeRootRank = cpl.communicationGet() -> unionCommLocCodeRootRanksGet();

  int rank;
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);

  int tag = 152;
  int unionCommRank;
  MPI_Status status;
  MPI_Comm_rank (unionComm, &unionCommRank);


  int CWP_Field_value_size = sizeof(int);
  if (unionCommRank == unionCommLocCodeRootRank) {
    MPI_Sendrecv (&(fieldNameIdx[0]),
                  fieldNameIdx.size(),
                  MPI_INT,
                  unionCommCplCodeRootRank,
                  tag,
                  &(fieldNameIdx_cpl[0]),
                  fieldNameIdx_cpl.size(),
                  MPI_INT,
                  unionCommCplCodeRootRank,
                  tag,
                  unionComm,
                  &status);
    tag++;
    fieldName_cpl.resize(fieldNameIdx_cpl[nb_field]);

    MPI_Sendrecv (&(fieldName[0]),
                  fieldName.size(),
                  MPI_CHAR,
                  unionCommCplCodeRootRank,
                  tag,
                  &(fieldName_cpl[0]),
                  fieldName_cpl.size(),
                  MPI_CHAR,
                  unionCommCplCodeRootRank,
                  tag,
                  unionComm,
                  &status);

    tag++;
    MPI_Sendrecv (&(fieldLocationV[0]),
                  sizeof(CWP_Dof_location_t) * fieldLocationV.size(),
                  MPI_BYTE,
                  unionCommCplCodeRootRank,
                  tag,
                  &(fieldLocationV_cpl[0]),
                  sizeof(CWP_Dof_location_t) * fieldLocationV_cpl.size(),
                  MPI_BYTE,
                  unionCommCplCodeRootRank,
                  tag,
                  unionComm,
                  &status);

    tag++;
    MPI_Sendrecv (&(fieldExch[0]),
                  sizeof(CWP_Field_exch_t) * fieldExch.size(),
                  MPI_BYTE,
                  unionCommCplCodeRootRank,
                  tag,
                  &(fieldExch_cpl[0]),
                  sizeof(CWP_Field_exch_t) * fieldExch_cpl.size(),
                  MPI_BYTE,
                  unionCommCplCodeRootRank,
                  tag,
                  unionComm,
                  &status);

  }


  int id     = cpl.localCodePropertiesGet()   -> idGet();
  int id_cpl = cpl.coupledCodePropertiesGet() -> idGet();
  bool both_local = cpl.spatialInterpGet(CWP_DOF_LOCATION_NODE)->_both_codes_are_local;
  if(both_local == 0 || (both_local == 1 && id < id_cpl) ){

    std::vector <int> tmp(3,0);
    std::vector <int> tmp_fieldNameIdx(fieldNameIdx.size(),0);
    std::vector <CWP_Dof_location_t> tmp_fieldLocationV(fieldLocationV.size());
    std::vector <CWP_Field_exch_t> tmp_fieldExch(fieldExch.size());

    string tmp_fieldName;
    if (id < id_cpl) {

      MPI_Bcast (&(fieldNameIdx_cpl[0]),
                 fieldNameIdx_cpl.size(),
                 MPI_INT,
                 unionCommLocCodeRootRank,
                 unionComm);

      MPI_Bcast (&(tmp_fieldNameIdx[0]),
                 tmp_fieldNameIdx.size(),
                 MPI_INT,
                 unionCommCplCodeRootRank,
                 unionComm);

      if (unionCommRank != unionCommLocCodeRootRank) {
        fieldName_cpl.resize(fieldNameIdx_cpl[nb_field],'r');
      }
      tmp_fieldName.resize(fieldNameIdx_cpl[nb_field],'n');

      MPI_Bcast (&(fieldName_cpl[0]),
                 fieldName_cpl.size(),
                 MPI_CHAR,
                 unionCommLocCodeRootRank,
                 unionComm);
      MPI_Bcast (&(tmp_fieldName[0]),
                 tmp_fieldName.size(),
                 MPI_CHAR,
                 unionCommCplCodeRootRank,
                 unionComm);

      MPI_Bcast (&(fieldLocationV_cpl[0]),
                 sizeof(CWP_Dof_location_t) * fieldLocationV_cpl.size(),
                 MPI_BYTE,
                 unionCommLocCodeRootRank,
                 unionComm);

      MPI_Bcast (&(tmp_fieldLocationV[0]),
                 sizeof(CWP_Dof_location_t) * tmp_fieldLocationV.size(),
                 MPI_BYTE,
                 unionCommCplCodeRootRank,
                 unionComm);

      MPI_Bcast (&(fieldExch_cpl[0]),
                 sizeof(CWP_Field_exch_t) * fieldExch_cpl.size(),
                 MPI_BYTE,
                 unionCommLocCodeRootRank,
                 unionComm);

      MPI_Bcast (&(tmp_fieldExch[0]),
                 sizeof(CWP_Field_exch_t) * tmp_fieldExch.size(),
                 MPI_BYTE,
                 unionCommCplCodeRootRank,
                 unionComm);

    }
      else{

         MPI_Bcast (&(tmp_fieldNameIdx[0]),
                    tmp_fieldNameIdx.size(),
                    MPI_INT,
                    unionCommCplCodeRootRank,
                    unionComm);
         MPI_Bcast (&(fieldNameIdx_cpl[0]),
                    fieldNameIdx_cpl.size(),
                    MPI_INT,
                    unionCommLocCodeRootRank,
                    unionComm);

         if (unionCommRank != unionCommLocCodeRootRank) {
           fieldName_cpl.resize(fieldNameIdx_cpl[nb_field]);
         }
         tmp_fieldName.resize(fieldNameIdx_cpl[nb_field],'n');

         MPI_Bcast (&(tmp_fieldName[0]),
                    tmp_fieldName.size(),
                    MPI_CHAR,
                    unionCommCplCodeRootRank,
                    unionComm);

         MPI_Bcast (&(fieldName_cpl[0]),
                    fieldName_cpl.size(),
                    MPI_CHAR,
                    unionCommLocCodeRootRank,
                    unionComm);

         MPI_Bcast (&(tmp_fieldLocationV[0]),
                    sizeof(CWP_Dof_location_t) * tmp_fieldLocationV.size(),
                    MPI_BYTE,
                    unionCommCplCodeRootRank,
                    unionComm);

         MPI_Bcast (&(fieldLocationV_cpl[0]),
                    sizeof(CWP_Dof_location_t) * fieldLocationV_cpl.size(),
                    MPI_BYTE,
                    unionCommLocCodeRootRank,
                    unionComm);

         MPI_Bcast (&(tmp_fieldExch[0]),
                    sizeof(CWP_Field_exch_t) * tmp_fieldExch.size(),
                    MPI_BYTE,
                    unionCommCplCodeRootRank,
                    unionComm);

         MPI_Bcast (&(fieldExch_cpl[0]),
                    sizeof(CWP_Field_exch_t) * fieldExch_cpl.size(),
                    MPI_BYTE,
                    unionCommLocCodeRootRank,
                    unionComm);
      }
  }
     else {

       cwipi::Coupling& cpl_cpl = _cpl_get(cpl.coupledCodePropertiesGet()->nameGet().c_str(),cpl_id);
       /* Only if the coupling exists */
       map <string, cwipi::Field *>* fields_cpl = cpl_cpl.fieldsGet();
       it = fields_cpl -> begin();
       fieldNameIdx_cpl.resize(0);
       fieldNameIdx_cpl.push_back(0);
       fieldLocationV_cpl.resize(0);
       fieldExch_cpl.resize(0);
       fieldName_cpl="";

       while (it != fields_cpl -> end()) {

        fieldNameIdx_cpl.push_back( fieldNameIdx_cpl[fieldNameIdx_cpl.size()-1] + it->first.size() );
        fieldName_cpl += it->first;
        cwipi::Field* field_cpl = it -> second;
        CWP_Dof_location_t fieldLocation = field_cpl -> typeGet();
        fieldLocationV_cpl.push_back(fieldLocation);
        CWP_Field_exch_t  exchangeType  = field_cpl -> exchangeTypeGet();
        fieldExch_cpl.push_back(exchangeType);

        it++;
      }
     }

     std::map<string,field_exch_type> field_cpl_map;
     for (int i = 0; i < nb_field; i++) {
       string field_name = fieldName_cpl.substr( fieldNameIdx_cpl[i], fieldNameIdx_cpl[i+1]-fieldNameIdx_cpl[i] );
       field_exch_type field_exch;
       field_exch.loc  = fieldLocationV_cpl[i];
       field_exch.exch = fieldExch_cpl     [i];
       field_cpl_map.insert( std::pair<string,field_exch_type>(field_name,field_exch) );
       //printf("field_name %s %s %i %i\n",field_name.c_str(),fieldName_cpl.c_str(),fieldNameIdx_cpl[i], fieldNameIdx_cpl[i+1]);
     }

    static const char *CWP_Dof_location_t_str[] = {"CWP_DOF_LOCATION_CELL_CENTER","CWP_DOF_LOCATION_NODE","CWP_DOF_LOCATION_USER"};
    static const char *CWP_Field_exch_t_str [] = {"CWP_FIELD_EXCH_SEND","CWP_FIELD_EXCH_RECV","CWP_FIELD_EXCH_SENDRECV"};
    std::vector<int> exchangeTypeByLocation(3,0);

     it = fields -> begin();
     while (it != fields -> end()) {
       if(it -> second ->  exchangeTypeGet() == CWP_FIELD_EXCH_SEND) {
         it -> second -> associatedCloudPointTypeSet(field_cpl_map[it->first].loc);
       }
       else if (it -> second ->  exchangeTypeGet() == CWP_FIELD_EXCH_RECV) {
         it -> second -> associatedCloudPointTypeSet( it -> second -> typeGet() );
       }
       else if (it -> second ->  exchangeTypeGet() == CWP_FIELD_EXCH_SENDRECV) {
         it -> second -> associatedCloudPointTypeSet( it -> second -> typeGet() );
       }
       else{
         PDM_error(__FILE__, __LINE__, 0, "Not correct exchange field value for this field.\n");
       }

       if (exchangeTypeByLocation[ static_cast<int>( it -> second -> associatedCloudPointTypeGet() )] == 0) {
         if (it -> second -> exchangeTypeGet()==CWP_FIELD_EXCH_SEND )
           exchangeTypeByLocation[ static_cast<int>( it -> second -> associatedCloudPointTypeGet() )] = 1;
         else if (it -> second -> exchangeTypeGet()==CWP_FIELD_EXCH_RECV)
           exchangeTypeByLocation[ static_cast<int>( it -> second -> associatedCloudPointTypeGet() )] = 2;
         else if (it -> second -> exchangeTypeGet()==CWP_FIELD_EXCH_SENDRECV)
           exchangeTypeByLocation[ static_cast<int>( it -> second -> associatedCloudPointTypeGet() )] = 3;
       }
       else {
         if ( exchangeTypeByLocation[ static_cast<int>( it -> second -> associatedCloudPointTypeGet() )] == 1 && it -> second -> exchangeTypeGet() != CWP_FIELD_EXCH_SEND){
            exchangeTypeByLocation[ static_cast<int>( it -> second -> associatedCloudPointTypeGet() )] = 3;
         }
         else if ( exchangeTypeByLocation[ static_cast<int>( it -> second -> associatedCloudPointTypeGet() )] == 2 && it -> second -> exchangeTypeGet() != CWP_FIELD_EXCH_RECV){
            exchangeTypeByLocation[ static_cast<int>( it -> second -> associatedCloudPointTypeGet() )] = 3;
         }
       }

       //printf(" %s typeGet() %s associatedCloudPointTypeGet %s rank %i exchangeTypeGet %i\n",it->first.c_str(),CWP_Dof_location_t_str[it -> second -> typeGet()],
       //CWP_Dof_location_t_str[static_cast<int>(it -> second -> associatedCloudPointTypeGet())],rank,it -> second ->  exchangeTypeGet());
       it++;
     }


      /*  Building of possible cloud points type vector */
      std::vector<CWP_Dof_location_t> locationV = {CWP_DOF_LOCATION_CELL_CENTER, CWP_DOF_LOCATION_NODE, CWP_DOF_LOCATION_USER};

      // Iteration over the possilbe cloud points type
      for(size_t i_location=0; i_location < locationV.size(); i_location++){
        CWP_Dof_location_t pointsCloudLocation = locationV[i_location];
        int spatialInterpComputeSend  = 0;
        int spatialInterpComputeRcv = 0;

        if (exchangeTypeByLocation[ static_cast<int>( pointsCloudLocation ) ] == 1) {
          spatialInterpComputeSend = 1;
        }
        else if (exchangeTypeByLocation[ static_cast<int>( pointsCloudLocation ) ] == 2) {
          spatialInterpComputeRcv = 1;
        }
        else if (exchangeTypeByLocation[ static_cast<int>( pointsCloudLocation ) ] == 3) {
          spatialInterpComputeSend = 1;
          spatialInterpComputeRcv  = 1;
        }
        else if (exchangeTypeByLocation[ static_cast<int>( pointsCloudLocation ) ] == 0) {
          spatialInterpComputeSend = 0;
          spatialInterpComputeRcv  = 0;
        }

        CWP_Field_exch_t exchange_type    ;
        CWP_Field_exch_t exchange_type_cpl;
        if (spatialInterpComputeRcv == 1 && spatialInterpComputeSend == 1 ) {
          if (id < id_cpl) {
             if (both_local == 1) {

              cwipi::Coupling& cpl_cpl = _cpl_get(cpl.coupledCodePropertiesGet() ->nameGet().c_str(),cpl_id);

              cpl.spatialInterpWeightsCompute(pointsCloudLocation, CWP_FIELD_EXCH_SEND);
              cpl_cpl.spatialInterpWeightsCompute(pointsCloudLocation,CWP_FIELD_EXCH_RECV);

              cpl.spatialInterpWeightsCompute(pointsCloudLocation, CWP_FIELD_EXCH_RECV);
              cpl_cpl.spatialInterpWeightsCompute(pointsCloudLocation,CWP_FIELD_EXCH_SEND);
             }
             else if (both_local == 0) {
              cpl.spatialInterpWeightsCompute(pointsCloudLocation, CWP_FIELD_EXCH_SEND);
              cpl.spatialInterpWeightsCompute(pointsCloudLocation, CWP_FIELD_EXCH_RECV);
            }
          }
          else {
            if (both_local == 0) {
              cpl.spatialInterpWeightsCompute(pointsCloudLocation, CWP_FIELD_EXCH_RECV);
              cpl.spatialInterpWeightsCompute(pointsCloudLocation, CWP_FIELD_EXCH_SEND);
            }
          }
         // printf("pointsCloudLocation %s rank %i %i %i id<id_cpl %i\n", CWP_Dof_location_t_str[static_cast<int>( pointsCloudLocation )],rank, spatialInterpComputeSend, spatialInterpComputeRcv,id<id_cpl );

        }
        else if (spatialInterpComputeRcv == 1 && spatialInterpComputeSend == 0) {
          exchange_type     = CWP_FIELD_EXCH_RECV ;
          exchange_type_cpl = CWP_FIELD_EXCH_SEND ;
          if(both_local == 1 && id < id_cpl) {
            cpl.spatialInterpWeightsCompute(pointsCloudLocation, exchange_type);
            cwipi::Coupling& cpl_cpl = _cpl_get(cpl.coupledCodePropertiesGet() ->nameGet().c_str(),cpl_id);
            cpl_cpl.spatialInterpWeightsCompute(pointsCloudLocation, exchange_type_cpl);
         }
          else if (both_local == 0) {
            cpl.spatialInterpWeightsCompute(pointsCloudLocation, exchange_type);
            //printf("pointsCloudLocation %s rank %i %i %i\n", CWP_Dof_location_t_str[static_cast<int>( pointsCloudLocation )],rank, spatialInterpComputeSend, spatialInterpComputeRcv );
          }
        }
        else if (spatialInterpComputeSend == 1 && spatialInterpComputeRcv == 0) {
          exchange_type     = CWP_FIELD_EXCH_SEND ;
          exchange_type_cpl = CWP_FIELD_EXCH_RECV ;
          if(both_local == 1 && id < id_cpl) {
            cpl.spatialInterpWeightsCompute(pointsCloudLocation, exchange_type);
            cwipi::Coupling& cpl_cpl = _cpl_get(cpl.coupledCodePropertiesGet() ->nameGet().c_str(),cpl_id);
            cpl_cpl.spatialInterpWeightsCompute(pointsCloudLocation, exchange_type_cpl);
          }
          else if (both_local == 0) {
            cpl.spatialInterpWeightsCompute(pointsCloudLocation, exchange_type);
           // printf("pointsCloudLocation %s rank %i %i %i\n", CWP_Dof_location_t_str[static_cast<int>( pointsCloudLocation )],rank, spatialInterpComputeSend, spatialInterpComputeRcv );
          }
        }
       if((both_local == 1 && id < id_cpl) || both_local == 0) MPI_Barrier(unionComm);
      } //end on location loop
  }



// void
// CWP_spatial_interp_update
// (
//  const char     *cpl_id,
//  const int       storage_id
// )
// {
// TODO
// }


// void
// CWP_spatial_interp_properties_set
// (
//  const char     *cpl_id,
//  const char     *fmt,
//                   ...
// )
// {
//   cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);

//   va_list pa;
//   va_start(pa, fmt);
//   cpl.spatialInterpPropertiesSet(fmt, &pa);
//   va_end(pa);
// }

/*----------------------------------------------------------------------------*
 * Functions about visualization                                              *
 *----------------------------------------------------------------------------*/


/**
 * \brief Enable visualization output.
 *
 * \param [in]  local_code_name  Local code name
 * \param [in]  cpl_id           Coupling identifier
 * \param [in]  freq             Output frequency
 * \param [in]  format           Output format to visualize exchanged fieldsDouble
 *                               on the coupled mesh. Choice between :
 *                               - "EnSight Gold"
 * \param [in]  format_option   Output options "opt1, opt2, ..."
 *                               - text : output text files
 *                               - binary : output binary files (default)
 */

void
CWP_Visu_set
(
 const char                 *local_code_name,
 const char                 *cpl_id,
 const int                   freq,
 const CWP_Visu_format_t     format,
 const char                 *format_option
)
{
  if(_is_coupled(local_code_name)){
    cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
    cpl.visuSet(freq, format, format_option);
  }
}

/*----------------------------------------------------------------------------*
 * Functions about User target points                                         *
 *----------------------------------------------------------------------------*/


/**
 * \brief Setting user target points.
 *
 * This function must be called if the degrees of freedom locations are
 * \ref CWP_DOF_LOCATION_USER
 *
 * \param [in]  local_code_name  Local code name
 * \param [in]  cpl_id           Coupling identifier
 * \param [in]  i_part           Current partition
 * \param [in]  n_pts            Number of points
 * \param [in]  coord            Coordinates (size = 3 * n_pts)
 *
 */

void
CWP_User_tgt_pts_set
(
 const char                 *local_code_name,
 const char                 *cpl_id,
 const int                   i_part,
 const int                   n_pts,
 double                      coord[]
)
{
  if(_is_coupled(local_code_name)){
    cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
    cpl.userTgtPtsSet(i_part, n_pts, coord);
  }
}

/*----------------------------------------------------------------------------*
 * Functions about Mesh                                                    *
 *----------------------------------------------------------------------------*/

/**
 * \brief Finalize interface mesh.
 *
 * This function computes the global numbers of mesh entities if they are
 * not provided.
 *
 * \param [in]  local_code_name  Local code name
 * \param [in]  cpl_id           Coupling identifier
 *
 */

void
CWP_Mesh_interf_finalize
(
 const char                 *local_code_name,
 const char                 *cpl_id
)
{
  if(_is_coupled(local_code_name)){
    cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
    cpl.meshFinalize();
  }
}



/**
 * \brief Set vertices.
 *
 * \param [in]  local_code_name  Local code name
 * \param [in]  cpl_id           Coupling identifier
 * \param [in]  i_part           Current partition
 * \param [in]  n_pts            Number of points
 * \param [in]  coord            Coordinates (size = 3 * \p n_pts)
 * \param [in]  global_num       Pointer to parent element number (or NULL)
 *
 */

void
CWP_Mesh_interf_vtx_set
(
 const char                 *local_code_name,
 const char                 *cpl_id,
 const int                   i_part,
 const int                   n_pts,
 double                      coord[],
 CWP_g_num_t                 global_num[]
)
{
  if(_is_coupled(local_code_name)){
    cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
    cpl.meshVtcsSet(i_part,
                    n_pts,
                    coord,
                    global_num);
  }
}

/**
 * \brief Add a connectivity block to the interface mesh.
 *
 * \param [in]  local_code_name  Local code name
 * \param [in]  cpl_id           Coupling identifier
 * \param [in]  block_type       Block type
 *
 * \return block identifier
 */

int
CWP_Mesh_interf_block_add
(
 const char           *local_code_name,
 const char           *cpl_id,
 const CWP_Block_t     block_type
)
{
  if(_is_coupled(local_code_name)){
    cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
    int block_id = cpl.meshBlockAdd(block_type);
    return block_id;
  }
}


/**
 * \brief Set a standard block to the interface mesh.
 *
 * This function adds a connectivity block to the interface mesh.
 * Definition of element connectivity is :
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
 * \param [in]  local_code_name  Local code name
 * \param [in]  cpl_id           Coupling identifier
 * \param [in]  i_part           Partition identifier
 * \param [in]  block_id         Block identifier
 * \param [in]  n_elts           Number of elements
 * \param [in]  connec           Connectivity (size = n_vertex_elt * n_elts)
 * \param [in]  global_num       Pointer to global element number (or NULL)
 *
 */

void
CWP_Mesh_interf_block_std_set
(
 const char        *local_code_name,
 const char        *cpl_id,
 const int          i_part,
 const int          block_id,
 const int          n_elts,
 int                connec[],
 CWP_g_num_t       global_num[]
)
{
  if(_is_coupled(local_code_name)){
    cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
    cpl.meshStdBlockSet(i_part,
                        block_id,
                        n_elts,
                        connec,
                        global_num);
  }
}


/*void
CWP_Mesh_interf_h_order_block_set
(
 const char        *local_code_name,
 const char        *cpl_id,
 const int          i_part,
 const int          block_id,
 const int          n_elts,
 const int          order,
 int                connec[],
 CWP_g_num_t        global_num[]
)
{
  cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
  cpl.meshHighOrderBlockSet(i_part,
                            block_id,
                            block_type,
                            n_elts,
                            order,
                            connec,
                            global_num);
}
*/

/**
 * \brief Set the connectivity of a polygon block in a interface mesh partition.
 *
 * \param [in]  local_code_name  Local code name
 * \param [in]  cpl_id           Coupling identifier
 * \param [in]  i_part           Current partition
 * \param [in]  block_id         Block identifier
 * \param [in]  n_elts           Number of elements
 * \param [in]  connec_idx       Connectivity index (\p connec_id[0] = 0 and
 *                               size = \p n_elts + 1)
 * \param [in]  connec           Connectivity (size = \p connec_idx[\p n_elts])
 * \param [in]  global_num       Pointer to global element number (or NULL)
 *
 */

void
CWP_Mesh_interf_f_poly_block_set
(
 const char             *local_code_name,
 const char             *cpl_id,
 const int               i_part,
 const int               block_id,
 const int               n_elts,
 int                     connec_idx[],
 int                     connec[],
 CWP_g_num_t             global_num[]
)
{

   if(_is_coupled(local_code_name)){
     cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
     cpl.meshFPolyBlockSet(i_part,
                           block_id,
                           n_elts,
                           connec_idx,
                           connec,
                           global_num);
  }
}


/**
 * \brief Adding a polyhedron connectivity block to the interface mesh.
 *
 * \param [in]  local_code_name   Local code name
 * \param [in]  cpl_id            Coupling identifier
 * \param [in]  i_part            Current partition
 * \param [in]  block_id          Block identifier
 * \param [in]  n_elts            Number of elements
 * \param [in]  connec_cells_idx  Polyhedron to face index
 *                                (\p src_poly_cell_face_idx[0] = 0 and
 *                                 size = \p n_elts + 1)
 * \param [in]  connec_cells      Polyhedron to face connectivity
 *                                (size = \p cell_face_idx[\p n_elts])
 * \param [in]  n_faces           Number of faces
 * \param [in]  connec_faces_idx  Polyhedron face to vertex index
 *                                (\p face_vertex_idx[0] = 0 and
 *                                 size = max(\p cell_face_connec) + 1)
 * \param [in]  connec_faces      Polyhedron face to vertex connectivity
 *                                (size = \p face_vertex_idx[\p n_elts])
 * \param [in]  global_num        Pointer to global element number (or NULL)
 *
 */

void
CWP_Mesh_interf_c_poly_block_set
(
 const char           *local_code_name,
 const char           *cpl_id,
 const int             i_part,
 const int             block_id,
 const int             n_elts,
 const int             n_faces,
 int                   connec_faces_idx[],
 int                   connec_faces[],
 int                   connec_cells_idx[],
 int                   connec_cells[],
 CWP_g_num_t           global_num[]
)
{
   cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
   cpl.meshCPolyBlockSet(i_part,
                         block_id,
                         n_elts,
                         n_faces,
                         connec_faces_idx,
                         connec_faces    ,
                         connec_cells_idx,
                         connec_cells    ,
                         global_num);
}


/**
 * \brief Define the interface mesh from a cell to face connectivity.
 *
 * \param [in]  local_code_name   Local code name
 * \param [in]  cpl_id            Coupling identifier
 * \param [in]  i_part            Current partition
 * \param [in]  n_cells           Number of cells
 * \param [in]  cell_face_idx     Polyhedron to face index
 *                                (\p src_poly_cell_face_idx[0] = 0 and
 *                                 size = \p n_elts + 1)
 * \param [in]  cell_face         Cell to face connectivity
 *                                (size = \p cell_face_idx[\p n_elts])
 * \param [in]  n_faces           Number of faces
 * \param [in]  face_vtx_idx      Polyhedron face to vertex index
 *                                (\p face_vtx_idx[0] = 0 and
 *                                 size = \p n_faces + 1)
 * \param [in]  face_vtx          Face to vertex connectivity
 *                                (size = \p face_vtx_idx[\p n_elts])
 * \param [in]  parent_num        Pointer to parent element number (or NULL)
 *
 */


void
CWP_Mesh_interf_from_cellface_set
(
 const char      *local_code_name,
 const char      *cpl_id,
 const int       i_part,
 const int       n_cells,
 int             cell_face_idx[],
 int             cell_face[],
 const int       n_faces,
 int             face_vtx_idx[],
 int             face_vtx[],
 CWP_g_num_t     parent_num[]
)
{
  cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
  cpl.meshFromCellFaceSet(i_part,
                          n_cells,
                          cell_face_idx,
                          cell_face,
                          n_faces,
                          face_vtx_idx,
                          face_vtx,
                          parent_num);
}


/**
 * \brief Define the surface interface mesh from a face to edge connectivity.
 *
 * \param [in]  local_code_name   Local code name
 * \param [in]  cpl_id            Coupling identifier
 * \param [in]  i_part            Current partition
 * \param [in]  n_faces           Number of cells
 * \param [in]  face_edge_idx     Polygon to edge index
 *                                (\p face_edge_idx[0] = 0 and
 *                                 size =  \p n_faces + 1)
 * \param [in]  face_edge         Face to edge connectivity
 *                                (size = \p face_edge_idx[\p n_faces])
 * \param [in]  n_edges           Number of faces
 * \param [in]  edge_vtx_idx      Polyhedron face to vertex index
 *                                (\p edge_vtx_idx[0] = 0 and
 *                                 size = \p n_edges + 1)
 * \param [in]  edge_vtx          Face to vertex connectivity
 *                                (size = \p edge_vtx_idx[\p n_edges])
 * \param [in]  parent_num        Pointer to parent element number (or NULL)
 *
 */

void
CWP_Mesh_interf_from_faceedge_set
(
 const char           *local_code_name,
 const char           *cpl_id,
 const int             i_part,
 const int             n_faces,
 int             face_edge_idx[],
 int             face_edge[],
 const int             n_edges,
 int             edge_vtx_idx[],
 int             edge_vtx[],
 CWP_g_num_t     parent_num[]
)
{
  cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
  cpl.meshFromFacesEdgeSet(i_part,
                        n_faces,
                        face_edge_idx,
                        face_edge,
                        n_edges,
                        edge_vtx_idx,
                        edge_vtx,
                        parent_num);
}


/**
 * \brief Delete interface mesh.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 *
 */

void
CWP_Mesh_interf_del
(
 const char *local_code_name,
 const char *cpl_id
)
{
   cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
   cpl.meshDel();
}

/*
void
CWP_Mesh_interf_shared_fvm_nodal
(
 const char   *local_code_name,
 const char   *cpl_id,
 const int     i_part,
 fvmc_nodal_t *fvmc_nodal
)
{
   cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
   cpl.fvmcNodalShared(i_part,
                       fvmc_nodal);
}
*/
/*----------------------------------------------------------------------------*
 * Functions about field                                                      *
 *----------------------------------------------------------------------------*/

/**
 *
 * \brief Create a new field.
 *
 * \param [in] local_code_name Local code name
 * \param [in]  cpl_id         Coupling identifier
 * \param [in]  field_id       Field id
 * \param [in]  data_type      Data type
 * \param [in]  storage        Storage type
 * \param [in]  n_component    Number of component
 * \param [in]  target_location Target location
 * \param [in]  exch_type      Exchange type
 * \param [in]  visu_status    Visualization status
 *
 */

void
CWP_Field_create
(
 const char                  *local_code_name,
 const char                  *cpl_id,
 const char                  *field_id,
 const CWP_Type_t             data_type,
 const CWP_Field_storage_t    storage,
 const int                    n_component,
 const CWP_Dof_location_t     target_location,
 const CWP_Field_exch_t       exch_type,
 const CWP_Status_t           visu_status)
{
  cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
  cpl.fieldCreate(field_id,
                  data_type,
                  storage,
                  n_component,
                  target_location,
                  exch_type,
                  visu_status);
}


/**
 *
 * \brief Setting of an user interpolation from location.
 *
 * This function takes into account an user interpolation function written with
 * void (*\ref CWP_Interp_from_location_t) interface.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 * \param [in] src_field_id    Source field id
 * \param [in] fct              Function
 *
 */

void
CWP_Interp_from_location_set
(
 const char                 *local_code_name,
 const char                 *cpl_id,
 const char                 *src_field_id,
  CWP_Interp_from_location_t fct
 )
{
  cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
  cpl.interpFromLocSet(src_field_id,fct);
}


/**
 *
 * \brief Set field data.
 *
 * \param [in] local_code_name   Local code name
 * \param [in] cpl_id            Coupling identifier
 * \param [in] field_id          Field identifier
 * \param [in] i_part            Current partition
 * \param [in] data              Storage array (Mapping)
 *
 */

void
CWP_Field_data_set
(
 const char           *local_code_name,
 const char           *cpl_id,
 const char           *field_id,
 const int             i_part,
 double                data[]
)
{
  if(_is_coupled(local_code_name)){
    cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
    cpl.fieldDataSet(field_id,i_part,data);
  }
}


/**
 *
 * \brief Get number of field components.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 * \param [in] field_id         Field identifier
 *
 * \return                      number of field components
 *
 */

int
CWP_Field_n_component_get
(
 const char                  *local_code_name,
 const char                  *cpl_id,
 const char                  *field_id
)
{
  cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
  return cpl.fieldNComponentGet(field_id);
}


/**
 *
 * \brief Get target degrees of freedom location.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 * \param [in] field_id         Field identifier
 *
 * \return                      Location of degrees of freedom
 *
 */

CWP_Dof_location_t
CWP_Field_target_dof_location_get
(
 const char                  *local_code_name,
 const char                  *cpl_id,
 const char                  *field_id
)
{
  cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
  return cpl.fieldTypeGet(field_id);
}

/*
 CWP_Type_t
 CWP_Field_type_get
 (
  const char                  *local_code_name,
  const char                  *cpl_id,
  const char                  *field_id
 )
 {
   cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
   return cpl.fieldTypeGet(field_id);
 }*/


/**
 *
 * \brief Get field storage type.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 * \param [in] field_id         Field identifier
 *
 * \return                      Field storage type
 */

CWP_Field_storage_t
CWP_Field_storage_get
(
 const char                  *local_code_name,
 const char                  *cpl_id,
 const char                  *field_id
)
{
  cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
  return cpl.fieldStorageGet(field_id);
}


/**
 * \brief Delete a field.
 *
 * \param [in] local_code_name Local code name
 * \param [in]  cpl_id         Coupling identifier
 * \param [in]  field_id       Field identifier
 *
 */

void
CWP_Field_del
(
 const char                  *local_code_name,
 const char                  *cpl_id,
 const char                  *field_id
)
{
  cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
  cpl.fieldDel(field_id);
}

/*----------------------------------------------------------------------------*
 * Functions about exchange                                                   *
 *----------------------------------------------------------------------------*/

/**
 * \brief Send a spatially interpolated field to the coupled code with
 *        nonblocking communications.
 *
 * This function is independant of \ref CWP_Time_exch_t mode. The user has to
 * manually check the consistency of the exchanges.
 *
 * \param [in]  local_code_name  Local code name
 * \param [in]  cpl_id          Coupling identifier
 * \param [in]  src_field_id    Source field id
 *
 *
 */

void
CWP_Field_issend
(
 const char        *local_code_name,
 const char        *cpl_id,
 const char        *src_field_id
)
{
  if(_is_coupled(local_code_name)){
    cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
    std::string referenceFieldID_str = const_cast<char*>(src_field_id);

    cpl.issend(referenceFieldID_str);
  }
}


/**
 *
 * \brief Receive a spatially interpolated field from the coupled code
 *        with nonblocking communications.
 *
 * This function is independant of \ref CWP_Time_exch_t mode. The user has to
 * manually check the consistency of the exchanges.
 *
 * \param [in] local_code_name  Local code name
 * \param [in]  cpl_id          Coupling identifier
 * \param [in]  tgt_field_id    Target field id
 *
 *
 */

void
CWP_Field_irecv
(
 const char        *local_code_name,
 const char        *cpl_id,
 const char        *tgt_field_id
)
{
  if(_is_coupled(local_code_name)){
    cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
    std::string targetFieldID_str = const_cast<char*>(tgt_field_id);
    cpl.irecv(targetFieldID_str);
  }
}


/**
 *
 * \brief Wait the end of an exchange related to request from \ref CWP_Field_issend.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 * \param [in] src_field_id     Source field id
 *
 */

void
CWP_Field_wait_issend
(
 const char  *local_code_name,
 const char  *cpl_id,
 const char  *src_field_id
)
{
  if(_is_coupled(local_code_name)){
    cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
    std::string srcFieldID_str = const_cast<char*>(src_field_id);
    cpl.waitIssend(srcFieldID_str);
  }
}


/**
 *
 * \brief Wait the end of an exchange related to request from \ref CWP_Field_irecv.
 *
 * This function waits the end of exchange related to request
 * from \ref CWP_Field_irecv
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 * \param [in] tgt_field_id     Target field id
 *
 */

void
CWP_Field_wait_irecv
(
 const char  *local_code_name,
 const char  *cpl_id,
 const char  *tgt_field_id
)
{
  if(_is_coupled(local_code_name)){
    cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
    std::string distantFieldID_str = const_cast<char*>(tgt_field_id);
    cpl.waitIrecv(distantFieldID_str);
  }
}


/*----------------------------------------------------------------------------*
 * Functions about user interpolation                                         *
 *----------------------------------------------------------------------------*/


// void
// CWP_Interp_from_loc_set
// (
//  const char                  *cpl_id,
//  CWP_Interp_from_location fct
// )
// {
//   cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
//   cpl.interpFromLocSet(fct);
// }


// void
// CWP_Interp_from_loc_set_f
// (
//  const char *cpl_id,
//  void       *fct
// )
// {
//   cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
//   cpl.interpFromLocSetF(fct);
// }


// void
// CWP_Interp_from_inter_set
// (
//  const char                  *cpl_id,
//  CWP_Interp_from_intersec_t fct
// )
// {
//   cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
//   cpl.interpFromInterSet(fct);
// }


// void
// CWP_Interp_from_inter_set_f
// (
//  const char *cpl_id,
//  void       *fct
// )
// {
//   cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
//   cpl.interpFromInterSetF(fct);
// }


// void
// CWP_Interp_from_closest_set
// (
//  const char                     *cpl_id,
//  CWP_Interp_from_closest_pts_t fct
// )
// {
//   cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
//   cpl.interpFromClosestSet(fct);
// }


// void
// CWP_Interp_from_closest_set_f
// (
//  const char *cpl_id,
//  void       *fct
// )
// {
//   cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
//   cpl.interpFromClosestSetF(fct);
// }

/*----------------------------------------------------------------------------*
 * Functions about current application control parameters                     *
 *----------------------------------------------------------------------------*/

/**
 *
 * \brief Add a new parameter and intialize it.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] param_name       Parameter name
 * \param [in] data_type        Parameter type
 * \param [in] initial_value    Initial value
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
    PDM_error(__FILE__, __LINE__, 0,
               "Not yet implemented data type\n");
  }
}



/**
 *
 * \brief Set a parameter.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] param_name       Parameter name
 * \param [in] data_type        Parameter type
 * \param [in] value            Value
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
    PDM_error(__FILE__, __LINE__, 0,
               "Not yet implemented data type\n");
  }
}



/**
 *
 * \brief Delete a parameter.
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
    PDM_error(__FILE__, __LINE__, 0,
               "Not yet implemented data type\n");
  }
}

/*----------------------------------------------------------------------------*
 * Functions about other application control parameters                       *
 *----------------------------------------------------------------------------*/


/**
 *
 * \brief Return the number of parameters for the code \p code_name.
 *
 * \param [in] code_name       Local or distant code name
 * \param [in] data_type       Parameter type,
 *
 * return  Number of parameters
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
    PDM_error(__FILE__, __LINE__, 0,
                "Not yet implemented data type\n");
  }
  return nParam;
}



/**
 *
 * \brief Return the list of parameters for the code \p code_name.
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
    PDM_error(__FILE__, __LINE__, 0,
               "Not yet implemented data type\n");
  }
}


/**
 *
 * \brief Is this \p code_name a parameter ?
 *
 * \param [in] code_name      Local or distant code name
 * \param [in] param_name     Parameter name
 * \param [in] data_type      Parameter type,
 *
 * return  1 : true / 0 : false
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
    PDM_error(__FILE__, __LINE__, 0,
               "Not yet implemented data type\n");
  }
  return isParam;
}


/**
 *
 * \brief Return the parameter value of \p param_name on \p code_name.
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
    PDM_error(__FILE__, __LINE__, 0,
              "Not yet implemented data type\n");
  }
  }
}


/**
 *
 * \brief Return the result of a reduce operation about a parameter.
 *
 * The parameter name has to be the same for all codes.
 *
 * \param [in]  op           Operation
 * \param [in]  param_name   Parameter name
 * \param [in]  data_type    Parameter type,
 * \param [out] res          Result
 * \param [in]  nCode        Number of codes
 * \param       ...          Codes name
 *
 */

void
CWP_Param_reduce
(
 const CWP_Op_t    op,
 const char       *param_name,
 const CWP_Type_t  data_type,
 void             *res,
 const int         nCode,
 ...
 )
{
  cwipi::CodePropertiesDB & properties =
    cwipi::CodePropertiesDB::getInstance();

  const string &nameStr = param_name;
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
    PDM_error(__FILE__, __LINE__, 0,
               "Not yet implemented data type\n");

    va_end(pa);
  }
}



/**
 *
 * \brief Lock access to local parameters from a distant code.
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
 * \brief Unlock access to local parameters from a distant code.
 *
 * \param [in]  code_name  Code to unlock
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



// Fonctions privee qui n'appartiennent pas a l'API qui doivent etre deplacees



void
CWP_surf_gen_init
(
  char* genName, int nx, int ny, int nPart, MPI_Comm* comm, double prop, double width, double randomVar
)
{
  if(_is_coupled(genName)){
    cwipi::surfMeshGeneratorDB &surfMeshDB =  cwipi::surfMeshGeneratorDB::getInstance();
    surfMeshDB.createMember(genName);
    cwipi::surfMeshGenerator* surfMesh = surfMeshDB.memberGet(genName);
    surfMesh -> init(nx,ny,nPart,comm,prop,width,randomVar);
  }
}

void
CWP_surf_gen_compute(char* genName)
{

  if(_is_coupled(genName)){
    cwipi::surfMeshGeneratorDB &surfMeshDB =  cwipi::surfMeshGeneratorDB::getInstance();
    cwipi::surfMeshGenerator* surfMesh = surfMeshDB.memberGet(genName);
    surfMesh -> computeMesh();
  }
}



void
CWP_surf_gen_by_block_get
( char* genName, int i_part,
  int* nVtx , double** coords, CWP_g_num_t** vtxGnum, int* nElts,
  int* nTri , int** eltsConnecTri , CWP_g_num_t** eltsGnumTri,
  int* nQuad, int** eltsConnecQuad, CWP_g_num_t** eltsGnumQuad,
  int* nPoly, int** eltsConnecPolyIndex, int** eltsConnecPoly, CWP_g_num_t** eltsGnumPoly
)
{


  if(_is_coupled(genName)){
    cwipi::surfMeshGeneratorDB &surfMeshDB =  cwipi::surfMeshGeneratorDB::getInstance();
    cwipi::surfMeshGenerator* surfMesh = surfMeshDB.memberGet(genName);
    *nVtx = surfMesh -> nVtxGet(i_part);
    *nElts = surfMesh -> nEltsGet(i_part);
    *coords = surfMesh -> coordsGet(i_part);
    *vtxGnum = surfMesh -> vtxGnumGet(i_part);

    *nTri = surfMesh -> nTriGet(i_part);
    *eltsConnecTri = surfMesh -> connecTriGet(i_part);
    *eltsGnumTri = surfMesh -> eltsGnumTriGet(i_part);

    *nQuad = surfMesh -> nQuadGet(i_part);
    *eltsConnecQuad = surfMesh -> connecQuadGet(i_part);
    *eltsGnumQuad = surfMesh -> eltsGnumQuadGet(i_part);

    *nPoly = surfMesh -> nPolyGet(i_part);
    *eltsConnecPoly = surfMesh -> connecPolyGet(i_part);
    *eltsConnecPolyIndex = surfMesh -> connecPolyIndexGet(i_part);
    *eltsGnumPoly = surfMesh -> eltsGnumPolyGet(i_part);
  }
  else{
    *nVtx = 0;
    *nElts = 0;
    *coords = NULL;
    *vtxGnum = NULL;

    *nTri = 0;
    *eltsConnecTri = NULL;
    *eltsGnumTri = NULL;

    *nQuad = 0;
    *eltsConnecQuad = NULL;
    *eltsGnumQuad = NULL;

    *nPoly = 0;
    *eltsConnecPoly = NULL;
    *eltsConnecPolyIndex = (int*)malloc(sizeof(int));
    (*eltsConnecPolyIndex)[0]=0;
    *eltsGnumPoly = NULL;
  }


}


void
CWP_surf_gen_one_connectivity_get
( char* genName, int i_part,
  int* nVtx , double** coords, CWP_g_num_t** vtxGnum,
  int* nElts, int** eltsConnecIndex, int** eltsConnec, CWP_g_num_t** eltsGnum
)
{

  if(_is_coupled(genName)){
    cwipi::surfMeshGeneratorDB &surfMeshDB =  cwipi::surfMeshGeneratorDB::getInstance();
    cwipi::surfMeshGenerator* surfMesh = surfMeshDB.memberGet(genName);

    *nVtx = surfMesh -> nVtxGet(i_part);
    *nElts = surfMesh -> nEltsGet(i_part);
    *coords = surfMesh -> coordsGet(i_part);
    *vtxGnum = surfMesh -> vtxGnumGet(i_part);

    *eltsConnec = surfMesh -> connecGet(i_part);
    *eltsConnecIndex = surfMesh -> connecIndexGet(i_part);
    *eltsGnum = surfMesh -> eltsGnumGet(i_part);
  }
  else{
    *nVtx = 0;
    *nElts = 0;
    *coords = NULL;
    *vtxGnum = NULL;

    *eltsConnec = NULL;
    *eltsConnecIndex = (int*)malloc(sizeof(int));
    (*eltsConnecIndex)[0]=0;
    *eltsGnum = NULL;
  }


}



void
CWP_surf_face_edge_get
( char* genName, int i_part,
  int* nVtx , double** coords, CWP_g_num_t** vtxGnum,
  int* nFace, int** faceEdgeIdx, int** faceEdge,
  int* nEdge, int** edgeVtxIdx, int** edgeVtx,
  CWP_g_num_t** faceLNToGN
)
{

  if(_is_coupled(genName)){
    cwipi::surfMeshGeneratorDB &surfMeshDB =  cwipi::surfMeshGeneratorDB::getInstance();
    cwipi::surfMeshGenerator* surfMesh = surfMeshDB.memberGet(genName);

    *nVtx = surfMesh -> nVtxGet(i_part);
    *nFace = surfMesh -> nFaceGet(i_part);
    *coords = surfMesh -> coordsGet(i_part);
    *vtxGnum = surfMesh -> vtxGnumGet(i_part);

    *nEdge = surfMesh -> nEdgeGet(i_part);
    *faceEdgeIdx = surfMesh -> faceEdgeIdxGet(i_part);
    *faceEdge = surfMesh -> faceEdgeGet(i_part);

    *edgeVtxIdx = surfMesh -> edgeVtxIdxGet(i_part);
    *edgeVtx = surfMesh -> edgeVtxGet(i_part);

    *faceLNToGN = surfMesh -> faceLNToGNGet(i_part);
  }
  else{
    *nVtx = 0;
    *coords = NULL;
    *vtxGnum = NULL;

    *faceEdgeIdx =(int*)malloc(sizeof(int));
    (*faceEdgeIdx)[0]=0;

    *faceEdge = NULL;

    *edgeVtxIdx = (int*)malloc(sizeof(int));
    (*edgeVtxIdx)[0]=0;

    *edgeVtx = NULL;

    *faceLNToGN = NULL;

  }


}






void
CWP_surf_gen_tri_field_get
( char* genName, int i_part,
  double** field
)
{

  if(_is_coupled(genName)){
    cwipi::surfMeshGeneratorDB &surfMeshDB =  cwipi::surfMeshGeneratorDB::getInstance();
    cwipi::surfMeshGenerator* surfMesh = surfMeshDB.memberGet(genName);
    *field = surfMesh -> specialFieldTriGet(i_part);
  }
  else{
    *field = NULL;
  }
}


void
CWP_surf_gen_quad_field_get
( char* genName, int i_part,
  double** field
)
{

  if(_is_coupled(genName)){
    cwipi::surfMeshGeneratorDB &surfMeshDB =  cwipi::surfMeshGeneratorDB::getInstance();
    cwipi::surfMeshGenerator* surfMesh = surfMeshDB.memberGet(genName);
    *field = surfMesh -> specialFieldQuadGet(i_part);
  }
  else{
    *field = NULL;
  }
}

void
CWP_surf_gen_poly_field_get
( char* genName, int i_part,
  double** field
)
{

  if(_is_coupled(genName)){
    cwipi::surfMeshGeneratorDB &surfMeshDB =  cwipi::surfMeshGeneratorDB::getInstance();
    cwipi::surfMeshGenerator* surfMesh = surfMeshDB.memberGet(genName);
    *field = surfMesh -> specialFieldPolyGet(i_part);
  }
  else{
    *field = NULL;
  }
}




 CWP_g_num_t*
 CWP_GlobalNumGet
 (
  const char  *local_code_name,
  const char  *cpl_id,
  const int    id_block,
  const int    i_part
 )
 {
   cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);

   return cpl.globalNumGet(id_block,i_part);

 }



MPI_Comm
CWP_Connectable_comm_get
(
  char* local_code_name
)
{

  if(_is_coupled(local_code_name)){
    cwipi::CodePropertiesDB & properties = cwipi::CodePropertiesDB::getInstance();
    const cwipi::CodeProperties & localCodeProperties = properties.codePropertiesGet(string(local_code_name));
    MPI_Comm connecComm = localCodeProperties.connectableCommGet();

    return connecComm;
  }

  return MPI_COMM_NULL;
}

/*-----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
