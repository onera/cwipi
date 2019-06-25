
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

#include "geometry.hxx"
#include "geomLocation.hxx"

 #include "mesh.hxx"
 #include "block.hxx"
 #include "blockStd.hxx"
 #include "blockFP.hxx"
 #include "blockCP.hxx"
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
 * This function create the MPI intra communicator for this code from
 * the MPI inter communicator that contains all code process. It is a
 * synchronization point between all codes
 *
 * \param [in]  global_comm       MPI global communicator
 * \param [in]  n_code            Number of codes on the current rank
 * \param [in]  code_names         Names of codes on the current rank (size = n_code)
 * \param [in]  is_coupled_rank   Is current rank used for coupling (size = n_code)
 * \param [in]  time_init         Time init (size = n_code)
 * \param [out] intra_comms        MPI intra communicators of each code
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

  PDM_printf("\ncwipi " CWIPI_VERSION " initializing\n");
  PDM_printf("------------------------\n\n");

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

  cwipi::Factory<cwipi::Geometry, CWP_Geom_t> &factoryGeom = 
    cwipi::Factory<cwipi::Geometry, CWP_Geom_t>::getInstance();

  factoryGeom.Register<cwipi::GeomLocation>(CWP_GEOM_LOCATION);
  // factoryGeom.Register<GeomIntersection>(CWP_GEOM_INTERSECTION);
  // factoryGeom.Register<GeomClosestPoint>(CWP_GEOM_CLOSEST_POINT);

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
    
  PDM_printf("CWP_Finalize\n");
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
 * \brief Number of codes known to CWIPI
 *
 * \return Number of codes
 *
 */

int
CWP_Codes_nb_get
(
)
{
  cwipi::CodePropertiesDB & properties =
    cwipi::CodePropertiesDB::getInstance();

  return properties.codesNbGet();
}


MPI_Comm  
CWP_Connectable_comm_get
(
  char* local_code_name
)
{
  cwipi::CodePropertiesDB & properties =
    cwipi::CodePropertiesDB::getInstance();


  const cwipi::CodeProperties & localCodeProperties = properties.codePropertiesGet(string(local_code_name));
  MPI_Comm connecComm = localCodeProperties.connectableCommGet();
  
  MPI_Comm globalComm = localCodeProperties.globalCommGet();
  int rank;
  MPI_Comm_rank(globalComm, &rank);
  
  //printf("CCCCCCC Après rank %i localColdeProperties %i \n",rank,localCodeProperties.idGet());

  return connecComm;
}


/**
 * \brief list of codes known to CWIPI
 *
 * \return Names list of codes
 *
 */

const char **
CWP_Codes_list_get
(
)
{
  cwipi::CodePropertiesDB & properties =
    cwipi::CodePropertiesDB::getInstance();
  return properties.codesListGet();
}


/**
 * \brief Number of codes known to CWIPI
 *
 * \return Number of local codes
 *
 */

int
CWP_Loc_codes_nb_get
(
)
{
  cwipi::CodePropertiesDB & properties =
    cwipi::CodePropertiesDB::getInstance();

  return properties.localCodesNbGet();
}


/**
 * \brief list of codes known to CWIPI
 *
 * \return Names list of local codes
 *
 */

const char **
CWP_Loc_codes_list_get
(
)
{
  cwipi::CodePropertiesDB & properties =
    cwipi::CodePropertiesDB::getInstance();
  return properties.localCodesListGet();
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
 * \brief Dump application properties.
 *
 * This function dump application properties.
 *
 */

void 
CWP_Properties_dump
(
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
 * \param [in]  local_code_name     Local code name
 * \param [in]  cpl_id              Coupling identifier
 * \param [in]  coupled_code_name   Distant or local coupled code name
 * \param [in]  comm_type           Communication type
 * \param [in]  geom_algo           Geometric algorithm
 * \param [in]  n_part              Number of interface partition 
 * \param [in]  displacement        Mesh displacement
 * \param [in]  recv_freq_type      Type of receiving frequency
 *
 */

void 
CWP_Cpl_create
(
 const char               *local_code_name,
 const char               *cpl_id,
 const char               *coupled_code_name,
 const CWP_Comm_t          comm_type, 
 const CWP_Geom_t          geom_algo,
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
  

 // printf("couplingDB.couplingCreate(proper\n");
  couplingDB.couplingCreate(properties.codePropertiesGet(local_application_str),
                            coupling_name_str,
                            properties.codePropertiesGet(coupled_application_str),
                            comm_type,
                            geom_algo,
                            n_part,
                            displacement,
                            recv_freq_type);
}

/**
 *
 * \brief Removing a coupling object
 *
 * This function delete a coupling object
 * 
 * \param [in] local_code_name   Local code name
 * \param [in] cpl_id            Coupling identifier
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



// void
// CWP_Exch
// (
//  const char *cpl_id
// )
// {
//   //TODO: Voir comment enchainer les appels, voir comment prendre en compte l'interpolation temporelle 
//   cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
//   cpl.exchange();
// }



int 
CWP_N_uncomputed_tgts_get
(
 const char *local_code_name,
 const char *cpl_id,
 const CWP_Field_value_t geometryLocation,
 const int  i_part
)
{
   cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);

   return cpl.nUncomputedTargetsGet(geometryLocation,i_part);
}



// const int *
// CWP_uncomputed_tgts_get
// (
//  const char *cpl_id
// )
// {
//   cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
//   return cpl.uncomputedTargetsGet();
// }


// int 
// CWP_n_computed_tgts_get
// (
//  const char *cpl_id
// )
// {
//   cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
//   return cpl.nComputedTargetsGet();
// }



// const int *
// CWP_computed_tgts_get
// (
//  const char *cpl_id
// )
// {
//   cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
//   return cpl.getLocatedPoint();
// }



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


 void 
 CWP_next_recv_time_set
 (const char     *local_code_name,
  const char                 *cpl_id,
  const double                next_time
 )
 {
   cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
   cpl.recvNextTimeSet(next_time);
 }

/*----------------------------------------------------------------------------*
 * Functions about geometry                                                   *
 *----------------------------------------------------------------------------*/


  void 
  CWP_Geom_compute(
                   const char        *cpl_id
                  )
  { 
    cwipi::CodePropertiesDB & properties =
    cwipi::CodePropertiesDB::getInstance();

    /* Get local codes list */
    int nb_loc_codes = properties.localCodesNbGet();
    const char** list_loc_codes = properties.localCodesListGet(); 

    /* Find the Geometry Location and exchange type of the fields DB. Which geomlocation and field_exch_type has to be computed ? */
    std::vector< std::map <CWP_Field_value_t, CWP_Field_exch_t> >* toComputeV = new std::vector< std::map <CWP_Field_value_t, CWP_Field_exch_t>  >(0);

    for(int i_codes=0;i_codes<nb_loc_codes;i_codes++){
      const char* loc_code_name = list_loc_codes[i_codes];
      
      /* Only if the coupling exists */
      if( _cpl_exist(loc_code_name,cpl_id) ) {
        cwipi::Coupling& cpl = _cpl_get(loc_code_name,cpl_id);
        std:map <std::string, cwipi::Field *>* fields = cpl.fieldsGet();
        std::map <CWP_Field_value_t, CWP_Field_exch_t> toCompute;

        std::map <std::string, cwipi::Field *>::iterator it = fields -> begin();
        while(it != fields -> end()){
          cwipi::Field* field = it -> second;
          CWP_Field_value_t fieldLocation = field -> typeGet();
          CWP_Field_exch_t  exchangeType  = field -> exchangeTypeGet();
       
          std::map <CWP_Field_value_t, CWP_Field_exch_t>::iterator finder = toCompute.find(fieldLocation);
          if(finder == toCompute.end()) {
            toCompute.insert({fieldLocation,exchangeType});
          }
          else {
            if(finder->second != exchangeType && finder->second != CWP_FIELD_EXCH_SENDRECV)
              toCompute.insert({fieldLocation,CWP_FIELD_EXCH_SENDRECV});
          }
          it++;
        }
        toComputeV -> push_back(toCompute);
      }//end if cpl_exist
    } //end for i_codes loop


    /*  Building of possible location vector */
    std::vector<CWP_Field_value_t> locationV;
    locationV.push_back(CWP_FIELD_VALUE_NODE);
    locationV.push_back(CWP_FIELD_VALUE_CELL_POINT);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    for(int i_location=0; i_location < locationV.size();i_location++){
      CWP_Field_value_t geometryLocation = locationV[i_location];
      for(int i_codes=0;i_codes<nb_loc_codes;i_codes++){   
        const char* loc_code_name = list_loc_codes[i_codes];

        /* Only if the coupling exists */
        if( _cpl_exist(loc_code_name,cpl_id) ) {
          cwipi::Coupling& cpl = _cpl_get(loc_code_name,cpl_id);
          std::map <CWP_Field_value_t, CWP_Field_exch_t> toCompute = (*toComputeV)[i_codes];
          std::map <CWP_Field_value_t, CWP_Field_exch_t>::iterator itToCompute = toCompute.begin();
          while(itToCompute != toCompute.end()){
            CWP_Field_value_t geometryLocationTest = itToCompute -> first;
            CWP_Field_exch_t exchange_type = itToCompute -> second;
            
            if(geometryLocationTest == geometryLocation) {
              if(geometryLocation == CWP_FIELD_VALUE_NODE) {
              if(exchange_type == CWP_FIELD_EXCH_SENDRECV )
                printf("localCodes %s rank %i CWP_FIELD_EXCH_SENDRECV CWP_FIELD_VALUE_NODE\n",loc_code_name,rank);
              else if(exchange_type == CWP_FIELD_EXCH_SEND)
                printf("localCodes %s rank %i CWP_FIELD_EXCH_SEND CWP_FIELD_VALUE_NODE\n",loc_code_name,rank);            
              else if(exchange_type == CWP_FIELD_EXCH_RECV)
                printf("localCodes %s rank %i CWP_FIELD_EXCH_RECV CWP_FIELD_VALUE_NODE\n",loc_code_name,rank);    
              }
              else if(geometryLocation == CWP_FIELD_VALUE_CELL_POINT) {
              if(exchange_type == CWP_FIELD_EXCH_SENDRECV )
                printf("localCodes %s rank %i CWP_FIELD_EXCH_SENDRECV CWP_FIELD_VALUE_CELL_POINT\n",loc_code_name,rank);
              else if(exchange_type == CWP_FIELD_EXCH_SEND)
                printf("localCodes %s rank %i CWP_FIELD_EXCH_SEND CWP_FIELD_VALUE_CELL_POINT\n",loc_code_name,rank);            
              else if(exchange_type == CWP_FIELD_EXCH_RECV)
                printf("localCodes %s rank %i CWP_FIELD_EXCH_RECV CWP_FIELD_VALUE_CELL_POINT\n",loc_code_name,rank);    
              }
              
              
              cpl.geomCompute(geometryLocation, exchange_type);
            }
            itToCompute++;
          }
        }//end if exist    
      } //end on i_codes loop  
      printf("rank %i End compute loop\n",rank);    
    } //end on location loop
    delete toComputeV;
  }



// void 
// CWP_geom_update
// (
//  const char     *cpl_id,
//  const int       storage_id
// )
// {
// TODO
// }


// void 
// CWP_geom_properties_set
// (
//  const char     *cpl_id,
//  const char     *fmt,
//                   ...
// )
// {
//   cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);

//   va_list pa;
//   va_start(pa, fmt);
//   cpl.geomPropertiesSet(fmt, &pa);
//   va_end(pa);
// }

/*----------------------------------------------------------------------------*
 * Functions about visualization                                              *
 *----------------------------------------------------------------------------*/


 void 
 CWP_Visu_set
 (const char                 *local_code_name,
  const char                 *cpl_id,
  const int                   freq,
  const CWP_Visu_format_t     format,
  const char                 *format_option
 )
 {
   cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
   cpl.visuSet(freq, format, format_option);
 }

/*----------------------------------------------------------------------------*
 * Functions about User target points                                         *
 *----------------------------------------------------------------------------*/


// void 
// CWP_user_tgt_pts_set
// (
//  const char                 *cpl_id,
//  const int                   n_pts,
//  double                      coord[]
// )
// {
//   cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
//   cpl.userTgtPtsSet(n_pts, coord);
// }

/*----------------------------------------------------------------------------*
 * Functions about Mesh                                                    *
 *----------------------------------------------------------------------------*/


  void 
  CWP_Mesh_interf_finalize 
  (const char                 *local_code_name,
   const char                 *cpl_id
  ) 
  {
    cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
    cpl.meshFinalize();
  }


 void 
 CWP_Mesh_interf_vtx_set
 (const char                 *local_code_name,
  const char                 *cpl_id,
  const int                   i_part,
  const int                   n_pts,
  double                      coord[],
  CWP_g_num_t                 global_num[]
 )
 {
   cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
   cpl.meshVtcsSet(i_part,
                   n_pts,
                   coord,
                   global_num);
 }

 int 
 CWP_Mesh_interf_block_add
 (
  const char           *local_code_name,
  const char           *cpl_id,
  const CWP_Block_t     block_type
 )
 {
   cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
   int block_id = cpl.meshBlockAdd(block_type);
   return block_id;
 }

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
   cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
   cpl.meshStdBlockSet(i_part,
                       block_id,
                       n_elts,
                       connec,
                       global_num);
 }


void 
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
 /* cpl.meshHighOrderBlockSet(i_part,
                            block_id,
                            block_type,
                            n_elts,
                            order,
                            connec,
                            global_num);*/
}


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
   cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
   cpl.meshFPolyBlockSet(i_part,
                         block_id,
                         n_elts,
                         connec_idx,
                         connec,
                         global_num);

}


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

/*----------------------------------------------------------------------------*
 * Functions about field                                                      *
 *----------------------------------------------------------------------------*/

 void
 CWP_Field_create
 (
  const char                  *codeName,
  const char                  *cpl_id,
  const char                  *field_id,
  const CWP_Type_t             data_type,
  const CWP_Field_storage_t    storage,
  const int                    n_component,
  const CWP_Field_value_t      value_target_proc_part_num,
  const CWP_Field_exch_t       exchange_type,
  const CWP_Status_t           visu_status)
 {
   cwipi::Coupling& cpl = _cpl_get(codeName,cpl_id);
   
   cpl.fieldCreate(field_id,
                   data_type,
                   storage,
                   n_component,
                   value_target_proc_part_num,
                   exchange_type, 
                   visu_status);
 }


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
   cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
   
   cpl.fieldDataSet(field_id,i_part,data);   
 }


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


 CWP_Field_value_t
 CWP_Field_target_proc_part_num_get
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


// CWP_error_t 
// CWP_sendrecv
// (
//  const char   *cpl_id,
//  const char   *src_field_id,
//  const char   *tgt_field_id,
//  int          *n_uncomputed_tgt
// )
// {
//   cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
//   const string &src_field_id_str = src_field_id;
//   const string &tgt_field_id_str = tgt_field_id;

//   return cpl.sendRecv(src_field_id_str,
//                       tgt_field_id_str,
//                       NULL,
//                       n_uncomputed_tgt);
// }



 void 
 CWP_Issend
 (const char        *local_code_name,
  const char        *cpl_id,
  const char        *referenceFieldID)
 {
   cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
   std::string referenceFieldID_str = const_cast<char*>(referenceFieldID);

   cpl.issend(referenceFieldID_str);
   
 }


 void 
 CWP_Irecv
 (const char        *local_code_name,
  const char        *cpl_id,
  const char        *targetFieldID)
 {
   cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);
   std::string targetFieldID_str = const_cast<char*>(targetFieldID);
   
   cpl.irecv(targetFieldID_str);
         
 }


 void 
 CWP_Wait_issend
 (
  const char  *local_code_name,
  const char  *cpl_id,
  const char  *src_field_id
 )
 {
   cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);

   std::string srcFieldID_str = const_cast<char*>(src_field_id);

   cpl.waitIssend(srcFieldID_str);
   
 }


 void 
 CWP_Wait_irecv
 (
  const char  *local_code_name,
  const char  *cpl_id,
  const char  *distant_field_id
 )
 {
   cwipi::Coupling& cpl = _cpl_get(local_code_name,cpl_id);

   MPI_Comm intraComm = cpl.localCodePropertiesGet()->intraCommGet();
   
   int size_comm = -1;
   MPI_Comm_size(intraComm,&size_comm);

   cwipi::CodeProperties* localCodeProperties = cpl.localCodePropertiesGet();
   MPI_Comm globalComm = localCodeProperties -> globalCommGet();
   
   int rank;
   MPI_Comm_rank(globalComm,&rank); 

   std::string distantFieldID_str = const_cast<char*>(distant_field_id);

   cpl.waitIrecv(distantFieldID_str);

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


/*----------------------------------------------------------------------------*
 * Functions about user interpolation                                         *
 *----------------------------------------------------------------------------*/


// void 
// CWP_Interp_from_loc_set
// (
//  const char                  *cpl_id,
//  CWP_Interp_from_target_proc_part_num_t fct
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
 * \brief Add a control parameter
 * 
 * Addition of a control parameter in the code properties.
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
    PDM_error(__FILE__, __LINE__, 0,
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
    PDM_error(__FILE__, __LINE__, 0,
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
    PDM_error(__FILE__, __LINE__, 0,
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
    PDM_error(__FILE__, __LINE__, 0,
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
    PDM_error(__FILE__, __LINE__, 0,
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
    PDM_error(__FILE__, __LINE__, 0,
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
    PDM_error(__FILE__, __LINE__, 0,
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
