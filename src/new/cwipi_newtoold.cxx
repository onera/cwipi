
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


#include <map>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cwp.h"
#include "cwipi_newtoold.h"



/*=============================================================================
 * Public function prototypes
 *============================================================================*/

void cwipi_init_new_to_old
(const MPI_Comm                           global_comm,
 const char                               *code_name,
 MPI_Comm                                 *local_comm
)
{

   const int           n_code = 1;
   const CWP_Status_t  is_coupled_rank = CWP_STATUS_ON;
   const double        time_init = 0.;

   CWP_Init(global_comm,
            n_code,
            (const char **) &(code_name),
            &is_coupled_rank,
            &time_init,
            local_comm);

}



void cwipi_create_coupling_new_to_old
( const char  *coupling_name,
  const cwipi_coupling_type_t coupling_type,
  const char  *codeCoupledName,
  const int    entitiesDim,
  const double tolerance,
  const cwipi_mesh_type_t mesh_type,
  const cwipi_solver_type_t solver_type,
  const int    output_frequency,
  const char  *output_format,
  const char  *output_format_option,
  ...
)
{

    const char** codeNameList = CWP_Loc_codes_list_get();
    const int nb_part = 1;
    const CWP_Comm_t comm_type = coupling_type_conv.at(coupling_type);
    const char* local_name = codeNameList[0];

    CWP_Cpl_create (local_name     ,                        
                    coupling_name  ,                       
                    codeCoupledName,
                    (CWP_Interface_t) entitiesDim    ,
                    comm_type      ,
                    CWP_SPATIAL_INTERP_FROM_LOCATION_DIST_CLOUD_SURF,
                    nb_part,                         
                    mesh_type_conv.at(mesh_type), (CWP_Time_exch_t) 1);
    
    CWP_Visu_format_t format;

    if(output_format == "EnSight Gold"){
      format = CWP_VISU_FORMAT_ENSIGHT;
    } 
    else{
      format = CWP_VISU_FORMAT_ENSIGHT;
    }
         
    const int freq = output_frequency;
    CWP_Visu_set(local_name,
                 coupling_name  ,                       
                 freq           ,
                 format,
                 output_format_option
                );   

}




