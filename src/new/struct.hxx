#ifndef __STRUCT_H__
#define __STRUCT_H__
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

/*
  This file is inspired from OpenPALM.
  OpenPALM is a free software under the GNU Lesser General Public License.
  See: https://www.cerfacs.fr/globc/PALM_WEB/
*/

/*----------------------------------------------------------------------------
 * Standard C++ library headers
 *----------------------------------------------------------------------------*/

#include <string>
#include <map>

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include "mpi.h"
#include "cwp.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Types definition
 *============================================================================*/

typedef struct t_code
{
  char                                **paramNames; // Param_list_get
  std::map<std::string, char *>         char_param_value;

  // std::map<std::string, p_server_cpl>   server_cpl;
} t_code, *p_code;

typedef struct t_cwp
{
  int          n_code_names;
  const char **code_names;     // Codes_list_get
  int          n_loc_code_names;
  const char **loc_code_names; // Loc_codes_list_get
  p_code       code;    // not implemented yet for n_code > 1
} t_cwp, *p_cwp;

typedef struct t_server_mpi
{
  MPI_Comm  global_comm;
  MPI_Comm *intra_comms;
} t_server_mpi, *p_server_mpi;

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __STRUCT_H__ */
