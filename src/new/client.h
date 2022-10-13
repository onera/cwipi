#ifndef __CLIENT_H__
#define __CLIENT_H__
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
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#ifdef WINDOWS

#else
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>
#include <arpa/inet.h>
#include <unistd.h>
#include <errno.h>
#endif

#include "cwp.h"


#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

#define CWP_CLIENTFLAG_VERBOSE    1
#define CWP_CLIENTFLAG_NOVERBOSE  0

/*============================================================================
 * Types definition
 *============================================================================*/

typedef struct t_client
{
  int server_port;
  int flags;
  int socket;
  int max_msg_size;
  int listen_socket;
  int connected_socket;
  int client_endianess;
  int server_endianess;
  char server_name[256];

}t_client,*p_client;

/*=============================================================================
 * Client CWIPI function interfaces
 *============================================================================*/

/**
 * \brief Initialize CWIPI.
 *
 * \param [in]  n_code         Number of codes on the current rank
 * \param [in]  code_names     Names of codes on the current rank (size = \p n_code)
 * \param [in]  is_active_rank Is current rank have to be used by CWIPI (size = \p n_code)
 * \param [in]  time_init      Initial time (size = \p n_code)
 * \param [out] intra_comms    MPI intra communicators of each code (size = \p n_code)
 *
 */

void
CWP_client_Init
(
  const int                n_code,
  const char             **code_names,
  const CWP_Status_t      *is_active_rank,
  const double            *time_init
);

/**
 *
 * \brief Finalize CWIPI.
 *
 */

void
CWP_client_Finalize
(
 void
);

/**
 *
 * \brief Param_lock CWIPI.
 *
 * \param [in]  code_name  Code to lock
 *
 */

void
CWP_client_Param_lock
(
const char *code_name
);

/**
 *
 * \brief Param_unlock CWIPI.
 *
 * \param [in]  code_name  Code to unlock
 *
 */


void
CWP_client_Param_unlock
(
const char *code_name
);

/**
 *
 * \brief Param_add CWIPI.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] param_name       Parameter name
 * \param [in] data_type        Parameter type
 * \param [in] initial_value    Initial value
 *
 */

void
CWP_client_Param_add
(
 const char        *local_code_name,
 const char        *param_name,
 const CWP_Type_t  data_type,
 void              *initial_value
);

/**
 *
 * \brief Param_get CWIPI.
 *
 * \param [in]  code_name  Local or distant code name
 * \param [in]  param_name Parameter name
 * \param [in]  data_type  Parameter type
 * \param [out] value      Parameter value
 *
 */

void
CWP_client_Param_get
(
 const char       *code_name,
 const char       *param_name,
 const CWP_Type_t  data_type,
 void             *value
);

/**
 * \brief Cpl_create CWIPI.
 *
 * \param [in]  local_code_name     Local code name
 * \param [in]  cpl_id              Coupling identifier
 * \param [in]  coupled_code_name   Distant or local coupled code name
 * \param [in]  comm_type           Communication type
 * \param [in]  spatial_interp      Spatial interpolation method
 * \param [in]  n_part              Number of interface partition
 * \param [in]  displacement        Mesh moving status
 * \param [in]  recv_freq_type      Type of receiving frequency
 *
 */

void
CWP_client_Cpl_create
(
 const char                *local_code_name,
 const char                *cpl_id,
 const char                *coupled_code_name,
 CWP_Interface_t            entities_dim,
 const CWP_Comm_t           comm_type,
 const CWP_Spatial_interp_t spatial_interp,
 const int                  n_part,
 const CWP_Dynamic_mesh_t   displacement,
 const CWP_Time_exch_t      recv_freq_type
);

/**
 *
 * \brief Cpl_del CWIPI.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 *
 */

void
CWP_client_Cpl_del
(
 const char *local_code_name,
 const char *cpl_id
);


/*=============================================================================
 * Client function interfaces
 *============================================================================*/

/* Connect to a server */

int
CWP_client_connect
(
 const char* server_name,
 int server_port,
 int flags
);

/* Disconnect */

int
CWP_client_disconnect
(
 void
);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __CLIENT_H__ */
