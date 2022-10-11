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
 * Public function interfaces
 *============================================================================*/

/* Connect to a server */

int
CWP_Connect
(
 const char* server_name,
 int server_port,
 int flags,
 p_client clt
);

/* Disconnect */

int
CWP_Disconnect
(
 p_client clt
);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __CLIENT_H__ */
