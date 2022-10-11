#ifndef __SERVER_H__
#define __SERVER_H__
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
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#ifdef WINDOWS

#else
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
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
/* server status */
#define CWP_SVRSTATE_WAITCONN      0
#define CWP_SVRSTATE_LISTENINGMSG  1
#define CWP_SVRSTATE_RECVPPUTDATA  2
#define CWP_SVRSTATE_SENDPGETDATA  3
#define CWP_SVRSTATE_TERMINATING   4

/* for debug */
#define CWP_SVRFLAG_VERBOSE    1
#define CWP_SVRFLAG_NOVERBOSE  0

/*============================================================================
 * Types definition
 *============================================================================*/

typedef struct t_server
{
  int port;
  int state;
  int flags;
  int max_msg_size;
  int listen_socket;
  int connected_socket;
  int client_endianess;
  int server_endianess;
  char host_name[256];
}t_server,*p_server;

/*=============================================================================
 * Public function interfaces
 *============================================================================*/

/* Create a server */

int
CWP_server_create
(
 int server_port,
 int flags,
 p_server svr
);

/* Kill a server */

int
CWP_server_kill
(
 p_server svr
);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __SERVER_H__ */
