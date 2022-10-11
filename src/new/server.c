/*
  This file is part of the CWIPI library.

  Copyright (C) 2017  ONERA

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
 *  System headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "server.h"
#include "message.h"
#include "transfer.h"
#include <pdm_error.h>
#include <pdm_mpi.h>
#include "pdm_logging.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/* Create a server */

int
CWP_server_create
(
 int server_port,
 int flags,
 p_server svr
)
{
  struct sockaddr_in server_addr;
  socklen_t d;
  int true=1;

  memset(svr,0,sizeof(t_server));
  svr->port             = server_port;
  svr->flags            = flags;
  svr->server_endianess = CWP_transfer_endian_machine();

  // retreive hostname
  if(gethostname(svr->host_name,sizeof(svr->host_name)) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "Could not get host name\n");
    return -1;
  }

  // verbose
  if (svr->flags & CWP_SVRFLAG_VERBOSE) {
    log_trace("CWP:Creating Server on %s port %i...\n", svr->host_name, svr->port);
  }

  // create socket (IPv4, binary data flux, unique protocol for AF_INET + SOCK_STREAM)
  svr->listen_socket = socket(AF_INET, SOCK_STREAM, 0);

  if (svr->listen_socket == -1) {
    PDM_error(__FILE__, __LINE__, 0, "Could not create socket to listen\n");
    return -1;
  }

  // set maximum message size
  getsockopt(svr->listen_socket, SOL_SOCKET, SO_RCVBUF, (char*)&svr->max_msg_size, &d);
  svr->max_msg_size = CWP_MSG_MAXMSGSIZE;

  // verbose
  if (svr->flags & CWP_SVRFLAG_VERBOSE) {
    log_trace("CWP:Max message size:%i\n",svr->max_msg_size);
  }

  // ensure restarted process can reconnect to same socket
  if (setsockopt(svr->listen_socket,SOL_SOCKET,SO_REUSEADDR,&true,sizeof(int)) == -1) {
    PDM_error(__FILE__, __LINE__, 0, "Setsockopt failed\n");
    return -1;
  }

  // fill data structure for ip adress + port
  server_addr.sin_family      = AF_INET;
  server_addr.sin_port        = htons(svr->port);
  server_addr.sin_addr.s_addr = INADDR_ANY;
  bzero(&(server_addr.sin_zero),8);

  // bind socket to ip adress + port
  if (bind(svr->listen_socket, (struct sockaddr *)&server_addr, sizeof(struct sockaddr)) == -1) {
    PDM_error(__FILE__, __LINE__, 0, "Unable to bind socket\n");
    return -1;
  }

  // wait for client to connect
  if (listen(svr->listen_socket, SOMAXCONN) == -1) {
    PDM_error(__FILE__, __LINE__, 0, "Unable to listen on socket\n");
    return -1;
  }

  svr->state = CWP_SVRSTATE_WAITCONN;

  // verbose
  if (svr->flags & CWP_SVRFLAG_VERBOSE) {
    log_trace("CWP:Server created on %s port %i\n", svr->host_name, svr->port);
  }

  return 0;
}

/* Kill a server */

int
CWP_server_kill
(
 p_server svr
)
{
  // verbose
  if (svr->flags & CWP_SVRFLAG_VERBOSE) {
    log_trace("CWP:Server shutting down\n");
  }

  // shutdown
#ifdef WINDOWS
  shutdown(svr->listen_socket,SD_BOTH);
  shutdown(svr->connected_socket,SD_BOTH);
#else
  shutdown(svr->listen_socket,SHUT_RDWR);
  shutdown(svr->connected_socket,SHUT_RDWR);
#endif

  memset(svr,0,sizeof(t_server));

  return 0;
}


#ifdef __cplusplus
}
#endif /* __cplusplus */
