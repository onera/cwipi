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
  Inspired from OpenPALM, a free software under the GNU Lesser General Public License
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
#include <pdm_error.h>
#include <pdm_mpi.h>

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
CWP_CreateServer
(
 int server_port,
 int flags,
 p_server svr
)
{
  struct sockaddr_in server_addr; // storage for IP adress and port
  socklen_t d; // length
  int true=1;

  memset(svr,0,sizeof(t_server));
  svr->port             = server_port;
  svr->flags            = flags;
  // svr->server_endianess = iplib_endian_machine(); // TO DO

  if(gethostname(svr->host_name,sizeof(svr->host_name)) != 0) { // IP adress
    log_trace("CWP:Warning could not get host name using loopback address 127.0.0.1");
    strcpy(svr->host_name, "127.0.0.1");
  }

  if (svr->flags & CWP_SVRFLAG_VERBOSE) {
    log_trace("CWP:Creating Server on %s port %i...\n", svr->host_name, svr->port);
  }

  // create communication point
  // AF_INET: processes on different hosts connected by IPv4
  // SOCK_STREAM: dialog support for binary data flux
  // 0: only one protocol type for combination of AF_INET + SOCK_STREAM
  svr->listen_socket = socket(AF_INET, SOCK_STREAM, 0);

  if (svr->listen_socket == -1) { // socket descriptor or -1 if failed
    PDM_error(__FILE__, __LINE__, 0, "Could not create socket to listen\n");
    return -1;
  }

  getsockopt(svr->listen_socket, SOL_SOCKET, SO_RCVBUF, (char*)&svr->max_msg_size, &d);
  // svr->max_msg_size = PALMONIP_MSG_MAXMSGSIZE;

  if (svr->flags & CWP_SVRFLAG_VERBOSE) {
    log_trace("CWP:Max message size:%i\n",svr->max_msg_size);
  }


  // if process would be killed SO_REUSEADDR ensure the restarted process can connect on the same socket (thus port) again
  if (setsockopt(svr->listen_socket,SOL_SOCKET,SO_REUSEADDR,&true,sizeof(int)) == -1) {
    PDM_error(__FILE__, __LINE__, 0, "Setsockopt failed\n");
    return -1;
  }

  server_addr.sin_family      = AF_INET;
  server_addr.sin_port        = htons(svr->port); // conversion from host to network byte order (short)
  server_addr.sin_addr.s_addr = INADDR_ANY;
  bzero(&(server_addr.sin_zero),8);

  // links sockets to IP adress and port
  if (bind(svr->listen_socket, (struct sockaddr *)&server_addr, sizeof(struct sockaddr)) == -1) {
    PDM_error(__FILE__, __LINE__, 0, "Unable to bind socket\n");
    return -1;
  }

  // waits for a client to connect
  if (listen(svr->listen_socket, SOMAXCONN) == -1) {
    PDM_error(__FILE__, __LINE__, 0, "Unable to listen on socket\n");
    return -1;
  }

  // svr->state = PALMONIP_SVRSTATE_WAITCONN;

  if (svr->flags & CWP_SVRFLAG_VERBOSE) {
    log_trace("CWP:Server created on %s port %i\n",svr->host_name,svr->port);
  }

  return 0;


}

/* Kill a server */

int
CWP_KillServer
(
 p_server svr
)
{

}


#ifdef __cplusplus
}
#endif /* __cplusplus */
