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

#include "client.h"
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

/* Connect to a server */

int
CWP_client_connect
(
 const char* server_name,
 int server_port,
 int flags,
 p_client clt
)
{
  struct hostent *host;
  struct sockaddr_in server_addr;
  socklen_t d;
  int il_cl_endian;

  memset(clt,0,sizeof(t_client));
  strncpy(clt->server_name,server_name,sizeof(clt->server_name));
  clt->server_port=server_port;
  clt->flags=flags;

  // verbose
  if (clt->flags & CWP_CLIENTFLAG_VERBOSE) {
    log_trace("CWP:Creating Client, connecting to %s:%i...\n",server_name,server_port);
  }

  host = (struct hostent *)gethostbyname(server_name);

  // create socket
  if ((clt->socket = socket(AF_INET, SOCK_STREAM, 0)) == -1) {
    PDM_error(__FILE__, __LINE__, 0, "Could not create Socket\n");
    return -1;
  }

  server_addr.sin_family = AF_INET;
  server_addr.sin_port = htons(server_port);
  server_addr.sin_addr = *((struct in_addr *)host->h_addr);
  bzero(&(server_addr.sin_zero),8);

  // connect
  while (connect(clt->socket, (struct sockaddr *)&server_addr,
                 sizeof(struct sockaddr)) == -1) {}

  // verbose
  if (clt->flags & CWP_CLIENTFLAG_VERBOSE) {
    log_trace("CWP:Client connected on %s port %i \n", server_name, server_port);
  }

  // set maximum message size
  getsockopt(clt->socket, SOL_SOCKET, SO_RCVBUF, (void*)&clt->max_msg_size, &d);
  clt->max_msg_size = CWP_MSG_MAXMSGSIZE;

  // exchange endianess // TO DO: create transfer.c and endianess
  clt->client_endianess = CWP_transfer_endian_machine();
  // il_cl_endian = htonl(clt->client_endianess);
  // transfer_writedata(clt->socket,clt->max_msg_size,
  //        (void*)&il_cl_endian,sizeof(int));
  // transfer_readdata(clt->socket,clt->max_msg_size,
  //          (void*)&clt->server_endianess,iLen);
  // clt->server_endianess = ntohl(clt->server_endianess);

  // verbose
  if (clt->flags & CWP_CLIENTFLAG_VERBOSE) {
    log_trace("CWP:client endian %i server indian %i\n",clt->client_endianess,clt->server_endianess);
  }

  return 0;

}

/* Disconnect */

int
CWP_client_disconnect
(
 p_client clt
)
{
  t_message msg;

  if (clt->flags & CWP_CLIENTFLAG_VERBOSE) {
    log_trace("CWP:Client shutting down\n");
  }

  // TO DO: send message to srv
  // NEWMESSAGE(msg,PALMONIP_MSG_DIE);

  Palm_On_Ip_ClientSendMsg(&msg);

  // shutdown
#ifdef WINDOWS
  shutdown(clt->socket,SD_BOTH);
#else
  shutdown(clt->socket,SHUT_RDWR);
#endif

  memset(clt,0,sizeof(t_client));

  return 0;
}

#ifdef __cplusplus
}
#endif /* __cplusplus */
