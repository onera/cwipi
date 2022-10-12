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
#include "cwp.h"
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

/*=============================================================================
 * Private function interfaces
 *============================================================================*/

// --> wrapper

void read_name(char *name) {
  int name_size;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &name_size, sizeof(int));
  name = malloc(name_size);
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) name, code_name_size);
}

/*=============================================================================
 * Server CWIPI function interfaces
 *============================================================================*/

void
CWP_server_Init
(
  p_server                 svr
)
{
  int      n_code;
  int      code_name_size;
  char         **code_names = NULL;
  CWP_Status_t  *is_active_rank = NULL;
  double        *time_init = NULL;

  // receive data
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  CWP_transfer_readdata(svr->connected_socket, svr->max_msg_size, &n_code, sizeof(int));
  code_names = malloc(sizeof(char *) * n_code);
  for (int i = 0; i < n_code; i++) {
    code_names[i] = NULL;
    read_name(code_names[i]);
  }
  is_active_rank = malloc(sizeof(int) * n_code);
  CWP_transfer_readdata(svr->connected_socket, svr->max_msg_size, is_active_rank, n_code * sizeof(int));
  time_init = malloc(sizeof(double) * n_code);
  CWP_transfer_readdata(svr->connected_socket, svr->max_msg_size, time_init, n_code * sizeof(double));

  // launch CWP_Init
  svr->intra_comms = malloc(sizeof(MPI_Comm) * n_code);
  CWP_Init(svr->global_comm,
           n_code,
           (const char **) code_names,
           is_active_rank,
           time_init,
           intra_comms); // stocker


  // verbose
  if (svr->flags & CWP_SVRFLAG_VERBOSE) {
    CWP_Properties_dump();
  }

  // free
  for (int i = 0; i < n_code; i++) {
    free(code_names[i]);
  }
  free(code_names);
  free(is_active_rank);
  free(time_init);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Finalize()
{
  CWP_Finalize();
}

void
CWP_server_Param_lock
(
 p_server                 svr
)
{
  char *code_name = NULL;
  read_name(code_name);

  CWP_Param_lock((const char *) code_name);

}

void
CWP_server_Param_unlock
(
 p_server                 svr
)
{
  char *code_name = NULL;
  read_name(code_name);

  CWP_Param_unlock((const char *) code_name);

}

void
CWP_server_Param_add
(
 p_server                 svr
)
{
  // read local code name
  char *local_code_name = NULL;
  read_name(local_code_name);

  // read param name
  char *param_name = NULL;
  read_name(param_name);

  // read initial value
  void initial_value;
  int data_type = -1;
  CWP_transfer_readdata(svr->connected_socket, svr->max_msg_size, &data_type, sizeof(int));

  switch (data_type) {

  case CWP_DOUBLE:
    double double_initial_value;
    CWP_transfer_readdata(svr->connected_socket, svr->max_msg_size, &double_initial_value, sizeof(double));
    initial_value = (void) double_initial_value;
    break;

  case CWP_INT:
    double int_initial_value;
    CWP_transfer_readdata(svr->connected_socket, svr->max_msg_size, &int_initial_value, sizeof(int));
    initial_value = (void) int_initial_value;
    break;

  case CWP_CHAR:
    char *char_initial_value = NULL;
    read_name(char_initial_value);
    initial_value = (void) char_initial_value;
    break;

  default:
    PDM_error(__FILE__, __LINE__, 0, "Received unknown CWP_Type_t %i\n", data_type);
  }

  CWP_Param_add(local_code_name,
                param_name,
                data_type,
                &initial_value);
}

/*============================================================================
 * Server function definitions
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

/* Message handler */

int
CWP_server_msg_handler
(
 p_server svr,
 p_message msg
)
{
  // verbose
  if (svr->flags & CWP_SVRFLAG_VERBOSE) {
    log_trace("CWP: Received msg->message_type: %d\n", msg->message_type);
  }

  switch(msg->message_type) {

  case CWP_MSG_DIE:
    svr->state=CWP_SVRSTATE_TERMINATING;

    // verbose
    if (svr->flags & CWP_SVRFLAG_VERBOSE) {
      log_trace("CWP: server recieved termination signal\n");
    }

    break;

  case CWP_MSG_CWP_INIT:

    // verbose
    if (svr->flags & CWP_SVRFLAG_VERBOSE) {
      log_trace("CWP: server received CWP_Init signal\n");
    }

    // launch
    CWP_server_Init(svr);

    break;

  case CWP_MSG_CWP_FINALIZE:

    // verbose
    if (svr->flags & CWP_SVRFLAG_VERBOSE) {
      log_trace("CWP: server received CWP_Finalize signal\n");
    }

    // launch
    CWP_server_Finalize();

    break;

  case CWP_MSG_CWP_PARAM_LOCK:

    // verbose
    if (svr->flags & CWP_SVRFLAG_VERBOSE) {
      log_trace("CWP: server received CWP_Param_lock signal\n");
    }

    // launch
    CWP_server_Param_lock(svr);

    break;

  case CWP_MSG_CWP_PARAM_UNLOCK:

    // verbose
    if (svr->flags & CWP_SVRFLAG_VERBOSE) {
      log_trace("CWP: server received CWP_Param_unlock signal\n");
    }

    // launch
    CWP_server_Param_unlock(svr);

    break;

  case CWP_MSG_CWP_PARAM_ADD:

    // verbose
    if (svr->flags & CWP_SVRFLAG_VERBOSE) {
      log_trace("CWP: server received CWP_Param_add signal\n");
    }

    // launch
    CWP_server_Param_add(svr);

    break;

  default:
    PDM_error(__FILE__, __LINE__, 0, "Received unknown message type %i, terminating server\n", msg->message_type);
    return -1;
  }

  return 0;
}

/* Run a server */

int
CWP_server_run
(
 p_server svr
)
{
  struct sockaddr_in client_addr;
  socklen_t sin_size;
  int il_sv_endian;
  int il_cl_endian;

  // accept client connexion
  svr->connected_socket = accept(svr->listen_socket, (struct sockaddr *)&client_addr,&sin_size);

  // verbose
  if (svr->flags & CWP_SVRFLAG_VERBOSE) {
    log_trace("CWP:got a connection from %s:%d\n",
    inet_ntoa(client_addr.sin_addr),ntohs(client_addr.sin_port));
  }

  // endianess
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  if(CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,
           &il_cl_endian,sizeof(int))!=0) {
    PDM_error(__FILE__, __LINE__, 0, "Client endian read failed\n");
    svr->state=CWP_SVRSTATE_LISTENINGMSG;
    return -1;
  }
  il_cl_endian =  ntohl(il_cl_endian);

  svr->client_endianess = il_cl_endian;

  il_sv_endian = htonl(svr->server_endianess);

  if(CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size,
           &il_sv_endian,sizeof(int))!=0) {
    PDM_error(__FILE__, __LINE__, 0, "Server endian send failed\n");
    svr->state=CWP_SVRSTATE_LISTENINGMSG;
    return -1;
  }

  svr->state=CWP_SVRSTATE_LISTENINGMSG;

  if (svr->flags & CWP_SVRFLAG_VERBOSE) {
    log_trace("Server : client endian %i server indian %i\n",svr->client_endianess,svr->server_endianess);
  }

  t_message msg;

  while (svr->state != CWP_SVRSTATE_TERMINATING) {

    CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,&msg, sizeof(t_message));

    if (CWP_server_msg_handler(svr,&msg) != 0) {
      PDM_error(__FILE__, __LINE__, 0, "Server message handling failed\n");
      return -1;
    }

  }

  return 0;
}


#ifdef __cplusplus
}
#endif /* __cplusplus */
