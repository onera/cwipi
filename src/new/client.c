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
#include "cwp.h"
#include <pdm_error.h>
#include <pdm_mpi.h>
#include "pdm_logging.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/* pointeur on client structure */

static t_client *clt;

/*=============================================================================
 * Private function interfaces
 *============================================================================*/

// --> wrapper

static void write_name(char * name) {
  int name_size = strlen(name);
  int endian_name_size = name_size;
  CWP_swap_endian_4bytes(&endian_name_size, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_name_size, sizeof(int));
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) name, name_size);
}

static void read_name(char *name) {
  int name_size;
  CWP_transfer_readdata(clt->socket,clt->max_msg_size,(void*) &name_size, sizeof(int));
  name = malloc(name_size);
  CWP_transfer_readdata(clt->socket,clt->max_msg_size,(void*) name, name_size);
}

// --> endianness

void  ip_swap_4bytes(char *f_bytes) {
  char a,b;
  if (clt->server_endianess == clt->client_endianess) {return;}
  a =  f_bytes[0];
  b =  f_bytes[1];
  f_bytes[0] = f_bytes[3];
  f_bytes[1] = f_bytes[2];
  f_bytes[2] = b;
  f_bytes[3] = a;
  return;
}

void  ip_swap_8bytes(char *cd_h_bytes) {
  char a,b,c,d;
  if (clt->server_endianess == clt->client_endianess) {return;}
  a =  cd_h_bytes[0];
  b =  cd_h_bytes[1];
  c =  cd_h_bytes[2];
  d =  cd_h_bytes[3];
  cd_h_bytes[0] = cd_h_bytes[7];
  cd_h_bytes[1] = cd_h_bytes[6];
  cd_h_bytes[2] = cd_h_bytes[5];
  cd_h_bytes[3] = cd_h_bytes[4];
  cd_h_bytes[4] = d;
  cd_h_bytes[5] = c;
  cd_h_bytes[6] = b;
  cd_h_bytes[7] = a;
  return;
}

void ip_swap_data_endian(char *data, const int datasize) {
  int i;
  if (clt->server_endianess == clt->client_endianess) {return;}

  for (i=0 ; i<datasize; i=i+4) {
    ip_swap_4bytes(&data[i]);
  }
  return;
}

/* convert message endian if client endianess different from server endianess */

int CWP_client_send_msg(p_message msg) {

  int data_size = sizeof(t_message);
  if (clt->server_endianess != clt->client_endianess) {
    ip_swap_4bytes((char *)&msg->message_type);
    ip_swap_4bytes((char *)&msg->flag);
    ip_swap_4bytes((char *)&msg->msg_tag);
    ip_swap_4bytes((char *)&msg->msg_time);
    ip_swap_4bytes((char *)&msg->data_size);
    ip_swap_4bytes((char *)&msg->data1);
    ip_swap_4bytes((char *)&msg->data2);
    ip_swap_4bytes((char *)&msg->data3);
  }
  return CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*)msg,data_size);
}

int CWP_swap_endian_4bytes(int *data,const  int datasize) {
  int i;
  for (i=0 ; i<datasize; i++) {
    ip_swap_4bytes((char*)&data[i]);
  }
  return 0;
}

int CWP_swap_endian_8bytes(double *data,const int datasize) {
  int i;
  for (i=0 ; i<datasize; i++) {
    ip_swap_8bytes((char*)&data[i]);
  }
  return 0;
}

/*=============================================================================
 * Client CWIPI function interfaces
 *============================================================================*/

void
CWP_client_Init
(
  const int                n_code,
  const char             **code_names,
  const CWP_Status_t      *is_active_rank,
  const double            *time_init
)
{
  t_message msg;

  // verbose
  if (clt->flags & CWP_CLIENTFLAG_VERBOSE) {
    log_trace("CWP:Client initiating CWP_Init\n");
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_INIT);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Init failed to send message header\n");
  }

  // endian swap
  int     endian_n_code         = n_code;
  int    *endian_is_active_rank = malloc(sizeof(int) * n_code);
  double *endian_time_init      = malloc(sizeof(double) * n_code);
  memcpy(endian_is_active_rank, is_active_rank, sizeof(int) * n_code);
  memcpy(endian_time_init, time_init, sizeof(double) * n_code);
  CWP_swap_endian_4bytes(&endian_n_code, 1);
  CWP_swap_endian_4bytes(endian_is_active_rank, 1);
  CWP_swap_endian_8bytes(endian_time_init, 1);

  // send arguments
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_n_code, sizeof(int));
  for (int i = 0; i < n_code; i++) {
    write_name(code_names[i]);
  }
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) endian_is_active_rank, n_code * sizeof(int));
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) endian_time_init, n_code * sizeof(double));

  // free
  free(endian_is_active_rank);
  free(endian_time_init);
}

void
CWP_client_Finalize()
{
  t_message msg;

  // verbose
  if (clt->flags & CWP_CLIENTFLAG_VERBOSE) {
    log_trace("CWP:Client initiating CWP_Finalize\n");
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_FINALIZE);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Finalize failed to send message header\n");
  }
}

void
CWP_client_Param_lock
(
const char *code_name
)
{
  t_message msg;

  // verbose
  if (clt->flags & CWP_CLIENTFLAG_VERBOSE) {
    log_trace("CWP:Client initiating CWP_Param_lock\n");
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_PARAM_LOCK);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Param_lock failed to send message header\n");
  }

  // send code name
  write_name(code_name);
}

void
CWP_client_Param_unlock
(
const char *code_name
)
{
  t_message msg;

  // verbose
  if (clt->flags & CWP_CLIENTFLAG_VERBOSE) {
    log_trace("CWP:Client initiating CWP_Param_unlock\n");
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_PARAM_UNLOCK);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Param_unlock failed to send message header\n");
  }

  // send code name
  write_name(code_name);
}

void
CWP_client_Param_add
(
 const char        *local_code_name,
 const char        *param_name,
 const CWP_Type_t  data_type,
 void              *initial_value
)
{
  t_message msg;

  // verbose
  if (clt->flags & CWP_CLIENTFLAG_VERBOSE) {
    log_trace("CWP:Client initiating CWP_Param_add\n");
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_PARAM_ADD);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Param_add failed to send message header\n");
  }

  // send local code name
  write_name(local_code_name);

  // send param name
  write_name(param_name);

  // send initial value
  int endian_data_type = data_type;
  CWP_swap_endian_4bytes(&endian_data_type, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_data_type, sizeof(int));

  switch (data_type) {

  case CWP_DOUBLE: ;
    double endian_double_initial_value = * ((double *) initial_value);
    CWP_swap_endian_8bytes(&endian_double_initial_value, 1);
    CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_double_initial_value, sizeof(double));
    break;

  case CWP_INT: ;
    int endian_int_initial_value = * ((int *) initial_value);
    CWP_swap_endian_4bytes(&endian_int_initial_value, 1);
    CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_int_initial_value, sizeof(int));
    break;

  case CWP_CHAR:
    write_name(initial_value);
    break;

  default:
    PDM_error(__FILE__, __LINE__, 0, "Unknown CWP_Type_t %i\n", data_type);
  }
}

void
CWP_client_Param_get
(
 const char       *code_name,
 const char       *param_name,
 const CWP_Type_t  data_type,
 void             *value
)
{
    t_message msg;

  // verbose
  if (clt->flags & CWP_CLIENTFLAG_VERBOSE) {
    log_trace("CWP:Client initiating CWP_Param_get\n");
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_PARAM_GET);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Param_get failed to send message header\n");
  }

  // send local code name
  write_name(code_name);

  // send param name
  write_name(param_name);

  // send value data type
  int endian_data_type = data_type;
  CWP_swap_endian_4bytes(&endian_data_type, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_data_type, sizeof(int));

  // receive value
  switch (data_type) {

  case CWP_DOUBLE:
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, value, sizeof(double));
    break;

  case CWP_INT:
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, value, sizeof(int));
    break;

  case CWP_CHAR: ;
    char *char_value = NULL;
    read_name(char_value);
    memcpy(value, char_value, strlen(char_value));
    break;

  default:
    PDM_error(__FILE__, __LINE__, 0, "Unknown CWP_Type_t %i\n", data_type);
  }

}

/*============================================================================
 * Client function definitions
 *============================================================================*/

/* Connect to a server */

int
CWP_client_connect
(
 const char* server_name,
 int server_port,
 int flags
)
{
  struct hostent *host;
  struct sockaddr_in server_addr;
  socklen_t d;
  int il_cl_endian;

  clt = malloc(sizeof(t_client));
  memset(clt,0,sizeof(t_client));
  strncpy(clt->server_name,server_name,sizeof(clt->server_name));
  clt->server_port=server_port;
  clt->flags=flags;

  // verbose
  if (clt->flags & CWP_CLIENTFLAG_VERBOSE) {
    log_trace("CWP:Creating Client, connecting to %s:%i...\n",server_name,server_port);
  }

  host = (struct hostent *) gethostbyname(server_name);

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

  // exchange endianess
  clt->client_endianess = CWP_transfer_endian_machine();
  il_cl_endian = htonl(clt->client_endianess);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,
         (void*)&il_cl_endian,sizeof(int));
  CWP_transfer_readdata(clt->socket,clt->max_msg_size,
           (void*)&clt->server_endianess,sizeof(int));
  clt->server_endianess = ntohl(clt->server_endianess);

  // verbose
  if (clt->flags & CWP_CLIENTFLAG_VERBOSE) {
    log_trace("CWP:client endian %i server indian %i\n",clt->client_endianess,clt->server_endianess);
  }

  return 0;

}

/* Disconnect */

int
CWP_client_disconnect
()
{
  t_message msg;

  if (clt->flags & CWP_CLIENTFLAG_VERBOSE) {
    log_trace("CWP:Client shutting down\n");
  }

  NEWMESSAGE(msg, CWP_MSG_DIE);

  CWP_client_send_msg(&msg);

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
