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
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>
#include <arpa/inet.h>
#include <unistd.h>
#include <errno.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "client.h"
#include "message.h"
#include "transfer.h"

#include "struct.hxx"

#include <cwp.h>
#include <pdm_error.h>
#include <pdm_mpi.h>
#include <pdm_logging.h>

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/* file struct definition */

static t_client *clt;
static t_cwp    *clt_cwp;

/*=============================================================================
 * Private function interfaces
 *============================================================================*/

// --> n_nodes

static int n_nodes_get(CWP_Block_t block_type) {
  int n_nodes = 0;

  switch(block_type) {

  case CWP_BLOCK_NODE: {
    n_nodes = 1;
    } break;

  case CWP_BLOCK_EDGE2: {
    n_nodes = 2;
    } break;

  case CWP_BLOCK_FACE_TRIA3: {
    n_nodes = 3;
    } break;

  case CWP_BLOCK_FACE_QUAD4: {
    n_nodes = 4;
    } break;

  case CWP_BLOCK_FACE_POLY: {
    // TO DO: what value ?
    PDM_error(__FILE__, __LINE__, 0, "Number of nodes unknown for CWP_BLOCK_FACE_POLY\n");
    return -1;
    } break;

  case CWP_BLOCK_CELL_TETRA4: {
    n_nodes = 4;
    } break;

  case CWP_BLOCK_CELL_HEXA8: {
    n_nodes = 8;
    } break;

  case CWP_BLOCK_CELL_PRISM6: {
    n_nodes = 6;
    } break;

  case CWP_BLOCK_CELL_PYRAM5: {
    n_nodes = 5;
    } break;

  case CWP_BLOCK_CELL_POLY: {
    // TO DO: what value ?
    PDM_error(__FILE__, __LINE__, 0, "Number of nodes unknown for CWP_BLOCK_CELL_POLY\n");
    return -1;
    } break;

  default:
    PDM_error(__FILE__, __LINE__, 0, "Unknown block type %d\n", block_type);
    return -1;
  }

  return n_nodes;
}

// --> endianness

static void ip_swap_4bytes(char *f_bytes) {
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

static void ip_swap_8bytes(char *cd_h_bytes) {
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

static void ip_swap_data_endian(char *data, const int datasize) {
  int i;
  if (clt->server_endianess == clt->client_endianess) {return;}

  for (i=0 ; i<datasize; i=i+4) {
    ip_swap_4bytes(&data[i]);
  }
  return;
}

/* convert message endian if client endianess different from server endianess */

static int CWP_client_send_msg(p_message msg) {

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

static int CWP_swap_endian_4bytes(int *data,const  int datasize) {
  int i;
  for (i=0 ; i<datasize; i++) {
    ip_swap_4bytes((char*)&data[i]);
  }
  return 0;
}

static int CWP_swap_endian_8bytes(double *data,const int datasize) {
  int i;
  for (i=0 ; i<datasize; i++) {
    ip_swap_8bytes((char*)&data[i]);
  }
  return 0;
}

// --> wrapper

static void write_name(const char * name) {
  int name_size = strlen(name) + 1;  // +1 for "\0"
  int endian_name_size = name_size;
  CWP_swap_endian_4bytes(&endian_name_size, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_name_size, sizeof(int));
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) name, name_size);
}

static void read_name(char **name) {
  int name_size;
  CWP_transfer_readdata(clt->socket,clt->max_msg_size,(void*) &name_size, sizeof(int));
  *name = (char *) realloc(*name, name_size);
  CWP_transfer_readdata(clt->socket,clt->max_msg_size,(void*) *name, name_size);
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
  int           endian_n_code         = n_code;
  CWP_Status_t *endian_is_active_rank = (CWP_Status_t *) malloc(sizeof(int) * n_code);
  double       *endian_time_init      = (double *) malloc(sizeof(double) * n_code);
  memcpy(endian_is_active_rank, is_active_rank, sizeof(CWP_Status_t) * n_code);
  memcpy(endian_time_init, time_init, sizeof(double) * n_code);
  CWP_swap_endian_4bytes(&endian_n_code, 1);
  CWP_swap_endian_4bytes((int *) endian_is_active_rank, n_code);
  CWP_swap_endian_8bytes(endian_time_init, n_code);

  // send arguments
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_n_code, sizeof(int));
  for (int i = 0; i < n_code; i++) {
    write_name(code_names[i]);
  }
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) endian_is_active_rank, n_code * sizeof(CWP_Status_t));
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
  CWP_Type_t endian_data_type = data_type;
  CWP_swap_endian_4bytes((int *) &endian_data_type, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_data_type, sizeof(CWP_Type_t));

  switch (data_type) {

  case CWP_DOUBLE: {
    double endian_double_initial_value = * ((double *) initial_value);
    CWP_swap_endian_8bytes(&endian_double_initial_value, 1);
    CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_double_initial_value, sizeof(double));
    } break;

  case CWP_INT: {
    int endian_int_initial_value = * ((int *) initial_value);
    CWP_swap_endian_4bytes(&endian_int_initial_value, 1);
    CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_int_initial_value, sizeof(int));
    } break;

  case CWP_CHAR: {
    char **p_char_value = (char **) initial_value;
    write_name(*p_char_value);
    } break;

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

  // send code name
  write_name(code_name);

  // send param name
  write_name(param_name);

  // send value data type
  CWP_Type_t endian_data_type = data_type;
  CWP_swap_endian_4bytes((int *) &endian_data_type, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_data_type, sizeof(CWP_Type_t));

  // receive value
  switch (data_type) {

  case CWP_DOUBLE: {
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, value, sizeof(double));
    } break;

  case CWP_INT: {
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, value, sizeof(int));
    } break;

  case CWP_CHAR: {
    int name_size;
    CWP_transfer_readdata(clt->socket,clt->max_msg_size,(void*) &name_size, sizeof(int));
    * (char **) value = (char *) malloc(sizeof(char) * name_size);
    CWP_transfer_readdata(clt->socket,clt->max_msg_size,(void*) * (char **) value, name_size);
    } break;

  default:
    PDM_error(__FILE__, __LINE__, 0, "Unknown CWP_Type_t %i\n", data_type);
  }

}

void
CWP_client_Param_set
(
 const char             *local_code_name,
 const char             *param_name,
 const CWP_Type_t        data_type,
 void                   *value
)
{
  t_message msg;

  // verbose
  if (clt->flags & CWP_CLIENTFLAG_VERBOSE) {
    log_trace("CWP:Client initiating CWP_Param_set\n");
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_PARAM_SET );

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Param_set failed to send message header\n");
  }

  // send local code name
  write_name(local_code_name);

  // send param name
  write_name(param_name);

  // send initial value
  CWP_Type_t endian_data_type = data_type;
  CWP_swap_endian_4bytes((int *) &endian_data_type, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_data_type, sizeof(CWP_Type_t));

  switch (data_type) {

  case CWP_DOUBLE: {
    double endian_double_value = * ((double *) value);
    CWP_swap_endian_8bytes(&endian_double_value, 1);
    CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_double_value, sizeof(double));
    } break;

  case CWP_INT: {
    int endian_int_value = * ((int *) value);
    CWP_swap_endian_4bytes(&endian_int_value, 1);
    CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_int_value, sizeof(int));
    } break;

  case CWP_CHAR: {
    write_name((char *) value);
    } break;

  default:
    PDM_error(__FILE__, __LINE__, 0, "Unknown CWP_Type_t %i\n", data_type);
  }
}

void
CWP_client_Param_del
(
 const char       *local_code_name,
 const char       *param_name,
 const CWP_Type_t  data_type
)
{
  t_message msg;

  // verbose
  if (clt->flags & CWP_CLIENTFLAG_VERBOSE) {
    log_trace("CWP:Client initiating CWP_Param_del\n");
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_PARAM_DEL);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Param_del failed to send message header\n");
  }

  // send local code name
  write_name(local_code_name);

  // send param name
  write_name(param_name);

  // send data type
  CWP_Type_t endian_data_type = data_type;
  CWP_swap_endian_4bytes((int *) &endian_data_type, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_data_type, sizeof(CWP_Type_t));
}

int
CWP_client_Param_n_get
(
 const char             *code_name,
 const CWP_Type_t        data_type
)
{
  t_message msg;

  // verbose
  if (clt->flags & CWP_CLIENTFLAG_VERBOSE) {
    log_trace("CWP:Client initiating CWP_Param_n_get\n");
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_PARAM_N_GET);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Param_n_get failed to send message header\n");
  }

  // send local code name
  write_name(code_name);

  // send data type
  CWP_Type_t endian_data_type = data_type;
  CWP_swap_endian_4bytes((int *) &endian_data_type, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size, (void*) &endian_data_type, sizeof(CWP_Type_t));

  // read n_param
  int n_param = -1;
  CWP_transfer_readdata(clt->socket, clt->max_msg_size, (void*) &n_param, sizeof(int));

  return n_param;
}

void
CWP_client_Param_list_get
(
 const char             *code_name,
 const CWP_Type_t        data_type,
 int                    *nParam,
 char                 ***paramNames
)
{
  t_message msg;

  // verbose
  if (clt->flags & CWP_CLIENTFLAG_VERBOSE) {
    log_trace("CWP:Client initiating CWP_Param_list_get\n");
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_PARAM_LIST_GET);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Param_list_get failed to send message header\n");
  }

  // send local code name
  write_name(code_name);

  // send data type
  CWP_Type_t endian_data_type = data_type;
  CWP_swap_endian_4bytes((int *) &endian_data_type, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size, (void*) &endian_data_type, sizeof(CWP_Type_t));

  // read n_param
  CWP_transfer_readdata(clt->socket, clt->max_msg_size, (void*) nParam, sizeof(int));

  // read param names
  *paramNames = (char **) malloc(sizeof(char *) * (*nParam));
  for (int i = 0; i < *nParam; i++) {
    (*paramNames)[i] = (char *) malloc(sizeof(char));
    read_name(&(*paramNames)[i]);
  }
}

int
CWP_client_Param_is
(
 const char             *code_name,
 const char             *param_name,
 const CWP_Type_t        data_type
)
{
  t_message msg;

  // verbose
  if (clt->flags & CWP_CLIENTFLAG_VERBOSE) {
    log_trace("CWP:Client initiating CWP_Param_is\n");
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_PARAM_IS);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Param_is failed to send message header\n");
  }

  // send local code name
  write_name(code_name);

  // send param name
  write_name(param_name);

  // send data type
  CWP_Type_t endian_data_type = data_type;
  CWP_swap_endian_4bytes((int *) &endian_data_type, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_data_type, sizeof(CWP_Type_t));

  // read bool
  int is_param = -1;
  CWP_transfer_readdata(clt->socket, clt->max_msg_size, (void*) &is_param, sizeof(int));

  return is_param;
}

void
CWP_client_Param_reduce
(
 const CWP_Op_t    op,
 const char       *param_name,
 const CWP_Type_t  data_type,
 void             *res,
 const int         nCode,
 ...
)
{
  t_message msg;

  // verbose
  if (clt->flags & CWP_CLIENTFLAG_VERBOSE) {
    log_trace("CWP:Client initiating CWP_Param_reduce\n");
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_PARAM_REDUCE);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Param_reduce failed to send message header\n");
  }

  // send operation type
  CWP_Op_t endian_op = op;
  CWP_swap_endian_4bytes((int *) &endian_op, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size, (void*) &endian_op, sizeof(CWP_Op_t));

  // send prameter name
  write_name(param_name);

  // send data type
  CWP_Type_t endian_data_type = data_type;
  CWP_swap_endian_4bytes((int *) &endian_data_type, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size, (void*) &endian_data_type, sizeof(CWP_Type_t));

  // send number of codes
  int endian_nCode = nCode;
  CWP_swap_endian_4bytes(&endian_nCode, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size, (void*) &endian_nCode, sizeof(int));

  // other variables TO DO: how to deal with endianness of va_list parameters ?
  va_list pa;
  va_start(pa, nCode);
  log_trace("CWP: CWP_Param_reduce va_list parameters sending not implemented yet\n");
  // CWP_transfer_writedata(clt->socket,clt->max_msg_size, (void*) &pa, sizeof(pa));

  // retreive result
  switch (data_type) {

  case CWP_DOUBLE: {
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, (void*) res, sizeof(double));
    } break;

  case CWP_INT: {
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, (void*) res, sizeof(int));
    } break;

  case CWP_CHAR: {
    res = (char *) malloc(sizeof(char));
    read_name((char **) &res);
    } break;

  default:
    PDM_error(__FILE__, __LINE__, 0, "Unknown CWP_Type_t %i\n", data_type);
  }
}

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
)
{
  t_message msg;

  // verbose
  if (clt->flags & CWP_CLIENTFLAG_VERBOSE) {
    log_trace("CWP:Client initiating CWP_Cpl_create\n");
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_CPL_CREATE);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Cpl_create failed to send message header\n");
  }

  // send code name
  write_name(local_code_name);

  // send coupling identifier
  write_name(cpl_id);

  // send coupled code name
  write_name(coupled_code_name);

  // send entities dimension
  CWP_Interface_t endian_entities_dim = entities_dim;
  CWP_swap_endian_4bytes((int *) &endian_entities_dim, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_entities_dim, sizeof(CWP_Interface_t));

  // send communication type
  CWP_Comm_t endian_comm_type = comm_type;
  CWP_swap_endian_4bytes((int *) &endian_comm_type, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_comm_type, sizeof(CWP_Comm_t));

  // send spatial interpolation
  CWP_Spatial_interp_t endian_spatial_interp = spatial_interp;
  CWP_swap_endian_4bytes((int *) &endian_spatial_interp, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_spatial_interp, sizeof(CWP_Spatial_interp_t));

  // send number of partitions
  int endian_n_part = n_part;
  CWP_swap_endian_4bytes(&endian_n_part, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_n_part, sizeof(int));

  // send displacement type
  CWP_Dynamic_mesh_t endian_displacement = displacement;
  CWP_swap_endian_4bytes((int *) &endian_displacement, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_displacement, sizeof(CWP_Dynamic_mesh_t));

  // send time exchange type
  CWP_Time_exch_t endian_recv_freq_type = recv_freq_type;
  CWP_swap_endian_4bytes((int *) &endian_recv_freq_type, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_recv_freq_type, sizeof(CWP_Time_exch_t));
}

void
CWP_client_Cpl_del
(
 const char *local_code_name,
 const char *cpl_id
)
{
  t_message msg;

  // verbose
  if (clt->flags & CWP_CLIENTFLAG_VERBOSE) {
    log_trace("CWP:Client initiating CWP_Cpl_del\n");
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_CPL_DEL);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Cpl_del failed to send message header\n");
  }

  // send code name
  write_name(local_code_name);

  // send coupling identifier
  write_name(cpl_id);
}

void
CWP_client_Properties_dump
(
void
)
{
  t_message msg;

  // verbose
  if (clt->flags & CWP_CLIENTFLAG_VERBOSE) {
    log_trace("CWP:Client initiating CWP_Properties_dump\n");
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_PROPERTIES_DUMP);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Properties_dump failed to send message header\n");
  }
}

void
CWP_client_Visu_set
(
 const char                 *local_code_name,
 const char                 *cpl_id,
 const int                   freq,
 const CWP_Visu_format_t     format,
 const char                 *format_option
)
{
  t_message msg;

  // verbose
  if (clt->flags & CWP_CLIENTFLAG_VERBOSE) {
    log_trace("CWP:Client initiating CWP_Visu_set\n");
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_VISU_SET);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Visu_set failed to send message header\n");
  }

  // send code name
  write_name(local_code_name);

  // send coupling identifier
  write_name(cpl_id);

  // send frequency
  int endian_freq = freq;
  CWP_swap_endian_4bytes(&endian_freq, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_freq, sizeof(int));

  // send format
  CWP_Visu_format_t endian_format = format;
  CWP_swap_endian_4bytes((int *) &endian_format, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_format, sizeof(CWP_Visu_format_t));

  // send format option
  write_name(format_option);
}

void
CWP_client_State_update
(
 const char* local_code_name,
 const CWP_State_t state
)
{
  t_message msg;

  // verbose
  if (clt->flags & CWP_CLIENTFLAG_VERBOSE) {
    log_trace("CWP:Client initiating CWP_State_update\n");
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_STATE_UPDATE);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_State_update failed to send message header\n");
  }

  // send code name
  write_name(local_code_name);

  // send state
  CWP_State_t endian_state = state;
  CWP_swap_endian_4bytes((int *) &endian_state, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_state, sizeof(CWP_State_t));
}

void
CWP_client_Time_update
(
 const char* local_code_name,
 const double current_time
)
{
  t_message msg;

  // verbose
  if (clt->flags & CWP_CLIENTFLAG_VERBOSE) {
    log_trace("CWP:Client initiating CWP_Time_update\n");
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_TIME_UPDATE);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Time_update failed to send message header\n");
  }

  // send code name
  write_name(local_code_name);

  // send time
  double endian_current_time = current_time;
  CWP_swap_endian_8bytes(&endian_current_time, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_current_time, sizeof(double));
}

void
CWP_client_Output_file_set
(
const char *output_filename
)
{
  t_message msg;

  // verbose
  if (clt->flags & CWP_CLIENTFLAG_VERBOSE) {
    log_trace("CWP:Client initiating CWP_Output_file_set\n");
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_OUTPUT_FILE_SET);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Output_file_set failed to send message header\n");
  }

  // send output file name or directory
  write_name(output_filename);
}

void
CWP_client_User_structure_set
(
 const char* local_code_name,
       void* user_structure
)
{
  PDM_UNUSED(local_code_name);
  PDM_UNUSED(user_structure);
  log_trace("CWP: CWP_User_structure_set not implemented in client/server mode\n");
}

void *
CWP_client_User_structure_get
(
 const char* local_code_name
)
{
  PDM_UNUSED(local_code_name);
  log_trace("CWP: CWP_User_structure_get not implemented in client/server mode\n");

  return 0;
}

CWP_State_t
CWP_client_State_get
(
 const char    *code_name
)
{
  t_message msg;

  // verbose
  if (clt->flags & CWP_CLIENTFLAG_VERBOSE) {
    log_trace("CWP:Client initiating CWP_State_get\n");
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_STATE_GET);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_State_get failed to send message header\n");
  }

  // send code name
  write_name(code_name);

  // read state
  CWP_State_t state = (CWP_State_t) -1;
  CWP_transfer_readdata(clt->socket, clt->max_msg_size, (void*) &state, sizeof(CWP_State_t));

  return state;
}

int
CWP_client_Codes_nb_get
(
 void
)
{
  t_message msg;

  // verbose
  if (clt->flags & CWP_CLIENTFLAG_VERBOSE) {
    log_trace("CWP:Client initiating CWP_Codes_nb_get\n");
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_CODES_NB_GET);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Codes_nb_get failed to send message header\n");
  }

  // read nb_codes
  clt_cwp->n_code_names = -1;
  CWP_transfer_readdata(clt->socket, clt->max_msg_size, (void*) &clt_cwp->n_code_names, sizeof(int));

  return clt_cwp->n_code_names;
}

const char **
CWP_client_Codes_list_get
(
void
)
{
  t_message msg;

  // verbose
  if (clt->flags & CWP_CLIENTFLAG_VERBOSE) {
    log_trace("CWP:Client initiating CWP_Codes_list_get\n");
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_CODES_LIST_GET);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Codes_list_get failed to send message header\n");
  }

  // read code names
  clt_cwp->n_code_names = -1;
  CWP_transfer_readdata(clt->socket, clt->max_msg_size, (void*) &clt_cwp->n_code_names, sizeof(int));
  clt_cwp->code_names = (const char **) malloc(sizeof(char *) * clt_cwp->n_code_names);
  for (int i = 0; i < clt_cwp->n_code_names; i++) {
    clt_cwp->code_names[i] = (const char *) malloc(sizeof(char));
    read_name((char **) &clt_cwp->code_names[i]);
  }

  return clt_cwp->code_names;
}

int
CWP_client_Loc_codes_nb_get
(
 void
)
{
  t_message msg;

  // verbose
  if (clt->flags & CWP_CLIENTFLAG_VERBOSE) {
    log_trace("CWP:Client initiating CWP_Loc_codes_nb_get\n");
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_LOC_CODES_NB_GET);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Loc_codes_nb_get failed to send message header\n");
  }

  // read nb_codes
  clt_cwp->n_loc_code_names = -1;
  CWP_transfer_readdata(clt->socket, clt->max_msg_size, (void*) &clt_cwp->n_loc_code_names, sizeof(int));

  return clt_cwp->n_loc_code_names;
}

const char **
CWP_client_Loc_codes_list_get
(
 void
)
{
  t_message msg;

  // verbose
  if (clt->flags & CWP_CLIENTFLAG_VERBOSE) {
    log_trace("CWP:Client initiating CWP_Loc_codes_list_get\n");
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_LOC_CODES_LIST_GET);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Loc_codes_list_get failed to send message header\n");
  }

  // read code names
  clt_cwp->n_loc_code_names = -1;
  CWP_transfer_readdata(clt->socket, clt->max_msg_size, (void*) &clt_cwp->n_loc_code_names, sizeof(int));
  clt_cwp->loc_code_names = (const char **) malloc(sizeof(char *) * clt_cwp->n_loc_code_names);
  for (int i = 0; i < clt_cwp->n_loc_code_names; i++) {
    clt_cwp->loc_code_names[i] = (const char *) malloc(sizeof(char));
    read_name((char **) &clt_cwp->loc_code_names[i]);
  }

  return clt_cwp->loc_code_names;
}

int
CWP_client_N_uncomputed_tgts_get
(
 const char *local_code_name,
 const char *cpl_id,
 const char *field_id,
 const int   i_part
)
{
  t_message msg;

  // verbose
  if (clt->flags & CWP_CLIENTFLAG_VERBOSE) {
    log_trace("CWP:Client initiating CWP_N_uncomputed_tgts_get\n");
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_N_UNCOMPUTED_TGTS_GET);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_N_uncomputed_tgts_get failed to send message header\n");
  }

  // send local code name
  write_name(local_code_name);

  // send coupling identifier
  write_name(cpl_id);

  // send field identifier
  write_name(field_id);

  // send i_part
  int endian_i_part = i_part;
  CWP_swap_endian_4bytes(&endian_i_part, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_i_part, sizeof(int));

  // read number of targets
  int nb_tgts = -1;
  CWP_transfer_readdata(clt->socket, clt->max_msg_size, (void*) &nb_tgts, sizeof(int));

  return nb_tgts;
}

const int *
CWP_client_Uncomputed_tgts_get
(
 const char *local_code_name,
 const char *cpl_id,
 const char *field_id,
 const int   i_part
)
{
  t_message msg;

  // verbose
  if (clt->flags & CWP_CLIENTFLAG_VERBOSE) {
    log_trace("CWP:Client initiating CWP_Uncomputed_tgts_get\n");
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_UNCOMPUTED_TGTS_GET);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Uncomputed_tgts_get failed to send message header\n");
  }

  // send local code name
  write_name(local_code_name);

  // send coupling identifier
  write_name(cpl_id);

  // send field identifier
  write_name(field_id);

  // send i_part
  int endian_i_part = i_part;
  CWP_swap_endian_4bytes(&endian_i_part, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_i_part, sizeof(int));

  // read number of targets
  int nb_tgts = -1;
  CWP_transfer_readdata(clt->socket, clt->max_msg_size, (void*) &nb_tgts, sizeof(int));
  int *tgts = (int *) malloc(sizeof(int) * nb_tgts);
  CWP_transfer_readdata(clt->socket, clt->max_msg_size, (void*) tgts, sizeof(int)*nb_tgts);

  return tgts;
}

int
CWP_client_N_computed_tgts_get
(
 const char *local_code_name,
 const char *cpl_id,
 const char *field_id,
 const int   i_part
)
{
  t_message msg;

  // verbose
  if (clt->flags & CWP_CLIENTFLAG_VERBOSE) {
    log_trace("CWP:Client initiating CWP_N_computed_tgts_get\n");
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_N_COMPUTED_TGTS_GET);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_N_computed_tgts_get failed to send message header\n");
  }

  // send local code name
  write_name(local_code_name);

  // send coupling identifier
  write_name(cpl_id);

  // send field identifier
  write_name(field_id);

  // send i_part
  int endian_i_part = i_part;
  CWP_swap_endian_4bytes(&endian_i_part, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_i_part, sizeof(int));

  // read number of targets
  int nb_tgts = -1;
  CWP_transfer_readdata(clt->socket, clt->max_msg_size, (void*) &nb_tgts, sizeof(int));

  return nb_tgts;
}

const int *
CWP_client_Computed_tgts_get
(
 const char *local_code_name,
 const char *cpl_id,
 const char *field_id,
 const int   i_part
)
{
  t_message msg;

  // verbose
  if (clt->flags & CWP_CLIENTFLAG_VERBOSE) {
    log_trace("CWP:Client initiating CWP_Computed_tgts_get\n");
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_COMPUTED_TGTS_GET);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Computed_tgts_get failed to send message header\n");
  }

  // send local code name
  write_name(local_code_name);

  // send coupling identifier
  write_name(cpl_id);

  // send field identifier
  write_name(field_id);

  // send i_part
  int endian_i_part = i_part;
  CWP_swap_endian_4bytes(&endian_i_part, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_i_part, sizeof(int));

  // read number of targets
  int nb_tgts = -1;
  CWP_transfer_readdata(clt->socket, clt->max_msg_size, (void*) &nb_tgts, sizeof(int));
  int *tgts = (int *) malloc(sizeof(int) * nb_tgts);
  CWP_transfer_readdata(clt->socket, clt->max_msg_size, (void*) tgts, sizeof(int)*nb_tgts);

  return tgts;
}

int
CWP_client_N_involved_srcs_get
(
 const char *local_code_name,
 const char *cpl_id,
 const char *field_id,
 const int   i_part
)
{
  t_message msg;

  // verbose
  if (clt->flags & CWP_CLIENTFLAG_VERBOSE) {
    log_trace("CWP:Client initiating CWP_N_involved_srcs_get\n");
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_N_INVOLVED_SRCS_GET);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_N_involved_srcs_get failed to send message header\n");
  }

  // send local code name
  write_name(local_code_name);

  // send coupling identifier
  write_name(cpl_id);

  // send field identifier
  write_name(field_id);

  // send i_part
  int endian_i_part = i_part;
  CWP_swap_endian_4bytes(&endian_i_part, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_i_part, sizeof(int));

  // read number of sources
  int nb_srcs = -1;
  CWP_transfer_readdata(clt->socket, clt->max_msg_size, (void*) &nb_srcs, sizeof(int));

  return nb_srcs;
}

const int *
CWP_client_Involved_srcs_get
(
 const char *local_code_name,
 const char *cpl_id,
 const char *field_id,
 const int   i_part
)
{
  t_message msg;

  // verbose
  if (clt->flags & CWP_CLIENTFLAG_VERBOSE) {
    log_trace("CWP:Client initiating CWP_Involved_srcs_get\n");
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_INVOLVED_SRCS_GET);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Involved_srcs_get failed to send message header\n");
  }

  // send local code name
  write_name(local_code_name);

  // send coupling identifier
  write_name(cpl_id);

  // send field identifier
  write_name(field_id);

  // send i_part
  int endian_i_part = i_part;
  CWP_swap_endian_4bytes(&endian_i_part, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_i_part, sizeof(int));

  // read number of sources
  int nb_srcs = -1;
  CWP_transfer_readdata(clt->socket, clt->max_msg_size, (void*) &nb_srcs, sizeof(int));
  int *srcs = (int *) malloc(sizeof(int) * nb_srcs);
  CWP_transfer_readdata(clt->socket, clt->max_msg_size, (void*) srcs, sizeof(int)*nb_srcs);

  return srcs;
}

void
CWP_client_Spatial_interp_weights_compute
(
 const char     *local_code_name,
 const char     *cpl_id
)
{
  t_message msg;

  // verbose
  if (clt->flags & CWP_CLIENTFLAG_VERBOSE) {
    log_trace("CWP:Client initiating CWP_Spatial_interp_weights_compute\n");
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_SPATIAL_INTERP_WEIGHTS_COMPUTE);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Spatial_interp_weights_compute failed to send message header\n");
  }

  // send local code name
  write_name(local_code_name);

  // send coupling identifier
  write_name(cpl_id);
}

void
CWP_client_Spatial_interp_property_set
(
 const char     *local_code_name,
 const char     *cpl_id,
 const char     *property_name,
 const char     *property_type,
 const char     *property_value
)
{
  t_message msg;

  // verbose
  if (clt->flags & CWP_CLIENTFLAG_VERBOSE) {
    log_trace("CWP:Client initiating CWP_Spatial_interp_property_set\n");
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_SPATIAL_INTERP_PROPERTY_SET);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Spatial_interp_property_set failed to send message header\n");
  }

  // send local code name
  write_name(local_code_name);

  // send coupling identifier
  write_name(cpl_id);

  // send property name
  write_name(property_name);

  // send property type
  write_name(property_type);

  // send property value
  write_name(property_value);
}

void
CWP_client_User_tgt_pts_set
(
 const char    *local_code_name,
 const char    *cpl_id,
 const int      i_part,
 const int      n_pts,
 double         coord[],
 CWP_g_num_t    global_num[]
)
{
  t_message msg;

  // verbose
  if (clt->flags & CWP_CLIENTFLAG_VERBOSE) {
    log_trace("CWP:Client initiating CWP_User_tgt_pts_set\n");
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_USER_TGT_PTS_SET);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_User_tgt_pts_set failed to send message header\n");
  }

  // send local code name
  write_name(local_code_name);

  // send coupling identifier
  write_name(cpl_id);

  // send i_part
  int endian_i_part = i_part;
  CWP_swap_endian_4bytes(&endian_i_part, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_i_part, sizeof(int));

  // send n_pts
  int endian_n_pts = n_pts;
  CWP_swap_endian_4bytes(&endian_n_pts, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_n_pts, sizeof(int));

  // send coord
  double *endian_coord = (double *) malloc(sizeof(double) * 3 * n_pts);
  memcpy(endian_coord, coord, sizeof(double) * 3 * n_pts);
  CWP_swap_endian_8bytes(endian_coord, 3 * n_pts);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) endian_coord, sizeof(double) * 3 * n_pts);

  // send global_num
  CWP_g_num_t *endian_global_num = (CWP_g_num_t *) malloc(sizeof(CWP_g_num_t) * n_pts);
  memcpy(endian_global_num, global_num, sizeof(CWP_g_num_t) * n_pts);
  if (sizeof(CWP_g_num_t) == 4) {
    CWP_swap_endian_4bytes((int *) endian_global_num, n_pts);
  }
  else  if (sizeof(CWP_g_num_t) == 8) {
    CWP_swap_endian_8bytes((double *) endian_global_num, n_pts);
  }
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) endian_global_num, sizeof(CWP_g_num_t) * n_pts);

  // free
  free(endian_coord);
  free(endian_global_num);
}

void
CWP_client_Mesh_interf_finalize
(
 const char           *local_code_name,
 const char           *cpl_id
)
{
  t_message msg;

  // verbose
  if (clt->flags & CWP_CLIENTFLAG_VERBOSE) {
    log_trace("CWP:Client initiating CWP_Mesh_interf_finalize\n");
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_MESH_INTERF_FINALIZE);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Mesh_interf_finalize failed to send message header\n");
  }

  // send local code name
  write_name(local_code_name);

  // send coupling identifier
  write_name(cpl_id);
}

void
CWP_client_Mesh_interf_vtx_set
(
 const char           *local_code_name,
 const char           *cpl_id,
 const int             i_part,
 const int             n_pts,
 double                coord[],
 CWP_g_num_t           global_num[]
)
{
  t_message msg;

  // verbose
  if (clt->flags & CWP_CLIENTFLAG_VERBOSE) {
    log_trace("CWP:Client initiating CWP_Mesh_interf_vtx_set\n");
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_MESH_INTERF_VTX_SET);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Mesh_interf_vtx_set failed to send message header\n");
  }

  // send local code name
  write_name(local_code_name);

  // send coupling identifier
  write_name(cpl_id);

  // send i_part
  int endian_i_part = i_part;
  CWP_swap_endian_4bytes(&endian_i_part, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_i_part, sizeof(int));

  // send n_pts
  int endian_n_pts = n_pts;
  CWP_swap_endian_4bytes(&endian_n_pts, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_n_pts, sizeof(int));

  // send coord
  double *endian_coord = (double *) malloc(sizeof(double) * 3 * n_pts);
  memcpy(endian_coord, coord, sizeof(double) * 3 * n_pts);
  CWP_swap_endian_8bytes(endian_coord, 3 * n_pts);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) endian_coord, sizeof(double) * 3 * n_pts);

  // send global_num
  CWP_g_num_t *endian_global_num = (CWP_g_num_t *) malloc(sizeof(CWP_g_num_t) * n_pts);
  memcpy(endian_global_num, global_num, sizeof(CWP_g_num_t) * n_pts);
  if (sizeof(CWP_g_num_t) == 4) {
    CWP_swap_endian_4bytes((int *) endian_global_num, n_pts);
  }
  else  if (sizeof(CWP_g_num_t) == 8) {
    CWP_swap_endian_8bytes((double *) endian_global_num, n_pts);
  }
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) endian_global_num, sizeof(CWP_g_num_t) * n_pts);

  // free
  free(endian_coord);
  free(endian_global_num);
}

int
CWP_client_Mesh_interf_block_add
(
 const char           *local_code_name,
 const char           *cpl_id,
 const CWP_Block_t     block_type
)
{
  t_message msg;

  // verbose
  if (clt->flags & CWP_CLIENTFLAG_VERBOSE) {
    log_trace("CWP:Client initiating CWP_Mesh_interf_block_add\n");
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_MESH_INTERF_BLOCK_ADD);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Mesh_interf_block_add failed to send message header\n");
  }

  // send local code name
  write_name(local_code_name);

  // send coupling identifier
  write_name(cpl_id);

  // send block type
  CWP_Block_t endian_block_type = block_type;
  CWP_swap_endian_4bytes((int *) &endian_block_type, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_block_type, sizeof(CWP_Block_t));

  // read block identifier
  int block_id = -2;
  CWP_transfer_readdata(clt->socket, clt->max_msg_size, (void*) &block_id, sizeof(int));

  return block_id;
}

void
CWP_client_Mesh_interf_block_std_set
(
 const char        *local_code_name,
 const char        *cpl_id,
 const int          i_part,
 const int          block_id,
 const CWP_Block_t  block_type,
 const int          n_elts,
 int                connec[],
 CWP_g_num_t        global_num[]
)
{
  t_message msg;

  // verbose
  if (clt->flags & CWP_CLIENTFLAG_VERBOSE) {
    log_trace("CWP:Client initiating CWP_Mesh_interf_block_std_set\n");
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_MESH_INTERF_BLOCK_STD_SET);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Mesh_interf_block_std_set failed to send message header\n");
  }

  // send local code name
  write_name(local_code_name);

  // send coupling identifier
  write_name(cpl_id);

  // send i_part
  int endian_i_part = i_part;
  CWP_swap_endian_4bytes(&endian_i_part, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_i_part, sizeof(int));

  // send block_id
  int endian_block_id = block_id;
  CWP_swap_endian_4bytes(&endian_block_id, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_block_id, sizeof(int));

  // send n_elts
  int endian_n_elts = n_elts;
  CWP_swap_endian_4bytes(&endian_n_elts, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_n_elts, sizeof(int));

  // send n_vtx_elt
  // TO DO: knowing the size of connec one could retreive n_vtx_elt also for POLY
  int n_vtx_elt = n_nodes_get(block_type);
  int endian_n_vtx_elt = n_vtx_elt;
  CWP_swap_endian_4bytes(&endian_n_vtx_elt, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_n_vtx_elt, sizeof(int));

  // send connectivity
  int *endian_connec = (int *) malloc(sizeof(int) * n_elts * n_vtx_elt);
  memcpy(endian_connec, connec, sizeof(int) * n_elts * n_vtx_elt);
  CWP_swap_endian_4bytes(endian_connec, n_elts * n_vtx_elt);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) endian_connec, sizeof(int) * n_elts * n_vtx_elt);

  // send global number
  CWP_g_num_t  *endian_global_num = (CWP_g_num_t  *) malloc(sizeof(CWP_g_num_t) * n_elts);
  memcpy(endian_global_num, global_num, sizeof(CWP_g_num_t) * n_elts);
  if (sizeof(CWP_g_num_t) == 4) {
    CWP_swap_endian_4bytes((int *) endian_global_num, n_elts);
  }
  else  if (sizeof(CWP_g_num_t) == 8) {
    CWP_swap_endian_8bytes((double *) endian_global_num, n_elts);
  }
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) endian_global_num, sizeof(CWP_g_num_t) * n_elts);

  // free
  free(endian_connec);
  free(endian_global_num);
}

void
CWP_client_Mesh_interf_f_poly_block_set
(
 const char             *local_code_name,
 const char             *cpl_id,
 const int               i_part,
 const int               block_id,
 const int               n_elts,
 int                     connec_idx[],
 int                     connec[],
 CWP_g_num_t             global_num[]
)
{
  t_message msg;

  // verbose
  if (clt->flags & CWP_CLIENTFLAG_VERBOSE) {
    log_trace("CWP:Client initiating CWP_Mesh_interf_f_poly_block_set\n");
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_MESH_INTERF_F_POLY_BLOCK_SET);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Mesh_interf_f_poly_block_set failed to send message header\n");
  }

  // send local code name
  write_name(local_code_name);

  // send coupling identifier
  write_name(cpl_id);

  // send i_part
  int endian_i_part = i_part;
  CWP_swap_endian_4bytes(&endian_i_part, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_i_part, sizeof(int));

  // send block_id
  int endian_block_id = block_id;
  CWP_swap_endian_4bytes(&endian_block_id, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_block_id, sizeof(int));

  // send n_elts
  int endian_n_elts = n_elts;
  CWP_swap_endian_4bytes(&endian_n_elts, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_n_elts, sizeof(int));

  // send connectivity index
  int *endian_connec_idx = (int *) malloc(sizeof(int) * (n_elts+1));
  memcpy(endian_connec_idx, connec_idx, sizeof(int) * (n_elts+1));
  CWP_swap_endian_4bytes(endian_connec_idx, n_elts + 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) endian_connec_idx, sizeof(int) * (n_elts+1));

  // send connectivity
  int *endian_connec = (int *) malloc(sizeof(int) * connec_idx[n_elts]);
  memcpy(endian_connec, connec, sizeof(int) * connec_idx[n_elts]);
  CWP_swap_endian_4bytes(endian_connec, connec_idx[n_elts]);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) endian_connec, sizeof(int) * connec_idx[n_elts]);

  // send global number
  CWP_g_num_t  *endian_global_num = (CWP_g_num_t  *) malloc(sizeof(CWP_g_num_t) * n_elts);
  memcpy(endian_global_num, global_num, sizeof(CWP_g_num_t) * n_elts);
  if (sizeof(CWP_g_num_t) == 4) {
    CWP_swap_endian_4bytes((int *) endian_global_num, n_elts);
  }
  else  if (sizeof(CWP_g_num_t) == 8) {
    CWP_swap_endian_8bytes((double *) endian_global_num, n_elts);
  }
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) endian_global_num, sizeof(CWP_g_num_t) * n_elts);

  // free
  free(endian_connec_idx);
  free(endian_connec);
  free(endian_global_num);
}

void
CWP_client_Mesh_interf_f_poly_block_get
(
 const char             *local_code_name,
 const char             *cpl_id,
 const int               i_part,
 const int               block_id,
 int                    *n_elts,
 int                   **connec_idx,
 int                   **connec,
 CWP_g_num_t           **global_num
)
{
  t_message msg;

  // verbose
  if (clt->flags & CWP_CLIENTFLAG_VERBOSE) {
    log_trace("CWP:Client initiating CWP_Mesh_interf_f_poly_block_get\n");
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_MESH_INTERF_F_POLY_BLOCK_GET);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Mesh_interf_f_poly_block_get failed to send message header\n");
  }

  // send local code name
  write_name(local_code_name);

  // send coupling identifier
  write_name(cpl_id);

  // send i_part
  int endian_i_part = i_part;
  CWP_swap_endian_4bytes(&endian_i_part, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_i_part, sizeof(int));

  // send block_id
  int endian_block_id = block_id;
  CWP_swap_endian_4bytes(&endian_block_id, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_block_id, sizeof(int));

  // read n_elts
  CWP_transfer_readdata(clt->socket, clt->max_msg_size, n_elts, sizeof(int));

  // read connectivity index
  CWP_transfer_readdata(clt->socket, clt->max_msg_size, *connec_idx, sizeof(int) * ((*n_elts)+1));

  // read connectivity
  CWP_transfer_readdata(clt->socket, clt->max_msg_size, *connec, sizeof(int) * (*connec_idx)[(*n_elts)]);

  // read global number
  CWP_transfer_readdata(clt->socket, clt->max_msg_size, *global_num, sizeof(CWP_g_num_t) * (*n_elts));
}

void
CWP_client_Mesh_interf_c_poly_block_set
(
 const char           *local_code_name,
 const char           *cpl_id,
 const int             i_part,
 const int             block_id,
 const int             n_elts,
 const int             n_faces,
 int                   connec_faces_idx[],
 int                   connec_faces[],
 int                   connec_cells_idx[],
 int                   connec_cells[],
 CWP_g_num_t           global_num[]
)
{
  t_message msg;

  // verbose
  if (clt->flags & CWP_CLIENTFLAG_VERBOSE) {
    log_trace("CWP:Client initiating CWP_Mesh_interf_c_poly_block_set\n");
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_MESH_INTERF_C_POLY_BLOCK_SET);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Mesh_interf_c_poly_block_set failed to send message header\n");
  }

  // send local code name
  write_name(local_code_name);

  // send coupling identifier
  write_name(cpl_id);

  // send i_part
  int endian_i_part = i_part;
  CWP_swap_endian_4bytes(&endian_i_part, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_i_part, sizeof(int));

  // send block_id
  int endian_block_id = block_id;
  CWP_swap_endian_4bytes(&endian_block_id, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_block_id, sizeof(int));

  // send n_elts
  int endian_n_elts = n_elts;
  CWP_swap_endian_4bytes(&endian_n_elts, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_n_elts, sizeof(int));

  // send n_faces
  int endian_n_faces = n_faces;
  CWP_swap_endian_4bytes(&endian_n_faces, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_n_faces, sizeof(int));

  // send connectivity of faces index
  int *endian_connec_faces_idx = (int *) malloc(sizeof(int) * (n_faces+1));
  memcpy(endian_connec_faces_idx, connec_faces_idx, sizeof(int) * (n_faces+1));
  CWP_swap_endian_4bytes(endian_connec_faces_idx, n_faces + 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) endian_connec_faces_idx, sizeof(int) * (n_faces+1));

  // send connec_facestivity faces
  int *endian_connec_faces = (int *) malloc(sizeof(int) * connec_faces_idx[n_faces]);
  memcpy(endian_connec_faces, connec_faces, sizeof(int) * connec_faces_idx[n_faces]);
  CWP_swap_endian_4bytes(endian_connec_faces, connec_faces_idx[n_faces]);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) endian_connec_faces, sizeof(int) * connec_faces_idx[n_faces]);

  // send connectivity of cells index
  int *endian_connec_cells_idx = (int *) malloc(sizeof(int) * (n_elts+1));
  memcpy(endian_connec_cells_idx, connec_cells_idx, sizeof(int) * (n_elts+1));
  CWP_swap_endian_4bytes(endian_connec_cells_idx, n_elts + 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) endian_connec_cells_idx, sizeof(int) * (n_elts+1));

  // send connec_cellstivity cells
  int *endian_connec_cells = (int *) malloc(sizeof(int) * connec_cells_idx[n_elts]);
  memcpy(endian_connec_cells, connec_cells, sizeof(int) * connec_cells_idx[n_elts]);
  CWP_swap_endian_4bytes(endian_connec_cells, connec_cells_idx[n_elts]);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) endian_connec_cells, sizeof(int) * connec_cells_idx[n_elts]);

  // send global number
  CWP_g_num_t  *endian_global_num = (CWP_g_num_t  *) malloc(sizeof(CWP_g_num_t) * n_elts);
  memcpy(endian_global_num, global_num, sizeof(CWP_g_num_t) * n_elts);
  if (sizeof(CWP_g_num_t) == 4) {
    CWP_swap_endian_4bytes((int *) endian_global_num, n_elts);
  }
  else  if (sizeof(CWP_g_num_t) == 8) {
    CWP_swap_endian_8bytes((double *) endian_global_num, n_elts);
  }
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) endian_global_num, sizeof(CWP_g_num_t) * n_elts);

  // free
  free(endian_connec_faces_idx);
  free(endian_connec_faces);
  free(endian_connec_cells_idx);
  free(endian_connec_cells);
  free(endian_global_num);
}

void
CWP_client_Mesh_interf_c_poly_block_get
(
 const char           *local_code_name,
 const char           *cpl_id,
 const int             i_part,
 const int             block_id,
 int                  *n_elts,
 int                  *n_faces,
 int                 **connec_faces_idx,
 int                 **connec_faces,
 int                 **connec_cells_idx,
 int                 **connec_cells,
 CWP_g_num_t         **global_num
)
{
  t_message msg;

  // verbose
  if (clt->flags & CWP_CLIENTFLAG_VERBOSE) {
    log_trace("CWP:Client initiating CWP_Mesh_interf_c_poly_block_get\n");
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_MESH_INTERF_C_POLY_BLOCK_GET);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Mesh_interf_c_poly_block_get failed to send message header\n");
  }

  // send local code name
  write_name(local_code_name);

  // send coupling identifier
  write_name(cpl_id);

  // send i_part
  int endian_i_part = i_part;
  CWP_swap_endian_4bytes(&endian_i_part, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_i_part, sizeof(int));

  // send block_id
  int endian_block_id = block_id;
  CWP_swap_endian_4bytes(&endian_block_id, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_block_id, sizeof(int));

  // read n_elts
  CWP_transfer_readdata(clt->socket, clt->max_msg_size, n_elts, sizeof(int));

  // read n_faces
  CWP_transfer_readdata(clt->socket, clt->max_msg_size, n_faces, sizeof(int));

  // read connectivity index
  CWP_transfer_readdata(clt->socket, clt->max_msg_size, *connec_faces_idx, sizeof(int) * ((*n_faces)+1));

  // read connectivity
  CWP_transfer_readdata(clt->socket, clt->max_msg_size, *connec_faces, sizeof(int) * (*connec_faces_idx)[(*n_faces)]);

  // read connectivity index
  CWP_transfer_readdata(clt->socket, clt->max_msg_size, *connec_cells_idx, sizeof(int) * ((*n_elts)+1));

  // read connectivity
  CWP_transfer_readdata(clt->socket, clt->max_msg_size, *connec_cells, sizeof(int) * (*connec_cells_idx)[(*n_elts)]);

  // read global number
  CWP_transfer_readdata(clt->socket, clt->max_msg_size, *global_num, sizeof(CWP_g_num_t) * (*n_elts));
}

void
CWP_client_Mesh_interf_del
(
 const char *local_code_name,
 const char *cpl_id
)
{
  t_message msg;

  // verbose
  if (clt->flags & CWP_CLIENTFLAG_VERBOSE) {
    log_trace("CWP:Client initiating CWP_Mesh_interf_del\n");
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_MESH_INTERF_DEL);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Mesh_interf_del failed to send message header\n");
  }

  // send local code name
  write_name(local_code_name);

  // send coupling identifier
  write_name(cpl_id);
}

void
CWP_client_Mesh_interf_from_cellface_set
(
 const char           *local_code_name,
 const char           *cpl_id,
 const int             i_part,
 const int             n_cells,
 int                   cell_face_idx[],
 int                   cell_face[],
 const int             n_faces,
 int                   face_vtx_idx[],
 int                   face_vtx[],
 CWP_g_num_t           global_num[]
)
{
  t_message msg;

  // verbose
  if (clt->flags & CWP_CLIENTFLAG_VERBOSE) {
    log_trace("CWP:Client initiating CWP_Mesh_interf_from_cellface_set\n");
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_MESH_INTERF_FROM_CELLFACE_SET);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Mesh_interf_from_cellface_set failed to send message header\n");
  }

  // send local code name
  write_name(local_code_name);

  // send coupling identifier
  write_name(cpl_id);

  // send i_part
  int endian_i_part = i_part;
  CWP_swap_endian_4bytes(&endian_i_part, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_i_part, sizeof(int));

  // send number of cells
  int endian_n_cells = n_cells;
  CWP_swap_endian_4bytes(&endian_n_cells, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_n_cells, sizeof(int));

  // send cell->face connectivity index
  int *endian_cell_face_idx = (int *) malloc(sizeof(int) * (n_cells+1));
  memcpy(endian_cell_face_idx, cell_face_idx, sizeof(int) * (n_cells+1));
  CWP_swap_endian_4bytes(endian_cell_face_idx, n_cells + 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) endian_cell_face_idx, sizeof(int) * (n_cells+1));

  // send cell->face connectivity
  int *endian_cell_face = (int *) malloc(sizeof(int) * cell_face_idx[n_cells]);
  memcpy(endian_cell_face, cell_face, sizeof(int) * cell_face_idx[n_cells]);
  CWP_swap_endian_4bytes(endian_cell_face, cell_face_idx[n_cells]);
  CWP_transfer_writedata(clt->socket, clt->max_msg_size, (void *) endian_cell_face, sizeof(int) * cell_face_idx[n_cells]);

  // send number of faces
  int endian_n_faces = n_faces;
  CWP_swap_endian_4bytes(&endian_n_faces, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_n_faces, sizeof(int));

  // send face->vertex connectivity index
  int *endian_face_vtx_idx = (int *) malloc(sizeof(int) * (n_faces+1));
  memcpy(endian_face_vtx_idx, face_vtx_idx, sizeof(int) * (n_faces+1));
  CWP_swap_endian_4bytes(endian_face_vtx_idx, n_faces + 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) endian_face_vtx_idx, sizeof(int) * (n_faces+1));

  // send face->vertex connectivity
  int *endian_face_vtx = (int *) malloc(sizeof(int) * face_vtx_idx[n_faces]);
  memcpy(endian_face_vtx, face_vtx, sizeof(int) * face_vtx_idx[n_faces]);
  CWP_swap_endian_4bytes(endian_face_vtx, face_vtx_idx[n_faces]);
  CWP_transfer_writedata(clt->socket, clt->max_msg_size, (void *) endian_face_vtx, sizeof(int) * face_vtx_idx[n_faces]);

  // send global number
  CWP_g_num_t  *endian_global_num = (CWP_g_num_t  *) malloc(sizeof(CWP_g_num_t) * n_cells);
  memcpy(endian_global_num, global_num, sizeof(CWP_g_num_t) * n_cells);
  if (sizeof(CWP_g_num_t) == 4) {
    CWP_swap_endian_4bytes((int *) endian_global_num, n_cells);
  }
  else  if (sizeof(CWP_g_num_t) == 8) {
    CWP_swap_endian_8bytes((double *) endian_global_num, n_cells);
  }
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) endian_global_num, sizeof(CWP_g_num_t) * n_cells);

  // free
  free(endian_cell_face_idx);
  free(endian_cell_face);
  free(endian_face_vtx_idx);
  free(endian_face_vtx);
  free(endian_global_num);
}

void
CWP_client_Mesh_interf_from_faceedge_set
(
 const char           *local_code_name,
 const char           *cpl_id,
 const int             i_part,
 const int             n_faces,
 int                   face_edge_idx[],
 int                   face_edge[],
 const int             n_edges,
 int                   edge_vtx_idx[],
 int                   edge_vtx[],
 CWP_g_num_t           global_num[]
)
{
  t_message msg;

  // verbose
  if (clt->flags & CWP_CLIENTFLAG_VERBOSE) {
    log_trace("CWP:Client initiating CWP_Mesh_interf_from_faceedge_set\n");
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_MESH_INTERF_FROM_FACEEDGE_SET);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Mesh_interf_from_faceedge_set failed to send message header\n");
  }

  // send local code name
  write_name(local_code_name);

  // send coupling identifier
  write_name(cpl_id);

  // send i_part
  int endian_i_part = i_part;
  CWP_swap_endian_4bytes(&endian_i_part, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_i_part, sizeof(int));

  // send number of faces
  int endian_n_faces = n_faces;
  CWP_swap_endian_4bytes(&endian_n_faces, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_n_faces, sizeof(int));

  // send face->vertex connectivity index
  int *endian_face_edge_idx = (int *) malloc(sizeof(int) * (n_faces+1));
  memcpy(endian_face_edge_idx, face_edge_idx, sizeof(int) * (n_faces+1));
  CWP_swap_endian_4bytes(endian_face_edge_idx, n_faces + 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) endian_face_edge_idx, sizeof(int) * (n_faces+1));

  // send face->vertex connectivity
  int *endian_face_edge = (int *) malloc(sizeof(int) * face_edge_idx[n_faces]);
  memcpy(endian_face_edge, face_edge, sizeof(int) * face_edge_idx[n_faces]);
  CWP_swap_endian_4bytes(endian_face_edge, face_edge_idx[n_faces]);
  CWP_transfer_writedata(clt->socket, clt->max_msg_size, (void *) endian_face_edge, sizeof(int) * face_edge_idx[n_faces]);

  // send number of edges
  int endian_n_edges = n_edges;
  CWP_swap_endian_4bytes(&endian_n_edges, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_n_edges, sizeof(int));

  // send edge->vertex connectivity index
  int *endian_edge_vtx_idx = (int *) malloc(sizeof(int) * (n_edges+1));
  memcpy(endian_edge_vtx_idx, edge_vtx_idx, sizeof(int) * (n_edges+1));
  CWP_swap_endian_4bytes(endian_edge_vtx_idx, n_edges + 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) endian_edge_vtx_idx, sizeof(int) * (n_edges+1));

  // send edge->vertex connectivity
  int *endian_edge_vtx = (int *) malloc(sizeof(int) * edge_vtx_idx[n_edges]);
  memcpy(endian_edge_vtx, edge_vtx, sizeof(int) * edge_vtx_idx[n_edges]);
  CWP_swap_endian_4bytes(endian_edge_vtx, edge_vtx_idx[n_edges]);
  CWP_transfer_writedata(clt->socket, clt->max_msg_size, (void *) endian_edge_vtx, sizeof(int) * edge_vtx_idx[n_edges]);

  // send global number
  CWP_g_num_t  *endian_global_num = (CWP_g_num_t  *) malloc(sizeof(CWP_g_num_t) * n_faces);
  memcpy(endian_global_num, global_num, sizeof(CWP_g_num_t) * n_faces);
  if (sizeof(CWP_g_num_t) == 4) {
    CWP_swap_endian_4bytes((int *) endian_global_num, n_faces);
  }
  else  if (sizeof(CWP_g_num_t) == 8) {
    CWP_swap_endian_8bytes((double *) endian_global_num, n_faces);
  }
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) endian_global_num, sizeof(CWP_g_num_t) * n_faces);

  // free
  free(endian_face_edge_idx);
  free(endian_face_edge);
  free(endian_edge_vtx_idx);
  free(endian_edge_vtx);
  free(endian_global_num);
}

void
CWP_client_Field_create
(
 const char                  *local_code_name,
 const char                  *cpl_id,
 const char                  *field_id,
 const CWP_Type_t             data_type,
 const CWP_Field_storage_t    storage,
 const int                    n_component,
 const CWP_Dof_location_t     target_location,
 const CWP_Field_exch_t       exch_type,
 const CWP_Status_t           visu_status
)
{
  t_message msg;

  // verbose
  if (clt->flags & CWP_CLIENTFLAG_VERBOSE) {
    log_trace("CWP:Client initiating CWP_Field_create\n");
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_FIELD_CREATE);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Field_create failed to send message header\n");
  }

  // send local code name
  write_name(local_code_name);

  // send coupling identifier
  write_name(cpl_id);

  // send field identifier
  write_name(field_id);

  // send data type
  CWP_Type_t endian_data_type = data_type;
  CWP_swap_endian_4bytes((int *) &endian_data_type, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_data_type, sizeof(CWP_Type_t));

  // send storage mode
  CWP_Field_storage_t endian_storage = storage;
  CWP_swap_endian_4bytes((int *) &endian_storage, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_storage, sizeof(CWP_Field_storage_t));

  // send number of components
  int endian_n_component = n_component;
  CWP_swap_endian_4bytes(&endian_n_component, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_n_component, sizeof(int));

  // send target location
  CWP_Dof_location_t endian_target_location = target_location;
  CWP_swap_endian_4bytes((int *) &endian_target_location, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_target_location, sizeof(CWP_Dof_location_t));

  // send exchange type
  CWP_Field_exch_t endian_exch_type = exch_type;
  CWP_swap_endian_4bytes((int *) &endian_exch_type, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_exch_type, sizeof(CWP_Field_exch_t));

  // send visualisation status
  CWP_Status_t endian_visu_status = visu_status;
  CWP_swap_endian_4bytes((int *) &endian_visu_status, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_visu_status, sizeof(CWP_Status_t));
}

void
CWP_client_Field_data_set
(
 const char              *local_code_name,
 const char              *cpl_id,
 const char              *field_id,
 const int                i_part,
 const CWP_Field_map_t    map_type,
 double                   data[],
 const int                n_component,
 int                      n_dof
)
{
  t_message msg;

  // verbose
  if (clt->flags & CWP_CLIENTFLAG_VERBOSE) {
    log_trace("CWP:Client initiating CWP_Field_data_set\n");
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_FIELD_DATA_SET);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Field_data_set failed to send message header\n");
  }

  // send local code name
  write_name(local_code_name);

  // send coupling identifier
  write_name(cpl_id);

  // send field identifier
  write_name(field_id);

  // send i_part
  int endian_i_part = i_part;
  CWP_swap_endian_4bytes(&endian_i_part, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_i_part, sizeof(int));

  // send map type
  CWP_Field_map_t endian_map_type = map_type;
  CWP_swap_endian_4bytes((int *) &endian_map_type, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_map_type, sizeof(CWP_Field_map_t));

  // send map with data
  int size = n_component * n_dof;
  int endian_size = size;
  CWP_swap_endian_4bytes(&endian_size, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_size, sizeof(int));
  double *endian_data = (double *) malloc(sizeof(double) * size);
  memcpy(endian_data, data, sizeof(double) * size);
  CWP_swap_endian_8bytes(endian_data, size);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) endian_data, sizeof(double) * size);

  // free
  free(endian_data);
}

int
CWP_client_Field_n_component_get
(
 const char      *local_code_name,
 const char      *cpl_id,
 const char      *field_id
)
{
  t_message msg;

  // verbose
  if (clt->flags & CWP_CLIENTFLAG_VERBOSE) {
    log_trace("CWP:Client initiating CWP_Field_n_component_get\n");
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_FIELD_N_COMPONENT_GET);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Field_n_component_get failed to send message header\n");
  }

  // send local code name
  write_name(local_code_name);

  // send coupling identifier
  write_name(cpl_id);

  // send field identifier
  write_name(field_id);

  // read number of components
  int n_components = -1;
  CWP_transfer_readdata(clt->socket, clt->max_msg_size, &n_components, sizeof(int));

  return n_components;
}

CWP_Dof_location_t
CWP_client_Field_target_dof_location_get
(
 const char      *local_code_name,
 const char      *cpl_id,
 const char      *field_id
)
{
  t_message msg;

  // verbose
  if (clt->flags & CWP_CLIENTFLAG_VERBOSE) {
    log_trace("CWP:Client initiating CWP_Field_target_dof_location_get\n");
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_FIELD_TARGET_DOF_LOCATION_GET);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Field_target_dof_location_get failed to send message header\n");
  }

  // send local code name
  write_name(local_code_name);

  // send coupling identifier
  write_name(cpl_id);

  // send field identifier
  write_name(field_id);

  // read location of target degrees of freedom
  CWP_Dof_location_t dof_location = (CWP_Dof_location_t) -1;
  CWP_transfer_readdata(clt->socket, clt->max_msg_size, &dof_location, sizeof(CWP_Dof_location_t));

  return dof_location;
}

CWP_Field_storage_t
CWP_client_Field_storage_get
(
 const char      *local_code_name,
 const char      *cpl_id         ,
 const char      *field_id
)
{
  t_message msg;

  // verbose
  if (clt->flags & CWP_CLIENTFLAG_VERBOSE) {
    log_trace("CWP:Client initiating CWP_Field_storage_get\n");
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_FIELD_STORAGE_GET);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Field_storage_get failed to send message header\n");
  }

  // send local code name
  write_name(local_code_name);

  // send coupling identifier
  write_name(cpl_id);

  // send field identifier
  write_name(field_id);

  // read field storage type
  CWP_Field_storage_t storage_type = (CWP_Field_storage_t) -1;
  CWP_transfer_readdata(clt->socket, clt->max_msg_size, &storage_type, sizeof(CWP_Field_storage_t));

  return storage_type;
}

void
CWP_client_Field_del
(
 const char      *local_code_name,
 const char      *cpl_id         ,
 const char      *field_id
)
{
  t_message msg;

  // verbose
  if (clt->flags & CWP_CLIENTFLAG_VERBOSE) {
    log_trace("CWP:Client initiating CWP_Field_del\n");
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_FIELD_DEL);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Field_del failed to send message header\n");
  }

  // send local code name
  write_name(local_code_name);

  // send coupling identifier
  write_name(cpl_id);

  // send field identifier
  write_name(field_id);
}

void
CWP_client_Field_issend
(
 const char     *local_code_name,
 const char     *cpl_id,
 const char     *src_field_id
)
{
  t_message msg;

  // verbose
  if (clt->flags & CWP_CLIENTFLAG_VERBOSE) {
    log_trace("CWP:Client initiating CWP_Field_issend\n");
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_FIELD_ISSEND);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Field_issend failed to send message header\n");
  }

  // send local code name
  write_name(local_code_name);

  // send coupling identifier
  write_name(cpl_id);

  // send source field identifier
  write_name(src_field_id);
}

void
CWP_client_Field_irecv
(
 const char        *local_code_name,
 const char        *cpl_id,
 const char        *tgt_field_id
)
{
  t_message msg;

  // verbose
  if (clt->flags & CWP_CLIENTFLAG_VERBOSE) {
    log_trace("CWP:Client initiating CWP_Field_irecv\n");
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_FIELD_IRECV);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Field_irecv failed to send message header\n");
  }

  // send local code name
  write_name(local_code_name);

  // send coupling identifier
  write_name(cpl_id);

  // send target field identifier
  write_name(tgt_field_id);
}

void
CWP_client_Field_wait_issend
(
 const char  *local_code_name,
 const char  *cpl_id,
 const char  *src_field_id
)
{
  t_message msg;

  // verbose
  if (clt->flags & CWP_CLIENTFLAG_VERBOSE) {
    log_trace("CWP:Client initiating CWP_Field_wait_issend\n");
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_FIELD_WAIT_ISSEND);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Field_wait_issend failed to send message header\n");
  }

  // send local code name
  write_name(local_code_name);

  // send coupling identifier
  write_name(cpl_id);

  // send source field identifier
  write_name(src_field_id);
}

void
CWP_client_Field_wait_irecv
(
 const char              *local_code_name,
 const char              *cpl_id,
 const char              *tgt_field_id,
 const int                i_part,
 const CWP_Field_map_t    map_type,
 const int                n_component,
 int                      n_dof,
 double                 **data
)
{
  t_message msg;

  // verbose
  if (clt->flags & CWP_CLIENTFLAG_VERBOSE) {
    log_trace("CWP:Client initiating CWP_Field_wait_irecv\n");
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_FIELD_WAIT_IRECV);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Field_wait_irecv failed to send message header\n");
  }

  // send local code name
  write_name(local_code_name);

  // send coupling identifier
  write_name(cpl_id);

  // send target field identifier
  write_name(tgt_field_id);

  // send data needed to retreive field
  // send i_part
  int endian_i_part = i_part;
  CWP_swap_endian_4bytes(&endian_i_part, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_i_part, sizeof(int));

  // send map type
  CWP_Field_map_t endian_map_type = map_type;
  CWP_swap_endian_4bytes((int *) &endian_map_type, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_map_type, sizeof(CWP_Field_map_t));

  // send size
  int size = n_component * n_dof;
  int endian_size = size;
  CWP_swap_endian_4bytes(&endian_size, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_size, sizeof(int));

  // read irecv Field_data
  CWP_transfer_readdata(clt->socket, clt->max_msg_size, *data, sizeof(double) * size);
}

void
CWP_client_Interp_from_location_unset
(
 const char                 *local_code_name,
 const char                 *cpl_id,
 const char                 *src_field_id
)
{
  t_message msg;

  // verbose
  if (clt->flags & CWP_CLIENTFLAG_VERBOSE) {
    log_trace("CWP:Client initiating CWP_Interp_from_location_unset\n");
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_INTERP_FROM_LOCATION_UNSET);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Interp_from_location_unset failed to send message header\n");
  }

  // send local code name
  write_name(local_code_name);

  // send coupling identifier
  write_name(cpl_id);

  // send source field identifier
  write_name(src_field_id);
}

void
CWP_client_Interp_from_location_set
(
 const char                 *local_code_name,
 const char                 *cpl_id,
 const char                 *src_field_id,
 CWP_Interp_from_location_t  fct // TO DO: deal with other than double data_type
)
{
  // t_message msg;

  // // verbose
  // if (clt->flags & CWP_CLIENTFLAG_VERBOSE) {
  //   log_trace("CWP:Client initiating CWP_Interp_from_location_set\n");
  // }

  // // create message
  // NEWMESSAGE(msg, CWP_MSG_CWP_INTERP_FROM_LOCATION_SET);

  // // send message
  // if (CWP_client_send_msg(&msg) != 0) {
  //   PDM_error(__FILE__, __LINE__, 0, "CWP_client_Interp_from_location_set failed to send message header\n");
  // }

  // // send local code name
  // write_name(local_code_name);

  // // send coupling identifier
  // write_name(cpl_id);

  // // send source field identifier
  // write_name(src_field_id);

  // // send struct fct
  // // --> interface_type
  // int endian_interface_type = fct->interface_type;
  // CWP_swap_endian_4bytes( (int *) &endian_interface_type, 1);
  // CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_interface_type, sizeof(int));

  // // --> code_name
  // write_name(fct->code_name);

  // // --> src_n_block
  // int endian_src_n_block = fct->src_n_block;
  // CWP_swap_endian_4bytes( (int *) &endian_src_n_block, 1);
  // CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_src_n_block, sizeof(int));

  // // --> src_blocks_type
  // CWP_Block_t *endian_src_blocks_type = (CWP_Block_t *) malloc(sizeof(CWP_Block_t) * fct->src_n_block);
  // memcpy(endian_src_blocks_type, fct->src_blocks_type, sizeof(CWP_Block_t) * fct->src_n_block);
  // CWP_swap_endian_4bytes( (int *) endian_src_blocks_type, fct->src_n_block);
  // CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) endian_src_blocks_type, sizeof(CWP_Block_t) * fct->src_n_block);

  // // --> src_i_part
  // int endian_src_i_part = fct->src_i_part;
  // CWP_swap_endian_4bytes( (int *) &endian_src_i_part, 1);
  // CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_src_i_part, sizeof(int));

  // // --> src_n_vtx
  // int endian_src_n_vtx = fct->src_n_vtx;
  // CWP_swap_endian_4bytes( (int *) &endian_src_n_vtx, 1);
  // CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_src_n_vtx, sizeof(int));

  // // --> src_vtx_coords
  // double *endian_src_vtx_coords = (double *) malloc(sizeof(double) *  3 * fct->src_n_vtx);
  // memcpy(endian_src_vtx_coords, fct->src_vtx_coords, sizeof(double) * 3 *  fct->src_n_vtx);
  // CWP_swap_endian_8bytes(endian_src_vtx_coords, 3 *  fct->src_n_vtx);
  // CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) endian_src_vtx_coords, sizeof(double) * 3 *  fct->src_n_vtx);

  // // --> src_vtx_global_num
  // CWP_g_num_t  *endian_src_vtx_global_num = (CWP_g_num_t  *) malloc(sizeof(CWP_g_num_t) * fct->src_n_vtx);
  // memcpy(endian_src_vtx_global_num, fct->src_vtx_global_num, sizeof(CWP_g_num_t) * fct->src_n_vtx);
  // CWP_swap_endian_8bytes(endian_src_vtx_global_num, fct->src_n_vtx); // TO DO 32 bit machine
  // CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) endian_src_vtx_global_num, sizeof(CWP_g_num_t) * fct->src_n_vtx);

  // // --> src_n_elts
  // int endian_src_n_elts = fct->src_n_elts;
  // CWP_swap_endian_4bytes( (int *) &endian_src_n_elts, 1);
  // CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_src_n_elts, sizeof(int));

  // // --> src_id_block
  // int *endian_src_id_block = (int *) malloc(sizeof(int) * fct->src_n_elts);
  // memcpy(endian_src_id_block, fct->src_id_block, sizeof(int) * fct->src_n_elts);
  // CWP_swap_endian_4bytes( (int *) endian_src_id_block, fct->src_n_elts);
  // CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) endian_src_id_block, sizeof(int) * fct->src_n_elts);

  // // --> src_elt_in_block
  // int *endian_src_elt_in_block = (int *) malloc(sizeof(int) * fct->src_n_elts);
  // memcpy(endian_src_elt_in_block, fct->src_elt_in_block, sizeof(int) * fct->src_n_elts);
  // CWP_swap_endian_4bytes( (int *) endian_src_elt_in_block, fct->src_n_elts);
  // CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) endian_src_elt_in_block, sizeof(int) * fct->src_n_elts);

  // // --> src_elt_vtx_idx
  // int *endian_src_elt_vtx_idx = (int *) malloc(sizeof(int) * (fct->src_n_elts + 1));
  // memcpy(endian_src_elt_vtx_idx, fct->src_elt_vtx_idx, sizeof(int) * (fct->src_n_elts + 1));
  // CWP_swap_endian_4bytes( (int *) endian_src_elt_vtx_idx, (fct->src_n_elts + 1));
  // CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) endian_src_elt_vtx_idx, sizeof(int) * (fct->src_n_elts + 1));

  // // --> src_elt_vtx
  // int *endian_src_elt_vtx = (int *) malloc(sizeof(int) * fct->src_elt_vtx_idx[fct->src_n_elts]);
  // memcpy(endian_src_elt_vtx, fct->src_elt_vtx, sizeof(int) * fct->src_elt_vtx_idx[fct->src_n_elts]);
  // CWP_swap_endian_4bytes( (int *) endian_src_elt_vtx, fct->src_elt_vtx_idx[fct->src_n_elts]);
  // CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) endian_src_elt_vtx, sizeof(int) * fct->src_elt_vtx_idx[fct->src_n_elts]);

  // // --> src_elts_global_num
  // CWP_g_num_t  *endian_src_elts_global_num = (CWP_g_num_t  *) malloc(sizeof(CWP_g_num_t) * fct->src_n_elts);
  // memcpy(endian_src_elts_global_num, fct->src_elts_global_num, sizeof(CWP_g_num_t) * fct->src_n_elts);
  // CWP_swap_endian_8bytes(endian_src_elts_global_num, fct->src_n_elts); // TO DO 32 bit machine
  // CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) endian_src_elts_global_num, sizeof(CWP_g_num_t) * fct->src_n_elts);

  // // --> tgt_n_pts
  // int endian_tgt_n_pts = fct->tgt_n_pts;
  // CWP_swap_endian_4bytes( (int *) &endian_tgt_n_pts, 1);
  // CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_tgt_n_pts, sizeof(int));

  // // --> tgt_pts_elt_idx
  // int *endian_tgt_pts_elt_idx = (int *) malloc(sizeof(int) * (fct->src_n_elts + 1));
  // memcpy(endian_tgt_pts_elt_idx, fct->tgt_pts_elt_idx, sizeof(int) * (fct->src_n_elts + 1));
  // CWP_swap_endian_4bytes( (int *) endian_tgt_pts_elt_idx, (fct->src_n_elts + 1));
  // CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) endian_tgt_pts_elt_idx, sizeof(int) * (fct->src_n_elts + 1));

  // // --> tgt_pts_coords
  // double *endian_tgt_pts_coords = (double *) malloc(sizeof(double) *  3 * fct->tgt_n_pts);
  // memcpy(endian_tgt_pts_coords, fct->tgt_pts_coords, sizeof(double) * 3 *  fct->tgt_n_pts);
  // CWP_swap_endian_8bytes(endian_tgt_pts_coords, 3 *  fct->tgt_n_pts);
  // CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) endian_tgt_pts_coords, sizeof(double) * 3 *  fct->tgt_n_pts);

  // // --> tgt_pts_dist
  // double *endian_tgt_pts_dist = (double *) malloc(sizeof(double) * fct->tgt_n_pts);
  // memcpy(endian_tgt_pts_dist, fct->tgt_pts_dist, sizeof(double) * fct->tgt_n_pts);
  // CWP_swap_endian_8bytes(endian_tgt_pts_dist, fct->tgt_n_pts);
  // CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) endian_tgt_pts_dist, sizeof(double) * fct->tgt_n_pts);

  // // --> tgt_pts_uvw
  // // TO DO: size = dim_interface * tgt_n_pts get dim_interface??
  // double *endian_tgt_pts_uvw = (double *) malloc(sizeof(double) * fct->tgt_n_pts);
  // memcpy(endian_tgt_pts_uvw, fct->tgt_pts_uvw, sizeof(double) * fct->tgt_n_pts);
  // CWP_swap_endian_8bytes(endian_tgt_pts_uvw, fct->tgt_n_pts);
  // CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) endian_tgt_pts_uvw, sizeof(double) * fct->tgt_n_pts);

  // // --> tgt_pts_weights_idx
  // int *endian_tgt_pts_weights_idx = (int *) malloc(sizeof(int) * (fct->tgt_n_pts + 1));
  // memcpy(endian_tgt_pts_weights_idx, fct->tgt_pts_weights_idx, sizeof(int) * (fct->tgt_n_pts + 1));
  // CWP_swap_endian_4bytes( (int *) endian_tgt_pts_weights_idx, (fct->tgt_n_pts + 1));
  // CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) endian_tgt_pts_weights_idx, sizeof(int) * (fct->tgt_n_pts + 1));

  // // --> tgt_pts_weights
  // double *endian_tgt_pts_weights = (double *) malloc(sizeof(double) * fct->tgt_pts_weights_idx[fct->tgt_n_pts]);
  // memcpy(endian_tgt_pts_weights, fct->tgt_pts_weights, sizeof(double) * fct->tgt_pts_weights_idx[fct->tgt_n_pts]);
  // CWP_swap_endian_8bytes(endian_tgt_pts_weights, fct->tgt_pts_weights_idx[fct->tgt_n_pts]);
  // CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) endian_tgt_pts_weights, sizeof(double) * fct->tgt_pts_weights_idx[fct->tgt_n_pts]);

  // // --> stride
  // int endian_stride = fct->stride;
  // CWP_swap_endian_4bytes( (int *) &endian_stride, 1);
  // CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_stride, sizeof(int));

  // // --> src_field_dof_location
  // CWP_Dof_location_t endian_src_field_dof_location = fct->src_field_dof_location;
  // CWP_swap_endian_4bytes( (int *) &endian_src_field_dof_location, 1);
  // CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_src_field_dof_location, sizeof(CWP_Dof_location_t));

  // // --> src_field
  // double *endian_src_field = (double *) malloc(sizeof(double) * fct->src_n_vtx * fct->stride);
  // memcpy(endian_src_field, fct->src_field, sizeof(double) * fct->src_n_vtx * fct->stride);
  // CWP_swap_endian_8bytes(endian_src_field, fct->src_n_vtx * fct->stride);
  // CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) endian_src_field, sizeof(double) * fct->src_n_vtx * fct->stride);

  // // --> tgt_field
  // double *endian_tgt_field = (double *) malloc(sizeof(double) * fct->tgt_n_pts * fct->stride);
  // memcpy(endian_tgt_field, fct->tgt_field, sizeof(double) * fct->tgt_n_pts * fct->stride);
  // CWP_swap_endian_8bytes(endian_tgt_field, fct->tgt_n_pts * fct->stride);
  // CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) endian_tgt_field, sizeof(double) * fct->tgt_n_pts * fct->stride);

  // // free
  // free((void *) local_code_name);
  // free((void *) cpl_id);
  // free((void *) src_field_id);
  // free(endian_src_blocks_type);
  // free(endian_src_vtx_coords);
  // free(endian_src_vtx_global_num);
  // free(endian_src_id_block);
  // free(endian_src_elt_in_block);
  // free(endian_src_elt_vtx_idx);
  // free(endian_src_elt_vtx);
  // free(endian_src_elts_global_num);
  // free(endian_tgt_pts_elt_idx);
  // free(endian_tgt_pts_coords);
  // free(endian_tgt_pts_dist);
  // free(endian_tgt_pts_uvw);
  // free(endian_tgt_pts_weights_idx);
  // free(endian_tgt_pts_weights);
  // free(endian_src_field);
  // free(endian_tgt_field);
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
  struct sockaddr_in *server_addr = (struct sockaddr_in *) malloc(sizeof(struct sockaddr_in));
  socklen_t max_msg_size_len = 0;
  int il_cl_endian;

  clt = (t_client *) malloc(sizeof(t_client));
  memset(clt,0,sizeof(t_client));
  strncpy(clt->server_name, server_name,sizeof(clt->server_name));
  clt->server_port=server_port;
  clt->flags=flags;

  // clt_cwp init
  clt_cwp = (p_cwp) malloc(sizeof(t_cwp));
  memset(clt_cwp, 0, sizeof(t_cwp));
  clt_cwp->code = (t_code *) malloc(sizeof(t_code));

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

  server_addr->sin_family = AF_INET;
  server_addr->sin_port = htons(server_port);
  server_addr->sin_addr = *((struct in_addr *)host->h_addr);
  bzero(&(server_addr->sin_zero),8);

  // connect
  while (connect(clt->socket, (struct sockaddr *) server_addr,
                 sizeof(struct sockaddr)) == -1) {}

  // verbose
  if (clt->flags & CWP_CLIENTFLAG_VERBOSE) {
    log_trace("CWP:Client connected on %s port %i \n", server_name, server_port);
  }

  // get maximum message size
  clt->max_msg_size = 0;
  max_msg_size_len  = sizeof(int);
  getsockopt(clt->socket, SOL_SOCKET, SO_RCVBUF, (void *)&clt->max_msg_size, &max_msg_size_len);
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
    log_trace("CWP:client endian %i server endian %i\n",clt->client_endianess,clt->server_endianess);
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

  // free
  if (clt_cwp->code_names != NULL) {
    for (int i = 0; i < clt_cwp->n_code_names; i++) {
      free((void *) clt_cwp->code_names[i]);
    }
    free((void **) clt_cwp->code_names);
  }

  if (clt_cwp->loc_code_names != NULL) {
    for (int i = 0; i < clt_cwp->n_loc_code_names; i++) {
      free((void *) clt_cwp->loc_code_names[i]);
    }
    free((void **) clt_cwp->loc_code_names);
  }

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
