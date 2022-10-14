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
  *name = realloc(*name, name_size);
  CWP_transfer_readdata(clt->socket,clt->max_msg_size,(void*) *name, name_size);
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

  // send code name
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
    char *char_value = malloc(sizeof(char));
    read_name(&char_value);
    memcpy(value, char_value, strlen(char_value));
    break;

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
  int endian_data_type = data_type;
  CWP_swap_endian_4bytes(&endian_data_type, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_data_type, sizeof(int));

  switch (data_type) {

  case CWP_DOUBLE: ;
    double endian_double_value = * ((double *) value);
    CWP_swap_endian_8bytes(&endian_double_value, 1);
    CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_double_value, sizeof(double));
    break;

  case CWP_INT: ;
    int endian_int_value = * ((int *) value);
    CWP_swap_endian_4bytes(&endian_int_value, 1);
    CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_int_value, sizeof(int));
    break;

  case CWP_CHAR:
    write_name(value);
    break;

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
  int endian_data_type = data_type;
  CWP_swap_endian_4bytes(&endian_data_type, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_data_type, sizeof(int));
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
  int endian_data_type = data_type;
  CWP_swap_endian_4bytes(&endian_data_type, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size, (void*) &endian_data_type, sizeof(int));

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
  int endian_data_type = data_type;
  CWP_swap_endian_4bytes(&endian_data_type, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size, (void*) &endian_data_type, sizeof(int));

  // read n_param
  CWP_transfer_readdata(clt->socket, clt->max_msg_size, (void*) nParam, sizeof(int));

  // read param names
  *paramNames = malloc(sizeof(*nParam));
  for (int i = 0; i < *nParam; i++) {
    (*paramNames)[i] = malloc(sizeof(char));
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
  int endian_data_type = data_type;
  CWP_swap_endian_4bytes(&endian_data_type, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_data_type, sizeof(int));

  // read bool
  int is_param = -1;
  CWP_transfer_readdata(clt->socket, clt->max_msg_size, (void*) &is_param, sizeof(int));

  return is_param;
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
  int endian_entities_dim = entities_dim;
  CWP_swap_endian_4bytes(&endian_entities_dim, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_entities_dim, sizeof(int));

  // send communication type
  int endian_comm_type = comm_type;
  CWP_swap_endian_4bytes(&endian_comm_type, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_comm_type, sizeof(int));

  // send spatial interpolation
  int endian_spatial_interp = spatial_interp;
  CWP_swap_endian_4bytes(&endian_spatial_interp, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_spatial_interp, sizeof(int));

  // send number of partitions
  int endian_n_part = n_part;
  CWP_swap_endian_4bytes(&endian_n_part, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_n_part, sizeof(int));

  // send displacement type
  int endian_displacement = displacement;
  CWP_swap_endian_4bytes(&endian_displacement, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_displacement, sizeof(int));

  // send time exchange type
  int endian_recv_freq_type = recv_freq_type;
  CWP_swap_endian_4bytes(&endian_recv_freq_type, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_recv_freq_type, sizeof(int));
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
  int endian_format = format;
  CWP_swap_endian_4bytes(&endian_format, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_format, sizeof(int));

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
  int endian_state = state;
  CWP_swap_endian_4bytes(&endian_state, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_state, sizeof(int));
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
  int state = -1;
  CWP_transfer_readdata(clt->socket, clt->max_msg_size, (void*) &state, sizeof(int));

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
  int nb_codes = -1;
  CWP_transfer_readdata(clt->socket, clt->max_msg_size, (void*) &nb_codes, sizeof(int));

  return nb_codes;
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
  int nb_codes = -1;
  CWP_transfer_readdata(clt->socket, clt->max_msg_size, (void*) &nb_codes, sizeof(int));
  char **code_names = malloc(nb_codes);
  for (int i = 0; i < nb_codes; i++) {
    code_names[i] = malloc(sizeof(char));
    read_name(&code_names[i]);
  }

  return code_names;
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
  int nb_local_codes = -1;
  CWP_transfer_readdata(clt->socket, clt->max_msg_size, (void*) &nb_local_codes, sizeof(int));

  return nb_local_codes;
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
  int nb_local_codes = -1;
  CWP_transfer_readdata(clt->socket, clt->max_msg_size, (void*) &nb_local_codes, sizeof(int));
  char **code_local_names = malloc(nb_local_codes);
  for (int i = 0; i < nb_local_codes; i++) {
    code_local_names[i] = malloc(sizeof(char));
    read_name(&code_local_names[i]);
  }

  return code_local_names;
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
