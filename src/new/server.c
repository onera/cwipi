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

static void read_name(char **name,  p_server svr) {
  int name_size;
  CWP_transfer_readdata(svr->connected_socket, svr->max_msg_size, (void*) &name_size, sizeof(int));
  *name = realloc(*name, name_size);
  CWP_transfer_readdata(svr->connected_socket, svr->max_msg_size, (void*) *name, name_size);
}

static void write_name(char * name, p_server svr) {
  int name_size = strlen(name) + 1; // +1 for "\0"
  int endian_name_size = name_size;
  // CWP_swap_endian_4bytes(&endian_name_size, 1); // TO DO: mandatory ?
  CWP_transfer_writedata(svr->connected_socket, svr->max_msg_size, (void*) &endian_name_size, sizeof(int));
  CWP_transfer_writedata(svr->connected_socket, svr->max_msg_size, (void*) name, name_size);
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
  int            n_code;
  char         **code_names = NULL;
  CWP_Status_t  *is_active_rank = NULL;
  double        *time_init = NULL;

  // receive data
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  CWP_transfer_readdata(svr->connected_socket, svr->max_msg_size, &n_code, sizeof(int));

  if (n_code > 1) {
    PDM_error(__FILE__, __LINE__, 0, "CWIPI client-server not implemented yet for n_code > 1\n");
  }

  code_names = malloc(sizeof(char *) * n_code);
  for (int i = 0; i < n_code; i++) {
    code_names[i] = malloc(sizeof(char));
    read_name(&code_names[i], svr);
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
           svr->intra_comms);

  // wait all ranks have receive msg
  MPI_Barrier(svr->intra_comms[0]);

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
CWP_server_Finalize
(
 p_server                 svr
)
{
  // wait all ranks have receive msg
  MPI_Barrier(svr->intra_comms[0]);

  CWP_Finalize();
  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Param_lock
(
 p_server                 svr
)
{
  // wait all ranks have receive msg
  MPI_Barrier(svr->intra_comms[0]);

  // launch
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *code_name = malloc(sizeof(char));
  read_name(&code_name, svr);

  CWP_Param_lock((const char *) code_name);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Param_unlock
(
 p_server                 svr
)
{
  // wait all ranks have receive msg
  MPI_Barrier(svr->intra_comms[0]);

  // launch
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *code_name = malloc(sizeof(char));
  read_name(&code_name, svr);

  CWP_Param_unlock((const char *) code_name);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Param_add
(
 p_server                 svr
)
{
  // wait all ranks have receive msg
  MPI_Barrier(svr->intra_comms[0]);

  // read local code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read param name
  char *param_name = malloc(sizeof(char));
  read_name(&param_name, svr);

  // read initial value
  int data_type = -1;
  CWP_transfer_readdata(svr->connected_socket, svr->max_msg_size, &data_type, sizeof(int));

  switch (data_type) {

  case CWP_DOUBLE: ;
    double double_initial_value;
    CWP_transfer_readdata(svr->connected_socket, svr->max_msg_size, &double_initial_value, sizeof(double));
    CWP_Param_add(local_code_name,
                  param_name,
                  data_type,
                  &double_initial_value);
    break;

  case CWP_INT: ;
    double int_initial_value;
    CWP_transfer_readdata(svr->connected_socket, svr->max_msg_size, &int_initial_value, sizeof(int));
    CWP_Param_add(local_code_name,
                  param_name,
                  data_type,
                  &int_initial_value);
    break;

  case CWP_CHAR: ;
    char *char_initial_value = malloc(sizeof(char));
    read_name(&char_initial_value, svr);
    CWP_Param_add(local_code_name,
                  param_name,
                  data_type,
                  &char_initial_value);
    break;

  default:
    PDM_error(__FILE__, __LINE__, 0, "Received unknown CWP_Type_t %i\n", data_type);
  }

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Param_get
(
 p_server                 svr
)
{
  // wait all ranks have receive msg
  MPI_Barrier(svr->intra_comms[0]);

  // read local code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read param name
  char *param_name = malloc(sizeof(char));
  read_name(&param_name, svr);

  // read value data_type
  int data_type = -1;
  CWP_transfer_readdata(svr->connected_socket, svr->max_msg_size, &data_type, sizeof(int));

  // send value
  svr->state=CWP_SVRSTATE_SENDPGETDATA;
  switch (data_type) {

  case CWP_DOUBLE: ;
    double double_value;
    CWP_Param_get(local_code_name,
                  param_name,
                  data_type,
                  &double_value);
    double endian_double_value = double_value;
    // CWP_swap_endian_8bytes(&endian_double_value, 1); TO DO: mandatory here?
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size,(void*) &endian_double_value, sizeof(double));
    break;

  case CWP_INT: ;
    int int_value;
    CWP_Param_get(local_code_name,
                  param_name,
                  data_type,
                  &int_value);
    int endian_int_value = int_value;
    // CWP_swap_endian_4bytes(&endian_int_value, 1); TO DO: mandatory here?
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size,(void*) &endian_int_value, sizeof(int));
    break;

  case CWP_CHAR: ;
    char *char_value = malloc(sizeof(char));
    CWP_Param_get(local_code_name,
                  param_name,
                  data_type,
                  &char_value);
    write_name(char_value, svr);
    break;

  default:
    PDM_error(__FILE__, __LINE__, 0, "Received unknown CWP_Type_t %i\n", data_type);
  }

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Param_set
(
 p_server                 svr
)
{
  // wait all ranks have receive msg
  MPI_Barrier(svr->intra_comms[0]);

  // read local code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read param name
  char *param_name = malloc(sizeof(char));
  read_name(&param_name, svr);

  // read initial value
  int data_type = -1;
  CWP_transfer_readdata(svr->connected_socket, svr->max_msg_size, &data_type, sizeof(int));

  switch (data_type) {

  case CWP_DOUBLE: ;
    double double_initial_value;
    CWP_transfer_readdata(svr->connected_socket, svr->max_msg_size, &double_initial_value, sizeof(double));
    CWP_Param_set(local_code_name,
                  param_name,
                  data_type,
                  &double_initial_value);
    break;

  case CWP_INT: ;
    double int_initial_value;
    CWP_transfer_readdata(svr->connected_socket, svr->max_msg_size, &int_initial_value, sizeof(int));
    CWP_Param_set(local_code_name,
                  param_name,
                  data_type,
                  &int_initial_value);
    break;

  case CWP_CHAR: ;
    char *char_initial_value = malloc(sizeof(char));
    read_name(&char_initial_value, svr);
    CWP_Param_set(local_code_name,
                  param_name,
                  data_type,
                  &char_initial_value);
    break;

  default:
    PDM_error(__FILE__, __LINE__, 0, "Received unknown CWP_Type_t %i\n", data_type);
  }

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Param_del
(
 p_server                 svr
)
{
  // wait all ranks have receive msg
  MPI_Barrier(svr->intra_comms[0]);

  // read local code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read param name
  char *param_name = malloc(sizeof(char));
  read_name(&param_name, svr);

  // read initial value
  int data_type = -1;
  CWP_transfer_readdata(svr->connected_socket, svr->max_msg_size, &data_type, sizeof(int));

  // launch
  CWP_Param_del(local_code_name,
                param_name,
                data_type);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Param_n_get
(
 p_server                 svr
)
{
  // wait all ranks have receive msg
  MPI_Barrier(svr->intra_comms[0]);

  // read local code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *code_name = malloc(sizeof(char));
  read_name(&code_name, svr);

  // read initial value
  int data_type = -1;
  CWP_transfer_readdata(svr->connected_socket, svr->max_msg_size, &data_type, sizeof(int));

  // launch
  int n_param = CWP_Param_n_get(code_name,
                                data_type);

  // send n_param
  svr->state=CWP_SVRSTATE_SENDPGETDATA;
  CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, (void*) &n_param, sizeof(int));

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Param_list_get
(
 p_server                 svr
)
{
  // wait all ranks have receive msg
  MPI_Barrier(svr->intra_comms[0]);

  // read local code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *code_name = malloc(sizeof(char));
  read_name(&code_name, svr);

  // read initial value
  int data_type = -1;
  CWP_transfer_readdata(svr->connected_socket, svr->max_msg_size, &data_type, sizeof(int));

  // launch
  int nParam = -1;
  char **paramNames = NULL;
  CWP_Param_list_get(code_name,
                     data_type,
                     &nParam,
                     &paramNames);

  // send nParam
  svr->state=CWP_SVRSTATE_SENDPGETDATA;
  CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, (void*) &nParam, sizeof(int));

  // send paramNames
  for (int i = 0; i < nParam; i++) {
    write_name(paramNames[i], svr);
  }

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Param_is
(
 p_server                 svr
)
{
  // wait all ranks have receive msg
  MPI_Barrier(svr->intra_comms[0]);

  // read local code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *code_name = malloc(sizeof(char));
  read_name(&code_name, svr);

  // read param name
  char *param_name = malloc(sizeof(char));
  read_name(&param_name, svr);

  // read initial value
  int data_type = -1;
  CWP_transfer_readdata(svr->connected_socket, svr->max_msg_size, &data_type, sizeof(int));

  // launch
  int is_param = CWP_Param_is(code_name,
                              param_name,
                              data_type);

  // send is_param
  svr->state=CWP_SVRSTATE_SENDPGETDATA;
  CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, (void*) &is_param, sizeof(int));

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Cpl_create
(
 p_server                 svr
)
{
  // wait all ranks have receive msg
  MPI_Barrier(svr->intra_comms[0]);

  // read code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = malloc(sizeof(char));
  read_name(&cpl_id, svr);

  // read coupled code name
  char *coupled_code_name = malloc(sizeof(char));
  read_name(&coupled_code_name, svr);

  // read entities dimension
  int entities_dim;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &entities_dim, sizeof(int));

  // read communication type
  int comm_type;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &comm_type, sizeof(int));

  // read spatial interpolation
  int spatial_interp;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &spatial_interp, sizeof(int));

  // read number of partitions
  int n_part;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &n_part, sizeof(int));

  // read displacement type
  int displacement;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &displacement, sizeof(int));

  // read time exchange type
  int recv_freq_type;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &recv_freq_type, sizeof(int));

  // launch
  CWP_Cpl_create(local_code_name,
                 cpl_id,
                 coupled_code_name,
                 entities_dim,
                 comm_type,
                 spatial_interp,
                 n_part,
                 displacement,
                 recv_freq_type);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Cpl_del
(
 p_server                 svr
)
{
  // wait all ranks have receive msg
  MPI_Barrier(svr->intra_comms[0]);

  // read code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = malloc(sizeof(char));
  read_name(&cpl_id, svr);

  // launch
  CWP_Cpl_del(local_code_name,
              cpl_id);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Properties_dump
(
  p_server                 svr
)
{
  // wait all ranks have receive msg
  MPI_Barrier(svr->intra_comms[0]);

  // launch
  CWP_Properties_dump();

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Visu_set
(
  p_server                 svr
)
{
  // wait all ranks have receive msg
  MPI_Barrier(svr->intra_comms[0]);

  // read code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = malloc(sizeof(char));
  read_name(&cpl_id, svr);

  // read frequency
  int freq;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &freq, sizeof(int));

  // read format
  int format;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &format, sizeof(int));

  // read format option
  char *format_option = malloc(sizeof(char));
  read_name(&format_option, svr);

  // launch
  CWP_Visu_set(local_code_name,
               cpl_id,
               freq,
               format,
               format_option);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_State_update
(
  p_server                 svr
)
{
  // wait all ranks have receive msg
  MPI_Barrier(svr->intra_comms[0]);

  // read code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read state
  int state;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &state, sizeof(int));

  // launch
  CWP_State_update(local_code_name,
                   state);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Time_update
(
  p_server                 svr
)
{
  // wait all ranks have receive msg
  MPI_Barrier(svr->intra_comms[0]);

  // read code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read current time
  double current_time;
  CWP_transfer_readdata(svr->connected_socket, svr->max_msg_size, &current_time, sizeof(double));

  // launch
  CWP_Time_update(local_code_name,
                  current_time);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_State_get
(
  p_server                 svr
)
{
  // wait all ranks have receive msg
  MPI_Barrier(svr->intra_comms[0]);

  // read code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *code_name = malloc(sizeof(char));
  read_name(&code_name, svr);

  // launch
  int state = CWP_State_get(code_name);

  // send state
  svr->state=CWP_SVRSTATE_SENDPGETDATA;
  CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, (void*) &state, sizeof(int));

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Codes_nb_get
(
  p_server                 svr
)
{
  // wait all ranks have receive msg
  MPI_Barrier(svr->intra_comms[0]);

  // launch
  int nb_codes = CWP_Codes_nb_get();

  // send number of codes
  svr->state=CWP_SVRSTATE_SENDPGETDATA;
  CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, (void*) &nb_codes, sizeof(int));

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Codes_list_get
(
  p_server                 svr
)
{
  // wait all ranks have receive msg
  MPI_Barrier(svr->intra_comms[0]);

  // get number of codes
  int nb_codes = CWP_Codes_nb_get();

  // send number of codes
  svr->state=CWP_SVRSTATE_SENDPGETDATA;
  CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, (void*) &nb_codes, sizeof(int));

  // launch
  const char **code_names = CWP_Loc_codes_list_get();
  for (int i = 0; i < nb_codes; i++) {
    write_name(code_names[i], svr);
  }

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Loc_codes_nb_get
(
  p_server                 svr
)
{
  // wait all ranks have receive msg
  MPI_Barrier(svr->intra_comms[0]);

  // launch
  int nb_local_codes = CWP_Loc_codes_nb_get();

  // send number of codes
  svr->state=CWP_SVRSTATE_SENDPGETDATA;
  CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, (void*) &nb_local_codes, sizeof(int));

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Loc_codes_list_get
(
  p_server                 svr
)
{
  // wait all ranks have receive msg
  MPI_Barrier(svr->intra_comms[0]);

  // get number of codes
  int nb_local_codes = CWP_Loc_codes_nb_get();

  // send number of codes
  svr->state=CWP_SVRSTATE_SENDPGETDATA;
  CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, (void*) &nb_local_codes, sizeof(int));

  // launch
  const char **local_code_names = CWP_Loc_codes_list_get();
  for (int i = 0; i < nb_local_codes; i++) {
    write_name(local_code_names[i], svr);
  }

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_N_uncomputed_tgts_get
(
  p_server                 svr
)
{
  // wait all ranks have receive msg
  MPI_Barrier(svr->intra_comms[0]);

  // read local code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = malloc(sizeof(char));
  read_name(&cpl_id, svr);

  // read field identifier
  char *field_id = malloc(sizeof(char));
  read_name(&field_id, svr);

  // read i_part
  int i_part;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &i_part, sizeof(int));

  // launch
  int nb_tgts = CWP_N_uncomputed_tgts_get(local_code_name,
                                          cpl_id,
                                          field_id,
                                          i_part);

  // send number of targets
  svr->state=CWP_SVRSTATE_SENDPGETDATA;
  CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, (void*) &nb_tgts, sizeof(int));

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Uncomputed_tgts_get
(
  p_server                 svr
)
{
  // wait all ranks have receive msg
  MPI_Barrier(svr->intra_comms[0]);

  // read local code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = malloc(sizeof(char));
  read_name(&cpl_id, svr);

  // read field identifier
  char *field_id = malloc(sizeof(char));
  read_name(&field_id, svr);

  // read i_part
  int i_part;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &i_part, sizeof(int));

  // launch
  int nb_tgts = CWP_N_uncomputed_tgts_get(local_code_name,
                                          cpl_id,
                                          field_id,
                                          i_part);

  int *tgts = CWP_Uncomputed_tgts_get(local_code_name,
                                          cpl_id,
                                          field_id,
                                          i_part);

  // send number of targets
  svr->state=CWP_SVRSTATE_SENDPGETDATA;
  CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, (void*) &nb_tgts, sizeof(int));
  CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, (void*) tgts, sizeof(int) * nb_tgts);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_N_computed_tgts_get
(
  p_server                 svr
)
{
  // wait all ranks have receive msg
  MPI_Barrier(svr->intra_comms[0]);

  // read local code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = malloc(sizeof(char));
  read_name(&cpl_id, svr);

  // read field identifier
  char *field_id = malloc(sizeof(char));
  read_name(&field_id, svr);

  // read i_part
  int i_part;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &i_part, sizeof(int));

  // launch
  int nb_tgts = CWP_N_computed_tgts_get(local_code_name,
                                          cpl_id,
                                          field_id,
                                          i_part);

  // send number of targets
  svr->state=CWP_SVRSTATE_SENDPGETDATA;
  CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, (void*) &nb_tgts, sizeof(int));

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Computed_tgts_get
(
  p_server                 svr
)
{
  // wait all ranks have receive msg
  MPI_Barrier(svr->intra_comms[0]);

  // read local code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = malloc(sizeof(char));
  read_name(&cpl_id, svr);

  // read field identifier
  char *field_id = malloc(sizeof(char));
  read_name(&field_id, svr);

  // read i_part
  int i_part;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &i_part, sizeof(int));

  // launch
  int nb_tgts = CWP_N_computed_tgts_get(local_code_name,
                                        cpl_id,
                                        field_id,
                                        i_part);

  int *tgts = CWP_Computed_tgts_get(local_code_name,
                                    cpl_id,
                                    field_id,
                                    i_part);

  // send number of targets
  svr->state=CWP_SVRSTATE_SENDPGETDATA;
  CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, (void*) &nb_tgts, sizeof(int));
  CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, (void*) tgts, sizeof(int) * nb_tgts);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_N_involved_srcs_get
(
  p_server                 svr
)
{
  // wait all ranks have receive msg
  MPI_Barrier(svr->intra_comms[0]);

  // read local code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = malloc(sizeof(char));
  read_name(&cpl_id, svr);

  // read field identifier
  char *field_id = malloc(sizeof(char));
  read_name(&field_id, svr);

  // read i_part
  int i_part;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &i_part, sizeof(int));

  // launch
  int nb_srcs = CWP_N_involved_srcs_get(local_code_name,
                                        cpl_id,
                                        field_id,
                                        i_part);

  // send number of sources
  svr->state=CWP_SVRSTATE_SENDPGETDATA;
  CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, (void*) &nb_srcs, sizeof(int));

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Involved_srcs_get
(
  p_server                 svr
)
{
  // wait all ranks have receive msg
  MPI_Barrier(svr->intra_comms[0]);

  // read local code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = malloc(sizeof(char));
  read_name(&cpl_id, svr);

  // read field identifier
  char *field_id = malloc(sizeof(char));
  read_name(&field_id, svr);

  // read i_part
  int i_part;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &i_part, sizeof(int));

  // launch
  int nb_srcs = CWP_N_involved_srcs_get(local_code_name,
                                        cpl_id,
                                        field_id,
                                        i_part);

  int *srcs = CWP_Involved_srcs_get(local_code_name,
                                    cpl_id,
                                    field_id,
                                    i_part);

  // send number of targets
  svr->state=CWP_SVRSTATE_SENDPGETDATA;
  CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, (void*) &nb_srcs, sizeof(int));
  CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, (void*) srcs, sizeof(int) * nb_srcs);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Spatial_interp_weights_compute
(
  p_server                 svr
)
{
  // wait all ranks have receive msg
  MPI_Barrier(svr->intra_comms[0]);

  // read local code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = malloc(sizeof(char));
  read_name(&cpl_id, svr);

  // launch
  CWP_Spatial_interp_weights_compute(local_code_name,
                                     cpl_id);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Spatial_interp_property_set
(
  p_server                 svr
)
{
  // wait all ranks have receive msg
  MPI_Barrier(svr->intra_comms[0]);

  // read local code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = malloc(sizeof(char));
  read_name(&cpl_id, svr);

  // read property name
  char *property_name = malloc(sizeof(char));
  read_name(&property_name, svr);

  // read property type
  char *property_type = malloc(sizeof(char));
  read_name(&property_type, svr);

  // read property value
  char *property_value = malloc(sizeof(char));
  read_name(&property_value, svr);

  // launch
  CWP_Spatial_interp_property_set(local_code_name,
                                  cpl_id,
                                  property_name,
                                  property_type,
                                  property_value);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_User_tgt_pts_set
(
  p_server                 svr
)
{
  // wait all ranks have receive msg
  MPI_Barrier(svr->intra_comms[0]);

  // read local code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = malloc(sizeof(char));
  read_name(&cpl_id, svr);

  // read i_part
  int i_part;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &i_part, sizeof(int));

  // read n_pts
  int n_pts;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &n_pts, sizeof(int));

  // read coord
  double *coord = malloc(sizeof(double) * n_pts);
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) coord, sizeof(double) * n_pts);

  // read global_num
  CWP_g_num_t *global_num = malloc(sizeof(CWP_g_num_t) * n_pts);
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) global_num, sizeof(CWP_g_num_t) * n_pts);

  // launch
  CWP_User_tgt_pts_set(local_code_name,
                       cpl_id,
                       i_part,
                       n_pts,
                       coord,
                       global_num);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Mesh_interf_finalize
(
  p_server                 svr
)
{
  // wait all ranks have receive msg
  MPI_Barrier(svr->intra_comms[0]);

  // read local code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = malloc(sizeof(char));
  read_name(&cpl_id, svr);

  // launch
  CWP_Mesh_interf_finalize(local_code_name,
                           cpl_id);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Mesh_interf_vtx_set
(
  p_server                 svr
)
{
  // wait all ranks have receive msg
  MPI_Barrier(svr->intra_comms[0]);

  // read local code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = malloc(sizeof(char));
  read_name(&cpl_id, svr);

  // read i_part
  int i_part;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &i_part, sizeof(int));

  // read n_pts
  int n_pts;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &n_pts, sizeof(int));

  // read coord
  double *coord = malloc(sizeof(double) * n_pts);
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) coord, sizeof(double) * n_pts);

  // read global_num
  CWP_g_num_t *global_num = malloc(sizeof(CWP_g_num_t) * n_pts);
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) global_num, sizeof(CWP_g_num_t) * n_pts);

  // launch
  CWP_Mesh_interf_vtx_set(local_code_name,
                          cpl_id,
                          i_part,
                          n_pts,
                          coord,
                          global_num);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Mesh_interf_block_add
(
  p_server                 svr
)
{
  // wait all ranks have receive msg
  MPI_Barrier(svr->intra_comms[0]);

  // read local code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = malloc(sizeof(char));
  read_name(&cpl_id, svr);

  // read block_type
  int block_type;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &block_type, sizeof(int));

  // launch
  int block_id = CWP_Mesh_interf_block_add(local_code_name,
                                           cpl_id,
                                           block_type);
  // send block identifier
  svr->state=CWP_SVRSTATE_SENDPGETDATA;
  CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, (void*) &block_id, sizeof(int));

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Mesh_interf_block_std_set
(
  p_server                 svr
)
{
  // wait all ranks have receive msg
  MPI_Barrier(svr->intra_comms[0]);

  // read local code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = malloc(sizeof(char));
  read_name(&cpl_id, svr);

  // read i_part
  int i_part;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &i_part, sizeof(int));

  // read block_id
  int block_id;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &block_id, sizeof(int));

  // read n_elts
  int n_elts;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &n_elts, sizeof(int));

  // read connectivity
  // TO DO: depends on choices at client side

  // read global number
  CWP_g_num_t *global_num = malloc(sizeof(CWP_g_num_t) * n_elts);
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) global_num, sizeof(CWP_g_num_t) * n_elts);

  // launch
  CWP_Mesh_interf_block_std_set(local_code_name,
                                cpl_id,
                                i_part,
                                block_id,
                                n_elts,
                                connec,
                                global_num);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

/*============================================================================
 * Server function definitions
 *============================================================================*/

/* Create a server */

int
CWP_server_create
(
 MPI_Comm global_comm,
 int server_port,
 int flags,
 p_server svr
)
{
  struct sockaddr_in server_addr;
  socklen_t d;
  int true=1;

  memset(svr,0,sizeof(t_server));
  svr->global_comm      = global_comm;
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
  // TO DO: memset or strncpy to avoid uninitialised byte(s)?
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
    CWP_server_Finalize(svr);

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

  case CWP_MSG_CWP_PARAM_GET:

    // verbose
    if (svr->flags & CWP_SVRFLAG_VERBOSE) {
      log_trace("CWP: server received CWP_Param_get signal\n");
    }

    // launch
    CWP_server_Param_get(svr);

    break;

  case CWP_MSG_CWP_PARAM_SET:

    // verbose
    if (svr->flags & CWP_SVRFLAG_VERBOSE) {
      log_trace("CWP: server received CWP_Param_set signal\n");
    }

    // launch
    CWP_server_Param_set(svr);

    break;

  case CWP_MSG_CWP_PARAM_DEL:

    // verbose
    if (svr->flags & CWP_SVRFLAG_VERBOSE) {
      log_trace("CWP: server received CWP_Param_del signal\n");
    }

    // launch
    CWP_server_Param_del(svr);

    break;

  case CWP_MSG_CWP_PARAM_N_GET:

    // verbose
    if (svr->flags & CWP_SVRFLAG_VERBOSE) {
      log_trace("CWP: server received CWP_Param_n_get signal\n");
    }

    // launch
    CWP_server_Param_n_get(svr);

    break;

  case CWP_MSG_CWP_PARAM_LIST_GET:

    // verbose
    if (svr->flags & CWP_SVRFLAG_VERBOSE) {
      log_trace("CWP: server received CWP_Param_list_get signal\n");
    }

    // launch
    CWP_server_Param_list_get(svr);

    break;

  case CWP_MSG_CWP_PARAM_IS:

    // verbose
    if (svr->flags & CWP_SVRFLAG_VERBOSE) {
      log_trace("CWP: server received CWP_Param_is signal\n");
    }

    // launch
    CWP_server_Param_is(svr);

    break;

  case CWP_MSG_CWP_CPL_CREATE:

    // verbose
    if (svr->flags & CWP_SVRFLAG_VERBOSE) {
      log_trace("CWP: server received CWP_Cpl_create signal\n");
    }

    // launch
    CWP_server_Cpl_create(svr);

    break;

  case CWP_MSG_CWP_CPL_DEL:

    // verbose
    if (svr->flags & CWP_SVRFLAG_VERBOSE) {
      log_trace("CWP: server received CWP_Cpl_del signal\n");
    }

    // launch
    CWP_server_Cpl_del(svr);

    break;

  case CWP_MSG_CWP_PROPERTIES_DUMP:

    // verbose
    if (svr->flags & CWP_SVRFLAG_VERBOSE) {
      log_trace("CWP: server received CWP_Properties_dump signal\n");
    }

    // launch
    CWP_server_Properties_dump(svr);

    break;

  case CWP_MSG_CWP_VISU_SET:

    // verbose
    if (svr->flags & CWP_SVRFLAG_VERBOSE) {
      log_trace("CWP: server received CWP_Visu_set signal\n");
    }

    // launch
    CWP_server_Visu_set(svr);

    break;

  case CWP_MSG_CWP_STATE_UPDATE:

    // verbose
    if (svr->flags & CWP_SVRFLAG_VERBOSE) {
      log_trace("CWP: server received CWP_State_update signal\n");
    }

    // launch
    CWP_server_State_update(svr);

    break;

  case CWP_MSG_CWP_TIME_UPDATE:

    // verbose
    if (svr->flags & CWP_SVRFLAG_VERBOSE) {
      log_trace("CWP: server received CWP_Time_update signal\n");
    }

    // launch
    CWP_server_Time_update(svr);

    break;

  case CWP_MSG_CWP_STATE_GET:

    // verbose
    if (svr->flags & CWP_SVRFLAG_VERBOSE) {
      log_trace("CWP: server received CWP_State_get signal\n");
    }

    // launch
    CWP_server_State_get(svr);

    break;

  case CWP_MSG_CWP_CODES_NB_GET:

    // verbose
    if (svr->flags & CWP_SVRFLAG_VERBOSE) {
      log_trace("CWP: server received CWP_Codes_nb_get signal\n");
    }

    // launch
    CWP_server_Codes_nb_get(svr);

    break;

  case CWP_MSG_CWP_CODES_LIST_GET:

    // verbose
    if (svr->flags & CWP_SVRFLAG_VERBOSE) {
      log_trace("CWP: server received CWP_Codes_list_get signal\n");
    }

    // launch
    CWP_server_Codes_list_get(svr);

    break;

  case CWP_MSG_CWP_LOC_CODES_NB_GET:

    // verbose
    if (svr->flags & CWP_SVRFLAG_VERBOSE) {
      log_trace("CWP: server received CWP_Loc_codes_nb_get signal\n");
    }

    // launch
    CWP_server_Loc_codes_nb_get(svr);

    break;

  case CWP_MSG_CWP_LOC_CODES_LIST_GET:

    // verbose
    if (svr->flags & CWP_SVRFLAG_VERBOSE) {
      log_trace("CWP: server received CWP_Loc_codes_list_get signal\n");
    }

    // launch
    CWP_server_Loc_codes_list_get(svr);

    break;

  case CWP_MSG_CWP_N_UNCOMPUTED_TGTS_GET:

    // verbose
    if (svr->flags & CWP_SVRFLAG_VERBOSE) {
      log_trace("CWP: server received CWP_N_uncomputed_tgts_get signal\n");
    }

    // launch
    CWP_server_N_uncomputed_tgts_get(svr);

    break;

  case CWP_MSG_CWP_UNCOMPUTED_TGTS_GET:

    // verbose
    if (svr->flags & CWP_SVRFLAG_VERBOSE) {
      log_trace("CWP: server received CWP_Uncomputed_tgts_get signal\n");
    }

    // launch
    CWP_server_Uncomputed_tgts_get(svr);

    break;

  case CWP_MSG_CWP_N_COMPUTED_TGTS_GET:

    // verbose
    if (svr->flags & CWP_SVRFLAG_VERBOSE) {
      log_trace("CWP: server received CWP_N_computed_tgts_get signal\n");
    }

    // launch
    CWP_server_N_computed_tgts_get(svr);

    break;

  case CWP_MSG_CWP_COMPUTED_TGTS_GET:

    // verbose
    if (svr->flags & CWP_SVRFLAG_VERBOSE) {
      log_trace("CWP: server received CWP_Computed_tgts_get signal\n");
    }

    // launch
    CWP_server_Computed_tgts_get(svr);

    break;

  case CWP_MSG_CWP_N_INVOLVED_SRCS_GET:

    // verbose
    if (svr->flags & CWP_SVRFLAG_VERBOSE) {
      log_trace("CWP: server received CWP_N_involved_srcs_get\n");
    }

    // launch
    CWP_server_N_involved_srcs_get(svr);

    break;

  case CWP_MSG_CWP_INVOLVED_SRCS_GET:

    // verbose
    if (svr->flags & CWP_SVRFLAG_VERBOSE) {
      log_trace("CWP: server received CWP_Involved_srcs_get signal\n");
    }

    // launch
    CWP_server_Involved_srcs_get(svr);

    break;

  case CWP_MSG_CWP_SPATIAL_INTERP_WEIGHTS_COMPUTE:

    // verbose
    if (svr->flags & CWP_SVRFLAG_VERBOSE) {
      log_trace("CWP: server received CWP_Spatial_interp_weights_compute\n");
    }

    // launch
    CWP_server_Spatial_interp_weights_compute(svr);

    break;

  case CWP_MSG_CWP_SPATIAL_INTERP_PROPERTY_SET:

    // verbose
    if (svr->flags & CWP_SVRFLAG_VERBOSE) {
      log_trace("CWP: server received CWP_User_tgt_pts_set signal\n");
    }

    // launch
    CWP_server_User_tgt_pts_set(svr);

    break;

  case CWP_MSG_CWP_USER_TGT_PTS_SET:

    // verbose
    if (svr->flags & CWP_SVRFLAG_VERBOSE) {
      log_trace("CWP: server received CWP_Spatial_interp_property_set signal\n");
    }

    // launch
    CWP_server_Spatial_interp_property_set(svr);

    break;

  case CWP_MSG_CWP_MESH_INTERF_FINALIZE:

    // verbose
    if (svr->flags & CWP_SVRFLAG_VERBOSE) {
      log_trace("CWP: server received CWP_Mesh_interf_finalize signal\n");
    }

    // launch
    CWP_server_Mesh_interf_finalize(svr);

    break;

  case CWP_MSG_CWP_MESH_INTERF_VTX_SET:

    // verbose
    if (svr->flags & CWP_SVRFLAG_VERBOSE) {
      log_trace("CWP: server received CWP_Mesh_interf_vtx_set signal\n");
    }

    // launch
    CWP_server_Mesh_interf_vtx_set(svr);

    break;

  case CWP_MSG_CWP_MESH_INTERF_BLOCK_ADD:

    // verbose
    if (svr->flags & CWP_SVRFLAG_VERBOSE) {
      log_trace("CWP: server received CWP_Mesh_interf_block_add signal\n");
    }

    // launch
    CWP_server_Mesh_interf_block_add(svr);

    break;

  case CWP_MSG_CWP_MESH_INTERF_BLOCK_STD_SET:

    // verbose
    if (svr->flags & CWP_SVRFLAG_VERBOSE) {
      log_trace("CWP: server received CWP_Mesh_interf_block_std_set signal\n");
    }

    // launch
    CWP_server_Mesh_interf_block_std_set(svr);

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
  int count = 0;

  while (svr->state != CWP_SVRSTATE_TERMINATING) {

    CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,&msg, sizeof(t_message));

    if (CWP_server_msg_handler(svr,&msg) != 0) {
      PDM_error(__FILE__, __LINE__, 0, "Server message handling failed\n");
      return -1;
    }

    count++;

  }

  return 0;
}


#ifdef __cplusplus
}
#endif /* __cplusplus */
