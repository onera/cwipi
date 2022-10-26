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

#include <unistd.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <errno.h>
#include <iostream>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "server.h"
#include "transfer.h"
#include "client.h"
#include "message.h"
#include "struct.hxx"

#include "cwp_priv.h"

#include <pdm_error.h>
#include <pdm_mpi.h>
#include <pdm_logging.h>
#include <pdm_mpi.h>

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/* file struct definition */

static t_server_mpi *svr_mpi;
static t_cwp        *svr_cwp;

/*=============================================================================
 * Private function interfaces
 *============================================================================*/

// --> wrapper

static void read_name(char **name,  p_server svr) {
  int name_size;
  CWP_transfer_readdata(svr->connected_socket, svr->max_msg_size, (void*) &name_size, sizeof(int));
  *name = (char *) realloc(*name, name_size);
  CWP_transfer_readdata(svr->connected_socket, svr->max_msg_size, (void*) *name, name_size);
}

static void write_name(char * name, p_server svr) {
  int name_size = strlen(name) + 1; // +1 for "\0"
  int endian_name_size = name_size;
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
  char         **code_names     = NULL;
  CWP_Status_t  *is_active_rank = NULL;
  double        *time_init      = NULL;

  // receive data
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  CWP_transfer_readdata(svr->connected_socket, svr->max_msg_size, &n_code, sizeof(int));

  // svr_cwp init
  svr_cwp = (p_cwp) malloc(sizeof(t_cwp));
  memset(svr_cwp, 0, sizeof(t_cwp)); // TO DO: find solution to remove and that maps behave fine...

  // mandatory for the map to work
  svr_cwp->char_param_value.clear();
  svr_cwp->coupling.clear();

  if (n_code > 1) {
    PDM_error(__FILE__, __LINE__, 0, "CWIPI client-server not implemented yet for n_code > 1\n");
  }

  code_names = (char **) malloc(sizeof(char *) * n_code);
  for (int i = 0; i < n_code; i++) {
    code_names[i] = (char *) malloc(sizeof(char));
    read_name(&code_names[i], svr);
  }
  is_active_rank = (CWP_Status_t  *) malloc(sizeof(int) * n_code);
  CWP_transfer_readdata(svr->connected_socket, svr->max_msg_size, is_active_rank, n_code * sizeof(CWP_Status_t));
  time_init = (double *) malloc(sizeof(double) * n_code);
  CWP_transfer_readdata(svr->connected_socket, svr->max_msg_size, time_init, n_code * sizeof(double));

  // launch CWP_Init
  svr_mpi->intra_comms = (MPI_Comm *) malloc(sizeof(MPI_Comm) * n_code);
  CWP_Init(svr_mpi->global_comm,
           n_code,
           (const char **) code_names,
           is_active_rank,
           time_init,
           svr_mpi->intra_comms);

  // wait all ranks have receive msg
  MPI_Barrier(svr_mpi->intra_comms[0]);

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
  MPI_Barrier(svr_mpi->intra_comms[0]);

  // free
  if (svr_mpi->intra_comms != NULL) free(svr_mpi->intra_comms);

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
  MPI_Barrier(svr_mpi->intra_comms[0]);

  // launch
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *code_name = (char *) malloc(sizeof(char));
  read_name(&code_name, svr);

  CWP_Param_lock((const char *) code_name);

  // free
  free(code_name);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Param_unlock
(
 p_server                 svr
)
{
  // wait all ranks have receive msg
  MPI_Barrier(svr_mpi->intra_comms[0]);

  // launch
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *code_name = (char *) malloc(sizeof(char));
  read_name(&code_name, svr);

  CWP_Param_unlock((const char *) code_name);

  // free
  free(code_name);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Param_add
(
 p_server                 svr
)
{
  // wait all ranks have receive msg
  MPI_Barrier(svr_mpi->intra_comms[0]);

  // read local code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read param name
  char *param_name = (char *) malloc(sizeof(char));
  read_name(&param_name, svr);

  // read data type
  CWP_Type_t data_type = (CWP_Type_t ) -1;
  CWP_transfer_readdata(svr->connected_socket, svr->max_msg_size, &data_type, sizeof(CWP_Type_t));

  switch (data_type) {

  case CWP_DOUBLE: {
    double double_initial_value;
    CWP_transfer_readdata(svr->connected_socket, svr->max_msg_size, &double_initial_value, sizeof(double));
    CWP_Param_add(local_code_name,
                  param_name,
                  data_type,
                  &double_initial_value);
    } break;

  case CWP_INT: {
    double int_initial_value;
    CWP_transfer_readdata(svr->connected_socket, svr->max_msg_size, &int_initial_value, sizeof(int));
    CWP_Param_add(local_code_name,
                  param_name,
                  data_type,
                  &int_initial_value);
    } break;

  case CWP_CHAR: {
    int name_size;
    CWP_transfer_readdata(svr->connected_socket, svr->max_msg_size, (void*) &name_size, sizeof(int));
    char *char_initial_value = (char *) malloc(name_size);
    CWP_transfer_readdata(svr->connected_socket, svr->max_msg_size, (void*) char_initial_value, name_size);
    std::string s(param_name);
    svr_cwp->char_param_value.insert(std::make_pair(s, char_initial_value));
    CWP_Param_add(local_code_name,
                  param_name,
                  data_type,
                  &svr_cwp->char_param_value[s]);
    } break;

  default:
    PDM_error(__FILE__, __LINE__, 0, "Received unknown CWP_Type_t %i\n", data_type);
  }

  // free
  free(local_code_name);
  free(param_name);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Param_get
(
 p_server                 svr
)
{
  // wait all ranks have receive msg
  MPI_Barrier(svr_mpi->intra_comms[0]);

  // read local code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read param name
  char *param_name = (char *) malloc(sizeof(char));
  read_name(&param_name, svr);

  // read value data_type
  CWP_Type_t data_type = (CWP_Type_t ) -1;
  CWP_transfer_readdata(svr->connected_socket, svr->max_msg_size, &data_type, sizeof(CWP_Type_t));

  // send value
  svr->state=CWP_SVRSTATE_SENDPGETDATA;
  switch (data_type) {

  case CWP_DOUBLE: {
    double double_value;

    CWP_Param_get(local_code_name,
                  param_name,
                  data_type,
                  &double_value);

    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size,(void*) &double_value, sizeof(double));
    } break;

  case CWP_INT: {
    int int_value;
    CWP_Param_get(local_code_name,
                  param_name,
                  data_type,
                  &int_value);

    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size,(void*) &int_value, sizeof(int));
    } break;

  case CWP_CHAR: {
    char *char_value = NULL;
    CWP_Param_get(local_code_name,
                  param_name,
                  data_type,
                  &char_value);
    write_name(char_value, svr);
    free(char_value);
    } break;

  default:
    PDM_error(__FILE__, __LINE__, 0, "Received unknown CWP_Type_t %i\n", data_type);
  }

  // free
  free(local_code_name);
  free(param_name);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Param_set
(
 p_server                 svr
)
{
  // wait all ranks have receive msg
  MPI_Barrier(svr_mpi->intra_comms[0]);

  // read local code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read param name
  char *param_name = (char *) malloc(sizeof(char));
  read_name(&param_name, svr);

  // read initial value
  CWP_Type_t data_type = (CWP_Type_t ) -1;
  CWP_transfer_readdata(svr->connected_socket, svr->max_msg_size, &data_type, sizeof(CWP_Type_t));

  switch (data_type) {

  case CWP_DOUBLE: {
    double double_initial_value;
    CWP_transfer_readdata(svr->connected_socket, svr->max_msg_size, &double_initial_value, sizeof(double));
    CWP_Param_set(local_code_name,
                  param_name,
                  data_type,
                  &double_initial_value);
    } break;

  case CWP_INT: {
    double int_initial_value;
    CWP_transfer_readdata(svr->connected_socket, svr->max_msg_size, &int_initial_value, sizeof(int));
    CWP_Param_set(local_code_name,
                  param_name,
                  data_type,
                  &int_initial_value);
    } break;

  case CWP_CHAR: {
    int name_size;
    CWP_transfer_readdata(svr->connected_socket, svr->max_msg_size, (void*) &name_size, sizeof(int));
    char *char_initial_value = (char *) malloc(name_size);
    CWP_transfer_readdata(svr->connected_socket, svr->max_msg_size, (void*) char_initial_value, name_size);
    std::string s(param_name);
    svr_cwp->char_param_value.insert(std::make_pair(s, char_initial_value));
    CWP_Param_set(local_code_name,
                  param_name,
                  data_type,
                  &svr_cwp->char_param_value[s]);
    } break;

  default:
    PDM_error(__FILE__, __LINE__, 0, "Received unknown CWP_Type_t %i\n", data_type);
  }

  // free
  free(local_code_name);
  free(param_name);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Param_del
(
 p_server                 svr
)
{
  // wait all ranks have receive msg
  MPI_Barrier(svr_mpi->intra_comms[0]);

  // read local code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read param name
  char *param_name = (char *) malloc(sizeof(char));
  read_name(&param_name, svr);

  // read initial value
  CWP_Type_t data_type = (CWP_Type_t ) -1;
  CWP_transfer_readdata(svr->connected_socket, svr->max_msg_size, &data_type, sizeof(CWP_Type_t));

  // launch
  CWP_Param_del(local_code_name,
                param_name,
                data_type);

  // free
  free(local_code_name);
  free(param_name);
  if (data_type == CWP_CHAR) {
    std::string s(param_name);
    if (svr_cwp->char_param_value[s] != NULL) free(svr_cwp->char_param_value[s]);
    svr_cwp->char_param_value.erase(s);
  }

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Param_n_get
(
 p_server                 svr
)
{
  // wait all ranks have receive msg
  MPI_Barrier(svr_mpi->intra_comms[0]);

  // read local code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *code_name = (char *) malloc(sizeof(char));
  read_name(&code_name, svr);

  // read initial value
  CWP_Type_t data_type = (CWP_Type_t ) -1;
  CWP_transfer_readdata(svr->connected_socket, svr->max_msg_size, &data_type, sizeof(CWP_Type_t));

  // launch
  int n_param = CWP_Param_n_get(code_name,
                                data_type);

  // send n_param
  svr->state=CWP_SVRSTATE_SENDPGETDATA;
  CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, (void*) &n_param, sizeof(int));

  // free
  free(code_name);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Param_list_get
(
 p_server                 svr
)
{
  // wait all ranks have receive msg
  MPI_Barrier(svr_mpi->intra_comms[0]);

  // read local code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *code_name = (char *) malloc(sizeof(char));
  read_name(&code_name, svr);

  // read data type
  CWP_Type_t data_type = (CWP_Type_t ) -1;
  CWP_transfer_readdata(svr->connected_socket, svr->max_msg_size, &data_type, sizeof(CWP_Type_t));

  // launch
  svr_cwp->n_param_names = -1;
  svr_cwp->param_names = NULL;
  CWP_Param_list_get(code_name,
                     data_type,
                     &svr_cwp->n_param_names,
                     &svr_cwp->param_names);

  // send nParam
  svr->state=CWP_SVRSTATE_SENDPGETDATA;
  CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, (void*) &svr_cwp->n_param_names, sizeof(int));

  // send paramNames
  for (int i = 0; i < svr_cwp->n_param_names; i++) {
    write_name(svr_cwp->param_names[i], svr);
  }

  // free
  free(code_name);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Param_is
(
 p_server                 svr
)
{
  // wait all ranks have receive msg
  MPI_Barrier(svr_mpi->intra_comms[0]);

  // read local code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *code_name = (char *) malloc(sizeof(char));
  read_name(&code_name, svr);

  // read param name
  char *param_name = (char *) malloc(sizeof(char));
  read_name(&param_name, svr);

  // read initial value
  CWP_Type_t data_type = (CWP_Type_t) -1;
  CWP_transfer_readdata(svr->connected_socket, svr->max_msg_size, &data_type, sizeof(CWP_Type_t));

  // launch
  int is_param = CWP_Param_is(code_name,
                              param_name,
                              data_type);

  // send is_param
  svr->state=CWP_SVRSTATE_SENDPGETDATA;
  CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, (void*) &is_param, sizeof(int));

  // free
  free(code_name);
  free(param_name);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Param_reduce
(
 p_server                 svr
)
{
  // wait all ranks have receive msg
  MPI_Barrier(svr_mpi->intra_comms[0]);

  // read operation
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  CWP_Op_t op = (CWP_Op_t) -1;
  CWP_transfer_readdata(svr->connected_socket, svr->max_msg_size, &op, sizeof(CWP_Op_t));

  // read param name
  char *param_name = (char *) malloc(sizeof(char));
  read_name(&param_name, svr);

  // read data type
  CWP_Type_t data_type = (CWP_Type_t) -1;
  CWP_transfer_readdata(svr->connected_socket, svr->max_msg_size, &data_type, sizeof(CWP_Type_t));

  // read number of codes
  int nCode = -1;
  CWP_transfer_readdata(svr->connected_socket, svr->max_msg_size, &nCode, sizeof(int));

  // read code names
  char **code_names = (char **) malloc(sizeof(char *) * nCode);
  for (int i = 0; i < nCode; i++) {
    code_names[i] = (char *) malloc(sizeof(char));
    read_name(&code_names[i], svr);
  }

  // send res
  svr->state=CWP_SVRSTATE_SENDPGETDATA;
  switch (data_type) {

  case CWP_DOUBLE: {
    double res = -1.;

    switch (nCode) {

      case 0: {
        CWP_Param_reduce(op,
                         param_name,
                         data_type,
                         (void *) &res,
                         nCode);
      } break;

      case 1: {
        CWP_Param_reduce(op,
                         param_name,
                         data_type,
                         (void *) &res,
                         nCode,
                         code_names[0]);
      } break;

      case 2: {
        CWP_Param_reduce(op,
                         param_name,
                         data_type,
                         (void *) &res,
                         nCode,
                         code_names[0],
                         code_names[1]);
      } break;

      case 3: {
        CWP_Param_reduce(op,
                         param_name,
                         data_type,
                         (void *) &res,
                         nCode,
                         code_names[0],
                         code_names[1],
                         code_names[2]);
      } break;

      case 4: {
        CWP_Param_reduce(op,
                         param_name,
                         data_type,
                         (void *) &res,
                         nCode,
                         code_names[0],
                         code_names[1],
                         code_names[2],
                         code_names[3]);
      } break;

      default:
        PDM_error(__FILE__, __LINE__, 0, "CWP_Param_reduce not implemented yet for more than 4 codes\n");
    }

    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, (void*) &res, sizeof(double));
    } break;

  case CWP_INT: {
    int res = -1;

    switch (nCode) {

      case 0: {
        CWP_Param_reduce(op,
                         param_name,
                         data_type,
                         (void *) &res,
                         nCode);
      } break;

      case 1: {
        CWP_Param_reduce(op,
                         param_name,
                         data_type,
                         (void *) &res,
                         nCode,
                         code_names[0]);
      } break;

      case 2: {
        CWP_Param_reduce(op,
                         param_name,
                         data_type,
                         (void *) &res,
                         nCode,
                         code_names[0],
                         code_names[1]);
      } break;

      case 3: {
        CWP_Param_reduce(op,
                         param_name,
                         data_type,
                         (void *) &res,
                         nCode,
                         code_names[0],
                         code_names[1],
                         code_names[2]);
      } break;

      case 4: {
        CWP_Param_reduce(op,
                         param_name,
                         data_type,
                         (void *) &res,
                         nCode,
                         code_names[0],
                         code_names[1],
                         code_names[2],
                         code_names[3]);
      } break;

      default:
        PDM_error(__FILE__, __LINE__, 0, "CWP_Param_reduce not implemented yet for more than 4 codes\n");
    }

    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, (void*) &res, sizeof(int));
    } break;

  default:
    PDM_error(__FILE__, __LINE__, 0, "Unknown CWP_Type_t %i or impossible on CWP_CHAR\n", data_type);
  }

  // free
  free(param_name);
  for (int i = 0; i < nCode; i++) {
    free(code_names[i]);
  }
  free(code_names);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Cpl_create
(
 p_server                 svr
)
{
  // wait all ranks have receive msg
  MPI_Barrier(svr_mpi->intra_comms[0]);

  // read code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = (char *) malloc(sizeof(char));
  read_name(&cpl_id, svr);

  // read coupled code name
  char *coupled_code_name = (char *) malloc(sizeof(char));
  read_name(&coupled_code_name, svr);

  // read entities dimension
  CWP_Interface_t entities_dim;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &entities_dim, sizeof(CWP_Interface_t));

  // read communication type
  CWP_Comm_t comm_type;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &comm_type, sizeof(CWP_Comm_t));

  // read spatial interpolation
  CWP_Spatial_interp_t spatial_interp;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &spatial_interp, sizeof(CWP_Spatial_interp_t));

  // read number of partitions
  int n_part;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &n_part, sizeof(int));

  // read displacement type
  CWP_Dynamic_mesh_t displacement;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &displacement, sizeof(CWP_Dynamic_mesh_t));

  // read time exchange type
  CWP_Time_exch_t recv_freq_type;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &recv_freq_type, sizeof(CWP_Time_exch_t));

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

  // create occurence in map
  std::string s(cpl_id);
  p_coupling coupling = (t_coupling *) malloc(sizeof(t_coupling));
  memset(coupling, 0, sizeof(t_coupling));
  svr_cwp->coupling.insert(std::make_pair(s, coupling));

  // mandatory
  coupling->field.clear();

  // free
  free(local_code_name);
  free(cpl_id);
  free(coupled_code_name);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Cpl_del
(
 p_server                 svr
)
{
  // wait all ranks have receive msg
  MPI_Barrier(svr_mpi->intra_comms[0]);

  // read code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = (char *) malloc(sizeof(char));
  read_name(&cpl_id, svr);

  // launch
  CWP_Cpl_del(local_code_name,
              cpl_id);

  // remove from map
  std::string s(cpl_id);

  p_coupling coupling = svr_cwp->coupling[s];

  if (coupling->vtx_coord != NULL) free(coupling->vtx_coord);
  if (coupling->vtx_global_num != NULL) free(coupling->vtx_global_num);
  if (coupling->connec_faces_idx != NULL) free(coupling->connec_faces_idx);
  if (coupling->connec_faces != NULL) free(coupling->connec_faces);
  if (coupling->connec_cells_idx != NULL) free(coupling->connec_cells_idx);
  if (coupling->connec_cells != NULL) free(coupling->connec_cells);
  if (coupling->cell_global_num != NULL) free(coupling->cell_global_num);
  if (coupling->connec_idx != NULL) free(coupling->connec_idx);
  if (coupling->connec != NULL) free(coupling->connec);
  if (coupling->elt_global_num != NULL) free(coupling->elt_global_num);
  if (coupling->face_edge_idx != NULL) free(coupling->face_edge_idx);
  if (coupling->face_edge != NULL) free(coupling->face_edge);
  if (coupling->edge_vtx_idx != NULL) free(coupling->edge_vtx_idx);
  if (coupling->edge_vtx != NULL) free(coupling->edge_vtx);
  if (coupling->face_global_num != NULL) free(coupling->face_global_num);
  if (coupling->std_connec != NULL) free(coupling->std_connec);
  if (coupling->std_global_num != NULL) free(coupling->std_global_num);
  if (coupling->usr_coord != NULL) free(coupling->usr_coord);
  if (coupling->usr_global_num != NULL) free(coupling->usr_global_num);
  if (coupling->property_name != NULL) free(coupling->property_name);
  if (coupling->property_type != NULL) free(coupling->property_type);
  if (coupling->property_value != NULL) free(coupling->property_value);

  if (!coupling->field.empty()) {
    for (auto const& x : coupling->field) {
      p_field field = coupling->field[x.first];
      if (x.second != NULL) {
        if (field->data != NULL) free(field->data);
      }
      free(x.second);
      coupling->field.erase(x.first);
    }
  }

  free(coupling);
  svr_cwp->coupling.erase(s);

  // free
  free(local_code_name);
  free(cpl_id);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Properties_dump
(
  p_server                 svr
)
{
  // wait all ranks have receive msg
  MPI_Barrier(svr_mpi->intra_comms[0]);

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
  MPI_Barrier(svr_mpi->intra_comms[0]);

  // read code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = (char *) malloc(sizeof(char));
  read_name(&cpl_id, svr);

  // read frequency
  int freq;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &freq, sizeof(int));

  // read format
  CWP_Visu_format_t format;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &format, sizeof(CWP_Visu_format_t));

  // read format option
  char *format_option = (char *) malloc(sizeof(char));
  read_name(&format_option, svr);

  // launch
  CWP_Visu_set(local_code_name,
               cpl_id,
               freq,
               format,
               format_option);

  // free
  free(local_code_name);
  free(cpl_id);
  free(format_option);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_State_update
(
  p_server                 svr
)
{
  // wait all ranks have receive msg
  MPI_Barrier(svr_mpi->intra_comms[0]);

  // read code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read state
  CWP_State_t state;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &state, sizeof(CWP_State_t));

  // launch
  CWP_State_update(local_code_name,
                   state);

  // free
  free(local_code_name);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Time_update
(
  p_server                 svr
)
{
  // wait all ranks have receive msg
  MPI_Barrier(svr_mpi->intra_comms[0]);

  // read code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read current time
  double current_time = - 1.0;
  CWP_transfer_readdata(svr->connected_socket, svr->max_msg_size, &current_time, sizeof(double));

  // launch
  CWP_Time_update(local_code_name,
                  current_time);

  // free
  free(local_code_name);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Output_file_set
(
  p_server                 svr
)
{
  // wait all ranks have receive msg
  MPI_Barrier(svr_mpi->intra_comms[0]);

  // read code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *output_filename = (char *) malloc(sizeof(char));
  read_name(&output_filename, svr);

  // create FILE *
  svr_cwp->output_file = NULL;

  svr_cwp->output_file = fopen(output_filename, "a+"); // TO DO: wich opening mode?

  // launch
  CWP_Output_file_set(svr_cwp->output_file);

  // free
  free(output_filename);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_User_structure_set
(
  p_server                 svr
)
{
  PDM_UNUSED(svr);
  log_trace("CWP: CWP_User_structure_set not implemented in client/server mode\n");
}

void
CWP_server_User_structure_get
(
  p_server                 svr
)
{
  PDM_UNUSED(svr);
  log_trace("CWP: CWP_User_structure_get not implemented in client/server mode\n");
}

void
CWP_server_State_get
(
  p_server                 svr
)
{
  // wait all ranks have receive msg
  MPI_Barrier(svr_mpi->intra_comms[0]);

  // read code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *code_name = (char *) malloc(sizeof(char));
  read_name(&code_name, svr);

  // launch
  int state = CWP_State_get(code_name);

  // send state
  svr->state=CWP_SVRSTATE_SENDPGETDATA;
  CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, (void*) &state, sizeof(int));

  // free
  free(code_name);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Codes_nb_get
(
  p_server                 svr
)
{
  // wait all ranks have receive msg
  MPI_Barrier(svr_mpi->intra_comms[0]);

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
  MPI_Barrier(svr_mpi->intra_comms[0]);

  // get number of codes
  int nb_codes = CWP_Codes_nb_get();

  // send number of codes
  svr->state=CWP_SVRSTATE_SENDPGETDATA;
  CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, (void*) &nb_codes, sizeof(int));

  // launch
  svr_cwp->code_names = CWP_Codes_list_get();
  for (int i = 0; i < nb_codes; i++) {
    write_name((char *) svr_cwp->code_names[i], svr);
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
  MPI_Barrier(svr_mpi->intra_comms[0]);

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
  MPI_Barrier(svr_mpi->intra_comms[0]);

  // get number of codes
  int nb_local_codes = CWP_Loc_codes_nb_get();

  // send number of codes
  svr->state=CWP_SVRSTATE_SENDPGETDATA;
  CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, (void*) &nb_local_codes, sizeof(int));

  // launch
  svr_cwp->loc_code_names = CWP_Loc_codes_list_get();
  for (int i = 0; i < nb_local_codes; i++) {
    write_name((char *) svr_cwp->loc_code_names[i], svr);
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
  MPI_Barrier(svr_mpi->intra_comms[0]);

  // read local code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = (char *) malloc(sizeof(char));
  read_name(&cpl_id, svr);

  // read field identifier
  char *field_id = (char *) malloc(sizeof(char));
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

  // free
  free(local_code_name);
  free(cpl_id);
  free(field_id);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Uncomputed_tgts_get
(
  p_server                 svr
)
{
  // wait all ranks have receive msg
  MPI_Barrier(svr_mpi->intra_comms[0]);

  // read local code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = (char *) malloc(sizeof(char));
  read_name(&cpl_id, svr);

  // read field identifier
  char *field_id = (char *) malloc(sizeof(char));
  read_name(&field_id, svr);

  // read i_part
  int i_part;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &i_part, sizeof(int));

  // launch
  int nb_tgts = CWP_N_uncomputed_tgts_get(local_code_name,
                                          cpl_id,
                                          field_id,
                                          i_part);

  const int *tgts = NULL;
  tgts = CWP_Uncomputed_tgts_get(local_code_name,
                                 cpl_id,
                                 field_id,
                                 i_part);

  // send number of targets
  svr->state=CWP_SVRSTATE_SENDPGETDATA;
  CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, (void*) &nb_tgts, sizeof(int));
  CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, (void*) tgts, sizeof(int) * nb_tgts);

  // free
  free(local_code_name);
  free(cpl_id);
  free(field_id);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_N_computed_tgts_get
(
  p_server                 svr
)
{
  // wait all ranks have receive msg
  MPI_Barrier(svr_mpi->intra_comms[0]);

  // read local code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = (char *) malloc(sizeof(char));
  read_name(&cpl_id, svr);

  // read field identifier
  char *field_id = (char *) malloc(sizeof(char));
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

  // free
  free(local_code_name);
  free(cpl_id);
  free(field_id);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Computed_tgts_get
(
  p_server                 svr
)
{
  // wait all ranks have receive msg
  MPI_Barrier(svr_mpi->intra_comms[0]);

  // read local code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = (char *) malloc(sizeof(char));
  read_name(&cpl_id, svr);

  // read field identifier
  char *field_id = (char *) malloc(sizeof(char));
  read_name(&field_id, svr);

  // read i_part
  int i_part;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &i_part, sizeof(int));

  // launch
  int nb_tgts = CWP_N_computed_tgts_get(local_code_name,
                                        cpl_id,
                                        field_id,
                                        i_part);

  const int *tgts = NULL;
  tgts = CWP_Computed_tgts_get(local_code_name,
                               cpl_id,
                               field_id,
                               i_part);

  // send number of targets
  svr->state=CWP_SVRSTATE_SENDPGETDATA;
  CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, (void*) &nb_tgts, sizeof(int));
  CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, (void*) tgts, sizeof(int) * nb_tgts);

  // free
  free(local_code_name);
  free(cpl_id);
  free(field_id);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_N_involved_srcs_get
(
  p_server                 svr
)
{
  // wait all ranks have receive msg
  MPI_Barrier(svr_mpi->intra_comms[0]);

  // read local code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = (char *) malloc(sizeof(char));
  read_name(&cpl_id, svr);

  // read field identifier
  char *field_id = (char *) malloc(sizeof(char));
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

  // free
  free(local_code_name);
  free(cpl_id);
  free(field_id);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Involved_srcs_get
(
  p_server                 svr
)
{
  // wait all ranks have receive msg
  MPI_Barrier(svr_mpi->intra_comms[0]);

  // read local code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = (char *) malloc(sizeof(char));
  read_name(&cpl_id, svr);

  // read field identifier
  char *field_id = (char *) malloc(sizeof(char));
  read_name(&field_id, svr);

  // read i_part
  int i_part;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &i_part, sizeof(int));

  // launch
  int nb_srcs = CWP_N_involved_srcs_get(local_code_name,
                                        cpl_id,
                                        field_id,
                                        i_part);

  const int *srcs = NULL;
  srcs = CWP_Involved_srcs_get(local_code_name,
                               cpl_id,
                               field_id,
                               i_part);

  // send number of targets
  svr->state=CWP_SVRSTATE_SENDPGETDATA;
  CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, (void*) &nb_srcs, sizeof(int));
  CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, (void*) srcs, sizeof(int) * nb_srcs);

  // free
  free(local_code_name);
  free(cpl_id);
  free(field_id);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Spatial_interp_weights_compute
(
  p_server                 svr
)
{
  // wait all ranks have receive msg
  MPI_Barrier(svr_mpi->intra_comms[0]);

  // read local code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = (char *) malloc(sizeof(char));
  read_name(&cpl_id, svr);

  // launch
  CWP_Spatial_interp_weights_compute(local_code_name,
                                     cpl_id);

  // free
  free(local_code_name);
  free(cpl_id);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Spatial_interp_property_set
(
  p_server                 svr
)
{
  // wait all ranks have receive msg
  MPI_Barrier(svr_mpi->intra_comms[0]);

  // read local code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = (char *) malloc(sizeof(char));
  read_name(&cpl_id, svr);

  // read property name
  std::string s(cpl_id);
  p_coupling coupling = svr_cwp->coupling[s];

  coupling->property_name = (char *) malloc(sizeof(char));
  read_name(&coupling->property_name, svr);

  // read property type
  coupling->property_type = (char *) malloc(sizeof(char));
  read_name(&coupling->property_type, svr);

  // read property value
  coupling->property_value = (char *) malloc(sizeof(char));
  read_name(&coupling->property_value, svr);

  // launch
  CWP_Spatial_interp_property_set(local_code_name,
                                  cpl_id,
                                  coupling->property_name,
                                  coupling->property_type,
                                  coupling->property_value);

  // free
  free(local_code_name);
  free(cpl_id);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_User_tgt_pts_set
(
  p_server                 svr
)
{
  // wait all ranks have receive msg
  MPI_Barrier(svr_mpi->intra_comms[0]);

  // read local code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = (char *) malloc(sizeof(char));
  read_name(&cpl_id, svr);

  // read i_part
  int i_part;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &i_part, sizeof(int));

  // read n_pts
  int n_pts;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &n_pts, sizeof(int));

  // read coord
  std::string s(cpl_id);
  p_coupling coupling = svr_cwp->coupling[s];

  coupling->usr_coord = (double *) malloc(sizeof(double) * 3 * n_pts);
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) coupling->usr_coord, sizeof(double) * 3 * n_pts);

  // read global_num
  int NULL_flag;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &NULL_flag, sizeof(int));

  // launch
  if (NULL_flag) {
    CWP_User_tgt_pts_set(local_code_name,
                         cpl_id,
                         i_part,
                         n_pts,
                         coupling->usr_coord,
                         NULL);
  }
  else {
    coupling->usr_global_num = (CWP_g_num_t *) malloc(sizeof(CWP_g_num_t) * n_pts);
    CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) coupling->usr_global_num, sizeof(CWP_g_num_t) * n_pts);

    CWP_User_tgt_pts_set(local_code_name,
                         cpl_id,
                         i_part,
                         n_pts,
                         coupling->usr_coord,
                         coupling->usr_global_num);
  }

  // free
  free(local_code_name);
  free(cpl_id);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Mesh_interf_finalize
(
  p_server                 svr
)
{
  // wait all ranks have receive msg
  MPI_Barrier(svr_mpi->intra_comms[0]);

  // read local code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = (char *) malloc(sizeof(char));
  read_name(&cpl_id, svr);

  // launch
  CWP_Mesh_interf_finalize(local_code_name,
                           cpl_id);

  // free
  free(local_code_name);
  free(cpl_id);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Mesh_interf_vtx_set
(
  p_server                 svr
)
{
  // wait all ranks have receive msg
  MPI_Barrier(svr_mpi->intra_comms[0]);

  // read local code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = (char *) malloc(sizeof(char));
  read_name(&cpl_id, svr);

  // read i_part
  int i_part;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &i_part, sizeof(int));

  // read n_pts
  int n_pts;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &n_pts, sizeof(int));

  // read coord
  std::string s(cpl_id);
  p_coupling coupling = svr_cwp->coupling[s];
  coupling->vtx_coord = (double *) malloc(sizeof(double) * 3 * n_pts);
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) coupling->vtx_coord, sizeof(double) * 3 * n_pts);

  // read global_num
  int NULL_flag;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &NULL_flag, sizeof(int));

  // launch
  if (NULL_flag) {
    CWP_Mesh_interf_vtx_set(local_code_name,
                            cpl_id,
                            i_part,
                            n_pts,
                            coupling->vtx_coord,
                            NULL);
  }
  else {
    coupling->vtx_global_num = (CWP_g_num_t *) malloc(sizeof(CWP_g_num_t) * n_pts);
    CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) coupling->vtx_global_num, sizeof(CWP_g_num_t) * n_pts);

    CWP_Mesh_interf_vtx_set(local_code_name,
                            cpl_id,
                            i_part,
                            n_pts,
                            coupling->vtx_coord,
                            coupling->vtx_global_num);
  }

  // free
  free(local_code_name);
  free(cpl_id);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Mesh_interf_block_add
(
  p_server                 svr
)
{
  // wait all ranks have receive msg
  MPI_Barrier(svr_mpi->intra_comms[0]);

  // read local code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = (char *) malloc(sizeof(char));
  read_name(&cpl_id, svr);

  // read block_type
  CWP_Block_t block_type;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &block_type, sizeof(CWP_Block_t));

  // launch
  int block_id = CWP_Mesh_interf_block_add(local_code_name,
                                           cpl_id,
                                           block_type);
  // send block identifier
  svr->state=CWP_SVRSTATE_SENDPGETDATA;
  CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, (void*) &block_id, sizeof(int));

  // free
  free(local_code_name);
  free(cpl_id);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Mesh_interf_block_std_set
(
  p_server                 svr
)
{
  // wait all ranks have receive msg
  MPI_Barrier(svr_mpi->intra_comms[0]);

  // read local code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = (char *) malloc(sizeof(char));
  read_name(&cpl_id, svr);

  // read i_part
  int i_part;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &i_part, sizeof(int));

  // read block_id
  int block_id;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &block_id, sizeof(int));;

  // read n_elts
  int n_elts;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &n_elts, sizeof(int));

  // read n_vtx_elt
  int n_vtx_elt;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &n_vtx_elt, sizeof(int));

  // read connectivity
  std::string s(cpl_id);
  p_coupling coupling = svr_cwp->coupling[s];

  coupling->std_connec = (int *) malloc(sizeof(int) * n_elts * n_vtx_elt);
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) coupling->std_connec, sizeof(int) * n_elts * n_vtx_elt);

  // read global number
  int NULL_flag;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &NULL_flag, sizeof(int));

  // launch
  if (NULL_flag) {
    CWP_Mesh_interf_block_std_set(local_code_name,
                                  cpl_id,
                                  i_part,
                                  block_id,
                                  n_elts,
                                  coupling->std_connec,
                                  NULL);
  }
  else {
    coupling->std_global_num = (CWP_g_num_t *) malloc(sizeof(CWP_g_num_t) * n_elts);
    CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) coupling->std_global_num, sizeof(CWP_g_num_t) * n_elts);

    CWP_Mesh_interf_block_std_set(local_code_name,
                                  cpl_id,
                                  i_part,
                                  block_id,
                                  n_elts,
                                  coupling->std_connec,
                                  coupling->std_global_num);
  }

  // free
  free(local_code_name);
  free(cpl_id);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Mesh_interf_f_poly_block_set
(
  p_server                 svr
)
{
  // wait all ranks have receive msg
  MPI_Barrier(svr_mpi->intra_comms[0]);

  // read local code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = (char *) malloc(sizeof(char));
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

  // read connectivity index
  std::string s(cpl_id);
  p_coupling coupling = svr_cwp->coupling[s];

  coupling->connec_idx = (int *) malloc(sizeof(int) * (n_elts+1));
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) coupling->connec_idx, sizeof(int) * (n_elts+1));

  // read connectivity
  coupling->connec = (int *) malloc(sizeof(int) * coupling->connec_idx[n_elts]);
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) coupling->connec, sizeof(int) * coupling->connec_idx[n_elts]);

  // read global number
  int NULL_flag;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &NULL_flag, sizeof(int));

  // launch
  if (NULL_flag) {
    CWP_Mesh_interf_f_poly_block_set(local_code_name,
                                     cpl_id,
                                     i_part,
                                     block_id,
                                     n_elts,
                                     coupling->connec_idx,
                                     coupling->connec,
                                     NULL);
  }
  else {
    coupling->elt_global_num = (CWP_g_num_t *) malloc(sizeof(CWP_g_num_t) * n_elts);
    CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) coupling->elt_global_num, sizeof(CWP_g_num_t) * n_elts);

    CWP_Mesh_interf_f_poly_block_set(local_code_name,
                                     cpl_id,
                                     i_part,
                                     block_id,
                                     n_elts,
                                     coupling->connec_idx,
                                     coupling->connec,
                                     coupling->elt_global_num);
  }

  // free
  free(local_code_name);
  free(cpl_id);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Mesh_interf_f_poly_block_get
(
  p_server                 svr
)
{
  // wait all ranks have receive msg
  MPI_Barrier(svr_mpi->intra_comms[0]);

  // read local code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = (char *) malloc(sizeof(char));
  read_name(&cpl_id, svr);

  // read i_part
  int i_part;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &i_part, sizeof(int));

  // read block_id
  int block_id;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &block_id, sizeof(int));

  // launch
  int          n_elts     = -1;
  int         *connec_idx = NULL;
  int         *connec     = NULL;
  CWP_g_num_t *global_num = NULL;
  CWP_Mesh_interf_f_poly_block_get(local_code_name,
                                   cpl_id,
                                   i_part,
                                   block_id,
                                   &n_elts,
                                   &connec_idx,
                                   &connec,
                                   &global_num);

  // send n_elts
  svr->state=CWP_SVRSTATE_SENDPGETDATA;
  CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, (void*) &n_elts, sizeof(int));

  // send connectivity index
  CWP_transfer_writedata(svr->connected_socket, svr->max_msg_size, connec_idx, sizeof(int) * (n_elts+1));

  // send connectivity
  CWP_transfer_writedata(svr->connected_socket, svr->max_msg_size, connec, sizeof(int) * connec_idx[n_elts]);

  // send global number
  int NULL_flag = 0;
  if (global_num == NULL) {
    NULL_flag = 1;
  }
  CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size,(void*) &NULL_flag, sizeof(int));

  if (!NULL_flag) {
    CWP_transfer_writedata(svr->connected_socket, svr->max_msg_size, global_num, sizeof(CWP_g_num_t) * n_elts);
  }

  // free
  free(local_code_name);
  free(cpl_id);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Mesh_interf_c_poly_block_set
(
  p_server                 svr
)
{
  // wait all ranks have receive msg
  MPI_Barrier(svr_mpi->intra_comms[0]);

  // read local code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = (char *) malloc(sizeof(char));
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

  // read n_faces
  int n_faces;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &n_faces, sizeof(int));

  // read connectivity faces index
  std::string s(cpl_id);
  p_coupling coupling = svr_cwp->coupling[s];

  coupling->connec_faces_idx = (int *) malloc(sizeof(int) * (n_faces+1));
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) coupling->connec_faces_idx, sizeof(int) * (n_faces+1));

  // read connectivity faces
  coupling->connec_faces = (int *) malloc(sizeof(int) * coupling->connec_faces_idx[n_faces]);
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) coupling->connec_faces, sizeof(int) * coupling->connec_faces_idx[n_faces]);

  // read connectivity index
  coupling->connec_cells_idx = (int *) malloc(sizeof(int) * (n_elts+1));
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) coupling->connec_cells_idx, sizeof(int) * (n_elts+1));

  // read connectivity
  coupling->connec_cells = (int *) malloc(sizeof(int) * coupling->connec_cells_idx[n_elts]);
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) coupling->connec_cells, sizeof(int) * coupling->connec_cells_idx[n_elts]);

  // read global number
  int NULL_flag;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &NULL_flag, sizeof(int));

  // launch
  if (NULL_flag) {
    CWP_Mesh_interf_c_poly_block_set(local_code_name,
                                     cpl_id,
                                     i_part,
                                     block_id,
                                     n_elts,
                                     n_faces,
                                     coupling->connec_faces_idx,
                                     coupling->connec_faces,
                                     coupling->connec_cells_idx,
                                     coupling->connec_cells,
                                     NULL);
  }
  else {
    coupling->cell_global_num = (CWP_g_num_t *) malloc(sizeof(CWP_g_num_t) * n_elts);
    CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) coupling->cell_global_num, sizeof(CWP_g_num_t) * n_elts);

    CWP_Mesh_interf_c_poly_block_set(local_code_name,
                                     cpl_id,
                                     i_part,
                                     block_id,
                                     n_elts,
                                     n_faces,
                                     coupling->connec_faces_idx,
                                     coupling->connec_faces,
                                     coupling->connec_cells_idx,
                                     coupling->connec_cells,
                                     coupling->cell_global_num);
  }

  // free
  free(local_code_name);
  free(cpl_id);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Mesh_interf_c_poly_block_get
(
  p_server                 svr
)
{
  // wait all ranks have receive msg
  MPI_Barrier(svr_mpi->intra_comms[0]);

  // read local code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = (char *) malloc(sizeof(char));
  read_name(&cpl_id, svr);

  // read i_part
  int i_part;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &i_part, sizeof(int));

  // read block_id
  int block_id;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &block_id, sizeof(int));

  // launch
  int          n_elts           = -1;
  int          n_faces          = -1;
  int         *connec_faces_idx = NULL;
  int         *connec_faces     = NULL;
  int         *connec_cells_idx = NULL;
  int         *connec_cells     = NULL;
  CWP_g_num_t *global_num       = NULL;
  CWP_Mesh_interf_c_poly_block_get(local_code_name,
                                   cpl_id,
                                   i_part,
                                   block_id,
                                   &n_elts,
                                   &n_faces,
                                   &connec_faces_idx,
                                   &connec_faces,
                                   &connec_cells_idx,
                                   &connec_cells,
                                   &global_num);

  // send n_elts
  svr->state=CWP_SVRSTATE_SENDPGETDATA;
  CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, (void*) &n_elts, sizeof(int));

  // send n_faces
  CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, (void*) &n_faces, sizeof(int));

  // send connectivity index
  CWP_transfer_writedata(svr->connected_socket, svr->max_msg_size, connec_faces_idx, sizeof(int) * (n_faces+1));

  // send connectivity
  CWP_transfer_writedata(svr->connected_socket, svr->max_msg_size, connec_faces, sizeof(int) * connec_faces_idx[n_faces]);

  // send connectivity index
  CWP_transfer_writedata(svr->connected_socket, svr->max_msg_size, connec_cells_idx, sizeof(int) * (n_elts+1));

  // send connectivity
  CWP_transfer_writedata(svr->connected_socket, svr->max_msg_size, connec_cells, sizeof(int) * connec_cells_idx[n_elts]);

  // send global number
  int NULL_flag = 0;
  if (global_num == NULL) {
    NULL_flag = 1;
  }

  CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size,(void*) &NULL_flag, sizeof(int));

  if (!NULL_flag) {
    CWP_transfer_writedata(svr->connected_socket, svr->max_msg_size, global_num, sizeof(CWP_g_num_t) * n_elts);
  }

  // free
  free(local_code_name);
  free(cpl_id);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Mesh_interf_del
(
  p_server                 svr
)
{
  // wait all ranks have receive msg
  MPI_Barrier(svr_mpi->intra_comms[0]);

  // read local code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = (char *) malloc(sizeof(char));
  read_name(&cpl_id, svr);

  // launch
  CWP_Mesh_interf_del(local_code_name,
                      cpl_id);

  // free
  free(local_code_name);
  free(cpl_id);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Mesh_interf_from_cellface_set
(
  p_server                 svr
)
{
  // wait all ranks have receive msg
  MPI_Barrier(svr_mpi->intra_comms[0]);

  // read local code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = (char *) malloc(sizeof(char));
  read_name(&cpl_id, svr);

  // read i_part
  int i_part;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &i_part, sizeof(int));

  // read n_cells
  int n_cells;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &n_cells, sizeof(int));

  // read connectivity cells index
  int *cell_face_idx = (int *) malloc(sizeof(int) * (n_cells+1));
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) cell_face_idx, sizeof(int) * (n_cells+1));

  // read connectivity cells
  int *cell_face = (int *) malloc(sizeof(int) * cell_face_idx[n_cells]);
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) cell_face, sizeof(int) * cell_face_idx[n_cells]);

  // read n_faces
  int n_faces;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &n_faces, sizeof(int));

  // read connectivity faces index
  int *face_vtx_idx = (int *) malloc(sizeof(int) * (n_faces+1));
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) face_vtx_idx, sizeof(int) * (n_faces+1));

  // read connectivity faces
  int *face_vtx = (int *) malloc(sizeof(int) * face_vtx_idx[n_faces]);
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) face_vtx, sizeof(int) * face_vtx_idx[n_faces]);

  // read global number
  int NULL_flag;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &NULL_flag, sizeof(int));

  // launch
  if (NULL_flag) {
    CWP_Mesh_interf_from_cellface_set(local_code_name,
                                      cpl_id,
                                      i_part,
                                      n_cells,
                                      cell_face_idx,
                                      cell_face,
                                      n_faces,
                                      face_vtx_idx,
                                      face_vtx,
                                      NULL);
  }
  else {
    CWP_g_num_t *global_num = (CWP_g_num_t *) malloc(sizeof(CWP_g_num_t) * n_cells);
    CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) global_num, sizeof(CWP_g_num_t) * n_cells);

    CWP_Mesh_interf_from_cellface_set(local_code_name,
                                      cpl_id,
                                      i_part,
                                      n_cells,
                                      cell_face_idx,
                                      cell_face,
                                      n_faces,
                                      face_vtx_idx,
                                      face_vtx,
                                      global_num);
  }

  // free
  free(local_code_name);
  free(cpl_id);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Mesh_interf_from_faceedge_set
(
  p_server                 svr
)
{
  // wait all ranks have receive msg
  MPI_Barrier(svr_mpi->intra_comms[0]);

  // read local code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = (char *) malloc(sizeof(char));
  read_name(&cpl_id, svr);

  // read i_part
  int i_part;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &i_part, sizeof(int));

  // read n_cells
  int n_faces;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &n_faces, sizeof(int));

  // read connectivity faces index
  std::string s(cpl_id);
  p_coupling coupling = svr_cwp->coupling[s];

  coupling->face_edge_idx = (int *) malloc(sizeof(int) * (n_faces+1));
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) coupling->face_edge_idx, sizeof(int) * (n_faces+1));

  // read connectivity faces
  coupling->face_edge = (int *) malloc(sizeof(int) * coupling->face_edge_idx[n_faces]);
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) coupling->face_edge, sizeof(int) * coupling->face_edge_idx[n_faces]);

  // read n_edges
  int n_edges;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &n_edges, sizeof(int));

  // read connectivity edges index
  coupling->edge_vtx_idx = (int *) malloc(sizeof(int) * (n_edges+1));
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) coupling->edge_vtx_idx, sizeof(int) * (n_edges+1));

  // read connectivity edges
  coupling->edge_vtx = (int *) malloc(sizeof(int) * coupling->edge_vtx_idx[n_edges]);
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) coupling->edge_vtx, sizeof(int) * coupling->edge_vtx_idx[n_edges]);

  // read global number
  int NULL_flag;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &NULL_flag, sizeof(int));

  // launch
  if (NULL_flag) {
    CWP_Mesh_interf_from_faceedge_set(local_code_name,
                                      cpl_id,
                                      i_part,
                                      n_faces,
                                      coupling->face_edge_idx,
                                      coupling->face_edge,
                                      n_edges,
                                      coupling->edge_vtx_idx,
                                      coupling->edge_vtx,
                                      NULL);
  }
  else {
    coupling->face_global_num = (CWP_g_num_t *) malloc(sizeof(CWP_g_num_t) * n_faces);
    CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) coupling->face_global_num, sizeof(CWP_g_num_t) * n_faces);

    CWP_Mesh_interf_from_faceedge_set(local_code_name,
                                      cpl_id,
                                      i_part,
                                      n_faces,
                                      coupling->face_edge_idx,
                                      coupling->face_edge,
                                      n_edges,
                                      coupling->edge_vtx_idx,
                                      coupling->edge_vtx,
                                      coupling->face_global_num);
  }

  // free
  free(local_code_name);
  free(cpl_id);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Field_create
(
  p_server                 svr
)
{
  // wait all ranks have receive msg
  MPI_Barrier(svr_mpi->intra_comms[0]);

  // read local code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = (char *) malloc(sizeof(char));
  read_name(&cpl_id, svr);

  // read field identifier
  char *field_id = (char *) malloc(sizeof(char));
  read_name(&field_id, svr);

  // read data_type
  CWP_Type_t data_type;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &data_type, sizeof(CWP_Type_t));

  // read storage
  CWP_Field_storage_t storage;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &storage, sizeof(CWP_Field_storage_t));

  // read n_component
  int n_component;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &n_component, sizeof(int));

  // read target_location
  CWP_Dof_location_t target_location;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &target_location, sizeof(CWP_Dof_location_t));

  // read exch_type
  CWP_Field_exch_t exch_type;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &exch_type, sizeof(CWP_Field_exch_t));

  // read visu_status
  CWP_Status_t visu_status;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &visu_status, sizeof(CWP_Status_t));

  // launch
  CWP_Field_create(local_code_name,
                   cpl_id,
                   field_id,
                   data_type,
                   storage,
                   n_component,
                   target_location,
                   exch_type,
                   visu_status);

  // create occurence in map
  std::string s1(cpl_id);
  p_coupling coupling = svr_cwp->coupling[s1];
  std::string s2(field_id);
  p_field field = (t_field *) malloc(sizeof(t_field));
  coupling->field.insert(std::make_pair(s2, field));

  // free
  free(local_code_name);
  free(cpl_id);
  free(field_id);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Field_data_set
(
  p_server                 svr
)
{
  // wait all ranks have receive msg
  MPI_Barrier(svr_mpi->intra_comms[0]);

  // read local code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = (char *) malloc(sizeof(char));
  read_name(&cpl_id, svr);

  // read field identifier
  char *field_id = (char *) malloc(sizeof(char));
  read_name(&field_id, svr);

  // read i_part
  int i_part;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &i_part, sizeof(int));

  // read map_type
  CWP_Field_map_t map_type;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &map_type, sizeof(CWP_Field_map_t));

  // read data array
  std::string s1(cpl_id);
  p_coupling coupling = svr_cwp->coupling[s1];
  std::string s2(field_id);
  p_field field = coupling->field[s2];

  int size;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &size, sizeof(int));
  field->data = (double *) malloc(sizeof(double) * size);
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) field->data, sizeof(double) * size);

  // launch
  CWP_Field_data_set(local_code_name,
                     cpl_id,
                     field_id,
                     i_part,
                     map_type,
                     field->data);

  // free
  free(local_code_name);
  free(cpl_id);
  free(field_id);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Field_n_component_get
(
  p_server                 svr
)
{
  // wait all ranks have receive msg
  MPI_Barrier(svr_mpi->intra_comms[0]);

  // read local code name
  svr->state = CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = (char *) malloc(sizeof(char));
  read_name(&cpl_id, svr);

  // read field identifier
  char *field_id = (char *) malloc(sizeof(char));
  read_name(&field_id, svr);

  // launch
  int n_component = CWP_Field_n_component_get(local_code_name,
                                              cpl_id,
                                              field_id);

  // send number of components
  svr->state = CWP_SVRSTATE_SENDPGETDATA;
  CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, (void*) &n_component, sizeof(int));

  // free
  free(local_code_name);
  free(cpl_id);
  free(field_id);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Field_target_dof_location_get
(
  p_server                 svr
)
{
  // wait all ranks have receive msg
  MPI_Barrier(svr_mpi->intra_comms[0]);

  // read local code name
  svr->state = CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = (char *) malloc(sizeof(char));
  read_name(&cpl_id, svr);

  // read field identifier
  char *field_id = (char *) malloc(sizeof(char));
  read_name(&field_id, svr);

  // launch
  int dof_location = CWP_Field_target_dof_location_get(local_code_name,
                                                       cpl_id,
                                                       field_id);

  // send number of components
  svr->state = CWP_SVRSTATE_SENDPGETDATA;
  CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, (void*) &dof_location, sizeof(int));

  // free
  free(local_code_name);
  free(cpl_id);
  free(field_id);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Field_storage_get
(
  p_server                 svr
)
{
  // wait all ranks have receive msg
  MPI_Barrier(svr_mpi->intra_comms[0]);

  // read local code name
  svr->state = CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = (char *) malloc(sizeof(char));
  read_name(&cpl_id, svr);

  // read field identifier
  char *field_id = (char *) malloc(sizeof(char));
  read_name(&field_id, svr);

  // launch
  int storage_type = CWP_Field_storage_get(local_code_name,
                                           cpl_id,
                                           field_id);

  // send number of components
  svr->state = CWP_SVRSTATE_SENDPGETDATA;
  CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, (void*) &storage_type, sizeof(int));

  // free
  free(local_code_name);
  free(cpl_id);
  free(field_id);

  svr->state = CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Field_del
(
  p_server                 svr
)
{
  // wait all ranks have receive msg
  MPI_Barrier(svr_mpi->intra_comms[0]);

  // read local code name
  svr->state = CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = (char *) malloc(sizeof(char));
  read_name(&cpl_id, svr);

  // read field identifier
  char *field_id = (char *) malloc(sizeof(char));
  read_name(&field_id, svr);

  // launch
  CWP_Field_del(local_code_name,
                cpl_id,
                field_id);

  // free
  std::string s1(cpl_id);
  p_coupling coupling = svr_cwp->coupling[s1];
  std::string s2(field_id);
  p_field field = coupling->field[s2];

  if (field->data != NULL) free(field->data);

  free(field);
  coupling->field.erase(s2);

  free(local_code_name);
  free(cpl_id);
  free(field_id);

  svr->state = CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Field_issend
(
  p_server                 svr
)
{
  // wait all ranks have receive msg
  MPI_Barrier(svr_mpi->intra_comms[0]);

  // read local code name
  svr->state = CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = (char *) malloc(sizeof(char));
  read_name(&cpl_id, svr);

  // read source field identifier
  char *src_field_id = (char *) malloc(sizeof(char));
  read_name(&src_field_id, svr);

  // launch
  CWP_Field_issend(local_code_name,
                   cpl_id,
                   src_field_id);

  // free
  free(local_code_name);
  free(cpl_id);
  free(src_field_id);

  svr->state = CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Field_irecv
(
  p_server                 svr
)
{
  // wait all ranks have receive msg
  MPI_Barrier(svr_mpi->intra_comms[0]);

  // read local code name
  svr->state = CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = (char *) malloc(sizeof(char));
  read_name(&cpl_id, svr);

  // read target field identifier
  char *tgt_field_id = (char *) malloc(sizeof(char));
  read_name(&tgt_field_id, svr);

  // launch
  CWP_Field_irecv(local_code_name,
                  cpl_id,
                  tgt_field_id);

  // free
  free(local_code_name);
  free(cpl_id);
  free(tgt_field_id);

  svr->state = CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Field_wait_issend
(
  p_server                 svr
)
{
  // wait all ranks have receive msg
  MPI_Barrier(svr_mpi->intra_comms[0]);

  // read local code name
  svr->state = CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = (char *) malloc(sizeof(char));
  read_name(&cpl_id, svr);

  // read source field identifier
  char *src_field_id = (char *) malloc(sizeof(char));
  read_name(&src_field_id, svr);

  // launch
  CWP_Field_wait_issend(local_code_name,
                        cpl_id,
                        src_field_id);

  // free
  free(local_code_name);
  free(cpl_id);
  free(src_field_id);

  svr->state = CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Field_wait_irecv
(
  p_server                 svr
)
{
  // wait all ranks have receive msg
  MPI_Barrier(svr_mpi->intra_comms[0]);

  // read local code name
  svr->state = CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = (char *) malloc(sizeof(char));
  read_name(&cpl_id, svr);

  // read target field identifier
  char *tgt_field_id = (char *) malloc(sizeof(char));
  read_name(&tgt_field_id, svr);

  // read i_part
  int i_part;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &i_part, sizeof(int));

  // read map_type
  CWP_Field_map_t map_type;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &map_type, sizeof(CWP_Field_map_t));

  // read size
  int size;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &size, sizeof(int));

  // launch
  CWP_Field_wait_irecv(local_code_name,
                       cpl_id,
                       tgt_field_id);

  // launch CWP_Field_data_get to retreive field and send to Client
  double *data = NULL;
  CWP_Field_data_get(local_code_name,
                     cpl_id,
                     tgt_field_id,
                     i_part,
                     map_type,
                     &data);

  // send data
  svr->state = CWP_SVRSTATE_SENDPGETDATA;
  CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, (void*) data, sizeof(double) * size);

  // free
  free(local_code_name);
  free(cpl_id);
  free(tgt_field_id);

  svr->state = CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Interp_from_location_unset
(
  p_server                 svr
)
{
  // wait all ranks have receive msg
  MPI_Barrier(svr_mpi->intra_comms[0]);

  // read local code name
  svr->state = CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = (char *) malloc(sizeof(char));
  read_name(&cpl_id, svr);

  // read source field identifier
  char *src_field_id = (char *) malloc(sizeof(char));
  read_name(&src_field_id, svr);

  // launch
  CWP_Interp_from_location_unset(local_code_name,
                                 cpl_id,
                                 src_field_id);

  // free
  free(local_code_name);
  free(cpl_id);
  free(src_field_id);

  svr->state = CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Interp_from_location_set
(
  p_server                 svr
)
{
  PDM_UNUSED(svr);
  // TO DO: create some standard interpolation function on server side for user to choose
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
  struct sockaddr_in *server_addr = (struct sockaddr_in *) malloc(sizeof(struct sockaddr_in));
  socklen_t max_msg_size_len;
  int vrai=1;

  memset(svr,0,sizeof(t_server));
  svr->port             = server_port;
  svr->flags            = flags;
  svr->server_endianess = CWP_transfer_endian_machine();

  svr_mpi               = (p_server_mpi) malloc(sizeof(t_server_mpi));
  svr_mpi->global_comm  = global_comm;

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

  // get maximum message size
  max_msg_size_len  = sizeof(int);
  svr->max_msg_size = 0;
  getsockopt(svr->listen_socket, SOL_SOCKET, SO_RCVBUF, (void *) &svr->max_msg_size, &max_msg_size_len);
  svr->max_msg_size = CWP_MSG_MAXMSGSIZE;

  // verbose
  if (svr->flags & CWP_SVRFLAG_VERBOSE) {
    log_trace("CWP:Max message size:%i\n", svr->max_msg_size);
  }

  // ensure restarted process can reconnect to same socket
  if (setsockopt(svr->listen_socket, SOL_SOCKET, SO_REUSEADDR, &vrai, sizeof(int)) == -1) {
    PDM_error(__FILE__, __LINE__, 0, "Setsockopt failed\n");
    return -1;
  }

  // fill data structure for ip adress + port
  server_addr->sin_family      = AF_INET;
  server_addr->sin_port        = htons(svr->port);
  server_addr->sin_addr.s_addr = INADDR_ANY;
  bzero(&(server_addr->sin_zero),8);

  // bind socket to ip adress + port
  if (bind(svr->listen_socket, (struct sockaddr *) server_addr, sizeof(struct sockaddr)) == -1) {
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

  // free
  free(server_addr);

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

  case CWP_MSG_CWP_PARAM_REDUCE:

    // verbose
    if (svr->flags & CWP_SVRFLAG_VERBOSE) {
      log_trace("CWP: server received CWP_Param_reduce signal\n");
    }

    // launch
    CWP_server_Param_reduce(svr);

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
      log_trace("CWP: server received CWP_Spatial_interp_property_set signal\n");
    }

    // launch
    CWP_server_Spatial_interp_property_set(svr);

    break;

  case CWP_MSG_CWP_USER_TGT_PTS_SET:

    // verbose
    if (svr->flags & CWP_SVRFLAG_VERBOSE) {
      log_trace("CWP: server received CWP_User_tgt_pts_set signal\n");
    }

    // launch
    CWP_server_User_tgt_pts_set(svr);

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


  case CWP_MSG_CWP_MESH_INTERF_F_POLY_BLOCK_SET:

    // verbose
    if (svr->flags & CWP_SVRFLAG_VERBOSE) {
      log_trace("CWP: server received CWP_Mesh_interf_f_poly_block_set signal\n");
    }

    // launch
    CWP_server_Mesh_interf_f_poly_block_set(svr);

    break;

  case CWP_MSG_CWP_MESH_INTERF_F_POLY_BLOCK_GET:

    // verbose
    if (svr->flags & CWP_SVRFLAG_VERBOSE) {
      log_trace("CWP: server received CWP_Mesh_interf_f_poly_block_get signal\n");
    }

    // launch
    CWP_server_Mesh_interf_f_poly_block_get(svr);

    break;

  case CWP_MSG_CWP_MESH_INTERF_C_POLY_BLOCK_SET:

    // verbose
    if (svr->flags & CWP_SVRFLAG_VERBOSE) {
      log_trace("CWP: server received CWP_Mesh_interf_c_poly_block_set signal\n");
    }

    // launch
    CWP_server_Mesh_interf_c_poly_block_set(svr);

    break;

  case CWP_MSG_CWP_MESH_INTERF_C_POLY_BLOCK_GET:

    // verbose
    if (svr->flags & CWP_SVRFLAG_VERBOSE) {
      log_trace("CWP: server received CWP_Mesh_interf_c_poly_block_get signal\n");
    }

    // launch
    CWP_server_Mesh_interf_c_poly_block_get(svr);

    break;

  case CWP_MSG_CWP_MESH_INTERF_DEL:

    // verbose
    if (svr->flags & CWP_SVRFLAG_VERBOSE) {
      log_trace("CWP: server received CWP_Mesh_interf_del signal\n");
    }

    // launch
    CWP_server_Mesh_interf_del(svr);

    break;

  case CWP_MSG_CWP_MESH_INTERF_FROM_CELLFACE_SET:

    // verbose
    if (svr->flags & CWP_SVRFLAG_VERBOSE) {
      log_trace("CWP: server received CWP_Mesh_interf_from_cellface_set signal\n");
    }

    // launch
    CWP_server_Mesh_interf_from_cellface_set(svr);

    break;

  case CWP_MSG_CWP_MESH_INTERF_FROM_FACEEDGE_SET:

    // verbose
    if (svr->flags & CWP_SVRFLAG_VERBOSE) {
      log_trace("CWP: server received CWP_Mesh_interf_from_faceedge_set signal\n");
    }

    // launch
    CWP_server_Mesh_interf_from_faceedge_set(svr);

    break;

  case CWP_MSG_CWP_FIELD_CREATE:

    // verbose
    if (svr->flags & CWP_SVRFLAG_VERBOSE) {
      log_trace("CWP: server received CWP_Field_create signal\n");
    }

    // launch
    CWP_server_Field_create(svr);

    break;

  case CWP_MSG_CWP_FIELD_DATA_SET:

    // verbose
    if (svr->flags & CWP_SVRFLAG_VERBOSE) {
      log_trace("CWP: server received CWP_Field_data_set signal\n");
    }

    // launch
    CWP_server_Field_data_set(svr);

    break;

  case CWP_MSG_CWP_FIELD_N_COMPONENT_GET:

    // verbose
    if (svr->flags & CWP_SVRFLAG_VERBOSE) {
      log_trace("CWP: server received CWP_Field_n_component_get signal\n");
    }

    // launch
    CWP_server_Field_n_component_get(svr);

    break;

  case CWP_MSG_CWP_FIELD_TARGET_DOF_LOCATION_GET:

    // verbose
    if (svr->flags & CWP_SVRFLAG_VERBOSE) {
      log_trace("CWP: server received CWP_Field_target_dof_location_get signal\n");
    }

    // launch
    CWP_server_Field_target_dof_location_get(svr);

    break;

  case CWP_MSG_CWP_FIELD_STORAGE_GET:

    // verbose
    if (svr->flags & CWP_SVRFLAG_VERBOSE) {
      log_trace("CWP: server received CWP_Field_storage_get signal\n");
    }

    // launch
    CWP_server_Field_storage_get(svr);

    break;

  case CWP_MSG_CWP_FIELD_DEL:

    // verbose
    if (svr->flags & CWP_SVRFLAG_VERBOSE) {
      log_trace("CWP: server received CWP_Field_del signal\n");
    }

    // launch
    CWP_server_Field_del(svr);

    break;

  case CWP_MSG_CWP_FIELD_ISSEND:

    // verbose
    if (svr->flags & CWP_SVRFLAG_VERBOSE) {
      log_trace("CWP: server received CWP_CWP_Field_issend signal\n");
    }

    // launch
    CWP_server_Field_issend(svr);

    break;

  case CWP_MSG_CWP_FIELD_IRECV:

    // verbose
    if (svr->flags & CWP_SVRFLAG_VERBOSE) {
      log_trace("CWP: server received CWP_Field_irecv signal\n");
    }

    // launch
    CWP_server_Field_irecv(svr);

    break;

  case CWP_MSG_CWP_FIELD_WAIT_ISSEND:

    // verbose
    if (svr->flags & CWP_SVRFLAG_VERBOSE) {
      log_trace("CWP: server received CWP_Field_wait_issend signal\n");
    }

    // launch
    CWP_server_Field_wait_issend(svr);

    break;

  case CWP_MSG_CWP_FIELD_WAIT_IRECV:

    // verbose
    if (svr->flags & CWP_SVRFLAG_VERBOSE) {
      log_trace("CWP: server received CWP_Field_wait_irecv signal\n");
    }

    // launch
    CWP_server_Field_wait_irecv(svr);

    break;

  case CWP_MSG_CWP_INTERP_FROM_LOCATION_UNSET:

    // verbose
    if (svr->flags & CWP_SVRFLAG_VERBOSE) {
      log_trace("CWP: server received CWP_Interp_from_location_unset signal\n");
    }

    // launch
    CWP_server_Interp_from_location_unset(svr);

    break;

  case CWP_MSG_CWP_INTERP_FROM_LOCATION_SET:

    // verbose
    if (svr->flags & CWP_SVRFLAG_VERBOSE) {
      log_trace("CWP: server received CWP_Interp_from_location_set signal\n");
    }

    // launch
    CWP_server_Interp_from_location_set(svr);

    break;

  case CWP_MSG_CWP_OUTPUT_FILE_SET:

    // verbose
    if (svr->flags & CWP_SVRFLAG_VERBOSE) {
      log_trace("CWP: server received CWP_Output_file_set signal\n");
    }

    // launch
    CWP_server_Output_file_set(svr);

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
  struct sockaddr_in *client_addr = (struct sockaddr_in *) malloc(sizeof(struct sockaddr_in));
  socklen_t client_addr_len = sizeof(struct sockaddr_in);
  int il_sv_endian;
  int il_cl_endian;

  // accept client connexion
  svr->connected_socket = accept(svr->listen_socket, (struct sockaddr *) client_addr, &client_addr_len);

  // verbose
  if (svr->flags & CWP_SVRFLAG_VERBOSE) {
    log_trace("CWP:got a connection from %s:%d\n",
    inet_ntoa(client_addr->sin_addr),ntohs(client_addr->sin_port));
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
    log_trace("Server : client endian %i server endian %i\n",svr->client_endianess,svr->server_endianess);
  }

  t_message msg;

  while (svr->state != CWP_SVRSTATE_TERMINATING) {

    CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,&msg, sizeof(t_message));

    if (CWP_server_msg_handler(svr,&msg) != 0) {
      PDM_error(__FILE__, __LINE__, 0, "Server message handling failed\n");
      return -1;
    }

  }

  // free
  free(client_addr);

  return 0;
}


#ifdef __cplusplus
}
#endif /* __cplusplus */
