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
#include <tuple>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "client.h"
#include "message.h"
#include "transfer.h"

#include "struct.hxx"

#include "cwp.h"
#include "cwp_priv.h"
#include "cwipi_config.h"

#include <pdm_error.h>
#include <pdm_mpi.h>
#include <pdm_io.h>
#include <pdm_version.h>
#include <pdm_printf.h>

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/* file struct definition */

static t_client *clt;
static t_field fields = t_field(); // could be removed since there is clt_cwp, but code needs to be adapted
static FILE* _cwipi_output_listing;
static t_cwp clt_cwp;

/*=============================================================================
 * Macro definitions
 *============================================================================*/

#define CWP_HEADER_SIZE    45

/*=============================================================================
 * Private function interfaces
 *============================================================================*/

// --> taken from cwp.cxx

static int
_cwipi_print_with_c
(
 const char     *const format,
 va_list         arg_ptr
)
{
  return vfprintf(_cwipi_output_listing, format, arg_ptr);
}

// --> n_vtx

static int n_vtx_get(CWP_Block_t block_type) {
  int n_vtx = 0;

  switch(block_type) {

  case CWP_BLOCK_NODE: {
    n_vtx = 1;
    } break;

  case CWP_BLOCK_EDGE2: {
    n_vtx = 2;
    } break;

  case CWP_BLOCK_FACE_TRIA3: {
    n_vtx = 3;
    } break;

  case CWP_BLOCK_FACE_QUAD4: {
    n_vtx = 4;
    } break;

  case CWP_BLOCK_FACE_POLY: {
    // TO DO: what value ?
    PDM_error(__FILE__, __LINE__, 0, "Number of nodes unknown for CWP_BLOCK_FACE_POLY\n");
    return -1;
    } break;

  case CWP_BLOCK_CELL_TETRA4: {
    n_vtx = 4;
    } break;

  case CWP_BLOCK_CELL_HEXA8: {
    n_vtx = 8;
    } break;

  case CWP_BLOCK_CELL_PRISM6: {
    n_vtx = 6;
    } break;

  case CWP_BLOCK_CELL_PYRAM5: {
    n_vtx = 5;
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

  return n_vtx;
}

// --> ho_n_vtx

static int ho_n_vtx_get(CWP_Block_t block_type, int order) {
switch (block_type) {
  case CWP_BLOCK_NODE :
    return 1;
    break;
  case CWP_BLOCK_EDGE2 :
    return 2;
    break;
  case CWP_BLOCK_FACE_TRIA3 :
    return 3;
    break;
  case CWP_BLOCK_FACE_QUAD4 :
    return 4;
    break;
  case CWP_BLOCK_CELL_TETRA4 :
    return 4;
    break;
  case CWP_BLOCK_CELL_PYRAM5 :
    return 5;
    break;
  case CWP_BLOCK_CELL_PRISM6 :
    return 6;
    break;
  case CWP_BLOCK_CELL_HEXA8 :
    return 8;
    break;

  case CWP_BLOCK_EDGEHO :
    return order + 1;
    break;
  case CWP_BLOCK_FACE_TRIAHO :
    return (order + 1) * (order + 2) / 2;
    break;
  case CWP_BLOCK_FACE_QUADHO :
    return (order + 1) * (order + 1);
    break;
  case CWP_BLOCK_CELL_TETRAHO :
    return (order + 1) * (order + 2) * (order + 3) / 6;
    break;
  case CWP_BLOCK_CELL_PYRAMHO :
    return (order + 1) * (order + 2) * (2*order + 3) / 6;
    break;
  case CWP_BLOCK_CELL_PRISMHO :
    return (order + 1) * (order + 1) * (order + 2) / 2;
    break;
  case CWP_BLOCK_CELL_HEXAHO :
    return (order + 1) * (order + 1) * (order + 1);
    break;
  default :
    PDM_error (__FILE__, __LINE__, 0, "Unknown order for Poly2D and Poly3D (block_type %d)\n", block_type);
  }
  return -1;
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

// --> verbose

static void verbose(t_message msg) {

  char *function = (char *) malloc(sizeof(char) * 99);

  switch (msg.message_type) {
  case CWP_MSG_CWP_INIT: {
    char name[] = "CWP_Init";
    strcpy(function, name);
    } break;

  case CWP_MSG_CWP_FINALIZE: {
    char name[] = "CWP_Finalize";
    strcpy(function, name);
    } break;

  case CWP_MSG_CWP_PARAM_LOCK: {
    char name[] = "CWP_Param_lock";
    strcpy(function, name);
    } break;

  case CWP_MSG_CWP_PARAM_UNLOCK: {
    char name[] = "CWP_Param_unlock";
    strcpy(function, name);
    } break;

  case CWP_MSG_CWP_PARAM_ADD: {
    char name[] = "CWP_Param_add";
    strcpy(function, name);
    } break;

  case CWP_MSG_CWP_PARAM_GET: {
    char name[] = "CWP_Param_get";
    strcpy(function, name);
    } break;

  case CWP_MSG_CWP_PARAM_SET: {
    char name[] = "CWP_Param_set";
    strcpy(function, name);
    } break;

  case CWP_MSG_CWP_PARAM_DEL: {
    char name[] = "CWP_Param_del";
    strcpy(function, name);
    } break;

  case CWP_MSG_CWP_PARAM_N_GET: {
    char name[] = "CWP_Param_n_get";
    strcpy(function, name);
    } break;

  case CWP_MSG_CWP_PARAM_LIST_GET: {
    char name[] = "CWP_Param_list_get";
    strcpy(function, name);
    } break;

  case CWP_MSG_CWP_PARAM_IS: {
    char name[] = "CWP_Param_is";
    strcpy(function, name);
    } break;

  case CWP_MSG_CWP_PARAM_REDUCE: {
    char name[] = "CWP_Param_reduce";
    strcpy(function, name);
    } break;

  case CWP_MSG_CWP_CPL_CREATE: {
    char name[] = "CWP_Cpl_create";
    strcpy(function, name);
    } break;

  case CWP_MSG_CWP_CPL_DEL: {
    char name[] = "CWP_Cpl_del";
    strcpy(function, name);
    } break;

  case CWP_MSG_CWP_PROPERTIES_DUMP: {
    char name[] = "CWP_Properties_dump";
    strcpy(function, name);
    } break;

  case CWP_MSG_CWP_VISU_SET: {
    char name[] = "CWP_Visu_set";
    strcpy(function, name);
    } break;

  case CWP_MSG_CWP_STATE_UPDATE: {
    char name[] = "CWP_State_update";
    strcpy(function, name);
    } break;

  case CWP_MSG_CWP_TIME_UPDATE: {
    char name[] = "CWP_Time_update";
    strcpy(function, name);
    } break;

  case CWP_MSG_CWP_STATE_GET: {
    char name[] = "CWP_State_get";
    strcpy(function, name);
    } break;

  case CWP_MSG_CWP_CODES_NB_GET: {
    char name[] = "CWP_Codes_nb_get";
    strcpy(function, name);
    } break;

  case CWP_MSG_CWP_CODES_LIST_GET: {
    char name[] = "CWP_Codes_list_get";
    strcpy(function, name);
    } break;

  case CWP_MSG_CWP_LOC_CODES_NB_GET: {
    char name[] = "CWP_Loc_codes_nb_get";
    strcpy(function, name);
    } break;

  case CWP_MSG_CWP_LOC_CODES_LIST_GET: {
    char name[] = "CWP_Loc_codes_list_get";
    strcpy(function, name);
    } break;

  case CWP_MSG_CWP_N_UNCOMPUTED_TGTS_GET: {
    char name[] = "CWP_N_uncomputed_tgts_get";
    strcpy(function, name);
    } break;

  case CWP_MSG_CWP_UNCOMPUTED_TGTS_GET: {
    char name[] = "CWP_Uncomputed_tgts_get";
    strcpy(function, name);
    } break;

  case CWP_MSG_CWP_N_COMPUTED_TGTS_GET: {
    char name[] = "CWP_N_computed_tgts_get";
    strcpy(function, name);
    } break;

  case CWP_MSG_CWP_COMPUTED_TGTS_GET: {
    char name[] = "CWP_Computed_tgts_get";
    strcpy(function, name);
    } break;

  case CWP_MSG_CWP_N_INVOLVED_SRCS_GET: {
    char name[] = "CWP_N_involved_srcs_get";
    strcpy(function, name);
    } break;

  case CWP_MSG_CWP_INVOLVED_SRCS_GET: {
    char name[] = "CWP_Involved_srcs_get";
    strcpy(function, name);
    } break;

  case CWP_MSG_CWP_SPATIAL_INTERP_WEIGHTS_COMPUTE: {
    char name[] = "CWP_Spatial_interp_weights_compute";
    strcpy(function, name);
    } break;

  case CWP_MSG_CWP_SPATIAL_INTERP_PROPERTY_SET: {
    char name[] = "CWP_Spatial_interp_property_set";
    strcpy(function, name);
    } break;

  case CWP_MSG_CWP_USER_TGT_PTS_SET: {
    char name[] = "CWP_User_tgt_pts_set";
    strcpy(function, name);
    } break;

  case CWP_MSG_CWP_MESH_INTERF_FINALIZE: {
    char name[] = "CWP_Mesh_interf_finalize";
    strcpy(function, name);
    } break;

  case CWP_MSG_CWP_MESH_INTERF_VTX_SET: {
    char name[] = "CWP_Mesh_interf_vtx_set";
    strcpy(function, name);
    } break;

  case CWP_MSG_CWP_MESH_INTERF_BLOCK_ADD: {
    char name[] = "CWP_Mesh_interf_block_add";
    strcpy(function, name);
    } break;

  case CWP_MSG_CWP_MESH_INTERF_BLOCK_STD_SET: {
    char name[] = "CWP_Mesh_interf_block_std_set";
    strcpy(function, name);
    } break;

  case CWP_MSG_CWP_STD_BLOCK_TYPE_GET: {
    char name[] = "CWP_std_block_type_get";
    strcpy(function, name);
    } break;

  case CWP_MSG_CWP_MESH_INTERF_BLOCK_STD_GET: {
    char name[] = "CWP_Mesh_interf_block_std_get";
    strcpy(function, name);
    } break;

  case CWP_MSG_CWP_MESH_INTERF_F_POLY_BLOCK_SET: {
    char name[] = "CWP_Mesh_interf_f_poly_block_set";
    strcpy(function, name);
    } break;

  case CWP_MSG_CWP_MESH_INTERF_F_POLY_BLOCK_GET: {
    char name[] = "CWP_Mesh_interf_f_poly_block_get";
    strcpy(function, name);
    } break;

  case CWP_MSG_CWP_MESH_INTERF_C_POLY_BLOCK_SET: {
    char name[] = "CWP_Mesh_interf_c_poly_block_set";
    strcpy(function, name);
    } break;

  case CWP_MSG_CWP_MESH_INTERF_C_POLY_BLOCK_GET: {
    char name[] = "CWP_Mesh_interf_c_poly_block_get";
    strcpy(function, name);
    } break;

  case CWP_MSG_CWP_MESH_INTERF_DEL: {
    char name[] = "CWP_Mesh_interf_del";
    strcpy(function, name);
    } break;

  case CWP_MSG_CWP_MESH_INTERF_FROM_CELLFACE_SET: {
    char name[] = "CWP_Mesh_interf_from_cellface_set";
    strcpy(function, name);
    } break;

  case CWP_MSG_CWP_MESH_INTERF_FROM_FACEEDGE_SET: {
    char name[] = "CWP_Mesh_interf_from_faceedge_set";
    strcpy(function, name);
    } break;

  case CWP_MSG_CWP_FIELD_CREATE: {
    char name[] = "CWP_Field_create";
    strcpy(function, name);
    } break;

  case CWP_MSG_CWP_FIELD_DATA_SET: {
    char name[] = "CWP_Field_data_set";
    strcpy(function, name);
    } break;

  case CWP_MSG_CWP_FIELD_N_COMPONENT_GET: {
    char name[] = "CWP_Field_n_component_get";
    strcpy(function, name);
    } break;

  case CWP_MSG_CWP_FIELD_TARGET_DOF_LOCATION_GET: {
    char name[] = "CWP_Field_target_dof_location_get";
    strcpy(function, name);
    } break;

  case CWP_MSG_CWP_FIELD_STORAGE_GET: {
    char name[] = "CWP_Field_storage_get";
    strcpy(function, name);
    } break;

  case CWP_MSG_CWP_FIELD_DEL: {
    char name[] = "CWP_Field_del";
    strcpy(function, name);
    } break;

  case CWP_MSG_CWP_FIELD_ISSEND: {
    char name[] = "CWP_Field_issend";
    strcpy(function, name);
    } break;

  case CWP_MSG_CWP_FIELD_IRECV: {
    char name[] = "CWP_Field_irecv";
    strcpy(function, name);
    } break;

  case CWP_MSG_CWP_FIELD_WAIT_ISSEND: {
    char name[] = "CWP_Field_wait_issend";
    strcpy(function, name);
    } break;

  case CWP_MSG_CWP_FIELD_WAIT_IRECV: {
    char name[] = "CWP_Field_wait_irecv";
    strcpy(function, name);
    } break;

  case CWP_MSG_CWP_INTERP_FUNCTION_UNSET: {
    char name[] = "CWP_Interp_function_unset";
    strcpy(function, name);
    } break;

  case CWP_MSG_CWP_INTERP_FUNCTION_SET: {
    char name[] = "CWP_Interp_function_set";
    strcpy(function, name);
    } break;

  case CWP_MSG_CWP_OUTPUT_FILE_SET: {
    char name[] = "CWP_Output_file_set";
    strcpy(function, name);
    } break;

  default:
    PDM_error(__FILE__, __LINE__, 0, "Unknown cwipi function %d\n", msg.message_type);
  }

  char *flag = (char *) malloc(sizeof(char) * 99);

  switch (msg.flag) {
  case CWP_SVR_BEGIN: {
    char name[] = "start of server function";
    strcpy(flag, name);
    } break;

  case CWP_SVR_LCH_BEGIN: {
    char name[] = "received cwipi function arguments";
    strcpy(flag, name);
    } break;

  case CWP_SVR_LCH_END: {
    char name[] = "cwipi function execution finished";
    strcpy(flag, name);
    } break;

  default:
    PDM_error(__FILE__, __LINE__, 0, "Unknown server flag %d\n", msg.flag);
  }

  // print
  PDM_printf("%s-CWP-SERVER: %s --> %s\n", clt->code_name, function, flag);
  PDM_printf_flush();

  // free
  free(function);
  free(flag);
}

/*============================================================================
 * Client function definitions
 *============================================================================*/

/* Connect to a server */

int
CWP_client_connect
(
 MPI_Comm  comm,
 const char* server_name,
 int server_port,
 int flags
)
{
  // struct hostent *host;
  // struct sockaddr_in *server_addr = (struct sockaddr_in *) malloc(sizeof(struct sockaddr_in));
  socklen_t max_msg_size_len = 0;
  int il_cl_endian;

  clt = (t_client *) malloc(sizeof(t_client));
  memset(clt,0,sizeof(t_client));
  strncpy(clt->server_name, server_name,sizeof(clt->server_name));
  clt->server_port = server_port;
  clt->flags = flags;
  clt->comm = comm;
  MPI_Comm_rank(clt->comm, &clt->i_rank);

  struct addrinfo *svr_info;
  // struct addrinfo hints = {};
  // hints.ai_family = AF_INET;
  // hints.ai_socktype = SOCK_STREAM;
  // hints.ai_protocol = IPPROTO_TCP;

  char port_str[16] = {};
  sprintf(port_str, "%d", server_port);

  getaddrinfo(server_name, port_str, NULL, &svr_info); // hint
  char *dst = (char *) malloc(sizeof(char) * INET_ADDRSTRLEN);
  inet_ntop(AF_INET, svr_info->ai_addr->sa_data, dst, INET_ADDRSTRLEN);

  // host = (struct hostent *) gethostbyname(server_name);
  // printf("host: %p\n", host);

  // create socket
  if ((clt->socket = socket(svr_info->ai_family, SOCK_STREAM, 0)) == -1) {
    PDM_error(__FILE__, __LINE__, 0, "Could not create Socket\n");
    return -1;
  }

  // connect
  while (connect(clt->socket, (struct sockaddr *) svr_info->ai_addr,
                 (int) svr_info->ai_addrlen) == -1) {}

  // verbose
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    printf("CWP-CLIENT: client connected to server on %s port %i\n", server_name, server_port);
    fflush(stdout);
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

  // send verbose
  int endian_flags = flags;
  CWP_swap_endian_4bytes(&endian_flags, 1);
  CWP_transfer_writedata(clt->socket, clt->max_msg_size, (void*) &endian_flags, sizeof(int));

  // free
  free(svr_info);

  return 0;

}

/* Disconnect */

int
CWP_client_disconnect
()
{
  t_message msg;

  MPI_Barrier(clt->comm);
  if ((clt->flags  & CWP_FLAG_VERBOSE) && (clt->i_rank == 0)) {
    printf("CWP-CLIENT: Client shutting down\n");
    fflush(stdout);
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

/*=============================================================================
 * Client CWIPI function interfaces
 *============================================================================*/

void
CWP_client_Init
(
  MPI_Comm                  comm,
  char                    *config,
  const int                n_code,
  const char             **code_names,
  const CWP_Status_t      *is_active_rank,
  const double            *time_init
)
{
  // one code par rank in client-server mode
  if (n_code != 1) {
    PDM_error(__FILE__, __LINE__, 0, "Expected 1 code per rank in client-server mode, got %d\n", n_code);
  }

  /*intra communicators */

  int my_rank;
  int total_rank;
  MPI_Comm_rank(comm, &my_rank);
  MPI_Comm_size(comm, &total_rank);

  // --> get local code name size
  int i_rank_size = strlen(code_names[0]) + 1;
  int *j_rank_size = NULL;

  if (my_rank == 0) {
    j_rank_size = (int *) malloc(sizeof(int) * total_rank);
  }

  MPI_Barrier(comm);

  MPI_Gather((const void *) &i_rank_size, 1, MPI_INT,
             (void *)       j_rank_size, 1, MPI_INT,
             0, comm);

  int total_size = 0;
  int *j_rank_idx = NULL;

  if (my_rank == 0) {
    j_rank_idx = (int *) malloc(sizeof(int) * total_rank);
    j_rank_idx[0] = 0;
    for (int i = 0; i < total_rank; i++) {
      total_size += j_rank_size[i];
      if (i > 0) {
        j_rank_idx[i] = j_rank_idx[i-1] + j_rank_size[i-1];
      }
    }
  }

  // --> get local code names
  char *j_rank_code_names = NULL;

  if (my_rank == 0) {
    j_rank_code_names = (char *) malloc(sizeof(char) * total_size);
  }

  MPI_Barrier(comm);

  MPI_Gatherv((const void *) code_names[0], i_rank_size, MPI_CHAR,
             (void *)       j_rank_code_names, j_rank_size, j_rank_idx, MPI_CHAR,
             0, comm);

  // --> create map
  char *key = NULL;
  int *value = NULL;
  int n_codes = 0;
  int index = 0;
  int total_n_codes_size = 0;

  if (my_rank == 0) {
    n_codes = 0;
    std::map<std::string, int> color;
    for (int i = 0; i < total_rank; i++) {
      std::string s(j_rank_code_names + index);
      if (color.find(s) == color.end()) {
        color[s] = n_codes;
        n_codes++;
      }
      index += j_rank_size[i];
    }

    key = (char *) malloc(sizeof(char) * total_size);
    value = (int *) malloc(sizeof(int) * n_codes);

    int iter = 0;
    for (auto const& x : color) {
      int size = x.first.length() + 1;
      memcpy(key + total_n_codes_size, x.first.c_str(), size);
      value[iter] = x.second;
      iter++;
      total_n_codes_size += size;
    }

    key = (char *) realloc((void *) key, total_n_codes_size);
  }

  MPI_Barrier(comm);

  // --> send code names vect associated color vect
  MPI_Bcast((void *) &n_codes, 1, MPI_INT, 0, comm);
  MPI_Bcast((void *) &total_n_codes_size, 1, MPI_INT, 0, comm);

  if (my_rank != 0) {
    key = (char *) malloc(sizeof(char) * total_n_codes_size);
    value = (int *) malloc(sizeof(int) * n_codes);
  }

  MPI_Barrier(comm);

  MPI_Bcast((void *) key, total_n_codes_size, MPI_CHAR, 0, comm);
  MPI_Bcast((void *) value, n_codes, MPI_INT, 0, comm);

  MPI_Barrier(comm);

  /* connect */

  // mpi
  int i_rank;
  int n_rank;

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&comm);
  PDM_MPI_Comm_rank(pdm_comm, &i_rank);
  PDM_MPI_Comm_size(pdm_comm, &n_rank);

  // read host_name and port from config file
  // --> open
  PDM_io_file_t *read = NULL;
  PDM_l_num_t    ierr;

  PDM_io_open(config,
              PDM_IO_FMT_BIN,
              PDM_IO_SUFF_MAN,
              "",
              PDM_IO_BACKUP_OFF,
              PDM_IO_KIND_MPI_SIMPLE,
              PDM_IO_MOD_READ,
              PDM_IO_NATIVE,
              pdm_comm,
              -1.,
              &read,
              &ierr);

  // --> global read offset in header
  char *buffer = (char *) malloc(CWP_HEADER_SIZE+1);
  for (int i = 0; i < CWP_HEADER_SIZE; i++) {
    buffer[i] = '\0';
  }

  PDM_io_global_read(read,
                     CWP_HEADER_SIZE * sizeof(char),
                     1,
                     buffer);

  char div_line[] = "\n";
  char *line = strtok(buffer, div_line);
  line = strtok(NULL, div_line);
  char *second_line = (char *) malloc(sizeof(char) * (strlen(line) + 1));
  memcpy(second_line, line, strlen(line) + 1);
  line = strtok(NULL, div_line);

  char div_word[] = " ";
  char *N = strtok(second_line, div_word);
  N = strtok(NULL, div_word);

  int nb_rank = atoi(N);

  char *word = strtok(line, div_word);
  word = strtok(NULL, div_word);

  int offset = atoi(word);

  // check number of ranks
  if (n_rank != nb_rank) {
    PDM_error(__FILE__, __LINE__, 0, "Client executing on %d and server on %d ranks. Should be on the same number of ranks\n", n_rank, nb_rank);
  }

  // --> block read hostname/port
  char *data = (char *) malloc(offset+1);

  for (int i = 0; i < offset+1; i++) {
    data[i] = '\0';
  }

  int one = 1;
  PDM_g_num_t i_rank_gnum = (PDM_g_num_t) (i_rank+1);

  PDM_io_par_interlaced_read(read,
                             PDM_STRIDE_CST_INTERLACED,
             (PDM_l_num_t *) &one,
               (PDM_l_num_t) offset,
               (PDM_l_num_t) one,
                             &i_rank_gnum,
                             data);

  char div[] = "/";
  char *str = strtok(data, div);
  char *server_name = (char *) malloc(strlen(str)+1);
  memcpy(server_name, str, strlen(str)+1);
  str = strtok(NULL, div);
  int server_port = atoi(str);

  // --> close
  PDM_io_close(read);
  PDM_io_free(read);

  // get verbose flag
  char* pPath;
  int flags = 0;
  pPath = getenv ("CWP_TCP_IP_VERBOSE");
  if (pPath!=NULL) {
    flags = atoi(pPath);
  }

  // connect
  if (CWP_client_connect(comm, server_name, server_port, flags) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "Client connexion failed\n");
  }

  // code name
  clt->code_name = (char *) malloc(sizeof(char) * (strlen(code_names[0])+1));
  memset(clt->code_name, 0, strlen(code_names[0])+1);
  memcpy(clt->code_name, code_names[0], strlen(code_names[0]));

  // free
  free(buffer);
  free(data);
  free(server_name);
  free(second_line);

  /* cwipi init */
  t_message msg;

  // verbose
  MPI_Barrier(clt->comm);
  if ((clt->flags  & CWP_FLAG_VERBOSE) && (clt->i_rank == 0)) {
    PDM_printf("%s-CWP-CLIENT: Client initiating CWP_Init\n", clt->code_name);
    PDM_printf_flush();
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_INIT);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Init failed to send message header\n");
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
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

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // free
  free(endian_is_active_rank);
  free(endian_time_init);

  // standard cwipi init message
  if (clt->i_rank == 0) {

    printf("\ncwipi " CWIPI_VERSION " initializing\n");
    printf("------------------------\n\n");

  }
}

void
CWP_client_Finalize()
{
  /* finalize */
  t_message msg;

  // verbose
  MPI_Barrier(clt->comm);
  if ((clt->flags  & CWP_FLAG_VERBOSE) && (clt->i_rank == 0)) {
    PDM_printf("%s-CWP-CLIENT: Client initiating CWP_Finalize\n", clt->code_name);
    PDM_printf_flush();
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_FINALIZE);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Finalize failed to send message header\n");
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  if (clt->code_name != NULL) free(clt->code_name);

  if (clt_cwp.code_names != NULL) {
    for (int i = 0; i < clt_cwp.n_code_names; i++) {
      if (clt_cwp.code_names[i] != NULL) free((void *) clt_cwp.code_names[i]);
    }
    free(clt_cwp.code_names);
  }


  if (clt_cwp.loc_code_names != NULL) {
    for (int i = 0; i < clt_cwp.n_loc_code_names; i++) {
      if (clt_cwp.loc_code_names[i] != NULL) free((void *) clt_cwp.loc_code_names[i]);
    }
    free(clt_cwp.loc_code_names);
  }

  if (clt_cwp.param_names != NULL) {
    for (int i = 0; i < clt_cwp.n_param_names; i++) {
      if (clt_cwp.param_names[i] != NULL) free(clt_cwp.param_names[i]);
    }
    free(clt_cwp.param_names);
  }

  if (!clt_cwp.char_param_value.empty()) {
    for (const auto& x : clt_cwp.char_param_value) {
      if (x.second != NULL) {
        free((void *) x.second);
      }
    }
    clt_cwp.char_param_value.clear();
  }

  /* disconnect */
  CWP_client_disconnect();
}

void
CWP_client_Param_lock
(
const char *code_name
)
{
  t_message msg;

  // verbose
  MPI_Barrier(clt->comm);
  if ((clt->flags  & CWP_FLAG_VERBOSE) && (clt->i_rank == 0)) {
    PDM_printf("%s-CWP-CLIENT: Client initiating CWP_Param_lock\n", clt->code_name);
    PDM_printf_flush();
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_PARAM_LOCK);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Param_lock failed to send message header\n");
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // send code name
  write_name(code_name);

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }
}

void
CWP_client_Param_unlock
(
const char *code_name
)
{
  t_message msg;

  // verbose
  MPI_Barrier(clt->comm);
  if ((clt->flags  & CWP_FLAG_VERBOSE) && (clt->i_rank == 0)) {
    PDM_printf("%s-CWP-CLIENT: Client initiating CWP_Param_unlock\n", clt->code_name);
    PDM_printf_flush();
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_PARAM_UNLOCK);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Param_unlock failed to send message header\n");
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // send code name
  write_name(code_name);

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }
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
  MPI_Barrier(clt->comm);
  if ((clt->flags  & CWP_FLAG_VERBOSE) && (clt->i_rank == 0)) {
    PDM_printf("%s-CWP-CLIENT: Client initiating CWP_Param_add\n", clt->code_name);
    PDM_printf_flush();
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_PARAM_ADD);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Param_add failed to send message header\n");
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
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

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
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
  MPI_Barrier(clt->comm);
  if ((clt->flags  & CWP_FLAG_VERBOSE) && (clt->i_rank == 0)) {
    PDM_printf("%s-CWP-CLIENT: Client initiating CWP_Param_get\n", clt->code_name);
    PDM_printf_flush();
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_PARAM_GET);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Param_get failed to send message header\n");
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // send code name
  write_name(code_name);

  // send param name
  write_name(param_name);

  // send value data type
  CWP_Type_t endian_data_type = data_type;
  CWP_swap_endian_4bytes((int *) &endian_data_type, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_data_type, sizeof(CWP_Type_t));

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // receive value
  switch (data_type) {

  case CWP_DOUBLE: {
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, value, sizeof(double));
    } break;

  case CWP_INT: {
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, value, sizeof(int));
    } break;

  case CWP_CHAR: {
    char *char_value = (char *) malloc(sizeof(char));
    read_name(&char_value);
    std::string s(param_name);
    clt_cwp.char_param_value.insert(std::make_pair(s, char_value));
    * (char **) value = (char *) clt_cwp.char_param_value[s];

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
  MPI_Barrier(clt->comm);
  if ((clt->flags  & CWP_FLAG_VERBOSE) && (clt->i_rank == 0)) {
    PDM_printf("%s-CWP-CLIENT: Client initiating CWP_Param_set\n", clt->code_name);
    PDM_printf_flush();
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_PARAM_SET );

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Param_set failed to send message header\n");
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
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

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
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
  MPI_Barrier(clt->comm);
  if ((clt->flags  & CWP_FLAG_VERBOSE) && (clt->i_rank == 0)) {
    PDM_printf("%s-CWP-CLIENT: Client initiating CWP_Param_del\n", clt->code_name);
    PDM_printf_flush();
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_PARAM_DEL);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Param_del failed to send message header\n");
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // send local code name
  write_name(local_code_name);

  // send param name
  write_name(param_name);

  // send data type
  CWP_Type_t endian_data_type = data_type;
  CWP_swap_endian_4bytes((int *) &endian_data_type, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_data_type, sizeof(CWP_Type_t));

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }
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
  MPI_Barrier(clt->comm);
  if ((clt->flags  & CWP_FLAG_VERBOSE) && (clt->i_rank == 0)) {
    PDM_printf("%s-CWP-CLIENT: Client initiating CWP_Param_n_get\n", clt->code_name);
    PDM_printf_flush();
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_PARAM_N_GET);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Param_n_get failed to send message header\n");
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // send local code name
  write_name(code_name);

  // send data type
  CWP_Type_t endian_data_type = data_type;
  CWP_swap_endian_4bytes((int *) &endian_data_type, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size, (void*) &endian_data_type, sizeof(CWP_Type_t));

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

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
  MPI_Barrier(clt->comm);
  if ((clt->flags  & CWP_FLAG_VERBOSE) && (clt->i_rank == 0)) {
    PDM_printf("%s-CWP-CLIENT: Client initiating CWP_Param_list_get\n", clt->code_name);
    PDM_printf_flush();
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_PARAM_LIST_GET);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Param_list_get failed to send message header\n");
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // send local code name
  write_name(code_name);

  // send data type
  CWP_Type_t endian_data_type = data_type;
  CWP_swap_endian_4bytes((int *) &endian_data_type, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size, (void*) &endian_data_type, sizeof(CWP_Type_t));

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // read n_param
  CWP_transfer_readdata(clt->socket, clt->max_msg_size, (void*) &(clt_cwp.n_param_names), sizeof(int));
  *nParam = clt_cwp.n_param_names;

  // read param names
  clt_cwp.param_names = (char **) malloc(sizeof(char *) * (clt_cwp.n_param_names));
  for (int i = 0; i < clt_cwp.n_param_names; i++) {
    (clt_cwp.param_names)[i] = (char *) malloc(sizeof(char));
    read_name(&((clt_cwp.param_names)[i]));
  }
  *paramNames = clt_cwp.param_names;
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
  MPI_Barrier(clt->comm);
  if ((clt->flags  & CWP_FLAG_VERBOSE) && (clt->i_rank == 0)) {
    PDM_printf("%s-CWP-CLIENT: Client initiating CWP_Param_is\n", clt->code_name);
    PDM_printf_flush();
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_PARAM_IS);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Param_is failed to send message header\n");
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // send local code name
  write_name(code_name);

  // send param name
  write_name(param_name);

  // send data type
  CWP_Type_t endian_data_type = data_type;
  CWP_swap_endian_4bytes((int *) &endian_data_type, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_data_type, sizeof(CWP_Type_t));

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

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
 const char      **code_names
)
{
  t_message msg;

  // verbose
  MPI_Barrier(clt->comm);
  if ((clt->flags  & CWP_FLAG_VERBOSE) && (clt->i_rank == 0)) {
    PDM_printf("%s-CWP-CLIENT: Client initiating CWP_Param_reduce\n", clt->code_name);
    PDM_printf_flush();
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_PARAM_REDUCE);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Param_reduce failed to send message header\n");
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
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

  // send code names
  for (int i = 0; i < nCode; i++) {
    write_name(code_names[i]);
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // retreive result
  switch (data_type) {

  case CWP_DOUBLE: {
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, (void*) res, sizeof(double));
    } break;

  case CWP_INT: {
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, (void*) res, sizeof(int));
    } break;

  default:
    PDM_error(__FILE__, __LINE__, 0, "Unknown CWP_Type_t %i or impossible on CWP_CHAR\n", data_type);
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
  MPI_Barrier(clt->comm);
  if ((clt->flags  & CWP_FLAG_VERBOSE) && (clt->i_rank == 0)) {
    PDM_printf("%s-CWP-CLIENT: Client initiating CWP_Cpl_create\n", clt->code_name);
    PDM_printf_flush();
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_CPL_CREATE);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Cpl_create failed to send message header\n");
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
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

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // create occurence in map
  std::string s(cpl_id);
  t_coupling coupling = t_coupling();
  clt_cwp.coupling.insert(std::make_pair(s, coupling));
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
  MPI_Barrier(clt->comm);
  if ((clt->flags  & CWP_FLAG_VERBOSE) && (clt->i_rank == 0)) {
    PDM_printf("%s-CWP-CLIENT: Client initiating CWP_Cpl_del\n", clt->code_name);
    PDM_printf_flush();
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_CPL_DEL);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Cpl_del failed to send message header\n");
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // send code name
  write_name(local_code_name);

  // send coupling identifier
  write_name(cpl_id);

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // free struct
  std::string s(cpl_id);

  if (clt_cwp.coupling[s].connec_faces_idx != NULL) free(clt_cwp.coupling[s].connec_faces_idx);
  if (clt_cwp.coupling[s].connec_faces != NULL    ) free(clt_cwp.coupling[s].connec_faces);
  if (clt_cwp.coupling[s].connec_cells_idx != NULL) free(clt_cwp.coupling[s].connec_cells_idx);
  if (clt_cwp.coupling[s].connec_cells != NULL    ) free(clt_cwp.coupling[s].connec_cells);
  if (clt_cwp.coupling[s].cell_global_num != NULL ) free(clt_cwp.coupling[s].cell_global_num);

  if (clt_cwp.coupling[s].connec_idx != NULL     ) free(clt_cwp.coupling[s].connec_idx);
  if (clt_cwp.coupling[s].connec != NULL         ) free(clt_cwp.coupling[s].connec);
  if (clt_cwp.coupling[s].elt_global_num != NULL ) free(clt_cwp.coupling[s].elt_global_num);

  if (clt_cwp.coupling[s].std_connec != NULL    ) free(clt_cwp.coupling[s].std_connec);
  if (clt_cwp.coupling[s].std_global_num != NULL) free(clt_cwp.coupling[s].std_global_num);

  if (clt_cwp.coupling[s].ho_std_connec != NULL    ) free(clt_cwp.coupling[s].ho_std_connec);
  if (clt_cwp.coupling[s].ho_std_global_num != NULL) free(clt_cwp.coupling[s].ho_std_global_num);

  if (!clt_cwp.coupling[s].field.empty()) {

    std::map<std::string, t_field>::iterator itr = clt_cwp.coupling[s].field.begin();
    while (itr != clt_cwp.coupling[s].field.end()) {
      if ((itr->second).srcs   != NULL) {
        free((itr->second).srcs);
        (itr->second).srcs = NULL;
      }
      if ((itr->second).c_tgts != NULL) {
        free((itr->second).c_tgts);
        (itr->second).c_tgts = NULL;
      }
      if ((itr->second).u_tgts != NULL) {
        free((itr->second).u_tgts);
        (itr->second).u_tgts = NULL;
      }
      if ((itr->second).data   != NULL) {
       free((itr->second).data);
       (itr->second).data = NULL;
      }
      itr = clt_cwp.coupling[s].field.erase(itr);
    }
  }

  clt_cwp.coupling.erase(s);
}

void
CWP_client_Properties_dump
(
void
)
{
  t_message msg;

  // verbose
  MPI_Barrier(clt->comm);
  if ((clt->flags  & CWP_FLAG_VERBOSE) && (clt->i_rank == 0)) {
    PDM_printf("%s-CWP-CLIENT: Client initiating CWP_Properties_dump\n", clt->code_name);
    PDM_printf_flush();
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_PROPERTIES_DUMP);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Properties_dump failed to send message header\n");
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // read
  int size = -1;
  CWP_transfer_readdata(clt->socket, clt->max_msg_size, (void*) &size, sizeof(int));
  char *recv = (char *) malloc(sizeof(char) * (size + 1));
  memset(recv, 0, (size + 1));
  CWP_transfer_readdata(clt->socket, clt->max_msg_size, (void*) recv, size);

  // out
  printf("%s\n", recv);

  // free
  free(recv);
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
  MPI_Barrier(clt->comm);
  if ((clt->flags  & CWP_FLAG_VERBOSE) && (clt->i_rank == 0)) {
    PDM_printf("%s-CWP-CLIENT: Client initiating CWP_Visu_set\n", clt->code_name);
    PDM_printf_flush();
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_VISU_SET);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Visu_set failed to send message header\n");
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
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

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }
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
  MPI_Barrier(clt->comm);
  if ((clt->flags  & CWP_FLAG_VERBOSE) && (clt->i_rank == 0)) {
    PDM_printf("%s-CWP-CLIENT: Client initiating CWP_State_update\n", clt->code_name);
    PDM_printf_flush();
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_STATE_UPDATE);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_State_update failed to send message header\n");
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // send code name
  write_name(local_code_name);

  // send state
  CWP_State_t endian_state = state;
  CWP_swap_endian_4bytes((int *) &endian_state, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_state, sizeof(CWP_State_t));

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }
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
  MPI_Barrier(clt->comm);
  if ((clt->flags  & CWP_FLAG_VERBOSE) && (clt->i_rank == 0)) {
    PDM_printf("%s-CWP-CLIENT: Client initiating CWP_Time_update\n", clt->code_name);
    PDM_printf_flush();
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_TIME_UPDATE);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Time_update failed to send message header\n");
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // send code name
  write_name(local_code_name);

  // send time
  double endian_current_time = current_time;
  CWP_swap_endian_8bytes(&endian_current_time, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_current_time, sizeof(double));

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }
}

void
CWP_client_Output_file_set
(
  FILE *output_file
)
{
  // set file for PDM_printf
  _cwipi_output_listing = output_file;
  PDM_printf_proxy_set(_cwipi_print_with_c);
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
  printf("CWP-CLIENT: CWP_User_structure_set not implemented in client/server mode\n");
}

void *
CWP_client_User_structure_get
(
 const char* local_code_name
)
{
  PDM_UNUSED(local_code_name);
  printf("CWP-CLIENT: CWP_User_structure_get not implemented in client/server mode\n");

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
  MPI_Barrier(clt->comm);
  if ((clt->flags  & CWP_FLAG_VERBOSE) && (clt->i_rank == 0)) {
    PDM_printf("%s-CWP-CLIENT: Client initiating CWP_State_get\n", clt->code_name);
    PDM_printf_flush();
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_STATE_GET);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_State_get failed to send message header\n");
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // send code name
  write_name(code_name);

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

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
  MPI_Barrier(clt->comm);
  if ((clt->flags  & CWP_FLAG_VERBOSE) && (clt->i_rank == 0)) {
    PDM_printf("%s-CWP-CLIENT: Client initiating CWP_Codes_nb_get\n", clt->code_name);
    PDM_printf_flush();
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_CODES_NB_GET);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Codes_nb_get failed to send message header\n");
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // read nb_codes
  int n_code_names = -1;
  CWP_transfer_readdata(clt->socket, clt->max_msg_size, (void*) &n_code_names, sizeof(int));

  return n_code_names;
}

const char **
CWP_client_Codes_list_get
(
void
)
{
  t_message msg;

  // verbose
  MPI_Barrier(clt->comm);
  if ((clt->flags  & CWP_FLAG_VERBOSE) && (clt->i_rank == 0)) {
    PDM_printf("%s-CWP-CLIENT: Client initiating CWP_Codes_list_get\n", clt->code_name);
    PDM_printf_flush();
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_CODES_LIST_GET);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Codes_list_get failed to send message header\n");
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // read code names
  clt_cwp.n_code_names = -1;
  CWP_transfer_readdata(clt->socket, clt->max_msg_size, (void*) &(clt_cwp.n_code_names), sizeof(int));
  clt_cwp.code_names = (char **) malloc(sizeof(char *) * clt_cwp.n_code_names);
  for (int i = 0; i < clt_cwp.n_code_names; i++) {
    (clt_cwp.code_names)[i] = (char *) malloc(sizeof(char));
    read_name((char **) &((clt_cwp.code_names)[i]));
  }

  return (const char**) clt_cwp.code_names;
}

int
CWP_client_Loc_codes_nb_get
(
 void
)
{
  t_message msg;

  // verbose
  MPI_Barrier(clt->comm);
  if ((clt->flags  & CWP_FLAG_VERBOSE) && (clt->i_rank == 0)) {
    PDM_printf("%s-CWP-CLIENT: Client initiating CWP_Loc_codes_nb_get\n", clt->code_name);
    PDM_printf_flush();
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_LOC_CODES_NB_GET);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Loc_codes_nb_get failed to send message header\n");
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // read nb_codes
  int n_loc_code_names = -1;
  CWP_transfer_readdata(clt->socket, clt->max_msg_size, (void*) &n_loc_code_names, sizeof(int));

  return n_loc_code_names;
}

const char **
CWP_client_Loc_codes_list_get
(
 void
)
{
  t_message msg;

  // verbose
  MPI_Barrier(clt->comm);
  if ((clt->flags  & CWP_FLAG_VERBOSE) && (clt->i_rank == 0)) {
    PDM_printf("%s-CWP-CLIENT: Client initiating CWP_Loc_codes_list_get\n", clt->code_name);
    PDM_printf_flush();
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_LOC_CODES_LIST_GET);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Loc_codes_list_get failed to send message header\n");
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // read code names
  clt_cwp.n_loc_code_names = -1;
  CWP_transfer_readdata(clt->socket, clt->max_msg_size, (void*) &(clt_cwp.n_loc_code_names), sizeof(int));
  clt_cwp.loc_code_names = (char **) malloc(sizeof(char *) * clt_cwp.n_loc_code_names);
  for (int i = 0; i < clt_cwp.n_loc_code_names; i++) {
    (clt_cwp.loc_code_names)[i] = (char *) malloc(sizeof(char));
    read_name((char **) &((clt_cwp.loc_code_names)[i]));
  }

  return (const char **) clt_cwp.loc_code_names;
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
  MPI_Barrier(clt->comm);
  if ((clt->flags  & CWP_FLAG_VERBOSE) && (clt->i_rank == 0)) {
    PDM_printf("%s-CWP-CLIENT: Client initiating CWP_N_uncomputed_tgts_get\n", clt->code_name);
    PDM_printf_flush();
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_N_UNCOMPUTED_TGTS_GET);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_N_uncomputed_tgts_get failed to send message header\n");
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
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

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

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
  MPI_Barrier(clt->comm);
  if ((clt->flags  & CWP_FLAG_VERBOSE) && (clt->i_rank == 0)) {
    PDM_printf("%s-CWP-CLIENT: Client initiating CWP_Uncomputed_tgts_get\n", clt->code_name);
    PDM_printf_flush();
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_UNCOMPUTED_TGTS_GET);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Uncomputed_tgts_get failed to send message header\n");
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
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

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // read number of targets
  int nb_tgts = -1;
  CWP_transfer_readdata(clt->socket, clt->max_msg_size, (void*) &nb_tgts, sizeof(int));

  std::string s1(cpl_id);
  std::string s2(field_id);

  clt_cwp.coupling[s1].field[s2].u_tgts = (int *) malloc(sizeof(int) * nb_tgts);
  CWP_transfer_readdata(clt->socket, clt->max_msg_size, (void*) clt_cwp.coupling[s1].field[s2].u_tgts, sizeof(int)*nb_tgts);

  return clt_cwp.coupling[s1].field[s2].u_tgts;
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
  MPI_Barrier(clt->comm);
  if ((clt->flags  & CWP_FLAG_VERBOSE) && (clt->i_rank == 0)) {
    PDM_printf("%s-CWP-CLIENT: Client initiating CWP_N_computed_tgts_get\n", clt->code_name);
    PDM_printf_flush();
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_N_COMPUTED_TGTS_GET);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_N_computed_tgts_get failed to send message header\n");
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
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

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

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
  MPI_Barrier(clt->comm);
  if ((clt->flags  & CWP_FLAG_VERBOSE) && (clt->i_rank == 0)) {
    PDM_printf("%s-CWP-CLIENT: Client initiating CWP_Computed_tgts_get\n", clt->code_name);
    PDM_printf_flush();
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_COMPUTED_TGTS_GET);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Computed_tgts_get failed to send message header\n");
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
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

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // read number of targets
  int nb_tgts = -1;
  CWP_transfer_readdata(clt->socket, clt->max_msg_size, (void*) &nb_tgts, sizeof(int));

  std::string s1(cpl_id);
  std::string s2(field_id);

  clt_cwp.coupling[s1].field[s2].c_tgts = (int *) malloc(sizeof(int) * nb_tgts);
  CWP_transfer_readdata(clt->socket, clt->max_msg_size, (void*) clt_cwp.coupling[s1].field[s2].c_tgts, sizeof(int)*nb_tgts);

  return clt_cwp.coupling[s1].field[s2].c_tgts;
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
  MPI_Barrier(clt->comm);
  if ((clt->flags  & CWP_FLAG_VERBOSE) && (clt->i_rank == 0)) {
    PDM_printf("%s-CWP-CLIENT: Client initiating CWP_N_involved_srcs_get\n", clt->code_name);
    PDM_printf_flush();
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_N_INVOLVED_SRCS_GET);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_N_involved_srcs_get failed to send message header\n");
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
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

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

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
  MPI_Barrier(clt->comm);
  if ((clt->flags  & CWP_FLAG_VERBOSE) && (clt->i_rank == 0)) {
    PDM_printf("%s-CWP-CLIENT: Client initiating CWP_Involved_srcs_get\n", clt->code_name);
    PDM_printf_flush();
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_INVOLVED_SRCS_GET);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Involved_srcs_get failed to send message header\n");
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
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

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // read number of sources
  int nb_srcs = -1;
  CWP_transfer_readdata(clt->socket, clt->max_msg_size, (void*) &nb_srcs, sizeof(int));

  std::string s1(cpl_id);
  std::string s2(field_id);

  clt_cwp.coupling[s1].field[s2].srcs = (int *) malloc(sizeof(int) * nb_srcs);
  CWP_transfer_readdata(clt->socket, clt->max_msg_size, (void*) clt_cwp.coupling[s1].field[s2].srcs, sizeof(int)*nb_srcs);

  return clt_cwp.coupling[s1].field[s2].srcs;
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
  MPI_Barrier(clt->comm);
  if ((clt->flags  & CWP_FLAG_VERBOSE) && (clt->i_rank == 0)) {
    PDM_printf("%s-CWP-CLIENT: Client initiating CWP_Spatial_interp_weights_compute\n", clt->code_name);
    PDM_printf_flush();
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_SPATIAL_INTERP_WEIGHTS_COMPUTE);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Spatial_interp_weights_compute failed to send message header\n");
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // send local code name
  write_name(local_code_name);

  // send coupling identifier
  write_name(cpl_id);

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }
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
  MPI_Barrier(clt->comm);
  if ((clt->flags  & CWP_FLAG_VERBOSE) && (clt->i_rank == 0)) {
    PDM_printf("%s-CWP-CLIENT: Client initiating CWP_Spatial_interp_property_set\n", clt->code_name);
    PDM_printf_flush();
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_SPATIAL_INTERP_PROPERTY_SET);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Spatial_interp_property_set failed to send message header\n");
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
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

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }
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
  MPI_Barrier(clt->comm);
  if ((clt->flags  & CWP_FLAG_VERBOSE) && (clt->i_rank == 0)) {
    PDM_printf("%s-CWP-CLIENT: Client initiating CWP_User_tgt_pts_set\n", clt->code_name);
    PDM_printf_flush();
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_USER_TGT_PTS_SET);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_User_tgt_pts_set failed to send message header\n");
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
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
  int NULL_flag = 0;
  if (global_num == NULL) {
    NULL_flag = 1;
  }
  CWP_swap_endian_4bytes(&NULL_flag, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &NULL_flag, sizeof(int));

  CWP_g_num_t *endian_global_num = (CWP_g_num_t *) malloc(sizeof(CWP_g_num_t) * n_pts);
  if (!NULL_flag) {
    memcpy(endian_global_num, global_num, sizeof(CWP_g_num_t) * n_pts);
    if (sizeof(CWP_g_num_t) == 4) {
      CWP_swap_endian_4bytes((int *) endian_global_num, n_pts);
    }
    else  if (sizeof(CWP_g_num_t) == 8) {
      CWP_swap_endian_8bytes((double *) endian_global_num, n_pts);
    }
    CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) endian_global_num, sizeof(CWP_g_num_t) * n_pts);
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // free
  free(endian_coord);
  if (endian_global_num != NULL) free(endian_global_num);
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
  MPI_Barrier(clt->comm);
  if ((clt->flags  & CWP_FLAG_VERBOSE) && (clt->i_rank == 0)) {
    PDM_printf("%s-CWP-CLIENT: Client initiating CWP_Mesh_interf_finalize\n", clt->code_name);
    PDM_printf_flush();
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_MESH_INTERF_FINALIZE);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Mesh_interf_finalize failed to send message header\n");
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // send local code name
  write_name(local_code_name);

  // send coupling identifier
  write_name(cpl_id);

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }
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
  MPI_Barrier(clt->comm);
  if ((clt->flags  & CWP_FLAG_VERBOSE) && (clt->i_rank == 0)) {
    PDM_printf("%s-CWP-CLIENT: Client initiating CWP_Mesh_interf_vtx_set\n", clt->code_name);
    PDM_printf_flush();
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_MESH_INTERF_VTX_SET);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Mesh_interf_vtx_set failed to send message header\n");
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
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
  int NULL_flag = 0;
  if (global_num == NULL) {
    NULL_flag = 1;
  }
  CWP_swap_endian_4bytes(&NULL_flag, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &NULL_flag, sizeof(int));

  CWP_g_num_t *endian_global_num = (CWP_g_num_t *) malloc(sizeof(CWP_g_num_t) * n_pts);
  if (!NULL_flag) {
    memcpy(endian_global_num, global_num, sizeof(CWP_g_num_t) * n_pts);
    if (sizeof(CWP_g_num_t) == 4) {
      CWP_swap_endian_4bytes((int *) endian_global_num, n_pts);
    }
    else  if (sizeof(CWP_g_num_t) == 8) {
      CWP_swap_endian_8bytes((double *) endian_global_num, n_pts);
    }
    CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) endian_global_num, sizeof(CWP_g_num_t) * n_pts);
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // free
  free(endian_coord);
  if (endian_global_num != NULL) free(endian_global_num);
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
  MPI_Barrier(clt->comm);
  if ((clt->flags  & CWP_FLAG_VERBOSE) && (clt->i_rank == 0)) {
    PDM_printf("%s-CWP-CLIENT: Client initiating CWP_Mesh_interf_block_add\n", clt->code_name);
    PDM_printf_flush();
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_MESH_INTERF_BLOCK_ADD);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Mesh_interf_block_add failed to send message header\n");
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // send local code name
  write_name(local_code_name);

  // send coupling identifier
  write_name(cpl_id);

  // send block type
  CWP_Block_t endian_block_type = block_type;
  CWP_swap_endian_4bytes((int *) &endian_block_type, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_block_type, sizeof(CWP_Block_t));

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

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
 const int          n_elts,
 int                connec[],
 CWP_g_num_t        global_num[]
)
{
  t_message msg;

  // verbose
  MPI_Barrier(clt->comm);
  if ((clt->flags  & CWP_FLAG_VERBOSE) && (clt->i_rank == 0)) {
    PDM_printf("%s-CWP-CLIENT: Client initiating CWP_Mesh_interf_block_std_set\n", clt->code_name);
    PDM_printf_flush();
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_MESH_INTERF_BLOCK_STD_SET);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Mesh_interf_block_std_set failed to send message header\n");
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
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
  CWP_Block_t block_type;
  CWP_transfer_readdata(clt->socket, clt->max_msg_size,(void *) &block_type, sizeof(int));
  int n_vtx_elt = n_vtx_get(block_type);
  int endian_n_vtx_elt = n_vtx_elt;
  CWP_swap_endian_4bytes(&endian_n_vtx_elt, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_n_vtx_elt, sizeof(int));

  // send connectivity
  int *endian_connec = (int *) malloc(sizeof(int) * n_elts * n_vtx_elt);
  memcpy(endian_connec, connec, sizeof(int) * n_elts * n_vtx_elt);
  CWP_swap_endian_4bytes(endian_connec, n_elts * n_vtx_elt);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) endian_connec, sizeof(int) * n_elts * n_vtx_elt);

  // send global number
  int NULL_flag = 0;
  if (global_num == NULL) {
    NULL_flag = 1;
  }
  CWP_swap_endian_4bytes(&NULL_flag, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &NULL_flag, sizeof(int));

  CWP_g_num_t  *endian_global_num = (CWP_g_num_t  *) malloc(sizeof(CWP_g_num_t) * n_elts);
  if (!NULL_flag) {
    memcpy(endian_global_num, global_num, sizeof(CWP_g_num_t) * n_elts);
    if (sizeof(CWP_g_num_t) == 4) {
      CWP_swap_endian_4bytes((int *) endian_global_num, n_elts);
    }
    else  if (sizeof(CWP_g_num_t) == 8) {
      CWP_swap_endian_8bytes((double *) endian_global_num, n_elts);
    }
    CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) endian_global_num, sizeof(CWP_g_num_t) * n_elts);
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // free
  free(endian_connec);
  if (endian_global_num != NULL) free(endian_global_num);
}

void
CWP_client_Mesh_interf_block_std_get
(
 const char        *local_code_name,
 const char        *cpl_id,
 const int          i_part,
 const int          block_id,
 int               *n_elts,
 int              **connec,
 CWP_g_num_t      **global_num
)
{
  t_message msg;

  // verbose
  MPI_Barrier(clt->comm);
  if ((clt->flags  & CWP_FLAG_VERBOSE) && (clt->i_rank == 0)) {
    PDM_printf("%s-CWP-CLIENT: Client initiating CWP_Mesh_interf_block_std_get\n", clt->code_name);
    PDM_printf_flush();
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_MESH_INTERF_BLOCK_STD_GET);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Mesh_interf_block_std_get failed to send message header\n");
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
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

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // read n_elts
  CWP_transfer_readdata(clt->socket, clt->max_msg_size, n_elts, sizeof(int));

  std::string s(cpl_id);

  // read connectivity
  int n_vtx_elt = -1;
  CWP_transfer_readdata(clt->socket, clt->max_msg_size, &n_vtx_elt, sizeof(int));
  clt_cwp.coupling[s].std_connec = (int *) malloc(sizeof(int) * (n_vtx_elt * (*n_elts)));
  CWP_transfer_readdata(clt->socket, clt->max_msg_size, clt_cwp.coupling[s].std_connec, sizeof(int) * (n_vtx_elt * (*n_elts)));
  *connec = clt_cwp.coupling[s].std_connec;

  // read global number
  int NULL_flag = -1;
  CWP_transfer_readdata(clt->socket, clt->max_msg_size, &NULL_flag, sizeof(int));

  if (!NULL_flag) {
    clt_cwp.coupling[s].std_global_num = (CWP_g_num_t *) malloc(sizeof(CWP_g_num_t) * (*n_elts));
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, clt_cwp.coupling[s].std_global_num, sizeof(CWP_g_num_t) * (*n_elts));
    *global_num = clt_cwp.coupling[s].std_global_num;
  } else {
    *global_num = NULL;
  }

}

CWP_Block_t
CWP_client_std_block_type_get
(
 const char             *local_code_name,
 const char             *cpl_id,
 const int               block_id
)
{
  t_message msg;

  // verbose
  MPI_Barrier(clt->comm);
  if ((clt->flags  & CWP_FLAG_VERBOSE) && (clt->i_rank == 0)) {
    PDM_printf("%s-CWP-CLIENT: Client initiating CWP_std_block_type_get\n", clt->code_name);
    PDM_printf_flush();
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_STD_BLOCK_TYPE_GET);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_std_block_type_get failed to send message header\n");
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // send local code name
  write_name(local_code_name);

  // send coupling identifier
  write_name(cpl_id);

  // send block_id
  int endian_block_id = block_id;
  CWP_swap_endian_4bytes(&endian_block_id, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_block_id, sizeof(int));

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // read block_type
  CWP_Block_t block_type;
  CWP_transfer_readdata(clt->socket, clt->max_msg_size,(void*) &block_type, sizeof(int));

  return block_type;
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
  MPI_Barrier(clt->comm);
  if ((clt->flags  & CWP_FLAG_VERBOSE) && (clt->i_rank == 0)) {
    PDM_printf("%s-CWP-CLIENT: Client initiating CWP_Mesh_interf_f_poly_block_set\n", clt->code_name);
    PDM_printf_flush();
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_MESH_INTERF_F_POLY_BLOCK_SET);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Mesh_interf_f_poly_block_set failed to send message header\n");
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
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
  int NULL_flag = 0;
  if (global_num == NULL) {
    NULL_flag = 1;
  }
  CWP_swap_endian_4bytes(&NULL_flag, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &NULL_flag, sizeof(int));

  CWP_g_num_t  *endian_global_num = (CWP_g_num_t  *) malloc(sizeof(CWP_g_num_t) * n_elts);
  if (!NULL_flag) {
    memcpy(endian_global_num, global_num, sizeof(CWP_g_num_t) * n_elts);
    if (sizeof(CWP_g_num_t) == 4) {
      CWP_swap_endian_4bytes((int *) endian_global_num, n_elts);
    }
    else  if (sizeof(CWP_g_num_t) == 8) {
      CWP_swap_endian_8bytes((double *) endian_global_num, n_elts);
    }
    CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) endian_global_num, sizeof(CWP_g_num_t) * n_elts);
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // free
  free(endian_connec_idx);
  free(endian_connec);
  if (endian_global_num != NULL) free(endian_global_num);
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
  MPI_Barrier(clt->comm);
  if ((clt->flags  & CWP_FLAG_VERBOSE) && (clt->i_rank == 0)) {
    PDM_printf("%s-CWP-CLIENT: Client initiating CWP_Mesh_interf_f_poly_block_get\n", clt->code_name);
    PDM_printf_flush();
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_MESH_INTERF_F_POLY_BLOCK_GET);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Mesh_interf_f_poly_block_get failed to send message header\n");
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
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

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // read n_elts
  CWP_transfer_readdata(clt->socket, clt->max_msg_size, n_elts, sizeof(int));

  std::string s(cpl_id);

  // read connectivity index
  clt_cwp.coupling[s].connec_idx = (int *) malloc(sizeof(int) * ((*n_elts)+1));
  CWP_transfer_readdata(clt->socket, clt->max_msg_size, clt_cwp.coupling[s].connec_idx, sizeof(int) * ((*n_elts)+1));
  *connec_idx = clt_cwp.coupling[s].connec_idx;

  // read connectivity
  clt_cwp.coupling[s].connec = (int *) malloc(sizeof(int) * (*connec_idx)[(*n_elts)]);
  CWP_transfer_readdata(clt->socket, clt->max_msg_size, clt_cwp.coupling[s].connec, sizeof(int) * (*connec_idx)[(*n_elts)]);
  *connec = clt_cwp.coupling[s].connec;

  // read global number
  int NULL_flag = -1;
  CWP_transfer_readdata(clt->socket, clt->max_msg_size, &NULL_flag, sizeof(int));

  if (!NULL_flag) {
    clt_cwp.coupling[s].elt_global_num = (CWP_g_num_t *) malloc(sizeof(CWP_g_num_t) * (*n_elts));
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, clt_cwp.coupling[s].elt_global_num, sizeof(CWP_g_num_t) * (*n_elts));
    *global_num = clt_cwp.coupling[s].elt_global_num;
  } else {
    *global_num = NULL;
  }
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
  MPI_Barrier(clt->comm);
  if ((clt->flags  & CWP_FLAG_VERBOSE) && (clt->i_rank == 0)) {
    PDM_printf("%s-CWP-CLIENT: Client initiating CWP_Mesh_interf_c_poly_block_set\n", clt->code_name);
    PDM_printf_flush();
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_MESH_INTERF_C_POLY_BLOCK_SET);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Mesh_interf_c_poly_block_set failed to send message header\n");
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
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
  int NULL_flag = 0;
  if (global_num == NULL) {
    NULL_flag = 1;
  }
  CWP_swap_endian_4bytes(&NULL_flag, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &NULL_flag, sizeof(int));

  CWP_g_num_t  *endian_global_num = (CWP_g_num_t  *) malloc(sizeof(CWP_g_num_t) * n_elts);
  if (!NULL_flag) {
    memcpy(endian_global_num, global_num, sizeof(CWP_g_num_t) * n_elts);
    if (sizeof(CWP_g_num_t) == 4) {
      CWP_swap_endian_4bytes((int *) endian_global_num, n_elts);
    }
    else  if (sizeof(CWP_g_num_t) == 8) {
      CWP_swap_endian_8bytes((double *) endian_global_num, n_elts);
    }
    CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) endian_global_num, sizeof(CWP_g_num_t) * n_elts);
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // free
  free(endian_connec_faces_idx);
  free(endian_connec_faces);
  free(endian_connec_cells_idx);
  free(endian_connec_cells);
  if (endian_global_num != NULL) free(endian_global_num);
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
  MPI_Barrier(clt->comm);
  if ((clt->flags  & CWP_FLAG_VERBOSE) && (clt->i_rank == 0)) {
    PDM_printf("%s-CWP-CLIENT: Client initiating CWP_Mesh_interf_c_poly_block_get\n", clt->code_name);
    PDM_printf_flush();
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_MESH_INTERF_C_POLY_BLOCK_GET);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Mesh_interf_c_poly_block_get failed to send message header\n");
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
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

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // read n_elts
  CWP_transfer_readdata(clt->socket, clt->max_msg_size, n_elts, sizeof(int));

  // read n_faces
  CWP_transfer_readdata(clt->socket, clt->max_msg_size, n_faces, sizeof(int));

  std::string s(cpl_id);

  // read connectivity index
  clt_cwp.coupling[s].connec_faces_idx = (int *) malloc(sizeof(int) * ((*n_faces)+1));
  CWP_transfer_readdata(clt->socket, clt->max_msg_size, clt_cwp.coupling[s].connec_faces_idx, sizeof(int) * ((*n_faces)+1));
  *connec_faces_idx = clt_cwp.coupling[s].connec_faces_idx;

  // read connectivity
  clt_cwp.coupling[s].connec_faces = (int *) malloc(sizeof(int) * (*connec_faces_idx)[(*n_faces)]);
  CWP_transfer_readdata(clt->socket, clt->max_msg_size, clt_cwp.coupling[s].connec_faces, sizeof(int) * (*connec_faces_idx)[(*n_faces)]);
  *connec_faces = clt_cwp.coupling[s].connec_faces;

  // read connectivity index
  clt_cwp.coupling[s].connec_cells_idx = (int *) malloc(sizeof(int) * ((*n_elts)+1));
  CWP_transfer_readdata(clt->socket, clt->max_msg_size, clt_cwp.coupling[s].connec_cells_idx, sizeof(int) * ((*n_elts)+1));
  *connec_cells_idx = clt_cwp.coupling[s].connec_cells_idx;

  // read connectivity
  clt_cwp.coupling[s].connec_cells = (int *) malloc(sizeof(int) * (*connec_cells_idx)[(*n_elts)]);
  CWP_transfer_readdata(clt->socket, clt->max_msg_size, clt_cwp.coupling[s].connec_cells, sizeof(int) * (*connec_cells_idx)[(*n_elts)]);
  *connec_cells = clt_cwp.coupling[s].connec_cells;

  // read global number
  int NULL_flag = -1;
  CWP_transfer_readdata(clt->socket, clt->max_msg_size, &NULL_flag, sizeof(int));

  if (!NULL_flag) {
    clt_cwp.coupling[s].cell_global_num = (CWP_g_num_t *) malloc(sizeof(CWP_g_num_t) * (*n_elts));
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, clt_cwp.coupling[s].cell_global_num, sizeof(CWP_g_num_t) * (*n_elts));
    *global_num = clt_cwp.coupling[s].cell_global_num;
  } else {
    *global_num = NULL;
  }

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
  MPI_Barrier(clt->comm);
  if ((clt->flags  & CWP_FLAG_VERBOSE) && (clt->i_rank == 0)) {
    PDM_printf("%s-CWP-CLIENT: Client initiating CWP_Mesh_interf_del\n", clt->code_name);
    PDM_printf_flush();
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_MESH_INTERF_DEL);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Mesh_interf_del failed to send message header\n");
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // send local code name
  write_name(local_code_name);

  // send coupling identifier
  write_name(cpl_id);

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }
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
  MPI_Barrier(clt->comm);
  if ((clt->flags  & CWP_FLAG_VERBOSE) && (clt->i_rank == 0)) {
    PDM_printf("%s-CWP-CLIENT: Client initiating CWP_Mesh_interf_from_cellface_set\n", clt->code_name);
    PDM_printf_flush();
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_MESH_INTERF_FROM_CELLFACE_SET);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Mesh_interf_from_cellface_set failed to send message header\n");
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
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
  int NULL_flag = 0;
  if (global_num == NULL) {
    NULL_flag = 1;
  }
  CWP_swap_endian_4bytes(&NULL_flag, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &NULL_flag, sizeof(int));

  CWP_g_num_t  *endian_global_num = (CWP_g_num_t  *) malloc(sizeof(CWP_g_num_t) * n_cells);
  if (!NULL_flag) {
    memcpy(endian_global_num, global_num, sizeof(CWP_g_num_t) * n_cells);
    if (sizeof(CWP_g_num_t) == 4) {
      CWP_swap_endian_4bytes((int *) endian_global_num, n_cells);
    }
    else  if (sizeof(CWP_g_num_t) == 8) {
      CWP_swap_endian_8bytes((double *) endian_global_num, n_cells);
    }
    CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) endian_global_num, sizeof(CWP_g_num_t) * n_cells);
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // free
  free(endian_cell_face_idx);
  free(endian_cell_face);
  free(endian_face_vtx_idx);
  free(endian_face_vtx);
  if (endian_global_num != NULL) free(endian_global_num);
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
  MPI_Barrier(clt->comm);
  if ((clt->flags  & CWP_FLAG_VERBOSE) && (clt->i_rank == 0)) {
    PDM_printf("%s-CWP-CLIENT: Client initiating CWP_Mesh_interf_from_faceedge_set\n", clt->code_name);
    PDM_printf_flush();
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_MESH_INTERF_FROM_FACEEDGE_SET);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Mesh_interf_from_faceedge_set failed to send message header\n");
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
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
  int *endian_edge_vtx = (int *) malloc(sizeof(int) * 2 * n_edges);
  memcpy(endian_edge_vtx, edge_vtx, sizeof(int) * 2 * n_edges);
  CWP_swap_endian_4bytes(endian_edge_vtx, 2 * n_edges);
  CWP_transfer_writedata(clt->socket, clt->max_msg_size, (void *) endian_edge_vtx, sizeof(int) * 2 * n_edges);

  // send global number
  int NULL_flag = 0;
  if (global_num == NULL) {
    NULL_flag = 1;
  }
  CWP_swap_endian_4bytes(&NULL_flag, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &NULL_flag, sizeof(int));

  CWP_g_num_t  *endian_global_num = (CWP_g_num_t  *) malloc(sizeof(CWP_g_num_t) * n_faces);
  if (!NULL_flag) {
    memcpy(endian_global_num, global_num, sizeof(CWP_g_num_t) * n_faces);
    if (sizeof(CWP_g_num_t) == 4) {
      CWP_swap_endian_4bytes((int *) endian_global_num, n_faces);
    }
    else  if (sizeof(CWP_g_num_t) == 8) {
      CWP_swap_endian_8bytes((double *) endian_global_num, n_faces);
    }
    CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) endian_global_num, sizeof(CWP_g_num_t) * n_faces);
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // free
  free(endian_face_edge_idx);
  free(endian_face_edge);
  free(endian_edge_vtx_idx);
  free(endian_edge_vtx);
  if (endian_global_num != NULL) free(endian_global_num);
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
  MPI_Barrier(clt->comm);
  if ((clt->flags  & CWP_FLAG_VERBOSE) && (clt->i_rank == 0)) {
    PDM_printf("%s-CWP-CLIENT: Client initiating CWP_Field_create\n", clt->code_name);
    PDM_printf_flush();
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_FIELD_CREATE);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Field_create failed to send message header\n");
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
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

  // add field
  std::string s1(local_code_name);
  std::string s2(cpl_id);
  std::string s3(field_id);

  fields.field_settings[std::make_tuple(s1, s2, s3)] = t_field_settings();
  fields.field_settings[std::make_tuple(s1, s2, s3)].n_component = n_component;

  t_coupling coupling = clt_cwp.coupling[s2];
  t_field field = t_field();
  coupling.field.insert(std::make_pair(s3, field));

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

}

void
CWP_client_Field_data_set
(
 const char              *local_code_name,
 const char              *cpl_id,
 const char              *field_id,
 const int                i_part,
 const CWP_Field_map_t    map_type,
 int                      n_entities,
 double                   data[]
)
{
  t_message msg;

  // verbose
  MPI_Barrier(clt->comm);
  if ((clt->flags  & CWP_FLAG_VERBOSE) && (clt->i_rank == 0)) {
    PDM_printf("%s-CWP-CLIENT: Client initiating CWP_Field_data_set\n", clt->code_name);
    PDM_printf_flush();
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_FIELD_DATA_SET);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Field_data_set failed to send message header\n");
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
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

  // complete field settings
  std::string s1(local_code_name);
  std::string s2(cpl_id);
  std::string s3(field_id);

  fields.field_settings[std::make_tuple(s1, s2, s3)].i_part     = i_part;
  fields.field_settings[std::make_tuple(s1, s2, s3)].n_entities = n_entities;
  fields.field_settings[std::make_tuple(s1, s2, s3)].map_type   = map_type;

  // send map with data
  int size = fields.field_settings[std::make_tuple(s1, s2, s3)].n_component * n_entities;
  int endian_size = size;
  CWP_swap_endian_4bytes(&endian_size, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_size, sizeof(int));
  double *endian_data = (double *) malloc(sizeof(double) * size);
  memcpy(endian_data, data, sizeof(double) * size);
  CWP_swap_endian_8bytes(endian_data, size);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) endian_data, sizeof(double) * size);

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

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
  MPI_Barrier(clt->comm);
  if ((clt->flags  & CWP_FLAG_VERBOSE) && (clt->i_rank == 0)) {
    PDM_printf("%s-CWP-CLIENT: Client initiating CWP_Field_n_component_get\n", clt->code_name);
    PDM_printf_flush();
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_FIELD_N_COMPONENT_GET);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Field_n_component_get failed to send message header\n");
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // send local code name
  write_name(local_code_name);

  // send coupling identifier
  write_name(cpl_id);

  // send field identifier
  write_name(field_id);

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

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
  MPI_Barrier(clt->comm);
  if ((clt->flags  & CWP_FLAG_VERBOSE) && (clt->i_rank == 0)) {
    PDM_printf("%s-CWP-CLIENT: Client initiating CWP_Field_target_dof_location_get\n", clt->code_name);
    PDM_printf_flush();
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_FIELD_TARGET_DOF_LOCATION_GET);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Field_target_dof_location_get failed to send message header\n");
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // send local code name
  write_name(local_code_name);

  // send coupling identifier
  write_name(cpl_id);

  // send field identifier
  write_name(field_id);

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

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
  MPI_Barrier(clt->comm);
  if ((clt->flags  & CWP_FLAG_VERBOSE) && (clt->i_rank == 0)) {
    PDM_printf("%s-CWP-CLIENT: Client initiating CWP_Field_storage_get\n", clt->code_name);
    PDM_printf_flush();
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_FIELD_STORAGE_GET);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Field_storage_get failed to send message header\n");
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // send local code name
  write_name(local_code_name);

  // send coupling identifier
  write_name(cpl_id);

  // send field identifier
  write_name(field_id);

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

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
  MPI_Barrier(clt->comm);
  if ((clt->flags  & CWP_FLAG_VERBOSE) && (clt->i_rank == 0)) {
    PDM_printf("%s-CWP-CLIENT: Client initiating CWP_Field_del\n", clt->code_name);
    PDM_printf_flush();
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_FIELD_DEL);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Field_del failed to send message header\n");
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // send local code name
  write_name(local_code_name);

  // send coupling identifier
  write_name(cpl_id);

  // send field identifier
  write_name(field_id);

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // free
  std::string s1(cpl_id);
  t_coupling coupling = clt_cwp.coupling[s1];
  std::string s2(field_id);
  t_field field = coupling.field[s2];

  if (field.data != NULL  ) free(field.data);
  field.data = NULL;
  if (field.u_tgts != NULL) free(field.u_tgts);
  field.u_tgts = NULL;
  if (field.c_tgts != NULL) free(field.c_tgts);
  field.u_tgts = NULL;
  if (field.srcs != NULL  ) free(field.srcs);
  field.srcs = NULL;

  coupling.field.erase(s2);
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
  MPI_Barrier(clt->comm);
  if ((clt->flags  & CWP_FLAG_VERBOSE) && (clt->i_rank == 0)) {
    PDM_printf("%s-CWP-CLIENT: Client initiating CWP_Field_issend\n", clt->code_name);
    PDM_printf_flush();
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_FIELD_ISSEND);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Field_issend failed to send message header\n");
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // send local code name
  write_name(local_code_name);

  // send coupling identifier
  write_name(cpl_id);

  // send source field identifier
  write_name(src_field_id);

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }
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
  MPI_Barrier(clt->comm);
  if ((clt->flags  & CWP_FLAG_VERBOSE) && (clt->i_rank == 0)) {
    PDM_printf("%s-CWP-CLIENT: Client initiating CWP_Field_irecv\n", clt->code_name);
    PDM_printf_flush();
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_FIELD_IRECV);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Field_irecv failed to send message header\n");
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // send local code name
  write_name(local_code_name);

  // send coupling identifier
  write_name(cpl_id);

  // send target field identifier
  write_name(tgt_field_id);

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }
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
  MPI_Barrier(clt->comm);
  if ((clt->flags  & CWP_FLAG_VERBOSE) && (clt->i_rank == 0)) {
    PDM_printf("%s-CWP-CLIENT: Client initiating CWP_Field_wait_issend\n", clt->code_name);
    PDM_printf_flush();
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_FIELD_WAIT_ISSEND);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Field_wait_issend failed to send message header\n");
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // send local code name
  write_name(local_code_name);

  // send coupling identifier
  write_name(cpl_id);

  // send source field identifier
  write_name(src_field_id);

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }
}

int
CWP_client_Field_wait_irecv
(
 const char              *local_code_name,
 const char              *cpl_id,
 const char              *tgt_field_id,
 double                 **data
)
{
  t_message msg;

  // verbose
  MPI_Barrier(clt->comm);
  if ((clt->flags  & CWP_FLAG_VERBOSE) && (clt->i_rank == 0)) {
    PDM_printf("%s-CWP-CLIENT: Client initiating CWP_Field_wait_irecv\n", clt->code_name);
  } // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_FIELD_WAIT_IRECV);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Field_wait_irecv failed to send message header\n");
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // send local code name
  write_name(local_code_name);

  // send coupling identifier
  write_name(cpl_id);

  // send target field identifier
  write_name(tgt_field_id);

  // send data needed to retreive field
  std::string s1(local_code_name);
  std::string s2(cpl_id);
  std::string s3(tgt_field_id);

  // send i_part
  int endian_i_part = fields.field_settings[std::make_tuple(s1, s2, s3)].i_part;
  CWP_swap_endian_4bytes(&endian_i_part, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_i_part, sizeof(int));

  // send map type
  CWP_Field_map_t endian_map_type = fields.field_settings[std::make_tuple(s1, s2, s3)].map_type;
  CWP_swap_endian_4bytes((int *) &endian_map_type, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_map_type, sizeof(CWP_Field_map_t));

  // send size
  int size = fields.field_settings[std::make_tuple(s1, s2, s3)].n_component * fields.field_settings[std::make_tuple(s1, s2, s3)].n_entities;
  int endian_size = size;
  CWP_swap_endian_4bytes(&endian_size, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_size, sizeof(int));

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // read irecv Field_data
  clt_cwp.coupling[s2].field[s3].data = (double *) malloc(sizeof(double) * size);
  CWP_transfer_readdata(clt->socket, clt->max_msg_size, clt_cwp.coupling[s2].field[s3].data, sizeof(double) * size);
  *data = clt_cwp.coupling[s2].field[s3].data;

  return size;
}

void
CWP_client_Interp_function_unset
(
 const char                 *local_code_name,
 const char                 *cpl_id,
 const char                 *src_field_id
)
{
  // t_message msg;

  // // verbose
  // if (clt->flags & CWP_FLAG_VERBOSE) {
  //   printf("CWP-CLIENT: Client initiating CWP_Interp_function_unset\n");
  // }

  // // create message
  // NEWMESSAGE(msg, CWP_MSG_CWP_INTERP_FUNCTION_UNSET);

  // // send message
  // if (CWP_client_send_msg(&msg) != 0) {
  //   PDM_error(__FILE__, __LINE__, 0, "CWP_client_Interp_function_unset failed to send message header\n");
  // }

  // // send local code name
  // write_name(local_code_name);

  // // send coupling identifier
  // write_name(cpl_id);

  // // send source field identifier
  // write_name(src_field_id);

 PDM_UNUSED(local_code_name);
 PDM_UNUSED(cpl_id);
 PDM_UNUSED(src_field_id);

 printf("CWP_client_Interp_function_unset not implemented yet\n");
}

void
CWP_client_Interp_function_set
(
 const char                 *local_code_name,
 const char                 *cpl_id,
 const char                 *src_field_id,
 CWP_Interp_function_t  fct
)
{
 PDM_UNUSED(local_code_name);
 PDM_UNUSED(cpl_id);
 PDM_UNUSED(src_field_id);
 PDM_UNUSED(fct);

 printf("CWP_client_Interp_function_set not implemented yet\n");
 /* TO DO: standard interpolation functions could be implemented
  *        at server side and user could choose out of them or
  *        implement some others
  */
}

void
CWP_client_Mesh_interf_block_ho_set
(
 const char        *local_code_name,
 const char        *cpl_id,
 const int          i_part,
 const int          block_id,
 const int          n_elts,
 const int          order,
 int                connec[],
 CWP_g_num_t        global_num[]
)
{
  t_message msg;

  // verbose
  MPI_Barrier(clt->comm);
  if ((clt->flags  & CWP_FLAG_VERBOSE) && (clt->i_rank == 0)) {
    PDM_printf("%s-CWP-CLIENT: Client initiating CWP_Mesh_interf_h_order_block_set\n", clt->code_name);
    PDM_printf_flush();
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_MESH_INTERF_BLOCK_HO_SET);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Mesh_interf_h_order_block_set failed to send message header\n");
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
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

  // send order
  int endian_order = order;
  CWP_swap_endian_4bytes(&endian_order, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_order, sizeof(int));

  // send n_vtx_elt
  CWP_Block_t block_type;
  CWP_transfer_readdata(clt->socket, clt->max_msg_size,(void *) &block_type, sizeof(int));
  int n_vtx_elt = ho_n_vtx_get(block_type, order);
  int endian_n_vtx_elt = n_vtx_elt;
  CWP_swap_endian_4bytes(&endian_n_vtx_elt, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_n_vtx_elt, sizeof(int));

  // send connectivity
  int *endian_connec = (int *) malloc(sizeof(int) * n_elts * n_vtx_elt);
  memcpy(endian_connec, connec, sizeof(int) * n_elts * n_vtx_elt);
  CWP_swap_endian_4bytes(endian_connec, n_elts * n_vtx_elt);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) endian_connec, sizeof(int) * n_elts * n_vtx_elt);

  // send global number
  int NULL_flag = 0;
  if (global_num == NULL) {
    NULL_flag = 1;
  }
  CWP_swap_endian_4bytes(&NULL_flag, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &NULL_flag, sizeof(int));

  CWP_g_num_t  *endian_global_num = (CWP_g_num_t  *) malloc(sizeof(CWP_g_num_t) * n_elts);
  if (!NULL_flag) {
    memcpy(endian_global_num, global_num, sizeof(CWP_g_num_t) * n_elts);
    if (sizeof(CWP_g_num_t) == 4) {
      CWP_swap_endian_4bytes((int *) endian_global_num, n_elts);
    }
    else  if (sizeof(CWP_g_num_t) == 8) {
      CWP_swap_endian_8bytes((double *) endian_global_num, n_elts);
    }
    CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) endian_global_num, sizeof(CWP_g_num_t) * n_elts);
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // free
  free(endian_connec);
  if (endian_global_num != NULL) free(endian_global_num);
}

void
CWP_client_Mesh_interf_block_ho_get
(
 const char        *local_code_name,
 const char        *cpl_id,
 const int          i_part,
 const int          block_id,
 int               *n_elts,
 int               *order,
 int              **connec,
 CWP_g_num_t      **global_num
)
{
  t_message msg;

  // verbose
  MPI_Barrier(clt->comm);
  if ((clt->flags  & CWP_FLAG_VERBOSE) && (clt->i_rank == 0)) {
    PDM_printf("%s-CWP-CLIENT: Client initiating CWP_Mesh_interf_block_ho_get\n", clt->code_name);
    PDM_printf_flush();
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_MESH_INTERF_BLOCK_HO_GET);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Mesh_interf_block_ho_get failed to send message header\n");
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
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

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // read n_elts
  CWP_transfer_readdata(clt->socket, clt->max_msg_size, n_elts, sizeof(int));

  // read order
  CWP_transfer_readdata(clt->socket, clt->max_msg_size, order, sizeof(int));

  std::string s(cpl_id);

  // read connectivity
  int n_vtx_elt = -1;
  CWP_transfer_readdata(clt->socket, clt->max_msg_size, &n_vtx_elt, sizeof(int));
  clt_cwp.coupling[s].ho_std_connec = (int *) malloc(sizeof(int) * (n_vtx_elt * (*n_elts)));
  CWP_transfer_readdata(clt->socket, clt->max_msg_size, clt_cwp.coupling[s].ho_std_connec, sizeof(int) * (n_vtx_elt * (*n_elts)));
  *connec = clt_cwp.coupling[s].ho_std_connec;

  // read global number
  int NULL_flag = -1;
  CWP_transfer_readdata(clt->socket, clt->max_msg_size, &NULL_flag, sizeof(int));

  if (!NULL_flag) {
    clt_cwp.coupling[s].ho_std_global_num = (CWP_g_num_t *) malloc(sizeof(CWP_g_num_t) * (*n_elts));
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, clt_cwp.coupling[s].ho_std_global_num, sizeof(CWP_g_num_t) * (*n_elts));
    *global_num = clt_cwp.coupling[s].ho_std_global_num;
  } else {
    *global_num = NULL;
  }
}

void
CWP_client_Mesh_interf_ho_ordering_from_IJK_set
(
 const char        *local_code_name,
 const char        *cpl_id,
 const CWP_Block_t  block_type,
 const int          order,
 const int          n_nodes,
 const int         *ijk_grid
)
{
  t_message msg;

  // verbose
  MPI_Barrier(clt->comm);
  if ((clt->flags  & CWP_FLAG_VERBOSE) && (clt->i_rank == 0)) {
    PDM_printf("%s-CWP-CLIENT: Client initiating CWP_Mesh_interf_ho_ordering_from_IJK_set\n", clt->code_name);
    PDM_printf_flush();
  }

  // create message
  NEWMESSAGE(msg, CWP_MSG_CWP_MESH_INTERF_HO_ORDERING_FROM_IJK_SET);

  // send message
  if (CWP_client_send_msg(&msg) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "CWP_client_Mesh_interf_ho_ordering_from_IJK_set failed to send message header\n");
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // send local code name
  write_name(local_code_name);

  // send coupling identifier
  write_name(cpl_id);

  // send block_type
  int endian_block_type = (int) block_type;
  CWP_swap_endian_4bytes(&endian_block_type, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_block_type, sizeof(CWP_Block_t));

  // send order
  int endian_order = order;
  CWP_swap_endian_4bytes(&endian_order, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_order, sizeof(int));

  // send n_nodes
  int endian_n_nodes = n_nodes;
  CWP_swap_endian_4bytes(&endian_n_nodes, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_n_nodes, sizeof(int));

  // send connectivity
  int elt_dim = PDM_Mesh_nodal_elt_dim_get(CWP_block_type_to_PDM_elt_type(block_type));
  int endian_size = n_nodes * elt_dim;
  CWP_swap_endian_4bytes(&endian_size, 1);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) &endian_size, sizeof(int));
  int *endian_ijk_grid = (int *) malloc(sizeof(int) * n_nodes * elt_dim);
  memcpy(endian_ijk_grid, ijk_grid, sizeof(int) * n_nodes * elt_dim);
  CWP_swap_endian_4bytes(endian_ijk_grid, n_nodes * elt_dim);
  CWP_transfer_writedata(clt->socket,clt->max_msg_size,(void*) endian_ijk_grid, sizeof(int) * n_nodes * elt_dim);

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // receive status msg
  MPI_Barrier(clt->comm);
  if (clt->flags  & CWP_FLAG_VERBOSE) {
    t_message message;
    CWP_transfer_readdata(clt->socket, clt->max_msg_size, &message, sizeof(t_message));
    if (clt->i_rank == 0) verbose(message);
  }

  // free
  free(endian_ijk_grid);
}

#ifdef __cplusplus
}
#endif /* __cplusplus */
