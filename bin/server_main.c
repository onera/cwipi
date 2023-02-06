/*
  This file is part of the CWIPI library.

  Copyright (C) 2011  ONERA

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

#include "client_server/server.h"
#include <pdm_error.h>
#include <pdm_io.h>
#include <pdm_mpi.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include "pdm_logging.h"
#include "pdm_printf.h"

#include "cwp_priv.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Util functions
 *============================================================================*/

static void
_usage(int exit_code)
{
  PDM_printf
    ("\n"
     "  Usage: \n\n"
     "  -c     Filename of the server configuration file.\n\n"
     "  -p     Begin and end port number for the server sockets port range.\n\n"
     "  -cn    Code identifier. \n\n"
     "  -h     This message.\n\n");

  exit(exit_code);
}

static void
_read_args
(
 int            argc,
 char         **argv,
 char         **config,     // filename for server ip adresses + ports
 int           *port_begin, // begin of port range
 int           *port_end,   // end of port range
 char         **code_name   // code name
)
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {

    if (strcmp(argv[i], "-h") == 0)
      _usage(EXIT_SUCCESS);

    else if (strcmp(argv[i], "-c") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *config = argv[i];
      }
    }

    else if (strcmp(argv[i], "-p") == 0) {

      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *port_begin = atoi(argv[i]);
      }

      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *port_end = atoi(argv[i]);
      }
    }

    else if (strcmp(argv[i], "-cn") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *code_name = argv[i];
      }
    }

    else
      _usage(EXIT_FAILURE);
    i++;

  }

}

static int
_substrcmp
(
 char* a,
 int s_a,
 char* b,
 int s_b
)
{
  char *sub_a = malloc(s_a);
  char *sub_b = malloc(s_b);
  strncpy(sub_a, a, s_a);
  strncpy(sub_b, b, s_b);
  int out = strcmp(sub_a, sub_b);
  free(sub_a);
  free(sub_b);
  return out;
}

/*=============================================================================
 * Main
 *============================================================================*/

int
main
(
 int argc,
 char *argv[]
)
{
  // default
  char *config     = NULL;
  int   port_begin = 49100;
  int   port_end   = 49150;
  char* code_name  = NULL;

  _read_args(argc,
             argv,
             &config,
             &port_begin,
             &port_end,
             &code_name);

  if (code_name == NULL) {
    PDM_error(__FILE__, __LINE__, 0, "Server must be launched with a non NULL code identifier.\n");
  }

  if (config == NULL) {
    config = (char *) "cwp_config_srv.txt";
  }

  // mpi
  int i_rank;
  int n_rank;

  MPI_Init(&argc, &argv);
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_rank(comm, &i_rank);
  MPI_Comm_size(comm, &n_rank);

  // create intracomm
  size_t  s_code_name = strlen(code_name) + 1;
  size_t *s_recv = malloc(sizeof(size_t) * n_rank);
  int *code_ids = malloc(sizeof(int) * n_rank);
  int *code_name_idx = malloc(sizeof(int) * n_rank);
  int n_code_name = 0;
  int s_total_recv = 0;
  int *s_recv_idx = malloc(sizeof(int) * n_rank);

  MPI_Allgather(&s_code_name,
                sizeof(size_t),
                MPI_UNSIGNED_CHAR,
                s_recv,
                n_rank * sizeof(size_t),
                MPI_UNSIGNED_CHAR,
                comm);

  for (int i = 0; i < n_rank; i++) {
    s_total_recv += s_recv[i];
    if (i == 0) {
      s_recv_idx[i] = 0;
    } else {
      s_recv_idx[i] = s_recv_idx[i-1] + s_recv[i];
    }
    code_ids[i] = -1;
  }

  char *code_names = malloc(s_total_recv);
  MPI_Allgatherv(code_name,
                 s_code_name,
                 MPI_CHAR,
                 code_names,
         (int *) s_recv,
                 s_recv_idx,
                 MPI_CHAR,
                 comm);

  for (int i = 0; i < i_rank + 1; i++) {
    for (int j = 0; j < n_code_name; j++) {
      int idx = code_name_idx[j];
      if (_substrcmp(&code_names[idx], (int) s_recv[idx], &code_names[i], (int) s_recv[i]) == 1) {
        code_ids[i] = code_ids[code_name_idx[j]];
      }
    }
    if (code_ids[i] == -1) {
      code_name_idx[n_code_name++] = i;
      code_ids[i] = i;
    }
  }

  printf("%s has id %d\n", code_name, code_ids[i_rank]);
  fflush(stdout);

  MPI_Comm intra_comm;
  MPI_Comm_split(comm, code_ids[i_rank], i_rank, &intra_comm);

  // free
  free(code_name_idx);
  free(code_ids);
  free(code_names);
  free(s_recv_idx);
  free(s_recv);

  int i_intra_rank;
  int n_intra_rank;
  MPI_Comm_rank(intra_comm, &i_intra_rank);
  MPI_Comm_size(intra_comm, &n_intra_rank);

  // port choice
  MPI_Comm comm_node;
  int i_rank_node;

  // shared comm split
  MPI_Comm_split_type(intra_comm, MPI_COMM_TYPE_SHARED, i_intra_rank, MPI_INFO_NULL, &comm_node);
  MPI_Comm_rank(comm_node, &i_rank_node);

  int server_port = port_begin + i_rank_node;

  // create server
  p_server svr = malloc(sizeof(t_server));
  if (CWP_server_create(comm, server_port, 0, svr) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "Server creation failed\n");
    return -1;
  }

  // write config file
  // --> retreive host_name size
  int  irank_host_name_size     = strlen(svr->host_name);
  int *all_jrank_host_name_size = malloc(sizeof(int) * n_intra_rank);

  MPI_Allgather(&irank_host_name_size,
                1,
                MPI_INT,
                all_jrank_host_name_size,
                1,
                MPI_INT,
                intra_comm);

  int max_host_name_size = -1;
  for (int i = 0; i < n_intra_rank; i++) {
    if (all_jrank_host_name_size[i] > max_host_name_size) {
      max_host_name_size = all_jrank_host_name_size[i];
    }
  }

  free(all_jrank_host_name_size);

  // --> create format and data string
  char format[99];
  sprintf(format,"%s%d.%ds/%s9.9d\n", "%", max_host_name_size, max_host_name_size, "%");

  char data[max_host_name_size + 12];
  sprintf(data, format, svr->host_name, svr->port);

  // --> open
  PDM_io_file_t *write = NULL;
  PDM_l_num_t    ierr;

  PDM_io_open(config,
              PDM_IO_FMT_BIN,
              PDM_IO_SUFF_MAN,
              "",
              PDM_IO_BACKUP_OFF,
              PDM_IO_KIND_MPI_SIMPLE, // PDM_IO_KIND_MPIIO_EO,
              PDM_IO_MOD_WRITE,
              PDM_IO_NATIVE,
              PDM_MPI_mpi_2_pdm_mpi_comm(&intra_comm),
              -1.,
              &write,
              &ierr);

  // --> global write of header
  char  buf[99];
  sprintf(buf, "FORMAT hostname/port\nN %10.10ld\nSIZE %5.5ld\n", (long) n_intra_rank, strlen(data)); // header_size = 45 char

  size_t s_buf =  strlen(buf);
  PDM_io_global_write(write,
        (PDM_l_num_t) sizeof(char),
        (PDM_l_num_t) s_buf,
                      buf);

  // --> block write of data
  int one = 1;
  PDM_g_num_t i_rank_gnum = (PDM_g_num_t) (i_intra_rank+1);

  PDM_io_par_interlaced_write(write,
                              PDM_STRIDE_CST_INTERLACED,
              (PDM_l_num_t *) &one,
                (PDM_l_num_t) strlen(data),
                (PDM_l_num_t) one,
                              &i_rank_gnum,
               (const void *) data);

  // --> close
  PDM_io_close(write);
  PDM_io_free(write);

  // verbose
  MPI_Barrier(comm);

  if (i_rank == 0) {
    printf("----------------------------------------------------------------------------\n");
    printf("All servers listening and cwipi config file created. You may connect clients\n");
    printf("----------------------------------------------------------------------------\n");
  }

  // accept
  CWP_server_run(svr);

  // shutdown server
  CWP_server_kill(svr);

  // free
  free(svr);

  // mpi finalize
  MPI_Comm_free(&intra_comm);
  MPI_Comm_free(&comm_node);
  MPI_Finalize();

  return 0;
}


#ifdef __cplusplus
}
#endif /* __cplusplus */
