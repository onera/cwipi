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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include <pdm_error.h>
#include <pdm_io.h>
#include <pdm_mpi.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <unistd.h>
#include "pdm_logging.h"
#include "pdm_printf.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*
 * Function definitions
 */

static void
_usage(int exit_code)
{
  PDM_printf
    ("\n"
     "  Usage: \n\n"
     "  -c     Filename of the server configuration file.\n\n"
     "  -p     Begin and end port number for the server sockets port range.\n\n"
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
 int           *port_end    // end of port range
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

    else
      _usage(EXIT_FAILURE);
    i++;

  }

}


/*
 * Main
 */

int main(int argc, char *argv[])
{

  // default
  char *config     = NULL;
  int   port_begin = 1024;
  int   port_end   = 49151;

  _read_args(argc,
             argv,
             &config,
             &port_begin,
             &port_end);

  if (config == NULL) {
    config = (char *) "cwp_config_srv.txt";
  }

  // mpi
  int i_rank;
  int n_rank;

  PDM_MPI_Init(&argc, &argv);
  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  // port choice (test numa?)
  int port = port_begin + i_rank;

  char *p =malloc(sizeof(int));
  sprintf(p,"%d",port_end);
  int max_port_size = strlen(p);
  free(p);

  // retreive host_name
  char *host_name = malloc(256);
  gethostname(host_name, 256);

  // determine max host_name size
  int  irank_host_name_size     = strlen(host_name);
  int *all_jrank_host_name_size = malloc(sizeof(int) * n_rank);

  PDM_MPI_Allgather(&irank_host_name_size,
                    1,
                    PDM_MPI_INT,
                    all_jrank_host_name_size,
                    1,
                    PDM_MPI_INT,
                    comm);

  int max_host_name_size = -1;
  for (int i = 0; i < n_rank; i++) {
    if (all_jrank_host_name_size[i] > max_host_name_size) {
      max_host_name_size = all_jrank_host_name_size[i];
    }
  }

  // create string: host_name/port\n "%10.10s/5.5d\n"
  char format0[2] = "%";

  char format1[2 * sizeof(int) + 4];
  sprintf(format1,"%d.%dd\n", max_port_size, max_port_size);

  char *format2 = malloc(2 * sizeof(int) + 5);
  format2 = format0;
  strcat(format2, format1);

  char format3[2 * sizeof(int) + 4];
  sprintf(format3,"%d.%ds/",max_host_name_size, max_host_name_size);

  char format5[2] = "%";

  char *format4 = malloc(2 * sizeof(int) + 5);
  format4 = format5;

  strcat(format4, format3);

  char *format = malloc(4 * sizeof(int) + 10);
  format = format4;

  strcat(format, format2);

  log_trace("%s", format);

  char *data = malloc(max_port_size + max_host_name_size);
  sprintf(data, format, host_name, port);

  log_trace("%s", data);

  // write with pdm_io
  // --> open

  PDM_io_file_t *unite = NULL;
  PDM_l_num_t              ierr;

  PDM_io_open(config,
              PDM_IO_FMT_BIN,
              PDM_IO_SUFF_MAN,
              "",
              PDM_IO_BACKUP_OFF,
              PDM_IO_KIND_MPIIO_EO,
              PDM_IO_MOD_WRITE,
              PDM_IO_NATIVE,
              comm,
              1.,
              &unite,
              &ierr);

  // --> global write: header and offset (= host_name_size + port_size + 2)

  char  buf[82];
  sprintf(buf, "FORMAT hostname/port\nSIZE %ld\n", strlen(data));

  size_t s_buf =  strlen(buf);
  PDM_io_global_write(unite,
        (PDM_l_num_t) sizeof(char),
        (PDM_l_num_t) s_buf,
                      buf);

  // --> par_block_write

  int one = 1;
  PDM_g_num_t debut_bloc = i_rank * strlen(data) + strlen(buf);

  PDM_io_par_block_write(unite,
                         PDM_STRIDE_CST_INTERLACED,
         (PDM_l_num_t *) &one, // n_composantes
           (PDM_l_num_t) sizeof(char) * strlen(data),
           (PDM_l_num_t) one,  // n_donnees
                         debut_bloc,
                         data);

  // --> close

  PDM_io_close(unite);

  // read with pdm_io
  // --> read data size

  int header_size = 0;

  FILE *f = fopen(config, "r");

  if (f == NULL) {
    PDM_error(__FILE__, __LINE__, 0, "Could not read file %s\n", config);
  }

  char line[999];

  int size = 0;
  int not_size = 1;

  while (1) {

    int stat = fscanf(f, "%s", line);

    if (stat == EOF) {
      break;
    }

    log_trace("%s\n", line);

    if (strstr(line, "SIZE") != NULL) {
      not_size = 0;
      header_size += strlen(line);
      // Get size
      fscanf(f, "%s", line);
      // fscanf(f, "%d", &size);
      header_size += strlen(line);
      size = atoi(line);
    }

    if (not_size) {
      header_size += strlen(line);
    }

  }

  header_size += 4; // for blanc and \n

  log_trace("size = %d\n", size);

  fclose(f);

  // --> open

  PDM_io_file_t *read = NULL;
  PDM_l_num_t    ierr1;

  PDM_io_open(config,
              PDM_IO_FMT_BIN,
              PDM_IO_SUFF_MAN,
              "",
              PDM_IO_BACKUP_OFF,
              PDM_IO_KIND_MPIIO_EO,
              PDM_IO_MOD_WRITE,
              PDM_IO_NATIVE,
              comm,
              1.,
              &read,
              &ierr1);


  // --> read data (hostname/port);

  int one1 = 1;
  PDM_g_num_t debut_bloc1 = i_rank * size + header_size;

  char *read_data = NULL;

  PDM_io_par_block_read(read,
                        PDM_STRIDE_CST_INTERLACED,
         (PDM_l_num_t *) &one1, // n_composantes
           (PDM_l_num_t) sizeof(char) * strlen(data),
           (PDM_l_num_t) one1,  // n_donnees
                         debut_bloc1,
                         read_data);

  // retreive hostname and port seperatly
  char d1[] = "/";
  char *p1 = strtok(read_data, d1);
  log_trace("%s\n", p1);
  p1 = strtok(NULL, d1);
  int port1 = atoi(p1);
  log_trace("%d\n", port1);

  // --> close

  PDM_io_close(read);

  // free
  free(all_jrank_host_name_size);

  PDM_MPI_Finalize();

  return 0;
}


#ifdef __cplusplus
}
#endif /* __cplusplus */
