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

  char p[sizeof(int)];
  sprintf(p,"%d",port_end);
  int max_port_size = strlen(p);

  // retreive host_name
  char host_name[256];
  gethostname(host_name, sizeof(host_name));

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
  char format0[4 * sizeof(int) + 6];
  sprintf(format0,"%d.%ds/%d.%dd\n",max_host_name_size, max_host_name_size, max_port_size, max_port_size);

  char format[4 * sizeof(int) + 7];
  format[] = "%";
  strcat(format, format0); // TO DO: does not work change !

  log_trace("format: %s\n", format);

  // write with pdm_io
  // --> open

  // --> global write: header and offset (= host_name_size + port_size + 2)

  // --> par_write

  // --> close

  // read with pdm_io

  // PDM_MPI_File f;
  // PDM_MPI_File_open(comm,
  //                   config,
  //                   PDM_MPI_MODE_WRONLY_CREATE,
  //                   &f);

  // PDM_MPI_Offset offset = (PDM_MPI_Offset) i_rank * 300;
  // int n_written_octets;
  // PDM_MPI_File_write_at_all(f,
  //                           offset,
  //                           data,
  //                           sizeof(char) * 300, // number of written octets
  //                           PDM_MPI_BYTE,
  //                           &n_written_octets);

  // if (n_written_octets != sizeof(char) * 300) {
  //   PDM_error(__FILE__, __LINE__, 0, "Incorrect number of octets written (%d/6)\n", n_written_octets);
  // }

  // PDM_MPI_File_close(&f);

  // // read i_rank-th couple ip adress + port

  // FILE *fl = fopen(config, "r");

  // if (fl == NULL) {
  //   PDM_error(__FILE__, __LINE__, 0, "Could not read file %s.txt\n", config);
  // }

  // int i    = 0;

  // char line[999];

  // while (1) {

  //   int stat = fscanf(fl, "%s", line);

  //   if (stat == EOF) {
  //     break;
  //   }

  //   if (i == i_rank) { // ntohs ?
  //     log_trace("line: %s\n", line);
  //     break;
  //   }

  //   i++;
  // }

  // fclose(fl);

  PDM_MPI_Finalize();

  return 0;
}


#ifdef __cplusplus
}
#endif /* __cplusplus */
