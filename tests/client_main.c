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

#include "client.h"
#include <pdm_error.h>
#include <pdm_io.h>
#include <pdm_mpi.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include "pdm_logging.h"
#include "pdm_printf.h"
#include "cwp.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

#define CWP_HEADER_SIZE    32

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
     "  -h     This message.\n\n");

  exit(exit_code);
}

static void
_read_args
(
 int            argc,
 char         **argv,
 char         **config  // filename for server ip adresses + ports
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

    else
      _usage(EXIT_FAILURE);
    i++;

  }

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

  _read_args(argc,
             argv,
             &config);

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
              comm,
              -1.,
              &read,
              &ierr);

  // --> global read offset in header
  char *buffer = malloc(CWP_HEADER_SIZE+1);
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

  char div_word[] = " ";
  char *word = strtok(line, div_word);
  word = strtok(NULL, div_word);

  int offset = atoi(word);

  // --> block read hostname/port
  char *data = malloc(offset+1);

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
  char *server_name = malloc(strlen(str)+1);
  memcpy(server_name, str, strlen(str)+1);
  str = strtok(NULL, div);
  int server_port = atoi(str);

  // --> close
  PDM_io_close(read);
  PDM_io_free(read);

  // connect
  if (CWP_client_connect(server_name, server_port, CWP_CLIENTFLAG_VERBOSE) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "Client connexion failed\n");
    return -1;
  }

  // CWP_Init
  int n_code = 0;
  const char **code_names = NULL;
  double *times_init = NULL;
  CWP_Status_t *is_coupled_rank = NULL;

  if (i_rank == 0) {
    n_code = 1;
    code_names = malloc(sizeof(char *) * n_code);
    code_names[0] = "code1";
    is_coupled_rank = malloc(sizeof(CWP_Status_t) * n_code);
    is_coupled_rank[0] = CWP_STATUS_ON;
  }

  if (i_rank == 1) {
    n_code = 1;
    code_names = malloc(sizeof(char *) * n_code);
    code_names[0] = "code2";
    is_coupled_rank = malloc(sizeof(CWP_Status_t) * n_code);
    is_coupled_rank[0] = CWP_STATUS_ON;
  }

  times_init = malloc(sizeof(double) * n_code);
  for (int i = 0 ; i < n_code ; i++) {
    times_init[i] = 0;
  }

  CWP_client_Init(n_code,
                  (const char **) code_names,
                  is_coupled_rank,
                  times_init);

  // CWP_User_structure_*
  CWP_client_User_structure_set("code1", NULL);
  CWP_client_User_structure_get("code1");

  // Outputfile
  CWP_client_Output_file_set("outputfile.txt");

  // CWP_Codes_*
  int n_codes = CWP_client_Codes_nb_get();
  printf("n_codes = %d\n", n_codes);
  char **codeNames = malloc(sizeof(char *));
  for (int i = 0; i < n_codes; i++) {
    codeNames[i] = malloc(99);
  }
  codeNames = CWP_client_Codes_list_get();
  for (int i = 0; i < n_codes; i++) {
    printf("i_rank: %d code_names[i] = %s\n", i_rank, codeNames[i]);
  }

  // CWP_Loc_codes_*
  int n_Loc_codes = CWP_client_Loc_codes_nb_get();
  printf("n_Loc_codes = %d\n", n_Loc_codes);
  char **LoccodeNames = malloc(sizeof(char *));
  for (int i = 0; i < n_Loc_codes; i++) {
    LoccodeNames[i] = malloc(99);
  }
  LoccodeNames = CWP_client_Loc_codes_list_get();
  for (int i = 0; i < n_Loc_codes; i++) {
    printf("i_rank: %d Loc_code_names[i] = %s\n", i_rank, LoccodeNames[i]);
  }

  // Properties_dump
  CWP_client_Properties_dump();

  // *_update
  if (i_rank == 0) {
    CWP_client_State_update("code1", CWP_INTERFACE_LINEAR);
    CWP_State_t state = CWP_client_State_get("code1");
    printf("state = %d\n", state);
    CWP_client_Time_update("code1", 0.1);
  }

  // CWP_Param_*
  if (i_rank == 0) {
    int toto = 42;
    CWP_client_Param_lock("code1");
    CWP_client_Param_add("code1", "toto", CWP_INT, &toto);
    CWP_client_Param_unlock("code1");
    double tata = 0.99;
    CWP_client_Param_lock("code1");
    CWP_client_Param_add("code1", "tata", CWP_DOUBLE, &tata);
    CWP_client_Param_unlock("code1");
    double tota = 0.55;
    CWP_client_Param_lock("code1");
    CWP_client_Param_set("code1", "tata", CWP_DOUBLE, &tota);
    CWP_client_Param_unlock("code1");
  }

  if (i_rank == 1) {
    CWP_client_Param_lock("code2");
    const char *A = "Bonjour code 1 !";
    CWP_client_Param_add("code2", "toto2", CWP_CHAR, &A);
    CWP_client_Param_unlock("code2");
  }

  PDM_MPI_Barrier(comm);

  double titi1;
  CWP_client_Param_get("code1", "tata", CWP_DOUBLE, &titi1);
  printf("i_rank: %d code 1 : tata : %f\n", i_rank, titi1);
  const char *titi2;
  CWP_client_Param_get("code2", "toto2", CWP_CHAR, &titi2);
  printf("i_rank: %d code 2 : toto2 : %s\n", i_rank, titi2);
  int titi;
  CWP_client_Param_get("code1", "toto", CWP_INT, &titi);
  printf("i_rank: %d code 1 : toto : %d\n", i_rank, titi);

  if (i_rank == 0) {
    CWP_client_Param_lock("code1");
    CWP_client_Param_del("code1", "toto", CWP_INT);
    CWP_client_Param_unlock("code1");
  }

  PDM_MPI_Barrier(comm);

  int code1_n_int = CWP_client_Param_n_get("code1", CWP_INT);
  printf("i_rank: %d code 1 : n_int_param : %d\n", i_rank, code1_n_int);
  int code1_n_double = CWP_client_Param_n_get("code1", CWP_DOUBLE);
  printf("i_rank: %d code 1 : n_double_param : %d\n", i_rank, code1_n_double);


  if (i_rank == 0) {
    double tatic = 107.52;
    CWP_client_Param_lock("code1");
    CWP_client_Param_add("code1", "tatic", CWP_DOUBLE, &tatic);
    CWP_client_Param_unlock("code1");
  }

  if (i_rank == 1) {
    double totic = 33.50;
    CWP_client_Param_lock("code2");
    CWP_client_Param_add("code2", "tatic", CWP_DOUBLE, &totic);
    CWP_client_Param_unlock("code2");
  }

  PDM_MPI_Barrier(comm);

  // reduce
  double res = 0;
  CWP_client_Param_reduce(CWP_OP_MAX, "tatic", CWP_DOUBLE, &res, 2);

  printf("i_rank: %d max(tatic) = res = %f\n", i_rank, res);

  char **param_names = NULL;

  CWP_client_Param_list_get("code1", CWP_DOUBLE, &code1_n_double, &param_names);
  for (int i = 0; i < code1_n_double; i++) {
    printf("i_rank: %d code 1 : param[%d] = %s\n", i_rank, i, param_names[i]);
  }

  int is_param = CWP_client_Param_is("code2", "toto2", CWP_CHAR);
  printf("i_rank: %d code 2 : toto2 is param = %d\n", i_rank, is_param);

  is_param = CWP_client_Param_is("code2", "tambouille", CWP_CHAR);
  printf("i_rank: %d code 2 : tambouille is param = %d\n", i_rank, is_param);

  char cpl_id1[] = "cpl1_code1_code2";
  CWP_Spatial_interp_t interp_method = CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE;

  if (i_rank == 0) {
    CWP_client_Cpl_create("code1",
                          cpl_id1,
                          "code2",
                          CWP_INTERFACE_VOLUME,
                          CWP_COMM_PAR_WITH_PART,
                          interp_method,
                          1,
                          CWP_DYNAMIC_MESH_STATIC,
                          CWP_TIME_EXCH_USER_CONTROLLED);
  }

  if (i_rank == 1) {
    CWP_client_Cpl_create("code2",
                          cpl_id1,
                          "code1",
                          CWP_INTERFACE_VOLUME,
                          CWP_COMM_PAR_WITH_PART,
                          interp_method,
                          1,
                          CWP_DYNAMIC_MESH_STATIC,
                          CWP_TIME_EXCH_USER_CONTROLLED);
  }

  if (i_rank == 0) {
    CWP_client_Cpl_del("code1", cpl_id1);
  }

  if (i_rank == 1) {
    CWP_client_Cpl_del("code2", cpl_id1);
  }

  // CWP_Finalize
  CWP_client_Finalize();

  // disconnect
  CWP_client_disconnect();

  PDM_MPI_Finalize();

  return 0;
}


#ifdef __cplusplus
}
#endif /* __cplusplus */
