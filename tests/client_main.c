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
#include <unistd.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "client_server/client.h"
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


  // mpi
  int i_rank;
  int n_rank;

  MPI_Init(&argc, &argv);
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_rank(comm, &i_rank);
  MPI_Comm_size(comm, &n_rank);

  if (config == NULL) {
    if (i_rank == 0) {
      config = (char *) "client_main_o/code1/cwp_config_srv.txt";
    }
    else {
      config = (char *) "client_main_o/code2/cwp_config_srv.txt";
    }
  }

  assert (n_rank == 2);

  // launch server

  if (i_rank == 0) {
    system("mkdir -p client_main_o/code1");
    system("mkdir -p client_main_o/code2");
    system("rm -f ./client_main_o/code1/cwp_config_srv.txt");
    system("rm -f ./client_main_o/code2/cwp_config_srv.txt");
    system("mpiexec -n 1 ../bin/cwp_server -cn code0 -p 49100 49100 -c \"client_main_o/code1/cwp_config_srv.txt\" : -n 1  ../bin/cwp_server -cn code1 -p 49101 49101 -c \"client_main_o/code2/cwp_config_srv.txt\" &");
  }

  while (access(config, R_OK) != 0) {
    sleep(1);
  }
  sleep(5);

  // CWP_Init
  int n_code = 0;
  const char **code_names = NULL;
  double *times_init = NULL;
  CWP_Status_t *is_coupled_rank = NULL;
  int id_code = 0;

  if (i_rank == 0) {
    n_code = 1;
    id_code = 0;
    code_names = malloc(sizeof(char *) * n_code);
    code_names[0] = "code1";
    is_coupled_rank = malloc(sizeof(CWP_Status_t) * n_code);
    is_coupled_rank[0] = CWP_STATUS_ON;
  }

  if (i_rank == 1) {
    n_code = 1;
    id_code = 1;
    code_names = malloc(sizeof(char *) * n_code);
    code_names[0] = "code2";
    is_coupled_rank = malloc(sizeof(CWP_Status_t) * n_code);
    is_coupled_rank[0] = CWP_STATUS_ON;
  }

  times_init = malloc(sizeof(double) * n_code);
  for (int i = 0 ; i < n_code ; i++) {
    times_init[i] = 0;
  }

  // Outputfile
  // if (i_rank == 0) {
  //   FILE *f = fopen("output_file_code1.txt", "w");
  //   CWP_client_Output_file_set(f);
  // }

  // if (i_rank == 1) {
  //   FILE *f = fopen("output_file_code2.txt", "w");
  //   CWP_client_Output_file_set(f);
  // }

  MPI_Comm intra_comm;
  MPI_Comm_split(comm, id_code, i_rank, &intra_comm); // smallest i_rank becomes rank 0 in intra_comm

  CWP_client_Init(intra_comm,
                  config,
                  n_code,
                  (const char **) code_names,
                  is_coupled_rank,
                  times_init);

  // // CWP_User_structure_*
  // CWP_client_User_structure_set("code1", NULL);
  // CWP_client_User_structure_get("code1");

  // // CWP_Codes_*
  // int n_codes = CWP_client_Codes_nb_get();
  // printf("n_codes = %d\n", n_codes);
  // char **codeNames = NULL;
  // codeNames = (char **) CWP_client_Codes_list_get();
  // for (int i = 0; i < n_codes; i++) {
  //   printf("i_rank: %d code_names[i] = %s\n", i_rank, codeNames[i]);
  // }

  // // CWP_Loc_codes_*
  // int n_Loc_codes = CWP_client_Loc_codes_nb_get();
  // printf("n_Loc_codes = %d\n", n_Loc_codes);
  // char **LoccodeNames = NULL;
  // LoccodeNames = (char **) CWP_client_Loc_codes_list_get();
  // for (int i = 0; i < n_Loc_codes; i++) {
  //   printf("i_rank: %d Loc_code_names[i] = %s\n", i_rank, LoccodeNames[i]);
  // }

  // // Properties_dump
  // CWP_client_Properties_dump();

  // // State_update
  // if (i_rank == 0) {
  //   CWP_client_State_update("code1", CWP_STATE_IN_PROGRESS);
  //   CWP_State_t state = CWP_client_State_get("code1");
  //   printf("state == CWP_STATE_IN_PROGRESS: %d\n", state == CWP_STATE_IN_PROGRESS);
  // }

  // // CWP_Param_*
  // if (i_rank == 0) {
  //   int toto = 42;
  //   CWP_client_Param_lock("code1");
  //   CWP_client_Param_add("code1", "toto", CWP_INT, &toto);
  //   double tata = 0.99;
  //   CWP_client_Param_add("code1", "tata", CWP_DOUBLE, &tata);
  //   double tota = 0.55;
  //   CWP_client_Param_set("code1", "tata", CWP_DOUBLE, &tota);
  //   CWP_client_Param_unlock("code1");
  // }

  // if (i_rank == 1) {
  //   CWP_client_Param_lock("code2");
  //   const char *A = "Bonjour code 1 !";
  //   CWP_client_Param_add("code2", "toto2", CWP_CHAR, &A);
  //   CWP_client_Param_unlock("code2");
  // }

  // MPI_Barrier(comm);

  // double titi1;
  // CWP_client_Param_get("code1", "tata", CWP_DOUBLE, &titi1);
  // printf("i_rank: %d code 1 : tata : %f\n", i_rank, titi1);
  // const char *titi2 = NULL;
  // CWP_client_Param_get("code2", "toto2", CWP_CHAR, &titi2);
  // printf("i_rank: %d code 2 : toto2 : %s\n", i_rank, titi2);
  // int titi;
  // CWP_client_Param_get("code1", "toto", CWP_INT, &titi);
  // printf("i_rank: %d code 1 : toto : %d\n", i_rank, titi);

  // int code1_n_double = -1;

  // if (i_rank == 0) {
  //   int code1_n_int = CWP_client_Param_n_get("code1", CWP_INT);
  //   printf("i_rank: %d code 1 : n_int_param : %d\n", i_rank, code1_n_int);
  //   code1_n_double = CWP_client_Param_n_get("code1", CWP_DOUBLE);
  //   printf("i_rank: %d code 1 : n_double_param : %d\n", i_rank, code1_n_double);
  // }

  // if (i_rank == 0) {
  //   double tatic = 107.52;
  //   CWP_client_Param_lock("code1");
  //   CWP_client_Param_add("code1", "tatic", CWP_DOUBLE, &tatic);
  //   CWP_client_Param_unlock("code1");
  // }

  // if (i_rank == 1) {
  //   double totic = 33.50;
  //   CWP_client_Param_lock("code2");
  //   CWP_client_Param_add("code2", "tatic", CWP_DOUBLE, &totic);
  //   CWP_client_Param_unlock("code2");
  // }

  // if (i_rank == 0) {
  //   CWP_client_Param_lock("code1");
  //   CWP_client_Param_del("code1", "toto", CWP_INT);
  //   CWP_client_Param_unlock("code1");
  // }

  // MPI_Barrier(comm);

  // double tita;
  // CWP_client_Param_get("code1", "tatic", CWP_DOUBLE, &tita);
  // printf("i_rank: %d code 1 : tatic : %f\n", i_rank, tita);
  // double tito;
  // CWP_client_Param_get("code2", "tatic", CWP_DOUBLE, &tito);
  // printf("i_rank: %d code 2 : tatic : %f\n", i_rank, tito);

  // // reduce
  // double res = 0;
  // char **g_code_names = malloc(sizeof(char *) * 2);
  // g_code_names[0] = malloc(sizeof(char) * 99);
  // g_code_names[1] = malloc(sizeof(char) * 99);
  // char test1[6] = "code1";
  // strcpy(g_code_names[0], test1);
  // char test2[6] = "code2";
  // strcpy(g_code_names[1], test2);
  // g_code_names[0][6] = '\0';
  // g_code_names[1][6] = '\0';

  // CWP_client_Param_reduce(CWP_OP_MAX, "tatic", CWP_DOUBLE, &res, 2, (const char **) g_code_names);

  // free(g_code_names[0]);
  // free(g_code_names[1]);
  // free(g_code_names);

  // printf("res = %f\n", res);

  // if (i_rank == 0) {
  //   char **param_names = NULL;
  //   CWP_client_Param_list_get("code1", CWP_DOUBLE, &code1_n_double, &param_names);

  //   for (int i = 0; i < code1_n_double; i++) {
  //     printf("i_rank: %d code 1 : param[%d] = %s\n", i_rank, i, param_names[i]);
  //   }

  // }

  // if (i_rank == 1) {
  //   int code2_n_char = -1;
  //   char **param_names = NULL;
  //   CWP_client_Param_list_get("code2", CWP_CHAR, &code2_n_char, &param_names);

  //   for (int i = 0; i < code2_n_char; i++) {
  //     printf("i_rank: %d code 1 : param[%d] = %s\n", i_rank, i, param_names[i]);
  //   }

  // }

  // if (i_rank == 1) {
  //   int is_param = CWP_client_Param_is("code2", "toto2", CWP_CHAR);
  //   printf("i_rank: %d code 2 : toto2 is param = %d\n", i_rank, is_param);

  //   is_param = CWP_client_Param_is("code2", "tambouille", CWP_CHAR);
  //   printf("i_rank: %d code 2 : tambouille is param = %d\n", i_rank, is_param);
  // }

  // MPI_Barrier(comm);

  char cpl_id1[] = "cpl1_code1_code2";
  CWP_Spatial_interp_t interp_method = CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE;

  if (i_rank == 0) {
    printf("CWP_client_Cpl_create rank 0\n");
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
    printf("CWP_client_Cpl_create rank 1\n");
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

  MPI_Barrier(comm);

  // // GLOBAL DATA exchange
  // const char *global_data_name  = "andouillette";
  // size_t  s_entity = sizeof(double);
  // int     stride   = 2;
  // int     n_entity = 2;

  // double *send_data = malloc(s_entity * stride * n_entity);
  // send_data[0] = 0.1;
  // send_data[1] = 0.2;
  // send_data[2] = 0.3;
  // send_data[3] = 0.4;

  // if (i_rank == 0) {
  //   CWP_client_Global_data_issend("code1",
  //                                 cpl_id1,
  //                                 global_data_name,
  //                                 s_entity,
  //                                 stride,
  //                                 n_entity,
  //                                 (void *) send_data);
  // }

  // double *recv_data = malloc(s_entity * stride * n_entity);
  // if (i_rank == 1) {
  //   CWP_client_Global_data_irecv("code2",
  //                                cpl_id1,
  //                                global_data_name,
  //                                s_entity,
  //                                stride,
  //                                n_entity,
  //                      (void **) &recv_data);
  // }

  // MPI_Barrier(comm);

  // if (i_rank == 0) {
  //   CWP_client_Global_data_wait_issend("code1",
  //                                      cpl_id1,
  //                                      global_data_name);
  // }

  // if (i_rank == 1) {
  //   CWP_client_Global_data_wait_irecv("code2",
  //                                     cpl_id1,
  //                                     global_data_name);
  // }

  // MPI_Barrier(comm);

  // for (int i = 0; i < 4; i++) {
  //   if (i_rank == 0) {
  //     printf("rank %d -- send[%d] : %f\n", i_rank, i, send_data[i]);
  //     fflush(stdout);
  //   } else {
  //     printf("rank %d -- recv[%d] : %f\n", i_rank, i, recv_data[i]);
  //     fflush(stdout);
  //   }
  // }

  // PART DATA exchange
  const char *part_data_name = "Fifi Brindacier";

  int n_part = 1;
  int n_elt = 12;
  int *n_elts = malloc(sizeof(int *) * n_part);
  for (int i_part = 0; i_part < n_part; i_part++) {
    n_elts[i_part] = n_elt;
  }
  CWP_g_num_t **gnum_elt = malloc(sizeof(CWP_g_num_t *) * n_part);
  for (int i_part = 0; i_part < n_part; i_part++) {
    gnum_elt[i_part] = malloc(sizeof(CWP_g_num_t) * n_elt);
  }
  for (int i_part = 0; i_part < n_part; i_part++) {
    for (int i = 0; i < n_elt; i++) {
      gnum_elt[i_part][i] = i + 1;
    }
  }

  if (i_rank == 0) {
    CWP_client_Part_data_create("code1",
                                cpl_id1,
                                part_data_name,
                                CWP_PARTDATA_SEND,
                                gnum_elt,
                                n_elts,
                                n_part);
  }

  if (i_rank == 1) {
    CWP_client_Part_data_create("code2",
                                cpl_id1,
                                part_data_name,
                                CWP_PARTDATA_RECV,
                                gnum_elt,
                                n_elts,
                                n_part);
  }

  MPI_Barrier(comm);

  double **part1_to_part2_data = NULL;
  int send_request = -1;
  double **part2_data = NULL;
  int recv_request = -1;
  int n_comp = 3;

  if (i_rank == 0) {

    part1_to_part2_data = malloc(sizeof(double *) * n_part);
    for (int i_part = 0; i_part < n_part; i_part++) {
      part1_to_part2_data[i_part] = malloc(sizeof(double) * n_elt * n_comp);
    }

    for (int i_part = 0; i_part < n_part; i_part++) {
      for (int i = 0; i < n_elt; i++) {
        for (int i_comp = 0; i_comp < n_comp; i_comp++) {
          part1_to_part2_data[i_part][3*i + i_comp] = 0.1*i + i_comp;
        }
      }
    }
    CWP_client_Part_data_issend("code1",
                                cpl_id1,
                                part_data_name,
                                sizeof(double),
                                n_comp,
                                (void **) part1_to_part2_data,
                                &send_request);
  }

  if (i_rank == 1) {

    part2_data = malloc(sizeof(double *) * n_part);
    for (int i_part = 0; i_part < n_part; i_part++) {
      part2_data[i_part] = malloc(sizeof(double) * n_elt * n_comp);
    }
    CWP_client_Part_data_irecv("code2",
                               cpl_id1,
                               part_data_name,
                               sizeof(double),
                               n_comp,
                               (void ***) &part2_data,
                               &recv_request);
  }

  MPI_Barrier(comm);

  if (i_rank == 0) {

    CWP_client_Part_data_wait_issend("code1",
                                     cpl_id1,
                                     part_data_name,
                                     send_request);
  }

  if (i_rank == 1) {

    CWP_client_Part_data_wait_irecv("code2",
                                    cpl_id1,
                                    part_data_name,
                                    recv_request);
  }

  MPI_Barrier(comm);

  // --> check
  for (int i_part = 0; i_part < n_part; i_part++) {
    for (int i = 0; i < n_elt; i++) {
      for (int i_comp = 0; i_comp < n_comp; i_comp++) {
        if (i_rank == 0) {
          printf("%d - %ld -> s[%d][%d][%d] : %f\n", i_rank, gnum_elt[i_part][i], i_part, i, i_comp, part1_to_part2_data[i_part][3*i + i_comp]);
        }
        else if (i_rank == 1) {
          printf("%d - %ld -> r[%d][%d][%d] : %f\n", i_rank, gnum_elt[i_part][i], i_part, i, i_comp, part2_data[i_part][3*i + i_comp]);
        }
      }
    }
  }

  MPI_Barrier(comm);

  if (i_rank == 0) {
    CWP_client_Part_data_del("code1",
                             cpl_id1,
                             part_data_name,
                             CWP_PARTDATA_SEND);
  }

  if (i_rank == 1) {
    CWP_client_Part_data_del("code2",
                             cpl_id1,
                             part_data_name,
                             CWP_PARTDATA_RECV);
  }

  if (i_rank == 0) {
    printf("CWP_client_Cpl_del rank 0\n");
    CWP_client_Cpl_del("code1", cpl_id1);
  }

  if (i_rank == 1) {
    printf("CWP_client_Cpl_del rank 1\n");
    CWP_client_Cpl_del("code2", cpl_id1);
  }

  // CWP_Finalize
  printf("CWP_client_Finalize\n");
  CWP_client_Finalize();

  // free
  printf("Free\n");
  if (times_init != NULL) free(times_init);
  if (code_names != NULL) free(code_names);
  if (is_coupled_rank != NULL) free(is_coupled_rank);

  printf("MPI_Finalize\n");

  MPI_Comm_free(&intra_comm);

  MPI_Finalize();

  exit(0);
  return 0;
}


#ifdef __cplusplus
}
#endif /* __cplusplus */
