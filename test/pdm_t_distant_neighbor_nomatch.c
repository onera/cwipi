#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "pdm.h"
#include "pdm_priv.h"

#include "pdm_mpi.h"
#include "pdm_config.h"
#include "pdm_distant_neighbor.h"
#include "pdm_points_merge.h"
#include "pdm_printf.h"
#include "pdm_logging.h"
#include "pdm_error.h"

/**
 *
 * \brief  Main
 *
 */

int
main
(
int argc,
char *argv[]
)
{
  int iRank;
  int nRank;

  PDM_MPI_Init (&argc, &argv);
  PDM_MPI_Comm_rank (PDM_MPI_COMM_WORLD, &iRank);
  PDM_MPI_Comm_size (PDM_MPI_COMM_WORLD, &nRank);

  // int point_list_j1[15] = {145, 153, 155, 163, 172, 180, 201, 206, 207, 211, 212, 215, 217, 222, 227};
  // int point_list_j2[15] = {136, 142, 150, 151, 160, 161, 170, 171, 188, 191, 194, 199, 205, 214, 221};

  /* Connection de join1 avec join2 */
  int connect_idx_j1[4] = {0, 1, 3, 4};

  int oppRank;
  int oppPart;
  if(nRank == 1){
    oppRank = 0;
    oppPart = 1;
  } else if (nRank == 2){
    oppRank = 1;
    oppPart = 0;
  }

  int connect_triplet_j1[12] = {// Fisrt
                                oppRank, oppPart, 0,
                                // Second
                                oppRank, oppPart, 1,
                                oppRank, oppPart, 0,
                                // Third
                                oppRank, oppPart, 1};


  /* Connection de join2 avec join1 */
  // const int n_faces_j2 = 3;
  int connect_idx_j2[3] = {0, 2, 4};
  int connect_triplet_j2[12] = {// First
                                0, 0, 0,
                                0, 0, 1,
                                // Second
                                0, 0, 1,
                                0, 0, 2};

  /*
   * Data
   */
  int point_list_j1[3] = {101, 201, 301};
  int point_list_j2[2] = {102, 202};

  // Test 1 :
  // int point_list_var_j1[6] = {101, 1010,
  //                             201, 2010,
  //                             301, 3010};
  // int point_list_var_stri_j1[3] = {2, 2, 2};

  // Test 2 :
  int point_list_var_j1[6] = {101,
                              201, 2010,
                              301, 3010};
  int point_list_var_stri_j1[3] = {1, 2, 2};

  // Test 3 :
  // int point_list_var_j1[6] = {101, 1010,
  //                             201,
  //                             301, 3010, 30100};
  // int point_list_var_stri_j1[3] = {2, 1, 3};

  // Test 1
  int point_list_var_j2[6] = {102, 1020, 10200,
                              202, 2020, 20200};
  int point_list_var_stri_j2[3] = {3, 3, 3};

  // Test 2
  // int point_list_var_j2[6] = {102,
  //                             202, 2020, 20200};
  // int point_list_var_stri_j2[2] = {1, 3};

  // Test 3
  // int point_list_var_j2[6] = {102, 1020, 10200,
  //                             202};
  // int point_list_var_stri_j2[2] = {3, 1};

  /*
   * Begin
   */

  int n_cloud;
  int *n_entity;
  int **candidates_idx;
  int **candidates_desc;
  if(nRank == 1){
    n_cloud = 2;
    candidates_idx  = (int **) malloc( n_cloud * sizeof(int**));
    candidates_desc = (int **) malloc( n_cloud * sizeof(int**));
    n_entity        = (int * ) malloc( n_cloud * sizeof(int* ));

    n_entity[0] = 3;
    candidates_idx[0]  = connect_idx_j1;
    candidates_desc[0] = connect_triplet_j1;

    n_entity[1] = 2;
    candidates_idx[1]  = connect_idx_j2;
    candidates_desc[1] = connect_triplet_j2;
  } else if ( nRank == 2){
    n_cloud = 1;
    candidates_idx  = (int **) malloc( n_cloud * sizeof(int**));
    candidates_desc = (int **) malloc( n_cloud * sizeof(int**));
    n_entity        = (int * ) malloc( n_cloud * sizeof(int* ));
    if(iRank == 0){
      n_entity[0]        = 3;
      candidates_idx[0]  = connect_idx_j1;
      candidates_desc[0] = connect_triplet_j1;
    } else if (iRank == 1){
      n_entity[0]        = 2;
      candidates_idx[0]  = connect_idx_j2;
      candidates_desc[0] = connect_triplet_j2;
    }

  } else {
    PDM_error(__FILE__, __LINE__, 0, "pdm_t_distant_neighbor error : Bad number of process for test cases \n");
  }

  /*
   *  Now we have the connection between the two cloud (in this case it's a perfect match )
   *  We need to find out the connection in the partition (not in the cloud )
   *  To do that we setup an exchange with distant neighbor with the local face of each partition
   */
  int pdn_id = PDM_distant_neighbor_create(PDM_MPI_COMM_WORLD,
                                           n_cloud,
                                           n_entity,
                                           candidates_idx,
                                           candidates_desc);

  /*
   *  SetUp exchange
   */
  int** send_entity_data = (int **) malloc( n_cloud * sizeof(int**));
  if(nRank == 1){
    send_entity_data[0] = point_list_j1;
    send_entity_data[1] = point_list_j2;
  } else if(nRank == 2){
    if(iRank == 0){
      send_entity_data[0] = point_list_j1;
    } else if(iRank == 1){
      send_entity_data[0] = point_list_j2;
    }
  }

  /*
   * Constant stride test
   */
  int stride = 1;
  int** recv_entity_data = NULL;
  PDM_distant_neighbor_exch(pdn_id,
                            sizeof(int),
                            PDM_STRIDE_CST,
                            stride,
                            NULL,
                            send_entity_data,
                            NULL,
                 (int***) &recv_entity_data);

  if(1 == 1){
    for(int ipart = 0; ipart < n_cloud; ipart++){
      int *_part_neighbor_idx  = candidates_idx[ipart];
      log_trace(" recv_entity_data[%d]::", ipart);
      for(int i_entity = 0; i_entity < _part_neighbor_idx[n_entity[ipart]]; i_entity++){
        log_trace("%d ", recv_entity_data[ipart][i_entity]);
      }
      log_trace("\n");
    }
  }

  /*
   * Variable stride test
   */
  int** send_entity_var_data = (int **) malloc( n_cloud * sizeof(int**));
  int** send_entity_var_stri = (int **) malloc( n_cloud * sizeof(int**));
  if(nRank == 1){
    send_entity_var_data[0] = point_list_var_j1;
    send_entity_var_stri[0] = point_list_var_stri_j1;
    send_entity_var_data[1] = point_list_var_j2;
    send_entity_var_stri[1] = point_list_var_stri_j2;
  } else if(nRank == 2){
    if(iRank == 0){
    send_entity_var_data[0] = point_list_var_j1;
    send_entity_var_stri[0] = point_list_var_stri_j1;
    } else if(iRank == 1){
    send_entity_var_data[0] = point_list_var_j2;
    send_entity_var_stri[0] = point_list_var_stri_j2;
    }
  }

  int** recv_entity_var_stri = NULL;
  int** recv_entity_var_data = NULL;
  PDM_distant_neighbor_exch(pdn_id,
                            sizeof(int),
                            PDM_STRIDE_VAR,
                            -1,
                            send_entity_var_stri,
                            send_entity_var_data,
                 (int***) &recv_entity_var_stri,
                 (int***) &recv_entity_var_data);

  /*
   * Free
   */
  PDM_distant_neighbor_free(pdn_id);
  free(candidates_idx);
  free(candidates_desc);
  free(n_entity);
  free(send_entity_data);
  free(send_entity_var_data);
  free(send_entity_var_stri);
  for(int i_cloud = 0; i_cloud < n_cloud; i_cloud++){
    free(recv_entity_data[i_cloud]);
    free(recv_entity_var_stri[i_cloud]);
    free(recv_entity_var_data[i_cloud]);
  }
  free(recv_entity_data);
  free(recv_entity_var_stri);
  free(recv_entity_var_data);
  PDM_MPI_Finalize();

  PDM_printf ("\nfin Test\n");

  return 0;

}

