#include <mpi.h>
#include <stdio.h>
#include <bft_mem.h>
#include <bft_printf.h>
#include "couplings.h"

int main 
( 
 int    argc,    /* Nombre d'arguments dans la ligne de commandes */
 char  *argv[]   /* Tableau des arguments de la ligne de commandes */
 )
{
  MPI_Init(&argc, &argv);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  FILE* outputFile = NULL;

  char* fileName = NULL;
  BFT_MALLOC(fileName, 16, char); 
  BFT_MALLOC(outputFile, 1, FILE);
  sprintf(fileName,"listing_code_C_%4.4d",rank); 

  outputFile = fopen(fileName,"w");

  MPI_Comm localComm;

  /* Initialisation
   * -------------- */

  couplings_init(MPI_COMM_WORLD, 
                 outputFile,
                 "CodeC",
                 &localComm);
  
  int currentRank;
  MPI_Comm_rank(localComm, &currentRank);


  bft_printf("\nDump apres initialisation\n");
  bft_printf("-------------------------\n");
  couplings_dump_application_properties();

  /* Test sur la transmission des parametres de controle
   * --------------------------------------------------- */

  /* Ajout de parametres */
  couplings_add_local_int_control_parameter("niterC", 0);
  couplings_add_local_double_control_parameter("physicalTimeC", 0.1);

  bft_printf("\nDump apres ajout de parametres\n");
  bft_printf("------------------------------\n");
  couplings_dump_application_properties();

  /* Modification des parametres */
  couplings_set_local_double_control_parameter("physicalTimeC", 0.2);
  couplings_set_local_int_control_parameter("niterC", 2);

  bft_printf("\nDump apres modification des parametres\n");
  bft_printf("--------------------------------------\n");
  couplings_dump_application_properties();

  /* Mise a jour des parametres de controle avec une autre application */
  couplings_synchronise_control_parameter("CodeFortran"); 
  bft_printf("\nDump apres synchronisation des parametres\n");
  bft_printf("-----------------------------------------\n");
  couplings_dump_application_properties();

  /* Recuperation des parametres */
  double physicalTime = couplings_get_local_double_control_parameter("physicalTimeC");
  int niter = couplings_get_local_int_control_parameter("niterC");

  bft_printf("\nAffichage apres recuperation des parametres\n");
  bft_printf("-------------------------------------------\n");
  bft_printf("get physicalTime : %12.5e\n", physicalTime);
  bft_printf("get niter : %i\n", niter);

  /* Supression des parametres */
  bft_printf("\nDump apres suppresion des parametres\n");
  bft_printf("------------------------------------\n");
  couplings_delete_local_double_control_parameter("physicalTimeC");
  couplings_delete_local_int_control_parameter("niterC");
  couplings_dump_application_properties();
  
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
  return 0;
}
