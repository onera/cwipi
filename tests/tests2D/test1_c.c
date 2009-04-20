#include <stdio.h>
#include <stddef.h>
#include <string.h>

#include <mpi.h>

// Bft pour allocation de la memoire et l'impression

#include <bft_mem.h>
#include <bft_error.h>
#include <bft_printf.h>

// Pour la construction du maillage

#include "creeMaillagePolygone2D.h"

// Pour le couplage

#include "couplings.h"

//
// Fonction d'interpolation bidon, juste pour voir si c'est bien pris
// en compte 

static void _dumpStatus(couplings_exchange_status_t status )
{
  switch(status) {
  case COUPLINGS_EXCHANGE_OK :
    bft_printf("Exchange Ok\n");
    break;
  case COUPLINGS_EXCHANGE_BAD_RECEIVING :
    bft_printf("Bad receiving\n");
    break;
  default :
    bft_error(__FILE__, __LINE__, 0,"bad exchange status\n");
  }
}

static void _dumpNotLocatedPoints(const char *coupling_id,
                                  const int nNotLocatedPoints)
{
  if ( nNotLocatedPoints > 0) {
    bft_printf("Points non localises :\n");
    const int* notLocatedPoints = couplings_get_not_located_points(coupling_id);
    for(int i = 0; i < nNotLocatedPoints; i++)
      bft_printf("%i ", notLocatedPoints[i]);
    bft_printf("\n");
  }
}

static void _interpolationBidon(const int entities_dim,
                                const int n_local_vertex,
                                const int n_local_element,
                                const int n_local_polhyedra,
                                const int n_distant_point,
                                const double local_coordinates[],
                                const int local_connectivity_index[],
                                const int local_connectivity[],
                                const int local_polyhedra_face_index[],
                                const int local_polyhedra_cell_to_face_connectivity[],
                                const int local_polyhedra_face_connectivity_index[],
                                const int local_polyhedra_face_connectivity[],
                                const double distant_points_coordinates[],
                                const int distant_points_location[],
                                const int distant_points_barycentric_coordinates_index[],
                                const double distant_points_barycentric_coordinates[],
                                const int stride,
                                const couplings_solver_type_t  solver_type,
                                const void *local_field,
                                void *distant_field)
{
  for (int i = 0; i < n_distant_point; i++)
    ((double *) distant_field)[i] = (double) i;
}



int main 
( 
 int    argc,    /* Nombre d'arguments dans la ligne de commandes */
 char  *argv[]   /* Tableau des arguments de la ligne de commandes */
 )
{

  FILE *outputFile;

  bft_mem_init("logmem.txt");

  MPI_Init(&argc, &argv);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);


  char* fileName = NULL;
  BFT_MALLOC(fileName, 21, char); 
  sprintf(fileName,"listing_code_C_%4.4d",rank); 

  outputFile = fopen(fileName,"w");
  BFT_FREE(fileName);

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
  couplings_synchronize_control_parameter("CodeFortran"); 
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

  /* -----------------------
   * Test couplage P1 <-> P1
   * ----------------------- */
    
  {

    /* Initialisation du couplage */
    bft_printf("Test 1 :  Test couplage P1 <-> P1\n");
    bft_printf("\n");
    
    couplings_create_coupling("test2D_1",         // Nom du couplage
                              "CodeFortran",                      // Code couplé
                              2,                            // Dimension des entités géométriques
                              0.1,                          // Tolérance géométrique
                              COUPLINGS_STATIC_MESH,        // Maillage statique
                              COUPLINGS_SOLVER_CELL_VERTEX, // Type de champs
                              1,                            // Frequence des post-traitement
                              "EnSight Gold",               // Format du post-traitement
                              "text");                      // Options de post-traitements
    
    /* Construction du maillage local (Decoupage par Metis si plusieurs procs) */
    
    int nVertex = 0;               // Nombre de sommets
    double *coords = NULL;         // Coordonnees des sommets
    int nElts = 0;                 // Nombre d'elements
    int *eltsConnecPointer = NULL; // Index par element dans la connectivite
    int *eltsConnec = NULL;        // Description de la connectivite
    
    const double xmin = -100;
    const double xmax =  100;
    const double ymin = -100;
    const double ymax =  100;
    const int    nx   = 16;
    const int    ny   = 16;
    const int   order = 1;
    
    creeMaillagePolygone2D(order,
                           localComm,
                           xmin,
                           xmax,
                           ymin,
                           ymax,
                           1,
                           nx,
                           ny,
                           &nVertex,
                           &coords,
                           &nElts,
                           &eltsConnecPointer,
                           &eltsConnec);
    
    bft_printf("   nombre de sommets : %i\n", nVertex);
    bft_printf("   nombre d'elements : %i\n", nElts);
    
    couplings_define_mesh("test2D_1",
                          nVertex,
                          nElts,
                          coords,
                          eltsConnecPointer,
                          eltsConnec);
    
    /* Envoi de la coordonnee X
       Reception de la coordonnee Y*/
    
    double* values = NULL;
    BFT_MALLOC(values, nVertex, double);
    for (int i = 0; i < nVertex; i++)
      values[i] = coords[3*i];
    
    double* localValues = NULL;
    BFT_MALLOC(localValues, nVertex, double);

    int nNotLocatedPoints;
    
    couplings_exchange_status_t status = couplings_exchange("test2D_1",
                                                            "echange1",
                                                            1,
                                                            1,     // n_step
                                                            0.1,   // physical_time
                                                            "cooX",
                                                            values,
                                                            "cooY",
                                                            localValues,
                                                            &nNotLocatedPoints);
    _dumpStatus(status);
    _dumpNotLocatedPoints("test2D_1", nNotLocatedPoints);
    /* Suppression de l'objet de couplage */
    
    couplings_delete_coupling("test2D_1");
    
    /* Liberation de la memoire */
    
    if (coords != NULL)
      BFT_FREE(coords);
    
    if (eltsConnec != NULL)
      BFT_FREE(eltsConnec);
    
    if (eltsConnecPointer != NULL)
      BFT_FREE(eltsConnecPointer);
    
    if (values != NULL)
      BFT_FREE(values);
    
    if (localValues != NULL)
      BFT_FREE(localValues);
    bft_printf("--------------------------------------------------------\n");
  }

  /* ------------------------------------
   * Test couplage P1 -> P0 puis P0 -> P1
   * ------------------------------------ */

  {   
    /* Initialisation du couplage */
    
    bft_printf("Test 2 : Test couplage P1 -> P0 puis P0 -> P1\n");
    bft_printf("\n");
    couplings_create_coupling("test2D_2",         // Nom du couplage
                              "CodeFortran",                      // Code couplé
                              2,                            // Dimension des entités géométriques
                              0.1,                          // Tolérance géométrique 
                              COUPLINGS_STATIC_MESH,        // Maillage statique
                              COUPLINGS_SOLVER_CELL_CENTER, // Type de champs
                              1,                            // Frequence des post-traitement
                              "EnSight Gold",               // Format du post-traitement
                              "text");                      // Options de post-traitements
    
    /* Construction du maillage local (Decoupage par Metis si plusieurs procs) */
    
    int nVertex = 0;               // Nombre de sommets
    double *coords = NULL;         // Coordonnees des sommets
    int nElts = 0;                 // Nombre d'elements
    int *eltsConnecPointer = NULL; // Index par element dans la connectivite
    int *eltsConnec = NULL;        // Description de la connectivite
    
    const double xmin = -100;
    const double xmax =  100;
    const double ymin = -100;
    const double ymax =  100;
    const int    nx   = 16;    
    const int    ny   = 16;    
    const int   order = 1;
    
    creeMaillagePolygone2D(order,
                           localComm,
                           xmin, 
                           xmax, 
                           ymin, 
                           ymax,
                           1,
                           nx,
                           ny,
                           &nVertex,
                           &coords,
                           &nElts,
                           &eltsConnecPointer,
                           &eltsConnec);
    
    bft_printf("   nombre de sommets : %i\n", nVertex);
    bft_printf("   nombre d'elements : %i\n", nElts);
    
    couplings_define_mesh("test2D_2", 
                          nVertex,
                          nElts,
                          coords,
                          eltsConnecPointer,
                          eltsConnec);
    
    /* Reception de la coordonnee Y
       Envoi de la coordonnee Y */
    
    double* localValues = NULL;
    BFT_MALLOC(localValues, nElts, double);

    /* Receive */
    
    int nNotLocatedPoints;
    
    couplings_exchange_status_t status = couplings_exchange("test2D_2",
                                                            "echange1",
                                                            1,     // stride
                                                            1,     // n_step
                                                            0.1,   // physical_time
                                                            NULL,
                                                            NULL,
                                                            "cooY",
                                                            localValues,
                                                            &nNotLocatedPoints);
    bft_printf("Send\n");
    _dumpStatus(status);
    _dumpNotLocatedPoints("test2D_2", nNotLocatedPoints);

    /* Send */

    status = couplings_exchange("test2D_2",
                                "echange2",
                                1, // stride
                                1,     // n_step
                                0.1,   // physical_time
                                "cooYY",
                                localValues,
                                NULL,
                                NULL,
                                &nNotLocatedPoints);
    bft_printf("Receive\n");
    _dumpStatus(status);
    _dumpNotLocatedPoints("test2D_2", nNotLocatedPoints);

    /* Suppression de l'objet de couplage */
    
    couplings_delete_coupling("test2D_2");
    
    /* Liberation de la memoire */
    
    if (coords != NULL)
      BFT_FREE(coords);
    
    if (eltsConnec != NULL)
      BFT_FREE(eltsConnec);
    
    if (eltsConnecPointer != NULL)
      BFT_FREE(eltsConnecPointer);
    
    if (localValues != NULL)
      BFT_FREE(localValues);
    bft_printf("--------------------------------------------------------\n");
    
  }


  /* -------------------------------------
   * Test de definition des points 
   * d'interpolation
   * ------------------------------------- */

  {
    /* Initialisation du couplage */
    bft_printf(" Test 3 : Test de definition des points d'interpolation\n");
    bft_printf("\n");

    couplings_create_coupling("test2D_3",         // Nom du couplage
                              "CodeFortran",                      // Code couplé
                              2,                            // Dimension des entités géométriques
                              0.1,                          // Tolérance géométrique
                              COUPLINGS_STATIC_MESH,        // Maillage statique
                              COUPLINGS_SOLVER_CELL_VERTEX, // Type de champs
                              1,                            // Frequence des post-traitement
                              "EnSight Gold",               // Format du post-traitement
                              "text");                      // Options de post-traitements
    
    /* Construction du maillage local (Decoupage par Metis si plusieurs procs) */
    
    int nVertex = 0;               // Nombre de sommets
    double *coords = NULL;         // Coordonnees des sommets
    int nElts = 0;                 // Nombre d'elements
    int *eltsConnecPointer = NULL; // Index par element dans la connectivite
    int *eltsConnec = NULL;        // Description de la connectivite
    
    const double xmin = -100;
    const double xmax =  100;
    const double ymin = -100;
    const double ymax =  100;
    const int    nx   = 16;
    const int    ny   = 16;
    const int   order = 1;
   
    creeMaillagePolygone2D(order,
                           localComm,
                           xmin,
                           xmax,
                           ymin,
                           ymax,
                           1,
                           nx,
                           ny,
                           &nVertex,
                           &coords,
                           &nElts,
                           &eltsConnecPointer,
                           &eltsConnec);

    bft_printf("   nombre de sommets : %i\n", nVertex);
    bft_printf("   nombre d'elements : %i\n", nElts);
    
    couplings_define_mesh("test2D_3",
                          nVertex,
                          nElts,
                          coords,
                          eltsConnecPointer,
                          eltsConnec);
    
    /* Definition des points a localiser */

    const int nptstolocate = 21;
    double *coordstolocate = NULL;
    BFT_MALLOC(coordstolocate, 3*nptstolocate, double); 

    coordstolocate[0] = -75.;
    coordstolocate[1] = -75.;
    coordstolocate[2] = 0.;

    coordstolocate[3] = -50.;
    coordstolocate[4] = -75.;
    coordstolocate[5] = 0.;

    coordstolocate[6] = -25.;
    coordstolocate[7] = -75.;
    coordstolocate[8] = 0.;

    coordstolocate[9] = 0.;
    coordstolocate[10] = -75;
    coordstolocate[11] = 0.;

    coordstolocate[12] = 25.;
    coordstolocate[13] = -75.;
    coordstolocate[14] = 0.;

    coordstolocate[15] = 50.;
    coordstolocate[16] = -75.;
    coordstolocate[17] = 0.;

    coordstolocate[18] = 75.;
    coordstolocate[19] = -75.;
    coordstolocate[20] = 0.;

    coordstolocate[21] = -75.;
    coordstolocate[22] = 25.;
    coordstolocate[23] = 0.;

    coordstolocate[24] = -50.;
    coordstolocate[25] = 25.;
    coordstolocate[26] = 0.;

    coordstolocate[27] = -25.;
    coordstolocate[28] = 25.;
    coordstolocate[29] = 0.;

    coordstolocate[30] = 0.;
    coordstolocate[31] = 25.;
    coordstolocate[32] = 0.;

    coordstolocate[33] = 25.;
    coordstolocate[34] = 25.;
    coordstolocate[35] = 0.;

    coordstolocate[36] = 50.;
    coordstolocate[37] = 25.;
    coordstolocate[38] = 0.;

    coordstolocate[39] = 75.;
    coordstolocate[40] = 25.;
    coordstolocate[41] = 0.;

    coordstolocate[42] = -75.;
    coordstolocate[43] = 50.;
    coordstolocate[44] = 0.;

    coordstolocate[45] = -50.;
    coordstolocate[46] = 50.;
    coordstolocate[47] = 0.;

    coordstolocate[48] = -25.;
    coordstolocate[49] = 50.;
    coordstolocate[50] = 0.;

    coordstolocate[51] = 0.;
    coordstolocate[52] = 50.;
    coordstolocate[53] = 0.;

    coordstolocate[54] = 25.;
    coordstolocate[55] = 50.;
    coordstolocate[56] = 0.;

    coordstolocate[57] = 50.;
    coordstolocate[58] = 50.;
    coordstolocate[59] = 0.;

    coordstolocate[60] = 75.;
    coordstolocate[61] = 50.;
    coordstolocate[62] = 0.;

    couplings_set_points_to_locate ("test2D_3", nptstolocate, coordstolocate);

    /* Envoi de la coordonnee X
       Reception de la coordonnee Y*/
    
    double* values = NULL;
    BFT_MALLOC(values, nVertex, double);
    for (int i = 0; i < nVertex; i++)
      values[i] = coords[3*i];
    
    double* localValues = NULL;
    BFT_MALLOC(localValues, nptstolocate, double);
    
    int nNotLocatedPoints;
    
    couplings_exchange_status_t status = couplings_exchange("test2D_3",
                                                            "echange1",
                                                            1,   //stride
                                                            1,     // n_step
                                                            0.1,   // physical_time
                                                            "cooX",
                                                            values,
                                                            "cooY",
                                                            localValues,
                                                            &nNotLocatedPoints);
    _dumpStatus(status);
    _dumpNotLocatedPoints("test2D_3", nNotLocatedPoints);

    bft_printf("valeurs recues test2D_3 :\n");
    for (int i = 0; i <  nptstolocate; i++)
      bft_printf("%f ",localValues[i]);
    bft_printf("\n");

    /* Suppression de l'objet de couplage */
    
    couplings_delete_coupling("test2D_3");
    
    /* Liberation de la memoire */
    
    if (coords != NULL)
      BFT_FREE(coords);
    
    if (eltsConnec != NULL)
      BFT_FREE(eltsConnec);
    
    if (eltsConnecPointer != NULL)
      BFT_FREE(eltsConnecPointer);
    
    if (values != NULL)
      BFT_FREE(values);
    
    if (localValues != NULL)
      BFT_FREE(localValues);

    if (coordstolocate != NULL)
      BFT_FREE(coordstolocate);
    bft_printf("--------------------------------------------------------\n");

  }

  /* -------------------------------------
   * Test de definition d'une fonction 
   * d'interpolation (callback)
   * ------------------------------------- */

  {

    /* Initialisation du couplage */
    bft_printf("Test 4 : Test de definition d'une fonction d'interpolation\n");
    bft_printf("\n");

    couplings_create_coupling("test2D_4",                   // Nom du couplage
                              "CodeFortran",                // Code couplé
                              2,                            // Dimension des entités géométriques
                              0.1,                          // Tolérance géométrique
                              COUPLINGS_STATIC_MESH,        // Maillage statique
                              COUPLINGS_SOLVER_CELL_VERTEX, // Type de champs
                              1,                            // Frequence des post-traitement
                              "EnSight Gold",               // Format du post-traitement
                              "text");                      // Options de post-traitements
    
    /* Definition d'une interpolation bidon, pour tester la fonctionnalité */

    couplings_set_interpolation_function("test2D_4", _interpolationBidon); 

    /* Construction du maillage local (Decoupage par Metis si plusieurs procs) */
     
    int nVertex = 0;               // Nombre de sommets
    double *coords = NULL;         // Coordonnees des sommets
    int nElts = 0;                 // Nombre d'elements
    int *eltsConnecPointer = NULL; // Index par element dans la connectivite
    int *eltsConnec = NULL;        // Description de la connectivite
    
    const double xmin = -100;
    const double xmax =  100;
    const double ymin = -100;
    const double ymax =  100;
    const int    nx   = 16;
    const int    ny   = 16;
    const int   order = 1;
    
    creeMaillagePolygone2D(order,
                           localComm,
                           xmin,
                           xmax,
                           ymin,
                           ymax,
                           1,
                           nx,
                           ny,
                           &nVertex,
                           &coords,
                           &nElts,
                           &eltsConnecPointer,
                           &eltsConnec);
    
     couplings_define_mesh("test2D_4",
                          nVertex,
                          nElts,
                          coords,
                          eltsConnecPointer,
                          eltsConnec);
    
    /* Envoi de la coordonnee X
       Reception de la coordonnee Y*/
    
    double* values = NULL;
    BFT_MALLOC(values, nVertex, double);
    for (int i = 0; i < nVertex; i++)
      values[i] = coords[3*i];
    
    double* localValues = NULL;
    BFT_MALLOC(localValues, nVertex, double);
    
    int nNotLocatedPoints;
    
    couplings_exchange_status_t status = couplings_exchange("test2D_4",
                                                            "echange1",
                                                            1, // stride
                                                            1,     // n_step
                                                            0.1,   // physical_time
                                                            "cooX",
                                                            values,
                                                            "cooY",
                                                            localValues,
                                                            &nNotLocatedPoints);
    _dumpStatus(status);
    _dumpNotLocatedPoints("test2D_4", nNotLocatedPoints);

    /* Suppression de l'objet de couplage */
    
    couplings_delete_coupling("test2D_4");
    
    /* Liberation de la memoire */
    
    if (coords != NULL)
      BFT_FREE(coords);
    
    if (eltsConnec != NULL)
      BFT_FREE(eltsConnec);
    
    if (eltsConnecPointer != NULL)
      BFT_FREE(eltsConnecPointer);
    
    if (values != NULL)
      BFT_FREE(values);
    
    if (localValues != NULL)
      BFT_FREE(localValues);
    bft_printf("--------------------------------------------------------\n");
  }

 
  /* -------------------------------------
   * test de la transmission d'un vecteur
   * ------------------------------------- */

  {
    /* Initialisation du couplage */
    bft_printf("Test 5 : test de la transmission d'un vecteur\n");
    bft_printf("\n");

    couplings_create_coupling("test2D_5",         // Nom du couplage
                              "CodeFortran",                      // Code couplé
                              2,                            // Dimension des entités géométriques
                              0.1,                          // Tolérance géométrique
                              COUPLINGS_STATIC_MESH,        // Maillage statique
                              COUPLINGS_SOLVER_CELL_VERTEX, // Type de champs
                              1,                            // Frequence des post-traitement
                              "EnSight Gold",               // Format du post-traitement
                              "text");                      // Options de post-traitements
    
    /* Construction du maillage local (Decoupage par Metis si plusieurs procs) */
    
    bft_printf("CodeA : construction du maillage\n");
    
    int nVertex = 0;               // Nombre de sommets
    double *coords = NULL;         // Coordonnees des sommets
    int nElts = 0;                 // Nombre d'elements
    int *eltsConnecPointer = NULL; // Index par element dans la connectivite
    int *eltsConnec = NULL;        // Description de la connectivite
    
    const double xmin = -100;
    const double xmax =  100;
    const double ymin = -100;
    const double ymax =  100;
    const int    nx   = 16;
    const int    ny   = 16;
    const int   order = 1;
    
    creeMaillagePolygone2D(order,
                           localComm,
                           xmin,
                           xmax,
                           ymin,
                           ymax,
                           1,
                           nx,
                           ny,
                           &nVertex,
                           &coords,
                           &nElts,
                           &eltsConnecPointer,
                           &eltsConnec);
    
    bft_printf("   nombre de sommets : %i\n", nVertex);
    bft_printf("   nombre d'elements : %i\n", nElts);
    
    couplings_define_mesh("test2D_5",
                          nVertex,
                          nElts,
                          coords,
                          eltsConnecPointer,
                          eltsConnec);
    
    /* Envoi de la coordonnee X
       Reception de la coordonnee Y*/
    
    double* values = NULL;
    BFT_MALLOC(values, nVertex, double);
    for (int i = 0; i < nVertex; i++)
      values[i] = coords[3*i];
    
    double* localValues = NULL;
    BFT_MALLOC(localValues, 3*nVertex, double);
    
    int nNotLocatedPoints;
    
    couplings_exchange_status_t status = couplings_exchange("test2D_5",
                                                            "echange1",
                                                            3, //stride
                                                            1,     // n_step
                                                            0.1,   // physical_time
                                                            "cooX",
                                                            coords,
                                                            "cooY",
                                                            localValues,
                                                            &nNotLocatedPoints);
    _dumpStatus(status);
    _dumpNotLocatedPoints("test2D_5", nNotLocatedPoints);

    /* Suppression de l'objet de couplage */
    
    couplings_delete_coupling("test2D_5");
    
    /* Liberation de la memoire */
    
    if (coords != NULL)
      BFT_FREE(coords);
    
    if (eltsConnec != NULL)
      BFT_FREE(eltsConnec);
    
    if (eltsConnecPointer != NULL)
      BFT_FREE(eltsConnecPointer);
    
    if (values != NULL)
      BFT_FREE(values);
    
    if (localValues != NULL)
      BFT_FREE(localValues);
    bft_printf("--------------------------------------------------------\n");
  }

  /* -------------------------------------
   * test des sorties d'erreur
   * ------------------------------------- */

  {
   
    /* Initialisation du couplage */
    bft_printf("Test 6 : Test des sorties d'erreur\n");
    bft_printf("\n");
    
    couplings_create_coupling("test2D_6",         // Nom du couplage
                              "CodeFortran",                      // Code couplé
                              2,                            // Dimension des entités géométriques
                              0.1,                          // Tolérance géométrique
                              COUPLINGS_STATIC_MESH,        // Maillage statique
                              COUPLINGS_SOLVER_CELL_VERTEX, // Type de champs
                              1,                            // Frequence des post-traitement
                              "EnSight Gold",               // Format du post-traitement
                              "text");                      // Options de post-traitements
    
    /* Construction du maillage local (Decoupage par Metis si plusieurs procs) */
    
    bft_printf("CodeA : construction du maillage\n");
    
    int nVertex = 0;               // Nombre de sommets
    double *coords = NULL;         // Coordonnees des sommets
    int nElts = 0;                 // Nombre d'elements
    int *eltsConnecPointer = NULL; // Index par element dans la connectivite
    int *eltsConnec = NULL;        // Description de la connectivite
    
    const double xmin = -150;
    const double xmax =  150;
    const double ymin = -150;
    const double ymax =  150;
    const int    nx   =  20;
    const int    ny   =  20;
    const int   order =  1;
    
    creeMaillagePolygone2D(order,
                           localComm,
                           xmin,
                           xmax,
                           ymin,
                           ymax,
                           1,
                           nx,
                           ny,
                           &nVertex,
                           &coords,
                           &nElts,
                           &eltsConnecPointer,
                           &eltsConnec);
    
    bft_printf("   nombre de sommets : %i\n", nVertex);
    bft_printf("   nombre d'elements : %i\n", nElts);
    
    couplings_define_mesh("test2D_6",
                          nVertex,
                          nElts,
                          coords,
                          eltsConnecPointer,
                          eltsConnec);
    
    /* Envoi de la coordonnee X
       Reception de la coordonnee Y*/
    
    double* values = NULL;
    BFT_MALLOC(values, nVertex, double);
    for (int i = 0; i < nVertex; i++)
      values[i] = coords[3*i];
    
    double* localValues = NULL;
    BFT_MALLOC(localValues, nVertex, double);
    
    int nNotLocatedPoints;
    
    couplings_exchange_status_t status = couplings_exchange("test2D_6",
                                                            "echange1",
                                                            1, // stride
                                                            1,     // n_step
                                                            0.1,   // physical_time
                                                            "cooX",
                                                            values,
                                                            "cooY",
                                                            localValues,
                                                            &nNotLocatedPoints);
    _dumpStatus(status);
    _dumpNotLocatedPoints("test2D_6", nNotLocatedPoints);

    //
    // Reception mais aucun envoi par code fortran
    // Controle de status

    status = couplings_exchange("test2D_6",
                                "echange2",
                                1,  // stride
                                1,     // n_step
                                0.1,   // physical_time
                                "cooX",
                                values,
                                "cooY",
                                localValues,
                                &nNotLocatedPoints);

    _dumpStatus(status);
    _dumpNotLocatedPoints("test2D_6", nNotLocatedPoints);
    /* Suppression de l'objet de couplage */
    
    couplings_delete_coupling("test2D_6");
    
    /* Liberation de la memoire */
    
    if (coords != NULL)
      BFT_FREE(coords);
    
    if (eltsConnec != NULL)
      BFT_FREE(eltsConnec);
    
    if (eltsConnecPointer != NULL)
      BFT_FREE(eltsConnecPointer);
    
    if (values != NULL)
      BFT_FREE(values);
    
    if (localValues != NULL)
      BFT_FREE(localValues);
    bft_printf("--------------------------------------------------------\n");
  }

  /* -------------------------
   * Test avec tri des elements
   * ----------------------- */
    
  {
    /* Initialisation du couplage */
    bft_printf("Test 7 :  Test tri du maillage\n");
    bft_printf("\n");
    
    couplings_create_coupling("test2D_7",         // Nom du couplage
                              "CodeFortran",                      // Code couplé
                              2,                            // Dimension des entités géométriques
                              0.1,                          // Tolérance géométrique
                              COUPLINGS_STATIC_MESH,        // Maillage statique
                              COUPLINGS_SOLVER_CELL_VERTEX, // Type de champs
                              1,                            // Frequence des post-traitement
                              "EnSight Gold",               // Format du post-traitement
                              "text");                      // Options de post-traitements
    
    /* Construction du maillage local (Decoupage par Metis si plusieurs procs) */
    
    int nVertex = 0;               // Nombre de sommets
    double *coords = NULL;         // Coordonnees des sommets
    int nElts = 0;                 // Nombre d'elements
    int *eltsConnecPointer = NULL; // Index par element dans la connectivite
    int *eltsConnec = NULL;        // Description de la connectivite
    
    const double xmin = -100;
    const double xmax =  100;
    const double ymin = -100;
    const double ymax =  100;
    const int    nx   = 20;
    const int    ny   = 20;
    const int   order = 1;
    
    creeMaillagePolygone2D(order,
                           localComm,
                           xmin,
                           xmax,
                           ymin,
                           ymax,
                           1,
                           nx,
                           ny,
                           &nVertex,
                           &coords,
                           &nElts,
                           &eltsConnecPointer,
                           &eltsConnec);
    
    bft_printf("   nombre de sommets : %i\n", nVertex);
    bft_printf("   nombre d'elements : %i\n", nElts);
    
    couplings_define_mesh("test2D_7",
                          nVertex,
                          nElts,
                          coords,
                          eltsConnecPointer,
                          eltsConnec);
    
    /* Envoi de la coordonnee X
       Reception de la coordonnee Y*/
    
    double* values = NULL;
    BFT_MALLOC(values, nVertex, double);
    for (int i = 0; i < nVertex; i++)
      values[i] = coords[3*i];
    
    double* localValues = NULL;
    BFT_MALLOC(localValues, nVertex, double);

    int nNotLocatedPoints;
    
    couplings_exchange_status_t status = couplings_exchange("test2D_7",
                                                            "echange1",
                                                            1, // stride
                                                            1,     // n_step
                                                            0.1,   // physical_time
                                                            "cooX",
                                                            values,
                                                            "cooY",
                                                            localValues,
                                                            &nNotLocatedPoints);
    _dumpStatus(status);
    _dumpNotLocatedPoints("test2D_7", nNotLocatedPoints);
    /* Suppression de l'objet de couplage */
    
    couplings_delete_coupling("test2D_7");
    
    /* Liberation de la memoire */
    
    if (coords != NULL)
      BFT_FREE(coords);
    
    if (eltsConnec != NULL)
      BFT_FREE(eltsConnec);
    
    if (eltsConnecPointer != NULL)
      BFT_FREE(eltsConnecPointer);
    
    if (values != NULL)
      BFT_FREE(values);
    
    if (localValues != NULL)
      BFT_FREE(localValues);
    bft_printf("--------------------------------------------------------\n");
  }

  /* -------------------------------------------------------
   * Test avec tri des elements pour un maillage cell-center
   * ------------------------------------------------------- */

  {   
    /* Initialisation du couplage */
    
    bft_printf("Test 8 : Test avec tri des elements pour un maillage cell-center");
    bft_printf("\n");
    couplings_create_coupling("test2D_8",         // Nom du couplage
                              "CodeFortran",                      // Code couplé
                              2,                            // Dimension des entités géométriques
                              0.1,                          // Tolérance géométrique 
                              COUPLINGS_STATIC_MESH,        // Maillage statique
                              COUPLINGS_SOLVER_CELL_CENTER, // Type de champs
                              1,                            // Frequence des post-traitement
                              "EnSight Gold",               // Format du post-traitement
                              "text");                      // Options de post-traitements
    
    /* Construction du maillage local (Decoupage par Metis si plusieurs procs) */
    
    int nVertex = 0;               // Nombre de sommets
    double *coords = NULL;         // Coordonnees des sommets
    int nElts = 0;                 // Nombre d'elements
    int *eltsConnecPointer = NULL; // Index par element dans la connectivite
    int *eltsConnec = NULL;        // Description de la connectivite
    
    const double xmin = -100;
    const double xmax =  100;
    const double ymin = -100;
    const double ymax =  100;
    const int    nx   = 16;
    const int    ny   = 16;
    const int   order = 1;
    
    creeMaillagePolygone2D(order,
                           localComm,
                           xmin, 
                           xmax, 
                           ymin, 
                           ymax,
                           1,
                           nx,
                           ny,
                           &nVertex,
                           &coords,
                           &nElts,
                           &eltsConnecPointer,
                           &eltsConnec);
    
    bft_printf("   nombre de sommets : %i\n", nVertex);
    bft_printf("   nombre d'elements : %i\n", nElts);
    
    couplings_define_mesh("test2D_8", 
                          nVertex,
                          nElts,
                          coords,
                          eltsConnecPointer,
                          eltsConnec);
    
    /* Reception de la coordonnee Y
       Envoi de la coordonnee Y */
    
    double* localValues = NULL;
    BFT_MALLOC(localValues, nElts, double);

    /* Receive */
    
    int nNotLocatedPoints;
    
    couplings_exchange_status_t status = couplings_exchange("test2D_8",
                                                            "echange1",
                                                            1, // stride
                                                            1,     // n_step
                                                            0.1,   // physical_time
                                                            NULL,
                                                            NULL,
                                                            "cooY",
                                                            localValues,
                                                            &nNotLocatedPoints);
    bft_printf("Send\n");
    _dumpStatus(status);
    _dumpNotLocatedPoints("test2D_8", nNotLocatedPoints);

    /* Send */

    status = couplings_exchange("test2D_8",
                                "echange2",
                                1, // stride
                                1,     // n_step
                                0.1,   // physical_time
                                "cooYY",
                                localValues,
                                NULL,
                                NULL,
                                &nNotLocatedPoints);
    bft_printf("Receive\n");
    _dumpStatus(status);
    _dumpNotLocatedPoints("test2D_8", nNotLocatedPoints);

    /* Suppression de l'objet de couplage */
    
    couplings_delete_coupling("test2D_8");
    
    /* Liberation de la memoire */
    
    if (coords != NULL)
      BFT_FREE(coords);
    
    if (eltsConnec != NULL)
      BFT_FREE(eltsConnec);
    
    if (eltsConnecPointer != NULL)
      BFT_FREE(eltsConnecPointer);
    
    if (localValues != NULL)
      BFT_FREE(localValues);
    bft_printf("--------------------------------------------------------\n");
    
  }

  /* ------------------------
   * Test simple localisation
   * ------------------------ */
    
  {

    /* Initialisation du couplage */
    bft_printf("Test 9 :  Test simple localisation\n");
    bft_printf("\n");
    
    couplings_create_coupling("test2D_9",         // Nom du couplage
                              "CodeFortran",                      // Code couplé
                              2,                            // Dimension des entités géométriques
                              0.1,                          // Tolérance géométrique
                              COUPLINGS_STATIC_MESH,        // Maillage statique
                              COUPLINGS_SOLVER_CELL_VERTEX, // Type de champs
                              1,                            // Frequence des post-traitement
                              "EnSight Gold",               // Format du post-traitement
                              "text");                      // Options de post-traitements
    
    /* Construction du maillage local (Decoupage par Metis si plusieurs procs) */
    
    int nVertex = 0;               // Nombre de sommets
    double *coords = NULL;         // Coordonnees des sommets
    int nElts = 0;                 // Nombre d'elements
    int *eltsConnecPointer = NULL; // Index par element dans la connectivite
    int *eltsConnec = NULL;        // Description de la connectivite
    
    const double xmin = -100;
    const double xmax =  100;
    const double ymin = -100;
    const double ymax =  100;
    const int    nx   = 16;
    const int    ny   = 16;
    const int   order = 1;
    
    creeMaillagePolygone2D(order,
                           localComm,
                           xmin,
                           xmax,
                           ymin,
                           ymax,
                           1,
                           nx,
                           ny,
                           &nVertex,
                           &coords,
                           &nElts,
                           &eltsConnecPointer,
                           &eltsConnec);
    
    bft_printf("   nombre de sommets : %i\n", nVertex);
    bft_printf("   nombre d'elements : %i\n", nElts);
    
    couplings_define_mesh("test2D_9",
                          nVertex,
                          nElts,
                          coords,
                          eltsConnecPointer,
                          eltsConnec);
    
    couplings_locate("test2D_9");

    const int n_located_distant_point = couplings_get_n_located_distant_points("test2D_9");

    const int *distant_location = couplings_get_distant_location("test2D_9");

    const int *distant_barycentric_coordinates_index = couplings_get_distant_barycentric_coordinates_index("test2D_9");

    const double* distant_barycentric_coordinates = couplings_get_distant_barycentric_coordinates("test2D_9");

    /* Suppression de l'objet de couplage */
    
    couplings_delete_coupling("test2D_9");
    
    /* Liberation de la memoire */
    
    if (coords != NULL)
      BFT_FREE(coords);
    
    if (eltsConnec != NULL)
      BFT_FREE(eltsConnec);
    
    if (eltsConnecPointer != NULL)
      BFT_FREE(eltsConnecPointer);
    
    bft_printf("--------------------------------------------------------\n");
  }

  /* Fin des communications MPI */
  /* -------------------------- */
  
  couplings_finalize();
  
  bft_mem_end();
  fclose(outputFile);
  return 0;
}
