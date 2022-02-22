#include <mpi.h>

#include "cwipi.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char *argv[]) {
  int mode;
  if (argc == 1) {
    mode = 0;
  }
  else if (argc == 2) {
    if (strcmp(argv[1], "-0") == 0) mode = 0;
    else if (strcmp(argv[1], "-1") == 0) mode = 1;
    else {
      mode = -1;
      printf("Unknown argument value\n");
      exit(0);
    }
  }
  else {
    mode = -1;
    printf("Wrong number of arguments\n");
    exit(0);
  }

  MPI_Init(&argc, &argv);

  int rank;
  int comm_world_size;

  MPI_Comm commcwipi = MPI_COMM_WORLD;
  MPI_Comm_rank(commcwipi, &rank);
  MPI_Comm_size(commcwipi, &comm_world_size);

  const char *coupling_name = "cpl12";

  const char *code_name, *code_coupled_name;
  if (rank == 0) {
    code_name = "avbp1";
    code_coupled_name = "avbp2";
  }
  else if (rank == 1) {
    code_name = "avbp2";
    code_coupled_name = "avbp1";
  }
  else {
    printf("Invalid rank\n");
    exit(0);
  }

  int nb_codes = 1;
  const char **code_names = (const char **) malloc(sizeof(char *) * nb_codes);
  double *time_init = (double *) malloc(sizeof(double) * nb_codes);
  MPI_Comm *cwp_comms = (MPI_Comm *) malloc(sizeof(MPI_Comm) * nb_codes);

  code_names[0] = code_name;
  time_init[0] = 0.;

  cwipi_init(MPI_COMM_WORLD, code_name, cwp_comms);
  printf("%d --- CWIPI initialised\n", rank);

  cwipi_create_coupling(coupling_name, CWIPI_COUPLING_PARALLEL_WITH_PARTITIONING, code_coupled_name, 2, 0.01, CWIPI_STATIC_MESH, CWIPI_SOLVER_CELL_VERTEX, 1, "EnSight Gold", "text");
  printf("%d --- Coupling created\n", rank);

  int nnode_cpl_s, ncell_cpl_s, nnode_cpl_d;
  int n_vtx_per_cell = 8;

  char fn[100];
  sprintf(fn, "avbp/avbp_dimensions_%d.out", rank);
  FILE *fil = fopen(fn, "r");
  fscanf(fil, "%d", &nnode_cpl_s);
  fscanf(fil, "%d", &ncell_cpl_s);
  fscanf(fil, "%d", &nnode_cpl_d);
  fclose(fil);

  int n_var, stride;
  if (mode == 0) {
    n_var = 6;
    stride = 4;
  }
  else {
    n_var = 1;
    stride = 1;
  }

  double *xyz_s = (double *) malloc(sizeof(double) * (3 * nnode_cpl_s));
  int *connind_s = (int *) malloc(sizeof(int) * (ncell_cpl_s + 1));
  int *elems_s = (int *) malloc(sizeof(int) * (n_vtx_per_cell * ncell_cpl_s));
  double *xyz_dest = (double *) malloc(sizeof(double) * (3 * nnode_cpl_d));

  double *sfields = (double *) malloc(sizeof(double) * (n_var * nnode_cpl_s));
  double *rfields = (double *) malloc(sizeof(double) * (n_var * nnode_cpl_d));

  // xyz_s
  sprintf(fn, "avbp/avbp_xyz_s_%d.out", rank);
  fil = fopen(fn, "r");
  for (int i = 0 ; i < 3 * nnode_cpl_s ; i++) {
    fscanf(fil, "%lf", &xyz_s[i]);
  }
  fclose(fil);
  // connind_s
  sprintf(fn, "avbp/avbp_connind_s_%d.out", rank);
  fil = fopen(fn, "r");
  for (int i = 0 ; i < 3 * ncell_cpl_s + 1 ; i++) {
    fscanf(fil, "%d", &connind_s[i]);
  }
  fclose(fil);
  // elems_s
  sprintf(fn, "avbp/avbp_elems_s_%d.out", rank);
  fil = fopen(fn, "r");
  for (int i = 0 ; i < 3 * n_vtx_per_cell * ncell_cpl_s ; i++) {
    fscanf(fil, "%d", &elems_s[i]);
  }
  fclose(fil);
  // xyz_dest
  sprintf(fn, "avbp/avbp_xyz_dest_%d.out", rank);
  fil = fopen(fn, "r");
  for (int i = 0 ; i < 3 * nnode_cpl_d ; i++) {
    fscanf(fil, "%lf", &xyz_dest[i]);
  }
  fclose(fil);

  // Fields
  for (int i = 0 ; i < n_var * nnode_cpl_d ; i++) {
    rfields[i] = 0.;
  }
  if (mode == 0) {
    // avbp_sfields
    sprintf(fn, "avbp/avbp_sfields_%d.out", rank);
    fil = fopen(fn, "r");
    for (int i = 0 ; i < n_var * nnode_cpl_s ; i++) {
      fscanf(fil, "%lf", &sfields[i]);
    }
    fclose(fil);
  }
  else {
    // Test valeurs non uniformes entrelacees
    for (int i = 0 ; i < n_var * nnode_cpl_s ; i++) {
      sfields[i] = 2. * (i + 1);
    }
  }

  cwipi_define_mesh(coupling_name, nnode_cpl_s, ncell_cpl_s, xyz_s, connind_s, elems_s);
  printf("%d --- Geometry set\n", rank);

  cwipi_set_points_to_locate(coupling_name, nnode_cpl_d, xyz_dest);
  printf("%d --- Points to locate set\n", rank);

  cwipi_locate(coupling_name);
  printf("%d --- Localisation done\n", rank);

  int n_not_located_points;
  int n_located_points;
  int n_distant_located_points;
  const int *located_points, *not_located_points;
  const int *distant_located_points;

  n_not_located_points = cwipi_get_n_not_located_points(coupling_name);
  n_located_points = cwipi_get_n_located_points(coupling_name);
  n_distant_located_points = cwipi_get_n_distant_points(coupling_name);
  located_points = cwipi_get_located_points(coupling_name);
  not_located_points = cwipi_get_not_located_points(coupling_name);

  distant_located_points = cwipi_get_distant_location(coupling_name);

  printf("%d --- n_not_located_points %d\n", rank, n_not_located_points);
  printf("%d --- n_located_points %d\n", rank, n_located_points);
  printf("%d --- n_distant_located_points %d\n", rank, n_distant_located_points);
  printf("not_located_points :");
  for (int i = 0 ; i < n_not_located_points ; ++i) printf("%d ", not_located_points[i]);
  printf("\nlocated_points: ");
  for (int i = 0 ; i < n_located_points ; ++i) printf("%d ", located_points[i]);
  printf("\ndistant_located_points: ");
  for (int i = 0 ; i < n_distant_located_points ; ++i) printf("%d ", distant_located_points[i]);
  printf("\n");

  cwipi_delete_coupling(coupling_name);
  printf("%d --- Coupling deleted\n", rank);

  cwipi_finalize();
  printf("%d --- CWIPI finalized\n", rank);

  free(code_names);
  free(time_init);
  free(cwp_comms);
  free(xyz_s);
  free(connind_s);
  free(elems_s);
  free(xyz_dest);
  free(sfields);
  free(rfields);

  MPI_Finalize();
}
