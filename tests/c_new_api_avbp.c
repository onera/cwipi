#include <mpi.h>

#include "cwp.h"
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
    int code_id;
    if (rank == 0) {
        code_name = "avbp1";
        code_id = 1;
        code_coupled_name = "avbp2";
    }
    else if (rank == 1) {
        code_name = "avbp2";
        code_id = 2;
        code_coupled_name = "avbp1";
    }
    else {
        printf("Invalid rank\n");
        exit(0);
    }

    int nb_codes = 1;
    const char **code_names = (const char **) malloc(sizeof(char *) * nb_codes);
    double *time_init = (double *) malloc(sizeof(double) * nb_codes);
    CWP_Status_t *statuses = (CWP_Status_t *) malloc(sizeof(CWP_Status_t) * nb_codes);
    MPI_Comm *cwp_comms = (MPI_Comm *) malloc(sizeof(MPI_Comm) * nb_codes);

    code_names[0] = code_name;
    time_init[0] = 0.;
    statuses[0] = CWP_STATUS_ON;

    CWP_Init(commcwipi, nb_codes, (const char **) code_names, statuses, time_init, cwp_comms);
//    MPI_Comm commcodepcw = cwp_comms[0];
    printf("%d --- CWIPI initialised\n", rank);

    CWP_Cpl_create(code_name, coupling_name, code_coupled_name,
                   CWP_INTERFACE_SURFACE, CWP_COMM_PAR_WITH_PART, CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE,
                   1, CWP_DYNAMIC_MESH_VARIABLE, CWP_TIME_EXCH_CPL_TIME_STEP);
    printf("%d --- Coupling created\n", rank);

    int nnode_cpl_s, ncell_cpl_s, nnode_cpl_d;
    int n_vtx_per_cell = 4;

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

    double *xyz_s     = (double *) malloc(sizeof(double) * (3 * nnode_cpl_s));
    int    *connind_s = (int *)    malloc(sizeof(int)    * (ncell_cpl_s + 1));
    int    *elems_s   = (int *)    malloc(sizeof(int)    * (n_vtx_per_cell * ncell_cpl_s));
    double *xyz_dest  = (double *) malloc(sizeof(double) * (3 * nnode_cpl_d));

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
    for (int i = 0 ; i <  n_var * nnode_cpl_d ; i++) {
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

    CWP_Mesh_interf_vtx_set(code_name, coupling_name, 0, nnode_cpl_s, xyz_s, NULL);
    int block_id = CWP_Mesh_interf_block_add(code_name, coupling_name, CWP_BLOCK_FACE_QUAD4);
    CWP_Mesh_interf_block_std_set(code_name, coupling_name, 0, block_id, ncell_cpl_s, elems_s, NULL);
    CWP_Mesh_interf_finalize(code_name, coupling_name);
    printf("%d --- Geometry set\n", rank);

    const char *field = "direct";
    const char *fieldr = "reverse";

    // field:  avbp1 -> avbp2
    // fieldr: avbp2-> avbp1
    if (code_id == 1) {
        CWP_Field_create(code_name, coupling_name, field, CWP_DOUBLE, CWP_FIELD_STORAGE_BLOCK, stride,
                         CWP_DOF_LOCATION_NODE, CWP_FIELD_EXCH_SEND, CWP_STATUS_ON);
        CWP_Field_data_set(code_name,
                           coupling_name,
                           field,
                           0,
                           CWP_FIELD_MAP_SOURCE,
                           sfields);

        CWP_Field_create(code_name, coupling_name, fieldr, CWP_DOUBLE, CWP_FIELD_STORAGE_BLOCK, stride,
                         CWP_DOF_LOCATION_USER, CWP_FIELD_EXCH_RECV, CWP_STATUS_ON);
        CWP_Field_data_set(code_name,
                           coupling_name,
                           fieldr,
                           0,
                           CWP_FIELD_MAP_TARGET,
                           rfields);
    }
    else if (code_id == 2) {
        CWP_Field_create(code_name, coupling_name, fieldr, CWP_DOUBLE, CWP_FIELD_STORAGE_BLOCK, stride,
                         CWP_DOF_LOCATION_NODE, CWP_FIELD_EXCH_SEND, CWP_STATUS_ON);
        CWP_Field_data_set(code_name,
                           coupling_name,
                           fieldr,
                           0,
                           CWP_FIELD_MAP_SOURCE,
                           sfields);

        CWP_Field_create(code_name, coupling_name, field, CWP_DOUBLE, CWP_FIELD_STORAGE_BLOCK, stride,
                         CWP_DOF_LOCATION_USER, CWP_FIELD_EXCH_RECV, CWP_STATUS_ON);
        CWP_Field_data_set(code_name,
                           coupling_name,
                           field,
                           0,
                           CWP_FIELD_MAP_TARGET,
                           rfields);
    }

    CWP_User_tgt_pts_set(code_name, coupling_name, 0, nnode_cpl_d, xyz_dest, NULL);
    printf("%d --- Points to locate set\n", rank);

    printf("%d --- Fields created\n", rank);

    CWP_Spatial_interp_weights_compute(code_name, coupling_name);
    printf("%d --- Localisation done\n", rank);

    int n_not_located_points = CWP_N_uncomputed_tgts_get(code_name, coupling_name, "fieldr", 0);
    int n_located_points = CWP_N_computed_tgts_get(code_name, coupling_name, "fieldr", 0);
    const int* located_points = CWP_Computed_tgts_get(code_name, coupling_name, "fieldr", 0);

    printf("n_not_located_points %d\n", n_not_located_points);
    printf("n_located_points %d\n", n_located_points);
    if (rank == 1) {
        for (int i = 0 ; i < n_located_points ; ++i) {
            printf("%d ", located_points[i]);
        }
        printf("\n");
    }

    if (code_id == 1) {
        CWP_Field_issend(code_name, coupling_name, field);
        CWP_Field_wait_issend(code_name, coupling_name, field);

        CWP_Field_irecv(code_name, coupling_name, fieldr);
        CWP_Field_wait_irecv(code_name, coupling_name, fieldr);
    }
    else if (code_id == 2) {
        CWP_Field_irecv(code_name, coupling_name, field);
        CWP_Field_wait_irecv(code_name, coupling_name, field);

        CWP_Field_issend(code_name, coupling_name, fieldr);
        CWP_Field_wait_issend(code_name, coupling_name, fieldr);
    }
    printf("%d --- Exchanges completed\n", rank);

  CWP_Cpl_del(code_name, coupling_name);
  CWP_Finalize();
  exit(0);

  sprintf(fn, "rfields_out_cnew_%d.out", rank);
    fil = fopen(fn, "w");
    for (int i = 0 ; i <  n_var * nnode_cpl_d ; ++i) {
        fprintf(fil, "%E        ", rfields[i]);
    }
    fclose(fil);
    printf(" --- Output written\n");

    CWP_Cpl_del(code_name, coupling_name);
    printf("%d --- Coupling deleted completed\n", rank);

    CWP_Finalize();
    printf("%d --- CWIPI finalized\n", rank);

    free(code_names);
    free(time_init);
    free(statuses);
    free(cwp_comms);
    free(xyz_s);
    free(connind_s);
    free(elems_s);
    free(xyz_dest);
    free(sfields);
    free(rfields);

    MPI_Finalize();
}
