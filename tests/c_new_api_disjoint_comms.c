#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include <pdm.h>
#include <pdm_dcube_gen.h>
#include <pdm_part.h>
#include <cwp.h>
#include <stdbool.h>

static void
create_dcube
        (
                PDM_MPI_Comm pdm_comm, PDM_g_num_t n_vtx_seg, double length, double xmin, double ymin, double zmin,
                int **n_vtx, int **n_faces, int **n_cells, double ***coord, int ***face_vtx_idx, int ***face_vtx, int ***cell_face_idx, int ***cell_face,
                PDM_g_num_t ***vtx_ln_to_gn, PDM_g_num_t ***cell_ln_to_gn
        ) {
    int n_part = 1;

    // Create PDM cube mesh
    PDM_dcube_t *dcube = PDM_dcube_gen_init(pdm_comm, n_vtx_seg, length, xmin, ymin, zmin, PDM_OWNERSHIP_UNGET_RESULT_IS_FREE);

    // Create mesh partitions
    int d_n_cell, d_n_face, d_n_vertices, s_face_vtx, s_face_group, n_face_group;
    int *d_face_vertex_idx = NULL, *d_face_group_idx = NULL;
    PDM_g_num_t *d_face_cell = NULL, *d_face_vertex = NULL, *d_face_group = NULL;
    double *d_vertex_coord = NULL;
    PDM_dcube_gen_dim_get(dcube, &n_face_group, &d_n_cell, &d_n_face, &d_n_vertices, &s_face_vtx, &s_face_group);
    PDM_dcube_gen_data_get(dcube, &d_face_cell, &d_face_vertex_idx, &d_face_vertex, &d_vertex_coord, &d_face_group_idx, &d_face_group);

    int ppart_id = 0;
    int *d_cell_part = (int *) malloc(sizeof(int) * d_n_cell);
    int have_dcell_part = 0;
    int *renum_properties_cell = NULL;
    int *renum_properties_face = NULL;
    int n_property_cell = 0;
    int n_property_face = 0;

    PDM_part_split_t method = PDM_PART_SPLIT_PTSCOTCH;
    PDM_part_create(&ppart_id, pdm_comm, method, "PDM_PART_RENUM_CELL_NONE", "PDM_PART_RENUM_FACE_NONE",
                    n_property_cell, renum_properties_cell, n_property_face, renum_properties_face, n_part,
                    d_n_cell, d_n_face, d_n_vertices, n_face_group, NULL, NULL, NULL, NULL, have_dcell_part, d_cell_part,
                    d_face_cell, d_face_vertex_idx, d_face_vertex, NULL, d_vertex_coord, NULL, d_face_group_idx, d_face_group);
    free(d_face_vertex_idx);
    free(d_face_group_idx);
    free(d_face_cell);
    free(d_face_vertex);
    free(d_face_group);
    free(d_vertex_coord);
    free(d_cell_part);

    // Get connectivity
    // Cells
    *n_cells = (int *) malloc(sizeof(int) * n_part);
    *cell_face_idx = (int **) malloc(sizeof(int *) * n_part);
    *cell_face = (int **) malloc(sizeof(int *) * n_part);
    *cell_ln_to_gn = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t *) * n_part);

    // Faces
    *n_faces = (int *) malloc(sizeof(int) * n_part);
    *face_vtx_idx = (int **) malloc(sizeof(int *) * n_part);
    *face_vtx = (int **) malloc(sizeof(int *) * n_part);

    // Vertices
    *n_vtx = (int *) malloc(sizeof(int) * n_part);
    *coord = (double **) malloc(sizeof(double *) * n_part);
    *vtx_ln_to_gn = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t *) * n_part);

    for (int i_part = 0 ; i_part < n_part ; i_part++) {
        int _n_cells, _n_faces, _n_face_part_bound, _n_vtx, _n_proc, _n_total_part, _s_cell_face, _s_face_vtx, _s_face_group, _n_edge_group2;
        PDM_part_part_dim_get(ppart_id, i_part, &_n_cells, &_n_faces, &_n_face_part_bound, &_n_vtx, &_n_proc, &_n_total_part,
                              &_s_cell_face, &_s_face_vtx, &_s_face_group, &_n_edge_group2);

        int *_cell_tag, *_cell_face_idx, *_cell_face, *_face_tag, *_face_cell, *_face_vtx_idx, *_face_vtx,
                *_face_part_bound_proc_idx, *_face_part_bound_part_idx, *_face_part_bound, *_vtx_tag, *_face_group_idx, *_face_group;
        double *_vtx_coords;
        PDM_g_num_t *_cell_ln_to_gn, *_face_ln_to_gn, *_vtx_ln_to_gn, *_face_group_ln_to_gn;
        PDM_part_part_val_get(ppart_id, i_part, &_cell_tag, &_cell_face_idx, &_cell_face, &_cell_ln_to_gn,
                              &_face_tag, &_face_cell, &_face_vtx_idx, &_face_vtx, &_face_ln_to_gn,
                              &_face_part_bound_proc_idx, &_face_part_bound_part_idx, &_face_part_bound, &_vtx_tag,
                              &_vtx_coords, &_vtx_ln_to_gn, &_face_group_idx, &_face_group, &_face_group_ln_to_gn);

        // Cells
        *n_cells[i_part] = _n_cells;
        *cell_face_idx[i_part] = (int *) malloc(sizeof(int) * (_n_cells + 1));
        *cell_face[i_part] = (int *) malloc(sizeof(int) * _s_cell_face);
        *cell_ln_to_gn[i_part] = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * _n_cells);

        memcpy(*cell_face_idx[i_part], _cell_face_idx, (_n_cells + 1) * sizeof(int));
        memcpy(*cell_face[i_part], _cell_face, _s_cell_face * sizeof(int));
        memcpy(*cell_ln_to_gn[i_part], _cell_ln_to_gn, _n_cells * sizeof(PDM_g_num_t));

        // Faces
        *n_faces[i_part] = _n_faces;
        *face_vtx_idx[i_part] = (int *) malloc(sizeof(int) * (_n_faces + 1));
        *face_vtx[i_part] = (int *) malloc(sizeof(int) * _s_face_vtx);

        memcpy(*face_vtx_idx[i_part], _face_vtx_idx, (_n_faces + 1) * sizeof(int));
        memcpy(*face_vtx[i_part], _face_vtx, _s_face_vtx * sizeof(int));

        // Vertices
        *n_vtx[i_part] = _n_vtx;
        *coord[i_part] = (double *) malloc(sizeof(double) * (3 * _n_vtx));
        *vtx_ln_to_gn[i_part] = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * _n_vtx);

        memcpy(*coord[i_part], _vtx_coords, 3 * _n_vtx * sizeof(double));
        memcpy(*vtx_ln_to_gn[i_part], _vtx_ln_to_gn, _n_vtx * sizeof(PDM_g_num_t));
    }

    // Free
    PDM_part_free(ppart_id);
    PDM_dcube_gen_free (dcube);
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    int rank;
    int comm_world_size;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_world_size);
    assert(comm_world_size > 0);

    // Input
//    bool cond_code1 = rank % 2 == 0;
//    bool cond_code2 = rank % 2 == 1;
    bool cond_code1 = rank == 0;
    bool cond_code2 = rank == 0 || rank == 1;
    bool cond_both = cond_code1 && cond_code2;

    int n_vtx_seg_code1 = 3, n_vtx_seg_code2 = 4;
    double x_min_code1 = 0., x_min_code2 = 0.2;
    double y_min_code1 = 0., y_min_code2 = 0.8;

    // Define the number of codes per rank
    int n_code;
    if (cond_both) n_code = 2;
    else if (cond_code1 || cond_code2) n_code = 1;
    else n_code = 0;

    int *code_id = (int *) malloc(n_code * sizeof(int));
    const char **code_names = (const char **) malloc(n_code * sizeof(char *));
    const char **coupled_code_names = (const char **) malloc(n_code * sizeof(char *));
    double *time_init = (double *) malloc(n_code * sizeof(double));
    PDM_g_num_t *n_vtx_seg = (PDM_g_num_t *) malloc(n_code * sizeof(PDM_g_num_t));
    double *x_min = (double *) malloc(n_code * sizeof(double));
    double *y_min = (double *) malloc(n_code * sizeof(double));
    CWP_Status_t *is_active_rank = (CWP_Status_t *) malloc(n_code * sizeof(CWP_Status_t));

    MPI_Comm *intra_comms = (MPI_Comm *) malloc(n_code * sizeof(MPI_Comm));

    // Define which rank works for which code
    if (cond_both) {
        code_id[0] = 1;
        code_id[1] = 2;
        code_names[0] = "code1";
        code_names[1] = "code2";
        coupled_code_names[0] = "code2";
        coupled_code_names[1] = "code1";
        time_init[0] = 0.;
        time_init[1] = 0.;
        is_active_rank[0] = CWP_STATUS_ON;
        is_active_rank[1] = CWP_STATUS_ON;
        n_vtx_seg[0] = n_vtx_seg_code1;
        n_vtx_seg[1] = n_vtx_seg_code2;
        x_min[0] = x_min_code1;
        x_min[1] = x_min_code2;
        y_min[0] = y_min_code1;
        y_min[1] = y_min_code2;
    }
    else if (cond_code1) {
        code_id[0] = 1;
        code_names[0] = "code1";
        coupled_code_names[0] = "code2";
        time_init[0] = 0.;
        is_active_rank[0] = CWP_STATUS_ON;
        n_vtx_seg[0] = n_vtx_seg_code1;
        x_min[0] = x_min_code1;
        y_min[0] = y_min_code1;
    }
    else if (cond_code2) {
        code_id[0] = 2;
        code_names[0] = "code2";
        coupled_code_names[0] = "code1";
        time_init[0] = 0.;
        is_active_rank[0] = CWP_STATUS_ON;
        n_vtx_seg[0] = n_vtx_seg_code2;
        x_min[0] = x_min_code2;
        y_min[0] = y_min_code2;
    }

    // Init cwipi
    CWP_Init(MPI_COMM_WORLD, n_code, (const char **) code_names, is_active_rank, time_init, intra_comms);
    printf("%d --- CWIPI initialized\n", rank);

    // Get the comm size and rank
    int *intra_comm_rank = (int *) malloc(n_code * sizeof(int));
    int *intra_comm_size = (int *) malloc(n_code * sizeof(int));
    for (int i_code = 0 ; i_code < n_code ; ++i_code) {
        MPI_Comm_rank(intra_comms[i_code], &intra_comm_rank[i_code]);
        MPI_Comm_size(intra_comms[i_code], &intra_comm_size[i_code]);
        assert(intra_comm_size[i_code] > 0);
    }

    // Gather ranks and master ranks
    int *ranks_on_code1 = (int *) malloc(comm_world_size * sizeof(int));
    int *ranks_on_code2 = (int *) malloc(comm_world_size * sizeof(int));
    int master_code1 = -1, master_code2 = -1, master_both = -1;

    int comm_nb = 0;
    if (cond_code1) {
        MPI_Allgather(&rank, 1, MPI_INT, ranks_on_code1, 1, MPI_INT, intra_comms[0]);
        MPI_Allreduce(&rank, &master_code1, 1, MPI_INT, MPI_MIN, intra_comms[0]);
    }
    if (cond_code2) {
        if (n_code == 1) comm_nb = 0;
        else if (n_code == 2) comm_nb = 1;
        MPI_Allgather(&rank, 1, MPI_INT, ranks_on_code2, 1, MPI_INT, intra_comms[comm_nb]);
        MPI_Allreduce(&rank, &master_code2, 1, MPI_INT, MPI_MIN, intra_comms[comm_nb]);
    }
    if (cond_both) {
        for (int i = 0 ; i < intra_comm_size[0] ; ++i) {
            for (int j = 0 ; j < intra_comm_size[1] ; ++j) {
                if (ranks_on_code1[i] == ranks_on_code2[j]) {
                    master_both = ranks_on_code1[i];
                    break;
                }
            }
        }
    }

    // Print the number and ranks for each code
    if (rank == master_code1) {
        printf("%d --- %d procs work for code %d (%.1f %%): ", rank, intra_comm_size[0], code_id[0], (double) intra_comm_size[0] / comm_world_size * 100);
        for (int i = 0 ; i < intra_comm_size[0] ; ++i) printf("%d ", ranks_on_code1[i]);
        printf("\n");
    }
    if (rank == master_code2) {
        if (n_code == 1) comm_nb = 0;
        else if (n_code == 2) comm_nb = 1;
        printf("%d --- %d procs work for code %d (%.1f %%): ", rank, intra_comm_size[comm_nb], code_id[comm_nb], (double) intra_comm_size[comm_nb] / comm_world_size * 100);
        for (int i = 0 ; i < intra_comm_size[comm_nb] ; ++i) printf("%d ", ranks_on_code2[i]);
        printf("\n");
    }
    if (rank == master_both) {
        int tmp_code1 = -1, tmp_code2 = -1;
        int nb_both_codes = 0;
        int *ranks_on_both = (int *) malloc(comm_world_size * sizeof(int));
        for (int i = 0 ; i < comm_world_size ; ++i) {
            for (int j = 0 ; j < intra_comm_size[0] ; ++j) {
                if (i == ranks_on_code1[j]) tmp_code1 = i;
            }
            for (int j = 0 ; j < intra_comm_size[1] ; ++j) {
                if (i == ranks_on_code2[j]) tmp_code2 = i;
            }
            if (tmp_code1 != -1 && tmp_code2 != -1) ranks_on_both[nb_both_codes++] = i;
            tmp_code1 = -1, tmp_code2 = -1;
        }
        printf("%d --- %d procs work for both codes (%.1f %%): ", rank, nb_both_codes, (double) nb_both_codes / comm_world_size * 100);
        for (int i = 0 ; i < nb_both_codes ; ++i) printf("%d ", ranks_on_both[i]);
        printf("\n");
    }

    // Create coupling and visu
    const char *cpl_name = "c_new_api_disjoint_comms";
    CWP_Spatial_interp_t interp_method = CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE;

    for (int i_code = 0 ; i_code < n_code ; ++i_code) {
        CWP_Cpl_create(code_names[i_code], cpl_name, coupled_code_names[i_code], CWP_INTERFACE_VOLUME, CWP_COMM_PAR_WITH_PART, interp_method, 1, CWP_DYNAMIC_MESH_STATIC, CWP_TIME_EXCH_EACH_TIME_STEP);
        printf("%d (%d, %s) --- Coupling created between %s and %s\n", rank, intra_comm_rank[i_code], code_names[i_code], code_names[i_code], coupled_code_names[i_code]);
    }

    for (int i_code = 0 ; i_code < n_code ; ++i_code) {
        CWP_Visu_set(code_names[i_code], cpl_name, 1, CWP_VISU_FORMAT_ENSIGHT, "text");
        printf("%d (%d, %s) --- Visu set\n", rank, intra_comm_rank[i_code], code_names[i_code]);
    }

    // Create PDM communicators
    PDM_MPI_Comm *pdm_intra_comms = (PDM_MPI_Comm *) malloc(n_code * sizeof(PDM_MPI_Comm));
    for (int i_code = 0 ; i_code < n_code ; ++i_code) {
        pdm_intra_comms[i_code] = PDM_MPI_mpi_2_pdm_mpi_comm((void *) &intra_comms[i_code]);
        printf("%d (%d, %s) --- PDM comm created\n", rank, intra_comm_rank[i_code], code_names[i_code]);
    }

    // Create geometry
    int **n_vtx, **n_faces, **n_cells;
    int ***face_vtx_idx, ***face_vtx;
    int ***cell_face, ***cell_face_idx;
    PDM_g_num_t ***vtx_ln_to_gn, ***cell_ln_to_gn;
    double ***coord;
    n_vtx = (int **) malloc(n_code * sizeof(int **));
    n_faces = (int **) malloc(n_code * sizeof(int **));
    n_cells = (int **) malloc(n_code * sizeof(int **));
    face_vtx_idx = (int ***) malloc(n_code * sizeof(int ***));
    face_vtx = (int ***) malloc(n_code * sizeof(int ***));
    cell_face_idx = (int ***) malloc(n_code * sizeof(int ***));
    cell_face = (int ***) malloc(n_code * sizeof(int ***));
    coord = (double ***) malloc(n_code * sizeof(double ***));
    vtx_ln_to_gn = (PDM_g_num_t ***) malloc(n_code * sizeof(PDM_g_num_t ***));
    cell_ln_to_gn = (PDM_g_num_t ***) malloc(n_code * sizeof(PDM_g_num_t ***));

    for (int i_code = 0 ; i_code < n_code ; ++i_code) {
        create_dcube(pdm_intra_comms[i_code], n_vtx_seg[i_code], 1., x_min[i_code], y_min[i_code], 0., &(n_vtx[i_code]), &n_faces[i_code], &n_cells[i_code], &coord[i_code], &face_vtx_idx[i_code], &face_vtx[i_code], &cell_face_idx[i_code], &cell_face[i_code], &vtx_ln_to_gn[i_code], &cell_ln_to_gn[i_code]);
        printf("%d (%d, %s) --- dcube created\n", rank, intra_comm_rank[i_code], code_names[i_code]);
    }

    // Set geometry
    for (int i_code = 0 ; i_code < n_code ; ++i_code) {
        CWP_Mesh_interf_vtx_set(code_names[i_code], cpl_name, 0, n_vtx[i_code][0], coord[i_code][0], vtx_ln_to_gn[i_code][0]);
        CWP_Mesh_interf_from_cellface_set(code_names[i_code], cpl_name, 0, n_cells[i_code][0], cell_face_idx[i_code][0], cell_face[i_code][0], n_faces[i_code][0], face_vtx_idx[i_code][0], face_vtx[i_code][0], cell_ln_to_gn[i_code][0]);
        CWP_Mesh_interf_finalize(code_names[i_code], cpl_name);
        printf("%d (%d, %s) --- Geometry set\n", rank, intra_comm_rank[i_code], code_names[i_code]);
    }

    // Create and initialise Fields: code1 -> code2
    const char *field_name = "coo";

    double **send_values = (double **) malloc(n_code * sizeof(double *));
    double **recv_values = (double **) malloc(n_code * sizeof(double *));
    for (int i_code = 0 ; i_code < n_code ; ++i_code) {
        send_values[i_code] = (double *) malloc(3 * n_vtx[i_code][0] * sizeof(double));
        recv_values[i_code] = (double *) malloc(3 * n_vtx[i_code][0] * sizeof(double));
        if (code_id[i_code] == 1) {
          for (int i = 0 ; i < 3 * n_vtx[i_code][0] ; i++) {
            send_values[i_code][i] = coord[i_code][0][i];
          }
        }
        if (code_id[i_code] == 2) {
          for (int i = 0 ; i < 3 * n_vtx[i_code][0] ; i++) {
            recv_values[i_code][i] = 0;
          }
        }
    }

    for (int i_code = 0 ; i_code < n_code ; ++i_code) {
        if (code_id[i_code] == 1) {
            CWP_Field_create(code_names[i_code], cpl_name, field_name, CWP_DOUBLE, CWP_FIELD_STORAGE_BLOCK, 3, CWP_DOF_LOCATION_NODE, CWP_FIELD_EXCH_SEND, CWP_STATUS_ON);
            CWP_Field_data_set(code_names[i_code], cpl_name, field_name, 0, CWP_FIELD_MAP_SOURCE, send_values[i_code]);
        }
        if (code_id[i_code] == 2) {
            CWP_Field_create(code_names[i_code], cpl_name, field_name, CWP_DOUBLE, CWP_FIELD_STORAGE_BLOCK, 3, CWP_DOF_LOCATION_NODE, CWP_FIELD_EXCH_RECV, CWP_STATUS_ON);
            CWP_Field_data_set(code_names[i_code], cpl_name, field_name, 0, CWP_FIELD_MAP_TARGET, recv_values[i_code]);
        }
        printf("%d (%d, %s) --- Field created and data set\n", rank, intra_comm_rank[i_code], code_names[i_code]);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // Compute weights
    for (int i_code = 0 ; i_code < n_code ; ++i_code) {
        CWP_Spatial_interp_weights_compute(code_names[i_code], cpl_name);
        printf("%d (%d, %s) --- Weights computed\n", rank, intra_comm_rank[i_code], code_names[i_code]);
    }

    int n_computed_tgts = 0, n_uncomputed_tgts = 0;
    const int *computed_tgts = NULL, *uncomputed_tgts = NULL;
    if (cond_code2) {
        n_computed_tgts = CWP_N_computed_tgts_get("code2", cpl_name, field_name, 0);
        n_uncomputed_tgts = CWP_N_uncomputed_tgts_get("code2", cpl_name, field_name, 0);
        computed_tgts = (int *) malloc(n_computed_tgts * sizeof(int));
        uncomputed_tgts = (int *) malloc(n_uncomputed_tgts * sizeof(int));
        computed_tgts = CWP_Computed_tgts_get("code2", cpl_name, field_name, 0);
        uncomputed_tgts = CWP_Uncomputed_tgts_get("code2", cpl_name, field_name, 0);
        printf("%d (%d, %s) --- n computed targets: %d\n", rank, intra_comm_rank[0], code_names[0], n_computed_tgts);
        printf("%d (%d, %s) --- n uncomputed targets: %d\n", rank, intra_comm_rank[0], code_names[0], n_uncomputed_tgts);
        if (n_computed_tgts != 0) {
            printf("%d (%d, %s) --- computed targets: ", rank, intra_comm_rank[0], code_names[0]);
            for (int i = 0 ; i < n_computed_tgts ; ++i) printf("%d ", computed_tgts[i]);
            printf("\n");
        }
        if (n_uncomputed_tgts != 0) {
            printf("%d (%d, %s) --- uncomputed targets: ", rank, intra_comm_rank[0], code_names[0]);
            for (int i = 0 ; i < n_uncomputed_tgts ; ++i) printf("%d ", uncomputed_tgts[i]);
            printf("\n");
        }
    }

    // Send and receive field
    for (int i_code = 0 ; i_code < n_code ; i_code++) {
        if (code_id[i_code] == 2) {
            CWP_Field_irecv(code_names[i_code], cpl_name, field_name);
            printf("%d (%d, %s) --- Received field\n", rank, intra_comm_rank[i_code], code_names[i_code]);
        }
        if (code_id[i_code] == 1) {
            CWP_Field_issend(code_names[i_code], cpl_name, field_name);
            printf("%d (%d, %s) --- Sent field\n", rank, intra_comm_rank[i_code], code_names[i_code]);
        }
    }

    for (int i_code = 0 ; i_code < n_code ; i_code++) {
        if (code_id[i_code] == 2) {
            CWP_Field_wait_irecv(code_names[i_code], cpl_name, field_name);
            printf("%d (%d, %s) --- Received field\n", rank, intra_comm_rank[i_code], code_names[i_code]);
        }
        if (code_id[i_code] == 1) {
            CWP_Field_wait_issend(code_names[i_code], cpl_name, field_name);
            printf("%d (%d, %s) --- Sent field\n", rank, intra_comm_rank[i_code], code_names[i_code]);
        }
    }

    for (int i_code = 0 ; i_code < n_code ; ++i_code) {
      if (code_id[i_code] == 2) {
        for (int i = 0 ; i < n_computed_tgts ; i++) {
            printf("%12.5e %12.5e\n", recv_values[i_code][3 * i], coord[i_code][0][3 * (computed_tgts[i] - 1)]);
//          assert(fabs(recv_values[i_code][3 * i] - coord[i_code][0][3 * (computed_tgts[i] - 1)]) < 1.e-4);
        }
      }
    }

  // Delete interf
    for (int i_code = 0 ; i_code < n_code ; i_code++) {
        CWP_Mesh_interf_del(code_names[i_code], cpl_name);
        printf("%d (%d, %s) --- Interface deleted\n", rank, intra_comm_rank[i_code], code_names[i_code]);
    }

    // Delete coupling
    for (int i_code = 0 ; i_code < n_code ; i_code++) {
        CWP_Cpl_del(code_names[i_code], cpl_name);
        printf("%d (%d, %s) --- Coupling deleted\n", rank, intra_comm_rank[i_code], code_names[i_code]);
    }


    //Finalize cwipi
    CWP_Finalize();
    printf("%d --- CWIPI finalized\n", rank);

    MPI_Finalize();
    exit(0);
}
