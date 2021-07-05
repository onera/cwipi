#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include <cwp.h>
#include <pdm_dcube_gen.h>
#include <pdm_part.h>

int main(int argc, char *argv[]) {
    // Init MPI
    MPI_Init(&argc, &argv);
    int rank;
    int comm_world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_world_size);
    assert(comm_world_size > 1);

    // Init CWIPI
    int n_code = 1;
    char **code_name = malloc(sizeof(char *) * n_code);
    char **coupled_code_name = malloc(sizeof(char *) * n_code);
    CWP_Status_t *is_active_rank = malloc(sizeof(CWP_Status_t) * n_code);
    double *time_init = malloc(sizeof(double) * n_code);
    MPI_Comm *intra_comm = malloc(sizeof(MPI_Comm) * n_code);

    is_active_rank[0] = CWP_STATUS_ON;
    time_init[0] = 0.;
    if (rank < comm_world_size / 2) {
        code_name[0] = "code1";
        coupled_code_name[0] = "code2";
    }
    else {
        code_name[0] = "code2";
        coupled_code_name[0] = "code1";
    }

    CWP_Init(MPI_COMM_WORLD, n_code, (const char **) code_name, is_active_rank, time_init, intra_comm);
    printf("%d: CWIPI Init OK, %s\n", rank, code_name[0]);

    // Create CWIPI coupling
    int n_part = 1;
    char *coupling_name = "c_new_api_dcube";
    CWP_Spatial_interp_t interp_method = CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE;
    CWP_Cpl_create(code_name[0], coupling_name, coupled_code_name[0], CWP_INTERFACE_VOLUME, CWP_COMM_PAR_WITH_PART, interp_method, n_part, CWP_DYNAMIC_MESH_STATIC, CWP_TIME_EXCH_EACH_TIME_STEP);

    // Setup visualisation
    CWP_Visu_set(code_name[0], coupling_name, 1, CWP_VISU_FORMAT_ENSIGHT, "text");

    // Create PDM cube mesh
    PDM_g_num_t n_vtx_seg = 3;
    double length = 1.;
    const double xmin = 0, ymin = 0., zmin = 0.;
    PDM_dcube_t *dcube = PDM_dcube_gen_init(PDM_MPI_COMM_WORLD, n_vtx_seg, length, xmin, ymin, zmin, PDM_OWNERSHIP_KEEP);

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
    PDM_part_create(&ppart_id, PDM_MPI_COMM_WORLD, method, "PDM_PART_RENUM_CELL_NONE", "PDM_PART_RENUM_FACE_NONE",
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
    int *n_cells = (int *) malloc(sizeof(int) * n_part);
    int **cell_face_idx = (int **) malloc(sizeof(int *) * n_part);
    int **cell_face = (int **) malloc(sizeof(int *) * n_part);
    PDM_g_num_t **cell_ln_to_gn = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t *) * n_part);

    // Faces
    int *n_faces = (int *) malloc(sizeof(int) * n_part);
    int **face_vtx_idx = (int **) malloc(sizeof(int *) * n_part);
    int **face_vtx = (int **) malloc(sizeof(int *) * n_part);

    // Vertices
    int *n_vtx = (int *) malloc(sizeof(int) * n_part);
    double **vtx_coord = (double **) malloc(sizeof(double *) * n_part);
    PDM_g_num_t **vtx_ln_to_gn = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t *) * n_part);

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
        n_cells[i_part] = _n_cells;
        cell_face_idx[i_part] = (int *) malloc(sizeof(int) * (_n_cells + 1));
        cell_face[i_part] = (int *) malloc(sizeof(int) * _s_cell_face);
        cell_ln_to_gn[i_part] = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * _n_cells);

        memcpy(cell_face_idx[i_part], _cell_face_idx, (_n_cells + 1) * sizeof(int));
        memcpy(cell_face[i_part], _cell_face, _s_cell_face * sizeof(int));
        memcpy(cell_ln_to_gn[i_part], _cell_ln_to_gn, _n_cells * sizeof(PDM_g_num_t));

        // Faces
        n_faces[i_part] = _n_faces;
        face_vtx_idx[i_part] = (int *) malloc(sizeof(int) * (_n_faces + 1));
        face_vtx[i_part] = (int *) malloc(sizeof(int) * _s_face_vtx);

        memcpy(face_vtx_idx[i_part], _face_vtx_idx, (_n_faces + 1) * sizeof(int));
        memcpy(face_vtx[i_part], _face_vtx, _s_face_vtx * sizeof(int));

        // Vertices
        n_vtx[i_part] = _n_vtx;
        vtx_coord[i_part] = (double *) malloc(sizeof(double) * (3 * _n_vtx));
        vtx_ln_to_gn[i_part] = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * _n_vtx);

        memcpy(vtx_coord[i_part], _vtx_coords, 3 * _n_vtx * sizeof(double));
        memcpy(vtx_ln_to_gn[i_part], _vtx_ln_to_gn, _n_vtx * sizeof(PDM_g_num_t));
    }

    // Free
    PDM_part_free(ppart_id);

    // Set CWIPI mesh
    CWP_Mesh_interf_vtx_set(code_name[0], coupling_name, 0, n_vtx[0], vtx_coord[0], NULL);
    CWP_Mesh_interf_from_cellface_set(code_name[0], coupling_name, 0, n_cells[0], cell_face_idx[0], cell_face[0], n_faces[0], face_vtx_idx[0], face_vtx[0], cell_ln_to_gn[0]);
    CWP_Mesh_interf_finalize(code_name[0], coupling_name);

    // Set fields
    double *send_val = NULL;
    double *recv_val = NULL;
    char *field_name = "cooX";

    if (strcmp(code_name[0], "code1")) {
        send_val = (double *) malloc(sizeof(double) * n_vtx[0]);
        for (int i = 0 ; i < n_vtx[0] ; i++) send_val[i] = vtx_coord[0][3 * i];
    }
    else recv_val = (double *) malloc(sizeof(double) * n_vtx[0]);

    CWP_Status_t visu_status = CWP_STATUS_OFF;
    MPI_Barrier(MPI_COMM_WORLD);

    if (strcmp(code_name[0], "code1")) {
        CWP_Field_create(code_name[0], coupling_name, field_name, CWP_DOUBLE, CWP_FIELD_STORAGE_BLOCK, 1, CWP_DOF_LOCATION_NODE, CWP_FIELD_EXCH_SEND, visu_status);
        CWP_Field_data_set(code_name[0], coupling_name, field_name, 0, send_val);
    }
    else {
        CWP_Field_create(code_name[0], coupling_name, field_name, CWP_DOUBLE, CWP_FIELD_STORAGE_BLOCK, 1, CWP_DOF_LOCATION_NODE, CWP_FIELD_EXCH_RECV, visu_status);
        CWP_Field_data_set(code_name[0], coupling_name, field_name, 0, recv_val);
    }

    // Compute weights
    CWP_Spatial_interp_weights_compute(code_name[0], coupling_name);

    // Finalize
    CWP_Mesh_interf_del(code_name[0], coupling_name);
    CWP_Cpl_del(code_name[0], coupling_name);

    CWP_Finalize();
    MPI_Finalize();
}
