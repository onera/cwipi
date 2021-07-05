#include <cstdlib>
#include <cstring>

#include "cwp.h"
#include <cassert>

using namespace std;

static char *CWP_fortran_to_c_string(const char *application_name_f, const int l_application_name_f) {
    char *application_name_c;
    int imin = 0;
    int imax = 0;

    while (imin < l_application_name_f && application_name_f[imin] == ' ') imin++;
    while (imax < l_application_name_f && application_name_f[l_application_name_f - imax - 1] == ' ') imax++;

    imax = l_application_name_f - imax - 1;

    assert(imax >= imin);

    if ((imax == l_application_name_f) || (imin == l_application_name_f)) {
        application_name_c = new char[1];
        application_name_c[0] = '\0';
    }
    else {
        int size = imax - imin + 2;
        application_name_c = new char[size];
        int index = 0;
        for (int k = imin ; k <= imax ; k++) application_name_c[index++] = application_name_f[k];
        application_name_c[index] = '\0';
    }

    return application_name_c;
}

void CWP_Init_cf(MPI_Fint f_global_comm, const int n_code, const char **f_code_names, const int *l_code_names, const int *is_active_rank, const double *time_init, MPI_Fint *f_intra_comms) {
    // Convert code names dealing with different size characters
    char **c_code_names = (char **) malloc(n_code * sizeof(char *));
    for (int i = 0 ; i < n_code ; i++) {
        c_code_names[i] = (char *) malloc((l_code_names[i] - 1) * sizeof(char));
        memcpy(c_code_names[i], *f_code_names + i * l_code_names[i], l_code_names[i]); // TODO This is not perfectly right for the length of the f_code_name which could vary
    }

    // Convert global MPI communicator
    MPI_Comm c_global_comm = MPI_Comm_f2c(f_global_comm);

    // Allocate local communicators in C
    auto *c_intra_comms = (MPI_Comm *) malloc(n_code * sizeof(MPI_Comm));

    CWP_Init(c_global_comm, n_code, (const char **) c_code_names, (CWP_Status_t *) is_active_rank, time_init, c_intra_comms);

    // Convert local communicators to Fortran
    for (int i = 0 ; i < n_code ; i++) f_intra_comms[i] = MPI_Comm_c2f(c_intra_comms[i]);

    free(c_code_names);
    free(c_intra_comms);
}

void CWP_Cpl_create_cf(const char *f_local_code_name, const int l_local_code_name, const char *f_cpl_id, const int l_cpl_id, const char *f_coupled_code_name, const int l_coupled_code_name, const CWP_Interface_t entities_dim, const CWP_Comm_t comm_type, const CWP_Spatial_interp_t spatial_interp, const int n_part, const CWP_Dynamic_mesh_t displacement, const CWP_Time_exch_t freq) {
    char *c_local_code_name, *c_cpl_id, *c_coupled_code_name;

    c_local_code_name = CWP_fortran_to_c_string(f_local_code_name, l_local_code_name);
    c_cpl_id = CWP_fortran_to_c_string(f_cpl_id, l_cpl_id);
    c_coupled_code_name = CWP_fortran_to_c_string(f_coupled_code_name, l_coupled_code_name);

    CWP_Cpl_create((const char *) c_local_code_name, (const char *) c_cpl_id, (const char *) c_coupled_code_name, entities_dim, comm_type, spatial_interp, n_part, displacement, freq);

    free(c_local_code_name);
    free(c_cpl_id);
    free(c_coupled_code_name);
}

void CWP_Cpl_del_cf(const char *f_local_code_name, const int l_local_code_name, const char *f_cpl_id, const int l_cpl_id) {
    char *c_local_code_name, *c_cpl_id;

    c_local_code_name = CWP_fortran_to_c_string(f_local_code_name, l_local_code_name);
    c_cpl_id = CWP_fortran_to_c_string(f_cpl_id, l_cpl_id);

    CWP_Cpl_del(c_local_code_name, c_cpl_id);

    free(c_local_code_name);
    free(c_cpl_id);
}

int CWP_N_uncomputed_tgts_get_cf(const char *f_local_code_name, const int l_local_code_name, const char *f_cpl_id, const int l_cpl_id, CWP_Dof_location_t target_location, int i_part) {
    char *c_local_code_name, *c_cpl_id;

    c_local_code_name = CWP_fortran_to_c_string(f_local_code_name, l_local_code_name);
    c_cpl_id = CWP_fortran_to_c_string(f_cpl_id, l_cpl_id);

    return CWP_N_uncomputed_tgts_get(c_local_code_name, c_cpl_id, target_location, i_part);
}

const int *CWP_Uncomputed_tgts_get_cf(const char *f_local_code_name, const int l_local_code_name, const char *f_cpl_id, const int l_cpl_id) {
    char *c_local_code_name, *c_cpl_id;

    c_local_code_name = CWP_fortran_to_c_string(f_local_code_name, l_local_code_name);
    c_cpl_id = CWP_fortran_to_c_string(f_cpl_id, l_cpl_id);

    return CWP_Uncomputed_tgts_get(c_local_code_name, c_cpl_id);
}

int CWP_N_computed_tgts_get_cf(const char *f_local_code_name, const int l_local_code_name, const char *f_cpl_id, const int l_cpl_id) {
    char *c_local_code_name, *c_cpl_id;

    c_local_code_name = CWP_fortran_to_c_string(f_local_code_name, l_local_code_name);
    c_cpl_id = CWP_fortran_to_c_string(f_cpl_id, l_cpl_id);

    return CWP_N_computed_tgts_get(c_local_code_name, c_cpl_id);
}

const int *CWP_Computed_tgts_get_cf(const char *f_local_code_name, const int l_local_code_name, const char *f_cpl_id, const int l_cpl_id) {
    char *c_local_code_name, *c_cpl_id;

    c_local_code_name = CWP_fortran_to_c_string(f_local_code_name, l_local_code_name);
    c_cpl_id = CWP_fortran_to_c_string(f_cpl_id, l_cpl_id);

    return CWP_Computed_tgts_get(c_local_code_name, c_cpl_id);
}

const double *CWP_Computed_tgts_dist_to_spatial_interp_get_cf(const char *f_local_code_name, const int l_local_code_name, const char *f_cpl_id, const int l_cpl_id) {
    char *c_local_code_name, *c_cpl_id;

    c_local_code_name = CWP_fortran_to_c_string(f_local_code_name, l_local_code_name);
    c_cpl_id = CWP_fortran_to_c_string(f_cpl_id, l_cpl_id);

    return CWP_Computed_tgts_dist_to_spatial_interp_get(c_local_code_name, c_cpl_id);
}

void CWP_Spatial_interp_weights_compute_cf(const char *f_local_code_name, const int l_local_code_name, const char *f_cpl_id, const int l_cpl_id) {
    char *c_local_code_name, *c_cpl_id;

    c_local_code_name = CWP_fortran_to_c_string(f_local_code_name, l_local_code_name);
    c_cpl_id = CWP_fortran_to_c_string(f_cpl_id, l_cpl_id);

    CWP_Spatial_interp_weights_compute(c_local_code_name, c_cpl_id);
}

void CWP_Visu_set_cf(const char *f_local_code_name, const int l_local_code_name, const char *f_cpl_id, const int l_cpl_id, int freq, CWP_Visu_format_t format, const char *f_format_option, const int l_format_option) {
    char *c_local_code_name, *c_cpl_id, *c_format_option;

    c_local_code_name = CWP_fortran_to_c_string(f_local_code_name, l_local_code_name);
    c_cpl_id = CWP_fortran_to_c_string(f_cpl_id, l_cpl_id);
    c_format_option = CWP_fortran_to_c_string(f_format_option, l_format_option);

    CWP_Visu_set(c_local_code_name, c_cpl_id, freq, format, c_format_option);
}

void CWP_Mesh_interf_finalize_cf(const char *f_local_code_name, int l_local_code_name, const char *f_cpl_id, int l_cpl_id) {
    char *c_local_code_name, *c_cpl_id;

    c_local_code_name = CWP_fortran_to_c_string(f_local_code_name, l_local_code_name);
    c_cpl_id = CWP_fortran_to_c_string(f_cpl_id, l_cpl_id);

    CWP_Mesh_interf_finalize(c_local_code_name, c_cpl_id);
}

void CWP_Mesh_interf_vtx_set_cf(const char *f_local_code_name, int l_local_code_name, const char *f_cpl_id, int l_cpl_id, int i_part, int n_pts, double coord[], CWP_g_num_t global_num[]) {
    char *c_local_code_name, *c_cpl_id;

    c_local_code_name = CWP_fortran_to_c_string(f_local_code_name, l_local_code_name);
    c_cpl_id = CWP_fortran_to_c_string(f_cpl_id, l_cpl_id);

    CWP_Mesh_interf_vtx_set(c_local_code_name, c_cpl_id, i_part, n_pts, coord, global_num);
}

int CWP_Mesh_interf_block_add_cf(const char *f_local_code_name, int l_local_code_name, const char *f_cpl_id, int l_cpl_id, CWP_Block_t block_type) {
    char *c_local_code_name, *c_cpl_id;

    c_local_code_name = CWP_fortran_to_c_string(f_local_code_name, l_local_code_name);
    c_cpl_id = CWP_fortran_to_c_string(f_cpl_id, l_cpl_id);

    return CWP_Mesh_interf_block_add(c_local_code_name, c_cpl_id, block_type);
}

void CWP_Mesh_interf_block_std_set_cf(const char *f_local_code_name, int l_local_code_name, const char *f_cpl_id, int l_cpl_id, int i_part, int block_id, int n_elts, int connec[], CWP_g_num_t global_num[]) {
    char *c_local_code_name, *c_cpl_id;

    c_local_code_name = CWP_fortran_to_c_string(f_local_code_name, l_local_code_name);
    c_cpl_id = CWP_fortran_to_c_string(f_cpl_id, l_cpl_id);

    CWP_Mesh_interf_block_std_set(c_local_code_name, c_cpl_id, i_part, block_id, n_elts, connec, global_num);
}

void CWP_Mesh_interf_f_poly_block_set_cf(const char *f_local_code_name, int l_local_code_name, const char *f_cpl_id, int l_cpl_id, int i_part, int block_id, int n_elts, int connec_idx[], int connec[], CWP_g_num_t global_num[]) {
    char *c_local_code_name, *c_cpl_id;

    c_local_code_name = CWP_fortran_to_c_string(f_local_code_name, l_local_code_name);
    c_cpl_id = CWP_fortran_to_c_string(f_cpl_id, l_cpl_id);

    CWP_Mesh_interf_f_poly_block_set(c_local_code_name, c_cpl_id, i_part, block_id, n_elts, connec_idx, connec, global_num);
}

void CWP_Mesh_interf_c_poly_block_set_cf(const char *f_local_code_name, int l_local_code_name, const char *f_cpl_id, int l_cpl_id, int i_part, int block_id, int n_elts, int n_faces,
                                         int connec_faces_idx[], int connec_faces[], int connec_cells_idx[], int connec_cells[], CWP_g_num_t global_num[]) {
    char *c_local_code_name, *c_cpl_id;

    c_local_code_name = CWP_fortran_to_c_string(f_local_code_name, l_local_code_name);
    c_cpl_id = CWP_fortran_to_c_string(f_cpl_id, l_cpl_id);

    CWP_Mesh_interf_c_poly_block_set(c_local_code_name, c_cpl_id, i_part, block_id, n_elts, n_faces, connec_faces_idx, connec_faces, connec_cells_idx, connec_cells, global_num);
}

void CWP_Mesh_interf_from_cellface_set_cf(const char *f_local_code_name, int l_local_code_name, const char *f_cpl_id, int l_cpl_id, int i_part, int n_cells,
                                          int cell_face_idx[], int cell_face[], int n_faces, int face_vtx_idx[], int face_vtx[], CWP_g_num_t parent_num[]) {
    char *c_local_code_name, *c_cpl_id;

    c_local_code_name = CWP_fortran_to_c_string(f_local_code_name, l_local_code_name);
    c_cpl_id = CWP_fortran_to_c_string(f_cpl_id, l_cpl_id);

    CWP_Mesh_interf_from_cellface_set(c_local_code_name, c_cpl_id, i_part, n_cells, cell_face_idx, cell_face, n_faces, face_vtx_idx, face_vtx, parent_num);
}

void CWP_Mesh_interf_from_faceedge_set_cf(const char *f_local_code_name, int l_local_code_name, const char *f_cpl_id, int l_cpl_id, int i_part, int n_faces,
                                          int face_edge_idx[], int face_edge[], int n_edges, int edge_vtx_idx[], int edge_vtx[], CWP_g_num_t parent_num[]) {
    char *c_local_code_name, *c_cpl_id;

    c_local_code_name = CWP_fortran_to_c_string(f_local_code_name, l_local_code_name);
    c_cpl_id = CWP_fortran_to_c_string(f_cpl_id, l_cpl_id);

    CWP_Mesh_interf_from_faceedge_set(c_local_code_name, c_cpl_id, i_part, n_faces, face_edge_idx, face_edge, n_edges, edge_vtx_idx, edge_vtx, parent_num);

}

void CWP_Field_create_cf(const char *f_local_code_name, int l_local_code_name, const char *f_cpl_id, int l_cpl_id, const char *f_field_id, int l_field_id, CWP_Type_t data_type, CWP_Field_storage_t storage, int n_component, CWP_Dof_location_t target_location, CWP_Field_exch_t exch_type, CWP_Status_t visu_status) {
    char *c_local_code_name, *c_cpl_id, *c_field_id;

    c_local_code_name = CWP_fortran_to_c_string(f_local_code_name, l_local_code_name);
    c_cpl_id = CWP_fortran_to_c_string(f_cpl_id, l_cpl_id);
    c_field_id = CWP_fortran_to_c_string(f_field_id, l_field_id);

    CWP_Field_create(c_local_code_name, c_cpl_id, c_field_id, data_type, storage, n_component, target_location, exch_type, visu_status);
}

void CWP_Field_data_set_cf(const char *f_local_code_name, int l_local_code_name, const char *f_cpl_id, int l_cpl_id, const char *f_field_id, int l_field_id, int i_part, double data[]) {
    char *c_local_code_name, *c_cpl_id, *c_field_id;

    c_local_code_name = CWP_fortran_to_c_string(f_local_code_name, l_local_code_name);
    c_cpl_id = CWP_fortran_to_c_string(f_cpl_id, l_cpl_id);
    c_field_id = CWP_fortran_to_c_string(f_field_id, l_field_id);

    CWP_Field_data_set(c_local_code_name, c_cpl_id, c_field_id, i_part, data);
}

void CWP_Field_issend_cf(const char *f_local_code_name, const int l_local_code_name, const char *f_cpl_id, const int l_cpl_id, const char *f_src_field_id, const int l_src_field_id) {
    char *c_local_code_name, *c_cpl_id, *c_src_field_id;

    c_local_code_name = CWP_fortran_to_c_string(f_local_code_name, l_local_code_name);
    c_cpl_id = CWP_fortran_to_c_string(f_cpl_id, l_cpl_id);
    c_src_field_id = CWP_fortran_to_c_string(f_src_field_id, l_src_field_id);

    CWP_Field_issend(c_local_code_name, c_cpl_id, c_src_field_id);

    free(c_local_code_name);
    free(c_cpl_id);
    free(c_src_field_id);
}

void CWP_Field_irecv_cf(const char *f_local_code_name, const int l_local_code_name, const char *f_cpl_id, const int l_cpl_id, const char *f_tgt_field_id, const int l_tgt_field_id) {
    char *c_local_code_name, *c_cpl_id, *c_tgt_field_id;

    c_local_code_name = CWP_fortran_to_c_string(f_local_code_name, l_local_code_name);
    c_cpl_id = CWP_fortran_to_c_string(f_cpl_id, l_cpl_id);
    c_tgt_field_id = CWP_fortran_to_c_string(f_tgt_field_id, l_tgt_field_id);

    CWP_Field_irecv(c_local_code_name, c_cpl_id, c_tgt_field_id);

    free(c_local_code_name);
    free(c_cpl_id);
    free(c_tgt_field_id);
}

void CWP_Field_wait_issend_cf(const char *f_local_code_name, const int l_local_code_name, const char *f_cpl_id, const int l_cpl_id, const char *f_src_field_id, const int l_src_field_id) {
    char *c_local_code_name, *c_cpl_id, *c_src_field_id;

    c_local_code_name = CWP_fortran_to_c_string(f_local_code_name, l_local_code_name);
    c_cpl_id = CWP_fortran_to_c_string(f_cpl_id, l_cpl_id);
    c_src_field_id = CWP_fortran_to_c_string(f_src_field_id, l_src_field_id);

    CWP_Field_wait_issend(c_local_code_name, c_cpl_id, c_src_field_id);

    free(c_local_code_name);
    free(c_cpl_id);
    free(c_src_field_id);
}

void CWP_Field_wait_irecv_cf(const char *f_local_code_name, const int l_local_code_name, const char *f_cpl_id, const int l_cpl_id, const char *f_tgt_field_id, const int l_tgt_field_id) {
    char *c_local_code_name, *c_cpl_id, *c_tgt_field_id;

    c_local_code_name = CWP_fortran_to_c_string(f_local_code_name, l_local_code_name);
    c_cpl_id = CWP_fortran_to_c_string(f_cpl_id, l_cpl_id);
    c_tgt_field_id = CWP_fortran_to_c_string(f_tgt_field_id, l_tgt_field_id);

    CWP_Field_wait_irecv(c_local_code_name, c_cpl_id, c_tgt_field_id);

    free(c_local_code_name);
    free(c_cpl_id);
    free(c_tgt_field_id);
}
