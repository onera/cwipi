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
#include "pdm_config.h"
#include "pdm_mpi.h"
#include "pdm_part.h"
#include "pdm_part_extension.h"
#include "pdm_dcube_gen.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_cellface_orient.h"

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/**
 *
 * \brief  Usage
 *
 */

static void
_usage(int exit_code)
{
  PDM_printf
    ("\n"
     "  Usage: \n\n"
     "  -n      <level>  Number of vertices on the cube side.\n\n"
     "  -l      <level>  Cube length.\n\n"
     "  -n_part <level>  Number of partitions par process.\n\n"
     "  -post            Ensight outputs (only if n_part == 1). \n\n"
     "  -parmetis        Call ParMETIS.\n\n"
     "  -pt-scocth       Call PT-Scotch.\n\n"
     "  -h               This message.\n\n");

  exit(exit_code);
}


/**
 *
 * \brief  Read arguments from the command line
 *
 * \param [in]      argc     Number of arguments
 * \param [in]      argv     Arguments
 * \param [inout]   n_vtx_seg  Number of vertices on the cube side
 * \param [inout]   length   Cube length
 * \param [inout]   n_part   Number of partitions par process
 * \param [inout]   post     Ensight outputs status
 * \param [inout]   method   Partitioner (1 ParMETIS, 2 Pt-Scotch)
 *
 */

static void
_read_args(int            argc,
           char         **argv,
           PDM_g_num_t  *n_vtx_seg,
           double        *length,
           int           *n_part,
	   int           *post,
	   int           *method)
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {

    if (strcmp(argv[i], "-h") == 0)
      _usage(EXIT_SUCCESS);

    else if (strcmp(argv[i], "-n") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        long _n_vtx_seg = atol(argv[i]);
        *n_vtx_seg = (PDM_g_num_t) _n_vtx_seg;
      }
    }
    else if (strcmp(argv[i], "-l") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *length = atof(argv[i]);
    }
    else if (strcmp(argv[i], "-n_part") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *n_part = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-post") == 0) {
      *post = 1;
    }
    else if (strcmp(argv[i], "-pt-scotch") == 0) {
      *method = 2;
    }
    else if (strcmp(argv[i], "-parmetis") == 0) {
      *method = 1;
    }
    else
      _usage(EXIT_FAILURE);
    i++;
  }
}


/**
 *
 * \brief  Main
 *
 */

int main(int argc, char *argv[])
{

  /*
   *  Set default values
   */

  PDM_g_num_t        n_vtx_seg = 10;
  double             length  = 1.;
  int                n_part   = 1;
  int                post    = 0;
#ifdef PDM_HAVE_PARMETIS
  PDM_part_split_t method  = PDM_PART_SPLIT_PARMETIS;
#else
#ifdef PDM_HAVE_PTSCOTCH
  PDM_part_split_t method  = PDM_PART_SPLIT_PTSCOTCH;
#endif
#endif

  /*
   *  Read args
   */

  _read_args(argc,
             argv,
             &n_vtx_seg,
             &length,
             &n_part,
             &post,
             (int *) &method);

  /*
   *  Init
   */

  struct timeval t_elaps_debut;

  int i_rank;
  int n_rank;

  PDM_MPI_Init(&argc, &argv);
  PDM_MPI_Comm_rank(PDM_MPI_COMM_WORLD, &i_rank);
  PDM_MPI_Comm_size(PDM_MPI_COMM_WORLD, &n_rank);

  int           dn_cell;
  int           dn_face;
  int           dn_vtx;
  int           n_face_group;
  PDM_g_num_t  *dface_cell      = NULL;
  int          *dface_vtx_idx   = NULL;
  PDM_g_num_t  *dface_vtx       = NULL;
  double       *dvtx_coord      = NULL;
  int          *dface_group_idx = NULL;
  PDM_g_num_t  *dface_group     = NULL;
  int           dface_vtxL;
  int           dFaceGroupL;

  /*
   *  Create distributed cube
   */

  PDM_MPI_Comm     comm = PDM_MPI_COMM_WORLD;

  PDM_dcube_t* dcube = PDM_dcube_gen_init(comm,
                                          n_vtx_seg,
                                          length,
                                          0.,
                                          0.,
                                          0.,
                                          PDM_OWNERSHIP_KEEP);

  PDM_dcube_gen_dim_get(dcube,
                         &n_face_group,
                         &dn_cell,
                         &dn_face,
                         &dn_vtx,
                         &dface_vtxL,
                         &dFaceGroupL);

  PDM_dcube_gen_data_get(dcube,
                          &dface_cell,
                          &dface_vtx_idx,
                          &dface_vtx,
                          &dvtx_coord,
                          &dface_group_idx,
                          &dface_group);

  if (0 == 1) {

    PDM_printf("[%i] n_face_group    : %i\n", i_rank, n_face_group);
    PDM_printf("[%i] dn_cell        : %i\n", i_rank, dn_cell);
    PDM_printf("[%i] dn_face        : %i\n", i_rank, dn_face);
    PDM_printf("[%i] dn_vtx         : %i\n", i_rank, dn_vtx);

    PDM_printf("[%i] dface_cell     : ", i_rank);
    for (int i = 0; i < 2 * dn_face; i++)
      PDM_printf(" "PDM_FMT_G_NUM, dface_cell[i]);
    PDM_printf("\n");

    PDM_printf("[%i] dface_vtx_idx   : ", i_rank);
    for (int i = 0; i < dn_face + 1; i++)
      PDM_printf(" %i", dface_vtx_idx[i]);
    PDM_printf("\n");

    PDM_printf("[%i] dface_vtx      : ", i_rank);
    for (int i = 0; i < dface_vtx_idx[dn_face]; i++)
      PDM_printf(" "PDM_FMT_G_NUM, dface_vtx[i]);
    PDM_printf("\n");

    PDM_printf("[%i] dvtx_coord     : ", i_rank);
    for (int i = 0; i < 3*dn_vtx; i++)
      PDM_printf(" %12.5e", dvtx_coord[i]);
    PDM_printf("\n");

    PDM_printf("[%i] dface_group_idx : ", i_rank);
    for (int i = 0; i < n_face_group + 1; i++)
      PDM_printf(" %i", dface_group_idx[i]);
    PDM_printf("\n");

    PDM_printf("[%i] dface_group    : ", i_rank);
    for (int i = 0; i < dface_group_idx[n_face_group]; i++)
      PDM_printf(" "PDM_FMT_G_NUM, dface_group[i]);
    PDM_printf("\n");

  }
  int ppart_id = 0;

  gettimeofday(&t_elaps_debut, NULL);

  /*
   *  Create mesh partitions
   */

  int have_dcell_part = 0;

  int *dcell_part = (int *) malloc(dn_cell*sizeof(int));
  int *renum_properties_cell = NULL;
  int *renum_properties_face = NULL;
  int n_property_cell = 0;
  int n_property_face = 0;

  PDM_part_create(&ppart_id,
                  comm,
                  method,
                  "PDM_PART_RENUM_CELL_NONE",
                  "PDM_PART_RENUM_FACE_NONE",
                  n_property_cell,
                  renum_properties_cell,
                  n_property_face,
                  renum_properties_face,
                  n_part,
                  dn_cell,
                  dn_face,
                  dn_vtx,
                  n_face_group,
                  NULL,
                  NULL,
                  NULL,
                  NULL,
                  have_dcell_part,
                  dcell_part,
                  dface_cell,
                  dface_vtx_idx,
                  dface_vtx,
                  NULL,
                  dvtx_coord,
                  NULL,
                  dface_group_idx,
                  dface_group);

  PDM_part_extension_t* part_ext = PDM_part_extension_create(1,
                                                             &n_part,
                                                             PDM_EXTEND_FROM_FACE,
                                                             1,
                                                             comm,
                                                             PDM_OWNERSHIP_KEEP);

  for (int i_part = 0; i_part < n_part; i_part++) {

    int n_cell;
    int n_face;
    int n_face_part_bound;
    int n_vtx;
    int n_proc;
    int n_total_part;
    int scell_face;
    int sface_vtx;
    int sface_group;
    int n_face_group2;

    PDM_part_part_dim_get(ppart_id,
                          i_part,
                          &n_cell,
                          &n_face,
                          &n_face_part_bound,
                          &n_vtx,
                          &n_proc,
                          &n_total_part,
                          &scell_face,
                          &sface_vtx,
                          &sface_group,
                          &n_face_group2);

    int          *cell_tag;
    int          *cell_face_idx;
    int          *cell_face;
    PDM_g_num_t  *cell_ln_to_gn;
    int          *face_tag;
    int          *face_cell;
    int          *face_vtx_idx;
    int          *face_vtx;
    PDM_g_num_t  *face_ln_to_gn;
    int          *face_part_bound_proc_idx;
    int          *face_part_bound_part_idx;
    int          *face_part_bound;
    int          *vtx_tag;
    double       *vtx;
    PDM_g_num_t  *vtx_ln_to_gn;
    int          *face_group_idx;
    int          *face_group;
    PDM_g_num_t  *face_group_ln_to_gn;

    PDM_part_part_val_get(ppart_id,
                          i_part,
                          &cell_tag,
                          &cell_face_idx,
                          &cell_face,
                          &cell_ln_to_gn,
                          &face_tag,
                          &face_cell,
                          &face_vtx_idx,
                          &face_vtx,
                          &face_ln_to_gn,
                          &face_part_bound_proc_idx,
                          &face_part_bound_part_idx,
                          &face_part_bound,
                          &vtx_tag,
                          &vtx,
                          &vtx_ln_to_gn,
                          &face_group_idx,
                          &face_group,
                          &face_group_ln_to_gn);

    PDM_cellface_orient (n_cell,
                         n_face,
                         n_vtx,
                         vtx,
                         cell_face_idx,
                         cell_face,
                         face_cell,
                         face_vtx_idx,
                         face_vtx);

    PDM_part_extension_set_part(part_ext, 0, i_part,
                                n_cell,
                                n_face,
                                n_face_part_bound,
                                n_face_group,
                                n_vtx,
                                cell_face_idx,
                                cell_face,
                                face_cell,
                                NULL, // face_edge_idx
                                NULL, // face_edge
                                face_vtx_idx,
                                face_vtx,
                                NULL, //edge_vtx
                                face_group_idx,
                                face_group,
                                NULL, // face_join_idx
                                NULL, // face_join
                                face_part_bound_proc_idx,
                                face_part_bound_part_idx,
                                face_part_bound,
                                cell_ln_to_gn,
                                face_ln_to_gn,
                                NULL, // edge_ln_to_gn
                                vtx_ln_to_gn,
                                face_group_ln_to_gn,
                                vtx);

    // PDM_log_trace_array_long (cell_ln_to_gn, n_cell, "cell_ln_to_gn::");
    if( 0 == 1) {
      PDM_printf("[%i] n_face_group     : %i\n", i_rank, n_face_group);
      PDM_printf("[%i] n_cell          : %i\n", i_rank, n_cell);
      PDM_printf("[%i] n_face          : %i\n", i_rank, n_face);
      PDM_printf("[%i] n_vtx           : %i\n", i_rank, n_vtx);
      PDM_printf("[%i] n_face_part_bound : %i\n", i_rank, n_face_part_bound);

      PDM_printf("[%i] cell_face     : ", i_rank);
      for (int i = 0; i < n_cell; i++) {
        for (int j = cell_face_idx[i]; j < cell_face_idx[i+1]; j++) {
          PDM_printf(" %i", cell_face[j]);
        }
        PDM_printf("\n");
      }

      PDM_printf("\n");

      PDM_printf("[%i]  face_part_bound    : ", i_rank);
      for (int i = 0; i < 4 * n_face_part_bound; i++)
        PDM_printf(" "PDM_FMT_G_NUM, face_part_bound[i]);
      PDM_printf("\n");

      PDM_printf("[%i]  cell_ln_to_gn    : ", i_rank);
      for (int i = 0; i < n_cell; i++)
        PDM_printf(" "PDM_FMT_G_NUM, cell_ln_to_gn[i]);
      PDM_printf("\n");

      PDM_printf("[%i] face_cell     : ", i_rank);
      for (int i = 0; i < 2 * n_face; i++)
        PDM_printf(" %i", face_cell[i]);
      PDM_printf("\n");

      PDM_printf("[%i] face_vtx      : ", i_rank);
      for (int i = 0; i < n_face; i++) {
        for (int j = face_vtx_idx[i]; j < face_vtx_idx[i+1]; j++) {
          PDM_printf(" %i", face_vtx[j]);
        }
        PDM_printf("\n");
      }

      PDM_printf("[%i]  face_ln_to_gn    : ", i_rank);
      for (int i = 0; i < n_face; i++)
        PDM_printf(" "PDM_FMT_G_NUM, face_ln_to_gn[i]);
      PDM_printf("\n");

      PDM_printf("[%i] vtx           : ", i_rank);
      for (int i = 0; i < 3 * n_vtx; i++)
        PDM_printf(" %12.5e", vtx[i]);
      PDM_printf("\n");

      PDM_printf("[%i] vtx_ln_to_gn     : ", i_rank);
      for (int i = 0; i <  n_vtx; i++)
        PDM_printf(" "PDM_FMT_G_NUM, vtx_ln_to_gn[i]);
      PDM_printf("\n");

      PDM_printf("[%i] face_group_idx : ", i_rank);
      for (int i = 0; i < n_face_group + 1; i++)
        PDM_printf(" %i", face_group_idx[i]);
      PDM_printf("\n");

      PDM_printf("[%i] face_group    : ", i_rank);
      for (int i = 0; i < n_face_group; i++) {
        for (int j = face_group_idx[i]; j < face_group_idx[i+1]; j++) {
          PDM_printf(" %i", face_group[j]);
        }
        PDM_printf("\n");
      }

      PDM_printf("[%i] face_group_ln_to_gn   : ", i_rank);
      for (int i = 0; i < n_face_group; i++) {
        for (int j = face_group_idx[i]; j < face_group_idx[i+1]; j++) {
          PDM_printf(" "PDM_FMT_G_NUM, face_group_ln_to_gn[j]);
        }
        PDM_printf("\n");
      }
    }
  }

  PDM_part_extension_compute(part_ext);

  for (int i_part = 0; i_part < n_part; i_part++) {

    double* vtx_coord_extended;
    PDM_g_num_t* border_vtx_ln_to_gn;
    PDM_g_num_t* border_cell_ln_to_gn;
    int n_vtx_extended = PDM_part_extension_coord_get(part_ext, 0, i_part, &vtx_coord_extended);
    int n_vtx_extended2 = PDM_part_extension_ln_to_gn_get(part_ext, 0, i_part, PDM_MESH_ENTITY_VERTEX, &border_vtx_ln_to_gn);
    int n_cell_extended = PDM_part_extension_ln_to_gn_get(part_ext, 0, i_part, PDM_MESH_ENTITY_CELL, &border_cell_ln_to_gn);
    assert(n_vtx_extended == n_vtx_extended2);
    if(0 == 1) {
      for(int i_vtx = 0; i_vtx < n_vtx_extended; ++i_vtx) {
        printf("[%i] vtx_coord_extended[%i] = %12.5e %12.5e %12.5e "PDM_FMT_G_NUM" \n", i_part, i_vtx, vtx_coord_extended[3*i_vtx], vtx_coord_extended[3*i_vtx+1], vtx_coord_extended[3*i_vtx+2], border_vtx_ln_to_gn[i_vtx]);
      }

      for(int i_cell = 0; i_cell < n_cell_extended; ++i_cell) {
        printf("[%i] border_cell_ln_to_gn[%i] = "PDM_FMT_G_NUM" \n", i_part, i_cell, border_cell_ln_to_gn[i_cell]);
      }
    }
  }

  PDM_part_extension_free(part_ext);


  free(dcell_part);
  PDM_part_free(ppart_id);

  PDM_dcube_gen_free(dcube);

  PDM_MPI_Finalize();

  return 0;
}
