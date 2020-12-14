#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "pdm.h"
#include "pdm_config.h"
#include "pdm_mpi.h"
#include "pdm_partitioning_algorithm.h"
#include "pdm_dmesh_partitioning.h"
#include "pdm_para_graph_dual.h"
#include "pdm_dcube_gen.h"
#include "pdm_printf.h"
#include "pdm_sort.h"
#include "pdm_distrib.h"
#include "pdm_error.h"
#include "pdm_priv.h"

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
 * \param [inout]   part_method Partitioner (1 ParMETIS, 2 Pt-Scotch)
 *
 */

static void
_read_args(int            argc,
           char         **argv,
           PDM_g_num_t  *n_vtx_seg,
           double        *length,
           int           *n_part,
           int           *post,
           int           *part_method)
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
      *part_method = 2;
    }
    else if (strcmp(argv[i], "-parmetis") == 0) {
      *part_method = 1;
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
  double             length    = 1.;
  int                n_part    = 1;
  int                post      = 0;
#ifdef PDM_HAVE_PARMETIS
  PDM_split_dual_t part_method  = PDM_SPLIT_DUAL_WITH_PARMETIS;
#else
#ifdef PDM_HAVE_PTSCOTCH
  PDM_split_dual_t part_method  = PDM_SPLIT_DUAL_WITH_PTSCOTCH;
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
     (int *) &part_method);

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
  PDM_g_num_t  *dface_cell = NULL;
  int          *dface_vtx_idx = NULL;
  PDM_g_num_t  *dface_vtx = NULL;
  double       *dvtx_coord = NULL;
  int          *dface_group_idx = NULL;
  PDM_g_num_t  *dface_group = NULL;
  int           dface_vtxL;
  int           dFaceGroupL;

  /*
   *  Create distributed cube
   */

  int          id;
  PDM_MPI_Comm     comm = PDM_MPI_COMM_WORLD;

  PDM_dcube_gen_init(&id,
                      comm,
                      n_vtx_seg,
                      length,
            		      0.,
		                  0.,
		                  0.);

  PDM_dcube_gen_dim_get(id,
                         &n_face_group,
                         &dn_cell,
                         &dn_face,
                         &dn_vtx,
                         &dface_vtxL,
                         &dFaceGroupL);

  PDM_dcube_gen_data_get(id,
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

  // printf(" PDM_PART_NONE::%d\n"     , PDM_PART_NONE);
  // printf(" PDM_PART_NULL::%d\n"     , PDM_PART_NULL);
  // printf(" PDM_PART_FACE_CELL::%d\n", PDM_PART_FACE_CELL);
  // printf(" PDM_PART_CELL_FACE::%d\n", PDM_PART_CELL_FACE);

  int flags = PDM_PART_FACE_CELL|PDM_PART_CELL_FACE;
  printf("PDM_HASFLAG(flags, PDM_PART_FACE_CELL) :: %d\n", PDM_HASFLAG(flags, PDM_PART_FACE_CELL) );
  printf("PDM_HASFLAG(flags, PDM_PART_CELL_FACE) :: %d\n", PDM_HASFLAG(flags, PDM_PART_CELL_FACE) );
  printf("PDM_HASFLAG(flags, PDM_PART_FACE_VTX) :: %d\n" , PDM_HASFLAG(flags, PDM_PART_FACE_VTX) );
  printf("x::PDM_HASFLAG(flags, PDM_PART_FACE_VTX) :: %x\n", PDM_PART_FACE_VTX);


  // PDM_dmesh_partitioning_part_get(1, 1, PDM_PART_FACE_CELL, NULL);
  // PDM_dmesh_partitioning_part_get(1, 1, PDM_PART_CELL_FACE, NULL);
  // PDM_dmesh_partitioning_part_get(1, 1, PDM_PART_FACE_VTX, NULL);
  // PDM_dmesh_partitioning_part_get(1, 1, PDM_PART_FACE_VTX|PDM_PART_FACE_CELL, NULL);

  // PDM_dmesh_partitioning_get(1, PDM_PART_FACE_CELL, NULL);
  // PDM_dmesh_partitioning_get(1, PDM_PART_CELL_FACE, NULL);
  // PDM_dmesh_partitioning_get(1, PDM_PART_FACE_VTX, NULL);
  // PDM_dmesh_partitioning_get(1, PDM_PART_FACE_VTX|PDM_PART_FACE_CELL, NULL);

  // int ppart_id = 0;

  gettimeofday(&t_elaps_debut, NULL);

  /*
   *  Create mesh partitions
   */
  PDM_g_num_t* cell_distribution = PDM_compute_entity_distribution(comm, dn_cell);
  PDM_g_num_t* face_distribution = PDM_compute_entity_distribution(comm, dn_face);
  PDM_g_num_t* part_distribution = PDM_compute_entity_distribution(comm, n_part );
  PDM_g_num_t* vtx_distribution  = PDM_compute_entity_distribution(comm, dn_vtx );

  // printf("part_distribution::\n");
  // for(int i_part = 0; i_part < n_rank+1; ++i_part){
  //   printf("%d ", part_distribution[i_part]);
  // }
  // printf("\n");

  /*
   * Compute dual graph
   */
  PDM_g_num_t* dual_graph_idx;
  PDM_g_num_t* dual_graph;
  int* dcell_face_idx;
  PDM_g_num_t* dcell_face;
  PDM_para_graph_dual_from_arc2node(comm,
                                    cell_distribution,
                                    face_distribution,
                                    dface_cell,
                    (PDM_g_num_t**) &dual_graph_idx,
                    (PDM_g_num_t**) &dual_graph,
                                    1,
                    (int        **) &dcell_face_idx,
                    (PDM_g_num_t**) &dcell_face);

  // Test graph creation from cell_face connectivity
  if (0 == 1) {

    if (0 == 1) {
      printf("dcell_face_idx :");
      for (int i =0; i < dn_cell+1; i++)
        printf(" %d", dcell_face_idx[i]);
      printf("\n");

      printf("dcell_face:: %d \n", dn_cell);
      for(int i = 0; i < dn_cell; ++i){
        printf("Local cell %d :", i);
        for(int j = dcell_face_idx[i]; j < dcell_face_idx[i+1]; ++j){
          printf(" "PDM_FMT_G_NUM"", dcell_face[j]);
        }
        PDM_printf("\n");
      }
    }

    free(dual_graph_idx);
    free(dual_graph);
    dual_graph_idx = NULL;
    dual_graph = NULL;

    PDM_para_graph_dual_from_node2arc(comm,
                                      cell_distribution,
                                      face_distribution,
                                      dcell_face_idx,
                                      dcell_face,
                      (PDM_g_num_t**) &dual_graph_idx,
                      (PDM_g_num_t**) &dual_graph);
  }

  // PDM_log_trace_array_long(dual_graph_idx, dn_cell+1              , "pdm_t_partitioning_dcube::dual_graph_idx::");
  // PDM_log_trace_array_long(dual_graph    , dual_graph_idx[dn_cell], "pdm_t_partitioning_dcube::dual_graph::");

  free(dual_graph_idx);
  free(dual_graph);
  // mpirun -np 2 ./paradigm/test/pdm_t_partitioning_dcube -n 23 -n_part 1 -parmetis
  PDM_para_graph_dual_from_combine_connectivity(comm,
                                                cell_distribution,
                                                face_distribution,
                                                vtx_distribution,
                                                dcell_face_idx,
                                                dcell_face,
                                                dface_vtx_idx,
                                                dface_vtx,
                               (PDM_g_num_t**) &dual_graph_idx,
                               (PDM_g_num_t**) &dual_graph);

  /*
   * Split it !!! CAUTION dn_cell can be different of the size of dual graph !!!
   */
  // printf("PDM_split_graph\n");
  int* cell_part    = (int *) malloc( sizeof(int) * dn_cell );
  int* dcell_weight = (int *) malloc( sizeof(int) * dn_cell );
  for(int i = 0; i < dn_cell; ++i){
    dcell_weight[i] = dual_graph_idx[i+1] - dual_graph_idx[i];
  }

  if( 0 == 1 ){
    printf("n_cell_block:: %d \n", dn_cell);
    for(int i = 0; i < dn_cell; ++i){
      printf(" dual_graph_idx = "PDM_FMT_G_NUM" ---> \n", dual_graph_idx[i]);
      for(int i_data = dual_graph_idx[i]; i_data < dual_graph_idx[i+1]; ++i_data){
        // printf("%d ", dual_graph[i_data]);
        printf("\t dual_graph[%d] = "PDM_FMT_G_NUM" \n", i_data, dual_graph[i_data]);
      }
      printf("\n");
    }
  }

  int tn_part = part_distribution[n_rank]-1;

  double *part_frac = NULL;
  if (0 == 0) {
  part_frac = (double *) malloc(sizeof(double) * tn_part );
  for (int i_part = 0; i_part < tn_part-1; i_part++)
  {
    if (i_part % 2 == 0) part_frac[i_part] = (double) 0.5*(1./tn_part);
    else                 part_frac[i_part] = (double) 1.5*(1./tn_part);
  }
  if (tn_part % 2 == 0) part_frac[tn_part-1] = (double) 1.5*(1./tn_part);
  else                  part_frac[tn_part-1] = (double) 1.*(1./tn_part);
  PDM_printf("Testing with heterogeneous part sizes\n");
  }

  PDM_para_graph_split (part_method,
                        cell_distribution,
                        dual_graph_idx,
                        dual_graph,
                        NULL,
                        NULL,
                        tn_part,
                        part_frac, // Or NULL for homogeneous parts
                        cell_part,
                        comm);

  if (0==1){
    printf("cell_part[%d]::", dn_cell);
    for(int i = 0; i < dn_cell; ++i){
      printf("%d ", cell_part[i]);
    }
    printf("\n");
  }
  if (part_frac != NULL)
    free(part_frac);

  /*
   * On dispose pour chaque cellule de la partition associé : il faut retrouver le
   *  cell_ln_to_gn
   *    Idée : faire un part_to_block avec ln_to_gn == cell_part et distribution imposé = dpart_proc
   *    On devrait recupérer n_part block de données avec le ln_to_gn !
   */
  PDM_g_num_t** pcell_ln_to_gn;
  int*  pn_cell;

  int n_res_part = PDM_part_assemble_partitions(comm,
                                                part_distribution,
                                                cell_distribution,
                                                cell_part,
                                    (int ** )  &pn_cell,
                            (PDM_g_num_t ***)  &pcell_ln_to_gn);

  /*
   * Tentative extented partition
   */
  // PDM_g_num_t** pcell_ln_to_gn_extented;
  // int*          pn_cell_extented;

  // PDM_extend_mesh(comm,
  //                 part_distribution,
  //                 cell_distribution,
  //                 cell_part,
  //                 n_res_part,
  //                 dual_graph_idx,
  //                 dual_graph,
  //                 pn_cell,
  //                 pcell_ln_to_gn,
  //                &pn_cell_extented,
  //                &pcell_ln_to_gn_extented);

  // pn_cell        = pn_cell_extented;
  // pcell_ln_to_gn = pcell_ln_to_gn_extented;
  /*
   *  At this stage we have the cell_ln_to_gn :
   *      --> We need to deduce the other if needed
   *           --> face_ln_to_gn
   *           --> vtx_ln_to_gn
   *  We can do by the relationship of each varibles (connectivity)
   *  For example face_ln_to_gn can be deduce with cell_face connectivity if we have cell_ln_to_gn
   */
  int** pcell_face_idx;
  int** pcell_face;
  int*  pn_faces;
  PDM_g_num_t** pface_ln_to_gn = NULL;

  PDM_part_dconnectivity_to_pconnectivity_sort(comm,
                                               cell_distribution,
                                               dcell_face_idx,
                                               dcell_face,
                                               n_res_part,
                                               pn_cell,
                      (const PDM_g_num_t ** )  pcell_ln_to_gn,
                            (int         ** ) &pn_faces,
                            (PDM_g_num_t ***) &pface_ln_to_gn,
                            (int         ***) &pcell_face_idx,
                            (int         ***) &pcell_face);

  int **pface_cell;
  PDM_part_reverse_pcellface(n_res_part,
                             pn_cell,
                             pn_faces,
            (const int **)   pcell_face_idx,
             (const int **)  pcell_face,
              (int    ***)  &pface_cell);

  if (0 == 1){
    for (int i_part=0; i_part < n_res_part; i_part++){
      PDM_printf("[%i] generated facecell part %i [%i]:", i_rank, i_part, pn_faces[i_part]);
      for (int iface=0 ; iface < pn_faces[i_part]; iface++)
        PDM_printf(" %d %d", pface_cell[i_part][2*iface], pface_cell[i_part][2*iface+1]);
      PDM_printf("\n");
    }
  }


  /*
   * Generate vtx
   */

  int** pface_vtx_idx;
  int** pface_vtx;
  int*  pn_vtx;
  PDM_g_num_t** pvtx_ln_to_gn = NULL;

  PDM_part_dconnectivity_to_pconnectivity_sort(comm,
                                               face_distribution,
                                               dface_vtx_idx,
                                               dface_vtx,
                                               n_res_part,
                                               pn_faces,
                        (const PDM_g_num_t **) pface_ln_to_gn,
                           (int         ** )  &pn_vtx,
                           (PDM_g_num_t ***)  &pvtx_ln_to_gn,
                           (int         ***)  &pface_vtx_idx,
                           (int         ***)  &pface_vtx);

  double **pvtx_coord = NULL;
  PDM_part_dcoordinates_to_pcoordinates(comm,
                                        n_part,
                                        vtx_distribution,
                                        dvtx_coord,
                                        pn_vtx,
                (const PDM_g_num_t **)  pvtx_ln_to_gn,
                          (double ***) &pvtx_coord);
  /*
   * On doit calculer le dcell_face car avec lui et le cell_ln_to_gn on retrouve facilement
   *   le dcell_face sur la partition donc on n'a plus qu'a trié pour avoir le face_ln_to_gn
   *      Une fois le face_ln_to_gn trouvé on doit avoir facilement
   *      le dface_cell mais en numéro absolu (attention au signe ) et la on a juste à chercher
   *      dans le cell_ln_to_gn l'index correspondant pour transformé le numero global en numéro local
   */

  /*
   *  L'astuce de tri pour le face_cell et le meme que pour le face_vtx -->
   *    permet d'avoir le vtx_ln_to_gn
   */

  /*
   *  Boundary condition (face group )
   */
  PDM_g_num_t** pface_group_ln_to_gn;
  int** pface_group;
  int** pface_group_idx;

  PDM_part_distgroup_to_partgroup(comm,
                                  face_distribution,
                                  n_face_group,
                                  dface_group_idx,
                                  dface_group,
                                  n_res_part,
                                  pn_faces,
           (const PDM_g_num_t **) pface_ln_to_gn,
                (int         ***) &pface_group_idx,
                (int         ***) &pface_group,
                (PDM_g_num_t ***) &pface_group_ln_to_gn);

  PDM_part_reorient_bound_faces(n_part,
                                pn_faces,
                                pface_cell,
                 (const int **) pcell_face_idx,
                      ( int **) pcell_face,
                 (const int **) pface_vtx_idx,
                       (int **) pface_vtx);

  /*
   * Graph communication build
   */
  int** pproc_face_bound_idx;
  int** ppart_face_bound_idx;
  int** pface_bound;

  int **face_is_bnd = (int **) malloc(n_part * sizeof(int*));
  for (int i_part = 0; i_part < n_res_part; i_part++) {
    face_is_bnd[i_part] = (int *) malloc(pn_faces[i_part]*sizeof(int));
    for (int i_face = 0; i_face < pn_faces[i_part]; i_face++){
      if (pface_cell[i_part][2*i_face+1] > 0)
        face_is_bnd[i_part][i_face] = 0;
      else
        face_is_bnd[i_part][i_face] = 1;
    }
  }

  PDM_part_generate_entity_graph_comm(comm,
                                      part_distribution,
                                      face_distribution,
                                      n_part,
                                      pn_faces,
              (const PDM_g_num_t **)  pface_ln_to_gn,
              (const int **)          face_is_bnd,
                           (int ***) &pproc_face_bound_idx,
                           (int ***) &ppart_face_bound_idx,
                           (int ***) &pface_bound);
  for (int i_part = 0; i_part < n_res_part; i_part++)
    free(face_is_bnd[i_part]);
  free(face_is_bnd);


  // Attention on veut garder l'orientation donc il y a un signe dans le face_cell / cell_face
  // Reflechir sur les connectivité d'edge également ...

  /*
   * Free
   */
  free(dual_graph_idx);
  free(dual_graph);
  free(cell_part);
  free(dcell_face);
  free(dcell_face_idx);
  free(dcell_weight);
  free(cell_distribution);
  free(face_distribution);
  free(part_distribution);
  free(vtx_distribution);
  for(int i_part = 0; i_part < n_res_part; ++i_part){
    free(pface_ln_to_gn[i_part]);
    free(pcell_ln_to_gn[i_part]);
    // free(pcell_ln_to_gn_extented[i_part]);
    free(pvtx_ln_to_gn[i_part]);
    free(pcell_face[i_part]);
    free(pcell_face_idx[i_part]);
    free(pface_vtx_idx[i_part]);
    free(pface_vtx[i_part]);
    free(pface_group_ln_to_gn[i_part]);
    free(pface_group[i_part]);
    free(pface_group_idx[i_part]);
    free(pproc_face_bound_idx[i_part]);
    free(ppart_face_bound_idx[i_part]);
    free(pface_bound[i_part]);
    free(pface_cell[i_part]);
    free(pvtx_coord[i_part]);
  }
  free(pcell_face);
  free(pcell_face_idx);
  free(pvtx_ln_to_gn);
  free(pface_vtx_idx);
  free(pface_vtx);
  free(pcell_ln_to_gn);
  // free(pcell_ln_to_gn_extented);
  free(pproc_face_bound_idx);
  free(ppart_face_bound_idx);
  free(pface_bound);
  free(pface_ln_to_gn);
  free(pn_cell);
  // free(pn_cell_extented);
  free(pn_faces);
  free(pn_vtx);
  free(pface_group_ln_to_gn);
  free(pface_group);
  free(pface_group_idx);
  free(pface_cell);
  free(pvtx_coord);



  PDM_dcube_gen_free(id);

  PDM_MPI_Finalize();

  return 0;
}


  // Unit test ...
  // int test_unique[10] = {10, 2, 1, 12, 31, 2, 31, 4, 5, 2};
  // // int ns = PDM_inplace_unique_long(test_unique, 4, 10);
  // // PDM_quick_sort_long(test_unique, 0, 9);
  // int ns = PDM_inplace_unique_long(test_unique, 0, 9);
  // printf("ns::%d\n", ns);
  // printf("test_unique::");
  // for(int i = 0; i < ns; ++i){
  //   printf("%d ", test_unique[i]);
  // }
  // printf("\n");

  // Unit test ...
  // int cell_cell_n[4] = {4, 5, 3, 2};
  // PDM_g_num_t dual_graph_test[14] = {31, 10, 5, 6,
  //                                    8, 8, 2, 1, 4,
  //                                    2, 3, 2,
  //                                    91, 92 };
  // int dual_graph_test_idx[5];

  // PDM_compress_connectivity(dual_graph_test, dual_graph_test_idx, cell_cell_n, 4);


  // if( 1 == 1 ){
  //   // printf("n_cell_block:: %d \n", n_cell_block);
  //   for(int i = 0; i < 4; ++i){
  //     printf(" dual_graph_test_idx = %d ---> \n", dual_graph_test_idx[i]);
  //     for(int i_data = dual_graph_test_idx[i]; i_data < dual_graph_test_idx[i+1]; ++i_data){
  //       // printf("%d ", _dual_graph[i_data]);
  //       printf("\t _dual_graph[%d] = %d \n", i_data, dual_graph_test[i_data]);
  //     }
  //     printf("\n");
  //   }
  // }
