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
#include "pdm_multipart.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_distrib.h"
#include "pdm_array.h"
#include "pdm_mesh_nodal.h"
#include "pdm_part_extension.h"
#include "pdm_dmesh_nodal_elements_utils.h"
#include "pdm_vtk.h"
#include "pdm_dmesh_nodal.h"
#include "pdm_logging.h"

#include "pdm_dcube_nodal_gen.h"
#include "pdm_domain_interface.h"
#include "pdm_dmesh_nodal_to_dmesh.h"

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
           PDM_g_num_t   *nx,
           PDM_g_num_t   *ny,
           PDM_g_num_t   *nz,
           int           *n_dom_i,
           int           *n_dom_j,
           int           *n_dom_k,
           int           *periodic_i,
           int           *periodic_j,
           int           *periodic_k,
           int           *t_elt,
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
        long _n = atol(argv[i]);
        *nx = (PDM_g_num_t) _n;
        *ny = (PDM_g_num_t) _n;
        *nz = (PDM_g_num_t) _n;
      }
    }
    else if (strcmp(argv[i], "-nx") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        long _n = atol(argv[i]);
        *nx = (PDM_g_num_t) _n;
      }
    }
    else if (strcmp(argv[i], "-ny") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        long _n = atol(argv[i]);
        *ny = (PDM_g_num_t) _n;
      }
    }
    else if (strcmp(argv[i], "-nz") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        long _n = atol(argv[i]);
        *nz = (PDM_g_num_t) _n;
      }
    }
    else if (strcmp(argv[i], "-ni") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *n_dom_i = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-nj") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *n_dom_j = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-nk") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *n_dom_k = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-pi") == 0) {
      *periodic_i = 1;
    }
    else if (strcmp(argv[i], "-pj") == 0) {
      *periodic_j = 1;
    }
    else if (strcmp(argv[i], "-pk") == 0) {
      *periodic_k = 1;
    }
    else if (strcmp(argv[i], "-t") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *t_elt = atoi(argv[i]);
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
// mpirun -np 1 ./test/pdm_t_join_domains -nx 3 -ny 3 -nz 2 -pi -t 8
// mpirun -np 1 ./test/pdm_t_join_domains -nx 3 -ny 3 -nz 2 -t 8 -ni 2 -pi 1
// mpirun -np 1 ./test/pdm_t_join_domains -nx 3 -ny 3 -nz 2 -t 8 -pi
// mpirun -np 1 ./test/pdm_t_join_domains -nx 3 -ny 3 -nz 2 -t 8 -pi -pj

int main
(
 int   argc,
 char *argv[]
 )
{
  /*
   *  Set default values
   */
  PDM_g_num_t          nx         = 10;
  PDM_g_num_t          ny         = 10;
  PDM_g_num_t          nz         = 10;
  int                  n_dom_i    = 1;
  int                  n_dom_j    = 1;
  int                  n_dom_k    = 1;
  int                  periodic_i = 0;
  int                  periodic_j = 0;
  int                  periodic_k = 0;
  double               length     = 1.;
  int                  n_part     = 1;
  int                  post       = 0;
  // PDM_Mesh_nodal_elt_t t_elt      = PDM_MESH_NODAL_TRIA3;
  PDM_Mesh_nodal_elt_t t_elt      = PDM_MESH_NODAL_HEXA8;
  // 2 -> tria
  // 3 -> quad
  // 5 -> tetra
  // 6 -> pyramid
  // 7 -> prism
  // 8 -> hexa

#ifdef PDM_HAVE_PARMETIS
  PDM_split_dual_t method  = PDM_SPLIT_DUAL_WITH_PARMETIS;
#else
#ifdef PDM_HAVE_PTSCOTCH
  PDM_split_dual_t method  = PDM_SPLIT_DUAL_WITH_PTSCOTCH;
#endif
#endif

  /*
   *  Read args
   */

  _read_args(argc,
             argv,
             &nx,
             &ny,
             &nz,
             &n_dom_i,
             &n_dom_j,
             &n_dom_k,
             &periodic_i,
             &periodic_j,
             &periodic_k,
             (int *) &t_elt,
             &length,
             &n_part,
             &post,
             (int *) &method);

  /*
   *  Init
   */
  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;
  int i_rank;
  int n_rank;

  PDM_MPI_Init(&argc, &argv);
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);



  /*
   *  Initialize structs
   */
  int dim = PDM_Mesh_nodal_elt_dim_get(t_elt);

  if (dim == 2) {
    n_dom_k    = 1;
    nz         = 1;
    periodic_k = 0;
  }

  // int n_interface =
  // n_dom_j*n_dom_k*(n_dom_i - 1 + periodic_i) +
  // n_dom_k*n_dom_i*(n_dom_j - 1 + periodic_j) +
  // n_dom_i*n_dom_j*(n_dom_k - 1 + periodic_k);

  int n_domain = n_dom_i * n_dom_j * n_dom_k;

  // PDM_dcube_nodal_t **dcube = (PDM_dcube_nodal_t **) malloc(sizeof(PDM_dcube_nodal_t *) * n_domain);
  PDM_dmesh_nodal_t **dmn   = (PDM_dmesh_nodal_t **) malloc(sizeof(PDM_dmesh_nodal_t *) * n_domain);
  PDM_dcube_nodal_t **dcube = NULL;
  PDM_domain_interface_t *dom_intrf = NULL;

  const int order = 1;

  PDM_dcube_nodal_cart_topo(comm,
                            n_dom_i,
                            n_dom_j,
                            n_dom_k,
                            periodic_i,
                            periodic_j,
                            periodic_k,
                            nx,
                            ny,
                            nz,
                            length,
                            0.,
                            0.,
                            0.,
                            t_elt,
                            order,
                            &dcube,
                            &dom_intrf,
                            PDM_OWNERSHIP_KEEP);
  // /* Check */
  // if (n_rank == 1) {

  //   double vect[4][3] = {
  //     {0.,             0.,             0.},
  //     {length*n_dom_i, 0.,             0.},
  //     {0,              length*n_dom_j, 0.},
  //     {0.,             0.,             length*n_dom_k}
  //   };

  //   double **all_coord = (double **) malloc(sizeof(double *) * n_domain);
  //   for (int i = 0; i < n_domain; i++) {
  //     all_coord[i] = PDM_DMesh_nodal_vtx_get(dmn[i]);
  //   }

  //   for (i_interface = 0; i_interface < n_interface; i_interface++) {

  //     int i_domain1 = interface_dom[i_interface][0];// - 1?
  //     int i_domain2 = interface_dom[i_interface][1];// - 1?

  //     for (int i = 0; i < interface_dn[i_interface]; i++) {

  //       PDM_g_num_t ivtx_1 = interface_ids[i_interface][2*i]   - 1;
  //       PDM_g_num_t ivtx_2 = interface_ids[i_interface][2*i+1] - 1;

  //       // if(ivtx_1 < 0 || ivtx_1 >= PDM_DMesh_nodal_n_vtx_get(dmn[i_domain1])) {
  //       //   printf("error : interface %d, domain 1:%d, pair %d, vtx "PDM_FMT_G_NUM"/%d\n",
  //       //          i_interface, i_domain1, i, ivtx_1+1,
  //       //          PDM_DMesh_nodal_n_vtx_get(dmn[i_domain1]));
  //       // }
  //       // if(ivtx_2 < 0 || ivtx_2 >= PDM_DMesh_nodal_n_vtx_get(dmn[i_domain2])) {
  //       //   printf("error : interface %d, domain 2:%d, pair %d, vtx "PDM_FMT_G_NUM"/%d\n",
  //       //          i_interface, i_domain2, i, ivtx_2+1,
  //       //          PDM_DMesh_nodal_n_vtx_get(dmn[i_domain2]));
  //       // }

  //       double *vtx_1 = all_coord[i_domain1] + 3*ivtx_1;
  //       double *vtx_2 = all_coord[i_domain2] + 3*ivtx_2;

  //       double dist2 = 0;
  //       for (int j = 0; j < 3; j++) {
  //         double delta = vtx_1[j] - (vtx_2[j] + vect[i_period[i_interface]][j]);
  //         dist2 += delta * delta;
  //       }

  //       // assert(dist2 < 1e-9);
  //       if (dist2 > 1e-9) {
  //         printf("!! interface %d, domains %d %d, vtx ("PDM_FMT_G_NUM": %f %f %f) ("PDM_FMT_G_NUM": %f %f %f), dist = %f\n",
  //                i_interface, i_domain1, i_domain2,
  //                ivtx_1+1,
  //                vtx_1[0], vtx_1[1], vtx_1[2],
  //                ivtx_2+1,
  //                vtx_2[0], vtx_2[1], vtx_2[2],
  //                sqrt(dist2));
  //       }
  //     }
  //   }
  //   free (all_coord);
  // }

  /*
   * Partitionnement
   */
  int* n_part_by_domain = (int *) malloc(n_domain * sizeof(int));
  for(int i = 0; i < n_domain; ++i) {
    n_part_by_domain[i] = n_part;
  }
  PDM_multipart_t *mpart_id = PDM_multipart_create(n_domain,
                                                   n_part_by_domain,
                                                   PDM_FALSE,
                                                   method,
                                                   PDM_PART_SIZE_HOMOGENEOUS,
                                                   NULL,
                                                   comm,
                                                   PDM_OWNERSHIP_KEEP);

  for (int i = 0; i < n_domain; i++) {
    dmn[i] = PDM_dcube_nodal_gen_dmesh_nodal_get(dcube[i]);
    PDM_dmesh_nodal_generate_distribution(dmn[i]);
    PDM_multipart_register_dmesh_nodal(mpart_id, i, dmn[i]);
  }

  PDM_multipart_run_ppart(mpart_id);



  /*
   * Create dmesh_nodal_to_dmesh to setup face and edge
   */
  // PDM_dmesh_nodal_to_dmesh_t* dmn_to_dm = PDM_dmesh_nodal_to_dmesh_create(n_domain, comm, PDM_OWNERSHIP_KEEP);
  // for (int i = 0; i < n_domain; i++) {
  //   PDM_dmesh_nodal_to_dmesh_add_dmesh_nodal(dmn_to_dm, i, dmn[i]);
  // }

  // PDM_dmesh_nodal_to_dmesh_compute(dmn_to_dm,
  //                                  PDM_DMESH_NODAL_TO_DMESH_TRANSFORM_TO_FACE,
  //                                  PDM_DMESH_NODAL_TO_DMESH_TRANSLATE_GROUP_TO_FACE);

  // PDM_dmesh_t **dm = malloc(n_domain * sizeof(PDM_dmesh_t *));

  // int          *dn_vtx        = malloc( n_domain * sizeof(int          ));
  // int          *dn_face       = malloc( n_domain * sizeof(int          ));
  // int         **dface_vtx_idx = malloc( n_domain * sizeof(int         *));
  // PDM_g_num_t **dface_vtx     = malloc( n_domain * sizeof(PDM_g_num_t *));

  // for (int i = 0; i < n_domain; i++) {
  //   PDM_dmesh_nodal_to_dmesh_get_dmesh(dmn_to_dm, 0, &dm[i]);

  //   int _dn_cell, _dn_face, _dn_edge, _dn_vtx, _n_bnd, _n_join;
  //   PDM_dmesh_dims_get(dm[i], &_dn_cell, &_dn_face, &_dn_edge, &_dn_vtx, &_n_bnd, &_n_join);

  //   dn_vtx [i] = _dn_vtx;
  //   dn_face[i] = PDM_dmesh_connectivity_get(dm[i],
  //                                           PDM_CONNECTIVITY_TYPE_FACE_VTX,
  //                                           &dface_vtx[i],
  //                                           &dface_vtx_idx[i],
  //                                           PDM_OWNERSHIP_KEEP);
  //   assert(dface_vtx_idx[i] != NULL);
  //   assert(dface_vtx[i] != NULL);
  // }

  // /*
  //  * Transform interface by vtx by interface by faces
  //  */
  // PDM_domain_interface_translate_vtx2face(dom_intrf,
  //                                         dn_vtx,
  //                                         dn_face,
  //                                         dface_vtx_idx,
  //                                         dface_vtx);


  /*
   *  Prepare pointer by domain and by part
   */
  int           *pn_n_part     = (int           *) malloc( n_domain * sizeof(int          *));
  int          **pn_vtx        = (int          **) malloc( n_domain * sizeof(int          *));
  PDM_g_num_t ***pvtx_ln_to_gn = (PDM_g_num_t ***) malloc( n_domain * sizeof(PDM_g_num_t **));
  for (int i_dom = 0; i_dom < n_domain; i_dom++) {
    pn_n_part    [i_dom] = n_part;
    pn_vtx       [i_dom] = (int          *) malloc( n_part * sizeof(int          ));
    pvtx_ln_to_gn[i_dom] = (PDM_g_num_t **) malloc( n_part * sizeof(PDM_g_num_t *));
    for (int i_part = 0; i_part < pn_n_part[i_dom]; i_part++) {
      pn_vtx[i_dom][i_part] = PDM_multipart_part_ln_to_gn_get(mpart_id,
                                                              i_dom,
                                                              i_part,
                                                              PDM_MESH_ENTITY_VERTEX,
                                                              &pvtx_ln_to_gn[i_dom][i_part],
                                                              PDM_OWNERSHIP_KEEP);
    }
  }

  /*
   * Deduction en partition du graphe entre domaine
   */
  PDM_part_domain_interface_t* pdi = PDM_domain_interface_to_part_domain_interface(dom_intrf,
                                                                                   pn_n_part,
                                                                                   NULL,
                                                                                   NULL,
                                                                                   pn_vtx,
                                                                                   NULL,
                                                                                   NULL,
                                                                                   pvtx_ln_to_gn);

  /*
   * Extension
   */
  int n_depth = 1;
  PDM_part_extension_t* part_ext = PDM_part_extension_create(n_domain,
                                                             &n_part,
                                                             PDM_EXTEND_FROM_VTX,
                                                             n_depth,
                                                             comm,
                                                             PDM_OWNERSHIP_KEEP);
  PDM_part_extension_part_domain_interface_shared_set(part_ext, pdi);
  int shift_part = 0;
  for (int i_dom = 0; i_dom < n_domain; i_dom++) {
    for (int i_part = 0; i_part < pn_n_part[i_dom]; i_part++){

      int  n_proc, tn_part;
      int  n_cell, n_face, n_vtx, n_bounds, n_joins, n_part_joins;
      int  scell_face, sface_vtx, sface_bound, sface_join;
      int  n_section;
      int *n_elt;

      PDM_multipart_part_dim_get(mpart_id,
                                 i_dom,
                                 i_part,
                                 &n_section,
                                 &n_elt,
                                 &n_cell,
                                 &n_face,
                                 &n_part_joins,
                                 &n_vtx,
                                 &n_proc,
                                 &tn_part,
                                 &scell_face,
                                 &sface_vtx,
                                 &sface_bound,
                                 &n_bounds,
                                 &sface_join,
                                 &n_joins);
      double       *vtx;
      int          *cell_face_idx, *cell_face, *face_cell, *face_vtx_idx, *face_vtx;
      int          *face_bound_idx, *face_bound, *face_join_idx, *face_join;
      int          *face_part_bound_proc_idx, *face_part_bound_part_idx, *face_part_bound;
      PDM_g_num_t  *cell_ln_to_gn, *face_ln_to_gn, *vtx_ln_to_gn, *face_bound_ln_to_gn, *face_join_ln_to_gn;
      int          *cell_tag, *face_tag, *vtx_tag;
      int         **elt_vtx_idx;
      int         **elt_vtx;
      PDM_g_num_t **elt_section_ln_to_gn;

      PDM_multipart_part_val_get(mpart_id, i_dom, i_part, &elt_vtx_idx, &elt_vtx, &elt_section_ln_to_gn,
                                 &cell_tag, &cell_face_idx, &cell_face, &cell_ln_to_gn,
                                 &face_tag, &face_cell, &face_vtx_idx, &face_vtx, &face_ln_to_gn,
                                 &face_part_bound_proc_idx, &face_part_bound_part_idx, &face_part_bound,
                                 &vtx_tag, &vtx, &vtx_ln_to_gn, &face_bound_idx, &face_bound,
                                 &face_bound_ln_to_gn, &face_join_idx, &face_join, &face_join_ln_to_gn);

      int *vtx_part_bound_proc_idx = NULL;
      int *vtx_part_bound_part_idx = NULL;
      int *vtx_part_bound          = NULL;
      PDM_multipart_part_graph_comm_vtx_data_get(mpart_id,
                                                 i_dom,
                                                 i_part,
                                                 &vtx_part_bound_proc_idx,
                                                 &vtx_part_bound_part_idx,
                                                 &vtx_part_bound);


      int *face_edge     = NULL;
      int *face_edge_idx = NULL;
      int n_face2 = PDM_multipart_part_connectivity_get(mpart_id,
                                                        i_dom,
                                                        i_part,
                                                        PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                                        &face_edge,
                                                        &face_edge_idx,
                                                        PDM_OWNERSHIP_KEEP);
      // assert(n_face == n_face2);
      int *edge_vtx     = NULL;
      int *edge_vtx_idx = NULL;
      int n_edge = PDM_multipart_part_connectivity_get(mpart_id,
                                                       i_dom,
                                                       i_part,
                                                       PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                                       &edge_vtx,
                                                       &edge_vtx_idx,
                                                       PDM_OWNERSHIP_KEEP);
      assert(edge_vtx_idx == NULL);
      PDM_g_num_t* edge_ln_to_gn = NULL;
      int n_edge2 = PDM_multipart_part_ln_to_gn_get(mpart_id,
                                                    i_dom,
                                                    i_part,
                                                    PDM_MESH_ENTITY_EDGE,
                                                    &edge_ln_to_gn,
                                                    PDM_OWNERSHIP_KEEP);

      PDM_UNUSED(n_face2);
      PDM_UNUSED(n_edge2);
      // assert(n_edge2 == n_edge);
      PDM_part_extension_set_part(part_ext, i_dom, i_part,
                                  n_cell,
                                  n_face,
                                  n_part_joins,
                                  n_bounds,
                                  n_edge,
                                  n_vtx,
                                  cell_face_idx,
                                  cell_face,
                                  face_cell,
                                  face_edge_idx,
                                  face_edge,
                                  face_vtx_idx,
                                  face_vtx,
                                  edge_vtx,
                                  face_bound_idx,
                                  face_bound,
                                  NULL, // face_join_idx
                                  NULL, // face_join
                                  face_part_bound_proc_idx,
                                  face_part_bound_part_idx,
                                  face_part_bound,
                                  vtx_part_bound_proc_idx,
                                  vtx_part_bound_part_idx,
                                  vtx_part_bound,
                                  cell_ln_to_gn,
                                  face_ln_to_gn,
                                  edge_ln_to_gn, // edge_ln_to_gn
                                  vtx_ln_to_gn,
                                  face_bound_ln_to_gn,
                                  vtx);

      /*
       *  Mini-Bricoloage
       */
      if(1 == 1) {
        printf("Flag 3 \n");
        PDM_part_mesh_nodal_elmts_t* pmne_vol = PDM_dmesh_nodal_to_part_mesh_nodal_elmts(dmn[shift_part+i_part],
                                                                                         PDM_GEOMETRY_KIND_VOLUMIC,
                                                                                         1, // n_part
                                                                                         &n_vtx,
                                                                                         &vtx_ln_to_gn,
                                                                                         &n_cell,
                                                                                         &cell_ln_to_gn,
                                                                                         NULL);
        int id_section = 0;
        PDM_Mesh_nodal_elt_t t_elt_loc = PDM_part_mesh_nodal_elmts_block_type_get(pmne_vol, id_section);
        int         *elmt_vtx                 = NULL;
        int         *parent_num               = NULL;
        PDM_g_num_t *numabs                   = NULL;
        PDM_g_num_t *parent_entitity_ln_to_gn = NULL;
        PDM_part_mesh_nodal_elmts_block_std_get(pmne_vol, id_section, 0, &elmt_vtx, &numabs, &parent_num, &parent_entitity_ln_to_gn);

        int* cell_num = (int *) malloc(n_cell * sizeof(int));
        for(int i = 0; i < n_cell; ++i) {
          cell_num[i] = i;
        }

        char filename[999];
        sprintf(filename, "out_volumic_%i_%i.vtk", i_part, i_rank);
        const char* field_name[] = {"cell_num", 0 };
        const int * field     [] = {cell_num};
        PDM_vtk_write_std_elements(filename,
                                   n_vtx,
                                   vtx,
                                   vtx_ln_to_gn,
                                   t_elt_loc,
                                   n_cell,
                                   elmt_vtx,
                                   cell_ln_to_gn,
                                   1,
                                   field_name,
                   (const int **)  field);

        free(cell_num);

        PDM_part_mesh_nodal_elmts_free(pmne_vol);
        printf("Flag 4 \n");
      }
    }
    shift_part += pn_n_part[i_dom];
  }


  double t1 = PDM_MPI_Wtime();
  PDM_part_extension_compute(part_ext);
  double dt = PDM_MPI_Wtime() - t1;
  printf("PDM_part_extension_compute : %12.5e \n", dt);
  t1 = PDM_MPI_Wtime();
  /*
   * Export current domain with elements
   *   - Interface with dmesh_nodal
   *   - Export vtk
   */
  shift_part = 0;
  for (int i_dom = 0; i_dom < n_domain; i_dom++){
    for (int i_part = 0; i_part < pn_n_part[i_dom]; i_part++) {

      int  n_proc, tn_part;
      int  n_cell, n_face, n_vtx, n_bounds, n_joins, n_part_joins;
      int  scell_face, sface_vtx, sface_bound, sface_join;
      int  n_section;
      int *n_elt;

      PDM_multipart_part_dim_get(mpart_id,
                                 i_dom,
                                 i_part,
                                 &n_section,
                                 &n_elt,
                                 &n_cell,
                                 &n_face,
                                 &n_part_joins,
                                 &n_vtx,
                                 &n_proc,
                                 &tn_part,
                                 &scell_face,
                                 &sface_vtx,
                                 &sface_bound,
                                 &n_bounds,
                                 &sface_join,
                                 &n_joins);
      double       *vtx;
      int          *cell_face_idx, *cell_face, *face_cell, *face_vtx_idx, *face_vtx;
      int          *face_bound_idx, *face_bound, *face_join_idx, *face_join;
      int          *face_part_bound_proc_idx, *face_part_bound_part_idx, *face_part_bound;
      PDM_g_num_t  *cell_ln_to_gn, *face_ln_to_gn, *vtx_ln_to_gn, *face_bound_ln_to_gn, *face_join_ln_to_gn;
      int          *cell_tag, *face_tag, *vtx_tag;
      int         **elt_vtx_idx;
      int         **elt_vtx;
      PDM_g_num_t **elt_section_ln_to_gn;

      PDM_multipart_part_val_get(mpart_id, i_dom, i_part, &elt_vtx_idx, &elt_vtx, &elt_section_ln_to_gn,
                                 &cell_tag, &cell_face_idx, &cell_face, &cell_ln_to_gn,
                                 &face_tag, &face_cell, &face_vtx_idx, &face_vtx, &face_ln_to_gn,
                                 &face_part_bound_proc_idx, &face_part_bound_part_idx, &face_part_bound,
                                 &vtx_tag, &vtx, &vtx_ln_to_gn, &face_bound_idx, &face_bound,
                                 &face_bound_ln_to_gn, &face_join_idx, &face_join, &face_join_ln_to_gn);


      double* vtx_coord_extended;
      PDM_g_num_t* border_vtx_ln_to_gn;
      PDM_g_num_t* border_cell_ln_to_gn;
      int n_vtx_extended = PDM_part_extension_coord_get(part_ext, i_dom, i_part, &vtx_coord_extended);
      int n_vtx_extended2 = PDM_part_extension_ln_to_gn_get(part_ext, i_dom, i_part, PDM_MESH_ENTITY_VERTEX, &border_vtx_ln_to_gn);
      int n_cell_extended = PDM_part_extension_ln_to_gn_get(part_ext, i_dom, i_part, PDM_MESH_ENTITY_CELL, &border_cell_ln_to_gn);
      assert(n_vtx_extended == n_vtx_extended2);
      if(0 == 1) {
        for(int i_vtx = 0; i_vtx < n_vtx_extended; ++i_vtx) {
          printf("[%i] vtx_coord_extended[%i] = %12.5e %12.5e %12.5e "PDM_FMT_G_NUM" \n", i_part, i_vtx, vtx_coord_extended[3*i_vtx], vtx_coord_extended[3*i_vtx+1], vtx_coord_extended[3*i_vtx+2], border_vtx_ln_to_gn[i_vtx]);
        }

        for(int i_cell = 0; i_cell < n_cell_extended; ++i_cell) {
          printf("[%i] border_cell_ln_to_gn[%i] = "PDM_FMT_G_NUM" \n", i_part, i_cell, border_cell_ln_to_gn[i_cell]);
        }
      }

      if(0 == 1) {
        PDM_log_trace_array_long(border_vtx_ln_to_gn , n_vtx_extended , "border_vtx_ln_to_gn :: ");
        PDM_log_trace_array_long(border_cell_ln_to_gn, n_cell_extended, "border_cell_ln_to_gn :: ");
      }


      // PDM_g_num_t* border_face_ln_to_gn;
      // PDM_g_num_t* border_edge_ln_to_gn;
      // int n_face_extended  = PDM_part_extension_ln_to_gn_get(part_ext, i_dom, i_part, PDM_MESH_ENTITY_FACE, &border_face_ln_to_gn);
      // int n_edge_extended  = PDM_part_extension_ln_to_gn_get(part_ext, i_dom, i_part, PDM_MESH_ENTITY_EDGE, &border_edge_ln_to_gn);

      // int *pface_edge_idx;
      // int *pface_edge;
      // int n_face_extended2 = PDM_part_extension_connectivity_get(part_ext, i_dom, i_part,
      //                                                            PDM_CONNECTIVITY_TYPE_FACE_EDGE,
      //                                                            &pface_edge,
      //                                                            &pface_edge_idx);

      // int *pedge_vtx_idx;
      // int *pedge_vtx;
      // int n_edge_extended2 = PDM_part_extension_connectivity_get(part_ext, i_dom, i_part,
      //                                                            PDM_CONNECTIVITY_TYPE_EDGE_VTX,
      //                                                            &pedge_vtx,
      //                                                            &pedge_vtx_idx);

      // PDM_g_num_t* edge_ln_to_gn = NULL;
      // int n_edge = PDM_multipart_part_ln_to_gn_get(mpart_id,
      //                                               i_dom,
      //                                               i_part,
      //                                               PDM_MESH_ENTITY_EDGE,
      //                                               &edge_ln_to_gn,
      //                                               PDM_OWNERSHIP_KEEP);

      // int n_edge_tot  = n_edge + n_edge_extended;
      // PDM_g_num_t *concat_edge_ln_to_gn  = malloc(    n_edge_tot  * sizeof(PDM_g_num_t));
      // for(int i_edge = 0; i_edge < n_edge; ++i_edge) {
      //   concat_edge_ln_to_gn[i_edge] = edge_ln_to_gn[i_edge];
      // }
      // for(int i_edge = 0; i_edge < n_edge_extended; ++i_edge) {
      //   concat_edge_ln_to_gn[n_edge+i_edge] = border_edge_ln_to_gn[i_edge];
      // }

      // PDM_log_trace_part_connectivity_gnum(pface_edge_idx,
      //                                      pface_edge,
      //                                      border_face_ln_to_gn,
      //                                      concat_edge_ln_to_gn,
      //                                      n_face_extended,
      //                                      "face_edge_extented");
      // free(concat_edge_ln_to_gn);

      int n_vtx_tot  = n_vtx + n_vtx_extended;
      int n_cell_tot = n_cell + n_cell_extended;

      PDM_g_num_t *concat_vtx_ln_to_gn  = malloc(    n_vtx_tot  * sizeof(PDM_g_num_t));
      PDM_g_num_t *concat_cell_ln_to_gn = malloc(    n_cell_tot * sizeof(PDM_g_num_t));
      int         *is_extend            = malloc(    n_cell_tot * sizeof(int        ));
      double      *concat_vtx_coord     = malloc(3 * n_vtx_tot  * sizeof(double     ));
      for(int i_vtx = 0; i_vtx < 3 * n_vtx_tot; ++i_vtx) {
        concat_vtx_coord[i_vtx] = -10000.;
      }

      for(int i_vtx = 0; i_vtx < n_vtx; ++i_vtx) {
        concat_vtx_ln_to_gn[i_vtx] = vtx_ln_to_gn[i_vtx];
        concat_vtx_coord[3*i_vtx  ] = vtx[3*i_vtx  ];
        concat_vtx_coord[3*i_vtx+1] = vtx[3*i_vtx+1];
        concat_vtx_coord[3*i_vtx+2] = vtx[3*i_vtx+2];
      }

      for(int i_cell = 0; i_cell < n_cell; ++i_cell) {
        concat_cell_ln_to_gn[i_cell] = cell_ln_to_gn[i_cell];
        is_extend           [i_cell] = 0;
      }

      for(int i_vtx = 0; i_vtx < n_vtx_extended; ++i_vtx) {
        concat_vtx_ln_to_gn[n_vtx+i_vtx] = border_vtx_ln_to_gn[i_vtx];
        concat_vtx_coord[3*(n_vtx+i_vtx)  ] = vtx_coord_extended[3*i_vtx  ];
        concat_vtx_coord[3*(n_vtx+i_vtx)+1] = vtx_coord_extended[3*i_vtx+1];
        concat_vtx_coord[3*(n_vtx+i_vtx)+2] = vtx_coord_extended[3*i_vtx+2];
      }
      for(int i_cell = 0; i_cell < n_cell_extended; ++i_cell) {
        concat_cell_ln_to_gn[n_cell+i_cell] = border_cell_ln_to_gn[i_cell];
        is_extend           [n_cell+i_cell] = 1;
      }

      // PDM_log_trace_part_connectivity_gnum(pedge_vtx_idx,
      //                                      pedge_vtx,
      //                                      border_edge_ln_to_gn,
      //                                      concat_vtx_ln_to_gn,
      //                                      n_edge_extended,
      //                                      "edge_vtx_extented");

      // PDM_log_trace_array_long(concat_vtx_ln_to_gn , n_vtx_tot , "concat_vtx_ln_to_gn :: ");
      // PDM_log_trace_array_long(concat_cell_ln_to_gn, n_cell_tot, "concat_cell_ln_to_gn :: ");

      printf("Flag 6 \n");
      PDM_part_mesh_nodal_elmts_t* pmne_vol = PDM_dmesh_nodal_to_part_mesh_nodal_elmts(dmn[shift_part+i_part],
                                                                                       PDM_GEOMETRY_KIND_VOLUMIC,
                                                                                       1, // n_part
                                                                                       &n_vtx_tot,
                                                                                       &concat_vtx_ln_to_gn,
                                                                                       &n_cell_tot,
                                                                                       &concat_cell_ln_to_gn,
                                                                                       NULL);

      int id_section = 0;
      PDM_Mesh_nodal_elt_t t_elt_loc = PDM_part_mesh_nodal_elmts_block_type_get(pmne_vol, id_section);
      int         *elmt_vtx                 = NULL;
      int         *parent_num               = NULL;
      PDM_g_num_t *numabs                   = NULL;
      PDM_g_num_t *parent_entitity_ln_to_gn = NULL;
      PDM_part_mesh_nodal_elmts_block_std_get(pmne_vol, id_section, 0, &elmt_vtx, &numabs, &parent_num, &parent_entitity_ln_to_gn);

      char filename[999];
      sprintf(filename, "out_volumic_extended_%i_%i.vtk", i_part, i_rank);

      const char* field_name[] = {"is_extend", 0 };
      const int * field     [] = {is_extend};
      PDM_vtk_write_std_elements(filename,
                                 n_vtx_tot,
                                 concat_vtx_coord,
                                 concat_vtx_ln_to_gn,
                                 t_elt_loc,
                                 n_cell_tot,
                                 elmt_vtx,
                                 concat_cell_ln_to_gn,
                                 1,
                                 field_name,
                (const int **)   field);

      PDM_part_mesh_nodal_elmts_free(pmne_vol);
      free(concat_vtx_ln_to_gn );
      free(concat_cell_ln_to_gn);
      free(concat_vtx_coord    );
      free(is_extend    );

    }
    shift_part += pn_n_part[i_dom];
  }
  dt = PDM_MPI_Wtime() - t1;
  printf("post-treatment : %12.5e \n", dt);

  PDM_part_extension_free(part_ext);


  /*
   *  Pour le debug : extration des faces avec le extract part + vtk
   */

  for (int i_dom = 0; i_dom < n_domain; i_dom++) {
    free(pn_vtx       [i_dom]);
    free(pvtx_ln_to_gn[i_dom]);
  }
  free(pn_vtx       );
  free(pvtx_ln_to_gn);
  free(pn_n_part);

  /*
   *  Free memory
   */
  // free(i_period);

  for (int i = 0; i < n_domain; i++) {
    PDM_dcube_nodal_gen_free(dcube[i]);
    // PDM_dcube_nodal_gen_free(dmn[i]);
  }
  PDM_multipart_free(mpart_id);

  PDM_part_domain_interface_free(pdi);
  PDM_UNUSED(pdi);

  // free(dm);
  // free(dn_vtx);
  // free(dn_face);
  // free(dface_vtx_idx);
  // free(dface_vtx);

  // PDM_dmesh_nodal_to_dmesh_free(dmn_to_dm);
  PDM_domain_interface_free(dom_intrf);
  free(n_part_by_domain);

  // for (int i = 0; i < n_interface; i++) {
  //   free(interface_ids[i]);
  //   free(interface_dom[i]);
  // }
  // free(interface_dn);
  // free(interface_ids);
  // free(interface_dom);
  free(dcube);
  free(dmn);

  if (i_rank == 0) {
    printf("-- End\n");
    fflush(stdout);
  }

  PDM_MPI_Finalize();

  return 0;
}
