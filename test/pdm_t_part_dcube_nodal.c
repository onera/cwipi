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
#include "pdm_dcube_nodal_gen.h"
#include "pdm_dmesh_nodal_to_dmesh.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_geom_elem.h"
#include "pdm_priv.h"
#include "pdm_predicate.h"

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



static double
cross(double a, double b, double c, double d)
{
  double w = d * c;
  double e = fma(-d, c,  w);
  double f = fma( a, b, -w);
  return f + e;
}


static long double
cross_ld(long double a, long double b, long double c, long double d)
{
  long double w = d * c;
  long double e = fmal(-d, c,  w);
  long double f = fmal( a, b, -w);
  return f + e;
}


static double
two_diff(double a, double b)
{
  double x      = (a - b);
  double bvirt  = (a - x);
  double avirt  = x + bvirt;
  double bround = bvirt - b;
  double around = a     - avirt;
  return around+bround;
}

static void
_robust_surface_vector2
(
 int    *face_vtx,
 double *vtx_coord,
 double *surf_vector
)
{
  int v1 = face_vtx[0] - 1;
  int v2 = face_vtx[1] - 1;
  int v3 = face_vtx[2] - 1;

  double u1[3], u2[3]; // b - a | c - a

  u1[0] = vtx_coord[3*v2    ] - vtx_coord[3*v1    ];
  u1[1] = vtx_coord[3*v2 + 1] - vtx_coord[3*v1 + 1];
  u1[2] = vtx_coord[3*v2 + 2] - vtx_coord[3*v1 + 2];
  u2[0] = vtx_coord[3*v3    ] - vtx_coord[3*v1    ];
  u2[1] = vtx_coord[3*v3 + 1] - vtx_coord[3*v1 + 1];
  u2[2] = vtx_coord[3*v3 + 2] - vtx_coord[3*v1 + 2];

  surf_vector[0] = cross(u1[1], u2[2], u1[2], u2[1] );
  surf_vector[1] = cross(u1[2], u2[0], u1[0], u2[2] );
  surf_vector[2] = cross(u1[0], u2[1], u1[1], u2[0] );

  // // double err1 = two_diff(vtx_coord[3*v2    ], vtx_coord[3*v1    ]);
  // // double err2 = two_diff(vtx_coord[3*v2 + 1], vtx_coord[3*v1 + 1]);
  // // double err3 = two_diff(vtx_coord[3*v2 + 2], vtx_coord[3*v1 + 2]);
  // // double err4 = two_diff(vtx_coord[3*v3    ], vtx_coord[3*v1    ]);
  // // double err5 = two_diff(vtx_coord[3*v3 + 1], vtx_coord[3*v1 + 1]);
  // // double err6 = two_diff(vtx_coord[3*v3 + 2], vtx_coord[3*v1 + 2]);

  // // printf("err1 = %20.16e | err2 = %20.16e | err3 = %20.16e | err4 = %20.16e | err5 = %20.16e | err6 = %20.16e \n", err1,
  // //        err2,
  // //        err3,
  // //        err4,
  // //        err5,
  // //        err6);



  // long double u1[3], u2[3]; // b - a | c - a

  // u1[0] = (vtx_coord[3*v2    ] - vtx_coord[3*v1    ]);
  // u1[1] = (vtx_coord[3*v2 + 1] - vtx_coord[3*v1 + 1]);
  // u1[2] = (vtx_coord[3*v2 + 2] - vtx_coord[3*v1 + 2]);
  // u2[0] = (vtx_coord[3*v3    ] - vtx_coord[3*v1    ]);
  // u2[1] = (vtx_coord[3*v3 + 1] - vtx_coord[3*v1 + 1]);
  // u2[2] = (vtx_coord[3*v3 + 2] - vtx_coord[3*v1 + 2]);
  // surf_vector[0] = cross_ld(u1[1], u2[2], u1[2], u2[1] );
  // surf_vector[1] = cross_ld(u1[2], u2[0], u1[0], u2[2] );
  // surf_vector[2] = cross_ld(u1[0], u2[1], u1[1], u2[0] );



}

static void
_robust_surface_vector2_ld
(
 int         *face_vtx,
 long double *vtx_coord,
 long double *surf_vector
)
{
  int v1 = face_vtx[0] - 1;
  int v2 = face_vtx[1] - 1;
  int v3 = face_vtx[2] - 1;

  long double u1[3], u2[3]; // b - a | c - a

  u1[0] = vtx_coord[3*v2    ] - vtx_coord[3*v1    ];
  u1[1] = vtx_coord[3*v2 + 1] - vtx_coord[3*v1 + 1];
  u1[2] = vtx_coord[3*v2 + 2] - vtx_coord[3*v1 + 2];
  u2[0] = vtx_coord[3*v3    ] - vtx_coord[3*v1    ];
  u2[1] = vtx_coord[3*v3 + 1] - vtx_coord[3*v1 + 1];
  u2[2] = vtx_coord[3*v3 + 2] - vtx_coord[3*v1 + 2];

  surf_vector[0] = cross_ld(u1[1], u2[2], u1[2], u2[1] );
  surf_vector[1] = cross_ld(u1[2], u2[0], u1[0], u2[2] );
  surf_vector[2] = cross_ld(u1[0], u2[1], u1[1], u2[0] );
}



static void
_robust_surface_vector
(
 int    *face_vtx,
 double *vtx_coord,
 double *surf_vector
 )
{
  double a[2], b[2], c[2];

  int v1 = face_vtx[0] - 1;
  int v2 = face_vtx[1] - 1;
  int v3 = face_vtx[2] - 1;

  a[0] = vtx_coord[3*v1 + 1];
  a[1] = vtx_coord[3*v1 + 2];
  b[0] = vtx_coord[3*v2 + 1];
  b[1] = vtx_coord[3*v2 + 2];
  c[0] = vtx_coord[3*v3 + 1];
  c[1] = vtx_coord[3*v3 + 2];
  surf_vector[0]  = PDM_predicate_orient2d (a, b, c);
  // surf_vector[0] += PDM_predicate_orient2d (b, c, a);
  // surf_vector[0] += PDM_predicate_orient2d (c, a, b);
  // surf_vector[0] /= 3;


  a[0] = vtx_coord[3*v1 + 2];
  a[1] = vtx_coord[3*v1 + 0];
  b[0] = vtx_coord[3*v2 + 2];
  b[1] = vtx_coord[3*v2 + 0];
  c[0] = vtx_coord[3*v3 + 2];
  c[1] = vtx_coord[3*v3 + 0];
  surf_vector[1] = PDM_predicate_orient2d (a, b, c);
  // surf_vector[1] += PDM_predicate_orient2d (b, c, a);
  // surf_vector[1] += PDM_predicate_orient2d (c, a, b);
  // surf_vector[1] /= 3;


  a[0] = vtx_coord[3*v1 + 0];
  a[1] = vtx_coord[3*v1 + 1];
  b[0] = vtx_coord[3*v2 + 0];
  b[1] = vtx_coord[3*v2 + 1];
  c[0] = vtx_coord[3*v3 + 0];
  c[1] = vtx_coord[3*v3 + 1];
  surf_vector[2] = PDM_predicate_orient2d (a, b, c);
  // surf_vector[2] += PDM_predicate_orient2d (b, c, a);
  // surf_vector[2] += PDM_predicate_orient2d (c, a, b);
  // surf_vector[2] /= 3;
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

  PDM_predicate_exactinit();

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

  // int           dn_cell;
  // int           dn_face;
  // int           dn_vtx;
  // int           n_face_group;
  // PDM_g_num_t  *delmt_cell = NULL;
  // double       *dvtx_coord = NULL;
  // int          *dface_group_idx = NULL;
  // PDM_g_num_t  *dface_group = NULL;
  // int           dface_vtxL;
  // int           dFaceGroupL;

  /*
   *  Create distributed cube
   */

  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;

  PDM_dcube_nodal_t* dcube = PDM_dcube_nodal_gen_create (comm,
                                                           n_vtx_seg,
                                                           n_vtx_seg,
                                                           n_vtx_seg,
                                                           length,
                                                           0.,
                                                           0.,
                                                           0.,
                                                           PDM_MESH_NODAL_TETRA4,
                                                           1,
                                                           PDM_OWNERSHIP_KEEP);
  PDM_dcube_nodal_gen_build (dcube);


  PDM_dmesh_nodal_t* dmn = PDM_dcube_nodal_gen_dmesh_nodal_get(dcube);
  PDM_dmesh_nodal_generate_distribution(dmn);

  /*
   *
   */
  double R[3][3] = {{0.9362934, -0.2896295, 0.1986693},
                    {0.3129918,  0.9447025, -0.0978434},
                    {-0.1593451,  0.1537920,  0.9751703}};

  PDM_g_num_t *vtx_distrib = PDM_dmesh_nodal_vtx_distrib_get(dmn);
  double      *dvtx_coord  = PDM_DMesh_nodal_vtx_get(dmn);
  int dn_vtx = vtx_distrib[i_rank+1] - vtx_distrib[i_rank];

  long double *dvtx_coord_prec = malloc(3 * dn_vtx * sizeof(long double));

  for(int i_vtx = 0; i_vtx < dn_vtx; ++i_vtx) {
    dvtx_coord_prec[3*i_vtx  ] = dvtx_coord[3*i_vtx  ];
    dvtx_coord_prec[3*i_vtx+1] = dvtx_coord[3*i_vtx+1];
    dvtx_coord_prec[3*i_vtx+2] = dvtx_coord[3*i_vtx+2];

    // dvtx_coord[3*i_vtx  ] *= 1.;
    dvtx_coord[3*i_vtx+1] *= 40.;
    dvtx_coord[3*i_vtx+2] *= 18400.;

    double x = dvtx_coord[3*i_vtx  ];
    double y = dvtx_coord[3*i_vtx+1];
    double z = dvtx_coord[3*i_vtx+2];
    //dvtx_coord[3*i_vtx+2] += 0.3 * cos(x * PDM_PI);
    for (int j = 0; j < 3; j++) {
      dvtx_coord[3*i_vtx+j] = R[j][0]*x + R[j][1]*y + R[j][2]*z;
    }

    // dvtx_coord_prec[3*i_vtx  ] *= 1.;
    dvtx_coord_prec[3*i_vtx+1] *= 40.;
    dvtx_coord_prec[3*i_vtx+2] *= 18400.;

    long double xp = dvtx_coord_prec[3*i_vtx  ];
    long double yp = dvtx_coord_prec[3*i_vtx+1];
    long double zp = dvtx_coord_prec[3*i_vtx+2];
    //dvtx_coord_prec[3*i_vtx+2] += 0.3 * cos(x * PDM_PI);
    for (int j = 0; j < 3; j++) {
      dvtx_coord_prec[3*i_vtx+j] = R[j][0]*xp + R[j][1]*yp + R[j][2]*zp;
    }

  }

  PDM_dmesh_nodal_dump_vtk(dmn, PDM_GEOMETRY_KIND_VOLUMIC , "out_volumic");
  PDM_dmesh_nodal_dump_vtk(dmn, PDM_GEOMETRY_KIND_SURFACIC, "out_surfacic");

  PDM_dmesh_nodal_to_dmesh_t* dmntodm = PDM_dmesh_nodal_to_dmesh_create(1, comm, PDM_OWNERSHIP_KEEP);

  PDM_dmesh_nodal_to_dmesh_add_dmesh_nodal(dmntodm, 0, dmn);

  PDM_dmesh_nodal_to_dmesh_set_post_treat_result(dmntodm, 1);

  PDM_dmesh_nodal_to_dmesh_compute(dmntodm,
                                   PDM_DMESH_NODAL_TO_DMESH_TRANSFORM_TO_FACE,
                                   PDM_DMESH_NODAL_TO_DMESH_TRANSLATE_GROUP_TO_FACE);
  // PDM_dmesh_nodal_to_dmesh_compute(dmntodm,
  //                                  PDM_DMESH_NODAL_TO_DMESH_TRANSFORM_TO_EDGE,
  //                                  PDM_DMESH_NODAL_TO_DMESH_TRANSLATE_GROUP_TO_EDGE);
  // PDM_dmesh_nodal_to_dmesh_transform_to_coherent_dmesh(dmntodm, 2);

  PDM_dmesh_t* dmesh = NULL;
  PDM_dmesh_nodal_to_dmesh_get_dmesh(dmntodm, 0, &dmesh);

  int         *dface_vtx_idx;
  PDM_g_num_t *tmp_dface_vtx;
  int dn_face = PDM_dmesh_connectivity_get(dmesh, PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                           &tmp_dface_vtx,
                                           &dface_vtx_idx,
                                           PDM_OWNERSHIP_KEEP);

  int         *dface_cell_idx = NULL;
  PDM_g_num_t *tmp_dface_cell     = NULL;
  dn_face = PDM_dmesh_connectivity_get(dmesh, PDM_CONNECTIVITY_TYPE_FACE_CELL,
                                       &tmp_dface_cell,
                                       &dface_cell_idx,
                                       PDM_OWNERSHIP_KEEP);
  assert(dface_cell_idx == NULL);

  int         *dcell_face_idx;
  PDM_g_num_t *dcell_face;
  int dn_cell = PDM_dmesh_connectivity_get(dmesh, PDM_CONNECTIVITY_TYPE_CELL_FACE,
                                           &dcell_face,
                                           &dcell_face_idx,
                                           PDM_OWNERSHIP_KEEP);

  // Cast to int
  int *pface_vtx  = (int * ) malloc( dface_vtx_idx[dn_face] * sizeof(int));
  int *pface_cell = (int * ) malloc( 2 * dn_face            * sizeof(int));

  for(int i = 0; i < dface_vtx_idx[dn_face]; ++i) {
    pface_vtx[i] = (int) tmp_dface_vtx[i];
  }
  for(int i = 0; i < 2*dn_face; ++i) {
    pface_cell[i] = (int) tmp_dface_cell[i];
  }


  double* check_closed_volume = malloc(3 * dn_cell * sizeof(double));
  double* surface_vector      = malloc(3 * dn_face * sizeof(double));

  long double* check_closed_volume_prec = malloc(3 * dn_cell * sizeof(long double));
  long double* surface_vector_prec      = malloc(3 * dn_face * sizeof(long double));
  for(int i_cell = 0; i_cell < dn_cell; ++i_cell ){
    check_closed_volume[3*i_cell  ] = 0.;
    check_closed_volume[3*i_cell+1] = 0.;
    check_closed_volume[3*i_cell+2] = 0.;

    check_closed_volume_prec[3*i_cell  ] = 0.;
    check_closed_volume_prec[3*i_cell+1] = 0.;
    check_closed_volume_prec[3*i_cell+2] = 0.;
  }

  for(int i_face = 0; i_face < dn_face; ++i_face) {

    int n_vtx_per_face = dface_vtx_idx[i_face+1] - dface_vtx_idx[i_face];
    assert(n_vtx_per_face == 3);

    #if 1
    _robust_surface_vector2 (&pface_vtx[dface_vtx_idx[i_face]],
                            dvtx_coord,
                            &surface_vector[3*i_face]);
    _robust_surface_vector2_ld (&pface_vtx[dface_vtx_idx[i_face]],
                                dvtx_coord_prec,
                                &surface_vector_prec[3*i_face]);
    #else
    PDM_geom_elem_tria_surface_vector(1,
                                      &pface_vtx[dface_vtx_idx[i_face]],
                                      dvtx_coord,
                                      &surface_vector[3*i_face],
                                      NULL,
                                      NULL);
    #endif

    int i_cell_l = pface_cell[2*i_face  ] - 1;
    int i_cell_r = pface_cell[2*i_face+1] - 1;

    check_closed_volume[3*i_cell_l  ] = (check_closed_volume[3*i_cell_l  ] + surface_vector[3*i_face  ]);
    check_closed_volume[3*i_cell_l+1] = (check_closed_volume[3*i_cell_l+1] + surface_vector[3*i_face+1]);
    check_closed_volume[3*i_cell_l+2] = (check_closed_volume[3*i_cell_l+2] + surface_vector[3*i_face+2]);

    if(i_cell_r >= 0) {
      check_closed_volume[3*i_cell_r  ] = (check_closed_volume[3*i_cell_r  ] - surface_vector[3*i_face  ]);
      check_closed_volume[3*i_cell_r+1] = (check_closed_volume[3*i_cell_r+1] - surface_vector[3*i_face+1]);
      check_closed_volume[3*i_cell_r+2] = (check_closed_volume[3*i_cell_r+2] - surface_vector[3*i_face+2]);
    }

    check_closed_volume_prec[3*i_cell_l  ] = (check_closed_volume_prec[3*i_cell_l  ] + surface_vector_prec[3*i_face  ]);
    check_closed_volume_prec[3*i_cell_l+1] = (check_closed_volume_prec[3*i_cell_l+1] + surface_vector_prec[3*i_face+1]);
    check_closed_volume_prec[3*i_cell_l+2] = (check_closed_volume_prec[3*i_cell_l+2] + surface_vector_prec[3*i_face+2]);

    if(i_cell_r >= 0) {
      check_closed_volume_prec[3*i_cell_r  ] = (check_closed_volume_prec[3*i_cell_r  ] - surface_vector_prec[3*i_face  ]);
      check_closed_volume_prec[3*i_cell_r+1] = (check_closed_volume_prec[3*i_cell_r+1] - surface_vector_prec[3*i_face+1]);
      check_closed_volume_prec[3*i_cell_r+2] = (check_closed_volume_prec[3*i_cell_r+2] - surface_vector_prec[3*i_face+2]);
    }

    printf(" surface_vector     [%i] = %20.16e,  %20.16e,  %20.16e\n", i_face, surface_vector[3*i_face  ], surface_vector[3*i_face+1], surface_vector[3*i_face+2] );
    printf(" surface_vector_prec[%i] = %20.16Le,  %20.16Le,  %20.16Le\n", i_face, surface_vector_prec[3*i_face  ], surface_vector_prec[3*i_face+1], surface_vector_prec[3*i_face+2] );

  }
  free (dvtx_coord_prec);

  for(int i_cell = 0; i_cell < dn_cell; ++i_cell ){
    double norm2 = sqrt(  check_closed_volume[3*i_cell  ]*check_closed_volume[3*i_cell  ]
                        + check_closed_volume[3*i_cell+1]*check_closed_volume[3*i_cell+1]
                        + check_closed_volume[3*i_cell+2]*check_closed_volume[3*i_cell+2]);

    printf("norm[%i] =  %20.16e ( %20.16e,  %20.16e,  %20.16e) \n", i_cell,
                                                                norm2, check_closed_volume[3*i_cell  ],
                                                                check_closed_volume[3*i_cell+1],
                                                                check_closed_volume[3*i_cell+2]);

    long double norm2_prec = sqrtl(  check_closed_volume_prec[3*i_cell  ]*check_closed_volume_prec[3*i_cell  ]
                                   + check_closed_volume_prec[3*i_cell+1]*check_closed_volume_prec[3*i_cell+1]
                                   + check_closed_volume_prec[3*i_cell+2]*check_closed_volume_prec[3*i_cell+2]);

    printf("norm2_prec[%i] =  %20.16Le ( %20.16Le,  %20.16Le,  %20.16Le) \n", i_cell,
                                                                norm2_prec, check_closed_volume_prec[3*i_cell  ],
                                                                check_closed_volume_prec[3*i_cell+1],
                                                                check_closed_volume_prec[3*i_cell+2]);

  }
  free (check_closed_volume_prec);


  for(int i_cell = 0; i_cell < dn_cell; ++i_cell ){
    int faces[4];
    faces[0] = dcell_face[dcell_face_idx[i_cell]  ];
    faces[1] = dcell_face[dcell_face_idx[i_cell]+1];
    faces[2] = dcell_face[dcell_face_idx[i_cell]+2];
    faces[3] = dcell_face[dcell_face_idx[i_cell]+3];

    printf(" -------- %i %i %i %i \n", faces[0], faces[1], faces[2], faces[3]);

    long double sx[4];
    long double sy[4];
    long double sz[4];

    sx[0] = PDM_SIGN(faces[0]) * surface_vector[3*(PDM_ABS(faces[0])-1)  ];
    sx[1] = PDM_SIGN(faces[1]) * surface_vector[3*(PDM_ABS(faces[1])-1)  ];
    sx[2] = PDM_SIGN(faces[2]) * surface_vector[3*(PDM_ABS(faces[2])-1)  ];
    sx[3] = PDM_SIGN(faces[3]) * surface_vector[3*(PDM_ABS(faces[3])-1)  ];

    sy[0] = PDM_SIGN(faces[0]) * surface_vector[3*(PDM_ABS(faces[0])-1)+1];
    sy[1] = PDM_SIGN(faces[1]) * surface_vector[3*(PDM_ABS(faces[1])-1)+1];
    sy[2] = PDM_SIGN(faces[2]) * surface_vector[3*(PDM_ABS(faces[2])-1)+1];
    sy[3] = PDM_SIGN(faces[3]) * surface_vector[3*(PDM_ABS(faces[3])-1)+1];

    sz[0] = PDM_SIGN(faces[0]) * surface_vector[3*(PDM_ABS(faces[0])-1)+2];
    sz[1] = PDM_SIGN(faces[1]) * surface_vector[3*(PDM_ABS(faces[1])-1)+2];
    sz[2] = PDM_SIGN(faces[2]) * surface_vector[3*(PDM_ABS(faces[2])-1)+2];
    sz[3] = PDM_SIGN(faces[3]) * surface_vector[3*(PDM_ABS(faces[3])-1)+2];

    long double sx_prec[4];
    long double sy_prec[4];
    long double sz_prec[4];

    sx_prec[0] = PDM_SIGN(faces[0]) * surface_vector_prec[3*(PDM_ABS(faces[0])-1)  ];
    sx_prec[1] = PDM_SIGN(faces[1]) * surface_vector_prec[3*(PDM_ABS(faces[1])-1)  ];
    sx_prec[2] = PDM_SIGN(faces[2]) * surface_vector_prec[3*(PDM_ABS(faces[2])-1)  ];
    sx_prec[3] = PDM_SIGN(faces[3]) * surface_vector_prec[3*(PDM_ABS(faces[3])-1)  ];

    sy_prec[0] = PDM_SIGN(faces[0]) * surface_vector_prec[3*(PDM_ABS(faces[0])-1)+1];
    sy_prec[1] = PDM_SIGN(faces[1]) * surface_vector_prec[3*(PDM_ABS(faces[1])-1)+1];
    sy_prec[2] = PDM_SIGN(faces[2]) * surface_vector_prec[3*(PDM_ABS(faces[2])-1)+1];
    sy_prec[3] = PDM_SIGN(faces[3]) * surface_vector_prec[3*(PDM_ABS(faces[3])-1)+1];

    sz_prec[0] = PDM_SIGN(faces[0]) * surface_vector_prec[3*(PDM_ABS(faces[0])-1)+2];
    sz_prec[1] = PDM_SIGN(faces[1]) * surface_vector_prec[3*(PDM_ABS(faces[1])-1)+2];
    sz_prec[2] = PDM_SIGN(faces[2]) * surface_vector_prec[3*(PDM_ABS(faces[2])-1)+2];
    sz_prec[3] = PDM_SIGN(faces[3]) * surface_vector_prec[3*(PDM_ABS(faces[3])-1)+2];

    double tmp_x = (sx[0] + sx[1] + sx[2] + sx[3]);
    double tmp_y = (sy[0] + sy[1] + sy[2] + sy[3]);
    double tmp_z = (sz[0] + sz[1] + sz[2] + sz[3]);

    long double tmp_x_prec = (sx_prec[0] + sx_prec[1] + sx_prec[2] + sx_prec[3]);
    long double tmp_y_prec = (sy_prec[0] + sy_prec[1] + sy_prec[2] + sy_prec[3]);
    long double tmp_z_prec = (sz_prec[0] + sz_prec[1] + sz_prec[2] + sz_prec[3]);

    long double check_prec = sqrtl( tmp_x_prec*tmp_x_prec + tmp_y_prec*tmp_y_prec + tmp_z_prec*tmp_z_prec );
    double      check      = sqrt( tmp_x*tmp_x + tmp_y*tmp_y + tmp_z*tmp_z );


    printf(" %i      - (%20.16e) - (%20.16e,  %20.16e,  %20.16e) \n", i_cell, check, tmp_x, tmp_y, tmp_z);
    printf(" %i prec - (%20.16Le) - (%20.16Le,  %20.16Le,  %20.16Le) \n", i_cell, check_prec, tmp_x_prec, tmp_y_prec, tmp_z_prec);


  }
  free (surface_vector_prec);
  /*printf("\n\n\n\n\n");


  for(int i_cell = 0; i_cell < dn_cell; ++i_cell ){
    printf("cell[%i] :\n", i_cell);
    int faces[4], sign[4], idx = 0;
    for (int j = dcell_face_idx[i_cell]; j < dcell_face_idx[i_cell+1]; j++) {
      int iface = (int) (dcell_face[j] - 1);
      printf("    (%20.16e,  %20.16e,  %20.16e) \n", i_cell,
                                                                norm2, check_closed_volume[3*i_cell  ],
                                                                check_closed_volume[3*i_cell+1],
                                                                check_closed_volume[3*i_cell+2]);
    }

    for (int j = 0; j < 3; j++) {
      check_closed_volume[3*i_cell+j] =
        surface_vector[3*faces[0]+j] +
        surface_vector[3*faces[1]+j] +
        surface_vector[3*faces[2]+j] +
        surface_vector[3*faces[3]+j];
    }

    double norm2 = sqrt(  check_closed_volume[3*i_cell  ]*check_closed_volume[3*i_cell  ]
                        + check_closed_volume[3*i_cell+1]*check_closed_volume[3*i_cell+1]
                        + check_closed_volume[3*i_cell+2]*check_closed_volume[3*i_cell+2]);

    printf("norm[%i] =  %20.16e ( %20.16e,  %20.16e,  %20.16e) \n", i_cell,
                                                                norm2, check_closed_volume[3*i_cell  ],
                                                                check_closed_volume[3*i_cell+1],
                                                                check_closed_volume[3*i_cell+2]);
                                                                }*/


  // PDM_mpi_win_shared_t* wins_recv = PDM_mpi_win_allocate_shared_create(n_data, sizeof(int), comm);
  // int* recv = (int *) PDM_mpi_win_allocate_shared_get(wins_recv);


  free(pface_cell);
  free(pface_vtx);

  free(check_closed_volume);
  free(surface_vector);

  PDM_dmesh_nodal_to_dmesh_free(dmntodm);
  gettimeofday(&t_elaps_debut, NULL);
  PDM_dcube_nodal_gen_free(dcube);

  PDM_MPI_Finalize();

  return 0;
}
