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
#include "pdm_dcube_nodal_gen.h"
#include "pdm_dmesh_nodal_to_dmesh.h"
#include "pdm_partitioning_algorithm.h"
#include "pdm_dconnectivity_transform.h"
#include "pdm_part_to_block.h"
#include "pdm_block_to_part.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_priv.h"
#include "pdm_vtk.h"
#include "pdm_logging.h"
#include "pdm_gnum.h"
#include "pdm_array.h"
#include "pdm_iso_surface.h"

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
           double        *length)
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
    else
      _usage(EXIT_FAILURE);
    i++;
  }
}

static const double _A = 1.072;
static const double _B = 1.044;
static const double _C = 0.286;
static const double _D = -0.042;
static const double _E = 0.067;
static const double _F = -0.025;


static
inline
double
_unit_circle
(
 double x,
 double y
)
{
  // return x * x + y * y - 0.125;
  // return _A*x*x + _B*x*y + _C*y*y + _D*x + _E*y + _F;
  // double a = (x*x + y*y - 1);
  // return a*a*a - x*x*(6*y*x*x + x*x - 2*y*y*y - 6*y);
  // return a*a*a - x*x*(x*x*(1 + 6*y) - y*(6 + 2*y*y));


  // return x*x*x*x*x*x - y*x*(x*x-y*y) + y*y*y*y*y*y;
  // return PDM_MAX(PDM_ABS(x), PDM_ABS(y)) - 0.25;
  // return PDM_MIN(PDM_ABS(x) + PDM_ABS(y) - 0.25, (x-0.1)*(x-0.1) + (y+0.2)*(y+0.2) - 0.07);

  double v1 = (x-0.23)*(x-0.23) + (y-0.28)*(y-0.28) - 0.03;
  double v2 = (x+0.23)*(x+0.23) + (y-0.28)*(y-0.28) - 0.03;
  double v3 = x*x + y*y - 0.1;
  return PDM_MIN(PDM_MIN(v1, v2), v3);
}


static
inline
void
_unit_circle_gradient
(
 double  x,
 double  y,
 double *df_dx,
 double *df_dy
)
{
  *df_dx = 2*x;
  *df_dy = 2*y;
  *df_dx = 2*_A*x + _B*y + _D;
  *df_dy = 2*_C*y + _B*x + _E;

  *df_dx = 6*x*x*x*x*x + (12*y*y - 24*y - 16)*x*x*x + (6*y*y*y*y + 4*y*y*y - 12*y*y + 12*y + 6)*x;
  // *df_dx = x*( (6 + 2*y*(6 + y*(6 + y*(2 + 3*y)))) + x*x*(16 + 12*y*(-2 + y)) + 6*x*x );
  double a = (x*x + y*y - 1);
  *df_dy = 6*y*a*a - x*x*(-6*y*y + 6*x*x - 6);
  // *df_dy = 6*(y*a*a - x*x*(-y*y + x*x - 1));


  *df_dx = 6*x*x*x*x*x - 3*y*x*x + y*y*y;
  *df_dx = 6*y*y*y*y*y + 3*x*y*y - x*x*x;


  // if (PDM_ABS(x) > PDM_ABS(y)) {
  //   *df_dx = PDM_SIGN(x);
  //   *df_dy = 0.;
  // } else {
  //   *df_dx = 0;
  //   *df_dy = PDM_SIGN(y);
  // }
  // if (PDM_ABS(x) + PDM_ABS(y) - 0.25 < (x-0.1)*(x-0.1) + (y+0.2)*(y+0.2) - 0.07) {
  //   *df_dx = PDM_SIGN(x);
  //   *df_dy = PDM_SIGN(y);
  // } else {
  //   *df_dx = 2*(x-0.1);
  //   *df_dy = 2*(y+0.2);
  // }

  double df_dx1 = 2*(x-0.23);
  double df_dy1 = 2*(y-0.28);
  double df_dx2 = 2*(x+0.23);
  double df_dy2 = 2*(y-0.28);
  double df_dx3 = 2*x;
  double df_dy3 = 2*y;

  double v1 = (x-0.23)*(x-0.23) + (y-0.2)*(y-0.2) - 0.03;
  double v2 = (x+0.23)*(x+0.23) + (y-0.2)*(y-0.2) - 0.03;
  double v3 = x*x + y*y - 0.1;

  if (v1 < v2) {
    if (v1 < v3) {
      *df_dx = df_dx1;
      *df_dy = df_dy1;
    } else {
      *df_dx = df_dx3;
      *df_dy = df_dy3;
    }
  } else {
    if (v2 < v3) {
      *df_dx = df_dx2;
      *df_dy = df_dy2;
    } else {
      *df_dx = df_dx3;
      *df_dy = df_dy3;
    }
  }
}


static void
_dump_vectors
(
 const char   *filename,
 const int     n_pts,
 const double  pts_coord[],
 const double  vector[]
 )
{
  FILE *f = fopen(filename, "w");

  fprintf(f, "# vtk DataFile Version 2.0\n");
  fprintf(f, "mesh\n");
  fprintf(f, "ASCII\n");
  fprintf(f, "DATASET UNSTRUCTURED_GRID\n");

  fprintf(f, "POINTS %d double\n", n_pts);
  for (int i = 0; i < n_pts; i++) {
    for (int j = 0; j < 3; j++) {
      fprintf(f, "%.20lf ", pts_coord[3*i+j]);
    }
    fprintf(f, "\n");
  }

  fprintf(f, "CELLS %d %d\n", n_pts, 2*n_pts);
  for (int i = 0; i < n_pts; i++) {
    fprintf(f, "1 %d\n", i);
  }

  fprintf(f, "CELL_TYPES %d\n", n_pts);
  for (int i = 0; i < n_pts; i++) {
    fprintf(f, "1\n");
  }

  fprintf(f, "POINT_DATA %d\n", n_pts);
  fprintf(f, "VECTORS vector double\n");
  for (int i = 0; i < n_pts; i++) {
    for (int j = 0; j < 3; j++) {
      fprintf(f, "%.20lf ", vector[3*i+j]);
    }
    fprintf(f, "\n");
  }

  fclose(f);
}


#if defined(PDM_HAVE_MKL) || defined(PDM_HAVE_LAPACK)
extern void dgesvd_(const char *jobu,
                    const char *jobvt,
                    int        *m,
                    int        *n,
                    double     *a,
                    int        *lda,
                    double     *s,
                    double     *u,
                    int        *ldu,
                    double     *vt,
                    int        *ldvt,
                    double     *work,
                    int        *lwork,
                    int        *info);


static void
_compute_least_square
(
 const int  m,
 const int  n,
 double    *u,
 double    *s,
 double    *v,
 double    *b,
 double    *x
 )
{
  /* Compute y := S^{-1} U^T b */
  double y[n];

  for (int i = 0; i < n; i++) {

    y[i] = 0.;

    if (PDM_ABS(s[i]) > 1.e-16) {
      for (int j = 0; j < m; j++) {
        y[i] += u[j + m*i] * b[j];
      }

      y[i] /= s[i];
    }
  }

  /* Compute x := V^T y */
  for (int i = 0; i < n; i++) {
    x[i] = 0.;

    for (int j = 0; j < m; j++) {
      x[i] += v[j + m*i] * y[j];
    }
  }
}
#endif


static
void
_iso_line
(
 PDM_MPI_Comm   comm,
 int            n_face,
 int            n_edge,
 int            n_vtx,
 int           *pface_edge_idx,
 int           *pface_edge,
 int           *pedge_vtx_idx,
 int           *pedge_vtx,
 PDM_g_num_t   *pface_ln_to_gn,
 PDM_g_num_t   *pedge_ln_to_gn,
 PDM_g_num_t   *pvtx_ln_to_gn,
 double        *pvtx_coord,
 double       **isoline_dvtx_coord,
 int           *isoline_dn_edge,
 int          **isoline_dedge_vtx_idx,
 int          **isoline_dedge_vtx
)
{
  int i_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);

  // PDM_UNUSED(n_vtx);
  // PDM_UNUSED(pedge_vtx_idx);
  // PDM_UNUSED(pedge_ln_to_gn);
  // PDM_UNUSED(pvtx_ln_to_gn);

  char filename[999];
  sprintf(filename, "out_equi_vtx_coord_%2.2d.vtk", i_rank);
  PDM_vtk_write_point_cloud(filename,
                            n_vtx,
                            pvtx_coord,
                            pvtx_ln_to_gn,
                            NULL);

  /*
   *  Tag edges that cross the iso-line,
   *  compute the intersection point
   *  and the normal at that point
   */
  int    *tag_edge    = PDM_array_zeros_int(n_edge);
  double *edge_coord  = (double *) malloc(sizeof(double) * n_edge * 3);
  double *edge_normal = (double *) malloc(sizeof(double) * n_edge * 3);

  for (int i = 0; i < n_edge; i++) {

    int i_vtx1 = pedge_vtx[2*i  ]-1;
    int i_vtx2 = pedge_vtx[2*i+1]-1;

    double x1 = pvtx_coord[3*i_vtx1  ];
    double y1 = pvtx_coord[3*i_vtx1+1];

    double x2 = pvtx_coord[3*i_vtx2  ];
    double y2 = pvtx_coord[3*i_vtx2+1];

    double val1 = _unit_circle(x1, y1);
    double val2 = _unit_circle(x2, y2);

    int sgn1 = PDM_SIGN(val1);
    int sgn2 = PDM_SIGN(val2);

    if (sgn1 != sgn2) {
      tag_edge[i] = 1;

      double grad1[2], grad2[2];
      _unit_circle_gradient(x1, y1, &grad1[0], &grad1[1]);
      _unit_circle_gradient(x2, y2, &grad2[0], &grad2[1]);

      double t = val1 / (val1 - val2);

      edge_coord[3*i  ] = (1. - t)*x1 + t*x2;
      edge_coord[3*i+1] = (1. - t)*y1 + t*y2;
      edge_coord[3*i+2] = 0.;

      double gx = (1. - t)*grad1[0] + t*grad2[0];
      double gy = (1. - t)*grad1[1] + t*grad2[1];
      double imag = 1. / sqrt(gx*gx + gy*gy);

      edge_normal[3*i  ] = gx * imag;
      edge_normal[3*i+1] = gy * imag;
      edge_normal[3*i+2] = 0.;
    }

    else {
      edge_coord[3*i  ] = 0.;
      edge_coord[3*i+1] = 0.;
      edge_coord[3*i+2] = 0.;

      edge_normal[3*i  ] = 0.;
      edge_normal[3*i+1] = 0.;
      edge_normal[3*i+2] = 0.;
    }

  }


  sprintf(filename, "edge_intersection_%2.2d.vtk", i_rank);
  _dump_vectors (filename,
                 n_edge,
                 edge_coord,
                 edge_normal);


  int n_face_edge_max = 0;
  for (int i = 0; i < n_face; i++) {
    n_face_edge_max = PDM_MAX(n_face_edge_max,
                              pface_edge_idx[i+1] - pface_edge_idx[i]);
  }

  int *pface_edge_inter_idx = (int *) malloc(sizeof(int) * (n_face + 1));
  pface_edge_inter_idx[0] = 0;
  PDM_g_num_t *pface_edge_inter_g_num = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * pface_edge_idx[n_face]);

  double *mat = (double *) malloc (sizeof(double) * n_face_edge_max * 2);
  double *rhs = (double *) malloc (sizeof(double) * n_face_edge_max);

  double *S = (double *) malloc(sizeof(double) * 2);
  double *U = (double *) malloc(sizeof(double) * n_face_edge_max * 2);
  double *V = (double *) malloc(sizeof(double) * n_face_edge_max * 2);

  *isoline_dvtx_coord = (double *) malloc(sizeof(double) * n_face * 3);
  for (int i = 0; i < n_face; i++) {

    pface_edge_inter_idx[i+1] = pface_edge_inter_idx[i];

    int n_tagged_edge = 0;
    for (int j = pface_edge_idx[i]; j < pface_edge_idx[i+1]; j++) {
      int iedge = PDM_ABS(pface_edge[j]) - 1;
      if (tag_edge[iedge]) {
        pface_edge_inter_g_num[pface_edge_inter_idx[i+1]++] = pedge_ln_to_gn[iedge];
        n_tagged_edge++;
      }
    }

    int k = 0;
    for (int j = pface_edge_idx[i]; j < pface_edge_idx[i+1]; j++) {
      int iedge = PDM_ABS(pface_edge[j]) - 1;

      if (tag_edge[iedge]) {

        // mat[2*k  ] = edge_normal[3*iedge  ];
        // mat[2*k+1] = edge_normal[3*iedge+1];
        mat[k              ] = edge_normal[3*iedge  ];
        mat[k+n_tagged_edge] = edge_normal[3*iedge+1];

        rhs[k] = PDM_DOT_PRODUCT(edge_normal + 3*iedge, edge_coord + 3*iedge);

        k++;
      }
    }

    assert(k >= 2);

    // log_trace("A =\n");
    // for (int j = 0; j < n_tagged_edge; j++) {
    //   log_trace("[%20.16f, %20.16f]\n", mat[j], mat[j+n_tagged_edge]);
    // }

    // log_trace("b = [%20.16f, %20.16f]\n", rhs[0], rhs[1]);

    /* Compute SVD of mat */
#if defined(PDM_HAVE_MKL) || defined(PDM_HAVE_LAPACK)
    int info = 0;
    int n_row = n_tagged_edge;
    int n_col = 2;
    int lwork = 15;//>= MAX(1,5*MIN(M,N))
    double work[15];
    dgesvd_("S",
            "S",
            &n_row,
            &n_col,
            mat,
            &n_row,
            S,
            U,
            &n_row,
            V,
            &n_row,
            work,
            &lwork,
            &info);
#else
    printf("Error : LAPACK or MKL are mandatory, recompile with them. \n");
    abort();
#endif

    // log_trace("S = [%20.16f, %20.16f]\n", S[0], S[1]);
    // log_trace("U =\n");
    // for (int j = 0; j < n_tagged_edge; j++) {
    //   log_trace("[%20.16f, %20.16f]\n", U[j], U[j+n_tagged_edge]);
    // }
    // log_trace("V  =\n");
    // for (int j = 0; j < n_tagged_edge; j++) {
    //   log_trace("[%20.16f, %20.16f]\n", V[j], V[j+n_tagged_edge]);
    // }

    /* Solve for iso-line vertex coordinates */
    double *sol = *isoline_dvtx_coord + 3*i;

    _compute_least_square (n_tagged_edge,
                           2,
                           U,
                           S,
                           V,
                           rhs,
                           sol);
    sol[2] = 0.;
    // log_trace("sol = [%20.16f, %20.16f]\n", sol[0], sol[1]);
  }
  free(mat);
  free(rhs);
  free(S);
  free(U);
  free(V);
  free(tag_edge);
  free(edge_coord);
  free(edge_normal);

  pface_edge_inter_g_num = realloc(pface_edge_inter_g_num, sizeof(PDM_g_num_t) * pface_edge_inter_idx[n_face]);

  PDM_log_trace_connectivity_long(pface_edge_inter_idx, pface_edge_inter_g_num, n_face, "pface_edge_inter : ");


  sprintf(filename, "isoline_dvtx_coord_%2.2d.vtk", i_rank);
  PDM_vtk_write_point_cloud(filename,
                            n_face,
                            *isoline_dvtx_coord,
                            pface_ln_to_gn,
                            NULL);





  int tn_face_edge_inter = pface_edge_inter_idx[n_face];

  PDM_part_to_block_t *ptb = PDM_part_to_block_create (PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                      PDM_PART_TO_BLOCK_POST_MERGE,
                                                      1.,
                                                      &pface_edge_inter_g_num,
                                                      NULL,
                                                      &tn_face_edge_inter,
                                                      1,
                                                      comm);

  PDM_g_num_t *pface_edge_inter_face_g_num = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * tn_face_edge_inter);
  for (int i = 0; i < n_face; i++) {
    for (int j = pface_edge_inter_idx[i]; j < pface_edge_inter_idx[i+1]; j++) {
      pface_edge_inter_face_g_num[j] = pface_ln_to_gn[i];
    }
  }

  double *pface_edge_inter_face_coord = (double *) malloc(sizeof(double) * tn_face_edge_inter * 3);
  for (int i = 0; i < n_face; i++) {
    for (int j = pface_edge_inter_idx[i]; j < pface_edge_inter_idx[i+1]; j++) {
      for (int k = 0; k < 3; k++) {
        pface_edge_inter_face_coord[3*j+k] = (*isoline_dvtx_coord)[3*i+k];
      }
    }
  }

  int *part_stride = PDM_array_const_int(tn_face_edge_inter, 1);

  *isoline_dn_edge = PDM_part_to_block_n_elt_block_get(ptb);

  int *dedge_isoline_vtx_n = NULL;
  PDM_part_to_block_exch (ptb,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_VAR_INTERLACED,
                          1,
                          &part_stride,
                (void **) &pface_edge_inter_face_g_num,
                          &dedge_isoline_vtx_n,
                (void **) isoline_dedge_vtx);

  PDM_log_trace_array_int(dedge_isoline_vtx_n, *isoline_dn_edge, "dedge_isoline_vtx_n : ");

  if (1) {
    int    *dedge_isoline_vtx_coord_n = NULL;
    double *dedge_isoline_vtx_coord   = NULL;
    PDM_part_to_block_exch (ptb,
                            3*sizeof(double),
                            PDM_STRIDE_VAR_INTERLACED,
                            1,
                            &part_stride,
                  (void **) &pface_edge_inter_face_coord,
                            &dedge_isoline_vtx_coord_n,
                  (void **) &dedge_isoline_vtx_coord);

    PDM_log_trace_array_int(dedge_isoline_vtx_coord_n, *isoline_dn_edge, "dedge_isoline_vtx_coord_n : ");

    sprintf(filename, "dedge_isoline_%2.2d.vtk", i_rank);
    PDM_vtk_write_lines (filename,
                         *isoline_dn_edge,
                         dedge_isoline_vtx_coord,
                         NULL,
                         NULL);
    free(dedge_isoline_vtx_coord_n);
    free(dedge_isoline_vtx_coord  );
  }

  free(part_stride);
  *isoline_dedge_vtx_idx = PDM_array_new_idx_from_sizes_int(dedge_isoline_vtx_n, *isoline_dn_edge);
  free(dedge_isoline_vtx_n);

  // PDM_log_trace_connectivity_long(*isoline_dedge_vtx_idx, *isoline_dedge_vtx, *isoline_dn_edge, "dedge_isoline_vtx : ");




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

  /*
   *  Read args
   */

  _read_args(argc,
             argv,
             &n_vtx_seg,
             &length);

  /*
   *  Init
   */

  int i_rank;
  int n_rank;

  PDM_MPI_Init(&argc, &argv);
  PDM_MPI_Comm_rank(PDM_MPI_COMM_WORLD, &i_rank);
  PDM_MPI_Comm_size(PDM_MPI_COMM_WORLD, &n_rank);

  /*
   *  Create distributed cube
   */

  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;

  PDM_dcube_nodal_t* dcube = PDM_dcube_nodal_gen_create (comm,
                                                         n_vtx_seg,
                                                         n_vtx_seg,
                                                         n_vtx_seg,
                                                         length,
                                                         -0.5*length,// -0.35,
                                                         -0.5*length,// -0.3,
                                                         0.,
                                                         PDM_MESH_NODAL_TRIA3,
                                                         // PDM_MESH_NODAL_QUAD4,
                                                         1,
                                                         PDM_OWNERSHIP_KEEP);
  PDM_dcube_nodal_gen_build (dcube);


  PDM_dmesh_nodal_t* dmn = PDM_dcube_nodal_gen_dmesh_nodal_get(dcube);
  PDM_dmesh_nodal_generate_distribution(dmn);

  PDM_g_num_t *vtx_distrib = PDM_dmesh_nodal_vtx_distrib_get(dmn);
  double      *dvtx_coord  = PDM_DMesh_nodal_vtx_get(dmn);
  int dn_vtx = vtx_distrib[i_rank+1] - vtx_distrib[i_rank];

  if(1 == 1) {
    // PDM_dmesh_nodal_dump_vtk(dmn, PDM_GEOMETRY_KIND_VOLUMIC , "out_volumic");
    PDM_dmesh_nodal_dump_vtk(dmn, PDM_GEOMETRY_KIND_SURFACIC, "out_surfacic");
  }

  PDM_dmesh_nodal_to_dmesh_t* dmntodm = PDM_dmesh_nodal_to_dmesh_create(1, comm, PDM_OWNERSHIP_KEEP);

  PDM_dmesh_nodal_to_dmesh_add_dmesh_nodal(dmntodm, 0, dmn);

  PDM_dmesh_nodal_to_dmesh_set_post_treat_result(dmntodm, 1);

  PDM_dmesh_nodal_to_dmesh_compute(dmntodm,
                                   PDM_DMESH_NODAL_TO_DMESH_TRANSFORM_TO_EDGE,
                                   PDM_DMESH_NODAL_TO_DMESH_TRANSLATE_GROUP_TO_EDGE);

  PDM_dmesh_t* dmesh = NULL;
  PDM_dmesh_nodal_to_dmesh_get_dmesh(dmntodm, 0, &dmesh);

  int         *dface_edge_idx;
  PDM_g_num_t *dface_edge;
  int dn_face = PDM_dmesh_connectivity_get(dmesh, PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                           &dface_edge,
                                           &dface_edge_idx,
                                           PDM_OWNERSHIP_KEEP);
  int         *dedge_vtx_idx;
  PDM_g_num_t *dedge_vtx;
  int dn_edge  = PDM_dmesh_connectivity_get(dmesh, PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                           &dedge_vtx,
                                           &dedge_vtx_idx,
                                           PDM_OWNERSHIP_KEEP);

  if(0 == 1) {
    PDM_log_trace_connectivity_long(dface_edge_idx, dface_edge, dn_face, "dface_edge ::");
    PDM_log_trace_connectivity_long(dedge_vtx_idx , dedge_vtx , dn_edge, "dedge_vtx  ::");
  }

  PDM_g_num_t* distrib_edge = NULL;
  PDM_dmesh_distrib_get(dmesh, PDM_MESH_ENTITY_EDGE  , &distrib_edge);
  assert(distrib_edge != NULL);

  PDM_g_num_t* distrib_face = NULL;
  PDM_dmesh_distrib_get(dmesh, PDM_MESH_ENTITY_FACE  , &distrib_face);
  assert(distrib_face != NULL);

  PDM_UNUSED(dn_face);
  PDM_UNUSED(dn_edge);
  PDM_UNUSED(dn_vtx);
  PDM_UNUSED(dvtx_coord);

  // Compute dfield and gradient field
  double *dfield          = (double *) malloc(     dn_vtx * sizeof(double));
  double *dgradient_field = (double *) malloc( 3 * dn_vtx * sizeof(double));

  for(int i = 0; i < dn_vtx; ++i) {

    double x1 = dvtx_coord[3*i  ];
    double y1 = dvtx_coord[3*i+1];
    dfield[i] = _unit_circle(x1, y1);

    _unit_circle_gradient(x1, y1, &dgradient_field[3*i], &dgradient_field[3*i+1]);

    dgradient_field[3*i+2] = 0;
  }

  PDM_iso_surface_t* isos = PDM_iso_surface_create(2, PDM_ISO_SURFACE_KIND_FIELD, 1, PDM_OWNERSHIP_KEEP, comm);
  // PDM_iso_surface_t* isos = PDM_iso_surface_create(2, PDM_ISO_SURFACE_KIND_PLANE, 1, PDM_OWNERSHIP_KEEP, comm);

  PDM_iso_surface_plane_equation_set(isos, 1., 0., 0., -0);
  // PDM_iso_surface_plane_equation_set(isos, 1., 0.5, 0.25, -0.0234);

  PDM_iso_surface_dconnectivity_set(isos,
                                    PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                    dface_edge,
                                    dface_edge_idx);
  PDM_iso_surface_dconnectivity_set(isos,
                                    PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                    dedge_vtx,
                                    dedge_vtx_idx);

  PDM_iso_surface_distrib_set(isos, PDM_MESH_ENTITY_FACE  , distrib_face);
  PDM_iso_surface_distrib_set(isos, PDM_MESH_ENTITY_EDGE  , distrib_edge);
  PDM_iso_surface_distrib_set(isos, PDM_MESH_ENTITY_VERTEX, vtx_distrib);

  PDM_iso_surface_dvtx_coord_set (isos, dvtx_coord     );
  PDM_iso_surface_dfield_set     (isos, dfield         );
  PDM_iso_surface_dgrad_field_set(isos, dgradient_field);

  PDM_iso_surface_compute(isos);

  PDM_iso_surface_free(isos);

  free(dfield);
  free(dgradient_field);

  PDM_dmesh_nodal_to_dmesh_free(dmntodm);
  PDM_dcube_nodal_gen_free(dcube);


  if (i_rank == 0) {
    printf("-- End\n");
    fflush(stdout);
  }

  PDM_MPI_Finalize();

  return 0;

  /*
   * Select gnum that contains iso-surface
   */

  PDM_g_num_t* edge_ln_to_gn = (PDM_g_num_t * ) malloc( dn_edge * sizeof(PDM_g_num_t));
  for(int i = 0; i < dn_edge; ++i) {
    edge_ln_to_gn[i] = distrib_edge[i_rank] + i + 1;
  }

  int          pn_vtx           = 0;
  PDM_g_num_t *pvtx_ln_to_gn  = NULL;
  int         *pedge_vtx_idx    = NULL;
  int         *pedge_vtx        = NULL;
  PDM_part_dconnectivity_to_pconnectivity_sort_single_part(comm,
                                                           distrib_edge,
                                                           dedge_vtx_idx,
                                                           dedge_vtx,
                                                           dn_edge,
                                     (const PDM_g_num_t *) edge_ln_to_gn,
                                                           &pn_vtx,
                                                           &pvtx_ln_to_gn,
                                                           &pedge_vtx_idx,
                                                           &pedge_vtx);

  double **tmp_pvtx_coord = NULL;
  PDM_part_dcoordinates_to_pcoordinates(comm,
                                        1,
                                        vtx_distrib,
                                        dvtx_coord,
                                        &pn_vtx,
                 (const PDM_g_num_t **) &pvtx_ln_to_gn,
                                        &tmp_pvtx_coord);
  double* pvtx_coord = tmp_pvtx_coord[0];
  free(tmp_pvtx_coord);

  /*
   * Select edge
   */
  int *dedge_tag = (int * ) malloc(dn_edge * sizeof(int));

  for(int i = 0; i < dn_edge; ++i) {

    int i_vtx1 = pedge_vtx[2*i  ]-1;
    int i_vtx2 = pedge_vtx[2*i+1]-1;

    dedge_tag[i] = 0;

    double x1 = pvtx_coord[3*i_vtx1  ];
    double y1 = pvtx_coord[3*i_vtx1+1];
    // int x3 = pvtx_coord[3*i_vtx1+2];

    double x2 = pvtx_coord[3*i_vtx2  ];
    double y2 = pvtx_coord[3*i_vtx2+1];

    double val1 = _unit_circle(x1, y1);
    double val2 = _unit_circle(x2, y2);

    int sgn1 = PDM_SIGN(val1);
    int sgn2 = PDM_SIGN(val2);

    if(sgn1 * sgn2 < 0) {
      dedge_tag[i] = 1;
    }
  }

  int* val = malloc(pn_vtx * sizeof(int));
  printf("pn_vtx = %i \n", pn_vtx);
  for(int i = 0; i < pn_vtx; ++i) {
    double x1 = pvtx_coord[3*i  ];
    double y1 = pvtx_coord[3*i+1];
    val[i]    = PDM_SIGN( _unit_circle(x1, y1) );
  }

  char filename[999];
  sprintf(filename, "out_edge_tag_%2.2d.vtk", i_rank);
  PDM_vtk_write_point_cloud(filename,
                            pn_vtx,
                            pvtx_coord,
                            NULL,
                            val);
  free(val);


  // PDM_log_trace_array_int(dedge_tag, dn_edge, "dedge_tag");
  free(pvtx_coord);

  /*
   * block_to_part on dface_edge
   */
  PDM_block_to_part_t* btp = PDM_block_to_part_create(distrib_edge,
                               (const PDM_g_num_t **) &dface_edge,
                                                      &dface_edge_idx[dn_face],
                                                      1,
                                                      comm);

  int strid_one = 1;
  int **tmp_dface_edge_tag = NULL;
  PDM_block_to_part_exch(btp,
                         sizeof(int),
                         PDM_STRIDE_CST_INTERLACED,
                         &strid_one,
            (void *  )   dedge_tag,
            (int  ***)   NULL,
            (void ***)  &tmp_dface_edge_tag);
  int *dface_edge_tag = tmp_dface_edge_tag[0];
  free(tmp_dface_edge_tag);
  free(dedge_tag);

  int         *dface_tag        = malloc(dn_face * sizeof(int        ));
  PDM_g_num_t *face_to_extract_gnum = malloc(dn_face * sizeof(PDM_g_num_t));
  int  n_face_tag = 0;
  for(int i = 0; i < dn_face; ++i) {
    dface_tag[i] = 0;
    for(int idx_face = dface_edge_idx[i]; idx_face < dface_edge_idx[i+1]; ++idx_face) {
      if(dface_edge_tag[idx_face] == 1) {
        dface_tag[i] = 1;
        face_to_extract_gnum[n_face_tag++] = distrib_face[i_rank] + i + 1;
        break;
      }
    }
  }

  if(0 == 1) {
    PDM_log_trace_array_int (dedge_tag, dn_face, "dface_tag");
    PDM_log_trace_array_long(face_to_extract_gnum, n_face_tag, "face_to_extract_gnum");
  }

  PDM_block_to_part_free(btp);

  free(edge_ln_to_gn);
  free(pvtx_ln_to_gn);
  free(pedge_vtx_idx);
  free(pedge_vtx);
  free(dface_edge_tag);

  /*
   *  Visu
   */
  int         *dface_vtx_idx = NULL;
  PDM_g_num_t *dface_vtx     = NULL;

  for(int i = 0; i < dface_edge_idx[dn_face]; ++i) {
    dface_edge[i] = PDM_ABS(dface_edge[i]);
  }

  PDM_deduce_combine_connectivity(comm,
                                  distrib_face,
                                  distrib_edge,
                                  dface_edge_idx,
                                  dface_edge,
                                  dedge_vtx_idx,
                                  dedge_vtx,
                                  1,
                                  &dface_vtx_idx,
                                  &dface_vtx);

  // PDM_log_trace_connectivity_long(dface_vtx_idx, dface_vtx, dn_face, "dface_edge ::");
  for(int i = 0; i < dface_vtx_idx[dn_face]; ++i) {
    dface_vtx[i] = PDM_ABS(dface_vtx[i]);
  }

  // PDM_log_trace_array_int(dface_vtx_idx, dn_face+1,"dface_vtx_idx :: " );
  // PDM_log_trace_array_int(dface_edge_idx, dn_face+1,"dface_edge_idx :: " );

  int          pn_extract_vtx           = 0;
  PDM_g_num_t *pn_extract_vtx_ln_to_gn  = NULL;
  int         *pface_vtx_idx    = NULL;
  int         *pface_vtx        = NULL;
  PDM_part_dconnectivity_to_pconnectivity_sort_single_part(comm,
                                                           distrib_face,
                                                           dface_vtx_idx,
                                                           dface_vtx,
                                                           n_face_tag,
                                     (const PDM_g_num_t *) face_to_extract_gnum,
                                                           &pn_extract_vtx,
                                                           &pn_extract_vtx_ln_to_gn,
                                                           &pface_vtx_idx,
                                                           &pface_vtx);

  double **tmp_pvtx_extract_coord = NULL;
  PDM_part_dcoordinates_to_pcoordinates(comm,
                                        1,
                                        vtx_distrib,
                                        dvtx_coord,
                                        &pn_extract_vtx,
                 (const PDM_g_num_t **) &pn_extract_vtx_ln_to_gn,
                                        &tmp_pvtx_extract_coord);
  double* pvtx_extract_coord = tmp_pvtx_extract_coord[0];
  free(tmp_pvtx_extract_coord);


  sprintf(filename, "out_mesh_%2.2d.vtk", i_rank);
  PDM_vtk_write_polydata(filename,
                         pn_extract_vtx,
                         pvtx_extract_coord,
                         pn_extract_vtx_ln_to_gn,//NULL,
                         n_face_tag,
                         pface_vtx_idx,
                         pface_vtx,
                         face_to_extract_gnum,
                         NULL);



  free(pn_extract_vtx_ln_to_gn);



  double* face_center = (double *) malloc( 3 * n_face_tag * sizeof(double));
  for(int i_face = 0; i_face < n_face_tag; ++i_face) {
    face_center[3*i_face  ] = 0.;
    face_center[3*i_face+1] = 0.;
    face_center[3*i_face+2] = 0.;
    int nvtx_on_face = pface_vtx_idx[i_face+1] - pface_vtx_idx[i_face];
    for(int idx_vtx = pface_vtx_idx[i_face]; idx_vtx < pface_vtx_idx[i_face+1]; ++idx_vtx) {
      int i_vtx = pface_vtx[idx_vtx]-1;
      face_center[3*i_face  ] += pvtx_extract_coord[3*i_vtx  ];
      face_center[3*i_face+1] += pvtx_extract_coord[3*i_vtx+1];
      face_center[3*i_face+2] += pvtx_extract_coord[3*i_vtx+2];
    }
    double inv = 1./(double) nvtx_on_face;
    face_center[3*i_face  ] = face_center[3*i_face  ] * inv;
    face_center[3*i_face+1] = face_center[3*i_face+1] * inv;
    face_center[3*i_face+2] = face_center[3*i_face+2] * inv;
  }

  sprintf(filename, "out_equi_face_coord_%2.2d.vtk", i_rank);
  PDM_vtk_write_point_cloud(filename,
                            n_face_tag,
                            face_center,
                            NULL,
                            NULL);


  free(pface_vtx_idx);
  free(pface_vtx);


  free(pvtx_extract_coord);

  /*
   * Rebuild partition that contains faces and reequilibrate
   */
  PDM_gen_gnum_t* gnum_equi = PDM_gnum_create(3, 1, PDM_FALSE, 0., comm, PDM_OWNERSHIP_USER);
  PDM_gnum_set_from_coords(gnum_equi, 0, n_face_tag, face_center, NULL);
  PDM_gnum_compute(gnum_equi);
  PDM_g_num_t* child_equi_face_gnum = PDM_gnum_get(gnum_equi, 0);
  PDM_gnum_free(gnum_equi);
  free(face_center);

  /*
   * Equilibrage avec le part_to_block
   */
  PDM_part_to_block_t *ptb = PDM_part_to_block_create (PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                       PDM_PART_TO_BLOCK_POST_CLEANUP,
                                                       1.,
                                                       &child_equi_face_gnum,
                                                       NULL,
                                                       &n_face_tag,
                                                       1,
                                                       comm);

  int n_face_equi = PDM_part_to_block_n_elt_block_get (ptb);
  PDM_g_num_t *block_face_g_num_child_equi = PDM_part_to_block_block_gnum_get (ptb);

  PDM_g_num_t *block_face_equi_parent_g_num = NULL;
  PDM_part_to_block_exch (ptb,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_CST_INTERLACED,
                          1,
                          NULL,
               (void **) &face_to_extract_gnum,
                          NULL,
               (void **) &block_face_equi_parent_g_num);

  if(0 == 1) {
    PDM_log_trace_array_long(block_face_equi_parent_g_num, n_face_equi, "block_face_equi_parent_g_num ::");
    PDM_log_trace_array_long(block_face_g_num_child_equi , n_face_equi, "block_face_g_num_child_equi  ::");
  }

  /*
   * Je prepare tout pour mon petit Bastien
   */
  int          pn_edge_equi        = 0;
  PDM_g_num_t *pequi_parent_edge_ln_to_gn = NULL;
  int         *pequi_face_edge_idx = NULL;
  int         *pequi_face_edge     = NULL;
  PDM_part_dconnectivity_to_pconnectivity_sort_single_part(comm,
                                                           distrib_face,
                                                           dface_edge_idx,
                                                           dface_edge,
                                                           n_face_equi,
                                     (const PDM_g_num_t *) block_face_equi_parent_g_num,
                                                           &pn_edge_equi,
                                                           &pequi_parent_edge_ln_to_gn,
                                                           &pequi_face_edge_idx,
                                                           &pequi_face_edge);

  PDM_gen_gnum_t* gnum_edge = PDM_gnum_create(3, 1, PDM_FALSE, 0., comm, PDM_OWNERSHIP_USER);
  PDM_gnum_set_from_parents(gnum_edge, 0, pn_edge_equi, pequi_parent_edge_ln_to_gn);
  PDM_gnum_compute(gnum_edge);
  PDM_g_num_t* pequi_edge_ln_to_gn = PDM_gnum_get(gnum_edge, 0);
  PDM_gnum_free(gnum_edge);


  int          pn_vtx_equi        = 0;
  PDM_g_num_t *pequi_parent_vtx_ln_to_gn = NULL;
  int         *pequi_edge_vtx_idx = NULL;
  int         *pequi_edge_vtx     = NULL;
  PDM_part_dconnectivity_to_pconnectivity_sort_single_part(comm,
                                                           distrib_edge,
                                                           dedge_vtx_idx,
                                                           dedge_vtx,
                                                           pn_edge_equi,
                                     (const PDM_g_num_t *) pequi_parent_edge_ln_to_gn,
                                                           &pn_vtx_equi,
                                                           &pequi_parent_vtx_ln_to_gn,
                                                           &pequi_edge_vtx_idx,
                                                           &pequi_edge_vtx);

  PDM_gen_gnum_t* gnum_vtx = PDM_gnum_create(3, 1, PDM_FALSE, 0., comm, PDM_OWNERSHIP_USER);
  PDM_gnum_set_from_parents(gnum_vtx, 0, pn_vtx_equi, pequi_parent_vtx_ln_to_gn);
  PDM_gnum_compute(gnum_vtx);
  PDM_g_num_t* pequi_vtx_ln_to_gn = PDM_gnum_get(gnum_vtx, 0);
  PDM_gnum_free(gnum_vtx);

  double **tmp_pequi_vtx_coord = NULL;
  PDM_part_dcoordinates_to_pcoordinates(comm,
                                        1,
                                        vtx_distrib,
                                        dvtx_coord,
                                        &pn_vtx_equi,
                 (const PDM_g_num_t **) &pequi_parent_vtx_ln_to_gn,
                                        &tmp_pequi_vtx_coord);
  double* pequi_vtx_coord = tmp_pequi_vtx_coord[0];
  free(tmp_pequi_vtx_coord);


  double *isoline_dvtx_coord    = NULL;
  int     isoline_dn_edge       = 0;
  int    *isoline_dedge_vtx_idx = NULL;
  int    *isoline_dedge_vtx     = NULL;
  _iso_line(comm,
            n_face_equi,
            pn_edge_equi,
            pn_vtx_equi,
            pequi_face_edge_idx,
            pequi_face_edge,
            pequi_edge_vtx_idx,
            pequi_edge_vtx,
            block_face_g_num_child_equi,
            pequi_edge_ln_to_gn,
            pequi_vtx_ln_to_gn,
            pequi_vtx_coord,
            &isoline_dvtx_coord,
            &isoline_dn_edge,
            &isoline_dedge_vtx_idx,
            &isoline_dedge_vtx);

  free(pequi_edge_ln_to_gn);
  free(pequi_vtx_ln_to_gn);


  free(pequi_vtx_coord);

  free(pequi_parent_vtx_ln_to_gn);
  free(pequi_edge_vtx_idx);
  free(pequi_edge_vtx);

  free(pequi_parent_edge_ln_to_gn);
  free(pequi_face_edge_idx);
  free(pequi_face_edge);


  free(block_face_equi_parent_g_num);

  if (isoline_dvtx_coord != NULL) {
    free(isoline_dvtx_coord);
  }
  if (isoline_dedge_vtx_idx != NULL) {
    free(isoline_dedge_vtx_idx);
  }
  if (isoline_dedge_vtx != NULL) {
    free(isoline_dedge_vtx);
  }


  PDM_part_to_block_free(ptb);
  free(dface_vtx_idx);
  free(dface_vtx);
  free(child_equi_face_gnum);


  free(face_to_extract_gnum);
  free(dface_tag);

  PDM_dmesh_nodal_to_dmesh_free(dmntodm);
  PDM_dcube_nodal_gen_free(dcube);


  if (i_rank == 0) {
    printf("-- End\n");
    fflush(stdout);
  }


  PDM_MPI_Finalize();

  return 0;
}
