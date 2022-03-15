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
#include "pdm_dmesh_nodal_priv.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_geom_elem.h"
#include "pdm_priv.h"
#include "pdm_dconnectivity_transform.h"
#include "pdm_partitioning_algorithm.h"
#include "pdm_vtk.h"
#include "pdm_logging.h"
#include "pdm_distrib.h"

/*============================================================================
 * Type definitions
 *============================================================================*/

#define _MIN(a,b) ((a) < (b) ? (a) : (b))
#define _MAX(a,b) ((a) > (b) ? (a) : (b))


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
           int           *order,
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
    else if (strcmp(argv[i], "-o") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *order = atoi(argv[i]);
      }
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



static void
_dmesh_nodal_dump_vtk
(
       PDM_dmesh_nodal_t   *dmn,
       int                  order,
       PDM_geometry_kind_t  geom_kind,
 const char                *filename_patter
)
{
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(dmn->comm, &i_rank);
  PDM_MPI_Comm_size(dmn->comm, &n_rank);

  int* sections_id = PDM_DMesh_nodal_sections_id_get(dmn, geom_kind);
  int n_section    = PDM_DMesh_nodal_n_section_get(dmn, geom_kind);


  const char *field_name = "group";
  int n_field = 0;
  _pdm_dmesh_nodal_elts_t *dmne = NULL;
  int *delt_group_idx = NULL;
  int *delt_group     = NULL;
  double **field = NULL;
  if (geom_kind == PDM_GEOMETRY_KIND_RIDGE) {
    dmne = dmn->ridge;
  } else if (geom_kind == PDM_GEOMETRY_KIND_SURFACIC && dmn->mesh_dimension == 3) {
    dmne = dmn->surfacic;
  }

  PDM_g_num_t *distrib_elt = NULL;
  if (dmne != NULL) {
    distrib_elt = PDM_compute_uniform_entity_distribution(dmn->comm,
                                                          dmne->n_g_elmts);
    PDM_log_trace_array_long(distrib_elt, n_rank+1, "distrib_elt : ");

    PDM_dgroup_entity_transpose(dmne->n_group_elmt,
                                dmne->dgroup_elmt_idx,
                                dmne->dgroup_elmt,
                (PDM_g_num_t *) distrib_elt,
                                &delt_group_idx,
                                &delt_group,
                                dmn->comm);
  }

  PDM_g_num_t shift = 0;
  for(int i_section = 0; i_section < n_section; ++i_section) {

    int id_section = sections_id[i_section];
    const PDM_g_num_t    *delmt_distribution = PDM_DMesh_nodal_distrib_section_get(dmn, geom_kind, id_section);
    int                   n_elt              = PDM_DMesh_nodal_section_n_elt_get  (dmn, geom_kind, id_section);
    PDM_g_num_t          *dconnec            = PDM_DMesh_nodal_section_std_get    (dmn, geom_kind, id_section);
    PDM_Mesh_nodal_elt_t  t_elt              = PDM_DMesh_nodal_section_type_get   (dmn, geom_kind, id_section);

    int         *dconnec_idx    = (int         * ) malloc( (n_elt+1) * sizeof(int        ));
    PDM_g_num_t *delmt_ln_to_gn = (PDM_g_num_t * ) malloc( (n_elt  ) * sizeof(PDM_g_num_t));

    int strid = PDM_Mesh_nodal_n_vtx_elt_get(t_elt, order);
    dconnec_idx[0] = 0;
    for(int i = 0; i < n_elt; ++i) {
      dconnec_idx[i+1] = dconnec_idx[i] + strid;
      delmt_ln_to_gn[i] = delmt_distribution[i_rank] + i + 1;
    }

    log_trace("section %d (%d) :\n", i_section, id_section);
    PDM_log_trace_array_long(delmt_ln_to_gn, n_elt, "  delmt_ln_to_gn : ");

    PDM_g_num_t *pvtx_ln_to_gn;
    int         *pcell_vtx_idx;
    int         *pcell_vtx;
    int          pn_vtx;
    PDM_part_dconnectivity_to_pconnectivity_sort_single_part(dmn->comm,
                                                             delmt_distribution,
                                                             dconnec_idx,
                                                             dconnec,
                                                             n_elt,
                                    (const PDM_g_num_t *)    delmt_ln_to_gn,
                                                            &pn_vtx,
                                                            &pvtx_ln_to_gn,
                                                            &pcell_vtx_idx,
                                                            &pcell_vtx);

    /*
     * Coordinates
     */
    PDM_g_num_t *vtx_distrib = PDM_dmesh_nodal_vtx_distrib_get(dmn);
    double      *dvtx_coord  = PDM_DMesh_nodal_vtx_get(dmn);
    // int          dn_vtx   = PDM_DMesh_nodal_n_vtx_get(dln->dmesh_nodal_in);
    // assert(dn_vtx == (vtx_distrib[i_rank+1]-vtx_distrib[i_rank]));
    double** tmp_pvtx_coord = NULL;
    PDM_part_dcoordinates_to_pcoordinates(dmn->comm,
                                          1,
                                          vtx_distrib,
                                          dvtx_coord,
                                          &pn_vtx,
                   (const PDM_g_num_t **) &pvtx_ln_to_gn,
                                          &tmp_pvtx_coord);

    double* pvtx_coord_out = tmp_pvtx_coord[0];

    /*
     * Groups
     */
    if (dmne != NULL) {
      for(int i = 0; i < n_elt; ++i) {
        delmt_ln_to_gn[i] += shift;
      }

      int **tmp_elt_group_idx = NULL;
      int **tmp_elt_group     = NULL;
      PDM_part_dentity_group_to_pentity_group(dmn->comm,
                                              1,
                                              distrib_elt,
                                              delt_group_idx,
                                              delt_group,
                                              &n_elt,
                      (const PDM_g_num_t **)  &delmt_ln_to_gn,
                                              &tmp_elt_group_idx,
                                              &tmp_elt_group);
      int *pelt_group_idx = tmp_elt_group_idx[0];
      int *pelt_group     = tmp_elt_group    [0];
      PDM_log_trace_connectivity_int(pelt_group_idx, pelt_group, n_elt, "pelt_group : ");
      free (tmp_elt_group_idx);
      free (tmp_elt_group);
      PDM_log_trace_array_long(delmt_ln_to_gn, n_elt, "  delmt_ln_to_gn (shifted) : ");

      n_field = 1;
      field = malloc (sizeof(double *) * n_field);
      field[0] = malloc (sizeof(double) * n_elt);
      for (int i = 0; i < n_elt; i++) {
        assert (pelt_group_idx[i+1] == pelt_group_idx[i] + 1);
        field[0][i] = (double) pelt_group[i];
      }
      free (pelt_group);
      free (pelt_group_idx);
    }

    /*
     *  Dump
     */
    char filename[999];
    sprintf(filename, "%s_section_%2.2d_%2.2d.vtk", filename_patter, i_section, i_rank);
    if (order == 1) {
      PDM_vtk_write_std_elements_double(filename,
                                        pn_vtx,
                                        pvtx_coord_out,
                                        pvtx_ln_to_gn,
                                        t_elt,
                                        n_elt,
                                        pcell_vtx,
                                        delmt_ln_to_gn,
                                        n_field,
                                        (const char   **) &field_name,
                                        (const double **) field);
    } else {
      PDM_vtk_write_std_elements_ho(filename,
                                    order,
                                    pn_vtx,
                                    pvtx_coord_out,
                                    pvtx_ln_to_gn,
                                    t_elt,
                                    n_elt,
                                    pcell_vtx,
                                    delmt_ln_to_gn,
                                    n_field,
                                    (const char   **) &field_name,
                                    (const double **) field);
    }
    free(tmp_pvtx_coord);
    free(pvtx_ln_to_gn);
    free(pcell_vtx_idx);
    free(pcell_vtx);

    free(dconnec_idx);
    free(delmt_ln_to_gn);

    free(pvtx_coord_out);

    shift += delmt_distribution[n_rank];

    if (dmne != NULL) {
      free (field[0]);
      free (field);
    }
  }

  if (dmne != NULL) {
    free (delt_group_idx);
    free (delt_group);
    free (distrib_elt);
  }
}





static inline void _swap_rows
(
 const int   row1,
 const int   row2,
 double    **A,
 double     *b,
 const int   n,
 const int   stride
 )
{
  if (row1 != row2) {
    for (int j = 0; j < n; j++) {
      double tmp = A[row1][j];
      A[row1][j] = A[row2][j];
      A[row2][j] = tmp;
    }

    for (int j = 0; j < stride; j++) {
      double tmp = b[stride*row1 + j];
      b[stride*row1 + j] = b[stride*row2 + j];
      b[stride*row1 + j] = tmp;
    }
  }
}

/**
 * Solve the linear system Ax = b using Gaussian elimination,
 * where A is an n*n matrix and b, x are n*stride matrices (stored linearly)
 *
 */

static int
_gauss_elim
(
 double    **A,
 double     *b,
 double     *x,
 const int   n,
 const int   stride,
 const int   inplace
 )
{
  double **_A = A;
  double  *_b = b;

  if (!inplace) {
    _A = malloc (sizeof(double *) * n);
    for (int i = 0; i < n; i++) {
      _A[i] = malloc (sizeof(double) * n);
      memcpy(_A[i], A[i], sizeof(double) * n);
    }

    _b = malloc (sizeof(double) * n * stride);
    memcpy(_b, b, sizeof(double) * n * stride);
  }


  for (int i = 0; i < n; i++) {

    /* Find pivot */
    double amax = PDM_ABS(_A[i][i]);
    int imax = i;
    for (int k = i+1; k < n; k++) {
      double aki = PDM_ABS(_A[k][i]);
      if (aki > amax) {
        amax = aki;
        imax = k;
      }
    }

    if (amax < 1.e-15) {
      /* matrix A is singular */
      return 0;
    }

    /* Swap rows i and imax */
    _swap_rows (i, imax, _A, _b, n, stride);

    /* Eliminate */
    double iamax = 1. / _A[i][i];

    for (int k = i+1; k < n; k++) {
      double r = _A[k][i] * iamax;
      for (int j = i+1; j < n; j++) {
        _A[k][j] -= r * _A[i][j];
      }
      _A[k][i] = 0.;

      for (int j = 0; j < stride; j++) {
        _b[stride*k + j] -= r * _b[stride*i + j];
      }
    }

  }


  /* Solve triangular system */
  memcpy(x, _b, sizeof(double) * n * stride);

  for (int i = n-1; i >= 0; i--) {

    for (int j = i+1; j < n; j++) {
      for (int k = 0; k < stride; k++) {
        x[stride*i + k] -= x[stride*j + k] * _A[i][j];
      }
    }

    double iai = 1. / _A[i][i];
    for (int k = 0; k < stride; k++) {
      x[stride*i + k] *= iai;
    }
  }


  /* Free memory */
  if (!inplace) {
    for (int i = 0; i < n; i++) {
      free (_A[i]);
    }
    free (_A);
    free (_b);
  }

  return 1;
}


static inline double _pow(const double x, const int p) {
  if (p == 0) return 1.;

  double y = x;
  for (int i = 1; i < p; i++) {
    y *= x;
  }
  return y;
}

static double **
_bezier_matrix_bar
(
 const int order
 )
{
  int n_nodes = PDM_Mesh_nodal_n_vtx_elt_get(PDM_MESH_NODAL_BARHO, order);

  double **b = malloc (sizeof(double *) * n_nodes);
  for (int i = 0; i < n_nodes; i++) {
    b[i] = malloc (sizeof(double) * n_nodes);
  }

  // compute binomial coefficients
  int coef[order/2+1];
  coef[0] = 1;
  for (int n = 2; n <= order; n++) {

    if (n%2 == 0) coef[n/2] = coef[n/2-1];

    for (int k = n/2; k > 0; k--) {
      coef[k] += coef[k-1];
    }
  }

  double in = 1. / (double) order;

  b[0][0] = 1.;
  for (int j = 1; j <= order; j++) {
    b[0][j] = 0.;
  }

  for (int j = 0; j <= order; j++) {
    int c = coef[PDM_MIN(j,order-j)];
    for (int i = 1; i <= order/2; i++) {
      double u = i * in;
      b[i][j] = c * _pow(u,j) * _pow(1. - u, order-j);
    }
  }

  for (int i = order/2+1; i <= order; i++) {
    for (int j = 0; j <= order; j++) {
      b[i][j] = b[order-i][order-j];
    }
  }

  return b;
}


static double **
_bezier_matrix_quad
(
 const int order
 )
{
  int n_nodes = PDM_Mesh_nodal_n_vtx_elt_get(PDM_MESH_NODAL_QUADHO, order);

  double **A = _bezier_matrix_bar(order);

  double **B = malloc (sizeof(double *) * n_nodes);
  for (int i = 0; i < n_nodes; i++) {
    B[i] = malloc (sizeof(double) * n_nodes);
  }

  for (int i = 0; i <= order; i++) {
    for (int j = 0; j <= order; j++) {
      int k = i + (order+1)*j;
      for (int ii = 0; ii <= order; ii++) {
        for (int jj = 0; jj <= order; jj++) {
          int l = ii + (order+1)*jj;
          B[k][l] = A[i][ii] * A[j][jj];
        }
      }
    }
  }
  for (int i = 0; i <= order; i++) {
    free (A[i]);
  }
  free (A);
  return B;
}


static double **
_bezier_matrix_tria
(
 const int order
 )
{
#define ij2idx(i, j) ((i) + (j)*(order + 1 - (j)))
  int n_nodes = PDM_Mesh_nodal_n_vtx_elt_get(PDM_MESH_NODAL_TRIAHO, order);

  double **b = malloc (sizeof(double *) * n_nodes);
  for (int i = 0; i < n_nodes; i++) {
    b[i] = malloc (sizeof(double) * n_nodes);
  }

  // compute trinomial coefficients
  const int n_coef = (order/2 + 1)*(order + 1 - order/2);
  int coef[n_coef];
  coef[0] = 1;
  for (int i = 1; i < n_coef; i++) {
    coef[i] = 0;
  }

  for (int n = 1; n <= order; n++) {
    for (int j = n/2; j >=0; j--) {
      int idx = ij2idx(n-j,j);
      for (int i = n-j; i >= j; i--) {
        if (i > 0) {
          if (i > j) {
            coef[idx] += coef[ij2idx(i-1,j)];
          } else {
            coef[idx] += coef[ij2idx(i,j-1)];
          }
        }

        if (j > 0) {
          coef[idx] += coef[ij2idx(i,j-1)];
        }
        idx--;
      }
    }
  }

  double in = 1. / (double) order;

  int icol = 0;
  for (int j = 0; j <= order; j++) {
    for (int i = 0; i <= order-j; i++) {

      int c = coef[ij2idx(PDM_MAX(i,j), PDM_MIN(i,j))];

      int irow = 0;
      for (int l = 0; l <= order; l++) {
        double v = l*in;
        for (int k = 0; k <= order-l; k++) {
          double u = k*in;
          b[irow][icol] = c * _pow(u,i) * _pow(v,j) * _pow(1 - u - v, order - i - j);
          irow++;
        }
      }
      icol++;
    }
  }

#undef ij2idx

  return b;
}




static void
_lagrange_to_bezier_bar
(
 const int  order,
 double    *lag,
 double    *bez
 )
{
  int n_nodes = PDM_Mesh_nodal_n_vtx_elt_get(PDM_MESH_NODAL_BARHO, order);

  if (order == 1) {
    memcpy (bez, lag, sizeof(double) * n_nodes * 3);
  }

  else if (order == 2) {

    for (int j = 0; j < 3; j++) {
      bez[j] = lag[j];
    }

    for (int j = 0; j < 3; j++) {
      bez[3+j] = -0.5*lag[j] + 2*lag[3+j] - 0.5*lag[6+j];
    }

    for (int j = 0; j < 3; j++) {
      bez[6+j] = lag[6+j];
    }
  }

  else if (order == 3) {

    double f833 = 5. / 6.;
    double f333 = 1. / 3.;

    for (int j = 0; j < 3; j++) {
      bez[j] = lag[j];
    }

    for (int j = 0; j < 3; j++) {
      bez[3+j] = -f833*lag[j] + 3*lag[3+j] - 1.5*lag[6+j] + f333*lag[9+j];
    }

    for (int j = 0; j < 3; j++) {
      bez[6+j] = f333*lag[j] - 1.5*lag[3+j] + 3*lag[6+j] - f833*lag[9+j];
    }

    for (int j = 0; j < 3; j++) {
      bez[9+j] = lag[9+j];
    }
  }

  else {
    double **B = _bezier_matrix_bar(order);

    _gauss_elim (B, lag, bez, n_nodes, 3, 0);

    if (0) {
      printf("B = \n");
      for (int i = 0; i < n_nodes; i++) {
        for (int j = 0; j < n_nodes; j++) {
          printf("%3.3f ", B[i][j]);
        }
        printf("\n");
      }

      printf("lag = \n");
      for (int i = 0; i < n_nodes; i++) {
        for (int j = 0; j < 3; j++) {
          printf("%3.3f ", lag[3*i+j]);
        }
        printf("\n");
      }

      printf("bez = \n");
      for (int i = 0; i < n_nodes; i++) {
        for (int j = 0; j < 3; j++) {
          printf("%3.3f ", bez[3*i+j]);
        }
        printf("\n");
      }
    }

    for (int i = 0; i < n_nodes; i++) {
      free (B[i]);
    }
    free (B);
  }
}

static void
_lagrange_to_bezier_tria
(
 const int  order,
 double    *lag,
 double    *bez
 )
{
  int n_nodes = PDM_Mesh_nodal_n_vtx_elt_get(PDM_MESH_NODAL_TRIAHO, order);

  if (order == 1) {
    memcpy (bez, lag, sizeof(double) * n_nodes * 3);
  }

  else if (order == 2) {
    for (int j = 0; j < 3; j++) {
      bez[j] = lag[j];
    }

    for (int j = 0; j < 3; j++) {
      bez[3+j] = -0.5*lag[j] + 2*lag[3+j] - 0.5*lag[6+j];
    }

    for (int j = 0; j < 3; j++) {
      bez[6+j] = lag[6+j];
    }

    for (int j = 0; j < 3; j++) {
      bez[9+j] = -0.5*lag[j] + 2*lag[9+j] - 0.5*lag[15+j];
    }

    for (int j = 0; j < 3; j++) {
      bez[12+j] = -0.5*lag[6+j] + 2*lag[12+j] - 0.5*lag[15+j];
    }

    for (int j = 0; j < 3; j++) {
      bez[15+j] = lag[15+j];
    }
  }

  else if (order == 3) {
    double f5_6 = 5. / 6.;
    double f1_3 = 1. / 3.;

    for (int j = 0; j < 3; j++) {
      bez[j] = lag[j];
    }

    for (int j = 0; j < 3; j++) {
      bez[3+j] = -f5_6*lag[j] + 3*lag[3+j] - 1.5*lag[6+j] + f1_3*lag[9+j];
    }

    for (int j = 0; j < 3; j++) {
      bez[6+j] = f1_3*lag[j] - 1.5*lag[3+j] + 3*lag[6+j] - f5_6*lag[9+j];
    }

    for (int j = 0; j < 3; j++) {
      bez[9+j] = lag[9+j];
    }

    for (int j = 0; j < 3; j++) {
      bez[12+j] = -f5_6*lag[j] + 3*lag[12+j] - 1.5*lag[21+j] + f1_3*lag[27+j];
    }

    for (int j = 0; j < 3; j++) {
      bez[15+j] = f1_3*lag[j] - 0.75*lag[3+j] - 0.75*lag[6+j] + f1_3*lag[9+j] - 0.75*lag[12+j] + 4.5*lag[15+j] - 0.75*lag[18+j] - 0.75*lag[21+j] - 0.75*lag[24+j] + f1_3*lag[27+j];
    }

    for (int j = 0; j < 3; j++) {
      bez[18+j] = -f5_6*lag[9+j] + 3*lag[18+j] - 1.5*lag[24+j] + f1_3*lag[27+j];
    }

    for (int j = 0; j < 3; j++) {
      bez[21+j] = f1_3*lag[j] - 1.5*lag[12+j] + 3*lag[21+j] - f5_6*lag[27+j];
    }

    for (int j = 0; j < 3; j++) {
      bez[24+j] = f1_3*lag[9+j] - 1.5*lag[18+j] + 3*lag[24+j] - f5_6*lag[27+j];
    }

    for (int j = 0; j < 3; j++) {
      bez[27+j] = lag[27+j];
    }
  }

  else {
    double **B = _bezier_matrix_tria(order);

    _gauss_elim (B, lag, bez, n_nodes, 3, 0);

    if (0) {
      printf("B = \n");
      for (int i = 0; i < n_nodes; i++) {
        for (int j = 0; j < n_nodes; j++) {
          printf("%3.3f ", B[i][j]);
        }
        printf("\n");
      }

      printf("lag = \n");
      for (int i = 0; i < n_nodes; i++) {
        for (int j = 0; j < 3; j++) {
          printf("%3.3f ", lag[3*i+j]);
        }
        printf("\n");
      }

      printf("bez = \n");
      for (int i = 0; i < n_nodes; i++) {
        for (int j = 0; j < 3; j++) {
          printf("%3.3f ", bez[3*i+j]);
        }
        printf("\n");
      }
    }

    for (int i = 0; i < n_nodes; i++) {
      free (B[i]);
    }
    free (B);
  }
}


static void
_lagrange_to_bezier_quad
(
 const int  order,
 double    *lag,
 double    *bez
 )
{
  int n_nodes = PDM_Mesh_nodal_n_vtx_elt_get(PDM_MESH_NODAL_QUADHO, order);

  if (order == 1) {
    memcpy (bez, lag, sizeof(double) * n_nodes * 3);
  }

  else if (order == 2) {
    for (int j = 0; j < 3; j++) {
      bez[j] = lag[j];
    }

    for (int j = 0; j < 3; j++) {
      bez[3+j] = -0.5*lag[j] + 2*lag[3+j] - 0.5*lag[6+j];
    }

    for (int j = 0; j < 3; j++) {
      bez[6+j] = lag[6+j];
    }

    for (int j = 0; j < 3; j++) {
      bez[9+j] = -0.5*lag[j] + 2*lag[9+j] - 0.5*lag[18+j];
    }

    for (int j = 0; j < 3; j++) {
      bez[12+j] = 0.25*lag[j] - lag[3+j] + 0.25*lag[6+j] - lag[9+j] + 4*lag[12+j] - lag[15+j] + 0.25*lag[18+j] - lag[21+j] + 0.25*lag[24+j];
    }

    for (int j = 0; j < 3; j++) {
      bez[15+j] = -0.5*lag[6+j] + 2*lag[15+j] - 0.5*lag[24+j];
    }

    for (int j = 0; j < 3; j++) {
      bez[18+j] = lag[18+j];
    }

    for (int j = 0; j < 3; j++) {
      bez[21+j] = -0.5*lag[18+j] + 2*lag[21+j] - 0.5*lag[24+j];
    }

    for (int j = 0; j < 3; j++) {
      bez[24+j] = lag[24+j];
    }
  }

  else if (order == 3) {

    double f833 = 5. / 6.;
    double f333 = 1. / 3.;
    double f694 = 25./ 36.;
    double f278 = 5. / 18.;
    double f111 = 1. / 9.;

    for (int j = 0; j < 3; j++) {
      bez[j] = lag[j];
    }

    for (int j = 0; j < 3; j++) {
      bez[3+j] = -f833*lag[j] + 3*lag[3+j] - 1.5*lag[6+j] + f333*lag[9+j];
    }

    for (int j = 0; j < 3; j++) {
      bez[6+j] = f333*lag[j] - 1.5*lag[3+j] + 3*lag[6+j] - f833*lag[9+j];
    }

    for (int j = 0; j < 3; j++) {
      bez[9+j] = lag[9+j];
    }

    for (int j = 0; j < 3; j++) {
      bez[12+j] = -f833*lag[j] + 3*lag[12+j] - 1.5*lag[24+j] + f333*lag[36+j];
    }

    for (int j = 0; j < 3; j++) {
      bez[15+j] = f694*lag[j] - 2.5*lag[3+j] + 1.25*lag[6+j] - f278*lag[9+j] - 2.5*lag[12+j] + 9*lag[15+j] - 4.5*lag[18+j] + lag[21+j] + 1.25*lag[24+j] - 4.5*lag[27+j] + 2.25*lag[30+j] - 0.5*lag[33+j] - f278*lag[36+j] + lag[39+j] - 0.5*lag[42+j] + f111*lag[45+j];
    }

    for (int j = 0; j < 3; j++) {
      bez[18+j] = -f278*lag[j] + 1.25*lag[3+j] - 2.5*lag[6+j] + f694*lag[9+j] + lag[12+j] - 4.5*lag[15+j] + 9*lag[18+j] - 2.5*lag[21+j] - 0.5*lag[24+j] + 2.25*lag[27+j] - 4.5*lag[30+j] + 1.25*lag[33+j] + f111*lag[36+j] - 0.5*lag[39+j] + lag[42+j] - f278*lag[45+j];
    }

    for (int j = 0; j < 3; j++) {
      bez[21+j] = -f833*lag[9+j] + 3*lag[21+j] - 1.5*lag[33+j] + f333*lag[45+j];
    }

    for (int j = 0; j < 3; j++) {
      bez[24+j] = f333*lag[j] -1.5*lag[12+j] + 3*lag[24+j] - f833*lag[36+j];
    }

    for (int j = 0; j < 3; j++) {
      bez[27+j] = -f278*lag[j] + lag[3+j] - 0.5*lag[6+j] + f111*lag[9+j] + 1.25*lag[12+j] - 4.5*lag[15+j] + 2.25*lag[18+j] - 0.5*lag[21+j] - 2.5*lag[24+j] + 9*lag[27+j] - 4.5*lag[30+j] + lag[33+j] + f694*lag[36+j] - 2.5*lag[39+j] + 1.25*lag[42+j] - f278*lag[45+j];
    }

    for (int j = 0; j < 3; j++) {
      bez[30+j] = f111*lag[j] - 0.5*lag[3+j] + lag[6+j] - f278*lag[9+j] - 0.5*lag[12+j] + 2.25*lag[15+j] - 4.5*lag[18+j] + 1.25*lag[21+j] + lag[24+j] - 4.5*lag[27+j] + 9*lag[30+j] - 2.5*lag[33+j] - f278*lag[36+j] + 1.25*lag[39+j] - 2.5*lag[42+j] + f694*lag[45+j];
    }

    for (int j = 0; j < 3; j++) {
      bez[33+j] = f333*lag[9+j] - 1.5*lag[21+j] + 3*lag[33+j] - f833*lag[45+j];
    }

    for (int j = 0; j < 3; j++) {
      bez[36+j] = lag[36+j];
    }

    for (int j = 0; j < 3; j++) {
      bez[39+j] = -f833*lag[36+j] + 3*lag[39+j] - 1.5*lag[42+j] + f333*lag[45+j];
    }

    for (int j = 0; j < 3; j++) {
      bez[42+j] = f333*lag[36+j] - 1.5*lag[39+j] + 3*lag[42+j] - f833*lag[45+j];
    }

    for (int j = 0; j < 3; j++) {
      bez[45+j] = lag[45+j];
    }
  }

  else {
    double **B = _bezier_matrix_quad(order);

    _gauss_elim (B, lag, bez, n_nodes, 3, 0);

    for (int i = 0; i < n_nodes; i++) {
      free (B[i]);
    }
    free (B);
  }
}



static void _check_tetra_face
(
 const char   *filename,
 const int     order,
 const double *coord
 )
{
  #define ij2idx(i, j) (i) + (order+1)*(j) - ((j)-1)*(j)/2

  int n_node = PDM_Mesh_nodal_n_vtx_elt_get(PDM_MESH_NODAL_TRIAHO, order);

  FILE *f = fopen(filename, "w");

  fprintf(f, "# vtk DataFile Version 2.0\n");
  fprintf(f, "mesh\n");
  fprintf(f, "ASCII\n");
  fprintf(f, "DATASET POLYDATA\n");

  fprintf(f, "POINTS %d double\n", n_node);
  for (int i = 0; i < n_node; i++) {
    for (int j = 0; j < 3; j++) {
      fprintf(f, "%.20lf ", coord[3*i+j]);
    }
    fprintf(f, "\n");
  }

  int n_tri = order*order;
  fprintf(f, "POLYGONS %d %d\n", n_tri, 4*n_tri);
  for (int j = 0; j < order; j++) {
    for (int i = 0; i < order-j; i++) {
      fprintf(f, "3 %d %d %d\n", ij2idx(i,j), ij2idx(i+1,j), ij2idx(i,j+1));
    }
  }

  for (int j = 0; j < order-1; j++) {
    for (int i = 0; i < order-j-1; i++) {
      fprintf(f, "3 %d %d %d\n", ij2idx(i,j+1), ij2idx(i+1,j), ij2idx(i+1,j+1));
    }
  }

  fclose(f);

  #undef ij2idx
}

static void
_lagrange_to_bezier_tetra
(
 const int  order,
 double    *lag,
 double    *bez
 )
{
  #define ijk2idx(i, j, k) ((i) + (j)*(order + 1 - (k)) - (j)*((j)-1)/2 + ((k)*((k)*((k) - 3*order - 6) + 3*order*(order + 4) + 11)) / 6)

  int n_nodes = PDM_Mesh_nodal_n_vtx_elt_get(PDM_MESH_NODAL_TETRAHO, order);

  if (order == 1) {
    memcpy (bez, lag, sizeof(double) * n_nodes * 3);
  }

  else {
    int n_nodes_tria = PDM_Mesh_nodal_n_vtx_elt_get(PDM_MESH_NODAL_TRIAHO, order);

    double *tria_lag = (double *) malloc(sizeof(double) * n_nodes_tria * 3);
    double *tria_bez = (double *) malloc(sizeof(double) * n_nodes_tria * 3);

    double **B = _bezier_matrix_tria(order);

    int idx;


    // hack for internal nodes
    memcpy (bez, lag, sizeof(double) * n_nodes * 3);

    // face w = 0
    // idx = 0;
    // for (int j = 0; j <= order; j++) {
    //   for (int i = 0; i <= order - j; i++) {
    //     int idx2 = ijk2idx(i,j,0);
    //     for (int l = 0; l < 3; l++) {
    //       tria_lag[idx++] = lag[3*idx2+l];
    //     }
    //   }
    // }
    if(0 == 1) {
      _check_tetra_face("face_w0.vtk", order, tria_lag);
    }
    _gauss_elim (B, lag, bez, n_nodes_tria, 3, 0);
    // _gauss_elim (B, tria_lag, tria_bez, n_nodes_tria, 3, 0);

    // idx = 0;
    // for (int j = 0; j <= order; j++) {
    //   for (int i = 0; i <= order - j; i++) {
    //     int idx2 = ijk2idx(i,j,0);
    //     for (int l = 0; l < 3; l++) {
    //       bez[3*idx2+l] = tria_bez[idx++];
    //     }
    //   }
    // }

    // face u = 0
    idx = 0;
    for (int k = 0; k <= order; k++) {
      for (int j = 0; j <= order - k; j++) {
        int idx2 = ijk2idx(0,j,k);
        for (int l = 0; l < 3; l++) {
          tria_lag[idx++] = lag[3*idx2+l];
        }
      }
    }

    _gauss_elim (B, tria_lag, tria_bez, n_nodes_tria, 3, 0);

    idx = 0;
    for (int k = 0; k <= order; k++) {
      for (int j = 0; j <= order - k; j++) {
        int idx2 = ijk2idx(0,j,k);
        for (int l = 0; l < 3; l++) {
          bez[3*idx2+l] = tria_bez[idx++];
        }
      }
    }


    // face v = 0
    idx = 0;
    for (int k = 0; k <= order; k++) {
      for (int i = 0; i <= order - k; i++) {
        int idx2 = ijk2idx(i,0,k);
        for (int l = 0; l < 3; l++) {
          tria_lag[idx++] = lag[3*idx2+l];
        }
      }
    }

    _gauss_elim (B, tria_lag, tria_bez, n_nodes_tria, 3, 0);

    idx = 0;
    for (int k = 0; k <= order; k++) {
      for (int i = 0; i <= order - k; i++) {
        int idx2 = ijk2idx(i,0,k);
        for (int l = 0; l < 3; l++) {
          bez[3*idx2+l] = tria_bez[idx++];
        }
      }
    }

    // face u+v+w = 1
    idx = 0;
    for (int k = 0; k <= order; k++) {
      for (int j = 0; j <= order - k; j++) {
        int i = order - j - k;
        int idx2 = ijk2idx(i,j,k);
        for (int l = 0; l < 3; l++) {
          tria_lag[idx++] = lag[3*idx2+l];
        }
      }
    }

    _gauss_elim (B, tria_lag, tria_bez, n_nodes_tria, 3, 0);

    idx = 0;
    for (int k = 0; k <= order; k++) {
      for (int j = 0; j <= order - k; j++) {
        int i = order - j - k;
        int idx2 = ijk2idx(i,j,k);
        for (int l = 0; l < 3; l++) {
          bez[3*idx2+l] = tria_bez[idx++];
        }
      }
    }

    for (int i = 0; i < n_nodes_tria; i++) {
      free (B[i]);
    }
    free (B);
    free (tria_lag);
    free (tria_bez);
  }
  #undef ijk2idx
}


static void
_lagrange_to_bezier_pyramid
(
 const int  order,
 double    *lag,
 double    *bez
 )
{
  #define ijk2idx(i, j, k) ((i) + (j)*(order+1-(k)) + ((k)*((k)*(2*(k) - 6*order - 9) + 6*order*(order + 3) + 13)) / 6)

  int n_nodes = PDM_Mesh_nodal_n_vtx_elt_get(PDM_MESH_NODAL_PYRAMIDHO, order);

  if (order == 1) {
    memcpy (bez, lag, sizeof(double) * n_nodes * 3);
  }

  else {

    int n_nodes_tria = PDM_Mesh_nodal_n_vtx_elt_get(PDM_MESH_NODAL_TRIAHO, order);
    int n_nodes_quad = PDM_Mesh_nodal_n_vtx_elt_get(PDM_MESH_NODAL_QUADHO, order);

    double *tria_lag = (double *) malloc(sizeof(double) * n_nodes_tria * 3);
    double *tria_bez = (double *) malloc(sizeof(double) * n_nodes_tria * 3);
    double *quad_lag = (double *) malloc(sizeof(double) * n_nodes_quad * 3);
    double *quad_bez = (double *) malloc(sizeof(double) * n_nodes_quad * 3);

    double **B_tria = _bezier_matrix_tria(order);
    double **B_quad = _bezier_matrix_quad(order);

    int idx;


    // hack for internal nodes
    memcpy (bez, lag, sizeof(double) * n_nodes * 3);


    // face w = 0
    _gauss_elim (B_quad, lag, bez, n_nodes_quad, 3, 0);


    // face u = 0
    idx = 0;
    for (int k = 0; k <= order; k++) {
      for (int j = 0; j <= order-k; j++) {
        int idx2 = ijk2idx(0,j,k);
        for (int l = 0; l < 3; l++) {
          tria_lag[idx++] = lag[3*idx2+l];
        }
      }
    }

    _gauss_elim (B_tria, tria_lag, tria_bez, n_nodes_tria, 3, 0);

    idx = 0;
    for (int k = 0; k <= order; k++) {
      for (int j = 0; j <= order-k; j++) {
        int idx2 = ijk2idx(0,j,k);
        for (int l = 0; l < 3; l++) {
          bez[3*idx2+l] = tria_bez[idx++];
        }
      }
    }


    // face u = 1-w
    idx = 0;
    for (int k = 0; k <= order; k++) {
      for (int j = 0; j <= order-k; j++) {
        int idx2 = ijk2idx(order-k,j,k);
        for (int l = 0; l < 3; l++) {
          tria_lag[idx++] = lag[3*idx2+l];
        }
      }
    }

    _gauss_elim (B_tria, tria_lag, tria_bez, n_nodes_tria, 3, 0);

    idx = 0;
    for (int k = 0; k <= order; k++) {
      for (int j = 0; j <= order-k; j++) {
        int idx2 = ijk2idx(order-k,j,k);
        for (int l = 0; l < 3; l++) {
          bez[3*idx2+l] = tria_bez[idx++];
        }
      }
    }


    // face v = 0
    idx = 0;
    for (int k = 0; k <= order; k++) {
      for (int i = 0; i <= order-k; i++) {
        int idx2 = ijk2idx(i,0,k);
        for (int l = 0; l < 3; l++) {
          tria_lag[idx++] = lag[3*idx2+l];
        }
      }
    }

    _gauss_elim (B_tria, tria_lag, tria_bez, n_nodes_tria, 3, 0);

    idx = 0;
    for (int k = 0; k <= order; k++) {
      for (int i = 0; i <= order-k; i++) {
        int idx2 = ijk2idx(i,0,k);
        for (int l = 0; l < 3; l++) {
          bez[3*idx2+l] = tria_bez[idx++];
        }
      }
    }


    // face v = 1-w
    idx = 0;
    for (int k = 0; k <= order; k++) {
      for (int i = 0; i <= order-k; i++) {
        int idx2 = ijk2idx(i,order-k,k);
        for (int l = 0; l < 3; l++) {
          tria_lag[idx++] = lag[3*idx2+l];
        }
      }
    }

    _gauss_elim (B_tria, tria_lag, tria_bez, n_nodes_tria, 3, 0);

    idx = 0;
    for (int k = 0; k <= order; k++) {
      for (int i = 0; i <= order-k; i++) {
        int idx2 = ijk2idx(i,order-k,k);
        for (int l = 0; l < 3; l++) {
          bez[3*idx2+l] = tria_bez[idx++];
        }
      }
    }



    for (int i = 0; i < n_nodes_tria; i++) {
      free (B_tria[i]);
    }
    free (B_tria);
    free (tria_lag);
    free (tria_bez);

    for (int i = 0; i < n_nodes_quad; i++) {
      free (B_quad[i]);
    }
    free (B_quad);
    free (quad_lag);
    free (quad_bez);
  }

  #undef ijk2idx
}



static void
_lagrange_to_bezier_prism
(
 const int  order,
 double    *lag,
 double    *bez
 )
{
  #define ijk2idx(i, j, k) ((i) + (j)*(order+1) - (j)*((j)-1)/2 + (k)*(order+1)*(order+2)/2)

  int n_nodes = PDM_Mesh_nodal_n_vtx_elt_get(PDM_MESH_NODAL_PRISMHO, order);

  if (order == 1) {
    memcpy (bez, lag, sizeof(double) * n_nodes * 3);
  }

  else {

    int n_nodes_tria = PDM_Mesh_nodal_n_vtx_elt_get(PDM_MESH_NODAL_TRIAHO, order);
    int n_nodes_quad = PDM_Mesh_nodal_n_vtx_elt_get(PDM_MESH_NODAL_QUADHO, order);

    double *tria_lag = (double *) malloc(sizeof(double) * n_nodes_tria * 3);
    double *tria_bez = (double *) malloc(sizeof(double) * n_nodes_tria * 3);
    double *quad_lag = (double *) malloc(sizeof(double) * n_nodes_quad * 3);
    double *quad_bez = (double *) malloc(sizeof(double) * n_nodes_quad * 3);

    double **B_tria = _bezier_matrix_tria(order);
    double **B_quad = _bezier_matrix_quad(order);

    int idx;

    // hack for internal nodes
    memcpy (bez, lag, sizeof(double) * n_nodes * 3);

    // face w = 0
    _gauss_elim (B_tria, lag, bez, n_nodes_tria, 3, 0);

    // face w = 1
    _gauss_elim (B_tria, &lag[order*n_nodes_tria], &bez[order*n_nodes_tria], n_nodes_tria, 3, 0);


    // face u = 0
    idx = 0;
    for (int k = 0; k <= order; k++) {
      for (int j = 0; j <= order; j++) {
        int idx2 = ijk2idx(0,j,k);
        for (int l = 0; l < 3; l++) {
          quad_lag[idx++] = lag[3*idx2+l];
        }
      }
    }

    _gauss_elim (B_quad, quad_lag, quad_bez, n_nodes_quad, 3, 0);

    idx = 0;
    for (int k = 0; k <= order; k++) {
      for (int j = 0; j <= order; j++) {
        int idx2 = ijk2idx(0,j,k);
        for (int l = 0; l < 3; l++) {
          bez[3*idx2+l] = quad_bez[idx++];
        }
      }
    }


    // face v = 0
    idx = 0;
    for (int k = 0; k <= order; k++) {
      for (int i = 0; i <= order; i++) {
        int idx2 = ijk2idx(i,0,k);
        for (int l = 0; l < 3; l++) {
          quad_lag[idx++] = lag[3*idx2+l];
        }
      }
    }

    _gauss_elim (B_quad, quad_lag, quad_bez, n_nodes_quad, 3, 0);

    idx = 0;
    for (int k = 0; k <= order; k++) {
      for (int i = 0; i <= order; i++) {
        int idx2 = ijk2idx(i,0,k);
        for (int l = 0; l < 3; l++) {
          bez[3*idx2+l] = quad_bez[idx++];
        }
      }
    }


    // face u+v = 1
    idx = 0;
    for (int k = 0; k <= order; k++) {
      for (int i = 0; i <= order; i++) {
        int idx2 = ijk2idx(i,order-i,k);
        for (int l = 0; l < 3; l++) {
          quad_lag[idx++] = lag[3*idx2+l];
        }
      }
    }

    _gauss_elim (B_quad, quad_lag, quad_bez, n_nodes_quad, 3, 0);

    idx = 0;
    for (int k = 0; k <= order; k++) {
      for (int i = 0; i <= order; i++) {
        int idx2 = ijk2idx(i,order-i,k);
        for (int l = 0; l < 3; l++) {
          bez[3*idx2+l] = quad_bez[idx++];
        }
      }
    }


    for (int i = 0; i < n_nodes_tria; i++) {
      free (B_tria[i]);
    }
    free (B_tria);
    free (tria_lag);
    free (tria_bez);

    for (int i = 0; i < n_nodes_quad; i++) {
      free (B_quad[i]);
    }
    free (B_quad);
    free (quad_lag);
    free (quad_bez);

  }

  #undef ijk2idx
}




static void
_lagrange_to_bezier_hexa
(
 const int  order,
 double    *lag,
 double    *bez
 )
{
  #define ijk2idx(i, j, k) ((i) + (order+1)*((j) + (order+1)*(k)))

  int n_nodes = PDM_Mesh_nodal_n_vtx_elt_get(PDM_MESH_NODAL_HEXAHO, order);

  if (order == 1) {
    memcpy (bez, lag, sizeof(double) * n_nodes * 3);
  }

  else {
    int n_nodes_quad = PDM_Mesh_nodal_n_vtx_elt_get(PDM_MESH_NODAL_QUADHO, order);

    double *quad_lag = (double *) malloc(sizeof(double) * n_nodes_quad * 3);
    double *quad_bez = (double *) malloc(sizeof(double) * n_nodes_quad * 3);

    double **B = _bezier_matrix_quad(order);

    int idx;


    // hack for internal nodes
    memcpy (bez, lag, sizeof(double) * n_nodes * 3);

    // face w = 0
    _gauss_elim (B, lag, bez, n_nodes_quad, 3, 0);

    // face w = 1
    _gauss_elim (B, &lag[order*n_nodes_quad], &bez[order*n_nodes_quad], n_nodes_quad, 3, 0);


    // face u = 0
    idx = 0;
    for (int k = 0; k <= order; k++) {
      for (int j = 0; j <= order; j++) {
        int idx2 = ijk2idx(0,j,k);
        for (int l = 0; l < 3; l++) {
          quad_lag[idx++] = lag[3*idx2+l];
        }
      }
    }

    _gauss_elim (B, quad_lag, quad_bez, n_nodes_quad, 3, 0);

    idx = 0;
    for (int k = 0; k <= order; k++) {
      for (int j = 0; j <= order; j++) {
        int idx2 = ijk2idx(0,j,k);
        for (int l = 0; l < 3; l++) {
          bez[3*idx2+l] = quad_bez[idx++];
        }
      }
    }


    // face u = 1
    idx = 0;
    for (int k = 0; k <= order; k++) {
      for (int j = 0; j <= order; j++) {
        int idx2 = ijk2idx(order,j,k);
        for (int l = 0; l < 3; l++) {
          quad_lag[idx++] = lag[3*idx2+l];
        }
      }
    }

    _gauss_elim (B, quad_lag, quad_bez, n_nodes_quad, 3, 0);

    idx = 0;
    for (int k = 0; k <= order; k++) {
      for (int j = 0; j <= order; j++) {
        int idx2 = ijk2idx(order,j,k);
        for (int l = 0; l < 3; l++) {
          bez[3*idx2+l] = quad_bez[idx++];
        }
      }
    }


    // face v = 0
    idx = 0;
    for (int i = 0; i <= order; i++) {
      for (int k = 0; k <= order; k++) {
        int idx2 = ijk2idx(i,0,k);
        for (int l = 0; l < 3; l++) {
          quad_lag[idx++] = lag[3*idx2+l];
        }
      }
    }

    _gauss_elim (B, quad_lag, quad_bez, n_nodes_quad, 3, 0);

    idx = 0;
    for (int i = 0; i <= order; i++) {
      for (int k = 0; k <= order; k++) {
        int idx2 = ijk2idx(i,0,k);
        for (int l = 0; l < 3; l++) {
          bez[3*idx2+l] = quad_bez[idx++];
        }
      }
    }


    // face v = 1
    idx = 0;
    for (int i = 0; i <= order; i++) {
      for (int k = 0; k <= order; k++) {
        int idx2 = ijk2idx(i,order,k);
        for (int l = 0; l < 3; l++) {
          quad_lag[idx++] = lag[3*idx2+l];
        }
      }
    }

    _gauss_elim (B, quad_lag, quad_bez, n_nodes_quad, 3, 0);

    idx = 0;
    for (int i = 0; i <= order; i++) {
      for (int k = 0; k <= order; k++) {
        int idx2 = ijk2idx(i,order,k);
        for (int l = 0; l < 3; l++) {
          bez[3*idx2+l] = quad_bez[idx++];
        }
      }
    }

    for (int i = 0; i < n_nodes_quad; i++) {
      free (B[i]);
    }
    free (B);
    free (quad_lag);
    free (quad_bez);

  }

  #undef ijk2idx
}


static void
_bezier_bounding_boxes
(
 PDM_dmesh_nodal_t   *dmn,
 int                  order,
 PDM_geometry_kind_t  geom_kind,
 const char          *filename_patter
 )
{
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(dmn->comm, &i_rank);
  PDM_MPI_Comm_size(dmn->comm, &n_rank);

  int *sections_id = PDM_DMesh_nodal_sections_id_get(dmn, geom_kind);
  int  n_section   = PDM_DMesh_nodal_n_section_get  (dmn, geom_kind);

  PDM_g_num_t shift = 0;
  for (int i_section = 0; i_section < n_section; ++i_section) {

    int id_section = sections_id[i_section];
    const PDM_g_num_t    *delmt_distribution = PDM_DMesh_nodal_distrib_section_get(dmn, geom_kind, id_section);
    int                   n_elt              = PDM_DMesh_nodal_section_n_elt_get  (dmn, geom_kind, id_section);
    PDM_g_num_t          *dconnec            = PDM_DMesh_nodal_section_std_get    (dmn, geom_kind, id_section);
    PDM_Mesh_nodal_elt_t  t_elt              = PDM_DMesh_nodal_section_type_get   (dmn, geom_kind, id_section);

    if (t_elt != PDM_MESH_NODAL_BAR2      &&
        t_elt != PDM_MESH_NODAL_TRIA3     &&
        t_elt != PDM_MESH_NODAL_QUAD4     &&
        t_elt != PDM_MESH_NODAL_BARHO     &&
        t_elt != PDM_MESH_NODAL_TRIAHO    &&
        t_elt != PDM_MESH_NODAL_QUADHO    &&
        t_elt != PDM_MESH_NODAL_TETRA4    &&
        t_elt != PDM_MESH_NODAL_TETRAHO   &&
        t_elt != PDM_MESH_NODAL_PYRAMID5  &&
        t_elt != PDM_MESH_NODAL_PYRAMIDHO &&
        t_elt != PDM_MESH_NODAL_PRISM6    &&
        t_elt != PDM_MESH_NODAL_PRISMHO   &&
        t_elt != PDM_MESH_NODAL_HEXA8     &&
        t_elt != PDM_MESH_NODAL_HEXAHO) continue;

    int         *dconnec_idx    = (int         * ) malloc( (n_elt+1) * sizeof(int        ));
    PDM_g_num_t *delmt_ln_to_gn = (PDM_g_num_t * ) malloc( (n_elt  ) * sizeof(PDM_g_num_t));

    int strid = PDM_Mesh_nodal_n_vtx_elt_get(t_elt, order);
    dconnec_idx[0] = 0;
    for(int i = 0; i < n_elt; ++i) {
      dconnec_idx[i+1] = dconnec_idx[i] + strid;
      delmt_ln_to_gn[i] = delmt_distribution[i_rank] + i + 1;
    }

    int *ijk_to_vtk = PDM_ho_ordering_ijk_to_user_get("PDM_HO_ORDERING_VTK",
                                                      t_elt,
                                                      order);


    PDM_g_num_t *pvtx_ln_to_gn;
    int         *pcell_vtx_idx;
    int         *pcell_vtx;
    int          pn_vtx;
    PDM_part_dconnectivity_to_pconnectivity_sort_single_part(dmn->comm,
                                                             delmt_distribution,
                                                             dconnec_idx,
                                                             dconnec,
                                                             n_elt,
                                    (const PDM_g_num_t *)    delmt_ln_to_gn,
                                                            &pn_vtx,
                                                            &pvtx_ln_to_gn,
                                                            &pcell_vtx_idx,
                                                            &pcell_vtx);

    /*
     * Coordinates
     */
    PDM_g_num_t *vtx_distrib = PDM_dmesh_nodal_vtx_distrib_get(dmn);
    double      *dvtx_coord  = PDM_DMesh_nodal_vtx_get(dmn);
    // int          dn_vtx   = PDM_DMesh_nodal_n_vtx_get(dln->dmesh_nodal_in);
    // assert(dn_vtx == (vtx_distrib[i_rank+1]-vtx_distrib[i_rank]));
    double** tmp_pvtx_coord = NULL;
    PDM_part_dcoordinates_to_pcoordinates(dmn->comm,
                                          1,
                                          vtx_distrib,
                                          dvtx_coord,
                                          &pn_vtx,
                   (const PDM_g_num_t **) &pvtx_ln_to_gn,
                                          &tmp_pvtx_coord);

    double* pvtx_coord_out = tmp_pvtx_coord[0];

    int n_nodes = PDM_Mesh_nodal_n_vtx_elt_get(t_elt, order);
    double *lagrange_coord = malloc (sizeof(double) * n_nodes * 3);
    double *bezier_coord   = malloc (sizeof(double) * n_nodes * 3);
    double *elt_coord      = malloc (sizeof(double) * n_elt * n_nodes * 3);
    int    *elt_vtx        = malloc (sizeof(int)    * n_elt * n_nodes);

    double *extents = malloc (sizeof(double) * n_elt * 6);
    int idx2 = 0;
    for (int i = 0; i < n_elt; i++) {
      double *_min = extents + 6*i;
      double *_max = _min + 3;

      for (int j = 0; j < 3; j++) {
        _min[j] =  1e30;
        _max[j] = -1e30;
      }

      int idx = 0;
      for (int k = pcell_vtx_idx[i]; k < pcell_vtx_idx[i+1]; k++) {
        int ivtx = pcell_vtx[k] - 1;
        for (int j = 0; j < 3; j++) {
          lagrange_coord[idx++] = pvtx_coord_out[3*ivtx + j];
        }
      }

      if (t_elt == PDM_MESH_NODAL_BAR2 ||
          t_elt == PDM_MESH_NODAL_BARHO) {
        _lagrange_to_bezier_bar (order, lagrange_coord, bezier_coord);
      }
      else if (t_elt == PDM_MESH_NODAL_TRIA3 ||
               t_elt == PDM_MESH_NODAL_TRIAHO) {
        _lagrange_to_bezier_tria (order, lagrange_coord, bezier_coord);
      }
      else if (t_elt == PDM_MESH_NODAL_QUAD4 ||
               t_elt == PDM_MESH_NODAL_QUADHO) {
        _lagrange_to_bezier_quad (order, lagrange_coord, bezier_coord);
      }
      else if (t_elt == PDM_MESH_NODAL_TETRA4 ||
               t_elt == PDM_MESH_NODAL_TETRAHO) {
        _lagrange_to_bezier_tetra (order, lagrange_coord, bezier_coord);
      }
      else if (t_elt == PDM_MESH_NODAL_PYRAMID5 ||
               t_elt == PDM_MESH_NODAL_PYRAMIDHO) {
        _lagrange_to_bezier_pyramid (order, lagrange_coord, bezier_coord);
      }
      else if (t_elt == PDM_MESH_NODAL_PRISM6 ||
               t_elt == PDM_MESH_NODAL_PRISMHO) {
        _lagrange_to_bezier_prism (order, lagrange_coord, bezier_coord);
      }
      else if (t_elt == PDM_MESH_NODAL_HEXA8 ||
               t_elt == PDM_MESH_NODAL_HEXAHO) {
        _lagrange_to_bezier_hexa (order, lagrange_coord, bezier_coord);
      }

      for (int k = 0; k < n_nodes; k++) {
        for (int j = 0; j < 3; j++) {
          elt_coord[3*idx2 + j] = bezier_coord[3*k + j];//lagrange_coord[3*k + j];//
          _min[j] = _MIN(_min[j], bezier_coord[3*k + j]);
          _max[j] = _MAX(_max[j], bezier_coord[3*k + j]);
        }
        elt_vtx[n_nodes*i + ijk_to_vtk[k]] = ++idx2;
      }
    }

    /*
     *  Dump
     */
    for(int i = 0; i < n_elt; ++i) {
      delmt_ln_to_gn[i] += shift;
    }

    char filename[999];
    sprintf(filename, "%s_bezier_section_%2.2d_%2.2d.vtk", filename_patter, i_section, i_rank);
    PDM_vtk_write_std_elements_ho(filename,
                                  order,
                                  n_elt * n_nodes,
                                  elt_coord,
                                  NULL,
                                  t_elt,
                                  n_elt,
                                  elt_vtx,
                                  delmt_ln_to_gn,
                                  0,
                                  NULL,
                                  NULL);
    free(elt_vtx);
    free(elt_coord);

    sprintf(filename, "%s_boxes_section_%2.2d_%2.2d.vtk", filename_patter, i_section, i_rank);
    PDM_vtk_write_boxes(filename,
                        n_elt,
                        extents,
                        delmt_ln_to_gn);
    free (extents);
    free (bezier_coord);
    free (lagrange_coord);

    free(tmp_pvtx_coord);
    free(pvtx_ln_to_gn);
    free(pcell_vtx_idx);
    free(pcell_vtx);

    free(dconnec_idx);
    free(delmt_ln_to_gn);

    free(pvtx_coord_out);

    shift += delmt_distribution[n_rank];
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
  PDM_g_num_t          nx     = 10;
  PDM_g_num_t          ny     = 10;
  PDM_g_num_t          nz     = 10;
  int                  order  = 1;
  double               length = 1.;
  int                  n_part = 1;
  int                  post   = 0;
  PDM_Mesh_nodal_elt_t t_elt  = PDM_MESH_NODAL_TRIA3;
  //  2 -> tria
  //  3 -> quad
  //  5 -> tetra
  //  6 -> pyramid
  //  7 -> prism
  //  8 -> hexa
  // 11 -> tria_ho
  // 12 -> quad_ho
  // 13 -> tetra_ho
  // 14 -> pyramid_ho
  // 15 -> prism_ho
  // 16 -> hexa_ho

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
             &nx,
             &ny,
             &nz,
             &order,
             (int *) &t_elt,
             &length,
             &n_part,
             &post,
             (int *) &method);

  if (t_elt == PDM_MESH_NODAL_TRIA3    ||
      t_elt == PDM_MESH_NODAL_QUAD4    ||
      t_elt == PDM_MESH_NODAL_TETRA4   ||
      t_elt == PDM_MESH_NODAL_PYRAMID5 ||
      t_elt == PDM_MESH_NODAL_PRISM6   ||
      t_elt == PDM_MESH_NODAL_HEXA8) {
    if (order != 1) {
      PDM_error(__FILE__, __LINE__, 0, "Invalid order %d for linear element type %d\n", order, (int) t_elt);
    }
  }

  /*
   *  Init
   */
  struct timeval t_elaps_debut;

  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;
  int i_rank;
  int n_rank;

  PDM_MPI_Init(&argc, &argv);
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  int dim = 3;
  if (t_elt == PDM_MESH_NODAL_TRIA3  ||
      t_elt == PDM_MESH_NODAL_QUAD4  ||
      t_elt == PDM_MESH_NODAL_TRIAHO ||
      t_elt == PDM_MESH_NODAL_QUADHO) {
    dim = 2;
  }

  if (order > 3) {
    int *ijk = NULL;

    for (PDM_Mesh_nodal_elt_t type = PDM_MESH_NODAL_BARHO;
         type <= PDM_MESH_NODAL_HEXAHO;
         type++) {

      if (type == PDM_MESH_NODAL_PYRAMIDHO) continue;

      ijk = PDM_vtk_lagrange_to_ijk(type, order);
      PDM_ho_ordering_user_to_ijk_add ("PDM_HO_ORDERING_VTK",
                                       type,
                                       order,
                                       PDM_Mesh_nodal_n_vtx_elt_get(type, order),
                                       ijk);
      free (ijk);
    }
  }

  /*
   *  Create distributed cube
   */
  PDM_dcube_nodal_t *dcube = PDM_dcube_nodal_gen_create(comm,
                                                        nx,
                                                        ny,
                                                        nz,
                                                        length,
                                                        0.,
                                                        0.,
                                                        0.,
                                                        t_elt,
                                                        order,
                                                        PDM_OWNERSHIP_KEEP);

  /*PDM_dcube_nodal_gen_ordering_set (dcube,
    "PDM_HO_ORDERING_VTK");*/

  PDM_dcube_nodal_gen_build (dcube);

  PDM_dmesh_nodal_t* dmn = PDM_dcube_nodal_gen_dmesh_nodal_get(dcube);

  PDM_dmesh_nodal_generate_distribution(dmn);


  /* Deform */
  double R[3][3] = {{0.9362934, -0.2896295, 0.1986693},
                    {0.3129918,  0.9447025, -0.0978434},
                    {-0.1593451,  0.1537920,  0.9751703}};

  PDM_g_num_t *vtx_distrib = PDM_dmesh_nodal_vtx_distrib_get(dmn);
  int dn_vtx = vtx_distrib[i_rank+1] - vtx_distrib[i_rank];
  double *dvtx_coord  = PDM_DMesh_nodal_vtx_get(dmn);
  double amplitude = 0.1;//0.07;
  double frequence = 4.;

  if (1) {
    for (int i = 0; i < dn_vtx; i++) {
      double x = (dvtx_coord[3*i    ] - 0.5) / length;
      double y = (dvtx_coord[3*i + 1] - 0.5) / length;
      double z = (dvtx_coord[3*i + 2] - 0.5) / length;

      //double scale = length * pow(2, order-1);

      if (dim == 2) {
        //dvtx_coord[3*i + 2] = scale * (pow(x, order) + pow(y, order));
        dvtx_coord[3*i + 2] = length * (x*x + y*y);
      } else {
        dvtx_coord[3*i    ] += amplitude*length*cos(frequence*y);
        dvtx_coord[3*i + 1] += amplitude*length*cos(frequence*z);
        dvtx_coord[3*i + 2] += amplitude*length*cos(frequence*x);
      }
    }

    if (1) {
      for (int i = 0; i < dn_vtx; i++) {
        double x = dvtx_coord[3*i  ];
        double y = dvtx_coord[3*i+1];
        double z = dvtx_coord[3*i+2];

        for (int j = 0; j < 3; j++) {
          dvtx_coord[3*i+j] = R[j][0]*x + R[j][1]*y + R[j][2]*z;
        }
      }
    }
  }

  if (t_elt > PDM_MESH_NODAL_HEXA8) {
    /* Bounding boxes */
    if (dim == 3) {
      _bezier_bounding_boxes(dmn, order, PDM_GEOMETRY_KIND_VOLUMIC, "out_volumic");
    }
    _bezier_bounding_boxes(dmn, order, PDM_GEOMETRY_KIND_SURFACIC, "out_surfacic");
    _bezier_bounding_boxes(dmn, order, PDM_GEOMETRY_KIND_RIDGE,    "out_ridge");

    /* Reorder */
    PDM_dmesh_nodal_reorder (dmn,
                             "PDM_HO_ORDERING_VTK",
                             order);
  }

  if (dim == 3) {
    _dmesh_nodal_dump_vtk(dmn, order, PDM_GEOMETRY_KIND_VOLUMIC, "out_volumic");
  }
  //PDM_dmesh_nodal_dump_vtk(dmn, PDM_GEOMETRY_KIND_VOLUMIC , "out_volumic");
  //PDM_dmesh_nodal_dump_vtk(dmn, PDM_GEOMETRY_KIND_SURFACIC, "out_surfacic");
  _dmesh_nodal_dump_vtk(dmn, order, PDM_GEOMETRY_KIND_SURFACIC, "out_surfacic");
  _dmesh_nodal_dump_vtk(dmn, order, PDM_GEOMETRY_KIND_RIDGE,    "out_ridge");
  _dmesh_nodal_dump_vtk(dmn, order, PDM_GEOMETRY_KIND_CORNER,   "out_corner");


  //PDM_dmesh_nodal_to_dmesh_free(dmntodm);
  gettimeofday(&t_elaps_debut, NULL);
  PDM_dcube_nodal_gen_free(dcube);


  PDM_MPI_Finalize();

  return 0;
}
