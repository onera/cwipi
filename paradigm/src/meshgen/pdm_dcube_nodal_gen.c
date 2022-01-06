#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "pdm_part.h"
#include "pdm_dmesh_nodal_priv.h"
#include "pdm_dcube_nodal_gen.h"
#include "pdm_dcube_nodal_gen_priv.h"
#include "pdm_ho_ordering.h"
#include "pdm_mpi.h"
#include "pdm_distrib.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_priv.h"
#include "pdm_logging.h"
#include "pdm_binary_search.h"

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Global variable
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/
static int
_get_n_sub_elt
(
 PDM_Mesh_nodal_elt_t t_elmt
)
{
  switch (t_elmt) {
    case PDM_MESH_NODAL_TRIA3:
    {
      return 2;
    }
    break;

    case PDM_MESH_NODAL_QUAD4:
    {
      return 1;
    }
    break;

    case PDM_MESH_NODAL_TETRA4:
    {
      return 5;
    }
    break;

    case PDM_MESH_NODAL_PYRAMID5:
    {
      return 3;
    }
    break;

    case PDM_MESH_NODAL_PRISM6:
    {
      return 2;
    }
    break;

    case PDM_MESH_NODAL_HEXA8:
    {
      return 1;
    }
    break;

    default :
      PDM_error(__FILE__, __LINE__, 0, "Unknown element type\n");
      break;

  }
  return -1;
}


static
PDM_g_num_t
_get_n_cell_abs
(
PDM_g_num_t          n_hexa,
PDM_Mesh_nodal_elt_t t_elmt
)
{
  return _get_n_sub_elt(t_elmt) * n_hexa;
}


static inline PDM_g_num_t
sub2ind
(
 PDM_dcube_nodal_t *dcube,
 PDM_g_num_t        icell,
 PDM_g_num_t        jcell,
 PDM_g_num_t        kcell,
 int                i,
 int                j,
 int                k
)
{
  return 1 + dcube->order*icell + i + (dcube->order*dcube->nx + 1)*(dcube->order*jcell + j + (dcube->order*dcube->ny + 1)*(dcube->order*kcell + k));
}



static
void
_generate_tria_surf
(
 PDM_dcube_nodal_t *dcube,
 PDM_dmesh_nodal_t *dmesh_nodal
)
{
  int n_rank, i_rank;
  PDM_MPI_Comm_rank(dcube->comm, &i_rank);
  PDM_MPI_Comm_size(dcube->comm, &n_rank);

  int order = dcube->order;

  int n_vtx_elt = PDM_Mesh_nodal_n_vtx_elt_get(PDM_MESH_NODAL_TRIA3, order);

  int dn_elt = 2 * dcube->dn_quad;

  /* Set up surface part */
  PDM_g_num_t *delt_vtx = (PDM_g_num_t *) malloc((n_vtx_elt * dn_elt) * sizeof(PDM_g_num_t));
  int idx = 0;
  for (int iquad = 0; iquad < dcube->dn_quad; iquad++) {
    PDM_g_num_t g = dcube->distrib_quad[i_rank] + iquad;

    PDM_g_num_t indi = g % dcube->nx;
    PDM_g_num_t indj = g / dcube->nx;

    for (int j = 0; j <= dcube->order; j++) {
      for (int i = 0; i <= dcube->order - j; i++) {
        delt_vtx[idx++] = sub2ind(dcube, indi, indj, 0, i, j, 0);
      }
    }

    for (int j = 0; j <= dcube->order; j++) {
      for (int i = 0; i <= dcube->order - j; i++) {
        delt_vtx[idx++] = sub2ind(dcube, indi+1, indj+1, 0, -i, -j, 0);
      }
    }

  }

  dmesh_nodal->surfacic->n_g_elmts = 2*dcube->distrib_quad[n_rank];

  int id_section = PDM_DMesh_nodal_elmts_section_add(dmesh_nodal->surfacic, PDM_MESH_NODAL_TRIA3);
  PDM_DMesh_nodal_elmts_section_std_set(dmesh_nodal->surfacic,
                                        id_section,
                                        dn_elt,
                                        delt_vtx,
                                        PDM_OWNERSHIP_KEEP);
}


static
void
_generate_quad_surf
(
 PDM_dcube_nodal_t *dcube,
 PDM_dmesh_nodal_t *dmesh_nodal
)
{
  int n_rank, i_rank;
  PDM_MPI_Comm_rank(dcube->comm, &i_rank);
  PDM_MPI_Comm_size(dcube->comm, &n_rank);

  int order = dcube->order;

  int n_vtx_elt = PDM_Mesh_nodal_n_vtx_elt_get(PDM_MESH_NODAL_QUAD4, order);

  int dn_elt = dcube->dn_quad;

  /* Set up surface part */
  PDM_g_num_t *delt_vtx = (PDM_g_num_t *) malloc((n_vtx_elt * dn_elt) * sizeof(PDM_g_num_t));
  int idx = 0;
  for (int iquad = 0; iquad < dcube->dn_quad; iquad++) {
    PDM_g_num_t g = dcube->distrib_quad[i_rank] + iquad;

    PDM_g_num_t indi = g % dcube->nx;
    PDM_g_num_t indj = g / dcube->nx;

    for (int j = 0; j <= dcube->order; j++) {
      for (int i = 0; i <= dcube->order; i++) {
        delt_vtx[idx++] = sub2ind(dcube, indi, indj, 0, i, j, 0);
      }
    }
  }

  dmesh_nodal->surfacic->n_g_elmts = dcube->distrib_quad[n_rank];

  int id_section = PDM_DMesh_nodal_elmts_section_add(dmesh_nodal->surfacic, PDM_MESH_NODAL_QUAD4);
  PDM_DMesh_nodal_elmts_section_std_set(dmesh_nodal->surfacic,
                                        id_section,
                                        dn_elt,
                                        delt_vtx,
                                        PDM_OWNERSHIP_KEEP);
}



static
void
_generate_tetra_vol
(
 PDM_dcube_nodal_t *dcube,
 PDM_dmesh_nodal_t *dmesh_nodal
)
{
  int n_rank, i_rank;
  PDM_MPI_Comm_rank(dcube->comm, &i_rank);
  PDM_MPI_Comm_size(dcube->comm, &n_rank);

  int order = dcube->order;

  int n_vtx_elt = PDM_Mesh_nodal_n_vtx_elt_get(PDM_MESH_NODAL_TETRA4, order);

  int dn_elt = 5 * dcube->dn_hexa;

  /* Set up volume part */
  PDM_g_num_t *delt_vtx = (PDM_g_num_t *) malloc((n_vtx_elt * dn_elt) * sizeof(PDM_g_num_t));
  int idx = 0;
  for (int ihexa = 0; ihexa < dcube->dn_hexa; ihexa++) {
    PDM_g_num_t g = dcube->distrib_hexa[i_rank] + ihexa;

    PDM_g_num_t indi = g % dcube->nx;
    PDM_g_num_t indj = ((g - indi) / dcube->nx) % dcube->ny;
    PDM_g_num_t indk = g / (dcube->nx * dcube->ny);

    if ((indi + indj + indk) % 2 == 0) {
      for (int k = 0; k <= dcube->order; k++) {
        for (int j = 0; j <= dcube->order - k; j++) {
          for (int i = 0; i <= dcube->order - j - k; i++) {
            delt_vtx[idx++] = sub2ind(dcube, indi, indj, indk, i, j, k);
          }
        }
      }

      for (int k = 0; k <= dcube->order; k++) {
        for (int j = 0; j <= dcube->order - k; j++) {
          for (int i = 0; i <= dcube->order - j - k; i++) {
            delt_vtx[idx++] = sub2ind(dcube, indi+1, indj, indk+1, -i, j, -k);
          }
        }
      }

      for (int k = 0; k <= dcube->order; k++) {
        for (int j = 0; j <= dcube->order - k; j++) {
          for (int i = 0; i <= dcube->order - j - k; i++) {
            delt_vtx[idx++] = sub2ind(dcube, indi, indj+1, indk+1, i, -j, -k);
          }
        }
      }

      for (int k = 0; k <= dcube->order; k++) {
        for (int j = 0; j <= dcube->order - k; j++) {
          for (int i = 0; i <= dcube->order - j - k; i++) {
            delt_vtx[idx++] = sub2ind(dcube, indi+1, indj+1, indk, -i, -j, k);
          }
        }
      }

      for (int k = 0; k <= dcube->order; k++) {
        for (int j = 0; j <= dcube->order - k; j++) {
          for (int i = 0; i <= dcube->order - j - k; i++) {
            delt_vtx[idx++] = sub2ind(dcube, indi, indj, indk, order-j-k, i+j, i+k);
          }
        }
      }
    }

    else {
      for (int k = 0; k <= dcube->order; k++) {
        for (int j = 0; j <= dcube->order - k; j++) {
          for (int i = 0; i <= dcube->order - j - k; i++) {
            delt_vtx[idx++] = sub2ind(dcube, indi, indj, indk, i+j+k, j, k);
          }
        }
      }

      for (int k = 0; k <= dcube->order; k++) {
        for (int j = 0; j <= dcube->order - k; j++) {
          for (int i = 0; i <= dcube->order - j - k; i++) {
            delt_vtx[idx++] = sub2ind(dcube, indi, indj, indk, order-j, order-k, i+j+k);
          }
        }
      }

      for (int k = 0; k <= dcube->order; k++) {
        for (int j = 0; j <= dcube->order - k; j++) {
          for (int i = 0; i <= dcube->order - j - k; i++) {
            delt_vtx[idx++] = sub2ind(dcube, indi, indj, indk, j, order-i, k);
          }
        }
      }

      for (int k = 0; k <= dcube->order; k++) {
        for (int j = 0; j <= dcube->order - k; j++) {
          for (int i = 0; i <= dcube->order - j - k; i++) {
            delt_vtx[idx++] = sub2ind(dcube, indi, indj, indk, i, j, i+j+k);
          }
        }
      }

      for (int k = 0; k <= dcube->order; k++) {
        for (int j = 0; j <= dcube->order - k; j++) {
          for (int i = 0; i <= dcube->order - j - k; i++) {
            delt_vtx[idx++] = sub2ind(dcube, indi, indj, indk, i+j, j+k, k+i);
          }
        }
      }
    }

  }


  dmesh_nodal->volumic->n_g_elmts = 5*dcube->distrib_hexa[n_rank];

  int id_section = PDM_DMesh_nodal_elmts_section_add(dmesh_nodal->volumic, PDM_MESH_NODAL_TETRA4);
  PDM_DMesh_nodal_elmts_section_std_set(dmesh_nodal->volumic,
                                        id_section,
                                        dn_elt,
                                        delt_vtx,
                                        PDM_OWNERSHIP_KEEP);
}


static
void
_generate_pyramid_vol
(
 PDM_dcube_nodal_t *dcube,
 PDM_dmesh_nodal_t *dmesh_nodal
)
{
  int n_rank, i_rank;
  PDM_MPI_Comm_rank(dcube->comm, &i_rank);
  PDM_MPI_Comm_size(dcube->comm, &n_rank);

  int order = dcube->order;

  int n_vtx_elt = PDM_Mesh_nodal_n_vtx_elt_get(PDM_MESH_NODAL_PYRAMID5, order);

  int dn_elt = 3 * dcube->dn_hexa;

  /* Set up volume part */
  PDM_g_num_t *delt_vtx = (PDM_g_num_t *) malloc((n_vtx_elt * dn_elt) * sizeof(PDM_g_num_t));
  int idx = 0;
  for (int ihexa = 0; ihexa < dcube->dn_hexa; ihexa++) {
    PDM_g_num_t g = dcube->distrib_hexa[i_rank] + ihexa;

    PDM_g_num_t indi = g % dcube->nx;
    PDM_g_num_t indj = ((g - indi) / dcube->nx) % dcube->ny;
    PDM_g_num_t indk = g / (dcube->nx * dcube->ny);

    int ind[3] = {indi, indj, indk};
    int u[3], v[3], w[3];

    if (indk%2 == 1) {
      if (indj%2 == 0) {
        if (indi%2 == 0) {
          u[0] =  1; u[1] =  0; u[2] =  0;
          v[0] =  0; v[1] =  1; v[2] =  0;
          w[0] =  0; w[1] =  0; w[2] =  1;
        } else {
          ind[2]++;
          u[0] =  0; u[1] =  0; u[2] = -1;
          v[0] =  0; v[1] =  1; v[2] =  0;
          w[0] =  1; w[1] =  0; w[2] =  0;
        }
      } else {
        if (indi%2 == 0) {
          ind[2]++;
          u[0] =  1; u[1] =  0; u[2] =  0;
          v[0] =  0; v[1] =  0; v[2] = -1;
          w[0] =  0; w[1] =  1; w[2] =  0;
        } else {
          ind[1]++;
          ind[2]++;
          u[0] =  0; u[1] = -1; u[2] =  0;
          v[0] =  0; v[1] =  0; v[2] = -1;
          w[0] =  1; w[1] =  0; w[2] =  0;
        }
      }
    } else {
      ind[0]++;
      if (indj%2 == 0) {
        if (indi%2 == 0) {
          u[0] =  0; u[1] =  0; u[2] =  1;
          v[0] =  0; v[1] =  1; v[2] =  0;
          w[0] = -1; w[1] =  0; w[2] =  0;
        } else {
          ind[2]++;
          u[0] = -1; u[1] =  0; u[2] =  0;
          v[0] =  0; v[1] =  1; v[2] =  0;
          w[0] =  0; w[1] =  0; w[2] = -1;
        }
      } else {
        if (indi%2 == 0) {
          ind[1]++;
          u[0] =  0; u[1] = -1; u[2] =  0;
          v[0] =  0; v[1] =  0; v[2] =  1;
          w[0] = -1; w[1] =  0; w[2] =  0;
        } else {
          u[0] = -1; u[1] =  0; u[2] =  0;
          v[0] =  0; v[1] =  0; v[2] =  1;
          w[0] =  0; w[1] =  1; w[2] =  0;
        }
      }
    }

    /*log_trace("indi,j,k = "PDM_FMT_G_NUM", "PDM_FMT_G_NUM", "PDM_FMT_G_NUM"\n", indi, indj, indk);
    log_trace(" u = (%d %d %d),  v = (%d %d %d),  w = (%d %d %d)\n",
              u[0], u[1], u[2],
              v[0], v[1], v[2],
              w[0], w[1], w[2]);*/

    // 1st sub-pyramid
    for (int k = 0; k <= dcube->order; k++) {
      for (int j = 0; j <= dcube->order - k; j++) {
        for (int i = 0; i <= dcube->order - k; i++) {
          delt_vtx[idx++] = sub2ind(dcube,
                                    ind[0] + u[0] + w[0],
                                    ind[1] + u[1] + w[1],
                                    ind[2] + u[2] + w[2],
                                    i*v[0] - j*w[0] - k*u[0],
                                    i*v[1] - j*w[1] - k*u[1],
                                    i*v[2] - j*w[2] - k*u[2]);
        }
      }
    }
    //PDM_log_trace_array_long(delt_vtx + idx - n_vtx_elt, n_vtx_elt, "  1st pyra: ");

    // 2nd sub-pyramid
    for (int k = 0; k <= dcube->order; k++) {
      for (int j = 0; j <= dcube->order - k; j++) {
        for (int i = 0; i <= dcube->order - k; i++) {
          delt_vtx[idx++] = sub2ind(dcube,
                                    ind[0],
                                    ind[1],
                                    ind[2],
                                    i*u[0] + j*v[0] + k*w[0],
                                    i*u[1] + j*v[1] + k*w[1],
                                    i*u[2] + j*v[2] + k*w[2]);
        }
      }
    }
    //PDM_log_trace_array_long(delt_vtx + idx - n_vtx_elt, n_vtx_elt, "  2nd pyra: ");

    // 3rd sub-pyramid
    for (int k = 0; k <= dcube->order; k++) {
      for (int j = 0; j <= dcube->order - k; j++) {
        for (int i = 0; i <= dcube->order - k; i++) {
          delt_vtx[idx++] = sub2ind(dcube,
                                    ind[0] + v[0] + w[0],
                                    ind[1] + v[1] + w[1],
                                    ind[2] + v[2] + w[2],
                                    -i*w[0] + j*u[0] - k*v[0],
                                    -i*w[1] + j*u[1] - k*v[1],
                                    -i*w[2] + j*u[2] - k*v[2]);
        }
      }
    }
    //PDM_log_trace_array_long(delt_vtx + idx - n_vtx_elt, n_vtx_elt, "  3rd pyra: ");
  }

  dmesh_nodal->volumic->n_g_elmts = 3*dcube->distrib_hexa[n_rank];

  int id_section = PDM_DMesh_nodal_elmts_section_add(dmesh_nodal->volumic, PDM_MESH_NODAL_PYRAMID5);
  PDM_DMesh_nodal_elmts_section_std_set(dmesh_nodal->volumic,
                                        id_section,
                                        dn_elt,
                                        delt_vtx,
                                        PDM_OWNERSHIP_KEEP);
}


static
void
_generate_prism_vol
(
 PDM_dcube_nodal_t *dcube,
 PDM_dmesh_nodal_t *dmesh_nodal
)
{
  int n_rank, i_rank;
  PDM_MPI_Comm_rank(dcube->comm, &i_rank);
  PDM_MPI_Comm_size(dcube->comm, &n_rank);

  int order = dcube->order;

  int n_vtx_elt = PDM_Mesh_nodal_n_vtx_elt_get(PDM_MESH_NODAL_PRISM6, order);

  int dn_elt = 2 * dcube->dn_hexa;

  /* Set up volume part */
  PDM_g_num_t *delt_vtx = (PDM_g_num_t *) malloc((n_vtx_elt * dn_elt) * sizeof(PDM_g_num_t));
  int idx = 0;
  for (int ihexa = 0; ihexa < dcube->dn_hexa; ihexa++) {
    PDM_g_num_t g = dcube->distrib_hexa[i_rank] + ihexa;

    PDM_g_num_t indi = g % dcube->nx;
    PDM_g_num_t indj = ((g - indi) / dcube->nx) % dcube->ny;
    PDM_g_num_t indk = g / (dcube->nx * dcube->ny);

    for (int k = 0; k <= dcube->order; k++) {
      for (int j = 0; j <= dcube->order; j++) {
        for (int i = 0; i <= dcube->order - j; i++) {
          delt_vtx[idx++] = sub2ind(dcube, indi, indj, indk, i, j, k);
        }
      }
    }

    for (int k = 0; k <= dcube->order; k++) {
      for (int j = 0; j <= dcube->order; j++) {
        for (int i = 0; i <= dcube->order - j; i++) {
          delt_vtx[idx++] = sub2ind(dcube, indi+1, indj+1, indk, -i, -j, k);
        }
      }
    }

  }

  dmesh_nodal->volumic->n_g_elmts = 2*dcube->distrib_hexa[n_rank];

  int id_section = PDM_DMesh_nodal_elmts_section_add(dmesh_nodal->volumic, PDM_MESH_NODAL_PRISM6);
  PDM_DMesh_nodal_elmts_section_std_set(dmesh_nodal->volumic,
                                        id_section,
                                        dn_elt,
                                        delt_vtx,
                                        PDM_OWNERSHIP_KEEP);
}


static
void
_generate_hexa_vol
(
 PDM_dcube_nodal_t *dcube,
 PDM_dmesh_nodal_t *dmesh_nodal
)
{
  int n_rank, i_rank;
  PDM_MPI_Comm_rank(dcube->comm, &i_rank);
  PDM_MPI_Comm_size(dcube->comm, &n_rank);

  int order = dcube->order;

  int n_vtx_elt = PDM_Mesh_nodal_n_vtx_elt_get(PDM_MESH_NODAL_HEXA8, order);

  int dn_elt = dcube->dn_hexa;

  /* Set up volume part */
  PDM_g_num_t *delt_vtx = (PDM_g_num_t *) malloc((n_vtx_elt * dn_elt) * sizeof(PDM_g_num_t));
  int idx = 0;
  for (int ihexa = 0; ihexa < dcube->dn_hexa; ihexa++) {
    PDM_g_num_t g = dcube->distrib_hexa[i_rank] + ihexa;

    PDM_g_num_t indi = g % dcube->nx;
    PDM_g_num_t indj = ((g - indi) / dcube->nx) % dcube->ny;
    PDM_g_num_t indk = g / (dcube->nx * dcube->ny);

    for (int k = 0; k <= dcube->order; k++) {
      for (int j = 0; j <= dcube->order; j++) {
        for (int i = 0; i <= dcube->order; i++) {
          delt_vtx[idx++] = sub2ind(dcube, indi, indj, indk, i, j, k);
        }
      }
    }

  }

  dmesh_nodal->volumic->n_g_elmts = dcube->distrib_hexa[n_rank];

  int id_section = PDM_DMesh_nodal_elmts_section_add(dmesh_nodal->volumic, PDM_MESH_NODAL_HEXA8);
  PDM_DMesh_nodal_elmts_section_std_set(dmesh_nodal->volumic,
                                        id_section,
                                        dn_elt,
                                        delt_vtx,
                                        PDM_OWNERSHIP_KEEP);
}


static
void
_generate_corners
(
 PDM_dcube_nodal_t *dcube,
 PDM_dmesh_nodal_t *dmesh_nodal
)
{
  int n_rank, i_rank;
  PDM_MPI_Comm_rank(dcube->comm, &i_rank);
  PDM_MPI_Comm_size(dcube->comm, &n_rank);

  int dim = dmesh_nodal->mesh_dimension;

  PDM_g_num_t gn_corner = 1 << dim;

  PDM_g_num_t *distrib_corner = PDM_compute_uniform_entity_distribution(dcube->comm, gn_corner);
  int dn_corner = (int) (distrib_corner[i_rank+1] - distrib_corner[i_rank]);

  PDM_g_num_t *delt_vtx = malloc (sizeof(PDM_g_num_t) * dn_corner);
  for (int icorner = 0; icorner < dn_corner; icorner++) {
    PDM_g_num_t g = distrib_corner[i_rank] + icorner;

    PDM_g_num_t i = g % 2;
    PDM_g_num_t k = g / 4;
    PDM_g_num_t j = (g - i - 4*k) / 2;

    delt_vtx[icorner] = sub2ind(dcube,
                                i*dcube->nx, j*dcube->ny, k*dcube->nz,
                                0, 0, 0);
  }
  free (distrib_corner);

  dmesh_nodal->corner->n_g_elmts = gn_corner;

  int id_section = PDM_DMesh_nodal_elmts_section_add(dmesh_nodal->corner, PDM_MESH_NODAL_POINT);
  PDM_DMesh_nodal_elmts_section_std_set(dmesh_nodal->corner,
                                        id_section,
                                        dn_corner,
                                        delt_vtx,
                                        PDM_OWNERSHIP_KEEP);
}


static
void
_generate_ridges
(
 PDM_dcube_nodal_t *dcube,
 PDM_dmesh_nodal_t *dmesh_nodal
)
{
  int n_rank, i_rank;
  PDM_MPI_Comm_rank(dcube->comm, &i_rank);
  PDM_MPI_Comm_size(dcube->comm, &n_rank);

  int order = dcube->order;

  int n_vtx_elt = PDM_Mesh_nodal_n_vtx_elt_get(PDM_MESH_NODAL_BAR2, order);

  int dn_elt = dcube->dn_bar;
  PDM_g_num_t gn_ridge = dcube->distrib_bar[n_rank];

  int dim = dmesh_nodal->mesh_dimension;
  int n_group = 4;
  if (dim == 3) {
    n_group = 12;
  }
  PDM_g_num_t group_idx[13] = {0};
  group_idx[1] = group_idx[0] + dcube->nx;
  group_idx[2] = group_idx[1] + dcube->nx;
  group_idx[3] = group_idx[2] + dcube->ny;
  group_idx[4] = group_idx[3] + dcube->ny;
  if (dim == 3) {
    group_idx[5] = group_idx[4] + dcube->nx;
    group_idx[6] = group_idx[5] + dcube->nx;
    group_idx[7] = group_idx[6] + dcube->ny;
    group_idx[8] = group_idx[7] + dcube->ny;
    for (int i = 9; i <= 12; i++) {
      group_idx[i] = group_idx[i-1] + dcube->nz;
    }
  }


  /* Set up ridge part */
  PDM_g_num_t *delt_vtx = (PDM_g_num_t *) malloc((n_vtx_elt * dn_elt) * sizeof(PDM_g_num_t));
  int idx = 0;
  for (int ielt = 0; ielt < dn_elt; ielt++) {

    PDM_g_num_t g = dcube->distrib_bar[i_rank] + ielt;
    PDM_g_num_t indi = 0, indj = 0, indk = 0;
    int dx = 0, dy = 0, dz = 0;

    if (g < group_idx[1]) {
      indi = g;
      indj = 0;
      indk = 0;
      dx = 1; dy = 0; dz = 0;
    }
    else if (g < group_idx[2]) {
      indi = g - group_idx[1] + 1;
      indj = dcube->ny;
      indk = 0;
      dx = -1; dy = 0; dz = 0;
    }
    else if (g < group_idx[3]) {
      indi = 0;
      indj = g - group_idx[2] + 1;
      indk = 0;
      dx = 0; dy = -1; dz = 0;
    }
    else if (g < group_idx[4]) {
      indi = dcube->nx;
      indj = g - group_idx[3];
      indk = 0;
      dx = 0; dy = 1; dz = 0;
    }

    else if (g < group_idx[5]) {
      indi = g - group_idx[4];
      indj = 0;
      indk = dcube->nz;
      dx = 1; dy = 0; dz = 0;
    }
    else if (g < group_idx[6]) {
      indi = g - group_idx[5];
      indj = dcube->ny;
      indk = dcube->nz;
      dx = 1; dy = 0; dz = 0;
    }
    else if (g < group_idx[7]) {
      indi = 0;
      indj = g - group_idx[6];
      indk = dcube->nz;
      dx = 0; dy = 1; dz = 0;
    }
    else if (g < group_idx[8]) {
      indi = dcube->nx;
      indj = g - group_idx[7];
      indk = dcube->nz;
      dx = 0; dy = 1; dz = 0;
    }

    else if (g < group_idx[9]) {
      indi = 0;
      indj = 0;
      indk = g - group_idx[8];
      dx = 0; dy = 0; dz = 1;
    }
    else if (g < group_idx[10]) {
      indi = dcube->nx;
      indj = 0;
      indk = g - group_idx[9];
      dx = 0; dy = 0; dz = 1;
    }
    else if (g < group_idx[11]) {
      indi = 0;
      indj = dcube->ny;
      indk = g - group_idx[10];
      dx = 0; dy = 0; dz = 1;
    }
    else if (g < group_idx[12]) {
      indi = dcube->nx;
      indj = dcube->ny;
      indk = g - group_idx[11];
      dx = 0; dy = 0; dz = 1;
    }

    for (int i = 0; i < n_vtx_elt; i++) {
      delt_vtx[idx++] = sub2ind(dcube, indi, indj, indk, i*dx, i*dy, i*dz);
    }
  }

  dmesh_nodal->ridge->n_g_elmts = gn_ridge;


  int id_section = PDM_DMesh_nodal_elmts_section_add(dmesh_nodal->ridge, PDM_MESH_NODAL_BAR2);
  PDM_DMesh_nodal_elmts_section_std_set(dmesh_nodal->ridge,
                                        id_section,
                                        dn_elt,
                                        delt_vtx,
                                        PDM_OWNERSHIP_KEEP);


  /* Groups */
  PDM_g_num_t *distrib[3];
  distrib[0] = PDM_compute_uniform_entity_distribution(dcube->comm, dcube->nx);
  distrib[1] = PDM_compute_uniform_entity_distribution(dcube->comm, dcube->ny);
  distrib[2] = NULL;
  if (dim == 3) {
    distrib[2] = PDM_compute_uniform_entity_distribution(dcube->comm, dcube->nz);
  }

  int dn[3] = {0};
  for (int i = 0; i < dim; i++) {
    dn[i] = (int) (distrib[i][i_rank+1] - distrib[i][i_rank]);
  }

  int *dgroup_elt_idx = malloc (sizeof(int) * (n_group + 1));
  dgroup_elt_idx[0] = 0;
  for (int k = 0; k < dim-1; k++) {
    for (int j = 0; j < 2; j++) {
      for (int i = 0; i < 2; i++) {
        int igroup = 4*k + 2*j + i;
        dgroup_elt_idx[igroup+1] = dgroup_elt_idx[igroup] + dn[j];
      }
    }
  }

  if (dim == 3) {
    for (int igroup = 9; igroup <= 12; igroup++) {
      dgroup_elt_idx[igroup] = dgroup_elt_idx[igroup-1] + dn[2];
    }
  }


  PDM_g_num_t *dgroup_elt = malloc (sizeof(PDM_g_num_t) * dgroup_elt_idx[n_group]);

  idx = 0;
  for (int k = 0; k < dim-1; k++) {
    for (int j = 0; j < 2; j++) {
      for (int i = 0; i < 2; i++) {
        int igroup = 4*k + 2*j + i;
        for (int ielt = 0; ielt < dn[j]; ielt++) {
          dgroup_elt[idx++] = group_idx[igroup] + distrib[j][i_rank] + ielt + 1;
        }
      }
    }
  }

  if (dim == 3) {
    for (int igroup = 8; igroup < 12; igroup++) {
      for (int ielt = 0; ielt < dn[2]; ielt++) {
        dgroup_elt[idx++] = group_idx[igroup] + distrib[2][i_rank] + ielt + 1;
      }
    }
  }

  if (0) {
    PDM_log_trace_connectivity_long(dgroup_elt_idx, dgroup_elt, n_group, "dgroup_elt : ");
  }

  PDM_DMesh_nodal_elmts_group_set(dmesh_nodal->ridge,
                                  n_group,
                                  dgroup_elt_idx,
                                  dgroup_elt,
                                  PDM_OWNERSHIP_KEEP);

  for (int i = 0; i < dim; i++) {
    free (distrib[i]);
  }
}


static inline int
_g_to_ijk_uv
(
 PDM_dcube_nodal_t *dcube,
 const PDM_g_num_t *group_idx,
 PDM_g_num_t        g,
 PDM_g_num_t       *ind,
 int               *u,
 int               *v
)
{
  PDM_g_num_t h;
  if (g < group_idx[1]) {
    h = g - group_idx[0];
    ind[0] = h % dcube->nx + 1;
    ind[1] = h / dcube->nx;
    ind[2] = 0;
    u[0] = -1; u[1] = 0; u[2] = 0;
    v[0] =  0; v[1] = 1; v[2] = 0;
    return 0;
  }
  else if (g < group_idx[2]) {
    h = g - group_idx[1];
    ind[0] = h % dcube->nx;
    ind[1] = h / dcube->nx;
    ind[2] = dcube->nz;
    u[0] = 1; u[1] = 0; u[2] = 0;
    v[0] = 0; v[1] = 1; v[2] = 0;
    return 1;
  }
  else if (g < group_idx[3]) {
    h = g - group_idx[2];
    ind[0] = 0;
    ind[1] = h % dcube->ny + 1;
    ind[2] = h / dcube->ny;
    u[0] = 0; u[1] = -1; u[2] = 0;
    v[0] = 0; v[1] =  0; v[2] = 1;
    return 2;
  }
  else if (g < group_idx[4]) {
    h = g - group_idx[3];
    ind[0] = dcube->nx;
    ind[1] = h % dcube->ny;
    ind[2] = h / dcube->ny;
    u[0] = 0; u[1] = 1; u[2] = 0;
    v[0] = 0; v[1] = 0; v[2] = 1;
    return 3;
  }
  else if (g < group_idx[5]) {
    h = g - group_idx[4];
    ind[0] = h / dcube->nz;
    ind[1] = 0;
    ind[2] = h % dcube->nz + 1;
    u[0] = 0; u[1] = 0; u[2] = -1;
    v[0] = 1; v[1] = 0; v[2] =  0;
    return 4;
  }
  else {//if (g < group_idx[6]) {
    h = g - group_idx[5];
    ind[0] = h / dcube->nz;
    ind[1] = dcube->ny;
    ind[2] = h % dcube->nz;
    u[0] = 0; u[1] = 0; u[2] = 1;
    v[0] = 1; v[1] = 0; v[2] = 0;
    return 5;
  }

  return -1;
}


static void
_set_surf_groups
(
  PDM_dcube_nodal_t    *dcube,
  PDM_dmesh_nodal_t    *dmesh_nodal,
  const int             n_group,
  const PDM_g_num_t    *group_idx,
  PDM_Mesh_nodal_elt_t  t_elt[3]
)
{
  int i_rank;
  PDM_MPI_Comm_rank(dcube->comm, &i_rank);

  int n_elt_quad[3];
  for (int i = 0; i < 3; i++) {
    n_elt_quad[i] = _get_n_sub_elt(t_elt[i]);
  }

  PDM_g_num_t *distrib[3];
  distrib[0] = PDM_compute_uniform_entity_distribution(dcube->comm,
                                                       n_elt_quad[0]*dcube->nx*dcube->ny);
  distrib[1] = PDM_compute_uniform_entity_distribution(dcube->comm,
                                                       n_elt_quad[1]*dcube->ny*dcube->nz);
  distrib[2] = PDM_compute_uniform_entity_distribution(dcube->comm,
                                                       n_elt_quad[2]*dcube->nz*dcube->nx);

  int dn[3];
  for (int i = 0; i < 3; i++) {
    dn[i] = (int) (distrib[i][i_rank+1] - distrib[i][i_rank]);
  }

  int *dgroup_elt_idx = malloc (sizeof(int) * (n_group + 1));
  dgroup_elt_idx[0] = 0;
  for (int j = 0; j < 3; j++) {
    for (int i = 0; i < 2; i++) {
      int igroup = 2*j + i;
      dgroup_elt_idx[igroup+1] = dgroup_elt_idx[igroup] + dn[j];
    }
  }

  PDM_g_num_t *dgroup_elt = malloc (sizeof(PDM_g_num_t) * dgroup_elt_idx[n_group]);
  int idx = 0;
  PDM_g_num_t shift = 0;
  for (int j = 0; j < 3; j++) {
    for (int i = 0; i < 2; i++) {
      int igroup = 2*j + i;
      for (int ielt = 0; ielt < dn[j]; ielt++) {
        dgroup_elt[idx++] = shift + distrib[j][i_rank] + ielt + 1;
      }
      shift += n_elt_quad[j]*(group_idx[igroup+1] - group_idx[igroup]);
    }
  }

  if (0) {
    PDM_log_trace_array_int(dgroup_elt_idx, n_group + 1, "dgroup_elt_idx : ");
    PDM_log_trace_connectivity_long(dgroup_elt_idx, dgroup_elt, n_group, "dgroup_elt : ");
  }

  PDM_DMesh_nodal_elmts_group_set(dmesh_nodal->surfacic,
                                  n_group,
                                  dgroup_elt_idx,
                                  dgroup_elt,
                                  PDM_OWNERSHIP_KEEP);

  for (int i = 0; i < 3; i++) {
    free (distrib[i]);
  }
}


static void
_generate_tetra_surf
(
 PDM_dcube_nodal_t *dcube,
 PDM_dmesh_nodal_t *dmesh_nodal
)
{
  int n_rank, i_rank;
  PDM_MPI_Comm_rank(dcube->comm, &i_rank);
  PDM_MPI_Comm_size(dcube->comm, &n_rank);

  int order = dcube->order;

  int n_vtx_elt = PDM_Mesh_nodal_n_vtx_elt_get(PDM_MESH_NODAL_TRIA3, order);

  dmesh_nodal->surfacic->n_g_elmts = 2*dcube->distrib_quad[n_rank];
  int dn_elt = 2 * dcube->dn_quad;

  const int n_group = 6;
  PDM_g_num_t group_idx[n_group+1];
  group_idx[0] = 0;
  group_idx[1] = group_idx[0] + dcube->nx*dcube->ny;
  group_idx[2] = group_idx[1] + dcube->nx*dcube->ny;
  group_idx[3] = group_idx[2] + dcube->ny*dcube->nz;
  group_idx[4] = group_idx[3] + dcube->ny*dcube->nz;
  group_idx[5] = group_idx[4] + dcube->nz*dcube->nx;
  group_idx[6] = group_idx[5] + dcube->nz*dcube->nx;

  /* Set up surface part */
  PDM_g_num_t *delt_vtx = (PDM_g_num_t *) malloc((n_vtx_elt * dn_elt) * sizeof(PDM_g_num_t));
  int idx = 0;
  for (int iquad = 0; iquad < dcube->dn_quad; iquad++) {

    PDM_g_num_t g = dcube->distrib_quad[i_rank] + iquad;

    PDM_g_num_t ind[3];
    int u[3], v[3];
    _g_to_ijk_uv (dcube, group_idx, g, ind, u, v);

    if ((ind[0] + ind[1] + ind[2]) % 2 == 0) {
      for (int j = 0; j <= dcube->order; j++) {
        for (int i = 0; i <= dcube->order - j; i++) {
          delt_vtx[idx++] = sub2ind(dcube,
                                    ind[0], ind[1], ind[2],
                                    i*u[0] + j*v[0],
                                    i*u[1] + j*v[1],
                                    i*u[2] + j*v[2]);
        }
      }

      for (int j = 0; j <= dcube->order; j++) {
        for (int i = 0; i <= dcube->order - j; i++) {
          delt_vtx[idx++] = sub2ind(dcube,
                                    ind[0], ind[1], ind[2],
                                    (order-i)*u[0] + (order-j)*v[0],
                                    (order-i)*u[1] + (order-j)*v[1],
                                    (order-i)*u[2] + (order-j)*v[2]);
        }
      }
    } else {
      for (int j = 0; j <= dcube->order; j++) {
        for (int i = 0; i <= dcube->order - j; i++) {
          delt_vtx[idx++] = sub2ind(dcube,
                                    ind[0], ind[1], ind[2],
                                    i*v[0] + (order-j)*u[0],
                                    i*v[1] + (order-j)*u[1],
                                    i*v[2] + (order-j)*u[2]);
        }
      }

      for (int j = 0; j <= dcube->order; j++) {
        for (int i = 0; i <= dcube->order - j; i++) {
          delt_vtx[idx++] = sub2ind(dcube,
                                    ind[0], ind[1], ind[2],
                                    (order-i)*v[0] + j*u[0],
                                    (order-i)*v[1] + j*u[1],
                                    (order-i)*v[2] + j*u[2]);
        }
      }
    }
  }


  int id_section = PDM_DMesh_nodal_elmts_section_add(dmesh_nodal->surfacic, PDM_MESH_NODAL_TRIA3);
  PDM_DMesh_nodal_elmts_section_std_set(dmesh_nodal->surfacic,
                                        id_section,
                                        dn_elt,
                                        delt_vtx,
                                        PDM_OWNERSHIP_KEEP);

  /* Groups */
  PDM_Mesh_nodal_elt_t t_elt_face[3] = {PDM_MESH_NODAL_TRIA3,
                                        PDM_MESH_NODAL_TRIA3,
                                        PDM_MESH_NODAL_TRIA3};
  _set_surf_groups (dcube,
                    dmesh_nodal,
                    n_group,
                    group_idx,
                    t_elt_face);
}



static void
_generate_pyramid_surf
(
 PDM_dcube_nodal_t *dcube,
 PDM_dmesh_nodal_t *dmesh_nodal
)
{
  int n_rank, i_rank;
  PDM_MPI_Comm_rank(dcube->comm, &i_rank);
  PDM_MPI_Comm_size(dcube->comm, &n_rank);

  int order = dcube->order;

  int n_vtx_tria = PDM_Mesh_nodal_n_vtx_elt_get(PDM_MESH_NODAL_TRIA3, order);
  int n_vtx_quad = PDM_Mesh_nodal_n_vtx_elt_get(PDM_MESH_NODAL_QUAD4, order);

  PDM_g_num_t gn_tria = 2*(dcube->nx*dcube->ny + dcube->ny*dcube->nz + dcube->nz*dcube->nx);
  PDM_g_num_t gn_quad = 0;

  if (dcube->nz%2 == 0) {
    gn_tria += 2*dcube->nx*dcube->ny;
  } else {
    gn_quad += dcube->nx*dcube->ny;
  }

  if (dcube->nx%2 == 0) {
    gn_tria += 2*dcube->ny*dcube->nz;
  } else {
    gn_quad += dcube->ny*dcube->nz;
  }

  if (dcube->ny%2 == 0) {
    gn_tria += 2*dcube->nz*dcube->nx;
  } else {
    gn_quad += dcube->nz*dcube->nx;
  }

  //log_trace("gn_tria = "PDM_FMT_G_NUM", gn_quad = "PDM_FMT_G_NUM"\n", gn_tria, gn_quad);

  dmesh_nodal->surfacic->n_g_elmts = gn_tria + gn_quad;

  PDM_g_num_t *distrib_tria = PDM_compute_uniform_entity_distribution(dcube->comm, gn_tria);
  PDM_g_num_t *distrib_quad = PDM_compute_uniform_entity_distribution(dcube->comm, gn_quad);

  int dn_tria = (int) (distrib_tria[i_rank+1] - distrib_tria[i_rank]);
  int dn_quad = (int) (distrib_quad[i_rank+1] - distrib_quad[i_rank]);

  const int n_group = 6;
  PDM_g_num_t group_idx[n_group+1];
  group_idx[0] = 0;
  group_idx[1] = group_idx[0] + dcube->nx*dcube->ny;
  group_idx[2] = group_idx[1] + dcube->nx*dcube->ny;
  group_idx[3] = group_idx[2] + dcube->ny*dcube->nz;
  group_idx[4] = group_idx[3] + dcube->ny*dcube->nz;
  group_idx[5] = group_idx[4] + dcube->nz*dcube->nx;
  group_idx[6] = group_idx[5] + dcube->nz*dcube->nx;

  int group_type[n_group];
  group_type[0] = 0;
  group_type[1] = (dcube->nz%2 == 1);
  group_type[2] = 0;
  group_type[3] = (dcube->nx%2 == 1);
  group_type[4] = 0;
  group_type[5] = (dcube->ny%2 == 1);

  /* Triangles */
  PDM_g_num_t *dtria_vtx = (PDM_g_num_t *) malloc((n_vtx_tria * dn_tria) * sizeof(PDM_g_num_t));
  int idx = 0;
  for (int itria = 0; itria < dn_tria; itria++) {

    PDM_g_num_t g = distrib_tria[i_rank] + itria;
    PDM_g_num_t g_quad = g / 2;

    if (g_quad >= group_idx[1] && dcube->nz%2 == 1) {
      g_quad += dcube->nx*dcube->ny;
    }
    if (g_quad >= group_idx[3] && dcube->nx%2 == 1) {
      g_quad += dcube->ny*dcube->nz;
    }

    PDM_g_num_t ind[3];
    int u[3], v[3];
    _g_to_ijk_uv (dcube, group_idx, g_quad, ind, u, v);


    if ((ind[0] + ind[1] + ind[2]) % 2 == 1) {
      if (g%2 == 0) {
        for (int j = 0; j <= dcube->order; j++) {
          for (int i = 0; i <= dcube->order - j; i++) {
            dtria_vtx[idx++] = sub2ind(dcube,
                                       ind[0], ind[1], ind[2],
                                       i*u[0] + j*v[0],
                                       i*u[1] + j*v[1],
                                       i*u[2] + j*v[2]);
          }
        }
      } else {
        for (int j = 0; j <= dcube->order; j++) {
          for (int i = 0; i <= dcube->order - j; i++) {
            dtria_vtx[idx++] = sub2ind(dcube,
                                       ind[0], ind[1], ind[2],
                                       (order-i)*u[0] + (order-j)*v[0],
                                       (order-i)*u[1] + (order-j)*v[1],
                                       (order-i)*u[2] + (order-j)*v[2]);
          }
        }
      }
    } else {
      if (g%2 == 0) {
        for (int j = 0; j <= dcube->order; j++) {
          for (int i = 0; i <= dcube->order - j; i++) {
            dtria_vtx[idx++] = sub2ind(dcube,
                                       ind[0], ind[1], ind[2],
                                       i*v[0] + (order-j)*u[0],
                                       i*v[1] + (order-j)*u[1],
                                       i*v[2] + (order-j)*u[2]);
          }
        }
      } else {
        for (int j = 0; j <= dcube->order; j++) {
          for (int i = 0; i <= dcube->order - j; i++) {
            dtria_vtx[idx++] = sub2ind(dcube,
                                       ind[0], ind[1], ind[2],
                                       (order-i)*v[0] + j*u[0],
                                       (order-i)*v[1] + j*u[1],
                                       (order-i)*v[2] + j*u[2]);
          }
        }
      }
    }
  }

  int id_tria = PDM_DMesh_nodal_elmts_section_add(dmesh_nodal->surfacic, PDM_MESH_NODAL_TRIA3);
  PDM_DMesh_nodal_elmts_section_std_set(dmesh_nodal->surfacic,
                                        id_tria,
                                        dn_tria,
                                        dtria_vtx,
                                        PDM_OWNERSHIP_KEEP);


  /* Quadrangles */
  PDM_g_num_t *dquad_vtx = (PDM_g_num_t *) malloc((n_vtx_quad * dn_quad) * sizeof(PDM_g_num_t));
  idx = 0;
  for (int iquad = 0; iquad < dn_quad; iquad++) {

    PDM_g_num_t g = group_idx[1] + distrib_quad[i_rank] + iquad;

    if (dcube->nz%2 == 0) {
      g += dcube->nx*dcube->ny;
    }
    if (g >= group_idx[2]) {
      g += dcube->ny*dcube->nz;
    }
    if (g >= group_idx[3] && dcube->nx%2 == 0) {
      g += dcube->ny*dcube->nz;
    }
    if (g >= group_idx[4]) {
      g += dcube->nz*dcube->nx;
    }

    PDM_g_num_t ind[3];
    int u[3], v[3];
    _g_to_ijk_uv (dcube, group_idx, g, ind, u, v);

    for (int j = 0; j <= dcube->order; j++) {
      for (int i = 0; i <= dcube->order; i++) {
        dquad_vtx[idx++] = sub2ind(dcube,
                                   ind[0], ind[1], ind[2],
                                   i*u[0] + j*v[0], i*u[1] + j*v[1], i*u[2] + j*v[2]);
      }
    }

  }

  int id_quad = PDM_DMesh_nodal_elmts_section_add(dmesh_nodal->surfacic, PDM_MESH_NODAL_QUAD4);
  PDM_DMesh_nodal_elmts_section_std_set(dmesh_nodal->surfacic,
                                        id_quad,
                                        dn_quad,
                                        dquad_vtx,
                                        PDM_OWNERSHIP_KEEP);

  free (distrib_tria);
  free (distrib_quad);


  /* Groups */
  PDM_g_num_t *distrib[n_group];
  PDM_g_num_t gn_elt_group;

  gn_elt_group = 2*dcube->nx*dcube->ny;
  distrib[0] = PDM_compute_uniform_entity_distribution(dcube->comm,
                                                       gn_elt_group);

  if (dcube->nz%2 == 1) {
    gn_elt_group /= 2;
  }
  distrib[1] = PDM_compute_uniform_entity_distribution(dcube->comm,
                                                       gn_elt_group);

  gn_elt_group = 2*dcube->ny*dcube->nz;
  distrib[2] = PDM_compute_uniform_entity_distribution(dcube->comm,
                                                       gn_elt_group);

  if (dcube->nx%2 == 1) {
    gn_elt_group /= 2;
  }
  distrib[3] = PDM_compute_uniform_entity_distribution(dcube->comm,
                                                       gn_elt_group);

  gn_elt_group = 2*dcube->nz*dcube->nx;
  distrib[4] = PDM_compute_uniform_entity_distribution(dcube->comm,
                                                       gn_elt_group);

  if (dcube->ny%2 == 1) {
    gn_elt_group /= 2;
  }
  distrib[5] = PDM_compute_uniform_entity_distribution(dcube->comm,
                                                       gn_elt_group);

  int dn[n_group];
  for (int i = 0; i < n_group; i++) {
    //log_trace("group %d ", i);
    //PDM_log_trace_array_long (distrib[i], n_rank+1, "distrib : ");
    dn[i] = (int) (distrib[i][i_rank+1] - distrib[i][i_rank]);
  }


  int *dgroup_elt_idx = malloc (sizeof(int) * (n_group + 1));
  dgroup_elt_idx[0] = 0;
  for (int i = 0; i < n_group; i++) {
    dgroup_elt_idx[i+1] = dgroup_elt_idx[i] + dn[i];
  }

  PDM_g_num_t *dgroup_elt = malloc (sizeof(PDM_g_num_t) * dgroup_elt_idx[n_group]);
  PDM_g_num_t shift[2] = {0, gn_tria};
  idx = 0;
  for (int i = 0; i < n_group; i++) {
    for (int ielt = 0; ielt < dn[i]; ielt++) {
      dgroup_elt[idx++] = shift[group_type[i]] + distrib[i][i_rank] + ielt + 1;
    }
    shift[group_type[i]] += distrib[i][n_rank];
  }

   if (0) {
    PDM_log_trace_array_int(dgroup_elt_idx, n_group + 1, "dgroup_elt_idx : ");
    PDM_log_trace_connectivity_long(dgroup_elt_idx, dgroup_elt, n_group, "dgroup_elt : ");
  }

  PDM_DMesh_nodal_elmts_group_set(dmesh_nodal->surfacic,
                                  n_group,
                                  dgroup_elt_idx,
                                  dgroup_elt,
                                  PDM_OWNERSHIP_KEEP);

  for (int i = 0; i < n_group; i++) {
    free (distrib[i]);
  }
}



static void
_generate_prism_surf
(
 PDM_dcube_nodal_t *dcube,
 PDM_dmesh_nodal_t *dmesh_nodal
)
{
  int n_rank, i_rank;
  PDM_MPI_Comm_rank(dcube->comm, &i_rank);
  PDM_MPI_Comm_size(dcube->comm, &n_rank);

  int order = dcube->order;

  int n_vtx_tria = PDM_Mesh_nodal_n_vtx_elt_get(PDM_MESH_NODAL_TRIA3, order);
  int n_vtx_quad = PDM_Mesh_nodal_n_vtx_elt_get(PDM_MESH_NODAL_QUAD4, order);

  PDM_g_num_t gn_tria = 4*dcube->nx*dcube->ny;
  PDM_g_num_t gn_quad = 2*(dcube->ny*dcube->nz + dcube->nz*dcube->nx);
  dmesh_nodal->surfacic->n_g_elmts = gn_tria + gn_quad;

  PDM_g_num_t *distrib_tria = PDM_compute_uniform_entity_distribution(dcube->comm, gn_tria);
  PDM_g_num_t *distrib_quad = PDM_compute_uniform_entity_distribution(dcube->comm, gn_quad);

  int dn_tria = (int) (distrib_tria[i_rank+1] - distrib_tria[i_rank]);
  int dn_quad = (int) (distrib_quad[i_rank+1] - distrib_quad[i_rank]);

  const int n_group = 6;
  PDM_g_num_t group_idx[n_group+1];
  group_idx[0] = 0;
  group_idx[1] = group_idx[0] + dcube->nx*dcube->ny;
  group_idx[2] = group_idx[1] + dcube->nx*dcube->ny;
  group_idx[3] = group_idx[2] + dcube->ny*dcube->nz;
  group_idx[4] = group_idx[3] + dcube->ny*dcube->nz;
  group_idx[5] = group_idx[4] + dcube->nz*dcube->nx;
  group_idx[6] = group_idx[5] + dcube->nz*dcube->nx;

  /* Set up surface part */
  /* Triangles */
  PDM_g_num_t *dtria_vtx = (PDM_g_num_t *) malloc((n_vtx_tria * dn_tria) * sizeof(PDM_g_num_t));
  int idx = 0;
  for (int itria = 0; itria < dn_tria; itria++) {

    PDM_g_num_t g = distrib_tria[i_rank] + itria;
    PDM_g_num_t g_quad = g / 2;

    PDM_g_num_t ind[3];
    int u[3], v[3];
    if (g_quad < group_idx[1]) {
      ind[0] = g_quad % dcube->nx;
      ind[1] = g_quad / dcube->nx;
      ind[2] = 0;
      u[0] = 0; u[1] = 1; u[2] = 0;
      v[0] = 1; v[1] = 0; v[2] = 0;
    }
    else {
      PDM_g_num_t h = g_quad - group_idx[1];
      ind[0] = h % dcube->nx;
      ind[1] = h / dcube->nx;
      ind[2] = dcube->nz;
      u[0] = 1; u[1] = 0; u[2] = 0;
      v[0] = 0; v[1] = 1; v[2] = 0;
    }

    if (g%2 == 0) {
      for (int j = 0; j <= dcube->order; j++) {
        for (int i = 0; i <= dcube->order - j; i++) {
          dtria_vtx[idx++] = sub2ind(dcube,
                                     ind[0], ind[1], ind[2],
                                     i*u[0] + j*v[0],
                                     i*u[1] + j*v[1],
                                     i*u[2] + j*v[2]);
        }
      }
    }
    else {
      for (int j = 0; j <= dcube->order; j++) {
        for (int i = 0; i <= dcube->order - j; i++) {
          dtria_vtx[idx++] = sub2ind(dcube,
                                     ind[0], ind[1], ind[2],
                                     (order-i)*u[0] + (order-j)*v[0],
                                     (order-i)*u[1] + (order-j)*v[1],
                                     (order-i)*u[2] + (order-j)*v[2]);
        }
      }
    }
  }

  int id_tria = PDM_DMesh_nodal_elmts_section_add(dmesh_nodal->surfacic, PDM_MESH_NODAL_TRIA3);
  PDM_DMesh_nodal_elmts_section_std_set(dmesh_nodal->surfacic,
                                        id_tria,
                                        dn_tria,
                                        dtria_vtx,
                                        PDM_OWNERSHIP_KEEP);

  /* Quadrangles */
  PDM_g_num_t *dquad_vtx = (PDM_g_num_t *) malloc((n_vtx_quad * dn_quad) * sizeof(PDM_g_num_t));
  idx = 0;
  for (int iquad = 0; iquad < dn_quad; iquad++) {

    PDM_g_num_t g = group_idx[2] + distrib_quad[i_rank] + iquad;

    PDM_g_num_t ind[3];
    int u[3], v[3];
    _g_to_ijk_uv (dcube, group_idx, g, ind, u, v);

    for (int j = 0; j <= dcube->order; j++) {
      for (int i = 0; i <= dcube->order; i++) {
        dquad_vtx[idx++] = sub2ind(dcube,
                                   ind[0], ind[1], ind[2],
                                   i*u[0] + j*v[0], i*u[1] + j*v[1], i*u[2] + j*v[2]);
      }
    }

  }

  int id_quad = PDM_DMesh_nodal_elmts_section_add(dmesh_nodal->surfacic, PDM_MESH_NODAL_QUAD4);
  PDM_DMesh_nodal_elmts_section_std_set(dmesh_nodal->surfacic,
                                        id_quad,
                                        dn_quad,
                                        dquad_vtx,
                                        PDM_OWNERSHIP_KEEP);

  free (distrib_tria);
  free (distrib_quad);

  /* Groups */
  PDM_Mesh_nodal_elt_t t_elt_face[3] = {PDM_MESH_NODAL_TRIA3,
                                        PDM_MESH_NODAL_QUAD4,
                                        PDM_MESH_NODAL_QUAD4};
  _set_surf_groups (dcube,
                    dmesh_nodal,
                    n_group,
                    group_idx,
                    t_elt_face);
}





static void
_generate_hexa_surf
(
 PDM_dcube_nodal_t *dcube,
 PDM_dmesh_nodal_t *dmesh_nodal
)
{
  int n_rank, i_rank;
  PDM_MPI_Comm_rank(dcube->comm, &i_rank);
  PDM_MPI_Comm_size(dcube->comm, &n_rank);

  int order = dcube->order;

  int n_vtx_elt = PDM_Mesh_nodal_n_vtx_elt_get(PDM_MESH_NODAL_QUAD4, order);

  dmesh_nodal->surfacic->n_g_elmts = dcube->distrib_quad[n_rank];
  int dn_elt = dcube->dn_quad;

  const int n_group = 6;
  PDM_g_num_t group_idx[n_group+1];
  group_idx[0] = 0;
  group_idx[1] = group_idx[0] + dcube->nx*dcube->ny;
  group_idx[2] = group_idx[1] + dcube->nx*dcube->ny;
  group_idx[3] = group_idx[2] + dcube->ny*dcube->nz;
  group_idx[4] = group_idx[3] + dcube->ny*dcube->nz;
  group_idx[5] = group_idx[4] + dcube->nz*dcube->nx;
  group_idx[6] = group_idx[5] + dcube->nz*dcube->nx;

  /* Set up surface part */
  PDM_g_num_t *delt_vtx = (PDM_g_num_t *) malloc((n_vtx_elt * dn_elt) * sizeof(PDM_g_num_t));
  int idx = 0;
  for (int iquad = 0; iquad < dcube->dn_quad; iquad++) {

    PDM_g_num_t g = dcube->distrib_quad[i_rank] + iquad;

    PDM_g_num_t ind[3];
    int u[3], v[3];
    _g_to_ijk_uv (dcube, group_idx, g, ind, u, v);

    for (int j = 0; j <= dcube->order; j++) {
      for (int i = 0; i <= dcube->order; i++) {
        delt_vtx[idx++] = sub2ind(dcube,
                                  ind[0], ind[1], ind[2],
                                  i*u[0] + j*v[0], i*u[1] + j*v[1], i*u[2] + j*v[2]);
      }
    }

  }



  int id_section = PDM_DMesh_nodal_elmts_section_add(dmesh_nodal->surfacic, PDM_MESH_NODAL_QUAD4);
  PDM_DMesh_nodal_elmts_section_std_set(dmesh_nodal->surfacic,
                                        id_section,
                                        dn_elt,
                                        delt_vtx,
                                        PDM_OWNERSHIP_KEEP);

  /* Groups */
  PDM_Mesh_nodal_elt_t t_elt_face[3] = {PDM_MESH_NODAL_QUAD4,
                                        PDM_MESH_NODAL_QUAD4,
                                        PDM_MESH_NODAL_QUAD4};
  _set_surf_groups (dcube,
                    dmesh_nodal,
                    n_group,
                    group_idx,
                    t_elt_face);
}

/*=============================================================================
 * Public function definitions
 *============================================================================*/


/**
 *
 * \brief Create a distributed cube
 *
 * \param [out]  id             dcube identifier
 * \param [in]   comm           Communicator
 * \param [in]   n_vtx_seg      Number of vertices in segments
 * \param [in]   length         Segment length
 * \param [in]   zero_x         Coordinates of the origin
 * \param [in]   zero_y         Coordinates of the origin
 * \param [in]   zero_z         Coordinates of the origin
 *
 */


PDM_dmesh_nodal_t *
PDM_dcube_nodal_gen_dmesh_nodal_get
(
 PDM_dcube_nodal_t *dcube
)
{
  return dcube->dmesh_nodal;
}







PDM_dcube_nodal_t *
PDM_dcube_nodal_gen_create
(
 PDM_MPI_Comm          comm,
 const PDM_g_num_t     nx,
 const PDM_g_num_t     ny,
 const PDM_g_num_t     nz,
 const double          length,
 const double          zero_x,
 const double          zero_y,
 const double          zero_z,
 PDM_Mesh_nodal_elt_t  t_elt,
 const int             order,
 PDM_ownership_t       owner
)
{
  int n_rank;
  int i_rank;

  PDM_MPI_Comm_size(comm, &n_rank);
  PDM_MPI_Comm_rank(comm, &i_rank);

  PDM_dcube_nodal_t *dcube = (PDM_dcube_nodal_t *) malloc(sizeof(PDM_dcube_nodal_t));

  dcube->comm   = comm;
  dcube->nx     = nx;
  dcube->ny     = ny;
  dcube->nz     = nz;
  dcube->length = length;
  dcube->zero_x = zero_x;
  dcube->zero_y = zero_y;
  dcube->zero_z = zero_z;
  dcube->t_elt  = t_elt;
  dcube->order  = order;
  dcube->owner  = owner;

  dcube->dmesh_nodal = NULL;
  dcube->ordering = NULL;
  dcube->distrib_bar  = NULL;
  dcube->distrib_quad = NULL;
  dcube->distrib_hexa = NULL;

  return dcube;
}


/**
 *
 * \brief Free a distributed cube
 *
 * \param [in]  id            dcube identifier
 *
 */

void
PDM_dcube_nodal_gen_free
(
 PDM_dcube_nodal_t *dcube
)
{
  if (dcube == NULL) {
    return;
  }

  if (dcube->distrib_bar  != NULL) free(dcube->distrib_bar);
  if (dcube->distrib_quad != NULL) free(dcube->distrib_quad);
  if (dcube->distrib_hexa != NULL) free(dcube->distrib_hexa);

  if (dcube->owner == PDM_OWNERSHIP_KEEP) {
    /* Si l'utilisateur fait le get il doit liberer le dmesh_nodal */
    PDM_DMesh_nodal_free(dcube->dmesh_nodal);
  }

  free (dcube);
}


void PDM_dcube_nodal_gen_ordering_set
(
 PDM_dcube_nodal_t *dcube,
 const char        *ordering
)
{
  dcube->ordering = (char *) ordering;
}


PDM_dmesh_nodal_t *
PDM_dcube_nodal_gen_build
(
 PDM_dcube_nodal_t *dcube
)
{
  double t1 = PDM_MPI_Wtime();

  int n_rank;
  int i_rank;

  PDM_MPI_Comm_size(dcube->comm, &n_rank);
  PDM_MPI_Comm_rank(dcube->comm, &i_rank);

  int dim = 3;
  if(dcube->t_elt == PDM_MESH_NODAL_TRIA3 || dcube->t_elt == PDM_MESH_NODAL_QUAD4 ) {
    dim = 2;
  }
  if(dcube->t_elt == PDM_MESH_NODAL_POINT || dcube->t_elt == PDM_MESH_NODAL_BAR2){
    PDM_error(__FILE__, __LINE__, 0, "Invalid t_elt for PDM_dcube_nodal_gen_init\n");
  }

  PDM_g_num_t n_vtx_x = dcube->order*dcube->nx + 1;
  PDM_g_num_t n_vtx_y = dcube->order*dcube->ny + 1;
  PDM_g_num_t n_vtx_z = dcube->order*dcube->nz + 1;


  PDM_g_num_t gn_vtx  = 0;
  PDM_g_num_t gn_bar  = 0;
  PDM_g_num_t gn_quad = 0;
  PDM_g_num_t gn_hexa = 0;
  PDM_g_num_t gn_elt  = 0;

  if (dim == 2) {
    gn_vtx  = n_vtx_x * n_vtx_y;
    gn_bar  = 2*(dcube->nx + dcube->ny);
    gn_quad = dcube->nx * dcube->ny;
    gn_elt  = gn_quad;
  } else {
    gn_vtx  = n_vtx_x * n_vtx_y * n_vtx_z;
    gn_bar  = 4*(dcube->nx + dcube->ny + dcube->nz);
    gn_quad = 2*(dcube->nx*dcube->ny + dcube->ny*dcube->nz + dcube->nz*dcube->nx);
    gn_hexa = dcube->nx * dcube->ny * dcube->nz;
    gn_elt  = gn_hexa;
  }

  //log_trace("gn_vtx = "PDM_FMT_G_NUM"\n", gn_vtx);


  /*
   * Create the dmesh_nodal that hold the resulting mesh
   */
  PDM_g_num_t gn_cell_abs = _get_n_cell_abs(gn_elt, dcube->t_elt);
  dcube->dmesh_nodal = PDM_DMesh_nodal_create(dcube->comm,
                                              dim,
                                              gn_vtx,
                                              gn_cell_abs,  /* n_cell */
                                              -1,           /* n_face */
                                              -1);          /* n_edge */

  PDM_g_num_t *distrib_vtx = PDM_compute_uniform_entity_distribution (dcube->comm, gn_vtx);
  int dn_vtx = (int) (distrib_vtx[i_rank+1] - distrib_vtx[i_rank]);

  double *dvtx_coord = malloc(sizeof(double) * dn_vtx * 3);
  PDM_DMesh_nodal_coord_set(dcube->dmesh_nodal,
                            dn_vtx,
                            dvtx_coord,
                            PDM_OWNERSHIP_KEEP); /* Le responsable de la mémoire est le dmesh_nodal */

  /*
   * Generate vertices
   */
  double step_x = dcube->length / (double) (n_vtx_x - 1);
  double step_y = dcube->length / (double) (n_vtx_y - 1);
  double step_z = dcube->length / (double) (n_vtx_z - 1);

  if (dim == 2) {
    for (int i_vtx = 0; i_vtx < dn_vtx; ++i_vtx) {

      PDM_g_num_t g_vtx = distrib_vtx[i_rank] + i_vtx;

      PDM_g_num_t indi = g_vtx % n_vtx_x;
      PDM_g_num_t indj = g_vtx / n_vtx_x;

      dvtx_coord[3 * i_vtx    ] = indi * step_x + dcube->zero_x;
      dvtx_coord[3 * i_vtx + 1] = indj * step_y + dcube->zero_y;
      dvtx_coord[3 * i_vtx + 2] = dcube->zero_z;
    }
  }
  else {
    for (int i_vtx = 0; i_vtx < dn_vtx; ++i_vtx) {

      PDM_g_num_t g_vtx = distrib_vtx[i_rank] + i_vtx;

      PDM_g_num_t indi = g_vtx % n_vtx_x;
      PDM_g_num_t indj = ((g_vtx - indi) / n_vtx_x) % n_vtx_y;
      PDM_g_num_t indk = g_vtx / (n_vtx_x * n_vtx_y);

      dvtx_coord[3 * i_vtx    ] = indi * step_x + dcube->zero_x;
      dvtx_coord[3 * i_vtx + 1] = indj * step_y + dcube->zero_y;
      dvtx_coord[3 * i_vtx + 2] = indk * step_z + dcube->zero_z;
    }
  }
  free (distrib_vtx);


  dcube->distrib_bar = PDM_compute_uniform_entity_distribution(dcube->comm, gn_bar);
  dcube->dn_bar = (int) (dcube->distrib_bar[i_rank+1] - dcube->distrib_bar[i_rank]);

  dcube->distrib_quad = PDM_compute_uniform_entity_distribution(dcube->comm, gn_quad);
  dcube->dn_quad = (int) (dcube->distrib_quad[i_rank+1] - dcube->distrib_quad[i_rank]);
  if (dim == 3) {
    dcube->distrib_hexa = PDM_compute_uniform_entity_distribution(dcube->comm, gn_hexa);
    dcube->dn_hexa = (int) (dcube->distrib_hexa[i_rank+1] - dcube->distrib_hexa[i_rank]);
  } else {
    dcube->distrib_hexa = NULL;
    dcube->dn_hexa = 0;
  }


  switch (dcube->t_elt) {
  case PDM_MESH_NODAL_TRIA3:
    {
      _generate_tria_surf (dcube, dcube->dmesh_nodal);
    }
    break;

  case PDM_MESH_NODAL_QUAD4:
    {
      _generate_quad_surf (dcube, dcube->dmesh_nodal);
    }
    break;

  case PDM_MESH_NODAL_TETRA4:
    {
      _generate_tetra_vol (dcube, dcube->dmesh_nodal);
      _generate_tetra_surf(dcube, dcube->dmesh_nodal);
    }
    break;

  case PDM_MESH_NODAL_PYRAMID5:
    {
      _generate_pyramid_vol (dcube, dcube->dmesh_nodal);
      _generate_pyramid_surf(dcube, dcube->dmesh_nodal);
    }
    break;

  case PDM_MESH_NODAL_PRISM6:
    {
      _generate_prism_vol (dcube, dcube->dmesh_nodal);
      _generate_prism_surf(dcube, dcube->dmesh_nodal);
    }
    break;

  case PDM_MESH_NODAL_HEXA8:
    {
      _generate_hexa_vol (dcube, dcube->dmesh_nodal);
      _generate_hexa_surf(dcube, dcube->dmesh_nodal);
    }
    break;

  default :
    PDM_error(__FILE__, __LINE__, 0, "Unknown element type\n");
    break;

  }

  _generate_corners(dcube, dcube->dmesh_nodal);

  _generate_ridges (dcube, dcube->dmesh_nodal);


  if (dcube->ordering != NULL) {
    PDM_dmesh_nodal_reorder (dcube->dmesh_nodal,
                             dcube->ordering,
                             dcube->order);
  }



  double t2 = PDM_MPI_Wtime();

  double delta_t = t2 - t1;
  double delta_max;
  double delta_min;

  PDM_MPI_Allreduce (&delta_t,
                     &delta_max,
                     1,
                     PDM_MPI_DOUBLE,
                     PDM_MPI_MAX,
                     dcube->comm);

  PDM_MPI_Allreduce (&delta_t,
                     &delta_min,
                     1,
                     PDM_MPI_DOUBLE,
                     PDM_MPI_MIN,
                     dcube->comm);

  if(i_rank == 0) {
    printf("[%i] PDM_dcube_nodal : duration min/max -> %12.5e %12.5e \n", n_rank, delta_min, delta_max);
  }


  return dcube->dmesh_nodal;
}
