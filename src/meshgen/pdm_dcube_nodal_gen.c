#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "pdm_part.h"
#include "pdm_dcube_nodal_gen.h"
#include "pdm_dcube_nodal_gen_priv.h"
#include "pdm_mpi.h"
#include "pdm_distrib.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_priv.h"


/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Global variable
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*
 * Generate IMIN plane into QUAD
 */
static
void
_decompose_into_quad_i
(
 PDM_dcube_nodal_t* dcube,
 PDM_dmesh_nodal_t* dmesh_nodal,
 int                plane
)
{
  int i_rank;
  PDM_MPI_Comm_rank(dcube->comm, &i_rank);

  PDM_g_num_t n_vtx_seg  = dcube->n_vtx_seg;
  int n_vtx_per_elmt_lim = 4;

  PDM_g_num_t* delmt_vtx = (PDM_g_num_t *) malloc( ( n_vtx_per_elmt_lim * dcube->dn_quad_seq_lim ) * sizeof(PDM_g_num_t));
  for(int i_quad = 0; i_quad < dcube->dn_quad_seq_lim; ++i_quad) { // QUAD

    /* Pour le remplissage du elmt_group */
    PDM_g_num_t g_quad_lim = dcube->distrib_quad_seg_lim[i_rank] + i_quad;

    PDM_g_num_t ipl_k = ( g_quad_lim                                    ) / ( dcube->n_g_hexa_cell_seg );
    PDM_g_num_t ipl_j = ( g_quad_lim - ipl_k * dcube->n_g_hexa_cell_seg );
    PDM_g_num_t ipl_i = plane; // Donc maillage en J,K

    delmt_vtx[n_vtx_per_elmt_lim * i_quad    ] = (ipl_i  ) + (ipl_j    ) * n_vtx_seg + ( ipl_k     ) * n_vtx_seg * n_vtx_seg + 1;
    delmt_vtx[n_vtx_per_elmt_lim * i_quad + 1] = (ipl_i  ) + (ipl_j    ) * n_vtx_seg + ( ipl_k + 1 ) * n_vtx_seg * n_vtx_seg + 1;
    delmt_vtx[n_vtx_per_elmt_lim * i_quad + 2] = (ipl_i  ) + (ipl_j + 1) * n_vtx_seg + ( ipl_k + 1 ) * n_vtx_seg * n_vtx_seg + 1;
    delmt_vtx[n_vtx_per_elmt_lim * i_quad + 3] = (ipl_i  ) + (ipl_j + 1) * n_vtx_seg + ( ipl_k     ) * n_vtx_seg * n_vtx_seg + 1;

    printf("IMIN :: ipl_j = "PDM_FMT_G_NUM", ipl_k = "PDM_FMT_G_NUM" \n", ipl_j, ipl_k);

  }

  int id_quad = PDM_DMesh_nodal_section_add(dmesh_nodal, PDM_MESH_NODAL_QUAD4);

  PDM_DMesh_nodal_section_std_set(dmesh_nodal,
                                  id_quad,
                                  dcube->dn_quad_seq_lim,
                                  delmt_vtx,
                                  dcube->owner_for_dmesh_nodal);
}

/*
 * Generate jmin plane into QUAD
 */
static
void
_decompose_into_quad_j
(
 PDM_dcube_nodal_t* dcube,
 PDM_dmesh_nodal_t* dmesh_nodal,
 int                plane
)
{
  int i_rank;
  PDM_MPI_Comm_rank(dcube->comm, &i_rank);

  PDM_g_num_t n_vtx_seg  = dcube->n_vtx_seg;
  int n_vtx_per_elmt_lim = 4;

  PDM_g_num_t* delmt_vtx = (PDM_g_num_t *) malloc( ( n_vtx_per_elmt_lim * dcube->dn_quad_seq_lim ) * sizeof(PDM_g_num_t));
  for(int i_quad = 0; i_quad < dcube->dn_quad_seq_lim; ++i_quad) { // QUAD

    PDM_g_num_t g_quad_lim = dcube->distrib_quad_seg_lim[i_rank] + i_quad;

    PDM_g_num_t ipl_k = ( g_quad_lim                                    ) / ( dcube->n_g_hexa_cell_seg );
    PDM_g_num_t ipl_j = plane; // Donc maillage en I,K
    PDM_g_num_t ipl_i = ( g_quad_lim - ipl_k * dcube->n_g_hexa_cell_seg );

    delmt_vtx[n_vtx_per_elmt_lim * i_quad    ] = (ipl_i  ) + (ipl_j    ) * n_vtx_seg + ( ipl_k     ) * n_vtx_seg * n_vtx_seg + 1;
    delmt_vtx[n_vtx_per_elmt_lim * i_quad + 1] = (ipl_i+1) + (ipl_j    ) * n_vtx_seg + ( ipl_k     ) * n_vtx_seg * n_vtx_seg + 1;
    delmt_vtx[n_vtx_per_elmt_lim * i_quad + 2] = (ipl_i+1) + (ipl_j    ) * n_vtx_seg + ( ipl_k + 1 ) * n_vtx_seg * n_vtx_seg + 1;
    delmt_vtx[n_vtx_per_elmt_lim * i_quad + 3] = (ipl_i  ) + (ipl_j    ) * n_vtx_seg + ( ipl_k + 1 ) * n_vtx_seg * n_vtx_seg + 1;

    printf("JMIN :: ipl_j = "PDM_FMT_G_NUM", ipl_k = "PDM_FMT_G_NUM" \n", ipl_j, ipl_k);

  }

  int id_quad = PDM_DMesh_nodal_section_add(dmesh_nodal, PDM_MESH_NODAL_QUAD4);

  PDM_DMesh_nodal_section_std_set(dmesh_nodal,
                                  id_quad,
                                  dcube->dn_quad_seq_lim,
                                  delmt_vtx,
                                  dcube->owner_for_dmesh_nodal);
}

/*
 * Generate jmin plane into QUAD
 */
static
void
_decompose_into_quad_k
(
 PDM_dcube_nodal_t* dcube,
 PDM_dmesh_nodal_t* dmesh_nodal,
 int                plane
)
{
  int i_rank;
  PDM_MPI_Comm_rank(dcube->comm, &i_rank);

  PDM_g_num_t n_vtx_seg  = dcube->n_vtx_seg;
  int n_vtx_per_elmt_lim = 4;

  PDM_g_num_t* delmt_vtx = (PDM_g_num_t *) malloc( ( n_vtx_per_elmt_lim * dcube->dn_quad_seq_lim ) * sizeof(PDM_g_num_t));
  for(int i_quad = 0; i_quad < dcube->dn_quad_seq_lim; ++i_quad) { // QUAD

    /* Pour le remplissage du elmt_group */
    PDM_g_num_t g_quad_lim = dcube->distrib_quad_seg_lim[i_rank] + i_quad;

    PDM_g_num_t ipl_k = plane; // Donc maillage en I,J
    PDM_g_num_t ipl_j = ( g_quad_lim                                    ) / ( dcube->n_g_hexa_cell_seg );
    PDM_g_num_t ipl_i = ( g_quad_lim - ipl_j * dcube->n_g_hexa_cell_seg );

    delmt_vtx[n_vtx_per_elmt_lim * i_quad    ] = (ipl_i  ) + (ipl_j    ) * n_vtx_seg + ( ipl_k ) * n_vtx_seg * n_vtx_seg + 1;
    delmt_vtx[n_vtx_per_elmt_lim * i_quad + 1] = (ipl_i  ) + (ipl_j + 1) * n_vtx_seg + ( ipl_k ) * n_vtx_seg * n_vtx_seg + 1;
    delmt_vtx[n_vtx_per_elmt_lim * i_quad + 2] = (ipl_i+1) + (ipl_j + 1) * n_vtx_seg + ( ipl_k ) * n_vtx_seg * n_vtx_seg + 1;
    delmt_vtx[n_vtx_per_elmt_lim * i_quad + 3] = (ipl_i+1) + (ipl_j    ) * n_vtx_seg + ( ipl_k ) * n_vtx_seg * n_vtx_seg + 1;

    printf("kMIN :: ipl_j = "PDM_FMT_G_NUM", ipl_k = "PDM_FMT_G_NUM" \n", ipl_j, ipl_k);

  }

  int id_quad = PDM_DMesh_nodal_section_add(dmesh_nodal, PDM_MESH_NODAL_QUAD4);

  PDM_DMesh_nodal_section_std_set(dmesh_nodal,
                                  id_quad,
                                  dcube->dn_quad_seq_lim,
                                  delmt_vtx,
                                  dcube->owner_for_dmesh_nodal);
}

/*
 * Generate jmin plane into TRI
 */
static
void
_decompose_into_tri_j
(
 PDM_dcube_nodal_t* dcube,
 PDM_dmesh_nodal_t* dmesh_nodal,
 int                plane
)
{
  int i_rank;
  PDM_MPI_Comm_rank(dcube->comm, &i_rank);

  PDM_g_num_t n_vtx_seg  = dcube->n_vtx_seg;
  int n_vtx_per_elmt_lim = 3;

  int dn_tri_seq_lim = 2 * dcube->dn_quad_seq_lim;

  PDM_g_num_t* delmt_vtx = (PDM_g_num_t *) malloc( ( n_vtx_per_elmt_lim * dn_tri_seq_lim ) * sizeof(PDM_g_num_t));
  for(int i_quad = 0; i_quad < dcube->dn_quad_seq_lim; ++i_quad) { // QUAD

    PDM_g_num_t g_quad_lim = dcube->distrib_quad_seg_lim[i_rank] + i_quad;

    PDM_g_num_t ipl_k = ( g_quad_lim                                    ) / ( dcube->n_g_hexa_cell_seg );
    PDM_g_num_t ipl_j = plane; // Donc maillage en I,K
    PDM_g_num_t ipl_i = ( g_quad_lim - ipl_k * dcube->n_g_hexa_cell_seg );

    PDM_g_num_t n = ipl_i+ipl_k;

    if( n % 2 == 0) {

      delmt_vtx[2 * n_vtx_per_elmt_lim * i_quad    ] = (ipl_i  ) + (ipl_j    ) * n_vtx_seg + ( ipl_k     ) * n_vtx_seg * n_vtx_seg + 1;
      delmt_vtx[2 * n_vtx_per_elmt_lim * i_quad + 1] = (ipl_i+1) + (ipl_j    ) * n_vtx_seg + ( ipl_k     ) * n_vtx_seg * n_vtx_seg + 1;
      delmt_vtx[2 * n_vtx_per_elmt_lim * i_quad + 2] = (ipl_i  ) + (ipl_j    ) * n_vtx_seg + ( ipl_k + 1 ) * n_vtx_seg * n_vtx_seg + 1;

      delmt_vtx[2 * n_vtx_per_elmt_lim * i_quad + 3] = (ipl_i+1) + (ipl_j    ) * n_vtx_seg + ( ipl_k     ) * n_vtx_seg * n_vtx_seg + 1;
      delmt_vtx[2 * n_vtx_per_elmt_lim * i_quad + 4] = (ipl_i+1) + (ipl_j    ) * n_vtx_seg + ( ipl_k + 1 ) * n_vtx_seg * n_vtx_seg + 1;
      delmt_vtx[2 * n_vtx_per_elmt_lim * i_quad + 5] = (ipl_i  ) + (ipl_j    ) * n_vtx_seg + ( ipl_k + 1 ) * n_vtx_seg * n_vtx_seg + 1;

    } else {

      delmt_vtx[2 * n_vtx_per_elmt_lim * i_quad    ] = (ipl_i  ) + (ipl_j    ) * n_vtx_seg + ( ipl_k     ) * n_vtx_seg * n_vtx_seg + 1;
      delmt_vtx[2 * n_vtx_per_elmt_lim * i_quad + 1] = (ipl_i+1) + (ipl_j    ) * n_vtx_seg + ( ipl_k     ) * n_vtx_seg * n_vtx_seg + 1;
      delmt_vtx[2 * n_vtx_per_elmt_lim * i_quad + 2] = (ipl_i+1) + (ipl_j    ) * n_vtx_seg + ( ipl_k + 1 ) * n_vtx_seg * n_vtx_seg + 1;

      delmt_vtx[2 * n_vtx_per_elmt_lim * i_quad + 3] = (ipl_i  ) + (ipl_j    ) * n_vtx_seg + ( ipl_k     ) * n_vtx_seg * n_vtx_seg + 1;
      delmt_vtx[2 * n_vtx_per_elmt_lim * i_quad + 4] = (ipl_i+1) + (ipl_j    ) * n_vtx_seg + ( ipl_k + 1 ) * n_vtx_seg * n_vtx_seg + 1;
      delmt_vtx[2 * n_vtx_per_elmt_lim * i_quad + 5] = (ipl_i  ) + (ipl_j    ) * n_vtx_seg + ( ipl_k + 1 ) * n_vtx_seg * n_vtx_seg + 1;

    }
    // printf("JMIN :: ipl_j = "PDM_FMT_G_NUM", ipl_k = "PDM_FMT_G_NUM" \n", ipl_j, ipl_k);

  }

  int id_tri = PDM_DMesh_nodal_section_add(dmesh_nodal, PDM_MESH_NODAL_TRIA3);

  PDM_DMesh_nodal_section_std_set(dmesh_nodal,
                                  id_tri,
                                  dn_tri_seq_lim,
                                  delmt_vtx,
                                  dcube->owner_for_dmesh_nodal);
}

/*
 * Generate jmin plane into TRI
 */
static
void
_decompose_into_tri_i
(
 PDM_dcube_nodal_t* dcube,
 PDM_dmesh_nodal_t* dmesh_nodal,
 int                plane
)
{
  int i_rank;
  PDM_MPI_Comm_rank(dcube->comm, &i_rank);

  PDM_g_num_t n_vtx_seg  = dcube->n_vtx_seg;
  int n_vtx_per_elmt_lim = 3;

  int dn_tri_seq_lim = 2 * dcube->dn_quad_seq_lim;

  PDM_g_num_t* delmt_vtx = (PDM_g_num_t *) malloc( ( n_vtx_per_elmt_lim * dn_tri_seq_lim ) * sizeof(PDM_g_num_t));
  for(int i_quad = 0; i_quad < dcube->dn_quad_seq_lim; ++i_quad) { // QUAD

    PDM_g_num_t g_quad_lim = dcube->distrib_quad_seg_lim[i_rank] + i_quad;

    PDM_g_num_t ipl_k = ( g_quad_lim                                    ) / ( dcube->n_g_hexa_cell_seg );
    PDM_g_num_t ipl_j = ( g_quad_lim - ipl_k * dcube->n_g_hexa_cell_seg );
    PDM_g_num_t ipl_i = plane; // Donc maillage en J,K

    PDM_g_num_t n = ipl_j+ipl_k;

    if( n % 2 == 0) {

      delmt_vtx[2 * n_vtx_per_elmt_lim * i_quad    ] = (ipl_i  ) + (ipl_j    ) * n_vtx_seg + ( ipl_k     ) * n_vtx_seg * n_vtx_seg + 1;
      delmt_vtx[2 * n_vtx_per_elmt_lim * i_quad + 1] = (ipl_i  ) + (ipl_j +1 ) * n_vtx_seg + ( ipl_k     ) * n_vtx_seg * n_vtx_seg + 1;
      delmt_vtx[2 * n_vtx_per_elmt_lim * i_quad + 2] = (ipl_i  ) + (ipl_j    ) * n_vtx_seg + ( ipl_k + 1 ) * n_vtx_seg * n_vtx_seg + 1;

      delmt_vtx[2 * n_vtx_per_elmt_lim * i_quad + 3] = (ipl_i  ) + (ipl_j +1 ) * n_vtx_seg + ( ipl_k     ) * n_vtx_seg * n_vtx_seg + 1;
      delmt_vtx[2 * n_vtx_per_elmt_lim * i_quad + 4] = (ipl_i  ) + (ipl_j +1 ) * n_vtx_seg + ( ipl_k + 1 ) * n_vtx_seg * n_vtx_seg + 1;
      delmt_vtx[2 * n_vtx_per_elmt_lim * i_quad + 5] = (ipl_i  ) + (ipl_j    ) * n_vtx_seg + ( ipl_k + 1 ) * n_vtx_seg * n_vtx_seg + 1;

    } else {

      delmt_vtx[2 * n_vtx_per_elmt_lim * i_quad    ] = (ipl_i  ) + (ipl_j    ) * n_vtx_seg + ( ipl_k     ) * n_vtx_seg * n_vtx_seg + 1;
      delmt_vtx[2 * n_vtx_per_elmt_lim * i_quad + 1] = (ipl_i  ) + (ipl_j +1 ) * n_vtx_seg + ( ipl_k     ) * n_vtx_seg * n_vtx_seg + 1;
      delmt_vtx[2 * n_vtx_per_elmt_lim * i_quad + 2] = (ipl_i  ) + (ipl_j +1 ) * n_vtx_seg + ( ipl_k + 1 ) * n_vtx_seg * n_vtx_seg + 1;

      delmt_vtx[2 * n_vtx_per_elmt_lim * i_quad + 3] = (ipl_i  ) + (ipl_j    ) * n_vtx_seg + ( ipl_k     ) * n_vtx_seg * n_vtx_seg + 1;
      delmt_vtx[2 * n_vtx_per_elmt_lim * i_quad + 4] = (ipl_i  ) + (ipl_j +1 ) * n_vtx_seg + ( ipl_k + 1 ) * n_vtx_seg * n_vtx_seg + 1;
      delmt_vtx[2 * n_vtx_per_elmt_lim * i_quad + 5] = (ipl_i  ) + (ipl_j    ) * n_vtx_seg + ( ipl_k + 1 ) * n_vtx_seg * n_vtx_seg + 1;


    }
    // printf("JMIN :: ipl_j = "PDM_FMT_G_NUM", ipl_k = "PDM_FMT_G_NUM" \n", ipl_j, ipl_k);

  }

  int id_tri = PDM_DMesh_nodal_section_add(dmesh_nodal, PDM_MESH_NODAL_TRIA3);

  PDM_DMesh_nodal_section_std_set(dmesh_nodal,
                                  id_tri,
                                  dn_tri_seq_lim,
                                  delmt_vtx,
                                  dcube->owner_for_dmesh_nodal);
}

/*
 * Generate jmin plane into TRI
 */
static
void
_decompose_into_tri_k
(
 PDM_dcube_nodal_t* dcube,
 PDM_dmesh_nodal_t* dmesh_nodal,
 int                plane
)
{
  int i_rank;
  PDM_MPI_Comm_rank(dcube->comm, &i_rank);

  PDM_g_num_t n_vtx_seg  = dcube->n_vtx_seg;
  int n_vtx_per_elmt_lim = 3;

  int dn_tri_seq_lim = 2 * dcube->dn_quad_seq_lim;

  PDM_g_num_t* delmt_vtx = (PDM_g_num_t *) malloc( ( n_vtx_per_elmt_lim * dn_tri_seq_lim ) * sizeof(PDM_g_num_t));
  for(int i_quad = 0; i_quad < dcube->dn_quad_seq_lim; ++i_quad) { // QUAD

    PDM_g_num_t g_quad_lim = dcube->distrib_quad_seg_lim[i_rank] + i_quad;

    PDM_g_num_t ipl_k = plane; // Donc maillage en I,J
    PDM_g_num_t ipl_j = ( g_quad_lim                                    ) / ( dcube->n_g_hexa_cell_seg );
    PDM_g_num_t ipl_i = ( g_quad_lim - ipl_j * dcube->n_g_hexa_cell_seg );

    PDM_g_num_t n = ipl_i+ipl_j;

    if( n % 2 == 0) {

      delmt_vtx[2 * n_vtx_per_elmt_lim * i_quad    ] = (ipl_i     ) + (ipl_j    ) * n_vtx_seg + ( ipl_k ) * n_vtx_seg * n_vtx_seg + 1;
      delmt_vtx[2 * n_vtx_per_elmt_lim * i_quad + 1] = (ipl_i     ) + (ipl_j +1 ) * n_vtx_seg + ( ipl_k ) * n_vtx_seg * n_vtx_seg + 1;
      delmt_vtx[2 * n_vtx_per_elmt_lim * i_quad + 2] = (ipl_i + 1 ) + (ipl_j    ) * n_vtx_seg + ( ipl_k ) * n_vtx_seg * n_vtx_seg + 1;

      delmt_vtx[2 * n_vtx_per_elmt_lim * i_quad + 3] = (ipl_i     ) + (ipl_j +1 ) * n_vtx_seg + ( ipl_k ) * n_vtx_seg * n_vtx_seg + 1;
      delmt_vtx[2 * n_vtx_per_elmt_lim * i_quad + 4] = (ipl_i + 1 ) + (ipl_j +1 ) * n_vtx_seg + ( ipl_k ) * n_vtx_seg * n_vtx_seg + 1;
      delmt_vtx[2 * n_vtx_per_elmt_lim * i_quad + 5] = (ipl_i + 1 ) + (ipl_j    ) * n_vtx_seg + ( ipl_k ) * n_vtx_seg * n_vtx_seg + 1;

    } else {

      delmt_vtx[2 * n_vtx_per_elmt_lim * i_quad    ] = (ipl_i     ) + (ipl_j    ) * n_vtx_seg + ( ipl_k ) * n_vtx_seg * n_vtx_seg + 1;
      delmt_vtx[2 * n_vtx_per_elmt_lim * i_quad + 1] = (ipl_i     ) + (ipl_j +1 ) * n_vtx_seg + ( ipl_k ) * n_vtx_seg * n_vtx_seg + 1;
      delmt_vtx[2 * n_vtx_per_elmt_lim * i_quad + 2] = (ipl_i + 1 ) + (ipl_j +1 ) * n_vtx_seg + ( ipl_k ) * n_vtx_seg * n_vtx_seg + 1;

      delmt_vtx[2 * n_vtx_per_elmt_lim * i_quad + 3] = (ipl_i     ) + (ipl_j    ) * n_vtx_seg + ( ipl_k ) * n_vtx_seg * n_vtx_seg + 1;
      delmt_vtx[2 * n_vtx_per_elmt_lim * i_quad + 4] = (ipl_i + 1 ) + (ipl_j +1 ) * n_vtx_seg + ( ipl_k ) * n_vtx_seg * n_vtx_seg + 1;
      delmt_vtx[2 * n_vtx_per_elmt_lim * i_quad + 5] = (ipl_i + 1 ) + (ipl_j    ) * n_vtx_seg + ( ipl_k ) * n_vtx_seg * n_vtx_seg + 1;


    }
    printf("JMIN :: ipl_j = "PDM_FMT_G_NUM", ipl_k = "PDM_FMT_G_NUM" \n", ipl_j, ipl_k);

  }

  int id_tri = PDM_DMesh_nodal_section_add(dmesh_nodal, PDM_MESH_NODAL_TRIA3);

  PDM_DMesh_nodal_section_std_set(dmesh_nodal,
                                  id_tri,
                                  dn_tri_seq_lim,
                                  delmt_vtx,
                                  dcube->owner_for_dmesh_nodal);
}


static
PDM_g_num_t
_get_n_cell_abs
(
PDM_g_num_t          n_hexa,
PDM_Mesh_nodal_elt_t t_elmt
)
{
  switch (t_elmt) {
    case PDM_MESH_NODAL_TETRA4    :
    {
      // Each hexa in split in 5 tetra and boundary
      return n_hexa * 5;
    }
    break;

    case PDM_MESH_NODAL_PRISM6    :
    {
      return n_hexa * 2;
    }
    break;

    case PDM_MESH_NODAL_HEXA8    :
    {
      return n_hexa;
    }
    break;

    default :
      PDM_error(__FILE__, __LINE__, 0, "Unknown element type\n");
      break;

  }
  return -1;
}

static
void
_generate_tetra_from_hexa
(
 PDM_dcube_nodal_t* dcube,
 PDM_dmesh_nodal_t* dmesh_nodal
)
{
  int i_rank;
  PDM_MPI_Comm_rank(dcube->comm, &i_rank);

  int n_vtx_per_elmt     = 4;

  PDM_g_num_t ngh       = dcube->n_g_hexa_cell_seg;
  PDM_g_num_t ngh2x     = ngh * ngh;
  PDM_g_num_t n_vtx_seg = dcube->n_vtx_seg;

  int dn_tetra_cell = 5 * dcube->dn_hexa_cell;

  /* Setup volumic part */
  PDM_g_num_t* delmt_vtx = (PDM_g_num_t *) malloc( ( n_vtx_per_elmt * dn_tetra_cell ) * sizeof(PDM_g_num_t));

  for(int i_cell = 0; i_cell < dcube->dn_hexa_cell; ++i_cell) { // HEXA

    PDM_g_num_t g_cell = dcube->distrib_hexa[i_rank] + i_cell;

    PDM_g_num_t indk = ( g_cell                             ) / ( ngh2x );
    PDM_g_num_t indj = ( g_cell - indk * ngh2x              ) / ( ngh   );
    PDM_g_num_t indi = ( g_cell - indk * ngh2x - indj * ngh );

    PDM_g_num_t n = indi+indj+indk;

    PDM_g_num_t ind1 = (indi     ) + (indj     ) * n_vtx_seg + ( indk     ) * n_vtx_seg * n_vtx_seg + 1; // A (  i,  j,k  )
    PDM_g_num_t ind2 = (indi + 1 ) + (indj     ) * n_vtx_seg + ( indk     ) * n_vtx_seg * n_vtx_seg + 1; // B (i+1,  j,k  )
    PDM_g_num_t ind3 = (indi + 1 ) + (indj + 1 ) * n_vtx_seg + ( indk     ) * n_vtx_seg * n_vtx_seg + 1; // C (i+1,j+1,k  )
    PDM_g_num_t ind4 = (indi     ) + (indj + 1 ) * n_vtx_seg + ( indk     ) * n_vtx_seg * n_vtx_seg + 1; // D (  i,j+1,k  )
    PDM_g_num_t ind5 = (indi     ) + (indj     ) * n_vtx_seg + ( indk + 1 ) * n_vtx_seg * n_vtx_seg + 1; // E (  i,  j,k+1)
    PDM_g_num_t ind6 = (indi + 1 ) + (indj     ) * n_vtx_seg + ( indk + 1 ) * n_vtx_seg * n_vtx_seg + 1; // F (i+1,  j,k+1)
    PDM_g_num_t ind7 = (indi + 1 ) + (indj + 1 ) * n_vtx_seg + ( indk + 1 ) * n_vtx_seg * n_vtx_seg + 1; // G (i+1,j+1,k+1)
    PDM_g_num_t ind8 = (indi     ) + (indj + 1 ) * n_vtx_seg + ( indk + 1 ) * n_vtx_seg * n_vtx_seg + 1; // H (  i,j+1,k+1)

    if( n % 2 == 0) {

      delmt_vtx[5 * n_vtx_per_elmt * i_cell     ] = ind1;
      delmt_vtx[5 * n_vtx_per_elmt * i_cell + 1 ] = ind2;
      delmt_vtx[5 * n_vtx_per_elmt * i_cell + 2 ] = ind4;
      delmt_vtx[5 * n_vtx_per_elmt * i_cell + 3 ] = ind5;

      delmt_vtx[5 * n_vtx_per_elmt * i_cell + 4 ] = ind2;
      delmt_vtx[5 * n_vtx_per_elmt * i_cell + 5 ] = ind3;
      delmt_vtx[5 * n_vtx_per_elmt * i_cell + 6 ] = ind4;
      delmt_vtx[5 * n_vtx_per_elmt * i_cell + 7 ] = ind7;

      delmt_vtx[5 * n_vtx_per_elmt * i_cell + 8 ] = ind4;
      delmt_vtx[5 * n_vtx_per_elmt * i_cell + 9 ] = ind5;
      delmt_vtx[5 * n_vtx_per_elmt * i_cell + 10] = ind7;
      delmt_vtx[5 * n_vtx_per_elmt * i_cell + 11] = ind8;

      delmt_vtx[5 * n_vtx_per_elmt * i_cell + 12] = ind2;
      delmt_vtx[5 * n_vtx_per_elmt * i_cell + 13] = ind5;
      delmt_vtx[5 * n_vtx_per_elmt * i_cell + 14] = ind6;
      delmt_vtx[5 * n_vtx_per_elmt * i_cell + 15] = ind7;

      delmt_vtx[5 * n_vtx_per_elmt * i_cell + 16] = ind2;
      delmt_vtx[5 * n_vtx_per_elmt * i_cell + 17] = ind4;
      delmt_vtx[5 * n_vtx_per_elmt * i_cell + 18] = ind5;
      delmt_vtx[5 * n_vtx_per_elmt * i_cell + 19] = ind7;

    } else {

      delmt_vtx[5 * n_vtx_per_elmt * i_cell     ] = ind1;
      delmt_vtx[5 * n_vtx_per_elmt * i_cell + 1 ] = ind3;
      delmt_vtx[5 * n_vtx_per_elmt * i_cell + 2 ] = ind4;
      delmt_vtx[5 * n_vtx_per_elmt * i_cell + 3 ] = ind8;

      delmt_vtx[5 * n_vtx_per_elmt * i_cell + 4 ] = ind1;
      delmt_vtx[5 * n_vtx_per_elmt * i_cell + 5 ] = ind2;
      delmt_vtx[5 * n_vtx_per_elmt * i_cell + 6 ] = ind3;
      delmt_vtx[5 * n_vtx_per_elmt * i_cell + 7 ] = ind6;

      delmt_vtx[5 * n_vtx_per_elmt * i_cell + 8 ] = ind3;
      delmt_vtx[5 * n_vtx_per_elmt * i_cell + 9 ] = ind6;
      delmt_vtx[5 * n_vtx_per_elmt * i_cell + 10] = ind7;
      delmt_vtx[5 * n_vtx_per_elmt * i_cell + 11] = ind8;

      delmt_vtx[5 * n_vtx_per_elmt * i_cell + 12] = ind6;
      delmt_vtx[5 * n_vtx_per_elmt * i_cell + 13] = ind5;
      delmt_vtx[5 * n_vtx_per_elmt * i_cell + 14] = ind8;
      delmt_vtx[5 * n_vtx_per_elmt * i_cell + 15] = ind1;

      delmt_vtx[5 * n_vtx_per_elmt * i_cell + 16] = ind6;
      delmt_vtx[5 * n_vtx_per_elmt * i_cell + 17] = ind3;
      delmt_vtx[5 * n_vtx_per_elmt * i_cell + 18] = ind1;
      delmt_vtx[5 * n_vtx_per_elmt * i_cell + 19] = ind8;

    }
    printf(" i_cell = %i -> %i \n", i_cell, 2*i_cell);

  }

  int id_tetra = PDM_DMesh_nodal_section_add(dmesh_nodal, PDM_MESH_NODAL_TETRA4);
  PDM_DMesh_nodal_section_std_set(dmesh_nodal,
                                  id_tetra,
                                  dn_tetra_cell,
                                  delmt_vtx,
                                  dcube->owner_for_dmesh_nodal);

  _decompose_into_tri_i(dcube, dmesh_nodal, 0);
  _decompose_into_tri_i(dcube, dmesh_nodal, dcube->n_g_hexa_cell_seg);
  _decompose_into_tri_j(dcube, dmesh_nodal, 0);
  _decompose_into_tri_j(dcube, dmesh_nodal, dcube->n_g_hexa_cell_seg);
  _decompose_into_tri_k(dcube, dmesh_nodal, 0);
  _decompose_into_tri_k(dcube, dmesh_nodal, dcube->n_g_hexa_cell_seg);

}


static
void
_generate_prism_from_hexa
(
 PDM_dcube_nodal_t* dcube,
 PDM_dmesh_nodal_t* dmesh_nodal
)
{
  int i_rank;
  PDM_MPI_Comm_rank(dcube->comm, &i_rank);

  int n_vtx_per_elmt     = 6;

  PDM_g_num_t ngh       = dcube->n_g_hexa_cell_seg;
  PDM_g_num_t ngh2x     = ngh * ngh;
  PDM_g_num_t n_vtx_seg = dcube->n_vtx_seg;

  int dn_prism_cell = 2 * dcube->dn_hexa_cell;

  /* Setup volumic part */
  PDM_g_num_t* delmt_vtx = (PDM_g_num_t *) malloc( ( n_vtx_per_elmt * dn_prism_cell ) * sizeof(PDM_g_num_t));

  for(int i_cell = 0; i_cell < dcube->dn_hexa_cell; ++i_cell) {

    PDM_g_num_t g_cell = dcube->distrib_hexa[i_rank] + i_cell;

    PDM_g_num_t indk = ( g_cell                             ) / ( ngh2x );
    PDM_g_num_t indj = ( g_cell - indk * ngh2x              ) / ( ngh   );
    PDM_g_num_t indi = ( g_cell - indk * ngh2x - indj * ngh );

    delmt_vtx[2 * n_vtx_per_elmt * i_cell     ] = (indi  ) + (indj    ) * n_vtx_seg + ( indk   ) * n_vtx_seg * n_vtx_seg + 1;
    delmt_vtx[2 * n_vtx_per_elmt * i_cell + 1 ] = (indi+1) + (indj+1  ) * n_vtx_seg + ( indk   ) * n_vtx_seg * n_vtx_seg + 1;
    delmt_vtx[2 * n_vtx_per_elmt * i_cell + 2 ] = (indi  ) + (indj+1  ) * n_vtx_seg + ( indk   ) * n_vtx_seg * n_vtx_seg + 1;
    delmt_vtx[2 * n_vtx_per_elmt * i_cell + 3 ] = (indi  ) + (indj    ) * n_vtx_seg + ( indk+1 ) * n_vtx_seg * n_vtx_seg + 1;
    delmt_vtx[2 * n_vtx_per_elmt * i_cell + 4 ] = (indi+1) + (indj+1  ) * n_vtx_seg + ( indk+1 ) * n_vtx_seg * n_vtx_seg + 1;
    delmt_vtx[2 * n_vtx_per_elmt * i_cell + 5 ] = (indi  ) + (indj+1  ) * n_vtx_seg + ( indk+1 ) * n_vtx_seg * n_vtx_seg + 1;

    delmt_vtx[2 * n_vtx_per_elmt * i_cell + 6 ] = (indi  ) + (indj    ) * n_vtx_seg + ( indk   ) * n_vtx_seg * n_vtx_seg + 1;
    delmt_vtx[2 * n_vtx_per_elmt * i_cell + 7 ] = (indi+1) + (indj    ) * n_vtx_seg + ( indk   ) * n_vtx_seg * n_vtx_seg + 1;
    delmt_vtx[2 * n_vtx_per_elmt * i_cell + 8 ] = (indi+1) + (indj+1  ) * n_vtx_seg + ( indk   ) * n_vtx_seg * n_vtx_seg + 1;
    delmt_vtx[2 * n_vtx_per_elmt * i_cell + 9 ] = (indi  ) + (indj    ) * n_vtx_seg + ( indk+1 ) * n_vtx_seg * n_vtx_seg + 1;
    delmt_vtx[2 * n_vtx_per_elmt * i_cell + 10] = (indi+1) + (indj    ) * n_vtx_seg + ( indk+1 ) * n_vtx_seg * n_vtx_seg + 1;
    delmt_vtx[2 * n_vtx_per_elmt * i_cell + 11] = (indi+1) + (indj+1  ) * n_vtx_seg + ( indk+1 ) * n_vtx_seg * n_vtx_seg + 1;

    // printf(" i_cell = %i -> %i \n", i_cell, 2*i_cell);

  }

  int id_prism = PDM_DMesh_nodal_section_add(dmesh_nodal, PDM_MESH_NODAL_PRISM6);
  PDM_DMesh_nodal_section_std_set(dmesh_nodal,
                                  id_prism,
                                  dn_prism_cell,
                                  delmt_vtx,
                                  dcube->owner_for_dmesh_nodal);

  _decompose_into_quad_i(dcube, dmesh_nodal, 0);
  _decompose_into_quad_i(dcube, dmesh_nodal, dcube->n_g_hexa_cell_seg);
  _decompose_into_quad_j(dcube, dmesh_nodal, 0);
  _decompose_into_quad_j(dcube, dmesh_nodal, dcube->n_g_hexa_cell_seg);
  _decompose_into_tri_k (dcube, dmesh_nodal, 0);
  _decompose_into_tri_k (dcube, dmesh_nodal, dcube->n_g_hexa_cell_seg);

}

static
void
_generate_hexa_from_hexa
(
 PDM_dcube_nodal_t* dcube,
 PDM_dmesh_nodal_t* dmesh_nodal
)
{
  int i_rank;
  PDM_MPI_Comm_rank(dcube->comm, &i_rank);

  int n_vtx_per_elmt     = 8;

  PDM_g_num_t ngh       = dcube->n_g_hexa_cell_seg;
  PDM_g_num_t ngh2x     = ngh * ngh;
  PDM_g_num_t n_vtx_seg = dcube->n_vtx_seg;

  /* Setup volumic part */
  PDM_g_num_t* delmt_vtx = (PDM_g_num_t *) malloc( ( n_vtx_per_elmt * dcube->dn_hexa_cell ) * sizeof(PDM_g_num_t));

  for(int i_cell = 0; i_cell < dcube->dn_hexa_cell; ++i_cell) { // HEXA

    PDM_g_num_t g_cell = dcube->distrib_hexa[i_rank] + i_cell;

    PDM_g_num_t indk = ( g_cell                             ) / ( ngh2x );
    PDM_g_num_t indj = ( g_cell - indk * ngh2x              ) / ( ngh   );
    PDM_g_num_t indi = ( g_cell - indk * ngh2x - indj * ngh );

    delmt_vtx[n_vtx_per_elmt * i_cell    ] = (indi  ) + (indj    ) * n_vtx_seg + ( indk   ) * n_vtx_seg * n_vtx_seg + 1;
    delmt_vtx[n_vtx_per_elmt * i_cell + 1] = (indi+1) + (indj    ) * n_vtx_seg + ( indk   ) * n_vtx_seg * n_vtx_seg + 1;
    delmt_vtx[n_vtx_per_elmt * i_cell + 2] = (indi+1) + (indj    ) * n_vtx_seg + ( indk+1 ) * n_vtx_seg * n_vtx_seg + 1;
    delmt_vtx[n_vtx_per_elmt * i_cell + 3] = (indi  ) + (indj    ) * n_vtx_seg + ( indk+1 ) * n_vtx_seg * n_vtx_seg + 1;
    delmt_vtx[n_vtx_per_elmt * i_cell + 4] = (indi  ) + (indj+1  ) * n_vtx_seg + ( indk   ) * n_vtx_seg * n_vtx_seg + 1;
    delmt_vtx[n_vtx_per_elmt * i_cell + 5] = (indi+1) + (indj+1  ) * n_vtx_seg + ( indk   ) * n_vtx_seg * n_vtx_seg + 1;
    delmt_vtx[n_vtx_per_elmt * i_cell + 6] = (indi+1) + (indj+1  ) * n_vtx_seg + ( indk+1 ) * n_vtx_seg * n_vtx_seg + 1;
    delmt_vtx[n_vtx_per_elmt * i_cell + 7] = (indi  ) + (indj+1  ) * n_vtx_seg + ( indk+1 ) * n_vtx_seg * n_vtx_seg + 1;

  }

  int id_hexa = PDM_DMesh_nodal_section_add(dmesh_nodal, PDM_MESH_NODAL_HEXA8);
  PDM_DMesh_nodal_section_std_set(dmesh_nodal,
                                  id_hexa,
                                  dcube->dn_hexa_cell,
                                  delmt_vtx,
                                  dcube->owner_for_dmesh_nodal);


  _decompose_into_quad_i(dcube, dmesh_nodal, 0);
  _decompose_into_quad_i(dcube, dmesh_nodal, dcube->n_g_hexa_cell_seg);
  _decompose_into_quad_j(dcube, dmesh_nodal, 0);
  _decompose_into_quad_j(dcube, dmesh_nodal, dcube->n_g_hexa_cell_seg);
  _decompose_into_quad_k(dcube, dmesh_nodal, 0);
  _decompose_into_quad_k(dcube, dmesh_nodal, dcube->n_g_hexa_cell_seg);

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

PDM_dcube_nodal_t*
PDM_dcube_nodal_gen_init
(
      PDM_MPI_Comm          comm,
const PDM_g_num_t           n_vtx_seg,
const double                length,
const double                zero_x,
const double                zero_y,
const double                zero_z,
      PDM_Mesh_nodal_elt_t  t_elt,
      PDM_ownership_t       owner
)
{
  int n_rank;
  int i_rank;

  PDM_MPI_Comm_size(comm, &n_rank);
  PDM_MPI_Comm_rank(comm, &i_rank);

  PDM_dcube_nodal_t *dcube = (PDM_dcube_nodal_t *) malloc(sizeof(PDM_dcube_nodal_t));

  /*
   * Build dcube structure
   */
  dcube->comm      = comm;
  dcube->n_vtx_seg = n_vtx_seg;
  dcube->length    = length;
  dcube->zero_x    = zero_x;
  dcube->zero_y    = zero_y;
  dcube->zero_z    = zero_z;
  dcube->t_elt     = t_elt;
  dcube->owner     = owner;

  // Si l'utilisateur prends les résulats, c'est à travers dmesh_nodal
  dcube->owner_for_dmesh_nodal = owner;

  PDM_g_num_t n_vtx        = n_vtx_seg  * n_vtx_seg * n_vtx_seg;
  dcube->n_g_hexa_cell_seg = n_vtx_seg - 1;
  dcube->n_g_hexa_cell     = dcube->n_g_hexa_cell_seg * dcube->n_g_hexa_cell_seg * dcube->n_g_hexa_cell_seg;

  PDM_g_num_t n_quad_seg_face = dcube->n_g_hexa_cell_seg * dcube->n_g_hexa_cell_seg;
  double step = length / (double) dcube->n_g_hexa_cell_seg;

  /*
   * Create the dmesh_nodal that hold the resulting mesh
   */
  PDM_g_num_t n_cell_abs = _get_n_cell_abs(dcube->n_g_hexa_cell, t_elt);
  dcube->dmesh_nodal = PDM_DMesh_nodal_create(dcube->comm,
                                              3,
                                              n_vtx,
                                              n_cell_abs,   /* n_cell */
                                              -1,           /* n_face */
                                              -1);          /* n_edge */

  dcube->dn_vtx = PDM_compute_uniform_dn_entity(dcube->comm, n_vtx);

  double* dvtx_coord = (double *) malloc( 3 * (dcube->dn_vtx ) * sizeof(double *));
  PDM_DMesh_nodal_coord_set(dcube->dmesh_nodal,
                            dcube->dn_vtx,
                            dvtx_coord,
                            dcube->owner_for_dmesh_nodal); /* Le responsable de la mémoire est le dmesh_nodal */

  PDM_g_num_t* distrib_vtx = PDM_dmesh_nodal_vtx_distrib_get(dcube->dmesh_nodal);

  PDM_g_num_t _dn_vtx = distrib_vtx[i_rank+1] - distrib_vtx[i_rank];
  dcube->dn_vtx       = (int) _dn_vtx;

  printf(" _dn_vtx = %i \n", _dn_vtx);
  /*
   * Generate vertex
   */
  for(int i_vtx = 0; i_vtx < _dn_vtx; ++i_vtx) {

    PDM_g_num_t g_vtx = distrib_vtx[i_rank] + i_vtx;

    PDM_g_num_t indk =  g_vtx                                 / ( n_vtx_seg * n_vtx_seg );
    PDM_g_num_t indj = (g_vtx - indk * n_vtx_seg * n_vtx_seg) / n_vtx_seg;
    PDM_g_num_t indi = (g_vtx - indk * n_vtx_seg * n_vtx_seg - indj * n_vtx_seg);

    // printf(" g_vtx = "PDM_FMT_G_NUM" -> ["PDM_FMT_G_NUM"/"PDM_FMT_G_NUM"/"PDM_FMT_G_NUM"] \n", g_vtx, indi, indj, indk);

    dvtx_coord[3 * i_vtx    ] = indi * step + zero_x;
    dvtx_coord[3 * i_vtx + 1] = indj * step + zero_y;
    dvtx_coord[3 * i_vtx + 2] = indk * step + zero_z;

  }

  /*
   * A Prevoir ici le rajout des vtx lié à la decomposition des hexa, ou raffinement (ie HO)
   */


  /*
   * Create the real hexa
   */
  dcube->distrib_hexa = PDM_compute_uniform_entity_distribution(dcube->comm, dcube->n_g_hexa_cell);

  PDM_g_num_t n_g_quad_lim = 6 * n_quad_seg_face;
  dcube->distrib_quad_lim     = PDM_compute_uniform_entity_distribution(dcube->comm, n_g_quad_lim);
  dcube->distrib_quad_seg_lim = PDM_compute_uniform_entity_distribution(dcube->comm, n_quad_seg_face);

  dcube->n_face_group = 6;

  PDM_g_num_t _dn_hexa_cell = dcube->distrib_hexa[i_rank+1] - dcube->distrib_hexa[i_rank];
  dcube->dn_hexa_cell       = (int) _dn_hexa_cell;

  PDM_g_num_t _dn_quad_lim  = dcube->distrib_quad_lim[i_rank+1] - dcube->distrib_quad_lim[i_rank];
  dcube->dn_quad_lim        = (int) _dn_quad_lim;

  PDM_g_num_t _dn_quad_seq_lim  = dcube->distrib_quad_seg_lim[i_rank+1] - dcube->distrib_quad_seg_lim[i_rank];
  dcube->dn_quad_seq_lim        = (int) _dn_quad_seq_lim;


  switch (t_elt) {
    case PDM_MESH_NODAL_TETRA4    :
    {
      // Each hexa in split in 5 tetra and boundary
      _generate_tetra_from_hexa(dcube, dcube->dmesh_nodal);
    }
    break;

    case PDM_MESH_NODAL_PRISM6    :
    {
      _generate_prism_from_hexa(dcube, dcube->dmesh_nodal);
    }
    break;

    case PDM_MESH_NODAL_HEXA8    :
    {
      _generate_hexa_from_hexa(dcube, dcube->dmesh_nodal);
    }
    break;

    default :
      PDM_error(__FILE__, __LINE__, 0, "Unknown element type\n");
      break;

  }

  return dcube;

}

/**
 *
 * \brief Return distributed cube size
 *
 * \param [in]   id            dcube identifier
 * \param [out]  n_face_group  Number of faces groups
 * \param [out]  dn_hexa_cell       Number of cells stored in this process
 * \param [out]  dn_face       Number of faces stored in this process
 * \param [out]  dn_vtx        Number of vertices stored in this process
 * \param [out]  sface_vtx     Length of dface_vtx array
 * \param [out]  sface_group   Length of dface_group array
 *
 */

void
PDM_dcube_nodal_gen_dim_get
(
 PDM_dcube_nodal_t  *dcube,
 int                *n_face_group,
 int                *dn_hexa_cell,
 int                *dn_face,
 int                *dn_vtx,
 int                *sface_vtx,
 int                *sface_group
)
{
  *n_face_group = dcube->n_face_group;
  *dn_hexa_cell = dcube->dn_hexa_cell;
  *dn_face      = -1;
  *dn_vtx       = dcube->dn_vtx;
  *sface_vtx    = 0; // dcube->dface_vtx_idx[dcube->dn_face];
  *sface_group  = 0; // dcube->dface_group_idx[dcube->n_face_group];
}

/**
 *
 * \brief Return distributed cube data
 *
 * \param [in]  id              dcube identifier
 * \param [out] dface_cell      Faces from cells connectivity (size = 2 * dn_face)
 * \param [out] dface_vtx_idx    Faces from vertices connectivity index (size = dn_face + 1)
 * \param [out] dface_vtx       Faces from vertices connectivity (size = sface_vtx)
 * \param [out] dvtx_coord      Vertices coordinates (size = 3 * dn_vtx)
 * \param [out] dface_group_idx Faces groups index (size = n_face_group + 1)
 * \param [out] dface_group     Faces groups (size = sface_group)
 *
 */
void
PDM_dcube_nodal_gen_data_get
(
 PDM_dcube_nodal_t  *dcube,
 PDM_g_num_t       **delmt_vtx,
 double            **dvtx_coord,
 int               **dface_group_idx,
 PDM_g_num_t       **dface_group
)
{
  *delmt_vtx       = dcube->delmt_vtx;
  *dvtx_coord      = dcube->dvtx_coord;
  *dface_group_idx = dcube->dface_group_idx;
  *dface_group     = dcube->dface_group;
}


PDM_dmesh_nodal_t*
PDM_dcube_nodal_gen_dmesh_nodal_get
(
 PDM_dcube_nodal_t  *dcube
)
{
  return dcube->dmesh_nodal;
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
PDM_dcube_nodal_t        *dcube
)
{
  free(dcube->distrib_hexa);
  free(dcube->distrib_quad_lim);
  free(dcube->distrib_quad_seg_lim);

  // if(dcube->owner == PDM_OWNERSHIP_KEEP) {
  //   if (dcube->dface_cell  != NULL)
  //     free(dcube->dface_cell);

  //   if (dcube->dface_vtx_idx  != NULL)
  //     free(dcube->dface_vtx_idx);

  //   if (dcube->dface_vtx  != NULL)
  //     free(dcube->dface_vtx);

  //   if (dcube->dvtx_coord  != NULL)
  //     free(dcube->dvtx_coord);

  //   if (dcube->dface_group_idx  != NULL)
  //     free(dcube->dface_group_idx);

  //   if (dcube->dface_group  != NULL)
  //     free(dcube->dface_group);
  // }

  if(dcube->owner == PDM_OWNERSHIP_KEEP) {
    /* Si l'utilisateur fait le get il doit liberer le dmesh_nodal */
    PDM_DMesh_nodal_free(dcube->dmesh_nodal, 0);
  }

  free(dcube);

}



// /*
//  * Generate jmax plane into QUAD
//  */
// static
// void
// _decompose_into_quad_jmax
// (
//  PDM_dcube_nodal_t* dcube,
//  PDM_dmesh_nodal_t* dmesh_nodal
// )
// {
//   int i_rank;
//   PDM_MPI_Comm_rank(dcube->comm, &i_rank);

//   PDM_g_num_t n_vtx_seg  = dcube->n_vtx_seg;
//   int n_vtx_per_elmt_lim = 4;
//   PDM_g_num_t ngh        = dcube->n_g_hexa_cell_seg;

//   PDM_g_num_t* delmt_vtx = (PDM_g_num_t *) malloc( ( n_vtx_per_elmt_lim * dcube->dn_quad_seq_lim ) * sizeof(PDM_g_num_t));
//   for(int i_quad = 0; i_quad < dcube->dn_quad_seq_lim; ++i_quad) { // QUAD

//     /* Pour le remplissage du elmt_group */
//     PDM_g_num_t g_quad_lim = dcube->distrib_quad_seg_lim[i_rank] + i_quad;

//     PDM_g_num_t ipl_k = ( g_quad_lim                                    ) / ( dcube->n_g_hexa_cell_seg );
//     PDM_g_num_t ipl_j = ngh; // Donc maillage en J,K
//     PDM_g_num_t ipl_i = ( g_quad_lim - ipl_k * dcube->n_g_hexa_cell_seg );

//     delmt_vtx[n_vtx_per_elmt_lim * i_quad    ] = (ipl_i  ) + (ipl_j    ) * n_vtx_seg + ( ipl_k     ) * n_vtx_seg * n_vtx_seg + 1;
//     delmt_vtx[n_vtx_per_elmt_lim * i_quad + 1] = (ipl_i  ) + (ipl_j    ) * n_vtx_seg + ( ipl_k + 1 ) * n_vtx_seg * n_vtx_seg + 1;
//     delmt_vtx[n_vtx_per_elmt_lim * i_quad + 2] = (ipl_i+1) + (ipl_j    ) * n_vtx_seg + ( ipl_k + 1 ) * n_vtx_seg * n_vtx_seg + 1;
//     delmt_vtx[n_vtx_per_elmt_lim * i_quad + 3] = (ipl_i+1) + (ipl_j    ) * n_vtx_seg + ( ipl_k     ) * n_vtx_seg * n_vtx_seg + 1;

//     printf("JMAX :: ipl_j = "PDM_FMT_G_NUM", ipl_k = "PDM_FMT_G_NUM" \n", ipl_j, ipl_k);

//   }

//   int id_quad = PDM_DMesh_nodal_section_add(dmesh_nodal, PDM_MESH_NODAL_QUAD4);

//   PDM_DMesh_nodal_section_std_set(dmesh_nodal,
//                                   id_quad,
//                                   dcube->dn_quad_seq_lim,
//                                   delmt_vtx,
//                                   dcube->owner_for_dmesh_nodal);
// }


// /*
//  * Generate IMAX plane into QUAD
//  */
// static
// void
// _decompose_into_quad_imax
// (
//  PDM_dcube_nodal_t* dcube,
//  PDM_dmesh_nodal_t* dmesh_nodal
// )
// {
//   int i_rank;
//   PDM_MPI_Comm_rank(dcube->comm, &i_rank);

//   PDM_g_num_t n_vtx_seg  = dcube->n_vtx_seg;
//   int n_vtx_per_elmt_lim = 4;
//   PDM_g_num_t ngh        = dcube->n_g_hexa_cell_seg;

//   PDM_g_num_t* delmt_vtx = (PDM_g_num_t *) malloc( ( n_vtx_per_elmt_lim * dcube->dn_quad_seq_lim ) * sizeof(PDM_g_num_t));
//   for(int i_quad = 0; i_quad < dcube->dn_quad_seq_lim; ++i_quad) { // QUAD

//     /* Pour le remplissage du elmt_group */
//     PDM_g_num_t g_quad_lim = dcube->distrib_quad_seg_lim[i_rank] + i_quad;

//     PDM_g_num_t ipl_k = ( g_quad_lim                                    ) / ( dcube->n_g_hexa_cell_seg );
//     PDM_g_num_t ipl_j = ( g_quad_lim - ipl_k * dcube->n_g_hexa_cell_seg );
//     PDM_g_num_t ipl_i = ngh; // Donc maillage en J,K

//     delmt_vtx[n_vtx_per_elmt_lim * i_quad    ] = (ipl_i  ) + (ipl_j    ) * n_vtx_seg + ( ipl_k     ) * n_vtx_seg * n_vtx_seg + 1;
//     delmt_vtx[n_vtx_per_elmt_lim * i_quad + 1] = (ipl_i  ) + (ipl_j + 1) * n_vtx_seg + ( ipl_k     ) * n_vtx_seg * n_vtx_seg + 1;
//     delmt_vtx[n_vtx_per_elmt_lim * i_quad + 2] = (ipl_i  ) + (ipl_j + 1) * n_vtx_seg + ( ipl_k + 1 ) * n_vtx_seg * n_vtx_seg + 1;
//     delmt_vtx[n_vtx_per_elmt_lim * i_quad + 3] = (ipl_i  ) + (ipl_j    ) * n_vtx_seg + ( ipl_k + 1 ) * n_vtx_seg * n_vtx_seg + 1;

//     printf("IMAX :: ipl_j = "PDM_FMT_G_NUM", ipl_k = "PDM_FMT_G_NUM" \n", ipl_j, ipl_k);

//   }

//   int id_quad = PDM_DMesh_nodal_section_add(dmesh_nodal, PDM_MESH_NODAL_QUAD4);

//   PDM_DMesh_nodal_section_std_set(dmesh_nodal,
//                                   id_quad,
//                                   dcube->dn_quad_seq_lim,
//                                   delmt_vtx,
//                                   dcube->owner_for_dmesh_nodal);
// }

// /*
//  * Generate jmax plane into QUAD
//  */
// static
// void
// _decompose_into_quad_kmax
// (
//  PDM_dcube_nodal_t* dcube,
//  PDM_dmesh_nodal_t* dmesh_nodal
// )
// {
//   int i_rank;
//   PDM_MPI_Comm_rank(dcube->comm, &i_rank);

//   PDM_g_num_t n_vtx_seg  = dcube->n_vtx_seg;
//   int n_vtx_per_elmt_lim = 4;
//   PDM_g_num_t ngh        = dcube->n_g_hexa_cell_seg;

//   PDM_g_num_t* delmt_vtx = (PDM_g_num_t *) malloc( ( n_vtx_per_elmt_lim * dcube->dn_quad_seq_lim ) * sizeof(PDM_g_num_t));
//   for(int i_quad = 0; i_quad < dcube->dn_quad_seq_lim; ++i_quad) { // QUAD

//     /* Pour le remplissage du elmt_group */
//     PDM_g_num_t g_quad_lim = dcube->distrib_quad_seg_lim[i_rank] + i_quad;

//     PDM_g_num_t ipl_k = ngh; // Donc maillage en I,J
//     PDM_g_num_t ipl_j = ( g_quad_lim                                    ) / ( dcube->n_g_hexa_cell_seg );
//     PDM_g_num_t ipl_i = ( g_quad_lim - ipl_j * dcube->n_g_hexa_cell_seg );

//     delmt_vtx[n_vtx_per_elmt_lim * i_quad    ] = (ipl_i  ) + (ipl_j    ) * n_vtx_seg + ( ipl_k ) * n_vtx_seg * n_vtx_seg + 1;
//     delmt_vtx[n_vtx_per_elmt_lim * i_quad + 1] = (ipl_i+1) + (ipl_j    ) * n_vtx_seg + ( ipl_k ) * n_vtx_seg * n_vtx_seg + 1;
//     delmt_vtx[n_vtx_per_elmt_lim * i_quad + 2] = (ipl_i+1) + (ipl_j + 1) * n_vtx_seg + ( ipl_k ) * n_vtx_seg * n_vtx_seg + 1;
//     delmt_vtx[n_vtx_per_elmt_lim * i_quad + 3] = (ipl_i  ) + (ipl_j + 1) * n_vtx_seg + ( ipl_k ) * n_vtx_seg * n_vtx_seg + 1;
//     // printf("KMAX :: ipl_j = "PDM_FMT_G_NUM", ipl_k = "PDM_FMT_G_NUM" \n", ipl_j, ipl_k);

//   }

//   int id_quad = PDM_DMesh_nodal_section_add(dmesh_nodal, PDM_MESH_NODAL_QUAD4);

//   PDM_DMesh_nodal_section_std_set(dmesh_nodal,
//                                   id_quad,
//                                   dcube->dn_quad_seq_lim,
//                                   delmt_vtx,
//                                   dcube->owner_for_dmesh_nodal);
// }

