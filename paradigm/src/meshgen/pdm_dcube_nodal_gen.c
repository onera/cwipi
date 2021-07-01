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
_decompose_into_bar_i
(
 PDM_dcube_nodal_t* dcube,
 PDM_dmesh_nodal_t* dmesh_nodal,
 int                plane,
 PDM_g_num_t        shift_elmt,
 PDM_g_num_t*       delmt_group
)
{
  int i_rank;
  PDM_MPI_Comm_rank(dcube->comm, &i_rank);

  PDM_g_num_t n_vtx_seg  = dcube->n_vtx_seg;
  int n_vtx_per_elmt_lim = 2;

  PDM_g_num_t* delmt_vtx = (PDM_g_num_t *) malloc( ( n_vtx_per_elmt_lim * dcube->dn_quad_seq_lim ) * sizeof(PDM_g_num_t));
  for(int i_quad = 0; i_quad < dcube->dn_quad_seq_lim; ++i_quad) { // QUAD

    /* Pour le remplissage du elmt_group */
    PDM_g_num_t g_quad_lim = dcube->distrib_quad_seg_lim[i_rank] + i_quad;
    PDM_g_num_t g_elmt     = shift_elmt + dcube->distrib_quad_seg_lim[i_rank] + i_quad; // Revient a donner le numero global de l'element (implicit par construction mais on a besoin de l'explicite pour les elmt_group)

    PDM_g_num_t ipl_j = g_quad_lim ;
    PDM_g_num_t ipl_i = plane; // Donc maillage en J,K

    delmt_vtx[n_vtx_per_elmt_lim * i_quad    ] = (ipl_i  ) + (ipl_j+1 ) * n_vtx_seg + 1;
    delmt_vtx[n_vtx_per_elmt_lim * i_quad + 1] = (ipl_i  ) + (ipl_j   ) * n_vtx_seg + 1;

    delmt_group[i_quad] = g_elmt;

    // printf(" delmt_group[%i] = %i (shift_elmt = %i)\n", i_quad, g_elmt, shift_elmt);

    // printf("_decompose_into_bar_i :: ipl_j = "PDM_FMT_G_NUM", ipl_k = "PDM_FMT_G_NUM" \n", delmt_vtx[n_vtx_per_elmt_lim * i_quad    ], ipl_j);

  }

  if(dcube->api_type == 0){
    int id_bar = PDM_DMesh_nodal_section_add(dmesh_nodal, PDM_MESH_NODAL_BAR2);
    PDM_DMesh_nodal_section_std_set(dmesh_nodal,
                                    id_bar,
                                    dcube->dn_quad_seq_lim,
                                    delmt_vtx,
                                    dcube->owner_for_dmesh_nodal);
  } else{
    int id_bar = PDM_DMesh_nodal_elmts_section_add(dmesh_nodal->ridge, PDM_MESH_NODAL_BAR2);
    printf("id_bar = %i \n", id_bar);
    PDM_DMesh_nodal_elmts_section_std_set(dmesh_nodal->ridge,
                                          id_bar,
                                          dcube->dn_quad_seq_lim,
                                          delmt_vtx,
                                          dcube->owner_for_dmesh_nodal);
  }

}

static
void
_decompose_into_bar_j
(
 PDM_dcube_nodal_t* dcube,
 PDM_dmesh_nodal_t* dmesh_nodal,
 int                plane,
 PDM_g_num_t        shift_elmt,
 PDM_g_num_t*       delmt_group
)
{
  int i_rank;
  PDM_MPI_Comm_rank(dcube->comm, &i_rank);

  PDM_g_num_t n_vtx_seg  = dcube->n_vtx_seg;
  int n_vtx_per_elmt_lim = 2;

  PDM_g_num_t* delmt_vtx = (PDM_g_num_t *) malloc( ( n_vtx_per_elmt_lim * dcube->dn_quad_seq_lim ) * sizeof(PDM_g_num_t));
  for(int i_quad = 0; i_quad < dcube->dn_quad_seq_lim; ++i_quad) { // QUAD

    /* Pour le remplissage du elmt_group */
    PDM_g_num_t g_quad_lim = dcube->distrib_quad_seg_lim[i_rank] + i_quad;
    PDM_g_num_t g_elmt     = shift_elmt + dcube->distrib_quad_seg_lim[i_rank] + i_quad; // Revient a donner le numero global de l'element (implicit par construction mais on a besoin de l'explicite pour les elmt_group)

    PDM_g_num_t ipl_j = plane;
    PDM_g_num_t ipl_i = g_quad_lim ; // Donc maillage en J,K

    delmt_vtx[n_vtx_per_elmt_lim * i_quad    ] = (ipl_i+1) + (ipl_j) * n_vtx_seg + 1;
    delmt_vtx[n_vtx_per_elmt_lim * i_quad + 1] = (ipl_i  ) + (ipl_j) * n_vtx_seg + 1;

    delmt_group[i_quad] = g_elmt;

    // printf(" delmt_group[%i] = %i (shift_elmt = %i)\n", i_quad, g_elmt, shift_elmt);

    // printf("IMIN :: ipl_j = "PDM_FMT_G_NUM", ipl_k = "PDM_FMT_G_NUM" \n", ipl_j, ipl_k);

  }

  if(dcube->api_type == 0){
      int id_bar = PDM_DMesh_nodal_section_add(dmesh_nodal, PDM_MESH_NODAL_BAR2);
      PDM_DMesh_nodal_section_std_set(dmesh_nodal,
                                      id_bar,
                                      dcube->dn_quad_seq_lim,
                                      delmt_vtx,
                                      dcube->owner_for_dmesh_nodal);
  } else{
    int id_bar = PDM_DMesh_nodal_elmts_section_add(dmesh_nodal->ridge, PDM_MESH_NODAL_BAR2);
    printf("id_bar = %i \n", id_bar);
    PDM_DMesh_nodal_elmts_section_std_set(dmesh_nodal->ridge,
                                          id_bar,
                                          dcube->dn_quad_seq_lim,
                                          delmt_vtx,
                                          dcube->owner_for_dmesh_nodal);
  }

}

/*
 * Generate IMIN plane into QUAD
 */
static
void
_decompose_into_quad_i
(
 PDM_dcube_nodal_t* dcube,
 PDM_dmesh_nodal_t* dmesh_nodal,
 int                plane,
 PDM_g_num_t        shift_elmt,
 PDM_g_num_t*       delmt_group
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
    PDM_g_num_t g_elmt     = shift_elmt + dcube->distrib_quad_seg_lim[i_rank] + i_quad; // Revient a donner le numero global de l'element (implicit par construction mais on a besoin de l'explicite pour les elmt_group)

    PDM_g_num_t ipl_k = ( g_quad_lim                                    ) / ( dcube->n_g_hexa_cell_seg );
    PDM_g_num_t ipl_j = ( g_quad_lim - ipl_k * dcube->n_g_hexa_cell_seg );
    PDM_g_num_t ipl_i = plane; // Donc maillage en J,K

    delmt_vtx[n_vtx_per_elmt_lim * i_quad    ] = (ipl_i  ) + (ipl_j    ) * n_vtx_seg + ( ipl_k     ) * n_vtx_seg * n_vtx_seg + 1;
    delmt_vtx[n_vtx_per_elmt_lim * i_quad + 1] = (ipl_i  ) + (ipl_j    ) * n_vtx_seg + ( ipl_k + 1 ) * n_vtx_seg * n_vtx_seg + 1;
    delmt_vtx[n_vtx_per_elmt_lim * i_quad + 2] = (ipl_i  ) + (ipl_j + 1) * n_vtx_seg + ( ipl_k + 1 ) * n_vtx_seg * n_vtx_seg + 1;
    delmt_vtx[n_vtx_per_elmt_lim * i_quad + 3] = (ipl_i  ) + (ipl_j + 1) * n_vtx_seg + ( ipl_k     ) * n_vtx_seg * n_vtx_seg + 1;

    delmt_group[i_quad] = g_elmt;

    // printf(" delmt_group[%i] = %i (shift_elmt = %i)\n", i_quad, g_elmt, shift_elmt);

    // printf("IMIN :: ipl_j = "PDM_FMT_G_NUM", ipl_k = "PDM_FMT_G_NUM" \n", ipl_j, ipl_k);

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
 int                plane,
 PDM_g_num_t        shift_elmt,
 PDM_g_num_t*       delmt_group
)
{
  int i_rank;
  PDM_MPI_Comm_rank(dcube->comm, &i_rank);

  PDM_g_num_t n_vtx_seg  = dcube->n_vtx_seg;
  int n_vtx_per_elmt_lim = 4;

  PDM_g_num_t* delmt_vtx = (PDM_g_num_t *) malloc( ( n_vtx_per_elmt_lim * dcube->dn_quad_seq_lim ) * sizeof(PDM_g_num_t));
  for(int i_quad = 0; i_quad < dcube->dn_quad_seq_lim; ++i_quad) { // QUAD

    PDM_g_num_t g_quad_lim = dcube->distrib_quad_seg_lim[i_rank] + i_quad;
    PDM_g_num_t g_elmt     = shift_elmt + dcube->distrib_quad_seg_lim[i_rank] + i_quad; // Revient a donner le numero global de l'element (implicit par construction mais on a besoin de l'explicite pour les elmt_group)


    PDM_g_num_t ipl_k = ( g_quad_lim                                    ) / ( dcube->n_g_hexa_cell_seg );
    PDM_g_num_t ipl_j = plane; // Donc maillage en I,K
    PDM_g_num_t ipl_i = ( g_quad_lim - ipl_k * dcube->n_g_hexa_cell_seg );

    delmt_vtx[n_vtx_per_elmt_lim * i_quad    ] = (ipl_i  ) + (ipl_j    ) * n_vtx_seg + ( ipl_k     ) * n_vtx_seg * n_vtx_seg + 1;
    delmt_vtx[n_vtx_per_elmt_lim * i_quad + 1] = (ipl_i+1) + (ipl_j    ) * n_vtx_seg + ( ipl_k     ) * n_vtx_seg * n_vtx_seg + 1;
    delmt_vtx[n_vtx_per_elmt_lim * i_quad + 2] = (ipl_i+1) + (ipl_j    ) * n_vtx_seg + ( ipl_k + 1 ) * n_vtx_seg * n_vtx_seg + 1;
    delmt_vtx[n_vtx_per_elmt_lim * i_quad + 3] = (ipl_i  ) + (ipl_j    ) * n_vtx_seg + ( ipl_k + 1 ) * n_vtx_seg * n_vtx_seg + 1;

    delmt_group[i_quad] = g_elmt;

    // printf("JMIN :: ipl_j = "PDM_FMT_G_NUM", ipl_k = "PDM_FMT_G_NUM" \n", ipl_j, ipl_k);

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
 int                plane,
 PDM_g_num_t        shift_elmt,
 PDM_g_num_t*       delmt_group
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
    PDM_g_num_t g_elmt     = shift_elmt + dcube->distrib_quad_seg_lim[i_rank] + i_quad; // Revient a donner le numero global de l'element (implicit par construction mais on a besoin de l'explicite pour les elmt_group)

    PDM_g_num_t ipl_k = plane; // Donc maillage en I,J
    PDM_g_num_t ipl_j = ( g_quad_lim                                    ) / ( dcube->n_g_hexa_cell_seg );
    PDM_g_num_t ipl_i = ( g_quad_lim - ipl_j * dcube->n_g_hexa_cell_seg );

    delmt_vtx[n_vtx_per_elmt_lim * i_quad    ] = (ipl_i  ) + (ipl_j    ) * n_vtx_seg + ( ipl_k ) * n_vtx_seg * n_vtx_seg + 1;
    delmt_vtx[n_vtx_per_elmt_lim * i_quad + 1] = (ipl_i  ) + (ipl_j + 1) * n_vtx_seg + ( ipl_k ) * n_vtx_seg * n_vtx_seg + 1;
    delmt_vtx[n_vtx_per_elmt_lim * i_quad + 2] = (ipl_i+1) + (ipl_j + 1) * n_vtx_seg + ( ipl_k ) * n_vtx_seg * n_vtx_seg + 1;
    delmt_vtx[n_vtx_per_elmt_lim * i_quad + 3] = (ipl_i+1) + (ipl_j    ) * n_vtx_seg + ( ipl_k ) * n_vtx_seg * n_vtx_seg + 1;

    delmt_group[i_quad] = g_elmt;
    // printf("kMIN :: ipl_j = "PDM_FMT_G_NUM", ipl_k = "PDM_FMT_G_NUM" \n", ipl_j, ipl_k);

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
_decompose_into_tri_i
(
 PDM_dcube_nodal_t* dcube,
 PDM_dmesh_nodal_t* dmesh_nodal,
 int                plane,
 PDM_g_num_t        shift_elmt,
 PDM_g_num_t*       delmt_group
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
    PDM_g_num_t g_elmt     = shift_elmt + 2 * (dcube->distrib_quad_seg_lim[i_rank] + i_quad); // Revient a donner le numero global de l'element (implicit par construction mais on a besoin de l'explicite pour les elmt_group)

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
    delmt_group[2*i_quad  ] = g_elmt;
    delmt_group[2*i_quad+1] = g_elmt+1;
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
_decompose_into_tri_j
(
 PDM_dcube_nodal_t* dcube,
 PDM_dmesh_nodal_t* dmesh_nodal,
 int                plane,
 PDM_g_num_t        shift_elmt,
 PDM_g_num_t*       delmt_group
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
    PDM_g_num_t g_elmt     = shift_elmt + 2 * (dcube->distrib_quad_seg_lim[i_rank] + i_quad); // Revient a donner le numero global de l'element (implicit par construction mais on a besoin de l'explicite pour les elmt_group)

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

    delmt_group[2*i_quad  ] = g_elmt;
    delmt_group[2*i_quad+1] = g_elmt+1;

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
 int                plane,
 PDM_g_num_t        shift_elmt,
 PDM_g_num_t*       delmt_group
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
    PDM_g_num_t g_elmt     = shift_elmt + 2 * (dcube->distrib_quad_seg_lim[i_rank] + i_quad); // Revient a donner le numero global de l'element (implicit par construction mais on a besoin de l'explicite pour les elmt_group)

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
    delmt_group[2*i_quad  ] = g_elmt;
    delmt_group[2*i_quad+1] = g_elmt+1;
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
_decompose_into_tri2_k
(
 PDM_dcube_nodal_t* dcube,
 PDM_dmesh_nodal_t* dmesh_nodal,
 int                plane,
 PDM_g_num_t        shift_elmt,
 PDM_g_num_t*       delmt_group
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
    PDM_g_num_t g_elmt     = shift_elmt + 2 * (dcube->distrib_quad_seg_lim[i_rank] + i_quad); // Revient a donner le numero global de l'element (implicit par construction mais on a besoin de l'explicite pour les elmt_group)

    PDM_g_num_t ipl_k = plane; // Donc maillage en I,J
    PDM_g_num_t ipl_j = ( g_quad_lim                                    ) / ( dcube->n_g_hexa_cell_seg );
    PDM_g_num_t ipl_i = ( g_quad_lim - ipl_j * dcube->n_g_hexa_cell_seg );

    delmt_vtx[2 * n_vtx_per_elmt_lim * i_quad    ] = (ipl_i     ) + (ipl_j    ) * n_vtx_seg + ( ipl_k ) * n_vtx_seg * n_vtx_seg + 1;
    delmt_vtx[2 * n_vtx_per_elmt_lim * i_quad + 1] = (ipl_i     ) + (ipl_j +1 ) * n_vtx_seg + ( ipl_k ) * n_vtx_seg * n_vtx_seg + 1;
    delmt_vtx[2 * n_vtx_per_elmt_lim * i_quad + 2] = (ipl_i + 1 ) + (ipl_j +1 ) * n_vtx_seg + ( ipl_k ) * n_vtx_seg * n_vtx_seg + 1;

    delmt_vtx[2 * n_vtx_per_elmt_lim * i_quad + 3] = (ipl_i     ) + (ipl_j    ) * n_vtx_seg + ( ipl_k ) * n_vtx_seg * n_vtx_seg + 1;
    delmt_vtx[2 * n_vtx_per_elmt_lim * i_quad + 4] = (ipl_i + 1 ) + (ipl_j +1 ) * n_vtx_seg + ( ipl_k ) * n_vtx_seg * n_vtx_seg + 1;
    delmt_vtx[2 * n_vtx_per_elmt_lim * i_quad + 5] = (ipl_i + 1 ) + (ipl_j    ) * n_vtx_seg + ( ipl_k ) * n_vtx_seg * n_vtx_seg + 1;

    // printf(" delmt_group[%i] = %i (shift_elmt = %i)\n", i_quad, g_elmt, shift_elmt);
    delmt_group[2*i_quad  ] = g_elmt;
    delmt_group[2*i_quad+1] = g_elmt+1;

    // printf("JMIN :: ipl_j = "PDM_FMT_G_NUM", ipl_k = "PDM_FMT_G_NUM" \n", ipl_j, ipl_k);

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
    case PDM_MESH_NODAL_TRIA3    :
    {
      // Each hexa in split in 5 tetra and boundary
      return n_hexa * 2;
    }
    break;

    case PDM_MESH_NODAL_QUAD4    :
    {
      // Each hexa in split in 5 tetra and boundary
      return n_hexa;
    }
    break;

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
_generate_tri_from_hexa
(
 PDM_dcube_nodal_t* dcube,
 PDM_dmesh_nodal_t* dmesh_nodal
)
{
  int n_rank, i_rank;
  PDM_MPI_Comm_rank(dcube->comm, &i_rank);
  PDM_MPI_Comm_size(dcube->comm, &n_rank);

  int n_vtx_per_elmt     = 3;

  PDM_g_num_t ngh       = dcube->n_g_hexa_cell_seg;
  PDM_g_num_t n_vtx_seg = dcube->n_vtx_seg;

  int dn_tri_cell = 2 * dcube->dn_hexa_cell;

  /* Setup volumic part */
  PDM_g_num_t* delmt_vtx = (PDM_g_num_t *) malloc( ( n_vtx_per_elmt * dn_tri_cell ) * sizeof(PDM_g_num_t));

  for(int i_cell = 0; i_cell < dcube->dn_hexa_cell; ++i_cell) { // HEXA

    PDM_g_num_t g_cell = dcube->distrib_hexa[i_rank] + i_cell;

    PDM_g_num_t indj = ( g_cell              ) / ( ngh   );
    PDM_g_num_t indi = ( g_cell - indj * ngh );

    delmt_vtx[2 * n_vtx_per_elmt * i_cell    ] = (indi  ) + (indj  ) * n_vtx_seg + 1;
    delmt_vtx[2 * n_vtx_per_elmt * i_cell + 1] = (indi+1) + (indj  ) * n_vtx_seg + 1;
    delmt_vtx[2 * n_vtx_per_elmt * i_cell + 2] = (indi  ) + (indj+1) * n_vtx_seg + 1;

    delmt_vtx[2 * n_vtx_per_elmt * i_cell + 3] = (indi+1) + (indj  ) * n_vtx_seg + 1;
    delmt_vtx[2 * n_vtx_per_elmt * i_cell + 4] = (indi+1) + (indj+1) * n_vtx_seg + 1;
    delmt_vtx[2 * n_vtx_per_elmt * i_cell + 5] = (indi  ) + (indj+1) * n_vtx_seg + 1;

  }

  int id_tri = PDM_DMesh_nodal_section_add(dmesh_nodal, PDM_MESH_NODAL_TRIA3);
  PDM_DMesh_nodal_section_std_set(dmesh_nodal,
                                  id_tri,
                                  dn_tri_cell,
                                  delmt_vtx,
                                  dcube->owner_for_dmesh_nodal);

  int         *delmt_group_idx = (int         *) malloc( (4 + 1                                   ) * sizeof(int        ));
  int dn_quad_face = dcube->distrib_quad_seg_lim[i_rank+1] - dcube->distrib_quad_seg_lim[i_rank];
  PDM_g_num_t *delmt_group     = (PDM_g_num_t *) malloc( (4 * dn_quad_face) * sizeof(PDM_g_num_t));

  delmt_group_idx[0] = 0;

  PDM_g_num_t shift_elmt = 2 * dcube->distrib_hexa[n_rank]+1; // Car les elements surfaciques se mettent après les volumiques et un tetra = 5 Hexa
  delmt_group_idx[1] = dn_quad_face;
  _decompose_into_bar_i(dcube, dmesh_nodal, 0, shift_elmt, &delmt_group[delmt_group_idx[0]]);

  shift_elmt += dcube->distrib_quad_seg_lim[n_rank];
  delmt_group_idx[2] = 2 * dn_quad_face;
  _decompose_into_bar_i(dcube, dmesh_nodal, dcube->n_g_hexa_cell_seg, shift_elmt, &delmt_group[delmt_group_idx[1]]);

  shift_elmt += dcube->distrib_quad_seg_lim[n_rank];
  delmt_group_idx[3] = 3 * dn_quad_face;
  _decompose_into_bar_j(dcube, dmesh_nodal, 0, shift_elmt, &delmt_group[delmt_group_idx[2]]);

  shift_elmt += dcube->distrib_quad_seg_lim[n_rank];
  delmt_group_idx[4] = 4 * dn_quad_face;
  _decompose_into_bar_j(dcube, dmesh_nodal, dcube->n_g_hexa_cell_seg, shift_elmt, &delmt_group[delmt_group_idx[3]]);

  // Rajout dans le dmesh_nodal
  PDM_DMesh_nodal_section_group_elmt_set(dmesh_nodal,
                                         4,
                                         delmt_group_idx,
                                         delmt_group,
                                         dcube->owner_for_dmesh_nodal);
}


static
void
_generate_quad_from_hexa
(
 PDM_dcube_nodal_t* dcube,
 PDM_dmesh_nodal_t* dmesh_nodal
)
{
  int n_rank, i_rank;
  PDM_MPI_Comm_rank(dcube->comm, &i_rank);
  PDM_MPI_Comm_size(dcube->comm, &n_rank);

  int n_vtx_per_elmt     = 4;

  PDM_g_num_t ngh       = dcube->n_g_hexa_cell_seg;
  PDM_g_num_t n_vtx_seg = dcube->n_vtx_seg;

  int dn_quad_cell = dcube->dn_hexa_cell;

  /* Setup volumic part */
  PDM_g_num_t* delmt_vtx = (PDM_g_num_t *) malloc( ( n_vtx_per_elmt * dcube->dn_hexa_cell ) * sizeof(PDM_g_num_t));

  for(int i_cell = 0; i_cell < dcube->dn_hexa_cell; ++i_cell) {

    PDM_g_num_t g_cell = dcube->distrib_hexa[i_rank] + i_cell;

    PDM_g_num_t indj = ( g_cell              ) / ( ngh   );
    PDM_g_num_t indi = ( g_cell - indj * ngh );

    delmt_vtx[n_vtx_per_elmt * i_cell    ] = (indi  ) + (indj  ) * n_vtx_seg + 1;
    delmt_vtx[n_vtx_per_elmt * i_cell + 1] = (indi+1) + (indj  ) * n_vtx_seg + 1;
    delmt_vtx[n_vtx_per_elmt * i_cell + 2] = (indi+1) + (indj+1) * n_vtx_seg + 1;
    delmt_vtx[n_vtx_per_elmt * i_cell + 3] = (indi  ) + (indj+1) * n_vtx_seg + 1;

    // printf(" i_cell = %i -> %i \n", i_cell, 2*i_cell);

  }

  printf(" ***** %i \n", dcube->api_type);
  if(dcube->api_type == 0){
    int id_quad = PDM_DMesh_nodal_section_add(dmesh_nodal, PDM_MESH_NODAL_QUAD4);
    PDM_DMesh_nodal_section_std_set(dmesh_nodal,
                                    id_quad,
                                    dn_quad_cell,
                                    delmt_vtx,
                                    dcube->owner_for_dmesh_nodal);
  } else {
    int id_quad = PDM_DMesh_nodal_elmts_section_add(dmesh_nodal->surfacic, PDM_MESH_NODAL_QUAD4);
    printf("id_quad = %i \n", id_quad);
    PDM_DMesh_nodal_elmts_section_std_set(dmesh_nodal->surfacic,
                                          id_quad,
                                          dn_quad_cell,
                                          delmt_vtx,
                                          dcube->owner_for_dmesh_nodal);
  }

  int         *delmt_group_idx = (int         *) malloc( (4 + 1                                   ) * sizeof(int        ));
  int dn_quad_face = dcube->distrib_quad_seg_lim[i_rank+1] - dcube->distrib_quad_seg_lim[i_rank];
  PDM_g_num_t *delmt_group     = (PDM_g_num_t *) malloc( (4 * dn_quad_face) * sizeof(PDM_g_num_t));

  delmt_group_idx[0] = 0;

  PDM_g_num_t shift_elmt = 0;
  if(dcube->api_type == 0) {
    shift_elmt = dcube->distrib_hexa[n_rank]+1; // Car les elements surfaciques se mettent après les volumiques et un tetra = 5 Hexa
  } else {
    shift_elmt = 1;
  }

  delmt_group_idx[1] = dn_quad_face;
  _decompose_into_bar_i(dcube, dmesh_nodal, 0, shift_elmt, &delmt_group[delmt_group_idx[0]]);

  shift_elmt += dcube->distrib_quad_seg_lim[n_rank];
  delmt_group_idx[2] = 2 * dn_quad_face;
  _decompose_into_bar_i(dcube, dmesh_nodal, dcube->n_g_hexa_cell_seg, shift_elmt, &delmt_group[delmt_group_idx[1]]);

  shift_elmt += dcube->distrib_quad_seg_lim[n_rank];
  delmt_group_idx[3] = 3 * dn_quad_face;
  _decompose_into_bar_j(dcube, dmesh_nodal, 0, shift_elmt, &delmt_group[delmt_group_idx[2]]);

  shift_elmt += dcube->distrib_quad_seg_lim[n_rank];
  delmt_group_idx[4] = 4 * dn_quad_face;
  _decompose_into_bar_j(dcube, dmesh_nodal, dcube->n_g_hexa_cell_seg, shift_elmt, &delmt_group[delmt_group_idx[3]]);

  // Rajout dans le dmesh_nodal
  if(dcube->api_type == 0){
    PDM_DMesh_nodal_section_group_elmt_set(dmesh_nodal,
                                           4,
                                           delmt_group_idx,
                                           delmt_group,
                                           dcube->owner_for_dmesh_nodal);
  } else {
    PDM_DMesh_nodal_elmts_group_set(dmesh_nodal->ridge,
                                    4,
                                    delmt_group_idx,
                                    delmt_group,
                                    dcube->owner_for_dmesh_nodal);
  }
}

static
void
_generate_tetra_from_hexa
(
 PDM_dcube_nodal_t* dcube,
 PDM_dmesh_nodal_t* dmesh_nodal
)
{
  int n_rank, i_rank;
  PDM_MPI_Comm_rank(dcube->comm, &i_rank);
  PDM_MPI_Comm_size(dcube->comm, &n_rank);

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
    // printf(" i_cell = %i -> %i \n", i_cell, 2*i_cell);

  }

  int id_tetra = PDM_DMesh_nodal_section_add(dmesh_nodal, PDM_MESH_NODAL_TETRA4);
  PDM_DMesh_nodal_section_std_set(dmesh_nodal,
                                  id_tetra,
                                  dn_tetra_cell,
                                  delmt_vtx,
                                  dcube->owner_for_dmesh_nodal);

  int         *delmt_group_idx = (int         *) malloc( (6 + 1                                   ) * sizeof(int        ));
  int dn_quad_face = dcube->distrib_quad_seg_lim[i_rank+1] - dcube->distrib_quad_seg_lim[i_rank];
  PDM_g_num_t *delmt_group     = (PDM_g_num_t *) malloc( (6 * 2 * dn_quad_face) * sizeof(PDM_g_num_t));

  delmt_group_idx[0] = 0;

  PDM_g_num_t shift_elmt = 5 * dcube->distrib_hexa[n_rank]+1; // Car les elements surfaciques se mettent après les volumiques et un tetra = 5 Hexa

  delmt_group_idx[1] = 2 * dn_quad_face;
  _decompose_into_tri_i(dcube, dmesh_nodal, 0, shift_elmt, &delmt_group[delmt_group_idx[0]]);

  shift_elmt += 2 * dcube->distrib_quad_seg_lim[n_rank];
  delmt_group_idx[2] = 4 * dn_quad_face;
  _decompose_into_tri_i(dcube, dmesh_nodal, dcube->n_g_hexa_cell_seg, shift_elmt, &delmt_group[delmt_group_idx[1]]);

  shift_elmt += 2 * dcube->distrib_quad_seg_lim[n_rank];
  delmt_group_idx[3] = 6 * dn_quad_face;
  _decompose_into_tri_j(dcube, dmesh_nodal, 0, shift_elmt, &delmt_group[delmt_group_idx[2]]);

  shift_elmt += 2 * dcube->distrib_quad_seg_lim[n_rank];
  delmt_group_idx[4] = 8 * dn_quad_face;
  _decompose_into_tri_j(dcube, dmesh_nodal, dcube->n_g_hexa_cell_seg, shift_elmt, &delmt_group[delmt_group_idx[3]]);

  shift_elmt += 2 * dcube->distrib_quad_seg_lim[n_rank];
  delmt_group_idx[5] = 10 * dn_quad_face;
  _decompose_into_tri_k(dcube, dmesh_nodal, 0, shift_elmt, &delmt_group[delmt_group_idx[4]]);

  shift_elmt += 2 * dcube->distrib_quad_seg_lim[n_rank];
  delmt_group_idx[6] = 12 * dn_quad_face;
  _decompose_into_tri_k(dcube, dmesh_nodal, dcube->n_g_hexa_cell_seg, shift_elmt, &delmt_group[delmt_group_idx[5]]);

  // Rajout dans le dmesh_nodal
  PDM_DMesh_nodal_section_group_elmt_set(dmesh_nodal,
                                         6,
                                         delmt_group_idx,
                                         delmt_group,
                                         dcube->owner_for_dmesh_nodal);
}


static
void
_generate_prism_from_hexa
(
 PDM_dcube_nodal_t* dcube,
 PDM_dmesh_nodal_t* dmesh_nodal
)
{
  int n_rank, i_rank;
  PDM_MPI_Comm_rank(dcube->comm, &i_rank);
  PDM_MPI_Comm_size(dcube->comm, &n_rank);

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

  // Prism = 4 * quad + 2 * tri
  int         *delmt_group_idx = (int         *) malloc( (6 + 1                                   ) * sizeof(int        ));
  int dn_quad_face = dcube->distrib_quad_seg_lim[i_rank+1] - dcube->distrib_quad_seg_lim[i_rank];
  PDM_g_num_t *delmt_group     = (PDM_g_num_t *) malloc( (4 * dn_quad_face + 2 * 2 * dn_quad_face) * sizeof(PDM_g_num_t));

  delmt_group_idx[0] = 0;

  PDM_g_num_t shift_elmt = 2 * dcube->distrib_hexa[n_rank]+1; // Car les elements surfaciques se mettent après les volumiques et un prism = 2 Hexa
  delmt_group_idx[1] = dn_quad_face;
  _decompose_into_quad_i(dcube, dmesh_nodal, 0, shift_elmt, &delmt_group[delmt_group_idx[0]]);

  shift_elmt += dcube->distrib_quad_seg_lim[n_rank];
  delmt_group_idx[2] = 2 * dn_quad_face;
  _decompose_into_quad_i(dcube, dmesh_nodal, dcube->n_g_hexa_cell_seg, shift_elmt, &delmt_group[delmt_group_idx[1]]);

  shift_elmt += dcube->distrib_quad_seg_lim[n_rank];
  delmt_group_idx[3] = 3 * dn_quad_face;
  _decompose_into_quad_j(dcube, dmesh_nodal, 0, shift_elmt, &delmt_group[delmt_group_idx[2]]);

  shift_elmt += dcube->distrib_quad_seg_lim[n_rank];
  delmt_group_idx[4] = 4 * dn_quad_face;
  _decompose_into_quad_j(dcube, dmesh_nodal, dcube->n_g_hexa_cell_seg, shift_elmt, &delmt_group[delmt_group_idx[3]]);

  shift_elmt += dcube->distrib_quad_seg_lim[n_rank];
  delmt_group_idx[5] = 4 * dn_quad_face + 2 * dn_quad_face;
  _decompose_into_tri2_k(dcube, dmesh_nodal, 0, shift_elmt, &delmt_group[delmt_group_idx[4]]);

  shift_elmt += 2 * dcube->distrib_quad_seg_lim[n_rank];
  delmt_group_idx[6] = 4 * dn_quad_face + 4 * dn_quad_face;
  _decompose_into_tri2_k(dcube, dmesh_nodal, dcube->n_g_hexa_cell_seg, shift_elmt, &delmt_group[delmt_group_idx[5]]);

  // Rajout dans le dmesh_nodal
  PDM_DMesh_nodal_section_group_elmt_set(dmesh_nodal,
                                         6,
                                         delmt_group_idx,
                                         delmt_group,
                                         dcube->owner_for_dmesh_nodal);

}

static
void
_generate_hexa_from_hexa
(
 PDM_dcube_nodal_t* dcube,
 PDM_dmesh_nodal_t* dmesh_nodal
)
{
  int n_rank, i_rank;
  PDM_MPI_Comm_rank(dcube->comm, &i_rank);
  PDM_MPI_Comm_size(dcube->comm, &n_rank);

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

  // Hexa = que des quad aux bords
  int         *delmt_group_idx = (int         *) malloc( (6 + 1                                   ) * sizeof(int        ));
  int dn_quad_face = dcube->distrib_quad_seg_lim[i_rank+1] - dcube->distrib_quad_seg_lim[i_rank];
  PDM_g_num_t *delmt_group     = (PDM_g_num_t *) malloc( (6 * dn_quad_face ) * sizeof(PDM_g_num_t));
  delmt_group_idx[0] = 0;

  int shift_elmt = dcube->distrib_hexa[n_rank]+1; // Car les elements surfaciques se mettent après les volumiques et on commence à 1
  delmt_group_idx[1] = dn_quad_face;
  _decompose_into_quad_i(dcube, dmesh_nodal, 0, shift_elmt, &delmt_group[delmt_group_idx[0]]);

  shift_elmt += dcube->distrib_quad_seg_lim[n_rank];
  delmt_group_idx[2] = 2 * dn_quad_face;
  _decompose_into_quad_i(dcube, dmesh_nodal, dcube->n_g_hexa_cell_seg, shift_elmt, &delmt_group[delmt_group_idx[1]]);

  shift_elmt += dcube->distrib_quad_seg_lim[n_rank];
  delmt_group_idx[3] = 3 * dn_quad_face;
  _decompose_into_quad_j(dcube, dmesh_nodal, 0, shift_elmt, &delmt_group[delmt_group_idx[2]]);

  shift_elmt += dcube->distrib_quad_seg_lim[n_rank];
  delmt_group_idx[4] = 4 * dn_quad_face;
  _decompose_into_quad_j(dcube, dmesh_nodal, dcube->n_g_hexa_cell_seg, shift_elmt, &delmt_group[delmt_group_idx[3]]);

  shift_elmt += dcube->distrib_quad_seg_lim[n_rank];
  delmt_group_idx[5] = 5 * dn_quad_face;
  _decompose_into_quad_k(dcube, dmesh_nodal, 0, shift_elmt, &delmt_group[delmt_group_idx[4]]);

  shift_elmt += dcube->distrib_quad_seg_lim[n_rank];
  delmt_group_idx[6] = 6 * dn_quad_face;
  _decompose_into_quad_k(dcube, dmesh_nodal, dcube->n_g_hexa_cell_seg, shift_elmt, &delmt_group[delmt_group_idx[5]]);

  // Rajout dans le dmesh_nodal
  PDM_DMesh_nodal_section_group_elmt_set(dmesh_nodal,
                                         6,
                                         delmt_group_idx,
                                         delmt_group,
                                         dcube->owner_for_dmesh_nodal);

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

  char *env_dmesh_nodal_api = getenv("PDM_NEW_DMESH_NODAL_API");
  int api_type = 0;
  if (env_dmesh_nodal_api != NULL) {
    if (atoi(env_dmesh_nodal_api) == 0) {
      api_type = 0;
    } else {
      api_type = 1;
    }
  }
  dcube->api_type = api_type;

  double t1 = PDM_MPI_Wtime();

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

  int dim = 3;
  if(t_elt == PDM_MESH_NODAL_TRIA3 || t_elt == PDM_MESH_NODAL_QUAD4 ) {
    dim = 2;
  }
  if(t_elt == PDM_MESH_NODAL_POINT || t_elt == PDM_MESH_NODAL_BAR2){
    PDM_error(__FILE__, __LINE__, 0, "Invalid t_elt for PDM_dcube_nodal_gen_init\n");
  }

  PDM_g_num_t n_vtx = -1;
  PDM_g_num_t n_quad_seg_face = -1;
  if(dim == 2) {
    n_vtx = n_vtx_seg  * n_vtx_seg;
    dcube->n_g_hexa_cell_seg = n_vtx_seg - 1;
    dcube->n_g_hexa_cell     = dcube->n_g_hexa_cell_seg * dcube->n_g_hexa_cell_seg;
    n_quad_seg_face          = dcube->n_g_hexa_cell_seg;
  } else {
    n_vtx = n_vtx_seg  * n_vtx_seg * n_vtx_seg;
    dcube->n_g_hexa_cell_seg = n_vtx_seg - 1;
    dcube->n_g_hexa_cell     = dcube->n_g_hexa_cell_seg * dcube->n_g_hexa_cell_seg * dcube->n_g_hexa_cell_seg;
    n_quad_seg_face          = dcube->n_g_hexa_cell_seg * dcube->n_g_hexa_cell_seg;
  }
  double step = length / (double) dcube->n_g_hexa_cell_seg;

  /*
   * Create the dmesh_nodal that hold the resulting mesh
   */
  PDM_g_num_t n_cell_abs = _get_n_cell_abs(dcube->n_g_hexa_cell, t_elt);
  dcube->dmesh_nodal = PDM_DMesh_nodal_create(dcube->comm,
                                              dim,
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

  // printf(" _dn_vtx = %i \n", _dn_vtx);
  /*
   * Generate vertex
   */
  if(dim == 2) {
    for(int i_vtx = 0; i_vtx < _dn_vtx; ++i_vtx) {

      PDM_g_num_t g_vtx = distrib_vtx[i_rank] + i_vtx;

      PDM_g_num_t indk = 0;
      PDM_g_num_t indj = (g_vtx ) / n_vtx_seg;
      PDM_g_num_t indi = (g_vtx  - indj * n_vtx_seg);

      // printf(" g_vtx = "PDM_FMT_G_NUM" -> ["PDM_FMT_G_NUM"/"PDM_FMT_G_NUM"/"PDM_FMT_G_NUM"] \n", g_vtx, indi, indj, indk);

      dvtx_coord[3 * i_vtx    ] = indi * step + zero_x;
      dvtx_coord[3 * i_vtx + 1] = indj * step + zero_y;
      dvtx_coord[3 * i_vtx + 2] = indk * step + zero_z;

    }
  } else {
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
  }

  /*
   * A Prevoir ici le rajout des vtx lié à la decomposition des hexa, ou raffinement (ie HO)
   */


  /*
   * Create the real hexa
   */
  dcube->distrib_hexa = PDM_compute_uniform_entity_distribution(dcube->comm, dcube->n_g_hexa_cell);

  PDM_g_num_t n_g_quad_lim = -1;
  if(dim == 2) {
    n_g_quad_lim = 4 * n_quad_seg_face;
    dcube->n_face_group = 4;
  } else {
    n_g_quad_lim = 6 * n_quad_seg_face;
    dcube->n_face_group = 6;
  }
  dcube->distrib_quad_lim     = PDM_compute_uniform_entity_distribution(dcube->comm, n_g_quad_lim);
  dcube->distrib_quad_seg_lim = PDM_compute_uniform_entity_distribution(dcube->comm, n_quad_seg_face);

  PDM_g_num_t _dn_hexa_cell = dcube->distrib_hexa[i_rank+1] - dcube->distrib_hexa[i_rank];
  dcube->dn_hexa_cell       = (int) _dn_hexa_cell;

  PDM_g_num_t _dn_quad_lim  = dcube->distrib_quad_lim[i_rank+1] - dcube->distrib_quad_lim[i_rank];
  dcube->dn_quad_lim        = (int) _dn_quad_lim;

  PDM_g_num_t _dn_quad_seq_lim  = dcube->distrib_quad_seg_lim[i_rank+1] - dcube->distrib_quad_seg_lim[i_rank];
  dcube->dn_quad_seq_lim        = (int) _dn_quad_seq_lim;

  PDM_DMesh_nodal_elmts_t* dmn_elmts_vol    = NULL;
  PDM_DMesh_nodal_elmts_t* dmn_elmts_surf   = NULL;
  PDM_DMesh_nodal_elmts_t* dmn_elmts_ridge  = NULL;
  PDM_DMesh_nodal_elmts_t* dmn_elmts_corner = NULL;

  // printf("dcube->n_g_hexa_cell = "PDM_FMT_G_NUM" \n", dcube->n_g_hexa_cell);
  // printf("n_g_quad_lim         = "PDM_FMT_G_NUM" \n", n_g_quad_lim);
  // printf("n_quad_seg_face      = "PDM_FMT_G_NUM" \n", n_quad_seg_face);
  if(dim == 2) {
    dmn_elmts_surf   = PDM_DMesh_nodal_elmts_create(dcube->comm, 2, dcube->n_g_hexa_cell );
    dmn_elmts_ridge  = PDM_DMesh_nodal_elmts_create(dcube->comm, 1, n_g_quad_lim);
    // dmn_elmts_corner = PDM_DMesh_nodal_elmts_create(dcube->comm, 0, n_ridge);
  } else {
    dmn_elmts_vol    = PDM_DMesh_nodal_elmts_create(dcube->comm, 3, dcube->n_g_hexa_cell );
    dmn_elmts_surf   = PDM_DMesh_nodal_elmts_create(dcube->comm, 2, n_g_quad_lim );
    // dmn_elmts_ridge  = PDM_DMesh_nodal_elmts_create(dcube->comm, 1, n_ridge);
    // dmn_elmts_corner = PDM_DMesh_nodal_elmts_create(dcube->comm, 0, n_ridge);
  }

  if(dmn_elmts_vol != NULL) {
    PDM_Mesh_nodal_add_desh_nodal_elmts(dcube->dmesh_nodal, dmn_elmts_vol   );
  }
  if(dmn_elmts_surf != NULL) {
    PDM_Mesh_nodal_add_desh_nodal_elmts(dcube->dmesh_nodal, dmn_elmts_surf  );
  }
  if(dmn_elmts_ridge != NULL) {
    PDM_Mesh_nodal_add_desh_nodal_elmts(dcube->dmesh_nodal, dmn_elmts_ridge );
  }
  if(dmn_elmts_corner != NULL) {
    PDM_Mesh_nodal_add_desh_nodal_elmts(dcube->dmesh_nodal, dmn_elmts_corner);
  }

  switch (t_elt) {
    case PDM_MESH_NODAL_TRIA3    :
    {
      // Each hexa in split in 5 tetra and boundary
      _generate_tri_from_hexa(dcube, dcube->dmesh_nodal);
    }
    break;

    case PDM_MESH_NODAL_QUAD4    :
    {
      // Each hexa in split in 5 tetra and boundary
      _generate_quad_from_hexa(dcube, dcube->dmesh_nodal);
    }
    break;

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

  return dcube;
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

  if(dcube->owner == PDM_OWNERSHIP_KEEP) {
    /* Si l'utilisateur fait le get il doit liberer le dmesh_nodal */
    PDM_DMesh_nodal_free(dcube->dmesh_nodal, 0);
  }

  free(dcube);

}




