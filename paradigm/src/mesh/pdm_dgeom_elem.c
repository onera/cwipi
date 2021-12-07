/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_priv.h"
#include "pdm_polygon.h"
#include "pdm_plane.h"
#include "pdm_dgeom_elem.h"
#include "pdm_block_to_part.h"
#include "pdm_plane.h"
#include "pdm_printf.h"
#include "pdm_error.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type
 *============================================================================*/

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Static function definitions
 *============================================================================*/

/*=============================================================================
 * Public function definitions
 *============================================================================*/

void
PDM_compute_center_from_descending_connectivity
(
  const int         *dentity1_entity2_idx,
  const PDM_g_num_t *dentity1_entity2,
  const int          dn_entity1,
  const PDM_g_num_t *dentity2_distrib,
  double            *dentity1_coord,
  double            *dentity2_coord,
  PDM_MPI_Comm       comm
)
{
  int *dentity1_entity2_sgn = malloc(dentity1_entity2_idx[dn_entity1] * sizeof(int));
  PDM_g_num_t *dentity1_entity2_abs = malloc(dentity1_entity2_idx[dn_entity1] * sizeof(PDM_g_num_t));
  for(int i = 0; i < dentity1_entity2_idx[dn_entity1]; ++i) {
    dentity1_entity2_sgn[i] = PDM_SIGN(dentity1_entity2[i]);
    dentity1_entity2_abs[i] = PDM_ABS (dentity1_entity2[i]);
  }
  PDM_block_to_part_t *btp_entity1_coord = PDM_block_to_part_create (dentity2_distrib,
                                              (const PDM_g_num_t **) &dentity1_entity2_abs,
                                                                     &dentity1_entity2_idx[dn_entity1],
                                                                     1,
                                                                     comm);
  /*for(int i = 0; i < dentity1_entity2_idx[dn_entity1]; ++i) {
    dentity1_entity2[i]     = dentity1_entity2[i]*dentity1_entity2_sgn[i];
    }*/
  free(dentity1_entity2_sgn);
  free(dentity1_entity2_abs);

  int strid_one = 1;
  double **tmp_entity1_entity2_coord;
  PDM_block_to_part_exch2 (btp_entity1_coord,
                           3 * sizeof(double),
                           PDM_STRIDE_CST,
                           &strid_one,
                  (void *) dentity2_coord,
                           NULL,
                (void ***) &tmp_entity1_entity2_coord);
  double *dentity1_entity2_coord = tmp_entity1_entity2_coord[0];
  free(tmp_entity1_entity2_coord);
  PDM_block_to_part_free(btp_entity1_coord);

  for(int i_entity1 = 0; i_entity1 < dn_entity1; ++i_entity1) {
    dentity1_coord[3*i_entity1  ] = 0.;
    dentity1_coord[3*i_entity1+1] = 0.;
    dentity1_coord[3*i_entity1+2] = 0.;
    int n_entity2_per_entity1 = dentity1_entity2_idx[i_entity1+1] - dentity1_entity2_idx[i_entity1];
    double inv = 1./n_entity2_per_entity1;
    for(int idx_entity2 = dentity1_entity2_idx[i_entity1]; idx_entity2 < dentity1_entity2_idx[i_entity1+1]; ++idx_entity2) {
      dentity1_coord[3*i_entity1  ] += dentity1_entity2_coord[3*idx_entity2  ];
      dentity1_coord[3*i_entity1+1] += dentity1_entity2_coord[3*idx_entity2+1];
      dentity1_coord[3*i_entity1+2] += dentity1_entity2_coord[3*idx_entity2+2];
    }
    dentity1_coord[3*i_entity1  ] = dentity1_coord[3*i_entity1  ]*inv;
    dentity1_coord[3*i_entity1+1] = dentity1_coord[3*i_entity1+1]*inv;
    dentity1_coord[3*i_entity1+2] = dentity1_coord[3*i_entity1+2]*inv;

  }
  free(dentity1_entity2_coord);

}


void
PDM_compute_dface_normal
(
  const int         *dface_vtx_idx,
  const PDM_g_num_t *dface_vtx,
  const int          dn_face,
  const PDM_g_num_t *dvtx_distrib,
  double            *dvtx_coord,
  double            *dface_normal,
  PDM_MPI_Comm       comm
)
{
  int *dface_vtx_sgn = malloc(dface_vtx_idx[dn_face] * sizeof(int));
  PDM_g_num_t *dface_vtx_abs = malloc(dface_vtx_idx[dn_face] * sizeof(PDM_g_num_t));
  for(int i = 0; i < dface_vtx_idx[dn_face]; ++i) {
    dface_vtx_sgn[i] = PDM_SIGN(dface_vtx[i]);
    dface_vtx_abs[i] = PDM_ABS (dface_vtx[i]);
  }
  PDM_block_to_part_t *btp_entity1_coord = PDM_block_to_part_create (dvtx_distrib,
                                              (const PDM_g_num_t **) &dface_vtx_abs,
                                                                     &dface_vtx_idx[dn_face],
                                                                     1,
                                                                     comm);
  free(dface_vtx_sgn);
  free(dface_vtx_abs);

  int strid_one = 1;
  double **tmp_face_vtx_coord;
  PDM_block_to_part_exch2 (btp_entity1_coord,
                           3 * sizeof(double),
                           PDM_STRIDE_CST,
                           &strid_one,
                  (void *) dvtx_coord,
                           NULL,
                (void ***) &tmp_face_vtx_coord);
  double *dface_vtx_coord = tmp_face_vtx_coord[0];
  free(tmp_face_vtx_coord);
  PDM_block_to_part_free(btp_entity1_coord);

  // double* _dface_vtx_ptr = dface_vtx_coord;
  for(int i_face = 0; i_face < dn_face; ++i_face) {

    dface_normal[3*i_face  ] = i_face;
    dface_normal[3*i_face+1] = i_face;
    dface_normal[3*i_face+2] = i_face;

    int n_vtx_per_face = dface_vtx_idx[i_face+1] - dface_vtx_idx[i_face];
    PDM_plane_normal(n_vtx_per_face, dface_vtx_coord + 3*dface_vtx_idx[i_face], &dface_normal[3*i_face  ]);

    // _dface_vtx_ptr += 3 * n_vtx_per_face;

  }
  free(dface_vtx_coord);
}


#ifdef __cplusplus
}
#endif /* __cplusplus */
