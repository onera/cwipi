/*============================================================================
 * Parallel partitioning
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  System headers
 *----------------------------------------------------------------------------*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"
#include "pdm_config.h"
#include "pdm_priv.h"
#include "pdm_timer.h"
#include "pdm_mpi.h"
#include "pdm_mpi_ext_dependencies.h"
#include "pdm_part_to_block.h"
#include "pdm_block_to_part.h"
#include "pdm_multi_block_to_part.h"

#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_sort.h"
#include "pdm_unique.h"
#include "pdm_binary_search.h"
#include "pdm_hash_tab.h"
#include "pdm_array.h"
#include "pdm_dmesh_nodal_priv.h"
#include "pdm_part_mesh_nodal_priv.h"
#include "pdm_part_mesh_nodal_elmts_priv.h"

#include "pdm_partitioning_nodal_algorithm.h"
#include "pdm_distrib.h"
#include "pdm_order.h"
// #include "pdm_para_graph_dual.h"

/*----------------------------------------------------------------------------
 *  Optional headers
 *----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Fortran function header
 *============================================================================*/

/*============================================================================
 * Local macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Global variable
 *============================================================================*/


/*============================================================================
 * Private function definitions
 *============================================================================*/

/*=============================================================================
 * Public function definitions
 *============================================================================*/


void
PDM_dmesh_nodal_elmts_to_part_mesh_nodal_elmts
(
 PDM_dmesh_nodal_elmts_t      *dmne,
 PDM_part_mesh_nodal_elmts_t  *pmne,
 int                           n_part,
 int                          *pn_elmt,
 PDM_g_num_t                 **elmt_ln_to_gn
)
{
  PDM_UNUSED(dmne);
  PDM_UNUSED(pmne);

  int n_section = dmne->n_section;
  PDM_g_num_t          **block_elmts_disbrib_idx = (PDM_g_num_t          ** ) malloc( n_section * sizeof(PDM_g_num_t          *));
  PDM_g_num_t          **block_elmts_connec      = (PDM_g_num_t          ** ) malloc( n_section * sizeof(PDM_g_num_t          *));
  int                  **block_elmts_n_vtx       = (int                  ** ) malloc( n_section * sizeof(int                  *));
  PDM_Mesh_nodal_elt_t **block_elmts_types       = (PDM_Mesh_nodal_elt_t ** ) malloc( n_section * sizeof(PDM_Mesh_nodal_elt_t *));
  int                  **stride_one              = (int                  ** ) malloc( n_section * sizeof(int                  *));
  int order = 1;
  for (int i_section = 0; i_section < n_section; i_section++) {
    int id_section = dmne->sections_id[i_section];
    block_elmts_disbrib_idx[i_section] = (PDM_g_num_t *) PDM_DMesh_nodal_elmts_distrib_section_get(dmne, id_section);

    PDM_Mesh_nodal_elt_t t_elt = PDM_DMesh_nodal_elmts_section_type_get(dmne, id_section);

    stride_one[i_section] = (int * ) malloc( 1 * sizeof(int));
    stride_one[i_section][0] = 1;

    switch (t_elt) {
      case PDM_MESH_NODAL_POINT:
      case PDM_MESH_NODAL_BAR2:
      case PDM_MESH_NODAL_TRIA3:
      case PDM_MESH_NODAL_QUAD4:
      case PDM_MESH_NODAL_TETRA4:
      case PDM_MESH_NODAL_PYRAMID5:
      case PDM_MESH_NODAL_PRISM6:
      case PDM_MESH_NODAL_HEXA8:
      {
        int n_elt           = PDM_DMesh_nodal_elmts_section_n_elt_get(dmne, id_section);
        block_elmts_n_vtx[i_section] = (int                  * ) malloc( n_elt * sizeof(int                 ));
        block_elmts_types[i_section] = (PDM_Mesh_nodal_elt_t * ) malloc( n_elt * sizeof(PDM_Mesh_nodal_elt_t));
        int n_vtx_per_elmt = PDM_Mesh_nodal_n_vertices_element (t_elt, order);
        for(int i = 0; i < n_elt; ++i) {
          block_elmts_n_vtx[i_section][i] = n_vtx_per_elmt;
          block_elmts_types[i_section][i] = t_elt;
        }
        block_elmts_connec[i_section] = PDM_DMesh_nodal_elmts_section_std_get(dmne, id_section);

        break;
      }
      case PDM_MESH_NODAL_POLY_2D:
      {
        PDM_error(__FILE__, __LINE__, 0, "Error PDM_sections_decompose_edges : Element type is not taking int account\n");
        break;
      }

      case PDM_MESH_NODAL_POLY_3D:
      {
        PDM_error(__FILE__, __LINE__, 0, "Error PDM_sections_decompose_edges : Element type is not taking int account\n");
        break;
      }

      default:
        PDM_error(__FILE__, __LINE__, 0, "Error PDM_sections_decompose_edges : Element type is not taking int account\n");
    }

  }

  // PDM_g_num_t* volumic_distrib =
  PDM_multi_block_to_part_t* mbtp = PDM_multi_block_to_part_create(dmne->section_distribution,
                                                                   n_section,
                                            (const PDM_g_num_t **) block_elmts_disbrib_idx,
                                            (const PDM_g_num_t **) elmt_ln_to_gn,
                                            (const int          *) pn_elmt,
                                                                   n_part,
                                                                   dmne->comm);

  free(block_elmts_disbrib_idx);

  /*
   * Exchange connectivity
   */
  int         **pelmts_stride;
  PDM_g_num_t **pelmts_connec;
  PDM_multi_block_to_part_exch2(mbtp,
                                sizeof(PDM_g_num_t),
                                PDM_STRIDE_VAR,
                                block_elmts_n_vtx,
                     (void ** ) block_elmts_connec,
                     (int  ***) &pelmts_stride,
                     (void ***) &pelmts_connec);

  /*
   * Exchange type of elements
   */
  PDM_Mesh_nodal_elt_t **pelmts_types;
  PDM_multi_block_to_part_exch2(mbtp,
                                sizeof(PDM_Mesh_nodal_elt_t),
                                PDM_STRIDE_CST,
                                stride_one,
                     (void ** ) block_elmts_types,
                                NULL,
                     (void ***) &pelmts_types);


  for (int i_section = 0; i_section < n_section; i_section++) {
    free(block_elmts_n_vtx[i_section]);
    free(block_elmts_types[i_section]);
    free(stride_one[i_section]);
  }
  free(block_elmts_n_vtx);
  free(block_elmts_types);
  free(block_elmts_connec);
  free(stride_one);

  PDM_multi_block_to_part_free(mbtp);

  /*
   * Rebuild the associate part_mesh_nodal
   *         --> CAUTION WE NEED TO REBUILD IN THE SAME ORDER FOR EACH PROC AND EACH PARTITION !!!!
   *         -->
   */
  abort(); // We need to setup the same skeloton over all procs !!!




  for(int i_part = 0; i_part < n_part; ++i_part) {
    free(pelmts_connec[i_part]);
    free(pelmts_stride[i_part]);
    free(pelmts_types [i_part]);
  }
  free(pelmts_connec);
  free(pelmts_stride);
  free(pelmts_types );


}
