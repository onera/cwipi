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
#include "pdm_logging.h"
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


PDM_part_mesh_nodal_elmts_t*
PDM_dmesh_nodal_elmts_to_part_mesh_nodal_elmts
(
 PDM_dmesh_nodal_elmts_t      *dmne,
 int                           n_part,
 int                          *pn_vtx,
 PDM_g_num_t                 **vtx_ln_to_gn,
 int                          *pn_elmt,
 PDM_g_num_t                 **elmt_ln_to_gn,
 PDM_g_num_t                 **pparent_entitity_ln_to_gn
)
{
  PDM_UNUSED(dmne);
  // PDM_UNUSED(pmne);
  PDM_part_mesh_nodal_elmts_t* pmne = PDM_part_mesh_nodal_elmts_create(dmne->mesh_dimension, n_part, dmne->comm);

  /*
   * A priori le vtx_ln_to_gn n'est pas trié
   */
  PDM_g_num_t **sorted_vtx_ln_to_gn = (PDM_g_num_t ** ) malloc( n_part * sizeof(PDM_g_num_t *));
  int         **vtx_order           = (int         ** ) malloc( n_part * sizeof(int         *));
  for(int i_part = 0; i_part < n_part; ++i_part) {
    sorted_vtx_ln_to_gn[i_part] = malloc(pn_vtx[i_part] * sizeof(PDM_g_num_t));
    vtx_order          [i_part] = malloc(pn_vtx[i_part] * sizeof(int        ));

    for(int i = 0; i < pn_vtx[i_part]; ++i) {
      sorted_vtx_ln_to_gn[i_part][i] = vtx_ln_to_gn[i_part][i];
      vtx_order          [i_part][i] = i;
    }
    PDM_sort_long(sorted_vtx_ln_to_gn[i_part], vtx_order[i_part], pn_vtx[i_part]);
    // PDM_log_trace_array_long(sorted_vtx_ln_to_gn[i_part] , pn_vtx[i_part] , "sorted_vtx_ln_to_gn :: ");
  }

  int n_section = dmne->n_section;
  PDM_g_num_t          **block_elmts_disbrib_idx = (PDM_g_num_t          ** ) malloc( n_section * sizeof(PDM_g_num_t          *));
  PDM_g_num_t          **block_elmts_connec      = (PDM_g_num_t          ** ) malloc( n_section * sizeof(PDM_g_num_t          *));
  int                  **block_elmts_n_vtx       = (int                  ** ) malloc( n_section * sizeof(int                  *));
  PDM_Mesh_nodal_elt_t **block_elmts_types       = (PDM_Mesh_nodal_elt_t ** ) malloc( n_section * sizeof(PDM_Mesh_nodal_elt_t *));
  int                  **stride_one              = (int                  ** ) malloc( n_section * sizeof(int                  *));
  int                  *pid_section              = (int                  * )  malloc( n_section * sizeof(int                   ));
  int order = 1;
  for (int i_section = 0; i_section < n_section; i_section++) {
    int id_section = dmne->sections_id[i_section];
    block_elmts_disbrib_idx[i_section] = (PDM_g_num_t *) PDM_DMesh_nodal_elmts_distrib_section_get(dmne, id_section);

    PDM_Mesh_nodal_elt_t t_elt = PDM_DMesh_nodal_elmts_section_type_get(dmne, id_section);

    stride_one[i_section] = (int * ) malloc( 1 * sizeof(int));
    stride_one[i_section][0] = 1;

    pid_section[i_section] = PDM_part_mesh_nodal_elmts_add(pmne, t_elt);

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
  int         **pelmts_stride = NULL;
  PDM_g_num_t **pelmts_connec = NULL;
  PDM_multi_block_to_part_exch2(mbtp,
                                sizeof(PDM_g_num_t),
                                PDM_STRIDE_VAR,
                                block_elmts_n_vtx,
                     (void ** ) block_elmts_connec,
                     (int  ***) &pelmts_stride,
                     (void ***) &pelmts_connec);

  int         **pelmts_stride_idx = malloc( n_part * sizeof(int *));
  for(int i_part = 0; i_part < n_part; ++i_part) {
    pelmts_stride_idx[i_part] = PDM_array_new_idx_from_sizes_int(pelmts_stride[i_part], pn_elmt[i_part]);
  }

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
   *  We don't need to exchange the section because we find it with binary search on section_distribution
   */
  int          *pelmt_by_section_n        = (int          *) malloc( (n_section+1) * sizeof(int          ));
  int         **connec                    = (int         **) malloc( (n_section+1) * sizeof(int         *));
  PDM_g_num_t **numabs                    = (PDM_g_num_t **) malloc( (n_section+1) * sizeof(PDM_g_num_t *));
  int         **parent_num                = (int         **) malloc( (n_section+1) * sizeof(int         *));
  PDM_g_num_t **sparent_entitity_ln_to_gn = (PDM_g_num_t **) malloc( (n_section+1) * sizeof(PDM_g_num_t *));

  for(int i_section = 0; i_section < n_section; ++i_section){
    sparent_entitity_ln_to_gn[i_section] = NULL;
  }

  for(int i_part = 0; i_part < n_part; ++i_part) {

    /* Reset */
    PDM_array_reset_int(pelmt_by_section_n, n_section, 0);

    for(int i_elmt = 0; i_elmt < pn_elmt[i_part]; ++i_elmt) {

      PDM_g_num_t g_num = elmt_ln_to_gn[i_part][i_elmt]-1;
      int i_section = PDM_binary_search_gap_long(g_num, dmne->section_distribution, n_section+1);

      /* We need to sort entry in each section */
      pelmt_by_section_n[i_section]++;
    }

    /* We allocate here and ownership if tranfert to PDM_part_mesh_nodal_elmts_t*/
    for(int i_section = 0; i_section < n_section; ++i_section){
      int n_elmt_in_section = pelmt_by_section_n[i_section];
      int id_section = pid_section[i_section];
      PDM_Mesh_nodal_elt_t t_elt = PDM_part_mesh_nodal_elmts_block_type_get(pmne, id_section);
      int n_vtx_per_elmt = PDM_Mesh_nodal_n_vertices_element (t_elt, order);
      connec    [i_section] = malloc( n_elmt_in_section * n_vtx_per_elmt * sizeof(int        ));
      numabs    [i_section] = malloc( n_elmt_in_section                  * sizeof(PDM_g_num_t));
      parent_num[i_section] = malloc( n_elmt_in_section                  * sizeof(int        ));

      if(pparent_entitity_ln_to_gn != NULL) {
        sparent_entitity_ln_to_gn[i_section] = malloc( n_elmt_in_section * sizeof(PDM_g_num_t));
      }

    }
    PDM_array_reset_int(pelmt_by_section_n, n_section, 0);

    /* For each section we rebuild the connectivity and the parent_num */
    for(int i_elmt = 0; i_elmt < pn_elmt[i_part]; ++i_elmt) {

      PDM_g_num_t g_num = elmt_ln_to_gn[i_part][i_elmt]-1;
      int i_section = PDM_binary_search_gap_long(g_num, dmne->section_distribution, n_section+1);
      int id_section = pid_section[i_section];

      PDM_Mesh_nodal_elt_t t_elt = PDM_part_mesh_nodal_elmts_block_type_get(pmne, id_section);
      int n_vtx_per_elmt = PDM_Mesh_nodal_n_vertices_element (t_elt, order);
      int idx_write      = pelmt_by_section_n[i_section]++;

      int idx_read_connec = pelmts_stride_idx[i_part][i_elmt];
      for(int i_vtx = 0; i_vtx < n_vtx_per_elmt; ++i_vtx){
        PDM_g_num_t vtx_g_num = pelmts_connec[i_part][idx_read_connec+i_vtx];
        int vtx_l_num = PDM_binary_search_long(vtx_g_num, sorted_vtx_ln_to_gn[i_part], pn_vtx[i_part]);
        assert(vtx_l_num != -1);
        connec[i_section][n_vtx_per_elmt*idx_write + i_vtx] = vtx_order[i_part][vtx_l_num]+1;
      }

      numabs    [i_section][idx_write] = g_num+1;
      parent_num[i_section][idx_write] = i_elmt+1;
      if(pparent_entitity_ln_to_gn != NULL) {
        sparent_entitity_ln_to_gn[i_section][idx_write] = pparent_entitity_ln_to_gn[i_part][i_elmt];
      }

    }

    /* Fill up structure */
    for(int i_section = 0; i_section < n_section; ++i_section){
      int n_elmt_in_section = pelmt_by_section_n[i_section];
      int id_section = pid_section[i_section];
      PDM_part_mesh_nodal_elmts_std_set(pmne,
                                        id_section,
                                        i_part,
                                        n_elmt_in_section,
                                        connec                   [i_section],
                                        numabs                   [i_section],
                                        parent_num               [i_section],
                                        sparent_entitity_ln_to_gn[i_section],
                                        PDM_OWNERSHIP_KEEP);

      connec                   [i_section] = NULL;
      numabs                   [i_section] = NULL;
      parent_num               [i_section] = NULL;
      sparent_entitity_ln_to_gn[i_section] = NULL;
    }
  }

  free(pelmt_by_section_n       );
  free(connec                   );
  free(parent_num               );
  free(numabs                   );
  free(sparent_entitity_ln_to_gn);

  for(int i_part = 0; i_part < n_part; ++i_part) {
    free(pelmts_connec      [i_part]);
    free(pelmts_stride      [i_part]);
    free(pelmts_types       [i_part]);
    free(pelmts_stride_idx  [i_part]);
    free(sorted_vtx_ln_to_gn[i_part]);
    free(vtx_order          [i_part]);
  }
  free(sorted_vtx_ln_to_gn);
  free(vtx_order);
  free(pelmts_connec);
  free(pelmts_stride);
  free(pelmts_types );
  free(pid_section  );
  free(pelmts_stride_idx);

  return pmne;
}



void
PDM_reverse_dparent_gnum
(
       PDM_g_num_t    *dparent_gnum,
       int            *dparent_sign,
       PDM_g_num_t    *parent_distrib,
       PDM_g_num_t    *delmt_child_distrib,
       int             n_part,
       int            *pn_parent,
       PDM_g_num_t   **pparent_gnum,
       int           **pn_child,
       PDM_g_num_t  ***pchild_gnum,
       PDM_g_num_t  ***pchild_parent_gnum,
       int          ***pchild_parent_sign,
 const PDM_MPI_Comm    comm
)
{
  // TODO : pchild_parent_gnum Keep or remove ?

  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  /* First part_to_block to map in parent block the child g_num */
  int dn_child_elmt = delmt_child_distrib[i_rank+1] - delmt_child_distrib[i_rank];
  PDM_part_to_block_t* ptb = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                      PDM_PART_TO_BLOCK_POST_MERGE,
                                                      1.,
                                                      &dparent_gnum,
                                                      NULL,
                                                      &dn_child_elmt,
                                                      1,
                                                      comm);

  int         *pblk_child_n    = (int         *) malloc( dn_child_elmt * sizeof(int        ));
  PDM_g_num_t *pblk_child_gnum = (PDM_g_num_t *) malloc( dn_child_elmt * sizeof(PDM_g_num_t));
  for(int i = 0; i < dn_child_elmt; ++i) {
    pblk_child_n   [i] = 1;
    pblk_child_gnum[i] = delmt_child_distrib[i_rank] + i + 1; // Donc le gnum de child ...
  }

  int         *blk_child_n = NULL;
  PDM_g_num_t *blk_child_gnum   = NULL;
  PDM_part_to_block_exch(ptb,
                         sizeof(PDM_g_num_t),
                         PDM_STRIDE_VAR,
                         -1,
                        &pblk_child_n,
             (void **)  &pblk_child_gnum,
                        &blk_child_n,
             (void **)  &blk_child_gnum);

  int *blk_dparent_sign_n = NULL;
  int *blk_dparent_sign = NULL;
  if(dparent_sign != NULL) {
    PDM_part_to_block_exch(ptb,
                           sizeof(int),
                           PDM_STRIDE_VAR,
                           -1,
                           &pblk_child_n,
                (void **)  &dparent_sign,
                           &blk_dparent_sign_n,
                (void **)  &blk_dparent_sign);
    free(blk_dparent_sign_n);
  }

  PDM_g_num_t n_g_parent = parent_distrib[n_rank]+1;
  PDM_g_num_t* block_distrib_tmp_idx = PDM_part_to_block_adapt_partial_block_to_block(ptb, &blk_child_n, n_g_parent);
  PDM_g_num_t* blk_dparent_gnum      = PDM_part_to_block_block_gnum_get(ptb);

  free(pblk_child_n   );
  free(pblk_child_gnum);

  /*
   * At this stage we have in each block of parent the global number of child
   * We need to resend into part
   */
  PDM_block_to_part_t* btp = PDM_block_to_part_create(block_distrib_tmp_idx,
                              (const PDM_g_num_t **)  pparent_gnum,
                                                      pn_parent,
                                                      n_part,
                                                      comm);

  int         **_pchild_n    = NULL;
  PDM_g_num_t **_pchild_gnum = NULL;
  PDM_block_to_part_exch2(btp,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_VAR,
                          blk_child_n,
                          blk_child_gnum,
                         &_pchild_n,
              (void ***) &_pchild_gnum);

  int         **_tmp_pchild_n    = NULL;
  PDM_g_num_t **_tmp_pchild_parent_gnum = NULL;
  PDM_block_to_part_exch2(btp,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_VAR,
                          blk_child_n,
                          blk_dparent_gnum,
                         &_tmp_pchild_n,
              (void ***) &_tmp_pchild_parent_gnum);

  int **_tmp_pchild_parent_sign = NULL;
  if(blk_dparent_sign != NULL) {
    for(int i_part = 0; i_part < n_part; ++i_part) {
      free(_tmp_pchild_n[i_part]);
    }
    free(_tmp_pchild_n);
    PDM_block_to_part_exch2(btp,
                            sizeof(int),
                            PDM_STRIDE_VAR,
                            blk_child_n,
                            blk_dparent_sign,
                           &_tmp_pchild_n,
                (void ***) &_tmp_pchild_parent_sign);
    free(blk_dparent_sign);
  }

  PDM_block_to_part_free(btp);
  PDM_part_to_block_free(ptb);
  free(blk_child_n);
  free(blk_child_gnum);
  free(block_distrib_tmp_idx);

  /*
   * At this stage we have for each partition the number AND the gnum of childs inside
   *          -> We sort pchild_gnum but normally it's unique
   *
   */
  *pn_child = (int *) malloc(n_part * sizeof(int));
  int* _pn_child = *pn_child;
  PDM_g_num_t **_pchild_parent_gnum = (PDM_g_num_t **) malloc(n_part * sizeof(PDM_g_num_t *));
  int **_pchild_parent_sign = NULL;
  if(blk_dparent_sign != NULL) {
    _pchild_parent_sign = (int **) malloc(n_part * sizeof(int *));
  }

  for(int i_part = 0; i_part < n_part; ++i_part) {

    int pn_child_tmp = 0;
    for(int i = 0; i < pn_parent[i_part]; ++i) {
      // printf("_pchild_n[i_part][%i] = %i \n",i, _pchild_n[i_part][i] );
      pn_child_tmp += _pchild_n[i_part][i];
      assert(_pchild_n[i_part][i] <= 1); // DOnc soit 0 soit 1
    }

    // PDM_log_trace_array_long(_pchild_gnum[i_part], pn_child_tmp, "_pchild_gnum :: ");

    int* unique_order = (int *) malloc( pn_child_tmp * sizeof(int));

    _pn_child   [i_part] = PDM_inplace_unique_long2(_pchild_gnum[i_part], unique_order, 0, pn_child_tmp-1);
    _pchild_gnum[i_part] = realloc(_pchild_gnum[i_part], _pn_child[i_part] * sizeof(PDM_g_num_t));

    _pchild_parent_gnum[i_part] = (PDM_g_num_t *) malloc( pn_child_tmp * sizeof(PDM_g_num_t));
    for(int i = 0; i < _pn_child[i_part]; ++i) {
      int idx_order = unique_order[i];
      PDM_g_num_t gnum = _tmp_pchild_parent_gnum[i_part][idx_order];
      _pchild_parent_gnum[i_part][i] = gnum;
    }

    if(blk_dparent_sign != NULL) {
      _pchild_parent_sign[i_part] = (int         *) malloc( pn_child_tmp * sizeof(int        ));
      for(int i = 0; i < _pn_child[i_part]; ++i) {
        int idx_order = unique_order[i];
        int sgn = _tmp_pchild_parent_sign[i_part][idx_order];
        _pchild_parent_sign[i_part][i] = sgn;
      }
      free(_tmp_pchild_parent_sign[i_part]);
    }
    free(_tmp_pchild_parent_gnum[i_part]);

    free(unique_order);
    free(_pchild_n[i_part]);
    free(_tmp_pchild_n[i_part]);
  }
  free(_pchild_n);
  free(_tmp_pchild_n);
  free(_tmp_pchild_parent_gnum);

  *pchild_gnum        = _pchild_gnum;
  *pchild_parent_gnum = _pchild_parent_gnum;
  if(dparent_sign != NULL) {
    *pchild_parent_sign = _pchild_parent_sign;
    free(_tmp_pchild_parent_sign);
  }
}
