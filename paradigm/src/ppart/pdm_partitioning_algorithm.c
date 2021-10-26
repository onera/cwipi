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

#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_sort.h"
#include "pdm_unique.h"
#include "pdm_binary_search.h"
#include "pdm_hash_tab.h"
#include "pdm_array.h"

#include "pdm_partitioning_algorithm.h"
#include "pdm_distrib.h"
#include "pdm_order.h"
// #include "pdm_para_graph_dual.h"
#include "pdm_logging.h"

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

static inline int _is_prime(int num)
{
  if ((num & 1)==0)
    return num == 2;
  else {
    for (int i = 3; i <= sqrt(num); i+=2) {
      if (num % i == 0)
        return 0;
    }
  }
  return 1;
}

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *  \brief Gather the entities splitted by the partitioner
 *   (usually cells) to their attributed partition, using the array mapping
 *   entities id to their assigned partition number.
 *   Each partition is hold by a (unique) process following the input partition
 *   distribution. The connection between partition members and original entities
 *   is made trought the local to global numbering computed by the function.
 *
 * \param [in]   comm                PDM_MPI communicator
 * \param [in]   part_distribution   Distribution of partitions over the processes (size=n_rank+1)
 * \param [in]   entity_distribution Distribution of entities over the processes (size=n_rank+1)
 * \param [in]   dentity_to_part     Id of assigned partition for each entity (size=dn_entity)
 * \param [out]  pn_entities         Number of entities in each partition (size n_part)
 * \param [out]  pentity_ln_to_gn    Array of local to global entity id for each partition (size n_part)
 *
 * \return       n_part              Number of partitions managed by this process
*/
int
PDM_part_assemble_partitions
(
 const PDM_MPI_Comm    comm,
       PDM_g_num_t    *part_distribution,
 const PDM_g_num_t    *entity_distribution,
 const int            *dentity_to_part,
       int           **pn_entity,
       PDM_g_num_t  ***pentity_ln_to_gn
)
{
  int i_rank;
  int n_rank;

  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  int dn_part   = part_distribution  [i_rank+1] - part_distribution  [i_rank];
  int dn_entity = entity_distribution[i_rank+1] - entity_distribution[i_rank];

  /*
   * On recréer un tableau pourtant dans scotch et metis le part est en int64 ...
   *     Il faudrait changer le ext_dependancies ..
   */
  PDM_g_num_t* dpart_ln_to_gn = (PDM_g_num_t * ) malloc( sizeof(PDM_g_num_t) * dn_entity );
  for(int i = 0; i < dn_entity; ++i){
    dpart_ln_to_gn[i] = (PDM_g_num_t) dentity_to_part[i] + 1;
  }

  /*
   * Each proc get all the entities affected to its partitions
   */
  PDM_part_to_block_t *ptb_partition =
   PDM_part_to_block_create2(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                             PDM_PART_TO_BLOCK_POST_MERGE,
                             1.,
                             &dpart_ln_to_gn,
                             part_distribution,
                             &dn_entity,
                             1,
                             comm);

  const int n_part_block = PDM_part_to_block_n_elt_block_get (ptb_partition);
  //printf("rank = %i | n_part_block = %i\n",i_rank,n_part_block);
  //printf("rank = %i | dn_part = %i\n",i_rank,dn_part);
  assert(n_part_block == dn_part || n_part_block==0);
  /*
   * Generate global numbering
   */
  int*             dentity_stri = (int         *) malloc( sizeof(int        ) * dn_entity );
  PDM_g_num_t* dentity_ln_to_gn = (PDM_g_num_t *) malloc( sizeof(PDM_g_num_t) * dn_entity );

  PDM_g_num_t shift_g = entity_distribution[i_rank] + 1;
  for(int i = 0; i < dn_entity; ++i){
    dentity_stri    [i] = 1;
    dentity_ln_to_gn[i] = (PDM_g_num_t) shift_g + i;
  }

  /*
   * Exchange
   */
  PDM_g_num_t* pentity_ln_to_gn_tmp = NULL;

  // C'est peut être un multiblock_to_part ???
  PDM_part_to_block_exch (ptb_partition,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_VAR,
                          1,
                          &dentity_stri,
                (void **) &dentity_ln_to_gn,
                          pn_entity,
                (void **) &pentity_ln_to_gn_tmp);

  free(dentity_stri);
  free(dentity_ln_to_gn);

  /* Reshape pentity_ln_to_gn */
  *pentity_ln_to_gn = (PDM_g_num_t **) malloc( sizeof(PDM_g_num_t *) * n_part_block);
  PDM_g_num_t **_pentity_ln_to_gn = *pentity_ln_to_gn;

  int offset = 0;
  for(int i_part = 0; i_part < n_part_block; ++i_part){

    int _pn_entity = (*pn_entity)[i_part];
    _pentity_ln_to_gn[i_part] = (PDM_g_num_t *) malloc( sizeof(PDM_g_num_t) * _pn_entity);

    for(int i_elmt = 0; i_elmt < _pn_entity; ++i_elmt){
      _pentity_ln_to_gn[i_part][i_elmt] = pentity_ln_to_gn_tmp[offset + i_elmt];
    }

    offset += _pn_entity;

    /* Panic verbose */
    if(0 == 1){
      printf("[%i] _pentity_ln_to_gn = ", i_rank);
      for(int i_data = 0; i_data < _pn_entity; ++i_data){
        printf(PDM_FMT_G_NUM" ", _pentity_ln_to_gn[i_part][i_data]);
      }
      printf("\n");
    }
  }

  if (n_part_block==0) { // in this case, the ln_to_gn is empty, but the size (i.e. 0) still needs to be there
    free(*pn_entity);
    *pn_entity = PDM_array_zeros_int(dn_part);
  }

  PDM_part_to_block_free (ptb_partition);

  free(pentity_ln_to_gn_tmp);
  free(dpart_ln_to_gn);

  return n_part_block;

}

/**
 *  \brief Construct the face->cell connectivity from the cell->face connectivity
 *   Since we assume that a face->cell is an array of two columns, this function
 *   can not be used to reverse other connectivities such as face->vtx.
 *
 *   This function respect the orientation : negative face in cell->face connectivity
 *   indicates that cell is right in the face->cell connectivity.
 *
 * \param [in]   n_part             Number of partitions
 * \param [in]   np_cell            Number of cells in each partition (size=n_part)
 * \param [in]   np_face            Number of faces in each partition (size=n_part)
 * \param [in]   pcell_face_idx     2d array of cell to face connectivity indexes
 *                                  (size = n_part*np_cell[i_part])
 * \param [in]   pcell_face         2d array of cell to face connectivity
 *                                  (size = n_part*pcell_face_idx[i_part][np_cell+1])
 * \param [out]  pface_face         2d array of face to cell connectivity
 *                                  (size = n_part*2*np_face[i_part])
*/
void
PDM_part_reverse_pcellface
(
  const int         n_part,
  const int        *np_cell,
  const int        *np_face,
  const int       **pcell_face_idx,
  const int       **pcell_face,
        int      ***pface_cell
)
{
  *pface_cell = (int **) malloc( n_part * sizeof(int *) );

  for (int i_part = 0; i_part < n_part; i_part++)
  {
    (*pface_cell)[i_part] = PDM_array_zeros_int(2 * np_face[i_part]);
    /* Shortcuts for this part */
    const int *_pcell_face_idx = pcell_face_idx[i_part];
    const int *_pcell_face     = pcell_face[i_part];
          int *_pface_cell     = (*pface_cell)[i_part];


    for (int i_cell = 0; i_cell < np_cell[i_part]; i_cell++){
      for (int i_data = _pcell_face_idx[i_cell]; i_data < _pcell_face_idx[i_cell+1]; i_data++){
        int l_face =  PDM_ABS (_pcell_face[i_data])-1;
        int sign    = PDM_SIGN(_pcell_face[i_data]);
        if (sign > 0) {
          _pface_cell[2*l_face  ] = i_cell+1;
        } else {
          _pface_cell[2*l_face+1] = i_cell+1;
        }
      }
    }
  }
}

/**
 *  \brief Reorient the boundary faces such that they have a outward normal for the boundary cell.
 *   This functions only uses topological information (face->cell connectivity) to determine
 *   if the face must be reoriented. Thus, input face->cell and cell->face connectivities must
 *   contain correct orientation information.
 *
 *
 * \param [in]    n_part             Number of partitions
 * \param [in]    pn_face            Number of faces in each partition (size=n_part)
 * \param [inout] pface_face         On each part, array of face to cell connectivity
 *                                   (size = n_part*2*np_face[i_part])
 * \param [in]    pcell_face_idx     On each part, array of cell to face connectivity indexes
 *                                   (size = n_part*np_cell[i_part])
 * \param [inout] pcell_face         On each part, array of cell to face connectivity
 *                                   (size = n_part*pcell_face_idx[i_part][np_cell+1])
 * \param [in]    pface_vtx_idx      On each part, array of face to vertex connectivity indexes
 *                                   (size = n_part*np_vtx[i_part])
 * \param [inout] pface_vtx          On each part, array of face to vertex connectivity
 *                                  (size = n_part*pface_vtx_idx[i_part][np_vtx+1])
*/
void
PDM_part_reorient_bound_faces
(
  const int         n_part,
  const int        *pn_face,
        int       **pface_cell,
  const int       **pcell_face_idx,
        int       **pcell_face,
  const int       **pface_vtx_idx,
        int       **pface_vtx,
        int       **pface_edge_idx,
        int       **pface_edge
)
{
  // printf("PDM_part_reorient_faces\n");

  for (int i_part = 0; i_part < n_part; i_part++){
    for (int i_face = 0; i_face < pn_face[i_part]; i_face++){
      if (pface_cell[i_part][2*i_face] == 0) {

        /* Swap connectivity in pface_cell */
        pface_cell[i_part][2*i_face]   = pface_cell[i_part][2*i_face+1];
        pface_cell[i_part][2*i_face+1] = 0;

        /* Remove sign in cell_face connectivity */
        int cell_lid = pface_cell[i_part][2*i_face]-1;
        for (int j = pcell_face_idx[i_part][cell_lid]; j < pcell_face_idx[i_part][cell_lid+1]; j++){
          if (pcell_face[i_part][j] == -1*(i_face+1)){
            pcell_face[i_part][j] = (i_face+1);
          }
        }

        /* Reverse face->vertex connectivity : start from same vtx but reverse array*/
        int offset     = pface_vtx_idx[i_part][i_face] + 1;
        int n_face_vtx = pface_vtx_idx[i_part][i_face+1] - pface_vtx_idx[i_part][i_face] - 1;
        int end_index  = n_face_vtx - 1;
        int tmp_swap;
        for (int j = 0; j < n_face_vtx/2; j++){
          tmp_swap = pface_vtx[i_part][offset+end_index];
          pface_vtx[i_part][offset+end_index] = pface_vtx[i_part][offset+j];
          pface_vtx[i_part][offset+j] = tmp_swap;
          end_index--;
        }

        /* Change sign of edge */
        if(pface_edge_idx != NULL) {
          for(int idx_edge = pface_edge_idx[i_part][i_face]; idx_edge < pface_edge_idx[i_part][i_face+1]; ++idx_edge) {
            pface_edge[i_part][idx_edge] = -pface_edge[i_part][idx_edge];
          }
        }

      }
    }
  }
}

/**
 *  \brief Recover partitioned entity groups (cell, face, vertex) from distributed
 *   entity groups. Return the list of local element id belonging to each group,
 *   and the position of those entities in the corresponding original (distributed) group.
 *
 *   This function is especially used to retrieve boundary conditions which are defined as
 *   face groups.
 *
 *   n_group is a global data that must be know by each process, even if they
 *   dont hold any group element.
 *
 * \param [in]   comm                PDM_MPI communicator
 * \param [in]   entity_distribution Distribution of entities over the processes (size=n_rank+1)
 * \param [in]   n_group             Number of groups defined for this entity
 * \param [in]   dgroup_idx          Number of distributed elements in each group (size=n_group+1)
 * \param [in]   dgroup              Global id of entities belonging to the groups (size=dgroup_idx[n_group])
 * \param [in]   n_part              Number of partitions
 * \param [in]   pn_entities         Number of entities in each partition (size=n_part)
 * \param [in]   pentity_ln_to_gn    Array of local to global entity id for each partition
 *                                   (size=n_part, each component size = pn_entities[i_part])
 * \param [out]  pgroup_idx          For each part, number of partitioned elements in each group
 *                                   (size = n_part, each component size = n_group+1)
 * \param [out]  pgroup              For each part, local id of entities belonging to the groups
 *                                   (size = n_part, each component size = pgroup_idx[n_group])
 * \param [out]  pgroup_ln_to_gn     For each part, position of entity in the original groups
 *                                   (size = n_part, each component size = pgroup_idx[n_group])
*/
void
PDM_part_distgroup_to_partgroup
(
 const PDM_MPI_Comm      comm,
 const PDM_g_num_t      *entity_distribution,
 const int               n_group,
 const int              *dgroup_idx,
 const PDM_g_num_t      *dgroup,
 const int               n_part,
 const int              *pn_entity,
 const PDM_g_num_t     **pentity_ln_to_gn,
       int            ***pgroup_idx,
       int            ***pgroup,
       PDM_g_num_t    ***pgroup_ln_to_gn
)
{
  // printf("PDM_part_distgroup_to_partgroup\n");
  int i_rank;
  int n_rank;

  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  int dn_entity = entity_distribution[i_rank+1] -  entity_distribution[i_rank];

  /*
   * Compute the groups distribution - Mandatory to have the link between rank and group
   */
  PDM_g_num_t** groups_distribution = (PDM_g_num_t **) malloc( n_group * sizeof(PDM_g_num_t *));

  for(int i_group = 0; i_group < n_group; ++i_group) {
    int dn_entity_group = dgroup_idx[i_group+1] - dgroup_idx[i_group];
    groups_distribution[i_group] = PDM_compute_entity_distribution(comm, dn_entity_group);
  }

  /*
   * Create exchange protocol
   *    - dgroup contains the absolute number of faces that contains the boundary conditions -> A kind of ln_to_gn
   */
  int dgroup_tot_size = dgroup_idx[n_group];
  PDM_part_to_block_t *ptb_group =
   PDM_part_to_block_create2(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                             PDM_PART_TO_BLOCK_POST_MERGE,
                             1.,
                             (PDM_g_num_t **) & dgroup,
                             entity_distribution,
                             &dgroup_tot_size,
                             1,
                             comm);

  /*
   * Prepare send
   *    Reminder :face_group have double indirection
   *       - For each group we have a list of face
   *       - All rank have a distribution of each group
   * part_data should contains : ( position in distributed array + group_id )
   */
  int*         part_stri = (int         *) malloc( sizeof(int        )     * dgroup_tot_size );
  PDM_g_num_t* part_data = (PDM_g_num_t *) malloc( sizeof(PDM_g_num_t) * 2 * dgroup_tot_size );

  int idx_data = 0;
  int idx_stri = 0;
  for(int i_group = 0; i_group < n_group; ++i_group) {
    int dn_entity_group = dgroup_idx[i_group+1] - dgroup_idx[i_group];
    for(int i_elmt = 0; i_elmt < dn_entity_group; ++i_elmt ) {
      part_data[idx_data++] = i_group;
      part_data[idx_data++] = i_elmt + groups_distribution[i_group][i_rank]+1; /* Numero dans le tableau distribue */
      part_stri[idx_stri++] = 2;
    }
  }

  /*
   * Exchange
   */
  const int n_entity_block = PDM_part_to_block_n_elt_block_get(ptb_group);
  PDM_g_num_t* blk_gnum  = PDM_part_to_block_block_gnum_get(ptb_group);

  int*         blk_stri = NULL;
  PDM_g_num_t* blk_data = NULL;
  PDM_part_to_block_exch (ptb_group,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_VAR,
                          1,
                          &part_stri,
                (void **) &part_data,
                          &blk_stri,
                (void **) &blk_data);

  /*
   * Post-treatement
   */
  if(1 == 0){
    int idx_elmt_d = 0;
    for(int i_elmt = 0; i_elmt < n_entity_block; ++i_elmt) {
      printf("[%d] - %d | "PDM_FMT_G_NUM" --> ", i_elmt, blk_stri[i_elmt], blk_gnum[i_elmt]);
      for(int i_data = 0; i_data < blk_stri[i_elmt]; ++i_data) {
        printf(PDM_FMT_G_NUM" ", blk_data[idx_elmt_d++]);
      }
      printf("\n");
    }
    printf("n_entity_block::%d\n", n_entity_block);
    printf("idx_entity_d  ::%d\n", idx_elmt_d/2);
  }

  /*
   * No choice we need to rebuild a proper blk_stri (but not blk_data)
   */
  int* blk_stri_full = PDM_array_zeros_int(dn_entity);

  /* On remet la stri au bonne endroit */
  // int idx_face = 0;
  for(int i_elmt = 0; i_elmt < n_entity_block; ++i_elmt) {
    int l_elmt = blk_gnum[i_elmt] - entity_distribution[i_rank] - 1;
    // printf("[%d] --> %d \n", i_elmt, l_elmt);
    blk_stri_full[l_elmt] = blk_stri[i_elmt];
  }
  free(blk_stri);



  /*
   *  Maintenant on peut faire un block_to_part sur le entity_ln_to_gn (le vrai) pour recuperer pour chaque partition
   *    le group
   */
  if(0 == 1){
    for(int i_part = 0; i_part < n_part; ++i_part) {
      printf("[%i][%d] pentity_ln_to_gn:: ", i_rank, i_part);
      for(int i_elmt = 0; i_elmt < pn_entity[i_part]; ++i_elmt){
        printf(PDM_FMT_G_NUM" ", pentity_ln_to_gn[i_part][i_elmt]);
      }
      printf("\n");
    }
  }

  /*
   * Prepare exchange protocol
   */
  PDM_block_to_part_t* btp = PDM_block_to_part_create(entity_distribution,
                               (const PDM_g_num_t **) pentity_ln_to_gn,
                                                      pn_entity,
                                                      n_part,
                                                      comm);

  /*
   * Exchange
   */
  int**         part_group_stri;
  PDM_g_num_t** part_group_data;
  PDM_block_to_part_exch2(btp,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_VAR,
                          blk_stri_full,
             (void *  )   blk_data,
             (int  ***)  &part_group_stri,
             (void ***)  &part_group_data);

  /*
   * Post-Treatment
   */
  if(0 == 1){
    for(int i_part = 0; i_part < n_part; ++i_part) {
      int idx_data_g = 0;
      printf("[%d] part_group_data --> ", i_part);
      for(int i_elmt = 0; i_elmt < pn_entity[i_part]; ++i_elmt){
        printf(" \t [%d] ", i_elmt);
        for(int i_data = 0; i_data < part_group_stri[i_part][i_elmt]; ++i_data) {
          printf(PDM_FMT_G_NUM" ", part_group_data[i_part][idx_data_g++]);
        }
        printf("\n");
      }
    }
  }

  /*
   * Donc sur chaque partition on a dans la numerotation des faces la liste des frontières
   * On doit maintenant recréer les tableaux identiques à pdm_part
   */
  *pgroup_ln_to_gn = (PDM_g_num_t **) malloc( n_part * sizeof(PDM_g_num_t *) );
  *pgroup          = (int         **) malloc( n_part * sizeof(int         *) );
  *pgroup_idx      = (int         **) malloc( n_part * sizeof(int         *) );

  /* Shortcut */
  PDM_g_num_t** _group_ln_to_gn = *pgroup_ln_to_gn;
  int**         _group          = *pgroup;
  int**         _group_idx      = *pgroup_idx;

  int* count_group = (int *) malloc( (n_group + 1) * sizeof(int));
  for(int i_part = 0; i_part < n_part; ++i_part) {

    /*
     * All partition allocate is own array of size n_group
     */
    _group_idx[i_part] = (int *) malloc( (n_group + 1) * sizeof(int) );
    for(int i_group = 0; i_group < n_group+1; ++i_group){
      _group_idx[i_part][i_group] = 0;
      count_group[i_group] = 0;
    }

    /*
     * First step to count
     */
    int idx_elmt = 0;
    for(int i_elmt = 0; i_elmt < pn_entity[i_part]; ++i_elmt){
      for(int j = 0; j < part_group_stri[i_part][i_elmt]; j+=2){
        PDM_g_num_t i_group = part_group_data[i_part][idx_elmt];
        _group_idx[i_part][i_group+1] += 1;
        idx_elmt += 2;
      }
    }

    /*
     *  Deduce idx and allocate
     */
    for(int i_group = 0; i_group < n_group; ++i_group){
      _group_idx[i_part][i_group+1] += _group_idx[i_part][i_group];
    }
    int pgroup_tot_size = _group_idx[i_part][n_group];

    _group_ln_to_gn[i_part] = (PDM_g_num_t * ) malloc( pgroup_tot_size * sizeof(PDM_g_num_t));
    _group[i_part]          = (int         * ) malloc( pgroup_tot_size * sizeof(int        ));

    /*
     * Panic verbose
     */
    if( 0 == 1){
      printf("pgroup_tot_size::%d\n", pgroup_tot_size );
      printf("[%d][%d] _group_idx::\n", i_rank, i_part );
      for(int i_group = 0; i_group < n_group+1; ++i_group){
        printf("%d ", _group_idx[i_part][i_group]);
      }
      printf("\n");
    }

    /*
     * Fill data
     */
    idx_elmt = 0;
    for(int i_elmt = 0; i_elmt < pn_entity[i_part]; ++i_elmt){
      for(int j = 0; j < part_group_stri[i_part][i_elmt]; j+=2){
        PDM_g_num_t i_group        = part_group_data[i_part][idx_elmt  ];
        PDM_g_num_t elmt_position  = part_group_data[i_part][idx_elmt+1];

        int idx = _group_idx[i_part][i_group] + count_group[i_group];

        _group         [i_part][idx] = i_elmt+1;
        _group_ln_to_gn[i_part][idx] = elmt_position;

        count_group[i_group] += 1;
        idx_elmt += 2;
      }
    }
  }


  /*
   * Panic verbose
   */
  if( 0 == 1 ){
    for(int i_part = 0; i_part < n_part; ++i_part) {
      printf("[%i][%d] Groups \n", i_rank, i_part);
      for(int i_group = 0; i_group < n_group; ++i_group){
        printf("\t [%d]_group::", i_group);
        for(int idx = _group_idx[i_part][i_group]; idx < _group_idx[i_part][i_group+1]; ++idx){
          printf("%d ", _group[i_part][idx]);
        }
        printf("\n");
        printf("\t [%d]_group_ln_to_gn::", i_group);
        for(int idx = _group_idx[i_part][i_group]; idx < _group_idx[i_part][i_group+1]; ++idx){
          printf(PDM_FMT_G_NUM" ", _group_ln_to_gn[i_part][idx]);
        }
        printf("\n");
      }
    }
  }

  free(count_group);
  free(blk_stri_full);
  free(part_data);
  free(part_stri);
  free(blk_data);

  for(int i_part = 0; i_part < n_part; ++i_part) {
    free(part_group_stri[i_part]);
    free(part_group_data[i_part]);
  }
  free(part_group_stri);
  free(part_group_data);
  PDM_part_to_block_free(ptb_group);
  PDM_block_to_part_free(btp);
  for(int i_group = 0; i_group < n_group; ++i_group) {
    free(groups_distribution[i_group]);
  }
  free(groups_distribution);
}


/**
 *  \brief Computes partition connectivities with children gn from block connectivities.
 *
 * \param [in]   comm                PDM_MPI communicator
 * \param [in]   entity_distribution Distribution of entities over the processes (size=n_rank+1)
 * \param [in]   dconnectivity_idx   Distributed connectivity indexes (size=dn_entity+1)
 * \param [in]   dconnectivity       Distributed connectivity (size=dconnectivity_idx[dn_entity])
 * \param [in]   n_part              Number of partitions
 * \param [in]   pn_entity           Number of entities in each partition (size=n_part)
 * \param [in]   pentity_ln_to_gn    Array of local to global entity id for each partition
 * \param [out]  pconnectivity_idx   For each part, partitioned connectivity indexes
 *                                   (size = n_part, each component size = pn_entity[i_part])
 * \param [out]  pconnectivity_abs   For each part, partitioned connectivity in global numbering (size = n_part,
 *                                   each component size = pconnectivity_idx[i_part][pn_entity[i_part]])
*/
static
void
_dconnectivity_to_pconnectivity_abs
(
 const PDM_MPI_Comm    comm,
 const PDM_g_num_t    *entity_distribution,
 const int            *dconnectivity_idx,
 const PDM_g_num_t    *dconnectivity,
 const int             n_part,
 const int            *pn_entity,
 const PDM_g_num_t   **pentity_ln_to_gn,
       int          ***pconnectivity_idx,
       PDM_g_num_t  ***pconnectivity_abs
)
{
  int i_rank;
  int n_rank;

  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  int dn_entity = entity_distribution[i_rank+1] - entity_distribution[i_rank];

  /*
   * Prepare exchange protocol
   */
  PDM_block_to_part_t* btp = PDM_block_to_part_create(entity_distribution,
                               (const PDM_g_num_t **) pentity_ln_to_gn,
                                                      pn_entity,
                                                      n_part,
                                                      comm);

  /*
   * Prepare data
   */
  int* blk_stri = (int *) malloc( sizeof(int) * dn_entity);
  for(int i_elmt = 0; i_elmt < dn_entity; ++i_elmt){
    blk_stri[i_elmt] = dconnectivity_idx[i_elmt+1] - dconnectivity_idx[i_elmt];
  }

  /*
   * Exchange
   */
  int**         pstride;
  PDM_block_to_part_exch2(btp,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_VAR,
                          blk_stri,
             (void *  )   dconnectivity,
             (int  ***)  &pstride,
             (void ***)   pconnectivity_abs);

  free(blk_stri);

  /*
   * Panic verbose
   */
  if(0 == 1){
    for(int i_part = 0; i_part < n_part; ++i_part){
      int idx_data = 0;
      printf("[%d] pconnectivity_abs:: \n", i_part);
      for(int i_elmt = 0; i_elmt < pn_entity[i_part]; ++i_elmt) {
        printf("[%d] --> ", pstride[i_part][i_elmt]);
        for(int i_data = 0; i_data < pstride[i_part][i_elmt]; ++i_data ){
          printf(PDM_FMT_G_NUM" ", *pconnectivity_abs[i_part][idx_data++] );
        }
        printf("\n");
      }
      printf("\n");
    }
  }

  *pconnectivity_idx = (int         **) malloc( n_part * sizeof(int         *) );
  int** _pconnectivity_idx       = *pconnectivity_idx;

  for(int i_part = 0; i_part < n_part; ++i_part){
    _pconnectivity_idx[i_part] = PDM_array_new_idx_from_sizes_int(pstride[i_part], pn_entity[i_part]);
  }

  // free
  PDM_block_to_part_free(btp);
  for(int i_part = 0; i_part < n_part; ++i_part) {
    free(pstride[i_part]);
  }
  free(pstride);
}

/**
 *  \brief Computes a local numbering from connectivities with children in gn
 *
 * \param [in]  pn_child              For each part, n_child over connectivity (hence, non-unique)
 * \param [in]  pconnectivity_abs   For each part, partitioned connectivity in global numbering (size = n_part,
 *                                   each component size = pconnectivity_idx[i_part][pn_entity[i_part]])
 * \param [out]  pn_child_entity     Number of (unique) child elements in each partition (size=n_part)
 * \param [out]  pchild_ln_to_gn     For each part, position of child entity in the original array
 *                                   (size = n_part, each component size = pn_child_entity[i_part])
 * \param [out]  unique_order_ptr   For each part, gn_to_ln permutation
*/
static
void
_create_pchild_local_num
(
  const int      n_part,
  int           *pn_child,
  PDM_g_num_t  **pconnectivity_abs,
  int          **pn_child_entity,
  PDM_g_num_t ***pchild_ln_to_gn,
  int         ***unique_order_ptr
)
{
  *pn_child_entity   = (int          *) malloc( n_part * sizeof(int          ) );
  *pchild_ln_to_gn   = (PDM_g_num_t **) malloc( n_part * sizeof(PDM_g_num_t *) );
  int*  _pn_child_entity         = *pn_child_entity;
  PDM_g_num_t** _pchild_ln_to_gn = *pchild_ln_to_gn;
  *unique_order_ptr = (int**) malloc(n_part * sizeof(int*));
  int** unique_order = *unique_order_ptr;
  for(int i_part = 0; i_part < n_part; ++i_part){
    /*
     * Save array
     */
    int n_child = pn_child[i_part];
    _pchild_ln_to_gn[i_part] = (PDM_g_num_t *) malloc( n_child * sizeof(PDM_g_num_t));

    for(int i_child = 0; i_child < n_child; ++i_child) {
      _pchild_ln_to_gn[i_part][i_child] = PDM_ABS(pconnectivity_abs[i_part][i_child]);
    }

    /*
     * Deduce ln_to_gn
     */
    // printf("Sort data between : 0 and %d \n", n_child);
    unique_order[i_part] = (int *) malloc( n_child* sizeof(int));
    int n_elmt_sort = PDM_inplace_unique_long2(_pchild_ln_to_gn[i_part], unique_order[i_part], 0, n_child-1);
    _pn_child_entity[i_part] = n_elmt_sort;

    //if(0 == 1 && i_rank == 0){
    //  printf("n_elmt_sort::%d\n", n_elmt_sort);
    //  printf("_pchild_ln_to_gn::");
    //  for(int i = 0; i < n_elmt_sort; ++i){
    //    printf(PDM_FMT_G_NUM" ", _pchild_ln_to_gn[i_part][i]);
    //  }
    //  printf("\n");
    //}

    /*
     * Realloc
     */
    _pchild_ln_to_gn[i_part] = (PDM_g_num_t *) realloc(_pchild_ln_to_gn[i_part], n_elmt_sort * sizeof(PDM_g_num_t) );

  }
}

/**
 *  \brief Computes a local numbering from connectivities with children in gn
 *
 * \param [in]  pconnectivity_abs   For each part, partitioned connectivity in global numbering (size = n_part,
 *                                   each component size = pconnectivity_idx[i_part][pn_entity[i_part]])
 * \param [in]  unique_order        For each part, gn_to_ln permutation
 * \param [out] pconnectivity       For each part, partitioned connectivity (size = n_part,
 *                                   each component size = pconnectivity_idx[i_part][pn_entity[i_part]])
*/
static
void
_pconnectivity_with_local_num
(
  const int      n_part,
  int           *pn_child,
  PDM_g_num_t  **pconnectivity_abs,
  int          **unique_order,
  int         ***pconnectivity
)
{
  *pconnectivity = (int**) malloc( n_part * sizeof(int*) );
  int** _pconnectivity = *pconnectivity;
  for(int i_part = 0; i_part < n_part; ++i_part){
    /*
     *  We need to regenerate the connectivity and pass it in local numbering
     *  This one is vectorisable if we remove PDM_binary_seach and idx_data++
     */
    // idx_data = 0;
    // for(int i_cell = 0; i_cell < n_elmts[i_part]; ++i_cell) {
    //   for(int i_data = 0; i_data < cell_stri[i_part][i_cell]; ++i_data ){
    //     int         g_sgn  = PDM_SIGN(pcell_face_tmp[i_part][idx_data]);
    //     PDM_g_num_t g_elmt = PDM_ABS (pcell_face_tmp[i_part][idx_data]);
    //     int l_elmt         = unique_order[idx_data];
    //     int l_elmt_old     = PDM_binary_search_long(g_elmt, _pchild_ln_to_gn[i_part], n_elmt_sort); /* In [0, n_elmt_sort-1] */

    //     // printf("[%d] - Search [%d] --> %d | %d \n", idx_data, g_elmt, l_elmt, l_elmt_old);

    //     /* Overwrite the pcell_face with local numbering and reput sign on it */
    //     _pconnectivity[i_part][idx_data++] = (l_elmt + 1) * g_sgn ;
    //   }
    // }

    int n_child = pn_child[i_part];
    _pconnectivity[i_part]   = (int*) malloc( n_child * sizeof(int        ));
    for(int idx = 0; idx < n_child; ++idx) {
      int g_sgn  = PDM_SIGN(pconnectivity_abs[i_part][idx]);
      int l_elmt = unique_order[i_part][idx];
      _pconnectivity[i_part][idx] = (l_elmt + 1) * g_sgn ;
    }
  }
}


void
PDM_part_multi_dconnectivity_to_pconnectivity_sort
(
 const PDM_MPI_Comm    comm,
 const int             n_part,
 const int             n_section,
 const int            *section_idx,
       PDM_g_num_t   **entity_distribution,
       int            *dconnectivity_idx,
       PDM_g_num_t    *dconnectivity,
       int           **pn_entity,
       PDM_g_num_t  ***pentity_ln_to_gn,
       int           **pn_child_entity,
       PDM_g_num_t  ***pchild_ln_to_gn,
       int         ****pconnectivity_idx,
       int         ****pconnectivity
)
{
  //printf("PDM_part_multi_dconnectivity_to_pconnectivity_sort\n");
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  *pconnectivity_idx = (int***) malloc(n_section * sizeof(int**));
  int* pn_vtx = (int*)malloc(n_part * sizeof(int));
  for (int i=0; i<n_part; ++i) pn_vtx[i]=0;
  PDM_g_num_t*** pconnectivity_abs = (PDM_g_num_t***) malloc(n_section * sizeof(PDM_g_num_t**));

  for (int i_section=0; i_section<n_section; ++i_section) {
    int* dconnectivity_section_idx = dconnectivity_idx + section_idx[i_section];
    int dn_entity = entity_distribution[i_section][i_rank+1] - entity_distribution[i_section][i_rank];
    int* dconnectivity_section_idx_at_0 = (int*)malloc((dn_entity+1) * sizeof(int));
    for (int i=0; i<dn_entity+1; ++i) {
      dconnectivity_section_idx_at_0[i] = dconnectivity_section_idx[i] - dconnectivity_section_idx[0];
    }
    PDM_g_num_t* dconnectivity_section = dconnectivity + dconnectivity_section_idx[0];

    // 0. create pconnectivity with global numbering
    _dconnectivity_to_pconnectivity_abs(
      comm,
      entity_distribution[i_section],
      dconnectivity_section_idx_at_0,
      dconnectivity_section,
      n_part,
      pn_entity[i_section],
      (const PDM_g_num_t**) pentity_ln_to_gn[i_section],
      &(*pconnectivity_idx)[i_section],&pconnectivity_abs[i_section]
    );
    free(dconnectivity_section_idx_at_0);

    for(int i_part = 0; i_part < n_part; ++i_part) {
      int n_elmts = pn_entity[i_section][i_part];
      pn_vtx[i_part] += (*pconnectivity_idx)[i_section][i_part][n_elmts];
    }
  }

  // 1. Create local numbering
  // Caution the recv connectivity can be negative
  PDM_g_num_t** pconnectivity_abs_cat = (PDM_g_num_t**)malloc( n_part *sizeof(PDM_g_num_t*));
  for (int i_part=0; i_part<n_part; ++i_part) {
    pconnectivity_abs_cat[i_part] = (PDM_g_num_t*)malloc( pn_vtx[i_part] *sizeof(PDM_g_num_t));
    int pos = 0;
    for (int i_section=0; i_section<n_section; ++i_section) {
      int n_elt = pn_entity[i_section][i_part];
      int n_vtx_section = (*pconnectivity_idx)[i_section][i_part][n_elt];
      for (int j=0; j<n_vtx_section; ++j) {
        pconnectivity_abs_cat[i_part][pos] = pconnectivity_abs[i_section][i_part][j];
        ++pos;
      }
    }
  }
  int** unique_order;
  _create_pchild_local_num(
    n_part,
    pn_vtx,
    pconnectivity_abs_cat,
    pn_child_entity,
    pchild_ln_to_gn,
    &unique_order
  );
  // 2. create pconnectivity with local numbering
  int** pconnectivity_cat = (int**)malloc(n_part * sizeof(int*));
  for (int i_part=0; i_part<n_part; ++i_part) {
    pconnectivity_cat[i_part] = (int*)malloc(pn_vtx[i_part] * sizeof(int));
  }
  _pconnectivity_with_local_num(
    n_part,
    pn_vtx,
    pconnectivity_abs_cat,
    unique_order,
    &pconnectivity_cat
  );
  int pos = 0;
  *pconnectivity = (int***)malloc(n_section * sizeof(int**));
  for (int i_section=0; i_section<n_section; ++i_section) {
    (*pconnectivity)[i_section] = (int**)malloc(n_part * sizeof(int*));
    assert(n_part==1); // TODOUX
    for(int i_part = 0; i_part < n_part; ++i_part) {
      int n_elmts = pn_entity[i_section][i_part];
      int n_vtx_section = (*pconnectivity_idx)[i_section][i_part][n_elmts];

      (*pconnectivity)[i_section][i_part] = (int*)malloc(n_vtx_section * sizeof(int));
      for (int i=0; i<n_vtx_section; ++i) {
        (*pconnectivity)[i_section][i_part][i] = pconnectivity_cat[i_part][pos];
        ++pos;
      }
    }
  }

  ///*
  // * Free
  // */
  for(int i_section = 0; i_section < n_section; ++i_section) {
    for(int i_part = 0; i_part < n_part; ++i_part) {
      free(pconnectivity_abs[i_section][i_part]);
    }
    free(pconnectivity_abs[i_section]);
  }
  free(pconnectivity_abs);

  for(int i_part = 0; i_part < n_part; ++i_part) {
    free(pconnectivity_abs_cat[i_part]);
    free(pconnectivity_cat[i_part]);
  }
  free(pconnectivity_abs_cat);
  free(pconnectivity_cat);

  free(pn_vtx);

  for(int i_part = 0; i_part < n_part; ++i_part) {
    free(unique_order[i_part]);
  }
  free(unique_order);
}

/**
 *  \brief Generated the partitioned connectivity (entity->child_elements) associated
 *   to the given distributed connectivity, using element distribution and element local
 *   to global numbering. In addition, return the partitioned number of unique child_element
 *   and the corresponding local to global numbering for the child elements.
 *
 *   The orientation data (ie negative index) present in the distributed connectivity,
 *   if any, are preserved, *meaning that boundary faces can be badly oriented on partitions*.
 *   See PDM_part_reorient_bound_faces function to correct this orientation.
 *
 * \param [in]   comm                PDM_MPI communicator
 * \param [in]   entity_distribution Distribution of entities over the processes (size=n_rank+1)
 * \param [in]   dconnectivity_idx   Distributed connectivity indexes (size=dn_entity+1)
 * \param [in]   dconnectivity       Distributed connectivity (size=dconnectivity_idx[dn_entity])
 * \param [in]   n_part              Number of partitions
 * \param [in]   pn_entity           Number of entities in each partition (size=n_part)
 * \param [in]   pentity_ln_to_gn    Array of local to global entity id for each partition
 *                                   (size=n_part, each component size = pn_entities[i_part])
 * \param [out]  pn_child_entity     Number of (unique) child elements in each partition (size=n_part)
 * \param [out]  pchild_ln_to_gn     For each part, position of child entity in the original array
 *                                   (size = n_part, each component size = pn_child_entity[i_part])
 * \param [out]  pconnectivity_idx   For each part, partitioned connectivity indexes
 *                                   (size = n_part, each component size = pn_entity[i_part])
 * \param [out]  pconnectivity       For each part, partitioned connectivity (size = n_part,
 *                                   each component size = pconnectivity_idx[i_part][pn_entity[i_part]])
*/
void
PDM_part_dconnectivity_to_pconnectivity_sort
(
 const PDM_MPI_Comm    comm,
 const PDM_g_num_t    *entity_distribution,
 const int            *dconnectivity_idx,
 const PDM_g_num_t    *dconnectivity,
 const int             n_part,
 const int            *pn_entity,
 const PDM_g_num_t   **pentity_ln_to_gn,
       int           **pn_child_entity,
       PDM_g_num_t  ***pchild_ln_to_gn,
       int          ***pconnectivity_idx,
       int          ***pconnectivity
)
{
  //printf("PDM_part_dconnectivity_to_pconnectivity\n");
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  // 0. create pconnectivity with global numbering
  PDM_g_num_t** pconnectivity_abs;
  _dconnectivity_to_pconnectivity_abs(
    comm,
    entity_distribution,
    dconnectivity_idx,
    dconnectivity,
    n_part,
    pn_entity,
    pentity_ln_to_gn,
    pconnectivity_idx,&pconnectivity_abs
  );
  int** _pconnectivity_idx       = *pconnectivity_idx;

  // 1. Create local numbering
  // Caution the recv connectivity can be negative
  int** unique_order;
  int* pn_child = (int*)malloc(n_part * sizeof(int));
  for(int i_part = 0; i_part < n_part; ++i_part) {
    int n_elmts = pn_entity[i_part];
    pn_child[i_part] = _pconnectivity_idx[i_part][n_elmts];
  }
  _create_pchild_local_num(
    n_part,
    pn_child,
    pconnectivity_abs,
    pn_child_entity,
    pchild_ln_to_gn,
    &unique_order
  );

  // 2. create pconnectivity with local numbering
  _pconnectivity_with_local_num(
    n_part,
    pn_child,
    pconnectivity_abs,
    unique_order,
    pconnectivity
  );

  /*
   * Panic verbose
   */
  if(0 == 1){
    for(int i_part = 0; i_part < n_part; ++i_part){
      int idx_data = 0;
      printf("[%d] pconnectivity:: \n", i_part);
      for(int i_cell = 0; i_cell < pn_entity[i_part]; ++i_cell) {
        int stride = _pconnectivity_idx[i_part][i_cell+1] - _pconnectivity_idx[i_part][i_cell];
        printf("[%d] --> ", stride);
        for(int i_data = 0; i_data < stride; ++i_data ){
          printf("%d ", (*pconnectivity)[i_part][idx_data++] );
        }
        printf("\n");
      }
      printf("\n");
    }
  }

  /*
   * Free
   */
  for(int i_part = 0; i_part < n_part; ++i_part) {
    free(pconnectivity_abs[i_part]);
    free(unique_order[i_part]);
  }
  free(pconnectivity_abs);
  free(unique_order);
  free(pn_child);
}


void
PDM_part_dconnectivity_to_pconnectivity_sort_single_part
(
 const PDM_MPI_Comm    comm,
 const PDM_g_num_t    *entity_distribution,
 const int            *dconnectivity_idx,
 const PDM_g_num_t    *dconnectivity,
 const int             pn_entity,
 const PDM_g_num_t    *pentity_ln_to_gn,
       int            *pn_child_entity,
       PDM_g_num_t   **pchild_ln_to_gn,
       int           **pconnectivity_idx,
       int           **pconnectivity
)
{
  int          *tmp_pn_child_entity   = NULL;
  PDM_g_num_t **tmp_pchild_ln_to_gn   = NULL;
  int         **tmp_pconnectivity_idx = NULL;
  int         **tmp_pconnectivity     = NULL;
  PDM_part_dconnectivity_to_pconnectivity_sort(comm,
                                               entity_distribution,
                                               dconnectivity_idx,
                                               dconnectivity,
                                               1,
                                               &pn_entity,
                                               &pentity_ln_to_gn,
                                               &tmp_pn_child_entity,
                                               &tmp_pchild_ln_to_gn,
                                               &tmp_pconnectivity_idx,
                                               &tmp_pconnectivity);


  *pn_child_entity   = tmp_pn_child_entity[0];
  *pchild_ln_to_gn   = tmp_pchild_ln_to_gn[0];
  *pconnectivity_idx = tmp_pconnectivity_idx[0];
  *pconnectivity     = tmp_pconnectivity[0];

  free(tmp_pn_child_entity);
  free(tmp_pchild_ln_to_gn);
  free(tmp_pconnectivity_idx);
  free(tmp_pconnectivity);
}

/**
 *  \brief Generated the partitioned connectivity (entity->child_elements) associated
 *   to the given distributed connectivity, using element distribution and element local
 *   to global numbering. In addition, return the partitioned number of unique child_element
 *   and the corresponding local to global numbering for the child elements.
 *
 *   The orientation data (ie negative index) present in the distributed connectivity,
 *   if any, are preserved, *meaning that boundary faces can be badly oriented on partitions*.
 *   See PDM_part_reorient_bound_faces function to correct this orientation.
 *
 * \param [in]   comm                PDM_MPI communicator
 * \param [in]   entity_distribution Distribution of entities over the processes (size=n_rank+1)
 * \param [in]   dconnectivity_idx   Distributed connectivity indexes (size=dn_entity+1)
 * \param [in]   dconnectivity       Distributed connectivity (size=dconnectivity_idx[dn_entity])
 * \param [in]   n_part              Number of partitions
 * \param [in]   pn_entity           Number of entities in each partition (size=n_part)
 * \param [in]   pentity_ln_to_gn    Array of local to global entity id for each partition
 *                                   (size=n_part, each component size = pn_entities[i_part])
 * \param [out]  pn_child_entity     Number of (unique) child elements in each partition (size=n_part)
 * \param [out]  pchild_ln_to_gn     For each part, position of child entity in the original array
 *                                   (size = n_part, each component size = pn_child_entity[i_part])
 * \param [out]  pconnectivity_idx   For each part, partitioned connectivity indexes
 *                                   (size = n_part, each component size = pn_entity[i_part])
 * \param [out]  pconnectivity       For each part, partitioned connectivity (size = n_part,
 *                                   each component size = pconnectivity_idx[i_part][pn_entity[i_part]])
*/
void
PDM_part_dconnectivity_to_pconnectivity_hash
(
 const PDM_MPI_Comm    comm,
 const PDM_g_num_t    *entity_distribution,
 const int            *dconnectivity_idx,
 const PDM_g_num_t    *dconnectivity,
 const int             n_part,
 const int            *pn_entity,
 const PDM_g_num_t   **pentity_ln_to_gn,
       int           **pn_child_entity,
       PDM_g_num_t  ***pchild_ln_to_gn,
       int          ***pconnectivity_idx,
       int          ***pconnectivity
)
{
  printf("PDM_part_dconnectivity_to_pconnectivity\n");
  int i_rank;
  int n_rank;

  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  int dn_entity = entity_distribution[i_rank+1] - entity_distribution[i_rank];

  /*
   * Prepare exchange protocol
   */
  PDM_block_to_part_t* btp = PDM_block_to_part_create(entity_distribution,
                               (const PDM_g_num_t **) pentity_ln_to_gn,
                                                      pn_entity,
                                                      n_part,
                                                      comm);

  /*
   * Prepare data
   */
  int* blk_stri = (int *) malloc( sizeof(int) * dn_entity);
  for(int i_elmt = 0; i_elmt < dn_entity; ++i_elmt){
    blk_stri[i_elmt] = dconnectivity_idx[i_elmt+1] - dconnectivity_idx[i_elmt];
  }

  /*
   * Exchange
   */
  int**         pstride;
  PDM_g_num_t** pconnectivity_tmp; /* We keep it in double precision because it contains global numbering */
  PDM_block_to_part_exch2(btp,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_VAR,
                          blk_stri,
             (void *  )   dconnectivity,
             (int  ***)  &pstride,
             (void ***)  &pconnectivity_tmp);

  free(blk_stri);

  /*
   * Panic verbose
   */
  if(0 == 1){
    for(int i_part = 0; i_part < n_part; ++i_part){
      int idx_data = 0;
      printf("[%d] pconnectivity_tmp:: \n", i_part);
      for(int i_elmt = 0; i_elmt < pn_entity[i_part]; ++i_elmt) {
        printf("[%d] --> ", pstride[i_part][i_elmt]);
        for(int i_data = 0; i_data < pstride[i_part][i_elmt]; ++i_data ){
          printf(PDM_FMT_G_NUM" ", pconnectivity_tmp[i_part][idx_data++] );
        }
        printf("\n");
      }
      printf("\n");
    }
  }


  /*
   * Post-treatment - Caution the recv connectivity can be negative
   */
  *pchild_ln_to_gn   = (PDM_g_num_t **) malloc( n_part * sizeof(PDM_g_num_t *) );
  *pconnectivity     = (int         **) malloc( n_part * sizeof(int         *) );
  *pconnectivity_idx = (int         **) malloc( n_part * sizeof(int         *) );
  *pn_child_entity   = (int          *) malloc( n_part * sizeof(int          ) );

  /* Shortcut */
  PDM_g_num_t** _pchild_ln_to_gn = *pchild_ln_to_gn;
  int** _pconnectivity           = *pconnectivity;
  int** _pconnectivity_idx       = *pconnectivity_idx;
  int*  _pn_child_entity         = *pn_child_entity;

  for(int i_part = 0; i_part < n_part; ++i_part){

    int n_elmts = pn_entity[i_part];
    /*
     *  First loop to count and setup part_strid_idx
     */
    _pconnectivity_idx[i_part] = PDM_array_new_idx_from_sizes_int(pstride[i_part], pn_entity[i_part]);



    /*
     * Save array
     */
    _pchild_ln_to_gn[i_part] = (PDM_g_num_t *) malloc( _pconnectivity_idx[i_part][n_elmts] * sizeof(PDM_g_num_t));
    _pconnectivity[i_part]   = (int *        ) malloc( _pconnectivity_idx[i_part][n_elmts] * sizeof(int        ));

    PDM_g_num_t *_keys_and_values = (PDM_g_num_t *) malloc( 2*_pconnectivity_idx[i_part][n_elmts] * sizeof(PDM_g_num_t));

    /* Create hash table */
    int keymax = (int) (1.3 * _pconnectivity_idx[i_part][n_elmts]); //next prime number following 1.3*nbofvalues
    while (!_is_prime(keymax))
      keymax++;
    //PDM_printf("[%i] nb val is %d, choose keymax = %d\n", i_rank, _pconnectivity_idx[i_part][n_elmts] ,keymax);
    PDM_hash_tab_t *hash_tab = PDM_hash_tab_create(PDM_HASH_TAB_KEY_INT, &keymax);

    /* Fill & use table */
    int idx_data = 0;
    int nbUnique = 0;
    for(int i_elmt = 0; i_elmt < n_elmts; ++i_elmt) {
      for(int i_data = 0; i_data < pstride[i_part][i_elmt]; ++i_data ){

        PDM_g_num_t key_value =  PDM_ABS(pconnectivity_tmp[i_part][idx_data]);
        int            g_sgn  = PDM_SIGN(pconnectivity_tmp[i_part][idx_data]);
        int hash_value = key_value % keymax;

        /* We will store in the tab the key (to maange collisions) and the position in
           lntogn list (to construct partioned connectivity) */
        PDM_g_num_t *_key_and_value = _keys_and_values + 2*idx_data;
        _key_and_value[0]   = key_value;
        _key_and_value[1]   = (PDM_g_num_t) nbUnique;

        int key_found_in_table = 0;
        int n_in_tab = PDM_hash_tab_n_data_get(hash_tab, &hash_value);

        /* If the slot exists in the tab (n_in_tab), check if the id is already registered,
           or if it is a hash collision. */
        if (n_in_tab > 0) {
          //if (i_rank == 0) PDM_printf("Check uniqueness/collision before adding %d at pos %d\n", key_value, hash_value);
          PDM_g_num_t **candidate_in_tab = (PDM_g_num_t **) PDM_hash_tab_data_get(hash_tab, &hash_value);
          for (int j = 0; j < n_in_tab; j++) {
            if (PDM_ABS(candidate_in_tab[j][0]) == key_value) {
              _pconnectivity[i_part][idx_data] = g_sgn*(candidate_in_tab[j][1] + 1);
              key_found_in_table++; //The key was trully in table
              //if (i_rank == 0) PDM_printf("  candidate key %d matches", candidate_in_tab[j][0]);
              break;
            }
            //else if (i_rank == 0) PDM_printf("  candidate key %d is a collision", candidate_in_tab[j][0]);
          }
          //if (i_rank == 0) PDM_printf("\n");
        }
        /* Slot is empty *or* we just have a collision -> add key in table */
        if (!key_found_in_table) {
          //if (i_rank == 0) PDM_printf("Add %d in table at pos %d\n", key_value, hash_value);
          PDM_hash_tab_data_add(hash_tab, (void *) &hash_value, _key_and_value);
          _pconnectivity[i_part][idx_data] = g_sgn*(nbUnique + 1);
          _pchild_ln_to_gn[i_part][nbUnique++] = PDM_ABS(pconnectivity_tmp[i_part][idx_data]);
        }
        idx_data++;
      }
    }

    PDM_hash_tab_free(hash_tab);
    free(_keys_and_values);

    _pn_child_entity[i_part] = nbUnique;

    if(0 == 1){
      //printf("n_elmt_sort::%d\n", n_elmt_sort);
      //printf("_pchild_ln_to_gn::");
      //for(int i = 0; i < nbUnique; ++i){
      //  printf(PDM_FMT_G_NUM" ", _pchild_ln_to_gn[i_part][i]);
      //}
      //printf("\n");
    }

    /*
     * Realloc
     */
    _pchild_ln_to_gn[i_part] = (PDM_g_num_t *) realloc(_pchild_ln_to_gn[i_part], nbUnique * sizeof(PDM_g_num_t) );
  }


  /*
   * Panic verbose
   */
  if(0 == 1){
    for(int i_part = 0; i_part < n_part; ++i_part){
      int idx_data = 0;
      printf("[%d] pconnectivity:: \n", i_part);
      for(int i_cell = 0; i_cell < pn_entity[i_part]; ++i_cell) {
        printf("[%d] --> ", pstride[i_part][i_cell]);
        for(int i_data = 0; i_data < pstride[i_part][i_cell]; ++i_data ){
          printf("%d ", (*pconnectivity)[i_part][idx_data++] );
        }
        printf("\n");
      }
      printf("\n");
    }
  }

  /*
   * Free
   */
  PDM_block_to_part_free(btp);
  for(int i_part = 0; i_part < n_part; ++i_part) {
    free(pstride[i_part]);
    free(pconnectivity_tmp[i_part]);
  }
  free(pstride);
  free(pconnectivity_tmp);
}


static int
min_sub_index
(
  int* array,
  int* sub_indices,
  int n_sub_indices,
  int i_rank
)
{
  PDM_UNUSED(i_rank);
  int min_idx = sub_indices[0];
  for (int i=1; i<n_sub_indices; ++i) {
    if (array[sub_indices[i]] < array[min_idx]) {
      min_idx = sub_indices[i];
    }
  }
  return min_idx;
}

/**
 *  \brief Generates the communication information at the partition interfaces for the
 *   given entity. The communication data associates to
 *   each partitioned entity belonging to an (internal) interface the 4-tuple
 *   (local id, opposite proc number, opposite part number on opposite proc, local id in the
 *   opposite partition).
 *   This list is sorted by opposite proc id, then by part id, and finally with respect
 *   to the entity global_id. Also return the stride indices pproc_bound_idx and ppart_bound_idx
 *   to acces the communication information for a given opposite proc id or global part id.
 *
 * \param [in]   comm                PDM_MPI communicator
 * \param [in]   part_distribution   Distribution of partitions over the processes (size=n_rank+1)
 * \param [in]   entity_distribution Distribution of entities over the processes (size=n_rank+1)
 * \param [in]   n_part              Number of partitions
 * \param [in]   pn_entity           Number of entities in each partition (size=n_part)
 * \param [in]   pentity_ln_to_gn    Array of local to global entity id for each partition
 *                                   (size=n_part, each component size = pn_entity[i_part])
 * \param [in]   pentity_hint        Can be used to indicate whether (1) or not (0) an entity is potentially
 *                                   shared with an other partition in order to minimize exchanged data
 *                                   (size=n_part, each component size = pn_entity[i_part]) or NULL
 * \param [out]  pproc_bound_idx     For each part, indexes of communication information related to the
 *                                   other procs (size=n_part, each component size=n_rank+1)
 * \param [out]  ppart_bound_idx     For each part, indexes of communication information related to the
 *                                   other (global id) parts (size=n_part, each component size=n_part_tot+1)
 * \param [out]  pentity_bound       For each part, communication information (see abobe) (size=n_part)
 */
void
PDM_part_generate_entity_graph_comm
(
 const PDM_MPI_Comm   comm,
 const PDM_g_num_t   *part_distribution,
 const PDM_g_num_t   *entity_distribution,
 const int            n_part,
 const int           *pn_entity,
 const PDM_g_num_t  **pentity_ln_to_gn,
 const int          **pentity_hint,
       int         ***pproc_bound_idx,
       int         ***ppart_bound_idx,
       int         ***pentity_bound,
       int         ***pentity_priority
)
{
  int i_rank;
  int n_rank;

  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  // printf(" PDM_part_assemble_partitions PDM_part_generate_entity_graph_comm flag 1 END\n");
  // PDM_MPI_Barrier(comm);

  int setup_priority = 0; // False
  if(pentity_priority == NULL){
    // printf(" pentity_priority not defined \n");
  } else {
    setup_priority = 1;
    assert(pentity_hint == NULL);
    // printf(" pentity_priority : %p \n", pentity_priority);
  }

  /*
   * First part : we put in data the triplet (iRank, iPart, iEntityLoc)
   */
  int** part_stri = (int ** ) malloc( n_part * sizeof(int *));
  int** part_data = (int ** ) malloc( n_part * sizeof(int *));
  for(int i_part = 0; i_part < n_part; ++i_part) {
    part_stri[i_part] = (int *) malloc( pn_entity[i_part] * sizeof(int));
    int idx_data = 0;
    if (pentity_hint != NULL) {
      //First pass to count
      for(int i_entity = 0; i_entity < pn_entity[i_part]; ++i_entity){
        if (pentity_hint[i_part][i_entity] == 1){
          idx_data++;
        }
      }
    } else {
      idx_data = pn_entity[i_part];
    }
    //Second pass to fill
    part_data[i_part] = (int *) malloc( 3 * idx_data * sizeof(int));
    if (pentity_hint != NULL){
      idx_data = 0;
      for(int i_entity = 0; i_entity < pn_entity[i_part]; ++i_entity) {
        if (pentity_hint[i_part][i_entity]==1){
          part_data[i_part][3*idx_data  ] = i_rank;
          part_data[i_part][3*idx_data+1] = i_part;
          part_data[i_part][3*idx_data+2] = i_entity;

          part_stri[i_part][i_entity] = 3;
          idx_data++;
        } else{
          part_stri[i_part][i_entity] = 0;
        }
      }
    } else {
      for(int i_entity = 0; i_entity < pn_entity[i_part]; ++i_entity) {
        part_data[i_part][3*i_entity  ] = i_rank;
        part_data[i_part][3*i_entity+1] = i_part;
        part_data[i_part][3*i_entity+2] = i_entity;

        part_stri[i_part][i_entity] = 3;
      }
    }
  }

  // printf(" PDM_part_assemble_partitions PDM_part_generate_entity_graph_comm flag 2 \n");
  // PDM_MPI_Barrier(comm);
  /*
   * Setup protocol exchange
   */
  PDM_part_to_block_t *ptb =
   PDM_part_to_block_create2(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                             PDM_PART_TO_BLOCK_POST_MERGE,
                             1.,
            (PDM_g_num_t **) pentity_ln_to_gn,
                             entity_distribution,
                             (int *) pn_entity,
                             n_part,
                             comm);
   const int n_entity_block = PDM_part_to_block_n_elt_block_get(ptb);

   /*
    * Exchange
    */
  int* blk_stri = NULL;
  int* blk_data = NULL;
  PDM_part_to_block_exch (ptb,
                          sizeof(int),
                          PDM_STRIDE_VAR,
                          1,
                          part_stri,
                (void **) part_data,
                          &blk_stri,
                (void **) &blk_data);
  // get the procs from which the data comes from
  int** proc_part_stri = (int ** ) malloc( n_part * sizeof(int *));
  int** proc_part_data = (int ** ) malloc( n_part * sizeof(int *));
  for(int i_part = 0; i_part < n_part; ++i_part) {
    proc_part_stri[i_part] = (int *) malloc( pn_entity[i_part] * sizeof(int));
    proc_part_data[i_part] = (int *) malloc( pn_entity[i_part] * sizeof(int));
    for (int i=0; i<pn_entity[i_part]; ++i) {
      proc_part_stri[i_part][i] = 1;
      proc_part_data[i_part][i] = i_rank;
    }
  }
  int* proc_blk_stri = NULL;
  int* proc_blk_data = NULL;
  PDM_part_to_block_exch (ptb,
                          sizeof(int),
                          PDM_STRIDE_VAR,
                          1,
                          proc_part_stri,
                (void **) proc_part_data,
                          &proc_blk_stri,
                (void **) &proc_blk_data);

  PDM_g_num_t* new_distrib = PDM_part_to_block_adapt_partial_block_to_block(ptb, &blk_stri, entity_distribution[n_rank]);
  free(new_distrib);
  new_distrib = PDM_part_to_block_adapt_partial_block_to_block(ptb, &proc_blk_stri, entity_distribution[n_rank]);
  // int dn_check = new_distrib[i_rank+1] - new_distrib[i_rank];
  // int dn_entity = entity_distribution[i_rank+1] - entity_distribution[i_rank];
  free(new_distrib);

  /*
   * Free
   */
  PDM_part_to_block_free(ptb);
  for(int i_part = 0; i_part < n_part; ++i_part) {
    free(part_stri[i_part]);
    free(part_data[i_part]);
    free(proc_part_stri[i_part]);
    free(proc_part_data[i_part]);
  }
  free(proc_part_stri);
  free(proc_part_data);
  free(part_stri);
  free(part_data);

  /*
   * Panic verbose
   */
  if(0 == 1){
    int idx_debug_data = 0;
    for(int i_block = 0; i_block < n_entity_block; ++i_block){
      printf("[%d]-blk_stri[%d]:: ", i_block, blk_stri[i_block]);
      for(int i_data = 0; i_data < blk_stri[i_block]; ++i_data){
        printf("%i ", blk_data[idx_debug_data++]);
      }
      printf("\n");
    }

  }

  // printf(" PDM_part_assemble_partitions PDM_part_generate_entity_graph_comm flag 3 \n");
  // PDM_MPI_Barrier(comm);
  /*
   * Post-treatment : For all non shared data we have blk_stri == 3
   *                  And for others we have n_shared * 3 data
   *                  In blk view we can easily remove all non shared data
   *                  We reexchange with block_to_part only the shared data
   * If priority is query an equilbrate choice is made up
   */
  int* blk_priority_data;
  if(setup_priority == 1){
    blk_priority_data = (int * ) malloc( n_entity_block * sizeof(int));
  }

  int idx_comp = 0;     /* Compressed index use to fill the buffer */
  int idx_data = 0;     /* Index in the block to post-treat        */

  int* n_owner_entity_by_rank = PDM_array_zeros_int(n_rank);
  int* proc_data_current = proc_blk_data;

  for(int i_block = 0; i_block < n_entity_block; ++i_block){

    /* Non shared data --> Compression */
    if(blk_stri[i_block] == 3){
      blk_stri[i_block] = 0;
      idx_data += 3;          /* Mv directly to the next */
      if(setup_priority == 1){
        blk_priority_data[i_block] = -1;
      }

    } else {
      for(int i_data = 0; i_data < blk_stri[i_block]; ++i_data){
        blk_data[idx_comp++] = blk_data[idx_data++];
      }
      /* Fill up entity priority */
      if(setup_priority == 1){
        int selected_owner_rank = min_sub_index(n_owner_entity_by_rank,proc_data_current,proc_blk_stri[i_block],i_rank);
        ++n_owner_entity_by_rank[selected_owner_rank];

        blk_priority_data[i_block] = selected_owner_rank;
      }
    }
    proc_data_current += proc_blk_stri[i_block];
  }
  free(n_owner_entity_by_rank);

  /*
   * Compress data
   */
  blk_data = (int *) realloc(blk_data, idx_comp * sizeof(int));

  /*
   * Panic verbose
   */
  if(0 == 1){
    int idx_debug_data = 0;
    for(int i_block = 0; i_block < n_entity_block; ++i_block){
      printf("[%d]-blk_stri[%d]:: ", i_block, blk_stri[i_block]);
      for(int i_data = 0; i_data < blk_stri[i_block]; ++i_data){
        printf("%i ", blk_data[idx_debug_data++]);
      }
      printf("\n");
    }
  }

  // int n_stride = 0;
  // for(int i_block = 0; i_block < n_entity_block; ++i_block){
  //   n_stride += blk_stri[i_block];
  // }


  // log_trace("n_stride = %i | n_entity_block = %i | dn_check = %i \n", n_stride, n_entity_block, dn_check);
  // if(dn_check != n_entity_block) {
  //   log_trace("STRANGE : n_stride = %i | n_entity_block = %i | dn_check = %i | dn_entity = %i \n", n_stride, n_entity_block, dn_check, dn_entity);
  // }

  // printf(" PDM_part_assemble_partitions PDM_part_generate_entity_graph_comm flag 4 \n");
  // PDM_MPI_Barrier(comm);
  /*
   * All data is now sort we cen resend to partition
   */
  PDM_block_to_part_t* btp = PDM_block_to_part_create(entity_distribution,
                               (const PDM_g_num_t **) pentity_ln_to_gn,
                                                      pn_entity,
                                                      n_part,
                                                      comm);

  // printf(" PDM_part_assemble_partitions PDM_part_generate_entity_graph_comm flag 4 - 1 \n");
  // PDM_MPI_Barrier(comm);
  PDM_block_to_part_exch2(btp,
                          sizeof(int),
                          PDM_STRIDE_VAR,
                          blk_stri,
             (void *  )   blk_data,
             (int  ***)  &part_stri,
             (void ***)  &part_data);
  // printf(" PDM_part_assemble_partitions PDM_part_generate_entity_graph_comm flag 4 - 2 \n");
  // PDM_MPI_Barrier(comm);

  if(setup_priority == 1 ){
    int stride_one = 1;
    PDM_block_to_part_exch2(btp,
                            sizeof(int),
                            PDM_STRIDE_CST,
                            &stride_one,
               (void *  )   blk_priority_data,
                            NULL,
               (void ***)   pentity_priority);
    free(blk_priority_data);
  }
  // printf(" PDM_part_assemble_partitions PDM_part_generate_entity_graph_comm flag 4 - 3 \n");
  // PDM_MPI_Barrier(comm);

  /*
   * Free
   */
  free(blk_data);
  free(blk_stri);
  free(proc_blk_data);
  free(proc_blk_stri);

  /*
   * Interface for pdm_part is  :
   *    - Face local number
   *    - Connected process
   *    - Connected partition on the connected process
   *    - Connected face local number in the connected partition
   */


  // printf(" PDM_part_assemble_partitions PDM_part_generate_entity_graph_comm flag 5 \n");
  // PDM_MPI_Barrier(comm);
  /* Allocate */
  *pproc_bound_idx   = (int **) malloc( ( n_part ) * sizeof(int * ));
  *ppart_bound_idx   = (int **) malloc( ( n_part ) * sizeof(int * ));
  *pentity_bound     = (int **) malloc( ( n_part ) * sizeof(int * ));

  /* Shortcut */
  int** _pproc_bound_idx   = *pproc_bound_idx;
  int** _ppart_bound_idx   = *ppart_bound_idx;
  int** _pentity_bound     = *pentity_bound;

  /*
   * Post-treatment in partition
   */
  for(int i_part = 0; i_part < n_part; ++i_part){

    /* Shortcut */
    int* _part_stri = part_stri[i_part];
    int* _part_data = part_data[i_part];

    /* First pass to count */
    int idx_part  = 0;
    int n_connect = 0;
    for(int i_entity = 0; i_entity < pn_entity[i_part]; ++i_entity) {
      n_connect += _part_stri[i_entity];
    }

    /* Rebuild in a specific array all the information to sort properly */
    n_connect = n_connect/3;
    PDM_g_num_t* connect_entity = (PDM_g_num_t * ) malloc( n_connect * 3 * sizeof(PDM_g_num_t ));
    int*         connect_info   = (int         * ) malloc( n_connect * 2 * sizeof(int         ));
    //printf("[%i][%i] n_connect::%d\n", i_rank, i_part, n_connect);

    idx_part  = 0;
    n_connect = 0;
    for(int i_entity = 0; i_entity < pn_entity[i_part]; ++i_entity) {
      for(int i_data = 0; i_data < _part_stri[i_entity]/3; ++i_data) {
        // printf(" idx_part = %d\n", idx_part);
        int opp_rank = _part_data[3*idx_part  ];
        int opp_part = _part_data[3*idx_part+1];

        /* On enleve de part_data les informations inutiles (car connu sur cette partition ) */
        int is_distant_bound = 1;
        if( opp_rank == i_rank ){
          if( opp_part == i_part ) {
            is_distant_bound = 0;
          }
        }

        if( is_distant_bound ) {
          connect_entity[3*n_connect  ] = opp_rank;
          connect_entity[3*n_connect+1] = opp_part;
          connect_entity[3*n_connect+2] = pentity_ln_to_gn[i_part][i_entity];

          connect_info[2*n_connect  ] = i_entity;
          connect_info[2*n_connect+1] = _part_data[3*idx_part+2];

          n_connect += 1;
        }
        idx_part += 1;
      }
    }
    // printf("n_connect_sort::%d\n", n_connect);

    /*
     * Sort
     */
    int* order = (int *) malloc( n_connect * sizeof(int) );
    PDM_order_gnum_s(connect_entity, 3, order, n_connect);

    /*
     * Panic verbose
     */
    if(0 == 1){
      printf("order::\n");
      for(int i = 0; i < n_connect; ++i){
        printf("%d ", order[i]);
      }
      printf("\n");
    }

    /* We need to recompute for each opposite part */
    _pproc_bound_idx[i_part]   = PDM_array_zeros_int(n_rank + 1);
    _ppart_bound_idx[i_part]   = PDM_array_zeros_int(part_distribution[n_rank] + 1);
    _pentity_bound[i_part]     = (int *) malloc( ( 4 * n_connect                ) * sizeof(int));


    /* Rebuild */
    for(int i = 0; i < n_connect; ++i){
      int idx = order[i];

      int opp_proc   = connect_entity[3*idx  ];
      int opp_part   = connect_entity[3*idx+1];

      // printf(" i::%d | idx::%d \n", i, idx);
      // printf(" opp_proc::%d | opp_part::%d \n", opp_proc, opp_part);
      int opp_part_g = opp_part + part_distribution[opp_proc];

      _pproc_bound_idx[i_part][opp_proc  +1] += 1;
      _ppart_bound_idx[i_part][opp_part_g+1] += 1;

      _pentity_bound[i_part][4*i  ] = connect_info[2*idx  ] + 1;
      _pentity_bound[i_part][4*i+3] = connect_info[2*idx+1] + 1;

      _pentity_bound[i_part][4*i+1] = opp_proc;
      _pentity_bound[i_part][4*i+2] = opp_part + 1;

    }

    /* Transform in idex */
    for(int i = 0; i < n_rank; ++i){
      _pproc_bound_idx[i_part][i+1] += _pproc_bound_idx[i_part][i];
    }
    for(int i = 0; i < part_distribution[n_rank]; ++i){
      _ppart_bound_idx[i_part][i+1] +=  _ppart_bound_idx[i_part][i];
    }

    free(order);
    free(connect_entity);
    free(connect_info);

  }

  // printf(" PDM_part_assemble_partitions PDM_part_generate_entity_graph_comm flag 6 \n");
  // PDM_MPI_Barrier(comm);
  /*
   * Panic Verbose
   */
  if(0 == 1){
    for(int i_part = 0; i_part < n_part; ++i_part){

      printf("[%d] _pproc_bound_idx::", i_part);
      for(int i = 0; i < n_rank+1; ++i){
        printf("%d ", _pproc_bound_idx[i_part][i]);
      }
      printf("\n");

      printf("[%d] _ppart_bound_idx::", i_part);
      for(int i = 0; i < part_distribution[n_rank]; ++i){
        printf("%d ", _ppart_bound_idx[i_part][i]);
      }
      printf("\n");
    }

  }

  /*
   * Free
   */
  for(int i_part = 0; i_part < n_part; ++i_part){
    free(part_stri[i_part]);
    free(part_data[i_part]);
  }
  free(part_stri);
  free(part_data);
  PDM_block_to_part_free(btp);
}

/**
 *  \brief Recover partitioned coordinates from distributed coordinates and
 *   vertex ln_to_gn indirection.
 *   This function basically calls PDM_block_to_part on to exchange vertex coordinates.
 *
 * \param [in]   comm                PDM_MPI communicator
 * \param [in]   n_part              Number of partitions
 * \param [in]   vertex_distribution Distribution of vertices over the processes (size=n_rank+1)
 * \param [in]   dvtx_coord          Coordinates of distributed vertices (size=3*dn_vtx)
 * \param [in]   pn_vtx              Number of vertices in each partition (size=n_part)
 * \param [in]   pvtx_ln_to_gn       For each part, position of vertices in the global numbering
 *                                   (size = n_part, each component size = pn_vtx[i_part])
 * \param [out]  pvtx_coord          Coordinates of partitioned vertices for each partition
 *                                   (size = n_part, each component size = 3*pn_vtx[i_part])
 */
void
PDM_part_dcoordinates_to_pcoordinates
(
  const PDM_MPI_Comm    comm,
  const int             n_part,
  const PDM_g_num_t    *vertex_distribution,
  const double         *dvtx_coord,
  const int            *pn_vtx,
  const PDM_g_num_t   **pvtx_ln_to_gn,
        double       ***pvtx_coord
)
{
  PDM_block_to_part_t* btp = PDM_block_to_part_create(vertex_distribution,
                               (const PDM_g_num_t **) pvtx_ln_to_gn,
                                                      pn_vtx,
                                                      n_part,
                                                      comm);

  int cst_stride = 3;
  int **pvtx_stride = NULL;
  PDM_block_to_part_exch2(btp,
                          sizeof(double),
                          PDM_STRIDE_CST,
                          &cst_stride,
             (void *  )   dvtx_coord,
             (int  ***)  &pvtx_stride,
             (void ***)   pvtx_coord);

  PDM_block_to_part_free(btp);
}

/**
 *  \brief Recover partitioned coordinates from distributed coordinates and
 *   vertex ln_to_gn indirection.
 *   This function basically calls PDM_block_to_part on to exchange vertex coordinates.
 *
 * \param [in]   comm                PDM_MPI communicator
 * \param [in]   n_part              Number of partitions
 * \param [in]   vertex_distribution Distribution of vertices over the processes (size=n_rank+1)
 * \param [in]   dvtx_coord          Coordinates of distributed vertices (size=3*dn_vtx)
 * \param [in]   pn_vtx              Number of vertices in each partition (size=n_part)
 * \param [in]   pvtx_ln_to_gn       For each part, position of vertices in the global numbering
 *                                   (size = n_part, each component size = pn_vtx[i_part])
 * \param [out]  pvtx_coord          Coordinates of partitioned vertices for each partition
 *                                   (size = n_part, each component size = 3*pn_vtx[i_part])
 */
void
PDM_part_dfield_to_pfield
(
  const PDM_MPI_Comm    comm,
  const int             n_part,
  size_t                s_data,
  const PDM_g_num_t    *field_distribution,
  const unsigned char  *dfield,
  const int            *pn_field,
  const PDM_g_num_t   **pfield_ln_to_gn,
        unsigned char ***pfield
)
{
  PDM_block_to_part_t* btp = PDM_block_to_part_create(field_distribution,
                               (const PDM_g_num_t **) pfield_ln_to_gn,
                                                      pn_field,
                                                      n_part,
                                                      comm);

  int cst_stride = 1;
  int **pfield_stride = NULL;
  PDM_block_to_part_exch2(btp,
                          s_data,
                          PDM_STRIDE_CST,
                          &cst_stride,
             (void *  )   dfield,
             (int  ***)  &pfield_stride,
             (void ***)   pfield);

  PDM_block_to_part_free(btp);
}

/**
 *  \brief Recover partitioned coordinates from distributed coordinates and
 *   vertex ln_to_gn indirection.
 *   This function basically calls PDM_block_to_part on to exchange vertex coordinates.
 *
 * \param [in]   comm                PDM_MPI communicator
 * \param [in]   n_part              Number of partitions
 * \param [in]   vertex_distribution Distribution of vertices over the processes (size=n_rank+1)
 * \param [in]   dvtx_coord          Coordinates of distributed vertices (size=3*dn_vtx)
 * \param [in]   pn_vtx              Number of vertices in each partition (size=n_part)
 * \param [in]   pvtx_ln_to_gn       For each part, position of vertices in the global numbering
 *                                   (size = n_part, each component size = pn_vtx[i_part])
 * \param [out]  pvtx_coord          Coordinates of partitioned vertices for each partition
 *                                   (size = n_part, each component size = 3*pn_vtx[i_part])
 */
void
PDM_part_dfield_to_pfield2
(
  const PDM_MPI_Comm     comm,
  const int              n_part,
  size_t                 s_data,
  PDM_stride_t           t_stride,
  const PDM_g_num_t     *field_distribution,
  const int             *dfield_stri,
  const unsigned char   *dfield,
  const int             *pn_field,
  const PDM_g_num_t    **pfield_ln_to_gn,
        int           ***pfield_stride,
        unsigned char ***pfield
)
{
  PDM_block_to_part_t* btp = PDM_block_to_part_create(field_distribution,
                               (const PDM_g_num_t **) pfield_ln_to_gn,
                                                      pn_field,
                                                      n_part,
                                                      comm);

  PDM_block_to_part_exch2(btp,
                          s_data,
                          t_stride,
             (int  *  )   dfield_stri,
             (void *  )   dfield,
             (int  ***)   pfield_stride,
             (void ***)   pfield);

  PDM_block_to_part_free(btp);
}



void
PDM_part_dentity_group_to_pentity_group
(
  const PDM_MPI_Comm     comm,
  const int              n_part,
  const PDM_g_num_t     *entity_distribution,
  const int             *dentity_group_idx,
  const int             *dentity_group,
  const int             *pn_entity,
  const PDM_g_num_t    **pentity_ln_to_gn,
  int                 ***pentity_group_idx,
  int                 ***pentity_group
)
{
  int i_rank;
  // int n_rank;

  PDM_MPI_Comm_rank(comm, &i_rank);
  // PDM_MPI_Comm_size(comm, &n_rank);

  int dn_entity = entity_distribution[i_rank+1] - entity_distribution[i_rank];

  int* dentity_group_n = (int * ) malloc( dn_entity * sizeof(int));
  for(int i = 0; i < dn_entity; ++i) {
    dentity_group_n[i] = dentity_group_idx[i+1] - dentity_group_idx[i];
  }
  int **pentity_group_n = NULL;
  PDM_part_dfield_to_pfield2(comm,
                             n_part,
                             sizeof(int),
                             PDM_STRIDE_VAR,
                             entity_distribution,
                             dentity_group_n,
         (unsigned char *  ) dentity_group,
                             pn_entity,
                             pentity_ln_to_gn,
                             &pentity_group_n,
         (unsigned char ***) pentity_group);


  int **_pentity_group_idx = malloc( (n_part) * sizeof(int *));
  for(int i_part = 0; i_part < n_part; ++i_part) {
    _pentity_group_idx[i_part] = malloc( (pn_entity[i_part]+1) * sizeof(int));
    _pentity_group_idx[i_part][0] = 0;
    //log_trace("pn_entity[i_part] = %i \n", pn_entity[i_part]);
    for(int i = 0; i < pn_entity[i_part]; ++i) {
      //log_trace(" [%i] --> %i \n", i, pentity_group_n[i_part][i]);
      _pentity_group_idx[i_part][i+1] = _pentity_group_idx[i_part][i] + pentity_group_n[i_part][i];
    }
    free(pentity_group_n[i_part]);
  }
  free(pentity_group_n);
  free(dentity_group_n);

  *pentity_group_idx = _pentity_group_idx;
}


// void
// PDM_part_multi_dfield_to_pfield
// (
//   const PDM_MPI_Comm      comm,
//   const int               n_part,
//   const int               n_field,
//   size_t                 *s_data,
//   const PDM_g_num_t      *field_distribution,
//   const unsigned char   **dfield,
//   const int              *pn_field,
//   const PDM_g_num_t     **pfield_ln_to_gn,
//         unsigned char ****pfield
// )
// {
// TO PDM_STRI_VAR
// }

/**
 *  \brief Extend an existing ln_to_gn from a connectivity
 *
 * \param [in]   comm                PDM_MPI communicator
 * \param [in]   part_distribution   Distribution of partitions over the processes (size=n_rank+1)
 * \param [in]   entity_distribution Distribution of entities over the processes (size=n_rank+1)
 * \param [in]   dentity_to_part     Id of assigned partition for each entity (size=dn_entity)
 * \param [out]  pn_entities         Number of entities in each partition (size n_part)
 * \param [out]  pentity_ln_to_gn    Array of local to global entity id for each partition (size n_part)
 *
 * \return       n_part              Number of partitions managed by this process
 */
void
PDM_extend_mesh
(
 const PDM_MPI_Comm    comm,
 const PDM_g_num_t    *part_distribution,
 const PDM_g_num_t    *entity_distribution,
 const int            *dentity_to_part,
 const int             n_part,
 const PDM_g_num_t    *dual_graph_idx,
 const PDM_g_num_t    *dual_graph,
 const int            *pn_entity,
       PDM_g_num_t   **pentity_ln_to_gn,
       int           **pn_entity_extented,
       PDM_g_num_t  ***pentity_ln_to_gn_extended
)
{
  PDM_UNUSED(comm);
  PDM_UNUSED(part_distribution);
  PDM_UNUSED(entity_distribution);
  PDM_UNUSED(dentity_to_part);
  PDM_UNUSED(n_part);
  PDM_UNUSED(dual_graph_idx);
  PDM_UNUSED(dual_graph);
  PDM_UNUSED(pn_entity);
  PDM_UNUSED(pentity_ln_to_gn);
  PDM_UNUSED(pn_entity_extented);
  PDM_UNUSED(pentity_ln_to_gn_extended);

  int i_rank;
  int n_rank;

  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  int dn_entity  = entity_distribution[i_rank+1]  -  entity_distribution[i_rank];

  /*
   * We search to extended the partition with the dual graph
   *     - We begin by get the dual_graph for each partition
   *
   */

  /*
   * Prepare exchange protocol
   */
  PDM_block_to_part_t* btp = PDM_block_to_part_create(entity_distribution,
                               (const PDM_g_num_t **) pentity_ln_to_gn,
                                                      pn_entity,
                                                      n_part,
                                                      comm);

  int* dual_graph_n = (int *) malloc( dn_entity * sizeof(int));
  for(int i = 0; i < dn_entity; ++i){
    dual_graph_n[i] = dual_graph_idx[i+1] - dual_graph_idx[i];
  }

  // PDM_log_trace_array_int(dual_graph_idx, dn_entity+1, "dual_graph_idx::");
  // PDM_log_trace_array_long(dual_graph, dual_graph_idx[dn_entity], "dual_graph::");

  /*
   * Exchange
   */
  int**         part_dual_graph_n;
  PDM_g_num_t** part_dual_graph;
  PDM_block_to_part_exch2(btp,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_VAR,
                          dual_graph_n,
             (void *  )   dual_graph,
             (int  ***)  &part_dual_graph_n,
             (void ***)  &part_dual_graph);
  free(dual_graph_n);

  int** part_dual_graph_idx = (int ** ) malloc( n_part * sizeof(int*));
  for(int i_part = 0; i_part < n_part; ++i_part) {
    part_dual_graph_idx[i_part] = PDM_array_new_idx_from_sizes_int(part_dual_graph_n[i_part], pn_entity[i_part]);
    // PDM_log_trace_array_int(part_dual_graph_n[i_part], pn_entity[i_part], "part_dual_graph_n[i_part]::");
    // PDM_log_trace_array_int(part_dual_graph_idx[i_part], pn_entity[i_part]+1, "part_dual_graph_idx[i_part]::");
    // PDM_log_trace_array_long(part_dual_graph[i_part], part_dual_graph_idx[i_part][pn_entity[i_part]], "part_dual_graph[i_part]::");
  }

  /*
   * Each partition have now for each cells the neighbour cells (via dual_graph )
   *     - Sort / Unique
   *     - We collect all boundary cells
   *     - we append to the older ln_to_gn all new entitiy
   */
  *pn_entity_extented = (int *) malloc( sizeof(int) * n_part);
  int* _pn_entity_extented = *pn_entity_extented;

  *pentity_ln_to_gn_extended = (PDM_g_num_t **) malloc( sizeof(PDM_g_num_t *) * n_part);
  PDM_g_num_t** _pentity_ln_to_gn_extended = *pentity_ln_to_gn_extended;
  for(int i_part = 0; i_part < n_part; ++i_part) {

    int new_size = PDM_inplace_unique_long(part_dual_graph[i_part], NULL, 0, part_dual_graph_idx[i_part][pn_entity[i_part]]-1);
    printf(" new_size         :: %i \n", new_size);
    printf(" pn_entity[i_part]:: %i \n", pn_entity[i_part]);

    _pentity_ln_to_gn_extended[i_part] = (PDM_g_num_t *) malloc( new_size * sizeof(PDM_g_num_t));

    for(int i_entity = 0; i_entity < pn_entity[i_part]; ++i_entity ) {
      _pentity_ln_to_gn_extended[i_part][i_entity] = pentity_ln_to_gn[i_part][i_entity];
    }

    PDM_sort_long(_pentity_ln_to_gn_extended[i_part], NULL, pn_entity[i_part]);

    int n_new_cell = pn_entity[i_part];
    for(int i_entity = 0; i_entity < new_size; ++i_entity) {
      /*
       * Si non existant on rajoute au bout
       */
      int ipos = PDM_binary_search_long(part_dual_graph[i_part][i_entity]+1,
                                        _pentity_ln_to_gn_extended[i_part],
                                        pn_entity[i_part]);

      // printf("ipos = %i | %i \n", ipos, part_dual_graph[i_part][i_entity]+1);

      // Append at end of the array
      if(ipos == -1) {
        _pentity_ln_to_gn_extended[i_part][n_new_cell++] = part_dual_graph[i_part][i_entity]+1;
      }
    }

    // /*
    //  * Compress
    //  */
    // for(int i_entity = 0; i_entity < pn_entity[i_part]; ++i_entity ) {
    //   _pentity_ln_to_gn_extended[i_part][i_entity] = pentity_ln_to_gn[i_part][i_entity];
    // }

    // printf("n_new_cell::%i\n", n_new_cell);
    // PDM_log_trace_array_long(pentity_ln_to_gn[i_part], pn_entity[i_part], "_pentity_ln_to_gn[i_part]::");
    // PDM_log_trace_array_long(_pentity_ln_to_gn_extended[i_part], n_new_cell, "_pentity_ln_to_gn_extended[i_part]::");

    // Realloc and setup size
    _pn_entity_extented[i_part] = n_new_cell;
    _pentity_ln_to_gn_extended[i_part] = (PDM_g_num_t *) realloc( _pentity_ln_to_gn_extended[i_part], n_new_cell * sizeof(PDM_g_num_t));

  }

  PDM_block_to_part_free(btp);
  for(int i_part = 0; i_part < n_part; ++i_part) {
    free(part_dual_graph_n[i_part]);
    free(part_dual_graph_idx[i_part]);
    free(part_dual_graph[i_part]);
  }
  free(part_dual_graph_n);
  free(part_dual_graph_idx);
  free(part_dual_graph);
}

void
PDM_setup_connectivity_idx
(
 int           dn_entity1,
 int           stride,
 PDM_g_num_t  *dentity1_dentity2,
 int         **dentity1_dentity2_idx,
 PDM_g_num_t **dentity1_dentity2_new
)
{
  int         *_dentity1_dentity2_idx = (int         *) malloc( (dn_entity1 + 1     ) * sizeof(int        ));
  PDM_g_num_t *_dentity1_dentity2_new = (PDM_g_num_t *) malloc( (stride * dn_entity1) * sizeof(PDM_g_num_t));

  _dentity1_dentity2_idx[0] = 0;
  for(int i = 0; i < dn_entity1; ++i) {
    _dentity1_dentity2_idx[i+1] = _dentity1_dentity2_idx[i];
    for(int is = 0; is < stride; ++is) {
      if(dentity1_dentity2[stride*i+is] != 0) {
        _dentity1_dentity2_new[_dentity1_dentity2_idx[i+1]++] = dentity1_dentity2[stride*i+is];
      }
    }
  }

  _dentity1_dentity2_new = realloc(_dentity1_dentity2_new, _dentity1_dentity2_idx[dn_entity1] * sizeof(PDM_g_num_t));

  *dentity1_dentity2_idx = _dentity1_dentity2_idx;
  *dentity1_dentity2_new = _dentity1_dentity2_new;
}



#ifdef __cplusplus
}
#endif /* __cplusplus */
