
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
#include "pdm_priv.h"
#include "pdm_mpi.h"
#include "pdm_dmesh_nodal.h"
#include "pdm_dmesh_nodal_priv.h"
#include "pdm_part_to_block.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_gnum.h"
#include "pdm_quick_sort.h"
#include "pdm_geom_elem.h"

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

/*----------------------------------------------------------------------------
 * Maximum number of sections depending of section type
 *----------------------------------------------------------------------------*/

/*typedef enum {

  PDM_BLOCK_ID_BLOCK_STD    = 0,
  PDM_BLOCK_ID_BLOCK_POLY2D = 1000000,
  PDM_BLOCK_ID_BLOCK_POLY3D = 2000000

  } PDM_section_id_section_t;*/ // --> moved to pdm_mesh_nodal.h

/*============================================================================
 * Global variable
 *============================================================================*/

/**
 * \brief Storage of mesh handles
 *
 */

static PDM_Handles_t *mesh_handles = NULL;

/*============================================================================
 * Private function definitions
 *============================================================================*/

/**
 * \def _find_pairs
 * Search common faces in a distribution
 *
 */
static int
_find_pairs
(
const int *idx_face,
const int *data,
const int  nFac,
      int  i_abs_face,  /* A passer ren reference ou a return */
      PDM_g_num_t *dface_vtx,
      int *dface_vtx_idx,
      PDM_g_num_t *dface_cell
)
{
  if(0 == 1){
    printf("_find_pairs : nFac : %d \n", nFac);
  }

  /*
   * Special Case -> To removed when Boundary is OK
   */
  // if(nFac == 1){
  //   // printf("Something strange append for in find_pairs - Check mesh ooo \n");

  //   /* Stokage en element externe */
  //   int iFacIdx = dface_vtx_idx[i_abs_face];

  //   int curFac = idx_face[0];
  //   int n_vtx1  = data[curFac+2];


  //   dface_cell[2*i_abs_face  ] = data[curFac];
  //   dface_cell[2*i_abs_face+1] = 0;

  //   // printf("----------: %d - %d - %d\n", curFac, n_vtx1, iFacIdx);


  //   for(int i_vtx=0; i_vtx < n_vtx1; i_vtx++){
  //     dface_vtx[iFacIdx+i_vtx] = data[curFac+3+i_vtx];
  //   }

  //   /* Update index */
  //   dface_vtx_idx[i_abs_face+1] = dface_vtx_idx[i_abs_face] + n_vtx1;
  //   i_abs_face++;

  //   return i_abs_face;
  //   // exit(1);
  // }

  /*
   * Make array of already treated face
   */
  int *AlreadyTreat = (int * ) malloc( sizeof(int *) * nFac);
  for(int i=0; i < nFac; i++){
    AlreadyTreat[i] = -1;
  }

  /*
   * Begin with the first face of the array
   */
  for(int iPos=0; iPos < nFac; iPos++){

    /** Get the current Idx in dface_vtx **/
    int iFacIdx = dface_vtx_idx[i_abs_face];

    /*
     * Verbose
     */
    // printf("---------------------------- : %d - %d - %d\n", iPos, idx_face[iPos ], data[idx_face[iPos ]]);

    /*
     * Face deja traité ??
     */
    if(AlreadyTreat[iPos] != 1){

      int curFac = idx_face[iPos ];
      int n_vtx1  = data[curFac+2];

      /** Prepare ElemCon **/
      int *ElmCon1 = (int *) malloc( sizeof(int *) * n_vtx1); /* To sort out of loop */
      int  tKey1   = 0;
      for(int i_vtx=0; i_vtx < n_vtx1; i_vtx++){
        ElmCon1[i_vtx] = data[curFac+3+i_vtx];
        tKey1        += data[curFac+3+i_vtx];
      }
      PDM_quick_sort_int(ElmCon1, 0, n_vtx1-1);

      if(0 == 1){
        for(int i_vtx=0; i_vtx<n_vtx1; i_vtx++){
          printf("ElmCon1[%d] : %d \n", i_vtx, ElmCon1[i_vtx]);
        }
      }

      /** Loop on candidates faces **/
      for(int iNex=iPos+1; iNex < nFac; iNex++){

        /*
         * Verbose
         */
        // printf("++++++++++++++++++++++++++++++++++ : %d - %d - %d \n", iNex, idx_face[iNex ], data[idx_face[iNex ]]);

        /*
         * On sait que la prochaine a traiter est forcement plus grande
         */
        int nexFac = idx_face[iNex ];
        int n_vtx2  = data[nexFac+2];

        // printf("+++++++++++++++++++ : %d - %d \n", n_vtx1, n_vtx2);
        /*
         * First sort - Number of Vertex is different -> Si Vtx different c'est pas la meme
         */
        if(n_vtx1 == n_vtx2){

          /*
           * Allocate memory and copy
           */
          int *ElmCon2 = (int *) malloc( sizeof(int *) * n_vtx1);
          int  tKey2   = 0;
          for(int i_vtx=0; i_vtx < n_vtx1; i_vtx++){
            ElmCon2[i_vtx] = data[nexFac+3+i_vtx];
            tKey2        += data[nexFac+3+i_vtx];
          }

          /*
           * Sort ElemCon2
           */
          PDM_quick_sort_int(ElmCon2, 0, n_vtx2-1);

          if(0 == 1){
            for(int i_vtx=0; i_vtx < n_vtx1; i_vtx++){
              printf("ElmCon2[%d] : %d \n", i_vtx, ElmCon2[i_vtx]);
            }
          }

          /** Compare **/
          int isSameFace = 1;
          for(int i_vtx=0; i_vtx < n_vtx1; i_vtx++){
             if(ElmCon1[i_vtx] != ElmCon2[i_vtx]){
                isSameFace = -1;
                break;
             }
          }

          /*
           * Verbose
           */
          if(0 == 1){
            printf("tKey : %d / %d \n", tKey1, tKey2);
            assert(tKey1 == tKey2);
          }

          /*
           * Fill the Face data if Face is the same
           */
          if(isSameFace == 1){

            /*
             * Stockage de parent elements
             *    -> TODOUX -> Manage sign if Entrante /Sortante
             *    Rendre un contenu explicit : - si sortant // + si entrant
             *    Il faut que le parent element d'une boundary soit le droit ...
             */
            // if(data[curFac] == -1){
            // if(data[curFac] > n_cell ){
            //   dface_cell[2*i_abs_face  ] = data[nexFac];
            //   // dface_cell[2*i_abs_face+1] = 0;//data[curFac];
            //   dface_cell[2*i_abs_face+1] = data[curFac];

            //   /*
            //    * Fill face_vtx connectivity array
            //    */
            //   for(int i_vtx=0; i_vtx < n_vtx1; i_vtx++){
            //     dface_vtx[iFacIdx+i_vtx] = data[nexFac+3+i_vtx];
            //   }
            //   printf("Impossible case now \n");
            //   exit(1);
            // }
            // else
            // {
              dface_cell[2*i_abs_face  ] = data[curFac];
              dface_cell[2*i_abs_face+1] = data[nexFac];

              /*
               * Fill face_vtx connectivity array
               */
              for(int i_vtx=0; i_vtx < n_vtx1; i_vtx++){
                dface_vtx[iFacIdx+i_vtx] = data[curFac+3+i_vtx];
              }
            // }

            /*
             * Stockage face_vtx connectivity Index
             */
            // printf("dface_vtx_idx est faux %d \n", n_vtx1);
            // dface_vtx_idx[++i_abs_face] = n_vtx1-1;
            // printf("dface_vtx_idx est faux [%d] -> %d \n", i_abs_face+1, dface_vtx_idx[i_abs_face+1]);
            dface_vtx_idx[i_abs_face+1] = dface_vtx_idx[i_abs_face] + n_vtx1;
            i_abs_face++;

            /*
             * Flags the two faces as treated
             */
            AlreadyTreat[iPos] = 1;
            AlreadyTreat[iNex] = 1;

          } /** End if (isSameFace) **/

          /** Free Cendidate ElemCon **/
          free(ElmCon2);

        } /* Enf if Same Number of Vertex */

      }
      /** Free current ElemCon **/
      free(ElmCon1);

    } /** End If alreadyTreat **/

    /* Boundary management **/
    if(AlreadyTreat[iPos] != 1){

      // printf("i_abs_face : %d \n", i_abs_face);
      // printf("iPos : %d \n", iPos);
      // printf("----------: %d - %d - %d\n", curFac, n_vtx1, iFacIdx);

      int iFacIdx2 = dface_vtx_idx[i_abs_face];
      int curFac  = idx_face[iPos];
      int n_vtx1   = data[curFac+2];

      dface_cell[2*i_abs_face  ] = data[curFac];
      dface_cell[2*i_abs_face+1] = 0;

      for(int i_vtx=0; i_vtx < n_vtx1; i_vtx++){
        dface_vtx[iFacIdx2+i_vtx] = data[curFac+3+i_vtx];
      }

      /* Update index */
      dface_vtx_idx[i_abs_face+1] = dface_vtx_idx[i_abs_face] + n_vtx1;
      i_abs_face++;

    }

  }

  /** Boundary management - Temporary **/


  /** Free **/
  free(AlreadyTreat);

  return i_abs_face;
}

static void
_make_absolute_face_numbering(PDM_DMesh_nodal_t* mesh)
{

  PDM_g_num_t n_faceProc = mesh->dn_face;
  PDM_g_num_t beg_NumAbs;

  PDM_MPI_Scan(&n_faceProc, &beg_NumAbs, 1, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, mesh->pdm_mpi_comm);
  beg_NumAbs -= n_faceProc;


  /** Compute the distribution of elements amont proc **/
  mesh->face_distrib = (PDM_g_num_t *) malloc((mesh->n_proc+1) * sizeof(PDM_g_num_t));
  PDM_g_num_t _dn_face = (PDM_g_num_t) mesh->dn_face;
  PDM_MPI_Allgather((void *) &_dn_face,
                    1,
                    PDM__PDM_MPI_G_NUM,
                    (void *) (&mesh->face_distrib[1]),
                    1,
                    PDM__PDM_MPI_G_NUM,
                    mesh->pdm_mpi_comm);

  // mesh->face_distrib[0] = 1;
  mesh->face_distrib[0] = 0;

  for (int i = 1; i < mesh->n_proc+1; i++) {
    mesh->face_distrib[i] +=  mesh->face_distrib[i-1];
  }

  if (0 == 1) {
    printf("beg_NumAbs::Face : "PDM_FMT_G_NUM" \n", beg_NumAbs);
    printf("mesh->face_distrib : "PDM_FMT_G_NUM,  mesh->face_distrib[0]);
    for (int i = 1; i < mesh->n_proc+1; i++) {
      printf(" "PDM_FMT_G_NUM, mesh->face_distrib[i]);
    }
    printf("\n");
  }
}

/**
 *
 * \brief Free vtx structure
 *
 * \param[inout]  vtx    Vertices
 *
 * \return        NULL
 *
 */

static
PDM_DMesh_nodal_vtx_t *
_vtx_free
(
 PDM_DMesh_nodal_vtx_t *vtx
)
{

  if (vtx != NULL) {

    if (vtx->distrib != NULL) {
      free (vtx->distrib);
      vtx->distrib = NULL;
    }

    free (vtx);
  }
  return NULL;
}


/**
 *
 * \brief Initialize a mesh
 *
 * \param [inout]  mesh        Mesh
 * \param [in]     n_part      Number of partitions
 */

static void
_mesh_init
(
PDM_DMesh_nodal_t *mesh,
const PDM_MPI_Comm comm,
      PDM_g_num_t  n_vtx,
      PDM_g_num_t  n_cell
)
{
  int n_proc;
  int i_proc;

  PDM_MPI_Comm_size (comm, &n_proc);
  PDM_MPI_Comm_rank (comm, &i_proc);

  mesh->n_proc                   = n_proc;
  mesh->i_proc                   = i_proc;

  mesh->n_som_abs                = n_vtx;
  mesh->n_cell_abs               = n_cell;

  mesh->vtx                      = malloc(sizeof(PDM_DMesh_nodal_vtx_t ));
  mesh->vtx->_coords             = NULL;
  mesh->vtx->n_vtx               = 0;

  mesh->sections_std             = NULL;
  mesh->sections_poly2d          = NULL;
  mesh->sections_poly3d          = NULL;

  mesh->pdm_mpi_comm             = comm;
  mesh->sections_id              = NULL;
  mesh->n_sections               = 0;

  mesh->n_dcell                  = -1;
  mesh->dcell_face                = NULL;
  mesh->dcell_face_idx             = NULL;
  mesh->cell_distrib             = NULL;

  mesh->dn_face                   = -1;
  mesh->_dface_vtx                = NULL;
  mesh->_dface_vtx_idx             = NULL;
  mesh->_dface_cell               = NULL;
  mesh->face_distrib             = NULL;

}

/**
 *
 * \brief Update sections identifier list
 *
 * \param [inout]  mesh        Mesh
 */

static void
_update_sections_id
(
PDM_DMesh_nodal_t *mesh
)
{
  int n_sections = 0;

  if (mesh->sections_std != NULL) {
    n_sections += PDM_Handles_n_get (mesh->sections_std);
  }

  if (mesh->sections_poly2d != NULL) {
    n_sections += PDM_Handles_n_get (mesh->sections_poly2d);
  }

  if (mesh->sections_poly3d != NULL) {
    n_sections += PDM_Handles_n_get (mesh->sections_poly3d);
  }

  if (mesh->n_sections < n_sections) {
    mesh->sections_id = (int *) realloc(mesh->sections_id, sizeof(int) * n_sections);
  }

  int k = 0;
  if (mesh->sections_std != NULL) {
    const int *id1 = PDM_Handles_idx_get (mesh->sections_std);
    int n = PDM_Handles_n_get (mesh->sections_std);
    for (int i = 0; i < n; i++) {
      mesh->sections_id[k++] = id1[i] + PDM_BLOCK_ID_BLOCK_STD;
    }
  }

  if (mesh->sections_poly2d != NULL) {
    const int *id1 = PDM_Handles_idx_get (mesh->sections_poly2d);
    int n = PDM_Handles_n_get (mesh->sections_poly2d);
    for (int i = 0; i < n; i++) {
      mesh->sections_id[k++] = id1[i] + PDM_BLOCK_ID_BLOCK_POLY2D;
    }
  }

  if (mesh->sections_poly3d != NULL) {
    const int *id1 = PDM_Handles_idx_get (mesh->sections_poly3d);
    int n = PDM_Handles_n_get (mesh->sections_poly3d);
    for (int i = 0; i < n; i++) {
      mesh->sections_id[k++] = id1[i] + PDM_BLOCK_ID_BLOCK_POLY3D;
    }
  }

  mesh->n_sections = n_sections;

}


/**
 *
 * \brief Free a standard section
 *
 * \param [inout]  _bloc_std    Standard section
 *
 * \return         Null
 *
 */

static
PDM_DMesh_nodal_section_std_t *
_section_std_free
(
PDM_DMesh_nodal_section_std_t *_section_std
)
{

  if (_section_std == NULL) {
    return NULL;
  }

  if (_section_std->distrib != NULL) {
    free (_section_std->distrib);
    _section_std->distrib = NULL;
  }

  free(_section_std);
  return NULL;
}


/**
 *
 * \brief Free a polygon section
 *
 * \param [inout]  _bloc_poly2d    Polygon section
 *
 * \return         Null
 *
 */

static
PDM_DMesh_nodal_section_poly2d_t *
_section_poly2d_free
(
PDM_DMesh_nodal_section_poly2d_t *_section_poly2d
)
{

  if (_section_poly2d == NULL) {
    return NULL;
  }

  if (_section_poly2d->distrib != NULL) {
    free (_section_poly2d->distrib);
    _section_poly2d->distrib = NULL;
  }


  free(_section_poly2d);

  return NULL;
}


/**
 *
 * \brief Free a polyhedron section
 *
 * \param [inout]  _section_poly3d    Polyhedron section
 *
 * \return         Null
 *
 */

static
PDM_DMesh_nodal_section_poly3d_t *
_section_poly3d_free
(
PDM_DMesh_nodal_section_poly3d_t *_section_poly3d
)
{
  if (_section_poly3d == NULL) {
    return NULL;
  }

  if (_section_poly3d->distrib != NULL) {
    free (_section_poly3d->distrib);
    _section_poly3d->distrib = NULL;
  }


  free(_section_poly3d);

  return NULL;

}



///**
// *
// * \brief _section_elt_faces_get
// *
// * \param [in]     mesh               Current mesh
// * \param [in]     id_section         Section identifier
// * \param [inout]  elt_face_vtx_idx   Index of element faces connectivity (preallocated)
// * \param [inout]  elt_face_vtx       Element faces connectivity (preallocated)
// *
// */

//static int
//_section_size_elt_faces_get
//(
//      PDM_DMesh_nodal_t *mesh,
//      int               *s_elt_face_vtx_idx,
//      int               *s_elt_face_vtx,
//      int               *s_elt_face_cell
//)
//{

//  int _s_elt_face_vtx_idx = 0;
//  int _s_elt_face_vtx = 0;


//  int n_sections_std = PDM_Handles_n_get (mesh->sections_std);
//  const int *list_ind = PDM_Handles_idx_get (mesh->sections_std);

//  for (int i = 0; i < n_sections_std; i++) {
//    PDM_DMesh_nodal_section_std_t *section =
//      (PDM_DMesh_nodal_section_std_t *) PDM_Handles_get (mesh->sections_std, list_ind[i]);
//    int n_face_elt = 0;
//    int n_sum_vtx_face = 0;

//    switch (section->t_elt) {
//    case PDM_MESH_NODAL_TRIA3:
//      n_face_elt = 3;
//      n_sum_vtx_face = 6;
//      break;
//    case PDM_MESH_NODAL_TETRA4:
//      n_face_elt = 4;
//      n_sum_vtx_face = 12;
//      break;
//    case PDM_MESH_NODAL_QUAD4:
//      n_face_elt = 4;
//      n_sum_vtx_face = 8;
//      break;
//    case PDM_MESH_NODAL_HEXA8:
//      n_face_elt = 6;
//      n_sum_vtx_face = 24;
//      break;
//    case PDM_MESH_NODAL_PYRAMID5:
//      n_face_elt = 5;
//      n_sum_vtx_face = 16;
//      break;
//    case PDM_MESH_NODAL_PRISM6:
//      n_face_elt = 5;
//      n_sum_vtx_face = 18;
//      break;
//    default:
//      PDM_error(__FILE__, __LINE__, 0, "Error _section_size_elt_faces_get : Element type is not taking int account\n");
//      exit(EXIT_FAILURE);
//    }

//    _s_elt_face_vtx_idx += section->n_elt * n_face_elt;
//    _s_elt_face_vtx += section->n_elt * n_sum_vtx_face;
//  }

//  int n_sections_poly2d = PDM_Handles_n_get (mesh->sections_poly2d);
//  list_ind = PDM_Handles_idx_get (mesh->sections_poly2d);

//  for (int i = 0; i < n_sections_poly2d; i++) {
//    PDM_DMesh_nodal_section_poly2d_t *section =
//      (PDM_DMesh_nodal_section_poly2d_t *) PDM_Handles_get (mesh->sections_poly2d, list_ind[i]);
//    _s_elt_face_vtx_idx += section->_connec_idx[section->n_elt];
//    _s_elt_face_vtx += 2 * section->_connec_idx[section->n_elt];
//  }

//  int n_sections_poly3d = PDM_Handles_n_get (mesh->sections_poly3d);
//  list_ind = PDM_Handles_idx_get (mesh->sections_poly3d);

//  for (int i = 0; i < n_sections_poly3d; i++) {
//    PDM_DMesh_nodal_section_poly3d_t *section =
//      (PDM_DMesh_nodal_section_poly3d_t *) PDM_Handles_get (mesh->sections_poly3d, list_ind[i]);
//    _s_elt_face_vtx_idx += section->n_face;
//    _s_elt_face_vtx += section->_facvtx[section->_facvtx_idx[section->n_face]];
//  }

//  *s_elt_face_cell =  _s_elt_face_vtx_idx;
//  *s_elt_face_vtx_idx = _s_elt_face_vtx_idx + 1;
//  *s_elt_face_vtx = _s_elt_face_vtx + 1;


//  return *s_elt_face_vtx - 1;

//}


///**
// *
// * \brief _section_elt_faces_get
// *
// * \param [in]     mesh               Current mesh
// * \param [in]     id_section         Section identifier
// * \param [inout]  elt_face_vtx_idx   Index of element faces connectivity (preallocated)
// * \param [inout]  elt_face_vtx       Element faces connectivity (preallocated)
// *
// */

//static void
//_section_elt_faces_add
//(
//      PDM_DMesh_nodal_t *mesh,
//const int                id_section,
//      int               *n_elt_current,
//      int               *n_face_current,
//      int               *elt_face_vtx_idx,
//      PDM_g_num_t       *elt_face_vtx,
//      PDM_g_num_t       *elt_face_cell
//)
//{
//  int _n_face_current = *n_face_current;

//  int         *_current_elt_face_vtx_idx = elt_face_vtx_idx + _n_face_current;
//  PDM_g_num_t *_current_elt_face_vtx     = elt_face_vtx + elt_face_vtx_idx[_n_face_current];
//  PDM_g_num_t *_current_elt_face_cell    = elt_face_cell +_n_face_current;

//  if (id_section < PDM_BLOCK_ID_BLOCK_POLY2D) {

//    PDM_DMesh_nodal_section_std_t *section =
//            (PDM_DMesh_nodal_section_std_t *) PDM_Handles_get (mesh->sections_std, id_section);

//    if (section == NULL) {
//      PDM_error (__FILE__, __LINE__, 0, "Bad section identifier\n");
//    }

//    switch (section->t_elt) {

//    case PDM_MESH_NODAL_TRIA3:
//      {
//        const int n_face_elt        = 3;
//        const int n_sum_vtx_face    = 6;
//        const int n_sum_vtx_elt     = 3;

//        for (int ielt = 0; ielt < section->n_elt; ielt++) {
//          for (int i_face = 0; i_face < n_face_elt; i_face++) {
//            _current_elt_face_vtx_idx[ielt * n_face_elt + i_face + 1] =
//              _current_elt_face_vtx_idx[ielt * n_face_elt + i_face] + 3;
//            _current_elt_face_cell[ielt * n_face_elt + i_face    ] = *n_elt_current + ielt + 1;
//          }

//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 0] = section->_connec[n_sum_vtx_elt * ielt    ];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 1] = section->_connec[n_sum_vtx_elt * ielt + 1];

//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 2] = section->_connec[n_sum_vtx_elt * ielt + 1];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 3] = section->_connec[n_sum_vtx_elt * ielt + 2];

//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 4] = section->_connec[n_sum_vtx_elt * ielt + 2];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 5] = section->_connec[n_sum_vtx_elt * ielt + 3];
//        }

//        *n_elt_current += section->n_elt;
//        *n_face_current += section->n_elt * n_face_elt;
//      }
//      break;

//    case PDM_MESH_NODAL_TETRA4:
//      {
//        const int n_face_elt        = 4;
//        const int n_sum_vtx_face    = 12;
//        const int n_sum_vtx_elt     = 4;

//        for (int ielt = 0; ielt < section->n_elt; ielt++) {
//          for (int i_face = 0; i_face < n_face_elt; i_face++) {
//            _current_elt_face_vtx_idx[ielt * n_face_elt + i_face + 1] =
//              _current_elt_face_vtx_idx[ielt * n_face_elt + i_face] + 3;
//            _current_elt_face_cell[ielt * n_face_elt + i_face    ] = *n_elt_current + ielt + 1;
//          }

//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 0]  = section->_connec[n_sum_vtx_elt * ielt    ];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 1]  = section->_connec[n_sum_vtx_elt * ielt + 1];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 2]  = section->_connec[n_sum_vtx_elt * ielt + 2];

//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 3]  = section->_connec[n_sum_vtx_elt * ielt    ];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 4]  = section->_connec[n_sum_vtx_elt * ielt + 3];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 5]  = section->_connec[n_sum_vtx_elt * ielt + 1];

//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 6]  = section->_connec[n_sum_vtx_elt * ielt    ];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 7]  = section->_connec[n_sum_vtx_elt * ielt + 2];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 8]  = section->_connec[n_sum_vtx_elt * ielt + 3];

//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 9]  = section->_connec[n_sum_vtx_elt * ielt + 1];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 10] = section->_connec[n_sum_vtx_elt * ielt + 3];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 11] = section->_connec[n_sum_vtx_elt * ielt + 2];

//        }

//        *n_elt_current += section->n_elt;
//        *n_face_current += section->n_elt * n_face_elt;

//    }
//      break;
//    case PDM_MESH_NODAL_QUAD4:
//      {
//        const int n_face_elt        = 4;
//        const int n_sum_vtx_face    = 8;
//        const int n_sum_vtx_elt     = 4;

//        for (int ielt = 0; ielt < section->n_elt; ielt++) {
//          for (int i_face = 0; i_face < n_face_elt; i_face++) {
//            _current_elt_face_vtx_idx[ielt * n_face_elt + i_face + 1] =
//              _current_elt_face_vtx_idx[ielt * n_face_elt + i_face] + 3;
//            _current_elt_face_cell[ielt * n_face_elt + i_face    ] = *n_elt_current + ielt + 1;
//          }

//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 0] = section->_connec[n_sum_vtx_elt * ielt    ];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 1] = section->_connec[n_sum_vtx_elt * ielt + 1];

//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 2] = section->_connec[n_sum_vtx_elt * ielt + 1];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 3] = section->_connec[n_sum_vtx_elt * ielt + 2];

//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 4] = section->_connec[n_sum_vtx_elt * ielt + 2];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 5] = section->_connec[n_sum_vtx_elt * ielt + 3];
//        }

//        *n_elt_current  += section->n_elt;
//        *n_face_current += section->n_elt * n_face_elt;

//      }
//      break;

//    case PDM_MESH_NODAL_HEXA8:
//      {
//        const int n_face_elt        = 6;
//        const int n_sum_vtx_face    = 24;
//        const int n_sum_vtx_elt     = 8;

//        for (int ielt = 0; ielt < section->n_elt; ielt++) {
//          for (int i_face = 0; i_face < n_face_elt; i_face++) {
//            _current_elt_face_vtx_idx[ielt * n_face_elt + i_face + 1] =
//              _current_elt_face_vtx_idx[ielt * n_face_elt + i_face] + 4;
//            _current_elt_face_cell[ielt * n_face_elt + i_face    ] = *n_elt_current + ielt + 1;
//          }

//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 0]  = section->_connec[n_sum_vtx_elt * ielt + 3];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 1]  = section->_connec[n_sum_vtx_elt * ielt + 2];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 2]  = section->_connec[n_sum_vtx_elt * ielt + 1];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 3]  = section->_connec[n_sum_vtx_elt * ielt    ];

//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 4]  = section->_connec[n_sum_vtx_elt * ielt + 6];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 5]  = section->_connec[n_sum_vtx_elt * ielt + 7];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 6]  = section->_connec[n_sum_vtx_elt * ielt + 4];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 7]  = section->_connec[n_sum_vtx_elt * ielt + 5];

//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 8]  = section->_connec[n_sum_vtx_elt * ielt + 4];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 9]  = section->_connec[n_sum_vtx_elt * ielt + 7];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 10] = section->_connec[n_sum_vtx_elt * ielt + 3];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 11] = section->_connec[n_sum_vtx_elt * ielt    ];

//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 12] = section->_connec[n_sum_vtx_elt * ielt + 7];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 13] = section->_connec[n_sum_vtx_elt * ielt + 6];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 14] = section->_connec[n_sum_vtx_elt * ielt + 2];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 15] = section->_connec[n_sum_vtx_elt * ielt + 3];

//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 16] = section->_connec[n_sum_vtx_elt * ielt + 2];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 17] = section->_connec[n_sum_vtx_elt * ielt + 6];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 18] = section->_connec[n_sum_vtx_elt * ielt + 5];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 19] = section->_connec[n_sum_vtx_elt * ielt + 1];

//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 20] = section->_connec[n_sum_vtx_elt * ielt + 1];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 21] = section->_connec[n_sum_vtx_elt * ielt + 5];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 22] = section->_connec[n_sum_vtx_elt * ielt + 4];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 23] = section->_connec[n_sum_vtx_elt * ielt + 0];
//        }

//        *n_elt_current  += section->n_elt;
//        *n_face_current += section->n_elt * n_face_elt;

//      }
//      break;

//    case PDM_MESH_NODAL_PYRAMID5:
//      {
//        const int n_face_elt        = 5;
//        const int n_sum_vtx_face    = 1*4 + 4*3;
//        const int n_sum_vtx_elt     = 5;

//        elt_face_vtx_idx[0] = 0;

//        for (int ielt = 0; ielt < section->n_elt; ielt++) {

//          for (int i_face = 0; i_face < n_face_elt; i_face++) {
//            _current_elt_face_cell[ielt * n_face_elt + i_face    ] = *n_elt_current + ielt + 1;
//          }

//          _current_elt_face_vtx_idx[ielt * n_face_elt + 1]  = elt_face_vtx_idx[ielt * n_face_elt    ] + 4;
//          _current_elt_face_vtx_idx[ielt * n_face_elt + 2]  = elt_face_vtx_idx[ielt * n_face_elt + 1] + 3;
//          _current_elt_face_vtx_idx[ielt * n_face_elt + 3]  = elt_face_vtx_idx[ielt * n_face_elt + 2] + 3;
//          _current_elt_face_vtx_idx[ielt * n_face_elt + 4]  = elt_face_vtx_idx[ielt * n_face_elt + 3] + 3;
//          _current_elt_face_vtx_idx[ielt * n_face_elt + 5]  = elt_face_vtx_idx[ielt * n_face_elt + 4] + 3;

//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 0]  = section->_connec[n_sum_vtx_elt * ielt + 3];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 1]  = section->_connec[n_sum_vtx_elt * ielt + 2];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 2]  = section->_connec[n_sum_vtx_elt * ielt + 1];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 3]  = section->_connec[n_sum_vtx_elt * ielt    ];

//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 4]  = section->_connec[n_sum_vtx_elt * ielt + 4];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 5]  = section->_connec[n_sum_vtx_elt * ielt    ];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 6]  = section->_connec[n_sum_vtx_elt * ielt + 1];

//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 7]  = section->_connec[n_sum_vtx_elt * ielt + 4];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 8]  = section->_connec[n_sum_vtx_elt * ielt + 1];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 9]  = section->_connec[n_sum_vtx_elt * ielt + 2];

//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 10] = section->_connec[n_sum_vtx_elt * ielt + 4];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 11] = section->_connec[n_sum_vtx_elt * ielt + 2];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 12] = section->_connec[n_sum_vtx_elt * ielt + 3];

//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 13] = section->_connec[n_sum_vtx_elt * ielt + 4];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 14] = section->_connec[n_sum_vtx_elt * ielt + 3];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 15] = section->_connec[n_sum_vtx_elt * ielt    ];

//        }

//        *n_elt_current  += section->n_elt;
//        *n_face_current += section->n_elt * n_face_elt;

//      }
//      break;

//    case PDM_MESH_NODAL_PRISM6:
//      {
//        const int n_face_elt        = 5;
//        const int n_sum_vtx_face    = 3*4 + 2*3;
//        const int n_sum_vtx_elt     = 6;

//        for (int ielt = 0; ielt < section->n_elt; ielt++) {

//          for (int i_face = 0; i_face < n_face_elt; i_face++) {
//            _current_elt_face_cell[ielt * n_face_elt + i_face    ] = *n_elt_current + ielt + 1;
//          }

//          _current_elt_face_vtx_idx[ielt * n_face_elt + 1]  = elt_face_vtx_idx[ielt * n_face_elt    ] + 3;
//          _current_elt_face_vtx_idx[ielt * n_face_elt + 2]  = elt_face_vtx_idx[ielt * n_face_elt + 1] + 3;
//          _current_elt_face_vtx_idx[ielt * n_face_elt + 3]  = elt_face_vtx_idx[ielt * n_face_elt + 2] + 4;
//          _current_elt_face_vtx_idx[ielt * n_face_elt + 4]  = elt_face_vtx_idx[ielt * n_face_elt + 3] + 4;
//          _current_elt_face_vtx_idx[ielt * n_face_elt + 5]  = elt_face_vtx_idx[ielt * n_face_elt + 4] + 4;

//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 0]  = section->_connec[n_sum_vtx_elt * ielt + 2];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 1]  = section->_connec[n_sum_vtx_elt * ielt + 1];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 2]  = section->_connec[n_sum_vtx_elt * ielt    ];

//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 3]  = section->_connec[n_sum_vtx_elt * ielt + 4];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 4]  = section->_connec[n_sum_vtx_elt * ielt + 5];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 5]  = section->_connec[n_sum_vtx_elt * ielt + 3];

//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 6]  = section->_connec[n_sum_vtx_elt * ielt + 5];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 7]  = section->_connec[n_sum_vtx_elt * ielt + 4];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 8]  = section->_connec[n_sum_vtx_elt * ielt + 1];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 9]  = section->_connec[n_sum_vtx_elt * ielt + 2];

//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 10] = section->_connec[n_sum_vtx_elt * ielt + 4];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 11] = section->_connec[n_sum_vtx_elt * ielt + 3];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 12] = section->_connec[n_sum_vtx_elt * ielt    ];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 13] = section->_connec[n_sum_vtx_elt * ielt + 1];

//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 14] = section->_connec[n_sum_vtx_elt * ielt + 3];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 15] = section->_connec[n_sum_vtx_elt * ielt + 5];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 16] = section->_connec[n_sum_vtx_elt * ielt + 2];
//          _current_elt_face_vtx[n_sum_vtx_face * ielt + 17] = section->_connec[n_sum_vtx_elt * ielt    ];

//        }

//        *n_elt_current  += section->n_elt;
//        *n_face_current += section->n_elt * n_face_elt;

//      }
//      break;
//    }
//  }

//  else if (id_section < PDM_BLOCK_ID_BLOCK_POLY3D) {

//    int _id = id_section - PDM_BLOCK_ID_BLOCK_POLY2D;
//    PDM_DMesh_nodal_section_poly2d_t *section =
//            (PDM_DMesh_nodal_section_poly2d_t *) PDM_Handles_get (mesh->sections_poly2d, _id);

//    if (section == NULL) {
//      PDM_error (__FILE__, __LINE__, 0, "Bad section identifier\n");
//    }

//    int idx = 0;

//    for (int ielt = 0; ielt < section->n_elt; ielt++) {
//      int n_face_elt = section->_connec_idx[section->n_elt];
//      *n_face_current += n_face_elt;
//      int idx2 = section->_connec_idx[ielt];
//      for (int i_face = 0; i_face < n_face_elt; i_face++) {
//        _current_elt_face_vtx_idx[idx + 1]  = _current_elt_face_vtx_idx[idx] + 2;
//        int inext = (i_face + 1) % n_face_elt;
//        _current_elt_face_vtx[2 * idx    ]  = section->_connec[idx2 + i_face];
//        _current_elt_face_vtx[2 * idx + 1]  = section->_connec[idx2 + inext];
//        _current_elt_face_cell[idx   ]  = *n_elt_current + ielt + 1;
//        idx += 1;
//      }
//    }

//    *n_elt_current += section->n_elt;

//  }

//  else {

//    PDM_error (__FILE__, __LINE__, 0, "PDM_BLOCK_ID_BLOCK_POLY3D : Not implemented yet\n");

//    //TODO: Compliqué car il faut faire des échanges Block_to_part pour recuperer les
//    // definitions des faces
//    // Il faut redupliquer les faces et les stocker comme pour les autres types de
//    // section
//    // Attention : Dans le cas d'un maillage avec une seule section poly3d, il ne faut
//    // rien faire.et ne pas passer dans cette fonction

//  }

//}


/*=============================================================================
 * Public function definitions
 *============================================================================*/

/**
 * \brief Create a Mesh nodal structure
 *
 * \param [in]   n_part   Number of partition on the current process
 * \param [in]   comm     MPI communicator
 *
 * \return       New mesh nodal handle
 *
 */

int
PDM_DMesh_nodal_create
(
const PDM_MPI_Comm comm,
      PDM_g_num_t  n_vtx,
      PDM_g_num_t  n_cell
)
{
  PDM_DMesh_nodal_t *mesh = (PDM_DMesh_nodal_t *) malloc (sizeof(PDM_DMesh_nodal_t));

   if (sizeof(int) != sizeof(PDM_g_num_t)) {


    printf("PDM_DMesh_nodal : Erreur : Cette fontion ne fonctionne pas en 64bit\n");
    exit(1);

  }

  _mesh_init (mesh, comm, n_vtx, n_cell);

  if (mesh_handles == NULL) {
    mesh_handles = PDM_Handles_create (4);
  }

  return PDM_Handles_store (mesh_handles, mesh);
}


/**
 * \brief Free a nodal mesh structure
 *
 * \param [in]  idx   Nodal mesh handle
 *
 * \return      NULL
 *
 */

void
PDM_DMesh_nodal_free
(
const int hdl
)
{

  PDM_DMesh_nodal_t * mesh = (PDM_DMesh_nodal_t *) PDM_Handles_get (mesh_handles, hdl);

  if (mesh != NULL) {

    if (mesh->sections_id != NULL) {
      free (mesh->sections_id);
    }

    mesh->sections_id = NULL;

     _vtx_free (mesh->vtx);

    /* free standard sections */

    if (mesh->sections_std != NULL) {
      int n_sections_std = PDM_Handles_n_get (mesh->sections_std);
      const int *list_ind = PDM_Handles_idx_get (mesh->sections_std);

      while (n_sections_std > 0) {
        PDM_DMesh_nodal_section_std_t *_bloc_std =
          (PDM_DMesh_nodal_section_std_t *) PDM_Handles_get (mesh->sections_std, list_ind[0]);
        _section_std_free(_bloc_std);
        PDM_Handles_handle_free (mesh->sections_std, list_ind[0], PDM_FALSE);
        n_sections_std = PDM_Handles_n_get (mesh->sections_std);
      }

      mesh->sections_std = PDM_Handles_free (mesh->sections_std);
    }

    /* Free polygon sections */

    if (mesh->sections_poly2d != NULL) {
      int n_sections_poly2d = PDM_Handles_n_get (mesh->sections_poly2d);
      const int *list_ind = PDM_Handles_idx_get (mesh->sections_poly2d);

      while (n_sections_poly2d > 0) {
        PDM_DMesh_nodal_section_poly2d_t *_bloc_poly2d =
          (PDM_DMesh_nodal_section_poly2d_t *) PDM_Handles_get (mesh->sections_poly2d, list_ind[0]);
        _section_poly2d_free(_bloc_poly2d);
        PDM_Handles_handle_free (mesh->sections_poly2d, list_ind[0], PDM_FALSE);
        n_sections_poly2d = PDM_Handles_n_get (mesh->sections_poly2d);
      }

      mesh->sections_poly2d = PDM_Handles_free (mesh->sections_poly2d);
    }

    /* Free polyhedron sections */

    if (mesh->sections_poly3d != NULL) {
      int n_sections_poly3d = PDM_Handles_n_get (mesh->sections_poly3d);
      const int *list_ind = PDM_Handles_idx_get (mesh->sections_poly3d);

      while (n_sections_poly3d > 0) {
        PDM_DMesh_nodal_section_poly3d_t *_bloc_poly3d =
          (PDM_DMesh_nodal_section_poly3d_t *) PDM_Handles_get (mesh->sections_poly3d, list_ind[0]);
        _section_poly3d_free(_bloc_poly3d);
        PDM_Handles_handle_free (mesh->sections_poly3d, list_ind[0], PDM_FALSE);
        n_sections_poly3d = PDM_Handles_n_get (mesh->sections_poly3d);
      }

      mesh->sections_poly3d = PDM_Handles_free (mesh->sections_poly3d);
    }

    if (mesh->sections_id != NULL) {
      free (mesh->sections_id);
    }

    if (mesh->dcell_face_idx != NULL) {
      free (mesh->dcell_face_idx);
    }

    if (mesh->dcell_face != NULL) {
      free (mesh->dcell_face);
    }

    if (mesh->cell_distrib != NULL) {
      free (mesh->cell_distrib);
    }

    if (mesh->_dface_vtx_idx != NULL) {
      free (mesh->_dface_vtx_idx);
    }

    if (mesh->_dface_vtx != NULL) {
      free (mesh->_dface_vtx);
    }

    if (mesh->face_distrib != NULL) {
      free (mesh->face_distrib);
    }

     free(mesh);

    PDM_Handles_handle_free (mesh_handles, hdl, PDM_FALSE);

    int n_mesh_array = PDM_Handles_n_get (mesh_handles);

    if (n_mesh_array == 0) {
      mesh_handles = PDM_Handles_free (mesh_handles);
    }
  }
}


/**
 * \brief Define partition vertices
 *
 * \param [in]  hdl       Distributed nodal mesh handle
 * \param [in]  n_vtx     Number of vertices
 * \param [in]  coords    Interlaced coordinates (size = 3 * \ref n_vtx)
 *
 */

void
PDM_DMesh_nodal_coord_set
(
 const int          hdl,
 const int          n_vtx,
 const PDM_real_t  *coords
)
{
  PDM_DMesh_nodal_t * mesh = (PDM_DMesh_nodal_t *) PDM_Handles_get (mesh_handles, hdl);

  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }

  PDM_DMesh_nodal_vtx_t *vtx = mesh->vtx;

  if (vtx->_coords != NULL) {
    PDM_error(__FILE__, __LINE__, 0, "Vertices are already defined\n");
  }

  /* Mapping memoire */

  vtx->n_vtx   = n_vtx;
  vtx->_coords = coords;

  vtx->distrib = (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) * (mesh->n_proc + 1));

  PDM_g_num_t *_distrib = vtx->distrib + 1;
  _distrib[0] = 0;
  PDM_g_num_t _n_vtx = n_vtx;

  PDM_MPI_Scan (&_n_vtx, _distrib, 1, PDM__PDM_MPI_G_NUM,
                PDM_MPI_SUM, mesh->pdm_mpi_comm);


}


/**
 * \brief  Return number of vertices
 *
 * \param [in]  hdl       Distributed nodal mesh handle
 *
 * \return  Number of vertices
 *
 */

int
PDM_DMesh_nodal_n_vtx_get
(
 const int          hdl
)
{
  PDM_DMesh_nodal_t * mesh =
          (PDM_DMesh_nodal_t *) PDM_Handles_get (mesh_handles, hdl);

  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }

  PDM_DMesh_nodal_vtx_t *vtx = mesh->vtx;

  return vtx->n_vtx;
}


/**
 * \brief  Return coordinates of vertices
 *
 * \param [in]  hdl       Distributed nodal mesh handle
 *
 * \return  Coordinates of vertices
 *
 */

const double *
PDM_DMesh_nodal_vtx_get
(
 const int          hdl
)
{
  PDM_DMesh_nodal_t * mesh =
          (PDM_DMesh_nodal_t *) PDM_Handles_get (mesh_handles, hdl);

  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }

  PDM_DMesh_nodal_vtx_t *vtx = mesh->vtx;

  return vtx->_coords;
}


/**
 * \brief  Return number of sections
 *
 * \param [in]  hdl            Nodal mesh handle
 *
 * \return  Number of sections
 *
 */

int
PDM_DMesh_nodal_n_sections_get
(
 const int   hdl
)
{
  PDM_DMesh_nodal_t * mesh =
          (PDM_DMesh_nodal_t *) PDM_Handles_get (mesh_handles, hdl);

  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }


  return mesh->n_sections;

}


/**
 * \brief  Return sections identifier
 *
 * \param [in]  hdl            Distributed nodal mesh handle
 *
 * \return  Blocks identifier
 *
 */

int *
PDM_DMesh_nodal_sections_id_get
(
const int   hdl
)
{
  PDM_DMesh_nodal_t * mesh =
          (PDM_DMesh_nodal_t *) PDM_Handles_get (mesh_handles, hdl);

  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }

  return mesh->sections_id;

}


/**
 * \brief  Return type of section
 *
 * \param [in]  hdl        Distributed nodal mesh handle
 * \param [in]  id_section   Block identifier
 *
 * \return  Type of section
 *
 */

PDM_Mesh_nodal_elt_t
PDM_DMesh_nodal_section_type_get
(
const int   hdl,
const int   id_section
)
{
  PDM_DMesh_nodal_t *mesh =
          (PDM_DMesh_nodal_t *) PDM_Handles_get (mesh_handles, hdl);

  PDM_Mesh_nodal_elt_t t_elt;

  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }

  if (id_section < PDM_BLOCK_ID_BLOCK_POLY2D) {

    t_elt = PDM_MESH_NODAL_POLY_3D;
    PDM_DMesh_nodal_section_std_t *section =
            (PDM_DMesh_nodal_section_std_t *) PDM_Handles_get (mesh->sections_std, id_section);

    if (section == NULL) {
      PDM_error (__FILE__, __LINE__, 0, "Bad section identifier\n");
    }

    t_elt = section->t_elt;
  }

  else if (id_section < PDM_BLOCK_ID_BLOCK_POLY3D) {

    t_elt = PDM_MESH_NODAL_POLY_2D;

  }

  else {

    t_elt = PDM_MESH_NODAL_POLY_3D;

  }

  return t_elt;

}


/**
 * \brief  Add a new section to the current mesh
 *
 * \param [in]  hdl            Distributed nodal mesh handle
 * \param [in]  st_free_data   Status of Release of the memory
 *                             when the section is destroyed
 * \param [in]  id_section       Block identifier
 *
 * \return Block identifier
 *
 */

int
PDM_DMesh_nodal_section_add
(
const int                    hdl,
const PDM_Mesh_nodal_elt_t   t_elt
)
{

  PDM_DMesh_nodal_t *mesh =
          (PDM_DMesh_nodal_t *) PDM_Handles_get (mesh_handles, hdl);

  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }

  int id_section = -1;

  switch (t_elt) {

  case PDM_MESH_NODAL_POINT    :
  case PDM_MESH_NODAL_BAR2     :
  case PDM_MESH_NODAL_TRIA3    :
  case PDM_MESH_NODAL_QUAD4    :
  case PDM_MESH_NODAL_TETRA4   :
  case PDM_MESH_NODAL_PYRAMID5 :
  case PDM_MESH_NODAL_PRISM6   :
  case PDM_MESH_NODAL_HEXA8    :
    {
      /* Mise a jour du tableau de stockage */

      if (mesh->sections_std == NULL) {
        mesh->sections_std = PDM_Handles_create (4);
      }

      /* Allocation du bloc */

      PDM_DMesh_nodal_section_std_t *section_std = (PDM_DMesh_nodal_section_std_t *) malloc(sizeof(PDM_DMesh_nodal_section_std_t));

      id_section = PDM_Handles_store (mesh->sections_std, section_std);

      /* Intialisation du bloc */

      section_std->t_elt = t_elt;

      section_std->n_elt      = -1;
      section_std->_connec    = NULL;
      section_std->distrib    = NULL;

      id_section += PDM_BLOCK_ID_BLOCK_STD;
      if (id_section >= PDM_BLOCK_ID_BLOCK_POLY2D) {
        PDM_error(__FILE__, __LINE__, 0, "The number of standard sections must be less than %d\n",
               PDM_BLOCK_ID_BLOCK_POLY2D);
        abort();
      }
    }

    break;

  case PDM_MESH_NODAL_POLY_2D  :
    {
      /* Mise a jour du tableau de stockage */

      if (mesh->sections_poly2d == NULL) {
        mesh->sections_poly2d = PDM_Handles_create (4);
      }

      /* Allocation du bloc */

      PDM_DMesh_nodal_section_poly2d_t *section_poly2d =
              (PDM_DMesh_nodal_section_poly2d_t *) malloc(sizeof(PDM_DMesh_nodal_section_poly2d_t));

      id_section = PDM_Handles_store (mesh->sections_poly2d, section_poly2d);

      /* Intialisation du bloc */

      section_poly2d->n_elt       = -1;
      section_poly2d->_connec_idx = NULL;
      section_poly2d->_connec     = NULL;
      section_poly2d->distrib     = NULL;

      id_section += PDM_BLOCK_ID_BLOCK_POLY2D;
      if (id_section >= PDM_BLOCK_ID_BLOCK_POLY3D) {
        PDM_error(__FILE__, __LINE__, 0, "The number of polygon sections must be less than %d\n",
               PDM_BLOCK_ID_BLOCK_POLY3D - PDM_BLOCK_ID_BLOCK_POLY2D);
      }
    }

    break;

  case PDM_MESH_NODAL_POLY_3D  :
    {
      /* Mise a jour du tableau de stockage */

      if (mesh->sections_poly3d == NULL) {
        mesh->sections_poly3d = PDM_Handles_create (4);
      }

      /* Allocation du bloc */

      PDM_DMesh_nodal_section_poly3d_t *section_poly3d =
              (PDM_DMesh_nodal_section_poly3d_t *) malloc(sizeof(PDM_DMesh_nodal_section_poly3d_t));

      id_section = PDM_Handles_store (mesh->sections_poly3d, section_poly3d);

      /* Intialisation du bloc */

      section_poly3d->n_elt        = -1;
      section_poly3d->n_face       = -1;
      section_poly3d->_facvtx_idx  = NULL;
      section_poly3d->_facvtx      = NULL;
      section_poly3d->_cellfac_idx = NULL;
      section_poly3d->_cellfac     = NULL;
      section_poly3d->distrib      = NULL;

      id_section += PDM_BLOCK_ID_BLOCK_POLY3D;

    }

    break;

  default :
    PDM_error(__FILE__, __LINE__, 0, "Unknown element type\n");
    break;

  }

  _update_sections_id (mesh);
  return id_section ;

}


/**
 * \brief Define a standard section
 *
 *  - PDM_MESH_NODAL_POINT :
 *
 *   1 x
 *
 *  - PDM_MESH_NODAL_BAR2 :
 *
 *   1 x-------x 2
 *
 *  - PDM_MESH_NODAL_TRIA3 :
 *
 *   1 x-------x 3
 *      \     /
 *       \   /
 *        \ /
 *         x 2
 *
 *  - PDM_MESH_NODAL_QUAD4 :
 *
 *      4 x-------x 3
 *       /       /
 *      /       /
 *   1 x-------x2
 *
 *   - PDM_MESH_NODAL_TETRA4 :
 *
 *         x 4
 *        /|\
 *       / | \
 *      /  |  \
 *   1 x- -|- -x 3
 *      \  |  /
 *       \ | /
 *        \|/
 *         x 2
 *
 *   - PDM_MESH_NODAL_PYRAMID5 :
 *
 *          5 x
 *           /|\
 *          //| \
 *         // |  \
 *      4 x/--|---x 3
 *       //   |  /
 *      //    | /
 *   1 x-------x 2
 *
 *  - PDM_MESH_NODAL_PRSIM6 :
 *
 *   4 x-------x 6
 *     |\     /|
 *     | \   / |
 *   1 x- \-/ -x 3
 *      \ 5x  /
 *       \ | /
 *        \|/
 *         x 2
 *
 *  - PDM_MESH_NODAL_HEXA8 :
 *
 *      8 x-------x 7
 *       /|      /|
 *      / |     / |
 *   5 x-------x6 |
 *     | 4x----|--x 3
 *     | /     | /
 *     |/      |/
 *   1 x-------x 2
 *
 * \param [in]  hdl            Distributed nodal mesh handle
 * \param [in]  id_section       Block identifier
 * \param [in]  n_elt          Number of elements
 * \param [in]  connect        Connectivity
 *
 */

void
PDM_DMesh_nodal_section_std_set
(
const int          hdl,
const int          id_section,
const int          n_elt,
      PDM_g_num_t *connec
)
{
  PDM_DMesh_nodal_t *mesh = (PDM_DMesh_nodal_t *) PDM_Handles_get (mesh_handles, hdl);

  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }

  int _id_section = id_section - PDM_BLOCK_ID_BLOCK_STD;

  PDM_printf(" PDM_Mesh_nodal_section_std_set - _id_section : %i  \n ", _id_section);
  PDM_printf(" PDM_Mesh_nodal_section_std_set - n_elt       : %i  \n ", n_elt);

  PDM_DMesh_nodal_section_std_t *section = (PDM_DMesh_nodal_section_std_t *)
     PDM_Handles_get (mesh->sections_std, _id_section);

  if (section == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad standard section identifier\n");
  }

  /* Mapping */

  section->n_elt   = n_elt;
  section->_connec = connec;

  section->distrib = (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) * (mesh->n_proc + 1));

  /* Creation of distribution */
  PDM_g_num_t begNumAbs;
  PDM_g_num_t _n_elt = n_elt;

  PDM_MPI_Scan (&_n_elt, &begNumAbs, 1, PDM__PDM_MPI_G_NUM,
                PDM_MPI_SUM, mesh->pdm_mpi_comm);
  begNumAbs -= _n_elt;

  PDM_MPI_Allgather((void *) &_n_elt,
                    1,
                    PDM__PDM_MPI_G_NUM,
                    (void *) (&section->distrib[1]),
                    1,
                    PDM__PDM_MPI_G_NUM,
                    mesh->pdm_mpi_comm);

  section->distrib[0] = 0;
  for (int i = 1; i < mesh->n_proc + 1; i++) {
    section->distrib[i] +=  section->distrib[i-1];
  }

}


/**
 * \brief Return standard section description
 *
 *  - PDM_MESH_NODAL_POINT :
 *
 *   1 x
 *
 *  - PDM_MESH_NODAL_BAR2 :
 *
 *   1 x-------x 2
 *
 *  - PDM_MESH_NODAL_TRIA3 :
 *
 *   1 x-------x 3
 *      \     /
 *       \   /
 *        \ /
 *         x 2
 *
 *  - PDM_MESH_NODAL_QUAD4 :
 *
 *      4 x-------x 3
 *       /       /
 *      /       /
 *   1 x-------x2
 *
 *   - PDM_MESH_NODAL_TETRA4 :
 *
 *         x 4
 *        /|\
 *       / | \
 *      /  |  \
 *   1 x- -|- -x 3
 *      \  |  /
 *       \ | /
 *        \|/
 *         x 2
 *
 *   - PDM_MESH_NODAL_PYRAMID5 :
 *
 *          5 x
 *           /|\
 *          //| \
 *         // |  \
 *      4 x/--|---x 3
 *       //   |  /
 *      //    | /
 *   1 x-------x 2
 *
 *  - PDM_MESH_NODAL_PRSIM6 :
 *
 *   4 x-------x 6
 *     |\     /|
 *     | \   / |
 *   1 x- \-/ -x 3
 *      \ 5x  /
 *       \ | /
 *        \|/
 *         x 2
 *
 *  - PDM_MESH_NODAL_HEXA8 :
 *
 *      8 x-------x 7
 *       /|      /|
 *      / |     / |
 *   5 x-------x6 |
 *     | 4x----|--x 3
 *     | /     | /
 *     |/      |/
 *   1 x-------x 2
 *
 * \param [in]  hdl            Distributed nodal mesh handle
 * \param [in]  id_section       Block identifier
 *
 * \return  connect           Connectivity
 *
 */

PDM_g_num_t *
PDM_DMesh_nodal_section_std_get
(
const int            hdl,
const int            id_section
)
{
  PDM_DMesh_nodal_t *mesh =
          (PDM_DMesh_nodal_t *) PDM_Handles_get (mesh_handles, hdl);

  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }

  int _id_section = id_section - PDM_BLOCK_ID_BLOCK_STD;

  const PDM_DMesh_nodal_section_std_t *section =
  (const PDM_DMesh_nodal_section_std_t *) PDM_Handles_get (mesh->sections_std, _id_section);

  if (section == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad standard section identifier\n");
  }

  return section->_connec;
}


/**
 * \brief Get number of section elements
 *
 * \param [in]  idx            Nodal mesh handle
 * \param [in]  id_section       Block identifier
 * \param [in]  id_part        Partition identifier
 *
 * \return      Number of elements
 *
 */


/**
 * \brief Get number of section elements
 *
 * \param [in]  hdl            Distributed nodal mesh handle
 * \param [in]  id_section       Block identifier
 *
 * \return      Number of elements
 *
 */

int
PDM_DMesh_nodal_section_n_elt_get
(
const int            hdl,
const int            id_section
)
{
  PDM_DMesh_nodal_t *mesh =
          (PDM_DMesh_nodal_t *) PDM_Handles_get (mesh_handles, hdl);

  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }

  int _id_section;

  if (id_section >= PDM_BLOCK_ID_BLOCK_POLY3D) {

    _id_section = id_section - PDM_BLOCK_ID_BLOCK_POLY3D;

    const PDM_DMesh_nodal_section_poly3d_t *section =
    (const PDM_DMesh_nodal_section_poly3d_t *)
     PDM_Handles_get (mesh->sections_poly3d, _id_section);

    if (section == NULL) {
      PDM_error (__FILE__, __LINE__, 0, "Bad standard section identifier\n");
    }

    return section->n_elt;
  }

  else if (id_section >= PDM_BLOCK_ID_BLOCK_POLY2D) {

    _id_section = id_section - PDM_BLOCK_ID_BLOCK_POLY2D;

    const PDM_DMesh_nodal_section_poly2d_t *section =
    (const PDM_DMesh_nodal_section_poly2d_t *)
     PDM_Handles_get (mesh->sections_poly2d, _id_section);

    if (section == NULL) {
      PDM_error (__FILE__, __LINE__, 0, "Bad polygon section identifier\n");
    }

    return section->n_elt;
  }

  else {

    _id_section = id_section - PDM_BLOCK_ID_BLOCK_STD;

    const PDM_DMesh_nodal_section_std_t *section =
    (const PDM_DMesh_nodal_section_std_t *)
     PDM_Handles_get (mesh->sections_std, _id_section);

    if (section == NULL) {
      PDM_error (__FILE__, __LINE__, 0, "Bad polyhedron section identifier\n");
    }

    return section->n_elt;
  }

}


/**
 * \brief Define a polygon section
 *
 * \param [in]  hdl            Distributed nodal mesh handle
 * \param [in]  id_section       Block identifier
 * \param [in]  n_elt          Number of elements
 * \param [in]  connect_idx    Connectivity index (size = \ref n_elt + 1)
 * \param [in]  connect        Connectivity (size = \ref connect_idx[\ref n_elt])
 *
 */

void
PDM_DMesh_nodal_section_poly2d_set
(
const int            hdl,
const int            id_section,
const PDM_l_num_t    n_elt,
      PDM_l_num_t   *connec_idx,
      PDM_g_num_t   *connec
)
{
  PDM_DMesh_nodal_t *mesh =
          (PDM_DMesh_nodal_t *) PDM_Handles_get (mesh_handles, hdl);

  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }

  int _id_section = id_section - PDM_BLOCK_ID_BLOCK_POLY2D;

  PDM_DMesh_nodal_section_poly2d_t *section =
          (PDM_DMesh_nodal_section_poly2d_t *) PDM_Handles_get (mesh->sections_poly2d, _id_section);

  if (section == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad standard section identifier\n");
  }

  /* Mapping */

  section->n_elt       = n_elt;
  section->_connec_idx = connec_idx;
  section->_connec     = connec;

  section->distrib = (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) * (mesh->n_proc + 1));

  /* Creation of distribution */
  PDM_g_num_t begNumAbs;
  PDM_g_num_t _n_elt = n_elt;

  PDM_MPI_Scan (&_n_elt, &begNumAbs, 1, PDM__PDM_MPI_G_NUM,
                PDM_MPI_SUM, mesh->pdm_mpi_comm);
  begNumAbs -= _n_elt;

  PDM_MPI_Allgather((void *) &_n_elt,
                    1,
                    PDM__PDM_MPI_G_NUM,
                    (void *) (&section->distrib[1]),
                    1,
                    PDM__PDM_MPI_G_NUM,
                    mesh->pdm_mpi_comm);

  section->distrib[0] = 0;
  for (int i = 1; i < mesh->n_proc + 1; i++) {
    section->distrib[i] +=  section->distrib[i-1];
  }

}


/**
 * \brief Return a polygon section description
 *
 * \param [in]  hdl            Distributed nodal mesh handle
 * \param [in]  id_section       Block identifier
 * \param [out] connect_idx    Connectivity index (size = \ref n_elt + 1)
 * \param [out] connect        Connectivity (size = \ref connect_idx[\ref n_elt])
 *
 */

void
PDM_DMesh_nodal_section_poly2d_get
(
 const int          hdl,
 const int          id_section,
       PDM_l_num_t  **connec_idx,
       PDM_g_num_t  **connec
)
{
  PDM_DMesh_nodal_t *mesh =
          (PDM_DMesh_nodal_t *) PDM_Handles_get (mesh_handles, hdl);

  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }

  int _id_section = id_section - PDM_BLOCK_ID_BLOCK_POLY2D;

  PDM_DMesh_nodal_section_poly2d_t *section =
          (PDM_DMesh_nodal_section_poly2d_t *) PDM_Handles_get (mesh->sections_poly2d, _id_section);

  if (section == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad standard section identifier\n");
  }

  *connec_idx = section->_connec_idx;
  *connec     = section->_connec;

}


/**
 * \brief Define a polyhedra section
 *
 * \param [in]  hdl            Distributed nodal mesh handle
 * \param [in]  id_section       Block identifier
 * \param [in]  n_elt          Number of polyhedra
 * \param [in]  n_face         Number of faces used to describe polyhedra
 * \param [in]  facvtx_idx     Index of face vertex connectivity
 * \param [in]  facvtx         Face vertex connectivity
 * \param [in]  cellfac_idx    Index of cell face connectivity
 * \param [in]  cellfac        Cell face connectivity
 *
 */

void
PDM_DMesh_nodal_section_poly3d_set
(
const int            hdl,
const int            id_section,
const PDM_l_num_t    n_elt,
const PDM_l_num_t    n_face,
      PDM_l_num_t   *facvtx_idx,
      PDM_g_num_t   *facvtx,
      PDM_l_num_t   *cellfac_idx,
      PDM_g_num_t   *cellfac
)
{
  PDM_DMesh_nodal_t *mesh = (PDM_DMesh_nodal_t *) PDM_Handles_get (mesh_handles, hdl);

  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }

  int _id_section = id_section - PDM_BLOCK_ID_BLOCK_POLY2D;

  PDM_DMesh_nodal_section_poly3d_t *section =
          (PDM_DMesh_nodal_section_poly3d_t *) PDM_Handles_get (mesh->sections_poly3d, _id_section);

  if (section == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad standard section identifier\n");
  }

  section->n_elt        = n_elt;
  section->n_face       = n_face;
  section->_facvtx_idx  = facvtx_idx;
  section->_facvtx      = facvtx;
  section->_cellfac_idx = cellfac_idx;
  section->_cellfac     = cellfac;

  section->distrib = (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) * (mesh->n_proc + 1));

  /* Creation of distribution */
  PDM_g_num_t begNumAbs;
  PDM_g_num_t _n_elt = n_elt;

  PDM_MPI_Scan (&_n_elt, &begNumAbs, 1, PDM__PDM_MPI_G_NUM,
                PDM_MPI_SUM, mesh->pdm_mpi_comm);
  begNumAbs -= _n_elt;

  PDM_MPI_Allgather((void *) &_n_elt,
                    1,
                    PDM__PDM_MPI_G_NUM,
                    (void *) (&section->distrib[1]),
                    1,
                    PDM__PDM_MPI_G_NUM,
                    mesh->pdm_mpi_comm);

  section->distrib[0] = 0;
  for (int i = 1; i < mesh->n_proc + 1; i++) {
    section->distrib[i] +=  section->distrib[i-1];
  }


}


/**
 * \brief Define a polyhedra section
 *
 * \param [in]  hdl            Distributed nodal mesh handle
 * \param [in]  id_section       Block identifier
 * \param [out]  n_face         Number of faces used to describe polyhedra
 * \param [out]  facvtx_idx     Index of face vertex connectivity
 * \param [out]  facvtx         Face vertex connectivity
 * \param [out]  cellfac_idx    Index of cell face connectivity
 * \param [out]  cellfac        Cell face connectivity
 *
 */

void
PDM_DMesh_nodal_section_poly3d_get
(
const int            hdl,
const int            id_section,
      PDM_l_num_t   *n_face,
      PDM_l_num_t  **facvtx_idx,
      PDM_g_num_t  **facvtx,
      PDM_l_num_t  **cellfac_idx,
      PDM_g_num_t  **cellfac
)
{
  PDM_DMesh_nodal_t *mesh =
          (PDM_DMesh_nodal_t *) PDM_Handles_get (mesh_handles, hdl);

  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }

  int _id_section = id_section - PDM_BLOCK_ID_BLOCK_POLY2D;

  PDM_DMesh_nodal_section_poly3d_t *section =
          (PDM_DMesh_nodal_section_poly3d_t *) PDM_Handles_get (mesh->sections_poly3d, _id_section);

  if (section == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad standard section identifier\n");
  }

  *n_face      = section->n_face;
  *facvtx_idx  = section->_facvtx_idx;
  *facvtx      = section->_facvtx;
  *cellfac_idx = section->_cellfac_idx;
  *cellfac     = section->_cellfac;

}


/**
 * \brief  Return total number of elements of a distributed mesh
 *
 * \param [in]  hdl       Distributed nodal mesh handle
 *
 * \return  Return number elements of a partition
 *
 */

PDM_g_num_t
PDM_DMesh_nodal_total_n_cell_get
(
const int  hdl
)
{
  PDM_DMesh_nodal_t * mesh =
          (PDM_DMesh_nodal_t *) PDM_Handles_get (mesh_handles, hdl);

  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }

  PDM_g_num_t total_n_cell = 0;

  int n_sections_std = PDM_Handles_n_get (mesh->sections_std);
  const int *list_ind = PDM_Handles_idx_get (mesh->sections_std);

  for (int i = 0; i < n_sections_std; i++) {
    PDM_DMesh_nodal_section_std_t *_section_std =
      (PDM_DMesh_nodal_section_std_t *) PDM_Handles_get (mesh->sections_std, list_ind[i]);
    total_n_cell += _section_std->distrib[mesh->n_proc];
  }

  int n_sections_poly2d = PDM_Handles_n_get (mesh->sections_poly2d);
  list_ind = PDM_Handles_idx_get (mesh->sections_poly2d);

  for (int i = 0; i < n_sections_poly2d; i++) {
    PDM_DMesh_nodal_section_poly2d_t *_section_poly2d =
      (PDM_DMesh_nodal_section_poly2d_t *) PDM_Handles_get (mesh->sections_poly2d, list_ind[i]);
    total_n_cell += _section_poly2d->distrib[mesh->n_proc];
  }

  int n_sections_poly3d = PDM_Handles_n_get (mesh->sections_poly3d);
  list_ind = PDM_Handles_idx_get (mesh->sections_poly3d);

  for (int i = 0; i < n_sections_poly3d; i++) {
    PDM_DMesh_nodal_section_poly3d_t *_section_poly3d =
      (PDM_DMesh_nodal_section_poly3d_t *) PDM_Handles_get (mesh->sections_poly3d, list_ind[i]);
    total_n_cell += _section_poly3d->distrib[mesh->n_proc];
  }

  return total_n_cell;

}


/**
 * \brief  Return total number of faces of a distributed mesh
 *
 * \param [in]  hdl       Distributed nodal mesh handle
 *
 * \return  Return total number of faces
 *
 */

PDM_g_num_t
PDM_DMesh_nodal_total_n_face_get
(
const int  hdl
)
{
  PDM_DMesh_nodal_t * mesh =
          (PDM_DMesh_nodal_t *) PDM_Handles_get (mesh_handles, hdl);

  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }

  // if (mesh->dcell_face == NULL) {
  //   PDM_error (__FILE__, __LINE__, 0, "Not implemented yet\n");
  // }

  return mesh->face_distrib[mesh->n_proc];

}


/**
 * \brief  Return total number of vertices of a distributed mesh
 *
 * \param [in]  hdl       Distributed nodal mesh handle
 *
 * \return  Return total number of vertices
 *
 */

PDM_g_num_t
PDM_DMesh_nodal_total_n_vtx_get
(
const int  hdl
)
{
  PDM_DMesh_nodal_t * mesh =
          (PDM_DMesh_nodal_t *) PDM_Handles_get (mesh_handles, hdl);

  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }

  PDM_DMesh_nodal_vtx_t *vtx = mesh->vtx;

  return vtx->distrib[mesh->n_proc];
}


/**
 * \brief  Compute cell->face connectivity
 *
 * \param [in]   hdl              Distributed nodal mesh handle
 *
 */

void
PDM_DMesh_nodal_cell_face_compute
(
const int   hdl
)
{


  PDM_printf("PDM_DMesh_nodal_cell_face_compute \n ");

  /* Get current structure to treat */
  PDM_DMesh_nodal_t *mesh = (PDM_DMesh_nodal_t *) PDM_Handles_get (mesh_handles, hdl);


  /* Verbose */
  if(1 == 0)
  {
    PDM_printf(" -------- \n ");
    PDM_printf(" n_proc     : %i \n ", mesh->n_proc);
    PDM_printf(" n_sections : %i \n ", mesh->n_sections);

    PDM_printf(" -------- \n ");

    // TODO : Sure ?? Check with Eric

    for(int i_section=0; i_section < mesh->n_sections; i_section++){

      PDM_DMesh_nodal_section_std_t *section_std = (PDM_DMesh_nodal_section_std_t *) PDM_Handles_get (mesh->sections_std, i_section);

      PDM_printf(" ------------------------------ \n ");
      PDM_printf(" i_section : %i  \n ", i_section);
      PDM_printf(" t_elt    : %i  \n ", section_std->t_elt);
      PDM_printf(" n_elt    : %i  \n ", section_std->n_elt);

      PDM_printf(" DitribElmt ... \n ");
      for(int iProc=0; iProc < mesh->n_proc + 1; iProc++)
        PDM_printf(PDM_FMT_G_NUM, section_std->distrib[iProc]);
      PDM_printf("\n");
    }

  }

  /* Creation of element distribution among all sections */
  mesh->section_distribution    = (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) * (mesh->n_sections + 1));
  mesh->section_distribution[0] = 0;

  for (int i_section = 0; i_section < mesh->n_sections; i_section++) {
    PDM_DMesh_nodal_section_std_t *section_std = (PDM_DMesh_nodal_section_std_t *) PDM_Handles_get (mesh->sections_std, i_section);
    mesh->section_distribution[i_section+1] = section_std->distrib[mesh->n_proc];
  }
  for (int i_section = 1; i_section < mesh->n_sections + 1; i_section++) {
    mesh->section_distribution[i_section] +=  mesh->section_distribution[i_section-1];
  }

  /* Verbose */
  if(1 == 0)
  {
    PDM_printf(" ------------------------------ \n ");
    for(int i_section=0; i_section < mesh->n_sections+1; i_section++){
      PDM_printf("%i ", mesh->section_distribution[i_section]);
    }
  }


  /* Build an approximate number of faces for current processor */
  PDM_g_num_t dnElemTot  = 0;
  for (int i_section = 0; i_section < mesh->n_sections; i_section++) {
    PDM_DMesh_nodal_section_std_t *section_std = (PDM_DMesh_nodal_section_std_t *) PDM_Handles_get (mesh->sections_std, i_section);
    dnElemTot += ( section_std->distrib[mesh->i_proc+1] - section_std->distrib[mesh->i_proc] );
  }
  int n_fac_approx = dnElemTot*8;

  if(1 == 0)
  {
    PDM_printf("n_fac_approx : %i  \n ", n_fac_approx);
  }


  /*
   * Allocate Disctribute Hash Key
   */
  PDM_g_num_t* ln_to_gn = (PDM_g_num_t *) malloc( sizeof(PDM_g_num_t) * n_fac_approx );


  /** Allocate (SurDim the part_data and part_stride ) **/
  int n_data_approx = dnElemTot*6*4*2;
  int* part_data = (int *) malloc( sizeof(int) * n_data_approx );
  int* part_stri = (int *) malloc( sizeof(int) * n_fac_approx  );

  /** Initialisation of some data **/
  int n_face  = 0;
  int n_data  = 0;

  /** Loop over Elements **/
  for (int i_section = 0; i_section < mesh->n_sections; i_section++) {

    /* Get current section */
    PDM_DMesh_nodal_section_std_t *section_std = (PDM_DMesh_nodal_section_std_t *) PDM_Handles_get (mesh->sections_std, i_section);

    /** Convert type to stride on elements **/
    int n_vtx_per_elmt = _get_size_of_element(section_std->t_elt);
    int n_face_per_elmt = _get_nbface_per_element(section_std->t_elt);

    /** Get the index wich compute the elements **/
    int*  n_vtx_per_face  = (int * ) malloc( n_face_per_elmt * sizeof(int * ) );
    int** tabFacVtx = (int **) malloc( n_face_per_elmt * sizeof(int **) );

    /** Fill up n_vtx_per_face and tabFacVtx **/
    _get_elmt_info(section_std->t_elt, n_vtx_per_face, tabFacVtx);

    // int BegE = section_std->distrib[mesh->i_proc  ]*n_vtx_per_elmt  ; //+ offset;
    // int EndE = section_std->distrib[mesh->i_proc+1]*n_vtx_per_elmt-1; //+ offset;

    /** Verbose **/
    // if(0 == 0){
    //   printf("BegE/EndE   : %d/%d     \n", BegE, EndE);
    // }

    /** Each section is composed of multiple elements **/
    for(int i_elmt = section_std->distrib[mesh->i_proc]; i_elmt < section_std->distrib[mesh->i_proc+1]; i_elmt++){

      // int i_offset = i_elmt*n_vtx_per_elmt-section_std->distrib[mesh->i_proc];
      int i_offset = (i_elmt-section_std->distrib[mesh->i_proc])*n_vtx_per_elmt;
      // printf("[%i] - i_offset : %i --- %i \n",mesh->i_proc,  i_offset, section_std->distrib[mesh->i_proc]);
      /** For each Elemt we decompose it to face  **/
      for(int i_face = 0; i_face < n_face_per_elmt; i_face++){

        /* Compute the key */
        int ikey = _compute_key(section_std->_connec, tabFacVtx[i_face], i_offset, n_vtx_per_face[i_face]);

        /* Build the ln_to_gn */
        // Attnetion ici c'est faux je pense !!!!
        if(n_face >= n_fac_approx){
          n_fac_approx = 2*n_fac_approx;
          ln_to_gn  = (PDM_g_num_t *) realloc(ln_to_gn   , sizeof(PDM_g_num_t) * n_fac_approx );
          part_stri = (int         *) realloc(part_stri, sizeof(int) * n_fac_approx );
        }
        if(n_data >= n_data_approx){
          n_data_approx = 2*n_data_approx;
          part_data = (int         *) realloc(part_data, sizeof(int) * n_data_approx );
        }

        ln_to_gn[n_face] = ikey;

        /** Prepare part_data - iElem iType n_vtx (i1, i2, i3, ... ) - **/
        part_data[n_data++] = i_elmt+1+mesh->section_distribution[i_section];
        part_data[n_data++] = 1;

        /* Setup n_vtx */
        part_data[n_data++] = n_vtx_per_face[i_face] ;

        for(int i_vtx=0; i_vtx < n_vtx_per_face[i_face]; i_vtx++){
          // printf("lVtxIdx[%d] : %d -> %d \n ", i_vtx, lVtxIdx[i_vtx],dCellVtxArr[lVtxIdx[i_vtx]]);
          // printf("lVtxIdx[%d] : %d -> %d \n ", BegE+tabFacVtx[i_face][i_vtx], fetch->_dElemArr[BegE+tabFacVtx[i_face][i_vtx]], ikey);
          // printf("lVtxIdx[%d] : %d -> %d \n ", i_vtx, section_std->_connec[i_offset+tabFacVtx[i_face][i_vtx]], ikey);
          part_data[n_data++] = section_std->_connec[i_offset+tabFacVtx[i_face][i_vtx]];
        }

        /** Prepare part_data **/
        part_stri[n_face] = 1+1+1+n_vtx_per_face[i_face];

        /** Go to next face **/
        n_face++;

      } /* End face loop */
    } /* End element loop */

    /** Free **/
    free(n_vtx_per_face);
    for(int i=0; i < n_face_per_elmt; i++){
      free(tabFacVtx[i]);
    }
    free(tabFacVtx);

  } /* End section loop */


  /*
   * Re-Allocate Array
   */
  ln_to_gn    = (PDM_g_num_t *) realloc((ln_to_gn   ), n_face * sizeof(PDM_g_num_t * ));
  part_stri = (int *) realloc((part_stri), n_face * sizeof(int * ));
  part_data = (int *) realloc((part_data), n_data * sizeof(int * ));

  /*
   * Verbose
   */
  if(0 == 1){
    printf("------------------------------------------- \n");
    printf("n_face       : %d \n", n_face );
    printf("n_data       : %d \n", n_data );
    printf("n_fac_approx  : %d \n", n_fac_approx );
    printf("------------------------------------------- \n");
    for (int i = 0; i < n_face; i++) {
      printf("[%d] - i/ikey  : %d - "PDM_FMT_G_NUM" \n", mesh->i_proc, i, ln_to_gn[i]);
    }
    printf("------------------------------------------- \n");
  }

  /*
   * Create PartToBlock Structure
   */
  PDM_part_to_block_t *ptb = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                      PDM_PART_TO_BLOCK_POST_MERGE,
                                                      1.,
                                                      &ln_to_gn,
                                                      NULL,
                                                      &n_face,
                                                      1,
                                                      mesh->pdm_mpi_comm);

  /*
   * Exchange data in one call ( better do multiple -> Eric ??)
   */
  int *blk_stri = NULL;
  int *blk_data = NULL;

  int data_size = PDM_part_to_block_exch(          ptb,
                                                  sizeof(PDM_g_num_t),
                                                  PDM_STRIDE_VAR,
                                                  -1,
                                                  &part_stri,
                                        (void **) &part_data,
                                                  &blk_stri,
                                        (void **) &blk_data);

  /*
   *  Get the size of the current process bloc
   */
  int blk_size = PDM_part_to_block_n_elt_block_get(ptb);

  /*
   * Verbose
   */
  if(0 == 1){
    printf("blk_size : %d\n", blk_size);
    for(int i = 0; i < blk_size; i++) {
      // printf("BlockData[%d]    : %d \n", i, blk_data[i]);
      printf("blk_stri[%d] : %d \n", i, blk_stri[i]);
    }
  }

  /*
   *  Creation of array of diplacement
   */
  int* blk_stri_idx = (int *) malloc( sizeof(int *) * (blk_size+1) ) ;
  blk_stri_idx[0] = 0;
  for (int i = 0; i < blk_size; i++) {
    blk_stri_idx[i+1] = blk_stri_idx[i] + blk_stri[i];
  }
  free(blk_stri);

  /*
   * Allocate Memory - face_vtx - face_cell
   */
  mesh->_dface_vtx     = (PDM_g_num_t *) malloc( sizeof(PDM_g_num_t *) * data_size/2); /* Not stupid at all */
  mesh->_dface_vtx_idx  = (int *) malloc( sizeof(int *) * data_size/2); /* Surdim as Hell */
  mesh->_dface_cell    = (PDM_g_num_t *) malloc( sizeof(PDM_g_num_t *) * data_size/2); /* Surdim as Hell */

  /*
   * Init AbsNumerotation
   */
  int i_abs_face = 0;
  mesh->_dface_vtx_idx[0] = 0;

  /*
   * Loop over blk_data : Each block correspond to a key
   *              Each key can identify the same faces
   *              Need to solve all conflict -> See _find_pairs
   */
  for(int i_blk=0; i_blk < blk_size; i_blk++){

    /*
     * Traitement des groupes de données
     * Reminder Blk = [ [iCell, iType, n_vtx, i1, i2, ...], ... ]
     */
    int nEntryFace = blk_stri_idx[i_blk+1] - blk_stri_idx[i_blk];

    /*
     * Verbose
     */
    if(0 == 1){
      printf("blk_stri[%d] : %d -> %d \n", i_blk, blk_stri_idx[i_blk], nEntryFace);
      for(int j=blk_stri_idx[i_blk]; j < blk_stri_idx[i_blk+1]; j++){
        printf("blk_data[%d] : %d \n", j, blk_data[j]);
      }
    }

    /*
     * First Pass - Identify face conflict
     */

    int nFac = 0;
    int i_fac_data = blk_stri_idx[i_blk];
    while(i_fac_data < blk_stri_idx[i_blk+1]){
      // printf("Number of Vextex for Face : %d \n", blk_data[i_fac_data+2]);
      i_fac_data += 3+blk_data[i_fac_data+2];
      nFac += 1;
    }

    /*
     *  Allocate array of displacement for all faces in conflict
     *         -> Set the adress of the block data of current faces
     *         -> Make a surdim before and not realloc ...
     */
    int *idx_face = (int *) malloc( sizeof(int *) * nFac + 1);

    /** New **/
    i_fac_data   = blk_stri_idx[i_blk];
    idx_face[0] = blk_stri_idx[i_blk];
    nFac       = 0;

    /** All Faces **/
    while(i_fac_data < blk_stri_idx[i_blk+1]){
      idx_face[nFac+1] = idx_face[nFac] + blk_data[i_fac_data+2] + 3;
      // printf("blk_data : %d \n ", blk_data[i_fac_data+2]);
      // printf("idx_face[%d] : %d \n ", nFac+1, idx_face[nFac+1]);
      i_fac_data     += blk_data[i_fac_data+2] + 3;
      nFac         += 1;
    }
    /* ************************************************** */


    /*
     * Verbose
     */
    if(0 == 1){
      for(int j=0; j < nFac+1; j++){
        printf("idx_face[%d] : %d \n", j, idx_face[j]);
        // printf("Deltaidx_face[%d] : %d \n", j, idx_face[j+1]-idx_face[j]);
      }
    }

    /*
     * Solve conflict between all faces and build connectivity array
     *   WATCH OUT -> i_abs_face is changed
     */
    i_abs_face = _find_pairs(idx_face,
                           blk_data,
                           nFac,
                           i_abs_face,
                           mesh->_dface_vtx,
                           mesh->_dface_vtx_idx,
                           mesh->_dface_cell);

    /*
     * Free
     */
    free(idx_face);

  }

  /*
   * Fill up fetch structure
   */
  mesh->dn_face = i_abs_face;

  /*
   * Realloc -> TODOUX
   */
  mesh->_dface_vtx_idx = (int *        ) realloc((mesh->_dface_vtx_idx), (mesh->dn_face + 1) * sizeof(int * ) );
  mesh->_dface_vtx    = (PDM_g_num_t *) realloc((mesh->_dface_vtx   ), mesh->_dface_vtx_idx[mesh->dn_face] * sizeof(PDM_g_num_t * ));
  mesh->_dface_cell   = (PDM_g_num_t *) realloc((mesh->_dface_cell  ), mesh->dn_face * 2 * sizeof(PDM_g_num_t * ));

  /*
   * Verbose
   */
  if(1 == 0){
    printf("fetch->dn_face : %d \n", mesh->dn_face);
    for(int i = 0; i < mesh->dn_face; i++) {
      // printf("BlockData[%d]    : %d \n", i, blk_data2[i]);
      printf("mesh->_dface_cell[%d] : "PDM_FMT_G_NUM"/"PDM_FMT_G_NUM" \n", i, mesh->_dface_cell[2*i], mesh->_dface_cell[2*i+1]);
    }

  }

  /*
   * Generate absolute numerotation of faces
   */
  _make_absolute_face_numbering(mesh);

  /*
   * Rebuild cell face
   */
  PDM_g_num_t *ln_to_gn_elem = (PDM_g_num_t *) malloc( sizeof(PDM_g_num_t) * 2 * mesh->dn_face );
  for (PDM_g_num_t i = mesh->face_distrib[mesh->i_proc]; i < mesh->face_distrib[mesh->i_proc+1]; i++) {
    int idx = (int) (i-mesh->face_distrib[mesh->i_proc]);
    ln_to_gn_elem[2*idx  ] = i+1;
    ln_to_gn_elem[2*idx+1] = i+1;
  }

  /*
   * Verbose
   */
  if(0 == 1){
    printf("------------------------------------------- \n");
    for (int i = mesh->face_distrib[mesh->i_proc]; i < mesh->face_distrib[mesh->i_proc+1]; i++) {
      // printf("[%d] - i/i_face  : %d - %d \n", mesh->i_proc, i, ln_to_gn_elem[i]);
      printf("[%d] - i/i_face  : %d - "PDM_FMT_G_NUM" / "PDM_FMT_G_NUM" \n", mesh->i_proc, i, ln_to_gn_elem[2*i], ln_to_gn_elem[2*i+1]);
    }
    printf("------------------------------------------- \n");
  }


  int n_face2 = 2*mesh->dn_face;

  int* part_stri2 = (int *) malloc( sizeof(int) * 2 * mesh->dn_face );

  /* Prepare dface_cell loc */
  PDM_g_num_t* dface_cell_tmp = (PDM_g_num_t *) malloc( sizeof(PDM_g_num_t) * 2 * mesh->dn_face );

  int idx_g = 0;
  for (int i = mesh->face_distrib[mesh->i_proc]; i < mesh->face_distrib[mesh->i_proc+1]; i++) {
    int idx = (int) (i-mesh->face_distrib[mesh->i_proc]);
    dface_cell_tmp[idx_g] = mesh->_dface_cell[2*idx];
    ln_to_gn_elem[idx_g] = i+1;
    part_stri2[idx_g] = 1;
    idx_g++;
    if(mesh->_dface_cell[2*idx+1] != 0){
      dface_cell_tmp[idx_g] = mesh->_dface_cell[2*idx+1];
      ln_to_gn_elem[idx_g  ] = i+1;
      part_stri2[idx_g] = 1;
      idx_g++;
    }
    else
    {
      dface_cell_tmp[idx_g] = mesh->_dface_cell[2*idx];
      ln_to_gn_elem[idx_g  ] = i+1;
      part_stri2[idx_g] = 0;
      idx_g++;
    }
  }

  n_face2 = idx_g;



  PDM_part_to_block_t *ptb2 = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                       PDM_PART_TO_BLOCK_POST_MERGE,
                                                       1.,
                                                       &dface_cell_tmp,
                                                       NULL,
                                                       &n_face2,
                                                       1,
                                                       mesh->pdm_mpi_comm);

  int *blk_stri2 = NULL;
  int *blk_data2 = NULL;

  /** Prepare part_stride for dface_cell - 2 for all **/
  // for(int i = 0; i < 2*mesh->dn_face ; i++){
  // for(int i = 0; i < idx_g ; i++){
  //   part_stri2[i] = 1;
  //   // part_stri2[i] = mesh->dn_face;
  // }

  PDM_part_to_block_exch(          ptb2,
                                   sizeof(PDM_g_num_t),
                                   PDM_STRIDE_VAR,
                                   -1,
                                   &part_stri2,
                         (void **) &ln_to_gn_elem,
                                   &blk_stri2,
                         (void **) &blk_data2);

  /*
   *  Get the size of the current process bloc
   */
  int delmt_tot = PDM_part_to_block_n_elt_block_get(ptb2); /* Volumic and surfacic */
  mesh->n_dcell = delmt_tot;

  /*
   * Verbose
   */
  if(0 == 1){
    printf("delmt_tot : %d\n", delmt_tot);
    for(int i = 0; i < delmt_tot; i++) {
      printf("BlockData[%d]    : %d \n", i, blk_data2[i]);
      printf("blk_stri2[%d] : %d \n", i, blk_stri2[i]);
    }
  }

  mesh->dcell_face_idx = (PDM_l_num_t * ) malloc( (delmt_tot + 1) * sizeof(PDM_l_num_t * ) );

  mesh->dcell_face_idx[0] = 0;
  for(int i = 0; i < delmt_tot; i++){
    mesh->dcell_face_idx[i+1] = mesh->dcell_face_idx[i] + blk_stri2[i];
  }

  // printf("delmt_tot : %d\n", mesh->dcell_face_idx[delmt_tot]);

  mesh->dcell_face = (PDM_g_num_t *) malloc(mesh->dcell_face_idx[delmt_tot] * sizeof(PDM_g_num_t * ));

  for(int i = 0; i < delmt_tot; i++){

    // int n_face_per_elmt = blk_stri2[i];
    int n_face_per_elmt = mesh->dcell_face_idx[i+1] - mesh->dcell_face_idx[i];
    // printf("n_face_per_elmt : %d\n", n_face_per_elmt);
    for(int iFac = 0; iFac < n_face_per_elmt; iFac++){
      mesh->dcell_face[mesh->dcell_face_idx[i]+iFac] = blk_data2[mesh->dcell_face_idx[i]+iFac];
    }
  }

  /*
   * Mask all face_cell with boundary
   */
  // for(int i = 0; i < mesh->dn_face; i++) {
  //    if(mesh->_dface_cell[2*i+1] > mesh->n_cell_abs){
  //       mesh->_dface_cell[2*i+1] = 0;
  //    }
  // }

  /* Free */
  PDM_part_to_block_free(ptb );
  PDM_part_to_block_free(ptb2);

  free(ln_to_gn);
  free(ln_to_gn_elem);
  free(dface_cell_tmp);
  free(part_data);
  free(part_stri);

  free(blk_data   );
  free(blk_stri_idx);
  free(blk_stri2  );
  free(blk_data2  );

  // PDM_error (__FILE__, __LINE__, 0, "Not implemented yet\n");
}


/**
 * \brief  Return cell \rightarrow face connectivity
 *
 * \param [in]   hdl              Distributed nodal mesh handle
 * \param [out]  dcell_face_idx   Index of distributed cell->face connectivity
 * \param [out]  dcell_face       Distributed cell->face connectivity
 *
 * \return     Number of cells on the current process
 *
 */

int
PDM_DMesh_nodal_cell_face_get
(
const int   hdl,
      int   **dcell_face_idx,
PDM_g_num_t **dcell_face
)
{
  PDM_DMesh_nodal_t * mesh =
          (PDM_DMesh_nodal_t *) PDM_Handles_get (mesh_handles, hdl);

  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }

  *dcell_face_idx = mesh->dcell_face_idx;
  *dcell_face = mesh->dcell_face;

  return mesh->n_dcell;

}

/**
 * \brief  Return face->cell connectivity
 *
 * \param [in]   hdl              Distributed nodal mesh handle
 * \param [out]  face_cell       Distributed face->cell connectivity
 *
 * \return     Number of cells on the current process
 *
 */
int
PDM_DMesh_nodal_face_cell_get
(
const int     hdl,
PDM_g_num_t **dface_cell
)
{
  PDM_DMesh_nodal_t * mesh =
          (PDM_DMesh_nodal_t *) PDM_Handles_get (mesh_handles, hdl);

  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }

  *dface_cell = mesh->_dface_cell;

  return mesh->dn_face;

}



/**
 * \brief  Return face \rightarrow vertex connectivity
 *
 * \param [in]   hdl              Distributed nodal mesh handle
 * \param [out]  dcell_face_idx   Index of distributed cell->face connectivity
 * \param [out]  dcell_face       Distributed cell->face connectivity
 *
 * \return     Number of faces on the current process
 *
 */

int
PDM_DMesh_nodal_face_vtx_get
(
const int   hdl,
      int   **_dface_vtx_idx,
PDM_g_num_t **_dface_vtx
)
{
  PDM_DMesh_nodal_t * mesh =
          (PDM_DMesh_nodal_t *) PDM_Handles_get (mesh_handles, hdl);

  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }

  *_dface_vtx_idx = mesh->_dface_vtx_idx;
  *_dface_vtx     = mesh->_dface_vtx;

  return mesh->dn_face;

}

/**
 * \brief  Return vertices distribution
 *
 * \param [in]  hdl  Distributed nodal mesh handle
 *
 * \return  A array of size \ref n_procs + 1
 *
 */

const PDM_g_num_t *
PDM_DMesh_nodal_distrib_vtx_get
(
 const int          hdl
)
{
  PDM_DMesh_nodal_t * mesh =
          (PDM_DMesh_nodal_t *) PDM_Handles_get (mesh_handles, hdl);

  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }

  PDM_DMesh_nodal_vtx_t *vtx = mesh->vtx;

  return vtx->distrib;

}


/**
 * \brief  Return section distribution
 *
 * \param [in]  hdl        Distributed nodal mesh handle
 * \param [in]  id_section   Block identifier
 *
 * \return  A array of size \ref n_procs + 1
 *
 */

const PDM_g_num_t *
PDM_DMesh_nodal_distrib_section_get
(
 const int   hdl,
 const int   id_section
)
{
  PDM_DMesh_nodal_t *mesh =
          (PDM_DMesh_nodal_t *) PDM_Handles_get (mesh_handles, hdl);

  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }

  int _id_section;

  if (id_section >= PDM_BLOCK_ID_BLOCK_POLY3D) {

    _id_section = id_section - PDM_BLOCK_ID_BLOCK_POLY3D;

    const PDM_DMesh_nodal_section_poly3d_t *section =
    (const PDM_DMesh_nodal_section_poly3d_t *)
     PDM_Handles_get (mesh->sections_poly3d, _id_section);

    if (section == NULL) {
      PDM_error (__FILE__, __LINE__, 0, "Bad standard section identifier\n");
    }

    return section->distrib;
  }

  else if (id_section >= PDM_BLOCK_ID_BLOCK_POLY2D) {

    _id_section = id_section - PDM_BLOCK_ID_BLOCK_POLY2D;

    const PDM_DMesh_nodal_section_poly2d_t *section =
    (const PDM_DMesh_nodal_section_poly2d_t *)
     PDM_Handles_get (mesh->sections_poly2d, _id_section);

    if (section == NULL) {
      PDM_error (__FILE__, __LINE__, 0, "Bad polygon section identifier\n");
    }

    return section->distrib;
  }

  else {

    _id_section = id_section - PDM_BLOCK_ID_BLOCK_STD;

    const PDM_DMesh_nodal_section_std_t *section =
    (const PDM_DMesh_nodal_section_std_t *)
     PDM_Handles_get (mesh->sections_std, _id_section);

    if (section == NULL) {
      PDM_error (__FILE__, __LINE__, 0, "Bad polyhedron section identifier\n");
    }

    return section->distrib;
  }
}


/**
 * \brief  Return cell distribution
 *
 * \param [in]  hdl  Distributed nodal mesh handle
 *
 * \return  A array of size \ref n_procs + 1
 *
 */

const PDM_g_num_t *
PDM_DMesh_nodal_distrib_cell_get
(
 const int  hdl
)
{
  PDM_DMesh_nodal_t * mesh =
          (PDM_DMesh_nodal_t *) PDM_Handles_get (mesh_handles, hdl);

  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }

  return mesh->cell_distrib;

}


/**
 * \brief  Return face distribution
 *
 * \param [in]  hdl  Distributed nodal mesh handle
 *
 * \return  A array of size \ref n_procs + 1
 *
 */

const PDM_g_num_t *
PDM_DMesh_nodal_distrib_face_get
(
 const int hdl
)
{
  PDM_DMesh_nodal_t * mesh =
          (PDM_DMesh_nodal_t *) PDM_Handles_get (mesh_handles, hdl);

  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }

  return mesh->face_distrib;

}


#ifdef __cplusplus
}
#endif /* __cplusplus */
