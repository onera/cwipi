/*----------------------------------------------------------------------------
 *  System headers
 *----------------------------------------------------------------------------*/
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_mpi.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_mesh_intersection_priv.h"
#include "pdm_mesh_intersection.h"
#include "pdm_timer.h"
#include "pdm_part_to_block.h"
#include "pdm_block_to_part.h"
#include "pdm_sort.h"
#include "pdm_array.h"
#include "pdm_logging.h"
#include "pdm_distrib.h"
#include "pdm_binary_search.h"
#include "pdm_part_mesh.h"
#include "pdm_part_mesh_priv.h"
#include "pdm_extract_part.h"
#include "pdm_extract_part_priv.h"
#include "pdm_part_connectivity_transform.h"
#include "pdm_vtk.h"
#include "pdm_box_priv.h"


#ifdef __cplusplus
extern "C"
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Local macro definitions
 *============================================================================*/

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Static function definitions
 *============================================================================*/

static
void
_export_vtk_2d
(
 const char               *pattern,
       PDM_extract_part_t *extrp_mesh
)
{
  int i_rank;
  PDM_MPI_Comm_rank(extrp_mesh->comm, &i_rank);

  for(int i_part = 0; i_part < extrp_mesh->n_part_out; ++i_part) {

    PDM_g_num_t *face_ln_to_gn = NULL;
    PDM_g_num_t *vtx_ln_to_gn  = NULL;
    int n_face = PDM_extract_part_ln_to_gn_get(extrp_mesh, i_part, PDM_MESH_ENTITY_FACE  , &face_ln_to_gn, PDM_OWNERSHIP_KEEP);
    int n_vtx  = PDM_extract_part_ln_to_gn_get(extrp_mesh, i_part, PDM_MESH_ENTITY_VERTEX, &vtx_ln_to_gn , PDM_OWNERSHIP_KEEP);

    double *vtx_coord = NULL;
    PDM_extract_part_vtx_coord_get(extrp_mesh, i_part, &vtx_coord, PDM_OWNERSHIP_KEEP);

    int  *face_edge     = NULL;
    int  *face_edge_idx = NULL;
    int  *edge_vtx      = NULL;
    int  *edge_vtx_idx  = NULL;
    PDM_extract_part_connectivity_get(extrp_mesh, i_part, PDM_CONNECTIVITY_TYPE_FACE_EDGE, &face_edge, &face_edge_idx, PDM_OWNERSHIP_KEEP);
    PDM_extract_part_connectivity_get(extrp_mesh, i_part, PDM_CONNECTIVITY_TYPE_EDGE_VTX , &edge_vtx , &edge_vtx_idx , PDM_OWNERSHIP_KEEP);

    int *face_vtx = NULL;
    PDM_compute_face_vtx_from_face_and_edge(n_face, face_edge_idx, face_edge, edge_vtx, &face_vtx);

    char filename[999];
    sprintf(filename, "%s_%i_%i.vtk", pattern, i_part, i_rank);
    PDM_vtk_write_polydata(filename,
                           n_vtx,
                           vtx_coord,
                           vtx_ln_to_gn,
                           n_face,
                           face_edge_idx,
                           face_vtx,
                           face_ln_to_gn,
                           NULL);


    free(face_vtx);
  }
}

static void
_export_vtk_3d
(
 const char               *name_chr,
       PDM_extract_part_t *extrp
)
{
  int i_rank;

  PDM_MPI_Comm_rank(extrp->comm, &i_rank);

  int          *pn_extract_cell        = (int          *) malloc(extrp->n_part_out * sizeof(int          ));
  int          *pn_extract_face        = (int          *) malloc(extrp->n_part_out * sizeof(int          ));
  int          *pn_extract_vtx         = (int          *) malloc(extrp->n_part_out * sizeof(int          ));
  int         **pextract_cell_face     = (int         **) malloc(extrp->n_part_out * sizeof(int         *));
  int         **pextract_cell_face_idx = (int         **) malloc(extrp->n_part_out * sizeof(int         *));
  int         **pextract_face_vtx      = (int         **) malloc(extrp->n_part_out * sizeof(int         *));
  int         **pextract_face_vtx_idx  = (int         **) malloc(extrp->n_part_out * sizeof(int         *));
  double      **pextract_vtx           = (double      **) malloc(extrp->n_part_out * sizeof(double      *));
  PDM_g_num_t **pextract_cell_ln_to_gn = (PDM_g_num_t **) malloc(extrp->n_part_out * sizeof(PDM_g_num_t *));
  PDM_g_num_t **pextract_face_ln_to_gn = (PDM_g_num_t **) malloc(extrp->n_part_out * sizeof(PDM_g_num_t *));
  PDM_g_num_t **pextract_vtx_ln_to_gn  = (PDM_g_num_t **) malloc(extrp->n_part_out * sizeof(PDM_g_num_t *));


  for(int i_part = 0; i_part < extrp->n_part_out; ++i_part) {

    pn_extract_cell[i_part] = PDM_extract_part_n_entity_get(extrp,
                                                            i_part,
                                                            PDM_MESH_ENTITY_CELL);

    pn_extract_face[i_part] = PDM_extract_part_n_entity_get(extrp,
                                                            i_part,
                                                            PDM_MESH_ENTITY_FACE);

    pn_extract_vtx[i_part] = PDM_extract_part_n_entity_get(extrp,
                                                           i_part,
                                                           PDM_MESH_ENTITY_VERTEX);

    PDM_extract_part_connectivity_get(extrp,
                                      i_part,
                                      PDM_CONNECTIVITY_TYPE_CELL_FACE,
                                      &pextract_cell_face[i_part],
                                      &pextract_cell_face_idx[i_part],
                                      PDM_OWNERSHIP_KEEP);

    PDM_extract_part_connectivity_get(extrp,
                                      i_part,
                                      PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                      &pextract_face_vtx[i_part],
                                      &pextract_face_vtx_idx[i_part],
                                      PDM_OWNERSHIP_KEEP);

    if(pextract_face_vtx[i_part] == NULL) {
      int *face_edge     = NULL;
      int *face_edge_idx = NULL;
      int *edge_vtx     = NULL;
      int *edge_vtx_idx = NULL;
      PDM_extract_part_connectivity_get(extrp,
                                        i_part,
                                        PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                        &face_edge,
                                        &face_edge_idx,
                                        PDM_OWNERSHIP_KEEP);
      PDM_extract_part_connectivity_get(extrp,
                                        i_part,
                                        PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                        &edge_vtx,
                                        &edge_vtx_idx,
                                        PDM_OWNERSHIP_KEEP);
      PDM_compute_face_vtx_from_face_and_edge(pn_extract_face[i_part], face_edge_idx, face_edge, edge_vtx, &pextract_face_vtx[i_part]);

      pextract_face_vtx_idx[i_part] = malloc((pn_extract_face[i_part]+1) * sizeof(int));
      for(int i = 0; i < pn_extract_face[i_part]+1; ++i) {
        pextract_face_vtx_idx[i_part][i] = face_edge_idx[i];
      }
    }


    PDM_extract_part_vtx_coord_get(extrp,
                                   i_part,
                                   &pextract_vtx[i_part],
                                   PDM_OWNERSHIP_KEEP);

    PDM_extract_part_ln_to_gn_get(extrp,
                                  i_part,
                                  PDM_MESH_ENTITY_CELL,
                                  &pextract_cell_ln_to_gn[i_part],
                                  PDM_OWNERSHIP_KEEP);

    PDM_extract_part_ln_to_gn_get(extrp,
                                  i_part,
                                  PDM_MESH_ENTITY_FACE,
                                  &pextract_face_ln_to_gn[i_part],
                                  PDM_OWNERSHIP_KEEP);

    PDM_extract_part_ln_to_gn_get(extrp,
                                  i_part,
                                  PDM_MESH_ENTITY_VERTEX,
                                  &pextract_vtx_ln_to_gn[i_part],
                                  PDM_OWNERSHIP_KEEP);

    // log_trace(" %s --> %i \n", name_chr, pn_extract_cell[i_part]);
    /* Vtk en légende */
    char filename[999];
    sprintf(filename, "%s_%3.3d_%3.3d.vtk", name_chr, i_part, i_rank);
    PDM_vtk_write_polydata(filename,
                           pn_extract_vtx[i_part],
                           pextract_vtx[i_part],
                           pextract_vtx_ln_to_gn[i_part],
                           pn_extract_face[i_part],
                           pextract_face_vtx_idx[i_part],
                           pextract_face_vtx[i_part],
                           pextract_face_ln_to_gn[i_part],
                           NULL);

  }

  // La visu concatene merge les valeurs donc on voit pas grand choses
  // _visu (name_chr,
  //        extrp->n_part_out,
  //        pn_extract_cell,
  //        pn_extract_face,
  //        pn_extract_vtx,
  //        pextract_cell_face_idx,
  //        pextract_cell_face,
  //        pextract_cell_ln_to_gn,
  //        pextract_face_vtx_idx,
  //        pextract_face_vtx,
  //        pextract_face_ln_to_gn,
  //        pextract_vtx,
  //        pextract_vtx_ln_to_gn);



  free(pn_extract_cell       );
  free(pn_extract_face       );
  free(pn_extract_vtx        );
  free(pextract_cell_face    );
  free(pextract_cell_face_idx);
  free(pextract_face_vtx     );
  free(pextract_face_vtx_idx );
  free(pextract_vtx          );
  free(pextract_cell_ln_to_gn);
  free(pextract_face_ln_to_gn);
  free(pextract_vtx_ln_to_gn );
}

static
void
_compute_extents_3d
(
 int      n_cell,
 int     *cell_face_idx,
 int     *cell_face,
 int     *face_vtx_idx,
 int     *face_vtx,
 double  *vtx_coord,
 double  *box_extents,
 double  *global_extents
)
{
  const double tolerance   = 1.e-12;
  const double eps_extents = 1.e-7;
  const int dim = 3;

  /* Loop over cell */
  for(int i_cell = 0; i_cell < n_cell; ++i_cell ) {

    double *_extents = box_extents + 6 * i_cell;

    /* Init */
    for (int i_dim = 0; i_dim < 3; i_dim++) {
      _extents[  i_dim] =  HUGE_VAL;
      _extents[3+i_dim] = -HUGE_VAL;
    }

    /* Loop over face and vtx to compute bbox */
    for(int idx_face = cell_face_idx[i_cell]; idx_face < cell_face_idx[i_cell+1]; ++idx_face) {

      int i_face = PDM_ABS(cell_face[idx_face])-1;

      for(int idx_vtx = face_vtx_idx[i_face]; idx_vtx < face_vtx_idx[i_face+1]; ++idx_vtx) {
        int i_vtx = face_vtx[idx_vtx]-1;

        for (int i_dim = 0; i_dim < 3; i_dim++) {
          double x = vtx_coord[3*i_vtx + i_dim];

          if (x < _extents[i_dim]) {
            _extents[i_dim] = x;
          }
          if (x > _extents[3+i_dim]) {
            _extents[3+i_dim] = x;
          }
        }
      }
    }

    double delta = 0.;
    for (int i_dim = 0; i_dim < 3; i_dim++) {
      double x = _extents[3+i_dim] - _extents[i_dim];

      if (delta < x) {
        delta = x;
      }
    }

    if (delta > eps_extents) {
      delta *= tolerance;
    } else {
      delta = eps_extents;
    }

    for (int i_dim = 0; i_dim < 3; i_dim++) {
      _extents[i_dim]   -= delta;
      _extents[3+i_dim] += delta;
    }

    for (int k = 0; k < dim; k++) {
      global_extents[k]       = PDM_MIN(_extents[k    ], global_extents[k    ]);
      global_extents[dim + k] = PDM_MAX(_extents[dim+k], global_extents[dim+k]);
    }

  } /* End loop cell */
}


static
void
_compute_extents_2d_from_face_vtx
(
 int      n_face,
 int     *face_vtx_idx,
 int     *face_vtx,
 double  *vtx_coord,
 double  *box_extents,
 double  *global_extents
)
{
  const double tolerance   = 1.e-12;
  const double eps_extents = 1.e-7;
  const int dim = 3;

  /* Loop over face */
  for(int i_face = 0; i_face < n_face; ++i_face ) {

    double *_extents = box_extents + 6 * i_face;

    /* Init */
    for (int i_dim = 0; i_dim < 3; i_dim++) {
      _extents[  i_dim] =  HUGE_VAL;
      _extents[3+i_dim] = -HUGE_VAL;
    }

    /* Loop over face and vtx to compute bbox */
    for(int idx_vtx = face_vtx_idx[i_face]; idx_vtx < face_vtx_idx[i_face+1]; ++idx_vtx) {
      int i_vtx = face_vtx[idx_vtx]-1;

      for (int i_dim = 0; i_dim < 3; i_dim++) {
        double x = vtx_coord[3*i_vtx + i_dim];

        if (x < _extents[i_dim]) {
          _extents[i_dim] = x;
        }
        if (x > _extents[3+i_dim]) {
          _extents[3+i_dim] = x;
        }
      }
    }

    double delta = 0.;
    for (int i_dim = 0; i_dim < 3; i_dim++) {
      double x = _extents[3+i_dim] - _extents[i_dim];

      if (delta < x) {
        delta = x;
      }
    }

    if (delta > eps_extents) {
      delta *= tolerance;
    } else {
      delta = eps_extents;
    }

    for (int i_dim = 0; i_dim < 3; i_dim++) {
      _extents[i_dim]   -= delta;
      _extents[3+i_dim] += delta;
    }

    for (int k = 0; k < dim; k++) {
      global_extents[k]       = PDM_MIN(_extents[k    ], global_extents[k    ]);
      global_extents[dim + k] = PDM_MAX(_extents[dim+k], global_extents[dim+k]);
    }

  } /* End loop cell */
}

static
void
_compute_part_mesh_extents
(
  PDM_part_mesh_t   *mesh,
  int                dim_mesh,
  double            *global_extents,
  double          ***extents_out
)
{
  int n_part = mesh->n_part;
  double **extents = malloc(n_part * sizeof(double *));
  if(dim_mesh == 3) {
    for(int i_part = 0; i_part < n_part; ++i_part) {

      int    *cell_face     = NULL;
      int    *cell_face_idx = NULL;
      int    *face_vtx_idx  = NULL;
      int    *face_vtx      = NULL;
      double *vtx_coord     = NULL;

      int n_cell = PDM_part_mesh_n_entity_get(mesh, i_part, PDM_MESH_ENTITY_CELL);

      PDM_part_mesh_connectivity_get(mesh, i_part, PDM_CONNECTIVITY_TYPE_CELL_FACE, &cell_face, &cell_face_idx, PDM_OWNERSHIP_USER);
      PDM_part_mesh_connectivity_get(mesh, i_part, PDM_CONNECTIVITY_TYPE_FACE_VTX , &face_vtx , &face_vtx_idx, PDM_OWNERSHIP_USER);

      PDM_part_mesh_vtx_coord_get(mesh, i_part, &vtx_coord, PDM_OWNERSHIP_USER); // Il faudrait un unchanged

      extents[i_part] = malloc(6 * n_cell * sizeof(double));
      if(face_vtx == NULL) {
        int    *face_edge_idx  = NULL;
        int    *face_edge      = NULL;
        int    *edge_vtx_idx   = NULL;
        int    *edge_vtx       = NULL;
        int n_face = PDM_part_mesh_n_entity_get(mesh, i_part, PDM_MESH_ENTITY_FACE);
        PDM_part_mesh_connectivity_get(mesh, i_part, PDM_CONNECTIVITY_TYPE_FACE_EDGE , &face_edge , &face_edge_idx, PDM_OWNERSHIP_USER);
        PDM_part_mesh_connectivity_get(mesh, i_part, PDM_CONNECTIVITY_TYPE_EDGE_VTX  , &edge_vtx  , &edge_vtx_idx , PDM_OWNERSHIP_USER);
        PDM_compute_face_vtx_from_face_and_edge(n_face, face_edge_idx, face_edge, edge_vtx, &face_vtx);
        _compute_extents_3d(n_cell, cell_face_idx, cell_face, face_edge_idx, face_vtx, vtx_coord, extents[i_part], global_extents);
        free(face_vtx);
      } else {
        _compute_extents_3d(n_cell, cell_face_idx, cell_face, face_vtx_idx, face_vtx, vtx_coord, extents[i_part], global_extents);
      }
    }
  } else if(dim_mesh == 2) {
    for(int i_part = 0; i_part < n_part; ++i_part) {

      int    *face_vtx      = NULL;
      int    *face_vtx_idx  = NULL;
      int    *face_edge_idx = NULL;
      int    *face_edge     = NULL;
      int    *edge_vtx_idx  = NULL;
      int    *edge_vtx      = NULL;
      double *vtx_coord     = NULL;

      // A gerer le cas mixte face_vtx ou face_edge + edge_vtx

      int n_face = PDM_part_mesh_n_entity_get(mesh, i_part, PDM_MESH_ENTITY_FACE);
      PDM_part_mesh_connectivity_get(mesh, i_part, PDM_CONNECTIVITY_TYPE_FACE_VTX , &face_vtx , &face_vtx_idx, PDM_OWNERSHIP_USER);

      PDM_part_mesh_vtx_coord_get(mesh, i_part, &vtx_coord, PDM_OWNERSHIP_USER); // Il faudrait un unchanged
      extents[i_part] = malloc(6 * n_face * sizeof(double));

      if(face_vtx != NULL) {
        _compute_extents_2d_from_face_vtx(n_face, face_vtx_idx, face_vtx, vtx_coord, extents[i_part], global_extents);
      } else {
        PDM_part_mesh_connectivity_get(mesh, i_part, PDM_CONNECTIVITY_TYPE_FACE_EDGE , &face_edge, &face_edge_idx, PDM_OWNERSHIP_USER);
        assert(face_edge != NULL);
        PDM_part_mesh_connectivity_get(mesh, i_part, PDM_CONNECTIVITY_TYPE_EDGE_VTX  , &edge_vtx, &edge_vtx_idx, PDM_OWNERSHIP_USER);
        assert(edge_vtx_idx == NULL);
        int n_edge = PDM_part_mesh_n_entity_get(mesh, i_part, PDM_MESH_ENTITY_EDGE);

        edge_vtx_idx = malloc((n_edge+1) * sizeof(int));
        for(int i_edge = 0; i_edge < n_edge+1; ++i_edge){
          edge_vtx_idx[i_edge] = 2 * i_edge;
        }
        _compute_extents_3d(n_face, face_edge_idx, face_edge, edge_vtx_idx, edge_vtx, vtx_coord, extents[i_part], global_extents);
        // _compute_extents_2d_from_face_edge(n_face, face_edge_idx, face_edge, edge_vtx, vtx_coord, extents[i_part], global_extents);
        free(edge_vtx_idx);
      }
    }

  } else {
    int    *edge_vtx_idx  = NULL;
    int    *edge_vtx      = NULL;
    double *vtx_coord     = NULL;
    for(int i_part = 0; i_part < n_part; ++i_part) {
      int n_edge = PDM_part_mesh_n_entity_get(mesh, i_part, PDM_MESH_ENTITY_EDGE);

      extents[i_part] = malloc(6 * n_edge * sizeof(double));
      PDM_part_mesh_vtx_coord_get(mesh, i_part, &vtx_coord, PDM_OWNERSHIP_USER); // Il faudrait un unchanged

      PDM_part_mesh_connectivity_get(mesh, i_part, PDM_CONNECTIVITY_TYPE_EDGE_VTX  , &edge_vtx, &edge_vtx_idx, PDM_OWNERSHIP_USER);
      assert(edge_vtx_idx == NULL);
      edge_vtx_idx = malloc((n_edge+1) * sizeof(int));
      for(int i_edge = 0; i_edge < n_edge+1; ++i_edge){
        edge_vtx_idx[i_edge] = 2 * i_edge;
      }

      _compute_extents_2d_from_face_vtx(n_edge, edge_vtx_idx, edge_vtx, vtx_coord, extents[i_part], global_extents);

      free(edge_vtx_idx);
    }
  }
  *extents_out = extents;
}

static
void
_select_elements_by_global_bbox
(
  PDM_part_mesh_t *mesh,
  int              dim_mesh,
  double         **box_extents,
  double          *g_mesh_global_extents,
  int            **n_extract_elmt_out,
  double        ***extract_box_extents_out,
  int           ***extract_elmt_init_location_out,
  PDM_g_num_t   ***extract_elmt_ln_to_gn_out
)
{
  int n_part = mesh->n_part;
  int i_rank;
  PDM_MPI_Comm_rank(mesh->comm, &i_rank);

  int          *n_extract_elmt             = malloc(n_part * sizeof(int         *));
  double      **extract_box_extents        = malloc(n_part * sizeof(double      *));
  int         **extract_elmt_init_location = malloc(n_part * sizeof(int         *));
  PDM_g_num_t **extract_elmt_ln_to_gn      = malloc(n_part * sizeof(PDM_g_num_t *));

  for(int i_part = 0; i_part < n_part; ++i_part) {

    int n_entity = 0;
    PDM_g_num_t* entity_ln_to_gn = NULL;
    if(dim_mesh == 3) {
      n_entity = PDM_part_mesh_n_entity_get(mesh, i_part, PDM_MESH_ENTITY_CELL);
      PDM_part_mesh_entity_ln_to_gn_get(mesh, i_part, PDM_MESH_ENTITY_CELL, &entity_ln_to_gn, PDM_OWNERSHIP_USER);
    } else if(dim_mesh == 2) {
      n_entity = PDM_part_mesh_n_entity_get(mesh, i_part, PDM_MESH_ENTITY_FACE);
      PDM_part_mesh_entity_ln_to_gn_get(mesh, i_part, PDM_MESH_ENTITY_FACE, &entity_ln_to_gn, PDM_OWNERSHIP_USER);
    } else if(dim_mesh == 1) {
      n_entity = PDM_part_mesh_n_entity_get(mesh, i_part, PDM_MESH_ENTITY_EDGE);
      PDM_part_mesh_entity_ln_to_gn_get(mesh, i_part, PDM_MESH_ENTITY_EDGE, &entity_ln_to_gn, PDM_OWNERSHIP_USER);
    }

    // char filename[999];
    // sprintf(filename, "titi_%i.vtk", n_entity);
    // PDM_vtk_write_boxes(filename,
    //                     n_entity,
    //                     box_extents[i_part],
    //                     NULL);

    n_extract_elmt[i_part] = 0;
    extract_box_extents       [i_part] = malloc(6 * n_entity * sizeof(double     ));
    extract_elmt_init_location[i_part] = malloc(3 * n_entity * sizeof(int        ));
    extract_elmt_ln_to_gn     [i_part] = malloc(    n_entity * sizeof(PDM_g_num_t));

    double *_box_extents = box_extents[i_part];

    for(int i = 0; i < n_entity; ++i) {

      double *box_min = _box_extents + 6*i;
      double *box_max = box_min + 3;

      int intersect = 1;
      for (int j = 0; j < 3; j++) {
        if (box_min[j] > g_mesh_global_extents[j+3] ||
            box_max[j] < g_mesh_global_extents[j  ]) {
          intersect = 0;
          break;
        }
      }

      if (intersect) {
        for (int j = 0; j < 6; j++) {
          extract_box_extents  [i_part][6*n_extract_elmt[i_part]+j] = _box_extents[6*i+j];
        }
        extract_elmt_init_location[i_part][3*n_extract_elmt[i_part]  ] = i_rank;
        extract_elmt_init_location[i_part][3*n_extract_elmt[i_part]+1] = i_part;
        extract_elmt_init_location[i_part][3*n_extract_elmt[i_part]+2] = i;

        extract_elmt_ln_to_gn[i_part][n_extract_elmt[i_part]] = entity_ln_to_gn[i];

        n_extract_elmt[i_part]++;
      }
    }
    extract_box_extents       [i_part] = realloc(extract_box_extents       [i_part], 6 * n_entity * sizeof(double     ));
    extract_elmt_init_location[i_part] = realloc(extract_elmt_init_location[i_part], 3 * n_entity * sizeof(int        ));
    extract_elmt_ln_to_gn     [i_part] = realloc(extract_elmt_ln_to_gn     [i_part],     n_entity * sizeof(PDM_g_num_t));

  }

  *n_extract_elmt_out             = n_extract_elmt;
  *extract_box_extents_out        = extract_box_extents;
  *extract_elmt_init_location_out = extract_elmt_init_location;
  *extract_elmt_ln_to_gn_out      = extract_elmt_ln_to_gn;
}


static
void
_redistrib_boxes
(
 PDM_MPI_Comm    comm,
 PDM_box_set_t  *boxes_mesh_a,
 PDM_box_set_t  *boxes_mesh_b,
 int            *box_a_to_box_b_idx,
 int            *box_a_to_box_b,
 int           **redistribute_box_a_to_box_b_idx,
 int           **redistribute_box_a_to_box_b
)
{

  int              n_elt_mesh_a    = PDM_box_set_get_size (boxes_mesh_a);
  int              n_elt_mesh_b    = PDM_box_set_get_size (boxes_mesh_b);

  PDM_g_num_t *gnum_elt_mesh_a = (PDM_g_num_t *) PDM_box_set_get_g_num (boxes_mesh_a);
  PDM_g_num_t *gnum_elt_mesh_b = (PDM_g_num_t *) PDM_box_set_get_g_num (boxes_mesh_b);

  /*****************************************************************************
   *                                                                           *
   *  Transfer intersection information from partitions to blocks              *
   * with PDM_part_to_block_exch function                                      *
   *                                                                           *
   *  Results :                                                                *
   *      - block_a_boxes_b_idx                                                *
   *      - block_a_boxes_b_gnum_data                                          *
   *                                                                           *
   ****************************************************************************/

  /*
   * Tentative Bruno :
   *   - Ponderate work
   *   - Extract only cell with job
   */
  double* weight = (double *) malloc( n_elt_mesh_a * sizeof(double));
  // PDM_g_num_t* extract_mesh_a_g_num = (PDM_g_num_t *) malloc( n_elt_mesh_a * sizeof(PDM_g_num_t));
  for (int i = 0; i < n_elt_mesh_a; i++) {
    weight[i] = box_a_to_box_b_idx[i+1] - box_a_to_box_b_idx[i];
  }

  // TODO : Geometric to better locality
  PDM_part_to_block_t *ptb_boxes_a = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                              PDM_PART_TO_BLOCK_POST_MERGE,
                                                              1.,
                                            (PDM_g_num_t **) &gnum_elt_mesh_a,
                                                             &weight,
                                                              &n_elt_mesh_a,
                                                              1,
                                                              comm);

  int n_elt_block_a = PDM_part_to_block_n_elt_block_get (ptb_boxes_a);
  free(weight);

  PDM_g_num_t *block_gnum_a = PDM_part_to_block_block_gnum_get (ptb_boxes_a);

  int *part_stride_a = (int *) malloc (sizeof(int) * n_elt_mesh_a);

  for (int i = 0; i < n_elt_mesh_a; i++) {
    part_stride_a[i] = box_a_to_box_b_idx[i+1] - box_a_to_box_b_idx[i];
  }

  /*
   * Exchange connectivity box_a_to_box_b
   */
  PDM_g_num_t *box_a_to_box_b_g_num = (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) *
                                                               box_a_to_box_b_idx[n_elt_mesh_a]);

  for (int k = 0; k < box_a_to_box_b_idx[n_elt_mesh_a]; k++) {
    box_a_to_box_b_g_num[k] =  gnum_elt_mesh_b[box_a_to_box_b[k]];
  }

  int         *block_a_boxes_b_stride;
  PDM_g_num_t *block_a_boxes_b_gnum_data;

  PDM_part_to_block_exch (ptb_boxes_a,
                         sizeof(PDM_g_num_t),
                         PDM_STRIDE_VAR_INTERLACED,
                         0,
                         &part_stride_a,
               (void **) &box_a_to_box_b_g_num,
                         &block_a_boxes_b_stride,
               (void **) &block_a_boxes_b_gnum_data);
  free(box_a_to_box_b_g_num);
  /*****************************************************************************
   *                                                                           *
   * Redistribute boxes_mesh_a intersections to ensure a good load balacing          *
   * in comm MPI communicator                                              *
   * This step removes intersections found many times on different ranks       *
   *                                                                           *
   * After this step, data are stored in a block with n_elt_block_a, block_gnum_a,
   * part_stride_a                                                              *
   *                                                                           *
   ****************************************************************************/
  int n_rank;
  PDM_MPI_Comm_size (comm, &n_rank);

  /*
   * Redistribute data boxes_mesh_a from blockB distribution with a PDM_box_distrib_t
   * structure
   * TODO: - Build a new PDM_box_distrib_t structure more simple
   *       - Hide PDM_box_distrib_t attributes
   */

  PDM_l_num_t *destination = PDM_part_to_block_destination_get (ptb_boxes_a);

  PDM_g_num_t n_g_elmt_mesh_a = PDM_box_set_get_global_size(boxes_mesh_a);

  PDM_box_distrib_t *distrib_a = PDM_box_distrib_create(n_elt_mesh_a,
                                                        n_g_elmt_mesh_a,
                                                        1, // Don't use in this case
                                                        comm);

  PDM_g_num_t n_g_elmt_mesh_b = PDM_box_set_get_global_size(boxes_mesh_b);

  PDM_box_distrib_t *distrib_b = PDM_box_distrib_create(n_elt_mesh_b,
                                                        n_g_elmt_mesh_b,
                                                        1, // Don't use in this case
                                                        comm);


  int *count_elts_a = (int *) malloc (sizeof(int) * n_rank);
  int *count_elts_b = (int *) malloc (sizeof(int) * n_rank);

  for (int i = 0; i < n_rank + 1; i++) {
    distrib_a->index[i] = 0;
    distrib_b->index[i] = 0;
  }

  for (int i = 0; i < n_rank; i++) {
    count_elts_a[i] = 0;
    count_elts_b[i] = 0;
  }

  for (int i = 0; i < n_elt_mesh_a; i++) {
    int t_rank = destination[i] + 1;

    if(part_stride_a[i] > 0 ) { // To see with Eric and Bastien --> I use it to extract only the intersect part
      distrib_a->index[t_rank]++;
      distrib_b->index[t_rank] += part_stride_a[i];
    }
  }

  for (int i = 0; i < n_rank; i++) {
    distrib_a->index[i+1] += distrib_a->index[i];
    distrib_b->index[i+1] += distrib_b->index[i];
  }

  distrib_a->list = (int *) malloc (sizeof(int) * distrib_a->index[n_rank]);
  distrib_b->list = (int *) malloc (sizeof(int) * distrib_b->index[n_rank]);

  for (int i = 0; i < n_elt_mesh_a; i++) {
    if(part_stride_a[i] > 0 ) { // To see with Eric and Bastien --> I use it to extract only the intersect part
      int t_rank = destination[i]; // EQU + 1; mais ce n est pas necessaire
      int idx_a = distrib_a->index[t_rank] + (count_elts_a[t_rank]++);
      distrib_a->list[idx_a] = i;
      int idx_b = distrib_b->index[t_rank] + count_elts_b[t_rank];
      count_elts_b[t_rank] += part_stride_a[i];
      int k=0;
      for (int j = box_a_to_box_b_idx[i]; j < box_a_to_box_b_idx[i+1]; j++) {
        distrib_b->list[idx_b+k++] = box_a_to_box_b[j];
      }
    }
  }


  free (part_stride_a);
  free (count_elts_a);
  free (count_elts_b);

  PDM_box_distrib_clean (distrib_a);
  PDM_box_distrib_clean (distrib_b);

  PDM_box_set_redistribute (distrib_a, boxes_mesh_a);
  PDM_box_set_redistribute (distrib_b, boxes_mesh_b);

  PDM_box_distrib_destroy (&distrib_a);
  PDM_box_distrib_destroy (&distrib_b);

  PDM_box_set_remove_duplicate (boxes_mesh_a);
  PDM_box_set_remove_duplicate (boxes_mesh_b);

  /*
   * All boxes are redistribute we need to update box_a_to_box_b array
   *    - Caution if morton / hilbert the array block_gnum_a IS NOT order
   */
  n_elt_mesh_a    = PDM_box_set_get_size (boxes_mesh_a); // Caution not the same of the fist call because redistibute
  n_elt_mesh_b    = PDM_box_set_get_size (boxes_mesh_b); // Caution not the same of the fist call because redistibute

  gnum_elt_mesh_a = (PDM_g_num_t *) PDM_box_set_get_g_num(boxes_mesh_a);
  gnum_elt_mesh_b = (PDM_g_num_t *) PDM_box_set_get_g_num(boxes_mesh_b);

  PDM_block_to_part_t *btp = PDM_block_to_part_create_from_sparse_block(block_gnum_a,
                                                                        n_elt_block_a,
                                                (const PDM_g_num_t **)  &gnum_elt_mesh_a,
                                                                        &n_elt_mesh_a,
                                                                        1,
                                                                        comm);

  PDM_part_to_block_free (ptb_boxes_a);
  int         **tmp_redistribute_box_a_to_box_b_n     = NULL;
  PDM_g_num_t **tmp_redistribute_box_a_to_box_b_g_num = NULL;
  PDM_block_to_part_exch(btp,
                         sizeof(PDM_g_num_t),
                         PDM_STRIDE_VAR_INTERLACED,
                         block_a_boxes_b_stride,
                         block_a_boxes_b_gnum_data,
                         &tmp_redistribute_box_a_to_box_b_n,
           (void ***)    &tmp_redistribute_box_a_to_box_b_g_num);
  free (block_a_boxes_b_stride);
  free (block_a_boxes_b_gnum_data);

  PDM_block_to_part_free(btp);

  int         *redistribute_box_a_to_box_b_n     = tmp_redistribute_box_a_to_box_b_n    [0];
  PDM_g_num_t *redistribute_box_a_to_box_b_g_num = tmp_redistribute_box_a_to_box_b_g_num[0];
  free(tmp_redistribute_box_a_to_box_b_n    );
  free(tmp_redistribute_box_a_to_box_b_g_num);


  /*
   * Translate in frame of B
   */

  int         *order              = (int         *) malloc(n_elt_mesh_b * sizeof(int        ));
  PDM_g_num_t *gnum_elt_mesh_b_cp = (PDM_g_num_t *) malloc(n_elt_mesh_b * sizeof(PDM_g_num_t));
  for(int i = 0; i < n_elt_mesh_b; ++i ) {
    order             [i] = i;
    gnum_elt_mesh_b_cp[i] = gnum_elt_mesh_b[i];
  }

  PDM_sort_long(gnum_elt_mesh_b_cp, order, n_elt_mesh_b);


  int *_redistribute_box_a_to_box_b_idx = (int *) malloc((n_elt_mesh_a+1) * sizeof(int));
  _redistribute_box_a_to_box_b_idx[0] = 0;
  int n_tot_connect = 0;
  for(int i = 0; i < n_elt_mesh_a; ++i) {
    n_tot_connect += redistribute_box_a_to_box_b_n[i];
    _redistribute_box_a_to_box_b_idx[i+1] = _redistribute_box_a_to_box_b_idx[i] + redistribute_box_a_to_box_b_n[i];
  }

  int *_redistribute_box_a_to_box_b = (int *) malloc( n_tot_connect * sizeof(int));

  for(int i = 0; i < n_tot_connect; ++i) {
    int pos = PDM_binary_search_long(redistribute_box_a_to_box_b_g_num[i], gnum_elt_mesh_b_cp, n_elt_mesh_b);
    _redistribute_box_a_to_box_b[i] = order[pos];
  }

  *redistribute_box_a_to_box_b_idx = _redistribute_box_a_to_box_b_idx;
  *redistribute_box_a_to_box_b     = _redistribute_box_a_to_box_b;

  free(order);
  free(gnum_elt_mesh_b_cp);
  free(redistribute_box_a_to_box_b_n    );
  free(redistribute_box_a_to_box_b_g_num);
}


static
PDM_extract_part_t*
_create_extract_part
(
 PDM_part_mesh_t *mesh,
 int              dim_mesh,
 PDM_box_set_t   *boxes_meshes
)
{
  int n_part_out = 1;
  PDM_extract_part_t* extrp_mesh = PDM_extract_part_create(dim_mesh,
                                                           mesh->n_part,
                                                           n_part_out,
                                                           PDM_EXTRACT_PART_KIND_FROM_TARGET,
                                                           PDM_SPLIT_DUAL_WITH_HILBERT, // Not used
                                                           PDM_FALSE,                   // compute_child_gnum
                                                           PDM_OWNERSHIP_KEEP,
                                                           mesh->comm);

  int              n_elt_mesh    = PDM_box_set_get_size (boxes_meshes);

  printf("n_elt_mesh = %i  \n", n_elt_mesh);

  PDM_g_num_t *gnum_elt_mesh = (PDM_g_num_t *) PDM_box_set_get_g_num (boxes_meshes);

  int *init_location_elt_mesh = (int  *) PDM_box_set_origin_get(boxes_meshes);


  for(int i_part = 0; i_part < mesh->n_part; ++i_part) {

    int n_cell = PDM_part_mesh_n_entity_get(mesh, i_part, PDM_MESH_ENTITY_CELL  );
    int n_face = PDM_part_mesh_n_entity_get(mesh, i_part, PDM_MESH_ENTITY_FACE  );
    int n_edge = PDM_part_mesh_n_entity_get(mesh, i_part, PDM_MESH_ENTITY_EDGE  );
    int n_vtx  = PDM_part_mesh_n_entity_get(mesh, i_part, PDM_MESH_ENTITY_VERTEX);

    PDM_g_num_t *cell_ln_to_gn = NULL;
    PDM_g_num_t *face_ln_to_gn = NULL;
    PDM_g_num_t *edge_ln_to_gn = NULL;
    PDM_g_num_t *vtx_ln_to_gn  = NULL;

    PDM_part_mesh_entity_ln_to_gn_get(mesh, i_part, PDM_MESH_ENTITY_CELL  , &cell_ln_to_gn, PDM_OWNERSHIP_USER);
    PDM_part_mesh_entity_ln_to_gn_get(mesh, i_part, PDM_MESH_ENTITY_FACE  , &face_ln_to_gn, PDM_OWNERSHIP_USER);
    PDM_part_mesh_entity_ln_to_gn_get(mesh, i_part, PDM_MESH_ENTITY_EDGE  , &edge_ln_to_gn, PDM_OWNERSHIP_USER);
    PDM_part_mesh_entity_ln_to_gn_get(mesh, i_part, PDM_MESH_ENTITY_VERTEX, &vtx_ln_to_gn , PDM_OWNERSHIP_USER);

    int *cell_face     = NULL;
    int *cell_face_idx = NULL;
    int *face_vtx      = NULL;
    int *face_vtx_idx  = NULL;
    int *face_edge     = NULL;
    int *face_edge_idx = NULL;
    int *edge_vtx      = NULL;
    int *edge_vtx_idx  = NULL;

    PDM_part_mesh_connectivity_get(mesh, i_part, PDM_CONNECTIVITY_TYPE_CELL_FACE, &cell_face, &cell_face_idx, PDM_OWNERSHIP_USER);
    PDM_part_mesh_connectivity_get(mesh, i_part, PDM_CONNECTIVITY_TYPE_FACE_VTX , &face_vtx , &face_vtx_idx , PDM_OWNERSHIP_USER);
    PDM_part_mesh_connectivity_get(mesh, i_part, PDM_CONNECTIVITY_TYPE_FACE_EDGE, &face_edge, &face_edge_idx, PDM_OWNERSHIP_USER);
    PDM_part_mesh_connectivity_get(mesh, i_part, PDM_CONNECTIVITY_TYPE_EDGE_VTX , &edge_vtx , &edge_vtx_idx , PDM_OWNERSHIP_USER);

    double *vtx_coord = NULL;
    PDM_part_mesh_vtx_coord_get(mesh, i_part, &vtx_coord, PDM_OWNERSHIP_USER);

    PDM_extract_part_part_set(extrp_mesh,
                              i_part,
                              n_cell,
                              n_face,
                              n_edge,
                              n_vtx ,
                              cell_face_idx,
                              cell_face,
                              face_edge_idx,
                              face_edge,
                              edge_vtx,
                              face_vtx_idx,
                              face_vtx,
                              cell_ln_to_gn,
                              face_ln_to_gn,
                              edge_ln_to_gn,
                              vtx_ln_to_gn,
                              vtx_coord);
  }




  /*  Setup target frame */
  PDM_extract_part_target_set(extrp_mesh, 0, n_elt_mesh, gnum_elt_mesh, init_location_elt_mesh);

  PDM_extract_part_compute(extrp_mesh);

  return extrp_mesh;
}

static
void
_mesh_intersection_vol_vol
(
 PDM_mesh_intersection_t *mi,
 PDM_extract_part_t      *extrp_mesh_a,
 PDM_extract_part_t      *extrp_mesh_b,
 int                     *redistribute_box_a_to_box_b_idx,
 int                     *redistribute_box_a_to_box_b
)
{
  PDM_UNUSED(mi);
  PDM_UNUSED(extrp_mesh_a);
  PDM_UNUSED(extrp_mesh_b);
  PDM_UNUSED(redistribute_box_a_to_box_b_idx);
  PDM_UNUSED(redistribute_box_a_to_box_b);

}

static
void
_mesh_intersection_vol_surf
(
 PDM_mesh_intersection_t *mi,
 PDM_extract_part_t      *extrp_mesh_a,
 PDM_extract_part_t      *extrp_mesh_b,
 int                     *redistribute_box_a_to_box_b_idx,
 int                     *redistribute_box_a_to_box_b
)
{
  PDM_UNUSED(mi);
  PDM_UNUSED(extrp_mesh_a);
  PDM_UNUSED(extrp_mesh_b);
  PDM_UNUSED(redistribute_box_a_to_box_b_idx);
  PDM_UNUSED(redistribute_box_a_to_box_b);

  /*
   * Panic vtk
   */
  if(0 == 1) {
    _export_vtk_3d("extrp_mesh_a", extrp_mesh_a);
    _export_vtk_2d("extrp_mesh_b", extrp_mesh_b);
  }

}


static
void
_mesh_intersection_surf_surf
(
 PDM_mesh_intersection_t *mi,
 PDM_extract_part_t      *extrp_mesh_a,
 PDM_extract_part_t      *extrp_mesh_b,
 int                     *redistribute_box_a_to_box_b_idx,
 int                     *redistribute_box_a_to_box_b
)
{
  PDM_UNUSED(mi);
  PDM_UNUSED(extrp_mesh_a);
  PDM_UNUSED(extrp_mesh_b);
  PDM_UNUSED(redistribute_box_a_to_box_b_idx);
  PDM_UNUSED(redistribute_box_a_to_box_b);

  /*
   * Panic vtk
   */
  if(0 == 1) {
    _export_vtk_2d("extrp_mesh_a", extrp_mesh_a);
    _export_vtk_2d("extrp_mesh_b", extrp_mesh_b);
  }

  /*
   * TODO : Faire la géometrie :D
   */


}


static
void
_mesh_intersection_surf_line
(
 PDM_mesh_intersection_t *mi,
 PDM_extract_part_t      *extrp_mesh_a,
 PDM_extract_part_t      *extrp_mesh_b,
 int                     *redistribute_box_a_to_box_b_idx,
 int                     *redistribute_box_a_to_box_b
)
{
  PDM_UNUSED(mi);
  PDM_UNUSED(extrp_mesh_a);
  PDM_UNUSED(extrp_mesh_b);
  PDM_UNUSED(redistribute_box_a_to_box_b_idx);
  PDM_UNUSED(redistribute_box_a_to_box_b);

}

/*=============================================================================
 * Public function definitions
 *============================================================================*/

PDM_mesh_intersection_t*
PDM_mesh_intersection_create
(
 const PDM_mesh_intersection_kind_t intersection_kind,
 const int                          dim_mesh_a,
 const int                          dim_mesh_b,
 const int                          n_part_mesh_a,
 const int                          n_part_mesh_b,
 const double                       project_coeff,
       PDM_MPI_Comm                 comm
)
{
  PDM_mesh_intersection_t *mi = (PDM_mesh_intersection_t *) malloc(sizeof(PDM_mesh_intersection_t));

  mi->comm = comm;
  mi->intersect_kind = intersection_kind;
  mi->n_part_mesh_a  = n_part_mesh_a;
  mi->n_part_mesh_b  = n_part_mesh_b;
  mi->dim_mesh_a     = dim_mesh_a;
  mi->dim_mesh_b     = dim_mesh_b;
  mi->project_coef   = project_coeff;

  mi->mesh_a = PDM_part_mesh_create(n_part_mesh_a, comm);
  mi->mesh_b = PDM_part_mesh_create(n_part_mesh_b, comm);

  return mi;
}

void
PDM_mesh_intersection_compute
(
  PDM_mesh_intersection_t  *mi
)
{
  /*
   * Compute extents of mesh_a and mesh_b
   */
  double **extents_mesh_a = NULL;
  double **extents_mesh_b = NULL;
  double global_extents[6] = { -HUGE_VAL, -HUGE_VAL, -HUGE_VAL, HUGE_VAL,  HUGE_VAL,  HUGE_VAL};
  double mesh_global_extents[2][6] = {{ HUGE_VAL,  HUGE_VAL,  HUGE_VAL,
                                       -HUGE_VAL, -HUGE_VAL, -HUGE_VAL},
                                      {HUGE_VAL,  HUGE_VAL,  HUGE_VAL,
                                       -HUGE_VAL, -HUGE_VAL, -HUGE_VAL}};
  double g_mesh_global_extents[2][6];
  double g_global_extents        [6];
  _compute_part_mesh_extents(mi->mesh_a, mi->dim_mesh_a, mesh_global_extents[0], &extents_mesh_a);
  _compute_part_mesh_extents(mi->mesh_b, mi->dim_mesh_b, mesh_global_extents[1], &extents_mesh_b);

  /*
   * Global extents exchange
   */
  const int dim = 3;
  for(int i_mesh = 0; i_mesh < 2; ++i_mesh) {
    PDM_MPI_Allreduce(mesh_global_extents[i_mesh], g_mesh_global_extents[i_mesh], dim,
                      PDM_MPI_DOUBLE, PDM_MPI_MIN, mi->comm);
    PDM_MPI_Allreduce(mesh_global_extents[i_mesh]+dim, g_mesh_global_extents[i_mesh]+dim, dim,
                      PDM_MPI_DOUBLE, PDM_MPI_MAX, mi->comm);
  }

  /* Union or intersection of global extents */
  for(int i_mesh = 0; i_mesh < 2; ++i_mesh) {
    for (int k = 0; k < 3; k++) {
      // Union
      // global_extents[k]     = PDM_MIN(mesh_global_extents[i_mesh][k  ], global_extents[k  ]);
      // global_extents[3 + k] = PDM_MAX(mesh_global_extents[i_mesh][3+k], global_extents[3+k]);
      // Intersection
      global_extents[k]     = PDM_MAX(g_mesh_global_extents[i_mesh][k  ], global_extents[k  ]);
      global_extents[3 + k] = PDM_MIN(g_mesh_global_extents[i_mesh][3+k], global_extents[3+k]);
    }
  }
  for(int i = 0; i < 6; ++i) {
    g_global_extents[i] = global_extents[i];
  }
  double max_range = -HUGE_VAL;
  double min_range =  HUGE_VAL;

  for (int k = 0; k < dim; k++) {
    max_range = PDM_MAX(max_range, (g_global_extents[dim+k] - g_global_extents[k]));
    min_range = PDM_MIN(min_range, (g_global_extents[dim+k] - g_global_extents[k]));
  }

  for (int k = 0; k < dim; k++) {
    g_global_extents[k]     += -max_range * 1.1e-3; // On casse la symetrie !
    g_global_extents[dim+k] +=  max_range * 1e-3;
  }


  /*
   * Extraction - In option ?
   */
  int           n_mesh = 2;
  int           n_part                    [n_mesh];
  int          *n_extract_elmt            [n_mesh];
  double      **extract_box_extents       [n_mesh];
  int         **extract_elmt_init_location[n_mesh];
  PDM_g_num_t **extract_elmt_ln_to_gn     [n_mesh];

  n_part[0] = mi->n_part_mesh_a;
  n_part[1] = mi->n_part_mesh_b;

  _select_elements_by_global_bbox(mi->mesh_a, mi->dim_mesh_a,
                                  extents_mesh_a, g_mesh_global_extents[1], // On enleve tout ce qui est en dehors de B
                                  &n_extract_elmt[0],
                                  &extract_box_extents[0],
                                  &extract_elmt_init_location[0],
                                  &extract_elmt_ln_to_gn[0]);
  _select_elements_by_global_bbox(mi->mesh_b, mi->dim_mesh_b,
                                  extents_mesh_b, g_mesh_global_extents[0], // On enleve tout ce qui est en dehors de A
                                  &n_extract_elmt[1],
                                  &extract_box_extents[1],
                                  &extract_elmt_init_location[1],
                                  &extract_elmt_ln_to_gn[1]);

  for(int i_part = 0; i_part < mi->n_part_mesh_a; ++i_part) {
    free(extents_mesh_a[i_part]);
  }
  for(int i_part = 0; i_part < mi->n_part_mesh_b; ++i_part) {
    free(extents_mesh_b[i_part]);
  }
  free(extents_mesh_a);
  free(extents_mesh_b);

  // Attention le dbtree fait  le init_location sauf que la il faut le forcer !!!!!!
  PDM_dbbtree_t *dbbtree_mesh_a = PDM_dbbtree_create (mi->comm, dim, g_global_extents);

  PDM_box_set_t  *boxes_mesh_a = PDM_dbbtree_boxes_set_with_init_location(dbbtree_mesh_a,
                                                                          mi->mesh_a->n_part,
                                                                          n_extract_elmt            [0],
                                                  (const int         **)  extract_elmt_init_location[0],
                                                  (const double      **)  extract_box_extents       [0],
                                                  (const PDM_g_num_t **)  extract_elmt_ln_to_gn     [0]);

  /*
   * Intersect with B
   */
  int *box_a_to_box_b_idx = NULL;
  int *box_a_to_box_b     = NULL;
  PDM_box_set_t  *boxes_mesh_b =  PDM_dbbtree_intersect_boxes_with_init_location_set(dbbtree_mesh_a,
                                                                                     mi->mesh_b->n_part,
                                                                                     n_extract_elmt            [1],
                                                              (const int         **) extract_elmt_init_location[1],
                                                              (const double      **) extract_box_extents       [1],
                                                              (const PDM_g_num_t **) extract_elmt_ln_to_gn     [1],
                                                                                     &box_a_to_box_b_idx,
                                                                                     &box_a_to_box_b);

  /* Free extraction */
  for(int i_mesh = 0; i_mesh < n_mesh; ++i_mesh) {
    for(int i_part = 0; i_part < n_part[i_mesh]; ++i_part) {
      free(extract_elmt_init_location[i_mesh][i_part]);
      free(extract_box_extents       [i_mesh][i_part]);
      free(extract_elmt_ln_to_gn     [i_mesh][i_part]);
    }
    free(n_extract_elmt            [i_mesh]);
    free(extract_elmt_init_location[i_mesh]);
    free(extract_box_extents       [i_mesh]);
    free(extract_elmt_ln_to_gn     [i_mesh]);
  }


  /*
   *  Redistrib all boxes (inplace) like overlay before extracting mesh
   */
  int *redistribute_box_a_to_box_b_idx = NULL;
  int *redistribute_box_a_to_box_b     = NULL;
  _redistrib_boxes(mi->comm,
                   boxes_mesh_a,
                   boxes_mesh_b,
                   box_a_to_box_b_idx,
                   box_a_to_box_b,
                   &redistribute_box_a_to_box_b_idx,
                   &redistribute_box_a_to_box_b);
  free(box_a_to_box_b_idx);
  free(box_a_to_box_b);


  /*
   * Extract part
   */
  PDM_extract_part_t* extrp_mesh_a = _create_extract_part(mi->mesh_a,
                                                          mi->dim_mesh_a,
                                                          boxes_mesh_a);
  PDM_extract_part_t* extrp_mesh_b = _create_extract_part(mi->mesh_b,
                                                          mi->dim_mesh_b,
                                                          boxes_mesh_b);

  PDM_dbbtree_free (dbbtree_mesh_a);
  PDM_box_set_destroy (&boxes_mesh_a);
  PDM_box_set_destroy (&boxes_mesh_b);

  /*
   * Geometry begin here ...
   */
  if(mi->dim_mesh_a == 3 && mi->dim_mesh_b == 3) {
    _mesh_intersection_vol_vol(mi,
                               extrp_mesh_a,
                               extrp_mesh_b,
                               redistribute_box_a_to_box_b_idx,
                               redistribute_box_a_to_box_b);
  } else if(mi->dim_mesh_a == 3 && mi->dim_mesh_b == 2) {
    // On suppose que l'utilisateur met A = Vol et B = Surf
    _mesh_intersection_vol_surf(mi,
                                extrp_mesh_a,
                                extrp_mesh_b,
                                redistribute_box_a_to_box_b_idx,
                                redistribute_box_a_to_box_b);
  } else if(mi->dim_mesh_a == 2 && mi->dim_mesh_b == 2) {
    _mesh_intersection_surf_surf(mi,
                                 extrp_mesh_a,
                                 extrp_mesh_b,
                                 redistribute_box_a_to_box_b_idx,
                                 redistribute_box_a_to_box_b);
  } else if(mi->dim_mesh_a == 2 && mi->dim_mesh_b == 1) {
    // On suppose que l'utilisateur met A = Vol et B = Surf
    _mesh_intersection_surf_line(mi,
                                 extrp_mesh_a,
                                 extrp_mesh_b,
                                 redistribute_box_a_to_box_b_idx,
                                 redistribute_box_a_to_box_b);
  } else {
    PDM_error(__FILE__, __LINE__, 0,
              "PDM_mesh_intersection_compute error : Cannot handle meshA with dim = %i and meshB = %i \n", mi->dim_mesh_a, mi->dim_mesh_b);
  }

  free(redistribute_box_a_to_box_b_idx);
  free(redistribute_box_a_to_box_b    );

  PDM_extract_part_free(extrp_mesh_a);
  PDM_extract_part_free(extrp_mesh_b);

}


void
PDM_mesh_intersection_part_set
(
  PDM_mesh_intersection_t  *mi,
  PDM_ol_mesh_t             i_mesh,
  int                       i_part,
  int                       n_cell,
  int                       n_face,
  int                       n_edge,
  int                       n_vtx,
  int                      *cell_face_idx,
  int                      *cell_face,
  int                      *face_edge_idx,
  int                      *face_edge,
  int                      *edge_vtx,
  int                      *face_vtx_idx,
  int                      *face_vtx,
  PDM_g_num_t              *cell_ln_to_gn,
  PDM_g_num_t              *face_ln_to_gn,
  PDM_g_num_t              *edge_ln_to_gn,
  PDM_g_num_t              *vtx_ln_to_gn,
  double                   *vtx_coord
)
{
  PDM_part_mesh_t* mesh = mi->mesh_a;
  if(i_mesh == PDM_OL_MESH_B) {
    mesh = mi->mesh_b;
  }

  PDM_part_mesh_n_entity_set(mesh, i_part, PDM_MESH_ENTITY_CELL  , n_cell);
  PDM_part_mesh_n_entity_set(mesh, i_part, PDM_MESH_ENTITY_FACE  , n_face);
  PDM_part_mesh_n_entity_set(mesh, i_part, PDM_MESH_ENTITY_EDGE  , n_edge);
  PDM_part_mesh_n_entity_set(mesh, i_part, PDM_MESH_ENTITY_VERTEX, n_vtx );

  PDM_part_mesh_connectivity_set(mesh, i_part, PDM_CONNECTIVITY_TYPE_CELL_FACE, cell_face, cell_face_idx, PDM_OWNERSHIP_USER);
  PDM_part_mesh_connectivity_set(mesh, i_part, PDM_CONNECTIVITY_TYPE_FACE_EDGE, face_edge, face_edge_idx, PDM_OWNERSHIP_USER);
  PDM_part_mesh_connectivity_set(mesh, i_part, PDM_CONNECTIVITY_TYPE_FACE_VTX , face_vtx , face_vtx_idx , PDM_OWNERSHIP_USER);
  PDM_part_mesh_connectivity_set(mesh, i_part, PDM_CONNECTIVITY_TYPE_EDGE_VTX , edge_vtx , NULL         , PDM_OWNERSHIP_USER);

  PDM_part_mesh_entity_ln_to_gn_set(mesh, i_part, PDM_MESH_ENTITY_CELL  , cell_ln_to_gn, PDM_OWNERSHIP_USER);
  PDM_part_mesh_entity_ln_to_gn_set(mesh, i_part, PDM_MESH_ENTITY_FACE  , face_ln_to_gn, PDM_OWNERSHIP_USER);
  PDM_part_mesh_entity_ln_to_gn_set(mesh, i_part, PDM_MESH_ENTITY_EDGE  , edge_ln_to_gn, PDM_OWNERSHIP_USER);
  PDM_part_mesh_entity_ln_to_gn_set(mesh, i_part, PDM_MESH_ENTITY_VERTEX, vtx_ln_to_gn , PDM_OWNERSHIP_USER);

  PDM_part_mesh_vtx_coord_set(mesh, i_part, vtx_coord, PDM_OWNERSHIP_USER);
}



void
PDM_mesh_intersection_free
(
 PDM_mesh_intersection_t* mi
)
{

  PDM_part_mesh_free(mi->mesh_a);
  PDM_part_mesh_free(mi->mesh_b);

  free(mi);
}


#ifdef __cplusplus
}
#endif /* __cplusplus */
