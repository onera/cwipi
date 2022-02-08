#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include <pdm_mpi.h>

#include "pdm.h"
#include "pdm_config.h"
#include "pdm_priv.h"
#include "pdm_part.h"
#include "pdm_poly_surf_gen.h"
#include "pdm_mesh_location.h"
#include "pdm_geom_elem.h"
#include "pdm_gnum.h"
#include "pdm_distrib.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_vtk.h"

/*============================================================================
 * Type definitions
 *============================================================================*/

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
     "  -t      <level>  Bounding boxes tolerance.\n\n"
     "  -n_part <level>  Number of partitions par process.\n\n"
     "  -p      <level>  Number of points to locate.\n\n"
     "  -octree          Use octree-based method.\n\n"
     "  -dbbree          Use dbbtree-based method.\n\n"
     "  -parmetis        Call ParMETIS.\n\n"
     "  -pt-scocth       Call PT-Scotch.\n\n"
     "  -h               This message.\n\n");

  exit(exit_code);
}

/**
 *
 * \brief  Read arguments from the command line
 *
 * \param [in]      argc       Number of arguments
 * \param [in]      argv       Arguments
 * \param [inout]   n_vtx_seg  Number of vertices on the cube side
 * \param [inout]   length     Cube length
 * \param [inout]   tolerance  Bounding boxes tolerance
 * \param [inout]   n_part     Number of partitions par process
 * \param [inout]   post       Ensight outputs status
 * \param [inout]   method     Partitioner (1 ParMETIS, 2 Pt-Scotch)
 *
 */

static void
_read_args(int                         argc,
           char                      **argv,
           PDM_g_num_t                *n_vtx_seg,
           double                     *lengthx,
           double                     *lengthy,
           double                     *xmin,
           double                     *ymin,
           double                     *tolerance,
           PDM_g_num_t                *n_pts)
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
        long _n_vtx_seg = atol(argv[i]);
        *n_vtx_seg = (PDM_g_num_t) _n_vtx_seg;
      }
    }

    else if (strcmp(argv[i], "-lx") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *lengthx = atof(argv[i]);
    }

    else if (strcmp(argv[i], "-ly") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *lengthy = atof(argv[i]);
    }

    else if (strcmp(argv[i], "-xmin") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *xmin = atof(argv[i]);
    }

    else if (strcmp(argv[i], "-ymin") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *ymin = atof(argv[i]);
    }

    else if (strcmp(argv[i], "-t") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *tolerance = atof(argv[i]);
    }

    else if (strcmp(argv[i], "-p") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *n_pts = atoi(argv[i]);
      }
    }

    else
      _usage(EXIT_FAILURE);
    i++;
  }
}



static void
_gen_point_cloud
(
 PDM_MPI_Comm        comm,
 const PDM_g_num_t   gn_pts,
 const double        xmin,
 const double        ymin,
 const double        lengthx,
 const double        lengthy,
 int                *n_pts,
 double            **pts_coord,
 PDM_g_num_t       **pts_g_num
 )
 {
  int i_rank, n_rank;
  PDM_MPI_Comm_rank (comm, &i_rank);
  PDM_MPI_Comm_size (comm, &n_rank);


  *n_pts = (int) (gn_pts / n_rank);
  if (i_rank < gn_pts % n_rank) {
    (*n_pts)++;
  }

  srand(0);

  PDM_g_num_t* distrib_pts = PDM_compute_entity_distribution(comm, (*n_pts));
  for(int i = 0; i < 2 * distrib_pts[i_rank]; ++i) {
    rand();
  }

  *pts_coord = malloc (sizeof(double) * (*n_pts) * 3);
  for (int i = 0; i < *n_pts; i++) {
    (*pts_coord)[3*i    ] = xmin + lengthx * (double) rand() / ((double) RAND_MAX);
    (*pts_coord)[3*i + 1] = ymin + lengthy * (double) rand() / ((double) RAND_MAX);
    (*pts_coord)[3*i + 2] = 0.;
  }

  double length = 0.5*(lengthx + lengthy);
  PDM_gen_gnum_t *gen_gnum = PDM_gnum_create (3, 1, PDM_FALSE, 1e-3, comm, PDM_OWNERSHIP_USER);
  double *char_length = malloc (sizeof(double) * (*n_pts));
  for (int i = 0; i < *n_pts; i++) {
    char_length[i] = length * 1e-6;
  }

  PDM_gnum_set_from_coords (gen_gnum, 0, *n_pts, *pts_coord, char_length);

  PDM_gnum_compute (gen_gnum);

  *pts_g_num = PDM_gnum_get (gen_gnum, 0);

  PDM_gnum_free (gen_gnum);
  free (char_length);

  free(distrib_pts);
 }


 static void
_get_connectivity
(
 PDM_part_t    *ppart,
 int            n_part,
 int          **nFace,
 int         ***faceEdgeIdx,
 int         ***faceEdge,
 int         ***faceVtxIdx,
 int         ***faceVtx,
 PDM_g_num_t ***faceLNToGN,
 int          **nEdge,
 int         ***edgeVtxIdx,
 int         ***edgeVtx,
 int          **nVtx,
 double      ***vtxCoord,
 PDM_g_num_t ***vtxLNToGN
 )
{
  *nFace = (int *) malloc(sizeof(int) * n_part);
  *faceEdgeIdx = (int **) malloc(sizeof(int *) * n_part);
  *faceEdge = (int **) malloc(sizeof(int *) * n_part);
  *faceVtxIdx = (int **) malloc(sizeof(int *) * n_part);
  *faceVtx = (int **) malloc(sizeof(int *) * n_part);
  *faceLNToGN = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t *) * n_part);

  *nEdge = (int *) malloc(sizeof(int) * n_part);
  *edgeVtxIdx = (int **) malloc(sizeof(int *) * n_part);
  *edgeVtx = (int **) malloc(sizeof(int *) * n_part);

  *nVtx = (int *) malloc(sizeof(int) * n_part);
  *vtxCoord = (double **) malloc(sizeof(double *) * n_part);
  *vtxLNToGN = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t *) * n_part);


  for (int ipart = 0; ipart < n_part; ipart++) {

    int _nFace;
    int _nEdge;
    int _nEdgePartBound;
    int _nVtx;
    int _nProc;
    int _nTPart;
    int _sFaceEdge;
    int _sEdgeVtx;
    int _sEdgeGroup;
    int _nEdgeGroup2;

    PDM_part_part_dim_get (ppart,
                           ipart,
                           &_nFace,
                           &_nEdge,
                           &_nEdgePartBound,
                           &_nVtx,
                           &_nProc,
                           &_nTPart,
                           &_sFaceEdge,
                           &_sEdgeVtx,
                           &_sEdgeGroup,
                           &_nEdgeGroup2);

    int         *_faceTag;
    int         *_faceEdgeIdx;
    int         *_faceEdge;
    PDM_g_num_t *_faceLNToGN;
    int         *_edgeTag;
    int         *_edgeFace;
    int         *_edgeVtxIdx;
    int         *_edgeVtx;
    PDM_g_num_t *_edgeLNToGN;
    int         *_edgePartBoundProcIdx;
    int         *_edgePartBoundPartIdx;
    int         *_edgePartBound;
    int         *_vtxTag;
    double      *_vtx;
    PDM_g_num_t *_vtxLNToGN;
    int         *_edgeGroupIdx;
    int         *_edgeGroup;
    PDM_g_num_t *_edgeGroupLNToGN;

    PDM_part_part_val_get (ppart,
                           ipart,
                           &_faceTag,
                           &_faceEdgeIdx,
                           &_faceEdge,
                           &_faceLNToGN,
                           &_edgeTag,
                           &_edgeFace,
                           &_edgeVtxIdx,
                           &_edgeVtx,
                           &_edgeLNToGN,
                           &_edgePartBoundProcIdx,
                           &_edgePartBoundPartIdx,
                           &_edgePartBound,
                           &_vtxTag,
                           &_vtx,
                           &_vtxLNToGN,
                           &_edgeGroupIdx,
                           &_edgeGroup,
                           &_edgeGroupLNToGN);

    /* Faces */
    (*nFace)[ipart] = _nFace;
    (*faceEdgeIdx)[ipart] = (int *) malloc(sizeof(int) * (_nFace + 1));
    (*faceEdge)[ipart] = (int *) malloc(sizeof(int) * _sFaceEdge);
    (*faceVtxIdx)[ipart] = (int *) malloc(sizeof(int) * (_nFace + 1));
    (*faceVtx)[ipart] = (int *) malloc(sizeof(int) * _sFaceEdge);
    (*faceLNToGN)[ipart] = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * _nFace);

    memcpy ((*faceEdgeIdx)[ipart], _faceEdgeIdx, (_nFace + 1) * sizeof(int));
    memcpy ((*faceEdge)[ipart], _faceEdge, _sFaceEdge * sizeof(int));
    memcpy ((*faceVtxIdx)[ipart], _faceEdgeIdx, (_nFace + 1) * sizeof(int));
    memcpy ((*faceLNToGN)[ipart], _faceLNToGN, _nFace * sizeof(PDM_g_num_t));

    /* Edges */
    (*nEdge)[ipart] = _nEdge;
    (*edgeVtxIdx) [ipart] = (int *) malloc(sizeof(int) * (_nEdge + 1));
    (*edgeVtx)[ipart] = (int *) malloc(sizeof(int) * _sEdgeVtx);

    memcpy ((*edgeVtxIdx)[ipart], _edgeVtxIdx, (_nEdge + 1) * sizeof(int));
    memcpy ((*edgeVtx)[ipart], _edgeVtx, _sEdgeVtx * sizeof(int));

    /* Vertices */
    (*nVtx)[ipart] = _nVtx;
    (*vtxCoord)[ipart] = (double *) malloc(sizeof(double) * (3 * _nVtx));
    (*vtxLNToGN)[ipart] = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * _nVtx);

    memcpy ((*vtxCoord)[ipart], _vtx, 3 *_nVtx * sizeof(double));
    memcpy ((*vtxLNToGN)[ipart], _vtxLNToGN, _nVtx * sizeof(PDM_g_num_t));


    /* Compute face-vtx connectivity */
    int *_faceVtx = (*faceVtx)[ipart];

    int *vtxEdgeIdx = (int *) malloc(sizeof(int) * (_nVtx + 1));

    for (int i = 0; i < _nVtx + 1; i++) {
      vtxEdgeIdx[i] = 0;
    }

    for (int i = 0; i < _nEdge; i++) {
      int ivtx1 = _edgeVtx[2*i];
      int ivtx2 = _edgeVtx[2*i + 1];

      vtxEdgeIdx[ivtx1] += 1;
      vtxEdgeIdx[ivtx2] += 1;
    }

    for (int i = 1; i < _nVtx + 1; i++) {
      vtxEdgeIdx[i] = vtxEdgeIdx[i] + vtxEdgeIdx[i-1];
    }

    int *vtxEdge = (int *) malloc(sizeof(int) * vtxEdgeIdx[_nVtx]);
    int *vtxEdgeN = (int *) malloc(sizeof(int) * _nVtx);
    for (int i = 0; i < _nVtx; i++) {
      vtxEdgeN[i] = 0;
    }

    for (int i = 0; i < _nEdge; i++) {
      int ivtx1 = _edgeVtx[2*i] - 1;
      int ivtx2 = _edgeVtx[2*i + 1] - 1;
      int iedge = i + 1;

      vtxEdge[vtxEdgeIdx[ivtx1] + vtxEdgeN[ivtx1]] = iedge;
      vtxEdge[vtxEdgeIdx[ivtx2] + vtxEdgeN[ivtx2]] = iedge;
      vtxEdgeN[ivtx1] += 1;
      vtxEdgeN[ivtx2] += 1;
    }
    free(vtxEdgeN);

    for (int i = 0; i < _nFace; i++) {
      int idx = _faceEdgeIdx[i];
      int __nEdge = _faceEdgeIdx[i+1] - idx;
      int *_edges = _faceEdge + idx;
      int *_vertices = _faceVtx + idx;

      int edge_cur = _edges[0];
      int vtx_deb =  _edgeVtx[2*(edge_cur - 1)];
      _vertices[0] = vtx_deb;
      int vtx_cur =  _edgeVtx[2*(edge_cur - 1) + 1];
      int idxVtx = 0;

      while (vtx_deb != vtx_cur) {
        _vertices[++idxVtx] = vtx_cur;
        int find_vtx = 0;

        for (int j = vtxEdgeIdx[vtx_cur - 1]; j <  vtxEdgeIdx[vtx_cur]; j++) {
          for (int k = 0; k < __nEdge; k++) {
            if ((_edges[k] == vtxEdge[j]) && (_edges[k] != edge_cur)) {
              edge_cur = _edges[k];
              if (_edgeVtx[2*(_edges[k]-1)] == vtx_cur) {
                vtx_cur = _edgeVtx[2*(_edges[k]-1) + 1];
              }
              else {
                vtx_cur = _edgeVtx[2*(_edges[k]-1)];
              }
              find_vtx = 1;
              break;
            }
          }
          if (find_vtx)
            break;
        }
        if (!find_vtx) {
          PDM_error(__FILE__, __LINE__, 0,"Error to compute vtxedge !!!!\n");
          abort();
        }
      }
    }

    free (vtxEdge);
    free (vtxEdgeIdx);

  }
}


static void
_create_split_mesh
(
 int                 activeRank,
 PDM_MPI_Comm        pdm_mpi_comm,
 double              xmin,
 double              ymin,
 PDM_g_num_t         nVtxSeg,
 double              length,
 int                 n_part,
 PDM_part_split_t    method,
 int                 haveRandom,
 int                 initRandom,
 PDM_g_num_t        *nGFace,
 PDM_g_num_t        *nGVtx,
 int               **nFace,
 int              ***faceEdgeIdx,
 int              ***faceEdge,
 int              ***faceVtxIdx,
 int              ***faceVtx,
 PDM_g_num_t      ***faceLNToGN,
 int               **nEdge,
 int              ***edgeVtxIdx,
 int              ***edgeVtx,
 int               **nVtx,
 double           ***vtxCoord,
 PDM_g_num_t      ***vtxLNToGN
 )
{
  int i_rank;
  int numProcs;

  if (activeRank) {

    PDM_MPI_Comm_rank (pdm_mpi_comm, &i_rank);
    PDM_MPI_Comm_size (pdm_mpi_comm, &numProcs);

    double       xmax = xmin + length;
    double       ymax = ymin + length;
    PDM_g_num_t  nx = nVtxSeg;
    PDM_g_num_t  ny = nVtxSeg;

    int          dNFace;
    int          dNVtx;
    int          dNEdge;
    int         *dFaceVtxIdx;
    PDM_g_num_t *dFaceVtx;
    double      *dVtxCoord;
    PDM_g_num_t *dFaceEdge;
    PDM_g_num_t *dEdgeVtx;
    PDM_g_num_t *dEdgeFace;
    int          nEdgeGroup;
    int         *dEdgeGroupIdx;
    PDM_g_num_t *dEdgeGroup;

    /*
     *  Create mesh
     */
    PDM_g_num_t nGEdge;

    PDM_poly_surf_gen (pdm_mpi_comm,
                       xmin,
                       xmax,
                       ymin,
                       ymax,
                       haveRandom,
                       initRandom,
                       nx,
                       ny,
                       nGFace,
                       nGVtx,
                       &nGEdge,
                       &dNVtx,
                       &dVtxCoord,
                       &dNFace,
                       &dFaceVtxIdx,
                       &dFaceVtx,
                       &dFaceEdge,
                       &dNEdge,
                       &dEdgeVtx,
                       &dEdgeFace,
                       &nEdgeGroup,
                       &dEdgeGroupIdx,
                       &dEdgeGroup);

    /*
     *  Create mesh partitions
     */
    int have_dCellPart = 0;

    int *dCellPart   = (int *) malloc (dNFace * sizeof(int));
    int *dEdgeVtxIdx = (int *) malloc ((dNEdge + 1) * sizeof(int));

    dEdgeVtxIdx[0] = 0;
    for (int i = 0; i < dNEdge; i++) {
      dEdgeVtxIdx[i+1] = 2 + dEdgeVtxIdx[i];
    }

    /*
     *  Split mesh
     */
    // int ppartId;

    int nPropertyCell = 0;
    int *renum_properties_cell = NULL;
    int nPropertyFace = 0;
    int *renum_properties_face = NULL;

    PDM_part_t *ppart = PDM_part_create (pdm_mpi_comm,
                                         method,
                                         "PDM_PART_RENUM_CELL_NONE",
                                         "PDM_PART_RENUM_FACE_NONE",
                                         nPropertyCell,
                                         renum_properties_cell,
                                         nPropertyFace,
                                         renum_properties_face,
                                         n_part,
                                         dNFace,
                                         dNEdge,
                                         dNVtx,
                                         nEdgeGroup,
                                         NULL,
                                         NULL,
                                         NULL,
                                         NULL,
                                         have_dCellPart,
                                         dCellPart,
                                         dEdgeFace,
                                         dEdgeVtxIdx,
                                         dEdgeVtx,
                                         NULL,
                                         dVtxCoord,
                                         NULL,
                                         dEdgeGroupIdx,
                                         dEdgeGroup);

    free (dCellPart);

    free (dVtxCoord);
    free (dFaceVtxIdx);
    free (dFaceVtx);
    free (dFaceEdge);
    free (dEdgeVtxIdx);
    free (dEdgeVtx);
    free (dEdgeFace);
    free (dEdgeGroupIdx);
    free (dEdgeGroup);

    _get_connectivity (ppart,
                       n_part,
                       nFace,
                       faceEdgeIdx,
                       faceEdge,
                       faceVtxIdx,
                       faceVtx,
                       faceLNToGN,
                       nEdge,
                       edgeVtxIdx,
                       edgeVtx,
                       nVtx,
                       vtxCoord,
                       vtxLNToGN);

    PDM_part_free (ppart);
  }

  else {
    *nFace = (int *) malloc(sizeof(int) * n_part);
    *faceEdgeIdx = (int **) malloc(sizeof(int *) * n_part);
    *faceEdge = (int **) malloc(sizeof(int *) * n_part);
    *faceVtxIdx = (int **) malloc(sizeof(int *) * n_part);
    *faceVtx = (int **) malloc(sizeof(int *) * n_part);
    *faceLNToGN = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t *) * n_part);

    *nEdge = (int *) malloc(sizeof(int) * n_part);
    *edgeVtxIdx = (int **) malloc(sizeof(int *) * n_part);
    *edgeVtx = (int **) malloc(sizeof(int *) * n_part);

    *nVtx = (int *) malloc(sizeof(int) * n_part);
    *vtxCoord = (double **) malloc(sizeof(double *) * n_part);
    *vtxLNToGN = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t *) * n_part);

    for (int ipart = 0; ipart < n_part; ipart++) {
      int _nFace = 0;
      int _sFaceEdge = 0;
      (*nFace)[ipart] = _nFace;
      (*faceEdgeIdx)[ipart] = (int *) malloc(sizeof(int) * (_nFace + 1));
      (*faceEdgeIdx)[ipart][0] = 0;
      (*faceEdge)[ipart] = (int *) malloc(sizeof(int) * _sFaceEdge);
      (*faceVtxIdx)[ipart] = (int *) malloc(sizeof(int) * (_nFace + 1));
      (*faceVtxIdx)[ipart][0] = 0;
      (*faceVtx)[ipart] = (int *) malloc(sizeof(int) * _sFaceEdge);
      (*faceLNToGN)[ipart] = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * _nFace);

      int _nEdge = 0;
      int _sEdgeVtx = 0;
      (*nEdge)[ipart] = _nEdge;
      (*edgeVtxIdx)[ipart] = (int *) malloc(sizeof(int) * (_nEdge + 1));
      (*edgeVtxIdx)[ipart][0] = 0;
      (*edgeVtx)[ipart] = (int *) malloc(sizeof(int) * _sEdgeVtx);

      int _nVtx = 0;
      (*nVtx)[ipart] = _nVtx;
      (*vtxCoord)[ipart] = (double *) malloc(sizeof(double) * (3 * _nVtx));
      (*vtxLNToGN)[ipart] = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * _nVtx);
    }
  }

  PDM_MPI_Bcast (nGFace, 1, PDM__PDM_MPI_G_NUM, 0, PDM_MPI_COMM_WORLD);
  PDM_MPI_Bcast (nGVtx, 1, PDM__PDM_MPI_G_NUM, 0, PDM_MPI_COMM_WORLD);
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
  PDM_g_num_t n_vtx_seg = 20;
  double      lengthx   = 0.6;
  double      lengthy   = 0.24;
  double      xmin      = 0.05;
  double      ymin      = 0.7;
  double      tolerance = 1e-6;
  int         n_part    = 1;
  PDM_g_num_t n_pts     = 40;

  PDM_mesh_location_method_t loc_method = PDM_MESH_LOCATION_OCTREE;
  PDM_part_split_t part_method  = PDM_PART_SPLIT_PARMETIS;


  /*
   *  Read args
   */

  _read_args (argc,
              argv,
              &n_vtx_seg,
              &lengthx,
              &lengthy,
              &xmin,
              &ymin,
              &tolerance,
              &n_pts);

  /*
   *  Init
   */

  int i_rank;
  int n_rank;

  PDM_MPI_Init (&argc, &argv);
  PDM_MPI_Comm_rank (PDM_MPI_COMM_WORLD, &i_rank);
  PDM_MPI_Comm_size (PDM_MPI_COMM_WORLD, &n_rank);

  char filename[999];


  /*
   *  Create partitionned surface mesh
   */

  if (i_rank == 0) {
    printf("-- Build surface mesh\n");
    fflush(stdout);
  }

  PDM_g_num_t   nGFace;
  PDM_g_num_t   nGVtx;
  int          *nFace       = NULL;
  PDM_g_num_t **faceLNToGN  = NULL;
  int         **faceEdgeIdx = NULL;
  int         **faceEdge    = NULL;
  int         **faceVtxIdx  = NULL;
  int         **faceVtx     = NULL;
  int          *nEdge       = NULL;
  int         **edgeVtxIdx  = NULL;
  int         **edgeVtx     = NULL;
  int          *nVtx        = NULL;
  double      **vtxCoord    = NULL;
  PDM_g_num_t **vtxLNToGN   = NULL;

  _create_split_mesh (1,
                      PDM_MPI_COMM_WORLD,
                      0.,
                      0.,
                      n_vtx_seg,
                      1.,
                      n_part,
                      part_method,
                      1,
                      0,
                      &nGFace,
                      &nGVtx,
                      &nFace,
                      &faceEdgeIdx,
                      &faceEdge,
                      &faceVtxIdx,
                      &faceVtx,
                      &faceLNToGN,
                      &nEdge,
                      &edgeVtxIdx,
                      &edgeVtx,
                      &nVtx,
                      &vtxCoord,
                      &vtxLNToGN);

  /* Write mesh */
  sprintf(filename, "00_part_mesh_%d.vtk", i_rank);
  PDM_vtk_write_polydata (filename,
                          nVtx[0],
                          vtxCoord[0],
                          vtxLNToGN[0],
                          nFace[0],
                          faceVtxIdx[0],
                          faceVtx[0],
                          faceLNToGN[0],
                          NULL);



  /************************
   *
   * Point cloud definition
   *
   ************************/
  if (i_rank == 0) {
    printf("-- Point cloud\n");
    fflush(stdout);
  }

  int n_pts_l;
  double      *pts_coord = NULL;
  PDM_g_num_t *pts_g_num = NULL;
  _gen_point_cloud (PDM_MPI_COMM_WORLD,
                    n_pts,
                    xmin,
                    ymin,
                    lengthx,
                    lengthy,
                    &n_pts_l,
                    &pts_coord,
                    &pts_g_num);


  /* Write cloud */
  sprintf(filename, "01_part_cloud_%d.vtk", i_rank);
  PDM_vtk_write_point_cloud (filename,
                             nVtx[0],
                             vtxCoord[0],
                             vtxLNToGN[0],
                             NULL);



  /************************
   *
   * Mesh location struct initializaiton
   *
   ************************/

  PDM_mesh_location_t* mesh_loc = PDM_mesh_location_create (PDM_MESH_NATURE_MESH_SETTED,//???
                                                            1,//const int n_point_cloud,
                                                            PDM_MPI_COMM_WORLD);

  /* Set point cloud(s) */
  PDM_mesh_location_n_part_cloud_set (mesh_loc,
                                      0,//i_point_cloud,
                                      1);//n_part

  PDM_mesh_location_cloud_set (mesh_loc,
                               0,//i_point_cloud,
                               0,//i_part,
                               n_pts_l,
                               pts_coord,
                               pts_g_num);

  PDM_mesh_location_mesh_global_data_set (mesh_loc,
                                          n_part);

  /* Set mesh */
  for (int ipart = 0; ipart < n_part; ipart++) {
    PDM_mesh_location_part_set_2d (mesh_loc,
                                   ipart,
                                   nFace[ipart],
                                   faceEdgeIdx[ipart],
                                   faceEdge[ipart],
                                   faceLNToGN[ipart],
                                   nEdge[ipart],
                                   edgeVtxIdx[ipart],
                                   edgeVtx[ipart],
                                   NULL,//edgeLNToGN[ipart],
                                   nVtx[ipart],
                                   vtxCoord[ipart],
                                   vtxLNToGN[ipart]);
  }




  /* Set location parameters */
  PDM_mesh_location_tolerance_set (mesh_loc,
                                   tolerance);

  PDM_mesh_location_method_set (mesh_loc,
                                loc_method);


  /*
   * Compute location
   */
  if (i_rank == 0) {
    printf("-- Locate\n");
    fflush(stdout);
  }

  PDM_mesh_location_compute (mesh_loc);

  PDM_mesh_location_dump_times (mesh_loc);



  /*
   * Finalize
   */
  PDM_mesh_location_free (mesh_loc, 0);

  for (int ipart = 0; ipart < n_part; ipart++) {
    free(faceEdgeIdx[ipart]);
    free(faceEdge[ipart]);
    free(faceVtxIdx[ipart]);
    free(faceVtx[ipart]);
    free(faceLNToGN[ipart]);
    free(edgeVtxIdx[ipart]);
    free(edgeVtx[ipart]);
    free(vtxCoord[ipart]);
    free(vtxLNToGN[ipart]);
  }
  free(faceVtxIdx);
  free(faceVtx);
  free(nFace);
  free(faceEdgeIdx);
  free(faceEdge);
  free(faceLNToGN);
  free(nEdge);
  free(edgeVtxIdx);
  free(edgeVtx);
  free(nVtx);
  free(vtxCoord);
  free(vtxLNToGN);

  free (pts_coord);
  free (pts_g_num);

  PDM_MPI_Finalize();

  if (i_rank == 0) {
    printf("-- End\n");
    fflush(stdout);
  }

  return 0;
}
