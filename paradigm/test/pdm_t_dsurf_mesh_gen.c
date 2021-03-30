#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_mpi.h"
#include "pdm_config.h"
#include "pdm_part.h"

#include "pdm_writer.h"
#include "pdm_geom_elem.h"
#include "pdm_part_to_block.h"
#include "pdm_poly_surf_gen.h"
#include "pdm_overlay.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_mpi_node_first_rank.h"
#include "pdm_priv.h"

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
_usage
(
int exit_code
)
{
  PDM_printf
    ("\n"
     "  Usage: \n\n"
     "  -nA       <value>  Cube A - Number of vertices on the side.\n\n"
     "  -lA       <value>  Cube A - length.\n\n"
     "  -xminA    <value>  Cube A - Origin X.\n\n"
     "  -yminA    <value>  Cube A - Origin Y.\n\n"
     "  -yminA    <value>  Cube A - Origin Y.\n\n"
     "  -n_partA  <level>  Cube A - Number of partitions par process.\n\n"
     "  -nB       <value>  Cube A - Number of vertices on the side (default : nA+4).\n\n"
     "  -lB       <value>  Cube A - length.\n\n"
     "  -xminB    <value>  Cube A - Origin X.\n\n"
     "  -yminB    <value>  Cube A - Origin Y.\n\n"
     "  -yminB    <value>  Cube A - Origin Y.\n\n"
     "  -n_partB  <level>  Cube A - Number of partitions par process.\n\n"
     "  -post              Ensight output.\n\n"
     "  -no_random         No random to define points coordinates.\n\n"
     "  -random_time_init  Intialize random with the current time.\n\n"
     "  -parmetis          Call ParMETIS.\n\n"
     "  -pt-scocth         Call PT-Scotch.\n\n"
     "  -n_proc_data       Number of processses where there are data.\n\n"
     "  -h                 This message.\n\n");

  exit (exit_code);
}


/**
 *
 * \brief  Read arguments from the command line
 *
 * \param [in]      argc     Number of arguments
 * \param [in]      argv     Arguments
 * \param [inout]   nVtxSeg  Number of vertices on the cube side
 * \param [inout]   length   Cube length
 * \param [inout]   n_part   Number of partitions par process
 * \param [inout]   post     Ensight outputs status
 * \param [inout]   method   Partitioner (1 ParMETIS, 2 Pt-Scotch)
 *
 */

static void
_read_args
(
 int            argc,
 char         **argv,
 PDM_g_num_t  *nVtxSegA,
 double        *lengthA,
 double        *xminA,
 double        *yminA,
 int           *n_partA,
 int           *post,
 int           *method,
 int           *haveRandom,
 int           *randomTimeInit,
 int           *randomMeshAInit,
 int           *nProcData
)
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {

    if (strcmp (argv[i], "-h") == 0)
      _usage(EXIT_SUCCESS);

    else if (strcmp (argv[i], "-nA") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        long _nVtxSeg = atol (argv[i]);
        *nVtxSegA = (PDM_g_num_t) _nVtxSeg;
      }
    }
    else if (strcmp (argv[i], "-lA") == 0) {
      i++;
      if (i >= argc)
        _usage (EXIT_FAILURE);
      else
        *lengthA = atof (argv[i]);
    }
    else if (strcmp (argv[i], "-xminA") == 0) {
      i++;
      if (i >= argc)
        _usage (EXIT_FAILURE);
      else
        *xminA = atof (argv[i]);
    }
    else if (strcmp (argv[i], "-yminA") == 0) {
      i++;
      if (i >= argc)
        _usage (EXIT_FAILURE);
      else
        *yminA = atof (argv[i]);
    }
    else if (strcmp (argv[i], "-n_partA") == 0) {
      i++;
      if (i >= argc)
        _usage (EXIT_FAILURE);
      else {
        *n_partA = atoi (argv[i]);
      }
    }
    else if (strcmp (argv[i], "-randomMeshAInit") == 0) {
      i++;
      if (i >= argc)
        _usage (EXIT_FAILURE);
      else {
        *randomMeshAInit = atoi (argv[i]);
      }
    }
    else if (strcmp (argv[i], "-no_random") == 0) {
      *haveRandom = 0;
    }
    else if (strcmp (argv[i], "-random_time_init") == 0) {
      *randomTimeInit = 1;
    }
    else if (strcmp (argv[i], "-post") == 0) {
      *post = 1;
    }
    else if (strcmp (argv[i], "-pt-scotch") == 0) {
      *method = 2;
    }
    else if (strcmp (argv[i], "-parmetis") == 0) {
      *method = 1;
    }
    else if (strcmp (argv[i], "-n_proc_data") == 0) {
      i++;
      if (i >= argc)
        _usage (EXIT_FAILURE);
      else {
        *nProcData = atoi (argv[i]);
      }
    }
    else
      _usage (EXIT_FAILURE);
    i++;
  }
}


/**
 *
 * \brief  Compute faceVtx connectivity
 *
 * \param [in]      ppartId  ppart identifier
 * \param [in]      n_part    Number of partitions
 *
 * \return          faceVtx connectivity for each partition of each mesh
 */

static void
_compute_faceVtx
(
 int           ppartId,
 int            n_part,
 int          **nFace,
 int         ***faceVtxIdx,
 int         ***faceVtx,
 PDM_g_num_t ***faceLNToGN,
 int          **nVtx,
 double      ***vtxCoord,
 PDM_g_num_t ***vtxLNToGN
)
{

  *nFace = (int *) malloc(sizeof(int) * n_part);
  *faceVtxIdx = (int **) malloc(sizeof(int *) * n_part);
  *faceVtx = (int **) malloc(sizeof(int *) * n_part);
  *faceLNToGN = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t *) * n_part);

  *nVtx = (int *) malloc(sizeof(int) * n_part);
  *vtxCoord = (double **) malloc(sizeof(double *) * n_part);
  *vtxLNToGN = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t *) * n_part);

  int id_ppart = ppartId;

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

    PDM_part_part_dim_get (id_ppart,
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

    int          *_faceTag;
    int          *_faceEdgeIdx;
    int          *_faceEdge;
    PDM_g_num_t *_faceLNToGN;
    int          *_edgeTag;
    int          *_edgeFace;
    int          *_edgeVtxIdx;
    int          *_edgeVtx;
    PDM_g_num_t *_edgeLNToGN;
    int          *_edgePartBoundProcIdx;
    int          *_edgePartBoundPartIdx;
    int          *_edgePartBound;
    int          *_vtxTag;
    double       *_vtx;
    PDM_g_num_t *_vtxLNToGN;
    int          *_edgeGroupIdx;
    int          *_edgeGroup;
    PDM_g_num_t  *_edgeGroupLNToGN;

    PDM_part_part_val_get (id_ppart,
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

    (*nFace)[ipart] = _nFace;
    (*faceVtxIdx) [ipart] = (int *) malloc(sizeof(int) * (_nFace + 1));
    (*faceVtx)[ipart] = (int *) malloc(sizeof(int) * _sFaceEdge);
    (*faceLNToGN)[ipart] = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * _nFace);

    memcpy ((*faceVtxIdx)[ipart], _faceEdgeIdx, (_nFace + 1) * sizeof(int));
    memcpy ((*faceLNToGN)[ipart], _faceLNToGN, _nFace * sizeof(PDM_g_num_t));

    (*nVtx)[ipart] = _nVtx;
    (*vtxCoord)[ipart] = (double *) malloc(sizeof(double) * (3 * _nVtx));
    (*vtxLNToGN)[ipart] = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * _nVtx);

    memcpy ((*vtxCoord)[ipart], _vtx, 3 *_nVtx * sizeof(double));
    memcpy ((*vtxLNToGN)[ipart], _vtxLNToGN, _nVtx * sizeof(PDM_g_num_t));

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


/**
 *
 * \brief  Create and split Mesh
 *
 * \param [in]      nVtxSeg  Number of arguments
 * \param [in]      length   Lenght of square
 * \param [in]      n_part    Number to obtain on this processus
 * \param [in]      post     mesh export status
 * \param [in]      method   Split method
 *
 * \return ppartId  ppart identifier
 *
 */

static void
_create_split_mesh
(
 int               activeRank,
 PDM_MPI_Comm      pdm_mpi_comm,
 double            xmin,
 double            ymin,
 PDM_g_num_t       nVtxSeg,
 double            length,
 int               n_part,
 PDM_part_split_t   method,
 int               haveRandom,
 int               initRandom,
 PDM_g_num_t      *nGFace,
 PDM_g_num_t      *nGVtx,
 int            **nFace,
 int            ***faceVtxIdx,
 int            ***faceVtx,
 PDM_g_num_t    ***faceLNToGN,
 int            **nVtx,
 double         ***vtxCoord,
 PDM_g_num_t    ***vtxLNToGN
)
{
  struct timeval t_elaps_debut;

  int i_rank;
  int numProcs;

  if (activeRank) {

    PDM_MPI_Comm_rank (pdm_mpi_comm, &i_rank);
    PDM_MPI_Comm_size (pdm_mpi_comm, &numProcs);

    double        xmax = xmin + length;
    double        ymax = ymin + length;
    PDM_g_num_t  nx = nVtxSeg;
    PDM_g_num_t  ny = nVtxSeg;

    int           dNFace;
    int           dNVtx;
    int           dNEdge;
    int          *dFaceVtxIdx;
    PDM_g_num_t *dFaceVtx;
    double       *dVtxCoord;
    PDM_g_num_t *dFaceEdge;
    PDM_g_num_t *dEdgeVtx;
    PDM_g_num_t *dEdgeFace;
    int           nEdgeGroup;
    int          *dEdgeGroupIdx;
    PDM_g_num_t   *dEdgeGroup;


    /*
     *  Create mesh i
     */

    gettimeofday(&t_elaps_debut, NULL);

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

    struct timeval t_elaps_fin;

    gettimeofday (&t_elaps_fin, NULL);

    long tranche_elapsed = (t_elaps_fin.tv_usec + 1000000 * t_elaps_fin.tv_sec) -
      (t_elaps_debut.tv_usec + 1000000 *
       t_elaps_debut.tv_sec);
    long tranche_elapsed_max = tranche_elapsed;
    double t_elapsed = (double) tranche_elapsed_max/1000000.;
    if (i_rank == 0)
      PDM_printf("[%d] Temps dans creeMaillagePolygone2D : %12.5e\n",
                 i_rank, t_elapsed);

    /*
     *  Create mesh partitions
     */

    int have_dCellPart = 0;

    int *dCellPart = (int *) malloc (dNFace*sizeof(int));
    int *dEdgeVtxIdx = (int *) malloc ((dNEdge+1)*sizeof(int));

    dEdgeVtxIdx[0] = 0;
    for (int i = 0; i < dNEdge; i++) {
      dEdgeVtxIdx[i+1] = 2 + dEdgeVtxIdx[i];
    }

    /*
     *  Split mesh i
     */

    int ppartId;

    int nPropertyCell = 0;
    int *renum_properties_cell = NULL;
    int nPropertyFace = 0;
    int *renum_properties_face = NULL;

    // PDM_part_create (&ppartId,
    //                  pdm_mpi_comm,
    //                  method,
    //                  "PDM_PART_RENUM_CELL_NONE",
    //                  "PDM_PART_RENUM_FACE_NONE",
    //                  nPropertyCell,
    //                  renum_properties_cell,
    //                  nPropertyFace,
    //                  renum_properties_face,
    //                  n_part,
    //                  dNFace,
    //                  dNEdge,
    //                  dNVtx,
    //                  nEdgeGroup,
    //                  NULL,
    //                  NULL,
    //                  NULL,
    //                  NULL,
    //                  have_dCellPart,
    //                  dCellPart,
    //                  dEdgeFace,
    //                  dEdgeVtxIdx,
    //                  dEdgeVtx,
    //                  NULL,
    //                  dVtxCoord,
    //                  NULL,
    //                  dEdgeGroupIdx,
    //                  dEdgeGroup);

    printf("dNFace = %i | dNEdge = %i | dNVtx = %i \n", dNFace, dNEdge, dNVtx);
    PDM_part_create (&ppartId,
                     pdm_mpi_comm,
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
                     dFaceVtxIdx,
                     dFaceEdge,
                     NULL,
                     NULL,
                     have_dCellPart,
                     dCellPart,
                     NULL,
                     dEdgeVtxIdx,
                     dEdgeVtx,
                     NULL,
                     dVtxCoord,
                     NULL,
                     dEdgeGroupIdx,
                     dEdgeGroup);
    free (dCellPart);

    double  *elapsed = NULL;
    double  *cpu = NULL;
    double  *cpu_user = NULL;
    double  *cpu_sys = NULL;

    PDM_part_time_get (ppartId,
                       &elapsed,
                       &cpu,
                       &cpu_user,
                       &cpu_sys);

    if (i_rank == 0)
      PDM_printf("[%d] Temps dans ppart : %12.5e\n",
                 i_rank,  elapsed[0]);

    /* Statistiques */

    int    cells_average;
    int    cells_median;
    double cells_std_deviation;
    int    cells_min;
    int    cells_max;
    int    bound_part_faces_average;
    int    bound_part_faces_median;
    double bound_part_faces_std_deviation;
    int    bound_part_faces_min;
    int    bound_part_faces_max;
    int    bound_part_faces_sum;

    PDM_part_stat_get (ppartId,
                       &cells_average,
                       &cells_median,
                       &cells_std_deviation,
                       &cells_min,
                       &cells_max,
                       &bound_part_faces_average,
                       &bound_part_faces_median,
                       &bound_part_faces_std_deviation,
                       &bound_part_faces_min,
                       &bound_part_faces_max,
                       &bound_part_faces_sum);

    if (i_rank == 0) {
      PDM_printf ("Statistics :\n");
      PDM_printf ("  - Number of cells :\n");
      PDM_printf ("       * average            : %i\n", cells_average);
      PDM_printf ("       * median             : %i\n", cells_median);
      PDM_printf ("       * standard deviation : %12.5e\n", cells_std_deviation);
      PDM_printf ("       * min                : %i\n", cells_min);
      PDM_printf ("       * max                : %i\n", cells_max);
      PDM_printf ("  - Number of faces exchanging with another partition :\n");
      PDM_printf ("       * average            : %i\n", bound_part_faces_average);
      PDM_printf ("       * median             : %i\n", bound_part_faces_median);
      PDM_printf ("       * standard deviation : %12.5e\n", bound_part_faces_std_deviation);
      PDM_printf ("       * min                : %i\n", bound_part_faces_min);
      PDM_printf ("       * max                : %i\n", bound_part_faces_max);
      PDM_printf ("       * total              : %i\n", bound_part_faces_sum);
    }

    free (dVtxCoord);
    free (dFaceVtxIdx);
    free (dFaceVtx);
    free (dFaceEdge);
    free (dEdgeVtxIdx);
    free (dEdgeVtx);
    free (dEdgeFace);
    free (dEdgeGroupIdx);
    free (dEdgeGroup);

    _compute_faceVtx (ppartId,
                      n_part,
                      nFace,
                      faceVtxIdx,
                      faceVtx,
                      faceLNToGN,
                      nVtx,
                      vtxCoord,
                      vtxLNToGN);

    PDM_part_free (ppartId);

  }
  else {
    *nFace = (int *) malloc(sizeof(int) * n_part);
    *faceVtxIdx = (int **) malloc(sizeof(int *) * n_part);
    *faceVtx = (int **) malloc(sizeof(int *) * n_part);
    *faceLNToGN = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t *) * n_part);

    *nVtx = (int *) malloc(sizeof(int) * n_part);
    *vtxCoord = (double **) malloc(sizeof(double *) * n_part);
    *vtxLNToGN = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t *) * n_part);

    for (int ipart = 0; ipart < n_part; ipart++) {
      int _nFace = 0;
      int _sFaceEdge = 0;
      (*nFace)[ipart] = _nFace;
      (*faceVtxIdx) [ipart] = (int *) malloc(sizeof(int) * (_nFace + 1));
      (*faceVtxIdx) [ipart][0] = 0;
      (*faceVtx)[ipart] = (int *) malloc(sizeof(int) * _sFaceEdge);
      (*faceLNToGN)[ipart] = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * _nFace);

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
 * \brief  Create and split Mesh
 *
 * \param [in]      argc     Number of arguments
 * \param [in]      argv     Arguments
 *
 */

static void
_export_ini_mesh
(
 const PDM_MPI_Comm pdm_mpi_comm,
 const int      n_part,
 int            *nFace,
 int            **faceVtxIdx,
 int            **faceVtx,
 PDM_g_num_t    **faceLNToGN,
 int            *nVtx,
 double         **vtxCoord,
 PDM_g_num_t    **vtxLNToGN
)
{

  int i_rank;
  int numProcs;

  PDM_MPI_Comm_rank (PDM_MPI_COMM_WORLD, &i_rank);
  PDM_MPI_Comm_size (PDM_MPI_COMM_WORLD, &numProcs);

  /*
   *  Export Mesh to Ensight
   */

  int id_cs;

  id_cs = PDM_writer_create ("Ensight",
                             PDM_WRITER_FMT_ASCII,
                             PDM_WRITER_TOPO_CONSTANTE,
                             PDM_WRITER_OFF,
                             "test_2d_surf_ens",
                             "mesh1",
                             pdm_mpi_comm,
                             PDM_IO_ACCES_MPI_SIMPLE,
                             1.,
                             NULL);

  /*
   * Creation des variables
   */

  int id_var_num_part;
  int id_var_coo_x;
  int id_var_coo_xyz;
  int id_geom;

  id_var_num_part = PDM_writer_var_create (id_cs,
                                           PDM_WRITER_OFF,
                                           PDM_WRITER_VAR_SCALAIRE,
                                           PDM_WRITER_VAR_ELEMENTS,
                                           "num_part");

  id_var_coo_x = PDM_writer_var_create (id_cs,
                                        PDM_WRITER_ON,
                                        PDM_WRITER_VAR_SCALAIRE,
                                        PDM_WRITER_VAR_SOMMETS,
                                        "coo_x");

  id_var_coo_xyz = PDM_writer_var_create (id_cs,
                                          PDM_WRITER_ON,
                                          PDM_WRITER_VAR_VECTEUR,
                                          PDM_WRITER_VAR_SOMMETS,
                                          "coo_xyz");

    /*
     * Creation de la geometrie
     */

    char nom_geom[6];
    strcpy (nom_geom,"mesh1");

    id_geom = PDM_writer_geom_create (id_cs,
                                             nom_geom,
                                             PDM_WRITER_OFF,
                                             PDM_WRITER_OFF,
                                             n_part);
    /*
     * Debut des ecritures
     */

    int *nsom_part  = (int *) malloc(sizeof(int) * n_part);

    int *n_part_procs = (int *) malloc(sizeof(int) * numProcs);

    PDM_MPI_Allgather ((void *) &n_part,      1, PDM_MPI_INT,
                       (void *) n_part_procs, 1, PDM_MPI_INT,
                       PDM_MPI_COMM_WORLD);

    int *debPartProcs = (int *) malloc(sizeof(int) * (numProcs + 1));

    debPartProcs[0] = 0;
    for (int i = 0; i < numProcs; i++) {
      debPartProcs[i+1] = debPartProcs[i] + n_part_procs[i];
    }

    free(n_part_procs);

    PDM_writer_step_beg (id_cs, 0.);

    int **_face_nb =  malloc(sizeof(int *) * n_part);
    int **_face_idx =  malloc(sizeof(int *) * n_part);

    for (int ipart = 0; ipart < n_part; ipart++) {

      PDM_writer_geom_coord_set (id_cs,
                                 id_geom,
                                 ipart,
                                 nVtx[ipart],
                                 vtxCoord[ipart],
                                 vtxLNToGN[ipart]);

      _face_nb[ipart] = malloc(sizeof(int) * nFace[ipart]);
      _face_idx[ipart] = malloc(sizeof(int) * nFace[ipart]);

      for (int j = 0; j < nFace[ipart]; j++) {
        _face_nb[ipart][j]  = faceVtxIdx[ipart][j+1] - faceVtxIdx[ipart][j];
        _face_idx[ipart][j] = faceVtxIdx[ipart][j] + 1;
      }

      PDM_writer_geom_faces_facesom_add (id_cs,
                                         id_geom,
                                         ipart,
                                         nFace[ipart],
                                         _face_idx[ipart],
                                         _face_nb[ipart],
                                         faceVtx[ipart],
                                         faceLNToGN[ipart]);
    }

    for (int ipart = 0; ipart < n_part; ipart++) {
      free (_face_nb[ipart]);
      free (_face_idx[ipart]);
    }

    free(_face_nb);
    free(_face_idx);

    PDM_writer_geom_write(id_cs,
                          id_geom);

    /* Creation des variables :
       - numero de partition
       - scalaire
       - vecteur
       - tenseur
    */

    PDM_real_t **val_num_part = (PDM_real_t **) malloc (sizeof(PDM_real_t *) * n_part);
    PDM_real_t **val_coo_x    = (PDM_real_t **) malloc (sizeof(PDM_real_t *) * n_part);
    PDM_real_t **val_coo_xyz  = (PDM_real_t **) malloc (sizeof(PDM_real_t *) * n_part);

    for (int ipart = 0; ipart < n_part; ipart++) {

      val_num_part[ipart] = (PDM_real_t *) malloc (sizeof(PDM_real_t) * nFace[ipart]);
      val_coo_x[ipart]    = (PDM_real_t *) malloc (sizeof(PDM_real_t) * nVtx[ipart]);
      val_coo_xyz[ipart]  = (PDM_real_t *) malloc (sizeof(PDM_real_t) * 3 * nVtx[ipart]);
      nsom_part[ipart]    = nVtx[ipart];

      for (int i = 0; i < nFace[ipart]; i++) {
        val_num_part[ipart][i] = ipart + 1 + debPartProcs[i_rank];
      }

      for (int i = 0; i < nVtx[ipart]; i++) {
        val_coo_x[ipart][i]       = vtxCoord[ipart][3*i];
        val_coo_xyz[ipart][3*i  ] = vtxCoord[ipart][3*i  ];
        val_coo_xyz[ipart][3*i+1] = vtxCoord[ipart][3*i+1];
        val_coo_xyz[ipart][3*i+2] = vtxCoord[ipart][3*i+2];
      }

      PDM_writer_var_set (id_cs,
                          id_var_num_part,
                          id_geom,
                          ipart,
                          val_num_part[ipart]);

      PDM_writer_var_set (id_cs,
                          id_var_coo_x,
                          id_geom,
                          ipart,
                          val_coo_x[ipart]);

      PDM_writer_var_set (id_cs,
                          id_var_coo_xyz,
                          id_geom,
                          ipart,
                          val_coo_xyz[ipart]);

    }

    PDM_writer_var_write (id_cs,
                          id_var_num_part);

    PDM_writer_var_free (id_cs,
                         id_var_num_part);

    PDM_writer_var_write (id_cs,
                          id_var_coo_x);

    PDM_writer_var_free (id_cs,
                         id_var_coo_x);

    PDM_writer_var_write (id_cs,
                          id_var_coo_xyz);

    PDM_writer_var_free (id_cs,
                         id_var_coo_xyz);

    for (int ipart = 0; ipart < n_part; ipart++) {
      free (val_num_part[ipart]);
      free (val_coo_x[ipart]);
      free (val_coo_xyz[ipart]);
    }

    free (val_num_part);
    free (val_coo_x);
    free (val_coo_xyz);
    free (nsom_part);

    PDM_writer_step_end (id_cs);
    PDM_writer_geom_data_free (id_cs,
                               id_geom);

    PDM_writer_geom_free (id_cs,
                          id_geom);
    PDM_writer_free (id_cs);

    free (debPartProcs);

}


/**
 *
 * \brief  Main
 *
 */

int
main
(
int argc,
char *argv[]
)
{
  PDM_MPI_Init (&argc, &argv);

  /*
   *  Set default values
   */

  PDM_g_num_t      n_vtx_segA = 4;
  double           lengthA  = 1.;
  double           xminA = 0.;
  double           yminA = 0.;
  int              n_partA   = 1;

  int              post    = 0;
  PDM_part_split_t method  = PDM_PART_SPLIT_PARMETIS;
  int              haveRandom = 1;
  int              randomTimeInit = 0;

  int              randomMeshAInit = -1;

  int              nProcData = -1;

  int              i_rank;
  int              numProcs;

  /*
   *  Read args
   */

  _read_args (argc,
              argv,
              &n_vtx_segA,
              &lengthA,
              &xminA,
              &yminA,
              &n_partA,
              &post,
              (int *) &method,
              &haveRandom,
              &randomTimeInit,
              &randomMeshAInit,
              &nProcData
              );

  PDM_MPI_Comm_rank (PDM_MPI_COMM_WORLD, &i_rank);
  PDM_MPI_Comm_size (PDM_MPI_COMM_WORLD, &numProcs);

  if (i_rank == 0) {
    PDM_printf ("%Parametres : \n");
    PDM_printf ("  - n_vtx_segA : %d\n", n_vtx_segA);
    PDM_printf ("  - lengthA : %f\n", lengthA);
    PDM_printf ("  - xminA : %d\n", xminA);
    ;
    PDM_printf ("  - n_partA : %d\n", n_partA);
    PDM_printf ("  - post : %d\n", post);
    PDM_printf ("  - method : %d\n", method);
    PDM_printf ("  - haveRandom : %d\n", haveRandom);
    PDM_printf ("  - randomTimeInit : %d\n", randomTimeInit);
  }

  /*
   *  Create a partitioned mesh
   */

  PDM_g_num_t      nGFace;
  PDM_g_num_t      nGVtx;
  int            *nFace;
  int            **faceVtxIdx;
  int            **faceVtx;
  PDM_g_num_t    **faceLNToGN;
  int            *nVtx;
  double         **vtxCoord;
  PDM_g_num_t    **vtxLNToGN;

  int initRandom = 0;

  int n_part = n_partA;

  PDM_MPI_Comm meshComm = PDM_MPI_COMM_WORLD;
  int activeRankMesh = 1;

  if (nProcData > 0 && nProcData < numProcs) {
    int rankInNode = PDM_io_mpi_node_rank (PDM_MPI_COMM_WORLD);

    int nNode = 0;
    int iNode = -1;
    int masterRank = 0;
    if (rankInNode == 0) {
      masterRank = 1;
    }

    int *rankInNodes = malloc(sizeof(int) * numProcs);

    PDM_MPI_Allreduce (&masterRank, &nNode, 1, PDM_MPI_INT, PDM_MPI_SUM, PDM_MPI_COMM_WORLD);
    PDM_MPI_Allgather (&rankInNode, 1, PDM_MPI_INT, rankInNodes, 1, PDM_MPI_INT, PDM_MPI_COMM_WORLD);

    activeRankMesh = 0;

    for (int i = 0; i < i_rank; i++) {
      if (rankInNodes[i] == 0) {
        iNode += 1;
      }
    }

    if (nProcData <= nNode) {
      if (iNode < nProcData && rankInNode == 0) {
        activeRankMesh = 1;
      }
    }

    else {

      if (rankInNode < (nProcData / nNode)) {
        activeRankMesh = 1;
      }
      if ((rankInNode == (nProcData / nNode)) && (iNode < (nProcData % nNode))) {
        activeRankMesh = 1;
      }

    }

    PDM_MPI_Comm_split(PDM_MPI_COMM_WORLD, activeRankMesh, i_rank, &meshComm);

    free (rankInNodes);
  }

  double xmin;
  double ymin;
  double length;
  PDM_g_num_t n_vtx_seg;

  n_vtx_seg = n_vtx_segA;
  length = lengthA;
  xmin = xminA;
  ymin = yminA;
  n_part = n_partA;

  if (randomTimeInit) {
    initRandom =  time( NULL );
  }

  if (randomMeshAInit != -1) {
    initRandom = randomMeshAInit;
  }

  _create_split_mesh (activeRankMesh,
                      meshComm,
                      xmin,
                      ymin,
                      n_vtx_seg,
                      length,
                      n_part,
                      method,
                      haveRandom,
                      initRandom,
                      &nGFace,
                      &nGVtx,
                      &nFace,
                      &faceVtxIdx,
                      &faceVtx,
                      &faceLNToGN,
                      &nVtx,
                      &vtxCoord,
                      &vtxLNToGN);

  if (activeRankMesh) {
    _export_ini_mesh (meshComm,
                      n_part,
                      nFace,
                      faceVtxIdx,
                      faceVtx,
                      faceLNToGN,
                      nVtx,
                      vtxCoord,
                      vtxLNToGN);
  }

  if (nProcData > 0 && nProcData < numProcs) {
    PDM_MPI_Comm_free(&meshComm);
  }

  if (1 == 0) {
    for (int ipart = 0; ipart < n_part; ipart++) {
      printf("vtxCoord :\n");
      for (int j = 0; j < nVtx[ipart]; j++) {
        printf(""PDM_FMT_G_NUM" : %12.5e %12.5e %12.5e\n", vtxLNToGN[ipart][j],
               vtxCoord[ipart][3*j],
               vtxCoord[ipart][3*j+1],
               vtxCoord[ipart][3*j+2]);
      }
      printf("faceVtx :\n");
      for (int j = 0; j < nFace[ipart]; j++) {
        printf(""PDM_FMT_G_NUM" :", faceLNToGN[ipart][j]);
        for (int k = faceVtxIdx[ipart][j]; k < faceVtxIdx[ipart][j+1]; k++) {
          printf(" %d", faceVtx[ipart][k]);
        }
          printf("\n");
      }
    }
  }

  PDM_MPI_Finalize ();

  return 0;

}
