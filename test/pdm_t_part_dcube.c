#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "pdm.h"
#include "pdm_mpi.h"
#include "pdm_part.h"
#include "pdm_dcube_gen.h"

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
  printf
    ("\n"
     "  Usage: \n\n"
     "  -n      <level>  Number of vertices on the cube side.\n\n"
     "  -l      <level>  Cube length.\n\n"
     "  -n_part <level>  Number of partitions par process.\n\n"
     "  -post            Ensight outputs (only if n_part == 1). \n\n"
     "  -parmetis        Call ParMETIS.\n\n"
     "  -pt-scocth       Call PT-Scotch.\n\n"
     "  -h               This message.\n\n");

  exit(exit_code);
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
_read_args(int            argc,
           char         **argv,
           PDM_g_num_t  *nVtxSeg,
           double        *length,
           int           *n_part,
	   int           *post,
	   int           *method)
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
        long _nVtxSeg = atol(argv[i]);
        *nVtxSeg = (PDM_g_num_t) _nVtxSeg;
      }
    }
    else if (strcmp(argv[i], "-l") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *length = atof(argv[i]);
    }
    else if (strcmp(argv[i], "-n_part") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *n_part = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-post") == 0) {
      *post = 1;
    }
    else if (strcmp(argv[i], "-pt-scotch") == 0) {
      *method = 2;
    }
    else if (strcmp(argv[i], "-parmetis") == 0) {
      *method = 1;
    }
    else
      _usage(EXIT_FAILURE);
    i++;
  }
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

  PDM_g_num_t        nVtxSeg = 10;
  double             length  = 1.;
  int                nPart   = 1;
  int                post    = 0;
  PDM_part_split_t   method  = PDM_PART_SPLIT_PTSCOTCH;

  /*
   *  Read args
   */

  _read_args(argc,
             argv,
             &nVtxSeg,
             &length,
             &nPart,
             &post,
             (int *) &method);

  /*
   *  Init
   */

  struct timeval t_elaps_debut;

  int myRank;
  int numProcs;

  PDM_MPI_Init(&argc, &argv);
  PDM_MPI_Comm_rank(PDM_MPI_COMM_WORLD, &myRank);
  PDM_MPI_Comm_size(PDM_MPI_COMM_WORLD, &numProcs);

  int           dNCell;
  int           dNFace;
  int           dNVtx;
  int           nFaceGroup;
  PDM_g_num_t *dFaceCell = NULL;
  int          *dFaceVtxIdx = NULL;
  PDM_g_num_t *dFaceVtx = NULL;
  double       *dVtxCoord = NULL;
  int          *dFaceGroupIdx = NULL;
  PDM_g_num_t *dFaceGroup = NULL;
  int           dFaceVtxL;
  int           dFaceGroupL;

  /*
   *  Create distributed cube
   */

  int          id;
  PDM_MPI_Comm     comm = PDM_MPI_COMM_WORLD;

  PDM_dcube_gen_init(&id,
                      comm,
                      nVtxSeg,
                      length);

  PDM_dcube_gen_dim_get(id,
                         &nFaceGroup,
                         &dNCell,
                         &dNFace,
                         &dNVtx,
                         &dFaceVtxL,
                         &dFaceGroupL);

  PDM_dcube_gen_data_get(id,
                          &dFaceCell,
                          &dFaceVtxIdx,
                          &dFaceVtx,
                          &dVtxCoord,
                          &dFaceGroupIdx,
                          &dFaceGroup);

  if (0 == 1) {

    printf("[%i] nFaceGroup    : %i\n", myRank, nFaceGroup);
    printf("[%i] dNCell        : %i\n", myRank, dNCell);
    printf("[%i] dNFace        : %i\n", myRank, dNFace);
    printf("[%i] dNVtx         : %i\n", myRank, dNVtx);

    printf("[%i] dFaceCell     : ", myRank);
    for (int i = 0; i < 2 * dNFace; i++)
      printf(" "PDM_FMT_G_NUM, dFaceCell[i]);
    printf("\n");

    printf("[%i] dFaceVtxIdx   : ", myRank);
    for (int i = 0; i < dNFace + 1; i++)
      printf(" %i", dFaceVtxIdx[i]);
    printf("\n");

    printf("[%i] dFaceVtx      : ", myRank);
    for (int i = 0; i < dFaceVtxIdx[dNFace]; i++)
      printf(" "PDM_FMT_G_NUM, dFaceVtx[i]);
    printf("\n");

    printf("[%i] dVtxCoord     : ", myRank);
    for (int i = 0; i < 3*dNVtx; i++)
      printf(" %12.5e", dVtxCoord[i]);
    printf("\n");

    printf("[%i] dFaceGroupIdx : ", myRank);
    for (int i = 0; i < nFaceGroup + 1; i++)
      printf(" %i", dFaceGroupIdx[i]);
    printf("\n");

    printf("[%i] dFaceGroup    : ", myRank);
    for (int i = 0; i < dFaceGroupIdx[nFaceGroup]; i++)
      printf(" "PDM_FMT_G_NUM, dFaceGroup[i]);
    printf("\n");

  }
  int ppartId = 0;

  gettimeofday(&t_elaps_debut, NULL);

  /*
   *  Create mesh partitions
   */

  int have_dCellPart = 0;

  int *dCellPart = (int *) malloc(dNCell*sizeof(int));

  PDM_part_create(&ppartId,
                  comm,
                  method,
                  PDM_PART_RENUM_CELL_NONE,
                  PDM_PART_RENUM_FACE_NONE,
                  nPart,
                  dNCell,
                  dNFace,
                  dNVtx,
                  nFaceGroup,
                  NULL,
                  NULL,
                  NULL,
                  NULL,
                  have_dCellPart,
                  dCellPart,
                  dFaceCell,
                  dFaceVtxIdx,
                  dFaceVtx,
                  NULL,
                  dVtxCoord,
                  NULL,
                  dFaceGroupIdx,
                  dFaceGroup);

  double  *elapsed = NULL;
  double  *cpu = NULL;
  double  *cpu_user = NULL;
  double  *cpu_sys = NULL;

  PDM_part_time_get(ppartId,
                 &elapsed,
                 &cpu,
                 &cpu_user,
                 &cpu_sys);

  printf("[%i]   - elapsed total                    : %12.5e\n", myRank, elapsed[0]);
  printf("[%i]   - elapsed building graph           : %12.5e\n", myRank, elapsed[1]);
  printf("[%i]   - elapsed splitting graph          : %12.5e\n", myRank, elapsed[2]);
  printf("[%i]   - elapsed building mesh partitions : %12.5e\n", myRank, elapsed[3]);

  printf("[%i]   - cpu total                        : %12.5e\n", myRank, cpu[0]);
  printf("[%i]   - cpu building graph               : %12.5e\n", myRank, cpu[1]);
  printf("[%i]   - cpu splitting graph              : %12.5e\n", myRank, cpu[2]);
  printf("[%i]   - cpu building mesh partitions     : %12.5e\n", myRank, cpu[3]);

  printf("[%i]   - cpu_user total                   : %12.5e\n", myRank, cpu_user[0]);
  printf("[%i]   - cpu_user building graph          : %12.5e\n", myRank, cpu_user[1]);
  printf("[%i]   - cpu_user splitting graph         : %12.5e\n", myRank, cpu_user[2]);
  printf("[%i]   - cpu_user building mesh partitions: %12.5e\n", myRank, cpu_user[3]);

  printf("[%i]   - cpu_sys total                    : %12.5e\n", myRank, cpu_sys[0]);
  printf("[%i]   - cpu_sys building graph           : %12.5e\n", myRank, cpu_sys[1]);
  printf("[%i]   - cpu_sys splitting graph          : %12.5e\n", myRank, cpu_sys[2]);
  printf("[%i]   - cpu_sys building mesh partitions : %12.5e\n", myRank, cpu_sys[3]);

  struct timeval t_elaps_fin;
  gettimeofday(&t_elaps_fin, NULL);

  long tranche_elapsed = (t_elaps_fin.tv_usec + 1000000 * t_elaps_fin.tv_sec) -
                         (t_elaps_debut.tv_usec + 1000000 *
                          t_elaps_debut.tv_sec);
  long tranche_elapsed_max = tranche_elapsed;
#if defined(__INTEL_COMPILER)
#pragma warning(push)
#pragma warning(disable:2259)
#endif
  double t_elapsed = (double) tranche_elapsed_max/1000000.;
#if defined(__INTEL_COMPILER)
#pragma warning(pop)
#endif

  printf("[%i]   - TEMPS DANS PART_CUBE  : %12.5e\n", myRank,  t_elapsed);

  if (0 == 1) {
    for (int ipart = 0; ipart < nPart; ipart++) {

      int nCell;
      int nFace;
      int nFacePartBound;
      int nVtx;
      int nProc;
      int nTPart;
      int sCellFace;
      int sFaceVtx;
      int sFaceGroup;

      PDM_part_part_dim_get(ppartId,
                         ipart,
                         &nCell,
                         &nFace,
                         &nFacePartBound,
                         &nVtx,
                         &nProc,
                         &nTPart,
                         &sCellFace,
                         &sFaceVtx,
                         &sFaceGroup);

      int          *cellTag;
      int          *cellFaceIdx;
      int          *cellFace;
      PDM_g_num_t *cellLNToGN;
      int          *faceTag;
      int          *faceCell;
      int          *faceVtxIdx;
      int          *faceVtx;
      PDM_g_num_t *faceLNToGN;
      int          *facePartBoundProcIdx;
      int          *facePartBoundPartIdx;
      int          *facePartBound;
      int          *vtxTag;
      double       *vtx;
      PDM_g_num_t *vtxLNToGN;
      int          *faceGroupIdx;
      int          *faceGroup;
      PDM_g_num_t *faceGroupLNToGN;

      PDM_part_part_val_get(ppartId,
                         ipart,
                         &cellTag,
                         &cellFaceIdx,
                         &cellFace,
                         &cellLNToGN,
                         &faceTag,
                         &faceCell,
                         &faceVtxIdx,
                         &faceVtx,
                         &faceLNToGN,
                         &facePartBoundProcIdx,
                         &facePartBoundPartIdx,
                         &facePartBound,
                         &vtxTag,
                         &vtx,
                         &vtxLNToGN,
                         &faceGroupIdx,
                         &faceGroup,
                         &faceGroupLNToGN);


      printf("[%i] nFaceGroup     : %i\n", myRank, nFaceGroup);
      printf("[%i] nCell          : %i\n", myRank, nCell);
      printf("[%i] nFace          : %i\n", myRank, nFace);
      printf("[%i] nVtx           : %i\n", myRank, nVtx);
      printf("[%i] nFacePartBound : %i\n", myRank, nFacePartBound);

      printf("[%i] cellFace     : ", myRank);
      for (int i = 0; i < nCell; i++) {
        for (int j = cellFaceIdx[i]; j < cellFaceIdx[i+1]; j++) {
          printf(" %i", cellFace[j]);
        }
        printf("\n");
      }

      printf("\n");

      printf("[%i]  cellLNToGN    : ", myRank);
      for (int i = 0; i < nCell; i++)
        printf(" "PDM_FMT_G_NUM, cellLNToGN[i]);
      printf("\n");

      printf("[%i] faceCell     : ", myRank);
      for (int i = 0; i < 2 * nFace; i++)
        printf(" %i", faceCell[i]);
      printf("\n");

      printf("[%i] faceVtx      : ", myRank);
      for (int i = 0; i < nFace; i++) {
        for (int j = faceVtxIdx[i]; j < faceVtxIdx[i+1]; j++) {
          printf(" %i", faceVtx[j]);
        }
        printf("\n");
      }

      printf("[%i]  faceLNToGN    : ", myRank);
      for (int i = 0; i < nFace; i++)
        printf(" "PDM_FMT_G_NUM, faceLNToGN[i]);
      printf("\n");

      printf("[%i] vtx           : ", myRank);
      for (int i = 0; i < 3 * nVtx; i++)
        printf(" %12.5e", vtx[i]);
      printf("\n");

      printf("[%i] vtxLNToGN     : ", myRank);
      for (int i = 0; i <  nVtx; i++)
        printf(" "PDM_FMT_G_NUM, vtxLNToGN[i]);
      printf("\n");

      printf("[%i] faceGroupIdx : ", myRank);
      for (int i = 0; i < nFaceGroup + 1; i++)
        printf(" %i", faceGroupIdx[i]);
      printf("\n");

      printf("[%i] faceGroup    : ", myRank);
      for (int i = 0; i < nFaceGroup; i++) {
        for (int j = faceGroupIdx[i]; j < faceGroupIdx[i+1]; j++) {
          printf(" %i", faceGroup[j]);
        }
        printf("\n");
      }

      printf("[%i] faceGroupLNToGN   : ", myRank);
      for (int i = 0; i < nFaceGroup; i++) {
        for (int j = faceGroupIdx[i]; j < faceGroupIdx[i+1]; j++) {
          printf(" "PDM_FMT_G_NUM, faceGroupLNToGN[j]);
        }
        printf("\n");
      }
    }
  }

  /* Calculs statistiques */

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

  PDM_part_stat_get(ppartId,
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

  if (myRank == 0) {
    printf("Statistics :\n");
    printf("  - Number of cells :\n");
    printf("       * average            : %i\n", cells_average);
    printf("       * median             : %i\n", cells_median);
    printf("       * standard deviation : %12.5e\n", cells_std_deviation);
    printf("       * min                : %i\n", cells_min);
    printf("       * max                : %i\n", cells_max);
    printf("  - Number of faces exchanging with another partition :\n");
    printf("       * average            : %i\n", bound_part_faces_average);
    printf("       * median             : %i\n", bound_part_faces_median);
    printf("       * standard deviation : %12.5e\n", bound_part_faces_std_deviation);
    printf("       * min                : %i\n", bound_part_faces_min);
    printf("       * max                : %i\n", bound_part_faces_max);
    printf("       * total              : %i\n", bound_part_faces_sum);
  }

  PDM_part_free(ppartId);

  PDM_dcube_gen_free(id);


  PDM_MPI_Finalize();

  return 0;
}
