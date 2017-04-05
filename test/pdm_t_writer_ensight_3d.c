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
#include "pdm_part.h"
#include "pdm_dcube_gen.h"

#include "pdm_writer.h"

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
      *method = PDM_PART_SPLIT_PTSCOTCH;
    }
    else if (strcmp(argv[i], "-parmetis") == 0) {
      *method = PDM_PART_SPLIT_PARMETIS;
    }
    else if (strcmp(argv[i], "-hilbert") == 0) {
      *method = PDM_PART_SPLIT_HILBERT;
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

  PDM_g_num_t  nVtxSeg = 10;
  double        length  = 1.;
  int           nPart   = 1;
  int           post    = 0;
  PDM_part_split_t method  = PDM_PART_SPLIT_PARMETIS;

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

  PDM_dcube_gen_init(&id,
                      PDM_MPI_COMM_WORLD,
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
  int ppartId = 0;

  gettimeofday(&t_elaps_debut, NULL);

  /*
   *  Create mesh partitions
   */

  int have_dCellPart = 0;

  int *dCellPart = (int *) malloc(dNCell*sizeof(int));

  PDM_part_create(&ppartId,
                  PDM_MPI_COMM_WORLD,
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

  free(dCellPart);

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
#pragma warning(push)
#pragma warning(disable:2259)
  double t_elapsed = (double) tranche_elapsed_max/1000000.;
#pragma warning(pop)  

  printf("[%i]   - TEMPS DANS PART_CUBE  : %12.5e\n", myRank,  t_elapsed);

  int id_cs = PDM_writer_create("Ensight",
                                PDM_WRITER_FMT_ASCII,
                                PDM_writer_TOPO_CONSTANTE,
                                PDM_writer_OFF,
                                "test_3d_ens",
                                "chrd3d",
                                PDM_MPI_COMM_WORLD,
                                PDM_IO_ACCES_MPI_SIMPLE,
                                1.,
                                NULL);

  /* Creation de la geometrie */

  int id_geom = PDM_writer_geom_create(id_cs,
                             "test3d_geom",
                             PDM_writer_OFF,
                             PDM_writer_OFF,
                             nPart);

  /* Creation des variables */

  int id_var_num_part = PDM_writer_var_create(id_cs,
                                    PDM_writer_OFF,
                                    PDM_WRITER_VAR_SCALAIRE,
                                    PDM_WRITER_VAR_ELEMENTS,
                                    "num_part");

  int id_var_coo_x = PDM_writer_var_create(id_cs,
                                 PDM_writer_ON,
                                 PDM_WRITER_VAR_SCALAIRE,
                                 PDM_WRITER_VAR_SOMMETS,
                                 "coo_x");

  int id_var_coo_xyz = PDM_writer_var_create(id_cs,
                                   PDM_writer_ON,
                                   PDM_WRITER_VAR_VECTEUR,
                                   PDM_WRITER_VAR_SOMMETS,
                                   "coo_xyz");

  /* Debut d'ecritures */

  PDM_writer_step_beg(id_cs, 0.);

  int id_bloc = PDM_writer_geom_bloc_add(id_cs,
                                 id_geom,
                                 PDM_writer_OFF,
                                 PDM_writer_HEXA8);

  PDM_l_num_t **connec = (PDM_l_num_t **) malloc(sizeof(PDM_l_num_t *) * nPart);

  PDM_real_t **val_num_part = (PDM_real_t **) malloc(sizeof(PDM_real_t *) * nPart);
  PDM_real_t **val_coo_x    = (PDM_real_t **) malloc(sizeof(PDM_real_t *) * nPart);
  PDM_real_t **val_coo_xyz  = (PDM_real_t **) malloc(sizeof(PDM_real_t *) * nPart);
  int *nsom_part  = (int *) malloc(sizeof(int) * nPart);

  int *nPartProcs = (int *) malloc(sizeof(int) * numProcs);

  PDM_MPI_Allgather((void *) &nPart,      1, PDM_MPI_INT,
                (void *) nPartProcs, 1, PDM_MPI_INT,
                PDM_MPI_COMM_WORLD);

  int *debPartProcs = (int *) malloc(sizeof(int) * (numProcs + 1));

  debPartProcs[0] = 0;
  for (int i = 0; i < numProcs; i++) {
    debPartProcs[i+1] = debPartProcs[i] + nPartProcs[i];
  }

  free(nPartProcs);

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

    assert(sizeof(PDM_g_num_t) == sizeof(PDM_g_num_t));

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

    
    val_num_part[ipart] = (PDM_real_t *) malloc(sizeof(PDM_real_t) * nCell);
    val_coo_x[ipart]    = (PDM_real_t *) malloc(sizeof(PDM_real_t) * nVtx);
    val_coo_xyz[ipart]  = (PDM_real_t *) malloc(sizeof(PDM_real_t) * 3 * nVtx);
    nsom_part[ipart]    = nVtx;

    for (int i = 0; i < nCell; i++) {
      val_num_part[ipart][i] = ipart + 1 + debPartProcs[myRank];
    }

    for (int i = 0; i < nVtx; i++) {
      val_coo_x[ipart][i]       = vtx[3*i];
      val_coo_xyz[ipart][3*i  ] = vtx[3*i  ];
      val_coo_xyz[ipart][3*i+1] = vtx[3*i+1];
      val_coo_xyz[ipart][3*i+2] = vtx[3*i+2];
    }

    PDM_writer_geom_coord_set(id_cs,
                      id_geom,
                      ipart,
                      nVtx,
                      vtx,
                      vtxLNToGN);

    /* Construction de la connectivite pour sortie graphique */

    int *tagSom = (int *) malloc(sizeof(int) * nVtx);
    for (int i = 0; i < nVtx; i++) {
      tagSom[i] = 0;
    }

    connec[ipart] = (int *) malloc(sizeof(int) * nCell * 8);

    for (int i = 0; i < nCell; i++) {

      int *_connec = connec[ipart] + 8 * i;

      /* On tag les sommets de la premi�re face + calcul du centre et de la normale */

      int ifac1 = cellFace[cellFaceIdx[i]] - 1;
      double c_1[3];
      int ifac_oppose = -1;
      double c_oppose[3];
      double n_1[3];
      int ifac_cote = -1;
      double n_oppose[3];
      int isom1 = faceVtx[faceVtxIdx[ifac1]] - 1;
      int isom2 = faceVtx[faceVtxIdx[ifac1]+1] - 1;

      int idx = 0;
      c_1[0] = 0;
      c_1[1] = 0;
      c_1[2] = 0;
      for (int k = faceVtxIdx[ifac1]; k < faceVtxIdx[ifac1+1]; k++) {
        int isom = faceVtx[k] - 1;
        c_1[0] += vtx[3*isom];
        c_1[1] += vtx[3*isom+1];
        c_1[2] += vtx[3*isom+2];
        tagSom[isom] = 1;
        _connec[idx++] = faceVtx[k];
      }
      c_1[0] *= 0.25;
      c_1[1] *= 0.25;
      c_1[2] *= 0.25;

      n_1[0] = 0;
      n_1[1] = 0;
      n_1[2] = 0;
      for (int k = faceVtxIdx[ifac1]; k < faceVtxIdx[ifac1+1]; k++) {
        int isom = faceVtx[k] - 1;
        int suiv;
        if (k == faceVtxIdx[ifac1+1] -1)
          suiv = faceVtxIdx[ifac1];
        else
          suiv = k + 1;
        int isom_suiv = faceVtx[suiv] - 1;

        double v1[3];
        double v2[3];
        v1[0] = vtx[3*isom]   -  c_1[0];
        v1[1] = vtx[3*isom+1] -  c_1[1];
        v1[2] = vtx[3*isom+2] -  c_1[2];

        v2[0] = vtx[3*isom_suiv]   -  c_1[0];
        v2[1] = vtx[3*isom_suiv+1] -  c_1[1];
        v2[2] = vtx[3*isom_suiv+2] -  c_1[2];

        n_1[0] += v1[1] * v2[2] - v1[2] * v2[1];
        n_1[1] += v1[2] * v2[0] - v1[0] * v2[2];
        n_1[2] += v1[0] * v2[1] - v1[1] * v2[0];
      }

      /* Recherche de la face oppos�e et de la face en contact de la premi�re arrete */

      for (int k = cellFaceIdx[i] + 1; k < cellFaceIdx[i+1]; k++) {
        int ifac = cellFace[k] - 1;
        int cpt1 = 0;
        int cpt2 = 0;
        for (int k1 = faceVtxIdx[ifac]; k1 < faceVtxIdx[ifac+1]; k1++) {
          int isom = faceVtx[k1] - 1;
          if (tagSom[isom] == 0)
            cpt1++;
          if ((tagSom[isom] == 1) && ((isom == isom1) || (isom == isom2)))
            cpt2++;
        }
        if (cpt1 == 4)
          ifac_oppose = ifac;
        if (cpt2 == 2)
          ifac_cote = ifac;
      }

      /* Calcul du centre de la normale de la face opposee */

      c_oppose[0] = 0.;
      c_oppose[1] = 0.;
      c_oppose[2] = 0.;
      for (int k = faceVtxIdx[ifac_oppose]; k < faceVtxIdx[ifac_oppose+1]; k++) {
        int isom = faceVtx[k] - 1;
        c_oppose[0] += vtx[3*isom];
        c_oppose[1] += vtx[3*isom+1];
        c_oppose[2] += vtx[3*isom+2];
        _connec[idx++] = faceVtx[k];
      }
      c_oppose[0] *= 0.25;
      c_oppose[1] *= 0.25;
      c_oppose[2] *= 0.25;

      n_oppose[0] = 0;
      n_oppose[1] = 0;
      n_oppose[2] = 0;

      for (int k = faceVtxIdx[ifac_oppose]; k < faceVtxIdx[ifac_oppose+1]; k++) {
        int isom = faceVtx[k] - 1;
        int suiv;
        if (k == faceVtxIdx[ifac_oppose+1] -1)
          suiv = faceVtxIdx[ifac_oppose];
        else
          suiv = k + 1;
        int isom_suiv = faceVtx[suiv] - 1;

        double v1[3];
        double v2[3];
        v1[0] = vtx[3*isom]   -  c_1[0];
        v1[1] = vtx[3*isom+1] -  c_1[1];
        v1[2] = vtx[3*isom+2] -  c_1[2];

        v2[0] = vtx[3*isom_suiv]   -  c_1[0];
        v2[1] = vtx[3*isom_suiv+1] -  c_1[1];
        v2[2] = vtx[3*isom_suiv+2] -  c_1[2];

        n_oppose[0] += v1[1] * v2[2] - v1[2] * v2[1];
        n_oppose[1] += v1[2] * v2[0] - v1[0] * v2[2];
        n_oppose[2] += v1[0] * v2[1] - v1[1] * v2[0];
      }

      /* Inversion eventuelle du sens de rotation de la face opposee et (ou) de la face courante */

      double v_c1_c_oppose[3];

      v_c1_c_oppose[0] = c_oppose[0] - c_1[0];
      v_c1_c_oppose[1] = c_oppose[1] - c_1[1];
      v_c1_c_oppose[2] = c_oppose[2] - c_1[2];

      double pscal1 = v_c1_c_oppose[0] * n_1[0]
                    + v_c1_c_oppose[1] * n_1[1]
                    + v_c1_c_oppose[2] * n_1[2];

      double pscal2 = v_c1_c_oppose[0] * n_oppose[0]
                    + v_c1_c_oppose[1] * n_oppose[1]
                    + v_c1_c_oppose[2] * n_oppose[2];

      if (pscal1 < 0) {
        int tmp = _connec[0];
        _connec[0] = _connec[1];
        _connec[1] = tmp;
        tmp = _connec[2];
        _connec[2] = _connec[3];
        _connec[3] = tmp;
      }

      if (pscal2 < 0) {
        int tmp = _connec[4];
        _connec[4] = _connec[5];
        _connec[5] = tmp;
        tmp = _connec[6];
        _connec[6] = _connec[7];
        _connec[7] = tmp;
      }

      /* Permutation eventuelle de la connectivite de la face sup�rieure */

      int premier_som_opp = -1;
      for (int k = faceVtxIdx[ifac_cote]; k < faceVtxIdx[ifac_cote + 1]; k++) {
        int isom = faceVtx[k];
        int pre;
        int suiv;
        int isom_pre;
        int isom_suiv;
        if (k == faceVtxIdx[ifac_cote+1] -1)
          suiv = faceVtxIdx[ifac_cote];
        else
          suiv = k + 1;
        if (k == faceVtxIdx[ifac_cote])
          pre = faceVtxIdx[ifac_cote+1] - 1;
        else
          pre = k - 1;
        isom_suiv = faceVtx[suiv];
        isom_pre = faceVtx[pre];
        if (isom == _connec[0]) {
          if (isom_suiv == _connec[1])
            premier_som_opp = isom_pre;
          else if (isom_pre == _connec[1])
            premier_som_opp = isom_suiv;
          else {
            printf("Erreur face opposee\n");
            abort();
          }
        }
      }

      int *_connec_face_oppose = _connec + 4;
      int copy_connec[4];

      int id1 = 0;
      for (int k = 0; k < 4; k++) {
        copy_connec[k] = _connec_face_oppose[k];
        if (premier_som_opp == _connec_face_oppose[k])
          id1 = k;
      }

      for (int k = 0; k < 4; k++) {
        _connec_face_oppose[k] = copy_connec[(id1 + k) % 4];
      }
    }

    free(tagSom);

    PDM_writer_geom_bloc_std_set(id_cs,
                         id_geom,
                         id_bloc,
                         ipart,
                         nCell,
                         connec[ipart],
                         cellLNToGN);
  }

  free(debPartProcs);

  PDM_writer_geom_write(id_cs,
              id_geom);

  /* Creation des variables :
      - numero de partition
      - scalaire
      - vecteur
      - tenseur
   */

  for (int ipart = 0; ipart < nPart; ipart++) {

    PDM_writer_var_set(id_cs,
               id_var_num_part,
               id_geom,
               ipart,
               val_num_part[ipart]);
  }

  PDM_writer_var_write(id_cs,
             id_var_num_part);

  PDM_writer_var_free(id_cs,
             id_var_num_part);

  for (int ipart = 0; ipart < nPart; ipart++) {
    free(val_num_part[ipart]);
  }
  free(val_num_part);

  for (int nstep = 0; nstep < 10; nstep++) {

    double tstep = nstep * 0.01;

    if (nstep > 0)
      PDM_writer_step_beg(id_cs, tstep);

    for (int ipart = 0; ipart < nPart; ipart++) {

      for (int i = 0; i < nsom_part[ipart]; i++) {
        val_coo_x[ipart][i]       = val_coo_x[ipart][i]       + val_coo_x[ipart][i]/length;
        val_coo_xyz[ipart][3*i]   = val_coo_xyz[ipart][3*i]   + val_coo_xyz[ipart][3*i]/length;
        val_coo_xyz[ipart][3*i+1] = val_coo_xyz[ipart][3*i+1] + val_coo_xyz[ipart][3*i+1]/length;
        val_coo_xyz[ipart][3*i+2] = val_coo_xyz[ipart][3*i+2] + val_coo_xyz[ipart][3*i+2]/length;
      }

      PDM_writer_var_set(id_cs,
                 id_var_coo_x,
                 id_geom,
                 ipart,
                 val_coo_x[ipart]);

      PDM_writer_var_set(id_cs,
                 id_var_coo_xyz,
                 id_geom,
                 ipart,
                 val_coo_xyz[ipart]);

    }

    PDM_writer_var_write(id_cs,
               id_var_coo_x);

    PDM_writer_var_write(id_cs,
               id_var_coo_xyz);

    PDM_writer_var_data_free(id_cs,
                    id_var_coo_x);

    PDM_writer_var_data_free(id_cs,
                    id_var_coo_xyz);

    PDM_writer_step_end(id_cs);
  }

  for (int ipart = 0; ipart < nPart; ipart++) {
    free(val_coo_x[ipart]);
    free(val_coo_xyz[ipart]);
  }
  free(val_coo_x);
  free(val_coo_xyz);
  free(nsom_part);

  PDM_writer_var_free(id_cs,
             id_var_coo_x);

  PDM_writer_var_free(id_cs,
             id_var_coo_xyz);

  /* Liberation memoire */

  PDM_writer_geom_data_free(id_cs,
                   id_geom);

  PDM_writer_geom_free(id_cs,
              id_geom);

  PDM_writer_free(id_cs);

  for (int ipart = 0; ipart < nPart; ipart++) {
    free(connec[ipart]);
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

  free(connec);

  PDM_MPI_Finalize();

  return 0;
}

