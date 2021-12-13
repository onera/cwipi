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
#include "pdm_part_to_block.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_gnum.h"
#include "pdm_para_octree.h"
#include "pdm_distrib.h"
#include "pdm_logging.h"

/*============================================================================
 * Macro definitions
 *============================================================================*/

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
     "  -n      <level>  Number of points (default : 10).\n\n"
     "  -radius <level>  Radius of domain (default : 10).\n\n"
     "  -local           Number of points is local (default : global).\n\n"
     "  -rand            Random definition of point coordinates (default : false).\n\n"
     "  -h               This message.\n\n");

  exit (exit_code);
}

/**
 *
 * \brief  Read arguments from the command line
 *
 * \param [in]    argc   Number of arguments
 * \param [in]    argv   Arguments
 * \param [inout] nPts   Number of points
 * \param [inout] ls     Low scalability
 * \param [inout] length Length of domains
 *
 */

static void
_read_args
(
 int            argc,
 char         **argv,
 PDM_g_num_t   *gn_box,
 double        *length
)
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {

    if (strcmp(argv[i], "-h") == 0)
      _usage(EXIT_SUCCESS);

    else if (strcmp(argv[i], "-n") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        long _n = atol(argv[i]);
        *gn_box = (PDM_g_num_t) _n;
      }
    }

    else if (strcmp(argv[i], "-l") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *length = atof(argv[i]);
      }
    }

    else {
      _usage(EXIT_FAILURE);
    }
    i++;
  }
}



static void
_random_boxes
(
 PDM_MPI_Comm        comm,
 const PDM_g_num_t   gn_box,
 double              length,
 int                *n_box,
 PDM_g_num_t       **box_g_num,
 double            **box_extents
 )
{
  int i_rank, n_rank;
  PDM_MPI_Comm_rank (comm, &i_rank);
  PDM_MPI_Comm_size (comm, &n_rank);

  *n_box = (int) (gn_box / n_rank);
  if (i_rank < gn_box % n_rank) {
    (*n_box)++;
  }

  PDM_g_num_t* distrib_box = PDM_compute_entity_distribution(comm, (*n_box));
  for(int i = 0; i < 6 * distrib_box[i_rank]; ++i) {
    rand();
  }

  *box_g_num = malloc (sizeof(PDM_g_num_t) * (*n_box));
  for (int i = 0; i < *n_box; i++) {
    (*box_g_num)[i] = distrib_box[i_rank] + i + 1;
  }

  *box_extents = malloc (sizeof(double) * (*n_box) * 6);
  for (int i = 0; i < *n_box; i++) {
    for (int j = 0; j < 3; j++) {
      double x1 = length * (double) rand() / ((double) RAND_MAX);
      double x2 = length * (double) rand() / ((double) RAND_MAX);

      (*box_extents)[6*i + j    ] = PDM_MIN (x1, x2);
      (*box_extents)[6*i + j + 3] = PDM_MAX (x1, x2);
    }
  }

  free (distrib_box);
}

// static void
// _read_points
// (
//  PDM_MPI_Comm   comm,
//  int           *n_pts,
//  double       **pts_coord,
//  PDM_g_num_t  **pts_g_num
//  )
// {
//   int i_rank;
//   PDM_MPI_Comm_rank (comm, &i_rank);
//   int n_rank;
//   PDM_MPI_Comm_size (comm, &n_rank);

//   double dpts_coord00[3*16] = {-2.3642203211784363e-01, 2.9211959918029606e-04, -2.3000000044703484e-02,
//                                -2.3647038638591766e-01, 9.2450733063742518e-04, -2.3000000044703484e-02,
//                                -2.3643817007541656e-01, 1.5251128934323788e-03, -2.3000000044703484e-02,
//                                -2.3638600111007690e-01, 2.1253244485706091e-03, -2.3000000044703484e-02,
//                                -2.3647220432758331e-01, 2.1739276126027107e-03, -2.3000000044703484e-02,
//                                -2.3640592396259308e-01, 2.7762399986386299e-03, -2.3000000044703484e-02,
//                                -2.3646549880504608e-01, 2.8099860064685345e-03, -2.3000000044703484e-02,
//                                -2.3644705116748810e-01, 3.4485808573663235e-03, -2.3000000044703484e-02,
//                                -2.3645696043968201e-01, 3.4542013891041279e-03, -2.3000000044703484e-02,
//                                -2.3639747500419617e-01, 3.4202726092189550e-03, -2.3000000044703484e-02,
//                                -2.3644654452800751e-01, 4.1058603674173355e-03, -2.3000000044703484e-02,
//                                -2.3639376461505890e-01, 6.0901041142642498e-03, -2.3000000044703484e-02,
//                                -2.3638451099395752e-01, 6.7791189067065716e-03, -2.2335803136229515e-02,
//                                -2.3622378706932068e-01, 6.6702286712825298e-03, -2.3000000044703484e-02,
//                                -2.3637551069259644e-01, 6.7612598650157452e-03, -2.3000000044703484e-02,
//                                -2.3638534545898438e-01, 6.7670890130102634e-03, -2.3000000044703484e-02 };

// double dpts_coord01[3*16] = {-2.3609909415245056e-01, 6.5949372947216034e-03, -2.3000000044703484e-02,
//                              -2.3600970208644867e-01, 6.5408772788941860e-03, -2.3000000044703484e-02,
//                              -2.3648110032081604e-01, -1.5137524087549536e-06, -7.7877985313534737e-03,
//                              -2.3647101223468781e-01, 2.2220288883545436e-05, -7.6369130983948708e-03,
//                              -2.3648110032081604e-01, -7.6107226050226018e-06, -7.5133782811462879e-03,
//                              -2.3617900907993317e-01, 1.2028258061036468e-03, -8.5132550448179245e-03,
//                              -2.3647975921630859e-01, 7.9786370042711496e-04, -8.4376204758882523e-03,
//                              -2.3618096113204956e-01, 7.0157845038920641e-04, -8.4335263818502426e-03,
//                              -2.3625497519969940e-01, 6.0786749236285686e-04, -8.3844512701034546e-03,
//                              -2.3648095130920410e-01, 2.7402592240832746e-04, -8.0809118226170540e-03,
//                              -2.3648107051849365e-01, 1.2792067718692124e-04, -8.2131624221801758e-03,
//                              -2.3648110032081604e-01, 3.4934419090859592e-05, -8.0979196354746819e-03,
//                              -2.3648102581501007e-01, 1.9569303549360484e-04, -8.1492755562067032e-03,
//                              -2.3648102581501007e-01, 1.9482655625324696e-04, -7.9802740365266800e-03,
//                              -2.3647095263004303e-01, 1.9181826792191714e-04, -7.9799881204962730e-03,
//                              -2.3644144833087921e-01, 1.8296993221156299e-04, -7.9791769385337830e-03 };

// double dpts_coord02[3*16] = {-2.3648110032081604e-01, 4.1268708628194872e-06, -7.9412525519728661e-03,
//                              -2.3648107051849365e-01, -1.2178508040960878e-04, -8.0305952578783035e-03,
//                              -2.3648110032081604e-01, -4.6144199586706236e-05, -7.9739456996321678e-03,
//                              -2.3648110032081604e-01, -3.2485015253769234e-05, -8.1585561856627464e-03,
//                              -2.3648099601268768e-01, -2.2748758783563972e-04, -8.1226900219917297e-03,
//                              -2.3646648228168488e-01, 1.4693035045638680e-03, -8.5044112056493759e-03,
//                              -2.3647652566432953e-01, 1.4726326335221529e-03, -8.5043208673596382e-03,
//                              -2.3585462570190430e-01, 7.8195352107286453e-03, -2.3000000044703484e-02,
//                              -2.3599010705947876e-01, 7.2134872898459435e-03, -2.3000000044703484e-02,
//                              -2.3614916205406189e-01, 7.3106153868138790e-03, -2.3000000044703484e-02,
//                              -2.3624655604362488e-01, 7.3699280619621277e-03, -2.3000000044703484e-02,
//                              -2.3596845567226410e-01, 7.8897932544350624e-03, -2.3000000044703484e-02,
//                              -2.3605753481388092e-01, 7.9447636380791664e-03, -2.3000000044703484e-02,
//                              -2.3612722754478455e-01, 7.9877413809299469e-03, -2.3000000044703484e-02,
//                              -2.3622444272041321e-01, 8.0475555732846260e-03, -2.3000000044703484e-02,
//                              -2.3618176579475403e-01, 8.0213267356157303e-03, -2.3000000044703484e-02 };

// double dpts_coord03[3*16] = {-2.3625783622264862e-01, 8.0680558457970619e-03, -2.3000000044703484e-02,
//                              -2.3634275794029236e-01, 8.1198299303650856e-03, -2.3000000044703484e-02,
//                              -2.3625710606575012e-01, 8.0802170559763908e-03, -2.2333409637212753e-02,
//                              -2.3631957173347473e-01, 8.1184189766645432e-03, -2.2333499044179916e-02,
//                              -2.3634184896945953e-01, 8.1319389864802361e-03, -2.2333530709147453e-02,
//                              -2.0289508998394012e-01, -1.3912524096667767e-02, 1.5854747965931892e-02,
//                              -2.0288951694965363e-01, -1.3889461755752563e-02, 1.5893017873167992e-02,
//                              -2.0288845896720886e-01, -1.3929520733654499e-02, 1.5790639445185661e-02,
//                              -2.0288468897342682e-01, -1.3925082981586456e-02, 1.5812428668141365e-02,
//                              -2.0288324356079102e-01, -1.3929874636232853e-02, 1.5746030956506729e-02,
//                              -2.3564243316650391e-01, -1.7665555700659752e-02, 1.3941368088126183e-02,
//                              -2.3564209043979645e-01, -1.7669880762696266e-02, 1.3985215686261654e-02,
//                              -2.3564238846302032e-01, -1.7665766179561615e-02, 1.4029084704816341e-02,
//                              -2.3564276099205017e-01, -1.7660545185208321e-02, 1.4050510711967945e-02,
//                              -2.3558072745800018e-01, -1.7660630866885185e-02, 1.4040173962712288e-02,
//                              -2.3558051884174347e-01, -1.7663745209574699e-02, 1.4018338173627853e-02 };

// double dpts_coord04[3*16] = {-2.3558054864406586e-01, -1.7663657665252686e-02, 1.3974273577332497e-02,
//                              -2.3558121919631958e-01, -1.7655113711953163e-02, 1.3931053690612316e-02,
//                              -2.3558044433593750e-01, -1.7664752900600433e-02, 1.3996304012835026e-02,
//                              -2.3558178544044495e-01, -1.7647657543420792e-02, 1.3910303823649883e-02,
//                              -2.3564855754375458e-01, -1.7584335058927536e-02, 1.3811076059937477e-02,
//                              -2.3564559221267700e-01, -1.7623897641897202e-02, 1.3850141316652298e-02,
//                              -2.3575571179389954e-01, -1.7678266391158104e-02, 1.3942795805633068e-02,
//                              -2.3572884500026703e-01, -1.7651015892624855e-02, 1.4076048508286476e-02,
//                              -2.3573297262191772e-01, -1.7595518380403519e-02, 1.4143785461783409e-02,
//                              -2.3581948876380920e-01, -1.7677305266261101e-02, 1.4131541363894939e-02,
//                              -2.3573718965053558e-01, -1.7538871616125107e-02, 1.4177283272147179e-02,
//                              -2.3575313389301300e-01, -1.7323430627584457e-02, 1.4161890372633934e-02,
//                              -2.3574590682983398e-01, -1.7421359196305275e-02, 1.4188742265105247e-02,
//                              -2.3570168018341064e-01, -1.7514932900667191e-02, 1.4191069640219212e-02,
//                              -2.3559935390949249e-01, -1.7408972606062889e-02, 1.4214887283742428e-02,
//                              -2.3565654456615448e-01, -1.7474954947829247e-02, 1.4206533320248127e-02 };

// double dpts_coord05[3*17] = {-2.3570434749126434e-01, -1.7478996887803078e-02, 1.4197986572980881e-02,
//                              -2.3566104471683502e-01, -1.7414171248674393e-02, 1.4203875325620174e-02,
//                              -2.3583920300006866e-01, -1.7412312328815460e-02, 1.4274570159614086e-02,
//                              -2.3584011197090149e-01, -1.7399985343217850e-02, 1.4360441826283932e-02,
//                              -2.3583328723907471e-01, -1.7492171376943588e-02, 1.4361641369760036e-02,
//                              -2.3582439124584198e-01, -1.7611725255846977e-02, 1.4334590174257755e-02,
//                              -2.3582325875759125e-01, -1.7626896500587463e-02, 1.4398392289876938e-02,
//                              -2.3581369221210480e-01, -1.7754409462213516e-02, 1.4430612325668335e-02,
//                              -2.3581764101982117e-01, -1.7701987177133560e-02, 1.4460830949246883e-02,
//                              -2.3581647872924805e-01, -1.7717454582452774e-02, 1.4210636727511883e-02,
//                              -2.3581843078136444e-01, -1.7691381275653839e-02, 1.4234273694455624e-02,
//                              -2.3581404983997345e-01, -1.7749680206179619e-02, 1.4244493097066879e-02,
//                              -2.3582871258258820e-01, -1.7553821206092834e-02, 1.4303647913038731e-02,
//                              -2.3582541942596436e-01, -1.7597883939743042e-02, 1.4289562590420246e-02,
//                              -2.3582576215267181e-01, -1.7593365162611008e-02, 1.4200356788933277e-02,
//                              -2.3582221567630768e-01, -1.7640773206949234e-02, 1.4488190412521362e-02,
//                              -2.3581092059612274e-01, -1.7791194841265678e-02, 1.4293002896010876e-02 };

//   if (i_rank == 0) {
//     *n_pts = 16;
//   } else if (i_rank == 1) {
//     *n_pts = 16;
//   } else if (i_rank == 2) {
//     *n_pts = 16;
//   } else if (i_rank == 3) {
//     *n_pts = 16;
//   } else if (i_rank == 4) {
//     *n_pts = 16;
//   } else if (i_rank == 5) {
//     *n_pts = 17;
//   } else  {
//     abort();
//   }

//   if (n_rank == 1) {
//     *n_pts = 16 + 16 + 16 + 16 + 16 + 17;
//   }




//   *pts_g_num = malloc (sizeof(PDM_g_num_t) * (*n_pts));
//   *pts_coord = malloc (sizeof(double     ) * (*n_pts) * 3);

//   if (n_rank == 1) {

//     memcpy(*pts_coord, dpts_coord00, sizeof(double) * 16 * 3);
//     memcpy(*pts_coord + 3*16, dpts_coord01, sizeof(double) * 16 * 3);
//     memcpy(*pts_coord + 3*(16 + 16), dpts_coord02, sizeof(double) * 16 * 3);
//     memcpy(*pts_coord + 3*(16 + 16 + 16), dpts_coord03, sizeof(double) * 16 * 3);
//     memcpy(*pts_coord + 3*(16 + 16 + 16 + 16), dpts_coord04, sizeof(double) * 16 * 3);
//     memcpy(*pts_coord + 3*(16 + 16 + 16 + 16 + 16), dpts_coord05, sizeof(double) * 17 * 3);

//   } else {

//     if (i_rank == 0) {
//       memcpy(*pts_coord, dpts_coord00, sizeof(double) * (*n_pts) * 3);
//     } else if (i_rank == 1) {
//       memcpy(*pts_coord, dpts_coord01, sizeof(double) * (*n_pts) * 3);
//     } else if (i_rank == 2) {
//       memcpy(*pts_coord, dpts_coord02, sizeof(double) * (*n_pts) * 3);
//     } else if (i_rank == 3) {
//       memcpy(*pts_coord, dpts_coord03, sizeof(double) * (*n_pts) * 3);
//     } else if (i_rank == 4) {
//       memcpy(*pts_coord, dpts_coord04, sizeof(double) * (*n_pts) * 3);
//     } else if (i_rank == 5) {
//       memcpy(*pts_coord, dpts_coord05, sizeof(double) * (*n_pts) * 3);
//     } else  {
//       abort();
//     }
//   }


//   PDM_g_num_t *distrib_pts = PDM_compute_entity_distribution(comm, (*n_pts));
//   for (int i = 0; i < *n_pts; i++) {
//     (*pts_g_num)[i] = distrib_pts[i_rank] + i + 1;
//   }
//   free (distrib_pts);
// }





static void
_read_points2
(
 PDM_MPI_Comm   comm,
 int           *n_pts,
 double       **pts_coord,
 PDM_g_num_t  **pts_g_num
 )
{
  int i_rank;
  PDM_MPI_Comm_rank (comm, &i_rank);
  int n_rank;
  PDM_MPI_Comm_size (comm, &n_rank);

  //assert(n_rank == 1);

  double coord[3*13] = {9.8115103639715859e-02, 5.0809236975297911e-01, 2.2837191775365356e-01,
                        2.0681046941169096e-01, 6.3192983972929873e-01, 3.6904037906743603e-01,
                        4.9462113149679321e-01, 9.1781326865209889e-01, 2.2456504077397521e-01,
                        7.7822140451903521e-01, 6.2575317692745158e-01, 6.4253503828893177e-01,
                        6.2835879749542045e-01, 7.8624259515490502e-01, 6.3375764697033798e-01,
                        7.8113486516807917e-01, 5.1134893888670441e-01, 9.2954303022452767e-01,
                        6.2741732725753330e-01, 6.2917677264622263e-01, 7.6557964275385237e-01,
                        5.7079123334995996e-01, 4.3588869922137291e-01, 9.7577328482911607e-01,
                        9.7776655109960886e-01, 2.9817392563827982e-01, 7.0606834709461230e-01,
                        5.0912678759504426e-01, 9.3618765055955733e-02, 7.7025472374644810e-01,
                        5.0207314246896328e-01, 2.4486308083630312e-01, 9.2353304398876290e-01,
                        3.6273379579779408e-01, 2.2136588524624978e-01, 6.4526300303417394e-01,
                        3.5429132562330523e-01, 3.6711760670743776e-01, 7.8287646871659744e-01};


  *n_pts = 13 / n_rank;
  if (i_rank < 13 % n_rank) {
    (*n_pts)++;
  }

  PDM_g_num_t *distrib_pts = PDM_compute_entity_distribution(comm, (*n_pts));
  PDM_log_trace_array_long (distrib_pts, n_rank+1, "distrib_pts : ");

  *pts_g_num = malloc (sizeof(PDM_g_num_t) * (*n_pts));
  *pts_coord = malloc (sizeof(double     ) * (*n_pts) * 3);

  memcpy(*pts_coord, coord + 3*distrib_pts[i_rank], sizeof(double) * (*n_pts) * 3);

  for (int i = 0; i < *n_pts; i++) {
    (*pts_g_num)[i] = distrib_pts[i_rank] + i + 1;
  }
  free (distrib_pts);
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief  Main
 *
 */

int
main (int argc, char *argv[])
{
  srand(0);

  PDM_MPI_Init (&argc, &argv);

  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;

  int i_rank, n_rank;
  PDM_MPI_Comm_rank (comm, &i_rank);
  PDM_MPI_Comm_size (comm, &n_rank);



  PDM_g_num_t gn_box = 10;
  double      length = 3.14;

  _read_args (argc,
              argv,
              &gn_box,
              &length);



  /*
   *  Random boxes
   */
  int n_box;
  PDM_g_num_t *box_g_num   = NULL;
  double      *box_extents = NULL;
  _random_boxes (comm,
                 gn_box,
                 length,
                 &n_box,
                 &box_g_num,
                 &box_extents);

  if (1) {
    for (int i = 0; i < n_box; i++) {
      log_trace("box "PDM_FMT_G_NUM" extents = %f %f %f  %f %f %f\n",
                box_g_num[i],
                box_extents[6*i + 0],
                box_extents[6*i + 1],
                box_extents[6*i + 2],
                box_extents[6*i + 3],
                box_extents[6*i + 4],
                box_extents[6*i + 5]);
    }
  }


  /*
   *  Random point
   */
  int n_pts = 0;
  double      *pts_coord = NULL;
  PDM_g_num_t *pts_g_num = NULL;
  if (0) {
    if (i_rank == 0) {
      n_pts = 1;
      pts_coord = malloc (sizeof(double) * n_pts * 3);
      pts_g_num = malloc (sizeof(PDM_g_num_t) * n_pts);

      for (int i = 0; i < 3; i++) {
        pts_coord[i] = 0.5 * length;
      }

      pts_g_num[0] = 1;
    }


  }

  else {
    _read_points2 (comm,
                  &n_pts,
                  &pts_coord,
                  &pts_g_num);
  }


  const int octree_depth_max = 31;
  const int octree_points_in_leaf_max = 1;
  const int octree_build_leaf_neighbours = 0;

  /* Create octree structure */
  PDM_para_octree_t *octree = PDM_para_octree_create (1,
                                                      octree_depth_max,
                                                      octree_points_in_leaf_max,
                                                      octree_build_leaf_neighbours,
                                                      comm);

  /* Set octree point cloud */
  PDM_para_octree_point_cloud_set (octree,
                                   0,
                                   n_pts,
                                   pts_coord,
                                   pts_g_num);

  /* Build parallel octree */
  PDM_para_octree_build (octree, NULL);



  log_trace("PDM_para_octree_build OK\n");
#if 0
  /* Find points inside boxes */
  int         *box_pts_idx   = NULL;
  PDM_g_num_t *box_pts_g_num = NULL;
  double      *box_pts_coord = NULL;
  PDM_para_octree_points_inside_boxes_with_copies (octree,
                                                   n_box,
                                                   box_extents,
                                                   box_g_num,
                                                   &box_pts_idx,
                                                   &box_pts_g_num,
                                                   &box_pts_coord);

  /* Free octree */
  PDM_para_octree_free (octree);

  // PDM_log_trace_connectivity_long (box_pts_idx,
  //                                  box_pts_g_num,
  //                                  n_box,
  //                                  "box_pts_g_num :");
  free (box_pts_idx);
  free (box_pts_g_num);
  free (box_pts_coord);

#endif



  free (box_extents);
  free (box_g_num);
  if (i_rank == 0) {
    free (pts_coord);
    free (pts_g_num);
  }

  PDM_MPI_Finalize();

  return 0;
}
