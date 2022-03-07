/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <assert.h>
#include <math.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_part_to_block.h"
#include "pdm_part_to_block_priv.h"
#include "pdm_mpi_node_first_rank.h"
#include "pdm.h"
#include "pdm_timer.h"
#include "pdm_priv.h"
#include "pdm_binary_search.h"
#include "pdm_sort.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_array.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type
 *============================================================================*/

/*=============================================================================
 * Static global variables
 *============================================================================*/

static const double  pdm_part_to_block_distrib_tol = 0.10;

/* Max. number of sub-iterations to get a well-balanced distribution */
static const int pdm_part_to_block_distrib_n_iter_max = 5;

static const int _sampling_factors[4] = {1, /* OD */
                                         2, /* 1D */
                                         2, /* 2D */
                                         4, /* 3D */};

/*=============================================================================
 * Static global variables
 *============================================================================*/

static  double t_elaps[3] = {0., 0., 0.};
static  double t_cpu[3] = {0., 0., 0.};
static  PDM_timer_t *t_timer[3] = {NULL, NULL, NULL};

static  int min_exch_rank[2] = {INT_MAX, INT_MAX};
static  int max_exch_rank[2] = {-1, -1};

static  unsigned long long exch_data[2] = {0, 0};

int n_ptb = 0;


/*=============================================================================
 * Static function definitions
 *============================================================================*/

static
void
_counting_sort_long
(
 PDM_part_to_block_t *ptb
)
{
  // On a un pb si la distrib est géant
  int n_elt_block_tot = ptb->data_distrib_index[ptb->i_rank+1] - ptb->data_distrib_index[ptb->i_rank];
  int* block_n = (int *) malloc(n_elt_block_tot * sizeof(int));

  PDM_g_num_t* block_gnum = malloc (sizeof(PDM_g_num_t) * ptb->tn_recv_data);

  for(int i = 0; i < n_elt_block_tot; ++i) {
    block_n[i] = 0;
  }

  for(int i = 0; i < ptb->tn_recv_data; ++i) {
    int lid = ptb->sorted_recv_gnum[i] - ptb->data_distrib_index[ptb->i_rank] - 1;
    block_n[lid] += 1;
    block_gnum[i] = ptb->sorted_recv_gnum[i];
  }

  int* block_idx = (int *) malloc( (n_elt_block_tot+1) * sizeof(int));
  block_idx[0] = 0;
  for(int i = 0; i < n_elt_block_tot; ++i) {
    block_idx[i+1] = block_idx[i] + block_n[i];
    block_n[i] = 0;
  }

  for(int i = 0; i < ptb->tn_recv_data; ++i) {
    int lid = block_gnum[i] - ptb->data_distrib_index[ptb->i_rank] - 1;
    int idx_write = block_idx[lid] + block_n[lid]++;
    ptb->sorted_recv_gnum[idx_write] = block_gnum[i];
    ptb->order[idx_write] = i;
  }

  free(block_n);
  free(block_idx);
  free(block_gnum);

}

/**
 * \brief   Evaluate a distribution array.
 *
 * \param [in]  n_ranges     Number of ranges in the distribution
 * \param [in]  distribution Number of elements associated to each range of the distribution
 * \param [in]  optim        Optimal count in each range
 *
 * \return  a fit associated to the distribution. If fit = 0, distribution is perfect.
 *
 */

static double
_evaluate_distribution(int          n_ranges,
                       double      *distribution,
                       double       optim)
{
  int  i;
  double  d_low = 0, d_up = 0, fit = 0;

  /*
     d_low is the max gap between the distribution count and the optimum when
     distribution is lower than optimum.
     d_up is the max gap between the distribution count and the optimum when
     distribution is greater than optimum.
  */

  for (i = 0; i < n_ranges; i++) {

    if (distribution[i] > optim)
      d_up = PDM_MAX(d_up, distribution[i] - optim);
    else
      d_low = PDM_MAX(d_low, optim - distribution[i]);

  }

  fit = (d_up + d_low) / optim;

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  if (cs_glob_rank_id <= 0)
    PDM_printf( "<DISTRIBUTION EVALUATION> optim: %g, fit: %g\n",
               optim, fit);
#endif

  return  fit;
}


/**
 * \brief Define a global distribution associated to a sampling array i.e. count
 * the number of elements in each range.
 *
 *   \param [in]    dim           2D or 3D
 *   \param [in]    n_ranks       number of ranks (= number of ranges)
 *   \param [in]    gsum_weight   global sum of all weightings
 *   \param [in]    n_codes       local number of Hilbert codes
 *   \param [in]    hilbert_codes local list of Hilbert codes to distribute
 *   \param [in]    weight        weighting related to each code
 *   \param [in]    order         ordering array
 *   \param [in]    sampling      sampling array
 *   \param [inout] c_freq        pointer to the cumulative frequency array
 *   \param [inout] g_distrib     pointer to a distribution array
 *   \param [in]    comm          mpi communicator
 */

static void
_define_rank_distrib( PDM_part_to_block_t *ptb,
                      int                 dim,
                      int                 n_ranks,
                      double              gsum_weight,
                      const PDM_g_num_t   sampling[],
                      double              cfreq[],
                      double              g_distrib[],
                      PDM_MPI_Comm        comm)
{
  int  id, rank_id;

  const int  sampling_factor = _sampling_factors[dim];
  const int  n_samples = sampling_factor * n_ranks;

  /* Initialization */

  double   *l_distrib = (double  *) malloc (n_samples * sizeof(double));

  for (id = 0; id < n_samples; id++) {
    l_distrib[id] = 0;
    g_distrib[id] = 0;
  }

  for (int i = 0; i < ptb->n_part; i++) {

    for (int j = 0; j < ptb->n_elt[i]; j++) {

      PDM_g_num_t _gnum_elt = PDM_ABS(ptb->gnum_elt[i][j]) - 1;

      int iSample = PDM_binary_search_gap_long (_gnum_elt,
                                                sampling,
                                                n_samples + 1);
      l_distrib[iSample] += ptb->weight[i][j];
    }
  }

  /* Define the global distribution */

  PDM_MPI_Allreduce(l_distrib, g_distrib, n_samples,
                    PDM_MPI_DOUBLE, PDM_MPI_SUM, comm);

  free(l_distrib);

  /* Define the cumulative frequency related to g_distribution */

  cfreq[0] = 0.;
  for (id = 0; id < n_samples; id++) {
    double _g_distrib  = (double)g_distrib[id];
    double _gsum_weight = (double)gsum_weight;
    cfreq[id+1] = cfreq[id] + _g_distrib/_gsum_weight;
  }
  cfreq[n_samples] = 1.0;

#if 0 && defined(DEBUG) && !defined(DEBUG) /* For debugging purpose only */

  if (cs_glob_rank_id <= 0) {

    FILE  *dbg_file = NULL;
    int  len;
    static int  loop_id1 = 0;

    len = strlen("DistribOutput_l.dat")+1+2;
    char  *rfilename = (char *) malloc (len * sizeof(char));
    sprintf(rfilename, "DistribOutput_l%02d.dat", loop_id1);

    loop_id1++;

    dbg_file = fopen(rfilename, "w");

    fprintf(dbg_file,
            "# Sample_id  |  OptCfreq  |  Cfreq  |  Sampling  |"
            "Global Distrib\n");
    for (i = 0; i < n_samples; i++)
      fprintf(dbg_file, "%8d %15.5f %15.10f %15.10f %10u\n",
              i, (double)i/(double)n_samples, cfreq[i],
              (double)(sampling[i]), distrib[i]);
    fprintf(dbg_file, "%8d %15.5f %15.10f %15.10f %10u\n",
            i, 1.0, 1.0, 1.0, 0);

    fclose(dbg_file);
    free(rfilename);

  }

#endif /* debugging output */

  /* Convert global distribution from n_samples to n_ranks */

  for (rank_id = 0; rank_id < n_ranks; rank_id++) {

    double   sum = 0.;
    int   shift = rank_id * sampling_factor;

    for (id = 0; id < sampling_factor; id++)
      sum += g_distrib[shift + id];
    g_distrib[rank_id] = sum;

  } /* End of loop on ranks */

#if 0 && defined(DEBUG) && !defined(NDEBUG) /* Sanity check in debug */
  {
    PDM_g_num_t   sum = 0;
    for (rank_id = 0; rank_id < n_ranks; rank_id++)
      sum += g_distrib[rank_id];

    if (sum != gsum_weight)
      PDM_error(__FILE__, __LINE__, 0,
                "Error while computing global distribution.\n"
                "sum = %u and gsum_weight = %u\n",
                sum, gsum_weight);
    exit(1);
  }
#endif /* sanity check */

}

/**
 * \brief Update a distribution associated to sampling to assume a well-balanced
 * distribution of the leaves of the tree.
 *
 *   \param [in]    dim      1D, 2D or 3D
 *   \param [in]    n_ranks  number of ranks (= number of ranges)
 *   \param [inout] c_freq   cumulative frequency array
 *   \param [inout] sampling pointer to pointer to a sampling array
 *   \param [in]    comm     mpi communicator
 */

static void
_update_sampling(int            dim,
                 int            n_ranks,
                 double         c_freq[],
                 PDM_g_num_t  *sampling[])
{
  int  i, j, next_id;
  double  target_freq, f_high, f_low, delta;
  PDM_g_num_t  s_low, s_high;

  PDM_g_num_t  *new_sampling = NULL, *_sampling = *sampling;

  const int  sampling_factor = _sampling_factors[dim];
  const int  n_samples = sampling_factor * n_ranks;
  const double  unit = 1/(double)n_samples;

  /* Compute new_sampling */

  new_sampling = ( PDM_g_num_t  *) malloc (sizeof(PDM_g_num_t) * (n_samples + 1));

  new_sampling[0] = _sampling[0];

  next_id = 1;

  for (i = 0; i < n_samples; i++) {

    target_freq = (i+1)*unit;

    /* Find the next id such as c_freq[next_id] >= target_freq */

    for (j = next_id; j < n_samples + 1; j++) {
      if (c_freq[j] >= target_freq) {
        next_id = j;
        break;
      }
    }

    /* Find new s such as new_s is equal to target_freq by
       a linear interpolation */

    f_low = c_freq[next_id-1];
    f_high = c_freq[next_id];

    s_low = _sampling[next_id-1];
    s_high = _sampling[next_id];

    if (f_high - f_low > 0) {
      delta = (target_freq - f_low) * (s_high - s_low) / (f_high - f_low);
      new_sampling[i+1] = (PDM_g_num_t) (s_low + delta);
    }
    else /* f_high = f_low */
      new_sampling[i+1] = (PDM_g_num_t) (s_low + 0.5 * (s_low + s_high));

#if 0 && defined(DEBUG) && !defined(NDEBUG)
    PDM_printf( " <_update_distrib> (rank: %d) delta: %g, target: %g,"
               " next_id: %d, f_low: %g, f_high: %g, s_low: %g, s_high: %g\n"
               "\t => new_sampling: %g\n",
               cs_glob_rank_id, delta, target_freq, next_id,
               f_low, f_high, s_low, s_high, new_sampling[i+1]);
#endif

  } /* End of loop on samples */

  new_sampling[n_samples] = _sampling[n_samples];

  free(_sampling);

  /* Return pointers */

  *sampling = new_sampling;
}


/**
 *
 * \brief  Define active ranks
 *
 * \param [inout]   ptb          Part to block structure
 *
 */

static void
_active_ranks
(
 PDM_part_to_block_t  *ptb
)
{

  assert (ptb->active_ranks == NULL);

  if (ptb->s_comm == 1) {
    ptb->is_my_rank_active = 1;
    ptb->n_active_ranks = 1;
    ptb->active_ranks = (int *) malloc(sizeof(int) * ptb->n_active_ranks);
    ptb->active_ranks[0] = ptb->i_rank;
  }

  else {

    switch (ptb->t_distrib) {

    case PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC : {
      ptb->is_my_rank_active = 1;
      ptb->n_active_ranks = ptb->s_comm;
      ptb->active_ranks   = (int *) malloc(sizeof(int) * ptb->n_active_ranks);
      for (int i = 0; i < ptb->n_active_ranks; i++) {
        ptb->active_ranks[i] = i;
      }
      break;
    }

    case PDM_PART_TO_BLOCK_DISTRIB_ONE_PROC_PER_NODE :
    case PDM_PART_TO_BLOCK_DISTRIB_PART_OF_NODE : {

      ptb->is_my_rank_active = 0;

      int rankInNode = PDM_io_mpi_node_rank(ptb->comm);
      if (rankInNode == 0) {
        ptb->is_my_rank_active = 1;
      }

      int *tag_active_ranks = (int *) malloc(sizeof(int) * ptb->s_comm);

      PDM_MPI_Allgather((void *) &ptb->is_my_rank_active, 1, PDM_MPI_INT,
                    (void *) tag_active_ranks, 1, PDM_MPI_INT,
                    ptb->comm);

      ptb->n_active_ranks = 0;
      for (int i = 0; i < ptb->s_comm; i++) {
        if (tag_active_ranks[i] == 1) {
          ptb->n_active_ranks += 1;
        }
      }

      ptb->active_ranks   = (int *) malloc(sizeof(int) * ptb->n_active_ranks);
      int n_active_ranks = 0;
      for (int i = 0; i < ptb->s_comm; i++) {
        if (tag_active_ranks[i] == 1) {
          ptb->active_ranks[n_active_ranks++] = i;
        }
      }

      if (ptb->t_distrib == PDM_PART_TO_BLOCK_DISTRIB_ONE_PROC_PER_NODE) {
        break;
      }

      int n_node = ptb->n_active_ranks;
      double part_active_node = PDM_MIN (1, ptb->part_active_node);
      part_active_node = PDM_MAX (0, part_active_node);
      ptb->n_active_ranks = (int) floor (n_node * part_active_node);
      ptb->n_active_ranks = PDM_MAX (1, ptb->n_active_ranks);

      n_active_ranks = 0;
      int coeff = n_node / ptb->n_active_ranks;
      for (int i = 0; i < n_node; i++) {
        if (i % coeff == 0) {
          ptb->active_ranks[n_active_ranks++] = ptb->active_ranks[i];
        }
        if (n_active_ranks == ptb->n_active_ranks) {
          break;
        }
      }

      assert (n_active_ranks == ptb->n_active_ranks);

      break;
    }

    default : {
      PDM_error(__FILE__, __LINE__, 0,"Error cs_part_to_bloc : unknown distribute type\n");
      abort();
    }
    }

    /* Dump */

    if (0 == 1) {
      PDM_printf("active ranks : ");
      for(int i = 0; i < ptb->n_active_ranks; i++)
        PDM_printf("%i ", ptb->active_ranks[i]);
      PDM_printf("\n");
    }
  }
}


/**
 *
 * \brief Distrib data
 *
 * \param [inout] ptb              Part to block structure
 * \param [in]    n_totalData      Total number of data
 * \param [out]   data_distrib_index Element global number
 * \param [out]   s_block_min       Local number of elements
 * \param [out]   s_block_max       Number of partition
 *
 */

static void
_distrib_data
(
 PDM_part_to_block_t *ptb,
 int                  user_distrib
)
{
  PDM_g_num_t _id_max = 0;
  PDM_g_num_t _id_max_max = 0;

  ptb->n_elt_proc= 0;
  //for (int i = 0; i < ptb->n_part; i++) {
  //  ptb->n_elt_proc+= ptb->n_elt[i];
  //  for (int j = 0; j < ptb->n_elt[i]; j++) {
  //    _id_max = _MAX (_id_max, ptb->gnum_elt[i][j]);
  //  }
  //}
  // fflush(stdout);

  for (int i = 0; i < ptb->n_part; i++) {

    ptb->n_elt_proc+= ptb->n_elt[i];
    if(user_distrib == 0) {
      for (int j = 0; j < ptb->n_elt[i]; j++) {
        PDM_g_num_t gnum = PDM_ABS(ptb->gnum_elt[i][j]);
        _id_max = PDM_MAX (_id_max, gnum);
      }
    }
  }



  if(user_distrib == 0) {

    PDM_MPI_Allreduce (&_id_max,
                       &_id_max_max,
                       1,
                       PDM__PDM_MPI_G_NUM,
                       PDM_MPI_MAX,
                       ptb->comm);

    if (ptb->weight == NULL) {
      PDM_g_num_t _n_rank_data = _id_max_max / ptb->n_active_ranks;
      PDM_g_num_t _rest = _id_max_max % ptb->n_active_ranks;

      int n_rank_data = (int) (_n_rank_data);
      int rest       = (int) (_rest);

      ptb->s_block_max = n_rank_data;
      ptb->s_block_min = n_rank_data;

      if (rest != 0) {
        ptb->s_block_max += 1;
      }

      for (int i = 0; i < ptb->s_comm + 1; i++) {
        ptb->data_distrib_index[i] = 0;
      }

      int k = 0;
      int idx = 0;
      for (int i = 0; i < ptb->s_comm; i++) {
        ptb->data_distrib_index[i+1] +=  ptb->data_distrib_index[i];
        if (idx < ptb->n_active_ranks) {
          if (ptb->active_ranks[idx] == i) {
            ptb->data_distrib_index[i+1] += n_rank_data;
            if (k < rest)
              ptb->data_distrib_index[i+1] += 1;
            k += 1;
            idx++;
          }

        }
      }

    }

    else {
      const int dim = 2;
      const int  n_active_ranks = ptb->n_active_ranks;
      const int  sampling_factor = _sampling_factors[dim];
      const int  n_samples = sampling_factor * n_active_ranks;

      double **weight = ptb->weight;
      PDM_MPI_Comm comm = ptb->comm;

      PDM_g_num_t *sampling = malloc(sizeof(PDM_g_num_t) * (n_samples + 1));

      double  lsum_weight = 0.;
      for (int i = 0; i < ptb->n_part; i++) {
        for (int j = 0; j < ptb->n_elt[i]; j++) {
          lsum_weight += weight[i][j];
        }
      }

      double  gsum_weight = 0.;
      PDM_MPI_Allreduce(&lsum_weight, &gsum_weight, 1,
                        PDM_MPI_DOUBLE, PDM_MPI_SUM, comm);

      double optim = gsum_weight / n_active_ranks;

      /* Define a naive sampling (uniform distribution) */

      PDM_g_num_t _n_sampleData = _id_max_max / n_samples;
      PDM_g_num_t _samplerest = _id_max_max % n_samples;

      sampling[0] = 0;
      int k = 0;
      for (int i = 0; i < n_samples; i++) {
        sampling[i+1] = sampling[i];
        sampling[i+1] += _n_sampleData;
        if (k < _samplerest) {
          sampling[i+1] += 1;
          k += 1;
        }
      }

      /* Define the distribution associated to the current sampling array */

      double *distrib = (double *) malloc (sizeof(double) * n_samples);
      double  *cfreq = (double *) malloc (sizeof(double) * (n_samples + 1));

      _define_rank_distrib(ptb,
                           dim,
                           n_active_ranks,
                           gsum_weight,
                           sampling,
                           cfreq,
                           distrib,
                           comm);

      /* Initialize best choice */

      double fit = _evaluate_distribution(n_active_ranks, distrib, optim);
      double best_fit = fit;

      PDM_g_num_t  *best_sampling =
        (PDM_g_num_t  *) malloc (sizeof(PDM_g_num_t) * (n_samples + 1));

      for (int i = 0; i < (n_samples + 1); i++) {
        best_sampling[i] = sampling[i];
      }

      /* Loop to get a better sampling array */

      for (int n_iters = 0;
           (   n_iters < pdm_part_to_block_distrib_n_iter_max
               && fit > pdm_part_to_block_distrib_tol);
           n_iters++)  {

        _update_sampling(dim, n_active_ranks, cfreq, &sampling);

        /* Compute the new distribution associated to the new sampling */

        _define_rank_distrib(ptb,
                             dim,
                             n_active_ranks,
                             gsum_weight,
                             sampling,
                             cfreq,
                             distrib,
                             comm);

        fit = _evaluate_distribution(n_active_ranks, distrib, optim);

        /* Save the best sampling array and its fit */

        if (fit < best_fit) {

          best_fit = fit;
          for (int i = 0; i < (n_samples + 1); i++) {
            best_sampling[i] = sampling[i];
          }
        }

      } /* End of while */

      free (distrib);
      free (cfreq);
      free (sampling);

      sampling = best_sampling;

      int *_active_ranks = ptb->active_ranks;

      PDM_g_num_t *rank_index = malloc (sizeof(PDM_g_num_t) * (n_active_ranks + 1));

      for (int i = 0; i < n_active_ranks + 1; i++) {
        int id = i * sampling_factor;
        rank_index[i] = sampling[id];
      }

      free (sampling);

      ptb->data_distrib_index[0] = 0;

      k = 0;
      for (int i = 0; i < n_active_ranks; i++) {
        int i_activeRank = _active_ranks[i];
        while (k < i_activeRank) {
          ptb->data_distrib_index[k+1] = ptb->data_distrib_index[k];
          k++;
        }
        ptb->data_distrib_index[k+1] = rank_index[i+1];
        k++;
      }
      while (k < ptb->s_comm) {
        ptb->data_distrib_index[k+1] = ptb->data_distrib_index[k];
        k++;
      }

      free (rank_index);

    }
  } // If User



  /* Affichage */

  if (1 == 0) {
    if (ptb->i_rank == 0) {
      PDM_printf("data_distrib_index : ");
      for(int i = 0; i < ptb->s_comm + 1; i++)
        PDM_printf(PDM_FMT_G_NUM" ", ptb->data_distrib_index[i]);
      PDM_printf("\n");
    }
  }

  ptb->n_send_data = (int *) malloc (sizeof(int) * ptb->s_comm);
  ptb->n_recv_data = (int *) malloc (sizeof(int) * ptb->s_comm);

  /* Pour chaque donnee le proc ou elle va etre envoyee */

  ptb->dest_proc = (int *) malloc (sizeof(int) * ptb->n_elt_proc);

  /* Calcul du nombre de donnees a envoyer a chaque procesus */

  for (int i = 0; i < ptb->s_comm; i++) {
    ptb->n_send_data[i] = 0;
  }

  int idx = -1;
  for (int i = 0; i < ptb->n_part; i++) {

    for (int j = 0; j < ptb->n_elt[i]; j++) {

      PDM_g_num_t _gnum_elt = PDM_ABS(ptb->gnum_elt[i][j]) - 1;

      int iproc = PDM_binary_search_gap_long (_gnum_elt,
                                              ptb->data_distrib_index,
                                              ptb->s_comm + 1);

      ptb->dest_proc[++idx] = iproc;
      // assert (ptb->dest_proc[idx] >= 0);
      ptb->n_send_data[iproc] += 1;
    }
  }

  PDM_MPI_Alltoall (ptb->n_send_data, 1, PDM_MPI_INT,
                    ptb->n_recv_data, 1, PDM_MPI_INT,
                    ptb->comm);

  ptb->i_send_data = (int *) malloc(sizeof(int) * ptb->s_comm);
  ptb->i_recv_data = (int *) malloc(sizeof(int) * ptb->s_comm);

  ptb->i_send_data[0] = 0;
  ptb->i_recv_data[0] = 0;
  for (int i = 1; i < ptb->s_comm; i++) {
    ptb->i_send_data[i] = ptb->i_send_data[i-1] + ptb->n_send_data[i-1];
    ptb->i_recv_data[i] = ptb->i_recv_data[i-1] + ptb->n_recv_data[i-1];
  }

  ptb->tn_recv_data = ptb->i_recv_data[ptb->s_comm - 1] +
                      ptb->n_recv_data[ptb->s_comm - 1];

  ptb->tn_send_data = ptb->i_send_data[ptb->s_comm - 1] +
                      ptb->n_send_data[ptb->s_comm - 1];

  PDM_g_num_t *send_gnum =
    (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) * ptb->tn_send_data) ;

  for (int i = 0; i < ptb->s_comm; i++)
    ptb->n_send_data[i] = 0;

  idx = 0;
  for (int i = 0; i < ptb->n_part; i++) {
    for (int j = 0; j < ptb->n_elt[i]; j++) {
      int iproc = ptb->dest_proc[idx];
      send_gnum[ptb->i_send_data[iproc] +
                ptb->n_send_data[iproc]] = PDM_ABS(ptb->gnum_elt[i][j]);
      idx++;
      ptb->n_send_data[iproc] += 1;
    }
  }

  ptb->sorted_recv_gnum =
    (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * ptb->tn_recv_data);

  PDM_MPI_Alltoallv(send_gnum,
                    ptb->n_send_data,
                    ptb->i_send_data,
                    PDM__PDM_MPI_G_NUM,
                    ptb->sorted_recv_gnum,
                    ptb->n_recv_data,
                    ptb->i_recv_data,
                    PDM__PDM_MPI_G_NUM,
                    ptb->comm);

  free(send_gnum);

  /*
   * Sort
   */
  if( ptb->t_post != PDM_PART_TO_BLOCK_POST_MERGE_UNIFORM)
  {
    ptb->order = malloc (sizeof(int) * ptb->tn_recv_data);

    char *env_var = NULL;
    int counting_sort = 0;
    env_var = getenv ("PDM_PTB_COUNTING_SORT");
    if (env_var != NULL) {
      counting_sort = atoi(env_var);
    }

    if(counting_sort == 0) {
      for (int i = 0; i < ptb->tn_recv_data; i++) {
        ptb->order[i] = i;
      }

      PDM_sort_long (ptb->sorted_recv_gnum,
                     ptb->order,
                     ptb->tn_recv_data);
    } else {
      if(ptb->i_rank == 0) {
        printf("_counting_sort_long \n");
      }
      _counting_sort_long(ptb);
    }

    // if(0 == 1) {
    //   PDM_log_trace_array_long(ptb->sorted_recv_gnum, ptb->tn_recv_data, "sorted_recv_gnum");
    //   PDM_log_trace_array_int (ptb->order           , ptb->tn_recv_data, "order");
    // }

  }

  ptb->n_elt_block = ptb->tn_recv_data;

  if( ptb->t_post == PDM_PART_TO_BLOCK_POST_MERGE_UNIFORM)
  {
    ptb->n_elt_block = ptb->data_distrib_index[ptb->i_rank+1] - ptb->data_distrib_index[ptb->i_rank];
    // printf("ptb->n_elt_block::%d \n", ptb->n_elt_block);
    // printf(" ptb->data_distrib_index[ptb->i_rank+1]::"PDM_FMT_G_NUM" \n", ptb->data_distrib_index[ptb->i_rank+1]);
    // printf(" ptb->data_distrib_index[ptb->i_rank]::"PDM_FMT_G_NUM" \n", ptb->data_distrib_index[ptb->i_rank]);
  }

  /*
   * Cleanup
   */

  if ( (ptb->t_post != PDM_PART_TO_BLOCK_POST_NOTHING) &&
       (ptb->t_post != PDM_PART_TO_BLOCK_POST_MERGE_UNIFORM) ) {

    int n_elt_block = 0;

    ptb->block_gnum = malloc (sizeof(PDM_g_num_t) * ptb->tn_recv_data);

    for (int i = 0; i < ptb->tn_recv_data; i++) {
      if (i == 0) {
        ptb->block_gnum[n_elt_block++] = ptb->sorted_recv_gnum[i];
      }
      else if (ptb->block_gnum[n_elt_block-1] != ptb->sorted_recv_gnum[i]) {
        ptb->block_gnum[n_elt_block++] = ptb->sorted_recv_gnum[i];
      }
    }
    ptb->n_elt_block = n_elt_block;
    ptb->block_gnum = realloc (ptb->block_gnum, sizeof(PDM_g_num_t) * ptb->n_elt_block);

    ptb->block_gnum_count = malloc (sizeof(int) * ptb->n_elt_block);
    n_elt_block = 0;
    if (ptb->tn_recv_data > 0)
      ptb->block_gnum_count[0] = 1;
    for (int i = 1; i < ptb->tn_recv_data; i++) {
      //If same than previous, juste increase stride counter
      if (ptb->sorted_recv_gnum[i-1] == ptb->sorted_recv_gnum[i]) {
        ptb->block_gnum_count[n_elt_block]++;
      }
      else { //Otherwise, count next
        ptb->block_gnum_count[++n_elt_block] = 1;
      }
    }
  }

  else {

    ptb->block_gnum = ptb->sorted_recv_gnum;

    ptb->block_gnum_count = PDM_array_const_int(ptb->tn_recv_data, 1);
  }
}



static void
_compute_global_weights
(
 PDM_part_to_block_t *ptb
 )
{
  /* Send local weights */
  int *send_count = PDM_array_zeros_int (ptb->s_comm);
  double *part_weight = malloc (sizeof(double) * ptb->tn_send_data);
  int idx = 0;
  for (int i = 0; i < ptb->n_part; i++) {
    for (int j = 0; j < ptb->n_elt[i]; j++) {
      int rank = ptb->dest_proc[idx++];
      part_weight[ptb->i_send_data[rank] + send_count[rank]++] = ptb->weight[i][j];
    }
  }

  double *recv_weight = malloc (sizeof(double) * ptb->tn_recv_data);
  PDM_MPI_Alltoallv (part_weight,
                     ptb->n_send_data,
                     ptb->i_send_data,
                     PDM_MPI_DOUBLE,
                     recv_weight,
                     ptb->n_recv_data,
                     ptb->i_recv_data,
                     PDM_MPI_DOUBLE,
                     ptb->comm);

  /* Sum received weights */
  double *block_weight = malloc (sizeof(double) * ptb->n_elt_block);
  for (int i = 0; i < ptb->n_elt_block; i++) {
    block_weight[i] = 0.;
  }

  idx = 0;
  for (int i = 0; i < ptb->tn_recv_data; i++) {
    int j = ptb->order[i];
    while (ptb->block_gnum[idx] < ptb->sorted_recv_gnum[i]) {
      idx++;
    }
    block_weight[idx] += recv_weight[j];
  }

  idx = 0;
  for (int i = 0; i < ptb->tn_recv_data; i++) {
    int j = ptb->order[i];
    while (ptb->block_gnum[idx] < ptb->sorted_recv_gnum[i]) {
      idx++;
    }
    recv_weight[j] = block_weight[idx];
  }
  free (block_weight);

  /* Send back global weights */
  PDM_MPI_Alltoallv (recv_weight,
                     ptb->n_recv_data,
                     ptb->i_recv_data,
                     PDM_MPI_DOUBLE,
                     part_weight,
                     ptb->n_send_data,
                     ptb->i_send_data,
                     PDM_MPI_DOUBLE,
                     ptb->comm);
  free (recv_weight);

  /* Store global weights */
  ptb->weight_g = malloc (sizeof(double *) * ptb->n_part);
  PDM_array_reset_int (send_count, ptb->s_comm, 0);
  idx = 0;
  for (int i = 0; i < ptb->n_part; i++) {
    ptb->weight_g[i] = malloc (sizeof(double) * ptb->n_elt[i]);
    for (int j = 0; j < ptb->n_elt[i]; j++) {
      int rank = ptb->dest_proc[idx++];
      ptb->weight_g[i][j] = part_weight[ptb->i_send_data[rank] + send_count[rank]++];
    }
  }
  free (send_count);
  free (part_weight);
}
/*=============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief Reset global statistic
 *
 */

void
PDM_part_to_block_global_statistic_reset
(
  void
)
{
  for (int i = 0; i < 3; i++) {
    t_elaps[i] = 0;
    t_cpu[i] = 0;
  }

  for (int i = 0; i < 2; i++) {
    min_exch_rank[i] = INT_MAX;
    max_exch_rank[i] = -1;
    exch_data[i] = 0;
  }
}


/**
 *
 * \brief Get global timer in part to block
 *
 * \param [in]   comm                 MPI communicator
 * \param [out]  min_exch_rank_send   Global min part of ranks used to send 
 * \param [out]  min_exch_rank_recv   Global min part of ranks used to receive 
 * \param [out]  max_exch_rank_send   Global max part of ranks used to send 
 * \param [out]  max_exch_rank_recv   Global max part of ranks used to receive
 * \param [out]  min_exch_data_send   Global min sent data for a rank 
 * \param [out]  min_exch_data_recv   Global min received data for a rank
 * \param [out]  max_exch_data_send   Global max sent data for a rank
 * \param [out]  max_exch_data_recv   Global max received data for a rank
 * 
 */

void
PDM_part_to_block_global_statistic_get
(
 PDM_MPI_Comm comm,
 int *min_exch_rank_send,
 int *min_exch_rank_recv,
 int *max_exch_rank_send,
 int *max_exch_rank_recv,
 unsigned long long *min_exch_data_send,
 unsigned long long *min_exch_data_recv,
 unsigned long long *max_exch_data_send,
 unsigned long long *max_exch_data_recv
)
{
  unsigned long long max_exch_data[2];
  unsigned long long min_exch_data[2];

  PDM_MPI_Allreduce (exch_data, min_exch_data, 2,
                     PDM_MPI_UNSIGNED_LONG_LONG, PDM_MPI_MIN, comm);
  
  PDM_MPI_Allreduce (exch_data, max_exch_data, 2,
                     PDM_MPI_UNSIGNED_LONG_LONG, PDM_MPI_MAX, comm);

  *min_exch_data_send = min_exch_data[0];
  *min_exch_data_recv = min_exch_data[1];
  *max_exch_data_send = max_exch_data[0];
  *max_exch_data_recv = max_exch_data[1];


  int max_max_exch_rank[2];
  int min_min_exch_rank[2];

  PDM_MPI_Allreduce (min_exch_rank, min_min_exch_rank, 2,
                     PDM_MPI_INT, PDM_MPI_MIN, comm);
  
  PDM_MPI_Allreduce (max_exch_rank, max_max_exch_rank, 2,
                     PDM_MPI_INT, PDM_MPI_MAX, comm);

  *min_exch_rank_send = min_min_exch_rank[0];
  *min_exch_rank_recv = min_min_exch_rank[1]; 
  *max_exch_rank_send = max_max_exch_rank[0]; 
  *max_exch_rank_recv = max_max_exch_rank[1];

}


/**
 *
 * \brief Get global timer in part to block
 *
 * \param [in]   comm              MPI communicator
 * \param [out]  min_elaps         Min elapsed time
 * \param [out]  max_elaps         Max elapsed time
 * \param [out]  min_cpu           Min cpu time
 * \param [out]  max_cpu           Max cpu time
 * \param [out]  min_elaps_create  Global min elapsed for create function
 * \param [out]  max_elaps_create  Global max elapsed for create function
 * \param [out]  min_cpu_create    Global min cpu for create function
 * \param [out]  max_cpu_create    Global max cpu for create function
 * \param [out]  min_elaps_create2 Global min elapsed for create2 function 
 * \param [out]  max_elaps_create2 Global max elapsed for create2 function
 * \param [out]  min_cpu_create2   Global min cpu for create2 function
 * \param [out]  max_cpu_create2   Global max cpu for create2 function
 * \param [out]  min_elaps_exch    Global min elapsed for exch function
 * \param [out]  max_elaps_exch    Global max elapsed for exch function
 * \param [out]  min_cpu_exch      Global min cpu for exch function
 * \param [out]  max_cpu_exch      Global max cpu for exch function
 * 
 */

void
PDM_part_to_block_global_timer_get
(
 PDM_MPI_Comm comm,
 double       *min_elaps_create,
 double       *max_elaps_create,
 double       *min_cpu_create,
 double       *max_cpu_create,
 double       *min_elaps_create2,
 double       *max_elaps_create2,
 double       *min_cpu_create2,
 double       *max_cpu_create2,
 double       *min_elaps_exch,
 double       *max_elaps_exch,
 double       *min_cpu_exch,
 double       *max_cpu_exch
)
{

  double min_elaps[3];
  double max_elaps[3];
  double min_cpu[3];
  double max_cpu[3];

  PDM_MPI_Allreduce (t_elaps, min_elaps, 3,
                     PDM_MPI_DOUBLE, PDM_MPI_MIN, comm);
  
  PDM_MPI_Allreduce (t_elaps, max_elaps, 3,
                     PDM_MPI_DOUBLE, PDM_MPI_MAX, comm);

  PDM_MPI_Allreduce (t_cpu, min_cpu, 3,
                     PDM_MPI_DOUBLE, PDM_MPI_MIN, comm);
  
  PDM_MPI_Allreduce (t_cpu, max_cpu, 3,
                     PDM_MPI_DOUBLE, PDM_MPI_MAX, comm);

  *min_elaps_create  = min_elaps[0];
  *max_elaps_create  = max_elaps[0];
  *min_cpu_create    = min_cpu[0];
  *max_cpu_create    = max_cpu[0];
  *min_elaps_create2 = min_elaps[1]; 
  *max_elaps_create2 = max_elaps[1]; 
  *min_cpu_create2   = min_cpu[1];
  *max_cpu_create2   = max_cpu[1];
  *min_elaps_exch    = min_elaps[2];
  *max_elaps_exch    = max_elaps[2];
  *min_cpu_exch      = min_cpu[2];
  *max_cpu_exch      = max_cpu[2];

}


/**
 *
 * \brief Create a partitioning to block redistribution
 *
 * \param [in]   t_distrib       Distribution type
 * \param [in]   t_post          Post processing type
 * \param [in]   part_active_node  Part of active nodes (\ref PDM_writer_BLOCK_DISTRIB_PART_OF_NODE mode)
 * \param [in]   gnum_elt        Element global number
 * \param [in]   weight          Weight of elements (or NULL)
 * \param [in]   n_elt           Local number of elements
 * \param [in]   n_part          Number of partition
 * \param [in]   comm            MPI communicator
 *
 * \return   Initialized PDM_part_to_block_t
 *
 */

PDM_part_to_block_t *
PDM_part_to_block_create
(
 PDM_part_to_block_distrib_t   t_distrib,
 PDM_part_to_block_post_t      t_post,
 double                        part_active_node,
 PDM_g_num_t                 **gnum_elt,
 double                      **weight,
 int                          *n_elt,
 int                           n_part,
 PDM_MPI_Comm                  comm
)
{

  if (n_ptb == 0) {
    t_timer[0] = PDM_timer_create ();
    t_timer[1] = PDM_timer_create ();
    t_timer[2] = PDM_timer_create ();
  }
  n_ptb++;

  double t1_elaps = PDM_timer_elapsed(t_timer[0]);
  double t1_cpu = PDM_timer_cpu(t_timer[0]);
  PDM_timer_resume(t_timer[0]);

  PDM_part_to_block_t *ptb =
    (PDM_part_to_block_t *) malloc (sizeof(PDM_part_to_block_t));

  if ((t_post != PDM_PART_TO_BLOCK_POST_MERGE) && (weight != NULL)) {
    PDM_error(__FILE__, __LINE__, 0,"PDM_part_to_block_create : weights are available only if PDM_PART_TO_BLOCK_POST_MERGE is selected\n");
  }

  ptb->t_distrib         = t_distrib;    /*!< Distribution type */
  ptb->t_post            = t_post;       /*!< Post processing type */
  ptb->n_active_ranks    = 0;            /*!< Number of active ranks */
  ptb->active_ranks      = NULL;         /*!< List of active ranks */
  ptb->comm              = comm;         /*!< MSG communicator */
  PDM_MPI_Comm_size (comm, &(ptb->s_comm));
  PDM_MPI_Comm_rank (comm, &(ptb->i_rank));
  ptb->is_my_rank_active   = 0;              /*!< Is active current rank */
  ptb->part_active_node    = part_active_node; /*!< Part of active nodes */

  ptb->n_part             = n_part;       /*!< Number of parts */
  ptb->n_elt              = n_elt;        /*!< Number of elements for any part */
  ptb->n_elt_proc         = 0;            /*!< Number of elements on the current rank */
  ptb->gnum_elt           = gnum_elt;     /*!< Global numbering of elements for any part */
  ptb->weight             = weight;
  ptb->dest_proc          = NULL;
  ptb->data_distrib_index =
    (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) * (ptb->s_comm + 1));   /*!< Data distribution on ranks */

  ptb->s_block_min   = INT_MAX;
  ptb->s_block_max   = 0;

  ptb->i_send_data   = NULL;  /*!< Data to send to other processes index (size = s_comm) */
  ptb->i_recv_data   = NULL;  /*!< Received Data from other processes index (size = s_comm) */
  ptb->n_send_data   = NULL;  /*!< Number of data to send to other processes (size = s_comm) */
  ptb->n_recv_data   = NULL;  /*!< Number of received Data from other processes (size = s_comm) */

  ptb->tn_send_data      = 0;     /*!< Total number of sended data */
  ptb->tn_recv_data      = 0;     /*!< Total number of received data */
  ptb->sorted_recv_gnum  = NULL;  /*!< Sorted recv global num */
  ptb->order             = NULL;  /*!< Order */
  ptb->n_elt_block       = 0;
  ptb->block_gnum        = NULL;  /*!< Global number of reveived data (size = tn_recv_data) */
  ptb->block_gnum_count  = NULL;  /*!< Number of occurence of reveived data (size = tn_recv_data) */

  ptb->weight_g          = NULL; /*!< Global weights of elements for any part */

  /* Asynchone */
  ptb->max_exch_request = 10;
  ptb->next_request     = 0;
  ptb->s_data           = (size_t          * ) malloc ( ptb->max_exch_request * sizeof(size_t           ) );
  ptb->t_stride         = (PDM_stride_t    * ) malloc ( ptb->max_exch_request * sizeof(PDM_stride_t     ) );
  ptb->cst_stride       = (int             * ) malloc ( ptb->max_exch_request * sizeof(int              ) );
  ptb->wait_status      = (int             * ) malloc ( ptb->max_exch_request * sizeof(int              ) );
  ptb->request_mpi      = (PDM_MPI_Request * ) malloc ( ptb->max_exch_request * sizeof(PDM_MPI_Request  ) );

  ptb->send_buffer      = (unsigned char  ** ) malloc ( ptb->max_exch_request * sizeof(unsigned char   *) );
  ptb->recv_buffer      = (unsigned char  ** ) malloc ( ptb->max_exch_request * sizeof(unsigned char   *) );
  ptb->recv_stride      = (int            ** ) malloc ( ptb->max_exch_request * sizeof(int             *) );
  ptb->n_send_buffer    = (int            ** ) malloc ( ptb->max_exch_request * sizeof(int             *) );
  ptb->i_send_buffer    = (int            ** ) malloc ( ptb->max_exch_request * sizeof(int             *) );
  ptb->n_recv_buffer    = (int            ** ) malloc ( ptb->max_exch_request * sizeof(int             *) );
  ptb->i_recv_buffer    = (int            ** ) malloc ( ptb->max_exch_request * sizeof(int             *) );

  for(int i_req = 0; i_req < ptb->max_exch_request; ++i_req) {
    ptb->send_buffer  [i_req] = NULL;
    ptb->recv_buffer  [i_req] = NULL;
    ptb->recv_stride  [i_req] = NULL;
    ptb->n_send_buffer[i_req] = NULL;
    ptb->i_send_buffer[i_req] = NULL;
    ptb->n_recv_buffer[i_req] = NULL;
    ptb->i_recv_buffer[i_req] = NULL;
  }
  /*
   * Active ranks definition
   */

  _active_ranks (ptb);

  /*
   * Data distribution definition
   */

  _distrib_data (ptb, 0);


  /*
   * Compute global weight for each element
   */
  if (ptb->weight != NULL && ptb->t_post == PDM_PART_TO_BLOCK_POST_MERGE) {
    _compute_global_weights (ptb);
  }

  int n_rank_recv = 0;
  int n_rank_send = 0;

  for (int i = 0; i < ptb->s_comm; i++) {
    if (ptb->i_rank != i && ptb->n_recv_data[i] > 0) {
      n_rank_recv += 1;
    }
    if (ptb->i_rank != i && ptb->n_send_data[i] > 0) {
      n_rank_send += 1;
    }
  }

  max_exch_rank[0] = PDM_MAX(max_exch_rank[0], n_rank_send); 
  max_exch_rank[1] = PDM_MAX(max_exch_rank[1], n_rank_recv); 
  min_exch_rank[0] = PDM_MIN(min_exch_rank[0], n_rank_send); 
  min_exch_rank[1] = PDM_MIN(min_exch_rank[1], n_rank_recv); 

  PDM_timer_hang_on(t_timer[0]);
  double t2_elaps = PDM_timer_elapsed(t_timer[0] );
  double t2_cpu = PDM_timer_cpu(t_timer[0]);

  t_elaps[0] += (t2_elaps - t1_elaps);
  t_cpu[0] += (t2_cpu - t1_cpu);

  return (PDM_part_to_block_t *) ptb;
}


/**
 *
 * \brief Create a partitioning to block redistribution
 *
 * \param [in]   t_distrib       Distribution type
 * \param [in]   t_post          Post processing type
 * \param [in]   part_active_node  Part of active nodes (\ref PDM_writer_BLOCK_DISTRIB_PART_OF_NODE mode)
 * \param [in]   gnum_elt        Element global number
 * \param [in]   weight          Weight of elements (or NULL)
 * \param [in]   n_elt           Local number of elements
 * \param [in]   n_part          Number of partition
 * \param [in]   comm            MPI communicator
 *
 * \return   Initialized PDM_part_to_block_t
 *
 */

PDM_part_to_block_t *
PDM_part_to_block_create_from_distrib
(
 PDM_part_to_block_distrib_t   t_distrib,
 PDM_part_to_block_post_t      t_post,
 double                         part_active_node,
 PDM_g_num_t                 **gnum_elt,
 const PDM_g_num_t            *data_distrib_index,
 int                          *n_elt,
 int                           n_part,
 PDM_MPI_Comm                  comm
)
{

  if (n_ptb == 0) {
    t_timer[0] = PDM_timer_create ();
    t_timer[1] = PDM_timer_create ();
    t_timer[2] = PDM_timer_create ();
  }
  n_ptb++;

  double t1_elaps = PDM_timer_elapsed(t_timer[1]);
  double t1_cpu = PDM_timer_cpu(t_timer[1]);
  PDM_timer_resume(t_timer[1]);

  PDM_part_to_block_t *ptb =
    (PDM_part_to_block_t *) malloc (sizeof(PDM_part_to_block_t));

  ptb->t_distrib         = t_distrib;    /*!< Distribution type */
  ptb->t_post            = t_post;       /*!< Post processing type */
  ptb->n_active_ranks    = 0;            /*!< Number of active ranks */
  ptb->active_ranks      = NULL;         /*!< List of active ranks */
  ptb->comm              = comm;         /*!< MSG communicator */
  PDM_MPI_Comm_size (comm, &(ptb->s_comm));
  PDM_MPI_Comm_rank (comm, &(ptb->i_rank));
  ptb->is_my_rank_active   = 0;                /*!< Is active current rank */
  ptb->part_active_node    = part_active_node; /*!< Part of active nodes */

  ptb->n_part             = n_part;       /*!< Number of parts */
  ptb->n_elt              = n_elt;        /*!< Number of elements for any part */
  ptb->n_elt_proc         = 0;            /*!< Number of elements on the current rank */
  ptb->gnum_elt           = gnum_elt;     /*!< Global numbering of elements for any part */
  ptb->weight             = NULL;
  ptb->dest_proc          = NULL;
  ptb->data_distrib_index =
    (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) * (ptb->s_comm + 1));   /*!< Data distribution on ranks */

  for(int i_rank = 0; i_rank < ptb->s_comm+1; i_rank++){
    ptb->data_distrib_index[i_rank] = data_distrib_index[i_rank];
  }

  ptb->s_block_min   = INT_MAX;
  ptb->s_block_max   = 0;

  ptb->i_send_data   = NULL;  /*!< Data to send to other processes index (size = s_comm) */
  ptb->i_recv_data   = NULL;  /*!< Received Data from other processes index (size = s_comm) */
  ptb->n_send_data   = NULL;  /*!< Number of data to send to other processes (size = s_comm) */
  ptb->n_recv_data   = NULL;  /*!< Number of received Data from other processes (size = s_comm) */

  ptb->tn_send_data     = 0;     /*!< Total number of sended data */
  ptb->tn_recv_data     = 0;     /*!< Total number of received data */
  ptb->sorted_recv_gnum = NULL;  /*!< Sorted recv global num */
  ptb->order            = NULL;  /*!< Order */
  ptb->n_elt_block      = 0;
  ptb->block_gnum       = NULL;  /*!< Global number of reveived data (size = tn_recv_data) */
  ptb->block_gnum_count = NULL;  /*!< Number of occurence of reveived data (size = tn_recv_data) */

  ptb->weight_g          = NULL; /*!< Global weights of elements for any part */

  /* Asynchone */
  ptb->max_exch_request = 10;
  ptb->next_request     = 0;
  ptb->s_data           = (size_t          * ) malloc ( ptb->max_exch_request * sizeof(size_t           ) );
  ptb->t_stride         = (PDM_stride_t    * ) malloc ( ptb->max_exch_request * sizeof(PDM_stride_t     ) );
  ptb->cst_stride       = (int             * ) malloc ( ptb->max_exch_request * sizeof(int              ) );
  ptb->wait_status      = (int             * ) malloc ( ptb->max_exch_request * sizeof(int              ) );
  ptb->request_mpi      = (PDM_MPI_Request * ) malloc ( ptb->max_exch_request * sizeof(PDM_MPI_Request  ) );

  ptb->send_buffer      = (unsigned char  ** ) malloc ( ptb->max_exch_request * sizeof(unsigned char   *) );
  ptb->recv_buffer      = (unsigned char  ** ) malloc ( ptb->max_exch_request * sizeof(unsigned char   *) );
  ptb->recv_stride      = (int            ** ) malloc ( ptb->max_exch_request * sizeof(int             *) );
  ptb->n_send_buffer    = (int            ** ) malloc ( ptb->max_exch_request * sizeof(int             *) );
  ptb->i_send_buffer    = (int            ** ) malloc ( ptb->max_exch_request * sizeof(int             *) );
  ptb->n_recv_buffer    = (int            ** ) malloc ( ptb->max_exch_request * sizeof(int             *) );
  ptb->i_recv_buffer    = (int            ** ) malloc ( ptb->max_exch_request * sizeof(int             *) );

  for(int i_req = 0; i_req < ptb->max_exch_request; ++i_req) {
    ptb->send_buffer  [i_req] = NULL;
    ptb->recv_buffer  [i_req] = NULL;
    ptb->recv_stride  [i_req] = NULL;
    ptb->n_send_buffer[i_req] = NULL;
    ptb->i_send_buffer[i_req] = NULL;
    ptb->n_recv_buffer[i_req] = NULL;
    ptb->i_recv_buffer[i_req] = NULL;
  }

  /*
   * Active ranks definition
   */
  _active_ranks (ptb);

  /*
   * Data distribution definition
   */
  // Rajouter un attribut de classe pour passer user ?
  _distrib_data (ptb, 1);

  int n_rank_recv = 0;
  int n_rank_send = 0;


  for (int i = 0; i < ptb->s_comm; i++) {
    if (ptb->i_rank != i && ptb->n_recv_data[i] > 0) {
      n_rank_recv += 1;
    }
    if (ptb->i_rank != i && ptb->n_send_data[i] > 0) {
      n_rank_send += 1;
    }
  }

  max_exch_rank[0] = PDM_MAX(max_exch_rank[0], n_rank_send); 
  max_exch_rank[1] = PDM_MAX(max_exch_rank[1], n_rank_recv); 
  min_exch_rank[0] = PDM_MIN(min_exch_rank[0], n_rank_send); 
  min_exch_rank[1] = PDM_MIN(min_exch_rank[1], n_rank_recv); 

  PDM_timer_hang_on(t_timer[1]);
  double t2_elaps = PDM_timer_elapsed(t_timer[1]);
  double t2_cpu = PDM_timer_cpu(t_timer[1]);

  t_elaps[1] += (t2_elaps - t1_elaps);
  t_cpu[1] += (t2_cpu - t1_cpu);

  return (PDM_part_to_block_t *) ptb;
}


/**
 *
 * \brief Return number of active ranks
 *
 * \param [in]   ptb          Part to block structure
 *
 * \return Number of active ranks
 *
 */

int
PDM_part_to_block_n_active_ranks_get
(
 PDM_part_to_block_t *ptb
)
{
  return ptb->n_active_ranks;
}


/**
 *
 * \brief Return if current rank is active
 *
 * \param [in]   ptb          Part to block structure
 *
 * \return  if current rank is active
 *
 */

int
PDM_part_to_block_is_active_rank
(
 PDM_part_to_block_t *ptb
 )
{
  return ptb->is_my_rank_active;
}


/**
 *
 * \brief Return active ranks
 *
 * \param [in]   ptb          Part to block structure
 *
 * \return  active ranks
 *
 */

int *
PDM_part_to_block_active_ranks_get
(
 PDM_part_to_block_t *ptb
)
{
  return ptb->active_ranks;
}


/**
 *
 * \brief Return number of element in the current process
 *
 * \param [in]   ptb          Part to block structure
 *
 * \return Number of element in the current process
 *
 */

int
PDM_part_to_block_n_elt_block_get
(
 PDM_part_to_block_t *ptb
)
{
  return ptb->n_elt_block;
}


/**
 *
 * \brief Return global numbers of element in the current process
 *
 * \param [in]   ptb          Part to block structure
 *
 * \return  Global numbers
 *
 */

PDM_g_num_t *
PDM_part_to_block_block_gnum_get
(
 PDM_part_to_block_t *ptb
)
{
  return ptb->block_gnum;
}

/**
 *
 * \brief Return numbers of occurence of each gnum element in the current process
 *
 * \param [in]   ptb          Part to block structure
 *
 * \return  Global numbers counter
 *
 */

int *
PDM_part_to_block_block_gnum_count_get
(
 PDM_part_to_block_t *ptb
)
{
  return ptb->block_gnum_count;
}


/**
 *
 * \brief Initialize a data exchange
 *
 * \param [in]   ptb          Part to block structure
 * \param [in]   s_data       Data size
 * \param [in]   t_stride     Stride type
 * \param [in]   var_stride   Variable stride (size = n_part) only for \ref PDM_writer_STRIDE_VAR
 * \param [in]   cst_stride   Stride only for \ref PDM_writer_STRIDE_CST
 *
 * \return       Size of highest block
 *
 */

int
PDM_part_to_block_exch
(
 PDM_part_to_block_t *ptb,
 size_t               s_data,
 PDM_stride_t         t_stride,
 int                  cst_stride,
 int                **part_stride,
 void               **part_data,
 int                **block_stride,
 void               **block_data
)
{

  double t1_elaps = PDM_timer_elapsed(t_timer[2]);
  double t1_cpu = PDM_timer_cpu(t_timer[2]);
  PDM_timer_resume(t_timer[2]);

  //if ( ptb->i_rank==0 ){
  //  printf("line %d PDM_part_to_block_exch s_data=%lu\n",__LINE__,s_data);
  //  fflush(stdout);
  //  printf("line %d PDM_part_to_block_exch t_stride=%d\n",__LINE__,t_stride);
  //  fflush(stdout);
  //  printf("line %d PDM_part_to_block_exch cst_stride=%d\n",__LINE__,cst_stride);
  //  fflush(stdout);
  //  printf("line %d PDM_part_to_block_exch ptb->n_elt[0]=%d\n",__LINE__,ptb->n_elt[0]);
  //  fflush(stdout);
  //  //for (int i = 0; i < ptb->n_part; i++) {
  //  //for (int j = 0; j < ptb->n_elt[0]; j++) {
  //  //  printf("line %d PDM_part_to_block_exch part_data[%d]==%d\n",__LINE__,j,(( int **) part_data)[0][j]);
  //  //  fflush(stdout);
  //  //}
  //  //}
  //}

  //*
  //*/


  if ((ptb->t_post == PDM_PART_TO_BLOCK_POST_MERGE) &&
      (t_stride ==  PDM_STRIDE_CST_INTERLACED)) {
    PDM_error(__FILE__, __LINE__, 0,"PDM_part_to_block_exch : PDM_writer_STRIDE_CST is not compatible PDM_writer_POST_MERGE post\n");
    abort ();
  }


  size_t *i_send_buffer = (size_t *) malloc (sizeof(size_t) * ptb->s_comm);
  size_t *i_recv_buffer = (size_t *) malloc (sizeof(size_t) * ptb->s_comm);
  int *n_send_buffer = (int *) malloc (sizeof(int) * ptb->s_comm);
  int *n_recv_buffer = (int *) malloc (sizeof(int) * ptb->s_comm);

  /*
   * Exchange Stride and build buffer properties
   */

  int *recv_stride = NULL;
  if (t_stride == PDM_STRIDE_VAR_INTERLACED) {

    for (int i = 0; i < ptb->s_comm; i++) {
      n_send_buffer[i] = 0;
      n_recv_buffer[i] = 0;
    }

    int *send_stride = (int *) malloc (sizeof(int) * ptb->tn_send_data);

    int idx = -1;
    for (int i = 0; i < ptb->n_part; i++) {
      for (int j = 0; j < ptb->n_elt[i]; j++) {
        int iproc = ptb->dest_proc[++idx];
        send_stride[ptb->i_send_data[iproc] + n_send_buffer[iproc]] = part_stride[i][j];
        n_send_buffer[iproc] += 1;
      }
    }

    recv_stride = (int *) malloc (sizeof(int) * ptb->tn_recv_data);
    assert (send_stride != NULL);
    assert (recv_stride != NULL);

    PDM_MPI_Alltoallv (send_stride,
                       ptb->n_send_data,
                       ptb->i_send_data,
                       PDM_MPI_INT,
                       recv_stride,
                       ptb->n_recv_data,
                       ptb->i_recv_data,
                       PDM_MPI_INT,
                       ptb->comm);

    /*
     * Build buffers
     */

    for (int i = 0; i < ptb->s_comm; i++) {
      int ibeg = ptb->i_send_data[i];
      int iend = ptb->i_send_data[i] + ptb->n_send_data[i];

      n_send_buffer[i] = 0;
      for (int k = ibeg; k < iend; k++)
        n_send_buffer[i] += send_stride[k];

      n_send_buffer[i] *= (int) s_data;

      if (i > 0) {
        i_send_buffer[i] = i_send_buffer[i-1] + n_send_buffer[i-1];
      }
      else {
        i_send_buffer[i] = 0;
      }

      ibeg = ptb->i_recv_data[i];
      iend = ptb->i_recv_data[i] + ptb->n_recv_data[i];

      n_recv_buffer[i] = 0;
      for (int k = ibeg; k < iend; k++)
        n_recv_buffer[i] += recv_stride[k];

      n_recv_buffer[i] *= (int) s_data;

      if (i > 0)
        i_recv_buffer[i] = i_recv_buffer[i-1] + n_recv_buffer[i-1];
      else
        i_recv_buffer[i] = 0;

    }

    free(send_stride);

  }

  else if (t_stride == PDM_STRIDE_CST_INTERLACED) {

    for (int i = 0; i < ptb->s_comm; i++) {

      i_send_buffer[i] = ptb->i_send_data[i] * cst_stride * (int) s_data;
      i_recv_buffer[i] = ptb->i_recv_data[i] * cst_stride * (int) s_data;

      n_send_buffer[i] = ptb->n_send_data[i] * cst_stride * (int) s_data;
      n_recv_buffer[i] = ptb->n_recv_data[i] * cst_stride * (int) s_data;

    }
  }

  size_t s_send_buffer = i_send_buffer[ptb->s_comm - 1] + n_send_buffer[ptb->s_comm -1];
  size_t s_recv_buffer = i_recv_buffer[ptb->s_comm - 1] + n_recv_buffer[ptb->s_comm -1];

  /*
   * Data exchange
   */


  unsigned char *send_buffer = (unsigned char *) malloc(sizeof(unsigned char) * s_send_buffer);
  unsigned char *recv_buffer = (unsigned char *) malloc(sizeof(unsigned char) * s_recv_buffer);

  assert (send_buffer != NULL);
  assert (recv_buffer != NULL);

  unsigned char **_part_data = (unsigned char **) part_data;

  for (int i = 0; i <  ptb->s_comm; i++) {
    n_send_buffer[i] = 0;
  }

  int idx = -1;
  for (int i = 0; i < ptb->n_part; i++) {

    size_t *i_part = NULL;
    if (t_stride == PDM_STRIDE_VAR_INTERLACED) {
      i_part = (size_t *) malloc (sizeof(size_t) * (ptb->n_elt[i] + 1));
      assert (i_part != NULL);

    i_part[0] = 0;
      for (int j = 1; j < ptb->n_elt[i] + 1; j++)
        i_part[j] = i_part[j-1] + ((size_t) part_stride[i][j-1] * s_data);
    }

    //if ( ptb->i_rank==0 ){
    //  printf("line %d PDM_part_to_block_exch  ptb->dest_proc[%d]=%d\n",__LINE__,idx,ptb->dest_proc[idx]);
    //  fflush(stdout);
    //  //printf("line %d PDM_part_to_block_exch s_octet_elt=%d\n",__LINE__,s_octet_elt);
    //  //fflush(stdout);
    //}


    for (int j = 0; j < ptb->n_elt[i]; j++) {
      int iproc = ptb->dest_proc[++idx];
      size_t s_octet_elt = 0;
      size_t i_part_elt = 0;

      if (t_stride == PDM_STRIDE_CST_INTERLACED) {
        s_octet_elt = cst_stride * (int) s_data;
        i_part_elt  = cst_stride * (int) s_data * j;
      }

      else if (t_stride == PDM_STRIDE_VAR_INTERLACED) {
        s_octet_elt = i_part[j+1] - i_part[j];
        i_part_elt  = i_part[j];
      }

      for (int k = 0; k < (int) s_octet_elt; k++) {
        send_buffer[i_send_buffer[iproc] + n_send_buffer[iproc]++] =
          _part_data[i][i_part_elt + k];
      }


    }

    if (i_part != NULL)
      free (i_part);
  }

  PDM_MPI_Alltoallv_l(send_buffer,
                      n_send_buffer,
                      i_send_buffer,
                      PDM_MPI_BYTE,
                      recv_buffer,
                      n_recv_buffer,
                      i_recv_buffer,
                      PDM_MPI_BYTE,
                      ptb->comm);

  for (int i = 0; i < ptb->s_comm; i++) {
    if (ptb->i_rank != i) {
      exch_data[1] += n_recv_buffer[i];
      exch_data[0] += n_send_buffer[i];
    }
  }

  free(send_buffer);
  free(n_send_buffer);
  free(i_send_buffer);
  free(n_recv_buffer);
  free(i_recv_buffer);

  unsigned char *_block_data = malloc(sizeof(unsigned char) * s_recv_buffer);
  assert(_block_data != NULL);

  *block_data = _block_data;
  int *i_recv_stride = NULL;
  int *i_block_stride = NULL;
  int s_block_data = ((int) sizeof(unsigned char) * s_recv_buffer) / (int) s_data;

  if (t_stride == PDM_STRIDE_VAR_INTERLACED) {
    *block_stride = NULL;
    int* _block_stride = NULL;
    if(ptb->tn_recv_data > 0){
      _block_stride = malloc(sizeof(int) * ptb->tn_recv_data);
    }
    *block_stride = _block_stride;
    for (int i = 0; i < ptb->tn_recv_data; i++) {
      _block_stride[i] = recv_stride[ptb->order[i]];
    }

    /*
     * Compute index in data
     */

    i_recv_stride = malloc (sizeof(int) * (ptb->tn_recv_data + 1));
    i_block_stride = malloc (sizeof(int) * (ptb->tn_recv_data + 1));

    i_recv_stride[0] = 0;
    i_block_stride[0] = 0;
    for (int i = 0; i < ptb->tn_recv_data; i++) {
      i_recv_stride[i+1]  = i_recv_stride[i] + recv_stride[i];
      i_block_stride[i+1] = i_block_stride[i] + _block_stride[i];
    }

    for (int i = 0; i < ptb->tn_recv_data; i++) {
      i_recv_stride[i+1]  *= (int) s_data;
      i_block_stride[i+1] *= (int) s_data;
    }

    /*
     * Sort Buffer
     */

    for (int i = 0; i < ptb->tn_recv_data; i++) {
      int old = ptb->order[i];
      int idOld = i_recv_stride[old];
      for (int k = i_block_stride[i]; k < i_block_stride[i+1]; k++) {
        _block_data[k] = recv_buffer[idOld++];
      }
    }

//    free (recv_buffer);
    free (recv_stride);
    free (i_recv_stride);

    /*
     * post processing
     */

    if (ptb->t_post != PDM_PART_TO_BLOCK_POST_NOTHING) {

      int idx1 = 0;
      int idx2 = 0;

      if (ptb->tn_recv_data == 1) {
        idx2 = i_block_stride[1];
      }

      for (int i = 1; i < ptb->tn_recv_data; i++) {
        if (i == 1) {
          idx2 = i_block_stride[1];
        }
        if (ptb->block_gnum[idx1] != ptb->sorted_recv_gnum[i]) {
          idx1 += 1;
          _block_stride[idx1] = _block_stride[i];
          if (ptb->t_post == PDM_PART_TO_BLOCK_POST_CLEANUP) {
            for (int k = i_block_stride[i]; k < i_block_stride[i+1]; k++) {
              _block_data[idx2++] = _block_data[k];
            }
          }
        }
        else {
          if (ptb->t_post == PDM_PART_TO_BLOCK_POST_MERGE) {
            _block_stride[idx1] += _block_stride[i];
          }
        }
      }

      /* Cleanup */

      if (ptb->t_post == PDM_PART_TO_BLOCK_POST_CLEANUP) {
        _block_data = realloc (_block_data, sizeof(unsigned char) * idx2);
        *block_data = _block_data;

        _block_stride = realloc (_block_stride, sizeof(int) * ptb->n_elt_block);

        *block_stride = _block_stride;
        s_block_data = idx2 / (int) s_data;

      }

    }

    free (i_block_stride);

  }

  else {

    /*
     * Sort Buffer
     */

    for (int i = 0; i < ptb->tn_recv_data; i++) {
      int n_octet = cst_stride * (int) s_data;
      int old = ptb->order[i];
      int idOld = old * n_octet;

      for (int k = i*n_octet; k < (i+1)*n_octet; k++) {
        _block_data[k] = recv_buffer[idOld++];
      }
    }

    /*
     * Post processing
     */

    if (ptb->t_post != PDM_PART_TO_BLOCK_POST_NOTHING) {
      int idx2 = 0;
      int idx1 = 0;

      assert (ptb->t_post != PDM_PART_TO_BLOCK_POST_MERGE);

      if (ptb->tn_recv_data == 1) {
        idx2 =  cst_stride * (int) s_data;
      }

      for (int i = 1; i < ptb->tn_recv_data; i++) {
        int n_octet = cst_stride * (int) s_data;
        if (i == 1) {
          idx2 = n_octet;
        }
        if (ptb->block_gnum[idx1] != ptb->sorted_recv_gnum[i]) {
          idx1 += 1;
          if (ptb->t_post == PDM_PART_TO_BLOCK_POST_CLEANUP) {
            int idx3 = i * cst_stride * (int) s_data;
            for (int k = 0; k < n_octet; k++) {
              _block_data[idx2++] = _block_data[idx3++];
            }
          }
        }
      }

      if (ptb->t_post == PDM_PART_TO_BLOCK_POST_CLEANUP) {
        _block_data = realloc (_block_data, sizeof(unsigned char) * idx2);
        *block_data = _block_data;
        s_block_data = idx2 / (int) s_data;
      }

    }
  }

  free (recv_buffer);

  PDM_timer_hang_on(t_timer[2]);
  double t2_elaps = PDM_timer_elapsed(t_timer[2]);
  double t2_cpu = PDM_timer_cpu(t_timer[2]);

  t_elaps[2] += (t2_elaps - t1_elaps);
  t_cpu[2] += (t2_cpu - t1_cpu);

  return s_block_data;

}

/**
 *
 * \brief Initialize a data exchange
 *
 * \param [in]   ptb          Part to block structure
 * \param [in]   s_data       Data size
 * \param [in]   t_stride     Stride type
 * \param [in]   var_stride   Variable stride (size = n_part) only for \ref PDM_writer_STRIDE_VAR
 * \param [in]   cst_stride   Stride only for \ref PDM_writer_STRIDE_CST
 *
 * \return       Size of highest block
 *
 */
int
PDM_part_to_block_async_exch
(
 PDM_part_to_block_t *ptb,
 size_t               s_data,
 PDM_stride_t         t_stride,
 int                  cst_stride,
 int                **part_stride,
 void               **part_data
)
{

  if ((ptb->t_post == PDM_PART_TO_BLOCK_POST_MERGE) &&
      (t_stride ==  PDM_STRIDE_CST_INTERLACED)) {
    PDM_error(__FILE__, __LINE__, 0,"PDM_part_to_block_exch : PDM_writer_STRIDE_CST is not compatible PDM_writer_POST_MERGE post\n");
    abort ();
  }

  /*
   * Take next id for message and buffer
   */
  int request_id = ptb->next_request++;

  ptb->i_send_buffer[request_id] = (int *) malloc (sizeof(int) * ptb->s_comm);
  ptb->i_recv_buffer[request_id] = (int *) malloc (sizeof(int) * ptb->s_comm);
  ptb->n_send_buffer[request_id] = (int *) malloc (sizeof(int) * ptb->s_comm);
  ptb->n_recv_buffer[request_id] = (int *) malloc (sizeof(int) * ptb->s_comm);

  /* Store additionnal information necessary for post-process */
  ptb->s_data    [request_id] = s_data;
  ptb->t_stride  [request_id] = t_stride;
  ptb->cst_stride[request_id] = cst_stride;

  /* Short cut */
  int* i_send_buffer = ptb->i_send_buffer[request_id];
  int* i_recv_buffer = ptb->i_recv_buffer[request_id];
  int* n_send_buffer = ptb->n_send_buffer[request_id];
  int* n_recv_buffer = ptb->n_recv_buffer[request_id];

  /*
   * Exchange Stride and build buffer properties
   */

  if (t_stride == PDM_STRIDE_VAR_INTERLACED) {

    for (int i = 0; i < ptb->s_comm; i++) {
      n_send_buffer[i] = 0;
      n_recv_buffer[i] = 0;
    }

    int *send_stride = (int *) malloc (sizeof(int) * ptb->tn_send_data);

    int idx = -1;
    for (int i = 0; i < ptb->n_part; i++) {
      for (int j = 0; j < ptb->n_elt[i]; j++) {
        int iproc = ptb->dest_proc[++idx];
        send_stride[ptb->i_send_data[iproc] + n_send_buffer[iproc]] = part_stride[i][j];
        n_send_buffer[iproc] += 1;
      }
    }

    ptb->recv_stride[request_id] = (int *) malloc (sizeof(int) * ptb->tn_recv_data);

    int* recv_stride = ptb->recv_stride[request_id];

    PDM_MPI_Alltoallv (send_stride,
                       ptb->n_send_data,
                       ptb->i_send_data,
                       PDM_MPI_INT,
                       recv_stride,
                       ptb->n_recv_data,
                       ptb->i_recv_data,
                       PDM_MPI_INT,
                       ptb->comm);

    /*
     * Build buffers
     */

    for (int i = 0; i < ptb->s_comm; i++) {
      int ibeg = ptb->i_send_data[i];
      int iend = ptb->i_send_data[i] + ptb->n_send_data[i];

      n_send_buffer[i] = 0;
      for (int k = ibeg; k < iend; k++)
        n_send_buffer[i] += send_stride[k];

      n_send_buffer[i] *= (int) s_data;

      if (i > 0) {
        i_send_buffer[i] = i_send_buffer[i-1] + n_send_buffer[i-1];
      }
      else {
        i_send_buffer[i] = 0;
      }

      ibeg = ptb->i_recv_data[i];
      iend = ptb->i_recv_data[i] + ptb->n_recv_data[i];

      n_recv_buffer[i] = 0;
      for (int k = ibeg; k < iend; k++)
        n_recv_buffer[i] += recv_stride[k];

      n_recv_buffer[i] *= (int) s_data;

      if (i > 0)
        i_recv_buffer[i] = i_recv_buffer[i-1] + n_recv_buffer[i-1];
      else
        i_recv_buffer[i] = 0;

    }

    free(send_stride);

  }

  else if (t_stride == PDM_STRIDE_CST_INTERLACED) {

    for (int i = 0; i < ptb->s_comm; i++) {

      i_send_buffer[i] = ptb->i_send_data[i] * cst_stride * (int) s_data;
      i_recv_buffer[i] = ptb->i_recv_data[i] * cst_stride * (int) s_data;

      n_send_buffer[i] = ptb->n_send_data[i] * cst_stride * (int) s_data;
      n_recv_buffer[i] = ptb->n_recv_data[i] * cst_stride * (int) s_data;

    }

  }

  int s_send_buffer = i_send_buffer[ptb->s_comm - 1] + n_send_buffer[ptb->s_comm -1];
  int s_recv_buffer = i_recv_buffer[ptb->s_comm - 1] + n_recv_buffer[ptb->s_comm -1];

  /*
   * Data exchange
   */
  ptb->send_buffer[request_id] = (unsigned char *) malloc(sizeof(unsigned char) * s_send_buffer);
  ptb->recv_buffer[request_id] = (unsigned char *) malloc(sizeof(unsigned char) * s_recv_buffer);

  /* Shortcut */
  unsigned char *send_buffer = ptb->send_buffer[request_id];
  unsigned char *recv_buffer = ptb->recv_buffer[request_id];

  unsigned char **_part_data = (unsigned char **) part_data;

  for (int i = 0; i <  ptb->s_comm; i++) {
    n_send_buffer[i] = 0;
  }

  int idx = -1;
  for (int i = 0; i < ptb->n_part; i++) {

    int *i_part = NULL;
    if (t_stride == PDM_STRIDE_VAR_INTERLACED) {
      i_part = (int *) malloc (sizeof(int) * (ptb->n_elt[i] + 1));

      i_part[0] = 0;
      for (int j = 1; j < ptb->n_elt[i] + 1; j++)
        i_part[j] = i_part[j-1] + (part_stride[i][j-1] * (int) s_data);
    }

    for (int j = 0; j < ptb->n_elt[i]; j++) {
      int iproc = ptb->dest_proc[++idx];
      int s_octet_elt = 0;
      int i_part_elt = 0;

      if (t_stride == PDM_STRIDE_CST_INTERLACED) {
        s_octet_elt = cst_stride * (int) s_data;
        i_part_elt  = cst_stride * (int) s_data * j;
      }

      else if (t_stride == PDM_STRIDE_VAR_INTERLACED) {
        s_octet_elt = i_part[j+1] - i_part[j];
        i_part_elt  = i_part[j];
      }

      for (int k = 0; k < s_octet_elt; k++) {
        send_buffer[i_send_buffer[iproc] + n_send_buffer[iproc]++] =
          _part_data[i][i_part_elt + k];
      }
    }

    if (i_part != NULL)
      free (i_part);
  }

  // PDM_MPI_Request* request_mpi =

  PDM_MPI_Ialltoallv(send_buffer,
                     n_send_buffer,
                     i_send_buffer,
                     PDM_MPI_BYTE,
                     recv_buffer,
                     n_recv_buffer,
                     i_recv_buffer,
                     PDM_MPI_BYTE,
                     ptb->comm,
                     &ptb->request_mpi[request_id]);

  ptb->wait_status[request_id] = 0;

  return request_id;
}


/**
 *
 * \brief Wait for an exchange
 *
 * \param [in]   ptb          Part to block structure
 * \param [in]   request_id   Internal id of the current exchange
 *
 */
void
PDM_part_to_block_async_wait
(
 PDM_part_to_block_t *ptb,
 int                  request_id
)
{
  // printf("PDM_part_to_block_async_wait::request_id::%d \n", request_id);

  assert(ptb->wait_status[request_id] == 0);

  int code = PDM_MPI_Wait(&ptb->request_mpi[request_id]);
  assert(code == PDM_MPI_SUCCESS);

  ptb->wait_status[request_id] = 1;

}

/**
 *
 * \brief Get the raw exchange buffer and stride and deallocate memory
 *
 * \param [in]   ptb          Part to block structure
 * \param [in]   request_id   Internal id of the current exchange
 * \param [out]  block_stride Block stride
 * \param [out]  block_data   Block data
 *
 */
int
PDM_part_to_block_asyn_get_raw
(
 PDM_part_to_block_t *ptb,
 int                  request_id,
 int                **block_stride,
 void               **block_data
)
{
  // printf("PDM_part_to_block_asyn_get_raw::request_id::%d \n", request_id);

  assert(ptb->wait_status[request_id] == 1);

  free(ptb->send_buffer  [request_id]);
  free(ptb->n_send_buffer[request_id]);
  free(ptb->i_send_buffer[request_id]);
  free(ptb->n_recv_buffer[request_id]);
  free(ptb->i_recv_buffer[request_id]);

  ptb->send_buffer  [request_id] = NULL;
  ptb->n_send_buffer[request_id] = NULL;
  ptb->i_send_buffer[request_id] = NULL;
  ptb->n_recv_buffer[request_id] = NULL;
  ptb->i_recv_buffer[request_id] = NULL;

  /* Mv pointer */
  *block_stride = ptb->recv_stride[request_id];
  *block_data   = ptb->recv_buffer[request_id];

  /* Nulliffy - User is now owner of the excahnge data */
  ptb->recv_stride[request_id] = NULL;
  ptb->recv_buffer[request_id] = NULL;

  ptb->wait_status[request_id] = 2;

  return ptb->tn_recv_data;

}


/**
 *
 * \brief Post-treatment of the recv buffer
 *
 * \param [in]   ptb          Part to block structure
 * \param [in]   request_id   Internal id of the current exchange
 * \param [out]  block_stride Block stride
 * \param [out]  block_data   Block data
 *
 */
void
PDM_part_to_block_asyn_post_treatment
(
 PDM_part_to_block_t *ptb,
 int                  request_id,
 int                **block_stride,
 void               **block_data
)
{
  PDM_UNUSED(block_stride);
  PDM_UNUSED(block_data);

  assert(ptb->wait_status[request_id] == 1);

  // size_t s_send_buffer = ptb->i_send_buffer[request_id][ptb->s_comm - 1] + ptb->n_send_buffer[request_id][ptb->s_comm -1];
  size_t s_recv_buffer = ptb->i_recv_buffer[request_id][ptb->s_comm - 1] + ptb->n_recv_buffer[request_id][ptb->s_comm -1];

  free(ptb->send_buffer  [request_id]);
  free(ptb->n_send_buffer[request_id]);
  free(ptb->i_send_buffer[request_id]);
  free(ptb->n_recv_buffer[request_id]);
  free(ptb->i_recv_buffer[request_id]);

  ptb->send_buffer  [request_id] = NULL;
  ptb->n_send_buffer[request_id] = NULL;
  ptb->i_send_buffer[request_id] = NULL;
  ptb->n_recv_buffer[request_id] = NULL;
  ptb->i_recv_buffer[request_id] = NULL;

  size_t s_data         = ptb->s_data    [request_id];
  PDM_stride_t t_stride = ptb->t_stride  [request_id];
  int cst_stride        = ptb->cst_stride[request_id];

  /* Classic post-treatment */
  unsigned char *_block_data = malloc(sizeof(unsigned char) * s_recv_buffer);
  assert(_block_data != NULL);

  int* recv_stride = ptb->recv_stride[request_id];
  unsigned char* recv_buffer = ptb->recv_buffer[request_id];

  *block_data = _block_data;
  *block_stride = NULL;
  int *i_recv_stride = NULL;
  int *i_block_stride = NULL;
  int s_block_data = ((int) sizeof(unsigned char) * s_recv_buffer) / (int) s_data;

  if (t_stride == PDM_STRIDE_VAR_INTERLACED) {
    int* _block_stride = NULL;
    if(ptb->tn_recv_data > 0){
      _block_stride = malloc(sizeof(int) * ptb->tn_recv_data);
    }
    *block_stride = _block_stride;
    for (int i = 0; i < ptb->tn_recv_data; i++) {
      _block_stride[i] = recv_stride[ptb->order[i]];
    }

    /*
     * Compute index in data
     */

    i_recv_stride = malloc (sizeof(int) * (ptb->tn_recv_data + 1));
    i_block_stride = malloc (sizeof(int) * (ptb->tn_recv_data + 1));

    i_recv_stride[0] = 0;
    i_block_stride[0] = 0;
    for (int i = 0; i < ptb->tn_recv_data; i++) {
      i_recv_stride[i+1]  = i_recv_stride[i] + recv_stride[i];
      i_block_stride[i+1] = i_block_stride[i] + _block_stride[i];
    }

    for (int i = 0; i < ptb->tn_recv_data; i++) {
      i_recv_stride[i+1]  *= (int) s_data;
      i_block_stride[i+1] *= (int) s_data;
    }

    /*
     * Sort Buffer
     */

    for (int i = 0; i < ptb->tn_recv_data; i++) {
      int old = ptb->order[i];
      int idOld = i_recv_stride[old];
      for (int k = i_block_stride[i]; k < i_block_stride[i+1]; k++) {
        _block_data[k] = recv_buffer[idOld++];
      }
    }

//    free (recv_buffer);
    // free (recv_stride);
    free (i_recv_stride);

    /*
     * post processing
     */

    if (ptb->t_post != PDM_PART_TO_BLOCK_POST_NOTHING) {

      int idx1 = 0;
      int idx2 = 0;

      if (ptb->tn_recv_data == 1) {
        idx2 = i_block_stride[1];
      }

      for (int i = 1; i < ptb->tn_recv_data; i++) {
        if (i == 1) {
          idx2 = i_block_stride[1];
        }
        if (ptb->block_gnum[idx1] != ptb->sorted_recv_gnum[i]) {
          idx1 += 1;
          _block_stride[idx1] = _block_stride[i];
          if (ptb->t_post == PDM_PART_TO_BLOCK_POST_CLEANUP) {
            for (int k = i_block_stride[i]; k < i_block_stride[i+1]; k++) {
              _block_data[idx2++] = _block_data[k];
            }
          }
        }
        else {
          if (ptb->t_post == PDM_PART_TO_BLOCK_POST_MERGE) {
            _block_stride[idx1] += _block_stride[i];
          }
        }
      }

      /* Cleanup */

      if (ptb->t_post == PDM_PART_TO_BLOCK_POST_CLEANUP) {
        _block_data = realloc (_block_data, sizeof(unsigned char) * idx2);
        *block_data = _block_data;

        _block_stride = realloc (_block_stride, sizeof(int) * ptb->n_elt_block);

        *block_stride = _block_stride;
        s_block_data = idx2 / (int) s_data;

      }

    }

    free (i_block_stride);

  }

  else {

    /*
     * Sort Buffer
     */

    for (int i = 0; i < ptb->tn_recv_data; i++) {
      int n_octet = cst_stride * (int) s_data;
      int old = ptb->order[i];
      int idOld = old * n_octet;

      for (int k = i*n_octet; k < (i+1)*n_octet; k++) {
        _block_data[k] = recv_buffer[idOld++];
      }
    }

    /*
     * Post processing
     */

    if (ptb->t_post != PDM_PART_TO_BLOCK_POST_NOTHING) {
      int idx2 = 0;
      int idx1 = 0;

      assert (ptb->t_post != PDM_PART_TO_BLOCK_POST_MERGE);

      if (ptb->tn_recv_data == 1) {
        idx2 =  cst_stride * (int) s_data;
      }

      for (int i = 1; i < ptb->tn_recv_data; i++) {
        int n_octet = cst_stride * (int) s_data;
        if (i == 1) {
          idx2 = n_octet;
        }
        if (ptb->block_gnum[idx1] != ptb->sorted_recv_gnum[i]) {
          idx1 += 1;
          if (ptb->t_post == PDM_PART_TO_BLOCK_POST_CLEANUP) {
            int idx3 = i * cst_stride * (int) s_data;
            for (int k = 0; k < n_octet; k++) {
              _block_data[idx2++] = _block_data[idx3++];
            }
          }
        }
      }

      if (ptb->t_post == PDM_PART_TO_BLOCK_POST_CLEANUP) {
        _block_data = realloc (_block_data, sizeof(unsigned char) * idx2);
        *block_data = _block_data;
        s_block_data = idx2 / (int) s_data;
      }

    }
  }

  /*
   * Free
   */
  free(ptb->recv_stride  [request_id]);
  free(ptb->recv_buffer  [request_id]);
  ptb->recv_stride[request_id] = NULL;
  ptb->recv_buffer[request_id] = NULL;

  ptb->wait_status[request_id] = 2;

  PDM_UNUSED(s_block_data);

}


/**
 *
 * \brief Free a part to block structure
 *
 * \param [inout] ptb         Part to block structure
 *
 * \return       NULL
 */

PDM_part_to_block_t *
PDM_part_to_block_free
(
 PDM_part_to_block_t *ptb
)
{

  if (ptb->active_ranks != NULL) {
    free (ptb->active_ranks);
    ptb->active_ranks = NULL;
  }
  if (ptb->dest_proc != NULL) {
    free (ptb->dest_proc);
    ptb->dest_proc = NULL;
  }
  if (ptb->data_distrib_index != NULL) {
    free (ptb->data_distrib_index);
    ptb->data_distrib_index = NULL;
  }
  if (ptb->i_send_data != NULL) {
    free (ptb->i_send_data);
    ptb->i_send_data = NULL;
  }
  if (ptb->i_recv_data != NULL) {
    free (ptb->i_recv_data);
    ptb->i_recv_data = NULL;
  }
  if (ptb->n_send_data != NULL) {
    free (ptb->n_send_data);
    ptb->n_send_data = NULL;
  }
  if (ptb->n_recv_data != NULL) {
    free (ptb->n_recv_data);
    ptb->n_recv_data = NULL;
  }
  if (ptb->sorted_recv_gnum != NULL) {
    free (ptb->sorted_recv_gnum);
    ptb->sorted_recv_gnum = NULL;
  }
  if (ptb->order != NULL) {
    free (ptb->order);
    ptb->order = NULL;
  }

  // if ((ptb->t_post != PDM_PART_TO_BLOCK_POST_NOTHING) && (ptb->block_gnum != NULL)) {
  //   free (ptb->block_gnum);
  //   ptb->block_gnum = NULL;
  // }
  if ((ptb->t_post != PDM_PART_TO_BLOCK_POST_NOTHING      ) &&
      (ptb->t_post != PDM_PART_TO_BLOCK_POST_MERGE_UNIFORM) && (ptb->block_gnum != NULL)) {
    free (ptb->block_gnum);
    ptb->block_gnum = NULL;
  }
  if (ptb->block_gnum_count != NULL) {
    free(ptb->block_gnum_count);
    ptb->block_gnum_count=NULL;
  }

  if (ptb->weight_g != NULL) {
    for (int i = 0; i < ptb->n_part; i++) {
      free (ptb->weight_g[i]);
    }
    free (ptb->weight_g);
  }

  /* This one check if all buffer has been correctly move or delete */
  for(int i_req = 0; i_req < ptb->max_exch_request; ++i_req) {
    assert(ptb->send_buffer  [i_req] == NULL);
    assert(ptb->recv_buffer  [i_req] == NULL);
    assert(ptb->recv_stride  [i_req] == NULL);
    assert(ptb->n_send_buffer[i_req] == NULL);
    assert(ptb->i_send_buffer[i_req] == NULL);
    assert(ptb->n_recv_buffer[i_req] == NULL);
    assert(ptb->i_recv_buffer[i_req] == NULL);
  }

  free(ptb->s_data       );
  free(ptb->t_stride     );
  free(ptb->cst_stride   );
  free(ptb->wait_status  );
  free(ptb->request_mpi  );
  free(ptb->send_buffer  );
  free(ptb->recv_buffer  );
  free(ptb->recv_stride  );
  free(ptb->n_send_buffer);
  free(ptb->i_send_buffer);
  free(ptb->n_recv_buffer);
  free(ptb->i_recv_buffer);

  free (ptb);

  n_ptb--;
  if (n_ptb == 0) {
    PDM_timer_free(t_timer[0]);
    PDM_timer_free(t_timer[1]);
    PDM_timer_free(t_timer[2]);
  }

  return NULL;
}


/**
 *
 * \brief Return block distribution index
 *
 * \param [in] ptb         Part to block structure
 *
 * \return  Distribution (size = communicator size + 1)
 */

PDM_g_num_t *
PDM_part_to_block_distrib_index_get
(
 PDM_part_to_block_t *ptb
)
{
  return ptb->data_distrib_index;
}


/**
 *
 * \brief Return processus destination
 *
 * \param [in] ptb         Part to block structure
 *
 * \return  Destination (size = sum of partition elements)
 */

PDM_l_num_t *
PDM_part_to_block_destination_get
(
 PDM_part_to_block_t *ptb
)
{
  return ptb->dest_proc;
}

PDM_g_num_t*
PDM_part_to_block_adapt_partial_block_to_block
(
 PDM_part_to_block_t  *ptb,
 int                 **block_n,
 PDM_g_num_t           n_g_block
)
{
  PDM_g_num_t *_block_distrib_idx = malloc (sizeof(PDM_g_num_t) * (ptb->s_comm + 1));

  for (int i = 0; i < ptb->s_comm + 1; i++) {
    _block_distrib_idx[i] = ptb->data_distrib_index[i];
  }
  int block_n_elt = PDM_part_to_block_n_elt_block_get (ptb);

  int block_n_elt_tot = _block_distrib_idx[ptb->i_rank+1] - _block_distrib_idx[ptb->i_rank];
  int* block_n_tmp = PDM_array_zeros_int(block_n_elt_tot);

  int* _block_n = *block_n;
  for (int i1 = 0; i1 < block_n_elt; i1++) {
    int i = (int) (ptb->block_gnum[i1] - _block_distrib_idx[ptb->i_rank] - 1);
    // printf(" ptb->block_gnum[%i] = %i --> %i (%i)\n", i1, ptb->block_gnum[i1], _block_distrib_idx[ptb->i_rank], i);
    block_n_tmp[i] = _block_n[i1];
  }

  *block_n = realloc(*block_n, sizeof(int) * block_n_elt_tot);
  _block_n = *block_n;
  for (int i1 = 0; i1 < block_n_elt_tot; i1++) {
    _block_n[i1] = block_n_tmp[i1];
  }

  free(block_n_tmp);

  PDM_g_num_t old_max = _block_distrib_idx[ptb->s_comm];
  PDM_g_num_t new_max = n_g_block;
  int diff_last = (int) (new_max - old_max);

  _block_distrib_idx[ptb->s_comm] = new_max;

  if (ptb->i_rank == (ptb->s_comm - 1)) {
    int new_size = block_n_elt_tot + diff_last;
    *block_n = realloc(*block_n, sizeof(int) * new_size);
    _block_n = *block_n;
    for (int i = block_n_elt_tot; i < new_size; i++) {
      _block_n[i] = 0;
    }
  }

  return _block_distrib_idx;
}





/**
 *
 * \brief Return global weights of element in the current process
 *
 * \param [in]   ptb          Part to block structure
 *
 * \return  Global weights
 *
 */

double **
PDM_part_to_block_global_weight_get
(
 PDM_part_to_block_t *ptb
)
{
  return ptb->weight_g;
}



/**
 *
 * \brief Get number of MPI ranks
 *
 * \param [in]   ptb          Part to block structure
 *
 * \return  Number of MPI ranks
 *
 */

int
PDM_part_to_block_n_ranks_get
(
 PDM_part_to_block_t *ptb
)
{
  return ptb->s_comm;
}


/**
 *
 * \brief Return total number of element in the current process (summed over all partitions)
 *
 * \param [in]   ptb          Part to block structure
 *
 * \return Total number of element in the current process
 *
 */

int
PDM_part_to_block_n_elt_proc_get
(
 PDM_part_to_block_t *ptb
 )
 {
  return ptb->n_elt_proc;
 }

#undef _MIN
#undef _MAX

#ifdef __cplusplus
}
#endif /* __cplusplus */
