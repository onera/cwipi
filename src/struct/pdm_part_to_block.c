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
#include "pdm_priv.h"
#include "pdm_binary_search.h"
#include "pdm_sort.h"
#include "pdm_printf.h"
#include "pdm_error.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

#define _MIN(a,b)   ((a) < (b) ?  (a) : (b))  /* Minimum of a et b */

#define _MAX(a,b)   ((a) > (b) ?  (a) : (b))  /* Maximum of a et b */

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
 * Static function definitions
 *============================================================================*/


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
      d_up = _MAX(d_up, distribution[i] - optim);
    else
      d_low = _MAX(d_low, optim - distribution[i]);

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
_define_rank_distrib( _pdm_part_to_block_t *ptb,
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

  /* gnum are supposed to be ordered */
  /* TODO: Refaire cette partie : les nombres ne sont pas tries */

  for (int i = 0; i < ptb->n_part; i++) {

    for (int j = 0; j < ptb->n_elt[i]; j++) {

      PDM_g_num_t _gnum_elt = ptb->gnum_elt[i][j] - 1;

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
#ifdef __INTEL_COMPILER
#pragma warning(push)
#pragma warning(disable:2259)
#endif
    double _g_distrib  = (double)g_distrib[id];
    double _gsum_weight = (double)gsum_weight;
#ifdef __INTEL_COMPILER
#pragma warning(pop)
#endif
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
 _pdm_part_to_block_t  *ptb
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
      double part_active_node = _MIN (1, ptb->part_active_node);
      part_active_node = _MAX (0, part_active_node);
      ptb->n_active_ranks = (int) floor (n_node * part_active_node);
      ptb->n_active_ranks = _MAX (1, ptb->n_active_ranks);

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
 _pdm_part_to_block_t *ptb,
 int                  user_distrib
)
{

  PDM_g_num_t _id_max = 0;
  PDM_g_num_t _id_max_max = 0;

  ptb->n_elt_proc= 0;
  for (int i = 0; i < ptb->n_part; i++) {
    ptb->n_elt_proc+= ptb->n_elt[i];
    for (int j = 0; j < ptb->n_elt[i]; j++) {
      _id_max = _MAX (_id_max, ptb->gnum_elt[i][j]);
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

      PDM_g_num_t _gnum_elt = ptb->gnum_elt[i][j] - 1;

      int iproc = PDM_binary_search_gap_long (_gnum_elt,
                                              ptb->data_distrib_index,
                                              ptb->s_comm + 1);

      ptb->dest_proc[++idx] = iproc;
      assert (ptb->dest_proc[idx] >= 0);
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
                ptb->n_send_data[iproc]] = ptb->gnum_elt[i][j];
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

  ptb->order = malloc (sizeof(int) * ptb->tn_recv_data);
  for (int i = 0; i < ptb->tn_recv_data; i++) {
    ptb->order[i] = i;
  }

  PDM_sort_long (ptb->sorted_recv_gnum,
                 ptb->order,
                 ptb->tn_recv_data);

  ptb->n_elt_block = ptb->tn_recv_data;

  /*
   * Cleanup
   */

  if (ptb->t_post != PDM_PART_TO_BLOCK_POST_NOTHING) {

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

  }

  else {

    ptb->block_gnum = ptb->sorted_recv_gnum;

  }
}

/*=============================================================================
 * Public function definitions
 *============================================================================*/


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
 * \return   Initialized cs_part_to_block
 *
 */

PDM_part_to_block_t *
PDM_part_to_block_create_cf
(
 PDM_part_to_block_distrib_t   t_distrib,
 PDM_part_to_block_post_t      t_post,
 double                        part_active_node,
 PDM_g_num_t                 **gnum_elt,
 double                      **weight,
 int                          *n_elt,
 int                           n_part,
 PDM_MPI_Fint                  fcomm
)
{
  const PDM_MPI_Comm _comm        = PDM_MPI_Comm_f2c(fcomm);
  return PDM_part_to_block_create (t_distrib, t_post, part_active_node,
                                   gnum_elt, weight, n_elt, n_part, _comm);
}

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

  _pdm_part_to_block_t *ptb =
    (_pdm_part_to_block_t *) malloc (sizeof(_pdm_part_to_block_t));

  if ((t_post != PDM_PART_TO_BLOCK_POST_MERGE) && (weight != NULL)) {
    PDM_error(__FILE__, __LINE__, 0,"PDM_part_to_block_create : weights are available only if PDM_PART_TO_BLOCK_POST_MERGE is selected\n");
  }

  ptb->t_distrib        = t_distrib;    /*!< Distribution type */
  ptb->t_post           = t_post;       /*!< Post processing type */
  ptb->n_active_ranks    = 0;            /*!< Number of active ranks */
  ptb->active_ranks      = NULL;         /*!< List of active ranks */
  ptb->comm             = comm;         /*!< MSG communicator */
  PDM_MPI_Comm_size (comm, &(ptb->s_comm));
  PDM_MPI_Comm_rank (comm, &(ptb->i_rank));
  ptb->is_my_rank_active   = 0;              /*!< Is active current rank */
  ptb->part_active_node   = part_active_node; /*!< Part of active nodes */

  ptb->n_part            = n_part;       /*!< Number of parts */
  ptb->n_elt             = n_elt;        /*!< Number of elements for any part */
  ptb->n_elt_proc        = 0;            /*!< Number of elements on the current rank */
  ptb->gnum_elt          = gnum_elt;     /*!< Global numbering of elements for any part */
  ptb->weight            = weight;
  ptb->dest_proc         = NULL;
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

  /*
   * Active ranks definition
   */

  _active_ranks (ptb);

  /*
   * Data distribution definition
   */

  _distrib_data (ptb, 0);

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
 * \return   Initialized cs_part_to_block
 *
 */

PDM_part_to_block_t *
PDM_part_to_block_create2_cf
(
 PDM_part_to_block_distrib_t   t_distrib,
 PDM_part_to_block_post_t      t_post,
 double                        part_active_node,
 PDM_g_num_t                 **gnum_elt,
 PDM_g_num_t                  *data_distrib_index,
 int                          *n_elt,
 int                           n_part,
 PDM_MPI_Fint                  fcomm
)
{
  const PDM_MPI_Comm _comm        = PDM_MPI_Comm_f2c(fcomm);
  return PDM_part_to_block_create2 (t_distrib, t_post, part_active_node,
                                    gnum_elt, data_distrib_index,
                                    n_elt, n_part, _comm);
}

PDM_part_to_block_t *
PDM_part_to_block_create2
(
 PDM_part_to_block_distrib_t   t_distrib,
 PDM_part_to_block_post_t      t_post,
 double                         part_active_node,
 PDM_g_num_t                 **gnum_elt,
 PDM_g_num_t                  *data_distrib_index,
 int                          *n_elt,
 int                           n_part,
 PDM_MPI_Comm                  comm
)
{

  _pdm_part_to_block_t *ptb =
    (_pdm_part_to_block_t *) malloc (sizeof(_pdm_part_to_block_t));

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

  /*
   * Active ranks definition
   */

  _active_ranks (ptb);

  /*
   * Data distribution definition
   */
  // Rajouter un attribut de classe pour passer user ?
  _distrib_data (ptb, 1);

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
  _pdm_part_to_block_t *_ptb = (_pdm_part_to_block_t *) ptb;
  return _ptb->n_active_ranks;
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
  _pdm_part_to_block_t *_ptb = (_pdm_part_to_block_t *) ptb;
  return _ptb->is_my_rank_active;
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
  _pdm_part_to_block_t *_ptb = (_pdm_part_to_block_t *) ptb;
  return _ptb->active_ranks;
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
  _pdm_part_to_block_t *_ptb = (_pdm_part_to_block_t *) ptb;
  return _ptb->n_elt_block;
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
  _pdm_part_to_block_t *_ptb = (_pdm_part_to_block_t *) ptb;
  return _ptb->block_gnum;
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
  _pdm_part_to_block_t *_ptb = (_pdm_part_to_block_t *) ptb;

  if ((_ptb->t_post == PDM_PART_TO_BLOCK_POST_MERGE) &&
      (t_stride ==  PDM_STRIDE_CST)) {
    PDM_error(__FILE__, __LINE__, 0,"PDM_part_to_block_exch : PDM_writer_STRIDE_CST is not compatible PDM_writer_POST_MERGE post\n");
    abort ();
  }

  size_t *i_send_buffer = (size_t *) malloc (sizeof(size_t) * _ptb->s_comm);
  size_t *i_recv_buffer = (size_t *) malloc (sizeof(size_t) * _ptb->s_comm);
  int *n_send_buffer = (int *) malloc (sizeof(int) * _ptb->s_comm);
  int *n_recv_buffer = (int *) malloc (sizeof(int) * _ptb->s_comm);

  /*
   * Exchange Stride and build buffer properties
   */

  int *recv_stride = NULL;
  if (t_stride == PDM_STRIDE_VAR) {

    for (int i = 0; i < _ptb->s_comm; i++) {
      n_send_buffer[i] = 0;
      n_recv_buffer[i] = 0;
    }

    int *send_stride = (int *) malloc (sizeof(int) * _ptb->tn_send_data);

    int idx = -1;
    for (int i = 0; i < _ptb->n_part; i++) {
      for (int j = 0; j < _ptb->n_elt[i]; j++) {
        int iproc = _ptb->dest_proc[++idx];
        send_stride[_ptb->i_send_data[iproc] + n_send_buffer[iproc]] = part_stride[i][j];
        n_send_buffer[iproc] += 1;
      }
    }

    recv_stride = (int *) malloc (sizeof(int) * _ptb->tn_recv_data);

    PDM_MPI_Alltoallv (send_stride,
                       _ptb->n_send_data,
                       _ptb->i_send_data,
                       PDM_MPI_INT,
                       recv_stride,
                       _ptb->n_recv_data,
                       _ptb->i_recv_data,
                       PDM_MPI_INT,
                       _ptb->comm);

    /*
     * Build buffers
     */

    for (int i = 0; i < _ptb->s_comm; i++) {
      int ibeg = _ptb->i_send_data[i];
      int iend = _ptb->i_send_data[i] + _ptb->n_send_data[i];

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

      ibeg = _ptb->i_recv_data[i];
      iend = _ptb->i_recv_data[i] + _ptb->n_recv_data[i];

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

  else if (t_stride == PDM_STRIDE_CST) {

    for (int i = 0; i < _ptb->s_comm; i++) {

      i_send_buffer[i] = _ptb->i_send_data[i] * cst_stride * (int) s_data;
      i_recv_buffer[i] = _ptb->i_recv_data[i] * cst_stride * (int) s_data;

      n_send_buffer[i] = _ptb->n_send_data[i] * cst_stride * (int) s_data;
      n_recv_buffer[i] = _ptb->n_recv_data[i] * cst_stride * (int) s_data;

    }
    fflush(stdout);
  }

  int s_send_buffer = i_send_buffer[_ptb->s_comm - 1] + n_send_buffer[_ptb->s_comm -1];
  int s_recv_buffer = i_recv_buffer[_ptb->s_comm - 1] + n_recv_buffer[_ptb->s_comm -1];

  /*
   * Data exchange
   */

  unsigned char *send_buffer = (unsigned char *) malloc(sizeof(unsigned char) * s_send_buffer);
  unsigned char *recv_buffer = (unsigned char *) malloc(sizeof(unsigned char) * s_recv_buffer);

  unsigned char **_part_data = (unsigned char **) part_data;

  for (int i = 0; i <  _ptb->s_comm; i++) {
    n_send_buffer[i] = 0;
  }

  int idx = -1;
  for (int i = 0; i < _ptb->n_part; i++) {

    int *i_part = NULL;
    if (t_stride == PDM_STRIDE_VAR) {
      i_part = (int *) malloc (sizeof(int) * (_ptb->n_elt[i] + 1));

      i_part[0] = 0;
      for (int j = 1; j < _ptb->n_elt[i] + 1; j++)
        i_part[j] = i_part[j-1] + (part_stride[i][j-1] * (int) s_data);
    }

    for (int j = 0; j < _ptb->n_elt[i]; j++) {
      int iproc = _ptb->dest_proc[++idx];
      int s_octet_elt = 0;
      int i_part_elt = 0;

      if (t_stride == PDM_STRIDE_CST) {
        s_octet_elt = cst_stride * (int) s_data;
        i_part_elt  = cst_stride * (int) s_data * j;
      }

      else if (t_stride == PDM_STRIDE_VAR) {
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

  PDM_MPI_Alltoallv_l(send_buffer,
                      n_send_buffer,
                      i_send_buffer,
                      PDM_MPI_BYTE,
                      recv_buffer,
                      n_recv_buffer,
                      i_recv_buffer,
                      PDM_MPI_BYTE,
                      _ptb->comm);

  free(send_buffer);
  free(n_send_buffer);
  free(i_send_buffer);
  free(n_recv_buffer);
  free(i_recv_buffer);

  unsigned char *_block_data = malloc(sizeof(unsigned char) * s_recv_buffer);
  *block_data = _block_data;
  *block_stride = NULL;
  int *i_recv_stride = NULL;
  int *i_block_stride = NULL;
  int s_block_data = ((int) sizeof(unsigned char) * s_recv_buffer) / (int) s_data;

  if (t_stride == PDM_STRIDE_VAR) {
    int* _block_stride = NULL;
    if(_ptb->tn_recv_data > 0){
      _block_stride = malloc(sizeof(int) * _ptb->tn_recv_data);
    }
    *block_stride = _block_stride;
    for (int i = 0; i < _ptb->tn_recv_data; i++) {
      _block_stride[i] = recv_stride[_ptb->order[i]];
    }

    /*
     * Compute index in data
     */

    i_recv_stride = malloc (sizeof(int) * (_ptb->tn_recv_data + 1));
    i_block_stride = malloc (sizeof(int) * (_ptb->tn_recv_data + 1));

    i_recv_stride[0] = 0;
    i_block_stride[0] = 0;
    for (int i = 0; i < _ptb->tn_recv_data; i++) {
      i_recv_stride[i+1]   = i_recv_stride[i] + recv_stride[i];
      i_block_stride[i+1] = i_block_stride[i] + _block_stride[i];
    }

    for (int i = 0; i < _ptb->tn_recv_data; i++) {
      i_recv_stride[i+1]   *= (int) s_data;
      i_block_stride[i+1] *= (int) s_data;
    }

    /*
     * Sort Buffer
     */

    for (int i = 0; i < _ptb->tn_recv_data; i++) {
      int old = _ptb->order[i];
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

    if (_ptb->t_post != PDM_PART_TO_BLOCK_POST_NOTHING) {

      int idx1 = 0;
      int idx2 = 0;

      if (_ptb->tn_recv_data == 1) {
        idx2 = i_block_stride[1];
      }

      for (int i = 1; i < _ptb->tn_recv_data; i++) {
        if (i == 1) {
          idx2 = i_block_stride[1];
        }
        if (_ptb->block_gnum[idx1] != _ptb->sorted_recv_gnum[i]) {
          idx1 += 1;
          _block_stride[idx1] = _block_stride[i];
          if (_ptb->t_post == PDM_PART_TO_BLOCK_POST_CLEANUP) {
            for (int k = i_block_stride[i]; k < i_block_stride[i+1]; k++) {
              _block_data[idx2++] = _block_data[k];
            }
          }
        }
        else {
          if (_ptb->t_post == PDM_PART_TO_BLOCK_POST_MERGE) {
            _block_stride[idx1] += _block_stride[i];
          }
        }
      }

      /* Cleanup */

      if (_ptb->t_post == PDM_PART_TO_BLOCK_POST_CLEANUP) {
        _block_data = realloc (_block_data, sizeof(unsigned char) * idx2);
        *block_data = _block_data;

        _block_stride = realloc (_block_stride, sizeof(int) * _ptb->n_elt_block);

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

    for (int i = 0; i < _ptb->tn_recv_data; i++) {
      int n_octet = cst_stride * (int) s_data;
      int old = _ptb->order[i];
      int idOld = old * n_octet;

      for (int k = i*n_octet; k < (i+1)*n_octet; k++) {
        _block_data[k] = recv_buffer[idOld++];
      }
    }

    /*
     * Post processing
     */

    if (_ptb->t_post != PDM_PART_TO_BLOCK_POST_NOTHING) {
      int idx2 = 0;
      int idx1 = 0;

      assert (_ptb->t_post != PDM_PART_TO_BLOCK_POST_MERGE);

      if (_ptb->tn_recv_data == 1) {
        idx2 =  cst_stride * (int) s_data;
      }

      for (int i = 1; i < _ptb->tn_recv_data; i++) {
        int n_octet = cst_stride * (int) s_data;
        if (i == 1) {
          idx2 = n_octet;
        }
        if (_ptb->block_gnum[idx1] != _ptb->sorted_recv_gnum[i]) {
          idx1 += 1;
          if (_ptb->t_post == PDM_PART_TO_BLOCK_POST_CLEANUP) {
            int idx3 = i * cst_stride * (int) s_data;
            for (int k = 0; k < n_octet; k++) {
              _block_data[idx2++] = _block_data[idx3++];
            }
          }
        }
      }

      if (_ptb->t_post == PDM_PART_TO_BLOCK_POST_CLEANUP) {
        _block_data = realloc (_block_data, sizeof(unsigned char) * idx2);
        *block_data = _block_data;
        s_block_data = idx2 / (int) s_data;
      }

    }
  }

  free (recv_buffer);
  return s_block_data;

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
  _pdm_part_to_block_t *_ptb = (_pdm_part_to_block_t *) ptb;

  if (_ptb->active_ranks != NULL) {
    free (_ptb->active_ranks);
    _ptb->active_ranks = NULL;
  }
  if (_ptb->dest_proc != NULL) {
    free (_ptb->dest_proc);
    _ptb->dest_proc = NULL;
  }
  if (_ptb->data_distrib_index != NULL) {
    free (_ptb->data_distrib_index);
    _ptb->data_distrib_index = NULL;
  }
  if (_ptb->i_send_data != NULL) {
    free (_ptb->i_send_data);
    _ptb->i_send_data = NULL;
  }
  if (_ptb->i_recv_data != NULL) {
    free (_ptb->i_recv_data);
    _ptb->i_recv_data = NULL;
  }
  if (_ptb->n_send_data != NULL) {
    free (_ptb->n_send_data);
    _ptb->n_send_data = NULL;
  }
  if (_ptb->n_recv_data != NULL) {
    free (_ptb->n_recv_data);
    _ptb->n_recv_data = NULL;
  }
  if (_ptb->sorted_recv_gnum != NULL) {
    free (_ptb->sorted_recv_gnum);
    _ptb->sorted_recv_gnum = NULL;
  }
  if (_ptb->order != NULL) {
    free (_ptb->order);
    _ptb->order = NULL;
  }

  if ((_ptb->t_post != PDM_PART_TO_BLOCK_POST_NOTHING) && (_ptb->block_gnum != NULL)) {
    free (_ptb->block_gnum);
    _ptb->block_gnum = NULL;
  }
  free (_ptb);
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
  _pdm_part_to_block_t *_ptb = (_pdm_part_to_block_t *) ptb;
  return _ptb->data_distrib_index;
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
  _pdm_part_to_block_t *_ptb = (_pdm_part_to_block_t *) ptb;
  return _ptb->dest_proc;
}

#undef _MIN
#undef _MAX

#ifdef __cplusplus
}
#endif /* __cplusplus */
