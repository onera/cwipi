/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_unique.h"
#include "pdm_quick_sort.h"
#include "pdm_radix_sort.h"

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

/*=============================================================================
 * Static function definitions
 *============================================================================*/


/**
 *
 * \brief Unique
 *
 * \param [inout]   a     Array to sort
 * \param [in]      l     First element
 * \param [in]      r     Last  element
 *
 */
int
PDM_inplace_unique_long
(
 PDM_g_num_t a[],
 int l,
 int r
)
{
  // PDM_quick_sort_long(a, l, r); /* Less optimal than PDM_sort_long */
  int array_size = r - l + 1;
  // printf("PDM_inplace_unique_long::array_size::%d\n", array_size);
  PDM_sort_long(&a[l], NULL, array_size);

  int new_size  = 1;
  int idx_write = l;
  PDM_g_num_t last_value = a[l];
  a[idx_write++] = last_value;
  for (int idx = l+1; idx <= r; idx++) {
    if(last_value != a[idx]){
      last_value = a[idx];
      a[idx_write++] = a[idx];
      new_size++;
    }
  }

  return new_size;
}


/**
 *
 * \brief Unique
 *
 * \param [inout]   a             Array to sort
 * \param [inout]   unique_order  Unique index in old numbering
 * \param [in]      l             First element
 * \param [in]      r             Last  element
 *
 */
int
PDM_inplace_unique_long2
(
 PDM_g_num_t a[],
 int unique_order[],
 int l,
 int r
)
{
  // printf("PDM_inplace_unique_long::a:: ");
  // for(int i = l; i < r; ++i){
  //   printf("%d ", a[i]);
  // }
  // printf("\n");


  // PDM_quick_sort_long(a, l, r); /* Less optimal than PDM_sort_long */
  int array_size = r - l + 1;
  // printf("PDM_inplace_unique_long::array_size::%d\n", array_size);
  int* order = (int *) malloc( (array_size) * sizeof(int));

  for(int i = 0; i < array_size; ++i){
    order[i] = i;
  }
  PDM_radix_sort_long(&a[l], order, array_size);
  // PDM_sort_long(&a[l], order, array_size);


  int first = a[l];
  for(int i = 1; i < array_size; ++i ) {
    if(a[i] < first){
      printf("The list is not sorted : a[%d] = %d > a[%d] = "PDM_FMT_G_NUM" \n", i-1, first, i, a[i]);
      // printf("Problem with list size : %d\n", itest);
      abort();
    }
    first = a[i];
  }

  // for(int i = 0; i < array_size; ++i){
  //   order[i] = i;
  // }
  // PDM_sort_long(&a[l], order, array_size);
  // PDM_quick_sort_long2(&a[l], 0, array_size, order);

  // printf("PDM_inplace_unique_long::a::sort:: ");
  // for(int i = l; i < r; ++i){
  //   printf("%d ", a[i]);
  // }
  // printf("\n");

  // printf("PDM_inplace_unique_long::a::order:: ");
  // for(int i = 0; i < array_size; ++i){
  //   printf("%d ", order[i]);
  //   unique_order[i] = -1;
  // }
  // printf("\n");

  int new_size  = 1;
  int idx_write = l;
  PDM_g_num_t last_value = a[l];
  int idx_save = l;
  unique_order[order[0]] = idx_save;
  a[idx_write++] = last_value;
  for (int idx = l+1; idx <= r; idx++) {
    if(last_value != a[idx]){
      last_value = a[idx];
      // printf(" order[%d] = %d\n", idx-l, order[idx-l]);
      idx_save = idx_write;
      unique_order[order[idx-l]] = idx_save;
      a[idx_write++] = a[idx];
      new_size++;
    }
    unique_order[order[idx-l]] = idx_save;
  }

  // printf("PDM_inplace_unique_long::a::unique_order:: ");
  // for(int i = 0; i < array_size; ++i){
  //   printf("%d ", unique_order[i]);
  // }
  // printf("\n");

  free(order);

  return new_size;
}


#ifdef __cplusplus
}
#endif /* __cplusplus */
