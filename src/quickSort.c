
#include <bft_mem.h>
#include "quickSort.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/* Fonctions privees */

int partitionner(int tableau[], int p, int r, int* index) {
  int pivot = tableau[p], i = p-1, j = r+1;
  int temp;
  while(1) {
    do
      j--;
    while(tableau[j] > pivot);
    do
      i++;
    while(tableau[i] < pivot);
    if(i<j) {
      temp = tableau[i];
      tableau[i] = tableau[j];
      tableau[j] = temp;
      if (index != NULL) {
        temp = index[i];
        index[i] = index[j];
        index[j] = temp;
      }
    }
    else
      return j;
  }
  return j;
}


/* Fonctions publiques */

void quickSort(int tableau[], int p, int r, int* index) {
  int q;
  int i;

  if (index != NULL) {
    for (i= 0; i <r+1; i++)
      index[i] = i+1;
  }

  if(p<r) {
    q = partitionner(tableau, p, r, index);
    quickSort(tableau, p, q, index);
    quickSort(tableau, q+1, r, index);
  }
  return;
}
 

#ifdef __cplusplus
}
#endif /* __cplusplus */

