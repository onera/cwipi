#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>


/* DSYEV prototype */
extern void dsyev( char* jobz, char* uplo, int* n, double* a, int* lda,
                   double* w, double* work, int* lwork, int* info );

static void matmul (double A[9], double x[3], double b[3]) {
  for (int i = 0; i < 3; i++) {
    b[i] = 0.;

    for (int j = 0; j < 3; j++) {
      b[i] += A[3*i + j] * x[j];
    }
  }
}

int main() {


  double A[9] = {1., 2., 3.,
                 2., 4., 1.,
                 3., 1., 5.};

  double B[9];
  memcpy(B, A, sizeof(double) * 9);

  double eig_val[3];
  double work[8];
  int info;

  int n = 3;
  int lwork = 3*n - 1;

  dsyev ("V",
         "U",
         &n,
         A,
         &n,
         eig_val,
         work,
         &lwork,
         &info);

  printf("info = %d\n", info);

  double *eig_vec = A;
  if (info == 0) {
    printf("eig_val = %f %f %f\n", eig_val[0], eig_val[1], eig_val[2]);
    printf("eig_vec =\n%f %f %f\n%f %f %f\n%f %f %f\n",
           eig_vec[0], eig_vec[1], eig_vec[2],
           eig_vec[3], eig_vec[4], eig_vec[5],
           eig_vec[6], eig_vec[7], eig_vec[8]);
  }


  for (int i = 0; i < 3; i++) {
    //double x[3] = {eig_vec[i], eig_vec[i+3], eig_vec[i+6]};
    double x[3] = {eig_vec[3*i], eig_vec[3*i+1], eig_vec[3*i+2]};
    double c[3];
    matmul(B, x, c);

    for (int j = 0; j < 3; j++) {
      c[j] -= eig_val[i]*x[j];
    }

    printf("check %d : %f %f %f\n", i, c[0], c[1], c[2]);
  }


  return 0;
}
