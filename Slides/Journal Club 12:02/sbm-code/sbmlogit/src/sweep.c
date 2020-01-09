#include <R.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>
#include <R_ext/Rdynload.h>

static int sweepaux (int n, int m, double *A, int k) {
  int i;
  double D = A[k * (n + 1)];
  if (D == 0) return 1; /* error */
  D = 1.0 / D;
  F77_NAME(dscal)(&m, &D, A + k, &n); /* A[k,] = A[k,] * D */
  for (i = 0; i < n; i++) {
    if (i != k) {
      double B = -A[i + k * n]; /* -A[i,k] */
      /* A[i,] = A[i,] + B * A[k,]: */
      F77_NAME(daxpy)(&m, &B, A + k, &n, A + i, &n);
      A[i + k * n] = B * D;
    }
  }
  A[k * (n + 1)] = D; /* A[k,k] = D */
  return 0;
}

SEXP sweep (SEXP sn, SEXP sm, SEXP sa, SEXP sind) {
  int n = INTEGER(sn)[0];
  int m = INTEGER(sm)[0];
  double *A = REAL(sa);
  int *ind = INTEGER(sind);
  int i, k = length(sind);
  for (i = 0; i < k; i++) sweepaux(n, m, A, ind[i]);
  return sa;
}

