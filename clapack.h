#ifndef _CLAPACK_H
#define _CLAPACK_H
#include "doublecomplex.h"

void zgetrf(int m, int n, doublecomplex* a, int lda, int* ipiv, int* info);
void zgetri(int n, doublecomplex* a, int lda, int* ipiv, int* info);
void zheev(char jobz, char uplo, int n, doublecomplex *a, int lda, double *w, int *info);
void dgetrf(int m, int n, double* a, int lda, int* ipiv, int* info);
void dgetri(int n, double* a, int lda, int* ipiv, int* info);

#endif
