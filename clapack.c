#include "clapack.h"
#include "doublecomplex.h"
#include <stdlib.h>
#include <stdio.h>

/* Declarations of Fortran LAPACK functions */
extern void zheev_(char*, char*, int*, doublecomplex*, int*, double*, doublecomplex*, int*, double*, int*);
extern void zgetrf_(int*, int*, doublecomplex*, int*, int*, int*);
extern void dgetrf_(int*, int*, double*, int*, int*, int*);
extern void zgetri_(int*, doublecomplex*, int*, int*, doublecomplex*, int*, int*);
extern void dgetri_(int*, double*, int*, int*, double*, int*, int*);

void zgetrf(int m, int n, doublecomplex* a, int lda, int* ipiv, int* info)
{
	zgetrf_(&m, &n, a, &lda, ipiv, info);
}

void dgetrf(int m, int n, double* a, int lda, int* ipiv, int* info)
{
	dgetrf_(&m, &n, a, &lda, ipiv, info);
}

void zgetri(int n, doublecomplex* a, int lda, int* ipiv, int* info)
{
	int lwork = -1;
	doublecomplex* work = malloc(sizeof(doublecomplex));
	if (NULL == work)
	{
		abort();
	}
	zgetri_(&n, a, &lda, ipiv, work, &lwork, info);
	lwork = (int) work[0].real;
	work = realloc(work, sizeof(doublecomplex) * lwork);
	if (NULL == work)
	{
		abort();
	}
	zgetri_(&n, a, &lda, ipiv, work, &lwork, info);
	free(work);
}

void dgetri(int n, double* a, int lda, int* ipiv, int* info)
{
	int lwork = -1;
	double* work = malloc(sizeof(double));
	if (NULL == work)
	{
		abort();
	}
	dgetri_(&n, a, &lda, ipiv, work, &lwork, info);
	lwork = (int) work[0];
	work = realloc(work, sizeof(double) * lwork);
	if (NULL == work)
	{
		abort();
	}
	dgetri_(&n, a, &lda, ipiv, work, &lwork, info);
	free(work);
}

void zheev(char jobz, char uplo, int n, doublecomplex *a, int lda, double *w, int *info)
{
	doublecomplex *work;
	double *rwork;
	int lwork = -1;

	work = malloc(sizeof(doublecomplex));
	if (NULL == work) {
		abort();
	}
	if (3*n - 2 >= 1) {
		rwork = malloc(sizeof(double)*(3*n - 2));
	}
	else {
		rwork = malloc(sizeof(double));
	}
	if (NULL == rwork) {
		abort();
	}

	zheev_(&jobz, &uplo, &n, a, &lda, w, work, &lwork, rwork, info);
	lwork = (int) work[0].real;
	work = realloc(work, sizeof(doublecomplex) * lwork);
	if (NULL == work) {
		abort();
	}
	zheev_(&jobz, &uplo, &n, a, &lda, w, work, &lwork, rwork, info);
	free(work);
	free(rwork);
}

