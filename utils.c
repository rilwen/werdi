#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "doublecomplex.h"
#include "clapack.h"

/* Put eigenvectors in matrix and eigenvalues in eigvals. */
int diagonalize_hermitian_matrix(doublecomplex * matrix, double * eigvals, int size)
{
    int info;
    zheev('V', 'U', size, &matrix[0], size, eigvals, &info);
    if (info != 0)
    {
        if (info < 0)
        {
            fprintf(stderr, "ZHEEV argument number %i has illegal value\n", -info);
        }
        else
        {
            fprintf(stderr, "The algorithm failed to converge: %i off-diagonal elements of an intermediate tridiagonal form failed to converge\n", info);
        }
    }
    return info;
}

void print_complex_matrix(FILE* out, doublecomplex* matrix, int dim)
{
	for (int row = 0; row < dim; ++row)
	{
		for (int col = 0; col < dim; ++col)
		{
			fprintf(out, "(%g, %g) ", matrix[col*dim + row].real, matrix[col*dim + row].imag);
		}
		fprintf(out, "\n");
	}
}

/* Multiply two complex numbers. */
doublecomplex cmplx_mul(doublecomplex a, doublecomplex b)
{
    doublecomplex m;
    m.real = a.real * b.real - a.imag * b.imag;
    m.imag = a.real * b.imag + a.imag * b.real;
    return m;
}

/* Complex conjugate a. */
doublecomplex cmplx_conj(doublecomplex a)
{
    doublecomplex m = {.real = a.real, .imag = - a.imag };
    return m;
}

/* Multiply a by hermitian matrix and put the result in b. */
void multiply_by_hermitian_matrix(const doublecomplex * matrix, const doublecomplex * a, int size, doublecomplex * b)
{    
    memset((void *) b, 0, size * sizeof(doublecomplex));
    for (int i = 0; i < size; i++)
    {
        doublecomplex ai = a[i];
        /* Now use the fact that matrix is hermitian and matrix[i][i].imag == 0 */
        b[i].real += ai.real * matrix[i * size + i].real;
        b[i].imag += ai.imag * matrix[i * size + i].real;
        for (int j = 0; j < i; j++)
        {
            doublecomplex mij = matrix[i * size + j];
            b[j].real += mij.real * ai.real - mij.imag * ai.imag;
            b[j].imag += mij.real * ai.imag + mij.imag * ai.real;
            /* Now use the fact that matrix is hermitian and matrix[i][j] = complex_conjugate(matrix[j][i]). */
            b[i].real += mij.real * a[j].real + mij.imag * a[j].imag;
            b[i].imag += mij.real * a[j].imag - mij.imag * a[j].real;
        }
    }
}

/* Scalar product of complex vectors <l|r> */
doublecomplex scalar_product(const doublecomplex * l, const doublecomplex * r, int size)
{
    doublecomplex sp = {.real = 0, .imag = 0};
    for (int n = 0; n < size; n++)
    {
        sp.real += l[n].real * r[n].real + l[n].imag * r[n].imag;
        sp.imag += l[n].real * r[n].imag - l[n].imag * r[n].real;
    }
    return sp;
}

/* Square of the l^2 norm. */
double norm2(const doublecomplex * v, int size)
{
    double nrm = 0;
    for (int n = 0; n < size; n++)
    {
        nrm += v[n].real * v[n].real + v[n].imag * v[n].imag;
    }
    return nrm;
}

/* l^2 norm of vector. */
double norm(const doublecomplex * v, int size)
{
    return sqrt(norm2(v, size));
}

/* Matrix element of hermitian M between vectors l and r <l|M|r> */
doublecomplex hermitian_matrix_element(const doublecomplex * matrix, const doublecomplex * l, const doublecomplex * r, const int size)
{
    doublecomplex sp = {.real = 0, .imag = 0};
    for (int i = size - 1; i >= 0; --i)
    {
        const double diagi = matrix[i * size + i].real;
        const doublecomplex li = l[i];
        const doublecomplex ri = r[i];
        sp.real += diagi * (li.real * ri.real + li.imag * ri.imag);
        sp.imag += diagi * (li.real * ri.imag - li.imag * ri.real);
        for (int j = i - 1; j >= 0; --j)
        {
            const doublecomplex mij = matrix[i * size + j]; // j-th row, i-th column
	    doublecomplex tmp = r[j];
            doublecomplex lr;
            lr.real = li.real * tmp.real + li.imag * tmp.imag;
            lr.imag = li.real * tmp.imag - li.imag * tmp.real;
            sp.real += mij.real * lr.real + mij.imag * lr.imag;
            sp.imag += mij.real * lr.imag - mij.imag * lr.real;
            /* Now make use of the hermicity. */
	    tmp = l[j];
            lr.real = tmp.real * ri.real + tmp.imag * ri.imag;
            lr.imag = tmp.real * ri.imag - tmp.imag * ri.real;
            sp.real += mij.real * lr.real - mij.imag * lr.imag;
            sp.imag += mij.real * lr.imag + mij.imag * lr.real;
        }
    }
    return sp;
}

/* Matrix element of complex M between vectors l and r <l|M|r> */
doublecomplex complex_matrix_element(const doublecomplex * matrix, const doublecomplex * l, const doublecomplex * r, int dim)
{
    doublecomplex sp = {.real = 0, .imag = 0};

    for (int i = dim - 1; i >= 0; --i)
    {
	    doublecomplex li = l[i];
	    for (int j = dim - 1; j >= 0; --j)
	    {
		    doublecomplex lr;
		    lr.real = li.real * r[j].real + li.imag * r[j].imag;
		    lr.imag = li.real * r[j].imag - li.imag * r[j].real;
		    doublecomplex mij = matrix[j * dim + i]; // i-th row, j-th column
		    sp.real += mij.real * lr.real - mij.imag * lr.imag;
		    sp.imag += mij.real * lr.imag + mij.imag * lr.real;
	    }
    }
    return sp;
}

/* Mean value of hermitian M for vector v: <v|M|v> */
double hermitian_matrix_mean_value(const doublecomplex * matrix, const doublecomplex * v, int size)
{
    double mv = 0;
    for (int i = 0; i < size; ++i)
    {
        doublecomplex vi = v[i];
        mv += matrix[i * size + i].real * (vi.real * vi.real + vi.imag * vi.imag);
        for (int j = 0; j < i; j++)
        {
            mv += 2 * (matrix[i * size + j].real * (vi.real * v[j].real + vi.imag * v[j].imag) + matrix[i * size + j].imag * (vi.real * v[j].imag - vi.imag * v[j].real));
        }
    }
    return mv;
}

/* A = A + c*B */
void add_matrix(doublecomplex* A, doublecomplex* B, double c, int dim)
{
	int nrelem = dim * dim;
	for (int i = 0; i < nrelem; i++)
	{
		A[i].real += c * B[i].real;
		A[i].imag += c * B[i].imag;
	}
}

/* A = A * c */
void multiply_matrix(doublecomplex* A, double c, int dim)
{
	int nrelem = dim * dim;
	for (int i = 0; i < nrelem; i++)
	{
		A[i].real *= c;
		A[i].imag *= c;
	}
}

/* Allocate memory, abort when unsucessful. */
void* alloc_memory(size_t size)
{
	if (size == 0)
	{
		fprintf(stderr, "Warning: zero memory allocated\n");
		return NULL;
	}
	void* ptr = malloc(size);
	if (ptr == NULL)
	{
		fprintf(stderr, "Memory allocation error while trying to reserve %ld bytes of memory\n", size);
		abort();
	}
	return ptr;
}

/* Create a complex matrix of dimension dim and set it to zero */
doublecomplex* initialize_cmplx_matrix(int dim)
{
        size_t size = dim * dim * sizeof(doublecomplex);
        doublecomplex* sz = alloc_memory(size);
        memset(sz, 0, size);
        return sz;
}

void whinge_and_die(const char* whinge)
{
	fprintf(stderr, "%s\n", whinge);
	abort();
}

void free_if_not_null(void* ptr)
{
	if (ptr != NULL)
	{
		free(ptr);
	}
}

void* copy_to_dynamic_array(const void* static_data, size_t size)
{
	void* dyn_arr = alloc_memory(size);
	memcpy(dyn_arr, static_data, size);
	return dyn_arr;
}

void invert_complex_matrix(doublecomplex* matrix, int dim)
{
	int info;
	int* pivots = alloc_memory(sizeof(int) * dim);
	zgetrf(dim, dim, matrix, dim, pivots, &info);
	if (info != 0)
	{
		fprintf(stderr, "ZGETRF error, INFO == %d\n", info);
	}
	zgetri(dim, matrix, dim, pivots, &info);
	if (info != 0)
	{
		fprintf(stderr, "ZGETRI error, INFO == %d\n", info);
	}
	free(pivots);
}

void invert_double_matrix(double* matrix, int dim)
{
	int info;
	int* pivots = alloc_memory(sizeof(int) * dim);
	dgetrf(dim, dim, matrix, dim, pivots, &info);
	if (info != 0)
	{
		fprintf(stderr, "DGETRF error, INFO == %d\n", info);
	}
	dgetri(dim, matrix, dim, pivots, &info);
	if (info != 0)
	{
		fprintf(stderr, "DGETRI error, INFO == %d\n", info);
	}
	free(pivots);
}

void matrix_product(const doublecomplex* a, const doublecomplex* b, doublecomplex* c, int rowsa, int colsa, int colsb)
{
	int rowsc = rowsa;
	int colsc = colsb;
	int rowsb = colsa;

	for (int rc = 0; rc < rowsc; rc++)
	{
		for (int cc = 0; cc < colsc; cc++)
		{
			c[cc*rowsc + rc] = cmplx(0, 0);
			for (int i = 0; i < colsa; i++)
			{
				c[cc*rowsc + rc].real += a[i*rowsa + rc].real * b[cc*rowsb + i].real - a[i*rowsa + rc].imag * b[cc*rowsb + i].imag;
				c[cc*rowsc + rc].imag += a[i*rowsa + rc].real * b[cc*rowsb + i].imag + a[i*rowsa + rc].imag * b[cc*rowsb + i].real;
				
			}
		}
	}
}

void square_matrix_product(const doublecomplex* a, const doublecomplex* b, doublecomplex* c, int dim)
{
	matrix_product(a, b, c, dim, dim, dim);
}

void hermitian_conjugate_matrix(doublecomplex* matrix, int dim)
{
	for (int i = 0; i < dim; i++)
	{
		matrix[i*dim + i].imag *= -1;
		for (int j = 0; j < i; j++)
		{
			doublecomplex tmp = matrix[i*dim + j];
			matrix[i*dim + j].real = matrix[j*dim + i].real;
			matrix[i*dim + j].imag = - matrix[j*dim + i].imag;

			matrix[j*dim + i].real = tmp.real;
			matrix[j*dim + i].imag = - tmp.imag;
		}
	}
}
