#ifndef _UTILS_H
#define _UTILS_H

#include <stdlib.h>
#include <stdio.h>
#include "doublecomplex.h"

/* Put eigenvectors in matrix and eigenvalues in eigvals. */
int diagonalize_hermitian_matrix(doublecomplex * matrix, double * eigvals, int size);

/* Multiply two complex numbers. */
doublecomplex cmplx_mul(doublecomplex a, doublecomplex b);

/* Complex conjugate a. */
doublecomplex cmplx_conj(doublecomplex a);

/* Multiply a by hermitian matrix and put the result in b. */
void multiply_by_hermitian_matrix(const doublecomplex * matrix, const doublecomplex * a, int size, doublecomplex * b);

/* Scalar product of complex vectors <l|r> */
doublecomplex scalar_product(const doublecomplex * l, const doublecomplex * r, int size);

/* Square of the l^2 norm. */
double norm2(const doublecomplex * v, int size);

/* l^2 norm of vector. */
double norm(const doublecomplex * v, int size);

/* Matrix element of hermitian M between vectors l and r <l|M|r> */
doublecomplex hermitian_matrix_element(const doublecomplex * matrix, const doublecomplex * l, const doublecomplex * r, const int size);

/* Matrix element of complex M between vectors l and r <l|M|r> */
doublecomplex complex_matrix_element(const doublecomplex * matrix, const doublecomplex * l, const doublecomplex * r, int dim);

/* Mean value of hermitian M for vector v: <v|M|v> */
double hermitian_matrix_mean_value(const doublecomplex * matrix, const doublecomplex * v, int size);

/* A = A + c*B */
void add_matrix(doublecomplex* A, doublecomplex* B, double c, int dim);

/* A = A * c */
void multiply_matrix(doublecomplex* A, double c, int dim);

/* Allocate memory, abort when unsucessful. */
void* alloc_memory(size_t size);

/* Create a complex matrix of dimension dim and set it to zero */
doublecomplex* initialize_cmplx_matrix(int dim);

/* Display an error message using British spelling. And die. */
void whinge_and_die(const char* whinge);

/* Free the pointer if it is not null */
void free_if_not_null(void* ptr);

/* Copy constents of a static array to dynamic array, freshly allocated */
void* copy_to_dynamic_array(const void* static_data, size_t size);

/* Invert a complex matrix in-place */
void invert_complex_matrix(doublecomplex* matrix, int dim);

/* Invert a double matrix in-place */
void invert_double_matrix(double* matrix, int dim);

/* Calculate C = A*B */
void square_matrix_product(const doublecomplex* a, const doublecomplex* b, doublecomplex* res, int dim);
void matrix_product(const doublecomplex* a, const doublecomplex* b, doublecomplex* res, int rowsa, int colsa, int colsb);

void hermitian_conjugate_matrix(doublecomplex* matrix, int dim);

void print_complex_matrix(FILE* out, doublecomplex* matrix, int dim);

#endif /* _UTILS_H */
