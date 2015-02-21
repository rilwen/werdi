#include <stdio.h>
#include <math.h>
#include "testing.h"

void test_double(double expected, double actual, double tolerance, const char* testname)
{
	if (fabs(actual - expected) > tolerance)
	{
		printf("Test %s FAILED:\n", testname);
		printf("Expected result %g with tolerance %g, but got %g instead\n", expected, tolerance, actual);
	}
	else
	{
		printf("Test %s SUCCEEDED\n", testname);
	}
}

void test_int(int expected, int actual, const char* testname)
{
	if (actual != expected)
	{
		printf("Test %s FAILED:\n", testname);
		printf("Expected result %d, but got %d instead\n", expected, actual);
	}
	else
	{
		printf("Test %s SUCCEEDED\n", testname);
	}
}

void test_doublecomplex(doublecomplex expected, doublecomplex actual, double tolerance, const char* testname)
{
	if (fabs(actual.real - expected.real) > tolerance || fabs(actual.imag - expected.imag) > tolerance)
	{
		printf("Test %s FAILED:\n", testname);
		printf("Expected result (%g,%g) with tolerance %g, but got (%g,%g) instead\n", expected.real, expected.imag, tolerance, actual.real, actual.imag);
	}
	else
	{
		printf("Test %s SUCCEEDED\n", testname);
	}
}

void test_success(int result, const char* testname)
{
	if (result != 0)
	{
		printf("Test %s SUCCEEDED\n", testname);
	}
	else
	{
		printf("Test %s FAILED\n", testname);
	}
}

int is_hermitian(doublecomplex * matrix, int dim)
{
	const double eps = 1E-14;
	
	for (int i = 0; i < dim; i++)
	{
		for (int j = 0; j <= i; j++)
		{
			double dr = fabs(matrix[i * dim + j].real - matrix[j * dim + i].real);
			double di = fabs(matrix[i * dim + j].imag + matrix[j * dim + i].imag);
			if (dr > eps || di > eps)
			{
				printf("Non-hermicity detected for element (%i, %i): dr == %g, di == %g\n", j, i, dr, di);
				return 0;
			}
		}
	}

	return 1;
}
