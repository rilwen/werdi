#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "utils.h"
#include "testing.h"
#include "kane.h"
#include "doublecomplex.h"


void test_hermitian_matrix_mean_value(void)
{
	doublecomplex w[2][2];
	
	w[0][0] = cmplx(1, 0);
	w[0][1] = cmplx(0, 0);
	w[1][0] = cmplx(0, 0);
	w[1][1] = cmplx(2, 0);

	doublecomplex v[2];
	v[0] = cmplx(1, 0);
	v[1] = cmplx(1, 0);

	double vwv = hermitian_matrix_mean_value(w, v, 2);
	test_double(3, vwv, 1E-16, "HMMV_DIAG_TEST");

	w[0][1] = cmplx(1, 1);
	w[1][0] = cmplx(1, -1);
	v[1] = cmplx(2, 0);
	vwv = hermitian_matrix_mean_value(w, v, 2);
	test_double(13, vwv, 1E-16, "HMMV_FULL_ONES_TEST");

	doublecomplex m[5][5];
	doublecomplex ev[5][5];
	const double dv = 0.01;;
	double val = 0.1;
	for (int i = 0; i < 5; i++)
	{
		for (int j = 0; j <= i; j++)
		{
			m[i][j].real = val;
			m[j][i].real = val;
			if (i != j)
			{
				m[i][j].imag = val;
				m[j][i].imag = -val;
			}
			else
			{
				m[i][j].imag = 0;
				m[j][i].imag = 0;
			}
			val += dv;
		}
	}
	memcpy(ev, m, 25 * sizeof(doublecomplex));
	double en[5];
	diagonalize_hermitian_matrix(ev, en, 5);
	for (int i = 0; i < 5; i++)
	{
		double mv = hermitian_matrix_mean_value(m, ev[i], 5);
		test_double(en[i], mv, 1E-8, "HMMV_EIGVECS_TEST");
	}
}

void test_hermitian_matrix_element(void)
{
	doublecomplex m[5][5];
	doublecomplex ev[5][5];
	double val = 0.1;
	const double dv = 0.01;
	
	for (int i = 0; i < 5; i++)
	{
		for (int j = 0; j <= i; j++)
		{
			m[i][j].real = val;
			m[j][i].real = val;
			if (i != j)
			{
				m[i][j].imag = val;
				m[j][i].imag = -val;
			}
			else
			{
				m[i][j].imag = 0;
				m[j][i].imag = 0;
			}
			val += dv;
		}
	}
	memcpy(ev, m, 25 * sizeof(doublecomplex));
	double en[5];
	diagonalize_hermitian_matrix(ev, en, 5);

	doublecomplex a[5];
	doublecomplex b[5];
	
	char test_name[64];
	
	for (int i = 0; i < 5; i++)
	{
		memset(a, 0, sizeof(doublecomplex) * 5);
		a[i].real = 1;
		a[i].imag = 0;
		snprintf(test_name, 64, "HMEL_TEST_DIAG_REAL_%d", i);
		doublecomplex mel = hermitian_matrix_element(m, a, a, 5);
		test_double(hermitian_matrix_mean_value(m, a, 5), mel.real, 1E-16, test_name);
		snprintf(test_name, 64, "HMEL_TEST_DIAG_IMAG_%d", i);
		test_double(0, mel.imag, 1E-16, test_name);
		
		mel = hermitian_matrix_element(m, ev[i], ev[i], 5);
		snprintf(test_name, 64, "HMEL_TEST_DIAG_EIG_REAL_%d", i);
		test_double(en[i], mel.real, 1E-8, test_name);
		snprintf(test_name, 64, "HMEL_TEST_DIAG_EIG_IMAG_%d", i);
		test_double(0, mel.imag, 1E-8, test_name);
		for (int j = 0; j < 5; j++)
		{
			memset(b, 0, sizeof(doublecomplex) * 5);
			b[j].real = 1;
			b[j].imag = 0;
			doublecomplex mel2 = hermitian_matrix_element(m, a, b, 5);
			snprintf(test_name, 64, "HMEL_TEST_REAL_%d_%d", i, j);
			test_double(m[j][i].real, mel2.real, 1E-16, test_name);
			snprintf(test_name, 64, "HMEL_TEST_IMAG_%d_%d", i, j);
			test_double(m[j][i].imag, mel2.imag, 1E-16, test_name);

			if (i != j)
			{
				mel2 = hermitian_matrix_element(m, ev[i], ev[j], 5);
				snprintf(test_name, 64, "HMEL_TEST_OFFDIAG_EIG_REAL_%d_%d", i, j);
				test_double(0, mel2.real, 1E-8, test_name);
				snprintf(test_name, 64, "HMEL_TEST_OFFDIAG_EIG_IMAG_%d_%d", i, j);
				test_double(0, mel2.imag, 1E-8, test_name);
			}
		}
	}
}

void setup_canonical_vector(doublecomplex v[], int length, int pos)
{
	memset(v, 0, sizeof(doublecomplex) * length);
	v[pos].real = 1;
	v[pos].imag = 0;
}

void test_complex_matrix_element(void)
{
	doublecomplex l[2];
	doublecomplex r[2];

	doublecomplex m[2][2];
	m[0][0].real = 0;
	m[0][0].imag = 0;
	m[0][0].real = 0;
	m[0][1].imag = 1;
	m[1][0].real = 1;
	m[1][0].imag = 0;
	m[1][1].real = 1;
	m[1][1].imag = 1;
	char test_name[128];

	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			setup_canonical_vector(l, 2, i);
			setup_canonical_vector(r, 2, j);
			doublecomplex v = complex_matrix_element(m, l, r, 2);
			snprintf(test_name, 128, "CMEL_TEST_%d_%d_REAL", i, j);
			test_double(m[j][i].real, v.real, 1E-10, test_name);
			snprintf(test_name, 128, "CMEL_TEST_%d_%d_IMAG", i, j);
			test_double(m[j][i].imag, v.imag, 1E-10, test_name);
		}
	}
}

void test_scalar_product_and_norm(void)
{
	doublecomplex l[3];
	doublecomplex r[3];

	l[0] = cmplx(0, 1);
	r[0] = cmplx(0, -1);
	l[1] = cmplx(1, 0);
	r[1] = cmplx(-1, 0);
	l[2] = cmplx(1, -1);
	r[2] = cmplx(0, 1);

	doublecomplex lr = scalar_product(l, r, 3);
	doublecomplex rl = scalar_product(r, l, 3);

	test_double(-3, lr.real, 1E-16, "SP_TEST_LR_REAL");
	test_double(1, lr.imag, 1E-16, "SP_TEST_LR_IMAG");
	test_double(0, lr.real - rl.real, 1E-16, "SP_TEST_CONJ_REAL");
	test_double(0, lr.imag + rl.imag, 1E-16, "SP_TEST_CONJ_IMAG");
	
	doublecomplex ll = scalar_product(l, l, 3);
	test_double(4, ll.real, 1E-16, "SP_TEST_LNORM_REAL");
	test_double(0, ll.imag, 1E-16, "SP_TEST_LNORM_IMAG");
	
	doublecomplex rr = scalar_product(r, r, 3);
	test_double(3, rr.real, 1E-16, "SP_TEST_RNORM_REAL");
	test_double(0, rr.imag, 1E-16, "SP_TEST_RNORM_IMAG");

	double nl = norm(l, 3);
	test_double(4, nl*nl, 1E-16, "NRM_TEST_L");
	double nr = norm(r, 3);
	test_double(3, nr*nr, 1E-12, "NRM_TEST_R");
}

void test_multiply_by_hermitian_matrix(void)
{
	doublecomplex a[3], b[3], m[3][3];
	double val = 0.1, dv = 0.01;
	
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j <= i; j++)
		{
			m[i][j].real = val;
			m[j][i].real = val;
			if (i != j)
			{
				m[i][j].imag = val;
				m[j][i].imag = -val;
			}
			else
			{
				m[i][j].imag = 0;
				m[j][i].imag = 0;
			}
			val += dv;
		}
	}

	for (int i = 0; i < 3; i++)
	{
		memset(a, 0, sizeof(doublecomplex) * 3);
		a[i] = cmplx(1, 0);
		multiply_by_hermitian_matrix(m, a, 3, b);
		for (int j = 0; j < 3; j++)
		{
			test_doublecomplex(m[i][j], b[j], 1E-16, "MBHM_TEST_ONE");
		}

		a[i] = cmplx(0, 1);
		multiply_by_hermitian_matrix(m, a, 3, b);
		for (int j = 0; j < 3; j++)
		{
			doublecomplex t;
			t.real = - m[i][j].imag;
			t.imag = m[i][j].real;
			test_doublecomplex(t, b[j], 1E-16, "MBHM_TEST_TWO");
		}
	}
}

void test_copy_to_dynamic_array(void)
{
	doublecomplex array[10];
	for (int i = 0; i < 10; i++)
	{
		array[i].real = i;
		array[i].imag = 0;
	}
	doublecomplex* copy = copy_to_dynamic_array(array, 10 * sizeof(doublecomplex));
	for (int i = 0; i < 10; i++)
	{
		test_doublecomplex(array[i], copy[i], 1E-16, "TEST_COPY");
	}
}

void test_diagonalize_hermitian_matrix(void)
{
	doublecomplex m[3][3];

	// Assume [column][row] indexing.
	m[0][0] = cmplx(1, 0);
	m[0][1] = cmplx(0, 1);
	m[0][2] = cmplx(0, 0);

	m[1][0] = cmplx(0, -1);
	m[1][1] = cmplx(2, 0);
	m[1][2] = cmplx(1, -1);

	m[2][0] = cmplx(0, 0);
	m[2][1] = cmplx(1, 1);
	m[2][2] = cmplx(-1, 0);

	double v[3];

	diagonalize_hermitian_matrix(m, v, 3);

	test_double(-1.618033988749895, v[0], 1E-12, "DIAG_V0");
	test_double(0.618033988749895, v[1], 1E-12, "DIAG_V1");
	test_double(3, v[2], 1E-12, "DIAG_V2");

	test_doublecomplex(cmplx(0.106913483023849, -0.106913483023849), m[0][0], 1E-12, "DIAG_VEC00");
	test_doublecomplex(cmplx(-0.279903132412072, -0.279903132412072), m[0][1], 1E-12, "DIAG_VEC01");
	test_doublecomplex(cmplx(0.905785563600590, 0), m[0][2], 1E-12, "DIAG_VEC02");
	test_doublecomplex(cmplx(-0.630603216165774, 0.630603216165774), m[1][0], 1E-12, "DIAG_VEC10");
	test_doublecomplex(cmplx(0.240868995160328, 0.240868995160329), m[1][1], 1E-12, "DIAG_VEC11");
	test_doublecomplex(cmplx(0.297730451690234, 0), m[1][2], 1E-12, "DIAG_VEC12");
	test_doublecomplex(cmplx(0.301511344577764, -0.301511344577764), m[2][0], 1E-12, "DIAG_VEC20");
	test_doublecomplex(cmplx(0.603022689155527, 0.603022689155527), m[2][1], 1E-12, "DIAG_VEC21");
	test_doublecomplex(cmplx(0.301511344577764, 0), m[2][2], 1E-12, "DIAG_VEC22");
}

int main(void)
{
	test_hermitian_matrix_mean_value();

	test_hermitian_matrix_element();

	test_scalar_product_and_norm();

	test_multiply_by_hermitian_matrix();

	test_complex_matrix_element();

	test_copy_to_dynamic_array();

	test_diagonalize_hermitian_matrix();

	return 0;
}

