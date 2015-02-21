#include <math.h>
#include "lattice.h"
#include "testing.h"

void test_reciprocal_vectors()
{
	double a[3][3], b[3][3];

	const double a0 = 0.5633;

	a[0][0] = a0 / 2;
	a[0][1] = a0 / 2;
	a[0][2] = 0;
	a[1][0] = a0 / 2;
	a[1][1] = 0;
	a[1][2] = a0 / 2;
	a[2][0] = 0;
	a[2][1] = a0 / 2;
	a[2][2] = a0 / 2;

	calculate_reciprocal_lattice_vectors(a, b);

	const double q = 11.15424;
	const double tol = 1E-5;

	test_double(q, b[0][0], tol, "B0_0");
	test_double(q, b[0][1], tol, "B0_1");
	test_double(-q, b[0][2], tol, "B0_2");
	test_double(q, b[1][0], tol, "B1_0");
	test_double(-q, b[1][1], tol, "B1_1");
	test_double(q, b[1][2], tol, "B1_2");
	test_double(-q, b[2][0], tol, "B2_0");
	test_double(q, b[2][1], tol, "B2_1");
	test_double(q, b[2][2], tol, "B2_2");
}

void test_decompose_cartesian_vectors()
{
	double b[3][3];

	b[0][0] = 1;
	b[0][1] = 2;
	b[0][2] = 0;
	b[1][0] = 0;
	b[1][1] = 1;
	b[1][2] = 0;
	b[2][0] = 1;
	b[2][1] = 1;
	b[2][2] = 1;

	double v1[3], v2[3], v3[3];

	decompose_cartesian_vectors(b, v1, v2, v3);

	const double tol = 1E-6;

	double c1[3], c2[3], c3[3];
	for (int i = 0; i < 3; ++i)
	{
		c1[i] = c2[i] = c3[i] = 0;
		for (int j = 0; j < 3; ++j)
		{
			c1[i] += v1[j] * b[j][i];
			c2[i] += v2[j] * b[j][i];
			c3[i] += v3[j] * b[j][i];
		}
	}
	test_double(1, c1[0], tol, "DCV1_0");
	test_double(0, c1[1], tol, "DCV1_1");
	test_double(0, c1[2], tol, "DCV1_2");
	test_double(0, c2[0], tol, "DCV2_0");
	test_double(1, c2[1], tol, "DCV2_1");
	test_double(0, c2[2], tol, "DCV2_2");
	test_double(0, c3[0], tol, "DCV3_0");
	test_double(0, c3[1], tol, "DCV3_1");
	test_double(1, c3[2], tol, "DCV3_2");
}

void test_volume_factors(void)
{
	double a[3][3] = { { 0.28165, 0.28165, 0 }, { 0.28165, 0, 0.28165 }, { 0, 0.28165, 0.28165 } };

	double vf_3d = calculate_volume_factor_3d(a);
	double vf_2d = calculate_volume_factor_2d(a);
	const double tol = 1E-8;

	test_double(1/sqrt(2), vf_3d, tol, "VF_3D");
	test_double(sqrt(3)/2, vf_2d, tol, "VF_2D");
}

int main(void)
{
	test_reciprocal_vectors();

	test_decompose_cartesian_vectors();

	test_volume_factors();

	return 0;
}

