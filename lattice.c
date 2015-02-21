#include <math.h>
#include "utils.h"
#include "constants.h"
#include "hamiltonian.h"

void calculate_reciprocal_lattice_vectors(const double a[3][3], double b[3][3])
{
	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			b[i][j] = a[i][j];				
		}
	}

	invert_double_matrix(&b[0][0], 3);
	
	for (int i = 0; i < 3; ++i)
	{
		b[i][i] *= 2 * M_PI;
		for (int j = 0; j < i; ++j)
		{
			double tmp = b[i][j];
			b[i][j] = 2 * M_PI * b[j][i];
			b[j][i] = 2 * M_PI * tmp;
		}
	}
}

void decompose_cartesian_vectors(const double basis[3][3], double kx[3], double ky[3], double kz[3])
{
	double tmp[3][3];
	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			tmp[i][j] = basis[i][j];
		}
	}

	invert_double_matrix(&tmp[0][0], 3);

	for (int i = 0; i < 3; ++i)
	{
		kx[i] = tmp[0][i];
		ky[i] = tmp[1][i];
		kz[i] = tmp[2][i];
	}
}

void calculate_primitive_cosines(const double a[3][3], const double b[3][3], double cosines[3])
{
	for (int i = 0; i < 3; ++i)
	{
		double sp = 0;
		double norm_a = 0;
		double norm_b = 0;

		for (int j = 0; j < 3; ++j)
		{
			sp += a[i][j]*b[i][j];
			norm_a += a[i][j]*a[i][j];
			norm_b += b[i][j]*b[i][j];
		}
		cosines[i] = sp / sqrt(norm_a * norm_b);
	}
}

static void vector_product(const double x[3], const double y[3], double vp[3])
{
	vp[0] = x[1]*y[2] - x[2]*y[1];
	vp[1] = x[2]*y[0] - x[0]*y[2];
	vp[2] = x[0]*y[1] - x[1]*y[0];
}

double calculate_volume_factor_3d(const double a[3][3])
{
	double a12[3];

	vector_product(a[0], a[1], a12);

	double volume = 0;
	double l[3];
	l[0] = l[1] = l[2] = 0;
	for (int i = 0; i < 3; ++i)
	{
		volume += a12[i]*a[2][i];
		l[0] += a[0][i]*a[0][i];
		l[1] += a[1][i]*a[1][i];
		l[2] += a[2][i]*a[2][i];
	}
	return fabs(volume) / sqrt(l[0]*l[1]*l[2]);
}

double calculate_volume_factor_2d(const double a[3][3])
{
	 double vp[3];
	 vector_product(a[0], a[1], vp);
	 return sqrt((vp[0]*vp[0] + vp[1]*vp[1] + vp[2]*vp[2]) / (a[0][0]*a[0][0] + a[0][1]*a[0][1] + a[0][2]*a[0][2]) / (a[1][0]*a[1][0] + a[1][1]*a[1][1] + a[1][2]*a[1][2]));
}

double calculate_ezz(double exx)
{
	return -2*exx*0.453;
}

/* hs->L and hs->exx must be set before calling this function. aSI_zs should be for zero strain and in nm */
void setup_geometry(ham_struct* hs, const double aSI_zs[3][3])
{
	double b[3][3]; // Reciprocal lattice vectors
	double aSI[3][3];

	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			aSI[i][j] = aSI_zs[i][j];
		}
	}
	// Apply strain.
	printf("exx == %g\n", hs->exx);
	double exx = hs->exx;
	double eyy = hs->exx;
	double ezz = calculate_ezz(exx);
	for (int i = 0; i < 3; ++i)
	{
		aSI[i][0] *= 1 + exx;
		aSI[i][1] *= 1 + eyy;
		aSI[i][2] *= 1 + ezz;
	}

	calculate_reciprocal_lattice_vectors(aSI, b);
	printf("Reciprocal lattice vectors\nb1 == [%g, %g, %g]\n", b[0][0], b[0][1], b[0][2]);
	printf("b2 == [%g, %g, %g]\n", b[1][0], b[1][1], b[1][2]);
	printf("b3 == [%g, %g, %g]\n", b[2][0], b[2][1], b[2][2]);
	printf("Primitive lattice vectors\na1 == [%g, %g, %g]\n", aSI[0][0], aSI[0][1], aSI[0][2]);
	printf("a2 == [%g, %g, %g]\n", aSI[1][0], aSI[1][1], aSI[1][2]);
	printf("a3 == [%g, %g, %g]\n", aSI[2][0], aSI[2][1], aSI[2][2]);

	// Transform b into our unit system
	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			b[i][j] *= hs->L;
		}
	}

	// Calculate lengths of reciprocal lattice vectors.
	for (int i = 0; i < 3; ++i)
	{
		hs->bvect_lengths[i] = 0;
		for (int j = 0; j < 3; ++j)
		{
			hs->bvect_lengths[i] += b[i][j]*b[i][j];
		}
		hs->bvect_lengths[i] = sqrt(hs->bvect_lengths[i]);
	}
	printf("Reciprocal lattice vectors\nb1 == [%g, %g, %g]\n", b[0][0], b[0][1], b[0][2]);
	printf("b2 == [%g, %g, %g]\n", b[1][0], b[1][1], b[1][2]);
	printf("b3 == [%g, %g, %g]\n", b[2][0], b[2][1], b[2][2]);

	printf("bvl == [%g, %g, %g]\n", hs->bvect_lengths[0], hs->bvect_lengths[1], hs->bvect_lengths[2]);
	// Express vectors (1,0,0),(0,1,0) and (0,0,1) in the b basis.
	decompose_cartesian_vectors(&b[0][0], hs->kx_vect, hs->ky_vect, hs->kz_vect);

	double cosines[3];
	calculate_primitive_cosines(&aSI[0][0], &b[0][0], cosines);
	hs->cosine_product_12 = cosines[0]*cosines[1];
	hs->cosine_product = hs->cosine_product_12*cosines[2];
	hs->volume_factor = calculate_volume_factor_3d(&aSI[0][0]);
	hs->volume_factor_12 = calculate_volume_factor_2d(&aSI[0][0]);
}

