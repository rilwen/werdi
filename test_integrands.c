#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <tb2.h>
#include <cuba.h>
#include "berry.h"
#include "hamiltonian.h"
#include "sixband.h"
#include "tightbinding.h"

/* Hamiltonian structure used to carry out calculations in model-independent manner (mostly) */
ham_struct* hamiltonian_structure;
berry_calculator* berry_calc;
double EF;

double (*integrand_for_cuba)(double, double, double);

double gaussian(double c1, double c2, double c3)
{
	double kx = c1 * hamiltonian_structure->bvect_lengths[0];
	double ky = c2 * hamiltonian_structure->bvect_lengths[1];
	double kz = c3 * hamiltonian_structure->bvect_lengths[2];
//	double kx = c1; double ky = c2; double kz = c3;
	double nrm2  = kx*kx + ky*ky + kz*kz;
	return exp( - nrm2 / 2);
}

double ahe(double c1, double c2, double c3)
{
	return ahe_integrand(berry_calc, c1, c2, c3, EF, 0, 0);
}

static void cuba_integrand(const int* ndim, const double x[], const int* ncomp, double f[])
{
	// x[] -- coordinates from (0,1) period each
	// f[0] -- value of the integrand for CUBA library

	assert(*ncomp == 1);

	f[0] = hamiltonian_structure->bril_zone_integrand(hamiltonian_structure, integrand_for_cuba, x);
}

static double calculate_integral(void)
// Calculates Xi in SI units for Delta given in [eV]
{
	double result;
	int nregions, neval, fail;
	double error, prob;

	
	Cuhre(hamiltonian_structure->periodic_bc_dim, 1, cuba_integrand, 1E-3, 1E-8, 0, 0, 2000000, 11, &nregions, &neval, &fail, &result, &error, &prob);

	switch (fail) {
		case 0:
			fprintf(stderr, "The desired accuracy was reached.\n");
			break;
		case -1:
			fprintf(stderr, "Dimension out of range\n");
			break;
		case 1:
			fprintf(stderr, "The desired accuracy was NOT reached\n");
			break;
	}

//	printf("error: %g\n", errorXi);
//	printf("prob for real: %g\n", probXi);
	double prefactor = hamiltonian_structure->sum2integral_prefactor(hamiltonian_structure);
	printf("Prefactor == %g\n", prefactor);
	printf("Absolute error == %g, relative error == %g\n", prefactor * error, fabs(error/result));

	return  prefactor * result;

}

int main(int argc, char* argv[])
{
	const double exx = 0.01;
	hamiltonian_structure = create_sixband_hamiltonian_structure(0.042, exx);
	const double Delta = 0.15;
	const double pSI = 0.1;

	
/*
	initialize_kane_spin_matrices();
	initialize_kane_constants(0);
	hamiltonian_structure = create_kane_hamiltonian_structure(0.042, 0);
*/
	integrand_for_cuba = gaussian;
	double integral = calculate_integral();
	printf("Gaussian == %g\n", integral);
	
	berry_calc = create_berry_calculator(hamiltonian_structure);

	integrand_for_cuba = ahe;
	double Vol = pow(hamiltonian_structure->L, 3);
	EF = find_single_EF(hamiltonian_structure, pSI*Vol, Delta, 40);
	hamiltonian_structure->D = Delta;
	printf("p[nm^-3] == %g\n", pSI);
	printf("Delta == %g\n", Delta);
	printf("exx == %g\n", exx);
	integral = calculate_integral();
	printf("AHE in 6band with zero temp simplfication == %g\n", integral);

	// Call this instead of free.
	hamiltonian_structure->destroy(hamiltonian_structure); // DIE!!!
	berry_calc->destroy(berry_calc);


	hamiltonian_structure = create_tightbinding_bulk_hamiltonian_structure(EXACT, GEOM_XY, JANCU, exx);
	berry_calc = create_berry_calculator(hamiltonian_structure);
	Vol = pow(hamiltonian_structure->L, 3);
	hamiltonian_structure->D = Delta;
	hamiltonian_structure->magdir[0] = 0;
	hamiltonian_structure->magdir[1] = 0;
	hamiltonian_structure->magdir[2] = 1;
	EF = find_single_EF(hamiltonian_structure, pSI * Vol, Delta, 40);
	integral = calculate_integral();
	printf("AHE in TB(Jancu,XY) with zero temp simplification and magnetization along [0,0,1] == %g\n", integral);
	hamiltonian_structure->magdir[0] = 1;
	hamiltonian_structure->magdir[1] = 0;
	hamiltonian_structure->magdir[2] = 0;
	EF = find_single_EF(hamiltonian_structure, pSI * Vol, Delta, 40);
	integral = calculate_integral();
	printf("AHE in TB(Jancu,XY) with zero temp simplification and magnetization along [1,0,0] == %g\n", integral);

	return 0;
}
