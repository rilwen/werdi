#include <stdio.h>
#include "testing.h"
#include "berry.h"
#include "hamiltonian.h"
#include "parabolic.h"
#include "sixband.h"
#include "tightbinding.h"
#include "utils.h"

int main(void)
{
	ham_struct* hs = create_parabolic_hamiltonian_structure();
	berry_calculator* bc = create_berry_calculator(hs);

	hs->D = 0.1; // Non-zero spin splitting.

	const double EF = 0.2;
	const double kBTemp = 0.05;

	const double cx = 1;
	const double cy = 0;
	const double cz = 0;

	for (double r = -1; r <= 1; r += 0.2)
	{
		double kx = r * cx;
		double ky = r * cy;
		double kz = r * cz;
		double curvs[PARABOLIC_HAM_DIM];
		berry_curvatures(bc, kx, ky, kz, curvs);
		for (int i = 0; i < hs->dim; i++)
		{
			test_double(0, curvs[i], 1E-16, "TEST_CURV_PARABOLIC");
			test_double(0, ahe_integrand(bc, kx, ky, kz, EF, 0, 0), 1E-16, "TEST_AHET0_PARABOLIC");
			test_double(0, ahe_integrand(bc, kx, ky, kz, EF, 0, kBTemp), 1E-16, "TEST_AHETNZ_PARABOLIC");
		}
	}


	bc->destroy(bc);
	hs->destroy(hs);

	hs = create_sixband_hamiltonian_structure(0.042, 0);
	hs->D = 0.15;
	bc = create_berry_calculator(hs);
	
	FILE* file = fopen("test_berry_results_6b.dat", "w");

	for (double r = -1; r <= 1; r += 0.00001)
	{
		double curvs[SIXBAND_HAM_DIM];
		double kx = r * cx;
		double ky = r * cy;
		double kz = r * cz;
		fprintf(file, "%g %g %g", kx, ky, kz);
		berry_curvatures(bc, kx, ky, kz, curvs);
		for (int i = 0; i < hs->dim; i++)
		{
			fprintf(file, " %g", curvs[i] * hs->L * hs->L);
		}
		fprintf(file, " %g %g", ahe_integrand(bc, kx, ky, kz, EF, 0, 0) * hs->L * hs->L, ahe_integrand(bc, kx, ky, kz, EF, 0, kBTemp) * hs->L * hs->L);
		fprintf(file, "\n");
	}
	
	fclose(file);
	
	printf("Saved curvatures and AHE integrand for 6-band model\n");
	fflush(stdout);

	return 0;

	hs = create_tightbinding_bulk_hamiltonian_structure(EXACT, GEOM_XY, JANCU, 0);
	hs->D = 0.17;
	bc = create_berry_calculator(hs);
	
	file = fopen("test_berry_results_tb_jancu_exact.dat", "w");
	scan_curvatures(bc, cx, cy, cz, file);
	fclose(file);

	bc->destroy(bc);
	hs->destroy(hs);

	printf("Saved curvatures and AHE integrand for TB XY Jancu model with exact differentiation\n");
	fflush(stdout);

	hs = create_tightbinding_bulk_hamiltonian_structure(FIVEPOINT, GEOM_XY, JANCU, 0);
	hs->D = 0.17;
	bc = create_berry_calculator(hs);
	
	file = fopen("test_berry_results_tb_jancu_fivepoint.dat", "w");
	scan_curvatures(bc, cx, cy, cz, file);
	fclose(file);
	
	bc->destroy(bc);
	hs->destroy(hs);

	printf("Saved curvatures and AHE integrand for TB XY Jancu model with 5-point differentiation\n");
	fflush(stdout);
	
	return 0;
}
