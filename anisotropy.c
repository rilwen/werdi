#include <math.h>
#include <stdio.h>
#include <assert.h>
#include "doublecomplex.h"
#include "constants.h"
#include "utils.h"
#include "hamiltonian.h"

double* total_energies(ham_struct* hs, double D, double pmin, double dp, int plen, int ncmax)
{
	const int en_cnt = calculate_en_cnt(hs, ncmax);

	double* energies = alloc_memory(en_cnt*sizeof(double));
	double* fermi_energies = alloc_memory(plen*sizeof(double));
	double* total_energies = alloc_memory(plen*sizeof(double));

	find_EF_preallocated(hs, pmin, dp, plen, D, 0, 1, ncmax, energies, fermi_energies);
	
	for (int i=0; i<plen; i++) {
		total_energies[i]=0;
		double ef = fermi_energies[i];
		for (int j=0; j<en_cnt; j++) {
			double en = energies[j];
			if (hs->inverted_fermi ? en < ef : en > ef)
				break;
			total_energies[i] += en;
		}
	}

	free(energies);
	free(fermi_energies);

	return total_energies;
}

double* unax_energies(ham_struct* hs, double D, double pmin, double dp, double pmax, const int ncmax)
{
	const int plen = count_values(pmin, pmax, dp);

	hs->magdir[0] = 1;
	hs->magdir[1] = 0;
	hs->magdir[2] = 0;
	const double* total_en_100 = total_energies(hs, D, pmin, dp, plen, ncmax);
	hs->magdir[0] = 0;
	hs->magdir[1] = 0;
	hs->magdir[2] = 1;
	const double* total_en_001 = total_energies(hs, D, pmin, dp, plen, ncmax);

	double* unax_en = alloc_memory(plen*sizeof(double));
	for (int i=0; i<plen; i++) {
		unax_en[i] = total_en_001[i] - total_en_100[i];
		printf("Anisotropy energy no %d == %g\n", i, unax_en[i]);
		fflush(stdout);
	}

	free(total_en_100);
	free(total_en_001);

	return unax_en;
}

