#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "hamiltonian.h"
#include "tightbinding.h"
#include "anisotropy.h"
#include "sixband.h"

int main(int argc, char* argv[])
{
	const double Delta = 0.15;
	const double pSI = 0.1;
	const double exx = 0.01;
	
	ham_struct* hs = create_tightbinding_bulk_hamiltonian_structure(EXACT, GEOM_XY, JANCU, exx);
//	ham_struct* hs = create_sixband_hamiltonian_structure(0.042, exx);
	const double Vol = hs->L*hs->L*hs->L;
	printf("Vol == %g\n", Vol);
	const double p = pSI * Vol;
	
	hs->D = Delta;
	const double pmax = 1.05*p;
	const int ncmax = 20;
	const double cstep=hs->find_cstep(hs, pmax, 2*ncmax);
	const double prefactor = hs->sum2integral_prefactor(hs) * pow(cstep, hs->periodic_bc_dim) / pow(hs->L, hs->periodic_bc_dim);
	printf("Prefactor == %g\n", prefactor);
	fflush(stdout);
	const double* un_en = unax_energies(hs, Delta, p, p, 1.05*p, ncmax);

	printf("Anizotropy energy [eV/nm^3] == %g\n", un_en[0] * prefactor);

	// Call this instead of free.
	hs->destroy(hs); // DIE!!!
	free(un_en);

	return 0;
}
