#include <stdio.h>
#include <math.h>
#include "tightbinding.h"
#include "hamiltonian.h"

int main(void)
{
	const double D = 0.15;
	const double pSI = 0.1;
	const double exx = 0.01;
	ham_struct* hs = create_tightbinding_bulk_hamiltonian_structure(EXACT, GEOM_XY, JANCU, exx);
	const double Vol = pow(hs->L, 3);
	const int nkmax = 50;
	const double p = pSI * Vol;
	hs->magdir[0] = 0;
	hs->magdir[1] = 0;
	hs->magdir[2] = 1;
	double ef_z = find_single_EF(hs, p, D, nkmax);
	hs->magdir[0] = 1;
	hs->magdir[1] = 0;
	hs->magdir[2] = 0;
	double ef_x = find_single_EF(hs, p, D, nkmax);
	printf("ef_z == %g, ef_x == %g\n", ef_z, ef_x);
	return 0;
}

