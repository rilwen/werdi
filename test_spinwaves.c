#include <stdio.h>
#include "hamiltonian.h"
#include "testing.h"
#include "sixband.h"
#include "spinwaves.h"
#include "doublecomplex.h"

int main(void)
{
	const double D = 0.1523;
	const double p = 0.1;
	double EF;

	ham_struct* hs = create_sixband_hamiltonian_structure(0.042, 0);
	hs->D = D;
	EF = find_single_EF(hs, p * hs->L * hs->L, D, 30);
	spinwave_calculator* sc = create_spinwave_calculator(hs);


	double qx, qy, qz;
	qx = 0.1;
	qy = qz = 0.2;
	double sw_c1, sw_c2, sw_c3;
	translate_into_inner_coords(hs, qx, qy, qz, &sw_c1, &sw_c2, &sw_c3);

	FILE* out = fopen("Epm_Epp_integrand_6b.dat", "w");
	for (double kx = -1; kx <= 1; kx += 0.01)
	{
		double c1, c2, c3;
		translate_into_inner_coords(hs, kx, 0, 0, &c1, &c2, &c3);
		doublecomplex epmint = Epm_integrand(sc, c1, c2, c3, sw_c1, sw_c2, sw_c3, EF, 0);
		doublecomplex eppint = Epp_integrand(sc, c1, c2, c3, sw_c1, sw_c2, sw_c3, EF, 0);
		fprintf(out, "%g %g %g %g %g %g %g\n", kx, 0.0, 0.0, epmint.real, epmint.imag, eppint.real, eppint.imag);
	}

	fclose(out);

	sc->destroy(sc);
	hs->destroy(hs);
}
