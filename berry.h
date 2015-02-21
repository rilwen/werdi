#ifndef _BERRY_H
#define _BERRY_H

#include <stdio.h>
#include "doublecomplex.h"
#include "hamiltonian.h"

typedef struct my_berry_calculator {
	const ham_struct* hs; /* Hamiltonian structure */
	doublecomplex* dHdkx; /* Work matrix to store dH/dkx in */
	doublecomplex* dHdky; /* Work matrix to store dH/dky in */
	doublecomplex* H; /* Work matrix to store the Hamiltonian (later, its eigenvectors) in */
	double* energies; /* Work vector to store the eigenstate energies */

	/* "Destructor" */
	void (*destroy) (struct my_berry_calculator* this);
} berry_calculator;

berry_calculator* create_berry_calculator(const ham_struct* hs);

/* Calculate Berry curvatures for selected bands, in program units */
void berry_curvatures(const berry_calculator* this, double c1, double c2, double c3, double* curvatures);

/* Calculate AHE integrand, in program units */
double ahe_integrand(const berry_calculator* this, double c1, double c2, double c3, double EF, double hG, double kBTemp);

void scan_curvatures(const berry_calculator* bc, double dx, double dy, double dz, FILE* out);

#endif
