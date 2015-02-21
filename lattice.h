#ifndef _LATTICE_H
#define _LATTICE_H

#include "hamiltonian.h"

void calculate_reciprocal_lattice_vectors(const double a[3][3], double b[3][3]);

void decompose_cartesian_vectors(const double basis[3][3], double kx[3], double ky[3], double kz[3]);

double calculate_volume_factor_3d(const double a[3][3]);

double calculate_volume_factor_2d(const double a[3][3]);

/* Function for calculating the zz strain for given xx strain. */
double calculate_ezz(double exx);

/* hs->L and hs->exx must be set before calling this function. aSI_zs should be for zero strain and in nm */
void setup_geometry(ham_struct* hs, const double aSI_zs[3][3]);

#endif
