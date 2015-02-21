#ifndef _SPINWAVES_H
#define _SPINWAVES_H

#include "hamiltonian.h"
#include "utils.h"

typedef struct my_spinwave_calculator
{
	ham_struct* hs;
	doublecomplex* eigs_k;
	doublecomplex* eigs_kplusq;
	double* energs_k;
	double* energs_kplusq;
	void (*destroy) (struct my_spinwave_calculator* this);
} spinwave_calculator;

spinwave_calculator* create_spinwave_calculator(ham_struct* hs);

/*
 * Spin wave spectrum is calculated from p bands only.
 *
 * sc - spinwave calculator object
 * c1, c2, c3 - reciprocal lattice coordinates of the point in which we calculate the integrand
 * sw_c1, sw_c2, sw_c3 - reciprocal lattice coordinates corresponding to the spinwave vector
 * fermi_energy
 * kBTemp
 */
doublecomplex Epm_integrand(spinwave_calculator* sc, double c1, double c2, double c3, double sw_c1, double sw_c2, double sw_c3, double fermi_energy, double kBTemp);
doublecomplex Epp_integrand(spinwave_calculator* sc, double c1, double c2, double c3, double sw_c1, double sw_c2, double sw_c3, double fermi_energy, double kBTemp);

/*
 * qDebay [nm^-1]
 */
double qDebaySI(double x);

#endif
