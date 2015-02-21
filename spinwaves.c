#include <math.h>
#include "spinwaves.h"
#include "hamiltonian.h"
#include "doublecomplex.h"
#include "fermi.h"
#include "constants.h"

/*
 * Spin wave spectrum is calculated from p bands only.
 *
 * sc - spinwave calculator object
 * c1, c2, c3 - reciprocal lattice coordinates of the point in which we calculate the integrand
 * sw_c1, sw_c2, sw_c3 - reciprocal lattice coordinates corresponding to the spinwave vector
 * plusminus - 0 for E_{+-}, 1 for E_{++}
 * fermi_energy
 * kBTemp
 */
static doublecomplex spinwave_integrand(spinwave_calculator* sc, double c1, double c2, double c3, double sw_c1, double sw_c2, double sw_c3, int plusminus, double fermi_energy, double kBTemp)
{
	ham_struct* hs = sc->hs;
	doublecomplex res = {.real = 0, .imag = 0};
	int zero_wave;
	if (sw_c1 == 0 && sw_c2 == 0 && sw_c3 == 0)
	{
		zero_wave = 1;
	}
	else
	{
		zero_wave = 0;
	}

	hs->gen_ham(hs, c1, c2, c3, sc->eigs_k);
	diagonalize_hermitian_matrix(sc->eigs_k, sc->energs_k, hs->dim);
	double* en_k = sc->energs_k;
	doublecomplex* ev_k = sc->eigs_k;

	const doublecomplex* second_op = plusminus ? hs->S_plus : hs->S_minus;

	double* en_kpq;
	doublecomplex* ev_kpq;
	if (! zero_wave)
	{
		hs->gen_ham(hs, c1 + sw_c1, c2 + sw_c2, c3 + sw_c3, sc->eigs_kplusq);
		diagonalize_hermitian_matrix(sc->eigs_kplusq, sc->energs_kplusq, hs->dim);
		en_kpq = sc->energs_kplusq;
		ev_kpq = sc->eigs_kplusq;
	}
	else
	{
		en_kpq = en_k;
		ev_kpq = ev_k;
	}

	int p_bands_end = hs->p_bands_start + hs->p_bands_cnt;
	for (int b1 = hs->p_bands_start; b1 < p_bands_end; b1++)
	{
		for (int b2 = hs->p_bands_start; b2 < p_bands_end; b2++)
		{
			if (zero_wave && b1 == b2)
			{
				/*
				 * For spin wave with zero wave vector, we skip elements with equal band number.
				 */
				continue;
			}

			double en1 = en_k[b1];
			double en2 = en_kpq[b2];

			double tmp;
			if (hs->inverted_fermi)
			{
				tmp = (inverted_fermi_dirac(en1, fermi_energy, kBTemp) - inverted_fermi_dirac(en2, fermi_energy, kBTemp)) / (en1 - en2);
				if (!isfinite(tmp))
				{
					tmp = inverted_fermi_dirac_derivative(en1, fermi_energy, kBTemp);
				}
			}
			else
			{
				tmp = (fermi_dirac(en1, fermi_energy, kBTemp) - fermi_dirac(en2, fermi_energy, kBTemp)) / (en1 - en2);
				if (!isfinite(tmp))
				{
					tmp = fermi_dirac_derivative(en1, fermi_energy, kBTemp);
				}
			}
			doublecomplex splus = complex_matrix_element(hs->S_plus, ev_k + b1 * hs->dim, ev_kpq + b2 * hs->dim, hs->dim);
			doublecomplex val2 = complex_matrix_element(second_op, ev_kpq + b2 * hs->dim, ev_k + b1 * hs->dim, hs->dim);
			res.real -= tmp * (splus.real * val2.real - splus.imag * val2.imag);
			res.imag -= tmp * (splus.real * val2.imag + splus.imag * val2.real);
		}
	}

	return res;
}

void destroy_spinwave_calculator(spinwave_calculator* this)
{
	free(this->energs_k);
	free(this->energs_kplusq);
	free(this->eigs_k);
	free(this->eigs_kplusq);
	free(this);
}

spinwave_calculator* create_spinwave_calculator(ham_struct* hs)
{
	spinwave_calculator* sc = alloc_memory(sizeof(spinwave_calculator));

	sc->hs = hs;
	sc->energs_k = alloc_memory(sizeof(double) * hs->dim);
	sc->energs_kplusq = alloc_memory(sizeof(double) * hs->dim);
	sc->eigs_k = alloc_memory(sizeof(doublecomplex) * hs->dim * hs->dim);
	sc->eigs_kplusq = alloc_memory(sizeof(doublecomplex) * hs->dim * hs->dim);
	sc->destroy = destroy_spinwave_calculator;

	return sc;
}

doublecomplex Epm_integrand(spinwave_calculator* sc, double c1, double c2, double c3, double sw_c1, double sw_c2, double sw_c3, double fermi_energy,
	double kBTemp)
{
	return spinwave_integrand(sc, c1, c2, c3, sw_c1, sw_c2, sw_c3, 0, fermi_energy, kBTemp);
}

doublecomplex Epp_integrand(spinwave_calculator* sc, double c1, double c2, double c3, double sw_c1, double sw_c2, double sw_c3, double fermi_energy,
	double kBTemp)
{
	return spinwave_integrand(sc, c1, c2, c3, sw_c1, sw_c2, sw_c3, 1, fermi_energy, kBTemp);
}

double qDebaySI(double x)
{
	return pow(24*x/M_PI, 1.0/3)*M_PI/a0SI;
}

