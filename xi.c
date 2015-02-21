#include "utils.h"
#include "fermi.h"
#include "hamiltonian.h"

double calculate_xi_integrand(const ham_struct* hs, const doublecomplex * eigenvectors, double* energies, double ef, double kBtemp)
{
	double xi = 0;
	const int size = hs->dim;
	const int limit = hs->p_bands_start + hs->p_bands_cnt;
	const doublecomplex* sx = hs->S_x;
	const doublecomplex* sy = hs->S_y;
	const doublecomplex* sz = hs->S_z;
	
	for (int n = hs->p_bands_start; n < limit; n++)
	{
		double fermi_factor;
		if (hs->inverted_fermi)
		{
			fermi_factor = inverted_fermi_dirac(energies[n], ef, kBtemp);
		}
		else
		{
			fermi_factor = fermi_dirac(energies[n], ef, kBtemp);
		}
		double mean_spin = 0;
		if (hs->magdir[0] != 0)
		{
			mean_spin += hs->magdir[0] * hermitian_matrix_mean_value(sx, &eigenvectors[n * size], size);
		}
		if (hs->magdir[1] != 0)
		{
			mean_spin += hs->magdir[1] * hermitian_matrix_mean_value(sy, &eigenvectors[n * size], size);
		}
		if (hs->magdir[2] != 0)
		{
			mean_spin += hs->magdir[2] * hermitian_matrix_mean_value(sz, &eigenvectors[n * size], size);
		}

		xi -= fermi_factor * mean_spin;
	}

	return xi;
}

double xi_integrand(const ham_struct* hs, double kx, double ky, double kz, double EF, double kBTemp)
{
	hs->gen_ham(hs, kx, ky, kz, hs->work_matrix);
	diagonalize_hermitian_matrix(hs->work_matrix, hs->work_vector, hs->dim);
	return calculate_xi_integrand(hs, hs->work_matrix, hs->work_vector, EF, kBTemp);
}
