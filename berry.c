#include <assert.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include "berry.h"
#include "hamiltonian.h"
#include "utils.h"
#include "fermi.h"

void destroy_berry_calculator(berry_calculator* this)
{
	free(this->dHdkx);
	free(this->dHdky);
	free(this->H);
	free(this->energies);
	free(this);
}

berry_calculator* create_berry_calculator(const ham_struct* hs)
{
	berry_calculator* bc = alloc_memory(sizeof(berry_calculator));
	bc->hs = hs;
	bc->dHdkx = alloc_memory(sizeof(doublecomplex) * hs->dim * hs->dim);
	bc->dHdky = alloc_memory(sizeof(doublecomplex) * hs->dim * hs->dim);
	bc->H = alloc_memory(sizeof(doublecomplex) * hs->dim * hs->dim);
	bc->energies = alloc_memory(sizeof(double) * hs->dim);
	bc->destroy = destroy_berry_calculator;
	return bc;
}

/*
Calculate the Berry curvature Omega_z = 2 Im <du/dky|du/dkx>
*/
static double berry_curvature(const doublecomplex* dHdkx, const doublecomplex* dHdky, const doublecomplex* eigvecs, const double* eigens, const int dim, const int band, double hG)
{
	double curvature = 0;

	const doublecomplex* band_eigvec = &eigvecs[band*dim];

	for (int i = 0; i < dim; ++i)
	{
		if (band == i)
		{
			continue;
		}

		double tmp1 = (eigens[band] - eigens[i]);
		double tmp2 = 1/(tmp1*tmp1 + hG*hG); // hbar Gamma = 0.2 eV
		 
		if (!isfinite(1/tmp1))
		{
			continue;
		}
		const doublecomplex* eigvec = &eigvecs[i*dim];
		doublecomplex ax = hermitian_matrix_element(dHdkx, eigvec, band_eigvec, dim);
		doublecomplex ay = hermitian_matrix_element(dHdky, band_eigvec, eigvec, dim);
		curvature += 2*(ay.real*ax.imag + ay.imag*ax.real)*tmp2 - 2*(ay.real*ax.real - ay.imag*ax.imag)*tmp2*hG/tmp1;
	}

	return curvature;
}

void berry_curvatures(const berry_calculator* bc, double c1, double c2, double c3, double* curvatures)
{
	const ham_struct* hs = bc->hs;
	int dim = hs->dim;

	hs->gen_ham_derivative(hs, c1, c2, c3, hs->kx_vect[0], hs->kx_vect[1], hs->kx_vect[2], bc->dHdkx);
	hs->gen_ham_derivative(hs, c1, c2, c3, hs->ky_vect[0], hs->ky_vect[1], hs->ky_vect[2], bc->dHdky);
	hs->gen_ham(hs, c1, c2, c3, bc->H);
	diagonalize_hermitian_matrix(bc->H, bc->energies, dim);

	for (int band = 0; band < dim; ++band)
	{
		curvatures[band] = berry_curvature(bc->dHdkx, bc->dHdky, bc->H, bc->energies, dim, band, 0);
	}
}

/*
static double ahe_integrand_t0(berry_calculator* this, double kx, double ky, double kz, double EF)
{
	double ahe_integrand = 0;
	int bnd_lim = this->hs->p_bands_start + this->hs->p_bands_cnt;
	double* energies = this->hs->work_vector;
	int dim = this->hs->dim;
	int inverted_fermi = this->hs->inverted_fermi;
	doublecomplex* evs = this->H;

	// Calculate dH/dkx
	// Since kx_vect describes a vector of length 1/this->L, in program units we have divide by length 1
	this->hs->gen_ham_derivative(this->hs, kx, ky, kz, this->hs->kx_vect[0], this->hs->kx_vect[1], this->hs->kx_vect[2], this->dHdkx);
	// Calculate dH/dky
	// Since ky_vect describes a vector of length 1/this->L, in program units we have divide by length 1
	this->hs->gen_ham_derivative(this->hs, kx, ky, kz, this->hs->ky_vect[0], this->hs->ky_vect[1], this->hs->ky_vect[2], this->dHdky);
	// Calculate and diagonalize the hamiltonian.
	this->hs->gen_ham(this->hs, kx, ky, kz, evs);
	diagonalize_hermitian_matrix(evs, energies, dim);
	
	for (int b1 = this->hs->p_bands_start; b1 < bnd_lim; b1++)
	{
		if (inverted_fermi ? energies[b1] > EF : energies[b1] < EF)
		{
			for (int b2 = 0; b2 < dim; b2++)
			{
				if (b2 >= this->hs->p_bands_start && b2 < bnd_lim && (inverted_fermi ? energies[b2] > EF : energies[b2] < EF))
				{
					continue;
				}
				else
				{
					doublecomplex cx = hermitian_matrix_element(this->dHdkx, evs + b2*dim, evs + b1*dim, dim);
					doublecomplex cy = hermitian_matrix_element(this->dHdky, evs + b2*dim, evs + b1*dim, dim);
					ahe_integrand += berry_curvature_contribution(cx, cy, energies[b1], energies[b2]);
				}
			}
		}
	}
	return ahe_integrand;
}
*/

double ahe_integrand(const berry_calculator *bc, double c1, double c2, double c3, double EF, double hG, double kBTemp)
{
	double ahe_integrand = 0;

	const ham_struct* hs = bc->hs;
	int dim = hs->dim;

	hs->gen_ham_derivative(hs, c1, c2, c3, hs->kx_vect[0], hs->kx_vect[1], hs->kx_vect[2], bc->dHdkx);
	hs->gen_ham_derivative(hs, c1, c2, c3, hs->ky_vect[0], hs->ky_vect[1], hs->ky_vect[2], bc->dHdky);
	hs->gen_ham(hs, c1, c2, c3, bc->H);
	diagonalize_hermitian_matrix(bc->H, bc->energies, dim);

	int b0 = hs->p_bands_start;
	int b1 = hs->p_bands_start + hs->p_bands_cnt;

	if (kBTemp != 0 || hG != 0)
	{
		for (int band = b0; band < b1; ++band)
		{
			double energy = bc->energies[band];
			double factor = hs->inverted_fermi ? inverted_fermi_dirac(energy, EF, kBTemp) : fermi_dirac(energy, EF, kBTemp);
			if (factor != 0)
			{
				ahe_integrand += factor * berry_curvature(bc->dHdkx, bc->dHdky, bc->H, bc->energies, dim, band, hG);
			}
		}
	}
	else
	{
		const double* energies = bc->energies;
		int inverted_fermi = hs->inverted_fermi;
		const doublecomplex* dHdkx = bc->dHdkx;
		const doublecomplex* dHdky = bc->dHdky;
		for (int band = b0; band < b1; ++band)
		{
			double en_band = energies[band];
			const doublecomplex* ev_band = bc->H + band*dim;
			if (inverted_fermi ? (en_band >= EF) : (en_band <= EF)) // Sum over occupied bands
			{
				for (int i = 0; i < dim; ++i)
				{
					double en_i = energies[i];
					if (! (i >= b0 && i < b1 && inverted_fermi ? (en_i >= EF) : (en_i <= EF))) // here sum only over un-occupied states
					{
						const doublecomplex* ev_i = bc->H + i*dim;
						double tmp1 = en_band - en_i;
						double tmp2 = 1 / (tmp1*tmp1 + hG*hG);
						if (isfinite(1/tmp1))
						{
							doublecomplex ax = hermitian_matrix_element(dHdkx, ev_i, ev_band, dim);
							doublecomplex ay = hermitian_matrix_element(dHdky, ev_band, ev_i, dim);
							ahe_integrand += 2*(ay.real*ax.imag + ay.imag*ax.real)*tmp2 - 2*(ay.real*ax.real - ay.imag*ax.imag)*tmp2*hG/tmp1;
						}
					}
				}
			}
		}
	}

	return ahe_integrand;
}

void scan_curvatures(const berry_calculator* bc, double dx, double dy, double dz, FILE* out)
{
	double kx, ky, kz;
	double c1, c2, c3;
	const double dr = 0.001;
	const ham_struct* hs = bc->hs;
	double* curvs = alloc_memory(sizeof(double) * hs->dim);
	for (double r = -1; r <= 1; r += dr)
	{
		kx = r * dx;
		ky = r * dy;
		kz = r * dz;
		translate_into_inner_coords(hs, kx, ky, kz, &c1, &c2, &c3);
		berry_curvatures(bc, c1, c2, c3, curvs);
		fprintf(out, "%g %g %g", kx / hs->L, ky / hs->L, kz / hs->L);
		for (int i = 0; i < hs->dim; i++)
		{
			fprintf(out, " %g", curvs[i] * hs->L * hs->L);
		}
		fprintf(out, "\n");
		fflush(out);
	}
	free(curvs);
}

