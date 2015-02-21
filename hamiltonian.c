#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include "hamiltonian.h"
#include "doublecomplex.h"
#include "utils.h"
#include "constants.h"

const double diff_eps = 1E-6; // Step size for numerical differentiation (for foot size 39)

void discretized_5pt_derivative(ham_struct * this, double kx, double ky, double kz, double cx, double cy, double cz, doublecomplex matrix[][])
{
	doublecomplex* work_matrix = this->work_matrix;

	this->gen_ham(this, kx, ky, kz, matrix); 

	kx += 2 * cx * diff_eps;
	ky += 2 * cy * diff_eps;
	kz += 2 * cz * diff_eps;
	this->gen_ham(this, kx, ky, kz, work_matrix);
	add_matrix(matrix, work_matrix, -1, this->dim);

	kx -= cx * diff_eps;
	ky -= cy * diff_eps;
	kz -= cz * diff_eps;
	this->gen_ham(this, kx, ky, kz, work_matrix);
	add_matrix(matrix, work_matrix, 8, this->dim);

	kx -= 2 * cx * diff_eps;
	ky -= 2 * cy * diff_eps;
	kz -= 2 * cz * diff_eps;
	this->gen_ham(this, kx, ky, kz, work_matrix);
	add_matrix(matrix, work_matrix, -8, this->dim);

	kx -= cx * diff_eps;
	ky -= cy * diff_eps;
	kz -= cz * diff_eps;
	this->gen_ham(this, kx, ky, kz, work_matrix);
	add_matrix(matrix, work_matrix, 1, this->dim);

	multiply_matrix(matrix, 1.0 / (12 * diff_eps), this->dim);
}

/*
all_energies: pointer to array of this->dim doubles
pband_energies: pointer to array of this->p_bands_cnt doubles

Selects the energies belonging to p-hole bands.
*/
static void select_pband_energies(ham_struct* this, double* restrict all_energies, double* restrict pband_energies)
{
	memcpy(pband_energies, all_energies + this->p_bands_start, sizeof(double) * this->p_bands_cnt);
}

/* It is assumed that energies aray has enough storage allocated. The storage can be
computed safely using the function calculate_en_cnt(this, nkmax) */
static double generate_energies(ham_struct* this, double cstep, int nkmax, double *energies)
{
	double extr_edgeen = this->inverted_fermi ? -1E20 : 1E20;
	double energy_minimum = this->inverted_fermi ? -1E20 : 1E20; // energy of the first occupied state
	double energy_minimum_c[3];
	int energy_minimum_band = -1;
	double* curr_pband_energies = energies;
	int n3_start, n3_end;
	if (this->periodic_bc_dim >= 3)
	{
		n3_start = - nkmax;
		n3_end = nkmax;
	}
	else
	{
		n3_start = 0;
		n3_end = 1;
	}
	const int pbc = this->p_bands_cnt;

	for (int nkx = - nkmax; nkx < nkmax; ++nkx)
	{
		double c1 = nkx * cstep;
		for (int nky = - nkmax; nky < nkmax; ++nky)
		{
			double c2 = nky * cstep;
			for (int nkz = n3_start; nkz < n3_end; ++nkz)
			{
				double c3 = nkz * cstep;
				this->gen_ham(this, c1, c2, c3, this->work_matrix);
				diagonalize_hermitian_matrix(this->work_matrix, this->work_vector, this->dim);
				select_pband_energies(this, this->work_vector, curr_pband_energies);

				for (int i = 0; i < pbc; ++i)
				{
					double cen = curr_pband_energies[i];
					if(this->inverted_fermi ? (cen > energy_minimum) : (cen < energy_minimum))
					{
						energy_minimum = cen;
						energy_minimum_c[0] = c1;
						energy_minimum_c[1] = c2;
						energy_minimum_c[2] = c3;
						energy_minimum_band = i + this->p_bands_start;
					}
				}					
				
				if (nkx == - nkmax || nkx == nkmax - 1 || nky == - nkmax || nky == nkmax - 1
					|| nkz == - nkmax || nkz == nkmax - 1)
				{
					// Look for the minimal or maximal energy on the edge.
					for (int i = 0; i < pbc; ++i)
					{
						double cen = curr_pband_energies[i];
						if (this->inverted_fermi ? (cen > extr_edgeen) : (cen < extr_edgeen))
						{
							extr_edgeen = cen;
						}
					}
				}

				curr_pband_energies += pbc;
			}
		}
	}

	printf("Energy \"minimum\" [eV] == %g found for c coordinates == (%g, %g, %g) and band %i\n", energy_minimum,
		energy_minimum_c[0], energy_minimum_c[1], energy_minimum_c[2], energy_minimum_band);
	fflush(stdout);
	
	return extr_edgeen;
}

static int cmp_double_asc(const void *p1, const void *p2)
{
	if (*((const double *)p1) < *((const double *)p2)) {
		return -1;
	} else if (*((const double *)p1) > *((const double *)p2)) {
		return 1;
	} else {
		return 0;
	}
}

static int cmp_double_desc(const void *p1, const void *p2)
{
	if (*((const double *)p1) > *((const double *)p2)) {
		return -1;
	} else if (*((const double *)p1) < *((const double *)p2)) {
		return 1;
	} else {
		return 0;
	}
}

/* Return the number of values in a for loop (for val = min; val <= max. val += step) */
int count_values(double min, double max, double step)
{
	return (int) round((max-min)/step) + 1;
}

/*
Factor by which we multiply integrals over d^3 c to convert them into limits of sums divided by V for V --> infty
*/
double sum2integral_prefactor_3d(ham_struct* this)
{
	return  this->bvect_lengths[0] * this->bvect_lengths[1] * this->bvect_lengths[2] * this->cosine_product /8/M_PI/M_PI/M_PI / this->volume_factor;
}

double sum2integral_prefactor_2d(ham_struct* this)
{
	return this->bvect_lengths[0] * this->bvect_lengths[1] * this->cosine_product_12 / 4 / M_PI / M_PI / this->volume_factor_12 / this->width;
}

/* Calculate the index of the energy corresponding to concentration p and k-vector step cstep */
static inline int concentration_to_index(ham_struct* this, double p, double cstep)
{
	int idx = (int) floor(p / pow(cstep, this->periodic_bc_dim) / this->sum2integral_prefactor(this)) - 1;
	return idx >= 0 ? idx : 0;
}

static double index_to_concentration(ham_struct* this, int idx, double cstep)
{
	if (idx < 0)
	{
		return 0;
	}
	return (idx + 1) * this->sum2integral_prefactor(this) * pow(cstep, this->periodic_bc_dim);
}

int calculate_en_cnt(ham_struct* this, int ncmax)
{
	return this->p_bands_cnt * pow(2 * ncmax, this->periodic_bc_dim);
}

/*
 * pmax - concentration which we expect to be the upper bound of the concentration we will obtain from the calculation. Some implementations
 * of the hamiltonian may not depend on this parameter.
 */
double find_concentration(ham_struct* this, double EF, int nkmax, double pmax)
{
	int en_cnt = calculate_en_cnt(this, nkmax);
	double cstep = this->find_cstep(this, pmax, 2*nkmax);

	double* energies = alloc_memory(sizeof(double) * en_cnt);

	generate_energies(this, cstep, nkmax, energies);
	int cnt = 0;
	for (int i = 0; i < en_cnt; i++)
	{
		cnt += this->inverted_fermi ? (energies[i] >= EF ? 1 : 0) : (energies[i] <= EF ? 1 : 0); // I felt so cool writing nested ?: operators :D		
	}
	free(energies);
	return index_to_concentration(this, cnt - 1, cstep);
}

/* It is assumed that EFtab has at least plen*Dlen*sizeof(double) size and energies has size at least calculate_en_cnt(this,nkmax)*sizeof(double */
void find_EF_preallocated(ham_struct* this, double pmin, double dp, int plen, double Dmin, double dD, int Dlen, int nkmax, double* energies, double* EFtab)
{
	assert(pmin > 0);
	assert(plen > 0);
	assert(dp > 0 || plen == 1);
	assert(Dlen > 0);
	assert(dD > 0 || Dlen == 1);
	
	const int en_cnt = calculate_en_cnt(this, nkmax);
	const double pmax = pmin + (plen - 1) * dp;
	const double cstep = this->find_cstep(this, pmax, 2 * nkmax);

	if (concentration_to_index(this, pmax, cstep) >= en_cnt)
	{
		fprintf(stderr, "Concentration too high too compute Fermi energies accurately\n");
		fprintf(stderr, "pmax == %g\n", pmax);
		fprintf(stderr, "nkmax == %d\n", nkmax);
		abort();
	}

	for (int nD = 0; nD < Dlen; nD++) {
		this->D = Dmin + nD * dD;
		double extr_edgeen = generate_energies(this, cstep, nkmax, energies);
		qsort(energies, en_cnt, sizeof(double), this->inverted_fermi ? cmp_double_desc : cmp_double_asc);

		for (int np = 0; np < plen; np++) {
			double p = pmin + np * dp;
			int pos = concentration_to_index(this, p, cstep);
			if (pos < 0) {
				pos = 0;
			}
			assert(pos < en_cnt);
			// check if the Fermi surface touches the lattice edge
			if (this->inverted_fermi ? energies[pos] <= extr_edgeen : energies[pos] >= extr_edgeen) {
				fprintf(stderr, "Fermi level touched the edge of the EF lattice!\n");
				fflush(stderr);
			}
			double currEF = energies[pos];
			EFtab[np*Dlen + nD] = currEF;
			if (en_cnt > 1)
			{
				double error;
				if (pos == 0 && en_cnt > 1)
				{
					error = fabs(energies[1] - energies[0]);
				}
				else if (pos == en_cnt - 1 && pos > 0)
				{
					error = fabs(energies[pos] - energies[pos - 1]);
				}
				else if (pos > 0 && pos < en_cnt - 1)
				{
					error = 0.5 * fabs(energies[pos + 1] - energies[pos - 1]);
				}
				printf("Obtained Fermi energy %g with approx. error %g\n", currEF, error);
			}
			else
			{
				printf("Obtained Fermi energy %g, unable to estimate error\n", currEF);
			}
			fflush(stdout);
		}
	}
}

double find_single_EF(ham_struct* this, double p, double D, int nkmax)
{
	int en_cnt = calculate_en_cnt(this, nkmax);
	double *en = alloc_memory(sizeof(double) * en_cnt);
	double EF = find_single_EF_preallocated(this, p, D, nkmax, en);
	free(en);
	return EF;
}


/* It is assumed that energies array has size at least calculate_en_cnt(this,nkmax)*sizeof(double */
double find_single_EF_preallocated(ham_struct* this, double p, double D, int nkmax, double* energies)
{
	double EF;

	find_EF_preallocated(this, p, 0, 1, D, 0, 1, nkmax, energies, &EF);

	return EF;
}

double* find_EF(ham_struct* this, double pmin, double dp, double pmax, double Dmin, double dD, double Dmax, int nkmax)
{
	int en_cnt = calculate_en_cnt(this, nkmax);
	double *en = alloc_memory(sizeof(double) * en_cnt);
	
	int Dlen = count_values(Dmin, Dmax, dD);
	int plen = count_values(pmin, pmax, dp);
	
	double *EFtab = alloc_memory(sizeof(double) * Dlen * plen);

	find_EF_preallocated(this, pmin, dp, plen, Dmin, dD, Dlen, nkmax, en, EFtab);

	free(en);

	return EFtab;
}

/* Create a vanilla ham_struct, ready for action! */
void setup_ham_struct(ham_struct* hs, double MnX, int dim, int periodic_bc_dim)
{
	hs->dim = dim;
	hs->work_matrix = alloc_memory(sizeof(doublecomplex) * dim * dim);
	hs->work_vector = alloc_memory(sizeof(double) * dim);
	hs->S_z = NULL;
	hs->S_x = NULL;
	hs->S_y = NULL;
	hs->S_minus = NULL;
	hs->S_plus = NULL;
	hs->magdir[0] = 0;
	hs->magdir[1] = 0;
	hs->magdir[2] = 1;
	hs->inverted_fermi = 0;
	hs->gen_ham = NULL;
	hs->destroy = destroy_ham_struct;
	hs->L = 0;
	hs->volume_factor = 0;
	hs->volume_factor_12 = 0;
	hs->MnX = MnX;
	hs->cosine_product = 0;
	hs->cosine_product_12 = 0;
	hs->width = 0;
	switch(periodic_bc_dim)
	{
		case 2:
			hs->sum2integral_prefactor = sum2integral_prefactor_2d;
			break;
		case 3:
			hs->sum2integral_prefactor = sum2integral_prefactor_3d;
			break;
		default:
			fprintf(stderr, "Unsupported dimension of periodic boundary constraints: %d\n", periodic_bc_dim);
	}
	hs->periodic_bc_dim = periodic_bc_dim;
	memset(hs->kx_vect, 0, sizeof(double) * 3);
	memset(hs->ky_vect, 0, sizeof(double) * 3);
	memset(hs->kz_vect, 0, sizeof(double) * 3);
	memset(hs->bvect_lengths, 0, sizeof(double) * 3);
}

/* Destroy a vanilla ham_struct */
void destroy_ham_struct(ham_struct* this)
{
	free_if_not_null(this->work_matrix);
	free_if_not_null(this->work_vector);
	free_if_not_null((void*) this->S_x);
	free_if_not_null((void*)this->S_y);
	free_if_not_null((void*)this->S_z);
	free_if_not_null((void*)this->S_minus);
	free_if_not_null((void*)this->S_plus);
	free(this);
}

/* Translate k-vector from program units and Cartesian coordinates to reciprocal basis coordinates */
void translate_into_inner_coords(const ham_struct* this, double kx, double ky, double kz, double* c1, double* c2, double* c3)
{
	(*c1) = this->kx_vect[0] * kx + this->ky_vect[0] * ky + this->kz_vect[0] * kz;
	(*c2) = this->kx_vect[1] * kx + this->ky_vect[1] * ky + this->kz_vect[1] * kz;
	(*c3) = this->kx_vect[2] * kx + this->ky_vect[2] * ky + this->kz_vect[2] * kz;
}

void scan_bands(ham_struct* this, double cx, double cy, double cz, FILE* output)
{
	double kx, ky, kz, r;
	double c1, c2, c3;
	const double dr = 0.001;
	for (r = -1; r <= 1; r+= dr)
	{
		kx = r * cx;
		ky = r * cy;
		kz = r * cz;
		translate_into_inner_coords(this, kx, ky, kz, &c1, &c2, &c3);
		this->gen_ham(this, c1, c2, c3, this->work_matrix);
		diagonalize_hermitian_matrix(this->work_matrix, this->work_vector, this->dim);
		fprintf(output, "%g %g %g", kx / this->L, ky / this->L, kz / this->L);
		for (int k = 0; k < this->dim; k++)
		{
			fprintf(output, " %g", this->work_vector[k]);
			fflush(output);
		}
		fprintf(output, "\n");
	}
}

void scan_bands_sz(ham_struct* this, double cx, double cy, double cz, FILE* output)
{
	double kx, ky, kz, r;
	double c1, c2, c3;
	const double dr = 0.001;
	for (r = -1; r <= 1; r+= dr)
	{
		kx = r * cx;
		ky = r * cy;
		kz = r * cz;
		translate_into_inner_coords(this, kx, ky, kz, &c1, &c2, &c3);
		this->gen_ham(this, c1, c2, c3, this->work_matrix);
		diagonalize_hermitian_matrix(this->work_matrix, this->work_vector, this->dim);
		fprintf(output, "%g %g %g", kx / this->L, ky / this->L, kz / this->L);
		for (int k = 0; k < this->dim; k++)
		{
			double msz = hermitian_matrix_mean_value(this->S_z, &this->work_matrix[k * this->dim], this->dim);
			fprintf(output, " %g", msz);
			fflush(output);
		}
		fprintf(output, "\n");
	}
}

double Delta0(ham_struct* this)
{
	return this->MnX * 4.0 / a0SI / a0SI / a0SI * BetaSI * Spin;
}

