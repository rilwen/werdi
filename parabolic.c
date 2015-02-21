#include <assert.h>
#include <math.h>
#include <string.h>
#include "constants.h"
#include "hamiltonian.h"
#include "parabolic.h"
#include "utils.h"
#include "kp.h"
#include "lattice.h"

doublecomplex parabolic_S_z[PARABOLIC_HAM_DIM][PARABOLIC_HAM_DIM];

static void initialize_S_z(void)
{
	memset(parabolic_S_z, 0, sizeof(doublecomplex) * PARABOLIC_HAM_DIM * PARABOLIC_HAM_DIM);
	assert(PARABOLIC_HAM_DIM % 2 == 0);
	int halfdim = PARABOLIC_HAM_DIM / 2;
	for (int i = 0; i < halfdim; i++)
	{
		parabolic_S_z[i][i] = cmplx(0.5, 0);
		parabolic_S_z[i + halfdim][i + halfdim] = cmplx(-0.5, 0);
	}
}

void calculate_parabolic_hamiltonian(ham_struct * this, double c1, double c2, double c3, doublecomplex* matrix)
{
	// k.p models have cubic cells
	double kx = c1*this->bvect_lengths[0];
	double ky = c2*this->bvect_lengths[1];
	double kz = c3*this->bvect_lengths[2];
	double c = 0.1, dc = 0.02;
	double Ekin = 0.5 * (kx*kx + ky*ky + kz*kz);

	memset(matrix, 0, PARABOLIC_HAM_DIM * PARABOLIC_HAM_DIM * sizeof(doublecomplex));
	for (int i = 0; i < PARABOLIC_HAM_DIM; i++)
	{
		matrix[i * PARABOLIC_HAM_DIM + i] = cmplx(c * Ekin, 0);
		c += dc;
		dc += dc;
	}

	for (int i = 0; i < PARABOLIC_HAM_DIM; i++)
	{
		for (int j = 0; j < PARABOLIC_HAM_DIM; j++)
		{
			matrix[i * PARABOLIC_HAM_DIM + j].real += this->D * parabolic_S_z[i][j].real;
			matrix[i * PARABOLIC_HAM_DIM + j].imag += this->D * parabolic_S_z[i][j].imag;
		}
	}
}

void calculate_parabolic_hamiltonian_derivative(ham_struct* this, double c1, double c2, double c3, double d1, double d2, double d3, doublecomplex dHkp[PARABOLIC_HAM_DIM][PARABOLIC_HAM_DIM])
{
	// k.p models have cubic cells
	double kx = c1*this->bvect_lengths[0];
	double ky = c2*this->bvect_lengths[1];
	double kz = c3*this->bvect_lengths[2];
	double cx = d1*this->bvect_lengths[0];
	double cy = d2*this->bvect_lengths[1];
	double cz = d3*this->bvect_lengths[2];
	double c = 0.1, dc = 0.02;
	
	double dEkindk = cx * kx + cy * ky + cz * kz;

	memset(dHkp, 0, PARABOLIC_HAM_DIM * PARABOLIC_HAM_DIM * sizeof(doublecomplex));
	for (int i = 0; i < PARABOLIC_HAM_DIM; i++)
	{
		dHkp[i][i] = cmplx(c * dEkindk, 0);
		c += dc;
		dc += dc;
	}
}

ham_struct* create_parabolic_hamiltonian_structure(void)
{
	initialize_S_z();
	ham_struct* hs = alloc_memory(sizeof(ham_struct));
	setup_ham_struct(hs, 0.05, PARABOLIC_HAM_DIM, 3);
	hs->gen_ham = calculate_parabolic_hamiltonian;
	hs->gen_ham_derivative = calculate_parabolic_hamiltonian_derivative;
	hs->p_bands_start = 0;
	hs->p_bands_cnt = PARABOLIC_HAM_DIM;
	size_t size = sizeof(doublecomplex) * PARABOLIC_HAM_DIM * PARABOLIC_HAM_DIM;
	hs->S_z = copy_to_dynamic_array(parabolic_S_z, size);
	hs->dim = PARABOLIC_HAM_DIM;
	hs->find_cstep = find_cstep_kp;
	hs->bril_zone_integrand = bril_zone_integrand_kp;
	hs->exx = 0;
	setup_kp_geometry(hs, DEFAULT_KP_L);

	return hs;
}
