#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include "doublecomplex.h"
#include "heig.h"
#include "constants.h"
#include "utils.h"
#include "hamiltonian.h"
#include "kp.h"
#include "lattice.h"

void heigf_(double* kx, double* ky, double* kz, double* C4_SI, double* exx, double* BG, double* thetaM, double* phiM, doublecomplex* H);

/* Calculate 6-band Hamiltonian */
void calculate_heig_hamiltonian(ham_struct* this, double c1, double c2, double c3, doublecomplex H[HEIG_HAM_DIM][HEIG_HAM_DIM])
{
	const double BG = - this->D / 6 / 27.2113845; // from eV to hartree; Delta = 6 BG
	const double C4_SI = -2.18e6;	
	const double thetaM = 0;
	const double phiM = 0;
	const double kx = c1 * HEIG_LENGTH_UNIT / DEFAULT_KP_L;
	const double ky = c2 * HEIG_LENGTH_UNIT / DEFAULT_KP_L;
	const double kz = c3 * HEIG_LENGTH_UNIT / DEFAULT_KP_L;
	heigf_(&kx, &ky, &kz, &C4_SI, &this->exx, &BG, &thetaM, &phiM, &H[0][0]);

	for (int i = 0; i < HEIG_HAM_DIM; ++i)
	{
		for (int j = 0; j < HEIG_HAM_DIM; ++j)
		{
			H[i][j].real *= 27.2113845;
			H[i][j].imag *= 27.2113845;
		}
	}
}

void setup_heig_spin_matrices(ham_struct* hs)
{
	doublecomplex m[HEIG_HAM_DIM][HEIG_HAM_DIM];
	const size_t size = sizeof(doublecomplex)*HEIG_HAM_DIM*HEIG_HAM_DIM;
	
	hs->S_x = NULL;
	hs->S_y = NULL;

	memset(m, 0, size);
	// [math column][math row]
	m[0][0] = cmplx(0.5, 0);
	m[1][1] = cmplx(1.0/6, 0);
	m[1][4] = cmplx(-sqrt(2.0)/3, 0);
	m[2][2] = cmplx(-1.0/6, 0);
	m[2][5] = cmplx(-sqrt(2.0)/3, 0);
	m[3][3] = cmplx(-0.5, 0);
	m[4][1] = cmplx(-sqrt(2.0)/3, 0);
	m[4][4] = cmplx(-1.0/6, 0);
	m[5][2] = cmplx(-sqrt(2.0)/3, 0);
	m[5][5] = cmplx(1.0/6, 0);
	hs->S_z = copy_to_dynamic_array(m, size);
	
	hs->S_minus = NULL;
	hs->S_plus = NULL;
}

/* Create the structure which generates 6-band hamiltonians */
ham_struct* create_heig_hamiltonian_structure(double MnX, double exx)
{
	ham_struct* hs = alloc_memory(sizeof(ham_struct));
	setup_ham_struct(hs, MnX, HEIG_HAM_DIM, 3);
	hs->gen_ham = calculate_heig_hamiltonian;
	hs->gen_ham_derivative = discretized_5pt_derivative;
	/* Six-band hamiltonian is p-band only */
	hs->p_bands_start = 0;
	hs->p_bands_cnt = HEIG_HAM_DIM;
	hs->find_cstep = find_cstep_kp;
	hs->inverted_fermi = 1;
	hs->bril_zone_integrand = bril_zone_integrand_kp;
	hs->exx = exx;
	setup_kp_geometry(hs, DEFAULT_KP_L);
	setup_heig_spin_matrices(hs);

	return hs;
}

