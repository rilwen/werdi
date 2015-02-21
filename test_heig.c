#include <stdio.h>
#include <stdlib.h>
#include "heig.h"
#include "hamiltonian.h"
#include "testing.h"
#include "doublecomplex.h"
#include "utils.h"

int main(void)
{
	ham_struct * hs = create_heig_hamiltonian_structure(0.05, 0.001);
	test_double(0.001, hs->exx, 0, "HEIG_EXX");
	test_int(3, hs->periodic_bc_dim, "HEIG_PBCD");

	doublecomplex matrix[HEIG_HAM_DIM][HEIG_HAM_DIM];

	hs->gen_ham(hs, 0.1, -0.2, 0.7, matrix);

	test_success(is_hermitian(matrix, HEIG_HAM_DIM), "HEIG_HAM_HERMICITY");

	hs->gen_ham_derivative(hs, 0.1, -0.2, 0.7, 1, 0, 0, matrix);
	test_success(is_hermitian(matrix, HEIG_HAM_DIM), "HEIG_HAM_DERIVATIVE_HERMICITY");

	double spin_evs[HEIG_HAM_DIM];

	diagonalize_hermitian_matrix(hs->S_z, spin_evs, HEIG_HAM_DIM);
	printf("Spin S_z matrix eigenvalues:\n");
	for (int i = 0; i < HEIG_HAM_DIM; i++)
	{
		printf("Eigenvalue number %d: %g\n", i, spin_evs[i]);
	}

/*
	diagonalize_hermitian_matrix(hs->S_x, spin_evs, HEIG_HAM_DIM);
	printf("Spin S_x matrix eigenvalues:\n");
	for (int i = 0; i < HEIG_HAM_DIM; i++)
	{
		printf("Eigenvalue number %d: %g\n", i, spin_evs[i]);
	}

	diagonalize_hermitian_matrix(hs->S_y, spin_evs, HEIG_HAM_DIM);
	printf("Spin S_y matrix eigenvalues:\n");
	for (int i = 0; i < HEIG_HAM_DIM; i++)
	{
		printf("Eigenvalue number %d: %g\n", i, spin_evs[i]);
	}
*/

	
	hs->destroy(hs);
}

