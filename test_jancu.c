#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "testing.h"
#include "constants.h"
#include "tightbinding.h"
#include "doublecomplex.h"
#include "utils.h"

int main(void)
{
	ham_struct * hs =  create_tightbinding_bulk_hamiltonian_structure(EXACT, GEOM_XY, JANCU, 0);
	test_int(3, hs->periodic_bc_dim, "TBJANCU_PBCD");
	test_double(0, hs->exx, 1E-12, "TBJANCU_EXX");
	test_double(1/sqrt(2), hs->volume_factor, 1E-6, "TBJANCU_VF");
	test_double(sqrt(3)/2, hs->volume_factor_12, 1E-6, "TBJANCU_VF_2D");
	test_double(2*sqrt(2) / 3 / sqrt(3), hs->cosine_product, 1E-6, "TBJANCU_CP");
	test_double(2.0/3, hs->cosine_product_12, 1E-6, "TBJANCU_CP_2D");
	test_double(0.5633 / 2 / M_PI / sqrt(3), hs->L, 1E-6, "TBJANCU_L");
	double aa = sqrt(3)/2;
	test_double(aa, hs->kx_vect[0], 1E-6, "TBJANCU_KX0");	
	test_double(aa, hs->kx_vect[1], 1E-6, "TBJANCU_KX1");	
	test_double(0, hs->kx_vect[2], 1E-6, "TBJANCU_KX2");	
	test_double(aa, hs->ky_vect[0], 1E-6, "TBJANCU_KY0");	
	test_double(0, hs->ky_vect[1], 1E-6, "TBJANCU_KY1");	
	test_double(aa, hs->ky_vect[2], 1E-6, "TBJANCU_KY2");	
	test_double(0, hs->kz_vect[0], 1E-6, "TBJANCU_KZ0");	
	test_double(aa, hs->kz_vect[1], 1E-6, "TBJANCU_KZ1");	
	test_double(aa, hs->kz_vect[2], 1E-6, "TBJANCU_KZ2");	
	test_double(1, hs->bvect_lengths[0], 1E-6, "TBJANCU_BL0");
	test_double(1, hs->bvect_lengths[1], 1E-6, "TBJANCU_BL1");
	test_double(1, hs->bvect_lengths[2], 1E-6, "TBJANCU_BL2");

	doublecomplex* matrix = hs->work_matrix;
	doublecomplex* matrix2 = initialize_cmplx_matrix(hs->dim);
	
	hs->gen_ham(hs, 0.1, -0.2, 0.7, matrix);

	test_success(is_hermitian(matrix, hs->dim), "TBJANCU_HAM_HERMICITY");

	hs->gen_ham_derivative(hs, 0.1, -0.2, 0.7, 1, 0, 0, matrix);
	test_success(is_hermitian(matrix, hs->dim), "TBJANCU_HAM_DERIVATIVE_HERMICITY");

	double* spin_evs = hs->work_vector;

	diagonalize_hermitian_matrix(hs->S_z, spin_evs, hs->dim);
	printf("Spin S_z matrix eigenvalues:\n");
	for (int i = 0; i < hs->dim; i++)
	{
		printf("Eigenvalue number %d: %g\n", i, spin_evs[i]);
	}

	diagonalize_hermitian_matrix(hs->S_x, spin_evs, hs->dim);
	printf("Spin S_x matrix eigenvalues:\n");
	for (int i = 0; i < hs->dim; i++)
	{
		printf("Eigenvalue number %d: %g\n", i, spin_evs[i]);
	}

	diagonalize_hermitian_matrix(hs->S_y, spin_evs, hs->dim);
	printf("Spin S_y matrix eigenvalues:\n");
	for (int i = 0; i < hs->dim; i++)
	{
		printf("Eigenvalue number %d: %g\n", i, spin_evs[i]);
	}

	hs->D = -1;
	hs->gen_ham(hs, 0.1, -0.2, 0.7, matrix);
	hs->magdir[2] = -hs->magdir[2];
	hs->D = 1;
	hs->gen_ham(hs, 0.1, -0.2, 0.7, matrix2);

	for (int i=0; i< hs->dim*hs->dim; i++){
		test_doublecomplex(matrix[i], matrix2[i], 1e-10, "TBJANCU_MAGDIR");	
	}
	
	
	hs->destroy(hs);
	free(matrix2);
}

