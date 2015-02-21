#include <stdio.h>
#include "kane.h"
#include "clapack.h"
#include "utils.h"
#include "testing.h"

int test_spin_matrix_Sz(void)
{
	doublecomplex Sz_copy[KANE_HAM_DIM][KANE_HAM_DIM];
	double ev[KANE_HAM_DIM];

	initialize_kane_spin_matrices();
	for (int i = 0; i < KANE_HAM_DIM; i++)
	{
		for (int j = 0; j < KANE_HAM_DIM; j++)
		{
			Sz_copy[i][j] = kane_spin_matrices.S_z[i][j];
		}
	}
	for (int i = 0; i < KANE_HAM_DIM; i++)
	{
		for (int j = 0; j < KANE_HAM_DIM; j++)
		{
			if (Sz_copy[i][j].real != Sz_copy[j][i].real || Sz_copy[i][j].imag != Sz_copy[j][i].imag)
			{
				printf("S_z matrix non-hermitian.\n");
				return 0;
			}
		}
	}
	diagonalize_hermitian_matrix(&Sz_copy[0][0], ev, KANE_HAM_DIM);
	for (int i = 0; i < KANE_HAM_DIM; i++)
	{
		printf("%g\n", ev[i]);
	}
	return 1;
}



int main(void)
{
	test_success(test_spin_matrix_Sz(), "TEST_SZ");
}
