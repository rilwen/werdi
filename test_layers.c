#include <stdio.h>

#include "berry.h"
#include "testing.h"
#include "tightbinding.h"
#include "hamiltonian.h"
#include "plane.h"

double conc_grad(double pos[3])
{
	return pos[0] + pos[1] - pos[2];
}

int main(void)
{
	ham_struct* hs = create_tightbinding_plane_hamiltonian_structure(FIVEPOINT, GEOM_Z, JANCU, 1, 0, 0, 0, NULL, 0);
	test_int(80, hs->dim, "NL1_DIM");
	plane_ham_struct* tbhs = (plane_ham_struct*) hs;
	test_int(1, tbhs->nlayers, "NL1_NL");
	test_int(2, hs->periodic_bc_dim, "NL1_PBCD");
	hs->gen_ham(hs, 0, 0, 0.1, hs->work_matrix);
	test_success(is_hermitian(hs->work_matrix, hs->dim), "NL1_HERM");
	hs->destroy(hs);

	hs = create_tightbinding_plane_hamiltonian_structure(FIVEPOINT, GEOM_Z, JANCU, 2, 0, 0, 0, NULL, 0);
	test_int(2 * 80, hs->dim, "NL2_DIM");
	test_int(2, hs->periodic_bc_dim, "NL2_PBCD");
	tbhs = (plane_ham_struct*) hs;
	test_int(2, tbhs->nlayers, "NL2_NL");
	hs->gen_ham(hs, 0, 0, 0.1, hs->work_matrix);
	test_success(is_hermitian(hs->work_matrix, hs->dim), "NL2_HERM");
	
	hs = create_tightbinding_plane_hamiltonian_structure(FIVEPOINT, GEOM_Z, JANCU, 2, 0, 0, 2, NULL, 0);
	hs->gen_ham(hs, 0, -0.2, 0.1, hs->work_matrix);
	test_success(is_hermitian(hs->work_matrix, hs->dim), "NL2EF_HERM");

	hs = create_tightbinding_plane_hamiltonian_structure(FIVEPOINT, GEOM_Z, JANCU, 2, 0, 0, 2, conc_grad, 0);
	hs->gen_ham(hs, 0, -0.2, 0.1, hs->work_matrix);
	test_success(is_hermitian(hs->work_matrix, hs->dim), "NL2EFCG_HERM");
/*	FILE * out;
	out = fopen("test_layers_scan_bands.dat", "w");
	scan_bands(hs, 1, 0, 0, out);
	fclose(out);
	printf("Saved scan of energy bands for 2 layers with GEOM_Z, JANCU, FIVEPOINT, Ez != 0 and custom conc gradient\n");
	fflush(stdout);*/
	hs->destroy(hs);

	hs = create_tightbinding_plane_hamiltonian_structure(FIVEPOINT, GEOM_Z, JANCU, 2, 0, 0, 0, NULL, 0);
	berry_calculator* bc = create_berry_calculator(hs);
	FILE* out = fopen("test_layers_scan_curvs_jancu_z_fivepoint_nl2.dat", "w");
	scan_curvatures(bc, 1, 0, 0, out);
	fclose(out);

	bc->destroy(bc);
	hs->destroy(hs);

	return 0;
}
