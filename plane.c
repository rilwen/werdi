#include "hamiltonian.h"
#include "utils.h"
#include "plane.h"

void destroy_plane_ham_struct(plane_ham_struct* this)
{
	free(this->single_plane_matrix);
	free(this->atom_positions);
	destroy_ham_struct(&this->super); // This will free all memory occupied by plane_ham_struct
}

void setup_plane_ham_struct(plane_ham_struct* this, int bulk_dim, int nlayers)
{
	this->super.dim = nlayers * bulk_dim;
	this->super.periodic_bc_dim = 2;
	this->single_plane_matrix = alloc_memory(bulk_dim * bulk_dim * sizeof(doublecomplex));
	this->nlayers = nlayers;
	this->atom_positions = alloc_memory(bulk_dim * 3 * sizeof(double));
	ham_struct* hs = (ham_struct*) this;
	hs->destroy = destroy_plane_ham_struct;
}
