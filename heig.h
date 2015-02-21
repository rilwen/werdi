#ifndef _HEIG_H
#define _HEIG_H

#include "hamiltonian.h"

#define HEIG_HAM_DIM 6
#define HEIG_LENGTH_UNIT 5.291772108E-2 // atomic units (in nm)
#define HEIG_ENERGY_UNIT 27.211384527.2113845 // atomic units (hartree, in eV)

ham_struct* create_heig_hamiltonian_structure(double MnX, double exx);

#endif
