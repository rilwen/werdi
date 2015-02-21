#ifndef _PLANE_H
#define _PLANE_H

#include "hamiltonian.h"

typedef struct my_plane_ham_struct {
        ham_struct super;
        int nlayers; // number of planes (0 for bulk)
        doublecomplex* single_plane_matrix; // work matrix for single plane or btw two single planes
        double electric_field[3]; // Electric field vector in program units
        double* atom_positions; // work vector for atomic positions
        double (*concentration_gradient)(double* pos); // function which calculates concentration gradient to be used for calculations of electric potential, if NULL then it is assumed to be always 0; it takes a pointer to 3-element arrya of coordinates (x,y,z)
} plane_ham_struct;

void destroy_plane_ham_struct(plane_ham_struct* phs);

void setup_plane_ham_struct(plane_ham_struct* this, int bulk_dim, int nlayers);

#endif
