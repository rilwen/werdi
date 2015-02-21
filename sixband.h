#ifndef _SIXBAND_H
#define _SIXBAND_H
#include "doublecomplex.h"
#include "hamiltonian.h"


/* Dimension of the Hamiltonian. */
#define SIXBAND_HAM_DIM 6

typedef struct
{
    doublecomplex S_x[SIXBAND_HAM_DIM][SIXBAND_HAM_DIM];
    doublecomplex S_y[SIXBAND_HAM_DIM][SIXBAND_HAM_DIM];
    doublecomplex S_z[SIXBAND_HAM_DIM][SIXBAND_HAM_DIM];
    doublecomplex S_plus[SIXBAND_HAM_DIM][SIXBAND_HAM_DIM];
    doublecomplex S_minus[SIXBAND_HAM_DIM][SIXBAND_HAM_DIM];
}
sixband_spin_matrices_struct;

/* Structure returned by the module */
typedef struct my_ham_struct_6band {
	ham_struct super;
	doublecomplex HBS[SIXBAND_HAM_DIM][SIXBAND_HAM_DIM];
} ham_struct_6band;

extern sixband_spin_matrices_struct sixband_spin_matrices;

/* Initialize S_z */
void initialize_sixband_spin_matrices(void);

/* Create the structure which generates Kane hamiltonians with exx strain */
ham_struct* create_sixband_hamiltonian_structure(double MnX, double exx);

#endif /* _SIXBAND_H */
