#ifndef _FOURBAND_H
#define _FOURBAND_H
#include "doublecomplex.h"
#include "hamiltonian.h"


/* Dimension of the Hamiltonian. */
#define FOURBAND_HAM_DIM 4

typedef struct
{
    doublecomplex S_x[FOURBAND_HAM_DIM][FOURBAND_HAM_DIM];
    doublecomplex S_y[FOURBAND_HAM_DIM][FOURBAND_HAM_DIM];
    doublecomplex S_z[FOURBAND_HAM_DIM][FOURBAND_HAM_DIM];
    doublecomplex S_plus[FOURBAND_HAM_DIM][FOURBAND_HAM_DIM];
    doublecomplex S_minus[FOURBAND_HAM_DIM][FOURBAND_HAM_DIM];
}
fourband_spin_matrices_struct;

/* Structure returned by the module */
typedef struct my_ham_struct_4band {
	ham_struct super;
	double exx;
	doublecomplex HBS[FOURBAND_HAM_DIM][FOURBAND_HAM_DIM];
} ham_struct_4band;

extern fourband_spin_matrices_struct fourband_spin_matrices;

/* Initialize S_z */
void initialize_fourband_spin_matrices(void);

/* Create the structure which generates Kane hamiltonians with exx strain */
ham_struct* create_fourband_hamiltonian_structure(double MnX, double exx);

#endif /* _FOURBAND_H */
