#ifndef _KANE_H
#define _KANE_H
#include "doublecomplex.h"
#include "hamiltonian.h"


// energy gap [eV]
#define Eg 1.519


/* Dimension of the Hamiltonian. */
#define KANE_HAM_DIM 8

/* Band at which p-hole bands start */
#define KANE_P_BANDS_START 0

/* Number of p-hole bands */
#define KANE_P_BANDS_CNT 6

typedef struct
{
    doublecomplex S_z[KANE_HAM_DIM][KANE_HAM_DIM];
}
kane_spin_matrices_struct;

extern kane_spin_matrices_struct kane_spin_matrices;

/* Initialize the constants in Kane Hamiltonian for a given value of free parameter A'. */
void initialize_kane_constants(double Ap);

/* Initialize S_z */
void initialize_kane_spin_matrices(void);

/* Create the structure which generates Kane hamiltonians:
MnX - Mn concentration
ExSI - electric field along x [V/nm]
*/
ham_struct * create_kane_hamiltonian_structure(double MnX, double ExSI);

/* Calculate the derivative of Kane Hamiltonian over kx */
void calculate_kane_hamiltonian_derivative_kx(double kx, double ky, double kz, doublecomplex H[KANE_HAM_DIM][KANE_HAM_DIM]);

/* Calculate the derivative of Kane Hamiltonian over ky */
void calculate_kane_hamiltonian_derivative_ky(double kx, double ky, double kz, doublecomplex H[KANE_HAM_DIM][KANE_HAM_DIM]);

/* Set B parameter to zero. */
void zero_B();

/* Reset parameter C0 */
void setC0(double value);

#endif /* _KANE_H */

