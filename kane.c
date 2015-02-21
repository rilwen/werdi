#include <math.h>
#include <stdio.h>
#include "doublecomplex.h"
#include "kane.h"
#include "constants.h"
#include "utils.h"
#include "hamiltonian.h"
#include "kp.h"

#define DEBUG 0

kane_spin_matrices_struct kane_spin_matrices;

struct
{
    /* Parameters from Ostromek PRB54 (1996), eqns. (1.4a)-(1.4j) */
    double Ap;
    double Ec;
    double P0;
    double B;
    double Ev;
    double g1;
    double g2;
    double g3;
    double C0;
}
kane_hamiltonian_constants;

/* Initialize the constants in Kane Hamiltonian for a given value of free parameter A'. */
/* A' is assumed to be in SI units [eV A^2] */
/* Taken from Ostromek PRB54 (1996) */
void initialize_kane_constants(double Ap)
{   
	const double L = 2.7604; /* length unit in angstrems */
	/*    double ApSI = Ap * 2.7604 * 2.7604;  A' in eV * Angstrem^2 */
	kane_hamiltonian_constants.Ap = Ap / L / L;

	double g1l = 6.6723 + 1.9027E-4 * Ap;
	double g2l = 1.8655 + 9.4762E-5 * Ap;
	double g3l = 2.6695 + 9.4762E-5 * Ap;
	double Ep = 22.827 - 0.42760 * Ap; 
	kane_hamiltonian_constants.B = (41.90 + 0.18978 * Ap) / L / L;
	kane_hamiltonian_constants.C0 = (0.1257 - 2.0881E-3 * Ap) / L;

	//    double Ep = 22.827 - 0.42760 * Ap; 

		//    kane_hamiltonian_constants.B = 41.90 / 2.7604 / 2.7604 + 0.18978 * Ap;
		//    kane_hamiltonian_constants.C0 = (0.12568 - 2.0881E-3 * ApSI) / 2.7604;

		/*    double g1l = 6.85 + 1.9027E-4 * Ap;
		      double g2l = 2.1 + 9.4762E-5 * Ap;
		      double g3l = 2.9 + 9.4762E-5 * Ap;*/
		//    double Ep = 22.827 - 0.42760 * Ap;
		//    kane_hamiltonian_constants.B = (41.90 + 0.18978 * Ap) / L / L;
		//    kane_hamiltonian_constants.B = 0;
		/*B accounts for the inversion asymmetry of the zinc-blende structure (Td
		  point group, Schoen.ies notation!.If B50, then the structure will be
		  that of diamond (Oh point group, Schoen.ies notation!*/

	//    kane_hamiltonian_constants.C0 = 0;

	double frac = Ep / (3 * Eg + DeltaSO);
	kane_hamiltonian_constants.g1 = g1l - frac;
	kane_hamiltonian_constants.g2 = g2l - 0.5 * frac;
	kane_hamiltonian_constants.g3 = g3l - 0.5 * frac;
	kane_hamiltonian_constants.P0 = sqrt(hbar * hbar * Ep / 2 / m0);
	kane_hamiltonian_constants.Ec = Eg;
	kane_hamiltonian_constants.Ev = 0;
	printf("B[SI] == %g, Ap[SI] = %g, gammal=%g %g %g\n gamma=%g %g %g\n frac=%g, Eg=%g\n", kane_hamiltonian_constants.B * L * L,
			Ap, g1l, g2l, g3l,
			kane_hamiltonian_constants.g1, kane_hamiltonian_constants.g2, kane_hamiltonian_constants.g3, frac, Eg);
}

void setC0(double value)
{
	kane_hamiltonian_constants.C0 = value;
}

/* Initialize S_z */
void initialize_kane_spin_matrices(void)
{
    const double sqrt2 =  1.4142135623730950488;
    
    /* We use Fortran addressing: [mathematical column][mathematical row] */
    kane_spin_matrices.S_z[0][0].real = -1;
    kane_spin_matrices.S_z[0][0].imag = 0;
    kane_spin_matrices.S_z[0][1].real = 0;
    kane_spin_matrices.S_z[0][1].imag = 0;
    kane_spin_matrices.S_z[0][2].real = 0;
    kane_spin_matrices.S_z[0][2].imag = 0;
    kane_spin_matrices.S_z[0][3].real = 0;
    kane_spin_matrices.S_z[0][3].imag = 0;
    kane_spin_matrices.S_z[0][4].real = 0;
    kane_spin_matrices.S_z[0][4].imag = 0;
    kane_spin_matrices.S_z[0][5].real = 0;
    kane_spin_matrices.S_z[0][5].imag = 0;
    kane_spin_matrices.S_z[0][6].real = 0;
    kane_spin_matrices.S_z[0][6].imag = 0;
    kane_spin_matrices.S_z[0][7].real = 0;
    kane_spin_matrices.S_z[0][7].imag = 0;
    
    kane_spin_matrices.S_z[1][1].real = 1;
    kane_spin_matrices.S_z[1][1].imag = 0;
    kane_spin_matrices.S_z[1][2].real = 0;
    kane_spin_matrices.S_z[1][2].imag = 0;
    kane_spin_matrices.S_z[1][3].real = 0;
    kane_spin_matrices.S_z[1][3].imag = 0;
    kane_spin_matrices.S_z[1][4].real = 0;
    kane_spin_matrices.S_z[1][4].imag = 0;
    kane_spin_matrices.S_z[1][5].real = 0;
    kane_spin_matrices.S_z[1][5].imag = 0;
    kane_spin_matrices.S_z[1][6].real = 0;
    kane_spin_matrices.S_z[1][6].imag = 0;
    kane_spin_matrices.S_z[1][7].real = 0;
    kane_spin_matrices.S_z[1][7].imag = 0;
    
    kane_spin_matrices.S_z[2][2].real = 1.0 / 3.0;
    kane_spin_matrices.S_z[2][2].imag = 0;
    kane_spin_matrices.S_z[2][3].real = 0;
    kane_spin_matrices.S_z[2][3].imag = 0;
    kane_spin_matrices.S_z[2][4].real = 0;
    kane_spin_matrices.S_z[2][4].imag = 0;
    kane_spin_matrices.S_z[2][5].real = 0;
    kane_spin_matrices.S_z[2][5].imag = 0;    
    kane_spin_matrices.S_z[2][6].real = 0;
    kane_spin_matrices.S_z[2][6].imag = 0;
    kane_spin_matrices.S_z[2][7].real = - 2.0 * sqrt2 / 3.0;
    kane_spin_matrices.S_z[2][7].imag = 0;
    
    kane_spin_matrices.S_z[3][3].real = 1;
    kane_spin_matrices.S_z[3][3].imag = 0;
    kane_spin_matrices.S_z[3][4].real = 0;
    kane_spin_matrices.S_z[3][4].imag = 0;
    kane_spin_matrices.S_z[3][5].real = 0;
    kane_spin_matrices.S_z[3][5].imag = 0;
    kane_spin_matrices.S_z[3][6].real = 0;
    kane_spin_matrices.S_z[3][6].imag = 0;
    kane_spin_matrices.S_z[3][7].real = 0;
    kane_spin_matrices.S_z[3][7].imag = 0;
    
    kane_spin_matrices.S_z[4][4].real = -1;
    kane_spin_matrices.S_z[4][4].imag = 0;
    kane_spin_matrices.S_z[4][5].real = 0;
    kane_spin_matrices.S_z[4][5].imag = 0;
    kane_spin_matrices.S_z[4][6].imag = 0;
    kane_spin_matrices.S_z[4][7].real = 0;
    kane_spin_matrices.S_z[4][7].imag = 0;
    
    kane_spin_matrices.S_z[5][5].real = - 1.0 / 3.0;
    kane_spin_matrices.S_z[5][5].imag = 0;
    kane_spin_matrices.S_z[5][6].real = - 2.0 * sqrt2 / 3.0;
    kane_spin_matrices.S_z[5][6].imag = 0;
    kane_spin_matrices.S_z[5][7].real = 0;
    kane_spin_matrices.S_z[5][7].imag = 0;
    
    kane_spin_matrices.S_z[6][6].real = 1.0 / 3.0;
    kane_spin_matrices.S_z[6][6].imag = 0;
    kane_spin_matrices.S_z[6][7].real = 0;
    kane_spin_matrices.S_z[6][7].imag = 0;
    
    kane_spin_matrices.S_z[7][7].real = - 1.0 / 3.0;
    kane_spin_matrices.S_z[7][7].imag = 0;
    
    /* Fill the upper half and scale. */
    for (int n = 0; n < KANE_HAM_DIM; n++)
    {
	    /* Divide by 2 to make eigenvalues +/- 0.5 */
	    kane_spin_matrices.S_z[n][n].real /= 2.0;
	    kane_spin_matrices.S_z[n][n].imag /= 2.0;
	    for (int m = 0; m < n; m++)
	    {
		    /* Divide by 2 to make eigenvalues +/- 0.5 */
		    kane_spin_matrices.S_z[m][n].real /= 2.0;
		    kane_spin_matrices.S_z[m][n].imag /= 2.0;
		    kane_spin_matrices.S_z[n][m].real = kane_spin_matrices.S_z[m][n].real;
		    kane_spin_matrices.S_z[n][m].imag = -kane_spin_matrices.S_z[m][n].imag;
	    }
    }
}

/* Calculate Kane Hamiltonian from Ostromek PRB54 (1996) */
void calculate_kane_hamiltonian_impl(double kx, double ky, double kz, double Delta, double Ex, doublecomplex H[KANE_HAM_DIM][KANE_HAM_DIM])
{
//    const double sxx = 0.20915456788495354460; // <s|x|px> in program units, == sqrt(3) * <r_s> / 3 for r_s = 0.1 nm; OBSOLETE!!!
    const double sxx = 2.22409; // Calculated by comparing commutators, see nAHE.c in newprogs
    const double sqrt2 =  1.4142135623730950488;
    const double sqrt3 =  1.73205080756887729353;
    const double sqrt6 =  2.4494897427831780982;
    double knormsq = kx * kx + ky * ky + kz * kz;
    double hbarsq2m = hbar * hbar / 2.0 / m0;
    double A = kane_hamiltonian_constants.Ec + (kane_hamiltonian_constants.Ap + hbarsq2m) * knormsq;
    double U = kane_hamiltonian_constants.P0 * kz /  sqrt3;
    double Vre = kane_hamiltonian_constants.P0 * kx / sqrt6;
    double Vim = - kane_hamiltonian_constants.P0 * ky / sqrt6;
    double Wim = kane_hamiltonian_constants.B * kx * ky / sqrt3;
    double Tre = kane_hamiltonian_constants.B * kz * kx / sqrt6;
    double Tim = kane_hamiltonian_constants.B * kz * ky / sqrt6;
    double P = - kane_hamiltonian_constants.Ev + kane_hamiltonian_constants.g1 * hbarsq2m * knormsq;
    double Q = kane_hamiltonian_constants.g2 * hbarsq2m * (knormsq - 3 * kz * kz);
    double Rre = - sqrt3 * hbarsq2m * kane_hamiltonian_constants.g2 * (kx * kx - ky * ky);
    double Rim = 2 * sqrt3 * hbarsq2m * kane_hamiltonian_constants.g3 * kx * ky;
    double Sre = 2 * sqrt3 * kane_hamiltonian_constants.g3 * hbarsq2m * kz * kx;
    double Sim = - 2 * sqrt3 * kane_hamiltonian_constants.g3 * hbarsq2m * kz * ky;
    double Z = kane_hamiltonian_constants.Ev - DeltaSO - kane_hamiltonian_constants.g1 * hbarsq2m * knormsq;
#if DEBUG == 1
    printf("Ap == %g\n", kane_hamiltonian_constants.Ap);
    printf("B == %g\n", kane_hamiltonian_constants.B);
    printf("Ev == %g\n", kane_hamiltonian_constants.Ev);
    printf("P0 == %g\n", kane_hamiltonian_constants.P0);
#endif

    /* We use Fortran addressing: [mathematical column][mathematical row] */
    H[0][0].real = A;
    H[0][0].imag = 0;
    H[0][1].real = 0;
    H[0][1].imag = 0;
    H[0][2].real = Tre + Vre;
    H[0][2].imag = Tim + Vim + Ex * sxx / sqrt6; // times - Ex for the electron hamiltonian
    H[0][3].real = 0;
    H[0][3].imag = 0;
    H[0][4].real = sqrt3 * (Vre - Tre);
    H[0][4].imag = sqrt3 * (Tim - Vim) + Ex * sxx /sqrt2;
    H[0][5].real = - sqrt2 * U;
    H[0][5].imag = - sqrt2 * Wim;
    H[0][6].real = - U;
    H[0][6].imag = - Wim;
    H[0][7].real = sqrt2 * (Tre + Vre);
    H[0][7].imag = sqrt2 * (Tim + Vim) + Ex * sxx / sqrt3;

    H[1][1].real = A;
    H[1][1].imag = 0;
    H[1][2].real = - sqrt2 * U;
    H[1][2].imag = - sqrt2 * Wim;
    H[1][3].real = - sqrt3 * (Tre + Vre);
    H[1][3].imag = - sqrt3 * (Tim + Vim) - Ex * sxx / sqrt2;
    H[1][4].real = 0;
    H[1][4].imag = 0;
    H[1][5].real = Tre - Vre;
    H[1][5].imag = Vim - Tim - Ex * sxx / sqrt6;
    H[1][6].real = sqrt2 * (Vre - Tre);
    H[1][6].imag = sqrt2 * (Tim - Vim) + Ex * sxx / sqrt3;
    H[1][7].real = U;
    H[1][7].imag = Wim;

    H[2][2].real = - P + Q;
    H[2][2].imag = 0;
    H[2][3].real = - Sre;
    H[2][3].imag = - Sim;
    H[2][4].real = Rre;
    H[2][4].imag = - Rim;
    H[2][5].real = 0;
    H[2][5].imag = 0;
    H[2][6].real = sqrt3 * Sre / sqrt2;
    H[2][6].imag = - sqrt3 * Sim / sqrt2;
    H[2][7].real = - sqrt2 * Q;
    H[2][7].imag = 0;

    H[3][3].real = - P - Q;
    H[3][3].imag = 0;
    H[3][4].real = 0;
    H[3][4].imag = 0;
    H[3][5].real = Rre;
    H[3][5].imag = - Rim;
    H[3][6].real = - sqrt2 * Rre;
    H[3][6].imag = sqrt2 * Rim;
    H[3][7].real = Sre / sqrt2;
    H[3][7].imag = - Sim / sqrt2;

    H[4][4].real = - P - Q;
    H[4][4].imag = 0;
    H[4][5].real = Sre;
    H[4][5].imag = Sim;
    H[4][6].real = Sre / sqrt2;
    H[4][6].imag = Sim / sqrt2;
    H[4][7].real = sqrt2 * Rre;
    H[4][7].imag = sqrt2 * Rim;

    H[5][5].real = - P + Q;
    H[5][5].imag = 0;
    H[5][6].real = sqrt2 * Q;
    H[5][6].imag = 0;
    H[5][7].real = sqrt3 * Sre / sqrt2;
    H[5][7].imag = sqrt3 * Sim / sqrt2;

    H[6][6].real = Z;
    H[6][6].imag = 0;
    H[6][7].real = 0;
    H[6][7].imag = 0;

    H[7][7].real = Z;
    H[7][7].imag = 0;

    /* Add spin-orbit coupling (to the lower half only) */
    H[0][2].real += - kane_hamiltonian_constants.C0 * kx / sqrt2;
    H[0][2].imag += kane_hamiltonian_constants.C0 * ky / sqrt2;
    H[0][3].real += 0;
    H[0][3].imag += 0;
    H[0][4].real += - sqrt3 * kane_hamiltonian_constants.C0 * kx / sqrt2;
    H[0][4].imag += - sqrt3 * kane_hamiltonian_constants.C0 * ky / sqrt2;
    H[0][5].real += kane_hamiltonian_constants.C0 * sqrt2 * kz;
    H[0][5].imag += 0;
    H[0][6].real += - 2 * kane_hamiltonian_constants.C0 * kz;
    H[0][6].imag += 0;
    H[0][7].real += 2 * kane_hamiltonian_constants.C0 * kx;
    H[0][7].imag += - 2 * kane_hamiltonian_constants.C0 * ky;

    H[1][2].real += kane_hamiltonian_constants.C0 * sqrt2 * kz;
    H[1][2].imag += 0;
    H[1][3].real += sqrt3 * kane_hamiltonian_constants.C0 * kx / sqrt2;
    H[1][3].imag += - sqrt3 * kane_hamiltonian_constants.C0 * ky / sqrt2;
    H[1][4].real += 0;
    H[1][4].imag += 0;
    H[1][5].real += kane_hamiltonian_constants.C0 * kx / sqrt2;
    H[1][5].imag += kane_hamiltonian_constants.C0 * ky / sqrt2;
    H[1][6].real += 2 * kane_hamiltonian_constants.C0 * kx;
    H[1][6].imag += 2 * kane_hamiltonian_constants.C0 * ky;
    H[1][7].real += 2 * kane_hamiltonian_constants.C0 * kz;
    H[1][7].imag += 0;

    /* Fill the upper half */
    for (int n = 0; n < KANE_HAM_DIM; n++)
    {
        for (int m = 0; m < n; m++)
        {
            H[n][m].real = H[m][n].real;
            H[n][m].imag = - H[m][n].imag;
        }
    }

    /* Add the Zeeman splitting */
    if (Delta != 0)
    {
        for (int m = 0; m < KANE_HAM_DIM; m++)
        {
            for (int n = 0; n < KANE_HAM_DIM; n++)
            {
                H[m][n].real -= Delta * kane_spin_matrices.S_z[m][n].real; // -Delta sz for holes spins
                H[m][n].imag -= Delta * kane_spin_matrices.S_z[m][n].imag;
            }
        }
    }

#if DEBUG == 1
    printf("Debug mode on\n");
    for (int i = 0; i < KANE_HAM_DIM; i++)
    {
        for (int j = 0; j <= i; j++)
        {
            double dr = H[i][j].real - H[j][i].real;
            double di = H[i][j].imag + H[j][i].imag;
            if (fabs(dr) > 1E-16 || fabs(di) > 1E-16) {
                printf("Hamiltonian non-hermicity detected for H[%d][%d]\n", i, j);
            }
        }
    }
#endif
}

/* Calculate the derivative of Kane Hamiltonian over k along direction (cx, cy, cz).*/
void calculate_kane_hamiltonian_derivative_impl(double kx, double ky, double kz, double cx, double cy, double cz, doublecomplex H[KANE_HAM_DIM][KANE_HAM_DIM])
{
    const double sqrt2 =  1.4142135623730950488;
    const double sqrt3 =  1.73205080756887729353;
    const double sqrt6 =  2.4494897427831780982;
    double knormsq = 2 * (cx * kx + cy * ky + cz * kz);
    double hbarsq2m = hbar * hbar / 2 / m0;
    double A = (kane_hamiltonian_constants.Ap + hbarsq2m) * knormsq;
    double U = kane_hamiltonian_constants.P0 * cz /  sqrt3;
    double Vre = kane_hamiltonian_constants.P0 * cx / sqrt6;
    double Vim = - kane_hamiltonian_constants.P0 * cy / sqrt6;
    double Wim = kane_hamiltonian_constants.B * (cx * ky + kx * cy) / sqrt3;
    double Tre = kane_hamiltonian_constants.B * (cz * kx + kz * cx) / sqrt6;
    double Tim = kane_hamiltonian_constants.B * (cz * ky + kz * cy) / sqrt6;
    double P = kane_hamiltonian_constants.g1 * hbarsq2m * knormsq;
    double Q = kane_hamiltonian_constants.g2 * hbarsq2m * (knormsq - 6 * cz * kz);
    double Rre = - sqrt3 * hbarsq2m * kane_hamiltonian_constants.g2 * 2 * (cx * kx - cy * ky);
    double Rim = 2 * sqrt3 * hbarsq2m * kane_hamiltonian_constants.g3 * (cx * ky + kx * cy);
    double Sre = 2 * sqrt3 * kane_hamiltonian_constants.g3 * hbarsq2m * (cz * kx + kz * cx);
    double Sim = - 2 * sqrt3 * kane_hamiltonian_constants.g3 * hbarsq2m * (cz * ky + kz * cy);
    double Z = - kane_hamiltonian_constants.g1 * hbarsq2m * knormsq;

    /* We use Fortran addressing: [mathematical column][mathematical row] */
    H[0][0].real = A;
    H[0][0].imag = 0;
    H[0][1].real = 0;
    H[0][1].imag = 0;
    H[0][2].real = Tre + Vre;
    H[0][2].imag = Tim + Vim;
    H[0][3].real = 0;
    H[0][3].imag = 0;
    H[0][4].real = sqrt3 * (Vre - Tre);
    H[0][4].imag = sqrt3 * (Tim - Vim);
    H[0][5].real = - sqrt2 * U;
    H[0][5].imag = - sqrt2 * Wim;
    H[0][6].real = - U;
    H[0][6].imag = - Wim;
    H[0][7].real = sqrt2 * (Tre + Vre);
    H[0][7].imag = sqrt2 * (Tim + Vim);
    
    H[1][1].real = A;
    H[1][1].imag = 0;
    H[1][2].real = - sqrt2 * U;
    H[1][2].imag = - sqrt2 * Wim;
    H[1][3].real = - sqrt3 * (Tre + Vre);
    H[1][3].imag = - sqrt3 * (Tim + Vim);
    H[1][4].real = 0;
    H[1][4].imag = 0;
    H[1][5].real = Tre - Vre;
    H[1][5].imag = Vim - Tim;
    H[1][6].real = sqrt2 * (Vre - Tre);
    H[1][6].imag = sqrt2 * (Tim - Vim);
    H[1][7].real = U;
    H[1][7].imag = Wim;
    
    H[2][2].real = - P + Q;
    H[2][2].imag = 0;
    H[2][3].real = - Sre;
    H[2][3].imag = - Sim;
    H[2][4].real = Rre;
    H[2][4].imag = - Rim;
    H[2][5].real = 0;
    H[2][5].imag = 0;
    H[2][6].real = sqrt3 * Sre / sqrt2;
    H[2][6].imag = - sqrt3 * Sim / sqrt2;
    H[2][7].real = - sqrt2 * Q;
    H[2][7].imag = 0;
    
    H[3][3].real = - P - Q;
    H[3][3].imag = 0;
    H[3][4].real = 0;
    H[3][4].imag = 0;
    H[3][5].real = Rre;
    H[3][5].imag = - Rim;
    H[3][6].real = - sqrt2 * Rre;
    H[3][6].imag = sqrt2 * Rim;
    H[3][7].real = Sre / sqrt2;
    H[3][7].imag = - Sim / sqrt2;
    
    H[4][4].real = - P - Q;
    H[4][4].imag = 0;
    H[4][5].real = Sre;
    H[4][5].imag = Sim;
    H[4][6].real = Sre / sqrt2;
    H[4][6].imag = Sim / sqrt2;
    H[4][7].real = sqrt2 * Rre;
    H[4][7].imag = sqrt2 * Rim;
    
    H[5][5].real = - P + Q;
    H[5][5].imag = 0;
    H[5][6].real = sqrt2 * Q;
    H[5][6].imag = 0;
    H[5][7].real = sqrt3 * Sre / sqrt2;
    H[5][7].imag = sqrt3 * Sim / sqrt2;
    
    H[6][6].real = Z;
    H[6][6].imag = 0;
    H[6][7].real = 0;
    H[6][7].imag = 0;
    
    H[7][7].real = Z;
    H[7][7].imag = 0;
    
    /* Add derivative of spin-orbit coupling (to the lower half only) */
    H[0][2].real += - kane_hamiltonian_constants.C0 * cx / sqrt2;
    H[0][2].imag += kane_hamiltonian_constants.C0 * cy / sqrt2;
    H[0][3].real += 0;
    H[0][3].imag += 0;
    H[0][4].real += - sqrt3 * kane_hamiltonian_constants.C0 * cx / sqrt2;
    H[0][4].imag += - sqrt3 * kane_hamiltonian_constants.C0 * cy / sqrt2;
    H[0][5].real += kane_hamiltonian_constants.C0 * sqrt2 * cz;
    H[0][5].imag += 0;
    H[0][6].real += - 2 * kane_hamiltonian_constants.C0 * cz;
    H[0][6].imag += 0;
    H[0][7].real += 2 * kane_hamiltonian_constants.C0 * cx;
    H[0][7].imag += - 2 * kane_hamiltonian_constants.C0 * cy;

    H[1][2].real += kane_hamiltonian_constants.C0 * sqrt2 * cz;
    H[1][2].imag += 0;
    H[1][3].real += sqrt3 * kane_hamiltonian_constants.C0 * cx / sqrt2;
    H[1][3].imag += - sqrt3 * kane_hamiltonian_constants.C0 * cy / sqrt2;
    H[1][4].real += 0;
    H[1][4].imag += 0;
    H[1][5].real += kane_hamiltonian_constants.C0 * cx / sqrt2;
    H[1][5].imag += kane_hamiltonian_constants.C0 * cy / sqrt2;
    H[1][6].real += 2 * kane_hamiltonian_constants.C0 * cx;
    H[1][6].imag += 2 * kane_hamiltonian_constants.C0 * cy;
    H[1][7].real += 2 * kane_hamiltonian_constants.C0 * cz;
    H[1][7].imag += 0;



    /* Fill the upper half */
    for (int n = 0; n < KANE_HAM_DIM; n++)
    {
        for (int m = 0; m < n; m++)
        {
            H[n][m].real = H[m][n].real;
            H[n][m].imag = - H[m][n].imag;
        }
    }
#if DEBUG == 1
    for (int i = 0; i < KANE_HAM_DIM; i++)
    {
        for (int j = 0; j <= i; j++)
        {
            double dr = H[i][j].real - H[j][i].real;
            double di = H[i][j].imag + H[j][i].imag;
            if (fabs(dr) > 1E-16 || fabs(di) > 1E-16) {
                printf("Hamiltonian derivative non-hermicity detected for H[%d][%d]\n", i, j);
            }
        }
    }
#endif
}

void zero_B(void)
{
    kane_hamiltonian_constants.B = 0;
}

/* Calculate the derivative of Kane Hamiltonian over kx */
void calculate_kane_hamiltonian_derivative_kx(double kx, double ky, double kz, doublecomplex H[KANE_HAM_DIM][KANE_HAM_DIM])
{
    calculate_kane_hamiltonian_derivative_impl(kx, ky, kz, 1, 0, 0, H);
}

/* Calculate the derivative of Kane Hamiltonian over ky */
void calculate_kane_hamiltonian_derivative_ky(double kx, double ky, double kz, doublecomplex H[KANE_HAM_DIM][KANE_HAM_DIM])
{
    calculate_kane_hamiltonian_derivative_impl(kx, ky, kz, 0, 1, 0, H);
}

void calculate_kane_hamiltonian(ham_struct* this, double kx, double ky, double kz, doublecomplex matrix[][])
{
	if (this->magdir[0]!=0 || this->magdir[1]!=0 || this->magdir[2]!=1){
		fprintf(stderr, "Magnetization direction not supported for this model.");
		abort();
	}

	calculate_kane_hamiltonian_impl(kx, ky, kz, this->D, this->Ex * this->L, matrix);
}

void calculate_kane_hamiltonian_derivative(ham_struct* this, double kx, double ky, double kz, double cx, double cy, double cz, doublecomplex matrix[][])
{
	calculate_kane_hamiltonian_derivative_impl(kx, ky, kz, cx, cy, cz, matrix);
}

/* Create the structure which generates Kane hamiltonians */
ham_struct* create_kane_hamiltonian_structure(double MnX, double ExSI)
{
	initialize_kane_spin_matrices();
	ham_struct * hs = alloc_memory(sizeof(ham_struct));
	setup_ham_struct(hs, MnX, KANE_HAM_DIM, 3);
	hs->gen_ham = calculate_kane_hamiltonian;
	hs->gen_ham_derivative = calculate_kane_hamiltonian_derivative;
	hs->p_bands_start = KANE_P_BANDS_START;
	hs->p_bands_cnt = KANE_P_BANDS_CNT;
	hs->dim = KANE_HAM_DIM;
	size_t size = sizeof(doublecomplex) * KANE_HAM_DIM * KANE_HAM_DIM;
	hs->S_z = copy_to_dynamic_array(kane_spin_matrices.S_z, size);
	hs->find_cstep = find_cstep_kp;
	hs->bril_zone_integrand = bril_zone_integrand_kp;
	setup_kp_geometry(hs, DEFAULT_KP_L);

	/*
	* Converted Ex from SI to program units in calculate_kane_hamiltonian.
	*/
	hs->Ex = ExSI;
	
	hs->inverted_fermi = 1;

	return hs;
}

