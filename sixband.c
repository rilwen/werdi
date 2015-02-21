#include <math.h>
#include <stdio.h>
#include <assert.h>
#include "doublecomplex.h"
#include "sixband.h"
#include "constants.h"
#include "utils.h"
#include "hamiltonian.h"
#include "kp.h"
#include "lattice.h"

#define DEBUG 0

sixband_spin_matrices_struct sixband_spin_matrices;

const double Dso=0.341;

/* Initialize S_x, S_y and S_z */
void initialize_sixband_spin_matrices(void)
{
	const double sqrt3 = 1.73205080756887729353;
	const double sqrt2 = 1.4142135623730950488;
	
	// addressing: sixband_spin_matrices.S_x[mathematical column][mathematical row]
	sixband_spin_matrices.S_x[0][0].real = 0.0;
	sixband_spin_matrices.S_x[0][1].real = 0.0;
	sixband_spin_matrices.S_x[0][2].real = 1.0/(2*sqrt3);
	sixband_spin_matrices.S_x[0][3].real = 0.0;
	sixband_spin_matrices.S_x[0][4].real = 1.0/(sqrt2*sqrt3);
	sixband_spin_matrices.S_x[0][5].real = 0.0;
	
	sixband_spin_matrices.S_x[1][0].real = 0.0;
	sixband_spin_matrices.S_x[1][1].real = 0.0;
	sixband_spin_matrices.S_x[1][2].real = 1.0/3;
	sixband_spin_matrices.S_x[1][3].real = 1.0/(2*sqrt3);
	sixband_spin_matrices.S_x[1][4].real = -1.0/(3*sqrt2);
	sixband_spin_matrices.S_x[1][5].real = 0.0;
	
	sixband_spin_matrices.S_x[2][0].real = 1.0/(2*sqrt3);
	sixband_spin_matrices.S_x[2][1].real = 1.0/3;
	sixband_spin_matrices.S_x[2][2].real = 0.0;
	sixband_spin_matrices.S_x[2][3].real = 0.0;
	sixband_spin_matrices.S_x[2][4].real = 0.0;
	sixband_spin_matrices.S_x[2][5].real = 1.0/(3*sqrt2);
	
	sixband_spin_matrices.S_x[3][0].real = 0.0;
	sixband_spin_matrices.S_x[3][1].real = 1.0/(2*sqrt3);
	sixband_spin_matrices.S_x[3][2].real = 0.0;
	sixband_spin_matrices.S_x[3][3].real = 0.0;
	sixband_spin_matrices.S_x[3][4].real = 0.0;
	sixband_spin_matrices.S_x[3][5].real = -1.0/(sqrt2*sqrt3);
	
	sixband_spin_matrices.S_x[4][0].real = 1.0/(sqrt2*sqrt3);
	sixband_spin_matrices.S_x[4][1].real = -1.0/(3*sqrt2);
	sixband_spin_matrices.S_x[4][2].real = 0.0;
	sixband_spin_matrices.S_x[4][3].real = 0.0;
	sixband_spin_matrices.S_x[4][4].real = 0.0;
	sixband_spin_matrices.S_x[4][5].real = -1.0/6;
	
	sixband_spin_matrices.S_x[5][0].real = 0.0;
	sixband_spin_matrices.S_x[5][1].real = 0.0;
	sixband_spin_matrices.S_x[5][2].real = 1.0/(3*sqrt2);
	sixband_spin_matrices.S_x[5][3].real = -1.0/(sqrt2*sqrt3);
	sixband_spin_matrices.S_x[5][4].real = -1.0/6;
	sixband_spin_matrices.S_x[5][5].real = 0.0;

	sixband_spin_matrices.S_y[0][0].imag = 0.0;
	sixband_spin_matrices.S_y[0][1].imag = 0.0;
	sixband_spin_matrices.S_y[0][2].imag = 1.0/(2*sqrt3);
	sixband_spin_matrices.S_y[0][3].imag = 0.0;
	sixband_spin_matrices.S_y[0][4].imag = 1.0/(sqrt2*sqrt3);
	sixband_spin_matrices.S_y[0][5].imag = 0.0;
	
	sixband_spin_matrices.S_y[1][0].imag = 0.0;
	sixband_spin_matrices.S_y[1][1].imag = 0.0;
	sixband_spin_matrices.S_y[1][2].imag = -1.0/3;
	sixband_spin_matrices.S_y[1][3].imag = 1.0/(2*sqrt3);
	sixband_spin_matrices.S_y[1][4].imag = 1.0/(3*sqrt2);
	sixband_spin_matrices.S_y[1][5].imag = 0.0;
	
	sixband_spin_matrices.S_y[2][0].imag = -1.0/(2*sqrt3);
	sixband_spin_matrices.S_y[2][1].imag = 1.0/3;
	sixband_spin_matrices.S_y[2][2].imag = 0.0;
	sixband_spin_matrices.S_y[2][3].imag = 0.0;
	sixband_spin_matrices.S_y[2][4].imag = 0.0;
	sixband_spin_matrices.S_y[2][5].imag = 1.0/(3*sqrt2);
	
	sixband_spin_matrices.S_y[3][0].imag = 0.0;
	sixband_spin_matrices.S_y[3][1].imag = -1.0/(2*sqrt3);
	sixband_spin_matrices.S_y[3][2].imag = 0.0;
	sixband_spin_matrices.S_y[3][3].imag = 0.0;
	sixband_spin_matrices.S_y[3][4].imag = 0.0;
	sixband_spin_matrices.S_y[3][5].imag = 1.0/(sqrt2*sqrt3);
	
	sixband_spin_matrices.S_y[4][0].imag = -1.0/(sqrt2*sqrt3);
	sixband_spin_matrices.S_y[4][1].imag = -1.0/(3*sqrt2);
	sixband_spin_matrices.S_y[4][2].imag = 0.0;
	sixband_spin_matrices.S_y[4][3].imag = 0.0;
	sixband_spin_matrices.S_y[4][4].imag = 0.0;
	sixband_spin_matrices.S_y[4][5].imag = -1.0/6;
	
	sixband_spin_matrices.S_y[5][0].imag = 0.0;
	sixband_spin_matrices.S_y[5][1].imag = 0.0;
	sixband_spin_matrices.S_y[5][2].imag = -1.0/(3*sqrt2);
	sixband_spin_matrices.S_y[5][3].imag = -1.0/(sqrt2*sqrt3);
	sixband_spin_matrices.S_y[5][4].imag = 1.0/6;
	sixband_spin_matrices.S_y[5][5].imag = 0.0;
	
	
	
	sixband_spin_matrices.S_z[0][0].real = 1.0/2;
	sixband_spin_matrices.S_z[0][1].real = 0.0;
	sixband_spin_matrices.S_z[0][2].real = 0.0;
	sixband_spin_matrices.S_z[0][3].real = 0.0;
	sixband_spin_matrices.S_z[0][4].real = 0.0;
	sixband_spin_matrices.S_z[0][5].real = 0.0;
	
	sixband_spin_matrices.S_z[1][0].real = 0.0;
	sixband_spin_matrices.S_z[1][1].real = -1.0/6;
	sixband_spin_matrices.S_z[1][2].real = 0.0;
	sixband_spin_matrices.S_z[1][3].real = 0.0;
	sixband_spin_matrices.S_z[1][4].real = 0.0;
	sixband_spin_matrices.S_z[1][5].real = -sqrt2/3.0;
	
	sixband_spin_matrices.S_z[2][0].real = 0.0;
	sixband_spin_matrices.S_z[2][1].real = 0.0;
	sixband_spin_matrices.S_z[2][2].real = 1.0/6;
	sixband_spin_matrices.S_z[2][3].real = 0.0;
	sixband_spin_matrices.S_z[2][4].real = -sqrt2/3.0;
	sixband_spin_matrices.S_z[2][5].real = 0.0;
	
	sixband_spin_matrices.S_z[3][0].real = 0.0;
	sixband_spin_matrices.S_z[3][1].real = 0.0;
	sixband_spin_matrices.S_z[3][2].real = 0.0;
	sixband_spin_matrices.S_z[3][3].real = -1.0/2;
	sixband_spin_matrices.S_z[3][4].real = 0.0;
	sixband_spin_matrices.S_z[3][5].real = 0.0;
	
	sixband_spin_matrices.S_z[4][0].real = 0.0;
	sixband_spin_matrices.S_z[4][1].real = 0.0;
	sixband_spin_matrices.S_z[4][2].real = -sqrt2/3.0;
	sixband_spin_matrices.S_z[4][3].real = 0.0;
	sixband_spin_matrices.S_z[4][4].real = -1.0/6;
	sixband_spin_matrices.S_z[4][5].real = 0.0;
	
	sixband_spin_matrices.S_z[5][0].real = 0.0;
	sixband_spin_matrices.S_z[5][1].real = -sqrt2/3.0;
	sixband_spin_matrices.S_z[5][2].real = 0.0;
	sixband_spin_matrices.S_z[5][3].real = 0.0;
	sixband_spin_matrices.S_z[5][4].real = 0.0;
	sixband_spin_matrices.S_z[5][5].real = 1.0/6;
	
	for (int m = 0; m < SIXBAND_HAM_DIM; m++) {
		for (int n = 0; n < SIXBAND_HAM_DIM; n++) {
			sixband_spin_matrices.S_x[m][n].imag = 0.0;
			sixband_spin_matrices.S_y[m][n].real = 0.0;
			sixband_spin_matrices.S_z[m][n].imag = 0.0;
			sixband_spin_matrices.S_plus[m][n].real = (sixband_spin_matrices.S_x[m][n].real - sixband_spin_matrices.S_y[m][n].imag)/sqrt(2);
			sixband_spin_matrices.S_plus[m][n].imag = (sixband_spin_matrices.S_x[m][n].imag + sixband_spin_matrices.S_y[m][n].real)/sqrt(2);
			sixband_spin_matrices.S_minus[m][n].real = (sixband_spin_matrices.S_x[m][n].real + sixband_spin_matrices.S_y[m][n].imag)/sqrt(2);
			sixband_spin_matrices.S_minus[m][n].imag = (sixband_spin_matrices.S_x[m][n].imag - sixband_spin_matrices.S_y[m][n].real)/sqrt(2);
		}
	}

}

static void calculate_Hkp(double kx, double ky, double kz, doublecomplex Hkp[SIXBAND_HAM_DIM][SIXBAND_HAM_DIM])
{
	const double sqrt2 = 1.4142135623730950488;
	const double sqrt3b2 = 1.2247448713915890491;

	double Hhh = (gamma1 + gamma2)*(kx*kx + ky*ky) + (gamma1 - 2*gamma2)*kz*kz;
	double Hlh = (gamma1 - gamma2)*(kx*kx + ky*ky) + (gamma1 + 2*gamma2)*kz*kz;
	double Hso = gamma1*(kx*kx + ky*ky + kz*kz) + 2*m0*Dso/(hbar*hbar);
	double bre = 2*sqrt(3)*gamma3*kx*kz;
	double bim = -2*sqrt(3)*gamma3*ky*kz;
	double cre = sqrt(3)*gamma2*(kx*kx - ky*ky);
	double cim = -2*sqrt(3)*gamma3*kx*ky;
	double d = -sqrt2*gamma2*(2*kz*kz - (kx*kx + ky*ky));


	/* Hkp is a 6x6 hermitian matrix */
	// addressing: H[mathematical column][mathematical row]

	Hkp[0][0].real = Hhh;
	Hkp[0][0].imag = 0;
	Hkp[0][1].real = -cre;
	Hkp[0][1].imag = cim;
	Hkp[0][2].real = -bre;
	Hkp[0][2].imag = bim;
	Hkp[0][3].real = 0;
	Hkp[0][3].imag = 0;
	Hkp[0][4].real = bre/sqrt2;
	Hkp[0][4].imag = -bim/sqrt2;
	Hkp[0][5].real = cre*sqrt2;
	Hkp[0][5].imag = -cim*sqrt2;
	
	Hkp[1][1].real = Hlh;
	Hkp[1][1].imag = 0;
	Hkp[1][2].real = 0;
	Hkp[1][2].imag = 0;
	Hkp[1][3].real = bre;
	Hkp[1][3].imag = -bim;
	Hkp[1][4].real = -bre*sqrt3b2;
	Hkp[1][4].imag = -bim*sqrt3b2;
	Hkp[1][5].real = -d;
	Hkp[1][5].imag = 0;
	
	Hkp[2][2].real = Hlh;
	Hkp[2][2].imag = 0;
	Hkp[2][3].real = -cre;
	Hkp[2][3].imag = cim;
	Hkp[2][4].real = d;
	Hkp[2][4].imag = 0;
	Hkp[2][5].real = -bre*sqrt3b2;
	Hkp[2][5].imag = bim*sqrt3b2;
	
	Hkp[3][3].real = Hhh;
	Hkp[3][3].imag = 0;
	Hkp[3][4].real = -cre*sqrt2;
	Hkp[3][4].imag = -cim*sqrt2;
	Hkp[3][5].real = bre/sqrt2;
	Hkp[3][5].imag = bim/sqrt2;
	
	Hkp[4][4].real = Hso;
	Hkp[4][4].imag = 0;
	Hkp[4][5].real = 0;
	Hkp[4][5].imag = 0;
	
	Hkp[5][5].real = Hso;
	Hkp[5][5].imag = 0;
	
	for (int m = SIXBAND_HAM_DIM - 1; m >= 0; --m) {
		for (int n = m; n < SIXBAND_HAM_DIM; ++n) {
			// sign chosen so that eigenvalues are positive
			Hkp[m][n].real *= hbar*hbar/(2*m0);
			Hkp[m][n].imag *= hbar*hbar/(2*m0);
#ifdef DEBUG
			assert(isfinite(Hkp[m][n].real));
			assert(isfinite(Hkp[m][n].imag));
#endif
		}
	}

	for (int i = SIXBAND_HAM_DIM - 1; i >= 0; --i) {
		for (int j = i - 1; j >= 0; --j) {
			Hkp[i][j].real = Hkp[j][i].real;
			Hkp[i][j].imag = - Hkp[j][i].imag;
		}
	}
}


static void calculate_Hkp_derivative(ham_struct* this, double c1, double c2, double c3, double d1, double d2, double d3, doublecomplex dHkp[SIXBAND_HAM_DIM][SIXBAND_HAM_DIM])
{
	const double sqrt2 = 1.4142135623730950488;
	const double sqrt3b2 = 1.2247448713915890491;

	// k.p models have cubic cells
	double kx = c1*this->bvect_lengths[0];
	double ky = c2*this->bvect_lengths[1];
	double kz = c3*this->bvect_lengths[2];
	double cx = d1*this->bvect_lengths[0];
	double cy = d2*this->bvect_lengths[1];
	double cz = d3*this->bvect_lengths[2];
	
	double Hhh = (gamma1 + gamma2)*2*(cx*kx + cy*ky) + 2*(gamma1 - 2*gamma2)*cz*kz;
	double Hlh = (gamma1 - gamma2)*2*(cx*kx + cy*ky) + 2*(gamma1 + 2*gamma2)*cz*kz;
	double Hso = 2*gamma1*(cx*kx + cy*ky + cz*kz);
	double bre = 2*sqrt(3)*gamma3*(cx*kz + kx*cz);
	double bim = -2*sqrt(3)*gamma3*(cy*kz + ky*cz);
	double cre = 2*sqrt(3)*gamma2*(cx*kx - cy*ky);
	double cim = -2*sqrt(3)*gamma3*(cx*ky + kx*cy);
	double d = -2*sqrt2*gamma2*(2*cz*kz - (cx*kx + cy*ky));


	/* dHkp is a 6x6 hermitian matrix */
	// addressing: H[mathematical column][mathematical row]

	dHkp[0][0].real = Hhh;
	dHkp[0][0].imag = 0;
	dHkp[0][1].real = -cre;
	dHkp[0][1].imag = cim;
	dHkp[0][2].real = -bre;
	dHkp[0][2].imag = bim;
	dHkp[0][3].real = 0;
	dHkp[0][3].imag = 0;
	dHkp[0][4].real = bre/sqrt2;
	dHkp[0][4].imag = -bim/sqrt2;
	dHkp[0][5].real = cre*sqrt2;
	dHkp[0][5].imag = -cim*sqrt2;
	
	dHkp[1][0].real = -cre;
	dHkp[1][0].imag = -cim;
	dHkp[1][1].real = Hlh;
	dHkp[1][1].imag = 0;
	dHkp[1][2].real = 0;
	dHkp[1][2].imag = 0;
	dHkp[1][3].real = bre;
	dHkp[1][3].imag = -bim;
	dHkp[1][4].real = -bre*sqrt3b2;
	dHkp[1][4].imag = -bim*sqrt3b2;
	dHkp[1][5].real = -d;
	dHkp[1][5].imag = 0;
	
	dHkp[2][0].real = -bre;
	dHkp[2][0].imag = -bim;
	dHkp[2][1].real = 0;
	dHkp[2][1].imag = 0;
	dHkp[2][2].real = Hlh;
	dHkp[2][2].imag = 0;
	dHkp[2][3].real = -cre;
	dHkp[2][3].imag = cim;
	dHkp[2][4].real = d;
	dHkp[2][4].imag = 0;
	dHkp[2][5].real = -bre*sqrt3b2;
	dHkp[2][5].imag = bim*sqrt3b2;
	
	dHkp[3][0].real = 0;
	dHkp[3][0].imag = 0;
	dHkp[3][1].real = bre;
	dHkp[3][1].imag = bim;
	dHkp[3][2].real = -cre;
	dHkp[3][2].imag = -cim;
	dHkp[3][3].real = Hhh;
	dHkp[3][3].imag = 0;
	dHkp[3][4].real = -cre*sqrt2;
	dHkp[3][4].imag = -cim*sqrt2;
	dHkp[3][5].real = bre/sqrt2;
	dHkp[3][5].imag = bim/sqrt2;
	
	dHkp[4][0].real = bre/sqrt2;
	dHkp[4][0].imag = bim/sqrt2;
	dHkp[4][1].real = -bre*sqrt3b2;
	dHkp[4][1].imag = bim*sqrt3b2;
	dHkp[4][2].real = d;
	dHkp[4][2].imag = 0;
	dHkp[4][3].real = -cre*sqrt2;
	dHkp[4][3].imag = cim*sqrt2;
	dHkp[4][4].real = Hso;
	dHkp[4][4].imag = 0;
	dHkp[4][5].real = 0;
	dHkp[4][5].imag = 0;
	
	dHkp[5][0].real = cre*sqrt2;
	dHkp[5][0].imag = cim*sqrt2;
	dHkp[5][1].real = -d;
	dHkp[5][1].imag = 0;
	dHkp[5][2].real = -bre*sqrt3b2;
	dHkp[5][2].imag = -bim*sqrt3b2;
	dHkp[5][3].real = bre/sqrt2;
	dHkp[5][3].imag = -bim/sqrt2;
	dHkp[5][4].real = 0;
	dHkp[5][4].imag = 0;
	dHkp[5][5].real = Hso;
	dHkp[5][5].imag = 0;
	
	for (int m = 0; m < SIXBAND_HAM_DIM; m++) {
		for (int n = 0; n < SIXBAND_HAM_DIM; n++) {
			// sign chosen so that eigenvalues are positive
			dHkp[m][n].real *= hbar*hbar/(2*m0);
			dHkp[m][n].imag *= hbar*hbar/(2*m0);
		}
	}

#if DEBUG == 1
    for (int i = 0; i < SIXBAND_HAM_DIM; i++)
    {
        for (int j = 0; j <= i; j++)
        {
            double dr = dHkp[i][j].real - dHkp[j][i].real;
            double di = dHkp[i][j].imag + dHkp[j][i].imag;
            if (fabs(dr) > 1E-16 || fabs(di) > 1E-16) {
                printf("Hamiltonian derivative non-hermicity detected for H[%d][%d]\n", i, j);
            }
        }
    }
#endif

}

static void generate_HBS(double exx, doublecomplex HBS[SIXBAND_HAM_DIM][SIXBAND_HAM_DIM])
{
	const double sqrt2 = 1.4142135623730950488;
	/* exx (=eyy) from Ohno exp
	   ezz = -2*exx*c12/c11, c12/c11 = 0.453 (Dietl et al. PRB63)
	   Q = ezz - (exx + eyy)/2
	   R = sqrt3*(exx - eyy)/2
	   */
	const double ezz = calculate_ezz(exx);
	const double Q = ezz - exx;// exx = eyy
	const double R = 0; // exx = eyy
	const double b = -1.7; //eV (Dietl et al. PRB63)

	/*  biaxial strain matrix */
	// addressing: HBS[mathematical column][mathematical row]

       HBS[0][0].real = -Q;
       HBS[0][0].imag = 0;
       HBS[0][1].real = -R;
       HBS[0][1].imag = 0;
       HBS[0][2].real = 0;
       HBS[0][2].imag = 0;
       HBS[0][3].real = 0;
       HBS[0][3].imag = 0;
       HBS[0][4].real = 0;
       HBS[0][4].imag = 0;
       HBS[0][5].real = sqrt2*R;
       HBS[0][5].imag = 0;
       
       HBS[1][0].real = -R;
       HBS[1][0].imag = 0;
       HBS[1][1].real = Q;
       HBS[1][1].imag = 0;
       HBS[1][2].real = 0;
       HBS[1][2].imag = 0;
       HBS[1][3].real = 0;
       HBS[1][3].imag = 0;
       HBS[1][4].real = 0;
       HBS[1][4].imag = 0;
       HBS[1][5].real = sqrt2*Q;
       HBS[1][5].imag = 0;
       
       HBS[2][0].real = 0;
       HBS[2][0].imag = 0;
       HBS[2][1].real = 0;
       HBS[2][1].imag = 0;
       HBS[2][2].real = Q;
       HBS[2][2].imag = 0;
       HBS[2][3].real = -R;
       HBS[2][3].imag = 0;
       HBS[2][4].real = -sqrt2*Q;
       HBS[2][4].imag = 0;
       HBS[2][5].real = 0;
       HBS[2][5].imag = 0;
       
       HBS[3][0].real = 0;
       HBS[3][0].imag = 0;
       HBS[3][1].real = 0;
       HBS[3][1].imag = 0;
       HBS[3][2].real = -R;
       HBS[3][2].imag = 0;
       HBS[3][3].real = -Q;
       HBS[3][3].imag = 0;
       HBS[3][4].real = -sqrt2*R;
       HBS[3][4].imag = 0;
       HBS[3][5].real = 0;
       HBS[3][5].imag = 0;
       
       HBS[4][0].real = 0;
       HBS[4][0].imag = 0;
       HBS[4][1].real = 0;
       HBS[4][1].imag = 0;
       HBS[4][2].real = -sqrt2*Q;
       HBS[4][2].imag = 0;
       HBS[4][3].real = -sqrt2*R;
       HBS[4][3].imag = 0;
       HBS[4][4].real = 0;
       HBS[4][4].imag = 0;
       HBS[4][5].real = 0;
       HBS[4][5].imag = 0;
       
       HBS[5][0].real = sqrt2*R;
       HBS[5][0].imag = 0;
       HBS[5][1].real = sqrt2*Q;
       HBS[5][1].imag = 0;
       HBS[5][2].real = 0;
       HBS[5][2].imag = 0;
       HBS[5][3].real = 0;
       HBS[5][3].imag = 0;
       HBS[5][4].real = 0;
       HBS[5][4].imag = 0;
       HBS[5][5].real = 0;
       HBS[5][5].imag = 0;
      
	for (int m = 0; m < 6; m++) {
		for (int n = 0; n < 6; n++) {
			HBS[m][n].real *= -b;
			HBS[m][n].imag *= -b;
		}
	}
}

/* Calculate 6-band Hamiltonian */
void calculate_sixband_hamiltonian(ham_struct_6band* this, double c1, double c2, double c3, doublecomplex H[SIXBAND_HAM_DIM][SIXBAND_HAM_DIM])
{
	const double D = this->super.D;
	const double* magdir = this->super.magdir;
	const double norm = sqrt(magdir[0]*magdir[0]+magdir[1]*magdir[1]+magdir[2]*magdir[2]);

	const ham_struct* hs = &this->super;
	// k.p models have cubic cells
	const double kx = c1*hs->bvect_lengths[0];
	const double ky = c2*hs->bvect_lengths[1];
	const double kz = c3*hs->bvect_lengths[2];
	
	calculate_Hkp(kx, ky, kz, H);

	for (int m = 0; m < SIXBAND_HAM_DIM; ++m) {
		for (int n = 0; n < SIXBAND_HAM_DIM; ++n) {
#ifdef DEBUG
			assert(isfinite(H[m][n].real));
			assert(isfinite(H[m][n].imag));
#endif
			// for possible other basis use, both real and imag parts of S added
			H[m][n].real = H[m][n].real + D/norm*(magdir[0]*sixband_spin_matrices.S_x[m][n].real + magdir[1]*sixband_spin_matrices.S_y[m][n].real + magdir[2]*sixband_spin_matrices.S_z[m][n].real) + this->HBS[m][n].real;
			H[m][n].imag = H[m][n].imag + D/norm*(magdir[0]*sixband_spin_matrices.S_x[m][n].imag + magdir[1]*sixband_spin_matrices.S_y[m][n].imag + magdir[2]*sixband_spin_matrices.S_z[m][n].imag) + this->HBS[m][n].imag;
#ifdef DEBUG
			assert(isfinite(H[m][n].real));
			assert(isfinite(H[m][n].imag));
#endif
		}
	}

#if DEBUG == 1
    for (int i = 0; i < SIXBAND_HAM_DIM; i++)
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

/* Create the structure which generates 6-band hamiltonians */
ham_struct* create_sixband_hamiltonian_structure(double MnX, double exx)
{
	initialize_sixband_spin_matrices();
	ham_struct_6band * hs6b = alloc_memory(sizeof(ham_struct_6band));
	ham_struct* hs = (ham_struct *) hs6b;
	setup_ham_struct(hs, MnX, SIXBAND_HAM_DIM, 3);
	hs->gen_ham = calculate_sixband_hamiltonian;
	assert(hs->gen_ham == hs6b->super.gen_ham);
	hs->gen_ham_derivative = calculate_Hkp_derivative;
	assert(hs->gen_ham_derivative == hs6b->super.gen_ham_derivative);
	assert(hs->work_matrix == hs6b->super.work_matrix);
	/* Six-band hamiltonian is p-band only */
	hs->p_bands_start = 0;
	hs->p_bands_cnt = SIXBAND_HAM_DIM;
	size_t size = sizeof(doublecomplex) * SIXBAND_HAM_DIM * SIXBAND_HAM_DIM;
	hs->S_z = copy_to_dynamic_array(sixband_spin_matrices.S_z, size);
	hs->S_x = copy_to_dynamic_array(sixband_spin_matrices.S_x, size);
	hs->S_y = copy_to_dynamic_array(sixband_spin_matrices.S_y, size);
	hs->S_minus = copy_to_dynamic_array(sixband_spin_matrices.S_minus, size);
	hs->S_plus = copy_to_dynamic_array(sixband_spin_matrices.S_plus, size);
	hs->find_cstep = find_cstep_kp;
	hs->inverted_fermi = 0;
	hs->bril_zone_integrand = bril_zone_integrand_kp;
	hs->exx = exx;
	setup_kp_geometry(hs, DEFAULT_KP_L);

	generate_HBS(exx, hs6b->HBS);
	
	return hs;
}

