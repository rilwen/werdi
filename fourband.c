#include <math.h>
#include <stdio.h>
#include <assert.h>
#include "doublecomplex.h"
#include "fourband.h"
#include "constants.h"
#include "utils.h"
#include "hamiltonian.h"
#include "kp.h"

#define DEBUG 0

fourband_spin_matrices_struct fourband_spin_matrices;

const double Dso=0.341;

/* Initialize S_x, S_y and S_z */
void initialize_fourband_spin_matrices(void)
{
	const double sqrt3 = 1.73205080756887729353;
	
	// addressing: fourband_spin_matrices.S_x[mathematical column][mathematical row]
	fourband_spin_matrices.S_x[0][0].real = 0.0;
	fourband_spin_matrices.S_x[0][1].real = 0.0;
	fourband_spin_matrices.S_x[0][2].real = 1.0/(2*sqrt3);
	fourband_spin_matrices.S_x[0][3].real = 0.0;
	
	fourband_spin_matrices.S_x[1][0].real = 0.0;
	fourband_spin_matrices.S_x[1][1].real = 0.0;
	fourband_spin_matrices.S_x[1][2].real = 1.0/3;
	fourband_spin_matrices.S_x[1][3].real = 1.0/(2*sqrt3);
	
	fourband_spin_matrices.S_x[2][0].real = 1.0/(2*sqrt3);
	fourband_spin_matrices.S_x[2][1].real = 1.0/3;
	fourband_spin_matrices.S_x[2][2].real = 0.0;
	fourband_spin_matrices.S_x[2][3].real = 0.0;
	
	fourband_spin_matrices.S_x[3][0].real = 0.0;
	fourband_spin_matrices.S_x[3][1].real = 1.0/(2*sqrt3);
	fourband_spin_matrices.S_x[3][2].real = 0.0;
	fourband_spin_matrices.S_x[3][3].real = 0.0;
	
	fourband_spin_matrices.S_y[0][0].imag = 0.0;
	fourband_spin_matrices.S_y[0][1].imag = 0.0;
	fourband_spin_matrices.S_y[0][2].imag = 1.0/(2*sqrt3);
	fourband_spin_matrices.S_y[0][3].imag = 0.0;
	
	fourband_spin_matrices.S_y[1][0].imag = 0.0;
	fourband_spin_matrices.S_y[1][1].imag = 0.0;
	fourband_spin_matrices.S_y[1][2].imag = -1.0/3;
	fourband_spin_matrices.S_y[1][3].imag = 1.0/(2*sqrt3);
	
	fourband_spin_matrices.S_y[2][0].imag = -1.0/(2*sqrt3);
	fourband_spin_matrices.S_y[2][1].imag = 1.0/3;
	fourband_spin_matrices.S_y[2][2].imag = 0.0;
	fourband_spin_matrices.S_y[2][3].imag = 0.0;
	
	fourband_spin_matrices.S_y[3][0].imag = 0.0;
	fourband_spin_matrices.S_y[3][1].imag = -1.0/(2*sqrt3);
	fourband_spin_matrices.S_y[3][2].imag = 0.0;
	fourband_spin_matrices.S_y[3][3].imag = 0.0;
	
	
	fourband_spin_matrices.S_z[0][0].real = 1.0/2;
	fourband_spin_matrices.S_z[0][1].real = 0.0;
	fourband_spin_matrices.S_z[0][2].real = 0.0;
	fourband_spin_matrices.S_z[0][3].real = 0.0;
	
	fourband_spin_matrices.S_z[1][0].real = 0.0;
	fourband_spin_matrices.S_z[1][1].real = -1.0/6;
	fourband_spin_matrices.S_z[1][2].real = 0.0;
	fourband_spin_matrices.S_z[1][3].real = 0.0;
	
	fourband_spin_matrices.S_z[2][0].real = 0.0;
	fourband_spin_matrices.S_z[2][1].real = 0.0;
	fourband_spin_matrices.S_z[2][2].real = 1.0/6;
	fourband_spin_matrices.S_z[2][3].real = 0.0;
	
	fourband_spin_matrices.S_z[3][0].real = 0.0;
	fourband_spin_matrices.S_z[3][1].real = 0.0;
	fourband_spin_matrices.S_z[3][2].real = 0.0;
	fourband_spin_matrices.S_z[3][3].real = -1.0/2;
	

	for (int m = 0; m < FOURBAND_HAM_DIM; m++) {
		for (int n = 0; n < FOURBAND_HAM_DIM; n++) {
			fourband_spin_matrices.S_x[m][n].imag = 0.0;
			fourband_spin_matrices.S_y[m][n].real = 0.0;
			fourband_spin_matrices.S_z[m][n].imag = 0.0;
			fourband_spin_matrices.S_plus[m][n].real = (fourband_spin_matrices.S_x[m][n].real - fourband_spin_matrices.S_y[m][n].imag)/sqrt(2);
			fourband_spin_matrices.S_plus[m][n].imag = (fourband_spin_matrices.S_x[m][n].imag + fourband_spin_matrices.S_y[m][n].real)/sqrt(2);
			fourband_spin_matrices.S_minus[m][n].real = (fourband_spin_matrices.S_x[m][n].real + fourband_spin_matrices.S_y[m][n].imag)/sqrt(2);
			fourband_spin_matrices.S_minus[m][n].imag = (fourband_spin_matrices.S_x[m][n].imag - fourband_spin_matrices.S_y[m][n].real)/sqrt(2);
		}
	}

}

static void calculate_Hkp(double kx, double ky, double kz, doublecomplex Hkp[FOURBAND_HAM_DIM][FOURBAND_HAM_DIM])
{
	
	double Hhh = (gamma1 + gamma2)*(kx*kx + ky*ky) + (gamma1 - 2*gamma2)*kz*kz;
	double Hlh = (gamma1 - gamma2)*(kx*kx + ky*ky) + (gamma1 + 2*gamma2)*kz*kz;
	double bre = 2*sqrt(3)*gamma3*kx*kz;
	double bim = -2*sqrt(3)*gamma3*ky*kz;
	double cre = sqrt(3)*gamma2*(kx*kx - ky*ky);
	double cim = -2*sqrt(3)*gamma3*kx*ky;


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
	
	Hkp[1][0].real = -cre;
	Hkp[1][0].imag = -cim;
	Hkp[1][1].real = Hlh;
	Hkp[1][1].imag = 0;
	Hkp[1][2].real = 0;
	Hkp[1][2].imag = 0;
	Hkp[1][3].real = bre;
	Hkp[1][3].imag = -bim;
	
	Hkp[2][0].real = -bre;
	Hkp[2][0].imag = -bim;
	Hkp[2][1].real = 0;
	Hkp[2][1].imag = 0;
	Hkp[2][2].real = Hlh;
	Hkp[2][2].imag = 0;
	Hkp[2][3].real = -cre;
	Hkp[2][3].imag = cim;
	
	Hkp[3][0].real = 0;
	Hkp[3][0].imag = 0;
	Hkp[3][1].real = bre;
	Hkp[3][1].imag = bim;
	Hkp[3][2].real = -cre;
	Hkp[3][2].imag = -cim;
	Hkp[3][3].real = Hhh;
	Hkp[3][3].imag = 0;
	
	for (int m = 0; m < FOURBAND_HAM_DIM; m++) {
		for (int n = 0; n < FOURBAND_HAM_DIM; n++) {
			// sign chosen so that eigenvalues are positive
			Hkp[m][n].real *= hbar*hbar/(2*m0);
			Hkp[m][n].imag *= hbar*hbar/(2*m0);
#ifdef DEBUG
			assert(isfinite(Hkp[m][n].real));
			assert(isfinite(Hkp[m][n].imag));
#endif
		}
	}
}


static void calculate_Hkp_derivative(ham_struct* this, double kx, double ky, double kz, double cx, double cy, double cz, doublecomplex dHkp[FOURBAND_HAM_DIM][FOURBAND_HAM_DIM])
{
	
	double Hhh = (gamma1 + gamma2)*2*(cx*kx + cy*ky) + 2*(gamma1 - 2*gamma2)*cz*kz;
	double Hlh = (gamma1 - gamma2)*2*(cx*kx + cy*ky) + 2*(gamma1 + 2*gamma2)*cz*kz;
	double bre = 2*sqrt(3)*gamma3*(cx*kz + kx*cz);
	double bim = -2*sqrt(3)*gamma3*(cy*kz + ky*cz);
	double cre = 2*sqrt(3)*gamma2*(cx*kx - cy*ky);
	double cim = -2*sqrt(3)*gamma3*(cx*ky + kx*cy);


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
	
	dHkp[1][0].real = -cre;
	dHkp[1][0].imag = -cim;
	dHkp[1][1].real = Hlh;
	dHkp[1][1].imag = 0;
	dHkp[1][2].real = 0;
	dHkp[1][2].imag = 0;
	dHkp[1][3].real = bre;
	dHkp[1][3].imag = -bim;
	
	dHkp[2][0].real = -bre;
	dHkp[2][0].imag = -bim;
	dHkp[2][1].real = 0;
	dHkp[2][1].imag = 0;
	dHkp[2][2].real = Hlh;
	dHkp[2][2].imag = 0;
	dHkp[2][3].real = -cre;
	dHkp[2][3].imag = cim;
	
	dHkp[3][0].real = 0;
	dHkp[3][0].imag = 0;
	dHkp[3][1].real = bre;
	dHkp[3][1].imag = bim;
	dHkp[3][2].real = -cre;
	dHkp[3][2].imag = -cim;
	dHkp[3][3].real = Hhh;
	dHkp[3][3].imag = 0;

	for (int m = 0; m < FOURBAND_HAM_DIM; m++) {
		for (int n = 0; n < FOURBAND_HAM_DIM; n++) {
			// sign chosen so that eigenvalues are positive
			dHkp[m][n].real *= hbar*hbar/(2*m0);
			dHkp[m][n].imag *= hbar*hbar/(2*m0);
		}
	}

#if DEBUG == 1
    for (int i = 0; i < FOURBAND_HAM_DIM; i++)
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

static void generate_HBS(double exx, doublecomplex HBS[FOURBAND_HAM_DIM][FOURBAND_HAM_DIM])
{
	/* exx (=eyy) from Ohno exp
	   ezz = -2*exx*c12/c11, c12/c11 = 0.453 (Dietl et al. PRB63)
	   Q = ezz - (exx + eyy)/2
	   R = sqrt3*(exx - eyy)/2
	   */
	const double ezz = -2*exx*0.453;
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
       
       HBS[1][0].real = -R;
       HBS[1][0].imag = 0;
       HBS[1][1].real = Q;
       HBS[1][1].imag = 0;
       HBS[1][2].real = 0;
       HBS[1][2].imag = 0;
       HBS[1][3].real = 0;
       HBS[1][3].imag = 0;
       
       HBS[2][0].real = 0;
       HBS[2][0].imag = 0;
       HBS[2][1].real = 0;
       HBS[2][1].imag = 0;
       HBS[2][2].real = Q;
       HBS[2][2].imag = 0;
       HBS[2][3].real = -R;
       HBS[2][3].imag = 0;
       
       HBS[3][0].real = 0;
       HBS[3][0].imag = 0;
       HBS[3][1].real = 0;
       HBS[3][1].imag = 0;
       HBS[3][2].real = -R;
       HBS[3][2].imag = 0;
       HBS[3][3].real = -Q;
       HBS[3][3].imag = 0;
       
	for (int m = 0; m < 6; m++) {
		for (int n = 0; n < 6; n++) {
			HBS[m][n].real *= -b;
			HBS[m][n].imag *= -b;
		}
	}
}

/* Calculate 6-band Hamiltonian */
void calculate_fourband_hamiltonian(ham_struct_4band* this, double kx, double ky, double kz, doublecomplex H[FOURBAND_HAM_DIM][FOURBAND_HAM_DIM])
{
	double D = this->super.D;
	double* magdir = this->super.magdir;
	double norm = sqrt(magdir[0]*magdir[0]+magdir[1]*magdir[1]+magdir[2]*magdir[2]);
	
	calculate_Hkp(kx, ky, kz, H);

	for (int m = 0; m < FOURBAND_HAM_DIM; m++) {
		for (int n = 0; n < FOURBAND_HAM_DIM; n++) {
#ifdef DEBUG
			assert(isfinite(H[m][n].real));
			assert(isfinite(H[m][n].imag));
#endif
			// for possible other basis use, both real and imag parts of S added
			H[m][n].real = H[m][n].real + D/norm*(magdir[0]*fourband_spin_matrices.S_x[m][n].real + magdir[1]*fourband_spin_matrices.S_y[m][n].real + magdir[2]*fourband_spin_matrices.S_z[m][n].real) + this->HBS[m][n].real;
			H[m][n].imag = H[m][n].imag + D/norm*(magdir[0]*fourband_spin_matrices.S_x[m][n].imag + magdir[1]*fourband_spin_matrices.S_y[m][n].imag + magdir[2]*fourband_spin_matrices.S_z[m][n].imag) + this->HBS[m][n].imag;
#ifdef DEBUG
			assert(isfinite(H[m][n].real));
			assert(isfinite(H[m][n].imag));
#endif
		}
	}

#if DEBUG == 1
    for (int i = 0; i < FOURBAND_HAM_DIM; i++)
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
ham_struct* create_fourband_hamiltonian_structure(double MnX, double exx)
{
	initialize_fourband_spin_matrices();
	ham_struct_4band * hs4b = alloc_memory(sizeof(ham_struct_4band));
	ham_struct* hs = (ham_struct *) hs4b;
	setup_ham_struct(hs, MnX, FOURBAND_HAM_DIM, 3);
	hs->gen_ham = calculate_fourband_hamiltonian;
	assert(hs->gen_ham == hs4b->super.gen_ham);
	hs->gen_ham_derivative = calculate_Hkp_derivative;
	assert(hs->gen_ham_derivative == hs4b->super.gen_ham_derivative);
	assert(hs->work_matrix == hs4b->super.work_matrix);
	/* Six-band hamiltonian is p-band only */
	hs->p_bands_start = 0;
	hs->p_bands_cnt = FOURBAND_HAM_DIM;
	size_t size = sizeof(doublecomplex) * FOURBAND_HAM_DIM * FOURBAND_HAM_DIM;
	hs->S_z = copy_to_dynamic_array(fourband_spin_matrices.S_z, size);
	hs->S_x = copy_to_dynamic_array(fourband_spin_matrices.S_x, size);
	hs->S_y = copy_to_dynamic_array(fourband_spin_matrices.S_y, size);
	hs->S_minus = copy_to_dynamic_array(fourband_spin_matrices.S_minus, size);
	hs->S_plus = copy_to_dynamic_array(fourband_spin_matrices.S_plus, size);
	hs->find_cstep = find_cstep_kp;
	hs->inverted_fermi = 0;
	hs->bril_zone_integrand = bril_zone_integrand_kp;
	setup_kp_geometry(hs, DEFAULT_KP_L);

	hs->periodic_bc_dim = 3;

	hs4b->exx = exx;
	generate_HBS(exx, hs4b->HBS);
	
	return hs;
}

