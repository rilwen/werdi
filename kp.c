#include <math.h>
#include "hamiltonian.h"
#include "constants.h"

/* k.p in 3D, in our units */
double find_cstep_kp(ham_struct* this, double pmax, int npts)
{
        return 2*M_PI*pow(pmax, 1.0/3)/(0.5*(npts))/pow(this->bvect_lengths[0]*this->bvect_lengths[1]*this->bvect_lengths[2], 1.0/3);
}

/* Integrand over Brillouin zone in k.p model, for 3D zone */
double bril_zone_integrand_kp(ham_struct* this, double (*int_fun) (double, double, double), const double* x)
{
        double tmp = x[0] * x[1] * x[2];

        // integrated function vanishes in infinity
        if (0 == tmp)
        {
                return 0;
        }

        tmp *= tmp;

        double c1 = (1 - x[0]) / x[0];
        double c2 = (1 - x[1]) / x[1];
        double c3 = (1 - x[2]) / x[2];

        if (fabs(c1) > momentum_cutoff)
        {
                return 0;
        }
        if (fabs(c2) > momentum_cutoff)
        {
                return 0;
        }
        if (fabs(c3) > momentum_cutoff)
        {
                return 0;
        }

        double val = 0;
        val += int_fun(c1, c2, c3);
        /* Visit other octants of the coordinate system as well as +,+,+ */
        val += int_fun(c1, c2, -c3);
        val += int_fun(c1, -c2, -c3);
        val += int_fun(c1, -c2, c3);
        val += int_fun(-c1, -c2, c3);
        val += int_fun(-c1, -c2, -c3);
        val += int_fun(-c1, c2, -c3);
        val += int_fun(-c1, c2, c3);

        return val / tmp;
}

void setup_kp_geometry(ham_struct* hs, double L)
{
	hs->L = L;
	hs->bvect_lengths[0] = 1;
	hs->bvect_lengths[1] = 1;
	hs->bvect_lengths[2] = 1;
	hs->volume_factor = 1;
	hs->volume_factor_12 = 1;
	hs->cosine_product = 1;
	hs->cosine_product_12 = 1;
	hs->kx_vect[0] = 1;
	hs->kx_vect[1] = 0;
	hs->kx_vect[2] = 0;
	hs->ky_vect[0] = 0;
	hs->ky_vect[1] = 1;
	hs->ky_vect[2] = 0;
	hs->kz_vect[0] = 0;
	hs->kz_vect[1] = 0;
	hs->kz_vect[2] = 1;
}

