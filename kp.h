#ifndef _KP_H
#define _KP_H

#include "hamiltonian.h"

/* default k.p length unit */
#define DEFAULT_KP_L 0.27604

/* k.p in 3D */
double find_cstep_kp(ham_struct* this, double pmax, int npts);

/* Integrand over Brillouin zone in k.p model, for 3D zone */
double bril_zone_integrand_kp(ham_struct* this, double (*int_fun) (double, double, double), const double* x);

/* Setup the k-space basis and related stuff which is convenient for k.p models */
void setup_kp_geometry(ham_struct* hs, double L);

#endif
