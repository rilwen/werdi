#ifndef _TIGHTBINDING_H
#define _TIGHTBINDING_H
#include "hamiltonian.h"

typedef enum { GEOM_XY, GEOM_Z } tb_geometry; // type of the geometry
typedef enum { JANCU, DICARLO } tb_parameterization; // type of the parameterization
typedef enum { EXACT, FIVEPOINT } tb_derivative; // type of the derivative
typedef enum { BULK, PLANE } tb_material; // type of the material

// Mn concentration for which we have input files.
extern const double TB_MN_X;

/* 
 * Construct a hamiltonian structure generating bulk TB hamiltonian.
 * derivative - method to calculate derivatives (EXACT or FIVEPOINT).
 * geometry - lattice geometry (GEOM_XY or GEOM_Z)
 * parameterization - JANCU or DICARLO
 * exx - xx strain
 */
ham_struct* create_tightbinding_bulk_hamiltonian_structure(tb_derivative derivative, tb_geometry geometry, tb_parameterization parameterization, double exx);

/* 
 * Construct a hamiltonian structure generating bulk TB hamiltonian.
 * derivative - method to calculate derivatives (FIVEPOINT).
 * geometry - lattice geometry (GEOM_XY or GEOM_Z)
 * parameterization - JANCU or DICARLO
 * exx - xx strain
 */
ham_struct* create_tightbinding_plane_hamiltonian_structure(tb_derivative derivative, tb_geometry geometry, tb_parameterization parameterization, int nlayers, double Ex, double Ey, double Ez, double (*concentration_gradient)(double*), double exx);

#endif
