#ifndef _CONSTANTS_H
#define _CONSTANTS_H

/*
 * new and good
 *
 * [S] = 1 hbar = 1.05457168E-34 kg*m^2/s		(action)
 * [m] = 1 m0 = 9.1093826E-31 kg			(mass)
 * [E] = 1 eV = 1.60217646E-19 J			(energy)
 * [T] = 1 K						(temperature)
 * [v] = sqrt([E]/[m]) = 4.1938E+5 m/s			(velocity)
 * [L] = [S]/sqrt([E]*[m]) = 2.7604E-10 m		(length)
 * [t] = [L]/[v] = 6.5821E-16 s				(time)
 * [Q] = 1e = 1.60217653E-19 C				(electric charge)
 * [Ex] = 3.622610520703641e+009 V/m			(electric field)
 * small fiels: E * a * e < 0.01 * EF (<< EF) -- 3.539823008849558e+006 V/m
 * a - lattice const., e - electrc charge, E - el. field
 *
 */

/* If not indicated otherwise, physical constants given in above units. */

/*spin-orbit coupling [eV]*/
#define DeltaSO 0.341
/* electron rest mass */
#define m0 1.0
/* Planck constant / 2 pi */
#define hbar 1.0
/* Boltzmann constant in eV/K */
#define kB 8.6173E-05
#define kappa 1.2
#define g0 2.0023193043718
//const double hartree_eV = 27.2113845;
//const double Dso = 0.34; // spin-orbit coupling [eV]
//const double D = 0.15; // band splitting [eV]
/* light speed in vacuum */
#define c_vac 7148.6
/* light speed in SI units */
#define c_SI  2997992458
//const double b = -0.07349864906726815021; /* = -2/hartree_eV; */
//const double d = -0.17639675776144356051; /* -4.8/hartree_eV; */
#define gamma1 6.85
#define gamma2 2.1
#define gamma3 2.9

#define BetaSI (1.2 * 0.054) // beta parameter in eV*nm^3, 1.2 factor included
#define a0SI 0.5633 // lattice parameter with no strain in nm
#define Spin 2.5 // Mn total spin

/* algorith constants */

// momentum cutoff (mainly applies minus and plus bands in Abolfath's PRB01 paper)
#define momentum_cutoff 1E10

#define M_PI 3.14159265358979323846264338327950288419716939937510

#endif
