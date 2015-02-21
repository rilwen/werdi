#ifndef _FERMI_H
#define _FERMI_H


/* Fermi-Dirac distribution */
double fermi_dirac(double E, double EF, double kBTemp);
/* Fermi-Dirac distribution for bands occupied from top to bottom. */
double inverted_fermi_dirac(double E, double EF, double kBTemp);
/* Derivative of the Fermi-Dirac distribution for bands occupied from top to bottom. */
double fermi_dirac_derivative(double E, double EF, double kBTemp);
/* Derivative of the Fermi-Dirac distribution over energy for bands occupied from top to bottom. */
double inverted_fermi_dirac_derivative(double E, double EF, double kBTemp);

#endif /* _FERMI_H */

