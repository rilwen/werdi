#include <math.h>
#include "fermi.h"

/* Fermi-Dirac distribution */
double fermi_dirac(double E, double EF, double kBTemp)
{
    if (kBTemp == 0)
    {
        return (E <= EF) ? 1 : 0;
    }
    else
    {
        return 1 / (exp((E - EF) / kBTemp) + 1);
    }
}

/* Fermi-Dirac distribution for bands occupied from top to bottom. */
double inverted_fermi_dirac(double E, double EF, double kBTemp)
{
    return fermi_dirac(EF, E, kBTemp);
}


/* Derivative of the Fermi-Dirac distribution over energy. */
double fermi_dirac_derivative(double E, double EF, double kBTemp)
{
    if (kBTemp == 0)
    {
        /* Simulate the Dirac delta. */
        if (E == EF)
        {
            return -1E20;
        }
        else
        {
            return 0;
        }
    }
    else
    {
        double tmp = exp((E - EF) / kBTemp);
        if (tmp == 0 || isinf(tmp))
        {
            return 0;
        }
        else
        {
            return - tmp / (kBTemp * (1 + tmp) * (1 + tmp));
        }
    }
}

/* Derivative of the Fermi-Dirac distribution over energy for bands occupied from top to bottom. */
double inverted_fermi_dirac_derivative(double E, double EF, double kBTemp)
{
    return fermi_dirac_derivative(EF, E, kBTemp);
}

