#ifndef _DOUBLECOMPLEX_H
#define _DOUBLECOMPLEX_H

typedef struct
{
    double real;
    double imag;
}
doublecomplex;

extern const doublecomplex CMPLX_ZERO; // ZERO

/* Reaturn a doublecomplex number */
doublecomplex cmplx(double real, double imag);

#endif /* _DOUBLECOMPLEX */

