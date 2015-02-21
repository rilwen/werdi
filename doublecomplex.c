#include "doublecomplex.h"

const doublecomplex CMPLX_ZERO = {.real = 0, .imag = 0};

doublecomplex cmplx(double real, double imag)
{
	doublecomplex z;
	z.real = real;
	z.imag = imag;
	return z;
}
