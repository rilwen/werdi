#include <math.h>
#include "secant.h"
#include "testing.h"

double square_test(double x, const void* a)
{
	const double* a_ptr = (const double*) a;
	return x * x - (*a_ptr);
}

double exp_test(double x, const void* a)
{
	const double* a_ptr = (const double*) a;
	return exp(x) - (*a_ptr);
}

int main(void)
{
	double a = 2;
	secant_results res = root_finder_secant(square_test, &a, 1.0, 1E-6, 1E-6, 20);
	test_double(sqrt(2), res.x, 1E-6, "SQUARE_ROOT");

	res = root_finder_secant(exp_test, &a, 1.0, 1E-6, 1E-6, 20);
	test_double(log(2), res.x, 1E-6, "LOG");
}
