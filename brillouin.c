#include <math.h>

double brillouin(double x, double S)
{
	 return ((2*S+1)/tanh((2*S+1)*x/2/S) - 1/tanh(x/2/S))/2/S;
}
