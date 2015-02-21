#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "secant.h"

const double slope_dx = 0.001;

static int find_nonzero_slope(double (*fun)(double, void*), void* params,
	double x, double fx, double* inv_slope, int* evalcnt, int maxevals)
{
	double x2;
	double fx2;
	int cnt = *evalcnt;
	
	x2 = x;
	fx2 = fx;
	do {
		x2 += slope_dx;
		fx2 = fun(x2, params);
		*inv_slope = (x - x2)/(fx - fx2);
		cnt++;
	} while (!isfinite(*inv_slope) && cnt < maxevals && isfinite(fx2));
	*evalcnt = cnt;

	if (!isfinite(fx2)) {
		return SECANT_FUN; // solved function returned NaN or +- Inf
	} else if (cnt >= maxevals) {
		return SECANT_EVALS; // too many function evaluations
	} else {
		return SECANT_OK;
	}
}

secant_results root_finder_secant(double (*fun)(double, void*), void* params,
	double x0, double tolx, double tolfun, int maxevals)
{
	double x; // current value of solution
	double px; // previous value of solution
	double fx; // value of fun(x, params)
	double pfx; // value of fun(px, params)
	double inv_slope; // inverse slope
	secant_results res;
	res.evalcnt = 1;

	x = x0;
	fx = fun(x, params);
	if (fabs(fx) < tolfun)
	{
		res.flag = SECANT_OK;
		res.x = x;
		res.fx = fx;
		return res;
	}
	if (!isfinite(fx)) {
		res.flag = SECANT_FUN;
		return res;
	}
	res.flag = find_nonzero_slope(fun, params, x, fx, &inv_slope, &res.evalcnt, maxevals);
	if (res.flag != SECANT_OK) {
		return res;
	}
	do {
		px = x;
		pfx = fx;
		x -= inv_slope*fx;
		fx = fun(x, params);
		res.evalcnt++;
		if (!isfinite(fx)) {
			res.flag = SECANT_FUN;
			return res;
		}
		inv_slope = (x - px)/(fx - pfx);
		if (!isfinite(inv_slope)) {
			res.flag = find_nonzero_slope(fun, params, x, fx, &inv_slope, &res.evalcnt, maxevals);
			if (res.flag != SECANT_OK) {
				return res;
			}
		}
	} while ((fabs(x - px) > tolx || fabs(fx - pfx) > tolfun) && res.evalcnt < maxevals);
	res.x = x;
	res.fx = fx;
	if (res.evalcnt >= maxevals) {
		res.flag = SECANT_EVALS;
	} else {
		res.flag = SECANT_OK;
	}
	return res;
}


