#ifndef _SECANT_H
#define _SECANT_H

#define SECANT_OK 0
#define SECANT_EVALS 1
#define SECANT_FUN 2

typedef struct {
	double x;
	double fx;
	int evalcnt;
	int flag;
} secant_results;

secant_results root_finder_secant(double (*fun)(double, void*), void* params,
	double x0, double tolx, double tolfun, int maxevals);

#endif
