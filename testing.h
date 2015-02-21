#ifndef _TESTING_H
#define _TESTING_H

#include "doublecomplex.h"

void test_double(double expected, double actual, double tolerance, const char* testname);
void test_int(int expected, int actual, const char* testname);
void test_doublecomplex(doublecomplex expected, doublecomplex actual, double tolerance, const char* testname);
void test_success(int result, const char* testname);
int is_hermitian(doublecomplex * matrix, int dim);

#endif /* _TESTING_H */
