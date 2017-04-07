#ifndef MY_TRIDIAG_EIG_H
#define MY_TRIDIAG_EIG_H

#include "matrix.h"

void tridiag_eig(matrix& X, matrix* out, double eps = 10e-16, bool checkErrors = false);
void check_tridiag_eig(matrix X, matrix* out);

#endif // MY_TRIDIAG_EIG_H
