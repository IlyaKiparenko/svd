#ifndef MY_BIDIAG_SVD_H
#define MY_BIDIAG_SVD_H

#include "matrix.h"

void bidiag_svd(matrix& A, matrix* out, bool isUpper, double eps = 10e-16, bool checkErrors = false);
void check_bidiag_svd(matrix A, matrix* out, bool isUpper);

#endif // MY_BIDIAG_SVD_H
