#ifndef MY_SVD_H
#define MY_SVD_H

#include "matrix.h"
#include <string>

void svd(matrix A, matrix* out, double eps = 10e-16, bool checkErrors = false);
void check_svd(matrix A, matrix* out, std::string name = "SVD");

#endif // MY_SVD_H
