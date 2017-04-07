#ifndef MY_BIDIAG_H_INCLUDED
#define MY_BIDIAG_H_INCLUDED

#include <cmath>
#include <stdio.h>
#include "matrix.h"





void bidiag(matrix& X, matrix* out, bool& isUpper, bool checkErrors = false);
void check_bidiag(matrix& a, matrix* out, bool isUpper);



#endif // MY_BIDIAG_H_INCLUDED
