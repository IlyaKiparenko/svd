#include "my_svd.h"

#include "matrix.h"
#include "my_bidiag.h"
#include "my_bidiag_svd.h"

void svd(matrix A, matrix* out, double eps, bool checkErrors) {
  matrix out1[3];
  bool isUpper;
  bidiag(A, out1, isUpper, checkErrors);
  matrix out2[3];
  bidiag_svd(out1[1], out2, isUpper, eps, checkErrors);
  out[0] = dot(out1[0], out2[0]);
  out[2] = dot(out1[2], out2[2]);
  out[1] = out2[1];

  if (checkErrors)
    check_svd(A, out);
}

void check_svd(matrix A, matrix* out, std::string name) {
  matrix S = cr_diag(out[1]);
  matrix A_ = dot(out[0], dot(S, out[2].T()));
  double error = 0.;
  for (int i = 0; i < A.w * A.h; i++)
    error += pow(A[i] - A_[i], 2);
  error = sqrt(error);
  printf("%s \n\tError = %e\n", name.c_str(), error);
}
