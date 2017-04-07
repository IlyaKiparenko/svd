#include "my_bidiag_svd.h"

#include "matrix.h"
#include "my_tridiag_eig.h"

void bidiag_svd(matrix& A, matrix* out, bool isUpper, double eps, bool checkErrors) {
  int w = A.w;
  matrix tridiag_1 = compact_tridiag(A);
  matrix out1[2];
  tridiag_eig(tridiag_1, out1, eps, checkErrors);

  matrix tridiag_2 = compact_tridiag_t(A);
  matrix out2[2];
  tridiag_eig(tridiag_2, out2, eps, checkErrors);

  for (int i = 0; i < w; i++)
    out1[0][i] = sqrt(out1[0][i]);

  matrix U = out1[1];
  matrix V = out2[1];

  // correct signs of singular vectors
  double* pA = A.data + A.stride + 1;
  double* pU = U.data, *pV = V.data;

  if (isUpper) {
    bool needChange = false;
    int shift = (w - 1) * A.stride;
    pU += shift;
    pV += shift;
    for (int i = 0; i < w; i++) {
      if (fabs(pU[0]) >= 1. / w) {
        needChange = ((A[w - 1]*pV[0])*(pU[0]*out1[0][i]) < 0.); // if signs are different
      } else {
        int k = 0;
        do {
          pU -= U.stride;
          k++;
        } while (fabs(pU[0]) < 1. / w);
        double t;
        pV -= (k)*V.stride;
        t = pV[V.stride]*A(w - k, 1) + pV[0]*A[w - k - 1];
        needChange = (t * (pU[0]*out1[0][i]) < 0.); // if signs are different
        pU = U.data + i + shift;
        pV = V.data + i + shift;

      }
      if (needChange)
        for (int j = 0; j > -w*w; j-=V.stride)
          pV[j] *= -1.;
      pU ++;
      pV ++;
    }
  } else {
    bool needChange = false;
    for (int i = 0; i < w; i++) {
      if (fabs(pV[0]) >= 1. / w) {
        needChange = ((A[0]*pU[0])*(pV[0]*out1[0][i]) < 0.); // if signs are different
      } else {
        int k = 0;
        do {
          pV += V.stride;
          k++;
        } while (fabs(pV[0]) < 1. / w);
        double t;
        pU += (k - 1)*U.stride;
        t = pU[0]*A(k, 1) + pU[U.stride]*A[k];
        needChange = (t * (pV[0]*out1[0][i]) < 0.); // if signs are different
        pU = U.data + i;
        pV = V.data + i;

      }
      if (needChange)
        for (int j = 0; j < w*w; j+=V.stride)
          pV[j] *= -1.;
      pU ++;
      pV ++;
    }
  }

  if (isUpper) {
    out[0] = U;
    out[2] = V;
  } else {
    out[0] = V;
    out[2] = U;
  }
  out[1] = out1[0];

  if (checkErrors)
    check_bidiag_svd(A, out, isUpper);
}

void check_bidiag_svd(matrix X, matrix* out, bool isUpper) {
  matrix A = cr_bidiag(X, isUpper);

  matrix G = cr_diag(out[1]);

  matrix A_ = dot(out[0], dot(G, out[2].T()));

  double error = 0.;
  for (int i = 0; i < A.w*A.h; i++)
    error += pow(A[i] - A_[i], 2);
  error = sqrt(error);
  printf("Bidiag Svd error = %e\n", error);
}
