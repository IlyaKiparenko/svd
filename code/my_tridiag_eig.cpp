#include "my_tridiag_eig.h"

#include "algorithm"

// source
// http://people.inf.ethz.ch/arbenz/ewp/Lnotes/chapter4.pdf (page 85-87)

// generate Givens rotation
void gen_r(double* data, int shift, double& c, double& s) {
  double n = sqrt(data[0]*data[0] + data[shift]*data[shift]);
  c = data[0] / n;
  s = data[shift] / n;
}

// generate Wilkinson's shift
double gen_shift(double* data, int shift) {
  double s = (data[0] - data[1+shift])/2;
  double p = fabs(s) + sqrt(s*s + data[1]*data[1]);
  return data[1+shift] - ((s > 0) ? 1. : -1)*data[1]*data[1] / p;
}

// apply Givens rotation matrix to A from left
void app_r(matrix A, int row, int col, int len, double c, double s) {
  double t;
  double* p1 = A.data + row*A.stride + col, *p2 = p1 + A.stride;
  for (int i = 0; i < len; i++) {
    t = p1[0];
    p1[0] = c*p1[0] +  s*p2[0];
    p2[0] = -s*t + c*p2[0];
    p1++;
    p2++;
  }
}

// apply Givens rotation matrix to A^T from left
void app_r_t(matrix A, int row, int col, int len, double c, double s) {
  double t;
  double* p1 = A.data + row*A.stride + col, *p2 = p1 + 1;
  for (int i = 0; i < len; i++) {
    t = p1[0];
    p1[0] = c*p1[0] + s*p2[0];
    p2[0] = -s*t + c*p2[0];
    p1 += A.stride;
    p2 += A.stride;
  }
}

// step of qr-algorithm with shift
void tridiag_qr_step(matrix A, double shift, int n, matrix koef) {
  double* koef2 = koef.data + koef.stride;
  double* p1 = A.data;
  double t;
  for (int i = 0, j = 0; j < n; i+= A.stride + 1, j++)
    A[i] -= shift;

  // A = QR
  for (int i = 0; i < n - 1; i++) {
    gen_r(p1, A.stride, koef[i], koef2[i]);
    app_r(A, i, i, 2, koef[i], koef2[i]);
    if (i < n - 2) {
      t = p1[2 + A.stride];
      p1[2] = t*koef2[i];
      p1 += A.stride;
      p1[2] = t*koef[i];
      p1++;
    }
  }
  for (int i = 0; i < n - 1; i++) {
    app_r_t(A, i, i, 2, koef[i], koef2[i]);
  }
  for (int i = 1, j = 1; j < n; i+=A.stride + 1, j++)
    A[i] = A[i + A.stride - 1];
  for (int i = 0, j = 0; j < n; i += A.stride + 1, j++)
    A[i] += shift;
}

void tridiag_eig(matrix& X, matrix* out, double eps, bool checkErrors) {
  matrix A = cr_tridiag(X);

  int w = A.w;
  matrix U(w, w, true);
  matrix S(w, 1, true);
  matrix b(w);
  for (int i = 0; i < w*w; i++)
    U[i] = 0.;
  for (int i = 0; i < w*w; i+= U.stride + 1)
    U[i] = 1.;

  int i = 0;
  double* p = A.data + (w-2)*A.stride + w - 2;
  matrix koef(w, 2);
  while (i < w - 1) {
    double shift = gen_shift(p, A.stride);
    tridiag_qr_step(A, shift, w - i, koef);
    for (int j = 0; j < w - i - 1; j++) {
      app_r(U, j, 0, w, koef[j], koef(j, 1));
    }
    if (fabs(p[1]) < eps) {
      i++;
      p -= A.stride + 1;
    }
  }

  for (int i = 0, j = 0; i < w; i++, j += A.stride + 1)
    S[i] = A[j];

  // sort in ascending order
  //int index[w];
  int* index = new int[w];
  for (int i = 0; i < w; i++)
    index[i] = i;
  std::sort(index, index + w, [&S](int i1, int i2) {return S[i1] < S[i2];});

  i = 0;
  int k = -1;
  double t = 0.;
  bool start = false;
  while (i < w) {
    if (i != index[i]) {

      if (!start) {
        t = S[i];
        b = U.slice(i);
        start = true;
        k = i;
      } else if (k == index[i]) {
        S[i] = t;
        U.slice(i) = b;
        start = false;
        index[i] = i;
        i = k + 1;
        continue;
      }

      S[i] = S[index[i]];
      U.slice(i) = U.slice(index[i]);
      int p = index[i];
      index[i] = i;
      i = p;
    } else
      i++;
  }

  delete [] index;
  out[0] = S;
  out[1] = U.T();

  if (checkErrors)
    check_tridiag_eig(X, out);
}

void check_tridiag_eig(matrix X, matrix* out) {
  matrix A = cr_tridiag(X);

  matrix G = cr_diag(out[0]);

  matrix A_ = dot(out[1], dot(G, out[1].T()));

  double error = 0.;
  for (int i = 0; i < A.w*A.h; i++)
    error += pow(A[i] - A_[i], 2);
  error = sqrt(error);
  printf("Tridiag Eig error = %e\n", error);
}

