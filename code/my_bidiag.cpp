#include <iostream>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <ctime>
using namespace std;

#include "my_bidiag.h"

// source
// http://www.sciencedirect.com/science/article/pii/S0024379504004276

// generate Householder reflection vector
void gen_hh(int shift, matrix x, double& beta, double& phi, matrix& v) {
  double alpha = 0.;
  for (int i = shift; i < x.w; i++) {
    alpha += x[i] * x[i];
    v[i] = x[i];
  }
  alpha = sqrt(alpha);
  if (alpha < 0.00001)
    beta = 0.;
  else
    beta = 1. / (alpha * (alpha + abs(v[shift])));
  if (v[shift] < 0.)
    alpha = -alpha;
  v[shift] += alpha;
  phi = -alpha;
}

// apply Householder reflection matrix to piece of A from left
void app_hh(int shift, int col, matrix& A, double beta, matrix& v) {
  double t;
  double* pA = A.data + col * A.stride;
  for (int j = col; j < A.h; j++) {
    t = 0.;
    for (int i = shift; i < A.w; i++)
      t += pA[i] * v[i];
    for (int i = shift; i < A.w; i++)
      pA[i] -= beta * t * v[i];
    pA += A.stride;
  }
}

// apply Householder reflection matrix to piece of A^T from left
void t_app_hh(int shift, int col, matrix& A, double beta, matrix& v) {
  double t;
  double* pA;

  for (int j = col; j < A.w; j++) {
    t = 0.;
    pA = A.data + j + shift * A.stride;
    for (int i = shift; i < A.h; i++) {
      t += pA[0] * v[i];
      pA += A.stride;
    }
    pA = A.data + j + shift * A.stride;
    for (int i = shift; i < A.h; i++) {
      pA[0] -= beta * t * v[i];
      pA += A.stride;
    }
  }
}

// calculate A^T * x (where x is vector)
void t_dot(matrix& A, matrix& x, matrix& y) {
  double t;
  double* p = A.data;
  int k;

  for (int j = 0; j < A.w; j++) {
    t = 0.;
    k = 0;
    for (int i = 0; i < A.h * A.w; i += A.stride) {
      t += p[i] * x[k];
      k++;
    }
    y[j] = t;
    p++;
  }
}

// calculate A * x (where x is vector)
void v_dot(matrix& A, matrix& x, matrix& y) {
  double t;
  double* pA = A.data;
  double* py = y.data + y.w - A.w;

  for (int j = 0; j < A.h; j++) {
    t = 0.;
    for (int i = 0; i < A.w; i++)
      t += pA[i] * x[i];
    py[j] = t;
    pA += A.stride;
  }
}

void bidiag(matrix& X, matrix* out, bool& isUpper, bool checkErrors) {
  isUpper = X.w < X.h;
  matrix A;
  if (isUpper)
    A = X.T();
  else
    A = X.copy();

  int m = A.w;
  int n = A.h;

  int l = m > n ? n : m;
  matrix U(m, l, true);
  matrix V(l, n, true);
  matrix B(l, 2, true); //a_i, b_i

  matrix v(m, 1);
  matrix z(m, 1);
  for (int i = 0; i < m; i++)
    z[i] = 0;
  double t = 0.;
  double beta;

  //Set V to Idenity
  for (int i = 0; i < V.w * V.h; i++)
    V[i] = 0.;
  for (int i = 0; i < V.w * V.h; i += V.stride + 1)
    V[i] = 1.;

  for (int i = 0; i < B.w * B.h; i++)
    B[i] = 0.;
  double* pA = A.data;
  double* pU = U.data;
  int k;
  for (k = 0; k < l - 2; k++) {
    if (fabs(B(k, 1)) < 0.000001)
      for (int i = 0; i < m; i++)
        pU[i] = pA[i];
    else
      for (int i = 0; i < m; i++)
        pU[i] = pA[i] - B(k, 1) * (pU - U.stride)[i];
    t = U.slice(k).norm();

    for (int i = 0; i < m; i++)
      pU[i] /= t;
    B(k, 0) = t;

    v_dot(A.slice(k + 1, n), U.slice(k), z);
    t = z.norm(k);

    if (!(t < 0.000001)) {
      gen_hh(k + 1, z, beta, B(k + 1, 1), v);
      t_app_hh(k + 1, 0 , A, beta, v);
      app_hh(k + 1, 0, V, beta, v);
    } else {
      B(k + 1, 1) = 0.;
    }

    pA += A.stride;
    pU += U.stride;
  }

  pU -= U.stride;

  for(int i = 0; i < 2; i++) {
    t = 0.;
    for (int i = 0; i < m; i++)
      t += pU[i] * pA[i];

    B(k, 1) = t;
    for (int i = 0; i < m; i++)
      (pU + U.stride)[i] = pA[i] - t * pU[i];

    t = U.slice(k).norm();
    B(k, 0) = t;
    pA += A.stride;
    pU += U.stride;
    for (int i = 0; i < m; i++)
      pU[i] /= t;
    k++;
  }

  out[1] = B;
  if (isUpper) {
    out[0] = U.T();
    out[2] = V;
  } else {
    out[0] = V;
    out[2] = U.T();
  }

  if (checkErrors)
    check_bidiag(X, out, isUpper);
}

void check_bidiag(matrix& a, matrix* out, bool isUpper) {
  matrix u = out[0];
  matrix v = out[2].T();

  matrix b = cr_bidiag(out[1], isUpper);

  matrix a_ = dot(u, dot(b, v));

  double error = 0.;

  for (int i = 0; i < a.w * a.h; i++)
    error += pow(a[i] - a_[i], 2);
  error = sqrt(error);
  printf("Bidiag Error = %e\n", error);
}




