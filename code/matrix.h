#ifndef MATRIX_H_INCLUDED
#define MATRIX_H_INCLUDED

#include <cmath>
#include <stdio.h>
#include <string>

class matrix;

#define printm(a) print(a, #a)

void print(matrix a, std::string name);

#define fprintm(a) fprint(a, #a)

void fprint(matrix& a, std::string name);

class matrix {
public:
  int w, h, stride;
  double* data;
  bool copyied;
  bool temp;
  matrix() {
    data = 0;
    copyied = false;
    temp = false;
  }
  matrix(int m, int n = 1, bool temp_ = false) {
    w = m;
    h = n;
    stride = m;
    data = new double[m * n];
    temp = temp_;
    copyied = false;
  }
  matrix(const matrix& b) {
    w = b.w;
    h = b.h;
    stride = b.stride;
    data = b.data;
    copyied = true;
    temp = false;
  }

  inline double& operator[](int i) {
    return data[i];
  }
  inline double& operator()(int i, int j) {
    return data[i + j * stride];
  }

  matrix& slice(int k1) {
    matrix* b = new matrix();
    b->w = w;
    b->h = 1;
    b->stride = w;
    b->data = data + k1 * stride;
    b->copyied = true;
    return *b;
  }

  matrix& slice(int k1, int k2) {
    matrix* b = new matrix();
    b->w = w;
    b->h = k2 - k1;
    b->stride = w;
    b->data = data + k1 * stride;
    b->copyied = true;
    return *b;
  }

  double norm(int shift = 0) {
    double t = 0.;
    for (int i = shift; i < w * h; i++)
      t += data[i] * data[i];
    return sqrtl(t);
  }

  void normalize() {
    double t = norm();
    for (int i = 0; i < w * h; i++)
      data[i] /= t;
  }

  void operator=(matrix b) {
    if (b.temp && data == 0) {
      // move
      w = b.w;
      h = b.h;
      copyied = false;
      b.copyied = true;
      stride = b.stride;
      data = b.data;
    } else {
      if (data == 0) {
        w = b.w;
        h = b.h;
        copyied = false;
        stride = b.stride;
        data = new double[w * h];
      }

      for (int i = 0; i < b.w * b.h; i++)
        data[i] = b.data[i];
    }
  }

  matrix& T() {
    matrix* b = new matrix(h, w);
    double* t = b->data;
    double* t2;
    for (int i = 0; i < w; i++) {
      t2 = data + i;
      for (int j = 0; j < h; j++) {
        t[j] = *t2;
        t2 += stride;
      }
      t += b->stride;
    }
    b->temp = true;
    return *b;
  }

  matrix copy() {
    matrix* b = new matrix(w, h);
    for (int i = 0; i < w*h; i++)
      b->data[i] = data[i];
    b->temp = true;
    return *b;
  }

  ~matrix() {
    if (!copyied)
      delete [] data;
  }
};

matrix cr_bidiag(matrix src, bool isUpper);
matrix cr_tridiag(matrix src);
matrix cr_diag(matrix s);

matrix compact_tridiag_t(matrix bidiag);
matrix compact_tridiag(matrix bidiag);

matrix dot(matrix a, matrix b);

#endif // MATRIX_H_INCLUDED
