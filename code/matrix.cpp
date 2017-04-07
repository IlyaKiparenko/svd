#include "matrix.h"

#define printm(a) print(a, #a)

void print(matrix a, std::string name = "Matrix") {
  printf("\t %s = \n", name.c_str());
  for (int j = 0; j < a.h; j++) {
    for (int i = 0; i < a.w; i++)
      printf("% .8f ", a(i, j));
    printf("\n");
  }
  printf("\n");
}

#define fprintm(a) fprint(a, #a)

void fprint(matrix& a, std::string name = "Matrix") {
  static bool first = true;
  FILE* f;
  if (first) {
    f = fopen("out.txt", "w");
    first = false;
  } else
    f = fopen("out.txt", "a");
  fprintf(f, "\t %s = \n", name.c_str());
  fprintf(f, "%i %i\n", a.w, a.h);
  for (int j = 0; j < a.h; j++) {
    for (int i = 0; i < a.w; i++)
      fprintf(f, "%.16f ", a(i, j));
    fprintf(f, "\n");
  }
  fprintf(f, "\n");

  fclose(f);
}

//create bidiagonal from compact form
matrix cr_bidiag(matrix src, bool isUpper) {
  int n = src.w;
  matrix* out = new matrix(n, n);
  double* pS = src.data;
  double* pO = out->data;
  for (int i = 0; i < n*n; i++)
    out->data[i] = 0.;
  for (int i = 0; i < n; i++) {
    pO[0] = pS[i];
    pO += out->stride + 1;
  }
  if (isUpper)
    pO = out->data + 1;
  else
    pO = out->data + out->stride;
  pS += src.stride;
  for (int i = 1; i < n; i++) {
    pO[0] = pS[i];
    pO += out->stride + 1;
  }
  out->temp = true;
  return *out;
}

//create tridiagonal from compact form
matrix cr_tridiag(matrix src) {
  int n = src.w;
  matrix* out = new matrix(n, n);
  double* pS = src.data;
  double* pO = out->data;
  for (int i = 0; i < n*n; i++)
    out->data[i] = 0.;
  for (int i = 0; i < n; i++) {
    pO[0] = pS[i];
    pO += out->stride + 1;
  }
  pO = out->data + 1;
  double* pO_2 = out->data + out->stride;
  pS += src.stride;
  for (int i = 1; i < n; i++) {
    pO_2[0] = pO[0] = pS[i];
    pO_2 += out->stride + 1;
    pO += out->stride + 1;
  }
  out->temp = true;
  return *out;
}

//create diagonal from compact form
matrix cr_diag(matrix s) {
  matrix* out = new matrix(s.w, s.w);
  double* p = out->data;
  int k = 0;
  for (int i = 0; i < s.w*s.w; i++)
    p[i] = 0.;
  for (int i = 0; i < s.w*s.w; i += out->stride + 1)
    p[i] = s[k++];
  out->temp = true;
  return *out;
}

//create tridiagonal compact form from bidiagonal compact form Tr = B*B^T
matrix compact_tridiag(matrix bidiag) {
  int w = bidiag.w;
  matrix* t  = new matrix(w, bidiag.h);

  double* second = bidiag.data + bidiag.stride + 1;
  double* second_t = t->data + t->stride + 1;
  for (int i = 0; i < w - 1; i++) {
    t->data[i] = pow(bidiag[i], 2) + pow(second[i], 2);
    second_t[i] = bidiag[i+1] * second[i];
  }
  t->data[w - 1] = pow(bidiag[w - 1], 2);
  t->temp = true;
  return *t;
}

//create tridiagonal compact form from bidiagonal compact form Tr = B^T*B
matrix compact_tridiag_t(matrix bidiag) {
  int w = bidiag.w;
  matrix* t  = new matrix(w, bidiag.h);

  double* second = bidiag.data + bidiag.stride + 1;
  double* second_t = t->data + t->stride + 1;
  for (int i = 1; i < w; i++) {
    t->data[i] = pow(bidiag[i], 2) + pow(second[i-1], 2);
    second_t[i-1] = bidiag[i-1] * second[i-1];
  }
  t->data[0] = pow(bidiag[0], 2);
  t->temp = true;
  return *t;
}

matrix dot(matrix a, matrix b) {
  matrix tB = b.T();
  matrix* c = new matrix(b.w, a.h);
  double t;
  double* t1 = a.data;
  double* t2 = tB.data;
  double* t3 = c->data;
  for (int i = 0; i < a.h; i++) {
    t2 = tB.data;
    for (int k = 0; k < b.w; k++) {
      t = 0.;
      for (int j = 0; j < a.w; j++)
        t += t1[j] * t2[j];
      t3[k] = t;
      t2 += tB.stride;
    }
    t1 += a.stride;
    t3 += c->stride;
  }
  c->temp = true;
  return *c;
}




