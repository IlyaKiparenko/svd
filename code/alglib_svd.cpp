#include "alglib_svd.h"

#include "alglib/LinAlg.h"
#include "alglib/ap.h"

void matrix_mult(alglib::real_2d_array a,
                 alglib::real_2d_array b,
                 alglib::real_2d_array& c,
                 bool FirstTransposed = false,
                 bool SecondTransposed = false) {
  int M = a.rows();
  int K = a.cols();
  int N = b.cols();
  c.setlength(M, N);
  alglib::rmatrixgemm(M, N, K,
                     1.0, a, 0, 0, FirstTransposed,
                     b, 0, 0, SecondTransposed,
                     0., c, 0, 0);
}

matrix algmatrix_to_matrix(alglib::real_2d_array a, bool transpose = false) {
  matrix *b;
  if (!transpose) {
    b = new matrix(a.rows(), a.cols(), true);
    for (int i = 0; i < b->w; i++)
      for (int j = 0; j < b->h; j++)
        b->data[i + j*b->w] = a[i][j];
  } else {
    b = new matrix(a.cols(), a.rows(), true);
    for (int i = 0; i < b->w; i++)
      for (int j = 0; j < b->h; j++)
        b->data[i + j*b->w] = a[j][i];
  }
  b->temp = true;
  return *b;
}

void alglib_svd(matrix X, matrix* out) {
  alglib::real_2d_array A;
  A.setcontent(X.w, X.h, X.data); // copy of X

  alglib::real_1d_array TauQ, TauP, TridiagonalMain, TridiagonalSecond;
  alglib::real_1d_array MainDiagonal, SecondDiagonal;


  alglib::real_2d_array Q, PT, U, VT, LeftSingular, RightSingular;

  // bidiagonal decomposition
  alglib::rmatrixbd(A, A.rows(), A.cols(), TauQ, TauP);

  int l = std::min(A.rows(), A.cols());
  alglib::rmatrixbdunpackpt(A, A.rows(), A.cols(), TauP, l, PT);
  alglib::rmatrixbdunpackq(A, A.rows(), A.cols(), TauQ, l, Q);

  bool isUpper;
  alglib::rmatrixbdunpackdiagonals(A, A.rows(), A.cols(),
                                   isUpper,
                                   MainDiagonal, SecondDiagonal);

  int w = MainDiagonal.length();

  TridiagonalMain.setlength(w);
  TridiagonalSecond.setlength(w);
  matrix Singular(w);

  // Tr = B*B^T
  int i;
  for (i = 0; i < w - 1; i++) {
    TridiagonalMain[i] = pow(MainDiagonal[i], 2) + pow(SecondDiagonal[i], 2);
    TridiagonalSecond[i] = MainDiagonal[i + 1]*SecondDiagonal[i];
  }
  TridiagonalMain[i] = pow(MainDiagonal[i], 2);

  //eig decomposition
  bool result = alglib::smatrixtdevd(TridiagonalMain, TridiagonalSecond,
                                     w, 2, LeftSingular);

  // Tr = B^T*B
  TridiagonalMain[0] = pow(MainDiagonal[0], 2);
  for (int i = 1; i < w; i++) {
    TridiagonalMain[i] = pow(MainDiagonal[i], 2) + pow(SecondDiagonal[i - 1], 2);
    TridiagonalSecond[i - 1] = MainDiagonal[i - 1]*SecondDiagonal[i - 1];
  }

  //eig decomposition
  result = alglib::smatrixtdevd(TridiagonalMain, TridiagonalSecond,
                                w, 2, RightSingular);

  // correct eigvectors sign, because they are singular vectors
  for (int i = 0; i < RightSingular.cols(); i++) {
    double t = MainDiagonal[0]*RightSingular[0][i];
    if (isUpper)
      t += SecondDiagonal[0]*RightSingular[1][i];
    if (t*LeftSingular[0][i] < 0.)
      for (int j = 0; j < RightSingular.rows(); j++)
        RightSingular[j][i] *= -1.;
  }

  for (int i = 0; i < w; i++)
    Singular[i] = pow(TridiagonalMain[i], 0.5);

  if (isUpper) {
    matrix_mult(Q, LeftSingular, U);
    matrix_mult(RightSingular, PT, VT, 1);
  } else {
    matrix_mult(Q, RightSingular, U);
    matrix_mult(LeftSingular, PT, VT, 1);
  }

  out[0] = algmatrix_to_matrix(U, true);
  out[1] = Singular;
  out[2] = algmatrix_to_matrix(VT);
}

void alglib_real_svd(matrix X, matrix* out) {
  alglib::real_2d_array U, VT, A;
  alglib::real_1d_array S;
  A.setcontent(X.w, X.h, X.data);
  alglib::rmatrixsvd(A, A.rows(), A.cols(), 1, 1, 2, S, U, VT);
  matrix s(S.length());
  for (int i = 0; i < s.w; i++)
    s[i] = S[i];
  s.temp = true;
  out[0] = algmatrix_to_matrix(U).T();
  out[1] = s;
  out[2] = algmatrix_to_matrix(VT);
}



