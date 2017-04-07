#include <iostream>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>

#include "my_timer.h"
#include "matrix.h"
#include "my_svd.h"
#include "alglib_svd.h"

int main(int argc, char** argv) {
  init_timer;
  int m, n;

  if (argc > 3) {
    m = atoi(argv[1]);
    n = atoi(argv[2]);
    int seed = atoi(argv[3]);
    srand(seed);
    printf("%i %i %i\n", m, n, seed);
  } else {
    m = 4;
    n = 4;
    srand(22);
  }

  matrix a(m, n);

  for (int j = 0; j < n; j++)
    for (int i = 0; i < m; i++)
      a(i, j) = rand() % 1000000 / 1000000. ;

  matrix out[3];
  timeit(svd(a, out, 10e-16, false));
  check_svd(a, out, "My svd");

  matrix out1[3];
  timeit(alglib_svd(a, out1));
  check_svd(a, out1, "My svd on alglib");

  matrix out2[3];
  timeit(alglib_real_svd(a, out2));
  check_svd(a, out2, "Real svd on alglib");
  return 0;
}
