#include <../inc/micro.hpp>
#include <../inc/csr.hpp>

int main (int argc, char *argv[])
{
  Problem problem (argc, argv);
  csr_matrix A;
  csr_vector x, dx, res;

  double eps [] = { 0.005, 0.0, 0.0 };
  int n = 10;

  // assembly res
  csr_alloc_A (A, n);
  csr_alloc_v (res, n);
  csr_alloc_v (dx, n);
  csr_alloc_v (x, n);

  // assembly A
  csr_assembly_A (A);

  // solve system  A * dx = -res

  // free memory
  csr_free_A (A);
  csr_free_v (res);
  csr_free_v (dx);
  csr_free_v (x);

  return 0;
}
