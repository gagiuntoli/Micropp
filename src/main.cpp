#include <iostream>
#include <ctime>
#include <cmath>
#include "micro.h"
#include "csr.h"
#include "ell.h"

using namespace std;

int main (int argc, char *argv[])
{
  try {

    csr_matrix A;
    csr_vector x, dx, res;
    Problem problem (argc, argv);
    clock_t start, end;
    double t_solve, t_assembly;

    double eps [] = { 0.005, 0.0, 0.0 };
    int nn = problem.nn;
    int dim = problem.dim;

    csr_alloc_A (A, nn*dim);
    csr_alloc_v (res, nn*dim);
    csr_alloc_v (dx, nn*dim);
    csr_alloc_v (x, nn*dim);
    csr_set_A (A, problem);

    cout << "Init micropp -> preparing ..." << endl;
    cout << "dim   = " << problem.dim << endl;
    cout << "nx    = " << problem.nx << endl;
    cout << "ny    = " << problem.ny << endl;
    cout << "nn    = " << problem.nn << endl;
    cout << "nelem = " << problem.nelem << endl;
    cout << "dx    = " << problem.dx << endl;
    cout << "dy    = " << problem.dy << endl;

    // assembly
    start = clock();
    csr_assembly_A   (A, problem);
    csr_assembly_res (res, problem);
    csr_assembly_res (dx, problem);
    end = clock();
    t_assembly = double(end - start) / CLOCKS_PER_SEC;

    // solve system  A * dx = -res
    start = clock();
    csr_cg (A, dx, res);
    end = clock();
    t_solve = double(end - start) / CLOCKS_PER_SEC;
    cout << "time assembly = "<<t_assembly<<endl;
    cout << "time solve    = "<<t_solve   <<endl;

    // free memory
    csr_free_A (A);
    csr_free_v (res);
    csr_free_v (dx);
    csr_free_v (x);

  } catch (int &e) {
    cerr << "Error : " << e << endl;
    return 1;
  }

  return 0;
}
