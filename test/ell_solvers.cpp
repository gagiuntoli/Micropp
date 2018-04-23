/* 
   Program to test ell.cpp library
*/

#include <iostream>
#include <iomanip> // print with format
#include <fstream>
#include "ell.h"

#define N 100

using namespace std;

int main(void)
{
  cout << "Testing ell.cpp ..." << endl;

  double *x = (double*)malloc(N*sizeof(double));
  double *b = (double*)malloc(N*sizeof(double));
  double m_e[4] = {1, -1, -1, 1};

  ell_matrix m;
  ell_solver solver;
  ell_init(&m, N, N, 4);
  int ix[2], iy[2];
  for (int i = 0 ; i < N-1 ; i++) {
    ix[0] = i; ix[1] = i+1;
    iy[0] = i; iy[1] = i+1;
    ell_add_vals(&m, ix, 2, iy, 2, m_e);
  }
  ell_set_zero_row(&m, 0  , 1);
  ell_set_zero_row(&m, N-1, 1);
  ell_set_zero_col(&m, 0  , 1);
  ell_set_zero_col(&m, N-1, 1);

//  cout << "matrix =" << endl;
//  ell_print_full(&m);

  for (int i = 0 ; i < N ; i++) {
    x[i] = 0.0;
    b[i] = 1.0;
  }
  b[0]   = 0.0;
  b[N-1] = 0.0;
  solver.max_its = 100;
  solver.min_tol = 1.0e-5;
  ell_solve_cg(&solver, &m, b, x);

  ofstream fl;
  fl.open ("x_b_cg.dat");
  fl << fixed;
  fl << setprecision(2);
  for (int i = 0 ; i < N ; i++)
    fl << setprecision(3) << x[i] << " " << b[i] << endl;
  fl.close();
  cout << "err = " << setw(10) << solver.err << " its = " << solver.its << endl;

  for (int i = 0 ; i < N ; i++)
    x[i] = 0.0;
  solver.max_its = 2000;
  solver.min_tol = 1.0e-5;
  ell_solve_jacobi(&solver, &m, b, x);

  fl.open ("x_b_jac.dat");
  fl << fixed;
  fl << setprecision(2);
  for (int i = 0 ; i < N ; i++)
    fl << setprecision(3) << x[i] << " " << b[i] << endl;
  fl.close();
  cout << "err = " << setw(10) << solver.err << " its = " << solver.its << endl;

  return 0;
}
