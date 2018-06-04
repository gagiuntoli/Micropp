#include <iostream>
#include <iomanip> // print with format
#include "micro.h"

#define MAX_TOL 1.0e-5

using namespace std;

void Problem::solve (void)
{
  ell_solver solver;
  solver.max_its = 2000;
  solver.min_tol = 1.0e-8;
  if (dim == 2)
    ell_solve_cgpd_2D (&solver, &A, dim, nx, ny, b, du);
  else if (dim == 3)
    ell_solve_cgpd_struct (&solver, &A, dim, dim, nn, b, du);
  //cout << "CG Its = " << solver.its << " Err = " << solver.err << endl;
}

void Problem::newtonRaphson (bool *non_linear)
{
  int its = 0;
  double tol;

  do {

    tol = Assembly_b(non_linear);

    if (tol < MAX_TOL) break;
    //cout << "NEWTON-R ITS = " << its << " TOL = " << tol << endl;

    Assembly_A();

    for (int i=0; i<nn*dim; i++)
      du[i] = 0.0;

    solve();

    for (int i=0; i<nn*dim; i++)
      u[i] = u[i] + du[i];

    its++;

  } while ((its < 20) && (tol > MAX_TOL));
  cout << "NEWTON-R ITS = " << its << " TOL = " << tol << endl;

}
