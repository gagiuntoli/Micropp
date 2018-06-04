#include <iostream>
#include <iomanip> // print with format
#include "micro.h"

#define CG_MAX_TOL 1.0e-8
#define CG_MAX_ITS 2000
#define NR_MAX_TOL 1.0e-5
#define NR_MAX_ITS 40

using namespace std;

void Problem::solve (void)
{
  ell_solver solver;
  solver.max_its = CG_MAX_ITS;
  solver.min_tol = CG_MAX_TOL;

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

    if (tol < NR_MAX_TOL) break;
    //cout << "NEWTON-R ITS = " << its << " TOL = " << tol << endl;

    Assembly_A();

    for (int i=0; i<nn*dim; i++)
      du[i] = 0.0;

    solve();

    for (int i=0; i<nn*dim; i++)
      u[i] = u[i] + du[i];

    its++;

  } while ((its < NR_MAX_ITS) && (tol > NR_MAX_TOL));
  cout << "NEWTON-R ITS = " << its << " TOL = " << tol << endl;

}
