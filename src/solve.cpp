#include <iostream>
#include <iomanip> // print with format
#include "micro.h"

using namespace std;

void Problem::solve (void)
{

  if (dim == 2) {
    ell_solve_cgpd_2D (&solver, &A, dim, nx, ny, b, du);
  } else if (dim == 3) {
    ell_solve_cgpd_struct (&solver, &A, dim, dim, nn, b, du);
  }
}

void Problem::newtonRaphson (double *vars_old, double *vars_new, bool *non_linear_flag)
{
  int its = 0;
  double tol;

  do {

    tol = Assembly_b(vars_old, vars_new, non_linear_flag);

    if (tol < NewRap_Tol) break;
    cout << "NewRap It =" << its << " Tol = " << tol << endl;

    Assembly_A(vars_old);

    for (int i=0; i<nn*dim; i++)
      du[i] = 0.0;

    if (dim == 2) {
      ell_solve_cgpd_2D (&solver, &A, dim, nx, ny, b, du);
    } else if (dim == 3) {
      ell_solve_cgpd_struct (&solver, &A, dim, dim, nn, b, du);
    }
    //cout << "CG Its = " << solver.its << " Err = " << solver.err << endl;

    for (int i=0; i<nn*dim; i++)
      u[i] = u[i] + du[i];

    its++;

  } while (its<NewRap_Its && tol>NewRap_Tol);
  cout << "NewRap It =" << its << " Tol = " << tol << endl;

}
