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

void Problem::newtonRaphson (void)
{
  int its = 0;
  double tol;

  do {

    tol = Assembly_b();

    if (flag_print_newton == true)
      cout << "NewRap It =" << its << " Tol = " << tol << endl;

    if (tol < NewRap_Tol) break;

    Assembly_A();

    for (int i=0; i<nn*dim; i++)
      du[i] = 0.0;

    if (dim == 2) {
      ell_solve_cgpd_2D (&solver, &A, dim, nx, ny, b, du);
    } else if (dim == 3) {
      ell_solve_cgpd_struct (&solver, &A, dim, dim, nn, b, du);
    }

    if (flag_print_solver == true)
      cout << "CG Its = " << solver.its << " Err = " << solver.err << endl;

    // u = u + du
    for (int i=0; i<nn*dim; i++)
      u[i] = u[i] + du[i];

    if (flag_print_u == true)
      for (int i=0; i<nn; i++)
	cout << setw(5) << u[i*dim] << " " << u[i*dim+1] << endl;

    if (flag_print_du == true)
      for (int i=0; i<nn; i++)
	cout << setw(5) << setprecision(4) <<  du[i*dim] << " " << du[i*dim+1] << endl;

    its++;

  } while (its<NewRap_Its && tol>NewRap_Tol);

}
