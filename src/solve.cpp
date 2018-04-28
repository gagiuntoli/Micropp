#include <iostream>
#include <iomanip> // print with format
#include "micro.h"

using namespace std;

void Problem::solve (void)
{

  ell_solve_cg (&solver, &A, b, du);

}

void Problem::newtonRaphson (void)
{
  int its = 0;
  double tol;

  do {

    tol = Assembly_b();
    cout << "NewRap It =" << its << " Tol = " << tol << endl;
    if (tol < NewRap_Tol) break;
    Assembly_A();

    for (int i=0; i<nn*dim; i++)
      du[i] = 0.0;
    ell_solve_cg (&solver, &A, b, du);
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
