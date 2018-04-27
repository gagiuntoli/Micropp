#include <iostream>
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

    ell_solve_cg (&solver, &A, b, du);

    // u = u + du
    for (int i=0; i<nn*dim; i++) {
      u[i] = u[i] + du[i];
    }

    its++;

  } while (its<NewRap_Its && tol>NewRap_Tol);

}
