#include "micro.h"

void Problem::solve (void)
{

  ell_solve_cg (&solver_ell, &A_ell, b, du);

}
