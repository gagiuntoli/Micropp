#include "micro.h"

void Problem::solve (void)
{

  ell_solve_cg (&solver, &A, b, du);

}
