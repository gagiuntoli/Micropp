/*
 *  MicroPP : 
 *  Finite element library to solve microstructural problems for composite materials.
 *
 *  Copyright (C) - 2018 - Guido Giuntoli <gagiuntoli@gmail.com>
 *  
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <iostream>
#include <iomanip> // print with format
#include "micro.h"

#define CG_MAX_TOL 1.0e-8
#define CG_MAX_ITS 2000
#define NR_MAX_TOL 1.0e-5
#define NR_MAX_ITS 40

using namespace std;

void Problem::solve()
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

  NR_norm = tol;
  NR_its = its;
  NR_non_linear = (*non_linear == true) ? 1:0;
}
