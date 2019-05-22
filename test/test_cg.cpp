/*
 *  This is a test example for MicroPP: a finite element library
 *  to solve microstructural problems for composite materials.
 *
 *  Copyright (C) - 2018 - Jimmy Aguilar Mena <kratsbinovish@gmail.com>
 *                         Guido Giuntoli <gagiuntoli@gmail.com>
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
#include <iomanip>

#include <ctime>
#include <cassert>

#include "micro.hpp"

using namespace std;

const double strain[6] = { 1., 2., 3., 1., 1., 1. };

int main (int argc, char *argv[])
{


	class test_t : public micropp<3> {

		public:
			test_t(const int size[3], const int micro_type, const double micro_params[5],
			       const struct material_base mat_params[2], const int solver)
				:micropp<3> (1, size, micro_type, micro_params, mat_params, NO_COUPLING, solver)
			{cout << solver <<endl;};

			~test_t() {};

			void assembly_and_solve(void)
			{

				const int ns[3] = { nx, ny, nz };
				const int nfield = dim;

				double *b = (double *) calloc(nndim, sizeof(double));
				double *du = (double *) calloc(nndim, sizeof(double));
				double *u = (double *) calloc(nndim, sizeof(double));

				double lerr, cg_err;
				memset(u, 0.0, nndim * sizeof(double));

				set_displ_bc(strain, u);
				lerr = assembly_rhs(u, nullptr, b);

				assembly_mat(&A, u, nullptr);
				if (solver == CGPILU)
					ell_ilu_factorization(&A);

				int cg_its = ell_solve(&A, b, du, &cg_err);

				for (int i = 0; i < nndim; ++i)
					u[i] += du[i];

				lerr = assembly_rhs(u, nullptr, b);
				cout
					<< "|RES| : " << lerr
					<< " CG_ITS : " << cg_its
					<< " CG_TOL : " << cg_err << endl;

				free(b);
				free(u);
				free(du);
			};

	};

	if (argc < 2) {
		/* argv[1] (n) : Problem size
		 * argv[2] (a) : Factor that says Ef = Em x a
		 */
		cerr << "Usage: " << argv[0] << " [n = 10] [solver = 0|1] [a = 1]" << endl;
		exit(1);
	}

	const int n = (argc > 1) ? atoi(argv[1]) : 10;
	const int solver = (argc > 2) ? atoi(argv[2]) : 0;
	const double a = (argc > 3) ? atoi(argv[3]) : 1.0;

	int size[3] = { n, n, n };

	int micro_type = 2;
	double micro_params[4] = { 1., 1., 1., 0.2 };

	material_base mat_params[2];
	material_set(&mat_params[0], 0, 1.0e8, 0.25, 5.0e4, 2.0e4, 0.0);
	material_set(&mat_params[1], 0, a * 1.0e8, 0.25, 5.0e4, 1.0e3, 0.0);

	test_t test(size, micro_type, micro_params, mat_params, solver);
	test.assembly_and_solve();

	return 0;
}
