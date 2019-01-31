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
			       const material_t mat_params[2])
				:micropp<3> (1, size, micro_type, micro_params, mat_params, NO_COUPLING)
			{};

			~test_t() {};

			void assembly_and_solve(void)
			{

				double lerr, cg_err;
				memset(u[0], 0.0, nndim * sizeof(double));

				set_displ_bc(strain, u[0]);
				lerr = assembly_rhs(u[0], nullptr, b[0]);

				assembly_mat(&A[0], u[0], nullptr);
				int cg_its = ell_solve_cgpd(&A[0], b[0], du[0], &cg_err);

				for (int i = 0; i < nndim; ++i)
					u[0][i] += du[0][i];

				lerr = assembly_rhs(u[0], nullptr, b[0]);
				cout
					<< "|RES| : " << lerr
					<< " CG_ITS : " << cg_its
					<< " CG_TOL : " << cg_err << endl;

			};

	};

	if (argc < 2) {
		/* argv[1] (n) : Problem size
		 * argv[2] (a) : Factor that says Ef = Em x a
		 */
		cerr << "Usage: " << argv[0] << " [n = 10] [a = 1]" << endl;
	}
	omp_set_dynamic(0);     // Explicitly disable dynamic teams
	omp_set_num_threads(1); // Use 4 threads for all consecutive parallel regions

	const int n = (argc > 1) ? atoi(argv[1]) : 10;
	const double a = (argc > 2) ? atoi(argv[2]) : 1.0;

	int size[3] = { n, n, n };

	int micro_type = 2;
	double micro_params[4] = { 1., 1., 1., 0.2 };

	double Em = 1.0e8;

	material_t mat_params[2];
	mat_params[0].set(Em, 0.25, 1.0e8, 1.0e4, 0);
	mat_params[1].set(Em * a, 0.25, 1.0e8, 1.0e4, 0);

	test_t test(size, micro_type, micro_params, mat_params);
	test.assembly_and_solve();

	return 0;
}
