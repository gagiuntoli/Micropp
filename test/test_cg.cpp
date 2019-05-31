/*
 *  This is a test example for MicroPP: a finite element library
 *  to solve microstructural problems for composite materials.
 *
 *  Copyright (C) - 2018 - Jimmy Aguilar Mena <kratsbinovish@gmail.com>
 *                         Guido Giuntoli <gagiuntoli@gmail.com>
 *                         Judicaël Grasset <judicael.grasset@stfc.ac.uk>
 *                         Alejandro Figueroa <afiguer7@maisonlive.gmu.edu>
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
			test_t(const micropp_params_t params)
				:micropp<3> (params){};

			~test_t() {};

			void assembly_and_solve(void)
			{

				const int ns[3] = { nx, ny, nz };
				const int nfield = dim;

				ell_matrix A;  // Jacobian
				ell_init(&A, nfield, dim, ns, CG_ABS_TOL, CG_REL_TOL, CG_MAX_ITS);
				double *b = (double *) calloc(nndim, sizeof(double));
				double *du = (double *) calloc(nndim, sizeof(double));
				double *u = (double *) calloc(nndim, sizeof(double));

				double lerr, cg_err;
				memset(u, 0.0, nndim * sizeof(double));

				set_displ_bc(strain, u);
				lerr = assembly_rhs(u, nullptr, b);

				assembly_mat(&A, u, nullptr);
				int cg_its = ell_solve_cgpd(&A, b, du, &cg_err);

				for (int i = 0; i < nndim; ++i)
					u[i] += du[i];

				lerr = assembly_rhs(u, nullptr, b);
				cout
					<< "|RES| : " << lerr
					<< " CG_ITS : " << cg_its
					<< " CG_TOL : " << cg_err << endl;

				ell_free(&A);
				free(b);
				free(u);
				free(du);
			};

	};

	if (argc < 2) {
		/* argv[1] (n) : Problem size
		 * argv[2] (a) : Factor that says Ef = Em x a
		 */
		cerr << "Usage: " << argv[0] << " [n = 10] [a = 1]" << endl;
	}

	const int n = (argc > 1) ? atoi(argv[1]) : 10;
	const double a = (argc > 2) ? atoi(argv[2]) : 1.0;

	micropp_params_t mic_params;

	mic_params.ngp = 1;
	mic_params.size[0] = n;
	mic_params.size[1] = n;
	mic_params.size[2] = n;
	mic_params.type = MIC_SPHERE;
	material_set(&mic_params.materials[0], 0, 1.0e7, 0.3, 0.0, 0.0, 0.0);
	material_set(&mic_params.materials[1], 0, a * 1.0e7, 0.3, 0.0, 0.0, 0.0);
	material_set(&mic_params.materials[2], 0, 1.0e7, 0.3, 0.0, 0.0, 0.0);

	mic_params.print();

	test_t test(mic_params);
	test.assembly_and_solve();

	return 0;
}
