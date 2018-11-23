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

#define D_EPS 5.0e-5

int main (int argc, char *argv[])
{
	const int dim = 3;
	if (argc < 5) {
		cerr << "Usage: " << argv[0] << " nx ny nz dir [steps]" << endl;
		return(1);
	}

	const int nx = atoi(argv[1]);
	const int ny = atoi(argv[2]);
	const int nz = atoi(argv[3]);
	const int dir = atoi(argv[4]);
	const int time_steps = (argc > 5 ? atoi(argv[5]) : 10);  // Optional value
	int size[3] = { nx, ny, nz };
	char filename[128];

	ofstream file;
	file.open("result.dat");

	int micro_type = MIC_QUAD_FIB_XZ;
	double micro_params[4] = { 1., 1., 1., .15 };

	material_t mat_params[2];
	mat_params[0].set(3.0e7, 0.25, 1.0e4, 1.0e5, 1);
	mat_params[1].set(3.0e9, 0.25, 1.0e5, 1.0e5, 0);

	micropp<3> micro(1, size, micro_type, micro_params, mat_params);
	micro.print_info();

	double sig[6], ctan[36];
	double eps[6] = { 0. };
	int sigma_solver_its[NR_MAX_ITS];
	double sigma_solver_err[NR_MAX_ITS];
	double sigma_nr_err[NR_MAX_ITS];

	for (int t = 0; t < time_steps; ++t) {

		cout << "time step = " << t << endl;

		if (t < time_steps / 2)
			eps[dir] += D_EPS;
		else
			eps[dir] -= D_EPS;

		micro.set_macro_strain(0, eps);
		micro.homogenize();
		micro.get_macro_stress(0, sig);
		micro.get_macro_ctan(0, ctan);
		int newton_its = micro.get_sigma_newton_its(0);
		int non_linear = micro.is_non_linear(0);
		micro.get_sigma_solver_its(0, sigma_solver_its);
		micro.get_sigma_solver_err(0, sigma_solver_err);
		micro.get_sigma_newton_err(0, sigma_nr_err);

		snprintf(filename, 128, "micropp_%d", t);
		micro.output(0, filename);

		micro.update_vars();

		cout << "non_linear = \t" << non_linear << "\tnewton its =\t" << newton_its << endl;
		for (int i = 0; i < newton_its + 1; ++i)
			cout    << " nr : " << i
				<< " nr_err : " << scientific << sigma_nr_err[i]
				<< " cg_its : " << sigma_solver_its[i]
				<< " cg_err : " << scientific << sigma_solver_err[i] << endl;

		cout << "eps =\t";
		for (int i = 0; i < 6; ++i)
			cout << setw(14) << eps[i] << "\t";
		cout << endl;

		cout << "sig =\t";
		for (int i = 0; i < 6; ++i)
			cout << setw(14) << sig[i] << "\t";
		cout << endl;

		cout << "ctan =\n";
		for (int i = 0; i < 6; ++i) {
			for (int j = 0; j < 6; ++j)
				cout << setw(14) << ctan[i * 6 + j] << "\t";
			cout << endl;
		}

		cout << endl;
		file << setw(14) << eps[dir] << "\t" << sig[dir] << "\t" << endl;
	}

	file.close();
	return 0;
}
