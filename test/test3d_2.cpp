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

int main (int argc, char *argv[])
{
	const int dim = 3;
	if (argc < 4) {
		cerr << "Usage: " << argv[0] << " nx ny nz [steps]" << endl;
		return(1);
	}

	const int size[3] = { atoi(argv[1]), atoi(argv[2]), atoi(argv[3]) };

	// Optional value
	const int time_steps = (argc > 4 ? atoi(argv[4]) : 10);

	// 2 materiales matriz y fibra (cilindro en z)
	const int micro_type = 3;
	const double d_eps = 0.01;
	const int dir = 2;

	double micro_params[4] = { 1.0, 1.0, 1.0, 0.15 };

	material_t mat_params[2];
	mat_params[0].set(1.0e6, 0.3, 5.0e4, 5.0e4, 1);
	mat_params[1].set(1.0e6, 0.3, 1.0e4, 0.0e-1, 0);

	micropp<3> micro(1, size, micro_type, micro_params, mat_params);

	double eps[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
	double sig[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
	for (int t = 0; t < time_steps; ++t) {
		cout << "time step = " << t << endl;
		if (t < 30)
			eps[dir] += d_eps;
		else if (t < 80)
			eps[dir] -= d_eps;
		else if (t < 130)
			eps[dir] += d_eps;
		else if (t < 250)
			eps[dir] -= d_eps;
		else
			eps[dir] += d_eps;

		micro.set_macro_strain(0, eps);
		micro.homogenize();

		micro.get_macro_stress(0, sig);

		micro.update_vars();

		cout << "eps = " << eps[dir] << endl;
		cout
			<< "sig = "
			<< sig[0] << " " << sig[1] << " " << sig[2] << " "
			<< sig[3] << " " << sig[4] << " " << sig[5]
			<< endl;

		cout << endl;
		micro.output (t, 0);
	}
	return 0;
}
