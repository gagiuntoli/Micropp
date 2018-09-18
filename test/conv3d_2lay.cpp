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

#define D_EPS 1.0e-3

int main (int argc, char *argv[])
{
	const int dim = 3;

	const int nmax = 11;
	const int dir = 3;
	const int time_steps = 80;

	material_t mat_params[2];
	mat_params[0].set(3.0e7, 0.25, 1.0e0, 2.0e5, 1);
	mat_params[1].set(3.0e7, 0.25, 1.0e5, 2.0e5, 0);
	double micro_params[5] = { 1., 1., 1., .5, 0. };
	int micro_type = 1; // 2 capas

	double sig[6];
	double res[nmax - 2][2];

	for (int n = 2; n < nmax; ++n) {

		cout << "Running n = " << n << endl;

		double eps[6] = { 0. };
		int size[3] = { n, n, n };
		micropp<3> micro(1, size, micro_type, micro_params, mat_params);

		for (int t = 0; t < time_steps; ++t) {

			eps[dir] += D_EPS;

			micro.set_macro_strain(0, eps);
			micro.homogenize();

			micro.write_info_files ();

			micro.update_vars();

		}
		micro.get_macro_stress(0, sig);

		res[n - 2][0] = n - 1;
		res[n - 2][1] = sig[dir];
	}

	ofstream file;
	file.open("result.dat");
	for (int n = 2; n < nmax; ++n)
		file << res[n - 2][0] << "\t" << res[n - 2][1] << endl;

	file.close();
	return 0;
}
