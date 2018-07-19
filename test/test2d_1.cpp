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
	if (argc < 3) {
		cerr << "Usage: " << argv[0] << " nx ny [steps]" << endl;
		return(1);
	}

	const int dim = 2;
	const int nx = atoi(argv[1]);
	const int ny = atoi(argv[2]);
	const int time_steps = (argc > 3 ? atoi(argv[3]) : 10);  // Optional value

	assert(nx > 1 && ny > 1);

	int size[3];
	size[0] = nx;
	size[1] = ny;

	int micro_type = 1; // 2 materiales matriz y fibra (3D esfera en matriz)
	double micro_params[5];
	micro_params[0] = 1.0; // lx
	micro_params[1] = 1.0; // ly
	micro_params[2] = 1.0; // lz
	micro_params[3] = 0.1; // grosor capa de abajo
	micro_params[4] = 0; // INV_MAX

	int mat_types[2]; // dos materiales lineales (type = 0)
	mat_types[0] = 0;
	mat_types[1] = 0;

	double mat_params[2 * MAX_MAT_PARAM];
	mat_params[0 * MAX_MAT_PARAM + 0] = 1.0e6; // E
	mat_params[0 * MAX_MAT_PARAM + 1] = 0.3;   // nu
	mat_params[0 * MAX_MAT_PARAM + 2] = 5.0e4; // Sy
	mat_params[0 * MAX_MAT_PARAM + 3] = 5.0e4; // Ka

 	mat_params[1 * MAX_MAT_PARAM + 0] = 1.0e6;
	mat_params[1 * MAX_MAT_PARAM + 1] = 0.3;
	mat_params[1 * MAX_MAT_PARAM + 2] = 1.0e4;
	mat_params[1 * MAX_MAT_PARAM + 3] = 0.0e-1;

	micropp<2> micro(1, size, micro_type, micro_params, mat_types, mat_params);

	int dir = 2;
	double sig[3], ctan[6];
	double eps[3] = { 0 };
	double d_eps = 0.01;

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
		micro.get_macro_ctan(0, ctan);

		micro.update_vars();
		micro.write_info_files ();

		cout << "eps =\t";
		for (int i = 0; i < 3; ++i)
			cout << setw(14) << eps[i] << "\t";
		cout << endl;

		cout << "sig =\t";
		for (int i = 0; i < 3; ++i)
			cout << setw(14) << sig[i] << "\t";
		cout << endl;

		cout << "ctan =\t";
		for (int i = 0; i < 3; ++i)
			cout << setw(14) << ctan[i] << "\t";
		cout << endl;

		micro.output (t, 0);
	}
	return 0;
}
