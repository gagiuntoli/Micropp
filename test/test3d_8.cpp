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

#include <cstring>
#include <ctime>
#include <cassert>

#include "micro.hpp"

using namespace std;

#define dim 3
#define nmaterials 2

int main(int argc, char **argv)
{

	if (argc < 4) {
		cerr << "Usage: " << argv[0] << " nx ny nz [ngp] [steps]" << endl;
		return(1);
	}

	const int nx = atoi(argv[1]);
	const int ny = atoi(argv[2]);
	const int nz = atoi(argv[3]);
	const int ngp = (argc > 4 ? atoi(argv[4]) : 2);
	const int time_steps = (argc > 5 ? atoi(argv[5]) : 10);  // Optional value

	assert(nx > 1 && ny > 1 && nz > 1 && ngp > 1 && time_steps > 0);

	int size[dim] = {nx, ny, nz};

	int micro_type = 1;	// 2 materiales matriz y fibra (3D esfera en matriz)

	double micro_params[5] = {1.0,		// lx
	                          1.0,		// ly
	                          1.0,		// lz
	                          0.1,		// Layer width
	                          1.0e-5};	// INV_MAX


	int mat_types[nmaterials] = {1, 0};	// dos materiales lineales (type = 0)

	double mat_params[2 * MAX_MAT_PARAM];
	mat_params[0 * MAX_MAT_PARAM + 0] = 1.0e6; // E
	mat_params[0 * MAX_MAT_PARAM + 1] = 0.3;   // nu
	mat_params[0 * MAX_MAT_PARAM + 2] = 5.0e4; // Sy
	mat_params[0 * MAX_MAT_PARAM + 3] = 5.0e4; // Ka

 	mat_params[1 * MAX_MAT_PARAM + 0] = 1.0e6;
	mat_params[1 * MAX_MAT_PARAM + 1] = 0.3;
	mat_params[1 * MAX_MAT_PARAM + 2] = 1.0e4;
	mat_params[1 * MAX_MAT_PARAM + 3] = 0.0e-1;

	int dir = 2;
	double d_eps = 0.01;
	double eps[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
	double sig[6], (*sig_test)[3];

    sig_test = (double (*)[3]) malloc (3 * ngp * sizeof(double));
	
	micropp_t micro(dim, ngp, size, micro_type, micro_params, mat_types, mat_params);

	for (int t = 0; t < time_steps; ++t) {
		cout << "Time step = " << t << endl;
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

		cout << "setting strains ..." << endl;
		for (int gp = 0; gp < ngp; ++gp) {
			micro.set_macro_strain(gp, eps);
			cout << "gp = " << gp << " eps = ";
			cout << scientific;

			for (int i = 0; i < 6; ++i)
				cout << setw(14) << eps[i] << " ";
			cout << endl;
		}

		double MacroStress[6], MacroCtan[36];
		cout << "Homogenizing ..." << endl;
		micro.homogenize();

		cout << "Getting stresses ..." << endl;
		for (int gp = 0; gp < ngp; ++gp) {
			micro.get_macro_stress(gp, sig);
			cout << "gp = " << gp << " sig = ";

			cout << scientific;
			for (int i = 0; i < 6; ++i)
				cout << setw(14) << sig[i] << " ";
			cout << endl;

            memcpy(sig_test[gp], sig, 3*sizeof(double));
		}

		for (int gp = 1; gp < ngp; ++gp) {
			cout << "Diff: \t";
			for (int i = 0; i < 3; ++i) {
				const double tmp = fabs(sig_test[gp][i] - sig_test[0][i]);
				cout << tmp << "\t";
				assert(tmp < 1.0e-6);
			}
			cout << endl;
		}

		micro.update_vars();
		micro.output (t, 1);
		micro.write_info_files();
		cout << endl;
	}
	return 0;
}
